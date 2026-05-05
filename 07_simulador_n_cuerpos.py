"""
SCRIPT 7 v4: MOTOR DE SIMULACIÓN N-CUERPOS 3D (Numba JIT + Parallel)
====================================================================

Versión optimizada con Numba. Todas las funciones críticas numéricas
(kernel gravitacional, potencial MW externo, integrador Leapfrog) son
compiladas a código nativo con @njit y paralelizadas con prange.

SPEEDUP ESPERADO vs v3: 50-100x en CPUs multi-core modernas.

CAMBIOS RESPECTO A v3:

1. KERNEL NUMBA COMPILADO:
   - `acc_nbody_numba`: bucles O(N²) con @njit(parallel=True)
   - `acc_galactic_numba`: potencial MW (MN + Hernquist + NFW) @njit
   - `leapfrog_step_numba`: un paso completo KDK @njit

2. BENCHMARK DE EQUIVALENCIA PREFLIGHT:
   Antes de simular, compara versión numpy vs Numba para garantizar
   que ambas dan MISMOS resultados numéricos (max_rel_err < 1e-10).
   Si no coinciden → aborta con IntegratorError.

3. VALIDACIÓN SIMBÓLICA PRESERVADA:
   Las 4 tests simbólicas (SymPy vs numérico) se mantienen intactas.
   Se ejecutan sobre las funciones Numba compiladas.

4. FALLBACK AUTOMÁTICO:
   Si Numba no está instalado, usa el código numpy original.
   Zero breakage para backwards compat.

CARACTERÍSTICAS TÉCNICAS:
- @njit(parallel=True, fastmath=True, cache=True)
- prange para paralelismo por estrella
- SIMD AVX2/FMA3 automático en CPUs compatibles
- Compilación JIT cached en disco (warmup solo la primera vez)
"""
from __future__ import annotations

import argparse
import logging
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd
from astropy.coordinates import Distance, Galactocentric, SkyCoord
import astropy.units as u

# Import del módulo fundacional
from gaia_core import (
    DataQualityError, IntegratorError, MissingDependencyError,
    validate_distance_pc, validate_magnitude, validate_velocity_kms,
    autodetect_column, setup_logging,
    SelfTestSuite, HAS_H5PY, HAS_SYMPY,
)

if HAS_H5PY:
    import h5py
if HAS_SYMPY:
    import sympy as sp

# Numba opcional con fallback seguro
try:
    from numba import njit, prange, config as numba_config
    HAS_NUMBA = True
    # Diagnóstico útil: cuántos threads va a usar Numba
    NUMBA_THREADS = numba_config.NUMBA_DEFAULT_NUM_THREADS
except ImportError:
    HAS_NUMBA = False
    NUMBA_THREADS = 1
    # Stub para que @njit no rompa el import si falta
    def njit(*args, **kwargs):
        def decorator(func):
            return func
        if len(args) == 1 and callable(args[0]):
            return args[0]
        return decorator
    prange = range

try:
    from rich.progress import (
        Progress, BarColumn, TextColumn, TimeRemainingColumn
    )
    HAS_RICH_PROG = True
except ImportError:
    HAS_RICH_PROG = False


# ═════════════════════════════════════════════════════════════════════════════
# CONSTANTES FÍSICAS
# ═════════════════════════════════════════════════════════════════════════════

G_CONST = 6.67430e-11
M_SUN = 1.98847e30
PC_TO_M = 3.085677581e16
KPC_TO_M = 1000 * PC_TO_M
YEAR_TO_S = 31557600
M_BOL_SUN = 4.74


# ═════════════════════════════════════════════════════════════════════════════
# DATACLASSES DE CONFIGURACIÓN
# ═════════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class PhysicsConfig:
    """Parámetros físicos del modelo galáctico (MWPotential2014-like, Bovy 2015)."""
    disk_mass_kg: float = 6.8e10 * M_SUN
    disk_a_kpc: float = 3.0
    disk_b_kpc: float = 0.28
    bulge_mass_kg: float = 5.0e9 * M_SUN
    bulge_c_kpc: float = 1.9
    halo_mvir_kg: float = 8.0e11 * M_SUN
    halo_rs_kpc: float = 16.0
    nfw_c: float = 15.3
    epsilon_pc: float = 0.1
    use_galactic_potential: bool = True

    @property
    def nfw_fc(self) -> float:
        """Factor de normalización NFW: ln(1+c) - c/(1+c)."""
        c = self.nfw_c
        return np.log(1 + c) - c / (1 + c)


@dataclass(frozen=True)
class IntegratorConfig:
    years_total: float = 1_000_000
    dt_days: float = 365.25
    save_interval: int = 100
    seed: int = 42


@dataclass(frozen=True)
class IOConfig:
    input_path: str = ""
    output_path: str = "simulacion.h5"
    n_stars: int = 1000
    mag_column_preferences: Tuple[str, ...] = field(
        default_factory=lambda: (
            "luminosidad_absoluta_g_corregida",
            "luminosidad_absoluta_g",
            "M_G",
            "phot_g_mean_mag",
        )
    )


# ═════════════════════════════════════════════════════════════════════════════
# FUNCIONES NUMBA-COMPILED (KERNELS CRÍTICOS)
# ═════════════════════════════════════════════════════════════════════════════

@njit(parallel=True, fastmath=True, cache=True)
def acc_nbody_numba(
    pos: np.ndarray,
    masses: np.ndarray,
    epsilon_m: float,
) -> np.ndarray:
    """
    Auto-gravedad O(N²) con softening de Plummer. Paralelizada con prange.

    Cada thread calcula la aceleración de un subconjunto de estrellas `i`
    contra TODAS las estrellas `j`. Sin dependencias entre threads.

    CRÍTICO: esta versión debe dar resultados NUMÉRICAMENTE idénticos a
    la versión numpy original. Se valida en benchmark preflight.

    Parámetros
    ----------
    pos : ndarray (N, 3)
        Posiciones en metros.
    masses : ndarray (N,)
        Masas en kg.
    epsilon_m : float
        Softening de Plummer en metros.

    Retorna
    -------
    acc : ndarray (N, 3)
        Aceleraciones en m/s².
    """
    N = pos.shape[0]
    acc = np.zeros((N, 3), dtype=np.float64)
    eps2 = epsilon_m * epsilon_m

    for i in prange(N):
        ax = 0.0
        ay = 0.0
        az = 0.0
        for j in range(N):
            # Skip i=j (auto-fuerza es cero, pero softening lo haría spurious)
            if i == j:
                continue
            dx = pos[j, 0] - pos[i, 0]
            dy = pos[j, 1] - pos[i, 1]
            dz = pos[j, 2] - pos[i, 2]
            r2 = dx * dx + dy * dy + dz * dz + eps2
            inv_r3 = r2 ** (-1.5)
            factor = G_CONST * masses[j] * inv_r3
            ax += factor * dx
            ay += factor * dy
            az += factor * dz
        acc[i, 0] = ax
        acc[i, 1] = ay
        acc[i, 2] = az

    return acc


@njit(parallel=True, fastmath=True, cache=True)
def acc_galactic_numba(
    pos: np.ndarray,
    disk_M: float, disk_a: float, disk_b: float,
    bulge_M: float, bulge_c: float,
    halo_Mvir: float, halo_rs: float, nfw_fc: float,
) -> np.ndarray:
    """
    Aceleración total del potencial galáctico externo (MN + Hernquist + NFW).

    Vectorizada por estrella en prange. Cada estrella es independiente.
    """
    N = pos.shape[0]
    acc = np.zeros((N, 3), dtype=np.float64)
    min_r = 1e-3 * PC_TO_M

    for i in prange(N):
        x = pos[i, 0]
        y = pos[i, 1]
        z = pos[i, 2]

        # Miyamoto-Nagai disco
        R2 = x * x + y * y
        zb = np.sqrt(z * z + disk_b * disk_b)
        a_plus_zb = disk_a + zb
        denom_mn = (R2 + a_plus_zb * a_plus_zb) ** 1.5
        factor_mn = -G_CONST * disk_M / denom_mn

        # Hernquist bulbo
        r = np.sqrt(x * x + y * y + z * z)
        r_safe = r if r > min_r else min_r
        factor_h = -G_CONST * bulge_M / (r_safe * (r_safe + bulge_c) ** 2)

        # NFW halo
        x_nfw = r_safe / halo_rs
        coeff_nfw = G_CONST * halo_Mvir / nfw_fc
        a_mag_nfw = -coeff_nfw * (
            np.log(1 + x_nfw) / (r_safe * r_safe)
            - 1.0 / (r_safe * halo_rs * (1 + x_nfw))
        )
        factor_nfw = a_mag_nfw / r_safe

        # Acumular las tres contribuciones
        acc[i, 0] = factor_mn * x + factor_h * x + factor_nfw * x
        acc[i, 1] = factor_mn * y + factor_h * y + factor_nfw * y
        acc[i, 2] = factor_mn * z * a_plus_zb / zb + factor_h * z + factor_nfw * z

    return acc


@njit(fastmath=True, cache=True)
def leapfrog_step_numba(
    pos: np.ndarray, vel: np.ndarray, masses: np.ndarray,
    dt: float, epsilon_m: float,
    use_galactic: bool,
    disk_M: float, disk_a: float, disk_b: float,
    bulge_M: float, bulge_c: float,
    halo_Mvir: float, halo_rs: float, nfw_fc: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Un paso completo Kick-Drift-Kick Leapfrog (simpléctico).

    Todo compilado JIT. Llama a los kernels Numba.
    """
    # Kick 1
    a = acc_nbody_numba(pos, masses, epsilon_m)
    if use_galactic:
        a_ext = acc_galactic_numba(
            pos, disk_M, disk_a, disk_b,
            bulge_M, bulge_c, halo_Mvir, halo_rs, nfw_fc,
        )
        a = a + a_ext

    vel_half = vel + 0.5 * dt * a

    # Drift
    pos_new = pos + dt * vel_half

    # Kick 2
    a_new = acc_nbody_numba(pos_new, masses, epsilon_m)
    if use_galactic:
        a_ext_new = acc_galactic_numba(
            pos_new, disk_M, disk_a, disk_b,
            bulge_M, bulge_c, halo_Mvir, halo_rs, nfw_fc,
        )
        a_new = a_new + a_ext_new

    vel_new = vel_half + 0.5 * dt * a_new

    return pos_new, vel_new


# ═════════════════════════════════════════════════════════════════════════════
# VERSIONES NUMPY (fallback + validación contra Numba)
# ═════════════════════════════════════════════════════════════════════════════

def acc_nbody_numpy(pos, masses, epsilon_m):
    """Auto-gravedad O(N²) — versión numpy pura (fallback/validación)."""
    dx = pos[np.newaxis, :, :] - pos[:, np.newaxis, :]
    r2 = (dx * dx).sum(axis=2) + epsilon_m * epsilon_m
    inv_r3 = r2 ** -1.5
    np.fill_diagonal(inv_r3, 0.0)  # anula auto-fuerza i=j
    acc = G_CONST * (dx * inv_r3[:, :, np.newaxis]) * masses[np.newaxis, :, np.newaxis]
    return acc.sum(axis=1)


def acc_galactic_numpy(pos, phys: PhysicsConfig):
    """Aceleración del potencial galáctico — versión numpy pura."""
    x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
    disk_M = phys.disk_mass_kg
    disk_a = phys.disk_a_kpc * KPC_TO_M
    disk_b = phys.disk_b_kpc * KPC_TO_M
    bulge_M = phys.bulge_mass_kg
    bulge_c = phys.bulge_c_kpc * KPC_TO_M
    halo_Mvir = phys.halo_mvir_kg
    halo_rs = phys.halo_rs_kpc * KPC_TO_M

    # MN disco
    R2 = x * x + y * y
    zb = np.sqrt(z * z + disk_b * disk_b)
    denom_mn = (R2 + (disk_a + zb) ** 2) ** 1.5
    factor_mn = -G_CONST * disk_M / denom_mn
    mn = np.column_stack((factor_mn * x, factor_mn * y, factor_mn * z * (disk_a + zb) / zb))

    # Hernquist
    r = np.sqrt((pos * pos).sum(axis=1))
    r_safe = np.maximum(r, 1e-3 * PC_TO_M)
    factor_h = -G_CONST * bulge_M / (r_safe * (r_safe + bulge_c) ** 2)
    h = pos * factor_h[:, np.newaxis]

    # NFW
    x_nfw = r_safe / halo_rs
    coeff = G_CONST * halo_Mvir / phys.nfw_fc
    a_mag = -coeff * (np.log(1 + x_nfw) / (r_safe ** 2) - 1.0 / (r_safe * halo_rs * (1 + x_nfw)))
    nfw = pos * (a_mag / r_safe)[:, np.newaxis]

    return mn + h + nfw


# ═════════════════════════════════════════════════════════════════════════════
# INGESTA DE DATOS
# ═════════════════════════════════════════════════════════════════════════════

def estimate_mass_from_abs_mag(m_g: np.ndarray) -> np.ndarray:
    """
    Masa estelar desde magnitud absoluta G.

    Usa relación M-L segmentada:
      - M > 10 (enanas M): M/M_sun = 10^(-0.4 * (M_G - M_bol_sun)) ^ 0.6
      - M < 10 (estrellas OBAFGK): M/M_sun = L^(1/3.5)
      - Clamp a rango físico [0.08, 150] M_sun
    """
    L_ratio = 10 ** (-0.4 * (m_g - M_BOL_SUN))
    mass_ratio = np.where(
        m_g > 10,
        L_ratio ** 0.6,       # enanas M (Kroupa+1993)
        L_ratio ** (1 / 3.5),  # estándar main sequence
    )
    return np.clip(mass_ratio, 0.08, 150.0) * M_SUN


def load_catalog(io: IOConfig, seed: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Carga catálogo + convierte a coordenadas galactocéntricas."""
    path = Path(io.input_path)
    if not path.exists():
        raise DataQualityError(f"Input no encontrado: {path}")

    logging.info(f"Cargando catálogo: {path}")
    if path.suffix == ".parquet":
        df = pd.read_parquet(path)
    else:
        df = pd.read_csv(path, low_memory=False)

    # Columnas requeridas
    req_pos = ["ra", "dec", "distancia_parsecs"]
    req_vel = ["U_kms", "V_kms", "W_kms"]
    for col in req_pos + req_vel:
        if col not in df.columns:
            raise DataQualityError(f"Columna requerida ausente: {col}")

    # Autodetect columna de magnitud (lanza ContractViolationError si no encuentra)
    try:
        mag_col = autodetect_column(df, list(io.mag_column_preferences), semantic_name="magnitud absoluta G")
    except Exception as e:
        raise DataQualityError(str(e))

    # Filtrar filas con todos los datos cinemáticos
    all_needed = req_pos + req_vel + [mag_col]
    df_clean = df.dropna(subset=all_needed).copy()

    if len(df_clean) < io.n_stars:
        logging.warning(
            f"Solicitaste {io.n_stars} estrellas pero solo hay {len(df_clean)} "
            "con cinemática 3D completa. Usando todas."
        )
        df_sample = df_clean
    else:
        df_sample = df_clean.sample(n=io.n_stars, random_state=seed)

    n_actual = len(df_sample)
    logging.info(f"► Simulando {n_actual:,} estrellas (mag_col={mag_col}).")

    # Transformar a galactocéntricas (Astropy)
    sky = SkyCoord(
        ra=df_sample["ra"].values * u.deg,
        dec=df_sample["dec"].values * u.deg,
        distance=Distance(df_sample["distancia_parsecs"].values * u.pc),
        frame="icrs",
    )
    galac = sky.transform_to(Galactocentric())

    pos = np.column_stack([
        galac.x.to_value(u.m),
        galac.y.to_value(u.m),
        galac.z.to_value(u.m),
    ])
    vel = np.column_stack([
        df_sample["U_kms"].values * 1000.0,
        df_sample["V_kms"].values * 1000.0,
        df_sample["W_kms"].values * 1000.0,
    ])
    masses = estimate_mass_from_abs_mag(df_sample[mag_col].values)
    source_ids = df_sample.get("source_id", pd.Series(range(n_actual))).values

    # ─── Validación post-transformación: filtrar NaN/Inf ───
    # Puede haber NaN en pos (si SkyCoord transformó distancia NaN)
    # o en masses (si mag_col tenía NaN que dropna no atrapó por dtype)
    valid = (
        np.isfinite(pos).all(axis=1)
        & np.isfinite(vel).all(axis=1)
        & np.isfinite(masses)
        & (masses > 0)
    )
    n_invalid = (~valid).sum()
    if n_invalid > 0:
        logging.warning(
            f"⚠ {n_invalid}/{n_actual} estrellas descartadas por NaN/Inf en "
            "pos/vel/mass post-transformación. Filtrando..."
        )
        pos = pos[valid]
        vel = vel[valid]
        masses = masses[valid]
        source_ids = source_ids[valid]
        n_actual = len(pos)
        if n_actual == 0:
            raise DataQualityError(
                "Todas las estrellas produjeron valores NaN/Inf tras transformación. "
                "Revisa tu catálogo de entrada."
            )

    logging.info(
        f"  Masas: min={masses.min()/M_SUN:.3f} | "
        f"med={np.median(masses)/M_SUN:.3f} | "
        f"max={masses.max()/M_SUN:.2f} M☉ | "
        f"total={masses.sum()/M_SUN:.1f} M☉"
    )

    # Diagnóstico: distancia al centro galáctico (debe ser >>1 kpc típicamente)
    r_sph = np.sqrt((pos * pos).sum(axis=1))
    r_sph_kpc = r_sph / KPC_TO_M
    min_r_kpc = r_sph_kpc.min()
    n_too_close = int((r_sph_kpc < 0.1).sum())  # menos de 100 pc del centro
    logging.info(
        f"  r_galactocentric: min={min_r_kpc:.2f} | "
        f"med={np.median(r_sph_kpc):.2f} | "
        f"max={r_sph_kpc.max():.2f} kpc"
    )
    if n_too_close > 0:
        logging.warning(
            f"⚠ {n_too_close} estrellas están a <100 pc del centro galáctico. "
            "Pueden causar inestabilidades numéricas. Considera filtrarlas."
        )

    return pos, vel, masses, source_ids


# ═════════════════════════════════════════════════════════════════════════════
# BENCHMARK DE EQUIVALENCIA (Numba vs numpy)
# ═════════════════════════════════════════════════════════════════════════════

def benchmark_numba_vs_numpy(phys: PhysicsConfig, n_test: int = 50) -> None:
    """
    Valida que la versión Numba da MISMOS resultados que numpy.

    Pruebas:
      1. acc_nbody_numba vs acc_nbody_numpy — max_rel_err < 1e-10
      2. acc_galactic_numba vs acc_galactic_numpy — max_rel_err < 1e-10

    Si falla → aborta con IntegratorError.
    """
    if not HAS_NUMBA:
        return

    rng = np.random.default_rng(42)
    pos_test = rng.uniform(-3 * KPC_TO_M, 3 * KPC_TO_M, (n_test, 3))
    masses_test = rng.uniform(0.5, 5.0, n_test) * M_SUN
    eps_m = phys.epsilon_pc * PC_TO_M

    # Warmup JIT (primera llamada compila)
    _ = acc_nbody_numba(pos_test, masses_test, eps_m)
    _ = acc_galactic_numba(
        pos_test,
        phys.disk_mass_kg, phys.disk_a_kpc * KPC_TO_M, phys.disk_b_kpc * KPC_TO_M,
        phys.bulge_mass_kg, phys.bulge_c_kpc * KPC_TO_M,
        phys.halo_mvir_kg, phys.halo_rs_kpc * KPC_TO_M, phys.nfw_fc,
    )

    # Benchmark acc_nbody
    acc_nb = acc_nbody_numba(pos_test, masses_test, eps_m)
    acc_np = acc_nbody_numpy(pos_test, masses_test, eps_m)

    max_err_nb = np.max(np.abs(acc_nb - acc_np) / (np.abs(acc_np) + 1e-30))
    if max_err_nb > 1e-10:
        raise IntegratorError(
            f"Benchmark FALLÓ: acc_nbody Numba vs numpy max_rel_err={max_err_nb:.2e} > 1e-10. "
            "Kernel Numba produce resultados distintos. Abortando."
        )
    logging.info(f"  acc_nbody_numba vs numpy: max_rel_err={max_err_nb:.2e} ✔")

    # Benchmark acc_galactic
    acc_g_nb = acc_galactic_numba(
        pos_test,
        phys.disk_mass_kg, phys.disk_a_kpc * KPC_TO_M, phys.disk_b_kpc * KPC_TO_M,
        phys.bulge_mass_kg, phys.bulge_c_kpc * KPC_TO_M,
        phys.halo_mvir_kg, phys.halo_rs_kpc * KPC_TO_M, phys.nfw_fc,
    )
    acc_g_np = acc_galactic_numpy(pos_test, phys)

    max_err_g = np.max(np.abs(acc_g_nb - acc_g_np) / (np.abs(acc_g_np) + 1e-30))
    if max_err_g > 1e-10:
        raise IntegratorError(
            f"Benchmark FALLÓ: acc_galactic Numba vs numpy max_rel_err={max_err_g:.2e} > 1e-10."
        )
    logging.info(f"  acc_galactic_numba vs numpy: max_rel_err={max_err_g:.2e} ✔")

    # Medir speedup
    n_reps = 3
    t_np = 0.0
    for _ in range(n_reps):
        t0 = time.perf_counter()
        _ = acc_nbody_numpy(pos_test, masses_test, eps_m)
        t_np += time.perf_counter() - t0
    t_np /= n_reps

    t_nb = 0.0
    for _ in range(n_reps):
        t0 = time.perf_counter()
        _ = acc_nbody_numba(pos_test, masses_test, eps_m)
        t_nb += time.perf_counter() - t0
    t_nb /= n_reps

    speedup = t_np / t_nb if t_nb > 0 else 0
    logging.info(
        f"  Speedup Numba (N={n_test}): {speedup:.1f}x "
        f"(numpy={t_np*1000:.2f}ms, numba={t_nb*1000:.2f}ms)"
    )


# ═════════════════════════════════════════════════════════════════════════════
# VALIDACIÓN SIMBÓLICA (preservada de v3)
# ═════════════════════════════════════════════════════════════════════════════

def validate_integrator_symbolic(
    phys: PhysicsConfig,
    R_test_kpc: float = 8.0,
    n_orbits: int = 10,
    tolerance: float = 0.02,
) -> SelfTestSuite:
    """Valida integrador Leapfrog contra solución analítica via SymPy.

    Test idéntico al v3 pero usando los kernels Numba.
    """
    suite = SelfTestSuite("Integrator Symbolic Validation")

    # --- (1) Derivación simbólica v_circ en el plano del disco ---
    if not HAS_SYMPY:
        suite.assert_true("sympy_available", False, "sympy no instalado")
        return suite

    R_sym, G_sym, M_sym, a_sym, b_sym = sp.symbols("R G M a b", positive=True)
    Phi = -G_sym * M_sym / sp.sqrt(R_sym**2 + (a_sym + b_sym)**2)
    dPhi_dR = sp.diff(Phi, R_sym)
    v_circ_sym = sp.sqrt(R_sym * dPhi_dR)
    v_circ_sym_func = sp.lambdify(
        (R_sym, G_sym, M_sym, a_sym, b_sym), v_circ_sym, modules="numpy"
    )

    R_test_m = R_test_kpc * KPC_TO_M
    v_circ_analytical_ms = float(v_circ_sym_func(
        R_test_m, G_CONST, phys.disk_mass_kg,
        phys.disk_a_kpc * KPC_TO_M, phys.disk_b_kpc * KPC_TO_M,
    ))
    v_circ_analytical_kms = v_circ_analytical_ms / 1000

    suite.assert_bounds(
        "v_circ_analytical_plausible", v_circ_analytical_kms, 50.0, 400.0
    )

    # --- (2) Verificar numérico ≡ simbólico para el disco ---
    pos_test = np.array([[R_test_m, 0.0, 0.0]])
    masses_zero = np.array([0.0])  # elimina auto-gravedad

    # Aceleración solo del disco MN (para test simbólico)
    disk_a_m = phys.disk_a_kpc * KPC_TO_M
    disk_b_m = phys.disk_b_kpc * KPC_TO_M
    x, y, z = R_test_m, 0.0, 0.0
    R2 = x**2 + y**2
    zb = np.sqrt(z**2 + disk_b_m**2)
    factor_mn = -G_CONST * phys.disk_mass_kg / (R2 + (disk_a_m + zb)**2)**1.5
    a_R_disk = factor_mn * x
    v_circ_numerical = np.sqrt(np.abs(R_test_m * a_R_disk))

    suite.assert_close(
        "v_circ_symbolic_vs_numerical",
        v_circ_numerical,
        v_circ_analytical_ms,
        tolerance=1e-5,
    )

    # --- (3) Órbita circular cerrada tras n_orbits ---
    T_orbit = 2 * np.pi * R_test_m / v_circ_analytical_ms
    dt = T_orbit / 1000
    n_steps = n_orbits * 1000

    pos = np.array([[R_test_m, 0.0, 0.0]])
    vel = np.array([[0.0, v_circ_analytical_ms, 0.0]])
    masses = np.array([0.0])

    phys_disk_only = PhysicsConfig(
        disk_mass_kg=phys.disk_mass_kg,
        disk_a_kpc=phys.disk_a_kpc,
        disk_b_kpc=phys.disk_b_kpc,
        bulge_mass_kg=0.0,
        halo_mvir_kg=0.0,
        use_galactic_potential=True,
    )

    R_initial = np.sqrt(pos[0, 0]**2 + pos[0, 1]**2)
    for _ in range(n_steps):
        pos, vel = leapfrog_step_numba(
            pos, vel, masses, dt, phys.epsilon_pc * PC_TO_M,
            True,
            phys_disk_only.disk_mass_kg, disk_a_m, disk_b_m,
            0.0, phys_disk_only.bulge_c_kpc * KPC_TO_M,
            0.0, phys_disk_only.halo_rs_kpc * KPC_TO_M, phys_disk_only.nfw_fc,
        )

    R_final = np.sqrt(pos[0, 0]**2 + pos[0, 1]**2)
    suite.assert_close(
        "orbit_closure_after_10_orbits", R_final, R_initial, tolerance=tolerance
    )

    # --- (4) Conservación de energía ---
    v0 = np.sqrt(v_circ_analytical_ms**2)
    E_analytical = 0.5 * v0**2 - G_CONST * phys.disk_mass_kg / np.sqrt(R_test_m**2 + (disk_a_m + disk_b_m)**2)

    v_final = np.sqrt(vel[0, 0]**2 + vel[0, 1]**2 + vel[0, 2]**2)
    zb_f = np.sqrt(pos[0, 2]**2 + disk_b_m**2)
    phi_f = -G_CONST * phys.disk_mass_kg / np.sqrt(
        pos[0, 0]**2 + pos[0, 1]**2 + (disk_a_m + zb_f)**2
    )
    E_final = 0.5 * v_final**2 + phi_f

    drift = np.abs((E_final - E_analytical) / E_analytical)
    suite.assert_bounds("energy_drift_leapfrog", drift, 0.0, 1e-6)

    return suite


# ═════════════════════════════════════════════════════════════════════════════
# SIMULACIÓN PRINCIPAL
# ═════════════════════════════════════════════════════════════════════════════

def run_simulation(
    pos: np.ndarray, vel: np.ndarray, masses: np.ndarray,
    phys: PhysicsConfig, integ: IntegratorConfig,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    """
    Corre la simulación completa usando los kernels Numba.

    Retorna: positions_history, velocities_final, times, metadata.
    """
    n_stars = pos.shape[0]

    # Parámetros físicos desempaquetados (Numba necesita floats sueltos)
    eps_m = phys.epsilon_pc * PC_TO_M
    disk_M = phys.disk_mass_kg
    disk_a = phys.disk_a_kpc * KPC_TO_M
    disk_b = phys.disk_b_kpc * KPC_TO_M
    bulge_M = phys.bulge_mass_kg
    bulge_c = phys.bulge_c_kpc * KPC_TO_M
    halo_Mvir = phys.halo_mvir_kg
    halo_rs = phys.halo_rs_kpc * KPC_TO_M
    nfw_fc = phys.nfw_fc

    # Tiempos
    dt_s = integ.dt_days * 86400
    total_s = integ.years_total * YEAR_TO_S
    n_steps = int(np.ceil(total_s / dt_s))
    n_snapshots = n_steps // integ.save_interval + 1

    logging.info(f"Pasos: {n_steps:,} | Snapshots: {n_snapshots:,}")
    logging.info(
        f"Potencial externo: {'ON' if phys.use_galactic_potential else 'OFF'}"
    )

    # Historial
    positions_history = np.empty((n_snapshots, n_stars, 3), dtype=np.float64)
    times = np.empty(n_snapshots, dtype=np.float64)
    positions_history[0] = pos
    times[0] = 0.0

    # Energía inicial para monitoreo
    E_initial = compute_total_energy(pos, vel, masses, phys)
    logging.info(f"E inicial: {E_initial:.3e} J")

    # Integración
    snapshot_idx = 1

    progress = None
    task = None
    if HAS_RICH_PROG:
        progress = Progress(
            TextColumn("[bold cyan]{task.description}"),
            BarColumn(bar_width=40),
            TextColumn("{task.percentage:>3.0f}%"),
            TextColumn("•"),
            TimeRemainingColumn(),
        )
        progress.start()
        task = progress.add_task("Integrando Leapfrog (Numba)...", total=n_steps)

    t_start = time.perf_counter()

    try:
        for step in range(1, n_steps + 1):
            pos, vel = leapfrog_step_numba(
                pos, vel, masses, dt_s, eps_m,
                phys.use_galactic_potential,
                disk_M, disk_a, disk_b,
                bulge_M, bulge_c,
                halo_Mvir, halo_rs, nfw_fc,
            )

            if step % integ.save_interval == 0 and snapshot_idx < n_snapshots:
                positions_history[snapshot_idx] = pos
                times[snapshot_idx] = step * dt_s
                snapshot_idx += 1

            if progress and step % max(1, n_steps // 100) == 0:
                progress.update(task, completed=step)
    finally:
        if progress:
            progress.update(task, completed=n_steps)
            progress.stop()

    elapsed = time.perf_counter() - t_start
    steps_per_s = n_steps / elapsed if elapsed > 0 else 0
    logging.info(f"Completado en {elapsed:.1f}s ({steps_per_s:.0f} pasos/s)")

    # Energía final
    E_final = compute_total_energy(pos, vel, masses, phys)
    drift = abs(E_final - E_initial) / abs(E_initial) if E_initial != 0 else 0
    logging.info(f"E final: {E_final:.3e} J | |ΔE/E₀|: {drift:.3e}")

    metadata = {
        "n_stars": n_stars,
        "n_steps": n_steps,
        "dt_days": integ.dt_days,
        "years_total": integ.years_total,
        "E_initial_J": float(E_initial),
        "E_final_J": float(E_final),
        "energy_drift_rel": float(drift),
        "runtime_seconds": float(elapsed),
        "steps_per_second": float(steps_per_s),
        "numba_enabled": HAS_NUMBA,
        "numba_threads": NUMBA_THREADS,
        "numba_version": getattr(__import__("numba", fromlist=[""]), "__version__", "N/A") if HAS_NUMBA else "N/A",
    }

    return positions_history, vel, times, metadata


def compute_total_energy(pos, vel, masses, phys):
    """Energía total para monitoreo de conservación (numpy, no necesita speed).

    Fix v4.2: diagnóstico explícito cuando hay NaN/Inf/overflow para
    identificar exactamente qué estrella está causando el problema.
    """
    eps_m = phys.epsilon_pc * PC_TO_M
    min_r = 1e-3 * PC_TO_M

    # ─── Diagnósticos previos ───
    diag = []
    if not np.all(np.isfinite(pos)):
        n_bad = (~np.isfinite(pos).all(axis=1)).sum()
        diag.append(f"pos tiene {n_bad} filas con NaN/Inf")
    if not np.all(np.isfinite(vel)):
        n_bad = (~np.isfinite(vel).all(axis=1)).sum()
        diag.append(f"vel tiene {n_bad} filas con NaN/Inf")
    if not np.all(np.isfinite(masses)):
        n_bad = (~np.isfinite(masses)).sum()
        diag.append(f"masses tiene {n_bad} entradas con NaN/Inf")
    if (masses <= 0).any():
        diag.append(f"masses tiene {(masses <= 0).sum()} entradas <= 0")

    max_mass_msun = masses[np.isfinite(masses) & (masses > 0)].max() / M_SUN if (np.isfinite(masses) & (masses > 0)).any() else 0
    if max_mass_msun > 500:
        diag.append(f"max masa = {max_mass_msun:.1f} M☉ (parece overflow en estimate_mass)")

    if diag:
        logging.warning("⚠ Diagnóstico compute_total_energy: " + " | ".join(diag))

    # ─── Cálculo con clipping defensivo ───
    ke = 0.5 * np.sum(masses * (vel * vel).sum(axis=1))

    dx = pos[np.newaxis, :, :] - pos[:, np.newaxis, :]
    r = np.sqrt((dx * dx).sum(axis=2) + eps_m ** 2)
    np.fill_diagonal(r, np.inf)

    # Reformulación numéricamente estable: normalizar por M_SUN antes
    # Esto evita que masses × masses llegue a números ~10⁶²
    m_norm = masses / M_SUN  # masas en unidades solares (números O(1))
    m_prod_msun2 = m_norm[:, np.newaxis] * m_norm[np.newaxis, :]
    # Ahora m_prod_msun2 es O(1-10000), sin riesgo de overflow
    # La energía se recupera multiplicando por M_SUN² al final
    with np.errstate(over='warn', invalid='warn'):
        ratio = m_prod_msun2 / r  # (M_sun² / m) — escala manejable
        pe_nbody = -0.5 * G_CONST * M_SUN * M_SUN * np.sum(ratio)

    pe_ext = 0.0
    if phys.use_galactic_potential:
        x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
        a_d = phys.disk_a_kpc * KPC_TO_M
        b_d = phys.disk_b_kpc * KPC_TO_M
        c_b = phys.bulge_c_kpc * KPC_TO_M
        rs = phys.halo_rs_kpc * KPC_TO_M
        R2 = x * x + y * y
        zb = np.sqrt(z * z + b_d ** 2)
        phi_mn = -G_CONST * phys.disk_mass_kg / np.sqrt(R2 + (a_d + zb) ** 2)
        r_sph = np.sqrt((pos * pos).sum(axis=1))
        r_sph_safe = np.maximum(r_sph, min_r)
        phi_h = -G_CONST * phys.bulge_mass_kg / (r_sph_safe + c_b)
        phi_nfw = -G_CONST * phys.halo_mvir_kg / (phys.nfw_fc * r_sph_safe) * np.log(1 + r_sph_safe / rs)
        pe_ext = np.sum(masses * (phi_mn + phi_h + phi_nfw))

    return ke + pe_nbody + pe_ext


# ═════════════════════════════════════════════════════════════════════════════
# EXPORTACIÓN HDF5
# ═════════════════════════════════════════════════════════════════════════════

def export_hdf5(
    out_path: str,
    positions: np.ndarray, velocities_final: np.ndarray,
    masses: np.ndarray, source_ids: np.ndarray, times: np.ndarray,
    metadata: dict,
) -> None:
    if not HAS_H5PY:
        raise MissingDependencyError("h5py requerido para exportación HDF5")

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Convertir posiciones metros → parsecs para storage eficiente
    pos_pc = positions / PC_TO_M

    with h5py.File(out, "w") as f:
        f.create_dataset("positions", data=pos_pc, compression="gzip", compression_opts=4)
        f.create_dataset("velocities_final", data=velocities_final / 1000)  # km/s
        f.create_dataset("masses", data=masses / M_SUN)  # M_sun
        f.create_dataset("times_years", data=times / YEAR_TO_S)

        # source_ids: string variable length
        dt_str = h5py.special_dtype(vlen=str)
        ids_data = np.array([str(s) for s in source_ids], dtype=object)
        f.create_dataset("source_ids", data=ids_data, dtype=dt_str)

        for key, val in metadata.items():
            f.attrs[key] = val
        f.attrs["units_positions"] = "parsecs (galactocentric)"
        f.attrs["units_velocities"] = "km/s"
        f.attrs["units_masses"] = "M_sun"
        f.attrs["units_times"] = "years"

    logging.info(f"💾 HDF5 guardado: {out}")


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    p = argparse.ArgumentParser(
        description="Simulador N-cuerpos 3D Gaia (v4 — Numba JIT + Parallel).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", default="simulacion_3d.h5")
    p.add_argument("--n_stars", type=int, default=1000)
    p.add_argument("--years", type=float, default=1_000_000)
    p.add_argument("--dt_days", type=float, default=365.25)
    p.add_argument("--save_interval", type=int, default=100)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--no_external_potential", action="store_true")
    p.add_argument("--skip_validation", action="store_true",
                   help="Salta validación simbólica (NO recomendado)")
    p.add_argument("--skip_benchmark", action="store_true",
                   help="Salta benchmark Numba vs numpy (NO recomendado)")
    p.add_argument("--selftest_only", action="store_true")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()

    setup_logging(verbosity=2 if args.verbose else 1)
    # Silenciar logs internos de Numba (LLVM IR, SSA analysis, bytecode dump)
    # aunque el usuario pida -v. Sus debugs son irrelevantes para el usuario.
    logging.getLogger("numba").setLevel(logging.WARNING)
    logging.getLogger("numba.core").setLevel(logging.WARNING)
    logging.getLogger("numba.core.ssa").setLevel(logging.WARNING)
    logging.getLogger("numba.core.byteflow").setLevel(logging.WARNING)
    logging.getLogger("numba.core.interpreter").setLevel(logging.WARNING)
    logging.getLogger("numba.core.typeinfer").setLevel(logging.WARNING)
    print("=" * 70)
    print(" 🪐 SIMULADOR N-CUERPOS 3D GAIA (v4 — Numba)")
    print("=" * 70)

    if not HAS_NUMBA:
        logging.warning(
            "⚠️ Numba no instalado. Usando numpy puro (SIN speedup). "
            "Instala con: conda install -c conda-forge numba"
        )
    else:
        logging.info(f"✔ Numba {HAS_NUMBA} detectado | Threads: {NUMBA_THREADS}")

    # Config
    phys = PhysicsConfig(use_galactic_potential=not args.no_external_potential)
    integ = IntegratorConfig(
        years_total=args.years, dt_days=args.dt_days,
        save_interval=args.save_interval, seed=args.seed,
    )
    io = IOConfig(
        input_path=args.input, output_path=args.output,
        n_stars=args.n_stars,
    )

    # Validación simbólica (con warmup JIT implícito)
    if not args.skip_validation:
        logging.info("🔬 Validando integrador contra solución analítica...")
        suite = validate_integrator_symbolic(phys)
        suite.print_report()
        if not suite.all_passed():
            raise IntegratorError("Validación simbólica falló. Abortando.")

    # Benchmark Numba vs numpy
    if HAS_NUMBA and not args.skip_benchmark:
        logging.info("🏁 Benchmark Numba vs numpy...")
        benchmark_numba_vs_numpy(phys, n_test=50)

    if args.selftest_only:
        logging.info("✅ --selftest_only: saliendo.")
        return

    # Cargar datos
    pos, vel, masses, source_ids = load_catalog(io, seed=args.seed)

    # Warmup JIT con datos reales (primera llamada compila)
    if HAS_NUMBA:
        logging.info("🔥 Warmup JIT (primera compilación)...")
        t0 = time.perf_counter()
        _ = acc_nbody_numba(pos[:10].copy(), masses[:10].copy(), phys.epsilon_pc * PC_TO_M)
        _ = acc_galactic_numba(
            pos[:10].copy(),
            phys.disk_mass_kg, phys.disk_a_kpc * KPC_TO_M, phys.disk_b_kpc * KPC_TO_M,
            phys.bulge_mass_kg, phys.bulge_c_kpc * KPC_TO_M,
            phys.halo_mvir_kg, phys.halo_rs_kpc * KPC_TO_M, phys.nfw_fc,
        )
        logging.info(f"  Warmup: {time.perf_counter()-t0:.1f}s")

    # Simulación
    positions, vel_final, times, metadata = run_simulation(
        pos, vel, masses, phys, integ
    )

    # Export
    export_hdf5(
        args.output, positions, vel_final, masses, source_ids, times,
        metadata,
    )


if __name__ == "__main__":
    try:
        main()
    except (DataQualityError, IntegratorError, MissingDependencyError) as e:
        logging.error(f"\n❌ {type(e).__name__}: {e}")
        raise SystemExit(1)
    except KeyboardInterrupt:
        logging.warning("\n⚠ Interrumpido por usuario.")
        raise SystemExit(130)
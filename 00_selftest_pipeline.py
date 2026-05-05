"""
SCRIPT 00: PIPELINE SELFTEST — El Guardián del Proyecto Gaia
============================================================

Corre una batería de asserts físicos end-to-end sobre TODO el pipeline
(catálogo sintético → enriquecimiento → UVW → clasificación → simulación).

Inspirado directamente en `_selftest` de vacuum_energy.py: cada test
afirma una propiedad física que DEBE cumplirse. Si alguno falla, algo se
rompió y no debe producirse resultado científico hasta fix.

USO:
    python 00_selftest_pipeline.py              # corre todo
    python 00_selftest_pipeline.py --fast       # sólo asserts rápidos
    python 00_selftest_pipeline.py --module 07  # sólo tests del script N

Tests incluidos:
  Bloque A — Constantes físicas y unidades
  Bloque B — Estimación de masas estelares (casos canónicos)
  Bloque C — Potencial galáctico (v_circ @ R_sol)
  Bloque D — Integrador Leapfrog (validación simbólica)
  Bloque E — Contrato de columnas entre scripts (sintético)
  Bloque F — Clasificación vectorizada (casos borde)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Asegurar que podemos importar los módulos del proyecto
sys.path.insert(0, str(Path(__file__).parent))

from gaia_core import (
    SelfTestSuite, setup_logging,
    validate_parallax, validate_distance_pc, validate_magnitude,
    validate_ruwe, validate_velocity_kms, validate_columns_contract,
    DataQualityError, ContractViolationError,
    HAS_SYMPY, HAS_H5PY,
)


def _load_module(name: str, filename: str):
    """Carga un módulo desde archivo registrándolo en sys.modules.

    El registro es necesario para que `@dataclass(frozen=True)` funcione
    correctamente en Python 3.12 (dataclasses internamente consulta
    sys.modules[cls.__module__] al resolver tipos).
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        name, Path(__file__).parent / filename
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module  # ← CRÍTICO: registrar antes de exec
    spec.loader.exec_module(module)
    return module


# ═════════════════════════════════════════════════════════════════════════════
# BLOQUE A — CONSTANTES FÍSICAS Y UNIDADES
# ═════════════════════════════════════════════════════════════════════════════

def test_physical_constants() -> SelfTestSuite:
    """Verifica que las constantes físicas y conversiones sean correctas."""
    s07 = _load_module("s07", "07_simulador_n_cuerpos.py")

    suite = SelfTestSuite("A — Constantes y Unidades")

    # 1 pc ≈ 3.086 × 10^16 m (valor IAU 2015)
    suite.assert_close(
        "PC_TO_M",
        s07.PC_TO_M, 3.0856775814913673e16,
        tolerance=1e-6
    )

    # 1 año juliano = 365.25 × 86400 s
    suite.assert_close(
        "YEAR_TO_S",
        s07.YEAR_TO_S, 365.25 * 86400,
        tolerance=1e-6
    )

    # Masa solar IAU 2015
    suite.assert_close(
        "M_SUN",
        s07.M_SUN, 1.98892e30,
        tolerance=1e-3
    )

    # Magnitud bolométrica solar (IAU 2015: M_bol_sun = 4.74)
    suite.assert_close(
        "M_BOL_SUN",
        s07.M_BOL_SUN, 4.74,
        tolerance=1e-3, relative=False
    )

    return suite


# ═════════════════════════════════════════════════════════════════════════════
# BLOQUE B — ESTIMACIÓN DE MASAS ESTELARES
# ═════════════════════════════════════════════════════════════════════════════

def test_mass_estimation() -> SelfTestSuite:
    """Casos canónicos de la relación Masa-Luminosidad segmentada."""
    s07 = _load_module("s07", "07_simulador_n_cuerpos.py")

    suite = SelfTestSuite("B — Masas Estelares (M-L Segmentada)")

    # Sol: M_G=4.67, T=5772K → debe dar ~1 M_sun
    m_sol = s07.estimate_mass_segmented(
        np.array([4.67]), teff=np.array([5772])
    )[0] / s07.M_SUN
    suite.assert_close("mass_Sun", m_sol, 1.0, tolerance=0.1)

    # Enana M: M_G=11.5, T=3200K → debe dar 0.08-0.3 M_sun
    m_dwarf = s07.estimate_mass_segmented(
        np.array([11.5]), teff=np.array([3200])
    )[0] / s07.M_SUN
    suite.assert_bounds("mass_M_dwarf", m_dwarf, lower=0.05, upper=0.35)

    # Tipo B (M_G=0, T=20000K) → debe dar 2-5 M_sun
    m_B = s07.estimate_mass_segmented(
        np.array([0.0]), teff=np.array([20000])
    )[0] / s07.M_SUN
    suite.assert_bounds("mass_B_type", m_B, lower=1.5, upper=5.0)

    # Sin teff (fallback): Sol debe seguir cercano a 1 M_sun (mayor error)
    m_sol_noteff = s07.estimate_mass_segmented(np.array([4.67]))[0] / s07.M_SUN
    suite.assert_bounds("mass_Sun_no_teff", m_sol_noteff, lower=0.7, upper=1.3)

    # Clamping: M_G muy alto (enana blanca tardía) no debe ir bajo 0.08
    m_low = s07.estimate_mass_segmented(np.array([20.0]))[0] / s07.M_SUN
    suite.assert_close("mass_clamp_low", m_low, 0.08, tolerance=1e-6)

    return suite


# ═════════════════════════════════════════════════════════════════════════════
# BLOQUE C — POTENCIAL GALÁCTICO
# ═════════════════════════════════════════════════════════════════════════════

def test_galactic_potential() -> SelfTestSuite:
    """Velocidad circular en R_sol debe estar en rango observacional."""
    s07 = _load_module("s07", "07_simulador_n_cuerpos.py")

    suite = SelfTestSuite("C — Potencial Galáctico MW")

    phys = s07.PhysicsConfig()

    # v_circ @ 8 kpc debe estar en rango observacional 200-250 km/s
    R_sol = 8.0 * s07.KPC_TO_M
    pos = np.array([[R_sol, 0.0, 0.0]])
    a_vec = s07.acc_galactic_numpy(pos, phys)
    a_rad = -a_vec[0, 0]
    v_circ = np.sqrt(a_rad * R_sol) / 1000  # km/s

    suite.assert_bounds(
        "v_circ_at_R_sol",
        v_circ, lower=200, upper=260
    )

    # Componentes individuales: los 3 deben contribuir aceleración hacia el centro
    phys_disk = s07.PhysicsConfig(
        bulge_mass_kg=0.0, halo_mvir_kg=0.0
    )
    phys_bulge = s07.PhysicsConfig(
        disk_mass_kg=0.0, halo_mvir_kg=0.0
    )
    phys_halo = s07.PhysicsConfig(
        disk_mass_kg=0.0, bulge_mass_kg=0.0
    )

    a_disk = -s07.acc_galactic_numpy(pos, phys_disk)[0, 0]
    a_bulge = -s07.acc_galactic_numpy(pos, phys_bulge)[0, 0]
    a_halo = -s07.acc_galactic_numpy(pos, phys_halo)[0, 0]

    suite.assert_true("disk_attractive_at_Rsol", a_disk > 0,
                       "aceleración hacia el centro")
    suite.assert_true("bulge_attractive_at_Rsol", a_bulge > 0)
    suite.assert_true("halo_attractive_at_Rsol", a_halo > 0)

    # A R_sol, el disco debe dominar (es un hecho observacional)
    suite.assert_true(
        "disk_dominates_at_Rsol",
        a_disk > a_halo,
        f"disk={a_disk:.2e} vs halo={a_halo:.2e}"
    )

    return suite


# ═════════════════════════════════════════════════════════════════════════════
# BLOQUE D — INTEGRADOR LEAPFROG (simbólico)
# ═════════════════════════════════════════════════════════════════════════════

def test_integrator() -> SelfTestSuite:
    """Delega en la validación simbólica interna del Script 07."""
    s07 = _load_module("s07", "07_simulador_n_cuerpos.py")

    phys = s07.PhysicsConfig()
    return s07.validate_integrator_symbolic(phys)


# ═════════════════════════════════════════════════════════════════════════════
# BLOQUE E — CONTRATO DE COLUMNAS ENTRE SCRIPTS (datos sintéticos)
# ═════════════════════════════════════════════════════════════════════════════

def _build_synthetic_dr3_catalog(N: int = 30, seed: int = 42) -> pd.DataFrame:
    """Catálogo Gaia DR3 sintético mínimo para testing del pipeline."""
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        'source_id': range(N),
        'ra': rng.uniform(0, 360, N),
        'dec': rng.uniform(-30, 30, N),
        'parallax': rng.uniform(1, 20, N),
        'parallax_error': rng.uniform(0.01, 0.1, N),
        'pmra': rng.uniform(-50, 50, N),
        'pmdec': rng.uniform(-50, 50, N),
        'ruwe': rng.uniform(0.9, 1.3, N),
        'phot_g_mean_mag': rng.uniform(8, 16, N),
        'phot_bp_mean_mag': rng.uniform(8, 16, N),
        'phot_rp_mean_mag': rng.uniform(8, 16, N),
        'phot_g_mean_flux': rng.uniform(1e5, 1e8, N),
        'phot_g_mean_flux_error': rng.uniform(10, 100, N),
        'bp_rp': rng.uniform(0.2, 2.5, N),
        'phot_bp_rp_excess_factor': rng.uniform(1.1, 1.4, N),
        'astrometric_chi2_al': rng.uniform(10, 100, N),
        'astrometric_n_good_obs_al': rng.integers(50, 200, N),
        'visibility_periods_used': rng.integers(10, 40, N),
        'radial_velocity': np.where(
            rng.random(N) > 0.3,
            rng.uniform(-200, 200, N),
            np.nan
        ),
        'teff_gspphot': rng.uniform(3000, 10000, N),
        'ag_gspphot': rng.uniform(0, 0.5, N),
        'phot_variable_flag': ['NOT_AVAILABLE'] * N,
    })


def test_pipeline_contracts() -> SelfTestSuite:
    """Flujo: datos DR3 sintéticos → Script 03 → asserts columnas → Script 06."""
    suite = SelfTestSuite("E — Contratos de Columnas (S03 → S06)")

    try:
        s03 = _load_module("s03", "03_pipeline_enriquecimiento.py")
    except Exception as e:
        suite.assert_true(
            "script_03_importable", False, f"no se pudo cargar: {e}"
        )
        return suite

    df = _build_synthetic_dr3_catalog(N=30)

    df = s03.compute_distances_and_magnitudes(df, s03.COLUMN_CONFIG, parallax_zp_mas=-0.017)
    df = s03.compute_colors(df, s03.COLUMN_CONFIG)
    df = s03.compute_kinematics(df, s03.COLUMN_CONFIG)
    df = s03.run_classifications(df, s03.COLUMN_CONFIG)

    # Contrato explícito: Script 06 requiere estas columnas
    cols_requeridas_s06 = [
        'clasificacion_espectral', 'luminosidad_absoluta_g',
        'color_bp_rp', 'V_radial', 'fase_evolutiva'
    ]

    try:
        validate_columns_contract(
            df, cols_requeridas_s06,
            producer="Script 03", consumer="Script 06"
        )
        suite.assert_true("contract_S03_to_S06", True)
    except ContractViolationError as e:
        suite.assert_true("contract_S03_to_S06", False, str(e))

    # V_3D debe ser NaN donde V_rad es NaN (NO inventar V=0)
    v3d_nan = df['velocidad_total_3d_kms'].isna().sum()
    vrad_nan = df['V_radial'].isna().sum()
    suite.assert_close(
        "V3D_honest_propagation_of_NaN",
        v3d_nan, vrad_nan,
        tolerance=0, relative=False
    )

    # Distancias deben ser positivas
    suite.assert_true(
        "distances_positive",
        (df['distancia_parsecs'] > 0).all()
    )

    return suite


# ═════════════════════════════════════════════════════════════════════════════
# BLOQUE F — CLASIFICACIÓN (casos borde)
# ═════════════════════════════════════════════════════════════════════════════

def test_classification_edge_cases() -> SelfTestSuite:
    """Casos borde de la clasificación: quasar, enana M, estrella tipo G."""
    suite = SelfTestSuite("F — Clasificación (Casos Borde)")

    try:
        s06 = _load_module("s06", "06_reporte_clasificacion.py")
    except Exception as e:
        suite.assert_true("script_06_importable", False, str(e))
        return suite

    # DataFrame canónico: 1 Sol, 1 enana M, 1 quasar, 1 gigante
    df = pd.DataFrame({
        'clasificacion_espectral': ['G', 'M', pd.NA, 'K'],
        'luminosidad_absoluta_g': [4.67, 11.5, -25.0, 0.0],
        'color_bp_rp': [0.8, 2.0, 0.3, 1.2],
        'V_radial': [10.0, -30.0, 5000.0, 50.0],
        'fase_evolutiva': [
            'Secuencia Principal (Tipo G/K)',
            'Secuencia Principal (Enana Roja/M)',
            'Indeterminado',
            'Gigante (Roja)'
        ],
    })

    tipo = s06.clasificar_objetos_vectorizado(df)

    # Sol → Estrella G (por letra espectral)
    suite.assert_true("sun_is_G_star", tipo.iloc[0] == 'Estrella G',
                     f"got: {tipo.iloc[0]}")

    # Enana M (M_G=11.5) → "Enana M/K tardía" por magnitud absoluta o "Estrella M"
    suite.assert_true(
        "M_dwarf_classified",
        tipo.iloc[1] in ('Enana M/K tardía', 'Estrella M'),
        f"got: {tipo.iloc[1]}"
    )

    # Quasar (M_G=-25) → debe clasificarse como Quasar/AGN
    suite.assert_true(
        "quasar_classified",
        'Quasar' in tipo.iloc[2] or 'AGN' in tipo.iloc[2],
        f"got: {tipo.iloc[2]}"
    )

    # No debe haber 'Desconocido' en ninguno (todos tienen datos)
    n_desconocidos = (tipo == 'Desconocido').sum()
    suite.assert_close(
        "zero_desconocidos",
        n_desconocidos, 0, tolerance=0, relative=False
    )

    return suite


# ═════════════════════════════════════════════════════════════════════════════
# BLOQUE VALIDADORES — Tests de los validators del gaia_core
# ═════════════════════════════════════════════════════════════════════════════

def test_validators() -> SelfTestSuite:
    """Los validators deben atrapar casos corruptos sin falsos positivos."""
    suite = SelfTestSuite("G — Validators de gaia_core")

    # Paralaje válido no debe disparar
    try:
        validate_parallax(np.array([1.0, 5.0, 100.0]))
        suite.assert_true("valid_parallax_ok", True)
    except DataQualityError:
        suite.assert_true("valid_parallax_ok", False, "falsopositivo")

    # Paralaje imposible DEBE disparar
    try:
        validate_parallax(np.array([1e6]))
        suite.assert_true("invalid_parallax_caught", False, "no detectó")
    except DataQualityError:
        suite.assert_true("invalid_parallax_caught", True)

    # Velocidad luz-absurda DEBE disparar
    try:
        validate_velocity_kms(np.array([1e5]))
        suite.assert_true("invalid_velocity_caught", False)
    except DataQualityError:
        suite.assert_true("invalid_velocity_caught", True)

    # Magnitud absoluta de quasar DEBE pasar (M_G=-25 es físico)
    try:
        validate_magnitude(np.array([-25.0]), abs_mag=True)
        suite.assert_true("quasar_magnitude_ok", True)
    except DataQualityError:
        suite.assert_true("quasar_magnitude_ok", False, "falsopositivo")

    # NaN en velocidad no debe disparar (es legítimo)
    try:
        validate_velocity_kms(np.array([10.0, np.nan, 20.0]))
        suite.assert_true("NaN_velocity_ok", True)
    except DataQualityError:
        suite.assert_true("NaN_velocity_ok", False)

    return suite


# ═════════════════════════════════════════════════════════════════════════════
# ORQUESTADOR PRINCIPAL
# ═════════════════════════════════════════════════════════════════════════════

def main():
    p = argparse.ArgumentParser(
        description="Pipeline SelfTest — corre asserts físicos end-to-end."
    )
    p.add_argument(
        "--fast", action="store_true",
        help="Sólo tests rápidos (saltar integrador simbólico)."
    )
    p.add_argument(
        "--module", type=str, default="all",
        choices=["all", "constants", "mass", "potential", "integrator",
                 "contracts", "classification", "validators"]
    )
    p.add_argument("-v", "--verbose", action="count", default=0)
    args = p.parse_args()

    setup_logging(verbosity=args.verbose)

    print("\n" + "═" * 70)
    print(" 🛡️  GAIA PIPELINE SELFTEST — Guardián del Proyecto")
    print("═" * 70)

    all_suites = []
    modules = {
        "constants": test_physical_constants,
        "mass": test_mass_estimation,
        "potential": test_galactic_potential,
        "integrator": test_integrator,
        "contracts": test_pipeline_contracts,
        "classification": test_classification_edge_cases,
        "validators": test_validators,
    }

    # Saltar integrator en modo --fast (es el más costoso)
    if args.fast and "integrator" in modules:
        del modules["integrator"]

    if args.module != "all":
        modules = {args.module: modules[args.module]}

    for name, test_func in modules.items():
        try:
            suite = test_func()
            suite.print_report()
            all_suites.append(suite)
        except Exception as e:
            print(f"\n❌ Suite '{name}' lanzó excepción: {e}")
            import traceback
            traceback.print_exc()

    # Resumen final
    total_passed = sum(s.summary()[0] for s in all_suites)
    total_tests = sum(s.summary()[1] for s in all_suites)
    all_ok = all(s.all_passed() for s in all_suites)

    print("\n" + "═" * 70)
    status_emoji = "✅" if all_ok else "❌"
    print(f" {status_emoji}  RESULTADO GLOBAL: {total_passed}/{total_tests} tests pasaron")
    print("═" * 70 + "\n")

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
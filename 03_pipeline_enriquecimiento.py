"""
SCRIPT 3 v3: ENRIQUECIMIENTO Y FILTRADO FOTOMÉTRICO/ASTROMÉTRICO GAIA DR3

Calcula parámetros astrofísicos y aplica filtros de calidad sobre catálogos
masivos de Gaia, soportando nativamente Parquet como formato I/O.

CAMBIOS RESPECTO A v2:

1. PARALLAX ZERO-POINT (Lindegren+2021):
   Gaia DR3 tiene un sesgo sistemático de ~−0.017 mas en las paralajes.
   Sin corregir → distancias sesgadas sistemáticamente (~10% a 1 kpc).
   Se expone como `--parallax_zp_mas` (default: -0.017). Impacto directo
   en `distancia_parsecs`, `M_G`, y todas las velocidades derivadas.

2. CONSOLA LIMPIA + LOG DETALLADO (patrón Script 02 v3):
   Warnings de Astropy y ruido de librerías → solo al archivo log.
   Consola muestra barra de progreso + resultados clave.

3. UNITS-SAFE (opcional con --units_safe):
   Velocidades vía Astropy u.mas/u.yr → u.km/u.s (matemáticamente
   equivalente al MAS_YR_TO_KM_S pero inmune a errores de unidades).

4. PARALELIZACIÓN OPCIONAL (--n_jobs, --chunk_size):
   Joblib paraleliza la transformación ICRS→galáctico por chunks para
   catálogos grandes (>500k filas). Transparente si joblib no está.

5. DOWNCAST DE FLOATS (--downcast_floats):
   float64→float32 donde sea seguro. Reduce tamaño Parquet ~50%.

6. SUMMARY JSON (--summary_json):
   KPIs del catálogo final (n_rows, cobertura de VR, medianas, etc.)
   útiles para CI, monitoreo, o comparación entre corridas.

7. INTEGRACIÓN CON GAIA_CORE:
   Custom exceptions (DataQualityError), validators físicos, y
   autodetect_column para retrocompat.
"""

import argparse
import json
import logging
import sys
import time
import warnings
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.utils.exceptions import AstropyWarning
import astropy.units as u

# Integración con módulo fundacional
try:
    from gaia_core import (
        DataQualityError, ContractViolationError,
        validate_parallax, validate_distance_pc, validate_magnitude,
        setup_logging as gc_setup_logging,
    )
    HAS_GAIA_CORE = True
except ImportError:
    HAS_GAIA_CORE = False
    # Fallback si gaia_core no está disponible: no romper el pipeline
    class DataQualityError(Exception): pass
    class ContractViolationError(Exception): pass

# Joblib opcional para paralelización
try:
    from joblib import Parallel, delayed
    HAS_JOBLIB = True
except ImportError:
    HAS_JOBLIB = False

try:
    from rich.console import Console
    from rich.progress import (
        Progress, BarColumn, TextColumn,
        TimeElapsedColumn, MofNCompleteColumn,
    )
    HAS_RICH = True
except ImportError:
    HAS_RICH = False


# ═════════════════════════════════════════════════════════════════════════════
# 1. CONSTANTES Y CONFIGURACIÓN
# ═════════════════════════════════════════════════════════════════════════════

MAS_YR_TO_KM_S: float = 4.74047e-3  # 1 mas/yr · 1 pc → km/s (CORREGIDO v3.1)
# Derivación: μ[rad/s] · d[m] / 1000 = V[km/s]
# = μ[mas/yr] · (π/648000000)/31557600 · d[pc] · 3.0857e16 / 1000
# = μ[mas/yr] · d[pc] · 4.74047e-3
# NOTA v3.1: v1/v2 tenían 4.74047 (sin e-3) — error por factor 1000 en V_tan.

# Gaia DR3 parallax zero-point (Lindegren+2021, A&A 649, A4)
# Valor global aproximado. Para precisión máxima, usar corrección dependiente
# de magnitud/color con `zero_point.get_zpt()` del pkg `gaiadr3-zeropoint`.
DEFAULT_PARALLAX_ZP_MAS: float = -0.017

# Mapeo de columnas Gaia DR3 → alias internos del pipeline.
# Si cambia el release (p.ej. DR4), modificar SOLO este dict.
COLUMN_CONFIG: Dict[str, str] = {
    "ra": "ra",
    "dec": "dec",
    "parallax": "parallax",
    "parallax_error": "parallax_error",
    "pmra": "pmra",
    "pmdec": "pmdec",
    "ruwe": "ruwe",
    "g_mag": "phot_g_mean_mag",
    "bp_mag": "phot_bp_mean_mag",
    "rp_mag": "phot_rp_mean_mag",
    "g_flux": "phot_g_mean_flux",
    "g_flux_error": "phot_g_mean_flux_error",
    "bp_flux": "phot_bp_mean_flux",
    "bp_flux_error": "phot_bp_mean_flux_error",
    "rp_flux": "phot_rp_mean_flux",
    "rp_flux_error": "phot_rp_mean_flux_error",
    "bp_rp_color": "bp_rp",
    "bp_rp_excess": "phot_bp_rp_excess_factor",
    "chi2_astro": "astrometric_chi2_al",
    "n_good_obs_astro": "astrometric_n_good_obs_al",
    "vis_periods": "visibility_periods_used",
    "rv": "radial_velocity",
    "teff": "teff_gspphot",
    "extinction_g": "ag_gspphot",
    "var_flag": "phot_variable_flag",
}

# Contrato explícito con el Script 06.
COLS_CONTRATO_S06 = [
    "clasificacion_espectral",
    "luminosidad_absoluta_g",
    "color_bp_rp",
    "V_radial",
    "fase_evolutiva",
]


# ═════════════════════════════════════════════════════════════════════════════
# 2. LOGGING DUAL (consola limpia + archivo detallado)
# ═════════════════════════════════════════════════════════════════════════════

def configure_logging(log_file: str, verbose: bool) -> Optional['Console']:
    """Logging dual: archivo captura todo (incluye warnings), consola mínimal.

    Retorna la instancia Console de Rich (para prints manuales) o None.
    """
    # Silenciar warnings Astropy en consola, redirigir a logging
    warnings.filterwarnings("ignore", category=AstropyWarning)
    logging.captureWarnings(True)

    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    if root.hasHandlers():
        root.handlers.clear()

    # Handler de archivo: nivel DEBUG (todo)
    log_path = Path(log_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(log_path, encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    ))
    root.addHandler(fh)

    # Handler de consola: INFO si verbose, WARNING si no
    console_level = logging.INFO if verbose else logging.WARNING
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(console_level)
    ch.setFormatter(logging.Formatter("%(message)s"))
    root.addHandler(ch)

    # Silenciar astropy en consola (manteniendo archivo intacto)
    for lg_name in ("astropy", "py.warnings"):
        logging.getLogger(lg_name).setLevel(logging.DEBUG)

    return Console() if HAS_RICH else None


# ═════════════════════════════════════════════════════════════════════════════
# 3. MÓDULO DE FÍSICA Y ASTROMETRÍA
# ═════════════════════════════════════════════════════════════════════════════

def compute_distances_and_magnitudes(
    df: pd.DataFrame, cfg: Dict[str, str], parallax_zp_mas: float
) -> pd.DataFrame:
    """Calcula distancia (1/π corregida) y magnitud absoluta G.

    FIX v3: aplica corrección de parallax zero-point (Lindegren+2021).
    π_corregida = π_raw − parallax_zp_mas
    (donde parallax_zp_mas ~ −0.017 mas en DR3)

    IMPORTANTE — Sesgo de Lutz-Kelker: el estimador 1/π es insesgado solo
    cuando S/N_parallax ≳ 10. Para S/N más bajos conviene usar distancias
    bayesianas (Bailer-Jones+2021 via `gaiadr3.geometric_distances`).
    """
    if cfg["parallax"] not in df.columns:
        return df

    parallax_raw = pd.to_numeric(df[cfg["parallax"]], errors="coerce")

    # Aplicar corrección de zero-point
    parallax_corr = parallax_raw - parallax_zp_mas

    df["parallax_corregido_mas"] = parallax_corr
    df["distancia_parsecs"] = np.where(
        parallax_corr > 0, 1000.0 / parallax_corr, np.nan
    )

    if cfg["g_mag"] in df.columns:
        g_mag = pd.to_numeric(df[cfg["g_mag"]], errors="coerce")
        dist = df["distancia_parsecs"]

        mask_valid = dist > 0
        df["luminosidad_absoluta_g_aparente"] = np.nan
        df.loc[mask_valid, "luminosidad_absoluta_g_aparente"] = (
            g_mag[mask_valid] - 5 * np.log10(dist[mask_valid]) + 5
        )

        if cfg["extinction_g"] in df.columns:
            extinction = pd.to_numeric(
                df[cfg["extinction_g"]], errors="coerce"
            ).fillna(0.0)
            df["luminosidad_absoluta_g_corregida"] = (
                df["luminosidad_absoluta_g_aparente"] - extinction
            )
        else:
            df["luminosidad_absoluta_g_corregida"] = df["luminosidad_absoluta_g_aparente"]

        # Contrato con Script 06
        df["luminosidad_absoluta_g"] = df["luminosidad_absoluta_g_corregida"]

    return df


def compute_colors(df: pd.DataFrame, cfg: Dict[str, str]) -> pd.DataFrame:
    """Asegura la existencia del color BP-RP."""
    if cfg["bp_rp_color"] in df.columns:
        df["color_bp_rp"] = pd.to_numeric(df[cfg["bp_rp_color"]], errors="coerce")
    elif cfg["bp_mag"] in df.columns and cfg["rp_mag"] in df.columns:
        df["color_bp_rp"] = (
            pd.to_numeric(df[cfg["bp_mag"]], errors="coerce")
            - pd.to_numeric(df[cfg["rp_mag"]], errors="coerce")
        )
    return df


def compute_kinematics(
    df: pd.DataFrame, cfg: Dict[str, str], units_safe: bool = False
) -> pd.DataFrame:
    """Velocidades tangenciales y 3D.

    Modos:
    - units_safe=False: fórmula compacta V_tan = 4.74047 · μ[mas/yr] · d[pc]
    - units_safe=True:  Astropy units (inmune a errores de unidades)

    Para velocidad radial: si no hay medida, V_3D = NaN (NO se inventa V_rad=0).
    """
    req_cols = [cfg["pmra"], cfg["pmdec"], "distancia_parsecs"]
    if not all(col in df.columns for col in req_cols):
        return df

    pmra = pd.to_numeric(df[cfg["pmra"]], errors="coerce")
    pmdec = pd.to_numeric(df[cfg["pmdec"]], errors="coerce")
    dist = df["distancia_parsecs"]

    if units_safe:
        # Units-safe con Astropy: inmune a errores de conversión
        mask = pmra.notna() & pmdec.notna() & dist.notna() & (dist > 0)
        v_tan_ra = np.full(len(df), np.nan)
        v_tan_dec = np.full(len(df), np.nan)

        if mask.any():
            pmra_q = pmra[mask].values * u.mas / u.yr
            pmdec_q = pmdec[mask].values * u.mas / u.yr
            dist_q = dist[mask].values * u.pc

            # V_tan = μ · d (con dimensionless_angles para convertir mas·pc/yr → km/s)
            v_ra = (pmra_q * dist_q).to(u.km / u.s, u.dimensionless_angles())
            v_dec = (pmdec_q * dist_q).to(u.km / u.s, u.dimensionless_angles())

            v_tan_ra[mask.values] = v_ra.value
            v_tan_dec[mask.values] = v_dec.value

        df["V_tan_ra_kms"] = v_tan_ra
        df["V_tan_dec_kms"] = v_tan_dec
    else:
        # Fórmula clásica (idéntica matemáticamente)
        df["V_tan_ra_kms"] = MAS_YR_TO_KM_S * pmra * dist
        df["V_tan_dec_kms"] = MAS_YR_TO_KM_S * pmdec * dist

    df["velocidad_tangencial_kms"] = np.sqrt(
        df["V_tan_ra_kms"] ** 2 + df["V_tan_dec_kms"] ** 2
    )

    # Velocidad radial (honesta: NaN si no hay medida)
    if cfg["rv"] in df.columns:
        v_rad = pd.to_numeric(df[cfg["rv"]], errors="coerce")
        df["V_radial_kms"] = v_rad
        df["V_radial"] = v_rad  # contrato S06
        df["velocidad_total_3d_kms"] = np.sqrt(
            df["V_tan_ra_kms"] ** 2 + df["V_tan_dec_kms"] ** 2 + v_rad ** 2
        )
    else:
        df["V_radial_kms"] = np.nan
        df["V_radial"] = np.nan
        df["velocidad_total_3d_kms"] = np.nan

    return df


def _coords_galactic_chunk(
    ra_deg: np.ndarray, dec_deg: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Helper pickle-able para paralelización."""
    coords = SkyCoord(ra=ra_deg * u.degree, dec=dec_deg * u.degree, frame='icrs')
    gal = coords.galactic
    return gal.l.degree, gal.b.degree


def compute_galactic_coords(
    df: pd.DataFrame, cfg: Dict[str, str],
    n_jobs: int = 1, chunk_size: int = 500_000
) -> pd.DataFrame:
    """Proyecta RA/DEC → (l, b) galácticas.

    Paralelización por chunks si n_jobs > 1 y joblib está disponible.
    Para catálogos <500k filas la paralelización no aporta (overhead > ganancia).
    """
    if cfg["ra"] not in df.columns or cfg["dec"] not in df.columns:
        return df

    ra = pd.to_numeric(df[cfg["ra"]], errors="coerce")
    dec = pd.to_numeric(df[cfg["dec"]], errors="coerce")
    mask = ra.notna() & dec.notna()

    df["l_galactic"] = np.nan
    df["b_galactic"] = np.nan

    if not mask.any():
        return df

    ra_vals = ra[mask].values
    dec_vals = dec[mask].values

    try:
        if n_jobs > 1 and HAS_JOBLIB and len(ra_vals) > chunk_size:
            # Paralelización por chunks
            slices = [
                slice(i, min(i + chunk_size, len(ra_vals)))
                for i in range(0, len(ra_vals), chunk_size)
            ]
            results = Parallel(n_jobs=n_jobs, prefer="processes")(
                delayed(_coords_galactic_chunk)(ra_vals[s], dec_vals[s])
                for s in slices
            )
            l_vals = np.concatenate([r[0] for r in results])
            b_vals = np.concatenate([r[1] for r in results])
        else:
            # Ruta serial (la común para catálogos <500k)
            l_vals, b_vals = _coords_galactic_chunk(ra_vals, dec_vals)

        df.loc[mask, "l_galactic"] = l_vals
        df.loc[mask, "b_galactic"] = b_vals
    except Exception as e:
        logging.error(f"Fallo calculando coordenadas galácticas: {e}")

    return df


# ═════════════════════════════════════════════════════════════════════════════
# 4. CLASIFICACIÓN VECTORIZADA
# ═════════════════════════════════════════════════════════════════════════════

def classify_spectral_type(teff: pd.Series) -> pd.Series:
    """Clasificación espectral por T_eff (Pecaut & Mamajek 2013)."""
    bins = [-np.inf, 3700, 5200, 6000, 7500, 10000, 30000, np.inf]
    labels = ["M", "K", "G", "F", "A", "B", "O"]
    return pd.cut(teff, bins=bins, labels=labels, right=True).astype("string")


def classify_evolutionary_phase(
    bp_rp: pd.Series, abs_mag_g: pd.Series
) -> pd.Series:
    """Fase evolutiva heurística en el CMD."""
    c = bp_rp.to_numpy()
    M = abs_mag_g.to_numpy()

    conds = [
        np.isnan(c) | np.isnan(M),
        (M > 10) & (c < 1.0),
        (M < 7 * c - 1) & (c > 1.0),
        (M < 7 * c - 1) & (c <= 1.0),
        (c > 1.5),
        (c > 0.8) & (c <= 1.5),
        (c > 0.4) & (c <= 0.8),
    ]
    choices = [
        "Indeterminado",
        "Posible Enana Blanca",
        "Gigante (Roja)",
        "Gigante/Subgigante (Azul/Blanca)",
        "Secuencia Principal (Enana Roja/M)",
        "Secuencia Principal (Tipo G/K)",
        "Secuencia Principal (Tipo F/A)",
    ]
    return pd.Series(
        np.select(conds, choices, default="Secuencia Principal (Tipo O/B)"),
        index=bp_rp.index, dtype="string"
    )


def run_classifications(df: pd.DataFrame, cfg: Dict[str, str]) -> pd.DataFrame:
    """Aplica todas las clasificaciones vectorizadas."""
    if cfg["teff"] in df.columns:
        teff = pd.to_numeric(df[cfg["teff"]], errors="coerce")
        df["clasificacion_espectral"] = classify_spectral_type(teff)

    if "color_bp_rp" in df.columns and "luminosidad_absoluta_g_corregida" in df.columns:
        df["fase_evolutiva_aprox"] = classify_evolutionary_phase(
            df["color_bp_rp"], df["luminosidad_absoluta_g_corregida"]
        )
        df["fase_evolutiva"] = df["fase_evolutiva_aprox"]  # contrato S06

    if cfg["var_flag"] in df.columns:
        df["indice_variabilidad"] = np.where(
            df[cfg["var_flag"]] == "VARIABLE",
            "Variable (Gaia)", "No Variable (Gaia)"
        )
    elif cfg["ruwe"] in df.columns:
        ruwe = pd.to_numeric(df[cfg["ruwe"]], errors="coerce")
        conds = [ruwe.isna(), ruwe > 1.4, ruwe > 1.0]
        choices = [np.nan, "Alta (RUWE > 1.4)", "Moderada (1.0 < RUWE <= 1.4)"]
        df["indice_variabilidad"] = np.select(
            conds, choices, default="Baja (RUWE <= 1.0)"
        )

    return df


# ═════════════════════════════════════════════════════════════════════════════
# 5. FILTROS (con logging mejorado de % de descarte)
# ═════════════════════════════════════════════════════════════════════════════

FilterType = Tuple[str, List[str], Callable[[pd.DataFrame], pd.Series]]


def get_quality_filters(
    args: argparse.Namespace, cfg: Dict[str, str]
) -> List[FilterType]:
    """Filtros de calidad estándar Gaia DR3."""
    return [
        (f"RUWE < {args.max_ruwe}",
         [cfg["ruwe"]],
         lambda d: pd.to_numeric(d[cfg["ruwe"]], errors="coerce") < args.max_ruwe),

        (f"Parallax S/N > {args.parallax_snr_min}",
         [cfg["parallax"], cfg["parallax_error"]],
         lambda d: (pd.to_numeric(d[cfg["parallax_error"]], errors="coerce") > 0) &
                   ((pd.to_numeric(d[cfg["parallax"]], errors="coerce") /
                     pd.to_numeric(d[cfg["parallax_error"]], errors="coerce")) > args.parallax_snr_min)),

        (f"Phot G S/N > {args.phot_g_snr_min}",
         [cfg["g_flux"], cfg["g_flux_error"]],
         lambda d: (pd.to_numeric(d[cfg["g_flux_error"]], errors="coerce") > 0) &
                   ((pd.to_numeric(d[cfg["g_flux"]], errors="coerce") /
                     pd.to_numeric(d[cfg["g_flux_error"]], errors="coerce")) > args.phot_g_snr_min)),

        ("BP/RP Excess Factor (Evans+2018)",
         [cfg["bp_rp_excess"], "color_bp_rp"],
         lambda d: pd.Series(True, index=d.index) if not args.use_excess_factor_filter else (
             (d["color_bp_rp"].notna()) &
             (pd.to_numeric(d[cfg["bp_rp_excess"]], errors="coerce") > 1.0 + 0.015 * d["color_bp_rp"] ** 2) &
             (pd.to_numeric(d[cfg["bp_rp_excess"]], errors="coerce") < 1.3 + 0.060 * d["color_bp_rp"] ** 2)
         )),

        ("Astrometric Chi2/DOF (Lindegren+2018)",
         [cfg["chi2_astro"], cfg["n_good_obs_astro"], cfg["g_mag"]],
         lambda d: pd.Series(True, index=d.index) if not args.use_astro_chi2_filter else (
             (pd.to_numeric(d[cfg["n_good_obs_astro"]], errors="coerce") > 5) &
             ((pd.to_numeric(d[cfg["chi2_astro"]], errors="coerce") /
               (pd.to_numeric(d[cfg["n_good_obs_astro"]], errors="coerce") - 5)) <
              1.44 * np.maximum(1.0, np.exp(-0.4 * (pd.to_numeric(d[cfg["g_mag"]], errors="coerce") - 19.5))))
         )),
    ]


def apply_filters(
    df: pd.DataFrame, filters: List[FilterType], step_name: str,
    console: Optional['Console'] = None
) -> pd.DataFrame:
    """Aplica filtros con reporte de % de descarte."""
    n_initial = len(df)
    logging.info(f"\n--- {step_name} ---")
    mask_global = pd.Series(True, index=df.index)

    for name, req_cols, condition_func in filters:
        missing_cols = [c for c in req_cols if c not in df.columns]
        if missing_cols:
            msg = f"Omitiendo filtro '{name}' (faltan: {missing_cols})"
            logging.warning(msg)
            continue

        mask_local = condition_func(df).reindex(df.index, fill_value=False)
        discarded = (~mask_local).sum()
        pct = 100 * discarded / len(df) if len(df) > 0 else 0
        logging.info(f"  {name}: descartó {discarded:,} ({pct:.1f}%)")
        mask_global &= mask_local

    df_filtered = df[mask_global].copy()
    surv_pct = 100 * len(df_filtered) / n_initial if n_initial > 0 else 0
    logging.info(
        f"► Tras {step_name}: {len(df_filtered):,}/{n_initial:,} "
        f"({surv_pct:.1f}% sobreviven)\n"
    )

    return df_filtered


# ═════════════════════════════════════════════════════════════════════════════
# 6. UTILIDADES: DOWNCAST, SUMMARY, CONTRATO
# ═════════════════════════════════════════════════════════════════════════════

def downcast_numeric(df: pd.DataFrame) -> pd.DataFrame:
    """Reduce precisión numérica cuando es seguro. float64→float32, int64→int32."""
    for c in df.select_dtypes(include=[np.number]).columns:
        dtype_str = str(df[c].dtype)
        if dtype_str == "float64":
            df[c] = pd.to_numeric(df[c], downcast="float")
        elif dtype_str == "int64":
            df[c] = pd.to_numeric(df[c], downcast="integer")
    return df


def build_summary(df: pd.DataFrame, args: argparse.Namespace) -> Dict:
    """KPIs del catálogo final."""
    has_v3d = "velocidad_total_3d_kms" in df.columns
    has_dist = "distancia_parsecs" in df.columns
    has_mg = "luminosidad_absoluta_g" in df.columns

    summary = {
        "n_rows": int(len(df)),
        "n_cols": int(len(df.columns)),
        "parallax_zp_mas_applied": args.parallax_zp_mas,
        "filters": {
            "max_ruwe": args.max_ruwe,
            "parallax_snr_min": args.parallax_snr_min,
            "phot_g_snr_min": args.phot_g_snr_min,
            "excess_factor_filter": args.use_excess_factor_filter,
            "astro_chi2_filter": args.use_astro_chi2_filter,
        },
        "coverage": {
            "frac_with_vrad": float(df["V_radial"].notna().mean()) if "V_radial" in df.columns else None,
            "frac_with_v3d": float(df["velocidad_total_3d_kms"].notna().mean()) if has_v3d else None,
            "frac_with_distance": float(df["distancia_parsecs"].notna().mean()) if has_dist else None,
            "frac_with_M_G": float(df["luminosidad_absoluta_g"].notna().mean()) if has_mg else None,
        },
        "distributions": {
            "distance_pc_median": float(df["distancia_parsecs"].median()) if has_dist else None,
            "distance_pc_p99": float(df["distancia_parsecs"].quantile(0.99)) if has_dist else None,
            "M_G_median": float(df["luminosidad_absoluta_g"].median()) if has_mg else None,
            "V_3D_kms_median": float(df["velocidad_total_3d_kms"].median()) if has_v3d else None,
        }
    }
    return summary


def self_check_contrato(df: pd.DataFrame) -> bool:
    """Verifica que existan las columnas requeridas por el Script 06."""
    faltantes = [c for c in COLS_CONTRATO_S06 if c not in df.columns]
    if faltantes:
        logging.warning(
            f"⚠️ CONTRATO S06 INCOMPLETO. Faltan columnas: {faltantes}. "
            "Script 06 clasificará estos objetos como 'Desconocido'."
        )
        return False
    logging.info("✔️ Contrato con Script 06 OK: todas las columnas presentes.")
    return True


# ═════════════════════════════════════════════════════════════════════════════
# 7. ORQUESTADOR
# ═════════════════════════════════════════════════════════════════════════════

def run_pipeline(
    args: argparse.Namespace, cfg: Dict[str, str],
    console: Optional['Console'] = None
) -> Optional[pd.DataFrame]:
    """Flujo principal del pipeline."""
    steps = [
        "Leyendo catálogo",
        "Distancias con corrección ZP",
        "Filtros de calidad",
        "Cinemática y coordenadas galácticas",
        "Clasificaciones",
        "Filtros físicos finales",
        "Guardando catálogo",
    ]

    # Warnings preliminares
    if args.parallax_snr_min < 5:
        logging.warning(
            f"parallax_snr_min={args.parallax_snr_min} < 5: sesgo de Lutz-Kelker "
            "significativo. Considera distancias bayesianas (Bailer-Jones+2021)."
        )

    if console:
        console.print(
            f"[bold cyan]🌟 Pipeline Gaia DR3 — Enriquecimiento v3[/bold cyan]\n"
            f"  Input:          [dim]{args.input}[/dim]\n"
            f"  Output:         [dim]{args.output}[/dim]\n"
            f"  Parallax ZP:    [yellow]{args.parallax_zp_mas} mas[/yellow] "
            f"(Lindegren+2021)\n"
            f"  Units-safe:     "
            f"[{'green' if args.units_safe else 'dim'}]{args.units_safe}[/]\n"
            f"  Paralelismo:    [yellow]n_jobs={args.n_jobs}[/yellow]\n"
            f"  Log detallado:  [dim]{args.log_file}[/dim]\n"
        )

    t0 = time.time()

    bar_cols = [
        TextColumn("[bold cyan]{task.description}"),
        BarColumn(bar_width=30),
        MofNCompleteColumn(),
        TextColumn("•"),
        TimeElapsedColumn(),
    ]

    if HAS_RICH and console:
        progress_ctx = Progress(*bar_cols, console=console)
    else:
        progress_ctx = Progress(*bar_cols) if HAS_RICH else None

    if progress_ctx:
        progress = progress_ctx
        progress.start()
        task = progress.add_task("Procesando", total=len(steps))
    else:
        progress = None
        task = None

    def _advance(desc: str):
        if progress:
            progress.update(task, description=desc, advance=1)

    try:
        # ─── Paso 1: Lectura ───
        if progress: progress.update(task, description=f"Paso 1: {steps[0]}")
        try:
            file_ext = Path(args.input).suffix.lower()
            if file_ext == '.parquet':
                df = pd.read_parquet(args.input)
            else:
                df = pd.read_csv(args.input, low_memory=False)
            logging.info(f"Cargado: {len(df):,} filas × {len(df.columns)} cols")
            _advance(f"Paso 2: {steps[1]}")
        except Exception as e:
            raise DataQualityError(f"Error leyendo archivo: {e}")

        # ─── Paso 2: Distancias + magnitudes con ZP correction ───
        df = compute_distances_and_magnitudes(df, cfg, args.parallax_zp_mas)
        df = compute_colors(df, cfg)
        _advance(f"Paso 3: {steps[2]}")

        # ─── Paso 3: Filtros de calidad ───
        quality_filters = get_quality_filters(args, cfg)
        df = apply_filters(df, quality_filters, "Filtros de Calidad", console)
        if df.empty:
            raise DataQualityError("Catálogo vacío tras filtros de calidad")
        _advance(f"Paso 4: {steps[3]}")

        # ─── Paso 4: Cinemática + coordenadas ───
        df = compute_kinematics(df, cfg, units_safe=args.units_safe)
        df = compute_galactic_coords(df, cfg, args.n_jobs, args.chunk_size)
        _advance(f"Paso 5: {steps[4]}")

        # ─── Paso 5: Clasificaciones ───
        df = run_classifications(df, cfg)
        _advance(f"Paso 6: {steps[5]}")

        # ─── Paso 6: Filtros físicos finales ───
        selection_filters = []
        if args.max_dist_pc is not None:
            selection_filters.append((
                f"Distancia ≤ {args.max_dist_pc} pc",
                ["distancia_parsecs"],
                lambda d: d["distancia_parsecs"] <= args.max_dist_pc
            ))
        if args.min_abs_b_gal is not None:
            selection_filters.append((
                f"|b_gal| ≥ {args.min_abs_b_gal}°",
                ["b_galactic"],
                lambda d: d["b_galactic"].abs() >= args.min_abs_b_gal
            ))

        if selection_filters:
            df = apply_filters(df, selection_filters, "Filtros Físicos", console)
            if df.empty:
                raise DataQualityError("Catálogo vacío tras filtros físicos")
        _advance(f"Paso 7: {steps[6]}")

        # ─── Self-check del contrato ───
        contrato_ok = self_check_contrato(df)

        # ─── Downcast opcional ───
        if args.downcast_floats:
            size_before = df.memory_usage(deep=True).sum() / 1e6
            df = downcast_numeric(df)
            size_after = df.memory_usage(deep=True).sum() / 1e6
            logging.info(
                f"Downcast: {size_before:.1f} MB → {size_after:.1f} MB "
                f"({100*(1-size_after/size_before):.1f}% reducción)"
            )

        # ─── Paso 7: Guardado ───
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)

        if args.format == 'parquet':
            df.to_parquet(
                out_path, engine='pyarrow', compression='snappy', index=False
            )
        else:
            df.to_csv(out_path, index=False, float_format='%.6f')

        logging.info(f"Guardado en: {out_path} ({args.format.upper()})")
        _advance("Completado")

        # ─── Summary JSON (opcional) ───
        if args.summary_json:
            summary = build_summary(df, args)
            summary_path = Path(args.summary_json)
            summary_path.parent.mkdir(parents=True, exist_ok=True)
            with open(summary_path, "w", encoding="utf-8") as f:
                json.dump(summary, f, indent=2, ensure_ascii=False)
            logging.info(f"Summary JSON: {summary_path}")

        # ─── Reporte final en consola ───
        elapsed = time.time() - t0
        if console:
            size_mb = out_path.stat().st_size / 1e6
            vr_cov = df["V_radial"].notna().mean() * 100 if "V_radial" in df.columns else 0
            console.print(
                f"\n[bold green]✅ Proceso completado[/bold green]\n"
                f"  Estrellas finales:   [yellow]{len(df):,}[/yellow]\n"
                f"  Cobertura V_radial:  [yellow]{vr_cov:.1f}%[/yellow]\n"
                f"  Tamaño archivo:      [dim]{size_mb:.1f} MB[/dim]\n"
                f"  Tiempo total:        [dim]{elapsed:.1f}s[/dim]\n"
                f"  Contrato S06:        "
                f"{'[green]✔ OK[/green]' if contrato_ok else '[red]✘ FALLÓ[/red]'}"
            )

    finally:
        if progress:
            progress.stop()

    return df


# ═════════════════════════════════════════════════════════════════════════════
# 8. CLI
# ═════════════════════════════════════════════════════════════════════════════

def get_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Pipeline Enriquecimiento Gaia DR3 (v3 — ZP + units-safe).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    g_io = p.add_argument_group("I/O")
    g_io.add_argument("--input", required=True, help="Catálogo DR3 (CSV/Parquet)")
    g_io.add_argument("--output", required=True, help="Destino del catálogo")
    g_io.add_argument("--format", choices=['csv', 'parquet'], default='parquet')
    g_io.add_argument("--log_file", default="03_enriquecimiento.log")
    g_io.add_argument("--summary_json", default=None,
                      help="Ruta opcional para KPIs en JSON")
    g_io.add_argument("-v", "--verbose", action="store_true")

    g_phys = p.add_argument_group("Física")
    g_phys.add_argument(
        "--parallax_zp_mas", type=float, default=DEFAULT_PARALLAX_ZP_MAS,
        help="Parallax zero-point (Lindegren+2021). Default: -0.017 mas"
    )
    g_phys.add_argument(
        "--units_safe", action="store_true",
        help="Usar Astropy units para velocidades (inmune a errores de unidades)"
    )

    g_perf = p.add_argument_group("Performance")
    g_perf.add_argument("--n_jobs", type=int, default=1,
                        help="Paralelismo para coord. galácticas (>1 requiere joblib)")
    g_perf.add_argument("--chunk_size", type=int, default=500000,
                        help="Chunks para paralelización")
    g_perf.add_argument("--downcast_floats", action="store_true",
                        help="float64→float32 para reducir tamaño")

    g_q = p.add_argument_group("Filtros de Calidad")
    g_q.add_argument("--max_ruwe", type=float, default=1.4)
    g_q.add_argument("--parallax_snr_min", type=float, default=10.0)
    g_q.add_argument("--phot_g_snr_min", type=float, default=50.0)
    g_q.add_argument("--use_excess_factor_filter", action="store_true")
    g_q.add_argument("--use_astro_chi2_filter", action="store_true")

    g_sel = p.add_argument_group("Filtros Físicos (Selección)")
    g_sel.add_argument("--max_dist_pc", type=float)
    g_sel.add_argument("--min_abs_b_gal", type=float)

    return p


def main() -> int:
    parser = get_parser()
    args = parser.parse_args()

    console = configure_logging(args.log_file, args.verbose)

    logging.info("=" * 80)
    logging.info("PIPELINE ENRIQUECIMIENTO GAIA DR3 v3")
    logging.info("=" * 80)
    logging.info(f"Args: {vars(args)}")

    try:
        final_df = run_pipeline(args, COLUMN_CONFIG, console)
        if final_df is not None:
            logging.info("\n✅ PROCESO COMPLETADO")
            return 0
        else:
            logging.error("\n❌ PROCESO TERMINÓ VACÍO")
            return 1
    except (DataQualityError, ContractViolationError) as e:
        logging.error(f"\n❌ {type(e).__name__}: {e}")
        if console:
            console.print(f"[bold red]❌ {type(e).__name__}:[/bold red] {e}")
        return 1
    except KeyboardInterrupt:
        logging.warning("\n⚠️ Interrumpido por usuario")
        return 130
    except Exception as e:
        logging.critical(f"\n❌ Error no controlado: {e}", exc_info=True)
        return 2


if __name__ == "__main__":
    sys.exit(main())
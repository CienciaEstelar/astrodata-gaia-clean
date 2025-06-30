# -*- coding: utf-8 -*-
"""
Procesa un catálogo CSV de Gaia, aplica filtros de calidad, calcula parámetros
físicos derivados y guarda un CSV listo para análisis científico.

Uso rápido
~~~~~~~~~~
python gaia_enrichment_pipeline.py \
  --input CatalogoPropio_Gaia.csv \
  --output Catalogo_Gaia_Enriquecido_v2.csv \
  --use_excess_factor_filter \
  --use_astro_chi2_filter

"""
# ─────────────────────────── 1. IMPORTS ────────────────────────────
import argparse
import logging
import os
from typing import Callable, List, Tuple

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from rich.progress import (
    Progress, BarColumn, TextColumn,
    TimeElapsedColumn, TimeRemainingColumn
)

# ─────────── 2. CONSTANTES Y CONFIG CENTRALIZADA ────────────
MAS_YR_TO_KM_S = 4.74047   # mas yr⁻¹ → km s⁻¹ pc⁻¹ (factor Vtan)

COLUMN_CONFIG = {
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
    "rv": "dr2_radial_velocity",
    "teff": "dr2_rv_template_teff",
    "extinction_g": "a_g_val",
    "var_flag": "phot_variable_flag",
}

# ─────────── 3. LOGGING ────────────
def setup_logging(log_file: str = "enriquecimiento.log", verbose: bool = False) -> None:
    """Configura logging dual (archivo + consola)."""
    log = logging.getLogger()
    log.setLevel(logging.DEBUG if verbose else logging.INFO)
    if log.hasHandlers():
        log.handlers.clear()

    file_fmt = logging.Formatter("%(asctime)s - %(levelname)-8s - %(message)s")
    console_fmt = logging.Formatter("%(message)s")

    fh = logging.FileHandler(log_file, mode="w", encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(file_fmt)

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG if verbose else logging.INFO)
    ch.setFormatter(console_fmt)

    log.addHandler(fh)
    log.addHandler(ch)

# ─────────── 4. CLASIFICACIONES VECTORIAL ────────────
def clasificacion_espectral_vectorizada(teff: pd.Series) -> pd.Series:
    bins = [-np.inf, 3700, 5200, 6000, 7500, 10000, 30000, np.inf]
    labels = ["M", "K", "G", "F", "A", "B", "O"]
    return pd.cut(teff, bins=bins, labels=labels, right=True).astype("string")

def clasificar_fase_evolutiva_vectorizado(df: pd.DataFrame,
                                          color_col: str,
                                          mg_col: str) -> pd.Series:
    c, M = df[color_col].to_numpy(), df[mg_col].to_numpy()
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
        "Indeterminado (Faltan Datos)",
        "Posible Enana Blanca",
        "Gigante (Roja)",
        "Gigante/Subgigante (Azul/Blanca)",
        "Secuencia Principal (Enana Roja/M)",
        "Secuencia Principal (Tipo G/K)",
        "Secuencia Principal (Tipo F/A)",
    ]
    return pd.Series(
        np.select(conds, choices, default="Secuencia Principal (Tipo O/B)"),
        index=df.index,
        dtype="string"
    )

def indice_variabilidad_vectorizado(ruwe: pd.Series) -> pd.Series:
    conds = [ruwe.isna(), ruwe > 1.4, ruwe > 1.0]
    choices = [np.nan, "Alta (RUWE > 1.4)", "Moderada (1.0 < RUWE ≤ 1.4)"]
    return pd.Series(
        np.select(conds, choices, default="Baja (RUWE ≤ 1.0)"),
        index=ruwe.index,
        dtype="string"
    )

# ─────────── 5. FILTRADO DE CALIDAD ────────────
Filtro = Tuple[str, List[str], Callable[[pd.DataFrame], pd.Series]]

def build_quality_mask(df: pd.DataFrame, filtros: List[Filtro]) -> pd.Series:
    mask = pd.Series(True, index=df.index)
    for nombre, req, fn in filtros:
        missing = [c for c in req if c not in df.columns]
        if missing:
            logging.warning(f"Filtro '{nombre}' omitido (faltan: {', '.join(missing)})")
            continue
        cond = fn(df).reindex(df.index, fill_value=False)
        logging.info(f"Filtro '{nombre}': descarta {(~cond).sum():,} filas")
        mask &= cond
    return mask

def crear_filtros_calidad(args, cfg) -> List[Filtro]:
    return [
        (f"RUWE < {args.max_ruwe}", [cfg["ruwe"]],
         lambda d: d[cfg["ruwe"]] < args.max_ruwe),

        (f"Parallax S/N > {args.parallax_snr_min}",
         [cfg["parallax"], cfg["parallax_error"]],
         lambda d: (d[cfg["parallax_error"]] > 0) &
                   ((d[cfg["parallax"]] /
                     d[cfg["parallax_error"]]) > args.parallax_snr_min)),

        (f"Phot G S/N > {args.phot_g_snr_min}",
         [cfg["g_flux"], cfg["g_flux_error"]],
         lambda d: (d[cfg["g_flux_error"]] > 0) &
                   ((d[cfg["g_flux"]] /
                     d[cfg["g_flux_error"]]) > args.phot_g_snr_min)),

        (f"Phot BP S/N > {args.phot_bp_snr_min}",
         [cfg["bp_flux"], cfg["bp_flux_error"]],
         lambda d: (d[cfg["bp_flux_error"]] > 0) &
                   ((d[cfg["bp_flux"]] /
                     d[cfg["bp_flux_error"]]) > args.phot_bp_snr_min)),

        (f"Phot RP S/N > {args.phot_rp_snr_min}",
         [cfg["rp_flux"], cfg["rp_flux_error"]],
         lambda d: (d[cfg["rp_flux_error"]] > 0) &
                   ((d[cfg["rp_flux"]] /
                     d[cfg["rp_flux_error"]]) > args.phot_rp_snr_min)),

        ("BP/RP Excess",
         [cfg["bp_rp_excess"], cfg["bp_rp_color"]],
         lambda d: pd.Series(True, d.index)
         if not args.use_excess_factor_filter else
         ((d[cfg["bp_rp_color"]].notna()) &
          (d[cfg["bp_rp_excess"]] >
           1.0 + 0.015 * d[cfg["bp_rp_color"]]**2) &
          (d[cfg["bp_rp_excess"]] <
           1.3 + 0.060 * d[cfg["bp_rp_color"]]**2))),

        ("Astro χ²/DOF",
         [cfg["chi2_astro"], cfg["n_good_obs_astro"], cfg["g_mag"]],
         lambda d: pd.Series(True, d.index)
         if not args.use_astro_chi2_filter else
         ((d[cfg["n_good_obs_astro"]] > 5) &
          ((d[cfg["chi2_astro"]] /
            (d[cfg["n_good_obs_astro"]] - 5)) <
           1.44 * np.maximum(1.0,
                             np.exp(-0.4 *
                                    (d[cfg["g_mag"]] - 19.5)))))),

        (f"Visibility periods > {args.vis_periods_min}",
         [cfg["vis_periods"]],
         lambda d: d[cfg["vis_periods"]] > args.vis_periods_min),
    ]

def apply_quality_filters(df: pd.DataFrame, args, cfg) -> pd.DataFrame:
    logging.info("\n" + "─"*28 + " Paso 2: Filtros de calidad " + "─"*28)
    if args.use_excess_factor_filter and cfg["bp_rp_color"] not in df.columns:
        df[cfg["bp_rp_color"]] = df[cfg["bp_mag"]] - df[cfg["rp_mag"]]
    mask = build_quality_mask(df, crear_filtros_calidad(args, cfg))
    df_out = df[mask].copy()
    logging.info(f"► Sobreviven {len(df_out):,} / {len(df):,} filas\n")
    return df_out

# ─────────── 6. CÁLCULOS DERIVADOS ────────────
def calculate_derived_parameters(df: pd.DataFrame, cfg) -> pd.DataFrame:
    logging.info("─"*25 + " Paso 3: Cálculos derivados " + "─"*24)

    # 1. Limpia a numérico todas las columnas “críticas” (ahora incluye ra y dec)
    for col in [cfg["parallax"], cfg["g_mag"], cfg["extinction_g"],
                cfg["pmra"], cfg["pmdec"], cfg["rv"], cfg["teff"],
                cfg["ruwe"], cfg["ra"], cfg["dec"]]:          # ← añadido ra/dec
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # 2. Distancia y magnitud absoluta
    df["dist_pc"] = np.where(df[cfg["parallax"]] > 0,
                             1000.0 / df[cfg["parallax"]], np.nan)

    df["M_G_app"] = df[cfg["g_mag"]] - 5 * np.log10(df["dist_pc"]) + 5
    df["M_G"] = df["M_G_app"] - df.get(cfg["extinction_g"], 0)

    # 3. Color BP-RP (si falta)
    if cfg["bp_rp_color"] not in df.columns:
        df[cfg["bp_rp_color"]] = (
            pd.to_numeric(df[cfg["bp_mag"]], errors="coerce") -
            pd.to_numeric(df[cfg["rp_mag"]], errors="coerce")
        )

    # 4. Velocidades
    df["V_tan_ra"]  = MAS_YR_TO_KM_S * df[cfg["pmra"]]  * df["dist_pc"]
    df["V_tan_dec"] = MAS_YR_TO_KM_S * df[cfg["pmdec"]] * df["dist_pc"]
    df["V_radial"]  = df[cfg["rv"]]
    df["V_3D"] = np.sqrt(df["V_tan_ra"]**2 + df["V_tan_dec"]**2 +
                         df["V_radial"].fillna(0)**2)

    # 5. Coordenadas galácticas
    if {cfg["ra"], cfg["dec"]} <= set(df.columns):
        mask = df[cfg["ra"]].notna() & df[cfg["dec"]].notna()
        if mask.any():
            icrs = SkyCoord(df.loc[mask, cfg["ra"]],  # ya son floats
                            df.loc[mask, cfg["dec"]],
                            unit="deg", frame="icrs")  # unidad explícita
            gal = icrs.galactic
            df.loc[mask, "l_gal"] = gal.l.degree
            df.loc[mask, "b_gal"] = gal.b.degree

    # 6. Clasificaciones rápidas
    df["tipo_espectral"]  = clasificacion_espectral_vectorizada(df[cfg["teff"]])
    df["fase_evolutiva"]  = clasificar_fase_evolutiva_vectorizado(
                                df, cfg["bp_rp_color"], "M_G")
    df["variabilidad"]    = indice_variabilidad_vectorizado(df[cfg["ruwe"]])

    return df


# ─────────── 7. FILTROS DE SELECCIÓN ────────────
def apply_selection_filters(df: pd.DataFrame, args) -> pd.DataFrame:
    logging.info("\n" + "─"*28 + " Paso 4: Filtros de selección " + "─"*29)
    initial = len(df)
    if args.max_dist_pc is not None and "dist_pc" in df.columns:
        df = df[df["dist_pc"] <= args.max_dist_pc].copy()
        logging.info(f"Distancia ≤ {args.max_dist_pc} pc → {len(df):,} filas")
    if args.min_abs_b_gal is not None and "b_gal" in df.columns:
        df = df[df["b_gal"].abs() >= args.min_abs_b_gal].copy()
        logging.info(f"|b_gal| ≥ {args.min_abs_b_gal}° → {len(df):,} filas")
    logging.info(f"► Final tras selección: {len(df):,}/{initial:,} filas\n")
    return df

# ─────────── 8. PIPELINE CON PROGRESO RICH ────────────
def run_enrichment_pipeline(args, cfg) -> pd.DataFrame | None:
    STEPS = (
        ("Leyendo CSV",                     "Paso 1/5"),
        ("Aplicando filtros de calidad",    "Paso 2/5"),
        ("Calculando parámetros derivados", "Paso 3/5"),
        ("Aplicando filtros de selección",  "Paso 4/5"),
        ("Guardando CSV final",             "Paso 5/5"),
    )

    bar_cols = [
        TextColumn("[progress.description]{task.description}"),
        BarColumn(bar_width=None),
        TextColumn("{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
    ]

    with Progress(*bar_cols, transient=False) as progress:
        task = progress.add_task("[cyan]Preparando…", total=len(STEPS))

        # ─ Paso 1
        progress.update(task, description=f"[cyan]{STEPS[0][1]}: {STEPS[0][0]}…")
        try:
            df = pd.read_csv(args.input, low_memory=False)
            original_n = len(df)
            progress.advance(task)
        except FileNotFoundError:
            logging.critical("Archivo de entrada no encontrado.")
            return None
        except Exception:
            logging.critical("Fallo al leer CSV:", exc_info=True)
            return None

        # ─ Paso 2
        progress.update(task, description=f"[cyan]{STEPS[1][1]}: {STEPS[1][0]}…")
        df = apply_quality_filters(df, args, cfg)
        if df.empty:
            logging.error("Todas las filas fueron descartadas por calidad.")
            return None
        progress.advance(task)

        # ─ Paso 3
        progress.update(task, description=f"[cyan]{STEPS[2][1]}: {STEPS[2][0]}…")
        df = calculate_derived_parameters(df, cfg)
        progress.advance(task)

        # ─ Paso 4
        progress.update(task, description=f"[cyan]{STEPS[3][1]}: {STEPS[3][0]}…")
        df = apply_selection_filters(df, args)
        if df.empty:
            logging.error("Todas las filas fueron descartadas por selección.")
            return None
        progress.advance(task)

        # ─ Paso 5
        progress.update(task, description=f"[green]{STEPS[4][1]}: {STEPS[4][0]}…")
        try:
            os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
            df.to_csv(args.output, index=False, float_format="%.8f")
            progress.advance(task)
        except Exception:
            logging.critical("No se pudo guardar el CSV final:", exc_info=True)
            return None

    return df

# ─────────── 9. CLI ────────────
def make_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Pipeline de enriquecimiento Gaia",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    g_files = p.add_argument_group("Archivos y logging")
    g_files.add_argument("--input", required=True, help="CSV de entrada de Gaia")
    g_files.add_argument("--output", required=True, help="CSV de salida enriquecido")
    g_files.add_argument("--log_file", default="enriquecimiento.log", help="Archivo de log")
    g_files.add_argument("-v", "--verbose", action="store_true", help="Log detallado (DEBUG)")

    g_q = p.add_argument_group("Filtros de calidad")
    g_q.add_argument("--max_ruwe", type=float, default=1.4)
    g_q.add_argument("--parallax_snr_min", type=float, default=10)
    g_q.add_argument("--phot_g_snr_min", type=float, default=50)
    g_q.add_argument("--phot_bp_snr_min", type=float, default=20)
    g_q.add_argument("--phot_rp_snr_min", type=float, default=20)
    g_q.add_argument("--use_excess_factor_filter", action="store_true")
    g_q.add_argument("--use_astro_chi2_filter", action="store_true")
    g_q.add_argument("--vis_periods_min", type=int, default=8)

    g_sel = p.add_argument_group("Filtros de selección")
    g_sel.add_argument("--max_dist_pc", type=float)
    g_sel.add_argument("--min_abs_b_gal", type=float)

    return p

# ─────────── 10. MAIN ────────────
def main() -> None:
    parser = make_parser()
    args = parser.parse_args()
    setup_logging(args.log_file, args.verbose)

    logging.info("="*80)
    logging.info("PIPELINE DE ENRIQUECIMIENTO GAIA — INICIO")
    logging.info(f"Argumentos: {vars(args)}")
    logging.info("="*80)

    try:
        final_df = run_enrichment_pipeline(args, COLUMN_CONFIG)
        if final_df is not None:
            logging.info("PROCESO COMPLETADO ✔️")
        else:
            logging.error("PROCESO TERMINÓ CON ERRORES ✖️")
    except KeyboardInterrupt:
        logging.warning("Proceso interrumpido por el usuario (Ctrl-C).")

    logging.info("="*80)

if __name__ == "__main__":
    main()

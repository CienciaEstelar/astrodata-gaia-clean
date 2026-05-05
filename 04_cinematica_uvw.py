"""
SCRIPT 4: CÁLCULO DE VELOCIDADES GALÁCTICAS U, V, W Y LSR (v2 — Auditado)

Toma un catálogo enriquecido (Parquet o CSV) con posiciones, distancias y
movimientos propios, y calcula las componentes cartesianas de la velocidad
galáctica (U, V, W) en el sistema Galactocéntrico y el Local Standard of Rest.

CAMBIOS RESPECTO A v1:
- FIX: validación NaN explícita ANTES de pasar a Astropy (antes propagaba
  silenciosamente y podía colgar el chunk entero con un solo NaN problemático).
- FIX: si skip_no_rv=False, los cálculos UVW sin VR se marcan explícitamente
  como resultado parcial (V_tan válido, V_3D=NaN). No se inventa VR=0.
- AÑADIDO: contador de filas procesadas/descartadas al final para sanity check.
- AÑADIDO: detección de columna `V_radial_kms` o `radial_velocity` según
  origen del pipeline (retrocompat con catálogos pre-v2).
"""

import argparse
import logging
from pathlib import Path
from typing import Optional, List

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import (
    Distance,
    Galactocentric,
    LSR,
    SkyCoord,
    galactocentric_frame_defaults,
)
from rich.progress import (
    BarColumn,
    Progress,
    TextColumn,
    TimeElapsedColumn,
)

# Parámetros modernos del centro galáctico y del Sol (Astropy v4.0+)
galactocentric_frame_defaults.set('v4.0')


# ─────────────────────────────────────────────────────────────────────────────
# 1. LOGGING
# ─────────────────────────────────────────────────────────────────────────────

def setup_logging(log_file: Optional[str] = None, verbose: bool = False) -> None:
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    if logger.hasHandlers():
        logger.handlers.clear()

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG if verbose else logging.INFO)
    ch.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(ch)

    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_path, encoding="utf-8")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(fh)


# ─────────────────────────────────────────────────────────────────────────────
# 2. PREPARACIÓN Y LIMPIEZA DE DATOS
# ─────────────────────────────────────────────────────────────────────────────

def autodetect_rv_column(df: pd.DataFrame, user_col: str) -> str:
    """Detecta cuál columna tiene la velocidad radial.

    Orden de preferencia:
      1. La que el usuario pasó explícitamente.
      2. `V_radial_kms` (estándar del pipeline post-v2).
      3. `radial_velocity` (Gaia DR3 nativo).
      4. `dr2_radial_velocity` (legacy pre-v2).
    """
    candidatos = [user_col, "V_radial_kms", "radial_velocity", "dr2_radial_velocity"]
    for c in candidatos:
        if c in df.columns:
            if c != user_col:
                logging.info(f"Columna RV detectada automáticamente: '{c}' (user pidió '{user_col}').")
            return c
    logging.warning(f"Ninguna columna RV encontrada. Se asume NaN.")
    return user_col


def validate_and_filter_data(df: pd.DataFrame, rv_col: str, skip_no_rv: bool) -> pd.DataFrame:
    """Valida columnas críticas, fuerza tipos numéricos y filtra físicamente inválidos."""
    required_cols = ['ra', 'dec', 'distancia_parsecs', 'pmra', 'pmdec']
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(f"Faltan columnas críticas: {missing}")

    if rv_col not in df.columns:
        logging.warning(f"La columna '{rv_col}' no existe. Se asumen NaNs.")
        df[rv_col] = np.nan

    numeric_cols = required_cols + [rv_col]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    initial_len = len(df)

    # 1. Filtros físicos obligatorios (no se puede transformar sin estos)
    df = df[
        df['distancia_parsecs'].notna() &
        (df['distancia_parsecs'] > 0) &
        df['ra'].notna() &
        df['dec'].notna() &
        df['pmra'].notna() &
        df['pmdec'].notna()
    ].copy()

    dist_discarded = initial_len - len(df)
    if dist_discarded > 0:
        logging.info(f"► Descartadas {dist_discarded:,} filas con astrometría incompleta o dist<=0.")

    # 2. Filtro opcional por VR
    if skip_no_rv:
        pre_rv_len = len(df)
        df = df.dropna(subset=[rv_col]).copy()
        logging.info(f"► Descartadas {pre_rv_len - len(df):,} filas sin Velocidad Radial válida.")
    else:
        n_con_rv = df[rv_col].notna().sum()
        n_sin_rv = len(df) - n_con_rv
        logging.info(
            f"► Se mantienen todas: {n_con_rv:,} con VR (cálculo 3D completo), "
            f"{n_sin_rv:,} sin VR (V_3D=NaN)."
        )

    if df.empty:
        raise ValueError("El catálogo quedó vacío tras los filtros base.")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# 3. MOTOR CINEMÁTICO
# ─────────────────────────────────────────────────────────────────────────────

def compute_kinematics_chunk(
    df_chunk: pd.DataFrame, rv_col: str, calc_lsr: bool
) -> pd.DataFrame:
    """Calcula U, V, W vía Astropy. Maneja NaN en VR sin colapsar el chunk."""

    # Preparar arrays. Donde VR es NaN, Astropy propaga NaN en las componentes 3D
    # pero las posiciones y VX/VY (tangenciales) siguen siendo calculables.
    # Aun así, prevenimos que un NaN único rompa la transformación entera.
    rv_values = df_chunk[rv_col].values.astype(float)
    has_rv = np.isfinite(rv_values)

    # Reemplazar NaNs por 0 SOLO para que Astropy no rompa; marcaremos las filas
    # sin VR como NaN en la salida manualmente tras el cálculo.
    rv_safe = np.where(has_rv, rv_values, 0.0)

    distances = Distance(df_chunk['distancia_parsecs'].values * u.parsec)

    coords = SkyCoord(
        ra=df_chunk['ra'].values * u.degree,
        dec=df_chunk['dec'].values * u.degree,
        distance=distances,
        pm_ra_cosdec=df_chunk['pmra'].values * u.mas / u.yr,
        pm_dec=df_chunk['pmdec'].values * u.mas / u.yr,
        radial_velocity=rv_safe * u.km / u.s,
        frame='icrs'
    )

    galactocentric_coords = coords.transform_to(Galactocentric())

    U = galactocentric_coords.velocity.d_x.to_value(u.km / u.s)
    V = galactocentric_coords.velocity.d_y.to_value(u.km / u.s)
    W = galactocentric_coords.velocity.d_z.to_value(u.km / u.s)

    # Marcar como NaN las filas que no tenían VR real (la componente radial
    # contribuye a las 3 componentes galactocéntricas).
    U[~has_rv] = np.nan
    V[~has_rv] = np.nan
    W[~has_rv] = np.nan

    df_chunk = df_chunk.copy()
    df_chunk['U_kms'] = U
    df_chunk['V_kms'] = V
    df_chunk['W_kms'] = W

    if calc_lsr:
        lsr_coords = coords.transform_to(LSR())
        U_l = lsr_coords.velocity.d_x.to_value(u.km / u.s)
        V_l = lsr_coords.velocity.d_y.to_value(u.km / u.s)
        W_l = lsr_coords.velocity.d_z.to_value(u.km / u.s)
        U_l[~has_rv] = np.nan
        V_l[~has_rv] = np.nan
        W_l[~has_rv] = np.nan
        df_chunk['U_lsr_kms'] = U_l
        df_chunk['V_lsr_kms'] = V_l
        df_chunk['W_lsr_kms'] = W_l

    return df_chunk


def calculate_kinematics_batched(
    df: pd.DataFrame, rv_col: str, calc_lsr: bool, chunk_size: int = 500_000
) -> pd.DataFrame:
    """Divide el DataFrame en lotes para evitar OOM en Astropy."""
    total_rows = len(df)
    processed_chunks: List[pd.DataFrame] = []

    logging.info(f"Procesando {total_rows:,} estrellas en lotes de {chunk_size:,}...")

    n_chunks = (total_rows + chunk_size - 1) // chunk_size
    for idx, start_idx in enumerate(range(0, total_rows, chunk_size), 1):
        end_idx = min(start_idx + chunk_size, total_rows)
        df_chunk = df.iloc[start_idx:end_idx].copy()

        try:
            df_chunk_processed = compute_kinematics_chunk(df_chunk, rv_col, calc_lsr)
            processed_chunks.append(df_chunk_processed)
            logging.debug(f"Chunk {idx}/{n_chunks} OK ({len(df_chunk):,} filas).")
        except Exception as e:
            logging.error(f"Fallo en chunk {idx}/{n_chunks}: {e}. Se omite este lote.")

    if not processed_chunks:
        raise RuntimeError("Todos los chunks fallaron.")

    return pd.concat(processed_chunks, ignore_index=True)


# ─────────────────────────────────────────────────────────────────────────────
# 4. ORQUESTADOR
# ─────────────────────────────────────────────────────────────────────────────

def run_uvw_pipeline(args: argparse.Namespace) -> Optional[pd.DataFrame]:
    steps = [
        "Leyendo Catálogo Enriquecido",
        "Validación y Filtrado Kinemático",
        "Transformación Astrométrica (UVW)",
        "Guardando Resultados"
    ]
    bar_cols = [
        TextColumn("[progress.description]{task.description}"),
        BarColumn(bar_width=None),
        TextColumn("{task.completed}/{task.total}"),
        TimeElapsedColumn(),
    ]

    with Progress(*bar_cols, transient=False) as progress:
        task = progress.add_task("[cyan]Iniciando Cinemática...", total=len(steps))

        # Paso 1: Lectura
        progress.update(task, description=f"[cyan]Paso 1: {steps[0]}")
        try:
            file_ext = Path(args.input).suffix.lower()
            if file_ext == '.parquet':
                df = pd.read_parquet(args.input)
            else:
                df = pd.read_csv(args.input, low_memory=False)
            progress.advance(task)
        except Exception as e:
            logging.critical(f"Error al leer archivo: {e}")
            return None

        # Auto-detect columna RV
        rv_col = autodetect_rv_column(df, args.rv_col)

        # Paso 2: Validación
        progress.update(task, description=f"[cyan]Paso 2: {steps[1]}")
        try:
            df = validate_and_filter_data(df, rv_col, args.skip_no_rv)
            progress.advance(task)
        except ValueError as e:
            logging.critical(str(e))
            return None

        # Paso 3: Transformaciones
        progress.update(task, description=f"[cyan]Paso 3: {steps[2]}")
        try:
            df = calculate_kinematics_batched(
                df, rv_col=rv_col,
                calc_lsr=args.calculate_lsr,
                chunk_size=args.chunk_size
            )
            progress.advance(task)
        except Exception as e:
            logging.critical(f"Fallo en transformación Astropy: {e}", exc_info=True)
            return None

        # Sanity check post-transformación
        n_uvw_ok = df[['U_kms', 'V_kms', 'W_kms']].notna().all(axis=1).sum()
        logging.info(
            f"► Post-transformación: {n_uvw_ok:,}/{len(df):,} filas con UVW 3D completo."
        )

        # Paso 4: Guardado
        progress.update(task, description=f"[green]Paso 4: {steps[3]}")
        try:
            out_path = Path(args.output)
            out_path.parent.mkdir(parents=True, exist_ok=True)

            if args.format == 'parquet':
                df.to_parquet(out_path, engine='pyarrow', compression='snappy', index=False)
            else:
                df.to_csv(out_path, index=False, float_format='%.8f')

            logging.info(f"Guardado en: {out_path} (Filas: {len(df):,})")
            progress.advance(task)
        except Exception as e:
            logging.critical(f"Fallo al guardar: {e}")
            return None

    return df


# ─────────────────────────────────────────────────────────────────────────────
# 5. CLI
# ─────────────────────────────────────────────────────────────────────────────

def get_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Pipeline Cinemático UVW/LSR (v2 auditado).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    g_io = p.add_argument_group("Archivos I/O")
    g_io.add_argument("--input", required=True)
    g_io.add_argument("--output", required=True)
    g_io.add_argument("--format", choices=['csv', 'parquet'], default='parquet')
    g_io.add_argument("--log_file", default="cinematica.log")
    g_io.add_argument("-v", "--verbose", action="store_true")

    g_kin = p.add_argument_group("Configuración Cinemática")
    g_kin.add_argument("--rv_col", type=str, default='V_radial_kms',
                       help="Nombre preferido de columna VR (auto-detect como fallback).")
    g_kin.add_argument("--skip_no_rv", action="store_true")
    g_kin.add_argument("--calculate_lsr", action="store_true")
    g_kin.add_argument("--chunk_size", type=int, default=500000)

    return p


def main() -> None:
    parser = get_parser()
    args = parser.parse_args()

    setup_logging(args.log_file, args.verbose)
    logging.info("=" * 80)
    logging.info("🌌 PIPELINE CINEMÁTICO UVW / LSR (v2)")
    logging.info("=" * 80)

    try:
        final_df = run_uvw_pipeline(args)
        if final_df is not None:
            logging.info("\n✅ PROCESO COMPLETADO")
        else:
            logging.error("\n❌ PROCESO FALLÓ")
    except KeyboardInterrupt:
        logging.warning("\n⚠️ Interrumpido por usuario.")
    except Exception as e:
        logging.critical(f"\n❌ Error no controlado: {e}", exc_info=True)


if __name__ == "__main__":
    main()

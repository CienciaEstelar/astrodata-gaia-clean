"""
SCRIPT 02: GAIA DR3 ADQL FOOTPRINT EXTRACTOR (v3 — UX Refinada)

Extrae el WCS de los FITS de SDSS y construye consultas ADQL
para descargar el catálogo exacto desde los servidores TAP de Gaia.

CAMBIOS RESPECTO A v2:
- Consola limpia: solo barra de progreso Rich + errores críticos.
- Warnings de Astropy (FITSFixedWarning) y logs de astroquery
  (Query finished, etc.) van SOLO al archivo log, no a pantalla.
- Control fino de verbosity: -v para mostrar más en consola.
- La barra Rich persiste correctamente incluso con | tee.
"""

import argparse
import logging
import random
import sys
import time
import warnings
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.exceptions import AstropyWarning
from astroquery.gaia import Gaia

try:
    from rich.console import Console
    from rich.progress import (
        BarColumn, Progress, TextColumn,
        TimeElapsedColumn, TimeRemainingColumn, MofNCompleteColumn,
    )
    HAS_RICH = True
except ImportError:
    HAS_RICH = False

# ─────────────────────────────────────────────────────────────────────────────
# COLUMNAS DR3
# ─────────────────────────────────────────────────────────────────────────────
GAIA_COLS = [
    "source_id", "solution_id", "ra", "dec", "parallax", "parallax_error",
    "pmra", "pmdec", "ruwe", "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag",
    "phot_g_mean_flux", "phot_g_mean_flux_error",
    "phot_bp_mean_flux", "phot_bp_mean_flux_error",
    "phot_rp_mean_flux", "phot_rp_mean_flux_error",
    "bp_rp", "phot_bp_rp_excess_factor",
    "astrometric_chi2_al", "astrometric_n_good_obs_al", "visibility_periods_used",
    "radial_velocity", "radial_velocity_error",
    "teff_gspphot", "logg_gspphot", "mh_gspphot",
    "ag_gspphot", "ebpminrp_gspphot",
    "phot_variable_flag"
]


# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURACIÓN DE LOGGING Y SILENCIADO DE WARNINGS
# ─────────────────────────────────────────────────────────────────────────────

def configure_logging_and_silence(log_file: str, verbose: bool) -> None:
    """Logging dual: consola minimal, archivo completo.

    Clave: los warnings de Astropy y los logs de astroquery se REDIRIGEN
    al archivo log en vez de ser impresos a consola.
    """
    # 1. Silenciar TODO warning de Astropy en la salida estándar
    warnings.filterwarnings("ignore", category=AstropyWarning)
    # Redirigir warnings al sistema de logging (que luego los manda al archivo)
    logging.captureWarnings(True)

    # 2. Logger raíz: nivel DEBUG (para que el archivo capture todo)
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    if root.hasHandlers():
        root.handlers.clear()

    # 3. Handler de archivo: TODO (DEBUG+) con timestamp
    log_path = Path(log_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(log_path, encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    ))
    root.addHandler(fh)

    # 4. Handler de consola: SOLO WARNING+ (o INFO si --verbose)
    console_level = logging.INFO if verbose else logging.WARNING
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(console_level)
    ch.setFormatter(logging.Formatter("%(message)s"))
    root.addHandler(ch)

    # 5. Silenciar específicamente astroquery en consola (deja archivo intacto)
    #    astroquery tiene su propio logger llamado "astroquery" que imprime
    #    "Query finished" y similares. Le subimos el nivel de CONSOLA.
    for logger_name in ("astroquery", "astroquery.gaia", "astroquery.utils.tap"):
        lg = logging.getLogger(logger_name)
        lg.setLevel(logging.DEBUG)   # al archivo va todo
        lg.propagate = True           # hereda handlers del root

    # 6. astropy también tiene su logger, mismo tratamiento
    logging.getLogger("astropy").setLevel(logging.DEBUG)


# ─────────────────────────────────────────────────────────────────────────────
# FUNCIONES CORE
# ─────────────────────────────────────────────────────────────────────────────

def wcs_to_box(fits_path: Path, margin_deg: float) -> Tuple[float, float, float, float]:
    """Extrae centro (RA, Dec) y dimensiones (w, h) del footprint del FITS."""
    with fits.open(fits_path, memmap=False) as hdul:
        h = hdul[0].header
        w = WCS(h)
        nx, ny = int(h.get("NAXIS1", 0)), int(h.get("NAXIS2", 0))

        pix = np.array([[1, 1], [nx, 1], [nx, ny], [1, ny]], dtype=float)
        sky = w.pixel_to_world(pix[:, 0], pix[:, 1])

        ra, dec = sky.ra.deg, sky.dec.deg
        ra_min, ra_max = float(np.min(ra)), float(np.max(ra))
        dec_min, dec_max = float(np.min(dec)), float(np.max(dec))

        ra_c = 0.5 * (ra_min + ra_max)
        dec_c = 0.5 * (dec_min + dec_max)

        cos_dec = max(np.cos(np.deg2rad(dec_c)), 1e-3)
        width = max((ra_max - ra_min) * cos_dec + 2 * margin_deg, 1e-4)
        height = max((dec_max - dec_min) + 2 * margin_deg, 1e-4)

        return ra_c, dec_c, width, height


def execute_tap_query(adql: str, retries: int, backoff: float, sleep_s: float) -> pd.DataFrame:
    """Ejecuta consulta TAP con reintentos exponenciales y jitter."""
    time.sleep(sleep_s)
    last_exc = None
    for attempt in range(max(1, retries)):
        try:
            Gaia.ROW_LIMIT = -1
            # verbose=False suprime el "Query finished" en stdout
            job = Gaia.launch_job_async(adql, dump_to_file=False, verbose=False)
            return job.get_results().to_pandas()
        except Exception as e:
            last_exc = e
            logging.warning(f"Intento {attempt+1}/{retries} falló: {e}")
            time.sleep(sleep_s * (backoff ** attempt) + random.random() * 0.2)
    raise RuntimeError(f"Fallo crítico conectando a Gaia TAP: {last_exc}")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="Script 02: Consulta ADQL a Gaia DR3 (v3 — UX refinada)."
    )
    p.add_argument("--fits-dir", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--band", default="r")
    p.add_argument("--margin-deg", type=float, default=0.02)
    p.add_argument("--limit", type=int, default=200000)
    p.add_argument("--log-file", default="02_consulta_gaia.log",
                   help="Archivo de log detallado (todos los warnings van aquí).")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Mostrar mensajes INFO en consola (default: solo WARNING+).")
    args = p.parse_args()

    # Configurar logging dual + silenciar warnings
    configure_logging_and_silence(args.log_file, args.verbose)

    # Silenciar también astroquery globalmente (por si algún módulo imprime)
    Gaia.ROW_LIMIT = -1

    # Buscar FITS
    fits_files = sorted(Path(args.fits_dir).glob(f"frame-{args.band}-*.fits"))
    if not fits_files:
        raise FileNotFoundError(f"No se encontraron FITS en {args.fits_dir}")

    console = Console() if HAS_RICH else None
    if console:
        console.print(
            f"[bold cyan]🌌 Gaia DR3 ADQL Extractor[/bold cyan]\n"
            f"  FITS encontrados: [green]{len(fits_files)}[/green] "
            f"(banda [yellow]{args.band}[/yellow])\n"
            f"  Log detallado:    [dim]{args.log_file}[/dim]\n"
        )

    chunks = []
    failures = []

    if HAS_RICH:
        progress = Progress(
            TextColumn("[bold cyan]{task.description}"),
            BarColumn(bar_width=40),
            MofNCompleteColumn(),
            TextColumn("•"),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TextColumn("•"),
            TimeElapsedColumn(),
            TextColumn("/"),
            TimeRemainingColumn(),
            TextColumn("• [yellow]{task.fields[stars]:,}[/yellow] ⭐"),
            console=console,
            refresh_per_second=4,
        )

        with progress:
            task = progress.add_task(
                "Consultando Gaia DR3",
                total=len(fits_files),
                stars=0,
            )
            total_stars = 0

            for fpath in fits_files:
                try:
                    ra_c, dec_c, w, h = wcs_to_box(fpath, args.margin_deg)
                    cols = ", ".join(GAIA_COLS)
                    adql = (
                        f"SELECT TOP {args.limit} {cols} "
                        f"FROM gaiadr3.gaia_source "
                        f"WHERE 1=CONTAINS(POINT('ICRS', ra, dec), "
                        f"BOX('ICRS', {ra_c:.8f}, {dec_c:.8f}, {w:.8f}, {h:.8f}))"
                    )
                    df = execute_tap_query(adql, retries=5, backoff=1.5, sleep_s=1.0)
                    if not df.empty:
                        chunks.append(df)
                        total_stars += len(df)
                    logging.info(f"OK {fpath.name}: {len(df):,} estrellas")
                except Exception as e:
                    failures.append((fpath.name, str(e)))
                    logging.error(f"FAIL {fpath.name}: {e}")

                progress.advance(task)
                progress.update(task, stars=total_stars)
    else:
        # Fallback sin Rich
        total_stars = 0
        for i, fpath in enumerate(fits_files, 1):
            try:
                ra_c, dec_c, w, h = wcs_to_box(fpath, args.margin_deg)
                cols = ", ".join(GAIA_COLS)
                adql = (
                    f"SELECT TOP {args.limit} {cols} "
                    f"FROM gaiadr3.gaia_source "
                    f"WHERE 1=CONTAINS(POINT('ICRS', ra, dec), "
                    f"BOX('ICRS', {ra_c:.8f}, {dec_c:.8f}, {w:.8f}, {h:.8f}))"
                )
                df = execute_tap_query(adql, retries=5, backoff=1.5, sleep_s=1.0)
                if not df.empty:
                    chunks.append(df)
                    total_stars += len(df)
                print(f"[{i}/{len(fits_files)}] OK: {len(df):,} estrellas | total: {total_stars:,}")
            except Exception as e:
                failures.append((fpath.name, str(e)))
                print(f"[{i}/{len(fits_files)}] FAIL: {fpath.name} - {e}")

    # Resultado final
    if console:
        console.print()  # línea en blanco

    if chunks:
        final_df = pd.concat(chunks, ignore_index=True)
        n_before_dedup = len(final_df)
        final_df.drop_duplicates("source_id", inplace=True)
        n_after = len(final_df)

        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path = out_path.with_suffix('.parquet')
        final_df.to_parquet(out_path, index=False)

        if console:
            console.print(
                f"[bold green]✅ Catálogo guardado[/bold green]\n"
                f"  Ruta:        [blue]{out_path}[/blue]\n"
                f"  Estrellas:   [yellow]{n_after:,}[/yellow] "
                f"(deduplicadas de {n_before_dedup:,})\n"
                f"  Tamaño:      [dim]{out_path.stat().st_size / 1e6:.1f} MB[/dim]"
            )
            if failures:
                console.print(
                    f"[yellow]⚠ {len(failures)} FITS fallaron[/yellow] "
                    f"(ver log: {args.log_file})"
                )
        else:
            print(f"✅ Éxito: {n_after:,} estrellas en {out_path}")
            if failures:
                print(f"⚠ {len(failures)} fallos (ver {args.log_file})")
    else:
        if console:
            console.print("[bold red]❌ No se obtuvieron resultados.[/bold red]")
        logging.warning("Ninguna query retornó estrellas.")


if __name__ == "__main__":
    main()
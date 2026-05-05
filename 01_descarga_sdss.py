"""
SCRIPT 01: SDSS MASS DOWNLOADER (v2 — Auditado)

Descarga masiva, paralela y resiliente de archivos FITS desde SDSS.
Implementa descompresión atómica (bz2 -> fits) y verificación opcional.

CAMBIOS RESPECTO A v1:
- FIX: timeout default subido de 15s a 60s (SDSS bajo carga puede demorar más).
- FIX: run y camcol ahora parametrizables (antes estaban hardcoded 001000/6).
- AÑADIDO: resumabilidad vía HTTP Range headers (si se cortó a mitad, retoma).
- AÑADIDO: backoff agresivo para errores de timeout (no solo 5xx).
- AÑADIDO: semáforo de concurrencia opcional para no abusar de SDSS (--rate).
"""

import argparse
import bz2
import logging
import os
import shutil
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from threading import Semaphore
from typing import Optional, Tuple

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

try:
    from astropy.io import fits
    HAS_ASTROPY = True
except ImportError:
    HAS_ASTROPY = False

try:
    from rich.progress import (
        BarColumn, Progress, TextColumn, TimeElapsedColumn, TimeRemainingColumn
    )
    HAS_RICH = True
except ImportError:
    HAS_RICH = False


def setup_logging(quiet: bool, log_file: Optional[str]) -> None:
    logger = logging.getLogger()
    logger.setLevel(logging.INFO if not quiet else logging.WARNING)
    if logger.hasHandlers():
        logger.handlers.clear()

    ch = logging.StreamHandler()
    ch.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(ch)

    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_path, encoding="utf-8")
        fh.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(fh)


def build_session(timeout: int, retries: int) -> requests.Session:
    sess = requests.Session()
    # Incluimos 408 (timeout) y 429 (rate limit) además de 5xx
    retry_strategy = Retry(
        total=retries,
        backoff_factor=1.5,  # agresivo: 1.5s, 3s, 6s, 12s...
        status_forcelist=[408, 429, 500, 502, 503, 504],
        allowed_methods=["GET", "HEAD"],
        raise_on_status=False,
    )
    adapter = HTTPAdapter(
        max_retries=retry_strategy,
        pool_connections=64,
        pool_maxsize=64
    )
    sess.mount("http://", adapter)
    sess.mount("https://", adapter)
    sess.request_timeout = timeout
    return sess


def _atomic_write(tmp_path: str, final_path: str) -> None:
    os.replace(tmp_path, final_path)


def decompress_bz2(bz2_path: str, fits_path: str) -> None:
    tmp_fits = fits_path + ".tmp"
    with bz2.BZ2File(bz2_path, "rb") as f_in, open(tmp_fits, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out, length=1024 * 1024)
    _atomic_write(tmp_fits, fits_path)


def verify_fits(fits_path: str) -> bool:
    if not HAS_ASTROPY:
        return True
    try:
        with fits.open(fits_path, memmap=False) as hdul:
            _ = hdul[0].header
            # Verificación más profunda: intentar acceder a los datos del primer HDU
            if hdul[0].data is None and len(hdul) > 1:
                _ = hdul[1].header
        return True
    except Exception as e:
        logging.warning(f"FITS corrupto {os.path.basename(fits_path)}: {e}")
        return False


def download_with_resume(
    session: requests.Session, url: str, dest_path: str, timeout: int
) -> None:
    """Descarga con soporte de resumado vía HTTP Range headers.

    Si ya existe `dest_path.tmp` parcial, retoma desde donde se cortó.
    """
    tmp_path = dest_path + ".tmp"
    resume_from = 0
    mode = "wb"
    headers = {}

    if os.path.exists(tmp_path):
        resume_from = os.path.getsize(tmp_path)
        if resume_from > 0:
            headers["Range"] = f"bytes={resume_from}-"
            mode = "ab"
            logging.debug(f"Retomando {os.path.basename(dest_path)} desde byte {resume_from}")

    with session.get(url, stream=True, timeout=timeout, headers=headers) as r:
        # Si el server no soporta Range, descarga desde cero
        if resume_from > 0 and r.status_code == 200:
            mode = "wb"
            logging.debug(f"Server ignoró Range; descargando {os.path.basename(dest_path)} completo.")
        r.raise_for_status()
        with open(tmp_path, mode) as f:
            for chunk in r.iter_content(chunk_size=256 * 1024):
                if chunk:
                    f.write(chunk)

    _atomic_write(tmp_path, dest_path)


def process_single_file(
    session: requests.Session,
    base_url: str,
    band: str,
    run: int,
    camcol: int,
    num: int,
    dest: str,
    verify: bool,
    keep_bz2: bool,
    timeout: int,
    semaphore: Optional[Semaphore],
) -> Tuple[str, bool, str]:

    bz2_name = f"frame-{band}-{run:06d}-{camcol}-{num:04d}.fits.bz2"
    fits_name = f"frame-{band}-{run:06d}-{camcol}-{num:04d}.fits"
    bz2_path = os.path.join(dest, bz2_name)
    fits_path = os.path.join(dest, fits_name)
    url = f"{base_url}/{bz2_name}"

    # Si ya existe el FITS final y verifica OK, skip
    if os.path.exists(fits_path):
        if verify and not verify_fits(fits_path):
            os.remove(fits_path)
        else:
            return fits_path, True, "skip"

    # Si tenemos bz2 completo pero no fits, descomprime y listo
    if not os.path.exists(fits_path) and os.path.exists(bz2_path):
        try:
            decompress_bz2(bz2_path, fits_path)
        except Exception:
            os.remove(bz2_path)

    # Si falta todo, descargar (con semáforo si aplica)
    if not os.path.exists(fits_path):
        if semaphore is not None:
            semaphore.acquire()
        try:
            download_with_resume(session, url, bz2_path, timeout)
            decompress_bz2(bz2_path, fits_path)
        except Exception as e:
            if os.path.exists(fits_path):
                os.remove(fits_path)
            return fits_path, False, str(e)
        finally:
            if semaphore is not None:
                semaphore.release()

    if verify and not verify_fits(fits_path):
        os.remove(fits_path)
        return fits_path, False, "verify_failed"

    if not keep_bz2 and os.path.exists(bz2_path):
        os.remove(bz2_path)

    return fits_path, True, "ok"


def main():
    p = argparse.ArgumentParser(description="Descargador concurrente de SDSS FITS (v2).")
    p.add_argument("--base-url", required=True)
    p.add_argument("--bands", nargs="+", default=["g", "i", "r", "u", "z"])
    p.add_argument("--start", type=int, required=True)
    p.add_argument("--end", type=int, required=True)
    p.add_argument("--run", type=int, default=1000,
                   help="SDSS run number (antes hardcoded 001000).")
    p.add_argument("--camcol", type=int, default=6,
                   help="SDSS camcol number (antes hardcoded 6).")
    p.add_argument("--dest", required=True)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--retries", type=int, default=5)
    p.add_argument("--timeout", type=int, default=60,
                   help="Timeout por request en segundos (antes era 15, demasiado agresivo).")
    p.add_argument("--rate", type=int, default=0,
                   help="Límite de descargas concurrentes (0 = sin límite; 4-6 recomendado para SDSS).")
    p.add_argument("--verify-fits", action="store_true")
    p.add_argument("--keep-bz2", action="store_true")
    p.add_argument("--quiet", action="store_true")
    p.add_argument("--log-file", default="sdss_download.log")
    args = p.parse_args()

    setup_logging(args.quiet, args.log_file)
    Path(args.dest).mkdir(parents=True, exist_ok=True)

    tasks = [(b, n) for b in args.bands for n in range(args.start, args.end + 1)]
    session = build_session(args.timeout, args.retries)
    semaphore = Semaphore(args.rate) if args.rate > 0 else None

    if semaphore:
        logging.info(f"Rate-limit activado: máx {args.rate} descargas concurrentes.")
    logging.info(f"Total tareas: {len(tasks)} | Run={args.run:06d} | Camcol={args.camcol}")

    failures = []
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {
            executor.submit(
                process_single_file, session, args.base_url, b, args.run, args.camcol,
                n, args.dest, args.verify_fits, args.keep_bz2, args.timeout, semaphore
            ): (b, n)
            for b, n in tasks
        }

        if HAS_RICH and not args.quiet:
            cols = [
                TextColumn("[cyan]{task.description}"),
                BarColumn(),
                TextColumn("{task.completed}/{task.total}"),
                TimeElapsedColumn(),
                TimeRemainingColumn(),
            ]
            with Progress(*cols) as prog:
                t = prog.add_task("Descargando...", total=len(tasks))
                for fut in as_completed(futures):
                    path, ok, msg = fut.result()
                    if not ok:
                        failures.append((futures[fut], msg))
                    prog.advance(t)
        else:
            for fut in as_completed(futures):
                path, ok, msg = fut.result()
                if not ok:
                    failures.append((futures[fut], msg))

    if failures:
        logging.warning(f"⚠️ {len(failures)} archivos fallaron. Muestra: {failures[:5]}")
    logging.info(f"✅ Descarga completa. Éxitos: {len(tasks) - len(failures)}/{len(tasks)}")


if __name__ == "__main__":
    main()

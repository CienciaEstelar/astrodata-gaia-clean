"""
SCRIPT 5: GENERADOR DE DIAGRAMA COLOR-MAGNITUD (CMD / HR) (v2 — Auditado)

Genera diagramas Hertzsprung-Russell (Color-Magnitud) a partir de catálogos
astronómicos enriquecidos. Soporta formatos nativos Parquet y CSV.

CAMBIOS RESPECTO A v1:
- AÑADIDO: modo `density` con hexbin + escala log (mucho mejor que scatter
  de 100k puntos rasterizados, que sobresatura y pierde estructura).
- AÑADIDO: overlay opcional de zonas canónicas del CMD (Main Sequence,
  Giant Branch, White Dwarf) para screening visual.
- FIX: muestreo aleatorio solo aplica en modo scatter; en hexbin usa TODOS
  los puntos (el hexbin ya reduce la carga visual).
- FIX: autodetect de columna de magnitud absoluta
  (`luminosidad_absoluta_g_corregida` o `luminosidad_absoluta_g`).
"""

import argparse
import logging
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ─────────────────────────────────────────────────────────────────────────────
# 1. LOGGING
# ─────────────────────────────────────────────────────────────────────────────

def setup_logging(verbose: bool = False) -> None:
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    if logger.hasHandlers():
        logger.handlers.clear()

    ch = logging.StreamHandler()
    ch.setFormatter(logging.Formatter("%(levelname)s - %(message)s"))
    logger.addHandler(ch)


# ─────────────────────────────────────────────────────────────────────────────
# 2. CARGA Y PREPARACIÓN
# ─────────────────────────────────────────────────────────────────────────────

def load_catalog(input_path: Path) -> pd.DataFrame:
    logging.info(f"Cargando datos desde: {input_path}")
    if not input_path.exists():
        raise FileNotFoundError(f"No se encontró: {input_path}")

    if input_path.suffix.lower() == '.parquet':
        df = pd.read_parquet(input_path)
    else:
        df = pd.read_csv(input_path, low_memory=False)

    logging.info(f"► {len(df):,} filas cargadas.")
    return df


def autodetect_mag_column(df: pd.DataFrame, user_col: str) -> str:
    """Detecta cuál columna de magnitud absoluta usar."""
    candidatos = [
        user_col,
        "luminosidad_absoluta_g_corregida",
        "luminosidad_absoluta_g",
        "luminosidad_absoluta_g_aparente",
    ]
    for c in candidatos:
        if c in df.columns:
            if c != user_col:
                logging.info(f"Magnitud detectada automáticamente: '{c}' (user pidió '{user_col}').")
            return c
    return user_col  # se validará como error después


def prepare_plot_data(
    df: pd.DataFrame, col_color: str, col_mag: str,
    limit_points: int, mode: str
) -> pd.DataFrame:
    if col_color not in df.columns:
        raise ValueError(f"Columna de color '{col_color}' no existe.")
    if col_mag not in df.columns:
        raise ValueError(f"Columna de magnitud '{col_mag}' no existe.")

    df_plot = df[[col_color, col_mag]].copy()

    initial_len = len(df_plot)
    df_plot.dropna(subset=[col_color, col_mag], inplace=True)
    final_len = len(df_plot)

    if initial_len - final_len > 0:
        logging.info(f"► Omitidas {initial_len - final_len:,} filas con NaN.")
    if final_len == 0:
        raise ValueError("No quedan estrellas válidas tras la limpieza.")

    # Muestreo aleatorio SOLO en modo scatter. En hexbin se usan todos.
    if mode == "scatter" and limit_points > 0 and final_len > limit_points:
        logging.info(f"Scatter mode: reduciendo de {final_len:,} a {limit_points:,} puntos.")
        df_plot = df_plot.sample(n=limit_points, random_state=42)
    elif mode == "hexbin":
        logging.info(f"Hexbin mode: usando {final_len:,} puntos completos.")

    return df_plot


# ─────────────────────────────────────────────────────────────────────────────
# 3. VISUALIZACIÓN
# ─────────────────────────────────────────────────────────────────────────────

def create_and_save_cmd(df: pd.DataFrame, args: argparse.Namespace) -> None:
    logging.info(f"Renderizando CMD (modo: {args.mode})...")

    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(10, 12))

    if args.mode == "hexbin":
        # Hexbin con escala log: revela estructura en zonas densas sin saturar
        hb = ax.hexbin(
            df[args.col_color].values,
            df[args.col_mag].values,
            gridsize=args.gridsize,
            bins='log',
            cmap='plasma',
            mincnt=1,
        )
        cb = fig.colorbar(hb, ax=ax, label='log₁₀(N estrellas)', shrink=0.7)
    else:  # scatter
        ax.scatter(
            df[args.col_color],
            df[args.col_mag],
            s=args.marker_size,
            alpha=args.alpha,
            edgecolors='none',
            rasterized=True
        )

    ax.set_xlabel(f"Color BP-RP ({args.col_color})", fontsize=12)
    ax.set_ylabel(f"Magnitud Absoluta ({args.col_mag})", fontsize=12)
    ax.set_title(args.title, fontsize=16, pad=15)

    # Convención: magnitudes más brillantes hacia arriba
    ax.invert_yaxis()

    # Límites
    if args.xlim:
        try:
            xmin, xmax = map(float, args.xlim.split(','))
            ax.set_xlim(xmin, xmax)
        except ValueError:
            logging.warning(f"--xlim inválido: '{args.xlim}'.")
    if args.ylim:
        try:
            ymin, ymax = map(float, args.ylim.split(','))
            ax.set_ylim(ymin, ymax)
        except ValueError:
            logging.warning(f"--ylim inválido: '{args.ylim}'.")

    # Overlay de zonas canónicas (opcional)
    if args.overlay_regions:
        add_canonical_regions(ax)

    ax.grid(True, linestyle='--', alpha=0.5)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    logging.info(f"Guardando en: {out_path}")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    logging.info("► Exportado.")

    if args.show_plot:
        plt.show()
    plt.close(fig)


def add_canonical_regions(ax) -> None:
    """Superpone zonas aproximadas del CMD Gaia (orientación, no astrofísicamente preciso)."""
    # Main Sequence: banda aproximada M_G = 7*(BP-RP) - 1 ± 2
    bp_rp = np.linspace(-0.5, 4.0, 100)
    ms_center = 7 * bp_rp - 1
    ax.fill_between(
        bp_rp, ms_center - 2, ms_center + 2,
        alpha=0.08, color='blue', label='Banda Secuencia Principal (aprox)'
    )
    # Rama Gigante: M_G ≈ 1 a -3, BP-RP > 1
    ax.fill_between(
        [1.0, 3.5], [-3, -3], [1, 1],
        alpha=0.08, color='red', label='Rama Gigante (aprox)'
    )
    # Enanas Blancas: M_G > 10 y BP-RP < 1
    ax.fill_between(
        [-0.5, 1.0], [10, 10], [16, 16],
        alpha=0.08, color='gray', label='Enanas Blancas (aprox)'
    )
    ax.legend(loc='lower left', fontsize=9, framealpha=0.85)


# ─────────────────────────────────────────────────────────────────────────────
# 4. CLI
# ─────────────────────────────────────────────────────────────────────────────

def get_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Generador de Diagramas Color-Magnitud (v2 — hexbin mode).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    g_io = p.add_argument_group("I/O")
    g_io.add_argument('--input', type=str, required=True)
    g_io.add_argument('--output', type=str, required=True)

    g_data = p.add_argument_group("Mapeo de Datos")
    g_data.add_argument('--col_color', type=str, default='color_bp_rp')
    g_data.add_argument('--col_mag', type=str, default='luminosidad_absoluta_g_corregida')
    g_data.add_argument('--limit_points', type=int, default=100000,
                        help='Límite en modo scatter; ignorado en modo hexbin.')

    g_plot = p.add_argument_group("Configuración Gráfica")
    g_plot.add_argument('--mode', choices=['scatter', 'hexbin'], default='hexbin',
                        help="hexbin es MUCHO mejor para catálogos > 50k.")
    g_plot.add_argument('--gridsize', type=int, default=120,
                        help="Resolución hexbin (mayor = más detalle).")
    g_plot.add_argument('--title', type=str, default='Diagrama Color-Magnitud (CMD)')
    g_plot.add_argument('--xlim', type=str, default=None)
    g_plot.add_argument('--ylim', type=str, default=None)
    g_plot.add_argument('--marker_size', type=float, default=0.5)
    g_plot.add_argument('--alpha', type=float, default=0.5)
    g_plot.add_argument('--overlay_regions', action='store_true',
                        help="Superponer zonas aproximadas (MS, Gigantes, EB).")
    g_plot.add_argument('--show_plot', action='store_true')

    g_sys = p.add_argument_group("Sistema")
    g_sys.add_argument('-v', '--verbose', action='store_true')

    return p


def main() -> None:
    parser = get_parser()
    args = parser.parse_args()

    setup_logging(verbose=args.verbose)
    logging.info("=" * 80)
    logging.info("🎨 GENERADOR DE CMD/HR (v2)")
    logging.info("=" * 80)

    try:
        input_path = Path(args.input)
        df_raw = load_catalog(input_path)

        # Auto-detect columna de magnitud
        args.col_mag = autodetect_mag_column(df_raw, args.col_mag)

        df_plot = prepare_plot_data(
            df_raw, args.col_color, args.col_mag, args.limit_points, args.mode
        )
        create_and_save_cmd(df_plot, args)
        logging.info("\n✅ PROCESO COMPLETADO")
    except KeyboardInterrupt:
        logging.warning("\n⚠️ Interrumpido por usuario.")
    except Exception as e:
        logging.error(f"\n❌ Error: {e}", exc_info=args.verbose)


if __name__ == "__main__":
    main()

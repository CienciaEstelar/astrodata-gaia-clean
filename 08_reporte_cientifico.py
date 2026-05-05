"""
SCRIPT 8: REPORTE CIENTÍFICO AUDITABLE (PDF multi-página)
==========================================================

Genera un reporte PDF profesional consolidando todos los resultados del
pipeline Gaia: catálogo depurado, cinemática UVW, clasificación, y
opcionalmente trayectorias de la simulación N-body del Script 07.

ESTRUCTURA DEL PDF:
-------------------
Página 1 — Portada con metadata del run (fecha, hash input, parámetros)
Página 2 — Estadísticas del catálogo (counts, coberturas, dists básicas)
Página 3 — Diagrama Color-Magnitud (CMD/HR) con zonas canónicas
Página 4 — Mapa galáctico (proyección Aitoff en coords l,b)
Página 5 — Distribuciones de velocidades UVW (3 paneles + medianas)
Página 6 — Test estadístico: KS vs Gaussiana LSR (rigor publicable)
Página 7 — Tabla de clasificación + gráfico de barras
Página 8 — Top candidatas hypervelocity (si hay datos de V_3D)
Página 9 — Simulación N-body: trayectorias + deriva de energía (opcional)
Página 10 — Apéndice: filtros aplicados, versiones, reproducibilidad

PRINCIPIOS DE DISEÑO:
- Rasterización inteligente: scatter plots densos usan rasterized=True
  para mantener PDFs < 10 MB aún con 100k estrellas.
- Estilo científico publicable: fuentes STIX, fondo blanco, grids sutiles.
- Reproducibilidad: hash SHA256 del input + versiones de librerías.
- Fail-loud: si falta una columna, dice qué y por qué.

INPUTS:
  --catalog       Parquet procesado (Script 03 o 06)
  --simulation    HDF5 del Script 07 (opcional)
  --output        Ruta del PDF final

PATRONES USADOS (de gaia_core):
- plot_context() — aplica estilo científico y lo restaura al salir
- DataQualityError / ContractViolationError — fail loud si falta algo
- validators físicos para sanity check antes de plotear
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import platform
import sys
import time
import warnings
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional, Tuple, Dict, Any, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
from scipy import stats as scipy_stats

from gaia_core import (
    GaiaPipelineError, DataQualityError, ContractViolationError,
    plot_context, setup_logging, HAS_H5PY,
)

if HAS_H5PY:
    import h5py

try:
    from rich.console import Console
    from rich.progress import Progress, BarColumn, TextColumn
    HAS_RICH = True
except ImportError:
    HAS_RICH = False


# ═════════════════════════════════════════════════════════════════════════════
# CONSTANTES FÍSICAS (Schönrich, Binney & Dehnen 2010 — valores canónicos LSR)
# ═════════════════════════════════════════════════════════════════════════════

# Local Standard of Rest (LSR) — movimiento peculiar del Sol respecto al LSR
# Schönrich+2010, MNRAS 403, 1829
SUN_UVW_LSR = {
    "U": 11.1,   # km/s, hacia el centro galáctico
    "V": 12.24,  # km/s, en dirección de rotación galáctica
    "W": 7.25,   # km/s, hacia el polo galáctico norte
}

# Dispersiones típicas del disco delgado para población de G/K dwarfs
# (Bensby+2003, Holmberg+2009)
DISK_DISPERSIONS = {
    "sigma_U": 35.0,  # km/s
    "sigma_V": 25.0,
    "sigma_W": 20.0,
}


# ═════════════════════════════════════════════════════════════════════════════
# DATACLASSES
# ═════════════════════════════════════════════════════════════════════════════

@dataclass
class ReportConfig:
    catalog_path: str
    simulation_path: Optional[str] = None
    output_path: str = "reporte_gaia.pdf"
    include_simulation: bool = False
    n_trajectory_samples: int = 30  # cuántas trayectorias plotear
    top_hypervelocity: int = 10      # cuántas HVS candidatas mostrar


@dataclass
class CatalogStats:
    """Estadísticas resumidas del catálogo para la portada."""
    n_stars: int
    n_cols: int
    coverage_vrad: float
    coverage_v3d: float
    coverage_uvw: float
    median_distance_pc: float
    median_M_G: float
    type_counts: Dict[str, int] = field(default_factory=dict)


# ═════════════════════════════════════════════════════════════════════════════
# HELPERS DE DATOS
# ═════════════════════════════════════════════════════════════════════════════

def compute_sha256(path: Path, chunk_size: int = 1 << 20) -> str:
    """Hash SHA256 del archivo para trazabilidad de reproducibilidad."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def load_catalog(path: Path) -> pd.DataFrame:
    """Carga el catálogo Parquet o CSV validando contratos mínimos."""
    if not path.exists():
        raise DataQualityError(f"No se encuentra el catálogo: {path}")

    if path.suffix == ".parquet":
        df = pd.read_parquet(path)
    elif path.suffix == ".csv":
        df = pd.read_csv(path, low_memory=False)
    else:
        raise DataQualityError(f"Formato no soportado: {path.suffix}")

    # Columnas mínimas para el reporte
    required = ["ra", "dec", "distancia_parsecs"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ContractViolationError(
            f"Catálogo no tiene columnas mínimas: {missing}. "
            f"Columnas disponibles: {list(df.columns)[:20]}..."
        )

    logging.info(f"Catálogo cargado: {len(df):,} estrellas × {len(df.columns)} columnas")
    return df


def compute_catalog_stats(df: pd.DataFrame) -> CatalogStats:
    """Calcula KPIs del catálogo para la portada."""
    n = len(df)

    coverage_vrad = float(df["V_radial"].notna().mean()) if "V_radial" in df.columns else 0.0
    coverage_v3d = float(df["velocidad_total_3d_kms"].notna().mean()) if "velocidad_total_3d_kms" in df.columns else 0.0

    has_uvw = all(c in df.columns for c in ["U_kms", "V_kms", "W_kms"])
    coverage_uvw = float(
        (df["U_kms"].notna() & df["V_kms"].notna() & df["W_kms"].notna()).mean()
    ) if has_uvw else 0.0

    median_dist = float(df["distancia_parsecs"].median())
    median_MG = float(df["luminosidad_absoluta_g"].median()) if "luminosidad_absoluta_g" in df.columns else np.nan

    type_counts = {}
    if "tipo_objeto" in df.columns:
        type_counts = df["tipo_objeto"].value_counts().to_dict()

    return CatalogStats(
        n_stars=n,
        n_cols=len(df.columns),
        coverage_vrad=coverage_vrad,
        coverage_v3d=coverage_v3d,
        coverage_uvw=coverage_uvw,
        median_distance_pc=median_dist,
        median_M_G=median_MG,
        type_counts=type_counts,
    )


# ═════════════════════════════════════════════════════════════════════════════
# PÁGINAS DEL PDF
# ═════════════════════════════════════════════════════════════════════════════

def page_cover(pdf: PdfPages, stats: CatalogStats,
               catalog_path: Path, sim_path: Optional[Path]) -> None:
    """Página 1: portada con metadata del run."""
    fig = plt.figure(figsize=(8.27, 11.69))  # A4
    fig.patch.set_facecolor('white')

    # Título
    fig.text(0.5, 0.92, "Reporte Científico — Proyecto Gaia",
             fontsize=22, weight="bold", ha="center", color="#1a1a2e")
    fig.text(0.5, 0.885, "Pipeline de Análisis Cinemático y Dinámico",
             fontsize=14, ha="center", style="italic", color="#16213e")

    # Separador
    fig.add_artist(plt.Line2D([0.15, 0.85], [0.855, 0.855],
                              color="#16213e", linewidth=1))

    # Metadata del run
    y = 0.81
    fig.text(0.15, y, "Fecha:", fontsize=10, weight="bold")
    fig.text(0.35, y, datetime.now().strftime("%Y-%m-%d %H:%M:%S"), fontsize=10)
    y -= 0.025

    fig.text(0.15, y, "Catálogo:", fontsize=10, weight="bold")
    fig.text(0.35, y, str(catalog_path.name), fontsize=9, family="monospace")
    y -= 0.025

    try:
        sha = compute_sha256(catalog_path)[:16]
        fig.text(0.15, y, "Hash SHA256:", fontsize=10, weight="bold")
        fig.text(0.35, y, f"{sha}...", fontsize=9, family="monospace")
        y -= 0.025
    except Exception:
        pass

    if sim_path and sim_path.exists():
        fig.text(0.15, y, "Simulación N-body:", fontsize=10, weight="bold")
        fig.text(0.35, y, str(sim_path.name), fontsize=9, family="monospace")
        y -= 0.025

    fig.text(0.15, y, "Sistema:", fontsize=10, weight="bold")
    fig.text(0.35, y, f"{platform.python_implementation()} {platform.python_version()} / {platform.system()}",
             fontsize=9, family="monospace")
    y -= 0.04

    # KPIs principales del catálogo
    fig.text(0.15, y, "Estadísticas del catálogo", fontsize=13,
             weight="bold", color="#1a1a2e")
    y -= 0.03

    kpis = [
        ("Estrellas totales", f"{stats.n_stars:,}"),
        ("Columnas", f"{stats.n_cols}"),
        ("Cobertura V_radial", f"{stats.coverage_vrad:.1%}"),
        ("Cobertura V_3D", f"{stats.coverage_v3d:.1%}"),
        ("Cobertura UVW", f"{stats.coverage_uvw:.1%}"),
        ("Distancia mediana", f"{stats.median_distance_pc:.0f} pc"),
        ("M_G mediana", f"{stats.median_M_G:.2f} mag"),
    ]
    for label, value in kpis:
        fig.text(0.18, y, label + ":", fontsize=10)
        fig.text(0.48, y, value, fontsize=10, family="monospace",
                 weight="bold", color="#0066cc")
        y -= 0.022

    # Footer
    fig.text(0.5, 0.05, "Generado por pipeline Gaia — Script 08",
             fontsize=9, ha="center", style="italic", color="#666")
    fig.text(0.5, 0.03, "Reproducibilidad científica garantizada por gaia_core",
             fontsize=8, ha="center", style="italic", color="#888")

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def page_catalog_stats(pdf: PdfPages, df: pd.DataFrame, stats: CatalogStats) -> None:
    """Página 2: estadísticas detalladas del catálogo."""
    fig, axes = plt.subplots(2, 2, figsize=(11.69, 8.27))
    fig.suptitle("Estadísticas del catálogo depurado", fontsize=14, weight="bold")

    # (1) Histograma de distancias
    ax = axes[0, 0]
    dist = df["distancia_parsecs"].dropna()
    ax.hist(dist, bins=60, color="#0066cc", alpha=0.7, edgecolor="black", linewidth=0.3)
    ax.axvline(np.median(dist), color="red", ls="--", lw=1.5,
               label=f"Mediana: {np.median(dist):.0f} pc")
    ax.set_xlabel("Distancia (pc)")
    ax.set_ylabel("N estrellas")
    ax.set_title("Distribución de distancias")
    ax.legend()

    # (2) Histograma de magnitudes absolutas (si existen)
    ax = axes[0, 1]
    if "luminosidad_absoluta_g" in df.columns:
        MG = df["luminosidad_absoluta_g"].dropna()
        ax.hist(MG, bins=60, color="#e66d00", alpha=0.7, edgecolor="black", linewidth=0.3)
        ax.axvline(np.median(MG), color="red", ls="--", lw=1.5,
                   label=f"Mediana: {np.median(MG):.2f}")
        ax.invert_xaxis()
        ax.set_xlabel("M_G (mag)")
        ax.set_ylabel("N estrellas")
        ax.set_title("Distribución de magnitudes absolutas")
        ax.legend()
    else:
        ax.text(0.5, 0.5, "M_G no disponible", ha="center", va="center", transform=ax.transAxes)

    # (3) Color BP-RP
    ax = axes[1, 0]
    if "color_bp_rp" in df.columns:
        c = df["color_bp_rp"].dropna()
        ax.hist(c, bins=60, color="#c00080", alpha=0.7, edgecolor="black", linewidth=0.3)
        ax.axvline(np.median(c), color="black", ls="--", lw=1.5,
                   label=f"Mediana: {np.median(c):.2f}")
        ax.set_xlabel("Color BP–RP")
        ax.set_ylabel("N estrellas")
        ax.set_title("Distribución de colores")
        ax.legend()

    # (4) Resumen en texto
    ax = axes[1, 1]
    ax.axis("off")
    lines = [
        f"TOTAL: {stats.n_stars:,} estrellas",
        f"",
        f"Cobertura por tipo de dato:",
        f"  • V_radial (Gaia): {stats.coverage_vrad:.1%}",
        f"  • V_3D completa:    {stats.coverage_v3d:.1%}",
        f"  • UVW galactocent.: {stats.coverage_uvw:.1%}",
        f"",
        f"Distancia mediana: {stats.median_distance_pc:.0f} pc",
        f"M_G mediana:       {stats.median_M_G:.2f}",
        f"",
        f"Tipos de objeto (top 5):",
    ]
    if stats.type_counts:
        top5 = sorted(stats.type_counts.items(), key=lambda x: -x[1])[:5]
        for tipo, n in top5:
            pct = 100 * n / stats.n_stars
            lines.append(f"  • {tipo}: {n:,} ({pct:.1f}%)")

    ax.text(0.05, 0.95, "\n".join(lines), fontsize=10, family="monospace",
            transform=ax.transAxes, verticalalignment="top")

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def page_cmd(pdf: PdfPages, df: pd.DataFrame) -> None:
    """Página 3: Diagrama Color-Magnitud (CMD)."""
    if "color_bp_rp" not in df.columns or "luminosidad_absoluta_g" not in df.columns:
        logging.warning("Sin BP-RP o M_G — saltando página CMD")
        return

    data = df[["color_bp_rp", "luminosidad_absoluta_g"]].dropna()
    if len(data) == 0:
        return

    fig, ax = plt.subplots(figsize=(8.27, 9.5))
    ax.hexbin(data["color_bp_rp"], data["luminosidad_absoluta_g"],
              gridsize=120, mincnt=1, bins="log", cmap="viridis")
    ax.invert_yaxis()
    ax.set_xlabel("Color BP − RP (mag)")
    ax.set_ylabel("Magnitud absoluta $M_G$ (mag)")
    ax.set_title(f"Diagrama Color-Magnitud — {len(data):,} estrellas\n"
                 f"Hexbin densidad log, Gaia DR3",
                 fontsize=13, weight="bold")

    # Overlay zonas canónicas
    ax.text(0.6, 4.5, "Secuencia\nprincipal", fontsize=10,
            ha="center", color="white", weight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="black", alpha=0.5))
    ax.text(2.5, 0, "Rama de\ngigantes", fontsize=10,
            ha="center", color="white", weight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="darkred", alpha=0.6))
    ax.text(0.0, 13, "Enanas\nblancas", fontsize=10,
            ha="center", color="white", weight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="navy", alpha=0.6))

    ax.set_xlim(-0.5, 4.5)
    ax.set_ylim(17, -5)

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def page_galactic_map(pdf: PdfPages, df: pd.DataFrame) -> None:
    """Página 4: mapa galáctico en proyección Aitoff."""
    if "l_galactic" not in df.columns or "b_galactic" not in df.columns:
        logging.warning("Sin l,b galácticas — saltando mapa Aitoff")
        return

    data = df[["l_galactic", "b_galactic"]].dropna()
    if len(data) == 0:
        return

    # Convertir a radianes y centrar longitudes en 0 (convención Aitoff)
    l_rad = np.deg2rad(((data["l_galactic"] + 180) % 360) - 180)
    b_rad = np.deg2rad(data["b_galactic"])

    fig = plt.figure(figsize=(11.69, 6.5))
    ax = fig.add_subplot(111, projection="aitoff")
    ax.scatter(l_rad, b_rad, s=0.5, alpha=0.3, c="#0066cc", rasterized=True)
    ax.grid(True, alpha=0.3)
    ax.set_title(f"Mapa galáctico (Aitoff) — {len(data):,} estrellas\n"
                 f"Coordenadas galácticas (l, b)",
                 fontsize=13, weight="bold", pad=20)
    ax.set_xlabel("Longitud galáctica l")
    ax.set_ylabel("Latitud galáctica b")

    pdf.savefig(fig, bbox_inches="tight", dpi=200)
    plt.close(fig)


def page_uvw_histograms(pdf: PdfPages, df: pd.DataFrame) -> None:
    """Página 5: distribuciones de velocidades UVW con overlay LSR."""
    cols = ["U_kms", "V_kms", "W_kms"]
    if not all(c in df.columns for c in cols):
        logging.warning("Sin UVW — saltando histogramas de velocidades")
        return

    data = df[cols].dropna()
    if len(data) == 0:
        return

    fig, axes = plt.subplots(3, 1, figsize=(8.27, 11.69), sharex=False)
    fig.suptitle(f"Distribuciones UVW galactocéntricas — {len(data):,} estrellas",
                 fontsize=14, weight="bold")

    labels = ["$U$ (hacia centro galáctico)", "$V$ (rotación galáctica)",
              "$W$ (hacia polo galáctico)"]
    colors = ["#0066cc", "#e66d00", "#c00080"]
    sun_vals = [SUN_UVW_LSR["U"], SUN_UVW_LSR["V"], SUN_UVW_LSR["W"]]
    sigma_disk = [DISK_DISPERSIONS["sigma_U"], DISK_DISPERSIONS["sigma_V"],
                  DISK_DISPERSIONS["sigma_W"]]

    for i, (col, label, color, sun_v, sigma) in enumerate(
        zip(cols, labels, colors, sun_vals, sigma_disk)
    ):
        ax = axes[i]
        v = data[col].values
        v_clipped = np.clip(v, -500, 500)  # limitar para histograma visible

        ax.hist(v_clipped, bins=80, color=color, alpha=0.7, density=True,
                edgecolor="black", linewidth=0.3, label=f"Datos (N={len(v)})")

        # Overlay Gaussiana LSR
        x_grid = np.linspace(-500, 500, 1000)
        gauss = scipy_stats.norm.pdf(x_grid, loc=sun_v, scale=sigma)
        ax.plot(x_grid, gauss, "k--", lw=2, alpha=0.8,
                label=f"Gaussiana LSR (Schönrich+2010)\n$\\mu={sun_v}$, $\\sigma={sigma}$")

        ax.axvline(np.median(v), color="red", ls=":", lw=1.5,
                   label=f"Mediana datos: {np.median(v):.1f}")
        ax.set_xlabel(f"{label} (km/s)")
        ax.set_ylabel("Densidad")
        ax.legend(loc="upper right", fontsize=8)
        ax.set_xlim(-300, 300)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def page_ks_test(pdf: PdfPages, df: pd.DataFrame) -> None:
    """Página 6: test KS para validar distribuciones UVW."""
    cols = ["U_kms", "V_kms", "W_kms"]
    if not all(c in df.columns for c in cols):
        return

    data = df[cols].dropna()
    if len(data) < 30:
        return

    fig, ax = plt.subplots(figsize=(8.27, 11.69))
    ax.axis("off")

    ax.text(0.5, 0.96, "Test Kolmogorov-Smirnov: UVW vs Gaussiana LSR",
            fontsize=14, weight="bold", ha="center", transform=ax.transAxes)
    ax.text(0.5, 0.92, "Hipótesis nula: la población sigue una Gaussiana con los\n"
                        "parámetros canónicos del disco galáctico (Schönrich+2010)",
            fontsize=10, ha="center", style="italic", transform=ax.transAxes)

    y = 0.82
    sun_vals = [SUN_UVW_LSR["U"], SUN_UVW_LSR["V"], SUN_UVW_LSR["W"]]
    sigma_disk = [DISK_DISPERSIONS["sigma_U"], DISK_DISPERSIONS["sigma_V"],
                  DISK_DISPERSIONS["sigma_W"]]

    results = []
    for col, mu, sigma in zip(cols, sun_vals, sigma_disk):
        v = data[col].values
        # Standardizar: si H0 es cierta, los datos están distribuidos N(mu, sigma)
        stat, pvalue = scipy_stats.kstest(v, "norm", args=(mu, sigma))
        results.append((col, mu, sigma, stat, pvalue))

    # Tabla de resultados
    header = f"{'Componente':<12} {'μ_LSR':>8} {'σ_disk':>8} {'KS_stat':>10} {'p-value':>14} {'Conclusión':<25}"
    ax.text(0.05, y, header, fontsize=10, family="monospace", weight="bold",
            transform=ax.transAxes)
    y -= 0.02
    ax.text(0.05, y, "─" * 80, fontsize=9, family="monospace", transform=ax.transAxes)
    y -= 0.03

    for col, mu, sigma, stat, pvalue in results:
        conclusion = "Rechaza H0 (α=0.05)" if pvalue < 0.05 else "No rechaza H0"
        row = f"{col:<12} {mu:>8.2f} {sigma:>8.1f} {stat:>10.4f} {pvalue:>14.3e} {conclusion:<25}"
        color = "#c00000" if pvalue < 0.05 else "#006600"
        ax.text(0.05, y, row, fontsize=10, family="monospace",
                color=color, transform=ax.transAxes)
        y -= 0.025

    y -= 0.03
    ax.text(0.05, y,
            "Interpretación física:",
            fontsize=11, weight="bold", transform=ax.transAxes)
    y -= 0.025

    interpretation = (
        "• En una muestra LOCAL del disco (vecindario solar) es normal que\n"
        "  los tests KS rechacen la H0 por tamaño grande de muestra + subpoblaciones.\n"
        "• Presencia de disco grueso, halo y moving groups causa desvíos.\n"
        "• Lo importante es que las MEDIANAS y DISPERSIONES estén en rango físico\n"
        "  (ver página anterior), no que la distribución sea Gaussiana perfecta.\n"
        "• Una desviación sistemática en μ puede indicar contaminación cinemática\n"
        "  (ej: pertenencia a un stream estelar, halo residual)."
    )
    ax.text(0.05, y, interpretation, fontsize=9, transform=ax.transAxes,
            verticalalignment="top")

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def page_classification(pdf: PdfPages, stats: CatalogStats) -> None:
    """Página 7: tabla y gráfico de clasificación."""
    if not stats.type_counts:
        return

    items = sorted(stats.type_counts.items(), key=lambda x: -x[1])
    types = [t for t, _ in items]
    counts = [c for _, c in items]
    pcts = [100 * c / stats.n_stars for c in counts]

    fig, (ax_bar, ax_table) = plt.subplots(2, 1, figsize=(8.27, 11.69),
                                            gridspec_kw={"height_ratios": [1.2, 1]})

    # Gráfico de barras horizontal
    y_pos = np.arange(len(types))
    bars = ax_bar.barh(y_pos, counts, color="#0066cc", alpha=0.8, edgecolor="black")
    ax_bar.set_yticks(y_pos)
    ax_bar.set_yticklabels(types, fontsize=10)
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("Número de estrellas")
    ax_bar.set_title("Clasificación de objetos", fontsize=13, weight="bold")

    # Añadir porcentaje al final de cada barra
    for bar, pct, cnt in zip(bars, pcts, counts):
        ax_bar.text(bar.get_width() + max(counts) * 0.01, bar.get_y() + bar.get_height() / 2,
                    f"{cnt:,} ({pct:.1f}%)", va="center", fontsize=9)

    # Tabla
    ax_table.axis("off")
    table_text = f"{'Tipo':<30} {'N':>10} {'%':>8}\n"
    table_text += "─" * 50 + "\n"
    for t, c, p in zip(types, counts, pcts):
        table_text += f"{t:<30} {c:>10,} {p:>7.2f}%\n"
    table_text += "─" * 50 + "\n"
    table_text += f"{'TOTAL':<30} {stats.n_stars:>10,} {100.0:>7.2f}%\n"

    ax_table.text(0.05, 0.95, table_text, fontsize=10, family="monospace",
                  verticalalignment="top", transform=ax_table.transAxes)

    plt.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def page_hypervelocity(pdf: PdfPages, df: pd.DataFrame, top_n: int = 10) -> None:
    """Página 8: top candidatas hypervelocity stars (HVS)."""
    if "velocidad_total_3d_kms" not in df.columns:
        return

    data = df.dropna(subset=["velocidad_total_3d_kms"]).copy()
    if len(data) == 0:
        return

    # Top N estrellas con mayor V_3D
    top = data.nlargest(top_n, "velocidad_total_3d_kms")

    fig, (ax_hist, ax_table) = plt.subplots(2, 1, figsize=(8.27, 11.69),
                                              gridspec_kw={"height_ratios": [1, 1.2]})

    # Histograma de V_3D con marcador de HVS
    v_clip = data["velocidad_total_3d_kms"].clip(upper=1500)
    ax_hist.hist(v_clip, bins=80, color="#0066cc", alpha=0.7, edgecolor="black", linewidth=0.3)
    ax_hist.axvline(500, color="orange", ls="--", lw=1.5,
                    label="Fast stars (500 km/s)")
    ax_hist.axvline(1000, color="red", ls="--", lw=1.5,
                    label="Hypervelocity stars (1000 km/s)")
    ax_hist.set_xlabel("$V_{3D}$ (km/s)")
    ax_hist.set_ylabel("N estrellas")
    ax_hist.set_title(f"Distribución de velocidades 3D — {len(data):,} estrellas",
                       fontsize=12, weight="bold")
    ax_hist.set_yscale("log")
    ax_hist.legend()

    # Tabla de top N
    ax_table.axis("off")
    ax_table.text(0.5, 0.98, f"Top {top_n} candidatas a Hypervelocity Stars (HVS)",
                  fontsize=13, weight="bold", ha="center", transform=ax_table.transAxes)

    header = f"{'#':<3} {'source_id':<22} {'RA':>8} {'Dec':>8} {'Dist (pc)':>10} {'V_3D (km/s)':>14}"
    header_y = 0.90
    ax_table.text(0.03, header_y, header, fontsize=9, family="monospace",
                  weight="bold", transform=ax_table.transAxes)
    ax_table.text(0.03, header_y - 0.022, "─" * 75, fontsize=9, family="monospace",
                  transform=ax_table.transAxes)

    y = header_y - 0.04
    for i, (_, row) in enumerate(top.iterrows(), 1):
        sid = str(row.get("source_id", "?"))[:20]
        ra = row.get("ra", np.nan)
        dec = row.get("dec", np.nan)
        dist = row.get("distancia_parsecs", np.nan)
        v3d = row.get("velocidad_total_3d_kms", np.nan)
        line = f"{i:<3} {sid:<22} {ra:>8.3f} {dec:>8.3f} {dist:>10.1f} {v3d:>14.1f}"
        color = "#c00000" if v3d > 1000 else ("#e66d00" if v3d > 500 else "black")
        ax_table.text(0.03, y, line, fontsize=9, family="monospace",
                      color=color, transform=ax_table.transAxes)
        y -= 0.025

    ax_table.text(0.03, y - 0.02,
                   "Nota: HVS reales requieren verificación con velocidades radiales\n"
                   "de alta precisión y cross-match con surveys espectroscópicos\n"
                   "(LAMOST, SDSS/APOGEE, Gaia RVS).",
                   fontsize=8, style="italic", color="#666",
                   transform=ax_table.transAxes, verticalalignment="top")

    plt.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def page_simulation(pdf: PdfPages, sim_path: Path, n_samples: int = 30) -> None:
    """Página 9: trayectorias de la simulación N-body."""
    if not HAS_H5PY:
        logging.warning("h5py no instalado, saltando página de simulación")
        return
    if not sim_path.exists():
        logging.warning(f"HDF5 no encontrado: {sim_path}")
        return

    with h5py.File(sim_path, "r") as f:
        positions = f["positions"][:]  # (n_snapshots, n_stars, 3) en pc
        times = f["times_years"][:]
        meta = dict(f.attrs)

    n_snap, n_stars, _ = positions.shape

    fig = plt.figure(figsize=(11.69, 8.27))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

    # Panel 1: XY del disco
    ax_xy = fig.add_subplot(gs[0, 0])
    # Pickear n_samples estrellas aleatorias
    rng = np.random.default_rng(42)
    sample_idx = rng.choice(n_stars, size=min(n_samples, n_stars), replace=False)
    for i in sample_idx:
        ax_xy.plot(positions[:, i, 0], positions[:, i, 1],
                   "-", alpha=0.4, lw=0.5, rasterized=True)
    # Marcar posiciones iniciales y finales
    ax_xy.scatter(positions[0, sample_idx, 0], positions[0, sample_idx, 1],
                  s=15, c="green", marker="o", zorder=5, label="t=0", edgecolors="black", linewidths=0.3)
    ax_xy.scatter(positions[-1, sample_idx, 0], positions[-1, sample_idx, 1],
                  s=15, c="red", marker="X", zorder=5, label=f"t={times[-1]/1e6:.2f} Myr",
                  edgecolors="black", linewidths=0.3)
    ax_xy.set_xlabel("X galáctico (pc)")
    ax_xy.set_ylabel("Y galáctico (pc)")
    ax_xy.set_title(f"Proyección XY del disco (n={n_samples} estrellas)",
                     fontsize=11, weight="bold")
    ax_xy.legend(loc="upper right", fontsize=9)
    ax_xy.set_aspect("equal")

    # Panel 2: XZ (perpendicular al disco)
    ax_xz = fig.add_subplot(gs[0, 1])
    for i in sample_idx:
        ax_xz.plot(positions[:, i, 0], positions[:, i, 2],
                   "-", alpha=0.4, lw=0.5, rasterized=True)
    ax_xz.set_xlabel("X galáctico (pc)")
    ax_xz.set_ylabel("Z galáctico (pc) — perpendicular al disco")
    ax_xz.set_title("Vista lateral del disco", fontsize=11, weight="bold")

    # Panel 3: distancia al centro galáctico vs tiempo
    ax_r = fig.add_subplot(gs[1, 0])
    r_t = np.linalg.norm(positions[:, sample_idx, :], axis=2)  # (n_snap, n_samples)
    for j in range(r_t.shape[1]):
        ax_r.plot(times / 1e6, r_t[:, j] / 1000, "-", alpha=0.4, lw=0.5)
    ax_r.set_xlabel("Tiempo (Myr)")
    ax_r.set_ylabel("Distancia al centro galáctico (kpc)")
    ax_r.set_title("Evolución radial", fontsize=11, weight="bold")

    # Panel 4: metadata de la simulación
    ax_meta = fig.add_subplot(gs[1, 1])
    ax_meta.axis("off")
    meta_lines = [
        "Simulación N-body (Script 07 v4)",
        "",
        f"  N estrellas:        {meta.get('n_stars', n_stars)}",
        f"  Tiempo simulado:    {meta.get('years_total', times[-1]):.2e} años",
        f"  dt:                 {meta.get('dt_days', '?')} días",
        f"  N pasos:            {meta.get('n_steps', '?'):,}",
        f"  N snapshots:        {n_snap}",
        "",
        f"  Integrador:         Leapfrog KDK (simpléctico)",
        f"  Potencial externo:  MN + Hernquist + NFW",
        "",
        f"  Energía inicial:    {meta.get('E_initial_J', 0):.3e} J",
        f"  Energía final:      {meta.get('E_final_J', 0):.3e} J",
        f"  |ΔE/E₀|:            {meta.get('energy_drift_rel', 0):.3e}",
        "",
        f"  Runtime:            {meta.get('runtime_seconds', 0):.1f} s",
        f"  Pasos/s:            {meta.get('steps_per_second', 0):.0f}",
        f"  Numba enabled:      {meta.get('numba_enabled', False)}",
        f"  Numba threads:      {meta.get('numba_threads', 1)}",
    ]
    ax_meta.text(0.05, 0.95, "\n".join(meta_lines), fontsize=9, family="monospace",
                 transform=ax_meta.transAxes, verticalalignment="top")

    fig.suptitle(f"Simulación dinámica: {n_stars} estrellas durante {times[-1]/1e6:.2f} Myr",
                 fontsize=13, weight="bold", y=0.995)

    pdf.savefig(fig, bbox_inches="tight", dpi=200)
    plt.close(fig)


def page_reproducibility(pdf: PdfPages, catalog_path: Path,
                           sim_path: Optional[Path], args: argparse.Namespace) -> None:
    """Página 10: apéndice con versiones + reproducibilidad."""
    fig, ax = plt.subplots(figsize=(8.27, 11.69))
    ax.axis("off")

    ax.text(0.5, 0.96, "Apéndice: Reproducibilidad",
            fontsize=14, weight="bold", ha="center", transform=ax.transAxes)

    y = 0.90
    ax.text(0.05, y, "Archivos de input:", fontsize=11, weight="bold",
            transform=ax.transAxes)
    y -= 0.025

    try:
        sha = compute_sha256(catalog_path)
        ax.text(0.05, y, f"  Catálogo: {catalog_path.name}", fontsize=9,
                family="monospace", transform=ax.transAxes)
        y -= 0.02
        ax.text(0.05, y, f"  SHA256:   {sha}", fontsize=8,
                family="monospace", transform=ax.transAxes)
        y -= 0.02
    except Exception as e:
        ax.text(0.05, y, f"  Error calculando hash: {e}", fontsize=9,
                transform=ax.transAxes)
        y -= 0.02

    if sim_path and sim_path.exists():
        try:
            sha_sim = compute_sha256(sim_path)
            ax.text(0.05, y, f"  Simulación: {sim_path.name}", fontsize=9,
                    family="monospace", transform=ax.transAxes)
            y -= 0.02
            ax.text(0.05, y, f"  SHA256:     {sha_sim}", fontsize=8,
                    family="monospace", transform=ax.transAxes)
            y -= 0.02
        except Exception:
            pass

    y -= 0.02
    ax.text(0.05, y, "Versiones de librerías:", fontsize=11, weight="bold",
            transform=ax.transAxes)
    y -= 0.025

    import numpy, pandas, matplotlib, scipy
    libs = [
        ("Python", platform.python_version()),
        ("numpy", numpy.__version__),
        ("pandas", pandas.__version__),
        ("matplotlib", matplotlib.__version__),
        ("scipy", scipy.__version__),
    ]
    try:
        import h5py as _h5
        libs.append(("h5py", _h5.__version__))
    except ImportError:
        pass
    try:
        import numba as _nb
        libs.append(("numba", _nb.__version__))
    except ImportError:
        pass
    try:
        import astropy
        libs.append(("astropy", astropy.__version__))
    except ImportError:
        pass

    for name, ver in libs:
        ax.text(0.1, y, f"{name}: {ver}", fontsize=9, family="monospace",
                transform=ax.transAxes)
        y -= 0.02

    y -= 0.03
    ax.text(0.05, y, "Comando para reproducir este reporte:", fontsize=11,
            weight="bold", transform=ax.transAxes)
    y -= 0.025

    cmd = (
        f"python 08_reporte_cientifico.py \\\n"
        f"  --catalog {catalog_path} \\\n"
    )
    if sim_path:
        cmd += f"  --simulation {sim_path} \\\n"
    cmd += f"  --output {args.output}"

    ax.text(0.05, y, cmd, fontsize=8, family="monospace",
            transform=ax.transAxes, verticalalignment="top")

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ═════════════════════════════════════════════════════════════════════════════
# ORQUESTADOR
# ═════════════════════════════════════════════════════════════════════════════

def generate_report(args: argparse.Namespace) -> None:
    """Genera el PDF completo."""
    catalog_path = Path(args.catalog)
    sim_path = Path(args.simulation) if args.simulation else None
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Cargar datos
    logging.info(f"📂 Cargando catálogo: {catalog_path}")
    df = load_catalog(catalog_path)
    stats = compute_catalog_stats(df)

    # Info básica
    logging.info(f"  → {stats.n_stars:,} estrellas | {stats.n_cols} columnas")
    logging.info(f"  → Cobertura V_radial: {stats.coverage_vrad:.1%}")
    logging.info(f"  → Cobertura UVW:       {stats.coverage_uvw:.1%}")

    # Generar PDF página por página
    with plot_context("scientific"):
        with PdfPages(output_path) as pdf:
            steps = [
                ("Portada", lambda: page_cover(pdf, stats, catalog_path, sim_path)),
                ("Estadísticas", lambda: page_catalog_stats(pdf, df, stats)),
                ("CMD", lambda: page_cmd(pdf, df)),
                ("Mapa galáctico", lambda: page_galactic_map(pdf, df)),
                ("Histogramas UVW", lambda: page_uvw_histograms(pdf, df)),
                ("Test KS", lambda: page_ks_test(pdf, df)),
                ("Clasificación", lambda: page_classification(pdf, stats)),
                ("Hypervelocity", lambda: page_hypervelocity(pdf, df, args.top_hypervelocity)),
            ]

            if sim_path and sim_path.exists():
                steps.append(
                    ("Simulación", lambda: page_simulation(pdf, sim_path, args.n_trajectory_samples))
                )

            steps.append(
                ("Reproducibilidad", lambda: page_reproducibility(pdf, catalog_path, sim_path, args))
            )

            for name, step_func in steps:
                try:
                    logging.info(f"📄 Generando página: {name}...")
                    step_func()
                except Exception as e:
                    logging.error(f"❌ Falló página '{name}': {e}", exc_info=True)

    # Tamaño del PDF generado
    size_mb = output_path.stat().st_size / 1e6
    logging.info(f"\n✅ Reporte generado: {output_path}")
    logging.info(f"   Tamaño: {size_mb:.2f} MB")
    logging.info(f"   Páginas: {len(steps)}")


# ═════════════════════════════════════════════════════════════════════════════
# CLI
# ═════════════════════════════════════════════════════════════════════════════

def main():
    p = argparse.ArgumentParser(
        description="Script 08: Reporte científico PDF del pipeline Gaia.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--catalog", required=True,
                   help="Catálogo procesado (Parquet del Script 03/04/06).")
    p.add_argument("--simulation", default=None,
                   help="HDF5 de simulación N-body (Script 07). Opcional.")
    p.add_argument("--output", default="reporte_gaia.pdf",
                   help="Ruta del PDF de salida.")
    p.add_argument("--n_trajectory_samples", type=int, default=30,
                   help="Trayectorias a plotear (si hay simulación).")
    p.add_argument("--top_hypervelocity", type=int, default=10,
                   help="Top N candidatas HVS a listar.")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()

    setup_logging(verbosity=2 if args.verbose else 1)

    print("=" * 70)
    print(" 📊 REPORTE CIENTÍFICO GAIA (Script 08)")
    print("=" * 70)

    try:
        generate_report(args)
    except (DataQualityError, ContractViolationError) as e:
        logging.error(f"❌ {type(e).__name__}: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        logging.warning("⚠ Interrumpido por usuario.")
        sys.exit(130)


if __name__ == "__main__":
    main()
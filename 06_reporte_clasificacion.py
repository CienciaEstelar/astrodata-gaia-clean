"""
SCRIPT 06: REPORTE Y CLASIFICACIÓN DE OBJETOS (v2 — Auditado)

Lee el catálogo depurado, clasifica los objetos basándose en criterios
astrofísicos 100% vectorizados y genera reportes profesionales.

CAMBIOS RESPECTO A v1:
- FIX CRÍTICO: las columnas que antes caían silenciosamente a NaN porque
  no coincidían con la salida del Script 03 ahora están alineadas. El
  Script 03 garantiza `luminosidad_absoluta_g`, `V_radial`, `fase_evolutiva`.
- FIX CRÍTICO: las reglas físicas estaban invertidas. M_G > 10 NO es Galaxia
  (es enana M / enana blanca); M_G ≈ -25 sí es Quasar; M_G ≈ -21 sí es Galaxia.
- AÑADIDO: validación explícita de contrato con Script 03. Si faltan
  columnas clave, el script falla RUIDOSAMENTE en vez de generar "Desconocido"
  masivamente sin avisar.
- AÑADIDO: reporte de columnas faltantes al inicio, antes de clasificar.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rich.console import Console
from rich.table import Table

console = Console()

# Columnas que el pipeline depurado DEBE proveer. Contrato con Script 03.
COLS_REQUERIDAS = [
    "clasificacion_espectral",
    "luminosidad_absoluta_g",
    "color_bp_rp",
    "V_radial",
    "fase_evolutiva",
]


def cargar_catalogo(ruta: str) -> pd.DataFrame:
    """Carga el catálogo detectando si es Parquet o CSV."""
    path = Path(ruta)
    if not path.is_file():
        console.print(f"[bold red]❌ Error:[/bold red] No se encuentra: {ruta}")
        sys.exit(1)

    try:
        if path.suffix.lower() == '.parquet':
            return pd.read_parquet(path)
        return pd.read_csv(path, low_memory=False)
    except Exception as e:
        console.print(f"[bold red]❌ Fallo al leer catálogo:[/bold red] {e}")
        sys.exit(1)


def validar_contrato(df: pd.DataFrame, estricto: bool) -> None:
    """Verifica que existan las columnas requeridas.

    Si `estricto=True`, aborta con código de error.
    Si `estricto=False`, emite warning y rellena con NaN (comportamiento legacy).
    """
    faltantes = [c for c in COLS_REQUERIDAS if c not in df.columns]
    if faltantes:
        console.print(
            f"[bold yellow]⚠️ Columnas faltantes respecto al contrato S03:[/bold yellow] {faltantes}"
        )
        if estricto:
            console.print(
                "[bold red]Abortando. Re-ejecuta el Script 03 (v2) para regenerar el catálogo.[/bold red]"
            )
            sys.exit(2)
        else:
            console.print("[yellow]Rellenando con NaN (modo permisivo).[/yellow]")
            for col in faltantes:
                df[col] = np.nan
    else:
        console.print("[green]✔️ Contrato con Script 03 verificado.[/green]")


def clasificar_objetos_vectorizado(df: pd.DataFrame) -> pd.Series:
    """Motor de clasificación vectorizado con reglas astrofísicamente correctas.

    ORDEN DE PRECEDENCIA (importante porque np.select usa la primera condición True):
      1. Keywords de objetos extragalácticos en 'clasificacion_espectral' (si vinieran
         de cross-match con catálogos externos tipo SDSS spectroscopic).
      2. Reglas por magnitud absoluta M_G (rangos físicamente sensatos).
      3. Reglas por letra espectral (MKGFABO + estrella genérica O).
      4. Reglas por fase evolutiva.
      5. Reglas por cinemática extrema.

    REFERENCIAS DE M_G (magnitud absoluta en banda G de Gaia):
      Quasar típico: M_G ≈ -25 a -28
      Galaxia L*:    M_G ≈ -21 a -23
      Supergigante:  M_G ≈ -6 a -10
      Sol:           M_G ≈ +4.67
      Enana M:       M_G ≈ +8 a +14
      Enana Blanca:  M_G ≈ +10 a +15
    """
    espectro = df['clasificacion_espectral'].astype(str).fillna('')
    fase = df['fase_evolutiva'].astype(str).fillna('').str.lower()
    v_rad_abs = df['V_radial'].abs()
    M_G = df['luminosidad_absoluta_g']

    # -- DEFINICIÓN DE REGLAS (orden de precedencia decreciente) --
    condiciones = [
        # 1. Keywords extragalácticas explícitas (si el catálogo las trae)
        espectro.str.contains('quasar|qso', case=False, na=False, regex=True),
        espectro.str.contains('blazar', case=False, na=False),
        espectro.str.contains('galaxy|galaxia|agn', case=False, na=False, regex=True),
        espectro.str.contains('nebula|nebulosa', case=False, na=False, regex=True),

        # 2. Por magnitud absoluta (FÍSICAMENTE CORRECTO AHORA)
        M_G < -22,                                  # Quasares/AGN luminosos
        (M_G >= -22) & (M_G < -15),                 # Galaxias o supergigantes muy luminosas
        (M_G >= -15) & (M_G < -6),                  # Supergigantes
        (M_G > 12),                                 # Enanas blancas / enanas M muy tardías
        (M_G > 8) & (M_G <= 12),                    # Enanas M / K tardías

        # 3. Por letra espectral (si sobreviven a las reglas anteriores)
        espectro.str.startswith('O', na=False),
        espectro.str.startswith('B', na=False),
        espectro.str.startswith('A', na=False),
        espectro.str.startswith('F', na=False),
        espectro.str.startswith('G', na=False),
        espectro.str.startswith('K', na=False),
        espectro.str.startswith('M', na=False),

        # 4. Por fase evolutiva (si el tipo espectral falla)
        fase.str.contains('gigante', na=False),
        fase.str.contains('enana blanca', na=False),
        fase.str.contains('secuencia principal', na=False),

        # 5. Por cinemática extrema (hypervelocity stars)
        v_rad_abs > 500,
    ]

    resultados = [
        # 1. Keywords
        'Quasar/AGN', 'Blazar', 'Galaxia', 'Nebulosa',
        # 2. Por M_G
        'Quasar/AGN candidato', 'Galaxia/Supergigante luminosa', 'Supergigante',
        'Enana Blanca candidata', 'Enana M/K tardía',
        # 3. Por letra espectral
        'Estrella O', 'Estrella B', 'Estrella A', 'Estrella F',
        'Estrella G', 'Estrella K', 'Estrella M',
        # 4. Por fase
        'Gigante', 'Enana Blanca', 'Secuencia Principal',
        # 5. Cinemática
        'Estrella de alta velocidad',
    ]

    assert len(condiciones) == len(resultados), (
        f"Desalineación condiciones/resultados: {len(condiciones)} vs {len(resultados)}"
    )

    return pd.Series(
        np.select(condiciones, resultados, default='Desconocido'),
        index=df.index
    )


def mostrar_reporte(df: pd.DataFrame) -> None:
    """Tabla de resumen con Rich."""
    conteo = df['tipo_objeto'].value_counts().reset_index()
    conteo.columns = ['Tipo', 'Cantidad']

    tabla = Table(
        title="Clasificación Gaia DR3 — Catálogo Depurado",
        header_style="bold cyan"
    )
    tabla.add_column("Tipo de Objeto", style="yellow")
    tabla.add_column("Cantidad", justify="right", style="green")
    tabla.add_column("Porcentaje", justify="right", style="magenta")

    total = len(df)
    for _, row in conteo.iterrows():
        porcentaje = (row['Cantidad'] / total) * 100
        tabla.add_row(
            str(row['Tipo']),
            f"{int(row['Cantidad']):,}",
            f"{porcentaje:.2f}%"
        )

    console.print(tabla)

    # Alerta sanitaria: si más del 50% cae en 'Desconocido' probablemente hay
    # un problema aguas arriba (columnas faltantes, cortes demasiado estrictos).
    desconocidos = conteo[conteo['Tipo'] == 'Desconocido']['Cantidad'].sum()
    if total > 0 and desconocidos / total > 0.5:
        console.print(
            f"[bold yellow]⚠️ {desconocidos:,} objetos ({desconocidos/total*100:.1f}%) "
            "quedaron sin clasificar. Revisa cobertura de teff_gspphot y paralaje.[/bold yellow]"
        )


def main():
    p = argparse.ArgumentParser(description="Clasificador de objetos Gaia DR3 (v2).")
    p.add_argument("--input", required=True,
                   help="Catálogo Parquet o CSV producido por el Script 03.")
    p.add_argument("--output", required=True,
                   help="Ruta de salida para el catálogo clasificado (Parquet o CSV).")
    p.add_argument(
        "--permisivo", action="store_true",
        help="No abortar si faltan columnas del contrato; rellenar con NaN."
    )
    args = p.parse_args()

    console.print("[bold cyan]🔭 Clasificador Gaia DR3 — Script 06[/bold cyan]")

    # 1. Cargar
    df = cargar_catalogo(args.input)
    console.print(f"► Cargadas {len(df):,} filas.")

    # 2. Validar contrato (aborta si faltan columnas clave, salvo --permisivo)
    validar_contrato(df, estricto=not args.permisivo)

    # 3. Clasificar
    df['tipo_objeto'] = clasificar_objetos_vectorizado(df)

    # 4. Reportar
    mostrar_reporte(df)

    # 5. Guardar
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.suffix == '.parquet':
        df.to_parquet(out_path, index=False)
    else:
        df.to_csv(out_path, index=False)

    console.print(f"\n[bold green]✔️ Finalizado.[/bold green] Catálogo en: {out_path}")


if __name__ == "__main__":
    main()
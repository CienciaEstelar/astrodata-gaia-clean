# Pipeline de Análisis Cinemático Gaia DR3

Pipeline científico completo para descarga, enriquecimiento, clasificación y simulación dinámica de estrellas usando datos del catálogo **Gaia DR3** de la ESA.

Desarrollado como proyecto de investigación en astrofísica computacional. Diseñado para ser ejecutado secuencialmente, desde la descarga de imágenes SDSS hasta la generación de un reporte PDF con estadísticas publicables.

---

## Requisitos previos

- **Python 3.10 o superior**
- Conexión a internet (Scripts 01 y 02 acceden a servidores externos de SDSS y Gaia)
- ~2 GB de espacio libre en disco para un catálogo de tamaño mediano
- Conocimientos recomendados: Python intermedio + nociones básicas de astrofísica (magnitudes, paralajes, diagramas HR)

---

## Instalación

```bash
# 1. Clonar o descargar el repositorio
cd GAIA-refactorizado

# 2. Crear entorno virtual (recomendado)
python -m venv venv
source venv/bin/activate        # Linux/Mac
# venv\Scripts\activate         # Windows

# 3. Instalar dependencias
pip install -r requirements.txt

# 4. (Opcional pero recomendado) Instalar aceleradores
pip install numba joblib sympy
```

> **Nota sobre Numba:** El Script 07 (simulador N-cuerpos) puede correr sin Numba usando numpy puro, pero será 50-100x más lento. Para simulaciones de más de 500 estrellas se recomienda instalarlo.

---

## Orden de ejecución

Los scripts están numerados. Deben correrse en orden: la salida de cada uno es la entrada del siguiente.

```
00 → 01 → 02 → 03 → 04 → 05 → 06 → 07 → 08
```

### Script 00 — Selftest del pipeline
```bash
python 00_selftest_pipeline.py
```
Corre una batería de tests físicos antes de procesar datos reales. Verifica constantes, masas estelares, potencial galáctico e integridad del integrador. Si algún test falla, algo está mal con la instalación.

---

### Script 01 — Descarga de imágenes SDSS
```bash
python 01_descarga_sdss.py \
  --base-url http://data.sdss.org/sas/dr17/eboss/photoObj/frames/301/1000/6 \
  --start 27 --end 30 \
  --dest data/FITS_BANDAS \
  --bands g r i
```
Descarga archivos FITS del Sloan Digital Sky Survey (SDSS) en paralelo con reintentos automáticos. Los FITS son las imágenes del cielo que definen la zona a consultar en Gaia.

---

### Script 02 — Consulta ADQL a Gaia DR3
```bash
python 02_consulta_gaia_adql.py \
  --fits-dir data/FITS_BANDAS \
  --output data/catalogo_gaia_raw.parquet
```
Extrae el footprint WCS de cada imagen FITS y consulta el servidor TAP de Gaia para obtener las estrellas dentro de ese campo. Genera el catálogo crudo en formato Parquet.

---

### Script 03 — Enriquecimiento del catálogo
```bash
python 03_pipeline_enriquecimiento.py \
  --input data/catalogo_gaia_raw.parquet \
  --output data/catalogo_enriquecido.parquet
```
Calcula distancias (con corrección de zero-point de paralaje, Lindegren+2021), magnitudes absolutas, velocidades tangenciales, coordenadas galácticas y clasificación espectral. Aplica filtros de calidad astrométrica y fotométrica estándar de Gaia DR3.

---

### Script 04 — Cinemática UVW
```bash
python 04_cinematica_uvw.py \
  --input data/catalogo_enriquecido.parquet \
  --output data/catalogo_con_uvw.parquet \
  --calculate_lsr
```
Transforma las velocidades observadas (movimiento propio + velocidad radial) a componentes galactocéntricas U, V, W usando Astropy. Incluye transformación al Local Standard of Rest (LSR).

---

### Script 05 — Diagrama Color-Magnitud (CMD/HR)
```bash
python 05_generador_cmd.py \
  --input data/catalogo_enriquecido.parquet \
  --output figuras/cmd.png \
  --mode hexbin \
  --overlay_regions
```
Genera el diagrama de Hertzsprung-Russell (CMD) con densidad hexagonal y overlay de regiones canónicas (Secuencia Principal, Rama de Gigantes, Enanas Blancas).

---

### Script 06 — Clasificación de objetos
```bash
python 06_reporte_clasificacion.py \
  --input data/catalogo_enriquecido.parquet \
  --output data/catalogo_clasificado.parquet
```
Clasifica cada objeto (estrella G, enana M, gigante, quasar, etc.) usando reglas vectorizadas basadas en magnitud absoluta, tipo espectral y cinemática.

---

### Script 07 — Simulador N-cuerpos
```bash
python 07_simulador_n_cuerpos.py \
  --input data/catalogo_con_uvw.parquet \
  --output data/simulacion.h5 \
  --n_stars 500 \
  --years 1000000
```
Simula la evolución dinámica de las estrellas en el potencial galáctico (disco Miyamoto-Nagai + bulbo Hernquist + halo NFW) usando un integrador Leapfrog simpléctico. Requiere que el catálogo de entrada tenga columnas UVW del Script 04.

---

### Script 08 — Reporte científico PDF
```bash
python 08_reporte_cientifico.py \
  --catalog data/catalogo_clasificado.parquet \
  --simulation data/simulacion.h5 \
  --output reporte_final.pdf
```
Genera un PDF de 10 páginas con: portada, estadísticas, CMD, mapa galáctico Aitoff, distribuciones UVW, test KS, tabla de clasificación, candidatas a hypervelocity stars, trayectorias de simulación y apéndice de reproducibilidad.

---

## Estructura del proyecto

```
.
├── 00_selftest_pipeline.py    # Tests de integridad física
├── 01_descarga_sdss.py        # Descarga FITS desde SDSS
├── 02_consulta_gaia_adql.py   # Consulta TAP a Gaia DR3
├── 03_pipeline_enriquecimiento.py  # Física + filtros de calidad
├── 04_cinematica_uvw.py       # Velocidades galactocéntricas UVW
├── 05_generador_cmd.py        # Diagrama Color-Magnitud
├── 06_reporte_clasificacion.py     # Clasificación de objetos
├── 07_simulador_n_cuerpos.py    # Simulación dinámica N-cuerpos
├── 08_reporte_cientifico.py   # Reporte PDF publicable
├── gaia_core.py               # Módulo fundacional (excepciones, validadores)
├── requirements.txt           # Dependencias
└── data/                      # Directorio de datos (generado al correr el pipeline)
```

---

## Nivel de dificultad

| Script | Nivel Python | Conceptos de astrofísica |
|---|---|---|
| 01, 05 | Intermedio | Básico |
| 02, 03, 06 | Intermedio-avanzado | Fotometría, paralajes, clasificación espectral |
| 04 | Avanzado | Cinemática galáctica, sistema LSR |
| 07 | Avanzado | Mecánica celeste, potenciales gravitacionales |
| 00, 08 | Avanzado | Transversal |

---

## Fuentes de datos externas

- **SDSS** (Script 01): [data.sdss.org](https://data.sdss.org) — imágenes fotométricas
- **Gaia DR3** (Script 02): servidor TAP público de la ESA [gea.esac.esa.int](https://gea.esac.esa.int/tap-server/tap)
- No se requieren credenciales ni API keys

---

## Referencias científicas principales

- Gaia Collaboration (2022) — *Gaia DR3 summary*, A&A 674, A1
- Lindegren et al. (2021) — *Parallax zero-point*, A&A 649, A4
- Schönrich, Binney & Dehnen (2010) — *Local kinematics and LSR*, MNRAS 403, 1829
- Bovy (2015) — *MWPotential2014*, ApJS 216, 29

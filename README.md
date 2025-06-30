# GAIA Enrichment Pipeline

![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green)
![Status](https://img.shields.io/badge/status-stable-brightgreen)

Pipeline CLI en Python para filtrar, limpiar y enriquecer catálogos Gaia DR2/EDR3/DR3.  
Aplica filtros fotométricos y astrométricos recomendados por la literatura y agrega parámetros físicos derivados como distancia, magnitud absoluta, velocidades 3D y clasificaciones estelares (tipo espectral, fase evolutiva, variabilidad).

---

## ✨ Características

- ✅ Filtros de calidad vectorizados (RUWE, BP/RP excess, χ²/dof, etc.)
- ✅ Parámetros derivados: `dist_pc`, `M_G`, `V_3D`, `l_gal`, `b_gal`
- ✅ Clasificaciones rápidas por temperatura y color
- ✅ Logging dual (archivo + consola)
- ✅ Barra de progreso Rich para seguimiento por etapas
- ✅ CLI completa vía `argparse`

---

## 🚀 Instalación

Requiere Python ≥ 3.9 y las siguientes bibliotecas:

```bash
pip install pandas numpy astropy rich


---

🧪 Uso rápido

python gaia_enrichment_pipeline.py \
  --input  CatalogoPropio_Gaia.csv \
  --output Catalogo_Gaia_Enriquecido_v2.csv \
  --use_excess_factor_filter \
  --use_astro_chi2_filter

> Opcional: agrega -v para ver logs detallados (nivel DEBUG).




---

⚙️ Argumentos principales

Flag	Descripción

--input	CSV original exportado desde Gaia Archive
--output	Nombre del CSV enriquecido que se generará
--use_excess_factor_filter	Activa filtro de fotometría contaminada
--use_astro_chi2_filter	Activa filtro astrométrico adicional (χ²/dof)
--max_dist_pc	(Opcional) Corta estrellas más lejanas
--min_abs_b_gal	(Opcional) Corta zonas cercanas al plano galáctico
--verbose o -v	Muestra logs extendidos mientras corre



---

📂 Ejemplo de salida

El archivo enriquecido incluirá columnas como:

dist_pc, M_G, bp_rp, V_3D

l_gal, b_gal (coordenadas galácticas)

tipo_espectral (OBAFGKM)

fase_evolutiva (Enana blanca, Secuencia principal, Gigante, etc.)

variabilidad (según RUWE)



---

📊 Visualización rápida (HR Diagram)

import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("Catalogo_Gaia_Enriquecido_v2.csv")
plt.scatter(df["bp_rp"], df["M_G"], s=1, c=df["bp_rp"], cmap="plasma_r")
plt.gca().invert_yaxis()
plt.xlabel("BP – RP")
plt.ylabel("M_G")
plt.title("Diagrama HR – Gaia")
plt.show()


---

📘 Referencias

Gaia Collaboration et al. (2018, 2021)

Lindegren et al. (2018) — RUWE y χ² filtering

Evans et al. (2018) — BP/RP Excess Factor

Bailer-Jones (2018) — Distance estimation from parallax



---

🪪 Licencia

Este proyecto está bajo licencia MIT. Eres libre de usarlo, adaptarlo y compartirlo con atribución.


---

🤝 Agradecimientos

Desarrollado con fines educativos y científicos.
Compatible con análisis de Astroinformática, evolución estelar y dinámica galáctica.


---

¿Te gustaría que también genere una versión en inglés para publicar en GitHub o dejarlo en español profesional?


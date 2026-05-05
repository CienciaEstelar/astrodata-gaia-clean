# Tutorial del Pipeline Gaia DR3
### De la descarga de datos estelares a la simulación dinámica — paso a paso

---

## BLOQUE 0 — ¿Qué vamos a hacer aquí?

Imaginate que tienes acceso a las coordenadas precisas, los colores y las velocidades de **más de mil millones de estrellas**. No estrellas inventadas en una simulación: estrellas reales de la Vía Láctea, observadas durante años por un satélite de la Agencia Espacial Europea llamado **Gaia**.

Eso es exactamente lo que existe. Y es exactamente lo que vamos a usar.

### Qué es Gaia DR3

Gaia es un telescopio espacial que lleva desde 2013 midiendo con una precisión increíble la posición, brillo y movimiento de cada estrella que puede ver. En 2022 publicó su tercer catálogo de datos, **Gaia DR3** ("Data Release 3"), con información de casi 1.500 millones de objetos.

El dato más valioso que mide Gaia es la **paralaje**: un pequeño desplazamiento angular que la estrella parece tener cuando la Tierra se mueve alrededor del Sol. Si conocés la paralaje, podés calcular la distancia. Y si conocés la distancia más el movimiento aparente en el cielo, podés calcular la velocidad real de la estrella en el espacio.

Esto es una revolución porque hasta Gaia, para la mayoría de las estrellas solo sabíamos su dirección. Ahora sabemos dónde están en 3D **y** a qué velocidad se mueven. Es como pasar de tener un mapa plano a tener un mapa con GPS en tiempo real.

### El problema científico que resuelve este pipeline

Una pregunta fundamental de la astronomía galáctica es: **¿cómo se mueven las estrellas cerca del Sol y qué nos dice eso sobre la estructura de la Vía Láctea?**

Para responderla necesitás hacer varias cosas encadenadas:
1. Descargar los datos de Gaia para una zona del cielo específica
2. Calcular distancias, magnitudes absolutas y velocidades reales
3. Clasificar qué tipo de objeto es cada estrella
4. Simular cómo esas estrellas van a evolucionar dinámicamente bajo la gravedad galáctica

Eso es exactamente lo que hace este pipeline, script por script.

### Qué vas a poder hacer al terminar

Cuando completes este tutorial vas a ser capaz de:

- Consultar el catálogo de Gaia directamente desde Python usando lenguaje ADQL
- Calcular distancias a estrellas a partir de paralajes con la corrección científica correcta
- Clasificar estrellas en el Diagrama Hertzsprung-Russell
- Calcular velocidades galácticas U, V, W en el sistema de referencia correcto
- Correr una simulación N-cuerpos de cientos de estrellas en el potencial gravitacional de la Vía Láctea
- Generar un reporte PDF con estadísticas y gráficos de nivel publicable

No hace falta que seas experto en astrofísica para llegar hasta ahí. Solo necesitás ganas de aprender y Python.

---

## BLOQUE 1 — Antes de empezar 🛠️

### Qué necesitás saber

| Necesario | No necesario |
|---|---|
| Python básico (funciones, listas, bucles) | Astrofísica avanzada |
| Saber instalar paquetes con `pip` | Álgebra lineal avanzada |
| Entender qué es un DataFrame de pandas | Experiencia previa con datos astronómicos |
| Saber leer un mensaje de error en Python | Conocer los formatos FITS o Parquet de antemano |

Si alguna vez usaste pandas para abrir un CSV y matplotlib para graficar algo, estás listo.

### Instalación

```bash
# 1. Clonar el repositorio
git clone https://github.com/juangalaz/pipeline-gaia-dr3.git
cd pipeline-gaia-dr3

# 2. Crear un entorno virtual (muy recomendado para no ensuciar tu Python global)
python -m venv venv
source venv/bin/activate        # En Linux/Mac
# venv\Scripts\activate         # En Windows

# 3. Instalar todo lo necesario
pip install -r requirements.txt
```

Si querés el máximo rendimiento en el simulador, también instalá esto:

```bash
pip install numba joblib sympy
```

Con `numba` el Script 07 va a correr unas 50-100 veces más rápido. Sin él funciona igual, solo más lento.

### Verificar que todo funciona: el selftest

Antes de tocar datos reales, el Script 00 corre una batería de tests físicos para asegurarse de que todo está instalado y funcionando correctamente. Es como hacer un chequeo médico antes de correr una maratón.

```bash
python 00_selftest_pipeline.py
```

Si todo está bien, vas a ver algo así:

```
══════════════════════════════════════════════════════════════════════
 🛡️  GAIA PIPELINE SELFTEST — Guardián del Proyecto
══════════════════════════════════════════════════════════════════════

═══ A — Constantes y Unidades [4/4 passed] ═══
  ✔ PASS  PC_TO_M (val=3.0856...e+16, expected=3.0856...e+16)
  ✔ PASS  YEAR_TO_S
  ✔ PASS  M_SUN
  ✔ PASS  M_BOL_SUN

═══ B — Masas Estelares (M-L Segmentada) [5/5 passed] ═══
  ✔ PASS  mass_Sun (val=0.98, expected=[0.9, 1.1])
  ✔ PASS  mass_M_dwarf
  ...

 ✅  RESULTADO GLOBAL: 28/28 tests pasaron
══════════════════════════════════════════════════════════════════════
```

Si algún test falla, el mensaje de error te va a decir exactamente qué está roto. Lo más común es que falte alguna librería: instalala con `pip install <nombre>` y volvé a correr.

> **Tip:** Si ves `❌ Suite 'X' lanzó excepción: ModuleNotFoundError`, significa que falta un paquete. Fijate cuál es en el mensaje y añadilo al comando de instalación.

---

## BLOQUE 2 — Recorrido script por script 🔭

### Script 00 — El guardián: selftest del pipeline

**Pregunta científica:** ¿Está todo configurado correctamente para producir resultados físicamente válidos?

**Input:** Nada (crea datos sintéticos internamente).
**Output:** Reporte de tests en consola. No genera archivos.

**Concepto clave — Tests físicos:**
Un test físico es diferente a un test de software normal. No chequea si una función devuelve el valor correcto de un mock: chequea que el resultado tiene sentido en el universo real. Por ejemplo, el test C verifica que la velocidad circular a 8 kpc del centro galáctico está entre 200 y 260 km/s, que es exactamente lo que miden los astrónomos con otras técnicas.

**Output esperado:**
```
═══ C — Potencial Galáctico MW [4/4 passed] ═══
  ✔ PASS  v_circ_at_R_sol (val=228.4, expected=[200, 260])
  ✔ PASS  disk_attractive_at_Rsol
  ✔ PASS  bulge_attractive_at_Rsol
  ✔ PASS  halo_attractive_at_Rsol
```

**Error común:**
```
❌ Suite 'constants' lanzó excepción: ModuleNotFoundError: No module named 'numba'
```
Solución: `pip install numba` — o simplemente ignorarlo si no querés Numba; el simulador tiene fallback automático.

---

### Script 01 — Descarga de imágenes SDSS

**Pregunta científica:** ¿De qué zona exacta del cielo queremos estudiar estrellas?

**Input:** URL del servidor SDSS, números de imagen (run, camcol, frame).
**Output:** Archivos `.fits` (imágenes astronómicas) en la carpeta de destino.

**Concepto clave — Formato FITS:**
FITS (Flexible Image Transport System) es el formato estándar de la astronomía, como el JPG de las fotos pero para imágenes científicas. Además de los datos de la imagen, cada archivo FITS tiene un "header" con metadatos: la posición en el cielo (RA/Dec), la escala de la imagen, la fecha de observación, etc. El Script 02 va a leer esos metadatos para saber exactamente qué zona del cielo consultar en Gaia.

```bash
python 01_descarga_sdss.py \
  --base-url http://data.sdss.org/sas/dr17/eboss/photoObj/frames/301/1000/6 \
  --start 27 --end 30 \
  --dest data/FITS_BANDAS \
  --bands g r \
  --threads 4
```

**Output esperado:**
```
Total tareas: 8 | Run=001000 | Camcol=6
Descargando... ████████████████████ 8/8 [0:01:23 / 0:00:00]
✅ Descarga completa. Éxitos: 8/8
```

**Error común:**
```
requests.exceptions.ConnectionError: HTTPSConnectionPool(host='data.sdss.org', ...)
```
Solución: verificá tu conexión a internet y que la URL sea accesible. SDSS a veces tiene mantención. Podés probar abrir la URL en el navegador primero.

---

### Script 02 — Consulta ADQL a Gaia DR3

**Pregunta científica:** ¿Qué estrellas hay exactamente en la zona del cielo que fotografiamos con SDSS?

**Input:** Las imágenes FITS descargadas en el paso anterior.
**Output:** Catálogo con miles de estrellas en formato Parquet.

**Concepto clave — ADQL y TAP:**
Para consultar datos astronómicos remotos existe un protocolo llamado **TAP** (Table Access Protocol) y un lenguaje de consulta llamado **ADQL** (Astronomical Data Query Language). ADQL es básicamente SQL con funciones especiales para trabajar con coordenadas celestes. Por ejemplo:

```sql
SELECT source_id, ra, dec, parallax, phot_g_mean_mag
FROM gaiadr3.gaia_source
WHERE 1=CONTAINS(
  POINT('ICRS', ra, dec),
  BOX('ICRS', 215.4, 52.7, 0.5, 0.5)
)
```

Esto le dice al servidor de Gaia: "dame todas las estrellas dentro de un cuadrado de 0.5° × 0.5° centrado en RA=215.4°, Dec=52.7°". El script construye esta consulta automáticamente para cada imagen FITS.

```bash
python 02_consulta_gaia_adql.py \
  --fits-dir data/FITS_BANDAS \
  --output data/catalogo_gaia_raw.parquet \
  --band r
```

**Output esperado:**
```
🌌 Gaia DR3 ADQL Extractor
  FITS encontrados: 4 (banda r)
  Log detallado: 02_consulta_gaia.log

Consultando Gaia DR3 ████████████████████ 4/4 • 100% • 0:02:15 / 0:00:00 • 12,847 ⭐

✅ Catálogo guardado
  Ruta:        data/catalogo_gaia_raw.parquet
  Estrellas:   12,847
  Tamaño:      8.3 MB
```

**Error común:**
```
RuntimeError: Fallo crítico conectando a Gaia TAP: ...
```
Solución: el servidor TAP de Gaia a veces está congestionado. Esperá unos minutos y volvé a intentar. El script tiene reintentos automáticos, pero si el servidor está caído del todo hay que esperar.

---

### Script 03 — Enriquecimiento del catálogo

**Pregunta científica:** ¿A qué distancia están estas estrellas y qué tan brillantes son en realidad?

**Input:** Catálogo crudo de Gaia (Parquet).
**Output:** Catálogo enriquecido con distancias, magnitudes absolutas, velocidades y clasificaciones.

**Concepto clave — Paralaje zero-point:**
La paralaje que mide Gaia tiene un pequeño sesgo sistemático de aproximadamente −0.017 mas (mili-arco-segundos) descubierto por Lindegren et al. en 2021. Si no lo corregís, las distancias calculadas tienen un error del ~1% a 100 pc y del ~10% a 1 kpc. Parece poco, pero cuando estás estudiando cinemática galáctica ese error acumulado arruina los resultados. El Script 03 aplica esta corrección automáticamente.

La fórmula de la distancia es simple:
```
distancia_pc = 1000 / (paralaje_mas - zero_point_mas)
```

```bash
python 03_pipeline_enriquecimiento.py \
  --input data/catalogo_gaia_raw.parquet \
  --output data/catalogo_enriquecido.parquet \
  --parallax_zp_mas -0.017
```

**Output esperado en consola:**
```
🌟 Pipeline Gaia DR3 — Enriquecimiento v3
  Input:         data/catalogo_gaia_raw.parquet
  Parallax ZP:   -0.017 mas (Lindegren+2021)

Procesando  ██████████████████████  7/7 • 0:00:14

✅ Proceso completado
  Estrellas finales:   9,341
  Cobertura V_radial:  34.2%
  Tamaño archivo:      12.7 MB
  Contrato S06:        ✔ OK
```

**Error común:** Ver que "Estrellas finales" es mucho menor que las originales (por ejemplo 2,000 de 12,000).
Explicación: los filtros de calidad astrométrica (RUWE < 1.4, S/N paralaje > 10, etc.) descartan estrellas con mediciones poco confiables. Es completamente normal: Gaia mide bien solo las estrellas con suficiente brillo y observaciones limpias. Si querés ser menos estricto, probá `--max_ruwe 2.0` o `--parallax_snr_min 5`.

---

### Script 04 — Cinemática UVW

**Pregunta científica:** ¿A qué velocidad se mueven estas estrellas en el espacio, en coordenadas galácticas?

**Input:** Catálogo enriquecido (con distancias y movimientos propios).
**Output:** Catálogo con columnas U_kms, V_kms, W_kms agregadas.

**Concepto clave — Velocidades UVW:**
Cuando vemos una estrella en el cielo, podemos medir dos cosas sobre su movimiento:
- **Movimiento propio** (pmra, pmdec): qué tan rápido se mueve en el plano del cielo, en arco-segundos por año.
- **Velocidad radial**: qué tan rápido se acerca o aleja de nosotros, en km/s (medida por efecto Doppler).

Combinando las tres componentes y conociendo la distancia, podemos calcular la velocidad 3D de la estrella en el sistema galactocéntrico:
- **U**: velocidad hacia/desde el centro galáctico
- **V**: velocidad en la dirección de la rotación galáctica
- **W**: velocidad perpendicular al disco galáctico

Esto lo hace Astropy con una transformación de coordenadas precisa.

```bash
python 04_cinematica_uvw.py \
  --input data/catalogo_enriquecido.parquet \
  --output data/catalogo_con_uvw.parquet \
  --calculate_lsr
```

**Output esperado:**
```
✅ PROCESO COMPLETADO
► Post-transformación: 3,198/9,341 filas con UVW 3D completo.
```
El número menor se debe a que solo las estrellas con velocidad radial medida pueden tener UVW completo. Las que solo tienen movimiento propio quedan con NaN en UVW, lo cual es honesto científicamente.

**Error común:**
```
ValueError: Faltan columnas críticas: ['distancia_parsecs']
```
Solución: asegurate de haber corrido el Script 03 primero y estar usando su archivo de salida como input.

---

### Script 05 — Diagrama Color-Magnitud (CMD)

**Pregunta científica:** ¿Qué tipos de estrellas hay en nuestra muestra?

**Input:** Catálogo enriquecido (necesita color BP-RP y magnitud absoluta G).
**Output:** Imagen PNG o PDF con el diagrama.

**Concepto clave — El Diagrama Hertzsprung-Russell:**
El CMD (Color-Magnitude Diagram) es el mapa más importante de la astronomía estelar. En el eje X va el color de la estrella (azul = caliente, rojo = frío) y en el eje Y va su brillo real (magnitud absoluta). Cuando graficás miles de estrellas, aparecen patrones clarísimos:

- **Secuencia Principal**: una banda diagonal donde viven la mayoría de las estrellas "normales" (como el Sol). Las estrellas pasan la mayor parte de su vida aquí.
- **Rama de Gigantes**: estrellas en la parte roja y brillante del gráfico. Son estrellas viejas que se expandieron.
- **Enanas Blancas**: abajo a la izquierda. Son los núcleos fríos de estrellas que ya murieron.

```bash
python 05_generador_cmd.py \
  --input data/catalogo_enriquecido.parquet \
  --output figuras/cmd_hexbin.png \
  --mode hexbin \
  --overlay_regions
```

**Output esperado:** Una imagen donde claramente se distingue la Secuencia Principal como una banda diagonal densa, con la Rama de Gigantes curvándose hacia arriba a la derecha.

**Error común:** El diagrama sale con todos los puntos amontonados en una esquina, sin estructura visible.
Solución: probablemente hay valores extremos (outliers) en los datos. Usá `--xlim -0.5,4 --ylim 15,-5` para fijar los límites del gráfico a la región físicamente relevante.

---

### Script 06 — Clasificación de objetos

**Pregunta científica:** ¿Qué tipo de objeto astronómico es cada entrada del catálogo?

**Input:** Catálogo enriquecido (con tipo espectral y magnitud absoluta).
**Output:** Catálogo con columna `tipo_objeto` y tabla de resumen en consola.

**Concepto clave — Clasificación vectorizada:**
En vez de recorrer el DataFrame fila por fila (lo cual sería lentísimo para 10.000 estrellas), el script usa `numpy.select` para aplicar todas las reglas de clasificación de una vez sobre arrays enteros. Es el equivalente de hacer 20 preguntas simultáneamente a todas las filas:

- Si M_G < −22 → Quasar/AGN candidato
- Si la letra espectral empieza con "G" → Estrella G
- Si M_G > 12 → Enana Blanca candidata
- ... y así sucesivamente, en orden de precedencia

```bash
python 06_reporte_clasificacion.py \
  --input data/catalogo_enriquecido.parquet \
  --output data/catalogo_clasificado.parquet
```

**Output esperado:**
```
✔️ Contrato con Script 03 verificado.

 Clasificación Gaia DR3 — Catálogo Depurado
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━┳━━━━━━━━━━━━┓
┃ Tipo de Objeto              ┃ Cantidad ┃ Porcentaje ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━╇━━━━━━━━━━━━┩
│ Estrella K                  │    3,241 │     34.70% │
│ Estrella G                  │    2,108 │     22.57% │
│ Estrella M                  │    1,893 │     20.27% │
│ Secuencia Principal         │      887 │      9.50% │
│ Enana M/K tardía            │      612 │      6.55% │
│ Gigante                     │      389 │      4.17% │
│ Enana Blanca candidata      │      211 │      2.26% │
└─────────────────────────────┴──────────┴────────────┘
```

**Error común:**
```
⚠️ Columnas faltantes respecto al contrato S03: ['fase_evolutiva']
Abortando. Re-ejecuta el Script 03 (v2) para regenerar el catálogo.
```
Solución: tu catálogo fue generado con una versión antigua del Script 03. Volvé a correr el Script 03 con los datos originales para regenerar el catálogo con todas las columnas requeridas.

---

### Script 07 — Simulador N-cuerpos

**Pregunta científica:** ¿Cómo van a evolucionar estas estrellas bajo la gravedad de la Vía Láctea en el próximo millón de años?

**Input:** Catálogo con velocidades UVW (salida del Script 04).
**Output:** Archivo HDF5 con las trayectorias completas de las estrellas.

**Concepto clave — Integrador Leapfrog:**
Para simular cómo se mueve una estrella bajo una fuerza gravitacional, necesitás un integrador numérico: un algoritmo que avanza la posición y velocidad en pequeños pasos de tiempo. El Script 07 usa el **Leapfrog KDK** (Kick-Drift-Kick), que tiene una propiedad especial: conserva energía perfectamente a largo plazo. Otros integradores como Runge-Kutta van acumulando pequeños errores que eventualmente desvían las órbitas. Leapfrog no.

El nombre "KDK" describe los tres pasos de cada iteración:
1. **Kick** (½ paso): actualizar velocidades con la aceleración actual
2. **Drift** (paso completo): mover las posiciones con las velocidades actualizadas
3. **Kick** (½ paso): corregir las velocidades con la nueva aceleración

```bash
python 07_simulador_n_cuerpos.py \
  --input data/catalogo_con_uvw.parquet \
  --output data/simulacion.h5 \
  --n_stars 200 \
  --years 1000000 \
  --dt_days 365
```

**Output esperado:**
```
🪐 SIMULADOR N-CUERPOS 3D GAIA (v4 — Numba)
✔ Numba detectado | Threads: 8
🔬 Validando integrador contra solución analítica...
  ✔ PASS  v_circ_analytical_plausible (val=197.3)
  ✔ PASS  orbit_closure_after_10_orbits
  ✔ PASS  energy_drift_leapfrog (val=3.2e-07)

Integrando Leapfrog (Numba)... ██████████████████████ 100% • 0:04:21
Completado en 261.3s (10,497 pasos/s)
E inicial: -2.847e+42 J | |ΔE/E₀|: 4.1e-07

💾 HDF5 guardado: data/simulacion.h5
```

**Error común:** La simulación corre muy lento (más de 30 minutos para 200 estrellas).
Solución: instalá Numba con `pip install numba`. La primera vez que corre, Numba compila el código (tarda ~1-2 minutos). Después es 50-100x más rápido.

---

### Script 08 — Reporte científico PDF

**Pregunta científica:** ¿Cómo presento estos resultados de forma clara y reproducible?

**Input:** Catálogo clasificado + (opcional) archivo HDF5 de simulación.
**Output:** PDF de 10 páginas con estadísticas, gráficos y metadata de reproducibilidad.

**Concepto clave — Reproducibilidad científica:**
El reporte incluye el **hash SHA256** del archivo de datos de entrada. Esto significa que cualquier persona que tenga el mismo archivo puede verificar que el reporte fue generado con exactamente esos datos y no otros. Es una práctica estándar en ciencia computacional.

```bash
python 08_reporte_cientifico.py \
  --catalog data/catalogo_clasificado.parquet \
  --simulation data/simulacion.h5 \
  --output reporte_gaia_final.pdf \
  -v
```

**Output esperado:**
```
======================================================================
 📊 REPORTE CIENTÍFICO GAIA (Script 08)
======================================================================
📂 Cargando catálogo: data/catalogo_clasificado.parquet
  → 9,341 estrellas | 28 columnas
  → Cobertura V_radial: 34.2%
  → Cobertura UVW: 34.1%
📄 Generando página: Portada...
📄 Generando página: Estadísticas...
📄 Generando página: CMD...
...
✅ Reporte generado: reporte_gaia_final.pdf
   Tamaño: 4.82 MB | Páginas: 10
```

**Error común:**
```
DataQualityError: No se encuentra el catálogo: data/catalogo_clasificado.parquet
```
Solución: verificá que el path al catálogo sea correcto. El Script 08 necesita el catálogo que produce el Script 06 (no el del 03 directamente).

---

## BLOQUE 3 — La física detrás del pipeline 🌌

### Qué es la paralaje y cómo la mide Gaia

Extendé el brazo y mirá tu dedo con el ojo derecho cerrado, luego con el izquierdo. El dedo parece moverse contra el fondo. Eso es paralaje: el desplazamiento aparente de un objeto cercano visto desde dos posiciones diferentes.

Gaia hace exactamente lo mismo pero con estrellas. Mientras la Tierra orbita alrededor del Sol, Gaia observa cada estrella desde posiciones separadas hasta 2 UA (la distancia Tierra-Sol). Las estrellas cercanas se "mueven" ligeramente respecto a las lejanas. Ese movimiento se mide en **mili-arco-segundos (mas)**.

La relación distancia-paralaje es:
```
d [parsecs] = 1 / π [arco-segundos] = 1000 / π [mas]
```

Una estrella a 100 parsecs tiene una paralaje de 10 mas. Una a 1000 parsecs tiene 1 mas. Gaia puede medir paralajes con precisión de 0.02 mas, lo que permite distancias confiables hasta ~5000 parsecs.

### Qué son las velocidades UVW y por qué importan

Las velocidades UVW describen el movimiento de una estrella en un sistema de coordenadas galáctico centrado en el Sol:

| Componente | Dirección | Valor típico (disco delgado) |
|---|---|---|
| **U** | Hacia el centro galáctico (positivo = adentro) | −20 a +20 km/s |
| **V** | En la dirección de rotación de la Galaxia | −30 a +10 km/s |
| **W** | Perpendicular al disco (hacia el polo norte galáctico) | −15 a +15 km/s |

Analizando estas velocidades podés distinguir diferentes poblaciones estelares: el **disco delgado** (estrellas jóvenes, velocidades bajas), el **disco grueso** (estrellas viejas, velocidades moderadas) y el **halo** (estrellas muy viejas, velocidades altísimas e irregulares). También podés detectar **corrientes estelares**: grupos de estrellas que se mueven juntas porque comparten un origen común.

### Cómo funciona un simulador N-cuerpos

Imaginate N bolas de billar en el espacio, pero en vez de colisionar se atraen mutuamente por gravedad. El problema es que no existe una fórmula exacta para calcular dónde van a estar después de un tiempo T cuando N > 2 (es el famoso "problema de los tres cuerpos").

La solución es numérica: dividir el tiempo en pasos muy pequeños (dt) y para cada paso:
1. Calcular la fuerza gravitacional de cada estrella sobre todas las demás (O(N²) operaciones)
2. Actualizar las velocidades y posiciones en consecuencia

Después de miles o millones de pasos, tenés la trayectoria completa. El Script 07 además incluye el **potencial galáctico externo** (la gravedad del disco, el bulbo y el halo de la Vía Láctea), que domina completamente el movimiento a escalas de kpc.

### Cómo leer el Diagrama Hertzsprung-Russell

El HR es uno de los gráficos más informativos de toda la astronomía. Aprender a "leerlo" te da inmediatamente información física de una estrella:

```
Magnitud
Absoluta
(más brillante)
    -5  ●  ●  ●    ← Supergigantes
     0     ● ●     ← Gigantes
    +5  ○ ● ●      ← Sol está aquí (G2V)
   +10    ● ●●     ← Enanas K/M
   +15  ○           ← Enanas blancas
         (más azul → más rojo)
         Color BP-RP →
```

- Estrellas **arriba**: más brillantes (gigantes, supergigantes)
- Estrellas **abajo**: más tenues (enanas)
- Estrellas **a la izquierda**: más calientes (azules, tipo O/B/A)
- Estrellas **a la derecha**: más frías (rojas, tipo K/M)
- La **diagonal** de arriba-izquierda a abajo-derecha es la Secuencia Principal: donde viven todas las estrellas que están quemando hidrógeno en su núcleo, incluido el Sol.

---

## BLOQUE 4 — Cómo extender el proyecto 🚀

Este pipeline es un punto de partida, no un límite. Acá van tres ideas concretas de proyectos que podés construir encima de lo que ya está.

### Idea 1 — Estudiar un tipo de estrella específico

**Objetivo:** En vez de descargar todo el campo de Gaia, filtrar la consulta ADQL para traer solo, por ejemplo, enanas blancas o estrellas de tipo A.

**Archivo a modificar:** `02_consulta_gaia_adql.py`

**Qué cambiar (~línea 224):**
```python
# Original — trae todo el campo
adql = (
    f"SELECT TOP {args.limit} {cols} "
    f"FROM gaiadr3.gaia_source "
    f"WHERE 1=CONTAINS(POINT('ICRS', ra, dec), "
    f"BOX('ICRS', {ra_c:.8f}, {dec_c:.8f}, {w:.8f}, {h:.8f}))"
)

# Modificado — solo estrellas azules calientes (BP-RP < 0.2, posibles tipo A/B)
adql = (
    f"SELECT TOP {args.limit} {cols} "
    f"FROM gaiadr3.gaia_source "
    f"WHERE 1=CONTAINS(POINT('ICRS', ra, dec), "
    f"BOX('ICRS', {ra_c:.8f}, {dec_c:.8f}, {w:.8f}, {h:.8f})) "
    f"AND bp_rp < 0.2 AND phot_g_mean_mag < 18"
)
```

**Nivel de dificultad:** Bajo. Solo hay que conocer las columnas de Gaia DR3 y las condiciones ADQL.

---

### Idea 2 — Simular el Sistema Solar como test

**Objetivo:** En vez de usar datos reales de Gaia, crear un "catálogo" artificial con los planetas del Sistema Solar y verificar que las órbitas son estables.

**Archivo a modificar:** `07_simulador_n_cuerpos.py`

**Qué cambiar (~línea 373, función `load_catalog`):**
```python
# En vez de leer el Parquet, retornar datos del Sistema Solar
def load_solar_system():
    """Posiciones y velocidades de los planetas en coordenadas galactocéntricas."""
    import astropy.constants as const

    # Posición en metros (distancias al Sol)
    pos = np.array([
        [0.0, 0.0, 0.0],                  # Sol
        [1.0 * AU_TO_M, 0.0, 0.0],        # Tierra (1 UA)
        [5.2 * AU_TO_M, 0.0, 0.0],        # Júpiter (5.2 UA)
    ])
    # Velocidad en m/s (velocidades orbitales)
    vel = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 29780.0, 0.0],   # Tierra: 29.78 km/s
        [0.0, 13070.0, 0.0],   # Júpiter: 13.07 km/s
    ])
    masses = np.array([M_SUN, 5.97e24, 1.90e27])
    return pos, vel, masses, np.array(['Sol', 'Tierra', 'Jupiter'])
```

**Nivel de dificultad:** Medio. Requiere conocer las unidades del simulador (metros, kg) y tener valores de referencia de los planetas.

---

### Idea 3 — Análisis de un cúmulo estelar

**Objetivo:** Modificar la consulta ADQL para apuntar a un cúmulo estelar conocido (por ejemplo, las Pléyades o el cúmulo de las Híades) y estudiar la cinemática de sus miembros.

**Archivos a modificar:** `02_consulta_gaia_adql.py` y `05_generador_cmd.py`

**En el Script 02 (~línea 165, argumento `--margin-deg`):**
Las Pléyades están en RA=56.75°, Dec=+24.12°. Podés crear un FITS sintético con ese WCS, o simplemente hardcodear la consulta directamente:

```python
# Consulta directa para las Pléyades (radio 2°)
adql = """
SELECT source_id, ra, dec, parallax, parallax_error,
       pmra, pmdec, phot_g_mean_mag, bp_rp, teff_gspphot
FROM gaiadr3.gaia_source
WHERE 1=CONTAINS(
  POINT('ICRS', ra, dec),
  CIRCLE('ICRS', 56.75, 24.12, 2.0)
)
AND parallax BETWEEN 5 AND 10
"""
```

El filtro `parallax BETWEEN 5 AND 10` (100-200 pc) selecciona solo las estrellas miembros del cúmulo y excluye las de fondo.

**En el Script 05**, el CMD del cúmulo va a mostrar una Secuencia Principal mucho más limpia y estrecha que la de un campo general, porque todas las estrellas del cúmulo tienen la misma edad y composición química.

**Nivel de dificultad:** Medio-alto. Requiere entender la relación entre paralaje y pertenencia a un cúmulo.

---

## BLOQUE 5 — Glosario 📖

### ADQL
Astronomical Data Query Language. Lenguaje de consulta similar a SQL diseñado para bases de datos astronómicas. Tiene funciones especiales como `CONTAINS`, `CIRCLE` y `BOX` para consultas espaciales en coordenadas celestes.

### Cinemática estelar
El estudio del movimiento de las estrellas sin considerar las fuerzas que lo causan. Describe posición, velocidad y trayectoria. La dinámica estelar sí considera las fuerzas (gravedad) y es lo que hace el Script 07.

### Clasificación espectral
Sistema de categorización de estrellas por temperatura superficial. Las clases principales son O, B, A, F, G, K, M en orden de más caliente (>30.000 K) a más fría (~3.000 K). El Sol es de tipo G2.

### Diagrama HR
Diagrama Hertzsprung-Russell. Gráfico que enfrenta la luminosidad o magnitud absoluta de una estrella contra su temperatura o color. Revela inmediatamente el tipo y estado evolutivo de la estrella.

### Enriquecimiento de datos
Proceso de tomar datos brutos (solo lo que mide el instrumento) y calcular cantidades físicas derivadas: distancias a partir de paralajes, magnitudes absolutas a partir de aparentes, velocidades en km/s a partir de movimientos propios en mas/año.

### Gaia DR3
Tercer catálogo de datos del satélite Gaia de la ESA (2022). Contiene astrometría, fotometría y espectroscopía de ~1.500 millones de objetos con una precisión sin precedentes en posiciones, paralajes y movimientos propios.

### Leapfrog
Integrador numérico simpléctico (conserva energía mecánica) usado para simular órbitas gravitacionales. Alterna medios pasos de velocidad y pasos completos de posición. Ideal para simulaciones largas porque no acumula errores de energía.

### Magnitud absoluta
La magnitud aparente que tendría una estrella si estuviera a exactamente 10 parsecs de distancia. Permite comparar el brillo real de estrellas sin importar su distancia al observador. El Sol tiene magnitud absoluta G ≈ +4.67.

### N-cuerpos
Problema físico de calcular la evolución de N partículas que interactúan mutuamente por gravedad. Para N > 2 no tiene solución analítica exacta y se resuelve numéricamente. La complejidad computacional naive es O(N²) por paso de tiempo.

### Numba JIT
Just-In-Time compiler para Python. Compila funciones Python marcadas con `@njit` a código máquina nativo la primera vez que se llaman, con velocidades comparables a C o Fortran. Esencial para el bucle O(N²) del simulador.

### Paralaje
Desplazamiento angular aparente de una estrella cercana visto desde dos posiciones diferentes de la Tierra (extremos opuestos de la órbita anual). Se mide en arco-segundos o mili-arco-segundos. Inversamente proporcional a la distancia: 1 parsec ↔ 1 arco-segundo de paralaje.

### Parsec
Unidad de distancia astronómica. Equivale a la distancia a la que una estrella tendría una paralaje de 1 arco-segundo. 1 parsec ≈ 3.26 años-luz ≈ 3.086 × 10¹⁶ metros. La estrella más cercana al Sol (Proxima Centauri) está a 1.3 parsecs.

### Secuencia Principal
La banda diagonal del Diagrama HR que contiene las estrellas que están en la etapa principal de su vida: quemando hidrógeno en su núcleo por fusión nuclear. El Sol lleva ~4.600 millones de años en esta etapa y seguirá otros ~5.000 millones de años más.

### TAP
Table Access Protocol. Protocolo estándar del IVOA (International Virtual Observatory Alliance) para acceder a bases de datos astronómicas remotas mediante consultas ADQL. Los servidores TAP de Gaia, SDSS, 2MASS y cientos de observatorios son accesibles con astroquery.

### UVW
Las tres componentes de la velocidad estelar en el sistema galactocéntrico local: U (hacia el centro galáctico), V (en la dirección de rotación de la Galaxia), W (perpendicular al disco). Son fundamentales para identificar a qué población estelar pertenece una estrella.

---

## BLOQUE 6 — Créditos y contexto 🙏

Este pipeline fue desarrollado por **Juan Galaz** como proyecto personal de investigación en astrofísica computacional. La idea era crear una herramienta completa que fuera desde el dato crudo del satélite Gaia hasta una simulación dinámica y un reporte publicable, sin depender de plataformas externas ni servidores propios.

Lo que empezó como un script de descarga creció hasta convertirse en ocho scripts interconectados con un módulo fundacional, un sistema de tests físicos y documentación pensada para que otra persona pueda usarlo sin que yo esté presente explicando qué hace cada función.

El proyecto fue donado a **Into Space Academia** para que estudiantes de astrofísica de todas las edades puedan explorar datos reales del universo, entender cómo se procesan científicamente y, si les da la curiosidad, modificarlo y mejorarlo.

Si encontrás un bug, querés sugerir algo, o simplemente querés contar qué proyecto hiciste con esto, bienvenido seas.

- 🌐 Web: [cienciaestelar.cl](https://cienciaestelar.cl)
- 💻 GitHub: [github.com/juangalaz](https://github.com/juangalaz)
- 🏫 Into Space Academia: educación en astronomía y astrofísica para Chile y Latinoamérica

---

*"Los datos de Gaia son públicos. Las herramientas para analizarlos también deberían serlo."*

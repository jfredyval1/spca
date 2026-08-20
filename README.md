# sPCA hidrogeoquímico — Cuenca Baja del Río Bogotá (CBRB)

Análisis de Componentes Principales **espacial** (sPCA, Jombart 2008) aplicado a los
muestreos hidrogeoquímicos del proyecto CBRB, con el **ACP clásico como control** en
todos los diagnósticos.

## Qué hace

El ACP clásico busca las combinaciones lineales de máxima **varianza** y no sabe nada
del espacio. El sPCA optimiza en cambio el criterio compuesto

```
C(x) = var(x) · I(x)        (I = índice de Moran sobre una red de vecindad)
```

de modo que sus autovalores pueden ser positivos (estructuras **globales**: gradientes,
parches, bloques) o negativos (estructuras **locales**: contraste entre vecinos).

La pregunta que ordena todo el documento es **qué gana el análisis al meter el espacio
en el criterio de optimización**: dónde el sPCA reproduce al ACP, dónde se aparta de él,
y cuánto de esa diferencia es propiedad de los datos y no de la red que se eligió para
mirarlos.

## Diseño de la comparación

Tres análisis sobre **exactamente la misma matriz** (`hq.spca.nb$tab`, ya centrada y
escalada, así que el ACP se ajusta con `center = FALSE, scale = FALSE`):

| Análisis | Red de conexión | Papel |
|---|---|---|
| ACP clásico (`dudi.pca`) | — | referencia sin espacio |
| sPCA Delaunay | `chooseCN(type = 1)` | red por defecto, solo geometría euclidiana |
| sPCA manual | grafo de Gabriel editado | fallas geológicas como barreras de flujo |

Cada diagnóstico —autovalores, `screeplot`, `global.rtest` / `local.rtest`,
descomposición `var × I`, cargas, clasificación jerárquica y coherencia espacial de los
grupos— se corre **con las dos redes**, para separar lo que aportan las fallas de lo que
ya estaba en los datos. Los signos de los ejes se alinean antes de restar cargas, porque
la orientación de un eje factorial es arbitraria.

## Los datos

60 puntos de agua subterránea (de 63 con química; se descartan 3 superficiales sin
geometría) y 12 variables: 10 concentraciones (Ca, Mg, Na, K, HCO₃, SO₄, Cl, Si, Mn, As),
pH y T.

- Las concentraciones son **datos composicionales**: se homogeneizan unidades a mg/L, se
  reemplazan los ceros por el método multiplicativo (Martín-Fernández et al., 2003) y se
  transforman con **CLR** (razón logarítmica centrada). Esto corrige la asimetría *y*
  elimina el efecto de dilución, de forma que ningún eje se lee como «mineralización
  total».
- pH y T entran sin transformar; la comparabilidad la resuelve `scale = TRUE`.
- El bloque CLR deja la matriz singular por construcción (rango 11 de 12).
- Coordenadas proyectadas a EPSG:9377 en kilómetros: la red se define por distancias
  euclidianas y no admite grados.

## Resultados principales

- **Dos componentes globales significativos con ambas redes** (manual: λ₁ = 1.903,
  λ₂ = 1.195; Delaunay: λ₁ = 1.541, λ₂ = 0.877). La dimensionalidad retenida no depende
  de la red, solo su nitidez.
- **Autocorrelación global significativa** (`global.rtest`: p = 0.0004 con Delaunay,
  p = 0.0007 con la red manual) y **estructura local descartada** (`local.rtest`:
  p = 0.89 y p = 0.63).
- **El espacio aporta de forma desigual.** Eje 1: el sPCA es casi idéntico al ACP
  (|r| = 0.98; +8.3 % de I por −3.5 % de varianza) — la varianza dominante *ya era*
  espacial. Eje 2: cede 10.4 % de varianza y gana **83.1 % de autocorrelación**
  (I: 0.299 → 0.548).
- **El eje 2 es una sustitución identificable**: el ACP lo apoya en As (carga 0.48), los
  dos sPCA lo sustituyen por HCO₃ (0.04 → 0.44), concentrando en un solo componente la
  firma silíceo-bicarbonatada que el ACP repartía entre dos ejes.
- **Clasificación (Ward, k = 5)**: 43.8 % de enlaces intragrupo frente a 20.7 % esperado
  por azar (p = 0.001) y 35.9 % de la partición del ACP; con Delaunay la ventaja casi
  desaparece (39.6 % vs 37.9 %), así que se reporta como condicionada a la red.

**Lectura de conjunto**: el aporte del sPCA aquí no está en producir coordenadas
distintas, sino en **someter a prueba la hipótesis espacial** —tests de permutación que
un ACP no puede dar— y en afinar un componente interpretable. Las fallas hacen la señal
más nítida, pero no la crean.

## Contenido del repositorio

```
cbrb_aplicacion_spca.qmd   documento principal: teoría + aplicación CBRB + mapa interactivo
cbrb_aplicacion_spca.html  render autocontenido (embed-resources)
spca_tutorial.qmd          tutorial de sPCA en español (Jombart 2008 / viñeta adegenet),
                           con los ejemplos genéticos originales
spca.qgz                   proyecto QGIS de apoyo (edición de la red de vecindad)
```

`data/`, `sql/` y `pdf/` no se versionan (ver `.gitignore`).

## Reproducir

Requiere R con: `adegenet`, `spdep`, `ade4`, `adespatial`, `sf`, `dplyr`, `RPostgres`,
`getPass`, `cluster`, `interp`, `leaflet`, `xml2`, más Quarto.

La química se lee de `data/acp_mediana.rds` (reconstruible desde la BDH con
`sql/acp_mediana.sql`) y las geometrías se consultan por conexión directa. Las
credenciales se leen de un archivo `.env` junto al `.qmd` y **no** se versionan:

```
IP_SERVER=...
PORT=...
DB_NAME=...
DB_USER=...
DB_PASS=...      # opcional; si falta, se pide por consola con getPass()
```

```sh
quarto render cbrb_aplicacion_spca.qmd
```

## Referencias

- Jombart, T., Devillard, S., Dufour, A.-B. & Pontier, D. (2008). *Revealing cryptic
  spatial patterns in genetic variability by a new multivariate method.* Heredity 101.
- Viñeta *adegenet* 2.1.6 — «spatial Principal Component Analysis».
- Martín-Fernández, J. A., Barceló-Vidal, C. & Pawlowsky-Glahn, V. (2003). *Dealing with
  zeros and missing values in compositional data sets.* Mathematical Geology 35.

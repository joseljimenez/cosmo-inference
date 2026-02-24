# Cosmological Parameter Inference Pipeline

Este pipeline end-to-end calcula la estadística \chi^2 total combinando distintos conjuntos de datos cosmológicos con el objetivo de evaluar la compatibilidad de un modelo teórico con las observaciones. En una etapa posterior, se utiliza un análisis MCMC para explorar el espacio de parámetros y obtener restricciones sobre los parámetros del modelo cosmológico mediante técnicas de ajuste estadístico. El 

### Datasets utilizados
    - Supernovas (Union3)
    - BAO (Baryon Acoustic Oscillations)
    - Cosmic Chronometers H(z)

### Como usar
- Instalar la librerías de python: numpy, astropy, pandas, scipy
- Descargar los datos
|`./data_download.sh`
- Evaluar los likelihoods:
`python evaluation.py`
- Para hacer el sampling MCMC, abrir y ejecutar el notebook
`emcee.ipynb`

### `data_download.sh`
Este fichero descarga datasets: Union3, BAO-DESI DR2 de repositorios oficiales

### cosmology.py
Contiene todas las ecuaciones de background cosmológico, ecuaciones de distancia implementadas mediante la ecuación de estado del Gas de Chaplygin.

### `likelihood.py`
Define y calcula funciones \chi^2 individuales para cada conjunto de datos (Union3, BAO, Cronómetros Cósmicos) y luego regresa la función \chi^2_total.


### Generalized Chaplygin Gas (GCG)

Ecuación de Estado

p = -A / ρ^α

Evolución de la densidad de energía:

ρ_GCG(z) = ρ0 [ As + (1 - As)(1+z)^{3(1+α)} ]^{1/(1+α)}

Parámetros libres:
- H0
- As
- α

ΛCDM se recupera en el límite α → 0.

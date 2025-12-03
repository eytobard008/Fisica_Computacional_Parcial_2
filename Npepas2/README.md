# N_Particulas

Simulador de partículas en una caja bidimensional con análisis de presión y distribución de velocidades.

## Estructura del Proyecto

```markdown
N_Particulas/
├── include/
│   └── Npepas.h           # Interfaz principal
├── src/
│   ├── main.cpp           # Programa principal
│   └── Npepas.cpp         # Implementación
├── scripts/               # Scripts gnuplot
├── results/               # Output de simulaciones
└── Makefile               # Sistema de compilación
├── documents/
│   ├── Kuramoto.tex
│   └── Kuramoto.pdf
└── README.md
```

## Características Principales

- Simulación de partículas en caja 2D (20m × 10m por defecto)
- Dos métodos de presión: por teoría cinética y por colisiones con paredes
- Dos modos de inicialización:
  * Malla regular: partículas organizadas en cuadrícula ordenada
  * Aleatorio: partículas colocadas aleatoriamente sin solapamiento
- Optimización opcional con malla espacial para detección de colisiones
- Análisis estadístico: distribución de velocidades Maxwell-Boltzmann
- Visualizaciones: animaciones y gráficas de análisis
- Métricas de rendimiento: conteo de cálculos para comparar métodos


## Compilación y Ejecución

```bash
# Compilar el proyecto
make

# Ejecutar la simulación (modo interactivo)
make run

# Ver animaciones y gráficas generadas
make view

# Limpiar ejecutable, resultados y scripts
make clean
```

## Parámetros de Configuración

### Parámetros solicitados al usuario:

- Radio de partículas (R) en metros
- Número de partículas (N)
- Tiempo máximo de simulación (tmax) en segundos
- Velocidad media inicial en m/s
- Dispersión de velocidad (desviación estándar) en m/s
- Optimización con malla (1=Sí, 0=No) para detección de colisiones


### Parámetros modificables en main.cpp:

```bash
double masa = 1.0;           // Masa de cada partícula (kg)
double W = 20.0, H = 10.0;   // Dimensiones de la caja (m)
double dt = 0.01;            // Paso de tiempo (s)
bool contarCalculos = true;  // Mostrar métricas de rendimiento
bool usarMalla = false;      // Disposición inicial (malla regular o aleatoria)
```


## Modos de Operación

### Disposición Inicial de Partículas
- **Malla regular**: Partículas organizadas en filas/columnas (modificar `usarMalla = true` en main.cpp)
- **Aleatoria**: Colocación aleatoria sin solapamiento (modificar `usarMalla = false` en main.cpp)

### Detección de Colisiones
- **Fuerza bruta**: Verifica todos los pares de partículas - seleccionar "No" (0) en optimización
- **Malla espacial**: Usa división espacial para reducir cálculos - seleccionar "Sí" (1) en optimización


## Visualizaciones Generadas

- `Npepas.gif`: Animación del sistema de partículas
- `histograma.gif`: Evolución de la distribución de velocidades
- `histograma_zoom.gif`: Vista detallada de la distribución
- `energia.png`: Gráfica de energía cinética total
- `presion.png`: Comparación de presiones teórica y medida

## Parámetros de Entrada

- Radio de las bolas (R)
- Número de partículas (N)
- Dimensiones de la caja (W, H)
- Tiempo máximo de simulación (tmax)
- Velocidad media inicial
- Desviación estándar de velocidad

## Límites Recomendados

El programa calcula automáticamente los límites recomendados basados en el empaquetamiento hexagonal máximo teórico. Para un funcionamiento óptimo, se recomienda:

- No exceder el 75% del empaquetamiento máximo
- Verificar que el diámetro de las partículas sea menor que las dimensiones de la caja

## Dependencias

- C++17
- Gnuplot (para visualizaciones)


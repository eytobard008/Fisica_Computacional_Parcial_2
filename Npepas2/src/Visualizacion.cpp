/**
 * @file Visualizacion.cpp
 * @brief Implementación de las funciones de visualización y gráficas
 * @author Esteban Tobar
 * @date 29 de octubre de 2025
 */

#include "../include/Npepas.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>

using namespace std;

// Función auxiliar para ejecutar comandos del sistema
inline void ejecutarComando(const string& comando) {
    int resultado = system(comando.c_str());
    (void)resultado;
}

// ==========================
// Generación de animación principal
// ==========================

void GenerarAnimacion(const string &archivoDat, const string &archivoPresion, double W, double H, int N, int pasos, double dt, double R, bool usarMallaOpt, double velocidadMedia, double dispersionVelocidad, const string &nombreGif) {
    // Leer datos de presión para mostrar en la animación
    vector<double> P_vel_values, P_col_values;
    ifstream presionFile(archivoPresion);
    string linea;

    getline(presionFile, linea);
    while (getline(presionFile, linea)) {
        stringstream ss(linea);
        double t, p_vel, p_col;
        if (ss >> t >> p_vel >> p_col) {
            P_vel_values.push_back(p_vel);
            P_col_values.push_back(p_col);
        }
    }
    presionFile.close();

    // Asegurar que tenemos datos para todos los pasos
    if (P_vel_values.size() != static_cast<size_t>(pasos)) {
        while (P_vel_values.size() < static_cast<size_t>(pasos)) {
            P_vel_values.push_back(0.0);
            P_col_values.push_back(0.0);
        }
    }

    // Crear script de gnuplot para la animación
    ofstream gnu("scripts/Npepas.gnu");
    if (!gnu.is_open()) {
        cerr << "Error creando script gnuplot\n";
        return;
    }

    gnu << "set terminal gif animate delay 5 size 800,500\n";
    gnu << "set output 'results/" << nombreGif << "'\n";
    gnu << "set xrange [0:" << W << "]\n";
    gnu << "set yrange [0:" << H << "]\n";
    gnu << "set size ratio " << H / W << "\n";
    gnu << "unset key\n";
    gnu << "set style fill solid\n";
    gnu << "set margins 0,0,0.1,0.1\n";

    gnu << "R = " << R << "\n";

    // Dibujar bordes de la caja
    gnu << "set border front lw 2\n";
    gnu << "set object 1 rectangle from 0,0 to " << W << "," << H << " front fillstyle empty border lc 'black' lw 3\n";

    // Dibujar líneas de malla si se usó optimización
    if (usarMallaOpt) {
        double tamanoCeldaVel = (velocidadMedia + dispersionVelocidad) * dt * 3.0;
        double tamanoCeldaRadio = 2.0 * R;
        double tamanoCelda = max(tamanoCeldaVel, tamanoCeldaRadio);
        
        for (double x = tamanoCelda; x < W; x += tamanoCelda) {
            gnu << "set arrow from " << x << ",0 to " << x << "," << H << " nohead lc 'gray' lw 0.5 front\n";
        }
        for (double y = tamanoCelda; y < H; y += tamanoCelda) {
            gnu << "set arrow from 0," << y << " to " << W << "," << y << " nohead lc 'gray' lw 0.5 front\n";
        }
    }

    // Almacenar valores de presión en arrays de gnuplot
    gnu << "array P_vel[" << pasos << "]\n";
    gnu << "array P_col[" << pasos << "]\n";
    for (int i = 0; i < pasos; i++) {
        gnu << "P_vel[" << i + 1 << "] = " << P_vel_values[i] << "\n";
        gnu << "P_col[" << i + 1 << "] = " << P_col_values[i] << "\n";
    }

    // Bucle para generar cada frame de la animación
    gnu << "do for [i=1:" << pasos << "] {\n";
    gnu << "  set label 1 sprintf('t = %.2f s', (i-1)*" << dt << ") at screen 0.02,0.98 front\n";
    gnu << "  set label 2 sprintf('P_{vel} = %.3f N/m', P_vel[i]) at screen 0.02,0.95 front\n";
    gnu << "  set label 3 sprintf('P_{col} = %.3f N/m', P_col[i]) at screen 0.02,0.92 front\n";
    gnu << "  \n";
    gnu << "  plot \\\n";

    // Graficar cada partícula como un círculo
    for (int j = 0; j < N; j++) {
        int col_x = 2 + 4 * j;
        int col_y = 3 + 4 * j;

        gnu << "    '" << archivoDat << "' every ::i-1::i-1 using " << col_x << ":" << col_y << ":(R) with circles lc 'black' notitle";

        if (j < N - 1) {
            gnu << ", \\\n";
        } else {
            gnu << "\n";
        }
    }

    gnu << "  unset label 1\n";
    gnu << "  unset label 2\n";
    gnu << "  unset label 3\n";
    gnu << "}\n";

    gnu.close();

    ejecutarComando("gnuplot scripts/Npepas.gnu");
}

// ==========================
// Generación de histogramas animados
// ==========================

void GenerarHistogramaAnimado(const string &archivoDat, int N, int pasos, double dt, double velocidadMedia, const string &nombreGif, const string &nombreGifZoom) {
    ifstream archivo(archivoDat);
    if (!archivo.is_open()) {
        cerr << "Error: no se pudo abrir " << archivoDat << " para generar histograma\n";
        return;
    }

    const int num_bins = 100;
    double min_velocidad = 0.0;
    double max_velocidad = velocidadMedia * 4;
    int max_frecuencia = 0;

    vector<double> velocidades(N);
    ofstream allHistograms("results/histograma_completo.dat");

    // Leer datos y calcular histogramas para cada frame
    for (int frame = 0; frame < pasos; frame++) {
        string linea;
        getline(archivo, linea);
        stringstream ss(linea);
        double t;
        ss >> t;

        // Extraer velocidades de todas las partículas
        for (int i = 0; i < N; i++) {
            double x, y, vx, vy;
            ss >> x >> y >> vx >> vy;
            velocidades[i] = sqrt(vx * vx + vy * vy);
        }

        // Calcular histograma
        vector<int> histograma(num_bins, 0);
        double bin_width = (max_velocidad - min_velocidad) / num_bins;

        for (int i = 0; i < N; i++) {
            double v = velocidades[i];
            int bin = static_cast<int>((v - min_velocidad) / bin_width);
            if (bin >= 0 && bin < num_bins) {
                histograma[bin]++;
                if (histograma[bin] > max_frecuencia) {
                    max_frecuencia = histograma[bin];
                }
            }
        }

        // Guardar histograma para este frame
        allHistograms << "# Frame " << frame << "\n";
        for (int i = 0; i < num_bins; i++) {
            double v = min_velocidad + (i + 0.5) * bin_width;
            allHistograms << v << " " << histograma[i] << "\n";
        }
        allHistograms << "\n\n";
    }

    allHistograms.close();
    archivo.close();

    // Calcular distribución teórica de Maxwell-Boltzmann
    double kT_over_m = 2.0 * velocidadMedia * velocidadMedia / M_PI;
    double bin_width = (max_velocidad - min_velocidad) / num_bins;

    // Encontrar valor máximo para el zoom
    double max_mb_value = 0.0;
    for (int i = 0; i < num_bins; i++) {
        double v = min_velocidad + (i + 0.5) * bin_width;
        double mb_val = (v / kT_over_m) * exp(-v*v/(2.0 * kT_over_m)) * N * bin_width;
        if (mb_val > max_mb_value) max_mb_value = mb_val;
    }

    double zoom_y_max = max_mb_value * 1.5;

    // Crear script para histograma animado completo
    ofstream gnu("scripts/histograma_animado.gnu");
    if (gnu.is_open()) {
        gnu << "set terminal gif animate delay 5 size 800,400\n";
        gnu << "set output 'results/" << nombreGif << "'\n";
        gnu << "set xlabel '|v| (m/s)'\n";
        gnu << "set ylabel 'Número de partículas'\n";
        gnu << "set grid\n";
        gnu << "set yrange [0:" << max_frecuencia << "]\n";
        gnu << "set xrange [" << min_velocidad << ":" << max_velocidad << "]\n";
        gnu << "set style fill solid\n";

        gnu << "kT_over_m = 2.0 * " << velocidadMedia << " * " << velocidadMedia << " / pi\n";
        gnu << "bin_width = " << bin_width << "\n";
        gnu << "maxwell(x) = (x / kT_over_m) * exp(-x*x/(2.0 * kT_over_m)) * " << N << " * bin_width\n";

        gnu << "do for [i=0:" << pasos - 1 << "] {\n";
        gnu << "  set title sprintf('t = %.2f s', i*" << dt << ")\n";
        gnu << "  plot 'results/histograma_completo.dat' index i with impulses lw 3 lc 'blue' title 'Distribución medida', \\\n";
        gnu << "       maxwell(x) with lines lw 2 lc 'red' title 'Maxwell-Boltzmann'\n";
        gnu << "}\n";
        gnu.close();
        
        ejecutarComando("gnuplot scripts/histograma_animado.gnu");
    }

    // Crear script para histograma animado con zoom
    ofstream gnuZoom("scripts/histograma_animado_zoom.gnu");
    if (gnuZoom.is_open()) {
        gnuZoom << "set terminal gif animate delay 5 size 800,400\n";
        gnuZoom << "set output 'results/" << nombreGifZoom << "'\n";
        gnuZoom << "set xlabel '|v| (m/s)'\n";
        gnuZoom << "set ylabel 'Número de partículas'\n";
        gnuZoom << "set grid\n";
        gnuZoom << "set yrange [0:" << zoom_y_max << "]\n";
        gnuZoom << "set xrange [0:" << max_velocidad * 0.6 << "]\n";
        gnuZoom << "set style fill solid\n";

        gnuZoom << "kT_over_m = 2.0 * " << velocidadMedia << " * " << velocidadMedia << " / pi\n";
        gnuZoom << "bin_width = " << bin_width << "\n";
        gnuZoom << "maxwell(x) = (x / kT_over_m) * exp(-x*x/(2.0 * kT_over_m)) * " << N << " * bin_width\n";

        gnuZoom << "do for [i=0:" << pasos - 1 << "] {\n";
        gnuZoom << "  set title sprintf('t = %.2f s (Zoom MB)', i*" << dt << ")\n";
        gnuZoom << "  plot 'results/histograma_completo.dat' index i with impulses lw 3 lc 'blue' title 'Distribución medida', \\\n";
        gnuZoom << "       maxwell(x) with lines lw 2 lc 'red' title 'Maxwell-Boltzmann'\n";
        gnuZoom << "}\n";
        gnuZoom.close();
        
        ejecutarComando("gnuplot scripts/histograma_animado_zoom.gnu");
    }
}

// ==========================
// Generación de gráficas estáticas
// ==========================

void GenerarGraficaEnergia(const string &archivoEnergia, const string &nombrePNG) {
    ofstream gnu("scripts/grafica_energia.gnu");
    if (gnu.is_open()) {
        gnu << "set terminal pngcairo size 800,400\n";
        gnu << "set output 'results/" << nombrePNG << "'\n";
        gnu << "set xlabel 'Tiempo (s)'\n";
        gnu << "set ylabel 'Energía cinética total (J)'\n";
        gnu << "set grid\n";
        gnu << "plot '" << archivoEnergia << "' with lines lw 2 title 'Energía cinética total'\n";
        gnu.close();
        
        ejecutarComando("gnuplot scripts/grafica_energia.gnu");
    }
}

void GenerarGraficaPresion(const string &archivoPresion, const string &nombrePNG) {
    ofstream gnu("scripts/grafica_presion.gnu");
    if (gnu.is_open()) {
        gnu << "set terminal pngcairo size 800,400\n";
        gnu << "set output 'results/" << nombrePNG << "'\n";
        gnu << "set xlabel 'Tiempo (s)'\n";
        gnu << "set ylabel 'Presión (N/m)'\n";
        gnu << "set grid\n";

        // Calcular promedio de presión por colisiones
        gnu << "stats '" << archivoPresion << "' using 3 nooutput\n";
        gnu << "promedio_col = STATS_mean\n";

        gnu << "set title sprintf('Presión: Teórica vs Medida (Promedio P_{col} = %.3f N/m)', promedio_col)\n";

        gnu << "plot '" << archivoPresion << "' using 1:2 with lines lw 2 title 'P por velocidad (teórica)', \\\n";
        gnu << "     '" << archivoPresion << "' using 1:3 with lines lw 2 title 'P por colisiones (medida)', \\\n";
        gnu << "     promedio_col with lines lw 2 lc 'black' dashtype 2 title sprintf('Promedio = %.3f N/m', promedio_col)\n";

        gnu.close();
        
        ejecutarComando("gnuplot scripts/grafica_presion.gnu");
    }
}

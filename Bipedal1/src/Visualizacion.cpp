#include "bipedal.h"
#include <fstream>
#include <cmath> 
#include <cstdlib>
#include <iostream>
#include <sys/stat.h>
#include <sstream>

using namespace std;

void crear_directorio(const string& path) {
    mkdir(path.c_str(), 0777);
}

void graficar_posiciones(const string& archivo_datos, 
                        const string& archivo_png) {
    crear_directorio("scripts");
    crear_directorio("results");
    
    ofstream script("scripts/posiciones.gnu");
    script << "set terminal pngcairo size 800,600 enhanced font 'Arial,12'\n";
    script << "set output '" << archivo_png << "'\n";
    script << "set xlabel 'Tiempo (s)'\n";
    script << "set ylabel 'Posición (m)'\n";
    script << "set title 'Posiciones de las cabezas'\n";
    script << "set grid\n";
    script << "plot '" << archivo_datos << "' u 1:2 w l lw 2 title 'Cabeza 1', \\\n";
    script << "     '" << archivo_datos << "' u 1:3 w l lw 2 title 'Cabeza 2'\n";
    script.close();
    
    int resultado = system("gnuplot scripts/posiciones.gnu");
    (void)resultado;
    cout << "Gráfica guardada: " << archivo_png << endl;
}

void graficar_velocidades(const string& archivo_datos, 
                         const string& archivo_png) {
    crear_directorio("scripts");
    crear_directorio("results");
    
    ofstream script("scripts/velocidades.gnu");
    script << "set terminal pngcairo size 800,600 enhanced font 'Arial,12'\n";
    script << "set output '" << archivo_png << "'\n";
    script << "set xlabel 'Tiempo (s)'\n";
    script << "set ylabel 'Velocidad (m/s)'\n";
    script << "set title 'Velocidades de las cabezas'\n";
    script << "set grid\n";
    script << "plot '" << archivo_datos << "' u 1:4 w l lw 2 title 'Cabeza 1', \\\n";
    script << "     '" << archivo_datos << "' u 1:5 w l lw 2 title 'Cabeza 2'\n";
    script.close();
    
    int resultado = system("gnuplot scripts/velocidades.gnu");
    (void)resultado;
    cout << "Gráfica guardada: " << archivo_png << endl;
}

void graficar_estados(const string& archivo_datos, 
                     const string& archivo_png) {
    crear_directorio("scripts");
    crear_directorio("results");
    
    ofstream script("scripts/estados.gnu");
    script << "set terminal pngcairo size 800,600 enhanced font 'Arial,12'\n";
    script << "set output '" << archivo_png << "'\n";
    script << "set xlabel 'Tiempo (s)'\n";
    script << "set ylabel 'Estado (1=ON, 0=OFF)'\n";
    script << "set title 'Estados de activación'\n";
    script << "set yrange [-0.1:1.1]\n";
    script << "set grid\n";
    script << "plot '" << archivo_datos << "' u 1:6 w steps lw 2 title 'Cabeza 1', \\\n";
    script << "     '" << archivo_datos << "' u 1:7 w steps lw 2 title 'Cabeza 2'\n";
    script.close();
    
    int resultado = system("gnuplot scripts/estados.gnu");
    (void)resultado;
    cout << "Gráfica guardada: " << archivo_png << endl;
}

void graficar_distancia(const string& archivo_datos, 
                       const string& archivo_png) {
    crear_directorio("scripts");
    crear_directorio("results");
    
    ofstream script("scripts/distancia.gnu");
    script << "set terminal pngcairo size 800,600 enhanced font 'Arial,12'\n";
    script << "set output '" << archivo_png << "'\n";
    script << "set xlabel 'Tiempo (s)'\n";
    script << "set ylabel 'Distancia (m)'\n";
    script << "set title 'Distancia entre cabezas'\n";
    script << "set grid\n";
    script << "plot '" << archivo_datos << "' u 1:(abs($3-$2)) w l lw 2 title '|x2 - x1|'\n";
    script.close();
    
    int resultado = system("gnuplot scripts/distancia.gnu");
    (void)resultado;
    cout << "Gráfica guardada: " << archivo_png << endl;
}

void animar_motor(const string& archivo_datos, 
                 const string& archivo_gif) {
    crear_directorio("scripts");
    crear_directorio("results");
    
    // Obtener rango de posiciones
    ifstream datos(archivo_datos);
    string linea;
    double x1, x2, min_x = 1e100, max_x = -1e100;
    
    getline(datos, linea);
    
    while (getline(datos, linea)) {
        stringstream ss(linea);
        double t;
        if (ss >> t >> x1 >> x2) {
            if (x1 < min_x) min_x = x1;
            if (x1 > max_x) max_x = x1;
            if (x2 < min_x) min_x = x2;
            if (x2 > max_x) max_x = x2;
        }
    }
    datos.close();
    
    // Añadir margen
    double margen = (max_x - min_x) * 0.1;
    if (margen < 1e-12) margen = 1e-9; // margen mínimo
    min_x -= margen;
    max_x += margen;
    
    // Crear script de animación
    ofstream script("scripts/animacion.gnu");
    script << "set terminal gif animate delay 10 size 800,400\n";
    script << "set output '" << archivo_gif << "'\n";
    script << "set xlabel 'Posición (m)'\n";
    script << "set ylabel 'Partícula'\n";
    script << "set title 'Motor Bipedal - Animación 1D'\n";
    script << "set yrange [0:3]\n";
    script << "set xrange [" << min_x << ":" << max_x << "]\n";
    script << "set grid\n";
    script << "set key off\n\n";
    
    // Dibujar líneas de referencia para el potencial
    double l0 = 8.3e-9;
    int num_periodos = int((max_x - min_x) / (2*l0)) + 2;
    double x_inicio = floor(min_x/(2*l0)) * (2*l0);
    for (int i = 0; i < num_periodos; i++) {
        double x = x_inicio + i * 2*l0;
        script << "set arrow from " << x << ",0 to " << x << ",3 nohead lc 'gray' dt 2\n";
    }
    
    // Animación (solo primeros 1000 frames)
    script << "do for [i=0:999] {\n";
    script << "  plot '" << archivo_datos << "' every ::i::i u 2:(1):(0.2e-9) w circles lc 'red' fill solid, \\\n";
    script << "       '" << archivo_datos << "' every ::i::i u 3:(2):(0.2e-9) w circles lc 'blue' fill solid\n";
    script << "}\n";
    script.close();
    
    int resultado = system("gnuplot scripts/animacion.gnu");
    (void)resultado;
    cout << "Animación guardada: " << archivo_gif << endl;
}

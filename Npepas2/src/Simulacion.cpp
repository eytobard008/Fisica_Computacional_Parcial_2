/**
 * @file Simulacion.cpp
 * @brief Implementación de la lógica de simulación y escritura de datos
 * @author Esteban Tobar
 * @date 29 de octubre de 2025
 */

#include "../include/Npepas.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include <ctime>

using namespace std;

// Variables globales para conteo de cálculos
extern long long choquesEntreParticulas;
extern long long calculosDistancias;
extern long long calculosMalla;
extern long long calculosVecindad;

// ==========================
// Funciones de entrada y validación
// ==========================

void SolicitarDatos(double& R, int& N, double& tmax, double& velocidadMedia, double& dispersionVelocidad, bool& usarMallaOpt) {
    double W = 20.0, H = 10.0;
    
    cout << "Tamaño de la caja: " << W << "m x " << H << "m\n";
    cout << "Radio de las bolas (m): ";
    cin >> R;
    
    MostrarLimites(R, W, H);
    
    cout << "Número de partículas N: ";
    cin >> N;
    cout << "Magnitud media de velocidad (m/s): ";
    cin >> velocidadMedia;
    cout << "Desviación estándar de la magnitud (m/s): ";
    cin >> dispersionVelocidad;
    cout << "Tiempo máximo tmax (s): ";
    cin >> tmax;
    cout << "¿Usar optimización con malla? (0=No, 1=Sí): ";
    cin >> usarMallaOpt;
}

bool VerificarDatos(double R, int N, double velocidadMedia, double dispersionVelocidad) {
    double W = 20.0, H = 10.0;
    
    if (velocidadMedia < 0) {
        cout << "Error: La velocidad media no puede ser negativa.\n";
        return false;
    }
    if (dispersionVelocidad < 0) {
        cout << "Error: La desviación estándar no puede ser negativa.\n";
        return false;
    }
    if (R <= 0 || N <= 0 || velocidadMedia <= 0) {
        cout << "Error: Parámetros deben ser positivos.\n";
        return false;
    }
    
    return VerificarAforo(R, N, W, H, false);
}

// ==========================
// Funciones de cálculo de aforo y límites
// ==========================

int CalcularAforoMaximo(double R, double W, double H) {
    // Calcula el número máximo de partículas que caben en disposición hexagonal
    // Esto evita solapamientos y maximiza el empaquetamiento
    if (2 * R > W || 2 * R > H) {
        return 0;
    }

    int nCols = static_cast<int>(W / (2 * R));
    int nRows = static_cast<int>(H / (2 * R));

    int total = 0;
    for (int i = 0; i < nRows; i++) {
        if (i % 2 == 0) {
            total += nCols;
        } else {
            total += (nCols - 1);
        }
    }

    return total;
}

bool VerificarAforo(double R, int N, double W, double H, bool usarMalla) {
    int aforoMaximo = CalcularAforoMaximo(R, W, H);

    if (aforoMaximo == 0) {
        cerr << "Error: El radio de las partículas es demasiado grande para la caja.\n";
        cerr << "El diámetro (2*R = " << 2 * R << " m) debe ser menor que ambos lados de la caja (W=" << W << " m, H=" << H << " m)\n";
        return false;
    }

    if (usarMalla && N > aforoMaximo) {
        cerr << "Error: No caben " << N << " partículas en malla regular. Máximo: " << aforoMaximo << "\n";
        return false;
    }

    if (!usarMalla && N > aforoMaximo * 0.75) {
        cout << "Advertencia: El número de partículas (" << N << ") supera el aforo recomendado\n";
    }

    return true;
}

void MostrarLimites(double R, double W, double H) {
    int aforoMaximo = CalcularAforoMaximo(R, W, H);
    int nMaxRecomendado = static_cast<int>(aforoMaximo * 0.75);
    
    cout << "Aforo máximo: " << aforoMaximo << " partículas\n";
    cout << "N máximo recomendado: " << nMaxRecomendado << " partículas\n";
}

// ==========================
// Funciones de inicialización del sistema
// ==========================

bool InicializarSistema(Bola *bolas, int N, double R, double W, double H, 
                       mt19937 &gen, double velocidadMedia, 
                       double dispersionVelocidad, double masa_particula, 
                       bool usarMalla) {
    if (!VerificarAforo(R, N, W, H, usarMalla)) {
        return false;
    }

    // Inicializar todas las partículas
    for (int i = 0; i < N; i++) {
        bolas[i].Inicie(0, 0, 0, 0, R, masa_particula);
    }

    if (usarMalla) {
        // Disposición en malla regular
        int nCols = static_cast<int>(W / (2 * R));
        int nRows = static_cast<int>(H / (2 * R));
        int contador = 0;

        for (int i = 0; i < nRows && contador < N; i++) {
            for (int j = 0; j < nCols && contador < N; j++) {
                double x = R + 2 * R * j;
                double y = R + 2 * R * i;
                bolas[contador].Inicie(x, y, 0, 0, R, masa_particula);
                contador++;
            }
        }
    } else {
        // Disposición aleatoria
        uniform_real_distribution<double> distX(R, W - R);
        uniform_real_distribution<double> distY(R, H - R);

        for (int i = 0; i < N; i++) {
            bool valido = false;
            int intentos = 0;
            const int max_intentos = 1000000;

            while (!valido && intentos < max_intentos) {
                double x = distX(gen);
                double y = distY(gen);
                bolas[i].Inicie(x, y, 0, 0, R, masa_particula);
                valido = true;
                
                for (int j = 0; j < i; j++) {
                    double dx = bolas[i].GetX() - bolas[j].GetX();
                    double dy = bolas[i].GetY() - bolas[j].GetY();
                    double distancia = sqrt(dx*dx + dy*dy);
                    if (distancia < 2.5 * R) {
                        valido = false;
                        break;
                    }
                }
                intentos++;
            }

            if (!valido) {
                cerr << "Error: No se pudo colocar la partícula " << i 
                     << " sin solapamiento después de " << max_intentos << " intentos.\n";
                return false;
            }
        }
    }

    // Asignar velocidades aleatorias
    normal_distribution<double> distMagnitud(velocidadMedia, dispersionVelocidad);
    uniform_real_distribution<double> distAngulo(0.0, 2.0 * M_PI);

    for (int i = 0; i < N; i++) {
        double magnitud = distMagnitud(gen);
        if (magnitud < 0) magnitud = -magnitud;

        double angulo = distAngulo(gen);
        double vx = magnitud * cos(angulo);
        double vy = magnitud * sin(angulo);

        bolas[i].Inicie(bolas[i].GetX(), bolas[i].GetY(), vx, vy, R, masa_particula);
    }

    return true;
}

// ==========================
// Funciones de cálculo de presión
// ==========================

double CalcularPresionPorVelocidad(Bola *bolas, int N, double W, double H, double masa) {
    // Calcula la presión usando la teoría cinética de gases:
    // P = (N * m * <v²>) / (2 * A) para 2D
    double sum_vx2 = 0.0, sum_vy2 = 0.0;

    for (int i = 0; i < N; i++) {
        double vx = bolas[i].GetVx();
        double vy = bolas[i].GetVy();
        sum_vx2 += vx * vx;
        sum_vy2 += vy * vy;
    }

    double avg_vx2 = sum_vx2 / N;
    double avg_vy2 = sum_vy2 / N;
    double A = W * H;

    double presion = (N * masa * (avg_vx2 + avg_vy2)) / (2.0 * A);
    return presion;
}

void CalcularPresionPorColisiones(Bola *bolas, int N, double W, double H, double dt, double masa, double& suma_vx_colisiones, double& suma_vy_colisiones, double& presion_colisiones) {
    // Calcula la presión basada en el momento transferido a las paredes
    // P = (Δp / Δt) / L para cada pared
    suma_vx_colisiones = 0.0;
    suma_vy_colisiones = 0.0;

    for (int i = 0; i < N; i++) {
        bolas[i].RebotarEnCaja(W, H, suma_vx_colisiones, suma_vy_colisiones);
    }

    double P_x = (suma_vx_colisiones * masa) / (H * dt);
    double P_y = (suma_vy_colisiones * masa) / (W * dt);
    presion_colisiones = (P_x + P_y) / 2.0;
}

// ==========================
// Función principal de simulación
// ==========================

void SimularSistema(Bola *bolas, int N, double W, double H, double dt, double tmax, 
                   const string &archivoSalida, const string &archivoEnergia, 
                   const string &archivoPresion, double masa, double velocidadMedia, 
                   double dispersionVelocidad, bool usarMallaOpt, bool contarCalculos,
                   double epsilon, double sigma) {
    // Abrir archivos de salida
    ofstream datos(archivoSalida);
    ofstream energiaFile(archivoEnergia);
    ofstream presionFile(archivoPresion);

    if (!datos.is_open() || !energiaFile.is_open() || !presionFile.is_open()) {
        cerr << "Error: no se pudieron abrir los archivos de salida\n";
        return;
    }

    // Escribir encabezados
    datos << "#" << setw(12) << "t";
    for (int i = 0; i < N; i++) {
        datos << setw(12) << "x" + to_string(i + 1) 
              << setw(12) << "y" + to_string(i + 1) 
              << setw(12) << "vx" + to_string(i + 1) 
              << setw(12) << "vy" + to_string(i + 1);
    }
    datos << "\n";

    energiaFile << "# t E_cinetica E_potencial E_total\n";
    presionFile << "# t P_cinetica P_virial P_total\n";

    int pasos = static_cast<int>(tmax / dt);
    datos << fixed << setprecision(6);
    energiaFile << fixed << setprecision(6);
    presionFile << fixed << setprecision(6);

    cout << "Simulando " << pasos << " pasos de tiempo (dt = " << dt << " s)...\n";

    // Calcular tamaño óptimo de celda para la malla espacial (PARA OPTIMIZACIÓN)
    double tamanoCeldaVel = (velocidadMedia + dispersionVelocidad) * dt * 3.0;
    double tamanoCeldaRadio = 2.5 * sigma;  // Usar sigma para cutoff
    double tamanoCelda = max(tamanoCeldaVel, tamanoCeldaRadio);
    
    // Variables para la malla espacial optimizada
    vector<Celda> malla;
    int celdasX, celdasY;
    vector<int> offsetsVecinos;
    vector<int> celdasUsadas;
    
    if (usarMallaOpt) {
        InicializarMalla(malla, celdasX, celdasY, offsetsVecinos, W, H, tamanoCelda);
        cout << "Optimización: usando " << celdasX << "×" << celdasY << " celdas (tamaño " << tamanoCelda << " m)\n";
    }

    // Parámetros Lennard-Jones
    double r_cutoff = 2.5 * sigma;
    double r_cutoff2 = r_cutoff * r_cutoff;

    // Inicializar fuerzas
    for (int i = 0; i < N; i++) {
        bolas[i].ReiniciarFuerza();
    }

    // Calcular fuerzas iniciales
    double energia_potencial = 0.0;
    double virial = 0.0;
    
    // Fuerzas entre partículas
    if (usarMallaOpt && contarCalculos) {
        calculosMalla = 0;
        calculosVecindad = 0;
        calculosDistancias = 0;
    }
    
    if (usarMallaOpt) {
        // Usar malla para calcular fuerzas
        ActualizarMalla(bolas, N, malla, celdasX, celdasY, celdasUsadas, tamanoCelda, W, H, contarCalculos);
        
        // Calcular fuerzas usando malla
        for (int indiceCelda : celdasUsadas) {
            Celda& celdaActual = malla[indiceCelda];
            int numParticulas = celdaActual.particulas.size();
            
            // Fuerzas dentro de la misma celda
            for (int i = 0; i < numParticulas; i++) {
                for (int j = i + 1; j < numParticulas; j++) {
                    int idx1 = celdaActual.particulas[i];
                    int idx2 = celdaActual.particulas[j];
                    
                    double dx = bolas[idx1].GetX() - bolas[idx2].GetX();
                    double dy = bolas[idx1].GetY() - bolas[idx2].GetY();
                    double r2 = dx*dx + dy*dy;
                    
                    if (contarCalculos) calculosDistancias++;
                    
                    if (r2 < r_cutoff2 && r2 > 1e-12) {
                        double r = sqrt(r2);
                        double sr6 = pow(sigma/r, 6);
                        double sr12 = sr6 * sr6;
                        
                        // Energía potencial
                        energia_potencial += 4.0 * epsilon * (sr12 - sr6);
                        
                        // Fuerza
                        double factor = 24.0 * epsilon * (2.0 * sr12 - sr6) / r2;
                        double fx = factor * dx;
                        double fy = factor * dy;
                        
                        bolas[idx1].SumarFuerza(fx, fy);
                        bolas[idx2].SumarFuerza(-fx, -fy);
                        
                        // Virial
                        virial += fx * dx + fy * dy;
                    }
                }
            }
            
            // Fuerzas con celdas vecinas
            for (int offset : offsetsVecinos) {
                int indiceVecino = indiceCelda + offset;
                if (offset == 0 || indiceVecino < 0 || indiceVecino >= celdasX*celdasY) continue;
                
                Celda& celdaVecina = malla[indiceVecino];
                int numParticulasVecina = celdaVecina.particulas.size();
                
                for (int i = 0; i < numParticulas; i++) {
                    for (int j = 0; j < numParticulasVecina; j++) {
                        int idx1 = celdaActual.particulas[i];
                        int idx2 = celdaVecina.particulas[j];
                        
                        if (idx1 >= idx2) continue;
                        
                        double dx = bolas[idx1].GetX() - bolas[idx2].GetX();
                        double dy = bolas[idx1].GetY() - bolas[idx2].GetY();
                        double r2 = dx*dx + dy*dy;
                        
                        if (contarCalculos) calculosDistancias++;
                        
                        if (r2 < r_cutoff2 && r2 > 1e-12) {
                            double r = sqrt(r2);
                            double sr6 = pow(sigma/r, 6);
                            double sr12 = sr6 * sr6;
                            
                            energia_potencial += 4.0 * epsilon * (sr12 - sr6);
                            double factor = 24.0 * epsilon * (2.0 * sr12 - sr6) / r2;
                            double fx = factor * dx;
                            double fy = factor * dy;
                            
                            bolas[idx1].SumarFuerza(fx, fy);
                            bolas[idx2].SumarFuerza(-fx, -fy);
                            virial += fx * dx + fy * dy;
                        }
                    }
                }
            }
        }
    } else {
        // Método directo O(N²)
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                double dx = bolas[i].GetX() - bolas[j].GetX();
                double dy = bolas[i].GetY() - bolas[j].GetY();
                double r2 = dx*dx + dy*dy;
                
                if (contarCalculos) calculosDistancias++;
                
                if (r2 < r_cutoff2 && r2 > 1e-12) {
                    double r = sqrt(r2);
                    double sr6 = pow(sigma/r, 6);
                    double sr12 = sr6 * sr6;
                    
                    energia_potencial += 4.0 * epsilon * (sr12 - sr6);
                    double factor = 24.0 * epsilon * (2.0 * sr12 - sr6) / r2;
                    double fx = factor * dx;
                    double fy = factor * dy;
                    
                    bolas[i].SumarFuerza(fx, fy);
                    bolas[j].SumarFuerza(-fx, -fy);
                    virial += fx * dx + fy * dy;
                }
            }
        }
    }
    
    // Fuerzas de paredes
    for (int i = 0; i < N; i++) {
        bolas[i].AgregarFuerzaParedesLJ(W, H, epsilon, sigma);
    }
    
    // Guardar fuerzas iniciales
    for (int i = 0; i < N; i++) {
        bolas[i].GuardarFuerzaAnterior();
    }

    auto start = chrono::high_resolution_clock::now();

    // Bucle principal de simulación con Velocity Verlet
    for (int step = 0; step < pasos; step++) {
        double t = step * dt;
        
        // --- PASO 1 de Velocity Verlet: actualizar posiciones ---
        for (int i = 0; i < N; i++) {
            bolas[i].MoverStep1(dt);
        }
        
        // --- Calcular nuevas fuerzas ---
        // Reiniciar fuerzas y energías
        for (int i = 0; i < N; i++) {
            bolas[i].ReiniciarFuerza();
        }
        energia_potencial = 0.0;
        virial = 0.0;
        
        if (usarMallaOpt) {
            // Actualizar malla
            ActualizarMalla(bolas, N, malla, celdasX, celdasY, celdasUsadas, tamanoCelda, W, H, contarCalculos);
            
            // Calcular fuerzas usando malla
            for (int indiceCelda : celdasUsadas) {
                Celda& celdaActual = malla[indiceCelda];
                int numParticulas = celdaActual.particulas.size();
                
                for (int i = 0; i < numParticulas; i++) {
                    for (int j = i + 1; j < numParticulas; j++) {
                        int idx1 = celdaActual.particulas[i];
                        int idx2 = celdaActual.particulas[j];
                        
                        double dx = bolas[idx1].GetX() - bolas[idx2].GetX();
                        double dy = bolas[idx1].GetY() - bolas[idx2].GetY();
                        double r2 = dx*dx + dy*dy;
                        
                        if (contarCalculos) calculosDistancias++;
                        
                        if (r2 < r_cutoff2 && r2 > 1e-12) {
                            double r = sqrt(r2);
                            double sr6 = pow(sigma/r, 6);
                            double sr12 = sr6 * sr6;
                            
                            energia_potencial += 4.0 * epsilon * (sr12 - sr6);
                            double factor = 24.0 * epsilon * (2.0 * sr12 - sr6) / r2;
                            double fx = factor * dx;
                            double fy = factor * dy;
                            
                            bolas[idx1].SumarFuerza(fx, fy);
                            bolas[idx2].SumarFuerza(-fx, -fy);
                            virial += fx * dx + fy * dy;
                        }
                    }
                }
                
                for (int offset : offsetsVecinos) {
                    int indiceVecino = indiceCelda + offset;
                    if (offset == 0 || indiceVecino < 0 || indiceVecino >= celdasX*celdasY) continue;
                    
                    Celda& celdaVecina = malla[indiceVecino];
                    int numParticulasVecina = celdaVecina.particulas.size();
                    
                    for (int i = 0; i < numParticulas; i++) {
                        for (int j = 0; j < numParticulasVecina; j++) {
                            int idx1 = celdaActual.particulas[i];
                            int idx2 = celdaVecina.particulas[j];
                            
                            if (idx1 >= idx2) continue;
                            
                            double dx = bolas[idx1].GetX() - bolas[idx2].GetX();
                            double dy = bolas[idx1].GetY() - bolas[idx2].GetY();
                            double r2 = dx*dx + dy*dy;
                            
                            if (contarCalculos) calculosDistancias++;
                            
                            if (r2 < r_cutoff2 && r2 > 1e-12) {
                                double r = sqrt(r2);
                                double sr6 = pow(sigma/r, 6);
                                double sr12 = sr6 * sr6;
                                
                                energia_potencial += 4.0 * epsilon * (sr12 - sr6);
                                double factor = 24.0 * epsilon * (2.0 * sr12 - sr6) / r2;
                                double fx = factor * dx;
                                double fy = factor * dy;
                                
                                bolas[idx1].SumarFuerza(fx, fy);
                                bolas[idx2].SumarFuerza(-fx, -fy);
                                virial += fx * dx + fy * dy;
                            }
                        }
                    }
                }
            }
        } else {
            // Método directo
            for (int i = 0; i < N; i++) {
                for (int j = i + 1; j < N; j++) {
                    double dx = bolas[i].GetX() - bolas[j].GetX();
                    double dy = bolas[i].GetY() - bolas[j].GetY();
                    double r2 = dx*dx + dy*dy;
                    
                    if (contarCalculos) calculosDistancias++;
                    
                    if (r2 < r_cutoff2 && r2 > 1e-12) {
                        double r = sqrt(r2);
                        double sr6 = pow(sigma/r, 6);
                        double sr12 = sr6 * sr6;
                        
                        energia_potencial += 4.0 * epsilon * (sr12 - sr6);
                        double factor = 24.0 * epsilon * (2.0 * sr12 - sr6) / r2;
                        double fx = factor * dx;
                        double fy = factor * dy;
                        
                        bolas[i].SumarFuerza(fx, fy);
                        bolas[j].SumarFuerza(-fx, -fy);
                        virial += fx * dx + fy * dy;
                    }
                }
            }
        }
        
        // Fuerzas de paredes
        for (int i = 0; i < N; i++) {
            bolas[i].AgregarFuerzaParedesLJ(W, H, epsilon, sigma);
        }
        
        // --- PASO 2 de Velocity Verlet: actualizar velocidades ---
        for (int i = 0; i < N; i++) {
            bolas[i].MoverStep2(dt);
        }
        
        // --- Calcular energías y presiones ---
        double energia_cinetica = 0.0;
        double sum_vx2 = 0.0, sum_vy2 = 0.0;
        
        for (int i = 0; i < N; i++) {
            double vx = bolas[i].GetVx();
            double vy = bolas[i].GetVy();
            energia_cinetica += 0.5 * masa * (vx*vx + vy*vy);
            sum_vx2 += vx * vx;
            sum_vy2 += vy * vy;
        }
        
        double energia_total = energia_cinetica + energia_potencial;
        
        // Presión: término cinético + término virial
        double area = W * H;
        double presion_cinetica = (N * masa * ((sum_vx2 + sum_vy2) / N)) / (2.0 * area);
        double presion_virial = virial / (2.0 * area);
        double presion_total = presion_cinetica + presion_virial;
        
        // --- Escribir datos del paso actual ---
        datos << setw(13) << t;
        for (int i = 0; i < N; i++) {
            datos << setw(12) << bolas[i].GetX() 
                  << setw(12) << bolas[i].GetY() 
                  << setw(12) << bolas[i].GetVx() 
                  << setw(12) << bolas[i].GetVy();
        }
        datos << "\n";
        
        energiaFile << t << " " << energia_cinetica << " " 
                   << energia_potencial << " " << energia_total << "\n";
        
        presionFile << t << " " << presion_cinetica << " " 
                   << presion_virial << " " << presion_total << "\n";
    }

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Tiempo de simulación: " << duration.count()/1000.0 << " segundos\n";

    datos.close();
    energiaFile.close();
    presionFile.close();
    
    cout << "Datos de simulación guardados en results/\n";
    
    // Mostrar estadísticas de cálculos si se solicitó
    if (contarCalculos) {
        long long totalCalculos = calculosDistancias + calculosMalla + calculosVecindad;
        cout << "\n=== CONTEO DE CÁLCULOS ===\n";
        cout << "Cálculos de distancia: " << calculosDistancias << endl;
        if (usarMallaOpt) {
            cout << "Operaciones de malla: " << calculosMalla << endl;
            cout << "Operaciones de vecindad: " << calculosVecindad << endl;
        }
        cout << "Total de cálculos: " << totalCalculos << endl;
    }
}

// ==========================
// Funciones de estadísticas
// ==========================

long long ContarChoquesPared(Bola *bolas, int N) {
    long long total = 0;
    for (int i = 0; i < N; i++) {
        total += bolas[i].GetChoquesPared();
    }
    return total;
}

void MostrarEstadisticas(Bola *bolas, int N, double dt, double tmax) {
    long long totalChoquesPared = ContarChoquesPared(bolas, N);
    double choquesPorSegundoPared = totalChoquesPared / tmax;
    double choquesPorSegundoParticulas = choquesEntreParticulas / tmax;

    cout << "\n=== ESTADÍSTICAS ===\n";
    cout << "Choques con paredes: " << totalChoquesPared << " (" << choquesPorSegundoPared << "/s)\n";
    cout << "Choques entre partículas: " << choquesEntreParticulas << " (" << choquesPorSegundoParticulas << "/s)\n";
    cout << "Total de choques: " << totalChoquesPared + choquesEntreParticulas << endl;
}

/**
 * @file Core.cpp
 * @brief Implementación de la física fundamental y métodos de la clase Bola
 * @author Esteban Tobar
 * @date 29 de octubre de 2025
 */

#include "../include/Npepas.h"
#include <cmath>
#include <vector>
#include <ctime>

using namespace std;

// Variables globales para conteo de cálculos
long long calculosDistancias = 0;
long long calculosMalla = 0;
long long calculosVecindad = 0;

// ==========================
// Implementación de la clase Bola
// ==========================

// ============ NUEVO: Métodos para Velocity Verlet ============

void Bola::Inicie(double x0, double y0, double vx0, double vy0, double R0, double m) {
    x = x0;
    y = y0;
    vx = vx0;
    vy = vy0;
    R = R0;
    masa = m;
    fx = fy = 0.0;
    fx_prev = fy_prev = 0.0;
    choquesPared = 0;
}

void Bola::MoverStep1(double dt) {
    // Paso 1: Actualizar posición usando fuerzas anteriores
    x += vx * dt + 0.5 * (fx_prev / masa) * dt * dt;
    y += vy * dt + 0.5 * (fy_prev / masa) * dt * dt;
}

void Bola::MoverStep2(double dt) {
    // Paso 2: Actualizar velocidad usando promedio de fuerzas
    vx += 0.5 * (fx_prev + fx) * dt / masa;
    vy += 0.5 * (fy_prev + fy) * dt / masa;
    
    // Guardar fuerzas actuales para el siguiente paso
    fx_prev = fx;
    fy_prev = fy;
}

void Bola::ReiniciarFuerza() {
    fx = 0.0;
    fy = 0.0;
}

void Bola::GuardarFuerzaAnterior() {
    fx_prev = fx;
    fy_prev = fy;
}

void Bola::SumarFuerza(double dfx, double dfy) {
    fx += dfx;
    fy += dfy;
}

// ============ NUEVO: Fuerza de Lennard-Jones ============

void Bola::CalcularFuerzaCon(Bola &otra, double epsilon, double sigma) {
    double dx = otra.x - x;
    double dy = otra.y - y;
    double r2 = dx*dx + dy*dy;
    
    // Evitar auto-interacción y distancia cero
    if (r2 < 1e-12) return;
    
    double r = sqrt(r2);
    
    // Calcular (sigma/r)^6 y (sigma/r)^12
    double sr6 = pow(sigma/r, 6);
    double sr12 = sr6 * sr6;
    
    // Fuerza de Lennard-Jones: F = 24*epsilon*(2*(sigma^12/r^13) - (sigma^6/r^7)) * (r_vector/r)
    // Simplificado: factor = 24*epsilon*(2*sr12 - sr6)/r2
    double factor = 24.0 * epsilon * (2.0 * sr12 - sr6) / r2;
    
    double fx_pair = factor * dx;
    double fy_pair = factor * dy;
    
    // Tercera ley de Newton
    SumarFuerza(fx_pair, fy_pair);
    otra.SumarFuerza(-fx_pair, -fy_pair);
}

// ============ NUEVO: Paredes con Lennard-Jones ============

void Bola::RebotarEnCajaLJ(double W, double H, double epsilon, double sigma, double& suma_vx_colisiones, double& suma_vy_colisiones) {
    // Pared izquierda (x = 0)
    if (x < 2.0 * sigma) {
        double dx = x;
        if (dx < 0.5 * sigma) dx = 0.5 * sigma; // Evitar división por cero
        double sr6 = pow(sigma/dx, 6);
        double sr12 = sr6 * sr6;
        double fx_wall = 24.0 * epsilon * (2.0 * sr12 - sr6) / (dx*dx) * dx;
        SumarFuerza(fx_wall, 0.0);
    }
    
    // Pared derecha (x = W)
    if (x > W - 2.0 * sigma) {
        double dx = W - x;
        if (dx < 0.5 * sigma) dx = 0.5 * sigma;
        double sr6 = pow(sigma/dx, 6);
        double sr12 = sr6 * sr6;
        double fx_wall = -24.0 * epsilon * (2.0 * sr12 - sr6) / (dx*dx) * dx;
        SumarFuerza(fx_wall, 0.0);
    }
    
    // Pared inferior (y = 0)
    if (y < 2.0 * sigma) {
        double dy = y;
        if (dy < 0.5 * sigma) dy = 0.5 * sigma;
        double sr6 = pow(sigma/dy, 6);
        double sr12 = sr6 * sr6;
        double fy_wall = 24.0 * epsilon * (2.0 * sr12 - sr6) / (dy*dy) * dy;
        SumarFuerza(0.0, fy_wall);
    }
    
    // Pared superior (y = H)
    if (y > H - 2.0 * sigma) {
        double dy = H - y;
        if (dy < 0.5 * sigma) dy = 0.5 * sigma;
        double sr6 = pow(sigma/dy, 6);
        double sr12 = sr6 * sr6;
        double fy_wall = -24.0 * epsilon * (2.0 * sr12 - sr6) / (dy*dy) * dy;
        SumarFuerza(0.0, fy_wall);
    }
}


// ==========================
// Implementación de malla espacial optimizada
// ==========================

void InicializarMalla(vector<Celda>& malla, int& celdasX, int& celdasY, vector<int>& offsetsVecinos, double W, double H, double tamanoCelda) {
    // Calcular número de celdas en cada dirección
    celdasX = static_cast<int>(ceil(W / tamanoCelda));
    celdasY = static_cast<int>(ceil(H / tamanoCelda));
    int totalCeldas = celdasX * celdasY;
    
    malla.resize(totalCeldas);
    
    // Precalcular offsets para los 9 vecinos (incluyendo celda central)
    // Esto evita calcular repetidamente los índices de celdas vecinas
    offsetsVecinos.clear();
    for (int dj = -1; dj <= 1; dj++) {
        for (int di = -1; di <= 1; di++) {
            offsetsVecinos.push_back(dj * celdasX + di);
        }
    }
}

void ActualizarMalla(Bola* bolas, int N, vector<Celda>& malla, int celdasX, int celdasY, vector<int>& celdasUsadas, double tamanoCelda, double W, double H, bool contarCalculos) {
    // Limpiar solo las celdas que fueron usadas en el frame anterior
    // Esto es más eficiente que limpiar todas las celdas
    for (int indice : celdasUsadas) {
        malla[indice].particulas.clear();
    }
    celdasUsadas.clear();

    // Asignar cada partícula a su celda correspondiente
    for (int i = 0; i < N; i++) {
        double x = bolas[i].GetX();
        double y = bolas[i].GetY();

        // Calcular índices de celda
        int celdaX = static_cast<int>(x / tamanoCelda);
        int celdaY = static_cast<int>(y / tamanoCelda);

        if (contarCalculos) {
            calculosMalla += 2; // Por los cálculos de celdaX y celdaY
        }

        // Verificar que la partícula esté dentro de los límites de la malla
        if (celdaX >= 0 && celdaX < celdasX && celdaY >= 0 && celdaY < celdasY) {
            int indice = celdaY * celdasX + celdaX;
            malla[indice].particulas.push_back(i);
            celdasUsadas.push_back(indice);
            
            if (contarCalculos) {
                calculosMalla += 2; // Por el cálculo del índice y seguimiento de celdas usadas
            }
        }
    }
}

void DetectarColisionesConMalla(Bola* bolas, int N, vector<Celda>& malla, int celdasX, int celdasY, const vector<int>& offsetsVecinos, const vector<int>& celdasUsadas, bool contarCalculos) {
    int totalCeldas = celdasX * celdasY;
    
    // Procesar solo las celdas que contienen partículas
    for (int indiceCelda : celdasUsadas) {
        Celda& celdaActual = malla[indiceCelda];
        int numParticulas = celdaActual.particulas.size();
        
        if (numParticulas == 0) continue;
        
        // Revisar colisiones entre partículas dentro de la misma celda
        // Se usa i+1 para evitar verificar el mismo par dos veces
        for (int i = 0; i < numParticulas; i++) {
            for (int j = i + 1; j < numParticulas; j++) {
                int idx1 = celdaActual.particulas[i];
                int idx2 = celdaActual.particulas[j];
                
                // Verificación rápida de colisión usando distancia al cuadrado
                double dx = bolas[idx1].GetX() - bolas[idx2].GetX();
                double dy = bolas[idx1].GetY() - bolas[idx2].GetY();
                double dist2 = dx * dx + dy * dy;
                double sumaR = bolas[idx1].GetR() + bolas[idx2].GetR();

                if (contarCalculos) {
                    calculosDistancias++;
                }

                if (dist2 <= sumaR * sumaR) {
                    bolas[idx1].Colisionar(bolas[idx2]);
                }
            }
        }
        
        // Revisar colisiones con partículas en celdas vecinas
        // Esto es necesario porque partículas en celdas adyacentes pueden colisionar
        for (int offset : offsetsVecinos) {
            int indiceVecino = indiceCelda + offset;
            
            // Verificar que la celda vecina esté dentro de los límites
            if (indiceVecino < 0 || indiceVecino >= totalCeldas) {
                if (contarCalculos) {
                    calculosVecindad++; // 1 cálculo por verificación de límites
                }
                continue;
            }
            
            // Verificación de offset cero
            if (contarCalculos) {
                calculosVecindad++; // Por la verificación offset == 0
            }
            
            // No comparar con los mismos
            if (offset == 0) continue;
            
            Celda& celdaVecina = malla[indiceVecino];
            int numParticulasVecina = celdaVecina.particulas.size();
            
            if (numParticulasVecina == 0) continue;
            
            // Comparar todas las partículas de la celda actual con todas las de la vecina
            for (int i = 0; i < numParticulas; i++) {
                for (int j = 0; j < numParticulasVecina; j++) {
                    int idx1 = celdaActual.particulas[i];
                    int idx2 = celdaVecina.particulas[j];
                    
                    // Asegurar que idx1 < idx2 para evitar verificar el mismo par dos veces
                    // cuando se procese la celda vecina
                    if (idx1 >= idx2) continue;
                    
                    // Verificación rápida de colisión
                    double dx = bolas[idx1].GetX() - bolas[idx2].GetX();
                    double dy = bolas[idx1].GetY() - bolas[idx2].GetY();
                    double dist2 = dx * dx + dy * dy;
                    double sumaR = bolas[idx1].GetR() + bolas[idx2].GetR();

                    if (contarCalculos) {
                        calculosDistancias++;
                    }

                    if (dist2 <= sumaR * sumaR) {
                        bolas[idx1].Colisionar(bolas[idx2]);
                    }
                }
            }
        }
        
        // Contar operaciones de gestión de vecindad por celda activa
        if (contarCalculos) {
            calculosVecindad += 9; // Por los 9 offsets procesados por cada celda activa
            calculosVecindad += numParticulas; // Por el procesamiento interno de la celda
        }
    }
}

double CalcularEnergiaPotencial(Bola *bolas, int N, double epsilon, double sigma) {
    double energia = 0.0;
    double r_cutoff = 2.5 * sigma;
    
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double dx = bolas[i].GetX() - bolas[j].GetX();
            double dy = bolas[i].GetY() - bolas[j].GetY();
            double r = sqrt(dx*dx + dy*dy);
            
            if (r < r_cutoff) {
                double sr6 = pow(sigma/r, 6);
                double sr12 = sr6 * sr6;
                energia += 4.0 * epsilon * (sr12 - sr6);
            }
        }
    }
    return energia;
}

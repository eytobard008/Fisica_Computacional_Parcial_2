/**
 * @file Npepas.h
 * @brief Interfaz principal para la simulación de N partículas en una caja 2D
 * @author Esteban Tobar
 * @date 29 de octubre de 2025
 */

#ifndef NPEPAS_H
#define NPEPAS_H

#include <string>
#include <random>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cmath>

/**
 * @brief Representa una partícula esférica en movimiento 2D
 */
class Bola {
private:
    double x, y;      ///< Coordenadas de posición (m)
    double vx, vy;    ///< Componentes de velocidad (m/s)
    double fx, fy;    ///< Fuerza actual (N)
    double fx_prev, fy_prev; ///< Fuerza del paso anterior (N)
    double R;         ///< Radio de la partícula (m)
    double masa;      ///< Masa de la partícula (kg)
    long long choquesPared; ///< Contador de choques con las paredes

public:
    void Inicie(double x0, double y0, double vx0, double vy0, double R0, double m);
    void MoverStep1(double dt);
    void MoverStep2(double dt);
    void ReiniciarFuerza();
    void SumarFuerza(double dfx, double dfy);
    void GuardarFuerzaAnterior();
    void AgregarFuerzaParedesLJ(double W, double H, double epsilon, double sigma);

    // Getters
    double GetX() const { return x; }
    double GetY() const { return y; }
    double GetVx() const { return vx; }
    double GetVy() const { return vy; }
    double GetR() const { return R; }
    double GetMasa() const { return masa; }
    double GetV() const { return std::sqrt(vx*vx + vy*vy); }
    long long GetChoquesPared() const { return choquesPared; }
    void CalcularFuerzaCon(Bola &otra, double epsilon, double sigma);
    void RebotarEnCajaLJ(double W, double H, double epsilon, double sigma, double& suma_vx_colisiones, double& suma_vy_colisiones);
    void CalcularFuerzaCon(Bola &otra, double epsilon, double sigma);
    void AgregarFuerzaParedesLJ(double W, double H, double epsilon, double sigma);
};

/**
 * @brief Estructura para celda en la malla espacial
 */
struct Celda {
    std::vector<int> particulas; ///< Índices de partículas en esta celda
};

// Funciones de configuración y validación
void SolicitarDatos(double& R, int& N, double& tmax, double& velocidadMedia, double& dispersionVelocidad, bool& usarMallaOpt);
bool VerificarDatos(double R, int N, double velocidadMedia, double dispersionVelocidad);
bool VerificarAforo(double R, int N, double W, double H, bool usarMalla);
void MostrarLimites(double R, double W, double H);
int CalcularAforoMaximo(double R, double W, double H);

// Funciones de inicialización del sistema
bool InicializarSistema(Bola *bolas, int N, double R, double W, double H, std::mt19937 &gen, 
                       double velocidadMedia, double dispersionVelocidad, double masa_particula, 
                       bool usarMalla = false);

// Funciones de simulación principal
void SimularSistema(Bola *bolas, int N, double W, double H, double dt, double tmax, 
                   const std::string &archivoSalida, const std::string &archivoEnergia, 
                   const std::string &archivoPresion, double masa, double velocidadMedia, 
                   double dispersionVelocidad, bool usarMalla, bool contarCalculos,
                   double epsilon, double sigma);

// Funciones de malla espacial optimizada
void InicializarMalla(std::vector<Celda>& malla, int& celdasX, int& celdasY, std::vector<int>& offsetsVecinos, double W, double H, double tamanoCelda);
void ActualizarMalla(Bola* bolas, int N, std::vector<Celda>& malla, int celdasX, int celdasY, std::vector<int>& celdasUsadas, double tamanoCelda, double W, double H, bool contarCalculos);

// Funciones de cálculo de presión
double CalcularPresionPorVelocidad(Bola *bolas, int N, double W, double H, double masa);

// Funciones de visualización y análisis
void GenerarAnimacion(const std::string &archivoDat, const std::string &archivoPresion, double W, double H, int N, int pasos, double dt, double R, bool usarMalla, double velocidadMedia, double dispersionVelocidad, const std::string &nombreGif);
void GenerarHistogramaAnimado(const std::string &archivoDat, int N, int pasos, double dt, double velocidadMedia, const std::string &nombreGif, const std::string &nombreGifZoom);
void GenerarGraficaEnergia(const std::string &archivoEnergia, const std::string &nombrePNG);
void GenerarGraficaPresion(const std::string &archivoPresion, const std::string &nombrePNG);

// Funciones de estadísticas y utilidades
void MostrarEstadisticas(Bola *bolas, int N, double dt, double tmax);
long long ContarChoquesPared(Bola *bolas, int N);

// Variables globales para métricas de rendimiento
extern long long calculosDistancias;
extern long long calculosMalla;
extern long long calculosVecindad;

#endif  // NPEPAS_H

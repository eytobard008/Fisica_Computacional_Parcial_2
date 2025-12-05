#include "bipedal.h"
#include <iostream>

using namespace std;

int main() {
    // Constantes físicas
    const double k_B = 1.380649e-23;
    const double T = 300.0;
    const double kBT = k_B * T;
    
    // Parámetros geométricos
    const double l0 = 8.3e-9;
    
    // Parámetros del potencial
    const double V0 = 20.0 * kBT;
    const double xM = 0.5 * l0;
    
    // Parámetros de fricción
    const double eta = 1e-3;
    const double R_cabeza = 2.5e-9;
    const double PI = 3.14159265358979323846;
    const double gamma = 6.0 * PI * eta * R_cabeza;
    
    // Tiempo característico
    const double tau = l0 * l0 * gamma / V0;
    
    // Parámetros de simulación
    const double dt = 0.001 * tau;
    const double t_on = 20.0 * tau;
    const double t_off = 20.0 * tau;
    const double duracion_total = 1000.0 * tau;
    
    // Resorte
    const double k_resorte = 10.0 * V0 / (l0 * l0);
    
    // Masa
    const double masa = 100.0 * 1.66e-27;
    
    cout << "SIMULACIÓN MOTOR BIPEDAL\n";
    cout << "Parámetros:\n";
    cout << "  l0 = " << l0 << " m\n";
    cout << "  V0 = " << V0 << " J\n";
    cout << "  tau = " << tau << " s\n";
    cout << "  dt = " << dt << " s\n";
    cout << "  t_on = t_off = " << t_on << " s\n";
    cout << "  gamma = " << gamma << " kg/s\n";
    cout << "  k_resorte = " << k_resorte << " N/m\n";
    cout << "  masa = " << masa << " kg\n\n";
    
    // Crear motor
    MotorBipedal motor(k_resorte, l0, xM, V0, gamma, kBT, t_on, t_off);
    motor.configurar(0.0, l0, 0.0, 0.0);
    
    // Ejecutar simulación
    simular_motor(duracion_total, dt, motor, "results/trayectorias.dat");
    
    // Generar gráficas
    graficar_posiciones("results/trayectorias.dat", "results/posiciones.png");
    graficar_velocidades("results/trayectorias.dat", "results/velocidades.png");
    graficar_estados("results/trayectorias.dat", "results/estados.png");
    graficar_distancia("results/trayectorias.dat", "results/distancia.png");
    animar_motor("results/trayectorias.dat", "results/animacion.gif");
    
    return 0;
}

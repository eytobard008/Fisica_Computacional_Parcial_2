#ifndef BIPEDAL_H
#define BIPEDAL_H

#include <vector>
#include <string>

/**
 * @brief Partícula en 1D
 */
class Particula {
public:
    double x;   ///< Posición (m)
    double v;   ///< Velocidad (m/s)
    double m;   ///< Masa (kg)
    
    Particula(double x0 = 0.0, double v0 = 0.0, double m0 = 1.0);
};

/**
 * @brief Motor molecular bipedal con flashing ratchet
 */
class MotorBipedal {
private:
    Particula cabeza1, cabeza2;
    double k;           ///< Constante del resorte (N/m)
    double l0_resorte;  ///< Longitud natural del resorte (m)
    double l0_pot;      ///< Mitad del período del potencial (m)
    double xM;          ///< Posición del máximo en [0, 2*l0_pot]
    double V0;          ///< Altura de la barrera (J)
    double gamma;       ///< Coeficiente de fricción (kg/s)
    double kBT;         ///< Energía térmica (J)
    double t;           ///< Tiempo actual (s)
    double t_on;        ///< Tiempo ON (s)
    double t_off;       ///< Tiempo OFF (s)
    bool estado1;       ///< true si cabeza1 siente potencial
    bool estado2;       ///< true si cabeza2 siente potencial
    
    /// Fuerza del potencial asimétrico en posición x
    double fuerza_potencial(double x, bool activo) const;
    
    /// Fuerza del resorte (sobre cabeza1, la de cabeza2 es opuesta)
    double fuerza_resorte() const;
    
    /// Ruido térmico (fuerza aleatoria para un paso dt)
    double ruido(double dt) const;
    
public:
    /**
     * @brief Constructor
     * @param k_resorte Constante del resorte (N/m)
     * @param l0 Longitud natural del resorte = mitad período potencial (m)
     * @param xM_pot Posición del máximo (m), típico 0.5*l0
     * @param V0_pot Altura de barrera (J)
     * @param gamma_val Coeficiente de fricción (kg/s)
     * @param kBT_val Energía térmica (J)
     * @param ton Tiempo ON (s)
     * @param toff Tiempo OFF (s)
     */
    MotorBipedal(double k_resorte, double l0, double xM_pot, double V0_pot,
                 double gamma_val, double kBT_val, double ton, double toff);
    
    /// Configura posiciones y velocidades iniciales
    void configurar(double x1, double x2, double v1 = 0.0, double v2 = 0.0);
    
    /// Avanza la simulación un paso dt (Velocity-Verlet + ruido)
    void avanzar(double dt);
    
    /// Actualiza estados ON/OFF según ciclo químico
    void actualizar_estado(double dt);
    
    // Getters
    double get_x1() const { return cabeza1.x; }
    double get_x2() const { return cabeza2.x; }
    double get_v1() const { return cabeza1.v; }
    double get_v2() const { return cabeza2.v; }
    bool get_estado1() const { return estado1; }
    bool get_estado2() const { return estado2; }
    double get_tiempo() const { return t; }
};

/**
 * @brief Ejecuta simulación completa y guarda datos
 */
void simular_motor(double duracion, double dt, MotorBipedal& motor,
                   const std::string& archivo);

/**
 * @brief Gráfica posiciones vs tiempo
 */
void graficar_posiciones(const std::string& archivo_datos,
                        const std::string& archivo_png);

/**
 * @brief Gráfica velocidades vs tiempo
 */
void graficar_velocidades(const std::string& archivo_datos,
                         const std::string& archivo_png);

/**
 * @brief Gráfica estados ON/OFF
 */
void graficar_estados(const std::string& archivo_datos,
                     const std::string& archivo_png);

/**
 * @brief Gráfica distancia entre cabezas
 */
void graficar_distancia(const std::string& archivo_datos,
                       const std::string& archivo_png);

/**
 * @brief Crea animación GIF del movimiento
 */
void animar_motor(const std::string& archivo_datos,
                 const std::string& archivo_gif);

#endif

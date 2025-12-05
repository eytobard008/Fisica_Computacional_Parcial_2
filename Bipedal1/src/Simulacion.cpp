#include "bipedal.h"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

void simular_motor(double duracion, double dt, MotorBipedal& motor,
                   const std::string& archivo) {
    ofstream file(archivo);
    if (!file.is_open()) {
        cerr << "Error: No se pudo abrir " << archivo << endl;
        return;
    }
    
    file << "# t(s) x1(m) x2(m) v1(m/s) v2(m/s) estado1 estado2\n";
    
    int pasos = static_cast<int>(duracion / dt);
    
    for (int i = 0; i <= pasos; i++) {
        file << motor.get_tiempo() << " "
             << motor.get_x1() << " "
             << motor.get_x2() << " "
             << motor.get_v1() << " "
             << motor.get_v2() << " "
             << motor.get_estado1() << " "
             << motor.get_estado2() << "\n";
        
        if (i < pasos) {
            motor.avanzar(dt);
        }
    }
    
    file.close();
    cout << "Simulación guardada en " << archivo << endl;
}

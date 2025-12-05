double MotorBipedal::fuerza_potencial(double x, bool activo, int cabeza) const {
    if (!activo) return 0.0;
    
    // Aplicar desplazamiento: cabeza2 siente V1(x - l0)
    double x_local = x;
    if (cabeza == 2) {
        x_local -= l0_pot;
    }
    
    // Potencial periódico con período 2*l0_pot
    double x_periodico = fmod(x_local, 2.0*l0_pot);
    if (x_periodico < 0) x_periodico += 2.0*l0_pot;
    
    // Fuerza del diente de sierra
    if (x_periodico < xM) {
        return -V0 / xM;  // Pendiente positiva
    } else {
        return V0 / (2.0*l0_pot - xM);  // Pendiente negativa
    }
}

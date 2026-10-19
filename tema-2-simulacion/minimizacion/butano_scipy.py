#!/usr/bin/env python3
"""Minimización del diedro del butano con scipy.optimize.minimize.

Mismo problema que descenso_1d.py con una molécula real y un minimizador de
biblioteca (BFGS por defecto). Potencial de Ryckaert-Bellemans (1975), el de
butano.py: V(ψ) = Σ Cₙ cosⁿψ con ψ = φ − 180°, en kJ/mol.

Uso:
    python3 butano_scipy.py
"""
import numpy as np
from scipy.optimize import minimize

# Ryckaert-Bellemans para el butano (kJ/mol)
C = [9.2789, 12.1557, -13.1200, -3.0597, 26.2403, -31.4950]

# energía del diedro φ (grados); φ llega en un vector
def energia(phi):
    psi = np.radians(phi[0] - 180)
    return sum(C[n] * np.cos(psi)**n for n in range(6))

# dos puntos de partida a cada lado de la barrera
for inicio in (40, 115, 125, 300):
    # BFGS ajusta φ hasta anular el gradiente
    r = minimize(energia, x0=[inicio])
    phi = r.x[0] % 360
    print(f"{inicio:3d}° -> {phi:5.1f}°  {r.fun:4.2f} kJ/mol")

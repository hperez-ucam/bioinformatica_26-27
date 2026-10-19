#!/usr/bin/env python3
"""Descenso por gradiente en una dimensión: el algoritmo de minimización más simple.

U(x) = x⁴ − 4x² + x es un pozo doble asimétrico: mínimo local en x ≈ 1,35 y
mínimo global en x ≈ −1,47, separados por un máximo en x ≈ 0,13. El algoritmo no
distingue uno de otro: termina en el mínimo de la cuenca de partida.

Uso:
    python3 descenso_1d.py 2.0      # termina en el mínimo local
    python3 descenso_1d.py -0.5     # termina en el mínimo global
"""
import sys

# energía del pozo doble y su derivada
def U(x):  return x**4 - 4*x**2 + x
def dU(x): return 4*x**3 - 8*x + 1

# punto de partida y tamaño del paso
x, paso = float(sys.argv[1]), 0.01
for i in range(10000):
    # fuerza: la derivada cambiada de signo
    F = -dU(x)
    # criterio de parada: fuerza bajo el umbral
    if abs(F) < 1e-6: break
    # avance hacia donde U disminuye
    x = x + paso * F
print(f"x = {x:.4f}  U = {U(x):.4f}  ({i} pasos)")

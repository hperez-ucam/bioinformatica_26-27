#!/usr/bin/env python3
"""linpack_mini.py - un LINPACK de bolsillo: la misma cuenta que usa el Top500.

Resuelve un sistema lineal A x = b de n ecuaciones con n incógnitas y dice
cuántas operaciones de coma flotante por segundo (FLOP/s) ha conseguido.
El Top500 hace exactamente esto (prueba HPL), con n de millones y en miles
de nodos; aquí, con n de unos miles y en un solo servidor.

Uso:
    python3 linpack_mini.py NUCLEOS [N]
    python3 linpack_mini.py 1          # un núcleo, n = 10000
    python3 linpack_mini.py 6 8000     # seis núcleos, n = 8000
"""
import os
import sys

# núcleos que se dejan usar a la biblioteca de álgebra lineal; hay que fijarlo
# ANTES de importar numpy, porque numpy lo lee al cargarse
nucleos = sys.argv[1] if len(sys.argv) > 1 else "1"
os.environ["OPENBLAS_NUM_THREADS"] = nucleos
os.environ["OMP_NUM_THREADS"] = nucleos

import time
import numpy as np

# tamaño del sistema: n ecuaciones con n incógnitas
n = int(sys.argv[2]) if len(sys.argv) > 2 else 10000
# matriz y vector aleatorios, siempre los mismos gracias a la semilla
rng = np.random.default_rng(42)
A = rng.random((n, n))
b = rng.random(n)

# reloj de pared antes y después de resolver
t0 = time.perf_counter()
x = np.linalg.solve(A, b)
t = time.perf_counter() - t0

# operaciones que cuenta HPL para resolver un sistema de n ecuaciones
flops = (2 / 3) * n**3 + 2 * n**2
# comprobación: el residuo debe ser minúsculo si la solución es correcta
residuo = np.linalg.norm(A @ x - b) / np.linalg.norm(b)

print(f"núcleos={nucleos}  n={n}  tiempo={t:.2f} s  "
      f"rendimiento={flops / t / 1e9:.1f} GFLOP/s  residuo={residuo:.1e}")

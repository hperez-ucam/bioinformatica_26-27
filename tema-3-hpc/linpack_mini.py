#!/usr/bin/env python3
"""linpack_mini.py - la prueba del Top500 (HPL) a pequeña escala.

Resuelve un sistema lineal A x = b de n ecuaciones con n incógnitas y calcula
cuántas operaciones de coma flotante por segundo (FLOP/s) ha hecho la máquina.
El Top500 hace esta misma prueba con n de millones y en miles de nodos.

Uso:
    python3 linpack_mini.py HILOS [N]
    python3 linpack_mini.py 1          # un hilo, n = 4000
    python3 linpack_mini.py 6 12000    # seis hilos, n = 12000 (2,2 GB de memoria)

La memoria crece con n al cuadrado: n = 4000 ocupa 0,3 GB; n = 8000, 1 GB;
n = 12000, 2,2 GB. Con n grande, lanzadlo de uno en uno (con la cola).
"""
import os
import sys

# hilos que puede usar la biblioteca de álgebra lineal; hay que fijarlo
# ANTES de importar numpy, porque numpy lo lee al cargarse
hilos = sys.argv[1] if len(sys.argv) > 1 else "1"
os.environ["OPENBLAS_NUM_THREADS"] = hilos
os.environ["OMP_NUM_THREADS"] = hilos

import time
import numpy as np

# tamaño del sistema: n ecuaciones con n incógnitas
n = int(sys.argv[2]) if len(sys.argv) > 2 else 4000

# memoria necesaria: la matriz (8 bytes por número) y la copia que hace la resolución
necesaria_gb = 2 * 8 * n * n / 1e9
# memoria disponible ahora: línea MemAvailable de /proc/meminfo, que viene en kB
with open("/proc/meminfo") as f:
    disponible_gb = int(next(l for l in f if l.startswith("MemAvailable")).split()[1]) / 1e6
# si no cabe con holgura, no empezar: sin memoria, el sistema mata procesos
if necesaria_gb > 0.5 * disponible_gb:
    sys.exit(f"n={n} necesita {necesaria_gb:.1f} GB y hay {disponible_gb:.1f} GB disponibles: "
             "usad un n menor o esperad vuestro turno en la cola")

# matriz y vector aleatorios, siempre los mismos gracias a la semilla
rng = np.random.default_rng(42)
A = rng.random((n, n))
b = rng.random(n)

# tiempo de reloj antes y después de resolver el sistema
t0 = time.perf_counter()
x = np.linalg.solve(A, b)
t = time.perf_counter() - t0

# operaciones que cuenta HPL para un sistema de n ecuaciones
flops = (2 / 3) * n**3 + (3 / 2) * n**2
# comprobación: el residuo es muy pequeño si la solución es correcta
residuo = np.linalg.norm(A @ x - b) / np.linalg.norm(b)

print(f"hilos={hilos}  n={n}  tiempo={t:.2f} s  "
      f"rendimiento={flops / t / 1e9:.1f} GFLOP/s  residuo={residuo:.1e}")

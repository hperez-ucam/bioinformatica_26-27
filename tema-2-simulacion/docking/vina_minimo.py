#!/usr/bin/env python3
"""Acoplamiento molecular mínimo con la API de Python de AutoDock Vina.

Proteasa del VIH-1 (PDB 1HVR) con su inhibidor XK263. Los ficheros de sistema/ los
genera preparar_minimo.py; el RMSD frente al ligando del cristal, calcular_rmsd.py.
La búsqueda es estocástica: con otra semilla o con otra exhaustividad puede devolver
otra pose y otra afinidad.

Uso (desde este directorio):
    python3 preparar_minimo.py 1HVR XK2   # una sola vez: escribe sistema/
    python3 vina_minimo.py 1 8       # semilla 1, exhaustividad 8
    python3 vina_minimo.py 777 8     # otra semilla
    python3 vina_minimo.py 1 32      # más búsqueda
"""
import sys
from pathlib import Path
from vina import Vina

# semilla y exhaustividad, desde la línea de órdenes
semilla, exh = int(sys.argv[1]), int(sys.argv[2])
v = Vina(sf_name="vina", seed=semilla, verbosity=0)
# receptor y ligando preparados, en formato PDBQT
v.set_receptor("sistema/1HVR_rec.pdbqt")
v.set_ligand_from_file("sistema/1HVR_lig.pdbqt")
# caja de 22 Å centrada en el ligando del cristal
centro = [-9.192, 15.906, 27.946]
v.compute_vina_maps(center=centro, box_size=[22] * 3)
# búsqueda estocástica: Monte Carlo + BFGS local
v.dock(exhaustiveness=exh, n_poses=5)
# afinidad (kcal/mol) de las cinco mejores poses
for i, e in enumerate(v.energies()[:, 0], 1):
    print(f"pose {i}: {e:6.2f} kcal/mol")
Path("salidas").mkdir(exist_ok=True)
salida = f"salidas/poses_s{semilla}_e{exh}.pdbqt"
v.write_poses(salida, n_poses=5, overwrite=True)

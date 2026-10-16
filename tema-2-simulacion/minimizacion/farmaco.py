#!/usr/bin/env python3
"""Minimización de un farmaco flexible, con el codigo más corto posible.

El butano tiene un grado de libertad y tres mínimos. Un farmaco tiene varios
enlaces que giran, y el número de mínimos crece tan deprisa que ya no se puede
dibujar la superficie: solo se puede muestrear.

La idea del guion es sencilla: generar muchos puntos de partida al azar,
minimizar cada uno, y ver cuantas energías distintas salen. Si una sola
minimización bastara, todas darian lo mismo.

Uso:
    python3 farmaco.py                         # ibuprofeno, 20 confórmeros
    python3 farmaco.py -n 50
    python3 farmaco.py --smiles "CC(=O)Oc1ccccc1C(=O)O" --nombre aspirina
"""

import argparse
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors

# RDKit avisa de cosas que aquí no importan
RDLogger.DisableLog("rdApp.*")

IBUPROFENO = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
RT = 8.314e-3 * 300      # kJ/mol a 300 K
KCAL = 4.184             # MMFF94 devuelve kcal/mol


def main():
    ap = argparse.ArgumentParser(description="Minimizacion de un farmaco flexible")
    ap.add_argument("--smiles", default=IBUPROFENO)
    ap.add_argument("--nombre", default="ibuprofeno")
    ap.add_argument("-n", type=int, default=20, help="número de puntos de partida")
    ap.add_argument("--semilla", type=int, default=2026)
    ap.add_argument("--sin-figura", action="store_true", help="no dibujar nada")
    args = ap.parse_args()

    # de texto SMILES a molécula
    mol = Chem.MolFromSmiles(args.smiles)
    # los hidrógenos hacen falta: sin ellos el campo de fuerza no vale
    mol = Chem.AddHs(mol)

    print(f"MOLÉCULA: {args.nombre}   {args.smiles}")
    print(f"  átomos (con H): {mol.GetNumAtoms()}")
    # cuenta los enlaces que pueden girar: son los grados de libertad de verdad
    print(f"  enlaces rotables: {Descriptors.NumRotatableBonds(Chem.RemoveHs(mol))}")

    # genera n geometrias 3D distintas. Es un método estocastico: con otra
    # semilla saldrian otros puntos de partida
    AllChem.EmbedMultipleConfs(mol, numConfs=args.n, randomSeed=args.semilla)
    # minimiza cada una con MMFF94; devuelve una pareja (0 si convergió, 1 si
    # agotó las iteraciones; energía en kcal/mol)
    resultados = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000)

    # se queda solo con las energías y las pasa a kJ/mol
    energias = np.array([e * KCAL for _, e in resultados])
    # si RDKit no ha podido generar ninguna geometría, avisar y salir
    if len(energias) == 0:
        raise SystemExit("RDKit no ha podido generar confórmeros para ese SMILES")
    # todo se mide respecto al mejor resultado obtenido
    relativas = np.sort(energias - energias.min())

    print(f"\nRESULTADO DE LAS {len(energias)} MINIMIZACIONES (kJ/mol sobre la mejor)")
    print("  " + "  ".join(f"{x:.1f}" for x in relativas))

    # dos resultados cuentan como el mismo mínimo si su energía difiere en
    # menos de medio kJ/mol
    distintas = len(np.unique(np.round(relativas / 0.5)))
    print(f"\n  energías distintas: {distintas} de {len(energias)} intentos")
    print(f"  dispersión: {relativas.max():.1f} kJ/mol = {relativas.max()/RT:.1f} RT")

    print("\nLo que hay que leer: una sola minimización habría devuelto una")
    print("cualquiera de estas energías, sin avisar de que había otras. Por eso")
    print("el cribado virtual parte de muchos confórmeros y no de uno.")

    if args.sin_figura:
        return

    import matplotlib
    # "Agg" dibuja a fichero sin abrir ninguna ventana: vale en un servidor sin pantalla
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 3.6), dpi=150)
    # un histograma de las energías relativas, en cajas de 0,5 kJ/mol
    ax.hist(relativas, bins=np.arange(0, relativas.max() + 1, 0.5), color="#8FA9BE", edgecolor="#004379")
    ax.set_xlabel("energía tras minimizar (kJ/mol respecto al mínimo más bajo)")
    ax.set_ylabel("confórmeros")
    fig.tight_layout()
    fig.savefig(f"{args.nombre}-minimos.png", dpi=150)
    print(f"\nfigura guardada en {args.nombre}-minimos.png")


if __name__ == "__main__":
    main()

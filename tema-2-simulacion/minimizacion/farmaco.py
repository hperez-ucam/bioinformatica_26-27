#!/usr/bin/env python3
"""Minimizacion de un farmaco flexible, con el codigo mas corto posible.

El butano tiene un grado de libertad y tres minimos. Un farmaco tiene varios
enlaces que giran, y el numero de minimos crece tan deprisa que ya no se puede
dibujar la superficie: solo se puede muestrear.

La idea del guion es sencilla: generar muchos puntos de partida al azar,
minimizar cada uno, y ver cuantas energias distintas salen. Si una sola
minimizacion bastara, todas darian lo mismo.

Uso:
    python3 farmaco.py                         # ibuprofeno, 20 confórmeros
    python3 farmaco.py -n 50
    python3 farmaco.py --smiles "CC(=O)Oc1ccccc1C(=O)O" --nombre aspirina
"""

import argparse
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors

# RDKit avisa de cosas que aqui no importan
RDLogger.DisableLog("rdApp.*")

IBUPROFENO = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
RT = 8.314e-3 * 300      # kJ/mol a 300 K
KCAL = 4.184             # MMFF94 devuelve kcal/mol


def main():
    ap = argparse.ArgumentParser(description="Minimizacion de un farmaco flexible")
    ap.add_argument("--smiles", default=IBUPROFENO)
    ap.add_argument("--nombre", default="ibuprofeno")
    ap.add_argument("-n", type=int, default=20, help="numero de puntos de partida")
    ap.add_argument("--semilla", type=int, default=2026)
    args = ap.parse_args()

    # de texto SMILES a molecula
    mol = Chem.MolFromSmiles(args.smiles)
    # los hidrogenos hacen falta: sin ellos el campo de fuerza no vale
    mol = Chem.AddHs(mol)

    print(f"MOLECULA: {args.nombre}   {args.smiles}")
    print(f"  atomos (con H): {mol.GetNumAtoms()}")
    # cuenta los enlaces que pueden girar: son los grados de libertad de verdad
    print(f"  enlaces rotables: {Descriptors.NumRotatableBonds(Chem.RemoveHs(mol))}")

    # genera n geometrias 3D distintas. Es un metodo estocastico: con otra
    # semilla saldrian otros puntos de partida
    AllChem.EmbedMultipleConfs(mol, numConfs=args.n, randomSeed=args.semilla)
    # minimiza cada una con MMFF94; devuelve una pareja (convergio, energia)
    resultados = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000)

    # se queda solo con las energias y las pasa a kJ/mol
    energias = np.array([e * KCAL for _, e in resultados])
    # todo se mide respecto al mejor resultado obtenido
    relativas = np.sort(energias - energias.min())

    print(f"\nRESULTADO DE LAS {args.n} MINIMIZACIONES (kJ/mol sobre la mejor)")
    print("  " + "  ".join(f"{x:.1f}" for x in relativas))

    # dos resultados cuentan como el mismo minimo si su energia difiere en
    # menos de medio kJ/mol
    distintas = len(np.unique(np.round(relativas / 0.5)))
    print(f"\n  energias distintas: {distintas} de {args.n} intentos")
    print(f"  dispersion: {relativas.max():.1f} kJ/mol = {relativas.max()/RT:.1f} RT")

    print("\nLo que hay que leer: una sola minimizacion habria devuelto una")
    print("cualquiera de estas energias, sin avisar de que habia otras. Por eso")
    print("el cribado virtual parte de muchos confórmeros y no de uno.")


if __name__ == "__main__":
    main()

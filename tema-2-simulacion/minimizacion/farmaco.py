#!/usr/bin/env python3
"""Minimizacion de energia de un farmaco flexible.

El butano tiene un grado de libertad y tres minimos. Un farmaco tiene varios
enlaces rotables y el numero de minimos crece de forma explosiva: ya no se puede
dibujar la superficie, solo muestrearla.

Lo que hace este guion:

  1. Construye la molecula a partir de su SMILES y le anade los hidrogenos.
  2. Genera N confórmeros de partida distintos (ETKDG, que es un metodo
     estocastico: cada semilla da un conjunto distinto).
  3. Minimiza cada uno con el campo de fuerza MMFF94.
  4. Agrupa los resultados por energia y por RMSD, y cuenta cuantos minimos
     *distintos* han aparecido.

La conclusion es la misma que en el butano, pero con consecuencias practicas:
una sola minimizacion no dice casi nada, porque solo informa del minimo que
tenia debajo el punto de partida. Por eso el cribado virtual parte de muchos
confórmeros y por eso el docking usa busqueda estocastica.

Uso:
    python3 farmaco.py                          # ibuprofeno, 20 confórmeros
    python3 farmaco.py --smiles "CC(=O)Oc1ccccc1C(=O)O" --nombre aspirina
    python3 farmaco.py -n 50 --sin-figura
"""

import argparse
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign, Descriptors

RDLogger.DisableLog("rdApp.*")

IBUPROFENO = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
RT = 8.314e-3 * 300            # kJ/mol a 300 K
KCAL = 4.184                   # MMFF94 devuelve kcal/mol


def preparar(smiles):
    """De SMILES a molecula 3D con hidrogenos."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SystemExit(f"SMILES no valido: {smiles}")
    return Chem.AddHs(mol)


def enlaces_rotables(mol):
    return Descriptors.NumRotatableBonds(Chem.RemoveHs(mol))


def generar_y_minimizar(mol, n, semilla=2026):
    """Genera n confórmeros y minimiza cada uno. Devuelve energias en kJ/mol.

    EmbedMultipleConfs es estocastico: con otra semilla saldrian otros puntos de
    partida y, por tanto, otros minimos. MMFFOptimizeMoleculeConfs hace el
    descenso hasta que la fuerza baja del umbral, igual que en el butano.
    """
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=n, randomSeed=semilla,
                                     pruneRmsThresh=-1.0)
    resultados = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000)
    energias = np.array([e * KCAL for convergido, e in resultados])
    convergidos = np.array([convergido == 0 for convergido, _ in resultados])
    return list(ids), energias, convergidos


def agrupar(mol, ids, energias, tol_e=0.5, tol_rmsd=0.5):
    """Agrupa confórmeros en minimos distintos: misma energia y mismo RMSD.

    Dos resultados son el mismo minimo si su energia difiere en menos de tol_e
    kJ/mol y se superponen por debajo de tol_rmsd angstrom.
    """
    orden = np.argsort(energias)
    grupos = []
    for i in orden:
        for g in grupos:
            j = g[0]
            if abs(energias[i] - energias[j]) < tol_e:
                rmsd = rdMolAlign.GetBestRMS(Chem.Mol(mol), Chem.Mol(mol),
                                             prbId=int(ids[i]), refId=int(ids[j]))
                if rmsd < tol_rmsd:
                    g.append(i)
                    break
        else:
            grupos.append([i])
    return grupos


def main():
    ap = argparse.ArgumentParser(description="Minimizacion de un farmaco flexible")
    ap.add_argument("--smiles", default=IBUPROFENO)
    ap.add_argument("--nombre", default="ibuprofeno")
    ap.add_argument("-n", type=int, default=20, help="numero de confórmeros de partida")
    ap.add_argument("--semilla", type=int, default=2026)
    ap.add_argument("--sin-figura", action="store_true")
    args = ap.parse_args()

    mol = preparar(args.smiles)
    print(f"MOLECULA: {args.nombre}   {args.smiles}")
    print(f"  atomos (con H): {mol.GetNumAtoms()}")
    print(f"  enlaces rotables: {enlaces_rotables(mol)}")

    ids, energias, convergidos = generar_y_minimizar(mol, args.n, args.semilla)
    print(f"\n{len(ids)} confórmeros generados; {convergidos.sum()} convergieron")

    rel = energias - energias.min()
    print(f"\nENERGIAS TRAS MINIMIZAR (kJ/mol, respecto al mejor)")
    print(f"  minimo        {rel.min():8.2f}")
    print(f"  maximo        {rel.max():8.2f}")
    print(f"  mediana       {np.median(rel):8.2f}")
    print(f"  dispersion    {rel.max() - rel.min():8.2f}   = {(rel.max() - rel.min()) / RT:.1f} RT")

    grupos = agrupar(mol, ids, energias)
    print(f"\nMINIMOS DISTINTOS ENCONTRADOS: {len(grupos)} de {len(ids)} intentos")
    print(f"{'#':>3} {'E rel (kJ/mol)':>15} {'veces':>7}  poblacion a 300 K")
    e_grupo = np.array([energias[g[0]] for g in grupos])
    pesos = np.exp(-(e_grupo - e_grupo.min()) / RT)
    pesos /= pesos.sum()
    for k, (g, p) in enumerate(zip(grupos, pesos), 1):
        if k > 10:
            print(f"    ... y {len(grupos) - 10} minimos mas")
            break
        print(f"{k:>3} {energias[g[0]] - energias.min():>15.2f} {len(g):>7}  {100 * p:5.1f} %")

    print(f"\nLo que hay que leer: {args.n} puntos de partida han dado {len(grupos)} minimos")
    print(f"distintos, repartidos en {rel.max() - rel.min():.1f} kJ/mol. Una sola")
    print("minimizacion habria devuelto uno cualquiera de ellos, sin avisar.")

    if args.sin_figura:
        return

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4), dpi=150)
    ax.hist(rel, bins=max(8, len(grupos)), color="#8FA9BE", edgecolor="#004379")
    ax.axvline(0, color="#9A7200", lw=2)
    ax.set_xlabel("energia tras minimizar, respecto al mejor (kJ/mol)")
    ax.set_ylabel("numero de confórmeros")
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{args.nombre}-minimos.png", dpi=150)
    print(f"\nfigura guardada en {args.nombre}-minimos.png")


if __name__ == "__main__":
    main()

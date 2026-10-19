#!/usr/bin/env python3
"""Minimización de una molécula flexible con RDKit desde varios puntos de partida.

Versión corta de farmaco.py. El ibuprofeno tiene cuatro enlaces rotables: se generan
al azar diez geometrías de partida (ETKDG), cada una se minimiza con el campo de
fuerza MMFF94 y se comparan las energías finales.

Uso:
    python3 ibuprofeno_rdkit.py
"""
from rdkit import Chem
from rdkit.Chem import AllChem

# ibuprofeno desde su SMILES, con hidrógenos
smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
# diez geometrías 3D de partida, al azar
AllChem.EmbedMultipleConfs(mol, 10, randomSeed=2026)
# MMFF94 minimiza cada una: (convergencia, kcal/mol)
res = AllChem.MMFFOptimizeMoleculeConfs(mol)
# energías en kJ/mol, ordenadas y relativas a la mejor
E = sorted(e * 4.184 for _, e in res)
dE = [e - E[0] for e in E]
print("kJ/mol sobre la mejor:")
# se imprimen de cinco en cinco
for i in range(0, len(dE), 5):
    print(" ".join(f"{x:4.1f}" for x in dE[i:i + 5]))
# energías distintas, con una tolerancia de 0,5 kJ/mol
print("distintas:", len({round(x * 2) for x in dE}), "de", len(dE))

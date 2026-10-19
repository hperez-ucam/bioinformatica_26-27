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
# energías relativas a la mejor, en kJ/mol
E = [e * 4.184 for _, e in res]
for e in sorted(E):
    print(f"{e - min(E):5.1f} kJ/mol")

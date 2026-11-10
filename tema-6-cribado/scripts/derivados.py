#!/usr/bin/env python3
# derivados.py - diseña derivados de un ligando de referencia cambiando un sustituyente.
#   python3 derivados.py > derivados.smi
#
# Parte del erlotinib escrito con un hueco {R} en la posición 3 del anillo de anilina
# y lo rellena con distintos sustituyentes. Es la vía de "diseño de derivados" del enunciado.
from rdkit import Chem

# erlotinib con un hueco {R} donde estaba el grupo etinilo (C#C)
plantilla = "COCCOc1cc2ncnc(Nc3cccc({R})c3)c2cc1OCCOC"
# nombre corto del derivado y sustituyente que se pone en el hueco
sustituyentes = {"H": "[H]", "F": "F", "Cl": "Cl", "Br": "Br", "Me": "C",
                 "OMe": "OC", "CF3": "C(F)(F)F", "NO2": "[N+](=O)[O-]",
                 "NH2": "N", "OH": "O"}

for nombre, grupo in sustituyentes.items():
    # se escribe el SMILES completo del derivado sustituyendo el hueco
    smiles = plantilla.format(R=grupo)
    # RDKit lee el SMILES: si es erróneo, MolFromSmiles devuelve None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SystemExit(f"SMILES no válido para {nombre}: {smiles}")
    # se imprime en formato "SMILES nombre", el que leen los demás scripts
    print(f"{Chem.MolToSmiles(mol)}\tderiv_{nombre}")

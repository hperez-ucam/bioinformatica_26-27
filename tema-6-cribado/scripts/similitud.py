#!/usr/bin/env python3
# similitud.py - LBVS: similitud de Tanimoto entre cada molécula y una o más referencias.
#   python3 similitud.py libreria.smi erlotinib [gefitinib ...] > ranking_lbvs.tsv
#
# Huella de Morgan de radio 2 (equivale a ECFP4) de 2048 bits. La similitud con varias
# referencias es la máxima de todas (la molécula se parece a la que más se le parezca).
import sys
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

# fichero de moléculas (SMILES, tabulador, nombre) y nombres de las referencias
fichero, referencias = sys.argv[1], sys.argv[2:]
# diccionario nombre -> huella de Morgan (vector de bits)
huellas, orden = {}, []
for linea in open(fichero):
    if not linea.strip() or linea.startswith("#"):
        continue
    smiles, nombre = linea.rstrip("\n").split("\t")[:2]
    # la huella codifica los entornos de cada átomo hasta 2 enlaces de distancia
    huellas[nombre] = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles), 2, nBits=2048)
    orden.append(nombre)

filas = []
for nombre in orden:
    # Tanimoto con cada referencia: bits comunes / bits en total
    sims = [DataStructs.TanimotoSimilarity(huellas[nombre], huellas[r])
            for r in referencias]
    filas.append((nombre, max(sims)))
# de mayor a menor similitud
filas.sort(key=lambda f: -f[1])
for nombre, s in filas:
    print(f"{nombre}\t{s:.3f}")

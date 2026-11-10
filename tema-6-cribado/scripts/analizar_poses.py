#!/usr/bin/env python3
# analizar_poses.py - dónde cae la mejor pose de cada ligando y qué residuos toca.
#   python3 analizar_poses.py sistema/1M17_rec.pdb sistema/1M17_lig_cristal.sdf resultados/
#
# Para la mejor pose (modo 1) de cada ligando: distancia de su centro al del ligando cristalográfico,
# residuos a menos de 4 Å, posibles enlaces de hidrógeno (N u O a menos de 3,3 Å) y cuánto
# se parecen sus contactos a los del ligando cristalográfico (índice de Jaccard).
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from meeko import PDBQTMolecule, RDKitMolCreate

receptor, cristal, directorio = sys.argv[1], sys.argv[2], Path(sys.argv[3])

# átomos pesados del receptor: coordenadas, elemento y residuo (nombre y número)
coord, elem, res = [], [], []
for l in open(receptor):
    if l.startswith("ATOM") and l[76:78].strip() != "H":
        coord.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        elem.append(l[76:78].strip())
        res.append(f"{l[17:20]}{int(l[22:26])}")
coord, elem, res = np.array(coord), np.array(elem), np.array(res)

def contactos(xyz, elementos):
    # distancias: átomos del ligando (filas) por átomos del receptor
    d = np.linalg.norm(xyz[:, None, :] - coord[None, :, :], axis=2)
    # residuos con algún átomo a menos de 4 Å del ligando
    cerca = sorted(set(res[(d < 4.0).any(axis=0)]),
                   key=lambda r: int(r[3:]))
    # N u O del ligando y del receptor a menos de 3,3 Å: posible enlace H
    polar_lig = np.array([e in ("N", "O") for e in elementos])
    polar_rec = np.isin(elem, ["N", "O"])
    cerca_polar = (d < 3.3) & polar_lig[:, None] & polar_rec[None, :]
    hb = sorted(set(res[cerca_polar.any(axis=0)]),
                key=lambda r: int(r[3:]))
    return cerca, hb

# ligando cristalográfico: referencia de contactos y de posición
ref = Chem.MolFromMolFile(cristal)
xyz_ref = ref.GetConformer().GetPositions()
c_ref, hb_ref = contactos(xyz_ref, [a.GetSymbol() for a in ref.GetAtoms()])
print(f"cristal\tcentro=0.0 Å\tcontactos={len(c_ref)}\tH-enlaces={','.join(hb_ref)}")
print("nombre\tE\tdist_centro_Å\tcontactos\tJaccard_cristal\tH-enlaces")

for fichero in sorted(directorio.glob("*_poses.pdbqt")):
    nombre = fichero.name.replace("_poses.pdbqt", "")
    # energía del modo 1 en las líneas «REMARK VINA RESULT»
    energia = [l.split()[3] for l in open(fichero) if l.startswith("REMARK VINA RESULT")][0]
    # meeko reconstruye la molécula; la primera conformación es el mejor modo
    mol = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule.from_file(str(fichero), skip_typing=True))[0]
    mol = Chem.RemoveHs(Chem.Mol(mol, confId=mol.GetConformers()[0].GetId()))
    xyz = mol.GetConformer().GetPositions()
    cerca, hb = contactos(xyz, [a.GetSymbol() for a in mol.GetAtoms()])
    # distancia entre el centro de la pose y el del ligando cristalográfico
    dist = np.linalg.norm(xyz.mean(axis=0) - xyz_ref.mean(axis=0))
    # Jaccard: residuos comunes / residuos en total entre pose y cristal
    jac = len(set(cerca) & set(c_ref)) / len(set(cerca) | set(c_ref))
    print(f"{nombre}\t{energia}\t{dist:.1f}\t{len(cerca)}\t{jac:.2f}\t{','.join(hb)}")

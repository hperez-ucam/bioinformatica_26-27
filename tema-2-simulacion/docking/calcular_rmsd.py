#!/usr/bin/env python3
# calcular_rmsd.py - RMSD de cada pose de Vina frente al ligando cristalográfico.
#   python calcular_rmsd.py sistema/1HVR_lig_cristal.sdf salidas/A_poses.pdbqt [...]
#
# RMSD sobre átomos pesados, en el sistema de referencia del receptor (sin
# superponer: la pose acoplada y la cristalográfica comparten coordenadas) y
# teniendo en cuenta la simetría molecular: si dos grupos son equivalentes
# (los dos naftilos o los dos bencilos de XK263), se toma la correspondencia
# de átomos que da el RMSD mínimo.
import sys
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from meeko import PDBQTMolecule, RDKitMolCreate

# Referencia: ligando del cristal, ya con órdenes de enlace (preparar_minimo.py).
ref = Chem.MolFromMolFile(sys.argv[1])

for fichero in sys.argv[2:]:
    # meeko reconstruye una molécula de RDKit por cada modo del .pdbqt de salida.
    poses = PDBQTMolecule.from_file(fichero, skip_typing=True)
    mol = RDKitMolCreate.from_pdbqt_mol(poses)[0]
    # La afinidad de cada modo (kcal/mol) está en las líneas «REMARK VINA RESULT».
    afinidades = [l.split()[3] for l in open(fichero) if l.startswith("REMARK VINA RESULT")]
    for i, conf in enumerate(mol.GetConformers(), start=1):
        pose = Chem.RemoveHs(Chem.Mol(mol, confId=conf.GetId()))
        # CalcRMS no superpone las moléculas y prueba todas las permutaciones por simetría.
        rmsd = rdMolAlign.CalcRMS(pose, ref)
        print(f"{fichero}\tmodo {i}\t{afinidades[i - 1]} kcal/mol\tRMSD = {rmsd:.2f} Å")

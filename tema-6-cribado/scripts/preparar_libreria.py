#!/usr/bin/env python3
# preparar_libreria.py - de un fichero de SMILES a un .pdbqt por molécula, con informe.
#   python3 preparar_libreria.py libreria.smi ligandos/
#
# El fichero de entrada tiene una molécula por línea: SMILES, tabulador y nombre.
# Para cada molécula: comprueba el SMILES, detecta duplicados y estereocentros sin
# definir, protona a pH 7,4 con Open Babel, genera una conformación 3D y la
# convierte a .pdbqt con meeko. Escribe ligandos/informe.tsv con lo que ha hecho.
import subprocess
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy

# fichero de entrada y directorio de salida, desde la línea de órdenes
entrada, salida = Path(sys.argv[1]), Path(sys.argv[2])
salida.mkdir(exist_ok=True)
# SMILES canónicos ya vistos, para detectar moléculas repetidas con otro nombre
vistos = {}
filas = ["nombre\testado\tPM\tatomos_pesados\ttorsiones\tcarga\testereocentros_sin_definir"]

for linea in entrada.read_text().splitlines():
    # se salta una línea vacía o un comentario
    if not linea.strip() or linea.startswith("#"):
        continue
    smiles, nombre = linea.split("\t")[:2]
    mol = Chem.MolFromSmiles(smiles)
    # un SMILES que RDKit no puede leer se anota y se continúa
    if mol is None:
        filas.append(f"{nombre}\tSMILES_no_valido\t\t\t\t\t")
        continue
    # el SMILES canónico es único por molécula: sirve para encontrar duplicados
    canonico = Chem.MolToSmiles(mol)
    if canonico in vistos:
        filas.append(f"{nombre}\tduplicado_de_{vistos[canonico]}\t\t\t\t\t")
        continue
    vistos[canonico] = nombre
    # estereocentros que el SMILES no define: habría que elegir isómero
    centros = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    sin_definir = sum(1 for _, q in centros if q == "?")
    # protonación a pH 7,4 con Open Babel (-p): aminas a catión, ácidos a anión
    orden = ["obabel", f"-:{canonico}", "-osmi", "-p", "7.4"]
    protonado = subprocess.run(orden, capture_output=True,
                               text=True).stdout.split()[0]
    mol = Chem.MolFromSmiles(protonado)
    # hidrógenos explícitos y 3D con ETKDG (semilla fija: repetible)
    mol3d = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol3d, randomSeed=42)
    # minimización con MMFF94 (Tema 2.3): lleva a un mínimo local
    AllChem.MMFFOptimizeMolecule(mol3d)
    # meeko: tipos AutoDock, cargas de Gasteiger y torsiones
    setup = MoleculePreparation().prepare(mol3d)[0]
    pdbqt, correcto, error = PDBQTWriterLegacy.write_string(setup)
    if not correcto:
        filas.append(f"{nombre}\tmeeko_falla\t\t\t\t\t")
        continue
    (salida / f"{nombre}.pdbqt").write_text(pdbqt)
    # torsiones activas: cada una abre y cierra una rama (BRANCH)
    torsiones = pdbqt.count("BRANCH") // 2
    carga = Chem.GetFormalCharge(mol)
    filas.append(f"{nombre}\tok\t{Descriptors.MolWt(mol):.1f}\t{mol.GetNumHeavyAtoms()}\t{torsiones}\t{carga:+d}\t{sin_definir}")

(salida / "informe.tsv").write_text("\n".join(filas) + "\n")
print("\n".join(filas))

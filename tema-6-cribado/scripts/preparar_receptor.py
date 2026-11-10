#!/usr/bin/env python3
# preparar_receptor.py - de un código PDB a receptor .pdbqt, ligando cristalográfico y caja.
#   python3 preparar_receptor.py 1M17 AQ4
#
# Es la parte de receptor de preparar_minimo.py (Tema 2), sin generar el ligando de partida.
# Escribe en sistema/: CODIGO_rec.pdbqt, CODIGO_lig_cristal.sdf y caja.txt
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

# código PDB y código de tres letras del ligando
pdb, lig = sys.argv[1].upper(), sys.argv[2].upper()
# directorio de salida, junto al script
destino = Path(__file__).resolve().parent / "sistema"
destino.mkdir(exist_ok=True)

# descarga de la entrada del PDB (formato PDB clásico)
url = f"https://files.rcsb.org/download/{pdb}.pdb"
texto = urllib.request.urlopen(url).read().decode()
(destino / f"{pdb}.pdb").write_text(texto)
lineas = texto.splitlines()

# receptor: solo las líneas ATOM (sin aguas, iones ni cofactores)
receptor = [l for l in lineas if l.startswith("ATOM")]
(destino / f"{pdb}_rec.pdb").write_text("\n".join(receptor) + "\nEND\n")
aguas = sum(l.startswith("HETATM") and l[17:20] == "HOH" for l in lineas)
print(f"receptor: {len(receptor)} átomos; aguas descartadas: {aguas}")

# ligando cristalográfico: líneas HETATM con su código en las columnas 18-20
bloque = "\n".join(l for l in lineas
                   if l.startswith("HETATM") and l[17:20] == lig)
cristal = Chem.MolFromPDBBlock(bloque, removeHs=False)
# un PDB no guarda órdenes de enlace: se toman del SMILES del diccionario
# de componentes químicos del PDB
url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{lig}"
ficha = json.load(urllib.request.urlopen(url))
smiles = ficha["rcsb_chem_comp_descriptor"]["SMILES_stereo"]
plantilla = Chem.MolFromSmiles(smiles)
cristal = AllChem.AssignBondOrdersFromTemplate(plantilla, cristal)
Chem.MolToMolFile(cristal, str(destino / f"{pdb}_lig_cristal.sdf"))
print(f"ligando {lig}: {cristal.GetNumAtoms()} átomos pesados, SMILES {smiles}")

# centro de la caja: media de las coordenadas del ligando del cristal
centro = cristal.GetConformer().GetPositions().mean(axis=0)
(destino / "caja.txt").write_text("%.3f %.3f %.3f\n" % tuple(centro))
print("centro de la caja: %.3f %.3f %.3f" % tuple(centro))

# receptor en .pdbqt con meeko: hidrógenos polares, cargas y tipos AutoDock.
# --default_altloc A: si un residuo tiene dos posiciones alternativas (en
# 1M17, 751 y 831) se queda con la A; sin la opción, meeko se detiene
subprocess.run(["mk_prepare_receptor.py", "--read_pdb",
                str(destino / f"{pdb}_rec.pdb"),
                "-o", str(destino / f"{pdb}_rec"), "-p",
                "--default_altloc", "A"], check=True)

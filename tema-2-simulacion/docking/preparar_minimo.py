#!/usr/bin/env python3
"""preparar_minimo.py - de un código PDB a los ficheros .pdbqt que usa vina_minimo.py.

Versión en Python de preparar.sh, con meeko y RDKit en lugar de AutoDockTools.
Escribe en el directorio sistema/:
    CODIGO.pdb             entrada original del PDB
    CODIGO_rec.pdb         receptor (solo líneas ATOM)
    CODIGO_rec.pdbqt       receptor con hidrógenos polares, cargas y tipos AutoDock
    CODIGO_lig_cristal.sdf ligando cristalográfico con órdenes de enlace (referencia del RMSD)
    CODIGO_lig.pdbqt       ligando que se acopla (conformación generada desde SMILES)
    caja.txt               centro geométrico del ligando cristalográfico

    python preparar_minimo.py 1HVR XK2
"""
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

pdb, lig = sys.argv[1].upper(), sys.argv[2].upper()
destino = Path(__file__).resolve().parent / "sistema"
destino.mkdir(exist_ok=True)

# 1. Descarga de la entrada del PDB (formato PDB clásico).
texto = urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb}.pdb").read().decode()
(destino / f"{pdb}.pdb").write_text(texto)
lineas = texto.splitlines()

# 2. Receptor: solo las líneas ATOM, sin aguas ni heteroátomos.
#    En 1HVR esto descarta también CSO 67 (cisteína oxidada, HETATM) de ambas
#    cadenas: queda un hueco en la cadena lejos del centro activo (> 15 Å).
receptor = [l for l in lineas if l.startswith("ATOM")]
(destino / f"{pdb}_rec.pdb").write_text("\n".join(receptor) + "\nEND\n")

# 3. Ligando cristalográfico: líneas HETATM con el código de tres letras indicado.
#    Las columnas 18-20 del formato PDB contienen el nombre del residuo.
bloque = "\n".join(l for l in lineas if l.startswith("HETATM") and l[17:20] == lig)
cristal = Chem.MolFromPDBBlock(bloque, removeHs=False)

# Un fichero PDB no guarda órdenes de enlace: se asignan a partir del SMILES
# del diccionario de componentes químicos del PDB (CCD).
url_ccd = f"https://data.rcsb.org/rest/v1/core/chemcomp/{lig}"
smiles = json.load(urllib.request.urlopen(url_ccd))["rcsb_chem_comp_descriptor"]["SMILES_stereo"]
plantilla = Chem.MolFromSmiles(smiles)
cristal = AllChem.AssignBondOrdersFromTemplate(plantilla, cristal)
Chem.MolToMolFile(cristal, str(destino / f"{pdb}_lig_cristal.sdf"))

# 4. Centro de la caja de búsqueda: media de las coordenadas del ligando cristalográfico.
centro = cristal.GetConformer().GetPositions().mean(axis=0)
(destino / "caja.txt").write_text("%.3f %.3f %.3f\n" % tuple(centro))
print(f"centro de la caja: {centro.round(3)}")

# 5. Ligando de partida: conformación 3D generada desde el SMILES (ETKDG), con
#    hidrógenos explícitos. Así el acoplamiento no parte de la pose cristalográfica.
partida = Chem.AddHs(plantilla)
AllChem.EmbedMolecule(partida, randomSeed=42)
AllChem.MMFFOptimizeMolecule(partida)

# meeko asigna tipos AutoDock, cargas de Gasteiger y el árbol de torsiones.
# rigid_macrocycles=True mantiene rígido el anillo de siete miembros (urea cíclica).
preparador = MoleculePreparation(rigid_macrocycles=True)
setup = preparador.prepare(partida)[0]
pdbqt_lig, correcto, error = PDBQTWriterLegacy.write_string(setup)
if not correcto:
    sys.exit(error)
(destino / f"{pdb}_lig.pdbqt").write_text(pdbqt_lig)
print(f"ligando: {cristal.GetNumAtoms()} átomos pesados, "
      f"{pdbqt_lig.count('BRANCH') // 2} torsiones activas")

# 6. Receptor .pdbqt con meeko: añade hidrógenos polares y tipos AutoDock
#    a partir de plantillas de aminoácidos.
subprocess.run(["mk_prepare_receptor.py", "--read_pdb", str(destino / f"{pdb}_rec.pdb"),
                "-o", str(destino / f"{pdb}_rec"), "-p"], check=True)

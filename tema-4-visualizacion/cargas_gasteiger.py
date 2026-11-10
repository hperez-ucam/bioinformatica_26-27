"""cargas_gasteiger.py — Tema 4, ejercicio E1 (v4, 10/11/2026).

Dibuja una molécula pequeña en 2D coloreada por sus cargas parciales de Gasteiger
(rojo: negativa; azul: positiva) e imprime los tres átomos más negativos.

Se ejecuta en el servidor de clase, donde python3 ya trae RDKit:
      cd ~/bioinformatica_26-27/tema-4-visualizacion
      python3 cargas_gasteiger.py "SMILES" nombre
Ejemplo (oleocantal):
      python3 cargas_gasteiger.py "C/C=C(/C=O)[C@@H](CC=O)CC(=O)OCCc1ccc(O)cc1" oleocantal
La figura queda en nombre_gasteiger.png; para verla, se copia al portátil con scp
(en el portátil, no en el servidor):
      scp USUARIO@31.70.148.22:bioinformatica_26-27/tema-4-visualizacion/oleocantal_gasteiger.png .
"""

# sys lee los argumentos de la línea de órdenes
import sys

# Chem lee el SMILES; AllChem calcula las cargas; rdMolDescriptors da la fórmula
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
# SimilarityMaps colorea cada átomo según un valor; rdMolDraw2D dibuja la molécula en un PNG
from rdkit.Chem.Draw import SimilarityMaps, rdMolDraw2D

# primer argumento: el SMILES; segundo: el nombre que llevará el fichero de salida
smiles, nombre = sys.argv[1], sys.argv[2]

# convertir el SMILES en una molécula; si está mal escrito, RDKit devuelve None
mol = Chem.MolFromSmiles(smiles)
if mol is None:
    sys.exit(f"SMILES no válido: {smiles}")

# fórmula molecular, para comprobar que la molécula es la que se pretendía
print(nombre, rdMolDescriptors.CalcMolFormula(mol))

# cargas de Gasteiger: método empírico y rápido
AllChem.ComputeGasteigerCharges(mol)
# una carga por átomo, en el orden de la molécula
cargas = [a.GetDoubleProp("_GasteigerCharge")
          for a in mol.GetAtoms()]

# átomos de la carga más negativa a la más positiva
orden = sorted(range(len(cargas)), key=lambda i: cargas[i])
# los tres primeros: número de átomo, elemento y carga
for i in orden[:3]:
    simbolo = mol.GetAtomWithIdx(i).GetSymbol()
    print(f"  átomo {i:2d}  {simbolo}  {cargas[i]:+.3f}")

# lienzo de 900 x 650 píxeles para el dibujo
lienzo = rdMolDraw2D.MolDraw2DCairo(900, 650)
# colorear cada átomo con su carga: bwr_r pone el negativo en rojo y el positivo en azul
SimilarityMaps.GetSimilarityMapFromWeights(mol, cargas, draw2d=lienzo,
                                           colorMap="bwr_r", contourLines=0)
lienzo.FinishDrawing()
# guardar la imagen como nombre_gasteiger.png
with open(f"{nombre}_gasteiger.png", "wb") as f:
    f.write(lienzo.GetDrawingText())
print(f"  figura: {nombre}_gasteiger.png")

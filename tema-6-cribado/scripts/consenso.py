#!/usr/bin/env python3
# consenso.py - une el ranking de docking (SBVS) y el de similitud (LBVS) en una tabla.
#   python3 consenso.py libreria.smi ligandos/informe.tsv resultados/ranking_sbvs.tsv ranking_lbvs.tsv
#
# Columnas: energía, eficiencia de ligando (LE = -energía / átomos pesados), puesto en cada
# ranking, suma de puestos (consenso) y propiedades de la regla de Lipinski.
import sys
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski

smi, informe, sbvs, lbvs = sys.argv[1:5]

def leer(fichero):
    # lee un tsv "nombre<TAB>valor" y devuelve un diccionario nombre -> float
    return {l.split("\t")[0]: float(l.split("\t")[1]) for l in open(fichero) if l.strip()}

energia = leer(sbvs)
similitud = leer(lbvs)
# átomos pesados de cada molécula, del informe de preparar_libreria.py (columna 4)
pesados = {l.split("\t")[0]: int(l.split("\t")[3]) for l in open(informe).read().splitlines()[1:]
           if l.split("\t")[1] == "ok"}
# SMILES de cada molécula, para las propiedades
smiles = {l.split("\t")[1].strip(): l.split("\t")[0] for l in open(smi) if l.strip() and not l.startswith("#")}

nombres = [n for n in energia if n in similitud and n in pesados]
# puesto 1 = mejor: menor energía (SBVS) y mayor similitud (LBVS)
p_sbvs = {n: i + 1 for i, n in enumerate(sorted(nombres, key=lambda n: energia[n]))}
p_lbvs = {n: i + 1 for i, n in enumerate(sorted(nombres, key=lambda n: -similitud[n]))}
# eficiencia de ligando: energía por átomo pesado, en positivo
le = {n: -energia[n] / pesados[n] for n in nombres}
p_le = {n: i + 1 for i, n in enumerate(sorted(nombres, key=lambda n: -le[n]))}

print("nombre\tE\tLE\tpuesto_E\tpuesto_LE\tTanimoto\tpuesto_sim\tsuma_E_sim\tPM\tlogP\tHBD\tHBA\tviolaciones_Lipinski")
for n in sorted(nombres, key=lambda n: p_sbvs[n] + p_lbvs[n]):
    m = Chem.MolFromSmiles(smiles[n])
    pm, logp = Descriptors.MolWt(m), Crippen.MolLogP(m)
    hbd, hba = Lipinski.NumHDonors(m), Lipinski.NumHAcceptors(m)
    # Lipinski: PM > 500, logP > 5, más de 5 donantes o más de 10 aceptores
    viol = (pm > 500) + (logp > 5) + (hbd > 5) + (hba > 10)
    print(f"{n}\t{energia[n]:.1f}\t{le[n]:.2f}\t{p_sbvs[n]}\t{p_le[n]}\t{similitud[n]:.2f}\t{p_lbvs[n]}\t{p_sbvs[n]+p_lbvs[n]}\t{pm:.0f}\t{logp:.1f}\t{hbd}\t{hba}\t{viol}")

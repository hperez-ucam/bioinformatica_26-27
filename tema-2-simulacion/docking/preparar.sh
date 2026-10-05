#!/bin/bash
# preparar.sh - de una entrada del PDB a los dos ficheros .pdbqt que pide Vina.
#
# Usa las mismas ordenes del Tema 1: wget, grep, grep -v, wc y redireccion.
# Lo unico nuevo son las dos ultimas ordenes, que anaden al fichero la
# informacion que el programa de docking necesita y que un PDB no trae:
# cargas parciales e hidrogenos polares.
#
#   ./preparar.sh 1HVR                 # descarga y prepara la entrada 1HVR
#
set -euo pipefail

PDB="${1:?uso: ./preparar.sh CODIGO_PDB}"
PDB="${PDB^^}"

echo "== 1. descargar =="
wget -q "https://files.rcsb.org/download/${PDB}.pdb" -O "${PDB}.pdb"
echo "   ${PDB}.pdb: $(wc -l < "${PDB}.pdb") lineas"

echo "== 2. separar receptor y ligando con grep =="
# El receptor son las lineas ATOM; las aguas y los heteroatomos sobran.
grep '^ATOM' "${PDB}.pdb" | grep -v 'HOH' > "${PDB}_rec.pdb"
echo "   receptor: $(wc -l < "${PDB}_rec.pdb") atomos"

# El ligando es un HETATM que no es agua ni ion. Primero vemos que hay:
echo "   heteroatomos presentes:"
grep '^HETATM' "${PDB}.pdb" | cut -c18-20 | sort | uniq -c | sort -rn | head -5

# El codigo del ligando se pasa como segundo argumento si no es el primero.
LIG="${2:-$(grep '^HETATM' "${PDB}.pdb" | grep -v 'HOH' | cut -c18-20 \
            | sort | uniq -c | sort -rn | head -1 | awk '{print $2}')}"
echo "   ligando elegido: ${LIG}"
grep '^HETATM' "${PDB}.pdb" | grep "${LIG}" > "${PDB}_lig.pdb"
echo "   ligando: $(wc -l < "${PDB}_lig.pdb") atomos"

echo "== 3. centro de la caja: el centro geometrico del ligando =="
awk '{x+=substr($0,31,8); y+=substr($0,39,8); z+=substr($0,47,8); n++}
     END {printf "   centro  x=%.2f  y=%.2f  z=%.2f  (%d atomos)\n", x/n, y/n, z/n, n}' \
     "${PDB}_lig.pdb" | tee caja.txt

echo "== 4. anadir hidrogenos y cargas: de .pdb a .pdbqt =="
# prepare_receptor y prepare_ligand vienen con AutoDockTools / ADFR Suite.
# -A hydrogens anade los hidrogenos que falten; -U nphs_lps quita los apolares
# y los pares libres, que el modelo de Vina no usa.
prepare_receptor -r "${PDB}_rec.pdb" -o "${PDB}_rec.pdbqt" -A hydrogens -U nphs_lps
prepare_ligand   -l "${PDB}_lig.pdb" -o "${PDB}_lig.pdbqt" -A hydrogens

echo "== 5. que ha cambiado =="
echo "   --- primera linea del .pdb:"
head -1 "${PDB}_rec.pdb"
echo "   --- primera linea del .pdbqt:"
grep -m1 '^ATOM' "${PDB}_rec.pdbqt"
echo
echo "   Las dos columnas nuevas del final son la carga parcial del atomo y su"
echo "   tipo segun AutoDock. Esa informacion no esta en un PDB: la anade el"
echo "   programa de preparacion a partir del campo de fuerza."

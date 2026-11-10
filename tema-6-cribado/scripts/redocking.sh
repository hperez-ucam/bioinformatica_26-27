#!/bin/bash
# redocking.sh - acopla de nuevo el ligando cristalográfico a su receptor y mide el RMSD.
#   ./redocking.sh sistema/1M17_rec.pdbqt sistema/1M17_lig_cristal.sdf ligandos/erlotinib.pdbqt "x y z"
#
# Es el control del cribado: si el método no recupera la pose del cristal, el ranking de
# una librería no es fiable. Prueba tres semillas con exhaustividad 8, una con 32 y una
# caja desplazada 15 Å (calcular_rmsd.py es el del Tema 2).
set -euo pipefail

# receptor, ligando cristalográfico con enlaces, ligando a acoplar y centro de la caja
receptor=$1; cristal=$2; ligando=$3; centro=$4
read -r cx cy cz <<< "$centro"
mkdir -p redocking

# una ejecución de Vina: $1 semilla, $2 exhaustividad,
# $3 desplazamiento de la caja en x (Å), $4 etiqueta
correr() {
  # centro en x desplazado (0 si la caja está bien puesta)
  local x
  x=$(echo "$cx + $3" | bc -l)
  vina --receptor "$receptor" --ligand "$ligando" \
       --center_x "$x" --center_y "$cy" --center_z "$cz" \
       --size_x 22 --size_y 22 --size_z 22 \
       --exhaustiveness "$2" --seed "$1" --cpu 6 --num_modes 5 \
       --out "redocking/$4.pdbqt" > "redocking/$4.log"
}

# tabla de resultados: etiqueta, tiempo y RMSD del mejor modo y del mejor de los cinco
echo -e "caso\ttiempo_s\tE_modo1\tRMSD_modo1\tmejor_RMSD_5_modos"
for caso in "s1_e8 1 8 0" "s7_e8 7 8 0" "s42_e8 42 8 0" "s1_e32 1 32 0" "caja_15A 1 8 15"; do
  # se separan etiqueta, semilla, exhaustividad y desplazamiento
  read -r etiqueta semilla exh desp <<< "$caso"
  t0=$(date +%s.%N)
  correr "$semilla" "$exh" "$desp" "$etiqueta"
  t1=$(date +%s.%N)
  # una línea por modo; se queda con la del modo 1 y con el mínimo de los cinco RMSD
  python3 calcular_rmsd.py "$cristal" "redocking/$etiqueta.pdbqt" 2>/dev/null \
    | awk -v caso="$etiqueta" -v t="$(echo "$t1 - $t0" | bc -l)" \
      'BEGIN{FS="\t"} NR==1{e=$3; r1=$4} {sub(/RMSD = /,"",$4); sub(/ Å/,"",$4); if(m==""||$4<m)m=$4}
       END{sub(/RMSD = /,"",r1); printf "%s\t%.0f\t%s\t%s\t%s Å\n", caso, t, e, r1, m}'
done

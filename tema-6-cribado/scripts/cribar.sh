#!/bin/bash
# cribar.sh - acoplamiento de todos los ligandos de un directorio con Vina, seis a la vez.
#   ./cribar.sh sistema/1M17_rec.pdbqt "x y z" ligandos/ resultados/ [semilla] [exhaustividad]
#
# Cada ligando se acopla con un solo hilo (--cpu 1) y xargs -P 6 reparte los ligandos entre
# seis procesos, uno por núcleo del servidor (Tema 3). Escribe en resultados/: las poses
# (nombre_poses.pdbqt), el registro (nombre.log) y ranking_sbvs.tsv ordenado por energía.
set -euo pipefail

# argumentos: receptor, centro de la caja, directorio de ligandos, directorio de resultados
receptor=$1; centro=$2; ligandos=$3; salida=$4
# semilla y exhaustividad: por defecto 1 y 8
semilla=${5:-1}; exh=${6:-8}
mkdir -p "$salida"
# lado de la caja en ångström (cubre el bolsillo del ATP con margen)
lado=22
# se separan las tres coordenadas del centro
read -r cx cy cz <<< "$centro"

# función que acopla un ligando; xargs la llama con su ruta ($1)
acoplar() {
  # nombre del ligando sin directorio ni extensión
  nombre=$(basename "$1" .pdbqt)
  vina --receptor "$receptor" --ligand "$1" \
       --center_x "$cx" --center_y "$cy" --center_z "$cz" \
       --size_x "$lado" --size_y "$lado" --size_z "$lado" \
       --exhaustiveness "$exh" --seed "$semilla" \
       --cpu 1 --num_modes 5 \
       --out "$salida/${nombre}_poses.pdbqt" \
       > "$salida/${nombre}.log"
}
# xargs usa subprocesos: se exporta la función y sus variables
export -f acoplar
export receptor cx cy cz lado exh semilla salida

inicio=$(date +%s)
# -P 6: seis procesos a la vez; -n 1: un ligando por proceso
ls "$ligandos"/*.pdbqt | xargs -P 6 -n 1 bash -c 'acoplar "$0"'
fin=$(date +%s)
n=$(ls "$ligandos"/*.pdbqt | wc -l)
echo "ligandos: $n  tiempo_total: $((fin - inicio)) s"

# tabla: nombre y mejor energía (primera línea de resultados del registro), de menor a mayor
for log in "$salida"/*.log; do
  nombre=$(basename "$log" .log)
  # la línea siguiente a "-----+" es el modo 1 y su segunda columna, la energía
  energia=$(grep -A2 -- '-----+' "$log" | awk 'NR==2{print $2}')
  echo -e "$nombre\t$energia"
done | sort -t$'\t' -k2,2n > "$salida/ranking_sbvs.tsv"

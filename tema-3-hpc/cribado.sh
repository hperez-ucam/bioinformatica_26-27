#!/bin/bash
# cribado.sh - seis acoplamientos moleculares independientes, primero uno detrás
# de otro y después los seis a la vez, para medir qué se gana al repartir tareas.
#
#   ./cribado.sh serie       # uno detrás de otro, cada uno con un hilo
#   ./cribado.sh paralelo    # los seis a la vez, cada uno con un hilo
#
# Las seis tareas son el mismo acoplamiento con seis semillas distintas. Son
# independientes (ninguna necesita el resultado de otra), como seis ligandos de
# un cribado virtual. Cada tarea anota cuánto ha durado.
#
# set -e: el script se detiene si una orden falla; -u: también si se usa una
# variable sin valor; -o pipefail: una tubería falla si falla cualquiera de sus partes
set -euo pipefail

# modo de ejecución: primer argumento; si no se da, serie
modo=${1:-serie}
# directorio de este script, para encontrar datos/ aunque se lance desde otro
aqui=$(cd "$(dirname "$0")" && pwd)
# directorio de salida de este modo, dentro de salidas/
salida="$aqui/salidas/cribado_$modo"
mkdir -p "$salida"

# una tarea: un acoplamiento con un solo hilo. Dentro de la función, $1 es el
# valor con que se la llama: la semilla
acoplar() {
  # tiempo de reloj al empezar esta tarea, en segundos con decimales
  local t0
  t0=$(date +%s.%N)
  vina --receptor "$aqui/datos/1HVR_rec.pdbqt" \
       --ligand "$aqui/datos/1HVR_lig.pdbqt" \
       --center_x -9.192 --center_y 15.906 --center_z 27.946 \
       --size_x 22 --size_y 22 --size_z 22 \
       --exhaustiveness 1 --seed "$1" --cpu 1 \
       --out "$salida/poses_s$1.pdbqt" > "$salida/log_s$1.txt"
  # duración de la tarea: tiempo al terminar menos tiempo al empezar
  awk -v a="$t0" -v b="$(date +%s.%N)" 'BEGIN{printf "%.1f\n", b - a}' \
    > "$salida/tiempo_s$1.txt"
}

echo "$(date +%H:%M:%S) en marcha: seis acoplamientos, modo $modo."
# tiempo de reloj al empezar, en segundos
inicio=$(date +%s)

if [ "$modo" = "serie" ]; then
  # bucle normal: una semilla no empieza hasta que acaba la anterior
  for semilla in 1 2 3 4 5 6; do
    acoplar "$semilla"
    echo "  semilla $semilla terminada"
  done
else
  # el & ejecuta cada acoplamiento en segundo plano y el bucle sigue sin
  # esperar: empiezan los seis a la vez
  for semilla in 1 2 3 4 5 6; do
    acoplar "$semilla" &
  done
  # wait espera a que terminen todas las tareas enviadas a segundo plano
  wait
fi

fin=$(date +%s)
echo "modo=$modo tareas=6 tiempo_total=$((fin - inicio))s"
# mejor afinidad y duración de cada semilla
for semilla in 1 2 3 4 5 6; do
  afinidad=$(grep -A3 'mode |' "$salida/log_s$semilla.txt" | tail -1 | awk '{print $2}')
  echo "semilla $semilla: $afinidad kcal/mol en $(cat "$salida/tiempo_s$semilla.txt") s"
done
# suma de las seis duraciones y duración de la tarea más larga
cat "$salida"/tiempo_s*.txt | awk '{s += $1; if ($1 > m) m = $1}
  END{printf "suma=%.1f s  mas_larga=%.1f s\n", s, m}'

#!/bin/bash
# cribado.sh - seis dockings independientes, primero uno detrás de otro y luego
# todos a la vez, para ver qué gana el paralelismo de tareas.
#
#   ./cribado.sh serie       # uno detrás de otro, cada uno con un núcleo
#   ./cribado.sh paralelo    # los seis a la vez, cada uno con un núcleo
#
# Las seis tareas son el mismo docking con seis semillas distintas: son
# independientes (ninguna necesita el resultado de otra), como seis ligandos
# de un cribado. Es el caso más fácil de paralelizar.
#
set -euo pipefail

modo=${1:-serie}
aqui=$(cd "$(dirname "$0")" && pwd)
salida="$aqui/salidas/cribado_$modo"
mkdir -p "$salida"

# una tarea: un docking con la semilla que se le pase y un solo núcleo
docking() {
  vina --receptor "$aqui/datos/1HVR_rec.pdbqt" --ligand "$aqui/datos/1HVR_lig.pdbqt" \
       --center_x -9.192 --center_y 15.906 --center_z 27.946 \
       --size_x 22 --size_y 22 --size_z 22 \
       --exhaustiveness 1 --seed "$1" --cpu 1 \
       --out "$salida/poses_s$1.pdbqt" > "$salida/log_s$1.txt"
}

# reloj de pared al empezar, en segundos
inicio=$(date +%s)

if [ "$modo" = "serie" ]; then
  # bucle normal: la siguiente semilla no empieza hasta que acaba la anterior
  for semilla in 1 2 3 4 5 6; do
    docking "$semilla"
  done
else
  # el & manda cada docking al fondo y sigue sin esperar: arrancan los seis
  for semilla in 1 2 3 4 5 6; do
    docking "$semilla" &
  done
  # wait espera a que terminen todos los que se mandaron al fondo
  wait
fi

fin=$(date +%s)
echo "modo=$modo tareas=6 tiempo_total=$((fin - inicio))s"
# mejor afinidad de cada semilla, para comprobar que son seis resultados
for semilla in 1 2 3 4 5 6; do
  echo "semilla $semilla: $(grep -A3 'mode |' "$salida/log_s$semilla.txt" | tail -1 | awk '{print $2}') kcal/mol"
done

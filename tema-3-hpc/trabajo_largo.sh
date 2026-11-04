#!/bin/bash
# trabajo_largo.sh - un cálculo largo y real (un docking con Vina) que puede usar
# uno o varios núcleos, y que mide su tiempo real y su tiempo de CPU.
#
#   ./trabajo_largo.sh 1        # con un núcleo
#   ./trabajo_largo.sh 4        # con cuatro núcleos
#   ./trabajo_largo.sh 6 12     # seis núcleos y exhaustividad 12 (más largo)
#
# La exhaustividad es cuántas búsquedas independientes hace Vina; Vina reparte
# esas búsquedas entre los núcleos, así que con exhaustividad 4 no usa más de
# 4 núcleos aunque se le pidan 6. Por defecto vale 4: con un núcleo tarda algo
# más de minuto y medio en el servidor vacío.
#
# Siempre es el MISMO trabajo: misma proteína (1HVR), mismo ligando, misma caja,
# misma semilla y misma exhaustividad. Lo único que cambia es cuántos núcleos
# puede usar Vina. Así, si el tiempo cambia, es por el reparto, no por la tarea.
#
# Al final escribe una línea con:
#   real  = tiempo de reloj, el que pasa en la pared desde que empieza hasta que acaba
#   user  = tiempo de CPU sumado de todos los núcleos, haciendo el cálculo
#   sys   = tiempo de CPU que el sistema operativo gasta por el programa (ficheros, memoria)
#   (user+sys)/real = cuántos núcleos ha tenido trabajando de media
#
set -euo pipefail

# número de núcleos, primer argumento; si no se da, uno
nucleos=${1:-1}
# exhaustividad, segundo argumento; por defecto 4
exh=${2:-4}
# directorio de este script, para encontrar los ficheros aunque se llame desde otro sitio
aqui=$(cd "$(dirname "$0")" && pwd)
# un directorio de salida propio para cada ejecución
salida="$aqui/salidas"
mkdir -p "$salida"

# carga media del servidor justo antes de empezar (primer número de uptime)
carga_antes=$(cut -d" " -f1 /proc/loadavg)

# /usr/bin/time mide el programa y escribe sus tiempos en un fichero aparte;
# %e = real, %U = user, %S = sys, todos en segundos
/usr/bin/time -f "%e %U %S" -o "$salida/tiempos_c${nucleos}_e$exh.txt" \
  vina --receptor "$aqui/datos/1HVR_rec.pdbqt" --ligand "$aqui/datos/1HVR_lig.pdbqt" \
       --center_x -9.192 --center_y 15.906 --center_z 27.946 \
       --size_x 22 --size_y 22 --size_z 22 \
       --exhaustiveness "$exh" --seed 1 --cpu "$nucleos" \
       --out "$salida/poses_c${nucleos}_e$exh.pdbqt" > "$salida/log_c${nucleos}_e$exh.txt"

# carga media al terminar
carga_despues=$(cut -d" " -f1 /proc/loadavg)
# leer los tres tiempos que ha dejado /usr/bin/time
read -r real user sys < "$salida/tiempos_c${nucleos}_e$exh.txt"
# afinidad de la mejor pose: primera línea de resultados de Vina, segunda columna
afinidad=$(grep -A3 "mode |" "$salida/log_c${nucleos}_e$exh.txt" | tail -1 | awk '{print $2}')
# núcleos que han trabajado de media: (user + sys) / real
efectivos=$(awk -v u="$user" -v s="$sys" -v r="$real" 'BEGIN{printf "%.1f", (u+s)/r}')

# una sola línea, fácil de copiar a la tabla común
linea="usuario=$USER nucleos=$nucleos exh=$exh real=${real}s user=${user}s sys=${sys}s nucleos_efectivos=$efectivos afinidad=$afinidad carga_antes=$carga_antes carga_despues=$carga_despues"
echo "$linea"
# si el profesor ha creado el directorio común /tmp/tema3, dejar allí la línea
# (un fichero por usuario, para que nadie pise a nadie)
if [ -d /tmp/tema3 ] && [ -w /tmp/tema3 ]; then
  echo "$(date +%H:%M:%S) $linea" >> "/tmp/tema3/$USER.txt"
fi

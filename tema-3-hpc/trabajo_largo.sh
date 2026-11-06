#!/bin/bash
# trabajo_largo.sh - un cálculo largo y real (un acoplamiento molecular con Vina)
# que puede usar uno o varios hilos y que mide su tiempo real y su tiempo de CPU.
#
#   ./trabajo_largo.sh 1        # con un hilo
#   ./trabajo_largo.sh 4        # con cuatro hilos
#   ./trabajo_largo.sh 6 12     # seis hilos y exhaustividad 12 (más largo)
#
# La exhaustividad es el número de búsquedas independientes que hace Vina. Vina
# reparte esas búsquedas entre los hilos: con exhaustividad 4 no usa más de
# 4 hilos aunque se le pidan 6. Aquí vale 4 si no se indica otra cosa.
#
# Siempre es el MISMO trabajo: misma proteína (1HVR), mismo ligando, misma caja,
# misma semilla y misma exhaustividad. Solo cambia cuántos hilos puede usar
# Vina. Si el tiempo cambia, es por el reparto, no por la tarea.
#
# Al final escribe una línea con estos campos (los tiempos, en segundos):
#   hilos      hilos pedidos a Vina (--cpu)
#   exh        exhaustividad
#   real       tiempo de reloj, desde que empieza hasta que acaba
#   user       tiempo de CPU del cálculo, sumado de todos los hilos
#   sys        tiempo de CPU que el sistema operativo gasta por el programa
#   ocupados   (user + sys) / real: hilos que han trabajado de media
#   afinidad   energía de la mejor pose, en kcal/mol
#   carga      carga media del servidor al terminar (procesos que piden CPU)
#
# set -e: el script se detiene si una orden falla; -u: también si se usa una
# variable sin valor; -o pipefail: una tubería falla si falla cualquiera de sus partes
set -euo pipefail

# número de hilos, primer argumento; si no se da, uno
hilos=${1:-1}
# exhaustividad, segundo argumento; si no se da, 4
exh=${2:-4}
# directorio de este script, para encontrar los datos aunque se lance desde otro
aqui=$(cd "$(dirname "$0")" && pwd)
# directorio de salida: salidas/, junto al script
salida="$aqui/salidas"
mkdir -p "$salida"
# nombre común de los ficheros de esta ejecución: llevan los hilos y la exhaustividad
base="$salida/h${hilos}_e${exh}"

# aviso al empezar: Vina no escribe nada en pantalla hasta que termina
echo "$(date +%H:%M:%S) en marcha: $hilos hilo(s), exhaustividad $exh. No escribe nada hasta terminar."

# /usr/bin/time mide el programa y escribe sus tiempos en un fichero aparte;
# %e = real, %U = user, %S = sys, los tres en segundos
/usr/bin/time -f "%e %U %S" -o "$base.tiempos" \
  vina --receptor "$aqui/datos/1HVR_rec.pdbqt" \
       --ligand "$aqui/datos/1HVR_lig.pdbqt" \
       --center_x -9.192 --center_y 15.906 --center_z 27.946 \
       --size_x 22 --size_y 22 --size_z 22 \
       --exhaustiveness "$exh" --seed 1 --cpu "$hilos" \
       --out "$base.pdbqt" > "$base.log"

# leer los tres tiempos que ha dejado /usr/bin/time
read -r real user sys < "$base.tiempos"
# afinidad de la mejor pose: primera fila de la tabla de Vina, segunda columna
afinidad=$(grep -A3 "mode |" "$base.log" | tail -1 | awk '{print $2}')
# hilos que han trabajado de media: (user + sys) / real
ocupados=$(awk -v u="$user" -v s="$sys" -v r="$real" 'BEGIN{printf "%.1f", (u+s)/r}')
# carga media del último minuto: primer número de /proc/loadavg
carga=$(cut -d" " -f1 /proc/loadavg)

# una sola línea con todos los datos de la medida
linea="usuario=$USER hilos=$hilos exh=$exh real=$real user=$user sys=$sys"
linea="$linea ocupados=$ocupados afinidad=$afinidad carga=$carga"
echo "$linea"
# la línea se guarda también en el directorio personal, con fecha y hora:
# la piden las tareas del tema
echo "$(date '+%d/%m %H:%M:%S') $linea" >> "$HOME/tema3_lineas.txt"
echo "Línea guardada en ~/tema3_lineas.txt"
# si el profesor ha creado el directorio común /tmp/tema3, la línea se deja
# también allí (un fichero por usuario, para que nadie escriba en el de otro)
if [ -d /tmp/tema3 ] && [ -w /tmp/tema3 ]; then
  echo "$(date +%H:%M:%S) $linea" >> "/tmp/tema3/$USER.txt"
fi

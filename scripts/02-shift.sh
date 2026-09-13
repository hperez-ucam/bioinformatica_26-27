#!/bin/bash
# shift descarta $1 y corre el resto una posición a la izquierda.
# Uso: ./02-shift.sh uno dos tres cuatro

echo "Script: $0"
echo "Al empezar, $# parámetros: $*"
echo "Ahora el primero es $1"
shift
echo "Tras un shift, el primero es $1 (era el segundo)"
shift
echo "Tras otro shift, el primero es $1 (era el tercero)"
echo "Y el segundo es $2 (era el cuarto)"

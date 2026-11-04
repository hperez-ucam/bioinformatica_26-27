#!/bin/bash
# muchos_ficheros.sh - los mismos datos guardados en un fichero o en 20 000,
# para ver que el disco y el sistema de ficheros también cuestan tiempo.
#
#   ./muchos_ficheros.sh
#
# Es el problema de un cribado mal organizado: un fichero pequeño por molécula.
# El cálculo es idéntico; lo que cambia es cuántas veces hay que abrir, escribir
# y cerrar un fichero.
#
set -euo pipefail

n=20000
aqui=$(cd "$(dirname "$0")" && pwd)
dir="$aqui/salidas/ficheros"
# empezar siempre de cero
rm -rf "$dir"
mkdir -p "$dir/muchos"

# 1) las n líneas en UN fichero; time mide el bloque entre llaves
echo "== escribir $n líneas en un solo fichero"
time { for i in $(seq "$n"); do echo "molecula_$i -7.5"; done > "$dir/uno.txt"; }

# 2) las mismas n líneas, una por fichero
echo "== escribir $n ficheros de una línea"
time { for i in $(seq "$n"); do echo "molecula_$i -7.5" > "$dir/muchos/m_$i.txt"; done; }

# 3) leerlo todo y contar líneas, en los dos casos
echo "== leer el fichero único"
time wc -l < "$dir/uno.txt"
echo "== leer los $n ficheros"
time cat "$dir"/muchos/m_*.txt | wc -l

# espacio que ocupa cada versión: cada fichero, por pequeño que sea, gasta un bloque
du -sh "$dir/uno.txt" "$dir/muchos"
# borrar los ficheros de prueba: 20 000 ficheros por alumno llenan el disco común
rm -rf "$dir"

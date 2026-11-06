#!/bin/bash
# muchos_ficheros.sh - los mismos datos guardados en un fichero o en 20 000,
# para medir que el disco y el sistema de ficheros también consumen tiempo.
#
#   ./muchos_ficheros.sh
#
# Es el problema de un cribado que guarda un fichero pequeño por molécula. El
# cálculo es idéntico; cambia cuántas veces se abre, se escribe y se cierra un
# fichero. Los ficheros de prueba se borran siempre al terminar.
#
# set -e: el script se detiene si una orden falla; -u: también si se usa una
# variable sin valor; -o pipefail: una tubería falla si falla cualquiera de sus partes
set -euo pipefail

# número de líneas (y de ficheros pequeños)
n=20000
# directorio de este script
aqui=$(cd "$(dirname "$0")" && pwd)
# directorio de trabajo de la prueba
dir="$aqui/salidas/ficheros"
# empezar siempre de cero
rm -rf "$dir"
mkdir -p "$dir/muchos"
# al salir, por la razón que sea (un error o Ctrl+C), borrar los ficheros de
# prueba: 20 000 ficheros por persona llenan el disco común
trap 'rm -rf "$dir"' EXIT

# 1) las n líneas en UN fichero; time mide el bloque entre llaves
echo "== escribir $n líneas en un solo fichero"
time { for i in $(seq "$n"); do echo "molecula_$i -7.5"; done > "$dir/uno.txt"; }

# 2) las mismas n líneas, una por fichero
echo "== escribir $n ficheros de una línea"
time { for i in $(seq "$n"); do echo "molecula_$i -7.5" > "$dir/muchos/m_$i.txt"; done; }

# 3) leerlo todo y contar las líneas, en los dos casos
echo "== leer el fichero único"
time wc -l < "$dir/uno.txt"
echo "== leer los $n ficheros"
# se entra en el directorio para que los 20 000 nombres sean cortos: con la
# ruta completa delante, la lista de nombres se acerca al tamaño máximo que
# admite una línea de órdenes (getconf ARG_MAX)
time ( cd "$dir/muchos" && cat m_*.txt | wc -l )

# 4) espacio en disco: cada fichero, por pequeño que sea, ocupa al menos un bloque
echo "== espacio en disco"
du -sh "$dir/uno.txt" "$dir/muchos"

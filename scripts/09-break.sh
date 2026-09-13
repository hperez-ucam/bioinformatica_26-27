#!/bin/bash
# break sale del bucle. Aquí para al encontrar f.
# Uso: ./09-break.sh a b f c

while [ "$#" -gt 0 ]; do
    if [ "$1" = "f" ]; then
        break
    fi
    echo "Parámetro: $1"
    shift
done

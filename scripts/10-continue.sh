#!/bin/bash
# continue salta a la siguiente vuelta. Aquí ignora f.
# Uso: ./10-continue.sh a b f c

while [ "$#" -gt 0 ]; do
    if [ "$1" = "f" ]; then
        shift
        continue
    fi
    echo "Parámetro: $1"
    shift
done

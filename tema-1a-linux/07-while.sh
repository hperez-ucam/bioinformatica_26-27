#!/bin/bash
# Mientras quede algún parámetro, lo muestra y hace shift.
# Uso: ./07-while.sh a dos tres

while [ -n "$1" ]; do
    echo "Parámetro: $1"
    shift
done

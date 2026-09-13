#!/bin/bash
# Ejemplo mediano: comprueba el fichero, cuenta líneas ATOM y avisa si no hay.
# Uso: ./00b-contar-atom.sh ejemplo.pdb

fichero=${1:-ejemplo.pdb}

if [ ! -f "$fichero" ]; then
    echo "No existe $fichero"
    exit 1
fi

n=$(grep -c '^ATOM' "$fichero" || true)
echo "ATOM en $fichero: $n"

if [ "$n" -eq 0 ]; then
    echo "Aviso: ninguna línea ATOM"
fi

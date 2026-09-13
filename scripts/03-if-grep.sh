#!/bin/bash
# Si el fichero contiene ATOM, lo dice. Si no, también.
# Uso: ./03-if-grep.sh ejemplo.pdb

fichero=${1:-ejemplo.pdb}

if grep -q ATOM "$fichero"; then
    echo "Encontrada la palabra ATOM en $fichero"
else
    echo "No encontrada la palabra ATOM en $fichero"
fi

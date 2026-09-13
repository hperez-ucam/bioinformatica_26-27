#!/bin/bash
# Tiene dos errores a propósito: ejecutadlo, leed el mensaje y corregidlo.
# Uso: ./12-errores.sh ejemplo.pdb   y después   ./12-errores.sh

fichero=$1

if [ -f "$fichero"]; then
    echo "Existe $fichero"
fi

if [ $fichero = "ejemplo.pdb" ]; then
    echo "Es el PDB de juguete"
fi

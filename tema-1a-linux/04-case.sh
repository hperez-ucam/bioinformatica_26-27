#!/bin/bash
# Clasifica el primer argumento por patrón (glob, no regex).
# Uso: ./04-case.sh pdb
#      ./04-case.sh fasta

case $1 in
    pdb) echo "Estructura 3D" ;;
    fasta|fa) echo "Secuencia" ;;
    *) echo "Extensión no reconocida" ;;
esac

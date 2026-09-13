#!/bin/bash
# Ejemplo complejo: recorre un lote, filtra por extensión y resume ATOM.
# Uso: ./00c-lote-pdb.sh ejemplo.pdb otro.pdb notas.txt

if [ "$#" -eq 0 ]; then
    echo "Uso: $0 fichero1 [fichero2 …]"
    exit 1
fi

for f in "$@"; do
    case "$f" in
        *.pdb|*.PDB)
            if [ ! -f "$f" ]; then
                echo "Falta: $f"
                continue
            fi
            n=$(grep -c '^ATOM' "$f" || true)
            echo "$f → $n ATOM"
            ;;
        *)
            echo "Ignorado (no PDB): $f"
            ;;
    esac
done

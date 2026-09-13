#!/bin/bash
# seq inicio paso fin, sustituido con $().
# Uso: ./06-for-seq.sh
# Valores: 0 5 10 15 20 25

for i in $(seq 0 5 25); do
    echo "i = $i"
done

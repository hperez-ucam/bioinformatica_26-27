#!/bin/bash
# until = mientras la expresión sea falsa.
# Uso: ./08-until.sh a dos tres

until [ -z "$1" ]; do
    echo "$1"
    shift
done

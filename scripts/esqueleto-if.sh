#!/bin/bash
# Esqueleto de if / elif / else / fi (plantilla de aula ejecutable).
# Uso: ./esqueleto-if.sh a   |   ./esqueleto-if.sh b   |   ./esqueleto-if.sh x

if [ "$1" = "a" ]; then
    echo "rama then"
elif [ "$1" = "b" ]; then
    echo "rama elif"
else
    echo "rama else"
fi

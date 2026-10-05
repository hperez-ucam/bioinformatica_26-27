#!/bin/bash
# correr-vina.sh - tres ejecuciones de AutoDock Vina sobre el mismo sistema,
# para ver que un metodo estocastico no devuelve siempre el mismo resultado.
#
#   ./correr-vina.sh 1HVR 2.5 3.1 10.4
#
# Las tres ejecuciones cambian un solo parametro cada vez:
#   A  exhaustiveness 8   semilla 1    (configuracion por defecto)
#   B  exhaustiveness 32  semilla 1    (mas busqueda, misma semilla)
#   C  exhaustiveness 8   semilla 777  (misma busqueda, otra semilla)
#
set -euo pipefail

PDB="${1:?uso: ./correr-vina.sh CODIGO_PDB centro_x centro_y centro_z}"
CX="${2:?falta centro_x}"; CY="${3:?falta centro_y}"; CZ="${4:?falta centro_z}"
TAM="${5:-22}"            # lado de la caja en angstrom

comun=(--receptor "${PDB}_rec.pdbqt" --ligand "${PDB}_lig.pdbqt"
       --center_x "$CX" --center_y "$CY" --center_z "$CZ"
       --size_x "$TAM" --size_y "$TAM" --size_z "$TAM"
       --num_modes 9)

ejecutar () {  # nombre exhaustividad semilla
    local nombre="$1" exh="$2" sem="$3"
    echo "== ${nombre}: exhaustiveness=${exh}  seed=${sem} =="
    /usr/bin/time -f "   tiempo: %e s" \
        vina "${comun[@]}" --exhaustiveness "$exh" --seed "$sem" \
             --out "salida_${nombre}.pdbqt" > "log_${nombre}.txt" 2>&1 || true
    # La tabla de resultados son las lineas que empiezan por un numero de modo.
    echo "   modo  afinidad(kcal/mol)  rmsd_l.b.  rmsd_u.b."
    grep -E '^\s+[0-9]+\s+-?[0-9]' "log_${nombre}.txt" | head -5 | sed 's/^/   /'
}

ejecutar A 8 1
ejecutar B 32 1
ejecutar C 8 777

echo
echo "== comparacion =="
for n in A B C; do
    mejor=$(grep -E '^\s+1\s+-?[0-9]' "log_${n}.txt" | awk '{print $2}')
    printf "   %s  mejor afinidad: %s kcal/mol\n" "$n" "${mejor:-sin resultado}"
done
echo
echo "   Si A y C difieren, es la semilla: la busqueda es estocastica."
echo "   Si B mejora a A, es que con exhaustiveness 8 la busqueda se quedaba"
echo "   corta. Lo que no se puede concluir de una sola ejecucion es nada."

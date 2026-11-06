#!/bin/bash
# tabla_clase.sh - resumen de las medidas de toda la clase: las líneas que
# trabajo_largo.sh deja en el directorio común /tmp/tema3.
#
#   ./tabla_clase.sh       # una fila por ronda (mismos hilos y exhaustividad)
#   ./tabla_clase.sh -t    # todas las líneas, una por medida
#
# set -u: el script se detiene si se usa una variable sin valor
set -u

# directorio común de la clase; lo crea el profesor
dir=/tmp/tema3
if ! ls "$dir"/*.txt > /dev/null 2>&1; then
  echo "Todavía no hay medidas en $dir."
  exit 0
fi

# ./tabla_clase.sh -t : todas las líneas, ordenadas por hora
if [ "${1:-}" = "-t" ]; then
  cat "$dir"/*.txt | sort
  exit 0
fi

echo "hilos  exh  trabajos  real_min  real_mediana  real_max  user_medio"
# de cada línea se sacan hilos, exhaustividad, real y user (campos clave=valor);
# se ordenan por ronda y por tiempo real
cat "$dir"/*.txt | awk '{
    for (i = 1; i <= NF; i++) { split($i, par, "="); dato[par[1]] = par[2] }
    print dato["hilos"], dato["exh"], dato["real"], dato["user"]
  }' | sort -k1,1n -k2,2n -k3,3n | awk '
  # al cambiar de ronda se escribe la fila de la anterior: mínimo, mediana y máximo
  function cierra() {
    if (n) printf "%5s  %3s  %8d  %8.1f  %12.1f  %8.1f  %10.1f\n", h, e, n, r[1], r[int((n + 1) / 2)], r[n], suma / n
  }
  { if ($1 != h || $2 != e) { cierra(); h = $1; e = $2; n = 0; suma = 0 }
    n++; r[n] = $3; suma += $4 }
  END { cierra() }'

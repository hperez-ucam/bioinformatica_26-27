#!/bin/bash
# cola.sh - cola mínima de trabajos para el servidor de clase, hecha con flock.
#
#   ./cola.sh ORDEN [ARGUMENTOS]   envía un trabajo: espera turno, lo ejecuta y avisa
#   ./cola.sh -l                   trabajos pendientes y en ejecución (como squeue)
#   ./cola.sh -a                   trabajos terminados, con sus tiempos (como sacct)
#   ./cola.sh -p N                 solo el profesor: prepara la cola con N turnos
#
# Como mucho N trabajos se ejecutan a la vez; los demás esperan en orden de
# llegada. Estados, con las abreviaturas de Slurm:
#   PD  pendiente (pending): el trabajo espera un turno libre
#   R   en ejecución (running)
#   CD  completado (completed); F si la orden termina con error (failed)
# Diferencias con Slurm: aquí se reparten turnos, no núcleos ni memoria; no hay
# límite de tiempo ni prioridades; y el trabajo se ejecuta en primer plano, como
# con srun: depende de la terminal desde la que se envía. Con sbatch, Slurm
# devuelve el control y guarda la salida del trabajo en un fichero.
#
# set -u: el script se detiene si se usa una variable sin valor; -o pipefail:
# una tubería falla si falla cualquiera de sus partes
set -uo pipefail

# directorio común de la cola; se puede cambiar con la variable COLA_DIR
cola=${COLA_DIR:-/tmp/tema3/cola}

# ---------------------------------------------------------------------------
# ./cola.sh -p N : preparar la cola (lo ejecuta el profesor una vez)
# ---------------------------------------------------------------------------
if [ "${1:-}" = "-p" ]; then
  # número de turnos: segundo argumento; por defecto 6, uno por núcleo físico
  n=${2:-6}
  # N tiene que ser un número entero mayor que cero
  case "$n" in
    ''|*[!0-9]*|0) echo "uso: ./cola.sh -p N   (N = número de turnos, 1 o más)"; exit 1 ;;
  esac
  # directorio de la cola y, dentro, el de las fichas de los trabajos
  mkdir -p "$cola/fichas"
  # /tmp/tema3 es el directorio común de la clase: todos pueden escribir en él y
  # cada uno solo borra lo suyo (permisos 1777, como /tmp)
  if [ "$cola" = /tmp/tema3/cola ]; then chmod 1777 /tmp/tema3; fi
  # en la cola solo escribe quien la prepara; los demás leen
  chmod 755 "$cola"
  # en fichas escribe toda la clase, y cada uno solo borra lo suyo (como en /tmp)
  chmod 1777 "$cola/fichas"
  # un fichero vacío por turno, turno_1 ... turno_N: sirven solo para poner cerrojos
  for i in $(seq 1 "$n"); do
    [ -e "$cola/turno_$i" ] || : > "$cola/turno_$i"
  done
  # si antes había más turnos, se quitan los que sobran
  for f in "$cola"/turno_*; do
    if [ "${f##*_}" -gt "$n" ]; then rm -f "$f"; fi
  done
  # se borran las fichas de trabajos que ya no existen (por ejemplo, tras un kill -9)
  for f in "$cola"/fichas/pd_* "$cola"/fichas/r_*; do
    if [ -e "$f" ] && [ ! -d "/proc/${f##*_}" ]; then rm -f "$f"; fi
  done
  # contador de trabajos: empieza en 0 la primera vez
  [ -e "$cola/contador" ] || echo 0 > "$cola/contador"
  # el contador lo actualiza cualquier usuario al enviar un trabajo
  chmod 666 "$cola/contador"
  echo "cola preparada en $cola con $n turnos"
  exit 0
fi

# sin la cola preparada no se puede seguir
if [ ! -d "$cola/fichas" ] || [ ! -e "$cola/turno_1" ]; then
  echo "La cola $cola no está preparada: avisad al profesor (./cola.sh -p N)."
  exit 1
fi

# ---------------------------------------------------------------------------
# ./cola.sh -l : trabajos pendientes y en ejecución de toda la clase
# ---------------------------------------------------------------------------
if [ "${1:-}" = "-l" ]; then
  echo "TRABAJO  USUARIO       ESTADO  DESDE     ORDEN"
  # primero las fichas de los pendientes (pd_) y después las de los que se ejecutan (r_)
  for f in "$cola"/fichas/pd_* "$cola"/fichas/r_*; do
    # el nombre de la ficha termina en el número de proceso (PID) del trabajo;
    # solo se enseña si ese proceso sigue existiendo
    if [ -d "/proc/${f##*_}" ]; then cat "$f"; fi
  done 2>/dev/null
  exit 0
fi

# ---------------------------------------------------------------------------
# ./cola.sh -a : trabajos ya terminados, con sus tiempos
# ---------------------------------------------------------------------------
if [ "${1:-}" = "-a" ]; then
  echo "TRABAJO  USUARIO       ESTADO  ESPERA  EJECUCION  TOTAL  ORDEN"
  # cada usuario apunta sus trabajos terminados en su fichero hecho_USUARIO
  cat "$cola"/fichas/hecho_* 2>/dev/null | sort -n
  exit 0
fi

# sin argumentos no hay ninguna orden que enviar
if [ $# -eq 0 ]; then
  echo "uso: ./cola.sh ORDEN [ARGUMENTOS]  |  ./cola.sh -l  |  ./cola.sh -a"
  exit 1
fi

# ---------------------------------------------------------------------------
# enviar un trabajo
# ---------------------------------------------------------------------------
# Ctrl+Z no detiene un trabajo de la cola: detenido, bloquearía a los que esperan detrás
trap '' TSTP
# la orden, tal como se escribió detrás de ./cola.sh, en una sola cadena para mostrarla
orden="$*"

# 1) número de trabajo: leer el contador, sumarle uno y guardarlo sin que otro trabajo lo cambie a la vez
# se abre el contador para lectura en el descriptor 8
exec 8< "$cola/contador"
# cerrojo exclusivo: si otro trabajo está cogiendo número, se espera aquí
flock 8
# el último número dado, más uno
id=$(( $(cat "$cola/contador") + 1 ))
# se guarda para el siguiente trabajo
echo "$id" > "$cola/contador"
# se suelta el cerrojo
flock -u 8
# y se cierra el descriptor
exec 8<&-

# 2) ficha del trabajo: su nombre lleva el estado, el número, el usuario y el PID ($$)
ficha=$(printf '%s/fichas/pd_%05d_%s_%s' "$cola" "$id" "$USER" "$$")
# al terminar, por la razón que sea (también con Ctrl+C), se borra la ficha
trap 'rm -f "$ficha"' EXIT
# hora de envío: en segundos para hacer cuentas
enviado=$(date +%s)
# contenido de la ficha: la línea que enseña ./cola.sh -l
printf '%7d  %-12s  %-6s  %s  %s\n' "$id" "$USER" PD "$(date +%H:%M:%S)" "$orden" > "$ficha"
# aviso en pantalla: estado PD
echo "[$(date +%H:%M:%S)] trabajo $id  PD  pendiente: $orden"

# ¿hay algún trabajo pendiente más antiguo que este y que siga vivo?
hay_alguien_delante() {
  local f
  # las fichas pd_ salen ordenadas por número de trabajo, del más antiguo al más nuevo
  for f in "$cola"/fichas/pd_*; do
    # al llegar a la ficha propia ya no queda nadie delante
    if [ "$f" = "$ficha" ]; then return 1; fi
    # una ficha más antigua cuenta si su proceso (el PID del final del nombre) sigue vivo
    if [ -d "/proc/${f##*_}" ]; then return 0; fi
  done
  return 1
}

# 3) esperar un turno: cada turno es un fichero y el turno es de quien tiene su cerrojo
turno=0
until [ "$turno" -gt 0 ]; do
  # solo el pendiente más antiguo intenta coger turno: orden de llegada
  if ! hay_alguien_delante; then
    for cerrojo in "$cola"/turno_*; do
      # se abre el fichero del turno en el descriptor 9
      exec 9< "$cerrojo"
      # flock -n intenta coger el cerrojo sin esperar; si lo tiene otro trabajo, falla
      if flock -n 9; then
        # turno conseguido: su número es el final del nombre del fichero
        turno=${cerrojo##*_}
        break
      fi
      # turno ocupado: se cierra el descriptor y se prueba el siguiente
      exec 9<&-
    done
  fi
  # sin turno todavía: esperar un segundo y volver a mirar
  if [ "$turno" -eq 0 ]; then sleep 1; fi
done

# 4) en ejecución: el cerrojo del turno sigue cogido mientras viva este script
empieza=$(date +%s)
# la ficha cambia de nombre (de pd_ a r_) y de contenido
nueva=$(printf '%s/fichas/r_%05d_%s_%s' "$cola" "$id" "$USER" "$$")
mv "$ficha" "$nueva"
ficha=$nueva
printf '%7d  %-12s  %-6s  %s  %s\n' "$id" "$USER" R "$(date +%H:%M:%S)" "$orden" > "$ficha"
# aviso en pantalla: estado R, con el turno y lo que se ha esperado
echo "[$(date +%H:%M:%S)] trabajo $id  R   en ejecución en el turno $turno, tras esperar $((empieza - enviado)) s"

# la orden se ejecuta tal como se escribió; se guarda su código de salida
"$@"
codigo=$?
# código 0 es completado (CD); cualquier otro, fallo (F)
if [ "$codigo" -eq 0 ]; then estado=CD; else estado=F; fi

# 5) terminado: tiempos de espera, de ejecución y total
fin=$(date +%s)
printf '[%s] trabajo %d  %-2s  terminado: espera %d s + ejecución %d s = %d s\n' "$(date +%H:%M:%S)" \
  "$id" "$estado" "$((empieza - enviado))" "$((fin - empieza))" "$((fin - enviado))"
# una línea en el registro propio, que es lo que enseña ./cola.sh -a
printf '%7d  %-12s  %-6s  %6d  %9d  %5d  %s\n' "$id" "$USER" "$estado" \
  "$((empieza - enviado))" "$((fin - empieza))" "$((fin - enviado))" "$orden" >> "$cola/fichas/hecho_$USER"
# al salir, el trap borra la ficha y el sistema suelta el cerrojo del turno
exit "$codigo"

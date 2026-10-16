#!/bin/bash
# md-dipeptido.sh - una dinámica molecular completa, de principio a fin, sobre
# el sistema más pequeno que tiene sentido simular: el dipeptido de alanina.
#
# 22 átomos de soluto y unas 650 moléculas de agua. Son unos 2 000 átomos en
# total, así que corre en CPU en un ordenador normal: la producción de 1 ns
# tarda del orden de media hora. No hace falta GPU ni cola de trabajos.
#
# Es el mismo flujo de siete etapas de la figura del 2.4, ejecutado paso a paso;
# las cabeceras van numeradas como las cajas de esa figura.
#
#   ./md-dipeptido.sh
#
# Se ejecuta desde la carpeta dinámica/, donde están minim.mdp y md.mdp; el
# PDB de partida lo escribe ../minimización/dipeptido.py, que debe estar al lado.
#
set -euo pipefail

# -- 0. comprobar que GROMACS esta disponible --------------------------------
# si no esta, el guion para aquí en vez de fallar a mitad
command -v gmx >/dev/null || { echo "No encuentro gmx. En el servidor de la asignatura deberia estar instalado: avisa al profesor."; exit 1; }

# -- 1. la molécula ----------------------------------------------------------
# la geometría de partida la escribe dipeptido.py, ya minimizada con MMFF94 y
# con los nombres de residuo y de átomo que espera el campo de fuerza. Asi la
# dinámica arranca exactamente de la misma molécula cuyo mapa calculamos
python3 ../minimizacion/dipeptido.py --pdb

# -- 2. topología: que campo de fuerza describe a cada átomo -----------------
# pdb2gmx lee el PDB y escribe la topología: a cada átomo le asigna su tipo y
# sus parámetros. El campo de fuerza y el modelo de agua se piden por nombre;
# -ignh descarta los hidrógenos del PDB (no los hay) y los reconstruye
gmx pdb2gmx -f dipeptido.pdb -o procesado.gro -p topol.top \
            -ff amber99sb-ildn -water tip3p -ignh

# -- 3. la caja y el disolvente --------------------------------------------------------------
# caja cubica con 1 nm de margen entre el soluto y el borde: menos margen y la
# molécula interaccionaria con su propia imagen periódica
gmx editconf -f procesado.gro -o caja.gro -c -d 1.0 -bt cubic

# -- 3b. el disolvente --------------------------------------------------------
# rellena la caja de agua; la topología se actualiza sola con el número de
# moléculas anadidas
gmx solvate -cp caja.gro -cs spc216.gro -o solvatado.gro -p topol.top

# -- 4. minimización ---------------------------------------------------------
# grompp junta estructura, topología y parámetros en un único fichero binario
gmx grompp -f minim.mdp -c solvatado.gro -p topol.top -o em.tpr
# mdrun ejecuta: aquí baja la energía hasta quitar los choques del solvatado
gmx mdrun -v -deffnm em

# -- 5. equilibrado a volumen constante ----------------------------------------------------------
# 50 ps a volumen constante para que el agua se ordene alrededor del soluto
# no hay restricciones de posición: con un soluto de 22 átomos y la caja ya a la
# densidad del agua el error es pequeño; en una proteína se añaden antes de equilibrar
cat > nvt-dipeptido.mdp <<'MDP'
integrator  = md          ; leap-frog
dt          = 0.002       ; 2 fs, posible porque los enlaces con H van restringidos
nsteps      = 25000       ; 25 000 x 2 fs = 50 ps
tcoupl      = V-rescale   ; termostato
tc-grps     = System
tau-t       = 0.1         ; ps, tiempo de acoplamiento
ref-t       = 300         ; K
constraints = h-bonds     ; LINCS sobre los enlaces con H
coulombtype = PME         ; electrostatica de largo alcance
rvdw        = 1.0         ; nm
rcoulomb    = 1.0         ; nm
gen-vel     = yes         ; velocidades iniciales de Maxwell-Boltzmann...
gen-temp    = 300         ; ...a esta temperatura
MDP
gmx grompp -f nvt-dipeptido.mdp -c em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt

# -- 5b. equilibrado a presión constante ------------------------------------
# 50 ps con barostato: la caja cambia de volumen hasta que el agua alcanza su
# densidad real. Es el paso que le falta al NVT para que la producción sea
# a T y P constantes, como el md.mdp de clase
cat > npt-dipeptido.mdp <<'MDP'
integrator  = md
dt          = 0.002
nsteps      = 25000
tcoupl      = V-rescale
tc-grps     = System
tau-t       = 0.1
ref-t       = 300
pcoupl      = C-rescale
tau-p       = 2.0
ref-p       = 1.0
compressibility = 4.5e-5
constraints = h-bonds
coulombtype = PME
rvdw        = 1.0
rcoulomb    = 1.0
continuation = yes
MDP
gmx grompp -f npt-dipeptido.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt

# -- 6. producción -----------------------------------------------------------
# 100 ps son 50 000 pasos de 2 fs, a T y P constantes como en el md.mdp de
# clase; es lo único que se analiza despues.
# El md.mdp de esta carpeta es el de una producción de 10 ns, para leerlo en
# clase; aquí se escribe uno más corto con otro nombre para no pisarlo
cat > md-dipeptido.mdp <<'MDP'
integrator  = md
dt          = 0.002
nsteps      = 500000
nstxout-compressed = 500   ; guarda coordenadas cada 1 ps
tcoupl      = V-rescale
tc-grps     = System
tau-t       = 0.1
ref-t       = 300
pcoupl      = C-rescale
tau-p       = 2.0
ref-p       = 1.0
compressibility = 4.5e-5
constraints = h-bonds
coulombtype = PME
rvdw        = 1.0
rcoulomb    = 1.0
continuation = yes        ; sigue del NPT: no regenera velocidades
MDP
gmx grompp -f md-dipeptido.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -v -deffnm md

# -- 7. análisis -------------------------------------------------------------
# saca los ángulos phi y psi a lo largo de la trayectoria: los mismos dos
# ángulos cuya superficie calculamos con dipeptido.py, ahora recorridos por
# una simulacion de verdad
gmx rama -s md.tpr -f md.xtc -o rama.xvg

echo
echo "Listo. rama.xvg tiene los pares (phi, psi) de toda la trayectoria."
echo "Compáralos con el mapa EN VACÍO de dipeptido.py: en agua el puente de"
echo "hidrógeno interno de C7eq se deshace y aparecen las regiones beta/PPII y"
echo "alfa del mapa de Ramachandran. El disolvente cambia el paisaje."

#!/bin/bash
# md-dipeptido.sh - una dinamica molecular completa, de principio a fin, sobre
# el sistema mas pequeno que tiene sentido simular: el dipeptido de alanina.
#
# 22 atomos de soluto y unas 650 moleculas de agua. Son unos 2000 atomos en
# total, asi que corre en CPU en un ordenador normal: la produccion de 100 ps
# tarda unos minutos. No hace falta GPU ni cola de trabajos.
#
# Es el mismo flujo de siete etapas de la diapositiva 2.4, ejecutado de verdad.
#
#   ./md-dipeptido.sh
#
set -euo pipefail

# -- 0. comprobar que GROMACS esta disponible --------------------------------
# si no esta, el guion para aqui en vez de fallar a mitad
command -v gmx >/dev/null || { echo "No encuentro gmx. Carga el modulo de GROMACS."; exit 1; }

# -- 1. la molecula ----------------------------------------------------------
# el dipeptido se escribe a mano en formato PDB: son pocos atomos y asi se ve
# que un PDB no es mas que una tabla de coordenadas
cat > dipeptido.pdb <<'PDB'
ATOM      1  CH3 ACE A   1       2.000   1.000   0.000  1.00  0.00           C
ATOM      2  C   ACE A   1       1.500   2.400   0.000  1.00  0.00           C
ATOM      3  O   ACE A   1       2.250   3.370   0.000  1.00  0.00           O
ATOM      4  N   ALA A   2       0.170   2.520   0.000  1.00  0.00           N
ATOM      5  CA  ALA A   2      -0.500   3.800   0.000  1.00  0.00           C
ATOM      6  CB  ALA A   2       0.000   4.600   1.200  1.00  0.00           C
ATOM      7  C   ALA A   2      -2.000   3.600   0.000  1.00  0.00           C
ATOM      8  O   ALA A   2      -2.750   4.570   0.000  1.00  0.00           O
ATOM      9  N   NME A   3      -2.500   2.350   0.000  1.00  0.00           N
ATOM     10  CH3 NME A   3      -3.950   2.100   0.000  1.00  0.00           C
TER
END
PDB

# -- 2. topologia: que campo de fuerza describe a cada atomo -----------------
# pdb2gmx lee el PDB y escribe la topologia. El 6 elige AMBER99SB-ILDN y el 1,
# el modelo de agua TIP3P; se pasan por la entrada estandar para no tener que
# contestar a mano
printf '6\n1\n' | gmx pdb2gmx -f dipeptido.pdb -o procesado.gro -p topol.top -ignh

# -- 3. la caja --------------------------------------------------------------
# caja cubica con 1 nm de margen entre el soluto y el borde: menos margen y la
# molecula interaccionaria con su propia imagen periodica
gmx editconf -f procesado.gro -o caja.gro -c -d 1.0 -bt cubic

# -- 4. el disolvente --------------------------------------------------------
# rellena la caja de agua; la topologia se actualiza sola con el numero de
# moleculas anadidas
gmx solvate -cp caja.gro -cs spc216.gro -o solvatado.gro -p topol.top

# -- 5. minimizacion ---------------------------------------------------------
# grompp junta estructura, topologia y parametros en un unico fichero binario
gmx grompp -f ../gromacs/minim.mdp -c solvatado.gro -p topol.top -o em.tpr -maxwarn 2
# mdrun ejecuta: aqui baja la energia hasta quitar los choques del solvatado
gmx mdrun -v -deffnm em

# -- 6. equilibrado ----------------------------------------------------------
# 50 ps a volumen constante para que el agua se ordene alrededor del soluto
cat > nvt.mdp <<'MDP'
integrator  = md
dt          = 0.002
nsteps      = 25000
tcoupl      = V-rescale
tc-grps     = System
tau-t       = 0.1
ref-t       = 300
constraints = h-bonds
coulombtype = PME
rvdw        = 1.0
rcoulomb    = 1.0
gen-vel     = yes
gen-temp    = 300
MDP
gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -deffnm nvt

# -- 7. produccion -----------------------------------------------------------
# 100 ps son 50 000 pasos de 2 fs; es lo unico que se analiza despues
cat > md.mdp <<'MDP'
integrator  = md
dt          = 0.002
nsteps      = 50000
nstxout-compressed = 500
tcoupl      = V-rescale
tc-grps     = System
tau-t       = 0.1
ref-t       = 300
constraints = h-bonds
coulombtype = PME
rvdw        = 1.0
rcoulomb    = 1.0
continuation = yes
MDP
gmx grompp -f md.mdp -c nvt.gro -p topol.top -o md.tpr -maxwarn 2
gmx mdrun -v -deffnm md

# -- 8. analisis -------------------------------------------------------------
# saca los angulos phi y psi a lo largo de la trayectoria: los mismos dos
# angulos cuya superficie calculamos con dipeptido.py, ahora recorridos por
# una simulacion de verdad
gmx rama -s md.tpr -f md.xtc -o rama.xvg

echo
echo "Listo. rama.xvg tiene los pares (phi, psi) de toda la trayectoria."
echo "Comparalos con el mapa que da dipeptido.py: la simulacion deberia"
echo "pasar la mayor parte del tiempo en los mismos minimos."

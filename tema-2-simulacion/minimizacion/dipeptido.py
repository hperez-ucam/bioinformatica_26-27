#!/usr/bin/env python3
"""El dipeptido de alanina: de una curva a una superficie.

En butano.py la molecula tenia UN grado de libertad y la energia cabia en una
curva. Aqui hay DOS, los angulos phi y psi del esqueleto peptidico, y la energia
ya es una superficie: el mapa de Ramachandran. Una proteina de verdad tiene
3N-6 grados de libertad y no se puede dibujar; este es el ultimo caso que si.

La molecula es ACE-ALA-NME, 22 atomos. Es el sistema de prueba clasico de la
dinamica molecular de biomoleculas.

Uso:
    python3 dipeptido.py              # mapa + minimizaciones + figura
    python3 dipeptido.py --paso 10    # mas rapido, mapa mas basto
"""

import argparse
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolTransforms

# RDKit avisa de cosas que aqui no importan; se silencia para que la salida
# se lea bien
RDLogger.DisableLog("rdApp.*")

# SMILES del dipeptido de alanina: acetilo - L-alanina - N-metilamida.
# El [C@@H] marca la quiralidad: sin ella RDKit elige un enantiomero al azar y
# medio mapa sale en espejo, con phi positivo. Las proteinas son de L-aminoacidos.
DIPEPTIDO = "CC(=O)N[C@@H](C)C(=O)NC"

# Indices de los atomos que definen cada angulo del esqueleto. Salen de contar
# los atomos pesados en el orden en que aparecen en el SMILES de arriba:
# 0=CH3(ACE) 1=C(ACE) 2=O 3=N 4=CA 5=CB 6=C 7=O 8=N(NME) 9=CH3(NME)
PHI = (1, 3, 4, 6)      # C(i-1) - N - CA - C
PSI = (3, 4, 6, 8)      # N - CA - C - N(i+1)

RT = 8.314e-3 * 300     # kJ/mol a 300 K
KCAL = 4.184            # MMFF94 da kcal/mol y aqui se trabaja en kJ/mol


def construir():
    """Devuelve el dipeptido en 3D, ya minimizado, listo para girar diedros."""
    # de texto SMILES a molecula; todavia sin hidrogenos y sin coordenadas
    mol = Chem.MolFromSmiles(DIPEPTIDO)
    # los hidrogenos hacen falta: sin ellos el campo de fuerza no vale
    mol = Chem.AddHs(mol)
    # genera unas coordenadas 3D de partida razonables
    AllChem.EmbedMolecule(mol, randomSeed=2026)
    # baja esa geometria inicial a un minimo, para partir de algo limpio
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def energia_relajada(mol, phi, psi, pasos=200):
    """Energia en kJ/mol con phi y psi fijos y el resto de la molecula relajado.

    Girar los dos angulos y medir la energia sin mas da numeros enormes: al
    rotar el esqueleto los hidrogenos chocan entre si y el termino repulsivo se
    dispara. Un mapa de Ramachandran de verdad se calcula **fijando phi y psi y
    dejando que todo lo demas se acomode**, que es lo que hace esta funcion.
    """
    # se trabaja sobre una copia para no estropear la geometria de referencia
    copia = Chem.Mol(mol)
    conf = copia.GetConformer()
    # se coloca el esqueleto en los dos angulos pedidos
    rdMolTransforms.SetDihedralDeg(conf, *PHI, float(phi))
    rdMolTransforms.SetDihedralDeg(conf, *PSI, float(psi))
    # parametros MMFF94 para esta molecula
    props = AllChem.MMFFGetMoleculeProperties(copia)
    ff = AllChem.MMFFGetMoleculeForceField(copia, props)
    # las dos restricciones que impiden que phi y psi se muevan al relajar:
    # se permite medio grado de juego y se penaliza con una constante alta
    ff.MMFFAddTorsionConstraint(*PHI, False, phi - 0.5, phi + 0.5, 1.0e4)
    ff.MMFFAddTorsionConstraint(*PSI, False, psi - 0.5, psi + 0.5, 1.0e4)
    # relaja el resto de grados de libertad: enlaces, angulos y los demas diedros
    ff.Minimize(maxIts=pasos)
    # se resta la energia de las restricciones para quedarse con la de la molecula
    limpio = AllChem.MMFFGetMoleculeForceField(copia, props)
    return limpio.CalcEnergy() * KCAL


def mapa(mol, paso=10):
    """Rejilla de energia sobre (phi, psi). Devuelve (phi, psi, energia)."""
    # los dos angulos recorren la circunferencia completa
    rejilla = np.arange(-180, 180, paso)
    # matriz vacia donde se van guardando las energias
    e = np.zeros((len(rejilla), len(rejilla)))
    for i, psi in enumerate(rejilla):
        for j, phi in enumerate(rejilla):
            e[i, j] = energia_relajada(mol, phi, psi)
    # todo se mide respecto al punto mas bajo del mapa
    return rejilla, rejilla, e - e.min()


def minimo_cercano(mol, phi0, psi0, paso=8.0, tol=0.05, max_iter=120):
    """Baja por la pendiente de la superficie relajada desde (phi0, psi0).

    Es el mismo descenso por gradiente del butano, pero en dos variables: la
    pendiente tiene ahora una componente en phi y otra en psi.
    """
    phi, psi = float(phi0), float(psi0)
    h = 5.0
    for n in range(1, max_iter + 1):
        # cuanto sube la energia al mover solo phi
        gphi = (energia_relajada(mol, phi + h, psi) - energia_relajada(mol, phi - h, psi)) / (2 * h)
        # cuanto sube al mover solo psi
        gpsi = (energia_relajada(mol, phi, psi + h) - energia_relajada(mol, phi, psi - h)) / (2 * h)
        # si ya no hay pendiente apreciable en ninguna direccion, se para
        if max(abs(gphi), abs(gpsi)) < tol:
            break
        # se avanza en sentido contrario a la pendiente, en las dos variables
        phi = (phi - paso * gphi + 180) % 360 - 180
        psi = (psi - paso * gpsi + 180) % 360 - 180
    return phi, psi, energia_relajada(mol, phi, psi), n


def region(phi, psi):
    """Nombre del minimo del dipeptido **en vacio** donde cae un punto.

    Ojo: en vacio los minimos no son los del mapa de Ramachandran que se ve en
    las proteinas. Sin agua alrededor, la molecula se pliega sobre si misma
    formando un puente de hidrogeno interno (los minimos C7). Es un buen
    recordatorio de que el disolvente no es decorado.
    """
    if -120 <= phi <= -40 and 40 <= psi <= 110:
        return "C7eq · puente de H interno"
    if 40 <= phi <= 120 and -110 <= psi <= -40:
        return "C7ax · el mismo, en espejo"
    if -180 <= phi <= -100 and (psi >= 140 or psi <= -140):
        return "C5 · extendida"
    if -130 <= phi <= -40 and -80 <= psi <= 10:
        return "alfa derecha"
    return "otra"


def main():
    ap = argparse.ArgumentParser(description="Superficie de energia del dipeptido de alanina")
    ap.add_argument("--paso", type=int, default=10, help="paso del mapa, en grados")
    ap.add_argument("--sin-figura", action="store_true")
    args = ap.parse_args()

    mol = construir()
    print(f"MOLECULA: dipeptido de alanina (ACE-ALA-NME), {mol.GetNumAtoms()} atomos")
    print(f"Dos grados de libertad: phi {PHI} y psi {PSI}")

    print(f"\nMAPA DE ENERGIA (paso de {args.paso} grados, con el resto relajado)")
    xs, ys, e = mapa(mol, args.paso)
    print(f"  rejilla de {len(xs)}x{len(ys)} = {e.size} puntos")
    print(f"  energia maxima sobre el minimo: {e.max():.1f} kJ/mol = {e.max()/RT:.0f} RT")
    # cuenta que fraccion del mapa esta por debajo de 2 RT: la zona realmente
    # accesible a temperatura ambiente
    accesible = 100 * (e < 2 * RT).sum() / e.size
    print(f"  fraccion del mapa por debajo de 2 RT: {accesible:.1f} %")

    print("\nMINIMIZACION DESDE DISTINTOS PUNTOS DE PARTIDA")
    print(f"{'inicio':>14} {'final':>16} {'E (kJ/mol)':>12} {'pasos':>6}   region")
    # el cero de energias es el punto mas bajo del mapa que se acaba de calcular
    suelo = min(energia_relajada(mol, float(x), float(y))
                for y in ys[::3] for x in xs[::3])
    for phi0, psi0 in [(-60, -45), (-135, 135), (60, 45), (-90, 0), (0, 0), (150, -150)]:
        phi, psi, ener, n = minimo_cercano(mol, phi0, psi0)
        print(f"  ({phi0:>4},{psi0:>5}) {f'({phi:7.1f},{psi:7.1f})':>16} "
              f"{ener - suelo:>12.1f} {n:>6}   {region(phi, psi)}")

    print("\nLo que hay que leer: puntos de partida distintos caen en regiones")
    print("distintas del mapa, igual que en el butano. La diferencia es que aqui")
    print("ya no hay tres minimos en una linea, sino cuencas en una superficie.")

    if args.sin_figura:
        return

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5.4, 4.8), dpi=150)
    # el mapa se dibuja recortado a 40 kJ/mol: por encima es todo pared
    m = ax.contourf(xs, ys, np.clip(e, 0, 40), levels=20, cmap="YlGnBu_r")
    fig.colorbar(m, ax=ax, label="energia sobre el minimo (kJ/mol)")
    ax.set_xlabel("phi (grados)")
    ax.set_ylabel("psi (grados)")
    ax.set_xticks(range(-180, 181, 90))
    ax.set_yticks(range(-180, 181, 90))
    fig.tight_layout()
    fig.savefig("dipeptido-ramachandran.png", dpi=150)
    print("\nfigura guardada en dipeptido-ramachandran.png")


if __name__ == "__main__":
    main()

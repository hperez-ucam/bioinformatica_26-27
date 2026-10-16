#!/usr/bin/env python3
"""El dipeptido de alanina: de una curva a una superficie.

En butano.py la molécula tenia UN grado de libertad y la energía cabia en una
curva. Aqui hay DOS, los ángulos phi y psi del esqueleto peptidico, y la energía
ya es una superficie: el mapa de Ramachandran. Una proteína real tiene
3N−6 grados de libertad y no se puede dibujar; este es el último caso que si.

La molécula es ACE-ALA-NME, 22 átomos. Es el sistema de prueba clásico de la
dinámica molecular de biomoleculas.

Uso:
    python3 dipeptido.py              # mapa + minimizaciones + figura
    python3 dipeptido.py --paso 20    # más rápido, mapa más basto
    python3 dipeptido.py --paso 5     # más fino, unos dos minutos
"""

import argparse
from pathlib import Path
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolTransforms

# RDKit avisa de cosas que aquí no importan; se silencia para que la salida
# se lea bien
RDLogger.DisableLog("rdApp.*")

# SMILES del dipeptido de alanina: acetilo - L-alanina - N-metilamida.
# El [C@@H] marca la quiralidad: sin ella RDKit elige un enantiomero al azar y
# medio mapa sale en espejo, con phi positivo. Las proteinas son de L-aminoacidos.
DIPEPTIDO = "CC(=O)N[C@@H](C)C(=O)NC"

# Indices de los átomos que definen cada ángulo del esqueleto. Salen de contar
# los átomos pesados en el orden en que aparecen en el SMILES de arriba:
# 0=CH3(ACE) 1=C(ACE) 2=O 3=N 4=CA 5=CB 6=C 7=O 8=N(NME) 9=CH3(NME)
PHI = (1, 3, 4, 6)      # C(i-1) - N - CA - C
PSI = (3, 4, 6, 8)      # N - CA - C - N(i+1)

RT = 8.314e-3 * 300     # kJ/mol a 300 K
KCAL = 4.184            # MMFF94 da kcal/mol y aquí se trabaja en kJ/mol


# Nombres de residuo y de átomo que espera el campo de fuerza AMBER99SB-ILDN de
# GROMACS para los átomos pesados del dipeptido, en el orden del SMILES.
# Los hidrógenos no hacen falta: pdb2gmx los reconstruye con -ignh.
NOMBRES_PDB = [
    ("ACE", 1, "CH3"), ("ACE", 1, "C"), ("ACE", 1, "O"),
    ("ALA", 2, "N"), ("ALA", 2, "CA"), ("ALA", 2, "CB"), ("ALA", 2, "C"), ("ALA", 2, "O"),
    ("NME", 3, "N"), ("NME", 3, "CH3"),
]


def escribir_pdb(mol, ruta="dipeptido.pdb"):
    """Escribe los átomos pesados del dipeptido en formato PDB, para GROMACS.

    Asi la geometría de partida de la dinámica es la misma que se ha minimizado
    aquí, y no una escrita a mano.
    """
    conf = mol.GetConformer()
    lineas = []
    for i, (res, num, nombre) in enumerate(NOMBRES_PDB):
        # posición del átomo i en angstrom (un objeto con .x, .y y .z)
        pos = conf.GetAtomPosition(i)
        x, y, z = pos.x, pos.y, pos.z
        # formato PDB de columnas fijas: ATOM, número, nombre, residuo, cadena,
        # número de residuo, coordenadas, ocupacion, factor B y elemento
        lineas.append(f"ATOM  {i + 1:5d} {nombre:<4s} {res:>3s} A{num:4d}    "
                      f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {nombre[0]:>2s}")
    lineas += ["TER", "END"]
    Path(ruta).write_text("\n".join(lineas) + "\n")
    return ruta


def construir():
    """Devuelve el dipeptido en 3D, ya minimizado, listo para girar diedros."""
    # de texto SMILES a molécula; todavia sin hidrógenos y sin coordenadas
    mol = Chem.MolFromSmiles(DIPEPTIDO)
    # los hidrógenos hacen falta: sin ellos el campo de fuerza no vale
    mol = Chem.AddHs(mol)
    # genera unas coordenadas 3D de partida razonables
    AllChem.EmbedMolecule(mol, randomSeed=2026)
    # baja esa geometría inicial a un mínimo, para partir de algo limpio
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def energia_relajada(mol, phi, psi, pasos=2000):
    """Energía en kJ/mol con phi y psi fijos y el resto de la molécula relajado.

    Girar los dos ángulos y medir la energía sin más da números enormes: al
    rotar el esqueleto los hidrógenos chocan entre si y el termino repulsivo se
    dispara. Un mapa de Ramachandran de verdad se calcula **fijando phi y psi y
    dejando que todo lo demas se acomode**, que es lo que hace esta función.
    """
    # se trabaja sobre una copia para no estropear la geometría de referencia
    copia = Chem.Mol(mol)
    conf = copia.GetConformer()
    # se coloca el esqueleto en los dos ángulos pedidos
    rdMolTransforms.SetDihedralDeg(conf, *PHI, float(phi))
    rdMolTransforms.SetDihedralDeg(conf, *PSI, float(psi))
    # parámetros MMFF94 para esta molécula
    props = AllChem.MMFFGetMoleculeProperties(copia)
    ff = AllChem.MMFFGetMoleculeForceField(copia, props)
    # las dos restricciones que impiden que phi y psi se muevan al relajar:
    # se permite medio grado de juego y se penaliza con una constante alta
    ff.MMFFAddTorsionConstraint(*PHI, False, phi - 0.5, phi + 0.5, 1.0e4)
    ff.MMFFAddTorsionConstraint(*PSI, False, psi - 0.5, psi + 0.5, 1.0e4)
    # relaja el resto de grados de libertad: enlaces, ángulos y los demas diedros
    ff.Minimize(maxIts=pasos)
    # se recalcula la energía con un campo sin restricciones: queda solo la de
    # la molécula, sin el termino de penalización
    limpio = AllChem.MMFFGetMoleculeForceField(copia, props)
    return limpio.CalcEnergy() * KCAL


def mapa(mol, paso=10):
    """Rejilla de energia sobre (phi, psi). Devuelve (phi, psi, energia)."""
    # los dos ángulos recorren la circunferencia completa
    rejilla = np.arange(-180, 180, paso)
    # matriz vacia donde se van guardando las energías
    e = np.zeros((len(rejilla), len(rejilla)))
    for i, psi in enumerate(rejilla):
        for j, phi in enumerate(rejilla):
            e[i, j] = energia_relajada(mol, phi, psi)
    # todo se mide respecto al punto más bajo del mapa
    return rejilla, rejilla, e - e.min()


def minimo_cercano(mol, phi0, psi0, paso=8.0, tol=0.05, max_iter=120):
    """Baja por la pendiente de la superficie relajada desde (phi0, psi0).

    Es el mismo descenso por gradiente del butano, pero en dos variables: la
    pendiente tiene ahora una componente en phi y otra en psi.
    """
    phi, psi = float(phi0), float(psi0)
    h = 5.0
    for n in range(1, max_iter + 1):
        # cuanto sube la energía al mover solo phi
        gphi = (energia_relajada(mol, phi + h, psi) - energia_relajada(mol, phi - h, psi)) / (2 * h)
        # cuanto sube al mover solo psi
        gpsi = (energia_relajada(mol, phi, psi + h) - energia_relajada(mol, phi, psi - h)) / (2 * h)
        # si ya no hay pendiente apreciable en ninguna dirección, se para
        if max(abs(gphi), abs(gpsi)) < tol:
            break
        # se avanza en sentido contrario a la pendiente, en las dos variables
        phi = (phi - paso * gphi + 180) % 360 - 180
        psi = (psi - paso * gpsi + 180) % 360 - 180
    return phi, psi, energia_relajada(mol, phi, psi), n


def region(phi, psi):
    """Nombre del mínimo del dipeptido **en vacio** donde cae un punto.

    Ojo: en vacio los mínimos no son los del mapa de Ramachandran que se ve en
    las proteinas. Sin agua alrededor, la molécula se pliega sobre si misma
    formando un puente de hidrógeno interno (los mínimos C7). Es un buen
    recordatorio de que el disolvente no es decorado.
    """
    if -120 <= phi <= -40 and 40 <= psi <= 110:
        return "C7eq · puente de hidrógeno interno"
    if 40 <= phi <= 120 and -110 <= psi <= -40:
        return "C7ax · el mismo, en espejo"
    if -180 <= phi <= -100 and (psi >= 140 or psi <= -140):
        return "C5 · extendida"
    if -130 <= phi <= -40 and -80 <= psi <= 10:
        return "alfa derecha"
    if 40 <= phi <= 100 and 0 <= psi <= 80:
        return "alfa izquierda"
    return "otra"


def main():
    ap = argparse.ArgumentParser(description="Superficie de energia del dipeptido de alanina")
    ap.add_argument("--paso", type=int, default=10, help="paso del mapa, en grados")
    ap.add_argument("--sin-figura", action="store_true", help="no dibujar nada")
    ap.add_argument("--pdb", action="store_true",
                    help="escribir dipeptido.pdb para la dinámica con GROMACS y salir")
    args = ap.parse_args()

    mol = construir()
    if args.pdb:
        print("escrito", escribir_pdb(mol))
        return
    print(f"MOLÉCULA: dipéptido de alanina (ACE-ALA-NME), {mol.GetNumAtoms()} átomos")
    print(f"Dos grados de libertad: phi {PHI} y psi {PSI}")

    print(f"\nMAPA DE ENERGIA (paso de {args.paso} grados, con el resto relajado)")
    xs, ys, e = mapa(mol, args.paso)
    print(f"  rejilla de {len(xs)}x{len(ys)} = {e.size} puntos")
    print(f"  energía maxima sobre el mínimo: {e.max():.1f} kJ/mol = {e.max()/RT:.0f} RT")
    # cuenta que fracción del mapa esta por debajo de 2 RT: la zona realmente
    # accesible a temperatura ambiente
    accesible = 100 * (e < 2 * RT).sum() / e.size
    print(f"  fracción del mapa por debajo de 2 RT: {accesible:.1f} %")

    print("\nMINIMIZACIÓN DESDE DISTINTOS PUNTOS DE PARTIDA")
    print(f"{'inicio':>14} {'final':>16} {'E (kJ/mol)':>12} {'pasos':>6}   región")
    # primero se localizan los mínimos; el cero de energías es el más bajo de
    # ellos (C7eq), no un punto de la rejilla, que nunca cae justo en el fondo
    resultados = [minimo_cercano(mol, phi0, psi0)
                  for phi0, psi0 in [(-60, -45), (-135, 135), (60, 45), (-90, 0), (0, 0), (150, -150)]]
    suelo = min(r[2] for r in resultados)
    for (phi0, psi0), (phi, psi, ener, n) in zip(
            [(-60, -45), (-135, 135), (60, 45), (-90, 0), (0, 0), (150, -150)], resultados):
        print(f"  ({phi0:>4},{psi0:>5}) {f'({phi:7.1f},{psi:7.1f})':>16} "
              f"{ener - suelo:>12.1f} {n:>6}   {region(phi, psi)}")

    print("\nLo que hay que leer: puntos de partida distintos caen en regiones")
    print("distintas del mapa, igual que en el butano. La diferencia es que aquí")
    print("ya no hay tres mínimos en una linea, sino cuencas en una superficie.")

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

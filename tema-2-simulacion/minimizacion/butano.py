#!/usr/bin/env python3
"""Minimizacion de energia de un confórmero: el butano.

El caso mas simple de los que vemos en el Tema 2. La molecula tiene un unico
grado de libertad interesante, el diedro C1-C2-C3-C4, asi que su superficie de
energia potencial es una curva que cabe en una pantalla. Todo lo que se dice en
clase sobre minimizar se puede comprobar aqui en diez segundos.

MODELO DE ENERGIA
-----------------
Potencial de Ryckaert-Bellemans para el butano de atomo unido, que es el que
lleva GROMACS en su ejemplo clasico:

    V(psi) = sum_{n=0..5} C_n cos^n(psi),   psi = phi - 180 grados

con los seis coeficientes en kJ/mol. Es un potencial *autocontenido*: ya
incluye el efecto de las interacciones entre los metilos terminales, de modo
que no hay que sumarle ningun termino de van der Waals aparte. Por eso el
diedro basta para describir la molecula.

Lo que el modelo reproduce bien: los tres minimos, cual es el global, y que las
barreras son mucho mayores que RT. Lo que no reproduce bien: la barrera de la
conformacion eclipsada sin sale en 44,8 kJ/mol frente a los 19-27 kJ/mol que
dan los experimentos. Un campo de fuerza es una aproximacion ajustada, y
conviene saber de que tamano es su error.

Referencia: Ryckaert, J.-P., & Bellemans, A. (1975). Molecular dynamics of
liquid n-butane near its boiling point. Chemical Physics Letters, 30(1),
123-125.

Uso:
    python3 butano.py                 # barrido + minimizaciones + figura
    python3 butano.py --sin-figura    # solo los numeros, sin matplotlib
"""

import argparse
import numpy as np

# --- parametros del modelo ---------------------------------------------------

# Ryckaert-Bellemans para el butano, kJ/mol
C = np.array([9.2789, 12.1557, -13.1200, -3.0597, 26.2403, -31.4950])

B = 0.153                       # nm, enlace C-C
THETA = np.radians(112.7)       # angulo C-C-C
RT = 8.314e-3 * 300             # kJ/mol a 300 K


def energia(phi):
    """Energia del diedro, en kJ/mol. Vale cero en la conformacion anti."""
    psi = phi - np.pi
    c = np.cos(psi)
    return sum(C[n] * c ** n for n in range(6))


def distancia_14(phi):
    """Distancia entre los carbonos terminales, en nm.

    Sale de construir las coordenadas de los cuatro carbonos y simplificar:

        r14^2 = b^2 [ (1 - 2 cos t)^2 + 2 sin^2(t) (1 - cos phi) ]

    con b la longitud de enlace y t el angulo C-C-C. Toda la dependencia del
    diedro esta en (1 - cos phi): separacion maxima en anti (phi = 180) y
    minima en la eclipsada sin (phi = 0). No entra en la energia; se calcula
    para ver que lo que encarece la eclipsada es el acercamiento de los metilos.
    """
    return B * np.sqrt((1 - 2 * np.cos(THETA)) ** 2
                       + 2 * np.sin(THETA) ** 2 * (1 - np.cos(phi)))


def posiciones(phi):
    """Coordenadas de los cuatro carbonos, en nm. Comprueba distancia_14."""
    c2 = np.array([0.0, 0.0, 0.0])
    c3 = np.array([B, 0.0, 0.0])
    c1 = np.array([B * np.cos(THETA), B * np.sin(THETA), 0.0])
    d = np.array([-np.cos(THETA),
                  np.sin(THETA) * np.cos(phi),
                  np.sin(THETA) * np.sin(phi)])
    return c1, c2, c3, c3 + B * d


def gradiente(phi, h=1e-5):
    """Derivada dE/dphi por diferencias finitas centradas.

    Un programa real deriva la energia analiticamente; aqui se hace por
    diferencias para que se vea que la fuerza no es mas que la pendiente
    cambiada de signo.
    """
    return (energia(phi + h) - energia(phi - h)) / (2 * h)


def minimizar(phi0, paso=0.01, tol=1e-4, max_iter=20000):
    """Descenso por gradiente desde phi0. Devuelve (phi, energia, pasos).

    En cada iteracion el diedro se mueve en sentido contrario a la pendiente.
    El criterio de parada es que la fuerza baje de un umbral, igual que en
    GROMACS se para cuando Fmax < emtol. El algoritmo no salta barreras: por eso
    el resultado depende del punto de partida.
    """
    phi = float(phi0)
    for n in range(1, max_iter + 1):
        g = gradiente(phi)
        if abs(g) < tol:
            return phi, energia(phi), n
        phi -= paso * g
    return phi, energia(phi), max_iter


def grados(phi):
    """Angulo en [0, 360)."""
    return np.degrees(phi) % 360.0


def extremos(n=36000):
    """Minimos y maximos del barrido, comparando cada punto con sus vecinos.

    El barrido es circular: el punto 0 grados y el 360 son el mismo, asi que la
    comparacion se hace con indices modulo n. Si no, la eclipsada sin (0 grados)
    se perderia por caer justo en el borde.
    """
    phi = np.linspace(0, 2 * np.pi, n, endpoint=False)
    e = energia(phi)
    mins, maxs = [], []
    for i in range(n):
        izq, der = e[(i - 1) % n], e[(i + 1) % n]
        if e[i] < izq and e[i] < der:
            mins.append((np.degrees(phi[i]), e[i]))
        if e[i] > izq and e[i] > der:
            maxs.append((np.degrees(phi[i]), e[i]))
    return mins, maxs


def nombre_conformero(g):
    """anti, gauche+ o gauche-, segun el diedro en grados."""
    g = g % 360
    if abs(g - 180) < 30:
        return "anti"
    return "gauche+" if g < 180 else "gauche-"


def main():
    ap = argparse.ArgumentParser(description="Minimizacion del diedro del butano")
    ap.add_argument("--sin-figura", action="store_true", help="no dibujar nada")
    args = ap.parse_args()

    mins, maxs = extremos()

    print("BARRIDO DEL DIEDRO C1-C2-C3-C4")
    print(f"{'':<9}{'diedro':>9} {'E (kJ/mol)':>12} {'E/RT':>7}   r(C1...C4)")
    for g, e in mins:
        print(f"{'minimo':<9}{g:>8.1f}° {e:>12.2f} {e / RT:>7.2f}   "
              f"{distancia_14(np.radians(g)):.3f} nm   {nombre_conformero(g)}")
    for g, e in maxs:
        print(f"{'maximo':<9}{g:>8.1f}° {e:>12.2f} {e / RT:>7.2f}   "
              f"{distancia_14(np.radians(g)):.3f} nm")

    print("\nPOBLACION A 300 K (Boltzmann sobre los tres minimos)")
    pesos = np.exp(-np.array([e for _, e in mins]) / RT)
    pesos /= pesos.sum()
    for (g, _), p in zip(mins, pesos):
        print(f"  {nombre_conformero(g):<8} {g:6.1f}°   {100 * p:5.1f} %")
    anti = sum(p for (g, _), p in zip(mins, pesos) if nombre_conformero(g) == "anti")
    print(f"  -> anti {100 * anti:.0f} %, gauche {100 * (1 - anti):.0f} %")

    print("\nMINIMIZACION DESDE DISTINTOS PUNTOS DE PARTIDA")
    print(f"{'inicial':>8} {'final':>9} {'E final':>9} {'pasos':>7}   llega a")
    for g0 in (10, 70, 115, 125, 175, 200, 245, 255, 300, 350):
        phi_f, e_f, n = minimizar(np.radians(g0))
        gf = grados(phi_f)
        print(f"{g0:>7}° {gf:>8.1f}° {e_f:>9.2f} {n:>7}   "
              f"{nombre_conformero(gf)}{' (global)' if nombre_conformero(gf) == 'anti' else ' (local)'}")

    print("\nLo que hay que leer en esta tabla: 115 y 125 grados estan separados")
    print("por 10 grados y acaban en minimos distintos, porque entre ellos pasa")
    print("la barrera de 120 grados. El algoritmo no elige el mejor minimo: se")
    print("queda en el que tiene debajo.")

    if args.sin_figura:
        return

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    g = np.linspace(0, 360, 1441)
    fig, ax = plt.subplots(figsize=(8, 4.2), dpi=150)
    ax.plot(g, energia(np.radians(g)), lw=2, color="#004379")
    for gm, em in mins:
        ax.plot(gm, em, "o", color="#9A7200", zorder=5)
    ax.axhline(RT, color="#C7CDD2", lw=1)
    ax.text(6, RT + 1.2, "RT a 300 K", fontsize=9, color="#565656")
    ax.set_xlabel("diedro C1-C2-C3-C4 (grados)")
    ax.set_ylabel("energia (kJ/mol)")
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 60))
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    fig.tight_layout()
    fig.savefig("butano-barrido.png", dpi=150)
    print("\nfigura guardada en butano-barrido.png")


if __name__ == "__main__":
    main()

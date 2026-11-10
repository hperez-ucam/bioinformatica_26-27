# Tema 4 · Visualización de biomoléculas

| Fichero | Ejercicio | Dónde se ejecuta |
|---|---|---|
| `cargas_gasteiger.py` | E1: cargas parciales de una molécula pequeña con RDKit | Servidor de clase |
| `e2_rmsd.pml` | E2: dos valores de RMSD para el mismo par de estructuras | PyMOL, en el portátil |

## E1, en el servidor

```bash
# actualizar el repositorio y entrar en el directorio del tema
cd ~/bioinformatica_26-27 && git pull
cd tema-4-visualizacion
# primero el SMILES, entre comillas; después el nombre del fichero de salida
python3 cargas_gasteiger.py "C/C=C(/C=O)[C@@H](CC=O)CC(=O)OCCc1ccc(O)cc1" oleocantal
```

La figura queda en `oleocantal_gasteiger.png`. Para verla, se copia al portátil desde una terminal
del portátil (PowerShell en Windows, Terminal en Mac), cambiando `USUARIO` por vuestro usuario:

```bash
scp USUARIO@31.70.148.22:bioinformatica_26-27/tema-4-visualizacion/oleocantal_gasteiger.png .
```

## E2, en PyMOL

En PyMOL: `File > Run Script…` y elegir `e2_rmsd.pml`, o escribir las órdenes una a una en la línea
de órdenes de PyMOL.

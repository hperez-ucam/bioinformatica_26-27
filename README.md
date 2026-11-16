# Bioinformática 26-27

Material de clase de Bioinformática. Grado en Biotecnología, 3.º curso, UCAM. Profesor: Horacio
Pérez Sánchez.

Aquí están los scripts que se escriben en clase, los ficheros con los que trabajan y los
resultados de los cálculos que se enseñan en las diapositivas. No hay que copiar nada a mano de
la pantalla.

Los enunciados, las tareas y las diapositivas están en el campus virtual.

## Cómo se usa

En el servidor de clase, la primera vez:

```bash
git clone https://github.com/hperez-ucam/bioinformatica_26-27.git
cd bioinformatica_26-27
```

Al empezar cada clase, para tener lo último:

```bash
cd ~/bioinformatica_26-27 && git pull
```

Los scripts de Python se ejecutan con `python3` y los de Bash con `./nombre.sh`. Si un script de
Bash da «Permission denied», se le da permiso una vez con `chmod +x nombre.sh`.

## Qué hay en cada directorio

| Directorio | Contenido | Tema |
|---|---|---|
| `tema-1a-linux/` | 17 scripts de Bash, de menos a más, y dos ficheros PDB para practicar | 1 |
| `tema-2-simulacion/` | Minimización de energía, campos de fuerza, dinámica molecular y acoplamiento | 2 |
| `tema-3-hpc/` | Cola de trabajos, cribado en paralelo, LINPACK y medida de tiempos | 3 |
| `tema-4-visualizacion/` | Cargas parciales con RDKit y RMSD en PyMOL. Tiene su propio README | 4 |
| `tema-6-cribado/` | El caso completo de cribado virtual: EGFR con erlotinib. Tiene su propio README | 6 |

Los Temas 5 y 7 no tienen directorio. Se trabajan con bases de datos y servidores web, no con
scripts propios.

## Tema 1A · Linux y Bash

Los scripts están numerados en el orden en que se ven en clase. Cada uno enseña una cosa y se
puede leer entero en un minuto.

| Script | Qué enseña |
|---|---|
| `00-clear-date.sh` | El script más corto posible: dos comandos seguidos |
| `00b-contar-atom.sh`, `00c-lote-pdb.sh` | Contar átomos en un PDB, y repetirlo sobre muchos ficheros |
| `01-parametros.sh`, `02-shift.sh` | Argumentos que se le pasan al script |
| `03-if-grep.sh`, `esqueleto-if.sh`, `04-case.sh` | Decidir: `if` y `case` |
| `05-for-lista.sh`, `06-for-seq.sh` | Repetir: `for` sobre una lista y sobre un rango |
| `07-while.sh`, `08-until.sh` | Repetir mientras se cumpla, o hasta que se cumpla |
| `09-break.sh`, `10-continue.sh` | Salir de un bucle y saltarse una vuelta |
| `11-asignar-mils.sh` | Un caso real: asignar valores a una lista |
| `12-errores.sh` | Los errores que más se cometen, cometidos a propósito |
| `13-crlf.sh` | Por qué un script escrito en Windows falla en Linux |

`1HVR.pdb` y `ejemplo.pdb` son los ficheros sobre los que se practica.

## Tema 2 · Simulación

| Directorio | Qué hay |
|---|---|
| `minimizacion/` | Seis scripts, de menos a más: una curva en una dimensión, el butano a mano y con `scipy`, el dipéptido de alanina, el ibuprofeno con RDKit y un fármaco |
| `campo-de-fuerza/` | Dos extractos del campo de fuerza de GROMACS: enlaces y no enlazantes, para leer qué hay dentro |
| `dinamica/` | Una dinámica molecular completa del dipéptido de alanina, con sus dos ficheros de parámetros |
| `docking/` | Acoplamiento con AutoDock Vina: preparar el sistema, lanzar el cálculo y medir el RMSD |

## Tema 3 · Computación de altas prestaciones

| Script | Qué hace |
|---|---|
| `cola.sh` | Una cola de trabajos para el servidor compartido, hecha con `flock`. Envía, lista y consulta, como hacen `squeue` y `sacct` |
| `cribado.sh` | Acopla muchos ligandos a la vez y compara el tiempo frente a hacerlo uno a uno |
| `linpack_mini.py` | La prueba del Top500 a pequeña escala: resuelve un sistema lineal y cuenta los FLOP/s |
| `trabajo_largo.sh` | Un trabajo que tarda, para ver la diferencia entre tiempo real y tiempo de CPU |
| `muchos_ficheros.sh` | Qué le pasa al disco cuando se escriben miles de ficheros pequeños |
| `tabla_clase.sh` | Reúne las medidas de toda la clase en una sola tabla |
| `dos_maquinas.md` | Guía del ejercicio E6: las mismas órdenes en el servidor y en Google Colab |

Lo que generan los scripts se guarda en `salidas/`, que no se versiona.

## Con qué se calculó

Todo lo que hay en `tema-6-cribado/resultados/` se calculó en el servidor de clase el
09/11/2026, con:

- **AutoDock Vina 1.2.5**
- **RDKit 2026.09.1**
- **meeko 0.8.0**
- **Open Babel 3.1.1**

La dinámica molecular del Tema 2 usa **GROMACS 2023.3**.

Los tiempos que aparecen en los ficheros son reales y sirven para calcular cuánto tardará lo
vuestro.

## Lo que no está aquí

- Las soluciones de las tareas y de los exámenes.
- Los enunciados. Mandan los del campus virtual.
- Las diapositivas y los apuntes. Están en el campus virtual.
- Los programas (GROMACS, AutoDock Vina, RDKit). Están instalados en el servidor de clase.

## La otra asignatura

Estructura e Ingeniería de Proteínas tiene su propio repositorio:
[hperez-ucam/eip_26-27](https://github.com/hperez-ucam/eip_26-27).

# Tema 6 · Cribado virtual

Material del Tema 6 de Bioinformática: los scripts, el sistema de estudio y los resultados del
caso que se usa en clase (EGFR con erlotinib, entrada **1M17** del PDB, 26 compuestos).

Todo lo que hay aquí se ha calculado en el servidor de clase el 09/11/2026 con AutoDock Vina
1.2.5, RDKit 2026.09.1, meeko 0.8.0 y Open Babel 3.1.1.

## Qué hay en cada directorio

| Directorio | Contenido |
|---|---|
| `scripts/` | Los nueve scripts del tema, comentados línea a línea |
| `sistema/` | Entrada 1M17 del PDB, receptor preparado (`.pdbqt`), ligando cristalizado y centro de la caja |
| `ligandos/` | Los 26 compuestos preparados para el acoplamiento, y el informe de preparación |
| `resultados/` | Poses calculadas y las cuatro tablas de resultados |
| `libreria.smi` | Los 26 compuestos en formato SMILES, con su nombre |
| `pymol/` | **Las seis sesiones de PyMOL de las diapositivas**, listas para abrir |

## Las sesiones de PyMOL

Son las mismas figuras que salen en las diapositivas. Se abren con `pymol pymol/<fichero>.pse`,
o desde PyMOL con **File › Open**, y se pueden girar, medir y cambiar.

| Fichero | Qué tiene montado |
|---|---|
| `1M17-entrada_v1.pse` | La entrada del PDB tal cual: proteína, ligando cristalizado y las 20 aguas |
| `1M17-receptor-preparado_v1.pse` | El mismo encuadre ya sin aguas ni ligando: lo que se acopla |
| `1M17-caja_v1.pse` | La caja de búsqueda bien puesta (verde) y desplazada 15 Å (roja) |
| `1M17-redocking-superposicion_v1.pse` | El modo 1 de la semilla 7 sobre el ligando del cristal |
| `1M17-redocking-dos-semillas_v1.pse` | Las poses de las semillas 1 y 7 frente al cristal |
| `1M17-poses-candidatos_v1.pse` | PD153035 e imatinib frente al cristal, con Met769 marcada |

## Para qué sirve cada script

| Script | Qué hace |
|---|---|
| `preparar_receptor.py` | De un código PDB al receptor en `.pdbqt`, el ligando cristalizado y la caja |
| `preparar_libreria.py` | De un fichero de SMILES a un `.pdbqt` por compuesto, con informe |
| `derivados.py` | Genera derivados de un ligando sustituyendo una posición |
| `redocking.sh` | Vuelve a acoplar el ligando cristalizado y mide el RMSD (validación) |
| `calcular_rmsd.py` | RMSD de cada pose frente al ligando del cristal (viene del Tema 2) |
| `cribar.sh` | Acopla todos los compuestos de un directorio, seis a la vez |
| `similitud.py` | Ranking por similitud de Tanimoto frente a una o varias referencias |
| `consenso.py` | Une los dos rankings y añade las propiedades de Lipinski |
| `analizar_poses.py` | Residuos en contacto y enlaces de hidrógeno de cada pose |
| `ver_pose.py` | Figura del sitio de unión con PyMOL |

## Cómo reproducir el caso completo

```bash
cd tema-6-cribado
# 1. receptor, ligando de referencia y caja de busqueda
python3 scripts/preparar_receptor.py 1M17 AQ4
# 2. los 26 compuestos, de SMILES a PDBQT
python3 scripts/preparar_libreria.py libreria.smi ligandos/
# 3. validacion: el ligando del cristal contra su propio receptor
./scripts/redocking.sh sistema/1M17_rec.pdbqt sistema/1M17_lig_cristal.sdf \
    ligandos/erlotinib.pdbqt "$(cat sistema/caja.txt)"
# 4. cribado de la libreria (unos 200 s en el servidor de clase)
./scripts/cribar.sh sistema/1M17_rec.pdbqt "$(cat sistema/caja.txt)" \
    ligandos/ resultados/ 1 8
# 5. cribado basado en ligandos
python3 scripts/similitud.py libreria.smi erlotinib > resultados/ranking_lbvs.tsv
# 6. integracion de los dos rankings
python3 scripts/consenso.py libreria.smi ligandos/informe.tsv \
    resultados/ranking_sbvs.tsv resultados/ranking_lbvs.tsv
# 7. contactos de cada pose
python3 scripts/analizar_poses.py sistema/1M17_rec.pdb \
    sistema/1M17_lig_cristal.sdf resultados/
```

Los resultados que ya están en `resultados/` se obtuvieron con la **semilla 1** y
**exhaustividad 8**. Con otra semilla los valores cambian: eso es parte de lo que se estudia en
el tema.

## Para el proyecto final

El proyecto se hace con **otro sistema**, elegido por cada uno. Estos ficheros sirven para
comprobar que el procedimiento funciona antes de aplicarlo al sistema propio, y como modelo de lo
que hay que entregar: qué se calcula, en qué orden y qué se anota de cada paso.

# Dos máquinas remotas: el servidor de clase y Colab

Ejercicio E6 del Tema 3. Mismas preguntas a dos máquinas que no son vuestro portátil.

## 1. En el servidor de clase (terminal, por ssh)

```bash
# núcleos que ve el sistema
nproc
# modelo de procesador, núcleos por zócalo e hilos por núcleo
lscpu | head -20
# memoria total y libre, en unidades legibles
free -h
# ¿hay GPU de NVIDIA? si no la hay, la orden no existe
nvidia-smi
# cuántas personas están conectadas ahora mismo
who | wc -l
```

## 2. En Google Colab (navegador, notebook)

Abrid un notebook nuevo en <https://colab.research.google.com>. En una celda, el signo `!`
delante hace que la línea se ejecute como orden de Linux **en la máquina de Google**, no en
vuestro portátil:

```
!nproc
!lscpu | head -20
!free -h
!nvidia-smi
```

Repetidlo después de *Entorno de ejecución → Cambiar tipo de entorno de ejecución → GPU*.

## 3. AlphaFold2 en Colab (ColabFold)

Notebook: <https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb>.
Lo único que se cambia es la secuencia (`query_sequence`). Con una proteína corta
(menos de 100 residuos), anotad cuánto tarda y qué parte del trabajo se hace fuera de Colab
(la búsqueda de secuencias parecidas, MSA, va a un servidor de MMseqs2).

## 4. Tabla que hay que entregar

| | Servidor de clase | Colab sin GPU | Colab con GPU |
|---|---|---|---|
| Núcleos (`nproc`) | | | |
| Memoria total | | | |
| GPU y su memoria | | | |
| Cómo se entra | ssh, terminal | navegador | navegador |
| ¿Cuánto dura la sesión? | | | |
| ¿Quién más la usa? | | | |

Y tres frases: qué ventaja tiene cada forma de trabajar, por qué AlphaFold2 completo no cabe
en el servidor de clase y qué no subiríais nunca a Colab (regla I0).

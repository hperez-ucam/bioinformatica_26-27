# Dos máquinas remotas: el servidor de clase y Colab

Guía del ejercicio E6 y de la Tarea 3D del Tema 3 (diapositivas 39 y 40). Tiene dos partes:

1. Las mismas órdenes en dos máquinas remotas: el servidor de clase, por `ssh`, y una máquina de
   Google Colab, desde el navegador.
2. La misma predicción de estructura en las dos máquinas, con ColabFold, que ejecuta AlphaFold2.

Hacen falta la cuenta del servidor de clase y una cuenta de Google. El enunciado de la Tarea 3D,
en el campus virtual, indica qué se entrega y qué órdenes se añaden a las de esta guía.

## 1. Las órdenes en el servidor de clase

En una terminal, conectados al servidor de clase por `ssh`:

```bash
# hilos (CPU lógicas) que ve el sistema operativo
nproc
# hilos por núcleo (Thread) y núcleos por zócalo (Core)
lscpu | grep -E 'Thread|Core'
# memoria total y memoria disponible para un programa nuevo (columna available)
free -h | head -2
# usuarios conectados y carga media del último minuto, de los últimos 5 y de los últimos 15
uptime
# GPU de NVIDIA: modelo y memoria; si la orden no existe, no hay GPU
nvidia-smi
```

Qué se copia de cada salida:

| Orden | Qué se copia |
|---|---|
| `nproc` | El número. Son hilos, no núcleos |
| `lscpu` | `Thread(s) per core` (hilos por núcleo) y `Core(s) per socket` (núcleos) |
| `free -h` | De la línea `Mem:`, las columnas `total` y `available`. La columna `free` no es la memoria disponible |
| `uptime` | El primer número de `load average`: la carga media del último minuto |
| `nvidia-smi` | El modelo de la GPU y su memoria total. Si la respuesta es `Command 'nvidia-smi' not found, but can be installed with: …`, la máquina no tiene GPU de NVIDIA y no hay que instalar nada |

## 2. Las mismas órdenes en Colab

1. Abrid en el navegador el *notebook* de ColabFold (Colab lo llama «cuaderno»):
   <https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb>
2. Iniciad sesión con la cuenta de Google.
3. Añadid una celda de código con el botón «+ Código» (en inglés, «+ Code») y escribid en ella las
   órdenes del apartado 1 con el signo `!` delante:

   ```
   # en una celda de Colab, el signo ! ejecuta la línea como orden de Linux en la máquina de Google
   !nproc
   !lscpu | grep -E 'Thread|Core'
   !free -h | head -2
   !uptime
   !nvidia-smi
   ```

4. Ejecutad la celda con su botón de ejecución o con Ctrl+Intro. La primera vez, Colab avisa de que
   el *notebook* no es de Google: aceptad el aviso para continuar.

Este *notebook* pide una GPU al abrirse. Se comprueba, y se cambia, en el menú «Entorno de
ejecución» → «Cambiar tipo de entorno de ejecución» (en inglés, «Runtime» → «Change runtime
type»). Colab no garantiza la GPU: si en ese momento no puede concederla, lo indica y permite
seguir sin ella. Anotad en ese caso que la sesión no tenía GPU.

## 3. Tabla de las dos máquinas

| Dato | Servidor de clase | Colab |
|---|---|---|
| Hilos (`nproc`) | | |
| Hilos por núcleo y núcleos (`lscpu`) | | |
| Memoria total y disponible (`free -h`) | | |
| Carga media del último minuto (`uptime`) | | |
| GPU y su memoria (`nvidia-smi`) | | |

## 4. AlphaFold2 con ColabFold, en Colab

La proteína es la misma para todos: la ubiquitina humana, de 76 residuos. Su secuencia es pública
y está en `datos/ubiquitina.fasta`:

```bash
# la primera línea, que empieza por >, es el nombre; la segunda es la secuencia
cat ~/bioinformatica_26-27/tema-3-hpc/datos/ubiquitina.fasta
```

En el mismo *notebook* del apartado 2:

1. En la primera celda, borrad la secuencia de ejemplo del campo `query_sequence` y pegad la de la
   ubiquitina: solo las letras, sin la línea que empieza por `>`.
2. En el campo `jobname`, escribid `ubiquitina`.
3. En «Advanced settings», elegid `1` en el campo `num_recycles` (por omisión vale `3`). Los demás
   campos no se cambian.
4. Menú «Entorno de ejecución» → «Ejecutar todas» (en inglés, «Runtime» → «Run all»).
5. Esperad sin cerrar la pestaña. El *notebook* instala ColabFold, obtiene el alineamiento múltiple
   de secuencias de la proteína, descarga los parámetros de AlphaFold2 y predice cinco modelos. La
   primera ejecución tarda varios minutos.
6. Al terminar, el navegador ofrece descargar un fichero `.zip` con los resultados. No hace falta
   para la tarea.

Qué se copia. La celda «Run Prediction» escribe un registro. De ese registro se copian:

- las cinco líneas que terminan en `took …s (1 recycles)`, una por modelo, con el tiempo de cada
  uno;
- la línea que contiene `rank_001`, con el pLDDT y el pTM del mejor modelo.

Para ver juntas esas líneas, añadid una celda al final del *notebook* y ejecutadla:

```
# líneas del registro con el tiempo de cada modelo y con el mejor modelo
!grep -E 'took|rank_001' */log.txt
```

## 5. La misma predicción en el servidor de clase

En el servidor de clase, ColabFold solo está instalado en la cuenta del profesor, que lanza esta
orden en clase. El servidor no tiene GPU.

```bash
# un modelo y un reciclado; el resultado queda en el directorio resultados/
colabfold_batch --num-models 1 --num-recycle 1 datos/ubiquitina.fasta resultados/
```

Registro de la ejecución del 06/11/2026 (ColabFold 1.6.3):

```
2026-11-06 11:05:58,099 Running colabfold 1.6.3
2026-11-06 11:05:58,107 WARNING: no GPU detected, will be using CPU
2026-11-06 11:05:58,186 Found 5 citations for tools or databases
2026-11-06 11:05:58,186 Query 1/1: ubiquitina (length 76)
2026-11-06 11:06:03,498 Setting max_seq=512, max_extra_seq=5120
2026-11-06 11:06:53,043 alphafold2_ptm_model_1_seed_000 recycle=0 pLDDT=94.6 pTM=0.813
2026-11-06 11:07:37,519 alphafold2_ptm_model_1_seed_000 recycle=1 pLDDT=94.9 pTM=0.815 tol=0.0911
2026-11-06 11:07:37,520 alphafold2_ptm_model_1_seed_000 took 90.7s (1 recycles)
2026-11-06 11:07:37,525 reranking models by 'plddt' metric
2026-11-06 11:07:37,525 rank_001_alphafold2_ptm_model_1_seed_000 pLDDT=94.9 pTM=0.815
2026-11-06 11:07:37,728 Done
```

Medido con `/usr/bin/time`: 101 s de tiempo real, 446 s de tiempo de CPU y 3,0 GB de memoria como
máximo.

## 6. Términos del registro

| Término | Significado |
|---|---|
| MSA | Alineamiento múltiple de secuencias (*multiple sequence alignment*): secuencias parecidas a la de la proteína, alineadas con ella. ColabFold lo obtiene antes de predecir |
| Modelo | AlphaFold2 tiene cinco redes entrenadas. Cada una da una predicción: `model_1` a `model_5` |
| Reciclado | Pasada adicional de la red sobre su propia predicción (`recycle`) |
| pLDDT | Confianza de la predicción, residuo a residuo, de 0 a 100. El registro da la media de la proteína; por encima de 90 se considera muy alta |
| pTM | Estimación de la calidad del plegamiento completo, de 0 a 1 |
| `rank_001` | El mejor de los modelos calculados, ordenados por pLDDT |

## Si algo falla

- Colab no se abre con la cuenta de la universidad: usad una cuenta personal de Google.
- Colab no concede GPU: repetid más tarde o seguid sin GPU. La predicción funciona, pero tarda
  más, y el registro lo indica con la línea `no GPU detected, will be using CPU`.
- La sesión de Colab se cierra por inactividad: volved a conectar, comprobad que `query_sequence`
  y `num_recycles` conservan vuestros valores y ejecutad todas las celdas otra vez.
- En el servidor de clase, `nvidia-smi` responde `Command 'nvidia-smi' not found, but can be
  installed with: …` y una lista de paquetes: es la respuesta esperada. No hay que instalar nada.

Guía escrita el 06/11/2026 para la versión 1.6.3 de ColabFold. La ejecución del servidor de clase
es real; los nombres de los campos son los del *notebook* en esa fecha y pueden cambiar en
versiones posteriores.

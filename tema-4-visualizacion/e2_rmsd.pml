# e2_rmsd.pml — ejercicio E2 del Tema 4 (v3): dos RMSD del mismo par. Uso: pymol -cq e2_rmsd.pml
set fetch_path, .
fetch 1ake, async=0
fetch 4ake, async=0
remove not chain A
remove solvent
align 1ake, 4ake
align 1ake, 4ake, cycles=0

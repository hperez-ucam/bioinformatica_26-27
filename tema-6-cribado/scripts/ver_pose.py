# ver_pose.py - figura del bolsillo con el ligando cristalográfico y la mejor pose de un ligando.
#   pymol -cq ver_pose.py -- PD153035 figura_PD153035.png
#
# Antes: obabel -ipdbqt resultados/NOMBRE_poses.pdbqt -f 1 -l 1 -O poses_pdb/NOMBRE.pdb
import sys
from pymol import cmd

# nombre del ligando y fichero de imagen, tras los dos guiones de la línea de órdenes
nombre, imagen = sys.argv[-2], sys.argv[-1]

# ligando cristalográfico (erlotinib) y receptor; el cristal se carga antes para medir distancias
cmd.load("sistema/1M17_lig_cristal.sdf", "cristal")
cmd.load("sistema/1M17_rec.pdb", "receptor")
cmd.load(f"poses_pdb/{nombre}.pdb", "pose")
# se parte de cero: nada visible, ni hidrógenos
cmd.hide("everything")
# lazo (cartoon) solo de los residuos a menos de 8 Å del cristal, semitransparente
cmd.select("cerca", "byres (receptor within 8 of cristal)")
cmd.show("cartoon", "cerca")
cmd.color("grey80", "receptor")
cmd.set("cartoon_transparency", 0.6)
# residuos del bolsillo (a menos de 4 Å del cristal) como bastones finos
cmd.select("bolsillo", "byres (receptor within 4 of cristal) and not hydro")
cmd.show("sticks", "bolsillo")
cmd.color("wheat", "bolsillo and elem C")
cmd.set("stick_radius", 0.08, "bolsillo")
# ligando cristalográfico en gris y pose en verde, más gruesos que el bolsillo
cmd.show("sticks", "(cristal or pose) and not hydro")
cmd.color("grey40", "cristal and elem C")
cmd.color("tv_green", "pose and elem C")
cmd.set("stick_radius", 0.22, "cristal or pose")
# bisagra (Met769), Thr766 y Lys721 con etiqueta
for resi, texto in (("769", "Met769"), ("766", "Thr766"), ("721", "Lys721")):
    cmd.show("sticks", f"receptor and resi {resi} and not hydro")
    cmd.color("salmon", f"receptor and resi {resi} and elem C")
    cmd.label(f"receptor and resi {resi} and name CA", f'"{texto}"')
# aspecto de la figura: fondo blanco, etiquetas negras y sin sombras
cmd.set("label_size", 30)
cmd.set("label_color", "black")
cmd.set("label_position", (0, 2.5, 3))
cmd.set("ray_shadow", 0)
cmd.bg_color("white")
# vista a lo largo del ligando cristalográfico y encuadre del bolsillo
cmd.orient("cristal")
cmd.zoom("cristal", 7)
cmd.png(imagen, width=1600, height=1100, dpi=200, ray=1)
cmd.quit()

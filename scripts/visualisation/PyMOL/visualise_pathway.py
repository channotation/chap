from pymol import cmd, stored

#from wobj import draw_wobj
import wobj as wobj


#cmd.bg_color(color = "white")
cmd.bg_color(color = "black")



file_structure = "output.pdb"


###############################################################################
# LOAD AND VISUALISE MOLECULE
###############################################################################

# load structure into object named structure:
cmd.load(file_structure, "structure")

# hide starting representation:
cmd.hide("all")

# protein pore:
sel_protein = "pol"
cmd.color("silver", sel_protein)
cmd.show("cartoon", sel_protein)

# pore lining (but not facing) residues:
sel_pore_lining = "q > 0.5 and b < 0.5"
cmd.show("sticks", sel_pore_lining)
cmd.color("orange", sel_pore_lining)

# pore facing residues:
sel_pore_facing = "b > 0.5"
cmd.show("sticks", sel_pore_facing)
cmd.color("yellow", sel_pore_facing)


###############################################################################
# ADD PORE SURFACE
###############################################################################


# read OBJ data from file:
obj = wobj.import_wobj("output.obj")

# draw all groups in the OBJ file:
wobj.draw_wobj(obj)


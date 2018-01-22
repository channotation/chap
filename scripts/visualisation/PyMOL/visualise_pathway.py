
from pymol import cmd, stored
from pymol.cgo import *


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


cmd.hide("all")


#obj = [ BEGIN, TRIANGLES ]
#obj.extend([COLOR, 1, 0, 0, VERTEX, 0, 0, 0])
#obj.extend([COLOR, 1, 1, 0, VERTEX, 1, 0, 0])
#obj.extend([COLOR, 1, 0, 1, VERTEX, 0, 1, 0])
#obj.extend([END])





#obj = [
#   BEGIN, LINES,
#   COLOR, 1.0, 1.0, 1.0,
 
#   VERTEX, 0.0, 0.0, 0.0,
#   VERTEX, 1.0, 0.0, 0.0,
# 
#   VERTEX, 0.0, 0.0, 0.0,
#   VERTEX, 0.0, 2.0, 0.0,
# 
#   VERTEX, 0.0, 0.0, 0.0,
#   VERTEX, 00, 0.0, 3.0,
# 
#   END
#   ]                                                                                            



def triangle ():
    name="test";
    obj=[ BEGIN, TRIANGLES ]
    obj.extend([COLOR, 1.0, 0.0, 0.0, VERTEX, 0.0, 0.0, 0.0])
    obj.extend([COLOR, 1.0, 0.0, 0.0, VERTEX, 2.0, 0.0, 0.0])
    obj.extend([COLOR, 1.0, 0.0, 0.0, VERTEX, 0.0, 2.0, 0.0])
    obj.extend([END])
    cmd.delete("test")
    cmd.load_cgo(obj,"test")

    print obj
    
    return
cmd.extend("triangle",triangle)


triangle()

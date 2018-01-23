###############################################################################
# SETUP
###############################################################################

# load libraries:
from pymol import cmd                       # PyMOL visualisation commands
from pymol import __script__                # location of this script
from os import path                         # ability to find pathname
import sys                                  # ability to append to import path
sys.path.append(path.dirname(__script__))   # find wobj.py in script directory
import argparse                             # command line argument parsing
import wobj as wobj                         # import and draw OBJ meshes

# parse command line arguments: 
parser = argparse.ArgumentParser()
parser.add_argument(
    "-structure",
    nargs = "?",
    const = "output.pdb",
    default = "output.pdb")
parser.add_argument(
    "-surface",
    nargs = "?",
    const = "output.obj",
    default = "output.obj")
parser.add_argument(
    "-property",
    nargs = "?",
    const = None,
    default = None)
args = parser.parse_args()


###############################################################################
# LOAD AND VISUALISE MOLECULE
###############################################################################

# load structure into object named structure:
cmd.load(args.structure, "structure")

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
obj = wobj.import_wobj(args.surface)

# draw all groups in the OBJ file:
wobj.draw_wobj(obj, args.property)


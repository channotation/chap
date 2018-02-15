# CHAP - The Channel Annotation Package
# 
# Copyright (c) 2016 - 2018 Gianni Klesse, Shanlin Rao, Mark S. P. Sansom, and 
# Stephen J. Tucker
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


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


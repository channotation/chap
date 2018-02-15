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

# try finding wobj library in current working directory:
set WOBJ_FILE "wobj.tcl"
if { [file exists $WOBJ_FILE] == 0 } {

    # try finding it in same path as this script:
    set WOBJ_FILE [file join [file dirname [info script]] wobj.tcl]
    if { [file exists $WOBJ_FILE] == 0 } {
        
        # raise error if neither location was successful:
        error "ERROR: Could not find wobj.tcl in working directory ([pwd]) \
               or script directory ([file dirname [info script]])."
    }
}

# source library with Wavefront OBJ parser:
source $WOBJ_FILE

# global settings:
axes location Off
color Display Background white
display projection Orthographic
display depthcue off
display ambientocclusion on
display rendermode GLSL

# have arguments been passed from command line?
if { [llength $argv] == 3 } {
    
    # get file names from input arguments:
    set FILE_STRUCTURE [lindex $argv 0]
    set FILE_PORE_SURFACE [lindex $argv 1]
    set PROPERTY [lindex $argv 2]

} else {

    # have arguments been set in console?
    if { [info exists FILE_STRUCTURE] == 0 } {

        # set default name:
        set FILE_STRUCTURE "output.pdb"
    }
    if { [info exists FILE_PORE_SURFACE] == 0 } {

        # set default name:
        set FILE_PORE_SURFACE "output.obj"
    }
    if { [info exists PROPERTY] == 0 } {

        # set default name:
        set PROPERTY ""
    }
}

# check that required files exist:
if { [file exists $FILE_STRUCTURE] == 0 } {
    
    error "ERROR: Structure file $FILE_STRUCTURE does not exist!"
}
if { [file exists $FILE_PORE_SURFACE] == 0 } {
    
    error "ERROR: Pore surface file $FILE_PORE_SURFACE does not exist!"
}



###############################################################################
# LOAD AND VISUALISE MOLECULE
###############################################################################

# load molecule:
mol delete top
if { [molinfo top] == -1 } {
    mol new $FILE_STRUCTURE
}

# delete default representation:
mol delrep 0 top

# protein:
mol addrep top
mol modselect 0 top (protein)
mol modstyle 0 top NewRibbons 0.3 10.0 4.1 0
mol modcolor 0 top ColorID 6
mol modmaterial 0 top AOEdgy

# pore lining (but not facing) residues:
mol addrep top
mol modselect 1 top (occupancy > 0.5 and beta <= 0.5)
mol modstyle 1 top Licorice 0.3 12 12
mol modcolor 1 top ColorID 3
mol modmaterial 1 top AOEdgy

# pore facing residues:
mol addrep top
mol modselect 2 top (beta > 0.5)
mol modstyle 2 top Licorice 0.3 12 12
mol modcolor 2 top ColorID 4
mol modmaterial 2 top AOEdgy


###############################################################################
# ADD PORE SURFACE
###############################################################################

# import an OBJ file:
set obj [WOBJ::import_wobj $FILE_PORE_SURFACE]

# draw OBJ mesh:
WOBJ::draw_wobj $obj $PROPERTY


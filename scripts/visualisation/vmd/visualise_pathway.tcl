###############################################################################
# SETUP
###############################################################################

# check for existence of WOBJ library:
if { [file exists "wobj.tcl"] == 0 } {
    error "ERROR: File wobj.tcl does not exist in working directory!"
}

# source library with Wavefront OBJ parser:
source wobj.tcl

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
mol modselect 1 top (occupancy > 0 and beta <= 0)
mol modstyle 1 top Licorice 0.3 12 12
mol modcolor 1 top ColorID 3
mol modmaterial 1 top AOEdgy

# pore facing residues:
mol addrep top
mol modselect 2 top (beta > 0)
mol modstyle 2 top Licorice 0.3 12 12
mol modcolor 2 top ColorID 32
mol modmaterial 2 top AOEdgy


###############################################################################
# ADD PORE SURFACE
###############################################################################

# import an OBJ file:
set obj [WOBJ::import_wavefront_obj $FILE_PORE_SURFACE]

# draw OBJ mesh:
WOBJ::draw_wavefront_obj $obj $PROPERTY


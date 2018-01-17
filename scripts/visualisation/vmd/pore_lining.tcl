###############################################################################
# SETUP
###############################################################################

# global settings:
axes location Off
color Display Background white
display projection Orthographic
display depthcue off
display ambientocclusion on
display rendermode GLSL

# have file names been passed as arguments?
if { [llength $argv] == 2 } {
    
    # get file names from input arguments:
    set FILE_STRUCTURE [lindex $argv 0]
    set FILE_PORE_SURFACE [lindex $argv 1]

} else {

    # have file names been set in console?
    if { [info exists FILE_STRUCTURE] == 0 } {

        # set default name:
        set FILE_STRUCTURE "output.pdb"
    }
    if { [info exists FILE_PORE_SURFACE] == 0 } {

        # set default name:
        set FILE_PORE_SURFACE "output.obj"
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


# Import OBJ Data:
# -----------------------------------------------------------------------------


source ~/repos/chap/scripts/visualisation/vmd/wobj.tcl

# import an OBJ file:
set filename $FILE_PORE_SURFACE
set obj [WOBJ::import_wavefront_obj $filename]


# Color Scale:
# -----------------------------------------------------------------------------
# Uncomment whichever color scale you prefer, rdbu, puor, and brbg are 
# divergent color scales, the others map between white and their respective
# saturated hue. Your can also create your own color scale by creating a list
# of three equal langth lists, where each of these lists contains the R, G, and
# B values of your colour table (these should be in [0,1] rather than [0, 255].

#set scale_colors [color_scale_greys]
set scale_colors [color_scale_reds]
#set scale_colors [color_scale_greens]
set scale_colors [color_scale_blues]
#set scale_colors [color_scale_oranges]
#set scale_colors [color_scale_purples]
#set scale_colors [color_scale_rdbu]
#set scale_colors [color_scale_puor]
#set scale_colors [color_scale_brbg]

# optionally invert color scale:
set scale_colors [revert_color_scale $scale_colors]


# PLOT GROUP
# -----------------------------------------------------------------------------
# TODO comment here

#set group_name "pathway_radius" 
#set group_name "pathway_avg_radius" 
#set group_name "pathway_avg_pl_hydrophobicity" 
set group_name "pathway_avg_pf_hydrophobicity" 
#set group_name "pathway_avg_energy" 
#set group_name "pathway_avg_density" 


# Draw Pore
# -----------------------------------------------------------------------------




# draw OBJ mesh:
WOBJ::draw_wavefront_obj $obj $group_name $scale_colors

rotate x by 90






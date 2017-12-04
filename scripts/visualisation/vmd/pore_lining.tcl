###############################################################################
# SETUP
###############################################################################

# check that all required variables are set:
if { [info exists FILE_STRUCTURE] == 0 } {
    puts "ERROR: Variable FILE_STRUCTURE is not set!"
    return
}
if { [info exists FILE_PORE_LINING] == 0 } {
    puts "ERROR: Variable FILE_PORE_LINING is not set!"
    return
}
if { [info exists FILE_PORE_SURFACE] == 0 } {
    puts "ERROR: Variable FILE_PORE_SURFACE is not set!"
    return
}

# global settings:
axes location Off
color Display Background white
display projection Orthographic
display depthcue off
display ambientocclusion on
display rendermode GLSL
color change rgb 7 1.000000 1.000000 1.000000

# load molecule:
set file_prot $FILE_STRUCTURE
mol delete top
if { [molinfo top] == -1 } {
    mol new $file_prot
}


###############################################################################
# BASIC PROTEIN
###############################################################################

# delete default representation:
mol delrep 0 top

# protein:
mol addrep top
mol modselect 0 top (protein)
mol modstyle 0 top NewRibbons 0.3 10.0 4.1 0
mol modcolor 0 top ColorID 6
mol modmaterial 0 top AOEdgy


###############################################################################
# PORE LINING RESIDUES
###############################################################################

# process for importing data:
proc import_pore_lining {filename} {

    # prepare lists as data containers:
    set res_id {}
    set pore_lining {}
    set pore_facing {}

    # parse csv file line by line:
    set linenumber 1;
    set csvfile [open $filename r]
    while { [gets $csvfile line] >= 0 } {

        # assume first line is header:
        if { $linenumber == 1 } {

            puts "$line"

        } else {

            # split line at whitespaces:
            set linespl [regexp -all -inline {\S+} $line]

            # add data to containers:
            # (NOTE: VMD index shifted by one wrt Gromacs)
            lappend res_id [expr [lindex $linespl 1] + 1]
            lappend pore_lining [lindex $linespl 5]
            lappend pore_facing [lindex $linespl 6]
        }

        # increment line counter:
        incr linenumber
    }
    close $csvfile

    # return data:
    return [list $res_id $pore_lining $pore_facing]
}


# read pore lining data from file:
set filename $FILE_PORE_LINING
set pore_lining_data [import_pore_lining $filename]


# seperate into several containers:
set res_id [lindex $pore_lining_data 0]
set pore_lining [lindex $pore_lining_data 1]
set pore_facing [lindex $pore_lining_data 2]


# loop over all residues:
foreach rid $res_id pl $pore_lining pf $pore_facing {
break

    # is residue pore lining?
    if { $pl == 1 } {
       
#        puts "pf = $pf"

        # is residue also pore facing:
        if { $pf == 1 } {

            set repnum [molinfo top get numreps]

            mol addrep top
#            mol modselect $repnum top (resid $rid and name CA)
            mol modselect $repnum top (resid $rid)
#            mol modstyle $repnum top VDW 0.5 12
            mol modstyle $repnum top Licorice
            mol modcolor $repnum top ColorID 32
            mol modmaterial $repnum top AOChalky
        
        } else {

            set repnum [molinfo top get numreps]

            mol addrep top
#            mol modselect $repnum top (resid $rid and name CA)
            mol modselect $repnum top (resid $rid)
#            mol modstyle $repnum top VDW 0.5 12
            mol modstyle $repnum top Licorice
            mol modcolor $repnum top ColorID 3
            mol modmaterial $repnum top AOChalky

        }
    }
}


###############################################################################
# PORE SURFACE IMPORT
###############################################################################

draw color blue
#draw material Transparent
# FIXME


###############################################################################
# PORE HOLE RESULT
###############################################################################


# Import OBJ Data:
# -----------------------------------------------------------------------------

# import the custom OBJ library:
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
#set scale_colors [color_scale_blues]
#set scale_colors [color_scale_oranges]
#set scale_colors [color_scale_purples]
#set scale_colors [color_scale_rdbu]
#set scale_colors [color_scale_puor]
#set scale_colors [color_scale_brbg]

# optionally invert color scale:
set scale_colors [revert_color_scale $scale_colors]


# Draw Pore
# -----------------------------------------------------------------------------

# draw OBJ mesh:
WOBJ::draw_wavefront_obj $obj $scale_colors








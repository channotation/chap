namespace eval WOBJ {

# global variables in this namespace:
set BASENAME obj.obj


# Process for importing OBJ file.
#
# INPUT:
# filename - name of obj file to import:
proc import_wavefront_obj {filename} {

    # prepare list of vertex coordinates:
    set vertices_x {}
    set vertices_y {}
    set vertices_z {}

    # prepare a list of vertex normals:
    set vertex_normals_x {}
    set vertex_normals_y {}
    set vertex_normals_z {}

    # prepare list of list of vertex indices:
    set vertex_idx_a {}
    set vertex_idx_b {}
    set vertex_idx_c {}

    # prepare list of vertex normal indices:
    set vertex_normal_idx_a {}
    set vertex_normal_idx_b {}
    set vertex_normal_idx_c {}

    # parse OBJ line by line:
    set linenumber 1;
    set objfile [open $filename r]
    while { [gets $objfile line] >= 0 } {
      
        # split line at space:
        set linespl [regexp -all -inline {\S+} $line]

        # check for vertex index:
        if { [lindex $linespl {0}] == "v"} {

            # add vertex to list:
            lappend vertices_x [lindex $linespl {1}]
            lappend vertices_y [lindex $linespl {2}] 
            lappend vertices_z [lindex $linespl {3}]
        }

        # check for vertex normal index:
        if { [lindex $linespl {0}] == "vn"} {

            # add vertex normal to list:
            lappend vertex_normals_x [lindex $linespl {1}]
            lappend vertex_normals_y [lindex $linespl {2}] 
            lappend vertex_normals_z [lindex $linespl {3}]
        }


        # check for face entries:
        if { [lindex $linespl 0] == "f"} {
      
            # sanity check:
            if { [expr [llength $linespl] > 4] } {

                puts "WARNING: All faces must be triangles!"
                puts "in line $linenumber:"
                puts $line
                break
            }

            # have vertex normals?
            if { [string match *//* $linespl] } {
                
                # add vertex to list:
                # NOTE: OBJ is one-based, Tcl is zero-based:
                lappend vertex_idx_a [expr [lindex [split [lindex $linespl {1}] "//"] 0] - 1]
                lappend vertex_idx_b [expr [lindex [split [lindex $linespl {2}] "//"] 0] - 1]
                lappend vertex_idx_c [expr [lindex [split [lindex $linespl {3}] "//"] 0] - 1]

                # add vertex normals to list:
                set vertspl_a [split [lindex $linespl {1}] "//"] 
                set vertspl_b [split [lindex $linespl {2}] "//"] 
                set vertspl_c [split [lindex $linespl {3}] "//"] 
                lappend vertex_normal_idx_a [expr [lindex $vertspl_a 2] - 1]
                lappend vertex_normal_idx_b [expr [lindex $vertspl_b 2] - 1]
                lappend vertex_normal_idx_c [expr [lindex $vertspl_c 2] - 1]

            } else {

                # add vertex to list:
                lappend vertex_idx_a [expr [lindex $linespl {1}] - 1]
                lappend vertex_idx_b [expr [lindex $linespl {2}] - 1]
                lappend vertex_idx_c [expr [lindex $linespl {3}] - 1]
            }
        }

        # increment line number counter:
        incr linenumber
    }
    close $objfile

    # return data:
    return [list $vertices_x $vertices_y $vertices_z \
                 $vertex_idx_a $vertex_idx_b $vertex_idx_c \
                 $vertex_normals_x $vertex_normals_y $vertex_normals_z \
                 $vertex_normal_idx_a $vertex_normal_idx_b $vertex_normal_idx_c]
}


# Import OBJ file for given frame. Is given a basename and will load file
# basename.framenumber.
#
# INPUT:
# basename - basename for file to be loaded
proc import_wavefront_obj_frame {basename} {

    # append frame number to file name:
    set molid [molinfo top]
    set frnr [molinfo $molid get frame]
    set filename "$basename.$frnr"

    puts "$filename"

    # load corresponding obj file:
    return [import_wavefront_obj $filename]
}


# Function that draws an OBJ object.
#
# INPUT:
# obj - list of lists containing vertices and faces
proc draw_wavefront_obj {obj} {

    # extract information from OBJ object:
    set vert_x [lindex $obj 0]
    set vert_y [lindex $obj 1]
    set vert_z [lindex $obj 2]
    set vertex_idx_a [lindex $obj 3]
    set vertex_idx_b [lindex $obj 4]
    set vertex_idx_c [lindex $obj 5]
    set vert_norm_x [lindex $obj 6]
    set vert_norm_y [lindex $obj 7]
    set vert_norm_z [lindex $obj 8]
    set vertex_norm_idx_a [lindex $obj 9]
    set vertex_norm_idx_b [lindex $obj 10]
    set vertex_norm_idx_c [lindex $obj 11]

    # sanity checks:
    if { [llength vertex_idx_a] != [llength vertex_idx_b] || \
         [llength vertex_idx_a] != [llength vertex_idx_c] } {

        error "Enequal length vertex index lists."
    }
    if { [llength vertex_norm_idx_a] != [llength vertex_norm_idx_b] || \
         [llength vertex_norm_idx_a] != [llength vertex_norm_idx_c] } {

        error "Enequal length vertex index lists."
    }

    set switch 1

    # have vertex normals?
    set num_vert_idx [llength $vertex_idx_a]
    set num_norm_idx [llength $vertex_norm_idx_a]
    if { $num_vert_idx == $num_norm_idx } {
    
        # loop over faces and draw them:
        foreach ia $vertex_idx_a ib $vertex_idx_b ic $vertex_idx_c {

            set switch [expr $switch * -1]
            if { $switch == 1 } {
                draw color red
            } else { 
                draw color blue
            }

            draw trinorm "[lindex $vert_x $ia] [lindex $vert_y $ia] [lindex $vert_z $ia]" \
                         "[lindex $vert_x $ib] [lindex $vert_y $ib] [lindex $vert_z $ib]" \
                         "[lindex $vert_x $ic] [lindex $vert_y $ic] [lindex $vert_z $ic]" \
                         "[lindex $vert_norm_x $ia] [lindex $vert_norm_y $ia] [lindex $vert_norm_z $ia]" \
                         "[lindex $vert_norm_x $ib] [lindex $vert_norm_y $ib] [lindex $vert_norm_z $ib]" \
                         "[lindex $vert_norm_x $ic] [lindex $vert_norm_y $ic] [lindex $vert_norm_z $ic]" 
        }

    } else {
        
        # loop over faces and draw them:
        foreach ia $vertex_idx_a ib $vertex_idx_b ic $vertex_idx_c {

            set switch [expr $switch * -1]
            if { $switch == 1 } {
                draw color red
            } else { 
                draw color blue
            }
            
            draw triangle "[lindex $vert_x $ia] [lindex $vert_y $ia] [lindex $vert_z $ia]" \
                          "[lindex $vert_x $ib] [lindex $vert_y $ib] [lindex $vert_z $ib]" \
                          "[lindex $vert_x $ic] [lindex $vert_y $ic] [lindex $vert_z $ic]" 
        }

    }
}


# Callback function that redraws OBJ when frame changes.
#
# INPUTS:
# args - callback arguments
# OBJ_GLOBAL_BASENAME - basename ob obj file, needs to be set as globally
proc draw_obj_frame_callback {args} {

    global WOBJ::BASENAME

    # load obj data:
    set obj [import_wavefront_obj_frame $BASENAME]

    # delete all current drawings:
    draw delete all
    
    # draw obj:
    draw_wavefront_obj $obj
}

}


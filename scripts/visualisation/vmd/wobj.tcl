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

    # prepare list of list of vertex indices:
    set faces_ia {}
    set faces_ib {}
    set faces_ic {}

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

        # check for face entries:
        if { [lindex $linespl 0] == "f"} {
           
            # sanity check:
            if { [expr [llength $linespl] > 4] } {

                puts "WARNING: All faces must be triangles!"
                puts "in line $linenumber:"
                puts $line
                break
            }

            # add vertex to list:
            # NOTE: OBJ is one-based, Tcl is zero-based:
            lappend faces_ia [expr [lindex $linespl {1}] - 1]
            lappend faces_ib [expr [lindex $linespl {2}] - 1]
            lappend faces_ic [expr [lindex $linespl {3}] - 1]
        }

        # increment line number counter:
        incr linenumber
    }
    close $objfile

    # return data:
    return [list $vertices_x $vertices_y $vertices_z $faces_ia $faces_ib $faces_ic]
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
    set faces_ia [lindex $obj 3]
    set faces_ib [lindex $obj 4]
    set faces_ic [lindex $obj 5]

    # loop over faces and draw them:
    foreach ia $faces_ia ib $faces_ib ic $faces_ic {

        draw triangle "[lindex $vert_x $ia] [lindex $vert_y $ia] [lindex $vert_z $ia]" \
                      "[lindex $vert_x $ib] [lindex $vert_y $ib] [lindex $vert_z $ib]" \
                      "[lindex $vert_x $ic] [lindex $vert_y $ic] [lindex $vert_z $ic]" 
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


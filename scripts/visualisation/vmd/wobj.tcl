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
    set vertices_w {}

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
            lappend vertices_w [lindex $linespl {4}]
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
    return [list $vertices_x $vertices_y $vertices_z $vertices_w \
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



# Returns the index of the color associated with a given scalar value in the 
# given color map.
#
# INPUT:
# scalar - scalar in the range [0, 1]
proc color_index_from_scalar {scalar color_lookup_table} {

    # sanity checks:
    if { $scalar < 0.0 } {
        error "Scalar input to color scale must be in unit interval."
    }
    if { $scalar > 1.0 } {
        error "Scalar input to color scale must be in unit interval."
    }

    # return index of associated color:
    set num_colors [llength $color_lookup_table]
    set interval_length [expr 1.0 / ($num_colors)]
    set lookup_table_index [expr int(floor($scalar / $interval_length))]
    return [lindex $color_lookup_table $lookup_table_index]
}


# Builds a color lookup table from a given set of colors. This is returned as
# list of color indices. The colors used in the lookup table are always the 
# last indeces of the VMD colors index set (i.e. 1056 - numColors to 1056).
#
# INPUT:
# color_rgb - list of three lists, containing R, G, and B values of colors.
proc setup_color_lookup_table { color_rgb } {

    # extract color lists:
    set col_r [lindex $color_rgb 0]
    set col_g [lindex $color_rgb 1]
    set col_b [lindex $color_rgb 2]

    # sanity checks:
    if { [llength $col_r] != [llength $col_b] } { 
        error "Inconsistent color definitions."
    }
    if { [llength $col_r] != [llength $col_g] } {
        error "Inconsistent color definitions."
    }
    if { [llength $col_r] > 100 } {
        error "Lookup table can not have more than 100 colors."
    }
    if { [llength $col_r] < 1 } { 
        error "Lookup table can not have less than 1 color."
    }

    # prepare table of indices:
    set num_colors [llength $col_r]
    set table {}
    for { set i 0 } { $i < $num_colors } { incr i } {
        lappend table [expr [colorinfo max] - $num_colors + $i]
    }

    # set table colours to given values:
    set table_size [llength $table]    
    for { set i 0 } { $i < $table_size } { incr i } {

        set r [lindex $col_r $i]
        set b [lindex $col_b $i]
        set g [lindex $col_g $i]
        color change rgb [lindex $table $i] $r $g $b
    }

    # return the indices:
    return $table;
}


# Function that draws an OBJ object.
#
# INPUT:
# obj - list of lists containing vertices and faces
proc draw_wavefront_obj {obj scale_colors} {

    # extract information from OBJ object:
    set vert_x [lindex $obj 0]
    set vert_y [lindex $obj 1]
    set vert_z [lindex $obj 2]
    set vert_w [lindex $obj 3]
    set vertex_idx_a [lindex $obj 4]
    set vertex_idx_b [lindex $obj 5]
    set vertex_idx_c [lindex $obj 6]
    set vert_norm_x [lindex $obj 7]
    set vert_norm_y [lindex $obj 8]
    set vert_norm_z [lindex $obj 9]
    set vertex_norm_idx_a [lindex $obj 10]
    set vertex_norm_idx_b [lindex $obj 11]
    set vertex_norm_idx_c [lindex $obj 12]

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

    # prepare the color lookup table:
    set color_lookup_table [setup_color_lookup_table $scale_colors]
    
    # have vertex normals?
    set num_vert_idx [llength $vertex_idx_a]
    set num_norm_idx [llength $vertex_norm_idx_a]
    if { $num_vert_idx == $num_norm_idx } {
    
        # loop over faces and draw them:
        foreach ia $vertex_idx_a ib $vertex_idx_b ic $vertex_idx_c {

            # get vertex weights:
            set wa [lindex $vert_w $ia]
            set wb [lindex $vert_w $ib]
            set wc [lindex $vert_w $ic]

            # calculate face weight:
            set face_weight [expr ($wa + $wb + $wc)/3.0]

            # change color index accordingly:
            draw color [color_index_from_scalar $face_weight $color_lookup_table]

            # draw triangle with normals:
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
#                draw color red
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


# Reverts order of color scale.
#
# INPUT:
# color_scale - a list of three lists containing R, G, and B values
proc revert_color_scale { color_scale } {
    
    set col_r [lreverse [lindex $color_scale 0]]
    set col_g [lreverse [lindex $color_scale 1]]
    set col_b [lreverse [lindex $color_scale 2]]

    puts $col_r
    puts $col_g
    puts $col_b

    return [list $col_r $col_g $col_b]
}


# Returns RGB values for grey scale.
#
# INPUT:
# none
proc color_scale_greys { } {

    set col_r {}
    set col_g {}
    set col_b {}
    
    set num_colors 100
    for { set i 0 } { $i < $num_colors } { incr i } {

        lappend col_r [expr 1.0/$num_colors*$i] 
        lappend col_g [expr 1.0/$num_colors*$i] 
        lappend col_b [expr 1.0/$num_colors*$i] 
    }

    # return as list of lists:
    set color_rgb [list $col_r $col_g $col_b]
    return $color_rgb
}


# Returns RGB values created from LAB space interpolation of the "Blues" 
# color scale of RColorBrewer.
#
# INPUT:
# none
proc color_scale_blues { } {

    set r {0.96863 0.95891 0.94963 0.94073 0.93220 0.92400 0.91609 0.90846 \
           0.90105 0.89385 0.88681 0.87991 0.87312 0.86639 0.85969 0.85297 \
           0.84618 0.83927 0.83220 0.82492 0.81737 0.80952 0.80131 0.79269 \
           0.78361 0.77402 0.76389 0.75324 0.74207 0.73042 0.71829 0.70570 \
           0.69268 0.67924 0.66538 0.65114 0.63651 0.62151 0.60615 0.59048 \
           0.57453 0.55835 0.54199 0.52550 0.50894 0.49237 0.47585 0.45946 \
           0.44328 0.42740 0.41193 0.39690 0.38233 0.36819 0.35448 0.34117 \
           0.32826 0.31571 0.30351 0.29163 0.28001 0.26863 0.25743 0.24637 \
           0.23542 0.22459 0.21386 0.20324 0.19272 0.18229 0.17196 0.16173 \
           0.15162 0.14165 0.13184 0.12222 0.11282 0.10363 0.09466 0.08593 \
           0.07746 0.06926 0.06139 0.05388 0.04680 0.04022 0.03443 0.02974 \
           0.02608 0.02339 0.02158 0.02060 0.02035 0.02076 0.02172 0.02315 
           0.02494 0.02698 0.02916 0.03137}
    set g {0.98431 0.97931 0.97431 0.96929 0.96426 0.95922 0.95417 0.94910 \
           0.94402 0.93893 0.93381 0.92867 0.92351 0.91833 0.91312 0.90792 \
           0.90271 0.89752 0.89235 0.88721 0.88212 0.87708 0.87211 0.86720 \
           0.86238 0.85765 0.85300 0.84838 0.84374 0.83903 0.83420 0.82921 \
           0.82400 0.81852 0.81274 0.80659 0.80005 0.79306 0.78560 0.77769 \
           0.76941 0.76079 0.75191 0.74282 0.73357 0.72423 0.71483 0.70545 \
           0.69612 0.68690 0.67785 0.66896 0.66020 0.65153 0.64290 0.63430 \
           0.62567 0.61697 0.60818 0.59926 0.59017 0.58088 0.57134 0.56155 \
           0.55154 0.54132 0.53094 0.52043 0.50982 0.49913 0.48842 0.47769 \
           0.46698 0.45633 0.44576 0.43530 0.42494 0.41468 0.40449 0.39436 \
           0.38427 0.37422 0.36418 0.35415 0.34412 0.33406 0.32397 0.31384 \
           0.30366 0.29343 0.28315 0.27281 0.26242 0.25197 0.24147 0.23092 \
           0.22032 0.20967 0.19897 0.18824}
    set b {1.00000 0.99671 0.99360 0.99065 0.98785 0.98520 0.98266 0.98025 \
           0.97793 0.97569 0.97354 0.97144 0.96939 0.96737 0.96537 0.96335 \
           0.96129 0.95915 0.95692 0.95457 0.95206 0.94938 0.94649 0.94336 \
           0.93998 0.93631 0.93235 0.92815 0.92375 0.91919 0.91453 0.90982 \
           0.90509 0.90039 0.89577 0.89128 0.88696 0.88285 0.87899 0.87534 \
           0.87188 0.86854 0.86529 0.86207 0.85886 0.85560 0.85224 0.84876 \
           0.84511 0.84124 0.83712 0.83275 0.82815 0.82335 0.81838 0.81326 \
           0.80804 0.80273 0.79735 0.79195 0.78654 0.78115 0.77581 0.77052 \
           0.76529 0.76007 0.75487 0.74965 0.74441 0.73911 0.73375 0.72830 \
           0.72275 0.71708 0.71128 0.70532 0.69917 0.69280 0.68615 0.67920 \
           0.67189 0.66419 0.65606 0.64746 0.63836 0.62870 0.61847 0.60762 \
           0.59612 0.58394 0.57106 0.55744 0.54307 0.52792 0.51196 0.49519 
           0.47758 0.45912 0.43980 0.41961}

    # return as list of lists:
    set color_rgb [list $r $g $b]
    return $color_rgb
}

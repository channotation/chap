namespace eval WOBJ {

# global variables in this namespace:
set BASENAME obj.obj


# Parses Wavefront MTL file and builds a dictionary of material properties.
#
# INPUT:
# filename - name of MTL file to import
proc import_wavefront_mtl {filename} {

	# handle empty file name:
	if { $filename == "" } {

		dict set colour_indeces undefined_material [colorinfo max]	
		return 
	}
	
    # parse MTL line by line:
    set linenumber 1;
    set mtlfile [open $filename r]
    while { [gets $mtlfile line] >= 0 } {

        # split line at space:
        set linespl [regexp -all -inline {\S+} $line]

		# check for new group:
		if { [lindex $linespl {0}] == "newmtl" } {
			
			set crnt_mtl_name [lindex $linespl {1}]
			puts $crnt_mtl_name
		}

		# ambient colour:
		if { [lindex $linespl {0}] == "Ka" } {
			
			dict set materials $crnt_mtl_name ka_r [lindex $linespl {1}]
			dict set materials $crnt_mtl_name ka_g [lindex $linespl {2}]
			dict set materials $crnt_mtl_name ka_b  [lindex $linespl {3}]
		}

		# diffuse colour:
		if { [lindex $linespl {0}] == "Kd" } {
			
			dict set materials $crnt_mtl_name kd_r [lindex $linespl {1}]
			dict set materials $crnt_mtl_name kd_g [lindex $linespl {2}]
			dict set materials $crnt_mtl_name kd_b  [lindex $linespl {3}]
		}

	}
	close $mtlfile

	# get number of materials:
	set num_materials [llength [dict keys $materials]]
	
	# sanity check:
	if { $num_materials > [expr [colorinfo max] - [colorinfo num]] } {
		error "Number of materials exceeds number of available slots!"
	}
	
	# loop over materials:
	set idx [colorinfo num]
	dict for {mtl def} $materials {

		# change color definition:
		dict with def {
			color change rgb $idx $ka_r $ka_g $ka_b
		}

	    # add to index dictionary:
		dict set colour_indeces $mtl $idx

		# increment colour index:
		incr idx
	}
	
	# return dictionary of colour indeces:
	return $colour_indeces
}


# Process for importing OBJ file.
#
# INPUT:
# filename - name of OBJ file to import
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

	# prepare list of groups:
	set group_names {}
	set crnt_group_name "undefined_group"

	# prepare list of materials:
	set mtl_lib ""
	set mtl_names {}
	set crnt_mtl_name "undefined_material"

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

		# check for material libraries:
		if { [lindex $linespl {0}] == "mtllib" } {
			set mtl_lib [lindex $linespl {1}]
		}

		# check for new group:
		if { [lindex $linespl {0}] == "g" } {
			
			set crnt_group_name [lindex $linespl {1}]
			puts $crnt_group_name
		}

		# check for new material:
		if { [lindex $linespl {0}] == "usemtl" } {
			set crnt_mtl_name [lindex $linespl {1}]
		}

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

			# add group and material label for each face:
			lappend group_names $crnt_group_name
			lappend mtl_names $crnt_mtl_name

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
                 $vertex_normal_idx_a $vertex_normal_idx_b $vertex_normal_idx_c \
				 $group_names $mtl_names $mtl_lib]
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
    if { [llength $col_r] > 104 } {
        error "Lookup table can not have more than 104 colors."
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
proc draw_wavefront_obj {obj group_name scale_colors} {



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
    set group_names [lindex $obj 13]
    set mtl_names [lindex $obj 14]
	set mtl_lib [lindex $obj 15]

    # sanity checks:
    if { [llength vertex_idx_a] != [llength vertex_idx_b] || \
         [llength vertex_idx_a] != [llength vertex_idx_c] } {

        error "Enequal length vertex index lists."
    }
    if { [llength vertex_norm_idx_a] != [llength vertex_norm_idx_b] || \
         [llength vertex_norm_idx_a] != [llength vertex_norm_idx_c] } {

        error "Enequal length vertex index lists."
    }

	# change colours to match material definitions:
	set mtl_color_indeces [WOBJ::import_wavefront_mtl $mtl_lib]

    # prepare the color lookup table:
    set color_lookup_table [setup_color_lookup_table $scale_colors]
    
    # have vertex normals?
    set num_vert_idx [llength $vertex_idx_a]
    set num_norm_idx [llength $vertex_norm_idx_a]
    if { $num_vert_idx == $num_norm_idx } {
    
        # loop over faces and draw them:
        foreach ia $vertex_idx_a ib $vertex_idx_b ic $vertex_idx_c grp $group_names mtl $mtl_names {

			# matching group name?
			if { $grp != $group_name } {
				continue
			}

			# change color according to material:
			draw color [dict get $mtl_color_indeces $mtl]

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

   			# change color according to material:
			draw color [dict get $mtl_color_indeces $mtl]
            
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


# Returns RGB values created from LAB space interpolation of the "Reds" 
# color scale of RColorBrewer.
#
# INPUT:
# none
proc color_scale_reds {} {
    set r {\
        1 0.99992 0.9998 0.99965 0.99945 0.99923 0.99896 0.99866\
        0.99832 0.99794 0.99752 0.99706 0.99656 0.99601 0.99542 0.99479\
        0.99413 0.99345 0.99276 0.99207 0.99138 0.99072 0.99008 0.9895\
        0.98897 0.98852 0.98815 0.98789 0.98772 0.98761 0.98757 0.98758\
        0.98762 0.98769 0.98778 0.98788 0.98799 0.98809 0.98818 0.98826\
        0.98832 0.98835 0.98834 0.98828 0.98817 0.98798 0.98771 0.98735\
        0.98689 0.98632 0.98562 0.98479 0.9838 0.98263 0.98121 0.97948\
        0.97738 0.97484 0.97182 0.96824 0.96405 0.95921 0.95365 0.94732\
        0.94016 0.93214 0.92328 0.9137 0.90348 0.89273 0.88155 0.87003\
        0.85826 0.84634 0.83436 0.8224 0.81057 0.79894 0.7876 0.7765\
        0.76558 0.75476 0.74394 0.73307 0.72205 0.71082 0.69929 0.68739\
        0.67506 0.66221 0.64878 0.6347 0.61993 0.60444 0.58819 0.57115\
        0.5533 0.5346 0.51504 0.4946 0.47327 0.45105 0.42793 0.40392\
    }
    set g {\
        0.96078 0.95787 0.95428 0.95003 0.94515 0.93968 0.93362 0.92702\
        0.91988 0.91225 0.90415 0.89559 0.88662 0.87724 0.86749 0.85739\
        0.84696 0.83624 0.82525 0.814 0.80253 0.79086 0.77901 0.767\
        0.75485 0.74259 0.73024 0.7178 0.70531 0.69277 0.68019 0.66759\
        0.655 0.64241 0.62985 0.61733 0.60487 0.59248 0.58018 0.56799\
        0.55591 0.54391 0.53197 0.52005 0.50814 0.49621 0.48421 0.47213\
        0.45992 0.44755 0.43499 0.42219 0.40911 0.39572 0.38206 0.36813\
        0.35396 0.33958 0.32502 0.31032 0.29552 0.28069 0.26589 0.2512\
        0.23672 0.22258 0.20885 0.19556 0.18277 0.1705 0.1588 0.1477\
        0.13724 0.12746 0.11838 0.11005 0.10248 0.09569 0.0897 0.08448\
        0.08 0.07621 0.07303 0.0704 0.06822 0.06639 0.06482 0.0634\
        0.06203 0.0606 0.05903 0.05723 0.05509 0.05255 0.04952 0.04592\
        0.04164 0.03662 0.03114 0.02527 0.0191 0.01274 0.00633 0\
    }
    set b {\
        0.94118 0.93682 0.93151 0.92528 0.91818 0.91026 0.90154 0.89207\
        0.8819 0.87106 0.85959 0.84755 0.83496 0.82187 0.80832 0.79437\
        0.78005 0.76542 0.75053 0.73542 0.72015 0.70475 0.68928 0.67378\
        0.6583 0.64287 0.62756 0.61238 0.59734 0.58245 0.56771 0.55312\
        0.53868 0.5244 0.51028 0.49631 0.48251 0.46887 0.4554 0.44209\
        0.42894 0.41595 0.40314 0.39049 0.37801 0.36571 0.35358 0.34163\
        0.32986 0.31829 0.3069 0.29571 0.28473 0.27396 0.26343 0.25315\
        0.24314 0.23343 0.22402 0.21494 0.20621 0.19784 0.18986 0.18228\
        0.17512 0.1684 0.16211 0.15624 0.15074 0.14561 0.14082 0.13633\
        0.13213 0.12819 0.12448 0.12098 0.11766 0.11449 0.11147 0.10856\
        0.10576 0.10308 0.10049 0.09799 0.09558 0.09325 0.091 0.08881\
        0.08669 0.08462 0.0826 0.08062 0.07867 0.07673 0.0748 0.07286\
        0.07088 0.06887 0.06679 0.06461 0.0621 0.05907 0.05542 0.05098\
    }
    set color_rgb [list $r $g $b]
    return $color_rgb
}


# Returns RGB values created from LAB space interpolation of the "Greens" 
# color scale of RColorBrewer.
#
# INPUT:
# none
proc color_scale_greens {} {
    set r {\
        0.96863 0.96522 0.96144 0.95728 0.95277 0.9479 0.94268 0.93712\
        0.93123 0.92502 0.91849 0.91166 0.90452 0.8971 0.88938 0.8814\
        0.87314 0.86463 0.85586 0.84686 0.83762 0.82816 0.81848 0.80859\
        0.7985 0.78822 0.77776 0.76712 0.7563 0.74531 0.73416 0.72284\
        0.71135 0.69971 0.68791 0.67595 0.66384 0.65159 0.63919 0.62665\
        0.61397 0.60114 0.58815 0.57499 0.56167 0.54815 0.53444 0.52052\
        0.50639 0.49201 0.47738 0.46247 0.44726 0.43175 0.416 0.40006\
        0.38399 0.36784 0.35169 0.33561 0.31969 0.30404 0.28878 0.27403\
        0.25997 0.24675 0.23445 0.22303 0.21243 0.2026 0.19345 0.1849\
        0.17683 0.16913 0.16167 0.15431 0.14688 0.13923 0.13116 0.12262\
        0.11361 0.10412 0.09413 0.08361 0.07249 0.06068 0.04803 0.0344\
        0.0217 0.01057 0.00107 0 0 0 0 0\
        0 0 0 0 0 0 0 0\
    }
    set g {\
        0.98824 0.98692 0.98545 0.98385 0.98211 0.98023 0.97821 0.97606\
        0.97377 0.97136 0.96881 0.96614 0.96334 0.96041 0.95736 0.9542\
        0.95092 0.94753 0.94403 0.94043 0.93672 0.93293 0.92904 0.92507\
        0.92101 0.91688 0.91267 0.90838 0.90402 0.89956 0.89502 0.89037\
        0.88562 0.88075 0.87577 0.87066 0.86541 0.86003 0.85451 0.84883\
        0.84301 0.83705 0.83096 0.82475 0.81843 0.81202 0.80552 0.79894\
        0.79229 0.78558 0.77883 0.77204 0.76521 0.75835 0.75143 0.74443\
        0.73732 0.73007 0.72268 0.71511 0.70735 0.69937 0.69115 0.68268\
        0.67394 0.66491 0.6556 0.64607 0.63636 0.6265 0.61654 0.60652\
        0.59649 0.58648 0.57653 0.56669 0.55698 0.54745 0.53813 0.52899\
        0.52 0.51111 0.50228 0.49346 0.48462 0.4757 0.46668 0.45751\
        0.44814 0.43856 0.4287 0.41855 0.40808 0.39726 0.38608 0.37452\
        0.36257 0.35021 0.33742 0.32419 0.31051 0.29637 0.28176 0.26667\
    }
    set b {\
        0.96078 0.95652 0.95186 0.94681 0.94138 0.93558 0.92943 0.92293\
        0.91611 0.90896 0.9015 0.89374 0.8857 0.87738 0.86879 0.85996\
        0.85089 0.84159 0.83209 0.8224 0.81252 0.80248 0.79229 0.78196\
        0.7715 0.76093 0.75027 0.73951 0.72867 0.71774 0.70671 0.6956\
        0.6844 0.6731 0.66172 0.65024 0.63867 0.627 0.61524 0.60339\
        0.59146 0.57949 0.56754 0.55564 0.54383 0.53217 0.52068 0.50943\
        0.49844 0.48777 0.47745 0.46754 0.45807 0.44904 0.44042 0.43216\
        0.42422 0.41656 0.40914 0.40191 0.39484 0.38789 0.381 0.37415\
        0.36729 0.36038 0.3534 0.34637 0.33926 0.3321 0.32487 0.31758\
        0.31022 0.30281 0.29533 0.28779 0.28018 0.27252 0.26478 0.25699\
        0.24917 0.24134 0.2335 0.22569 0.21792 0.21022 0.2026 0.19508\
        0.1877 0.18047 0.17342 0.16656 0.15993 0.15355 0.14742 0.14156\
        0.136 0.13074 0.1258 0.12118 0.11689 0.11291 0.10925 0.10588\
    }
    set color_rgb [list $r $g $b]
    return $color_rgb
}


# Returns RGB values created from LAB space interpolation of the "Blues" 
# color scale of RColorBrewer.
#
# INPUT:
# none
proc color_scale_blues {} {
    set r {\
        0.96863 0.95928 0.95033 0.94175 0.9335 0.92557 0.91791 0.91051\
        0.90333 0.89635 0.88953 0.88285 0.87627 0.86978 0.86333 0.85689\
        0.85041 0.84385 0.83716 0.83031 0.82326 0.81595 0.80835 0.80041\
        0.79209 0.78334 0.77412 0.7644 0.75419 0.74351 0.73237 0.7208\
        0.7088 0.69639 0.68359 0.67041 0.65685 0.64294 0.62869 0.6141\
        0.59919 0.584 0.56858 0.55296 0.5372 0.52133 0.5054 0.48948\
        0.47362 0.45788 0.44235 0.4271 0.41222 0.39776 0.38372 0.37008\
        0.35684 0.34398 0.33148 0.31933 0.30751 0.29598 0.28472 0.27369\
        0.26285 0.25215 0.24157 0.2311 0.22072 0.21045 0.20027 0.19018\
        0.18017 0.17026 0.16045 0.15075 0.14117 0.13174 0.1225 0.11345\
        0.1046 0.09595 0.08752 0.07933 0.07138 0.06372 0.05639 0.04943\
        0.04291 0.03693 0.03193 0.02791 0.02482 0.02258 0.02115 0.02045\
        0.0204 0.02094 0.02197 0.02341 0.02517 0.02715 0.02925 0.03137\
    }
    set g {\
        0.98431 0.97951 0.97469 0.96987 0.96504 0.9602 0.95535 0.95048\
        0.9456 0.94071 0.9358 0.93087 0.92592 0.92095 0.91595 0.91095\
        0.90594 0.90094 0.89596 0.891 0.88607 0.88119 0.87635 0.87158\
        0.86687 0.86224 0.85769 0.85322 0.84878 0.84433 0.83981 0.8352\
        0.83044 0.82549 0.82031 0.81485 0.80908 0.80296 0.79645 0.7895\
        0.78212 0.77436 0.76626 0.75789 0.74928 0.7405 0.73158 0.72259\
        0.71355 0.70454 0.69558 0.68673 0.67802 0.66947 0.66104 0.6527\
        0.64441 0.63614 0.62785 0.61951 0.6111 0.60257 0.5939 0.58505\
        0.576 0.56672 0.5572 0.54749 0.53761 0.52759 0.51745 0.50723\
        0.49695 0.48664 0.47634 0.46605 0.45582 0.44566 0.4356 0.42565\
        0.41577 0.40597 0.39622 0.38652 0.37685 0.3672 0.35756 0.34792\
        0.33826 0.32858 0.31886 0.3091 0.2993 0.28944 0.27954 0.26959\
        0.25958 0.24953 0.23943 0.22928 0.21908 0.20884 0.19856 0.18824\
    }
    set b {\
        1 0.99683 0.99383 0.99098 0.98828 0.9857 0.98324 0.98089\
        0.97864 0.97647 0.97437 0.97233 0.97034 0.96839 0.96646 0.96453\
        0.96257 0.96057 0.95849 0.95632 0.95402 0.95158 0.94897 0.94616\
        0.94314 0.93988 0.93635 0.93255 0.92852 0.92431 0.91995 0.91549\
        0.91096 0.90642 0.90189 0.89742 0.89305 0.88882 0.88478 0.88094\
        0.87733 0.87391 0.87063 0.86746 0.86435 0.86126 0.85817 0.85502\
        0.85178 0.84842 0.84489 0.84116 0.83721 0.83301 0.8286 0.82401\
        0.81926 0.81437 0.80937 0.80428 0.79913 0.79394 0.78874 0.78355\
        0.77839 0.77329 0.76823 0.76321 0.75821 0.7532 0.74818 0.74313\
        0.73802 0.73286 0.72761 0.72226 0.7168 0.71122 0.7055 0.6996\
        0.69349 0.68714 0.6805 0.67355 0.66625 0.65855 0.65044 0.64187\
        0.6328 0.62321 0.61307 0.60235 0.591 0.57903 0.56638 0.55306\
        0.53903 0.52427 0.50877 0.49251 0.47548 0.45765 0.43903 0.41961\
    }
    set color_rgb [list $r $g $b]
    return $color_rgb
}

# Returns RGB values created from LAB space interpolation of the "Oranges" 
# color scale of RColorBrewer.
#
# INPUT:
# none
proc color_scale_oranges {} {
    set r {\
        1 0.99994 0.99982 0.99964 0.99941 0.99914 0.99882 0.99848\
        0.9981 0.99771 0.9973 0.99688 0.99645 0.99603 0.9956 0.99519\
        0.99479 0.9944 0.99404 0.99371 0.9934 0.99312 0.99287 0.99264\
        0.99245 0.99227 0.99212 0.99198 0.99186 0.99175 0.99164 0.99155\
        0.99147 0.99142 0.9914 0.99143 0.99151 0.99168 0.99194 0.99231\
        0.9928 0.99335 0.99393 0.9945 0.99499 0.99538 0.99562 0.99566\
        0.99545 0.99496 0.99413 0.99292 0.99128 0.98922 0.98675 0.98391\
        0.98073 0.97725 0.97349 0.96948 0.96526 0.96086 0.95632 0.95165\
        0.9469 0.94208 0.93717 0.93209 0.92674 0.92106 0.91496 0.90836\
        0.90118 0.89335 0.88478 0.87541 0.86516 0.85394 0.84171 0.82852\
        0.81448 0.79972 0.78436 0.7685 0.75226 0.73575 0.71908 0.70237\
        0.68571 0.66921 0.65299 0.63712 0.62169 0.60676 0.59239 0.57864\
        0.56557 0.55324 0.54171 0.53103 0.52128 0.51249 0.50472 0.49804\
    }
    set g {\
        0.96078 0.95659 0.95236 0.94809 0.94377 0.9394 0.93496 0.93044\
        0.92584 0.92115 0.91636 0.91146 0.90645 0.90131 0.89604 0.8906\
        0.88498 0.87916 0.8731 0.86679 0.8602 0.85331 0.84609 0.83852\
        0.83057 0.82222 0.81345 0.80425 0.79467 0.78477 0.77458 0.76416\
        0.75356 0.74283 0.73202 0.72118 0.71035 0.69959 0.68893 0.67844\
        0.66812 0.65795 0.6479 0.63794 0.62804 0.61817 0.60831 0.59842\
        0.58849 0.57847 0.56836 0.55812 0.54772 0.53717 0.52649 0.51568\
        0.50477 0.49378 0.48272 0.47162 0.46048 0.44933 0.43817 0.42703\
        0.41592 0.40486 0.39385 0.38292 0.37208 0.36137 0.3508 0.34042\
        0.33026 0.32035 0.31076 0.30153 0.29271 0.28436 0.27653 0.26921\
        0.26239 0.25602 0.25008 0.24452 0.23928 0.23434 0.22963 0.22512\
        0.22075 0.21648 0.21229 0.20812 0.20396 0.19978 0.19557 0.19131\
        0.18697 0.18255 0.17803 0.17337 0.16857 0.16359 0.15839 0.15294\
    }
    set b {\
        0.92157 0.91502 0.90806 0.9007 0.89294 0.88479 0.87626 0.86735\
        0.85807 0.84844 0.83845 0.82812 0.81745 0.80645 0.79512 0.78346\
        0.77147 0.75912 0.74642 0.73336 0.71994 0.70614 0.69196 0.67739\
        0.66243 0.64708 0.63132 0.61517 0.59869 0.58195 0.56502 0.54795\
        0.5308 0.51364 0.49652 0.47951 0.46266 0.44603 0.42967 0.41364\
        0.39796 0.38262 0.36758 0.35282 0.33832 0.32405 0.30999 0.29611\
        0.2824 0.26882 0.25535 0.24196 0.22863 0.21536 0.20215 0.18903\
        0.17602 0.16314 0.15041 0.13785 0.12548 0.11334 0.10145 0.08985\
        0.07862 0.06781 0.05748 0.04761 0.03821 0.03001 0.02331 0.01792\
        0.01365 0.01035 0.00786 0.00607 0.00483 0.00406 0.00364 0.00351\
        0.00364 0.00397 0.00449 0.00515 0.00592 0.00678 0.00771 0.00867\
        0.00966 0.01066 0.01164 0.0126 0.01352 0.01438 0.01516 0.01586\
        0.01646 0.01696 0.01733 0.01747 0.01737 0.01704 0.01648 0.01569\
    }
    set color_rgb [list $r $g $b]
    return $color_rgb
}


# Returns RGB values created from LAB space interpolation of the "Purples" 
# color scale of RColorBrewer.
#
# INPUT:
# none
proc color_scale_purples {} {
    set r {\
        0.98824 0.98521 0.98203 0.9787 0.97521 0.97157 0.96777 0.96382\
        0.9597 0.95542 0.95097 0.94636 0.94158 0.93663 0.9315 0.92619\
        0.92069 0.91498 0.90907 0.90293 0.89655 0.88994 0.88307 0.87594\
        0.86854 0.86086 0.85288 0.84461 0.83608 0.82732 0.81836 0.80923\
        0.79995 0.79057 0.78111 0.7716 0.76208 0.75258 0.74312 0.73375\
        0.72447 0.71527 0.70613 0.69705 0.68799 0.67896 0.66993 0.66089\
        0.65182 0.6427 0.63353 0.62427 0.61492 0.60548 0.59598 0.58644\
        0.57691 0.56741 0.55798 0.54866 0.53948 0.53048 0.5217 0.51318\
        0.50496 0.49707 0.48951 0.48225 0.47523 0.46843 0.4618 0.4553\
        0.44891 0.44257 0.43628 0.42998 0.42366 0.41729 0.41085 0.40434\
        0.39777 0.39114 0.38447 0.37775 0.371 0.36423 0.35744 0.35063\
        0.34383 0.33704 0.33026 0.32351 0.31679 0.31011 0.30347 0.29689\
        0.29036 0.28391 0.27752 0.27123 0.26502 0.25892 0.25293 0.24706\
    }
    set g {\
        0.98431 0.98025 0.97618 0.97211 0.96802 0.96391 0.95976 0.95556\
        0.9513 0.94698 0.94258 0.93809 0.93351 0.92882 0.92401 0.91908\
        0.91401 0.90878 0.9034 0.89783 0.89209 0.88615 0.88 0.87363\
        0.86703 0.86019 0.85311 0.84576 0.83815 0.83029 0.82216 0.81377\
        0.80512 0.79621 0.78704 0.7776 0.7679 0.75793 0.7477 0.73721\
        0.72647 0.71555 0.7045 0.69338 0.68224 0.67115 0.66015 0.64931\
        0.63868 0.62832 0.61827 0.6086 0.59935 0.59052 0.58199 0.57369\
        0.56553 0.55741 0.54924 0.54093 0.5324 0.52355 0.5143 0.50456\
        0.49424 0.48325 0.47161 0.45939 0.44666 0.43351 0.41999 0.40618\
        0.39215 0.37797 0.36371 0.34944 0.33522 0.32114 0.30725 0.29357\
        0.2801 0.26682 0.25374 0.24083 0.22809 0.21551 0.20307 0.19077\
        0.17858 0.16647 0.15444 0.14245 0.13048 0.11851 0.10652 0.09446\
        0.08228 0.0699 0.0572 0.04398 0.03054 0.01872 0.00858 0\
    }
    set b {\
        0.99216 0.98975 0.98736 0.98498 0.9826 0.98023 0.97784 0.97544\
        0.97302 0.97058 0.96811 0.9656 0.96305 0.96046 0.95781 0.95511\
        0.95235 0.94953 0.94664 0.94368 0.94065 0.93754 0.93434 0.93106\
        0.92769 0.92423 0.92067 0.917 0.91322 0.90931 0.90526 0.90106\
        0.89669 0.89215 0.88741 0.88248 0.87734 0.87197 0.86637 0.86052\
        0.85445 0.84819 0.8418 0.83534 0.82884 0.82237 0.81598 0.8097\
        0.8036 0.79772 0.79211 0.78682 0.7819 0.77734 0.77306 0.76901\
        0.76512 0.76131 0.75752 0.75368 0.74973 0.7456 0.74121 0.73651\
        0.73143 0.7259 0.71994 0.71359 0.70692 0.69997 0.6928 0.68546\
        0.67802 0.67051 0.66301 0.65555 0.64819 0.64099 0.63399 0.62719\
        0.62059 0.61416 0.60788 0.60175 0.59575 0.58985 0.58405 0.57833\
        0.57267 0.56706 0.56148 0.55592 0.55038 0.54486 0.53934 0.53384\
        0.52836 0.52288 0.51741 0.51195 0.5065 0.50106 0.49563 0.4902\
    }
    set color_rgb [list $r $g $b]
    return $color_rgb
}


# Returns RGB values created from LAB space interpolation of the "RdBu" 
# color scale of RColorBrewer. This is a divergent color scale.
#
# INPUT:
# none
proc color_scale_rdbu {} {
    set r {\
        0.69804 0.70822 0.71865 0.72929 0.7401 0.75107 0.76216 0.77334\
        0.78458 0.79585 0.80711 0.81834 0.82951 0.84059 0.85155 0.86234\
        0.87293 0.88328 0.89333 0.90305 0.91241 0.92135 0.92983 0.93782\
        0.94527 0.95213 0.95836 0.96394 0.96889 0.97324 0.97702 0.98026\
        0.98301 0.98531 0.9872 0.98872 0.98995 0.99093 0.99172 0.9924\
        0.99298 0.99338 0.99353 0.99336 0.9928 0.99177 0.9902 0.98802\
        0.98515 0.98151 0.97705 0.97168 0.96533 0.95799 0.94972 0.94057\
        0.93058 0.91981 0.90828 0.89605 0.88315 0.86961 0.85547 0.84075\
        0.82548 0.80966 0.79331 0.77643 0.75902 0.74107 0.72258 0.70356\
        0.684 0.66391 0.64328 0.62212 0.60042 0.57819 0.55543 0.5322\
        0.50857 0.48462 0.4604 0.436 0.41148 0.38691 0.36235 0.33787\
        0.31354 0.28945 0.26569 0.24236 0.21962 0.19763 0.17667 0.15709\
        0.13944 0.12443 0.11299 0.10607 0.10439 0.10813 0.11676 0.12941\
    }
    set g {\
        0.09412 0.11843 0.1415 0.16384 0.18573 0.20736 0.22884 0.25025\
        0.27165 0.29306 0.31451 0.33601 0.35757 0.37917 0.40082 0.42248\
        0.44412 0.46573 0.48725 0.50867 0.52993 0.55099 0.57183 0.59238\
        0.6126 0.63245 0.65187 0.67084 0.68935 0.70741 0.72503 0.7422\
        0.75893 0.77523 0.79109 0.80651 0.82151 0.83607 0.85021 0.86391\
        0.87716 0.88986 0.90193 0.91329 0.92385 0.93352 0.94222 0.94985\
        0.95632 0.96156 0.96546 0.96794 0.96893 0.96843 0.96658 0.96349\
        0.95929 0.95412 0.9481 0.94136 0.93402 0.92621 0.91805 0.90968\
        0.90121 0.89277 0.88435 0.8759 0.86734 0.85863 0.84969 0.84047\
        0.83089 0.8209 0.81045 0.79947 0.78791 0.77571 0.76282 0.74929\
        0.73518 0.72057 0.70554 0.69015 0.67449 0.65861 0.64259 0.6265\
        0.6104 0.59436 0.57845 0.56272 0.54721 0.53197 0.51702 0.50239\
        0.48811 0.47422 0.46073 0.44766 0.43504 0.42289 0.4112 0.4\
    }
    set b {\
        0.16863 0.17651 0.1848 0.1935 0.20261 0.21214 0.22208 0.23245\
        0.24325 0.25447 0.26611 0.27818 0.29068 0.3036 0.31695 0.33073\
        0.34497 0.35967 0.37485 0.39051 0.40668 0.42335 0.44055 0.45827\
        0.47653 0.49534 0.5147 0.53459 0.55494 0.57569 0.59675 0.61805\
        0.63952 0.66107 0.68264 0.70415 0.72552 0.74668 0.76754 0.78803\
        0.80806 0.82752 0.8463 0.86428 0.88135 0.89739 0.91229 0.92593\
        0.9382 0.94897 0.95814 0.96559 0.9712 0.97501 0.97716 0.97782\
        0.97714 0.97529 0.97243 0.96873 0.96435 0.95944 0.95417 0.94871\
        0.94322 0.93783 0.9326 0.92748 0.92242 0.91738 0.91231 0.90717\
        0.90191 0.89649 0.89086 0.88497 0.8788 0.87228 0.86538 0.85811\
        0.85052 0.84265 0.83453 0.82621 0.81773 0.80913 0.80045 0.79174\
        0.78302 0.77436 0.76577 0.75731 0.749 0.74087 0.73293 0.72522\
        0.71776 0.71057 0.70368 0.69711 0.69088 0.68502 0.67956 0.67451\
    }
    set color_rgb [list $r $g $b]
    return $color_rgb
}


# Returns RGB values created from LAB space interpolation of the "PuOr"
# color scale of RColorBrewer. This is a divergent color scale.
#
# INPUT:
# none
proc color_scale_puor {} {
    set r {\
        0.70196 0.7165 0.73095 0.74529 0.75953 0.77365 0.78763 0.80147\
        0.81512 0.82858 0.84182 0.85482 0.86754 0.87997 0.89206 0.90377\
        0.91506 0.92585 0.93612 0.9458 0.95484 0.9632 0.97083 0.97769\
        0.98371 0.98888 0.99313 0.99648 0.99899 1 1 1\
        1 1 1 1 0.99892 0.99776 0.99667 0.99577\
        0.99506 0.99446 0.99383 0.99309 0.99212 0.99083 0.98913 0.98693\
        0.98415 0.9807 0.9765 0.97148 0.96555 0.95873 0.95112 0.94283\
        0.93394 0.92455 0.91473 0.90456 0.89412 0.88349 0.87273 0.86191\
        0.8511 0.84035 0.82967 0.81899 0.80827 0.79745 0.78648 0.77531\
        0.76389 0.75217 0.74012 0.72767 0.7148 0.70145 0.68761 0.67329\
        0.65858 0.64353 0.62823 0.61273 0.5971 0.58141 0.5657 0.55003\
        0.53447 0.51906 0.50385 0.48887 0.47416 0.45972 0.44557 0.4317\
        0.41812 0.40481 0.39176 0.37896 0.36636 0.35394 0.34164 0.32941\
    }
    set g {\
        0.3451 0.35347 0.36277 0.37295 0.38396 0.39575 0.40826 0.42145\
        0.43526 0.44964 0.46454 0.47991 0.49569 0.51184 0.5283 0.54501\
        0.56189 0.57889 0.59592 0.61291 0.62981 0.64654 0.66304 0.67922\
        0.69503 0.71039 0.72523 0.73952 0.75328 0.76657 0.77942 0.79187\
        0.80397 0.81576 0.82729 0.83859 0.84972 0.86072 0.87163 0.88251\
        0.89331 0.90393 0.91421 0.92403 0.93324 0.9417 0.94927 0.95581\
        0.96117 0.96522 0.9678 0.96878 0.96802 0.96552 0.9614 0.95578\
        0.94881 0.94061 0.93131 0.92105 0.90995 0.89816 0.8858 0.87299\
        0.85988 0.84657 0.83312 0.81952 0.80577 0.79186 0.77779 0.76355\
        0.74914 0.73456 0.7198 0.70486 0.68975 0.67444 0.65896 0.64327\
        0.62738 0.61126 0.5949 0.57827 0.56138 0.5442 0.52671 0.5089\
        0.49075 0.47225 0.45337 0.43409 0.4144 0.39426 0.37366 0.35254\
        0.33086 0.30856 0.28556 0.26173 0.2369 0.21083 0.18308 0.15294\
    }
    set b {\
        0.02353 0.00614 0 0 0 0 0 0\
        0 0 0.00664 0.02932 0.05642 0.08148 0.10538 0.1288\
        0.1521 0.17544 0.19894 0.22264 0.24659 0.27078 0.29519 0.31981\
        0.34459 0.3695 0.39449 0.41951 0.44458 0.46969 0.49484 0.52004\
        0.54529 0.57059 0.59593 0.62132 0.64675 0.67223 0.69775 0.72331\
        0.74883 0.7741 0.79891 0.82303 0.84626 0.86836 0.8891 0.90826\
        0.9256 0.9409 0.95391 0.9644 0.97215 0.97718 0.97972 0.98002\
        0.97835 0.97495 0.9701 0.96404 0.95703 0.94933 0.94119 0.93286\
        0.9246 0.91664 0.90903 0.90169 0.89454 0.8875 0.8805 0.87345\
        0.86627 0.85889 0.85124 0.84323 0.83479 0.82585 0.81635 0.8063\
        0.79576 0.7848 0.77346 0.76181 0.7499 0.73779 0.72553 0.71317\
        0.70077 0.68838 0.67604 0.66382 0.65172 0.63978 0.62802 0.61646\
        0.60512 0.59401 0.58316 0.57258 0.56229 0.55231 0.54265 0.53333\
    }
    set color_rgb [list $r $g $b]
    return $color_rgb
}


# Returns RGB values created from LAB space interpolation of the "BrBG" 
# color scale of RColorBrewer. This is a divergent color scale.
#
# INPUT:
# none
proc color_scale_brbg {} {
    set r {\
        0.54902 0.56677 0.58415 0.60119 0.61788 0.63423 0.65022 0.66584\
        0.68107 0.69589 0.71028 0.72421 0.73766 0.75061 0.76304 0.77494\
        0.78632 0.79717 0.80752 0.81739 0.82679 0.83578 0.84439 0.8527\
        0.86076 0.86865 0.87646 0.88425 0.89201 0.89973 0.90737 0.91491\
        0.92232 0.92955 0.93656 0.94331 0.94974 0.9558 0.96143 0.96657\
        0.97116 0.97511 0.97834 0.98079 0.98238 0.98305 0.98273 0.98136\
        0.97889 0.97526 0.97042 0.96432 0.95692 0.94822 0.93832 0.92727\
        0.91515 0.90204 0.88798 0.87306 0.85733 0.84084 0.82367 0.80585\
        0.78744 0.76849 0.74902 0.72906 0.70863 0.68774 0.66642 0.6447\
        0.62261 0.60017 0.57743 0.55443 0.53121 0.50783 0.48434 0.46079\
        0.43725 0.41375 0.39034 0.36707 0.34396 0.32105 0.29839 0.27598\
        0.25385 0.23203 0.21051 0.18931 0.16843 0.14788 0.12765 0.10776\
        0.08818 0.06893 0.05002 0.03194 0.01831 0.00942 0.00478 0.00392\
    }
    set g {\
        0.31765 0.32559 0.3349 0.34554 0.35743 0.37052 0.38472 0.39995\
        0.41614 0.43319 0.45102 0.46955 0.48869 0.50837 0.52849 0.54894\
        0.56962 0.59042 0.61123 0.63193 0.65243 0.67261 0.69236 0.71159\
        0.73016 0.74798 0.76494 0.78096 0.79607 0.81031 0.82371 0.8363\
        0.84814 0.85924 0.86964 0.87939 0.88852 0.89706 0.90506 0.91256\
        0.91957 0.92608 0.93208 0.93755 0.94249 0.94686 0.95067 0.9539\
        0.95653 0.95854 0.95993 0.96066 0.96074 0.96016 0.95896 0.95716\
        0.95481 0.95192 0.94855 0.94471 0.94045 0.93579 0.93078 0.92544\
        0.91982 0.91395 0.9078 0.90133 0.89448 0.8872 0.87943 0.87112\
        0.86221 0.85268 0.84245 0.83149 0.81975 0.8072 0.79378 0.77956\
        0.76461 0.74904 0.73291 0.71633 0.69937 0.68213 0.66467 0.64708\
        0.62944 0.61184 0.59433 0.57701 0.55991 0.5431 0.5266 0.51048\
        0.49478 0.47954 0.46479 0.45059 0.43697 0.42397 0.41164 0.4\
    }
    set b {\
        0.03922 0.03216 0.02859 0.0286 0.03243 0.04049 0.0521 0.0661\
        0.08192 0.0992 0.11769 0.13725 0.15778 0.17919 0.20143 0.2244\
        0.248 0.27215 0.29674 0.32167 0.34684 0.37213 0.39743 0.42261\
        0.44754 0.4721 0.49615 0.51961 0.54252 0.56492 0.58687 0.60842\
        0.62961 0.6505 0.67113 0.69156 0.71184 0.73201 0.75213 0.77225\
        0.79233 0.81218 0.83164 0.8505 0.86859 0.88572 0.9017 0.91635\
        0.92948 0.9409 0.95044 0.95789 0.9631 0.96604 0.96691 0.96588\
        0.96314 0.95887 0.95326 0.94649 0.93875 0.93022 0.92108 0.91153\
        0.90174 0.89188 0.882 0.87206 0.86202 0.85183 0.84145 0.83085\
        0.81998 0.80879 0.79725 0.78533 0.77297 0.76015 0.74682 0.73302\
        0.71878 0.70416 0.6892 0.67394 0.65844 0.64274 0.62689 0.61092\
        0.59488 0.57883 0.56279 0.5468 0.53091 0.51513 0.49949 0.48401\
        0.46871 0.45362 0.43876 0.42415 0.40981 0.39577 0.38203 0.36863\
    }
    set color_rgb [list $r $g $b]
    return $color_rgb
}


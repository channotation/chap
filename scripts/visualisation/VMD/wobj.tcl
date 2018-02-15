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


namespace eval WOBJ {

# global variables in this namespace:
set BASENAME obj.obj
set FALLBACK_COLOR 8

# Parses Wavefront MTL file and builds a dictionary of material properties.
#
# INPUT:
# filename - name of MTL file to import
proc import_wmtl {filename} {

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
proc import_wobj {filename} {

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
            puts "loading OBJ group $crnt_group_name"
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
proc draw_wobj {obj group_name} {

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

    # does group name exist?
    set unique_group_names [lsort -unique $group_names]
    set use_group_colours [lsearch -exact $unique_group_names $group_name]

    # make fallback colour accessible:
    global FALLBACK_COLOR

    # change colours to match material definitions:
	set mtl_color_indeces [WOBJ::import_wmtl $mtl_lib]
    
    # have vertex normals?
    set num_vert_idx [llength $vertex_idx_a]
    set num_norm_idx [llength $vertex_norm_idx_a]
    if { $num_vert_idx == $num_norm_idx } {
    
        # loop over faces and draw them:
        foreach ia $vertex_idx_a ib $vertex_idx_b ic $vertex_idx_c grp $group_names mtl $mtl_names {

			# matching group name?
			if { $grp != $group_name && $use_group_colours >= 0 } {
				continue
			}

            # use color scale or solid color?
            if { $use_group_colours >= 0 } {
                # change color according to material:
			    draw color [dict get $mtl_color_indeces $mtl]
            } else {
                draw color $WOBJ::FALLBACK_COLOR 
            }

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
    draw_wobj $obj
}

}


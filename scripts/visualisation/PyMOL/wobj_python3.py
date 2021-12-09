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


from pymol import cmd
from pymol.cgo import *


class WavefrontMtlImporter:
    """ Reads material specifications from Wavefront MTL files.

    Currently only supports reading ambient (Ka), diffuse (Kd), and specular
    (Ks) colours from each material definition introduces with newmtl.
    """

    # dictionary of materials:
    __materials = dict()

    # most recently read material name:
    __crnt_material = ""


    def read(self, filename):

        # read data from MTL file:
        with open(filename, "r") as f:

            # file is parsed line by line:
            for line in f: self.__parse_line(line)

        # return materials dictionary:
        return self.__materials


    def __parse_line(self, line):
        """ Parses one line from an MTL file.

        Currently only supports newmtl, Ka, Kd, and Ks records. Depending on
        which of these records in encountered, the corresponding parse function
        is called to further process the line. If an unknown record is
        encountered, this method will throw an exception.

        Args:
            line: a string read from an MTL file
        """

        # split line at space:
        split_line = line.split(" ")

        # what type of record does the lien contain?
        if split_line[0] == "newmtl": self.__parse_newmtl_line(split_line)
        elif split_line[0] == "Ka": self.__parse_ka_line(split_line)
        elif split_line[0] == "Kd": self.__parse_kd_line(split_line)
        elif split_line[0] == "Ks": self.__parse_ks_line(split_line)
        else: raise Exception("Encountered unknown record type: " + line)


    def __parse_newmtl_line(self, split_line):
        """ Parses a newmtl record read from a MTL file.

        This will create a dictionary for the new material and also change the
        current matrial member variable accordingly.

        Args:
            split_line: list of strings, MTL line split at spaces
        """

        # change current material:
        self.__crnt_material = split_line[1].rstrip()

        # add current material to dictionary:
        self.__materials[self.__crnt_material] = dict()


    def __parse_ka_line(self, split_line):
        """ Parses an ambient colour record read from a MTL file.

        This will add an ambient colour definition to the current material.

        Args:
            split_line: list of strings, MTL line split at spaces
        """

        # add ambient colour value to materials dictionary:
        self.__materials[self.__crnt_material]["Ka"] = list(map(float, split_line[1:4]))


    def __parse_kd_line(self, split_line):
        """ Parses a diffuse colour record read from a MTL file.

        This will add an diffuse colour definition to the current material.

        Args:
            split_line: list of strings, MTL line split at spaces
        """

        # add ambient colour value to materials dictionary:
        self.__materials[self.__crnt_material]["Kd"] = list(map(float, split_line[1:4]))


    def __parse_ks_line(self, split_line):
        """ Parses a specular colour record read from a MTL file.

        This will add an specular colour definition to the current material.

        Args:
            split_line: list of strings, MTL line split at spaces
        """

        # add ambient colour value to materials dictionary:
        self.__materials[self.__crnt_material]["Ks"] = list(map(float, split_line[1:4]))



class WavefrontObjImporter:
    """ Reads scene data from Wavefront OBJ files.

    Currently supports reading of vertices and vertex normals, but not texture
    vertices. Expects every face to have a vertex and vertex normal index
    specification.
    """

    # dictionary of faces in groups:
    __groups = dict()

    # dictionary of materials read from MTL file:
    __materials = dict()

    # lists for vertices and vertex normals:
    __vertices = []
    __normals = []

    # keep track of object, group, and material associated with a face:
    __crnt_object = ""
    __crnt_group = ""
    __crnt_material = ""


    def read(self, filename):
        """ Reads data from a Wavefront OBJ file.

        Args:
            filename: name of the OBJ file to read

        Returns:
            dict: a dictionary of Compiled Graphics Objects
        """

        # read data from OBJ file:
        with open(filename, "r") as f:

            # file is parsed line by line:
            for line in f: self.__parse_line(line)

        # create a dictionary for the CGO objects to be returned:
        cgo = dict()

        # loop over groups:
        for key in self.__groups:

            # add faces to CGO object:
            cgo[key] = self.__create_cgo_object(key)

        # return dictionary of cgo objects:
        return cgo


    def __parse_line(self, line):
        """ Parses a single line read from an OBJ file.

        This will split the line string at spaces and check which type of
        record the line contains. It will then call the appropriate line
        parsing function.

        Args:
            line: string representing line from OBJ file
        """

        # split string at space:
        split_line = line.split(" ")

        # determine type of record:
        if split_line[0] == "#": return           # ignore comments
        elif split_line[0] == "mtllib": self.__parse_mtllib_line(split_line)
        elif split_line[0] == "o": self.__parse_o_line(split_line)
        elif split_line[0] == "g": self.__parse_g_line(split_line)
        elif split_line[0] == "usemtl": self.__parse_usemtl_line(split_line)
        elif split_line[0] == "v": self.__parse_v_line(split_line)
        elif split_line[0] == "vn": self.__parse_vn_line(split_line)
        elif split_line[0] == "f": self.__parse_f_line(split_line)
        elif split_line[0] == "\n": return        # ignore empty lines
        else: raise Exception("Encountered unknorn record type: " + line)


    def __parse_mtllib_line(self, split_line):
        """ Parses a mtllib record read from an OBJ file.

        This will use the WavefrontMtlImporter class to parse the MTL file
        specified in the record. The resulting materials are then added to an
        internal dictionary.

        Args:
            split_line: list of strings created by splitting line at spaces
        """

        mtl_importer = WavefrontMtlImporter()
        self.__materials = mtl_importer.read(split_line[1].rstrip())


    def __parse_o_line(self, split_line):
        """ Parses an o record read from an OBJ file.

        This will simply change the current object variable.

        Args:
            split_line: list of strings created by splitting line at spaces
        """

        # new current object:
        __crnt_object = split_line[1]


    def __parse_usemtl_line(self, split_line):
        """ Parses an usemtl record read from an OBJ file.

        This will change the current material to the value contained in the
        given record.

        Args:
            split_line: list of strings created by splitting line at spaces
        """

        # new current material:
        self.__crnt_material = split_line[1].rstrip()


    def __parse_g_line(self, split_line):
        """ Parses a g record read from an OBJ file.

        This changes the current group name and also creates a dictionary
        entry for this group name. Inside this dictionary, an empty list of
        faces is created.

        Args:
            split_line: list of strings created by splitting line at spaces
        """

        # new current group:
        self.__crnt_group = split_line[1].rstrip()

        # add group to group dictionary:
        self.__groups[self.__crnt_group] = dict()
        self.__groups[self.__crnt_group]["faces"] = []


    def __parse_v_line(self, split_line):
        """ Parses a v record read from an OBJ file.

        This adds a new 3-vector of floats to the internal list of vertices.

        Args:
            split_line: list of strings created by splitting line at spaces
        """

        # append vertex 3-vector to list of vertices:
        self.__vertices.append(list(map(float, split_line[1:4])))


    def __parse_vn_line(self, split_line):
        """ Parses a vn record read from an OBJ file.

        This adds a new 3-vector of floats to the internal list of vertex
        normals.

        Args:
            split_line: list of strings created by splitting line at spaces
        """

        # append vertex normal 3-vector to list of vertices:
        self.__normals.append(list(map(float, split_line[1:4])))


    def __parse_f_line(self, split_line):
        """ Parses a f record read from an OBJ file.

        This adds a new entry to the list of faces for the current group. A
        face entry consists of two 3-vectors of vertex indices and vertex
        normal indices as well as the name of the current material.

        Args:
            split_line: list of strings created by splitting line at spaces
        """

        # split at slash to separate out vertex, texture and normal index:
        all_indices = list(map(lambda x: x.split("/"), split_line[1:4]))

        # create a face with material and index into vertex and normal lists:
        face = dict()
        face["material"] = self.__crnt_material
        face["vert_idx"] = list(map(int, list(map(lambda x: x[0], all_indices))) )
        face["norm_idx"] = list(map(int, list(map(lambda x: x[2], all_indices))) )

        # append face to list of faces for the current group:
        self.__groups[self.__crnt_group]["faces"].append(face)


    def __create_cgo_object(self, groupname):
        """ Creates a Compiled Graphics Object from the fully parsed OBJ file.

        This function is called after all records have been read from an OBJ
        file and associated MTL file. It goes through the list of all faces
        for the given group and creates a CGO TRIANGLE entry for each.

        Args:
            groupname: string, must be a key in the groups dictionary

        Returns:
            CGO: a Compiled Graphics Object compatible with PyMOL
        """

        # create a CGO object:
        cgo_object = [ BEGIN, TRIANGLES ]

        # loop over faces in group:
        for face in self.__groups[groupname]["faces"]:

            # obtain vertices (note index shift wrt OBJ file):
            vertex_a = self.__vertices[face["vert_idx"][0] - 1]
            vertex_b = self.__vertices[face["vert_idx"][1] - 1]
            vertex_c = self.__vertices[face["vert_idx"][2] - 1]

            # obtain normals (note index shift wrt OBJ file):
            normal_a = self.__normals[face["norm_idx"][0] - 1]
            normal_b = self.__normals[face["norm_idx"][1] - 1]
            normal_c = self.__normals[face["norm_idx"][2] - 1]

            # obtain ambient colour from material:
            ka = self.__materials[face["material"]]["Ka"]

            # add triangular vertices to CGO object:
            cgo_object.extend([
                COLOR, ka[0], ka[1], ka[2],
                NORMAL, normal_a[0], normal_a[1], normal_a[2],
                VERTEX, vertex_a[0], vertex_a[1], vertex_a[2]])
            cgo_object.extend([
                COLOR, ka[0], ka[1], ka[2],
                NORMAL, normal_b[0], normal_b[1], normal_b[2],
                VERTEX, vertex_b[0], vertex_b[1], vertex_b[2]])
            cgo_object.extend([
                COLOR, ka[0], ka[1], ka[2],
                NORMAL, normal_c[0], normal_c[1], normal_c[2],
                VERTEX, vertex_c[0], vertex_c[1], vertex_c[2]])


        # end this CGO object:
        cgo_object.extend([END])

        # return CGO obtject:
        return cgo_object


def import_wobj(filename):
    """ Reads data from OBJ file and returns a dictionary of CGO objects.

    Args:
        filename: Name of the OBJ file to load.

    Returns:
        dict: Dictionary of CGO objects where the keys are OBJ group names.
    """

    # create OBJ importer:
    obj_importer = WavefrontObjImporter()

    # import data from OBJ file:
    return obj_importer.read(filename)


def draw_wobj(obj, groupname = None):
    """ Draws a named Compiled Graphics Object (CGO) from a given dictionary.

    Takes a dictionary of Compiled Graphics Objects and a string as input and
    uses PyMOL functions to draw the CGO associated with a given string. If
    the name is none (the default) all CGOs in the dictionary will be drawn.
    This can then be manually hidden or shown using the PyMOL GUI.

    Args:
        obj: Dictionary of names and Compiled Graphics Objects
        groupname: Name of the Compiled Graphics Object to draw
    """

    # delete existing objects under this name:
    cmd.delete(groupname)

    # has a group name been specified?
    if groupname != None:

        # load the requested CGO object:
        cmd.load_cgo(obj[groupname], groupname)

    else:

        # load all CGO objects in dictionary:
        for key in obj.keys(): cmd.load_cgo(obj[key], key)

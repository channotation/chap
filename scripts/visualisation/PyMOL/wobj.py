from pymol import cmd
from pymol.cgo import *



class WavefrontMtlImporter:

    # dictionary of materials:
    materials = dict()

    # most recently read material name:
    crnt_material = ""


    def read(self, filename):

        # read data from MTL file:
        with open(filename, "r") as f:

            # file is parsed line by line:
            for line in f: self.parse_line(line)
        
        # return materials dictionary:
        return self.materials


    def parse_line(self, line):

        # split line at space:
        split_line = line.split(" ")
    
        # what type of record does the lien contain?
        if split_line[0] == "newmtl": self.parse_newmtl_line(split_line)
        elif split_line[0] == "Ka": self.parse_ka_line(split_line)
        elif split_line[0] == "Kd": self.parse_kd_line(split_line)
        elif split_line[0] == "Ks": self.parse_ks_line(split_line)
        else: raise Exception("Encountered unknown record type: " + line)


    def parse_newmtl_line(self, split_line):

        # change current material:
        self.crnt_material = split_line[1].rstrip()

        # add current material to dictionary:
        self.materials[self.crnt_material] = dict()


    def parse_ka_line(self, split_line):

        # add ambient colour value to materials dictionary:
        self.materials[self.crnt_material]["Ka"] = map(float, split_line[1:4])


    def parse_kd_line(self, split_line):

        # add ambient colour value to materials dictionary:
        self.materials[self.crnt_material]["Kd"] = map(float, split_line[1:4])


    def parse_ks_line(self, split_line):

        # add ambient colour value to materials dictionary:
        self.materials[self.crnt_material]["Ks"] = map(float, split_line[1:4])




class WavefrontObjImporter:
   
    # dictionary of faces in groups:
    groups = dict()

    # dictionary of materials read from MTL file:
    materials = dict()

    # lists for vertices and vertex normals:
    vertices = []
    normals = []

    # faces represented as lists of vertex, texture, and vertex normal indices: 
    vertex_indices = []
    vertex_normal_indices = []


    # keep track of object, group, and material associated with a face:
    crnt_object = ""
    crnt_group = ""
    crnt_material = ""

    
    def read(self, filename):

        # read data from OBJ file:
        with open(filename, "r") as f:

            # file is parsed line by line:
            for line in f: self.parse_line(line)

        # create a dictionary for the CGO objects to be returned:
        cgo = dict()

        # loop over groups:
        for key in self.groups:
       
            # add faces to CGO object:
            cgo[key] = self.create_cgo_object(key)


        # return dictionary of cgo objects:
        return cgo





    
    
    
    def parse_line(self, line):

        # split string at space:
        split_line = line.split(" ")

        # determine type of record:
        if split_line[0] == "#": return           # ignore comments
        elif split_line[0] == "mtllib": self.parse_mtllib_line(split_line) 
        elif split_line[0] == "o": self.parse_o_line(split_line)
        elif split_line[0] == "g": self.parse_g_line(split_line)
        elif split_line[0] == "usemtl": self.parse_usemtl_line(split_line)
        elif split_line[0] == "v": self.parse_v_line(split_line)
        elif split_line[0] == "vn": self.parse_vn_line(split_line)
        elif split_line[0] == "f": self.parse_f_line(split_line)
        elif split_line[0] == "\n": return        # ignore empty lines
        else: raise Exception("Encountered unknorn record type: " + line)
        
    
    def parse_mtllib_line(self, split_line):

        mtl_importer = WavefrontMtlImporter()
        self.materials = mtl_importer.read(split_line[1].rstrip())

     
    def parse_o_line(self, split_line):
        
        # new current object:
        crnt_object = split_line[1]
   

    def parse_usemtl_line(self, split_line):
        
        # new current material:
        self.crnt_material = split_line[1].rstrip()
    

    def parse_g_line(self, split_line):
    
        # new current group:
        self.crnt_group = split_line[1].rstrip()
        
        # add group to group dictionary:
        self.groups[self.crnt_group] = dict()
        self.groups[self.crnt_group]["faces"] = []



    def parse_v_line(self, split_line):
    
        # append vertex 3-vector to list of vertices:
        self.vertices.append(map(float, split_line[1:4]))

    
    def parse_vn_line(self, split_line):
    
        # append vertex normal 3-vector to list of vertices:
        self.normals.append(map(float, split_line[1:4]))
    

    def parse_f_line(self, split_line):


        # split at slash to separate out vertex, texture and normal index: 
        all_indices = map(lambda x: x.split("/"), split_line[1:4])

        # create a face with material and index into vertex and normal lists:
        face = dict()
        face["material"] = self.crnt_material
        face["vert_idx"] = map(int, map(lambda x: x[0], all_indices)) 
        face["norm_idx"] = map(int, map(lambda x: x[2], all_indices)) 
 
        # append face to list of faces for the current group:
        self.groups[self.crnt_group]["faces"].append(face)


        # append each index type to individual list:
        # TODO remove
        self.vertex_indices.append(
            map(int, map(lambda x: x[0],
            all_indices)) )
        self.vertex_normal_indices.append(
            map(int, map(lambda x: x[2],
            all_indices)) )


    def create_cgo_object(self, groupname):

        # create a CGO object:
        cgo_object = [ BEGIN, TRIANGLES ]

        # loop over faces in group:
        for face in self.groups[groupname]["faces"]:

            # obtain vertices (note index shift wrt OBJ file):
            vertex_a = self.vertices[face["vert_idx"][0] - 1]
            vertex_b = self.vertices[face["vert_idx"][1] - 1]
            vertex_c = self.vertices[face["vert_idx"][2] - 1]
            
            # obtain normals (note index shift wrt OBJ file):
            normal_a = self.normals[face["norm_idx"][0] - 1]
            normal_b = self.normals[face["norm_idx"][1] - 1]
            normal_c = self.normals[face["norm_idx"][2] - 1]

            # obtain ambient colour from material:
            ka = self.materials[face["material"]]["Ka"]

            # add triangular vertices to CGO object:
            cgo_object.extend([
                COLOR, ka[0], ka[1], ka[2],
                VERTEX, vertex_a[0], vertex_a[1], vertex_a[2],
                NORMAL, normal_a[0], normal_a[1], normal_a[2]])
            cgo_object.extend([
                COLOR, ka[0], ka[1], ka[2],
                VERTEX, vertex_b[0], vertex_b[1], vertex_b[2],
                NORMAL, normal_b[0], normal_b[1], normal_b[2]])
            cgo_object.extend([
                COLOR, ka[0], ka[1], ka[2],
                VERTEX, vertex_c[0], vertex_c[1], vertex_c[2],
                NORMAL, normal_c[0], normal_c[1], normal_c[2]])


        # end this CGO object:
        cgo_object.extend([END])

        # return CGO obtject:
        return cgo_object
 
        


importer = WavefrontObjImporter()

cmd.hide("all")


cgo = importer.read("output.obj")

print cgo.keys()

cmd.delete("test")
cmd.load_cgo(cgo["avg_pl_hydrophobicity"], "test")




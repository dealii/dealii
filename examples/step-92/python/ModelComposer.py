#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""
ModelComposer  --- version 1.0, October 2024 ---

(C) 2023-2024, Stephan Voss, stvoss@gmx.de, Neunkirchen am Brand, Germany

 * ---------------------------------------------------------------------
 *
 * This is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found at
 * https://www.gnu.org/licenses/old-licenses/lgpl-2.1
 * https://www.gnu.org/licenses/licenses/lgpl-3.0
 *
 * ---------------------------------------------------------------------


The "MeshComposer_pkl.py" is designed for
 1) converting step-files into serialized (pickle-) objects and store it in temporary files
 2) import the serialized (pickle-) objects into model components of class type "Component"
 3) manage a set of components in a "Model" class object
 
The "Model" class object is finally used to merge all the components correctly into a single model.
Components can be cut out of others or added.
     
Finally, the merged objects can be meshed. "Model" objects function "mesh()" utilizes gmsh that task. 


Actually, an "occ_converter" object is needed for each component of a model and
a step file may only represent a single component.

NOTE:
    This is actually a workaround and shall be further developed in the future such that step-files containing 
    complete assemblies and component-IDs should be handled at once 
    instead of having to create separate step-files for each model component
    (as it is in this example and the used classes in it's actual version in October 2024) 

"""
import math             # ... basic math ...obvious ;-)
import os               # used for system and file system operations
import gmsh             # gmsh is utilized for meshing
import occ_converter    # importing, converting and exporting step-files into serialized (pickle-) objects and/or gmsh .geo files


try:
    from multiprocessing import cpu_count
    
except:
    exit()

"""
Component class

used for
 1) converting step-files into serialized (pickle-) objects and store it in temporary files
 2) import the serialized (pickle-) objects into model components of class type "Component"
 3) get geometry description and provide interface to "Model" object for retrieving that geometry and joining it into a single model for meshing
 
 Component objects are organized in a tree.
 They can manage a list of objects and the objects in that list refer to this object as its parent. 
"""
class Component():
    
    parent = None           # parent Component of this object
    MeshSize = 1            # element size for the mesh (of this component)
    
    Name = ''               # clear taxt Component name
    PhysicalID = None       # (2, ID) is a (2D) surface, (3, ID) or simply ID is representing a (3D) volume
    entities = []           # list of entities (gmsh)
    filename = ''           # file name (path) for this Components model file (step)
    action = ''             # 

    MeshSizeList = []
    
    verbosity = -1
    
    """
    Dictionary of "child" objects in the sub-tree structure of this object:
    """
    objects = {}


    # object initialization
    def __init__(self,Name,filename='',PhysicalID=None,MeshSize=1,action='add',_verbosity=-1):
        self.clear()

        self.Name = Name
        self.PhysicalID = PhysicalID if hasattr(PhysicalID, "__len__") else None if PhysicalID==None else (3,PhysicalID)
        self.MeshSize=MeshSize
        self.filename = filename
        self.action = action.lower()
        self.verbosity=_verbosity
        
        
    # reset/clear this object
    def clear(self):
        self.objects = {}
        self.entities = []
        self.MeshSizeList = []

    
    # remove a child-Component    
    def remove(self,objectkey):
        
        if self.objects.get(objectkey,None)!=None:
            self.objects.pop(objectkey)
        
    
    # return the root element of the Component object tree (the root has no parent)
    def get_root(self):
        l_object = self
        while l_object.parent!=None:
            l_object = l_object.parent
            
        return l_object
    
    
    # return the gmsh .geo file object where this Components contributes to
    def get_GeoFile(self):
        l_root = self.get_root()
        if l_root!=None:
            if hasattr(l_root,'GeoFile'):
                return l_root.GeoFile
            
        return None

    
    """
    add another Component "object" into this objects child list "self.objects"
    
    Consider that element to be added (merged) into the volume of this object
    """
    def add(self,object):
        self.objects[object.Name] = object
        object.parent = self 


    """
    add another Component "object" into this objects child list "self.objects"
    
    ....BUT - NOTE:
    
    Consider that element to be cut out from the volume of this object
    """
    def cut(self,object):
        object.action = 'cut'
        self.add(object)


    """
    set verbosity level - in regard to logging debug-information
    """    
    def set_verbosity(self,verbosity):
        self.verbosity = verbosity
        for key, object in self.objects.items():
            object.set_verbosity(verbosity)
            

    """
    serialize the model specified by step-file "_filename"
    and export a pickle-object into a file (self.pickle_filename)
    
    if _force == False: only convert if pickle-file does not exist
    if _force == True: force conversion and generation of a pickle-file
    """
    def step2pickle(self,_filename,_force=False):
        self.pickle_filename, _ = os.path.splitext(_filename)
        self.pickle_filename += '.pkl'

        if _force or not os.path.exists(self.pickle_filename):
            occ_conv_tmp = occ_converter.occ_converter()
            occ_conv_tmp.occ_import(_filename)
            ofile = open(self.pickle_filename,'wb')
            occ_conv_tmp.pickle_export(_output_file=ofile)
            ofile.close()
            del occ_conv_tmp
            
        self.occ_conv = occ_converter.occ_converter()

        #self.occ_conv.occ_import(_filename)

        #self.occ_conv.offset_all_entities(_offsets=[1000,2000,3000,4000])
        infile = open(self.pickle_filename,'rb')
        self.occ_conv.pickle_import(infile)
        if self.PhysicalID!=None:
            try:
                _ = iter(self.PhysicalID)
                self.occ_conv.physical = (self.PhysicalID[0],self.Name,self.PhysicalID[1])
                #print('Moin 1',self.occ_conv.physical)
            except:
                self.occ_conv.physical = (3,self.Name,self.PhysicalID)
                #print('Moin 2',self.occ_conv.physical)
            
        infile.close()
        

    """
    recursively import all models of the entire tree (see dict "objects" of child-elements)
    and serialize them all
    """
    def import_occ(self):
        if self.filename!=None and self.filename.__len__()>0:
            if self.verbosity>10:
                print('model \''+self.Name+'\' importing \''+self.filename+'\'')
                
            self.step2pickle(self.filename) # import geometry
            #self.entities = gmsh.model.occ.importShapes(self.filename,format=self.action)
            
        for key,item in self.objects.items():
            item.import_occ()

    
    """
    find a child-object by its clear-text Name
    """
    def find_object(self,Name):
        if Name==self.Name:
            return self
        
        for key,item in self.objects.items():
            object = item.find(Name)
            if object!=None:
                return object
        
        return None
    
    
    """
    recursively collect all entity IDs and all Component objects of all child elements and return them in two lists
    """
    def get_entities_recursively(self,action='add'):
        l_action=action.lower()
        ret_entities = self.entities[:]
        ret_objects = [self] * self.entities.__len__()
        for key, object in self.objects.items():
            if object.action!=l_action:
                continue
            
            tmp_entities, tmp_objects = object.get_entities_recursively()
            for ix in range(tmp_entities.__len__()):
                ret_entities.append(tmp_entities[ix])
                ret_objects.append(tmp_objects[ix])
                
        return ret_entities, ret_objects
    
    
    """
    recursively do fragmentation of this (self) object and its sub-tree Components (self.objects)
    and write the gmsh .geo file 
    """
    def fragment(self,_indent=0):
        self.serialize_entities()
        self.fragment_entities(action='add',_indent=_indent)
        
        l_geofile = self.get_GeoFile()
        l_groups = self.get_geo_file_levels()
        for group in l_groups:
            for line in group:
                print(line,file=l_geofile)
                

    """
    recursively serialize all entity IDs in order to be unique numbers
    
    This may be not the case after importing them using import_occ() and step2pickle() 
    """
    def serialize_entities(self,_offsets=[0,0,0,0]):
        
        l_offsets = _offsets
        
        if hasattr(self,'occ_conv'):

            self.occ_conv.offset_all_entities(_offsets=l_offsets)
            l_offsets = [x for x in self.occ_conv.calc_entities_min_max()]

        for _, object in self.objects.items():
            l_offsets = object.serialize_entities(_offsets=l_offsets)

        return l_offsets
    
    
    """
    recursively get the gmsh .geo file entities grouped by its levels (0...6 
    = Points, Lines, Curves, Surface Loops, ...)
    """
    def get_geo_file_levels(self):
        
        l_groups = [[] for _ in range(7)]
        
        if hasattr(self,'geo_lines'):
            for group_ix in range(self.geo_lines.__len__()):
                l_groups[group_ix].extend(self.geo_lines[group_ix])

        for _, object in self.objects.items():
            tmp_groups = object.get_geo_file_levels()
            for group_ix in range(tmp_groups.__len__()):
                l_groups[group_ix].extend(tmp_groups[group_ix])

        if hasattr(self,'geo_physicals'):
            l_groups[6].extend(self.geo_physicals)
            
        return l_groups
        
        
    """
    recursively do fragmentation of this (self) objects entities
    and its sub-tree Components (self.objects) entities
    """
    def fragment_entities(self,action='add',_indent=0,_export=True):
        l_action=self.action.lower()
        #print('fragment '+' '*_indent+' '+self.action+' \''+ self.Name+'\'')

        l_surfaces = []
        l_interfaces = []
        l_add = []
        l_cut = []

        for _, object in self.objects.items():
            
            object.fragment_entities(action=l_action,_indent = _indent+4,_export=False)
            if hasattr(self,'occ_conv') and hasattr(object,'occ_conv'):
                
                if object.action.lower()=='cut':
                    self.occ_conv.merge_entities(object.occ_conv)
                    
                    tmp_occ_conv = occ_converter.occ_converter()
                    tmp_occ_conv.merge(self.occ_conv)
                    tmp_occ_conv.merge(object.occ_conv)
                    tmp_surfaces, tmp_interfaces = tmp_occ_conv.get_ext_surface_and_interface_curves()

                    tmp_ents = []
                    for key,surface in object.occ_conv.surfacedict.items():
                        tmp_ents.extend([int(math.fabs(x)) for x in object.occ_conv.surfacedict[key]])
                        
                    l_surfaces.extend(tmp_ents)
                    if tmp_interfaces.__len__()>0 and self.action.lower()=='add':
                        l_add.extend(tmp_ents)
                        l_cut.extend(tmp_interfaces)
                    else:
                        l_cut.extend(tmp_interfaces)
                    
                    l_interfaces.extend(tmp_interfaces)
                        
                    object.geo_physicals = object.occ_conv.get_geo_lines(_physical_only=True)[-1]

                    object.occ_conv.surfacedict.clear()
                else:
                    for key,_ in object.occ_conv.surfacedict.items():
                        l_surfaces.extend([int(math.fabs(x)) for x in object.occ_conv.surfacedict[key]])

            
            object.geo_lines = object.occ_conv.get_geo_lines()
            if hasattr(object,'add_cut_in_surface_loop'):
                tmp_add, tmp_cut = object.add_cut_in_surface_loop
                    
                tmp_occ_conv = occ_converter.occ_converter()
                tmp_occ_conv.merge(object.occ_conv)
                for skey,sval in tmp_occ_conv.surfacedict.items():
                    tmp_occ_conv.surfacedict[skey].extend(tmp_add)
                    tmp_occ_conv.surfacedict[skey] = list(set(tmp_occ_conv.surfacedict[skey]))

                for skey,sval in tmp_occ_conv.surfacedict.items():
                    tmp_occ_conv.surfacedict[skey] = [x for x in list(set(sval)) if int(math.fabs(x)) not in tmp_cut]
                     
                tmp_geo_lines = tmp_occ_conv.get_geo_lines()
                object.geo_lines = [object.geo_lines[ix] if ix!=4 else tmp_geo_lines[ix] for ix in range(object.geo_lines.__len__())]

            
        if hasattr(self,'occ_conv'):
            l_cut = list(set(l_cut))
            l_add = [x for x in list(set(l_add)) if int(math.fabs(x)) not in l_cut]
            l_surfaces.extend(l_add)
            l_surfaces = [x for x in list(set(l_surfaces)) if int(math.fabs(x)) not in l_cut]
            if l_add.__len__()>0 or l_cut.__len__()>0:
                self.add_cut_in_surface_loop = [l_cut,l_add]
            
            tmp_surfaces = {}
            for skey,sval in self.occ_conv.surfacedict.items():
                sval.extend(l_surfaces)
                tmp_surfaces[skey] =  [x for x in list(set(sval)) if int(math.fabs(x)) not in l_cut]
                self.occ_conv.surfacedict[skey] = tmp_surfaces[skey]



    """
    callback-function for gmsh mesher for returning local mesh element size 
    """
    def MeshSizeCallback(self, entity_dim, entity_tag,x, y, z, lc):
        
        ret_val = None
        
        for key, object in self.objects.items():
            tmpsize = object.MeshSizeCallback(entity_dim, entity_tag, x,y,z,lc)
            if tmpsize!=None:
                ret_val = tmpsize
                break

        if (entity_dim,entity_tag) in self.MeshSizeList:
            if ret_val==None or (self.MeshSize!=None and self.MeshSize<ret_val):
                ret_val = self.MeshSize
                        
        return ret_val


"""
Model class

A "Model" class object shall be the root of a Component object tree.
This root is finally used to merge (and fragment) all the components correctly into a single model.
Components can be cut out of others or added.
"""
class Model(Component):
    
    MeshPartitions = 1          # number of mesh partitions
    MeshFilename = ''           # .msh file name for final mesh output file
    SubdivisionAlgorithm = 2    # type of mesh - actually "all hexahedra"
    supported_SubdivisionAlgorithms = {'none' : 0,
                             'quadrangles' : 1,
                             'all quadrangles' : 1,
                             'tetrahedra' : 1,
                             'hexahedra' : 2,
                             'all hexahedra' : 2,
                             'barycentric': 3, }
    
    n_threads_per_worker = None
    MeshScalingFactor = 1.0
    
    """
    Model initialization
    """
    def __init__(self,
                 Name,                      # clear text name of this model
                 MeshFilename,              # .msh file name for final mesh output file
                 MeshSize=1,                # standard/default mesh element size if not defined for each sub-Component
                 MeshPartitions=1,          # number of mesh partitions
                 _MeshScalingFactor=1.0,    # scaling factor for coordinates (e.g. scale from m to mm)
                 MeshSubdivisionAlgorithm='all hexahedra',  # mesh element type
                 _n_threads_per_worker=None):   # number of threads for meshing process
        
        super(Model,self).__init__(Name,MeshSize=MeshSize)
        
        self.MeshFilename = MeshFilename
        self.GeoFileName, _ = os.path.splitext(MeshFilename)
        self.GeoFileName += '.geo'
        self.n_threads_per_worker = _n_threads_per_worker
        self.MeshScalingFactor=_MeshScalingFactor
        
        self.set_N_partitions(MeshPartitions)
        self.SubdivisionAlgorithm = self.supported_SubdivisionAlgorithms.get(MeshSubdivisionAlgorithm,self.SubdivisionAlgorithm)
        
      
    # set number of mesh partitions to generate  
    def set_N_partitions(self,MeshPartitions):
        self.MeshPartitions = MeshPartitions
        
    
    """
    callback-function for gmsh mesher for returning local mesh element size 
    """
    def MeshSizeCallback(self, entity_dim, entity_tag,x, y, z, lc):
        tmp_size = super(Model,self).MeshSizeCallback(entity_dim,entity_tag,x,y,z,lc);
        if tmp_size!=None:
            return tmp_size    
        
        return self.MeshSize
        
    
    """
    export the complete .geo file for this model and its sub-Component tree
     (as needed to import in gmsh)
    """
    def GeoExport(self):
        self.import_occ()
        print(self.GeoFileName)
        self.GeoFile = open(self.GeoFileName,'w')
        self.fragment()
        self.GeoFile.close()
    
    
    """
    set options for gmsh
    """
    def mesh_set_options(self):
        #gmsh.option.setNumber("General.Terminal", 0)
        
        l_cpu_count = cpu_count()
        l_OCCparallel = 0
        
        l_thread_count = l_cpu_count if self.n_threads_per_worker==None else self.n_threads_per_worker
        l_OCCparallel = 1 if l_thread_count == l_cpu_count else 0
        
        algorithms = [ [6, 3 ], [6,1]]
        #print(gmsh.option.getString('General.BuildInfo'))
        #print(gmsh.option.getString('General.BuildOptions'))
        gmsh.option.setNumber("General.InitialModule", 0) # automatic mode
        gmsh.option.setNumber("General.NumThreads",l_thread_count)#max([int(0.5+0.5*cpu_count()),cpu_count()-4]))
        gmsh.option.setNumber("Geometry.OCCParallel",l_OCCparallel)
        #gmsh.option.setNumber("Geometry.OCCFastUnbind",2)
        #gmsh.option.setNumber("Geometry.OCCBooleanSimplify",2)
        #gmsh.option.setNumber("Mesh.MaxNumThreads1D",0)
        #gmsh.option.setNumber("Mesh.MaxNumThreads2D",0)
        #gmsh.option.setNumber("Mesh.MaxNumThreads3D",0)
        #gmsh.option.setNumber("Mesh.MaxNumThreads1D",32)#max([int(0.5+0.5*cpu_count()),cpu_count()-4]))
        #gmsh.option.setNumber("Mesh.MaxNumThreads2D",32)#max([int(0.5+0.5*cpu_count()),cpu_count()-4]))
        #gmsh.option.setNumber("Mesh.MaxNumThreads3D",32)#max([int(0.5+0.5*cpu_count()),cpu_count()-4]))
        gmsh.option.setNumber("General.Verbosity", self.verbosity) # 0: silent except for fatal errors, 1: +errors, 2: +warnings, 3: +direct, 4: +information, 5: +status, 99: +debug
        gmsh.option.setNumber("General.Axes", 3) # 0: none, 1: simple axes, 2: box, 3: full grid, 4: open grid, 5: ruler

        gmsh.model.add(self.Name)
        gmsh.logger.start()
    
        #gmsh.option.setNumber("Mesh.MeshSizeMin", 0.2*self.MeshScalingFactor)
        #gmsh.option.setNumber("Mesh.MeshSizeMax", 13*self.MeshScalingFactor)
        #gmsh.option.setNumber("Mesh.Format",1)
        #gmsh.option.setNumber("Mesh.MshFileVersion",4.1)
        gmsh.option.setNumber("Mesh.Binary", 1 )
    
        gmsh.option.setNumber("Mesh.MaxRetries",100)
        #gmsh.option.setNumber("Geometry.Tolerance",1e-12)
        #gmsh.option.setNumber("Geometry.MatchMeshTolerance",1e-12)

        
        # mesh.Algorithm in [3, 5, 6, 8]
        # mesh.Algorithm3D in [1, 3, 10]
        #gmsh.option.setNumber("Mesh.Algorithm", 3) # 3: initial mesh 6: Frontal-Delaunay, 8: frontal Delaunay for Quads, 9: Packing of Parallelograms
        #gmsh.option.setNumber("Mesh.Algorithm3D", 3 ) # 1: delauney, 3: initial mesh, 9: R-tree, 11: Quasi-structured
        #gmsh.option.setNumber("Mesh.Algorithm", 9)#6) # 3: initial mesh 6: Frontal-Delaunay, 8: frontal Delaunay for Quads, 9: Packing of Parallelograms
        #gmsh.option.setNumber("Mesh.Algorithm3D", 9)#1 ) # 1: delauney, 3: initial mesh, 9: R-tree, 11: Quasi-structured
    
        if self.SubdivisionAlgorithm>1:
            gmsh.option.setNumber("Mesh.Algorithm", 8) # 1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 
            #                                            # 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms, 11: Quasi-structured Quad
            gmsh.option.setNumber("Mesh.Algorithm3D", 10 ) # 1: Delaunay, 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT

        #gmsh.option.setNumber("Mesh.Algorithm", 3) # 1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 
        #gmsh.option.setNumber("Mesh.Algorithm3D", 3 ) # 1: Delaunay, 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT

        #gmsh.option.setNumber("Mesh.CompoundClassify",0)
        #gmsh.option.setNumber("Mesh.MeshScalingFactor", 1 )
        #gmsh.option.setNumber("Mesh.Hexahedra", 1 )
        #gmsh.option.setNumber("Mesh.HighOrderDistCAD", 0 )
        #gmsh.option.setNumber("Mesh.HighOrderIterMax", 100 )
        #gmsh.option.setNumber("Mesh.HighOrderNumLayers", 6 )
        #gmsh.option.setNumber("Mesh.HighOrderOptimize", 1)
        #gmsh.option.setNumber("Mesh.HighOrderPassMax", 25 )
        #gmsh.option.setNumber("Mesh.HighOrderPeriodic", 0 )
        #gmsh.option.setNumber("Mesh.HighOrderPoissonRatio", 0.33 )
        #gmsh.option.setNumber("Mesh.HighOrderSavePeriodic", 0 )
        #gmsh.option.setNumber("Mesh.HighOrderPrimSurfMesh", 0 )
        #gmsh.option.setNumber("Mesh.HighOrderThresholdMin", 0.1 )
        #gmsh.option.setNumber("Mesh.HighOrderThresholdMax", 2 )
        #gmsh.option.setNumber("Mesh.HighOrderFastCurvingNewAlgo", 1 )
        #gmsh.option.setNumber("Mesh.HighOrderCurveOuterBL", 0 )
        #gmsh.option.setNumber("Mesh.HighOrderMaxRho", 0.3 )
        #gmsh.option.setNumber("Mesh.HighOrderMaxAngle", 0.1745329277777778 )
        #gmsh.option.setNumber("Mesh.HighOrderMaxInnerAngle", 0.5235987833333333 )
        #gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1 )
        gmsh.option.setNumber("Mesh.Optimize", 1 )
        #gmsh.option.setNumber("Geometry.ScalingFactor", self.MeshScalingFactor )
        gmsh.option.setNumber("Mesh.ScalingFactor", self.MeshScalingFactor )
        #gmsh.option.setNumber("Mesh.HighOrderThresholdMin",0.1*self.MeshScalingFactor )
        #gmsh.option.setNumber("Mesh.HighOrderThresholdMax",2*self.MeshScalingFactor )
        #gmsh.option.setNumber("Mesh.HighOrderMaxRho",0.3*self.MeshScalingFactor )
        #gmsh.option.setNumber("Mesh.ToleranceInitialDelaunay",1e-12*self.MeshScalingFactor )
        #gmsh.option.setNumber("Mesh.ToleranceReferenceElement",1e-6*self.MeshScalingFactor )
        #gmsh.option.setNumber("Mesh.AnisoMax",1e+33 )
        #gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0 )
        #gmsh.option.setNumber("Mesh.MeshSizeFromCurvatureIsotropic", 0 )
        #gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1 )
        #gmsh.option.setNumber("Mesh.MeshSizeFromParametricPoints", 1 )
        #gmsh.option.setNumber("Mesh.RandomSeed", 1 )
        #gmsh.option.setNumber("Mesh.RecombineAll", 1 )
        #gmsh.option.setNumber("Mesh.AngleToleranceFacetOverlap",0.05)
        #gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1 ) # 0: simple, 1: blossom, 2: simple full-quad, 3: blossom full quad
        #gmsh.option.setNumber("Mesh.ReadGroupsOfElements", 1 )
        #gmsh.option.setNumber("Mesh.RecombineOptimizeTopology", 5 )
        #gmsh.option.setNumber("Mesh.Recombine3DAll", 1 )
        #gmsh.option.setNumber("Mesh.Recombine3DLevel", 0 )
        gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", self.SubdivisionAlgorithm ) # 2 -> all hexahedra
        
        #gmsh.option.setNumber("Mesh.RecombineAll", 1)
        #gmsh.option.setNumber("Mesh.CharacteristicLengthFactor", 1.455)#0.250)#0.305)
        #gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 2) # 2 <- deal.ii needs all hexahedra
        gmsh.option.setNumber("Mesh.Smoothing", 10 )
        
        #gmsh.option.setString("Geometry.OCCTargetUnit","M")
        
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature",0)#12) 
        gmsh.option.setNumber("Mesh.MinimumCircleNodes",12)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvatureIsotropic",0)
        gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 0)#1 )
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)#2)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        
        gmsh.option.setNumber("Mesh.RefineSteps",0) # for MeshAdapt algorithm only
        
        #gmsh.option.setNumber("Geometry.OCCImportLabels", 1) # import colors from STEP
        #gmsh.option.setNumber("Geometry.Color,Points", _pColor)
        #gmsh.option.setNumber("Geometry.Color,Curves", _cColor)
        #gmsh.option.setNumber("Geometry.Color,Surfaces", _sColor)
        #gmsh.option.setNumber("Geometry.Color,Volumes", _vColor)

        #gmsh.option.setNumber("Mesh.SaveAll",0)
        
        if self.MeshPartitions>1:
            gmsh.option.setNumber("Mesh.PartitionCreateTopology",1)
            gmsh.option.setNumber("Mesh.PartitionCreatePhysicals",1)
            gmsh.option.setNumber("Mesh.PartitionCreateGhostCells",1)
            gmsh.option.setNumber("Mesh.NbPartitions",self.MeshPartitions)
            gmsh.option.setNumber("Mesh.PartitionSplitMeshFiles",1)
            gmsh.option.setNumber("Mesh.PartitionOldStyleMsh2",1)
            gmsh.option.setNumber("Mesh.MetisAlgorithm",1) # 1 scheitert bei make_grid, 2 scheitert bei distribute_dofs
                              
        gmsh.model.mesh.setSizeCallback(self.MeshSizeCallback)

    
    """
    do the meshing
    """
    def mesh(self,
             GUI_popup=False,   # show in graphical user-interface??
             mesh=True):        # mesh, if True, otherwise: only show geometry
        
        self.GeoExport()

        gmsh.initialize(interruptible=False)
    
        gmsh.clear()

        self.mesh_set_options()
        
        if self.verbosity>10:
            print('IMPORT')

        gmsh.merge(self.GeoFileName)

        if self.verbosity>10:
            print('MESHING')
        try:
            if mesh:
                #pass
                gmsh.model.mesh.generate(3)
        except:
            pass
    
        if self.MeshPartitions>1:
            gmsh.model.mesh.partition(self.MeshPartitions)
        #ent = gmsh.model.occ.getEntities(2)
        #print(lent,ent)
        
        gmsh.write(self.MeshFilename)
        log = gmsh.logger.get()
        tmpfarray = self.MeshFilename.split('.')[:-1]
        tmpfarray.append('log')
        logfilename = '.'.join(tmpfarray)
        lfile=open(logfilename,'a')
        lfile.writelines(map(lambda x: x + '\n', log))
        lfile.close()
        print("Logger has recorded " + str(len(log)) + " lines")
        gmsh.logger.stop()
        
        if GUI_popup:
            gmsh.fltk.run()
        
        gmsh.finalize()


    
if __name__=='__main__':
    
    pass

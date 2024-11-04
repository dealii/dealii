#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""
conductors-example

A (set of) minimum example(s) to use gmsh for meshing a model and importing in deal.II based solver
 --- example version 1.0, October 2024 ---

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

A simple example to generate a gmsh Mesh-File out of a set of step-files representing model components

NOTE:
    This is actually a workaround and shall be further developed in the future such that step-files containing 
    complete assemblies and component-IDs should be handled at once 
    instead of having to create separate step-files for each model component
    (as it is in this example and the used classes in it's actual version in October 2024) 

The model can consist of
- 1 or 2 conductors with 2 terminals each
- optionally a solid wall in between two conductors (and none, if only one conductor is modeled)

Each (of the 1 or 2) conductor(s) can be defined in C-shaped geometry or as a simple straight bar.

Shape of the conductor(s) is controlled by setting boolean variable "g_C_shape" - see below.
The number of conductor(s) is controlled by setting boolean variable "g_2_conductors" - see below.
Integration of an optional solid wall between 2 conductors is controlled by setting boolean variable "g_solidW" - see below.

Optionally, the solver can invoked - see function "exec_solver()" and variables "g_solvercmd", "g_prmfile" 
"""
import cadquery as cq                       # needed for model- and step-file generation
import os, sys                              # needed for basic file-operations
from shutil import rmtree                   # needed for basic file-operations
from occ_converter import occ_converter     # converting step-files into serialized (pickle-) objects (incl. file-import/export)

"""
global parameters section:

TO DO: change this to your Linux-user's sub-folder in /home path
"""
g_userhome = 'stephan'

g_C_shape       = True         # if True: make C-shaped conductor(s) with two terminals each 
                                # if False: make just straight bar conductor(s) with two terminals each
g_2_conductors  = True         # if True: make 2 conductor of same shape, otherwise: make only a single conductor
g_solidW        = True         # for 2 conductors only: if True: make solid wall between the two conductors, otherwise: no additional solid
g_meshtype      = 'hexahedra'   # 'hexahedra' or 'tetrahedra'
g_solvercmd     = './Maxwell'   # /bin/sh command for calling the solver
g_prmfile       = 'param_H.prm' # path to parameter file to be used by solver
g_MeshSize      = 0.2           # Mesh element size


outp_type='step'                # 'step' is actually the only available option as only step-files can be handled in this version
outp_fileext = '.'+outp_type

hostsys = sys.platform.lower()

"""
on Windows systems, this example is only running using WSL
"""
Def_LINUX = 'linux' in hostsys or 'unix' in hostsys

g_home = ''
    
# the folder where all model 
g_output_folder = os.path.join(g_home,'tmp/')
g_model_folder = os.path.join(g_output_folder,'model/')

g_scale = 1.0e3 # scale to milli-meter

g_gmsh_options = ['-3', '-a', '-algo del3d', '-clmin 100e-3', '-clmax 300e-3']
gmeshsizes = {}

from ModelComposer import Component, Model
    

"""
main modeling section (.....basically, the example code):
"""    
if __name__=='__main__':

    # 1.) create a new CadQuery assembly
    assy = cq.Assembly()
    models = []     # collect all 
    
    # clean up and re-create output-folders:
    if os.path.exists(g_model_folder):
        rmtree(g_model_folder)

    if not os.path.exists(g_model_folder):
        os.makedirs(g_model_folder)
    
    # model geometry
    terminal_dx = 1.0e-1                                    # terminals thickness
    C_shape_gap_dz = 0.5                                    # z-gap in C-shape conductos 
    shape_dz = 0.2 if g_C_shape else 0.07                   # conductive height (z-coordinate)
    shape_dl = shape_dz                                     # temporary value to cut out a box of correct size for generating the C-shaped conductor(s)
    shape_dy = shape_dz                                     # conductive depth (y-coordinate)
    shapes_distance_dy = 0.8 if g_2_conductors else 0.0     # if 2 conductors: y-distance between conductor y-centers
    shape_dx = 1                                            # x-extent of conductor (not necessarily the conductive length)
    boundary_dx = 1 if g_C_shape else 1.001*2*terminal_dx   # distance from conductor to boundary
    shift_x = 0.4999*boundary_dx-terminal_dx                # shift (C-shaped) conductor so that both electrodes almost touch the boundary 
    theight = 2*shape_dz+C_shape_gap_dz                     # for C-shaped conductor: total height (z)
    
    solid_dy = 0.1*(shapes_distance_dy - shape_dy) if g_2_conductors else 0.01      # y-size of a solid wall between 2 conductor shapes
    solid_dx = 0.8*(shape_dx+boundary_dx)                   # x-size of a solid wall between 2 conductor shapes
    
    # values needed for estimation of ohmic resistance of a single conductor
    l_R = 0         # sum up the conductive length 
    A_R = 0         # sum up conductive cross section
    kappa = 56.0e6  # electric conductivity
    voltage = 2*1.0 # total voltage between two terminals (for calculating the current - potentially needed for Neumann boundary conditions)
    
    if g_C_shape:
        # generate C-shaped conductor (without terminals)
        ubx2 = cq.Workplane("XY").box((shape_dx-shape_dl)*g_scale,shape_dy*g_scale,C_shape_gap_dz*g_scale).translate((0.5*shape_dl*g_scale,0,0))
        ubx = cq.Workplane("XY").box(shape_dx*g_scale,shape_dy*g_scale,theight*g_scale)
        ubx_1 = ubx.cut(ubx2).translate((shift_x*g_scale,-0.5*shapes_distance_dy*g_scale,0))
        ubx_2 = ubx.cut(ubx2).translate((shift_x*g_scale,0.5*shapes_distance_dy*g_scale,0))
        # estimate resistive length and cross section
        l_R = (2*(shape_dx-0.5*shape_dl) + theight-theight)*g_scale
        A_R = shape_dy*g_scale*shape_dz*g_scale
    else:
        # generate bar-shaped conductor (without terminals)
        ubx = cq.Workplane("XY").box(shape_dx*g_scale,shape_dy*g_scale,shape_dy*g_scale)
        ubx_1 = ubx.translate((0,-0.5*shapes_distance_dy*g_scale,0))
        ubx_2 = ubx.translate((0,0.5*shapes_distance_dy*g_scale,0))
        # estimate resistive length and cross section
        l_R = shape_dx*g_scale
        A_R = shape_dy*g_scale*shape_dy*g_scale
        
    # create step-file for conductor 1
    tmp_name = 'conductor1'
    ubx_1.solids().tag(tmp_name)
    tmpassy = cq.Assembly(name=tmp_name).add(ubx_1,color=cq.Color('green'),metadata ={'ID':801,'name':tmp_name})
    tmpassy.save(os.path.join(g_model_folder,tmp_name+outp_fileext))
    models.append(tmp_name)

    # create step-file for conductor 2
    tmp_name = 'conductor2'
    ubx_2.solids().tag(tmp_name)
    tmpassy = cq.Assembly(name=tmp_name).add(ubx_2,color=cq.Color('green'),metadata ={'ID':901,'name':tmp_name})
    tmpassy.save(os.path.join(g_model_folder,tmp_name+outp_fileext))
    models.append(tmp_name)

    if g_C_shape:
        # generate terminal 1 shape for C-shaped conductor
        el = cq.Workplane("XY").box(terminal_dx*g_scale,shape_dy*g_scale,shape_dz*g_scale)
        el1 = el.translate(((shift_x+0.5*shape_dx+0.5*terminal_dx)*g_scale,0.0,0.5*(C_shape_gap_dz+shape_dz)*g_scale))
        el_11 = el1.translate((0,-0.5*shapes_distance_dy*g_scale,0))
    else:
        # generate terminal 1 shape for bar-shaped conductor
        el = cq.Workplane("XY").box(terminal_dx*g_scale,shape_dy*g_scale,shape_dy*g_scale)
        el1 = el.translate(((0.5*shape_dx+0.5*terminal_dx)*g_scale,0.0,0.0))
        el_11 = el1.translate((0,-0.5*shapes_distance_dy*g_scale,0))
        
    # create step-file for terminal 1 from conductor 1
    tmp_name = 'terminal11'
    el_11.solids().tag(tmp_name)
    tmpassy = cq.Assembly(name=tmp_name).add(el_11,color=cq.Color('cyan'),metadata ={'ID':802,'name':tmp_name})
    tmpassy.save(os.path.join(g_model_folder,tmp_name+outp_fileext))
    models.append(tmp_name)

    # create step-file for terminal 1 from conductor 2
    el_21 = el1.translate((0,0.5*shapes_distance_dy*g_scale,0))
    tmp_name = 'terminal21'
    el_21.solids().tag(tmp_name)
    tmpassy = cq.Assembly(name=tmp_name).add(el_21,color=cq.Color('cyan'),metadata ={'ID':902,'name':tmp_name})
    tmpassy.save(os.path.join(g_model_folder,tmp_name+outp_fileext))
    models.append(tmp_name)

    if g_C_shape:
        # generate terminal 2 shape for C-shaped conductor
        el2 = el.translate(((shift_x+0.5*shape_dx+0.5*terminal_dx)*g_scale,0.0,-0.5*(C_shape_gap_dz+shape_dz)*g_scale))
        el_12 = el2.translate((0,-0.5*shapes_distance_dy*g_scale,0))
    else:
        # generate terminal 2 shape for bar-shaped conductor
        el2 = el.translate((-(0.5*shape_dx+0.5*terminal_dx)*g_scale,0.0,0.0))
        el_12 = el2.translate((0,-0.5*shapes_distance_dy*g_scale,0))

    # create step-file for terminal 2 from conductor 1
    tmp_name = 'terminal12'
    el_12.solids().tag(tmp_name)
    tmpassy = cq.Assembly(name=tmp_name,).add(el_12,color=cq.Color('cyan'),metadata ={'ID':803,'name':tmp_name})
    tmpassy.save(os.path.join(g_model_folder,tmp_name+outp_fileext))
    models.append(tmp_name)

    # create step-file for terminal 2 from conductor 2
    el_22 = el2.translate((0,0.5*shapes_distance_dy*g_scale,0))
    tmp_name = 'terminal22'
    el_22.solids().tag(tmp_name)
    tmpassy = cq.Assembly(name=tmp_name,).add(el_22,color=cq.Color('cyan'),metadata ={'ID':903,'name':tmp_name})
    tmpassy.save(os.path.join(g_model_folder,tmp_name+outp_fileext))
    models.append(tmp_name)

    # create step-file for simulation domain
    ebx = cq.Workplane("XY").box((shape_dx+boundary_dx)*g_scale,(shape_dx+boundary_dx)*g_scale,(shape_dx+boundary_dx)*g_scale)#.translate((-0.5,0,0))
    tmp_name = 'domain'
    ebx.solids().tag(tmp_name)
    tmpassy = cq.Assembly(name=tmp_name,).add(ebx,color=cq.Color('cyan'),metadata ={'ID':804,'name':tmp_name})
    tmpassy.save(os.path.join(g_model_folder,tmp_name+outp_fileext))
    models.append(tmp_name)

    # create step-file for solid wall
    epsbx = cq.Workplane("XY").box(solid_dx*g_scale,solid_dy*g_scale,solid_dx*g_scale)#.translate((-0.5,0,0))
    tmp_name = 'solidinsulation'
    epsbx.solids().tag(tmp_name)
    tmpassy = cq.Assembly(name=tmp_name,).add(epsbx,color=cq.Color('cyan'),metadata ={'ID':806,'name':tmp_name})
    tmpassy.save(os.path.join(g_model_folder,tmp_name+outp_fileext))
    models.append(tmp_name)
    
    
    """
    now compose all these step-files of the model components into a single mesh
    
    For this we use the "Model" and "Component" classes from "ModelComposer.py"
    
    The "ModelComposer.py" is designed for
     1) converting step-files into serialized (pickle-) objects and store it in temporary files
     2) import the serialized (pickle-) objects into model components of class type "Component"
     3) manage a set of components in a "Model" class object
     
     The "Model" class object is finally used to merge all the components correctly into a single model.
     Components can be cut out of others or added.
     
     Finally, the merged object can be meshed. gmsh is utilized for that task. 
    """
    t_models = [os.path.join(g_model_folder,x+outp_fileext) for x in models]

    # final mesh-filename
    gmeshfilename = os.path.join(g_model_folder,'testmodel.msh')

    MeshSize = g_MeshSize 
    l_MS_scale = g_scale

    # A) create a "Model" object for managing the components as a model and finally for meshing
    TestObject = Model('TestObject',gmeshfilename,MeshPartitions=1,MeshSize=MeshSize*l_MS_scale,
                        _MeshScalingFactor=1.0/g_scale,MeshSubdivisionAlgorithm=g_meshtype)
    
    TestObject.set_verbosity(2)

    # A.1) Create a "Component" object for the simulation domain
    #      All elements within this domain have to be added or cut-out from this Component!
    TO_domain = Component('domain',t_models[6],PhysicalID=0,MeshSize=MeshSize*l_MS_scale)

    # A.2) Create two "Component"s for the conductor(s) (Maybe not all of them will be added/cut-out from the domain - depending on "g_2_conductors" flag)
    physical_conductor1 = 10
    physical_conductor2 = 20
    physical_solid1 = 30

    TO_conductor1 = Component('conductor1',t_models[0],PhysicalID=physical_conductor1,MeshSize=MeshSize*l_MS_scale)
    TO_conductor2 = Component('conductor2',t_models[1],PhysicalID=physical_conductor2,MeshSize=MeshSize*l_MS_scale)
    
    # A.3) Create a "Component" for the solid wall (Maybe this will not even be added/cut-out from the domain - depending on "g_solidW" flag)
    TO_solid1 = Component('solid1',t_models[7],PhysicalID=physical_solid1,MeshSize=MeshSize*l_MS_scale)

    TO_domain.add(TO_conductor1)        # add 1st. conductor to simulation domain
    if g_2_conductors:
        TO_domain.add(TO_conductor2)    # add 2nd. conductor to simulation domain
    
        if g_solidW:
            TO_domain.add(TO_solid1)    # add solid wall to simulation domain

    # add the simulation domain object to the "Model" object
    TestObject.add(TO_domain)

    # A.4) For applying boundary conditions (either Dirichlet or Neumann), the terminals have to be cut-out from the connected conductors:
    # A.4.1) cut-out the terminals from conductor 1
    TO_conductor1.cut(Component('conductor1_terminal_1',t_models[2],PhysicalID=(2,11),MeshSize=MeshSize*l_MS_scale))
    TO_conductor1.cut(Component('conductor1_terminal_2',t_models[4],PhysicalID=(2,12),MeshSize=MeshSize*l_MS_scale))
    if g_2_conductors:
        # A.4.2) cut-out the terminals from conductor 2
        TO_conductor2.cut(Component('conductor2_terminal_1',t_models[3],PhysicalID=(2,21),MeshSize=MeshSize*l_MS_scale))
        TO_conductor2.cut(Component('conductor2_terminal_2',t_models[5],PhysicalID=(2,22),MeshSize=MeshSize*l_MS_scale))

    do_meshing=True # set this to False, if you only want to see the model geometry in gmsh-GUI
    
    # A.5) finalize model-composition (and do the meshing by using gmsh)
    TestObject.mesh(GUI_popup=True,mesh=do_meshing)

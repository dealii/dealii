====================================================
MESH CONVERSION TOOL
v. 1.2

Use: Convert an ABAQUS .inp file to an AVS UCD file.

Author: Jean-Paul Pelteret
        jppelteret.uct@gmail.com
        modified by: Timo Heister, heister@clemson.edu
===================================================


This tool is intended to support converting ABAQUS .inp files into a
format that deal.II can read. To use this tool, you need to configure
deal.II using the "-DDEAL_II_COMPONENT_MESH_CONVERTER=ON" argument
when calling 'cmake'. The resulting executable will then be placed in
the installation directory.



Cubit features captured by ABAQUS .inp files:
---------------------------------------------
1. Multiple material-id's can be defined in the mesh. This is done by specifying blocksets in the pre-processor.
2. Arbitrary surface boundaries can be defined in the mesh. This is done by specifying sidesets in the pre-processor. In particular, boundaries are not confined to just surfaces (in 3d) - individual element faces can be added to the sideset as well. This is useful when a boundary condition is to be applied on a complex shape boundary that is difficult to define using "surfaces" alone. Similar can be done in 2d.

How to use:
-----------
1. It is necessary to save the mesh as an Abaqus input file .inp. In Cubit, this is done as follows:
1.1. Go to "Analysis setup mode" by clicking on the disc icon in the toolbar on the right.
1.2. Select "Export Mesh" under "Operation" by clicking on the necessary icon in the toolbar on the right.
1.3. Select an output file. I have found that in Cubit version 11.0 and 12.0 that it is necessary to click on the browse button and type it in the dialogue that pops up.
1.4. Select the dimension to output in.
1.5. Tick the overwrite box.
1.6. If using Cubit v12.0 onwards, uncheck the box "Export using Cubit ID's". The conversion process will encounter errors if box if left checked.
1.7. Click apply.

2. Converting the mesh
IMPORTANT NOTE: From v12.0 onwards, the format of the Abaqus file that Cubit outputted changed. This has been accounted for but an requires a flag to be passed to the conversion program.
                I describe outputs from before v12.0 as "Abaqus OLD" and v12.0 an on (up to v12.1 at this point) as "Abaqus NEW".
2.1. Compile the program using 'make' or another suitable tool
2.2. Run the program with the following command line arguments:
       './convert_mesh <spatial_dimension>  /path/to/input_file.inp /path/to/output_file.ucd'

       The first input argument is the spatial dimension of the input file, the second is the path to the read-in Abaqus .inp file, and the third is the name of the file to which you wish to write the output AVS .ucd file to.
       An example of the correct program usage is:
e.g.  './convert_mesh 3 mesh/3d/test_in.inp mesh/3d/test_out.ucd'


Notes:
------
1. A few example / test meshes are provided in the directory "mesh". The .inp files that have been generated are all of the Abaqus OLD format except for one - CC.cub has been converted into both the "new" and "old" formats and both subsequently converted into deal.II readable meshes.

Copyright:
----------
This program is distributed under the GNU GPL v2.0 copyright. Details can be found at:
http://www.gnu.org/licenses/gpl-2.0.html

This copyright is extended to the example mesh files distributed with this program, namely
./mesh/2d/2d_test.*
./mesh/2d/quad.*
./mesh/3d/CC.*
./mesh/3d/CC_cubit_new.*
./mesh/3d/CC_cubit_old.*
./mesh/3d/test_cube_1.*
./mesh/3d/test_cube_pave_1.*
./mesh/3d/test_cube_two_materials.*

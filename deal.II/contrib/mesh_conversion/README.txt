====================================================
MESH CONVERSION TOOL
v. 1.1

Use: Convert an ABAQUS .inp file to an AVS UCD file.

Author: Jean-Paul Pelteret
        jppelteret.uct@gmail.com
===================================================

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
1.6. Click apply.
1.7. IMPORTANT NOTE: From v12.0 onwards, the format of the Abaqus file that Cubit outputted changed. This has been accounted for but an requires a flag to be passed to the conversion program.
                     I describe outputs from before v12.0 as "Abaqus OLD" and v12.0 an on (up to v12.1 at this point) as "Abaqus NEW".

2. Converting the mesh
2.1. Compile the program using 'make' or another suitable tool
2.2. Run the program with the following command line arguments:
       './convert_mesh <spatial_dimension> <ABAQUS input_file_type> /path/to/input_file.inp /path/to/output_file.ucd'

       The first input argument is the spatial dimension of the input file, the second is the type of input file (0 for Abaqus OLD, 1 for Abaqus NEW), the third is the path to the read-in Abaqus .inp file, and the fourth is the name of the file to which you wish to write the output AVS .ucd file to.
       An example of the correct program usage is: 
e.g.  './convert_mesh 3 0 mesh/3d/test_in.inp mesh/3d/test_out.ucd'


Notes:
------
1. This tool was made with the specific intention of the output file being used as a mesh for deal.II.
2. It has been tested with the deal.II 6.2-pre subversion as distributed on 12 Jan 2009.
3. Testing...
   Abaqus OLD: It has been tested for both 2d and 3d meshes, although more thoroughly for the latter. The 2d grids that were tested were drawn in the X-Y plane.
   Abaqus NEW: Has been tested only for 3d meshes (on a complex mesh, which is not provided as an example, but this hopefully covers all the possible problems one may encounter).
4. A few example / test meshes are provided in the directory "mesh". These are all of the Abaqus OLD format.

Copyright:
----------
This program is distributed under the GNU GPL v2.0 copyright. Details can be found at:
http://www.gnu.org/licenses/gpl-2.0.html

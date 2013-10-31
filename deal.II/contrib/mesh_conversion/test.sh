#!/bin/bash
# ======================================================
# MESH CONVERSION TOOL
# v. 1.2
#
# This program is distributed under the terms of the
# GNU GPL v2.0 license.
#
# Author: Jean-Paul Pelteret
#         jppelteret.uct@gmail.com
#         modified by: Timo Heister, heister@clemson.edu
# ======================================================

mesh=mesh/2d/2d_test
./convert_mesh 2 $mesh.inp test/`basename $mesh.ucd` >/dev/null
diff -q $mesh.ucd test/`basename $mesh.ucd`

mesh=mesh/2d/quad
./convert_mesh 2 $mesh.inp test/`basename $mesh.ucd` >/dev/null
diff -q $mesh.ucd test/`basename $mesh.ucd`

mesh=mesh/3d/test_cube_1
./convert_mesh 3 $mesh.inp test/`basename $mesh.ucd` >/dev/null
diff -q $mesh.ucd test/`basename $mesh.ucd`

mesh=mesh/3d/test_cube_two_materials
./convert_mesh 3 $mesh.inp test/`basename $mesh.ucd` >/dev/null
diff -q $mesh.ucd test/`basename $mesh.ucd`

mesh=mesh/3d/CC_cubit_old
./convert_mesh 3 $mesh.inp test/`basename $mesh.ucd` >/dev/null
diff -q $mesh.ucd test/`basename $mesh.ucd`

mesh=mesh/3d/CC_cubit_new
./convert_mesh 3 $mesh.inp test/`basename $mesh.ucd` >/dev/null
diff -q $mesh.ucd test/`basename $mesh.ucd`

mesh=mesh/3d/test_cube_pave_1
./convert_mesh 3 $mesh.inp test/`basename $mesh.ucd` >/dev/null
diff -q $mesh.ucd test/`basename $mesh.ucd`

mesh=mesh/3d/other_simple
./convert_mesh 3 $mesh.inp test/`basename $mesh.ucd` >/dev/null
diff -q $mesh.ucd test/`basename $mesh.ucd`

# mesh=mesh/3d/other
# ./convert_mesh 3 $mesh.inp test/`basename $mesh.ucd` >/dev/null
# diff -q $mesh.ucd test/`basename $mesh.ucd`

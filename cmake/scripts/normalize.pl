## ---------------------------------------------------------------------
##
## Copyright (C) 2001 - 2015 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#
# Remove insignificant volatile data from output files of tests
#
# Data affected:
#  JobID line (containing date)
#  line number of exceptions
#  start and final residual in iterations
#  small doubles
#

# Convert windows to unix line endings. This is necessary to be able to run
# the testsuite on windows (using cygwin's diff/perl)
s/\r$//;

# Remove JobID

s/JobId.*//;

# Remove Input File Name:

s/# Input file name:.*//;

# Several date and time strings

s/%%Creation Date:.*//;
s/\"created\".*//;
s/# Time =.*//;
s/# Date =.*//;
s/^\s+Time =.*//;
s/^\s+Date =.*//;
s/Time tag:.*//g;
s/by the deal.II library on.*//;

# Exceptions

s/line <\d+> of file <.*\//file </;

# See if we have a -0.0... (not followed by any other digit) and replace it
# by the same number without the negative sign
s/-0\.(0+)(?!\d)/0.\1/g;

# remove deal.II debug output
s/^DEAL.*::_.*\n//g;

# Normalize version string by replacing (for example) 'written by
# deal.II 8.1.0-pre' by written by 'written by deal.II x.y.z'
s/written by deal\.II \d+\.\d+\.\d+(-pre|-rc\d*|)/written by deal.II x.y.z/;


# different p4est versions output different text in VTU output. For
# example, we get these kinds of differences:
# ***************
# *** 6,14 ****
#       <PPoints>
#         <PDataArray type="Float32" Name="Position" NumberOfComponents="3" form# at="ascii"/>
#       </PPoints>
# !     <PCellData Scalars="mpirank,treeid">
# !       <PDataArray type="Int32" Name="mpirank" format="ascii"/>
#         <PDataArray type="Int32" Name="treeid" format="ascii"/>
#       </PCellData>
#       <PPointData>
#       </PPointData>
# --- 6,15 ----
#       <PPoints>
#         <PDataArray type="Float32" Name="Position" NumberOfComponents="3" form# at="ascii"/>
#       </PPoints>
# !     <PCellData Scalars="treeid,level,mpirank">
#         <PDataArray type="Int32" Name="treeid" format="ascii"/>
# +       <PDataArray type="UInt8" Name="level" format="ascii"/>
# +       <PDataArray type="Int32" Name="mpirank" format="ascii"/>
#       </PCellData>
#       <PPointData>
#       </PPointData>
#
# To deal with these issues, we simply delete these lines
s/.*<PCellData Scalars.*\n//g;
s/.*<PDataArray type.*(mpirank|level).*\n//g;

#
# Different boost versions output output the opening bracket for json
# output on a new line. Thus always transform
#     "label": {
#
# into
#     "label":
#     {
#
s/^(\s*)(".*":) \{$/\1\2\n\1\{/;

######################################################################
# $Id$
#
# Copyright (C) 2001, the deal.II authors
#
# Remove insignificant volatile data from output files of tests
#
# Data affected:
#  JobID line (containing date)
#  line number of exceptions
#  start and final residual in iterations
#  small doubles
######################################################################

# Remove JobID

s/JobId.*//;

# Several data and time strings

s/%%Creation Date:.*//;
s/\"created\".*//;
s/# Time =.*//;
s/# Date =.*//; 
s/^\s+Time =.*//;
s/^\s+Date =.*//;
s/Time tag:.*//g;

# Exceptions

s/line <\d+> of file <.*\//file </;

# Make small exponentials zero

s/-?\d?\.\d+e-[123456789]\d+/0.00/g;

# All doubles have two decimals
#s/(\.\d\d)\d+/\1/g;

# Zeroes are zero

s/-0\.00/0.00/g;

# Residual values

s/value.*//;

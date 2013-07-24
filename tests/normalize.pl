######################################################################
# $Id$
#
# Copyright (C) 2001, 2003, 2005, 2010, 2011, 2012, 2013, the deal.II authors
#
# Remove insignificant volatile data from output files of tests
#
# Data affected:
#  JobID line (containing date)
#  line number of exceptions
#  start and final residual in iterations
#  small doubles
######################################################################

# Remove absolute path names. Get the path to the testsuite from
# the location of this script
$D = $0;
$D =~ s!tests/normalize.pl!!;
s!$D!DEAL_II_PATH/!g;

# Remove JobID

s/JobId.*//;

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

# Make small exponentials zero

#s/-?\d?\.\d+e-[123456789]\d+/0.00/g;

# See if we have a -0.0... (not followed by any other digit) and replace it
# by the same number without the negative sign
s/-0\.(0+)(?!\d)/0.\1/g;

# Residual values

#s/value.*//;
#s/with residual.*//;


# remove deal.II debug output
s/^DEAL.*::_.*\n//g;

# Normalize version string by replacing (for example) 'written by
# deal.II 8.1.pre' by written by 'written by deal.II x.y.z'
s/written by deal\.II \d+\.\d+\.(pre|\d+)/written by deal.II x.y.z/;

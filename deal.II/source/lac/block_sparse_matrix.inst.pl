######################################################################
# $Id$
######################################################################
#
# Copyright (C) 2011 by the deal.II authors
#
######################################################################

use strict;
use warnings;

use dealiitemplates;

my $scalar_functions = <<'EOT'
template class BlockSparseMatrix<S1>;
EOT
    ;

######################################################################
# End of definitions, now come the loops
######################################################################

foreach my $r1 (@real_scalars)
{
    my $t = $scalar_functions;
    $t =~ s/S1/$r1/g;
    print $t;
}

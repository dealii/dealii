######################################################################
# $Id$
######################################################################
# Postprocess logstream output and create a data file for gnuplot
######################################################################

use strict;

my $step;     # The iteration step in the adaptive loop
my @dofs;     # The number of degrees of freedom in each step
my @error;    # The error of the solution
my @estimate; # The a posteriori error estimate
my @steps;    # The number of multigrid iteration steps
while(<>)
{
    $step = $1 if m/DEAL::Step\s*(\d+)/;
    $dofs[$step] = $1 if m/DEAL::DoFHandler\s*(\d+)/;
    $error[$step] = $1 if m/DEAL::Error\s*(\S+)/;
    $estimate[$step] = $1 if m/DEAL::Estimate\s*(\S+)/;
    $steps[$step] = $1 if m/DEAL:\w+::Convergence step\s*(\S+)/;
}

for (my $i=0;$i<=$step;++$i)
{
    printf "%-3d\t%-7d\t%g\t%g\t%d\n", $i, $dofs[$i], $error[$i], $estimate[$i], $steps[$i];
}

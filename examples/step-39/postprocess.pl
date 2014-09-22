######################################################################
######################################################################
# Postprocess logstream output and create a data file for gnuplot
######################################################################

use strict;

my $step;     # The iteration step in the adaptive loop
my @dofs;     # The number of degrees of freedom in each step
my @error;    # The energy error of the solution
my @l2error;  # The L2-error of the solution
my @estimate; # The a posteriori error estimate
my @steps;    # The number of multigrid iteration steps

my $energy_order = 0.;
my $l2_order = 0.;

while(<>)
{
    $step = $1 if m/DEAL::Step\s*(\d+)/;
    $dofs[$step] = $1 if m/DEAL::DoFHandler\s*(\d+)/;
    $error[$step] = $1 if m/DEAL::energy-error:\s*(\S+)/;
    $l2error[$step] = $1 if m/DEAL::L2-error:\s*(\S+)/;
    $estimate[$step] = $1 if m/DEAL::Estimate\s*(\S+)/;
    $steps[$step] = $1 if m/DEAL:\w+::Convergence step\s*(\S+)/;
}

print '#step dofs error estimate l2error iterations efficiency order l2order', "\n";

for (my $i=0;$i<=$step;++$i)
{
    if ($i>0)
    {
	my $hlog = -1/2.* log($dofs[$i-1]/$dofs[$i]);
	$energy_order = log($error[$i-1]/$error[$i]) / $hlog;
	$l2_order = log($l2error[$i-1]/$l2error[$i]) / $hlog;
    }
    my $eff = $error[$i]/$estimate[$i];
    printf "%-3d\t%-7d\t%e\t%e\t%e\t%d\t%f\t%f\t%f\n", $i, $dofs[$i], $error[$i], $estimate[$i], $l2error[$i],
    $steps[$i], $eff, $energy_order, $l2_order;
}

######################################################################
# $Id$
# Copyright 2005 by the deal.II authors
#
######################################################################
# Extracts xfig and eps files from the output of
#  bits/grid_generator_01.exe
######################################################################
# call:
#  perl thisfile outpufile
######################################################################

# Just open it in order to be able to close it
open OUT, ">0d-void.eps";

while (<>)
{
    if (m/DEAL:(\d)d-GridTest::(.+)/)
    {
	my $name = $2;
	$name .= ".eps" if ($1==1);
	$name .= ".fig" if ($1==2);
	$name .= ".eps" if ($1==3);
	close OUT;
	open OUT, ">$1d-$name\n";
	next;
    }
    print OUT;
}

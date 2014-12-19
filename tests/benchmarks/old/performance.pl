#
# Copyright (C) 2001 by the deal.II authors
#
#
# Make gnuplot table from output of performance
#
# Usage:
#  perl performance.pl performance.log
#
# order of values per line is
#  1. dofs per cell
#  2. Time for generating quadrature points
#  3. Time for generating shape values
#  4. Time for generating shape gradients
#  5. Time for generating shape second derivatives
#  6. Time for generating typical mix for matrix generation
#  7. Time for generating all of the above

while (<>)
{
    if (/(\d+\.\d+):DEAL:(\d+):([^:]+)::End/)
    {
	printf "%3d", $2 if ($3 eq "points");
	printf "\t%f", $1;
	print "\n" if ($3 eq "all---");
    }
}

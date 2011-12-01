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

my $vector_functions = <<'EOT'

template void ConstraintMatrix::condense<V1 >(const V1 &, V1&) const;
template void ConstraintMatrix::condense<V1 >(V1 &vec) const;
template void ConstraintMatrix::set_zero<V1 >(V1&) const;
template void ConstraintMatrix::distribute_local_to_global<V1 > (
    const Vector<double>&, const std::vector<unsigned int> &, V1&, const FullMatrix<double>&) const;
template void ConstraintMatrix::distribute<V1 >(const V1 &, V1&) const;
template void ConstraintMatrix::distribute<V1 >(V1 &) const;
EOT
    ;

my $scalar_vector_functions = <<'EOT'

EOT
    ;

my $scalar_functions = <<'EOT'

template void ConstraintMatrix::condense<S1>(const SparseMatrix<S1>&, SparseMatrix<S1> &) const;
template void ConstraintMatrix::condense<S1>(SparseMatrix<S1>&) const;
template void ConstraintMatrix::condense<S1>(BlockSparseMatrix<S1>&) const;
EOT
    ;

my $scalar_scalar_functions = <<'EOT'
template void ConstraintMatrix::condense<S1,Vector<S2> >(SparseMatrix<S1>&, Vector<S2>&) const;
template void ConstraintMatrix::condense<S1,BlockVector<S2> >(BlockSparseMatrix<S1>&, BlockVector<S2>&) const;
template void ConstraintMatrix::condense<S1,Vector<S2> >(
    const SparseMatrix<S1>&, const Vector<S2>&, SparseMatrix<S1> &, Vector<S2>&) const;
template void ConstraintMatrix::condense<S1,BlockVector<S2> >(
    const SparseMatrix<S1>&, const BlockVector<S2>&, SparseMatrix<S1> &, BlockVector<S2>&) const;
EOT
    ;

######################################################################
# End of definitions, now come the loops
######################################################################

foreach my $v (@sequential_vectors)
{
    my $t = $vector_functions;
    $t =~ s/V1/$v/g;
    print $t;
    
}

foreach my $r1 (@real_scalars)
{
    my $t = $scalar_functions;
    $t =~ s/S1/$r1/g;
    print $t;
    
    foreach my $r2 (@real_scalars)
    {
	my $t = $scalar_scalar_functions;
	$t =~ s/S1/$r1/g;
	$t =~ s/S2/$r2/g;
	print $t
    }
}

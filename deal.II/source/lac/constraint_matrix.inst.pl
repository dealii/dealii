######################################################################
# $Id$
######################################################################
#
# Copyright (C) 2011, 2012 by the deal.II authors
#
# This file is subject to QPL and may not be  distributed
# without copyright and license information. Please refer
# to the file deal.II/doc/license.html for the  text  and
# further information on this license.
#
######################################################################

use strict;
use warnings;
use dealiitemplates;

my $vector_functions = <<'EOT'

template void ConstraintMatrix::condense<V1 >(const V1 &, V1&) const;
template void ConstraintMatrix::condense<V1 >(V1 &vec) const;
template void ConstraintMatrix::distribute_local_to_global<V1 > (
    const Vector<double>&, const std::vector<types::global_dof_index> &, V1&, const FullMatrix<double>&) const;
template void ConstraintMatrix::distribute<V1 >(const V1 &, V1&) const;
template void ConstraintMatrix::distribute<V1 >(V1 &) const;
EOT
    ;

my $vector_functions_also_parallel = <<'EOT'

template void ConstraintMatrix::set_zero<V1 >(V1&) const;

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

multisubst($vector_functions, ['V1'], \@sequential_vectors);
multisubst($vector_functions, ['V1'], \@deal_parallel_vectors);
multisubst($vector_functions_also_parallel, ['V1'], \@sequential_vectors);
multisubst($vector_functions_also_parallel, ['V1'], \@parallel_vectors);
multisubst($scalar_functions, ['S1'], \@real_scalars);
multisubst($scalar_scalar_functions, ['S1','S2'], \@real_scalars, \@real_scalars);

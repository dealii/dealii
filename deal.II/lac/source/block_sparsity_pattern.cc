//----------------------------  block_block_sparsity_pattern.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  block_block_sparsity_pattern.cc  ---------------------------

#include <lac/block_sparsity_pattern.templates.h>

// explicit instantiations
template class BlockSparsityPattern<2,2>;
template class BlockSparsityPattern<2,3>;

template class BlockSparsityPattern<3,2>;
template class BlockSparsityPattern<3,3>;

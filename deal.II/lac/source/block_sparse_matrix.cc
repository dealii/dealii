//----------------------------  block_block_sparse_matrix.cc  ---------------------------
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
//----------------------------  block_block_sparse_matrix.cc  ---------------------------

#include <lac/block_sparse_matrix.templates.h>

// explicit instantiations
template class BlockSparseMatrix<double,2,2>;
template class BlockSparseMatrix<double,2,3>;

template class BlockSparseMatrix<double,3,2>;
template class BlockSparseMatrix<double,3,3>;



template class BlockSparseMatrix<float,2,2>;
template class BlockSparseMatrix<float,2,3>;

template class BlockSparseMatrix<float,3,2>;
template class BlockSparseMatrix<float,3,3>;


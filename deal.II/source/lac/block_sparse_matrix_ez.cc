//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <lac/block_sparse_matrix_ez.h>
#include <lac/block_sparse_matrix_ez.templates.h>

DEAL_II_NAMESPACE_OPEN

// explicit instantiations
template class BlockSparseMatrixEZ<double>;
template class BlockSparseMatrixEZ<float>;

DEAL_II_NAMESPACE_CLOSE

//----------------------------  mg_base.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mg_base.cc  ---------------------------


#include <multigrid/mg_base.h>
#include <lac/vector.h>
#include <lac/block_vector.h>


template <class VECTOR>
MGTransfer<VECTOR>::~MGTransfer()
{}


template <class VECTOR>
MGMatrixBase<VECTOR>::~MGMatrixBase()
{}


template <class VECTOR>
MGCoarseGrid<VECTOR>::~MGCoarseGrid()
{}


// Explicit instantiations

template class MGTransfer<Vector<double> >;
template class MGTransfer<Vector<float> >;
template class MGTransfer<BlockVector<double> >;
template class MGTransfer<BlockVector<float> >;

template class MGMatrixBase<Vector<double> >;
template class MGMatrixBase<Vector<float> >;
template class MGMatrixBase<BlockVector<double> >;
template class MGMatrixBase<BlockVector<float> >;

template class MGCoarseGrid<Vector<double> >;
template class MGCoarseGrid<Vector<float> >;
template class MGCoarseGrid<BlockVector<double> >;
template class MGCoarseGrid<BlockVector<float> >;

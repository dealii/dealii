//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <multigrid/mg_base.h>
#include <lac/vector.h>
#include <lac/block_vector.h>


template <class VECTOR>
MGTransferBase<VECTOR>::~MGTransferBase()
{}


template <class VECTOR>
MGMatrixBase<VECTOR>::~MGMatrixBase()
{}


template <class VECTOR>
MGSmootherBase<VECTOR>::~MGSmootherBase()
{}


template <class VECTOR>
MGCoarseGridBase<VECTOR>::~MGCoarseGridBase()
{}


// Explicit instantiations

template class MGTransferBase<Vector<double> >;
template class MGTransferBase<Vector<float> >;
template class MGTransferBase<BlockVector<double> >;
template class MGTransferBase<BlockVector<float> >;

template class MGMatrixBase<Vector<double> >;
template class MGMatrixBase<Vector<float> >;
template class MGMatrixBase<BlockVector<double> >;
template class MGMatrixBase<BlockVector<float> >;

template class MGSmootherBase<Vector<float> >;
template class MGSmootherBase<Vector<double> >;
template class MGSmootherBase<BlockVector<float> >;
template class MGSmootherBase<BlockVector<double> >;

template class MGCoarseGridBase<Vector<double> >;
template class MGCoarseGridBase<Vector<float> >;
template class MGCoarseGridBase<BlockVector<double> >;
template class MGCoarseGridBase<BlockVector<float> >;

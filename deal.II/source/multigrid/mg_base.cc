//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2005, 2006, 2010, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/parallel_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/multigrid/mg_base.h>


DEAL_II_NAMESPACE_OPEN


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

#include "mg_base.inst"

DEAL_II_NAMESPACE_CLOSE

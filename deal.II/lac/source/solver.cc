//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <lac/solver.h>
#include <lac/vector_memory.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/petsc_vector.h>
#include <lac/petsc_block_vector.h>


// a few instantiations of static members. Hope that we catch all that
// are required

template <class VECTOR>
PrimitiveVectorMemory<VECTOR>
Solver<VECTOR>::static_vector_memory;

template class Solver<Vector<double> >;
template class Solver<BlockVector<double> >;

template class Solver<Vector<float> >;
template class Solver<BlockVector<float> >;

#ifdef DEAL_II_USE_PETSC
template class Solver<PETScWrappers::Vector>;
template class Solver<PETScWrappers::BlockVector>;
#endif

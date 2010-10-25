//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/vector_memory.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/petsc_vector.h>
#include <lac/petsc_block_vector.h>
#include <lac/petsc_parallel_vector.h>
#include <lac/petsc_parallel_block_vector.h>
#include <lac/trilinos_vector.h>
#include <lac/trilinos_block_vector.h>

DEAL_II_NAMESPACE_OPEN


template <typename VECTOR>
typename GrowingVectorMemory<VECTOR>::Pool GrowingVectorMemory<VECTOR>::pool;

template <typename VECTOR>
Threads::ThreadMutex GrowingVectorMemory<VECTOR>::mutex;


// -------------------------------------------------------------
// explicit instantiations

#include "vector_memory.inst"

//TODO: Fold this into the list of vectors to be instantiated
#ifdef DEAL_II_USE_PETSC
    template class VectorMemory<PETScWrappers::Vector>;
    template class GrowingVectorMemory<PETScWrappers::Vector>;

    template class VectorMemory<PETScWrappers::BlockVector>;
    template class GrowingVectorMemory<PETScWrappers::BlockVector>;

    template class VectorMemory<PETScWrappers::MPI::Vector>;
    template class GrowingVectorMemory<PETScWrappers::MPI::Vector>;

    template class VectorMemory<PETScWrappers::MPI::BlockVector>;
    template class GrowingVectorMemory<PETScWrappers::MPI::BlockVector>;
#endif

#ifdef DEAL_II_USE_TRILINOS
    template class VectorMemory<TrilinosWrappers::Vector>;
    template class GrowingVectorMemory<TrilinosWrappers::Vector>;

    template class VectorMemory<TrilinosWrappers::BlockVector>;
    template class GrowingVectorMemory<TrilinosWrappers::BlockVector>;

    template class VectorMemory<TrilinosWrappers::MPI::Vector>;
    template class GrowingVectorMemory<TrilinosWrappers::MPI::Vector>;

    template class VectorMemory<TrilinosWrappers::MPI::BlockVector>;
    template class GrowingVectorMemory<TrilinosWrappers::MPI::BlockVector>;
#endif

DEAL_II_NAMESPACE_CLOSE

//----------------------------  petsc_parallel_block_vector.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  petsc_parallel_block_vector.cc  ---------------------------


#include <lac/petsc_parallel_block_vector.h>
#include <lac/petsc_block_vector.h>

  
#ifdef DEAL_II_USE_PETSC


namespace PETScWrappers
{
  namespace MPI
  {
    BlockVector &
    BlockVector::operator = (const PETScWrappers::BlockVector &v)
    {
      Assert (v.get_block_indices() == this->get_block_indices(),
              ExcNonMatchingBlockVectors());

      for (unsigned int i=0; i<this->n_blocks(); ++i)
        this->block(i) = v.block(i);
      
      return *this;
    }
  }

}

#else
// On gcc2.95 on Alpha OSF1, the native assembler does not like empty
// files, so provide some dummy code
namespace { void dummy () {} }
#endif // DEAL_II_USE_PETSC

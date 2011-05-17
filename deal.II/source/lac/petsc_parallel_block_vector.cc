//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005, 2006, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <deal.II/lac/petsc_parallel_block_vector.h>
  
#ifdef DEAL_II_USE_PETSC

#  include <deal.II/lac/petsc_block_vector.h>

DEAL_II_NAMESPACE_OPEN

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

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_PETSC

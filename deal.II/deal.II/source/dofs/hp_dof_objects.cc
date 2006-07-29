//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/memory_consumption.h>
#include <dofs/hp_dof_objects.h>
#include <dofs/hp_dof_levels.h>
#include <dofs/hp_dof_handler.h>


//TODO: Several of the functions in this file are identical. some template trickery should make it possible to merge them

namespace internal
{
  namespace hp
  {
    template <int structdim>
    unsigned int
    DoFObjects<structdim>::memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (dofs) +
              MemoryConsumption::memory_consumption (dof_offsets));
    }


				     // explicit instantiations
    template
    unsigned int
    DoFObjects<1>::
    memory_consumption () const;
    
#if deal_II_dimension >= 2

    template
    unsigned int
    DoFObjects<2>::
    memory_consumption () const;
    
#endif
    
    
#if deal_II_dimension >= 3

    template
    unsigned int
    DoFObjects<3>::
    memory_consumption () const;
    
#endif
  }
}

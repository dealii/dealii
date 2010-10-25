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
#include <hp/dof_objects.h>

DEAL_II_NAMESPACE_OPEN

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
    
    template
    unsigned int
    DoFObjects<2>::
    memory_consumption () const;
    
    template
    unsigned int
    DoFObjects<3>::
    memory_consumption () const;
  }
}

DEAL_II_NAMESPACE_CLOSE

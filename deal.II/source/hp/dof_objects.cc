//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <deal.II/base/memory_consumption.h>
#include <deal.II/hp/dof_objects.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace hp
  {
    template <int structdim>
    std::size_t
    DoFObjects<structdim>::memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (dofs) +
              MemoryConsumption::memory_consumption (dof_offsets));
    }


				     // explicit instantiations
    template
    std::size_t
    DoFObjects<1>::memory_consumption () const;
    
    template
    std::size_t
    DoFObjects<2>::memory_consumption () const;
    
    template
    std::size_t
    DoFObjects<3>::memory_consumption () const;
  }
}

DEAL_II_NAMESPACE_CLOSE

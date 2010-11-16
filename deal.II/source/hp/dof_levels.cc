//----------------------------  hp_dof_levels.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  hp_dof_levels.cc  ------------------------


#include <base/memory_consumption.h>
#include <hp/dof_levels.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace hp
  {

    unsigned int
    DoFLevel<1>::memory_consumption () const
    {
      return (DoFLevel<0>::memory_consumption() +
              MemoryConsumption::memory_consumption (lines));
    }
    


    unsigned int
    DoFLevel<2>::memory_consumption () const
    {
      return (DoFLevel<0>::memory_consumption () +
              MemoryConsumption::memory_consumption (quads));
    }



    unsigned int
    DoFLevel<3>::memory_consumption () const
    {
      return (DoFLevel<0>::memory_consumption () +
              MemoryConsumption::memory_consumption (hexes));
    }



    unsigned int
    DoFLevel<0>::memory_consumption () const
    {
      return MemoryConsumption::memory_consumption (active_fe_indices);
    }
    
  }
}

DEAL_II_NAMESPACE_CLOSE

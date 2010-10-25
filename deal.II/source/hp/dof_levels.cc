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

#if deal_II_dimension == 1

    unsigned int
    DoFLevel<1>::memory_consumption () const
    {
      return (DoFLevel<0>::memory_consumption() +
              MemoryConsumption::memory_consumption (lines));
    }
    
#endif


#if deal_II_dimension == 2

    unsigned int
    DoFLevel<2>::memory_consumption () const
    {
      return (DoFLevel<0>::memory_consumption () +
              MemoryConsumption::memory_consumption (quads));
    }

#endif


#if deal_II_dimension == 3

    unsigned int
    DoFLevel<3>::memory_consumption () const
    {
      return (DoFLevel<0>::memory_consumption () +
              MemoryConsumption::memory_consumption (hexes));
    }

#endif

  }
}

DEAL_II_NAMESPACE_CLOSE

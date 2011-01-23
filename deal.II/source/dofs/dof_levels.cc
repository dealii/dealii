//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2005, 2006, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/memory_consumption.h>
#include <dofs/dof_levels.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace DoFHandler
  {
    
    std::size_t
    DoFLevel<0>::memory_consumption () const
    {
      return MemoryConsumption::memory_consumption (cell_dof_indices_cache);
    }


    
    std::size_t
    DoFLevel<1>::memory_consumption () const
    {
      return (DoFLevel<0>::memory_consumption () +
              MemoryConsumption::memory_consumption (lines));
    }

    
    
    std::size_t
    DoFLevel<2>::memory_consumption () const
    {
      return (DoFLevel<0>::memory_consumption () +
              MemoryConsumption::memory_consumption (quads));
    }

    

    std::size_t
    DoFLevel<3>::memory_consumption () const
    {
      return (DoFLevel<0>::memory_consumption () +
              MemoryConsumption::memory_consumption (hexes));
    }
    
  }
  
}

DEAL_II_NAMESPACE_CLOSE

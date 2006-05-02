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
#include <dofs/hp_dof_levels.h>


namespace internal
{
  namespace hp
  {
    unsigned int
    DoFLevel<1>::memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (line_dofs) +
              MemoryConsumption::memory_consumption (line_dof_offsets));
    }



    unsigned int
    DoFLevel<2>::memory_consumption () const
    {
      return (DoFLevel<1>::memory_consumption () +
              MemoryConsumption::memory_consumption (quad_dofs) +
              MemoryConsumption::memory_consumption (quad_dof_offsets));
    }



    unsigned int
    DoFLevel<3>::memory_consumption () const
    {
      return (DoFLevel<2>::memory_consumption () +
              MemoryConsumption::memory_consumption (hex_dofs) +
              MemoryConsumption::memory_consumption (hex_dof_offsets));
    }
  }
}

//----------------------------  hp_dof_levels.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  hp_dof_levels.cc  ---------------------------


#include <base/memory_consumption.h>
#include <dofs/hp_dof_levels.h>


namespace internal
{
  
  unsigned int
  hpDoFLevel<1>::memory_consumption () const
  {
    return MemoryConsumption::memory_consumption (line_dofs);
  }



  unsigned int
  hpDoFLevel<2>::memory_consumption () const
  {
    return (hpDoFLevel<1>::memory_consumption () +
            MemoryConsumption::memory_consumption (quad_dofs));
  }



  unsigned int
  hpDoFLevel<3>::memory_consumption () const
  {
    return (hpDoFLevel<2>::memory_consumption () +
            MemoryConsumption::memory_consumption (hex_dofs));
  }
}

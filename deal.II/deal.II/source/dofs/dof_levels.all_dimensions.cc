//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/exceptions.h>
#include <base/memory_consumption.h>
#include <dofs/dof_levels.h>


namespace internal
{
  namespace DoFHandler
  {
    
    unsigned int
    DoFLevel<0>::memory_consumption () const
    {
      return MemoryConsumption::memory_consumption (cell_dof_indices_cache);
    }

  }
  
}

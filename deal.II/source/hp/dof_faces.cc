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
#include <hp/dof_faces.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace hp
  {

    unsigned int
    DoFFaces<1>::memory_consumption () const
    {
      return 0;
    }



    unsigned int
    DoFFaces<2>::memory_consumption () const
    {
      return MemoryConsumption::memory_consumption (lines);
    }



    unsigned int
    DoFFaces<3>::memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (quads) +
	      MemoryConsumption::memory_consumption (quads) );
    }

  }
}

DEAL_II_NAMESPACE_CLOSE

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
#include <deal.II/hp/dof_faces.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace hp
  {

    std::size_t
    DoFFaces<1>::memory_consumption () const
    {
      return 0;
    }



    std::size_t
    DoFFaces<2>::memory_consumption () const
    {
      return MemoryConsumption::memory_consumption (lines);
    }



    std::size_t
    DoFFaces<3>::memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (quads) +
	      MemoryConsumption::memory_consumption (quads) );
    }

  }
}

DEAL_II_NAMESPACE_CLOSE

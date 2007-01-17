//---------------------------------------------------------------------------
//    $Id: dof_objects.cc 14220 2006-12-04 23:49:48Z bangerth $
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


#include <base/exceptions.h>
#include <base/memory_consumption.h>
#include <dofs/dof_objects.h>
#include <dofs/dof_handler.h>
#include <fe/fe.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace DoFHandler
  {
    template <int dim>
    unsigned int
    DoFObjects<dim>::memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (dofs));
    }


    // explicit instantiations
    template unsigned int DoFObjects<1>::memory_consumption () const;
    template unsigned int DoFObjects<2>::memory_consumption () const;
    template unsigned int DoFObjects<3>::memory_consumption () const;
  }
}

DEAL_II_NAMESPACE_CLOSE

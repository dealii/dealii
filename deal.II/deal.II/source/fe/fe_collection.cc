//----------------------------  fe_collection.cc  ---------------------------
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
//----------------------------  fe_collection.cc  ---------------------------

#include <base/memory_consumption.h>
#include <fe/fe_collection.h>


template <int dim>
unsigned int
FECollection<dim>::memory_consumption () const
{
  return (sizeof(*this) +
          MemoryConsumption::memory_consumption (finite_elements));
}




// explicit instantiations
template class FECollection<deal_II_dimension>;

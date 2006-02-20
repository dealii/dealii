//----------------------------  fe_collection.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_collection.cc  ---------------------------

#include <base/memory_consumption.h>
#include <fe/fe_collection.h>


namespace hp
{
  template <int dim>
  void FECollection<dim>::push_back (const FiniteElement<dim> &new_fe)
  {
                                     // check that the new element has the right
                                     // number of components. only check with
                                     // the first element, since all the other
                                     // elements have already passed the test
                                     // against the first element
    if (finite_elements.size() != 0)
      Assert (new_fe.n_components() == finite_elements[0]->n_components(),
              ExcMessage ("All elements inside a collection need to have the "
                          "same number of vector components!"));
  
    finite_elements
      .push_back (boost::shared_ptr<const FiniteElement<dim> >(new_fe.clone()));
  }



  template <int dim>
  unsigned int
  FECollection<dim>::memory_consumption () const
  {
    unsigned int mem
      = (sizeof(*this) +
         MemoryConsumption::memory_consumption (finite_elements));
    for (unsigned int i=0; i<finite_elements.size(); ++i)
      mem += finite_elements[i]->memory_consumption();

    return mem;
  }


// explicit instantiations
  template class FECollection<deal_II_dimension>;
  
}

//----------------------------  q_collection.cc  ----------------------------
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
//----------------------------  q_collection.cc  ----------------------------

#include <base/memory_consumption.h>
#include <fe/q_collection.h>



namespace hp
{
  
  template <int dim>
  QCollection<dim>::QCollection ()
                  :
                  single_quadrature (false)
  {
  }


  template <int dim>
  QCollection<dim>::QCollection (const Quadrature<dim> &quadrature)
                  :
                  single_quadrature (true)
  {
    quadratures
      .push_back (boost::shared_ptr<const Quadrature<dim> >(new Quadrature<dim>(quadrature)));
  }

  

  template <int dim>
  inline
  const Quadrature<dim> &
  QCollection<dim>::operator[] (const unsigned int index) const
  {
                                     // if we have only a single quadrature
                                     // that was given during construction,
                                     // return this one. otherwise pick out
                                     // the correct one
    if (single_quadrature)
      return *quadratures[0];
    else
      {
	Assert (index < quadratures.size (),
		ExcIndexRange (index, 0, quadratures.size ()));
	return *quadratures[index];
      }
  }


  template <int dim>
  unsigned int
  QCollection<dim>::memory_consumption () const
  {
    return (sizeof(*this) +
	    MemoryConsumption::memory_consumption (quadratures));
  }


  template <int dim>
  unsigned int
  QCollection<dim>::
  push_back (const Quadrature<dim> &new_quadrature)
  {
                                     // A QCollection, which was initialized
                                     // as single QCollection cannot
                                     // administrate other Quadratures.
    Assert (!single_quadrature, ExcNotInitialized ());

    quadratures
      .push_back (boost::shared_ptr<const Quadrature<dim> >(new Quadrature<dim>(new_quadrature)));
    return quadratures.size ();
  }


// explicit instantiations
  template class QCollection<deal_II_dimension>;
#if deal_II_dimension >= 2
  template class QCollection<deal_II_dimension-1>;
#endif

  
}

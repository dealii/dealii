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

DEAL_II_NAMESPACE_OPEN


namespace hp
{
  
  template <int dim>
  QCollection<dim>::QCollection ()
  {}


  template <int dim>
  QCollection<dim>::QCollection (const Quadrature<dim> &quadrature)
  {
    quadratures
      .push_back (boost::shared_ptr<const Quadrature<dim> >(new Quadrature<dim>(quadrature)));
  }

  

  template <int dim>
  QCollection<dim>::
  QCollection (const QCollection<dim> &q_collection)
                  :
                  Subscriptor (),
                                                   // copy the array
                                                   // of shared
                                                   // pointers. nothing
                                                   // bad should
                                                   // happen -- they
                                                   // simply all point
                                                   // to the same
                                                   // objects, and the
                                                   // last one to die
                                                   // will delete the
                                                   // mappings
                  quadratures (q_collection.quadratures)
  {}



  template <int dim>
  unsigned int
  QCollection<dim>::memory_consumption () const
  {
    return (sizeof(*this) +
	    MemoryConsumption::memory_consumption (quadratures));
  }


  template <int dim>
  void
  QCollection<dim>::push_back (const Quadrature<dim> &new_quadrature)
  {
    quadratures
      .push_back (boost::shared_ptr<const Quadrature<dim> >(new Quadrature<dim>(new_quadrature)));
  }


// explicit instantiations
  template class QCollection<deal_II_dimension>;
#if deal_II_dimension >= 2
  template class QCollection<deal_II_dimension-1>;
#endif

  
}
DEAL_II_NAMESPACE_CLOSE

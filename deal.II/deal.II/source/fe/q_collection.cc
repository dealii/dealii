//----------------------------  q_collection.cc  ----------------------------
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
//----------------------------  q_collection.cc  ----------------------------

#include <base/memory_consumption.h>
#include <fe/q_collection.h>


template <int dim>
QCollection<dim>::QCollection () :
    single_quadrature (false)
{
}


template <int dim>
QCollection<dim>::QCollection (const Quadrature<dim> &quadrature) :
    single_quadrature (true)
{
    quadratures.push_back (&quadrature);
}


template <int dim>
inline
const Quadrature<dim> &
QCollection<dim>::get_quadrature (const unsigned int active_fe_index) const
{
    SmartPointer<const Quadrature<dim> > quad;

    if (single_quadrature)
	quad = quadratures[0];
    else
    {
	Assert (active_fe_index < quadratures.size (),
		ExcIndexRange (active_fe_index, 0, quadratures.size ()));
	quad = quadratures[active_fe_index];
    }

    return *quad;
}


template <int dim>
unsigned int
QCollection<dim>::memory_consumption () const
{
    return (sizeof(*this) +
	    MemoryConsumption::memory_consumption (quadratures));
}


template <int dim>
unsigned int QCollection<dim>::
add_quadrature (const Quadrature<dim> &new_quadrature)
{
    // A QCollection, which was initialized as single QCollection cannot
    // administrate other Quadratures.
    Assert (!single_quadrature,
	    ExcNotInitialized ());
    quadratures.push_back (&new_quadrature);
    return (quadratures.size ());
}


// explicit instantiations
template class QCollection<deal_II_dimension>;
#if deal_II_dimension >= 2
template class QCollection<deal_II_dimension-1>;
#endif

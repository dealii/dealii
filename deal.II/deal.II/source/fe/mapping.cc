//------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//------------------------------------------------------------------------

#include <grid/tria.h>
#include <fe/mapping.h>


#if (deal_II_dimension == 1)

template<>
const unsigned int Mapping<1>::normal_directions[2] =
{
  1, 0
};

#endif

#if (deal_II_dimension == 2)

template<>
const unsigned int Mapping<2>::normal_directions[4] =
{
  2, 0, 3, 1
};

#endif

#if (deal_II_dimension == 3)

template<>
const unsigned int Mapping<3>::normal_directions[6] =
{
  3, 2, 5, 0, 4, 1
};

#endif


template <int dim>
Mapping<dim>::~Mapping ()
{}


/*------------------------------ InternalDataBase ------------------------------*/


template <int dim>
Mapping<dim>::InternalDataBase::InternalDataBase ():
		update_flags(update_default),
		update_once(update_default),
		update_each(update_default),
		first_cell(true)
{}



template <int dim>
Mapping<dim>::InternalDataBase::~InternalDataBase ()
{}



template <int dim>
unsigned int
Mapping<dim>::InternalDataBase::memory_consumption () const
{
  return sizeof(*this);
}


/*------------------------------ InternalData ------------------------------*/



template class Mapping<deal_II_dimension>;




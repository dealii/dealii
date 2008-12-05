//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <grid/tria.h>
#include <fe/mapping.h>

DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim>
Mapping<dim, spacedim>::~Mapping ()
{}


/*------------------------------ InternalDataBase ------------------------------*/


template <int dim, int spacedim>
Mapping<dim, spacedim>::InternalDataBase::InternalDataBase ():
		update_flags(update_default),
		update_once(update_default),
		update_each(update_default),
		first_cell(true)
{}



template <int dim, int spacedim>
Mapping<dim, spacedim>::InternalDataBase::~InternalDataBase ()
{}



template <int dim, int spacedim>
unsigned int
Mapping<dim, spacedim>::InternalDataBase::memory_consumption () const
{
  return sizeof(*this);
}


/*------------------------------ InternalData ------------------------------*/



template class Mapping<deal_II_dimension>;

#if deal_II_dimension != 3
template class Mapping<deal_II_dimension,deal_II_dimension+1>;
#endif

DEAL_II_NAMESPACE_CLOSE

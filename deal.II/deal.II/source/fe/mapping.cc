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

DEAL_II_NAMESPACE_CLOSE

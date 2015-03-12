// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include <deal.II/grid/tria.h>
#include <deal.II/fe/mapping.h>

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
std::size_t
Mapping<dim, spacedim>::InternalDataBase::memory_consumption () const
{
  return sizeof(*this);
}


/*------------------------------ InternalData ------------------------------*/



// explicit instantiations
#include "mapping.inst"


DEAL_II_NAMESPACE_CLOSE

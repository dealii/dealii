// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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


// This function is deprecated and has been replaced by transform above
template<int dim, int spacedim>
void
Mapping<dim,spacedim>::transform_covariant (
  const VectorSlice<const std::vector<Tensor<1,dim> > > input,
  const unsigned int                 offset,
  VectorSlice<std::vector<Tensor<1,spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data) const
{
  Assert (offset == 0, ExcInternalError());

  transform(input, output, mapping_data, mapping_covariant);
}



// This function is deprecated and has been replaced by transform above
template <int dim, int spacedim>
void
Mapping<dim, spacedim>::transform_covariant (
  const VectorSlice<const std::vector<DerivativeForm<1, dim,spacedim> > > input,
  const unsigned int                 offset,
  VectorSlice<std::vector<Tensor<2,spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data) const
{
  Assert (offset == 0, ExcInternalError());

  transform(input, output, mapping_data, mapping_covariant);
}



// This function is deprecated and has been replaced by transform above
template<int dim, int spacedim>
void
Mapping<dim,spacedim>::transform_contravariant (
  const VectorSlice<const std::vector<Tensor<1,dim> > > input,
  const unsigned int                 offset,
  VectorSlice<std::vector<Tensor<1,spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data) const
{
  Assert (offset == 0, ExcInternalError());

  transform(input, output, mapping_data, mapping_contravariant);
}



// This function is deprecated and has been replaced by transform above
template<int dim, int spacedim>
void
Mapping<dim,spacedim>::transform_contravariant (
  const VectorSlice<const std::vector<DerivativeForm<1, dim,spacedim> > > input,
  const unsigned int                 offset,
  VectorSlice<std::vector<Tensor<2,spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data) const
{
  Assert (offset == 0, ExcInternalError());

  transform(input, output, mapping_data, mapping_contravariant);
}


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

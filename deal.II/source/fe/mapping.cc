//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2005, 2006, 2009, 2011 by the deal.II authors
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
    const VectorSlice<const std::vector<Tensor<2,dim> > > input,
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
  const VectorSlice<const std::vector<Tensor<2,dim> > > input,
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

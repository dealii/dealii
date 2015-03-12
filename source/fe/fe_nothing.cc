// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2014 by the deal.II authors
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


#include <deal.II/fe/fe_nothing.h>

DEAL_II_NAMESPACE_OPEN

namespace
{
  const char *
  zero_dof_message = "This element has no shape functions.";
}




template <int dim, int spacedim>
FE_Nothing<dim,spacedim>::FE_Nothing (const unsigned int n_components)
  :
  FiniteElement<dim,spacedim>
  (FiniteElementData<dim>(std::vector<unsigned>(dim+1,0),
                          n_components, 0,
                          FiniteElementData<dim>::unknown),
   std::vector<bool>(),
   std::vector<ComponentMask>() )
{
// in most other elements we have to set up all sorts of stuff
// here. there isn't much that we have to do here; in particular,
// we can simply leave the restriction and prolongation matrices
// empty since their proper size is in fact zero given that the
// element here has no degrees of freedom
}


template <int dim, int spacedim>
FiniteElement<dim,spacedim> *
FE_Nothing<dim,spacedim>::clone() const
{
  return new FE_Nothing<dim,spacedim>(*this);
}



template <int dim, int spacedim>
std::string
FE_Nothing<dim,spacedim>::get_name () const
{
  std::ostringstream namebuf;
  namebuf << "FE_Nothing<" << dim << ">(";
  if (this->n_components() > 1)
    namebuf << this->n_components();
  namebuf << ")";
  return namebuf.str();
}



template <int dim, int spacedim>
UpdateFlags
FE_Nothing<dim,spacedim>::update_once (const UpdateFlags /*flags*/) const
{
  return update_default;
}



template <int dim, int spacedim>
UpdateFlags
FE_Nothing<dim,spacedim>::update_each (const UpdateFlags /*flags*/) const
{
  return update_default;
}



template <int dim, int spacedim>
double
FE_Nothing<dim,spacedim>::shape_value (const unsigned int /*i*/,
                                       const Point<dim> & /*p*/) const
{
  Assert(false,ExcMessage(zero_dof_message));
  return 0;
}



template <int dim, int spacedim>
typename Mapping<dim,spacedim>::InternalDataBase *
FE_Nothing<dim,spacedim>::get_data (const UpdateFlags  /*flags*/,
                                    const Mapping<dim,spacedim> & /*mapping*/,
                                    const Quadrature<dim> & /*quadrature*/) const
{
  // Create a default data object.  Normally we would then
  // need to resize things to hold the appropriate numbers
  // of dofs, but in this case all data fields are empty.
  typename Mapping<dim,spacedim>::InternalDataBase *data
    = new typename FiniteElement<dim,spacedim>::InternalDataBase();
  return data;
}



template <int dim, int spacedim>
void
FE_Nothing<dim,spacedim>::
fill_fe_values (const Mapping<dim,spacedim> & /*mapping*/,
                const typename Triangulation<dim,spacedim>::cell_iterator & /*cell*/,
                const Quadrature<dim> & /*quadrature*/,
                typename Mapping<dim,spacedim>::InternalDataBase & /*mapping_data*/,
                typename Mapping<dim,spacedim>::InternalDataBase & /*fedata*/,
                FEValuesData<dim,spacedim> & /*data*/,
                CellSimilarity::Similarity & /*cell_similarity*/) const
{
  // leave data fields empty
}



template <int dim, int spacedim>
void
FE_Nothing<dim,spacedim>::
fill_fe_face_values (const Mapping<dim,spacedim> & /*mapping*/,
                     const typename Triangulation<dim,spacedim>::cell_iterator & /*cell*/,
                     const unsigned int /*face*/,
                     const Quadrature<dim-1> & /*quadrature*/,
                     typename Mapping<dim,spacedim>::InternalDataBase & /*mapping_data*/,
                     typename Mapping<dim,spacedim>::InternalDataBase & /*fedata*/,
                     FEValuesData<dim,spacedim> & /*data*/) const
{
  // leave data fields empty
}

template <int dim, int spacedim>
void
FE_Nothing<dim,spacedim>::
fill_fe_subface_values (const Mapping<dim,spacedim> & /*mapping*/,
                        const typename Triangulation<dim,spacedim>::cell_iterator & /*cell*/,
                        const unsigned int /*face*/,
                        const unsigned int /*subface*/,
                        const Quadrature<dim-1> & /*quadrature*/,
                        typename Mapping<dim,spacedim>::InternalDataBase & /*mapping_data*/,
                        typename Mapping<dim,spacedim>::InternalDataBase & /*fedata*/,
                        FEValuesData<dim,spacedim> & /*data*/) const
{
  // leave data fields empty
}


template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_Nothing<dim,spacedim> ::
compare_for_face_domination (const FiniteElement<dim,spacedim> &) const
{
  return FiniteElementDomination::no_requirements;
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_Nothing<dim,spacedim> ::
hp_vertex_dof_identities (const FiniteElement<dim,spacedim> &/*fe_other*/) const
{
  // the FE_Nothing has no
  // degrees of freedom, so there
  // are no equivalencies to be
  // recorded
  return std::vector<std::pair<unsigned int, unsigned int> > ();
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_Nothing<dim,spacedim> ::
hp_line_dof_identities (const FiniteElement<dim,spacedim> &/*fe_other*/) const
{
  // the FE_Nothing has no
  // degrees of freedom, so there
  // are no equivalencies to be
  // recorded
  return std::vector<std::pair<unsigned int, unsigned int> > ();
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_Nothing<dim,spacedim> ::
hp_quad_dof_identities (const FiniteElement<dim,spacedim> &/*fe_other*/) const
{
  // the FE_Nothing has no
  // degrees of freedom, so there
  // are no equivalencies to be
  // recorded
  return std::vector<std::pair<unsigned int, unsigned int> > ();
}


template <int dim, int spacedim>
bool
FE_Nothing<dim,spacedim> ::
hp_constraints_are_implemented () const
{
  return true;
}


template <int dim, int spacedim>
void
FE_Nothing<dim,spacedim>::
get_face_interpolation_matrix (const FiniteElement<dim,spacedim> &/*source_fe*/,
                               FullMatrix<double>       &interpolation_matrix) const
{
  // since this element has no face dofs, the
  // interpolation matrix is necessarily empty

  Assert (interpolation_matrix.m() == 0,
          ExcDimensionMismatch (interpolation_matrix.m(),
                                0));
  Assert (interpolation_matrix.n() == 0,
          ExcDimensionMismatch (interpolation_matrix.m(),
                                0));
}


template <int dim, int spacedim>
void
FE_Nothing<dim,spacedim>::
get_subface_interpolation_matrix (const FiniteElement<dim,spacedim> & /*source_fe*/,
                                  const unsigned int /*index*/,
                                  FullMatrix<double>  &interpolation_matrix) const
{
  // since this element has no face dofs, the
  // interpolation matrix is necessarily empty

  Assert (interpolation_matrix.m() == 0,
          ExcDimensionMismatch (interpolation_matrix.m(),
                                0));
  Assert (interpolation_matrix.n() == 0,
          ExcDimensionMismatch (interpolation_matrix.m(),
                                0));
}



// explicit instantiations
#include "fe_nothing.inst"


DEAL_II_NAMESPACE_CLOSE


// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/fe/fe_nothing.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim>
FE_Nothing<dim, spacedim>::FE_Nothing(const ReferenceCell &type,
                                      const unsigned int   n_components,
                                      const bool           dominate)
  : FiniteElement<dim, spacedim>(
      FiniteElementData<dim>(std::vector<unsigned>(dim + 1, 0),
                             type,
                             n_components,
                             0,
                             FiniteElementData<dim>::unknown),
      std::vector<bool>(),
      std::vector<ComponentMask>())
  , dominate(dominate)
{
  Assert(n_components >= 1,
         ExcMessage("A finite element needs to have at least one "
                    "vector component."));

  // in most other elements we have to set up all sorts of stuff
  // here. there isn't much that we have to do here; in particular,
  // we can simply leave the restriction and prolongation matrices
  // empty since their proper size is in fact zero given that the
  // element here has no degrees of freedom
}



template <int dim, int spacedim>
FE_Nothing<dim, spacedim>::FE_Nothing(const unsigned int n_components,
                                      const bool         dominate)
  : FE_Nothing<dim, spacedim>(ReferenceCells::get_hypercube<dim>(),
                              n_components,
                              dominate)
{}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_Nothing<dim, spacedim>::clone() const
{
  return std::make_unique<FE_Nothing<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
std::string
FE_Nothing<dim, spacedim>::get_name() const
{
  std::ostringstream namebuf;
  namebuf << "FE_Nothing<" << Utilities::dim_string(dim, spacedim) << ">(";

  std::vector<std::string> name_components;
  if (this->reference_cell() != ReferenceCells::get_hypercube<dim>())
    name_components.push_back(this->reference_cell().to_string());
  if (this->n_components() > 1)
    name_components.push_back(std::to_string(this->n_components()));
  if (dominate)
    name_components.emplace_back("dominating");

  for (const std::string &comp : name_components)
    {
      namebuf << comp;
      if (comp != name_components.back())
        namebuf << ", ";
    }
  namebuf << ")";

  return namebuf.str();
}



template <int dim, int spacedim>
UpdateFlags
FE_Nothing<dim, spacedim>::requires_update_flags(const UpdateFlags flags) const
{
  return flags;
}



template <int dim, int spacedim>
double
FE_Nothing<dim, spacedim>::shape_value(const unsigned int /*i*/,
                                       const Point<dim> & /*p*/) const
{
  Assert(false, ExcMessage("This element has no shape functions."));
  return 0;
}



template <int dim, int spacedim>
std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
FE_Nothing<dim, spacedim>::get_data(
  const UpdateFlags /*update_flags*/,
  const Mapping<dim, spacedim> & /*mapping*/,
  const Quadrature<dim> & /*quadrature*/,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    & /*output_data*/) const
{
  // Create a default data object.  Normally we would then
  // need to resize things to hold the appropriate numbers
  // of dofs, but in this case all data fields are empty.
  return std::make_unique<
    typename FiniteElement<dim, spacedim>::InternalDataBase>();
}



template <int dim, int spacedim>
void
FE_Nothing<dim, spacedim>::fill_fe_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &,
  const CellSimilarity::Similarity,
  const Quadrature<dim> &,
  const Mapping<dim, spacedim> &,
  const typename Mapping<dim, spacedim>::InternalDataBase &,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim> &,
  const typename FiniteElement<dim, spacedim>::InternalDataBase &,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &) const
{
  // leave data fields empty
}



template <int dim, int spacedim>
void
FE_Nothing<dim, spacedim>::fill_fe_face_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &,
  const unsigned int,
  const hp::QCollection<dim - 1> &,
  const Mapping<dim, spacedim> &,
  const typename Mapping<dim, spacedim>::InternalDataBase &,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim> &,
  const typename FiniteElement<dim, spacedim>::InternalDataBase &,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &) const
{
  // leave data fields empty
}



template <int dim, int spacedim>
void
FE_Nothing<dim, spacedim>::fill_fe_subface_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &,
  const unsigned int,
  const unsigned int,
  const Quadrature<dim - 1> &,
  const Mapping<dim, spacedim> &,
  const typename Mapping<dim, spacedim>::InternalDataBase &,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim> &,
  const typename FiniteElement<dim, spacedim>::InternalDataBase &,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &) const
{
  // leave data fields empty
}



template <int dim, int spacedim>
bool
FE_Nothing<dim, spacedim>::is_dominating() const
{
  return dominate;
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_Nothing<dim, spacedim>::compare_for_domination(
  const FiniteElement<dim, spacedim> &fe,
  const unsigned int                  codim) const
{
  Assert(codim <= dim, ExcImpossibleInDim(dim));

  if (!dominate)
    // if FE_Nothing does not dominate, there are no requirements
    return FiniteElementDomination::no_requirements;
  else if (dynamic_cast<const FE_Nothing<dim> *>(&fe) != nullptr)
    // if it does and the other is FE_Nothing, either can dominate
    return FiniteElementDomination::either_element_can_dominate;
  else
    // otherwise we dominate whatever FE is provided
    return FiniteElementDomination::this_element_dominates;
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_Nothing<dim, spacedim>::hp_vertex_dof_identities(
  const FiniteElement<dim, spacedim> & /*fe_other*/) const
{
  // the FE_Nothing has no
  // degrees of freedom, so there
  // are no equivalencies to be
  // recorded
  return std::vector<std::pair<unsigned int, unsigned int>>();
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_Nothing<dim, spacedim>::hp_line_dof_identities(
  const FiniteElement<dim, spacedim> & /*fe_other*/) const
{
  // the FE_Nothing has no
  // degrees of freedom, so there
  // are no equivalencies to be
  // recorded
  return std::vector<std::pair<unsigned int, unsigned int>>();
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_Nothing<dim, spacedim>::hp_quad_dof_identities(
  const FiniteElement<dim, spacedim> & /*fe_other*/,
  const unsigned int) const
{
  // the FE_Nothing has no
  // degrees of freedom, so there
  // are no equivalencies to be
  // recorded
  return std::vector<std::pair<unsigned int, unsigned int>>();
}


template <int dim, int spacedim>
bool
FE_Nothing<dim, spacedim>::hp_constraints_are_implemented() const
{
  return true;
}



template <int dim, int spacedim>
void
FE_Nothing<dim, spacedim>::get_interpolation_matrix(
  const FiniteElement<dim, spacedim> & /*source_fe*/,
  FullMatrix<double> &interpolation_matrix) const
{
  // Since this element has no dofs,
  // the interpolation matrix is necessarily empty.
  Assert(interpolation_matrix.m() == 0,
         ExcDimensionMismatch(interpolation_matrix.m(), 0));
  Assert(interpolation_matrix.n() == 0,
         ExcDimensionMismatch(interpolation_matrix.n(), 0));
}



template <int dim, int spacedim>
void
FE_Nothing<dim, spacedim>::get_face_interpolation_matrix(
  const FiniteElement<dim, spacedim> & /*source_fe*/,
  FullMatrix<double> &interpolation_matrix,
  const unsigned int) const
{
  // since this element has no face dofs, the
  // interpolation matrix is necessarily empty
  Assert(interpolation_matrix.m() == 0,
         ExcDimensionMismatch(interpolation_matrix.m(), 0));
  Assert(interpolation_matrix.n() == 0,
         ExcDimensionMismatch(interpolation_matrix.m(), 0));
}



template <int dim, int spacedim>
void
FE_Nothing<dim, spacedim>::get_subface_interpolation_matrix(
  const FiniteElement<dim, spacedim> & /*source_fe*/,
  const unsigned int /*index*/,
  FullMatrix<double> &interpolation_matrix,
  const unsigned int) const
{
  // since this element has no face dofs, the
  // interpolation matrix is necessarily empty
  Assert(interpolation_matrix.m() == 0,
         ExcDimensionMismatch(interpolation_matrix.m(), 0));
  Assert(interpolation_matrix.n() == 0,
         ExcDimensionMismatch(interpolation_matrix.m(), 0));
}



template <int dim, int spacedim>
std::pair<Table<2, bool>, std::vector<unsigned int>>
FE_Nothing<dim, spacedim>::get_constant_modes() const
{
  // since this element has no dofs, there are no constant modes
  return {Table<2, bool>{}, std::vector<unsigned int>{}};
}



// explicit instantiations
#include "fe/fe_nothing.inst"


DEAL_II_NAMESPACE_CLOSE

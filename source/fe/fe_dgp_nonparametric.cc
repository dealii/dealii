// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_dgp_nonparametric.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <memory>
#include <sstream>


DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
FE_DGPNonparametric<dim, spacedim>::FE_DGPNonparametric(
  const unsigned int degree)
  : FiniteElement<dim, spacedim>(
      FiniteElementData<dim>(get_dpo_vector(degree),
                             1,
                             degree,
                             FiniteElementData<dim>::L2),
      std::vector<bool>(
        FiniteElementData<dim>(get_dpo_vector(degree), 1, degree)
          .n_dofs_per_cell(),
        true),
      std::vector<ComponentMask>(
        FiniteElementData<dim>(get_dpo_vector(degree), 1, degree)
          .n_dofs_per_cell(),
        ComponentMask(std::vector<bool>(1, true))))
  , polynomial_space(Polynomials::Legendre::generate_complete_basis(degree))
{
  const unsigned int n_dofs = this->n_dofs_per_cell();
  for (const unsigned int ref_case :
       RefinementCase<dim>::all_refinement_cases())
    if (ref_case != RefinementCase<dim>::no_refinement)
      {
        if (dim != 2 && ref_case != RefinementCase<dim>::isotropic_refinement)
          // do nothing, as anisotropic
          // refinement is not
          // implemented so far
          continue;

        const unsigned int nc = this->reference_cell().template n_children<dim>(
          RefinementCase<dim>(ref_case));
        for (unsigned int i = 0; i < nc; ++i)
          {
            this->prolongation[ref_case - 1][i].reinit(n_dofs, n_dofs);
            // Fill prolongation matrices with
            // embedding operators
            for (unsigned int j = 0; j < n_dofs; ++j)
              this->prolongation[ref_case - 1][i](j, j) = 1.;
          }
      }

  // restriction can be defined
  // through projection for
  // discontinuous elements, but is
  // presently not implemented for DGPNonparametric
  // elements.
  //
  // if it were, then the following
  // snippet would be the right code
  //    if ((degree < Matrices::n_projection_matrices) &&
  //        (Matrices::projection_matrices[degree] != 0))
  //      {
  //        restriction[0].fill (Matrices::projection_matrices[degree]);
  //      }
  //    else
  //                                   // matrix undefined, set size to zero
  //      for (unsigned int i=0;i<GeometryInfo<dim>::max_children_per_cell;++i)
  //        restriction[i].reinit(0, 0);
  // since not implemented, set to
  // "empty". however, that is done in the
  // default constructor already, so do nothing
  //  for (unsigned int i=0;i<GeometryInfo<dim>::max_children_per_cell;++i)
  //    this->restriction[i].reinit(0, 0);

  // note further, that these
  // elements have neither support
  // nor face-support points, so
  // leave these fields empty
}



template <int dim, int spacedim>
std::string
FE_DGPNonparametric<dim, spacedim>::get_name() const
{
  // note that the
  // FETools::get_fe_by_name
  // function depends on the
  // particular format of the string
  // this function returns, so they
  // have to be kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_DGPNonparametric<" << Utilities::dim_string(dim, spacedim)
          << ">(" << this->degree << ")";

  return namebuf.str();
}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_DGPNonparametric<dim, spacedim>::clone() const
{
  return std::make_unique<FE_DGPNonparametric<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
double
FE_DGPNonparametric<dim, spacedim>::shape_value(const unsigned int i,
                                                const Point<dim>  &p) const
{
  (void)i;
  (void)p;
  AssertIndexRange(i, this->n_dofs_per_cell());
  AssertThrow(false,
              (typename FiniteElement<dim>::ExcUnitShapeValuesDoNotExist()));
  return 0;
}



template <int dim, int spacedim>
double
FE_DGPNonparametric<dim, spacedim>::shape_value_component(
  const unsigned int i,
  const Point<dim>  &p,
  const unsigned int component) const
{
  (void)i;
  (void)p;
  (void)component;
  AssertIndexRange(i, this->n_dofs_per_cell());
  AssertIndexRange(component, 1);
  AssertThrow(false,
              (typename FiniteElement<dim>::ExcUnitShapeValuesDoNotExist()));
  return 0;
}



template <int dim, int spacedim>
Tensor<1, dim>
FE_DGPNonparametric<dim, spacedim>::shape_grad(const unsigned int i,
                                               const Point<dim>  &p) const
{
  (void)i;
  (void)p;
  AssertIndexRange(i, this->n_dofs_per_cell());
  AssertThrow(false,
              (typename FiniteElement<dim>::ExcUnitShapeValuesDoNotExist()));
  return Tensor<1, dim>();
}


template <int dim, int spacedim>
Tensor<1, dim>
FE_DGPNonparametric<dim, spacedim>::shape_grad_component(
  const unsigned int i,
  const Point<dim>  &p,
  const unsigned int component) const
{
  (void)i;
  (void)p;
  (void)component;
  AssertIndexRange(i, this->n_dofs_per_cell());
  AssertIndexRange(component, 1);
  AssertThrow(false,
              (typename FiniteElement<dim>::ExcUnitShapeValuesDoNotExist()));
  return Tensor<1, dim>();
}



template <int dim, int spacedim>
Tensor<2, dim>
FE_DGPNonparametric<dim, spacedim>::shape_grad_grad(const unsigned int i,
                                                    const Point<dim>  &p) const
{
  (void)i;
  (void)p;
  AssertIndexRange(i, this->n_dofs_per_cell());
  AssertThrow(false,
              (typename FiniteElement<dim>::ExcUnitShapeValuesDoNotExist()));
  return Tensor<2, dim>();
}



template <int dim, int spacedim>
Tensor<2, dim>
FE_DGPNonparametric<dim, spacedim>::shape_grad_grad_component(
  const unsigned int i,
  const Point<dim>  &p,
  const unsigned int component) const
{
  (void)i;
  (void)p;
  (void)component;
  AssertIndexRange(i, this->n_dofs_per_cell());
  AssertIndexRange(component, 1);
  AssertThrow(false,
              (typename FiniteElement<dim>::ExcUnitShapeValuesDoNotExist()));
  return Tensor<2, dim>();
}


//---------------------------------------------------------------------------
// Auxiliary functions
//---------------------------------------------------------------------------


template <int dim, int spacedim>
std::vector<unsigned int>
FE_DGPNonparametric<dim, spacedim>::get_dpo_vector(const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim + 1, static_cast<unsigned int>(0));
  dpo[dim] = deg + 1;
  for (unsigned int i = 1; i < dim; ++i)
    {
      dpo[dim] *= deg + 1 + i;
      dpo[dim] /= i + 1;
    }
  return dpo;
}



template <int dim, int spacedim>
UpdateFlags
FE_DGPNonparametric<dim, spacedim>::requires_update_flags(
  const UpdateFlags flags) const
{
  UpdateFlags out = flags;

  if (flags & (update_values | update_gradients | update_hessians))
    out |= update_quadrature_points;

  return out;
}



//---------------------------------------------------------------------------
// Data field initialization
//---------------------------------------------------------------------------

template <int dim, int spacedim>
std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
FE_DGPNonparametric<dim, spacedim>::get_data(
  const UpdateFlags update_flags,
  const Mapping<dim, spacedim> &,
  const Quadrature<dim> &,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    & /*output_data*/) const
{
  // generate a new data object
  auto data_ptr =
    std::make_unique<typename FiniteElement<dim, spacedim>::InternalDataBase>();
  data_ptr->update_each = requires_update_flags(update_flags);

  // other than that, there is nothing we can add here as discussed
  // in the general documentation of this class

  return data_ptr;
}



//---------------------------------------------------------------------------
// Fill data of FEValues
//---------------------------------------------------------------------------

template <int dim, int spacedim>
void
FE_DGPNonparametric<dim, spacedim>::fill_fe_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &,
  const CellSimilarity::Similarity,
  const Quadrature<dim> &,
  const Mapping<dim, spacedim> &,
  const typename Mapping<dim, spacedim>::InternalDataBase &,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
                                                                &mapping_data,
  const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &output_data) const
{
  Assert(fe_internal.update_each & update_quadrature_points,
         ExcInternalError());

  const unsigned int n_q_points = mapping_data.quadrature_points.size();

  std::vector<double> values(
    (fe_internal.update_each & update_values) ? this->n_dofs_per_cell() : 0);
  std::vector<Tensor<1, dim>> grads(
    (fe_internal.update_each & update_gradients) ? this->n_dofs_per_cell() : 0);
  std::vector<Tensor<2, dim>> grad_grads(
    (fe_internal.update_each & update_hessians) ? this->n_dofs_per_cell() : 0);
  std::vector<Tensor<3, dim>> empty_vector_of_3rd_order_tensors;
  std::vector<Tensor<4, dim>> empty_vector_of_4th_order_tensors;

  if (fe_internal.update_each & (update_values | update_gradients))
    for (unsigned int i = 0; i < n_q_points; ++i)
      {
        polynomial_space.evaluate(mapping_data.quadrature_points[i],
                                  values,
                                  grads,
                                  grad_grads,
                                  empty_vector_of_3rd_order_tensors,
                                  empty_vector_of_4th_order_tensors);

        if (fe_internal.update_each & update_values)
          for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
            output_data.shape_values[k][i] = values[k];

        if (fe_internal.update_each & update_gradients)
          for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
            output_data.shape_gradients[k][i] = grads[k];

        if (fe_internal.update_each & update_hessians)
          for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
            output_data.shape_hessians[k][i] = grad_grads[k];
      }
}



template <int dim, int spacedim>
void
FE_DGPNonparametric<dim, spacedim>::fill_fe_face_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &,
  const unsigned int,
  const hp::QCollection<dim - 1> &,
  const Mapping<dim, spacedim> &,
  const typename Mapping<dim, spacedim>::InternalDataBase &,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
                                                                &mapping_data,
  const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &output_data) const
{
  Assert(fe_internal.update_each & update_quadrature_points,
         ExcInternalError());

  const unsigned int n_q_points = mapping_data.quadrature_points.size();

  std::vector<double> values(
    (fe_internal.update_each & update_values) ? this->n_dofs_per_cell() : 0);
  std::vector<Tensor<1, dim>> grads(
    (fe_internal.update_each & update_gradients) ? this->n_dofs_per_cell() : 0);
  std::vector<Tensor<2, dim>> grad_grads(
    (fe_internal.update_each & update_hessians) ? this->n_dofs_per_cell() : 0);
  std::vector<Tensor<3, dim>> empty_vector_of_3rd_order_tensors;
  std::vector<Tensor<4, dim>> empty_vector_of_4th_order_tensors;

  if (fe_internal.update_each & (update_values | update_gradients))
    for (unsigned int i = 0; i < n_q_points; ++i)
      {
        polynomial_space.evaluate(mapping_data.quadrature_points[i],
                                  values,
                                  grads,
                                  grad_grads,
                                  empty_vector_of_3rd_order_tensors,
                                  empty_vector_of_4th_order_tensors);

        if (fe_internal.update_each & update_values)
          for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
            output_data.shape_values[k][i] = values[k];

        if (fe_internal.update_each & update_gradients)
          for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
            output_data.shape_gradients[k][i] = grads[k];

        if (fe_internal.update_each & update_hessians)
          for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
            output_data.shape_hessians[k][i] = grad_grads[k];
      }
}



template <int dim, int spacedim>
void
FE_DGPNonparametric<dim, spacedim>::fill_fe_subface_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &,
  const unsigned int,
  const unsigned int,
  const Quadrature<dim - 1> &,
  const Mapping<dim, spacedim> &,
  const typename Mapping<dim, spacedim>::InternalDataBase &,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
                                                                &mapping_data,
  const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &output_data) const
{
  Assert(fe_internal.update_each & update_quadrature_points,
         ExcInternalError());

  const unsigned int n_q_points = mapping_data.quadrature_points.size();

  std::vector<double> values(
    (fe_internal.update_each & update_values) ? this->n_dofs_per_cell() : 0);
  std::vector<Tensor<1, dim>> grads(
    (fe_internal.update_each & update_gradients) ? this->n_dofs_per_cell() : 0);
  std::vector<Tensor<2, dim>> grad_grads(
    (fe_internal.update_each & update_hessians) ? this->n_dofs_per_cell() : 0);
  std::vector<Tensor<3, dim>> empty_vector_of_3rd_order_tensors;
  std::vector<Tensor<4, dim>> empty_vector_of_4th_order_tensors;

  if (fe_internal.update_each & (update_values | update_gradients))
    for (unsigned int i = 0; i < n_q_points; ++i)
      {
        polynomial_space.evaluate(mapping_data.quadrature_points[i],
                                  values,
                                  grads,
                                  grad_grads,
                                  empty_vector_of_3rd_order_tensors,
                                  empty_vector_of_4th_order_tensors);

        if (fe_internal.update_each & update_values)
          for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
            output_data.shape_values[k][i] = values[k];

        if (fe_internal.update_each & update_gradients)
          for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
            output_data.shape_gradients[k][i] = grads[k];

        if (fe_internal.update_each & update_hessians)
          for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
            output_data.shape_hessians[k][i] = grad_grads[k];
      }
}



template <int dim, int spacedim>
void
FE_DGPNonparametric<dim, spacedim>::get_face_interpolation_matrix(
  const FiniteElement<dim, spacedim> &x_source_fe,
  FullMatrix<double>                 &interpolation_matrix,
  const unsigned int) const
{
  // this is only implemented, if the source
  // FE is also a DGPNonparametric element. in that case,
  // both elements have no dofs on their
  // faces and the face interpolation matrix
  // is necessarily empty -- i.e. there isn't
  // much we need to do here.
  (void)interpolation_matrix;
  AssertThrow(
    (x_source_fe.get_name().find("FE_DGPNonparametric<") == 0) ||
      (dynamic_cast<const FE_DGPNonparametric<dim, spacedim> *>(&x_source_fe) !=
       nullptr),
    (typename FiniteElement<dim, spacedim>::ExcInterpolationNotImplemented()));

  Assert(interpolation_matrix.m() == 0,
         ExcDimensionMismatch(interpolation_matrix.m(), 0));
  Assert(interpolation_matrix.n() == 0,
         ExcDimensionMismatch(interpolation_matrix.n(), 0));
}



template <int dim, int spacedim>
void
FE_DGPNonparametric<dim, spacedim>::get_subface_interpolation_matrix(
  const FiniteElement<dim, spacedim> &x_source_fe,
  const unsigned int,
  FullMatrix<double> &interpolation_matrix,
  const unsigned int) const
{
  // this is only implemented, if the source
  // FE is also a DGPNonparametric element. in that case,
  // both elements have no dofs on their
  // faces and the face interpolation matrix
  // is necessarily empty -- i.e. there isn't
  // much we need to do here.
  (void)interpolation_matrix;
  AssertThrow(
    (x_source_fe.get_name().find("FE_DGPNonparametric<") == 0) ||
      (dynamic_cast<const FE_DGPNonparametric<dim, spacedim> *>(&x_source_fe) !=
       nullptr),
    (typename FiniteElement<dim, spacedim>::ExcInterpolationNotImplemented()));

  Assert(interpolation_matrix.m() == 0,
         ExcDimensionMismatch(interpolation_matrix.m(), 0));
  Assert(interpolation_matrix.n() == 0,
         ExcDimensionMismatch(interpolation_matrix.n(), 0));
}



template <int dim, int spacedim>
bool
FE_DGPNonparametric<dim, spacedim>::hp_constraints_are_implemented() const
{
  return true;
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_DGPNonparametric<dim, spacedim>::hp_vertex_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  // there are no such constraints for DGPNonparametric
  // elements at all
  if (dynamic_cast<const FE_DGPNonparametric<dim, spacedim> *>(&fe_other) !=
      nullptr)
    return std::vector<std::pair<unsigned int, unsigned int>>();
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_DGPNonparametric<dim, spacedim>::hp_line_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  // there are no such constraints for DGPNonparametric
  // elements at all
  if (dynamic_cast<const FE_DGPNonparametric<dim, spacedim> *>(&fe_other) !=
      nullptr)
    return std::vector<std::pair<unsigned int, unsigned int>>();
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_DGPNonparametric<dim, spacedim>::hp_quad_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int) const
{
  // there are no such constraints for DGPNonparametric
  // elements at all
  if (dynamic_cast<const FE_DGPNonparametric<dim, spacedim> *>(&fe_other) !=
      nullptr)
    return std::vector<std::pair<unsigned int, unsigned int>>();
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_DGPNonparametric<dim, spacedim>::compare_for_domination(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  codim) const
{
  Assert(codim <= dim, ExcImpossibleInDim(dim));

  // vertex/line/face domination
  // ---------------------------
  if (codim > 0)
    // this is a discontinuous element, so by definition there will
    // be no constraints wherever this element comes together with
    // any other kind of element
    return FiniteElementDomination::no_requirements;

  // cell domination
  // ---------------
  if (const FE_DGPNonparametric<dim, spacedim> *fe_nonparametric_other =
        dynamic_cast<const FE_DGPNonparametric<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_nonparametric_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_nonparametric_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_Nothing<dim> *fe_nothing =
             dynamic_cast<const FE_Nothing<dim> *>(&fe_other))
    {
      if (fe_nothing->is_dominating())
        return FiniteElementDomination::other_element_dominates;
      else
        // the FE_Nothing has no degrees of freedom and it is typically used
        // in a context where we don't require any continuity along the
        // interface
        return FiniteElementDomination::no_requirements;
    }

  DEAL_II_NOT_IMPLEMENTED();
  return FiniteElementDomination::neither_element_dominates;
}



template <int dim, int spacedim>
bool
FE_DGPNonparametric<dim, spacedim>::has_support_on_face(
  const unsigned int,
  const unsigned int) const
{
  return true;
}



template <int dim, int spacedim>
std::size_t
FE_DGPNonparametric<dim, spacedim>::memory_consumption() const
{
  DEAL_II_NOT_IMPLEMENTED();
  return 0;
}



template <int dim, int spacedim>
unsigned int
FE_DGPNonparametric<dim, spacedim>::get_degree() const
{
  return this->degree;
}



// explicit instantiations
#include "fe/fe_dgp_nonparametric.inst"


DEAL_II_NAMESPACE_CLOSE

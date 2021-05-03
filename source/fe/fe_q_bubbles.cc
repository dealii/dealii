// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q_bubbles.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <memory>
#include <sstream>
#include <vector>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace FE_Q_Bubbles
  {
    namespace
    {
      template <int dim, int spacedim>
      inline void
      compute_embedding_matrices(
        const dealii::FE_Q_Bubbles<dim, spacedim> &   fe,
        std::vector<std::vector<FullMatrix<double>>> &matrices,
        const bool                                    isotropic_only)
      {
        const unsigned int dpc    = fe.n_dofs_per_cell();
        const unsigned int degree = fe.degree;

        // Initialize quadrature formula on fine cells
        std::unique_ptr<Quadrature<dim>> q_fine;
        Quadrature<1>                    q_dummy(std::vector<Point<1>>(1),
                              std::vector<double>(1, 1.));
        switch (dim)
          {
            case 1:
              if (spacedim == 1)
                q_fine = std::make_unique<QGauss<dim>>(degree + 1);
              else if (spacedim == 2)
                q_fine =
                  std::make_unique<QAnisotropic<dim>>(QGauss<1>(degree + 1),
                                                      q_dummy);
              else
                q_fine =
                  std::make_unique<QAnisotropic<dim>>(QGauss<1>(degree + 1),
                                                      q_dummy,
                                                      q_dummy);
              break;
            case 2:
              if (spacedim == 2)
                q_fine = std::make_unique<QGauss<dim>>(degree + 1);
              else
                q_fine =
                  std::make_unique<QAnisotropic<dim>>(QGauss<1>(degree + 1),
                                                      QGauss<1>(degree + 1),
                                                      q_dummy);
              break;
            case 3:
              q_fine = std::make_unique<QGauss<dim>>(degree + 1);
              break;
            default:
              Assert(false, ExcInternalError());
          }

        Assert(q_fine.get() != nullptr, ExcInternalError());
        const unsigned int nq = q_fine->size();

        // loop over all possible refinement cases
        unsigned int ref_case = (isotropic_only) ?
                                  RefinementCase<dim>::isotropic_refinement :
                                  RefinementCase<dim>::cut_x;
        for (; ref_case <= RefinementCase<dim>::isotropic_refinement;
             ++ref_case)
          {
            const unsigned int nc =
              GeometryInfo<dim>::n_children(RefinementCase<dim>(ref_case));

            for (unsigned int i = 0; i < nc; ++i)
              {
                Assert(matrices[ref_case - 1][i].n() == dpc,
                       ExcDimensionMismatch(matrices[ref_case - 1][i].n(),
                                            dpc));
                Assert(matrices[ref_case - 1][i].m() == dpc,
                       ExcDimensionMismatch(matrices[ref_case - 1][i].m(),
                                            dpc));
              }

            // create a respective refinement on the triangulation
            dealii::Triangulation<dim, spacedim> tr;
            GridGenerator::hyper_cube(tr, 0, 1);
            tr.begin_active()->set_refine_flag(RefinementCase<dim>(ref_case));
            tr.execute_coarsening_and_refinement();

            dealii::DoFHandler<dim, spacedim> dh(tr);
            dh.distribute_dofs(fe);

            dealii::FEValues<dim, spacedim> fine(get_default_linear_mapping(tr),
                                                 fe,
                                                 *q_fine,
                                                 update_quadrature_points |
                                                   update_JxW_values |
                                                   update_values);

            const unsigned int n_dofs = dh.n_dofs();

            FullMatrix<double> fine_mass(n_dofs);
            FullMatrix<double> coarse_rhs_matrix(n_dofs, dpc);

            std::vector<std::vector<types::global_dof_index>> child_ldi(
              nc, std::vector<types::global_dof_index>(fe.n_dofs_per_cell()));

            // now create the mass matrix and all the right_hand sides
            unsigned int                                           child_no = 0;
            typename dealii::DoFHandler<dim>::active_cell_iterator cell =
              dh.begin_active();
            for (; cell != dh.end(); ++cell, ++child_no)
              {
                fine.reinit(cell);
                cell->get_dof_indices(child_ldi[child_no]);

                for (unsigned int q = 0; q < nq; ++q)
                  for (unsigned int i = 0; i < dpc; ++i)
                    for (unsigned int j = 0; j < dpc; ++j)
                      {
                        const unsigned int gdi = child_ldi[child_no][i];
                        const unsigned int gdj = child_ldi[child_no][j];
                        fine_mass(gdi, gdj) += fine.shape_value(i, q) *
                                               fine.shape_value(j, q) *
                                               fine.JxW(q);
                        Point<dim> quad_tmp;
                        for (unsigned int k = 0; k < dim; ++k)
                          quad_tmp(k) = fine.quadrature_point(q)(k);
                        coarse_rhs_matrix(gdi, j) +=
                          fine.shape_value(i, q) * fe.shape_value(j, quad_tmp) *
                          fine.JxW(q);
                      }
              }

            // now solve for all right-hand sides simultaneously
            dealii::FullMatrix<double> solution(n_dofs, dpc);
            fine_mass.gauss_jordan();
            fine_mass.mmult(solution, coarse_rhs_matrix);

            // and distribute to the fine cell matrices
            for (unsigned int child_no = 0; child_no < nc; ++child_no)
              for (unsigned int i = 0; i < dpc; ++i)
                for (unsigned int j = 0; j < dpc; ++j)
                  {
                    const unsigned int gdi = child_ldi[child_no][i];
                    // remove small entries
                    if (std::fabs(solution(gdi, j)) > 1.e-12)
                      matrices[ref_case - 1][child_no](i, j) = solution(gdi, j);
                  }
          }
      }
    } // namespace
  }   // namespace FE_Q_Bubbles
} // namespace internal


template <int dim, int spacedim>
FE_Q_Bubbles<dim, spacedim>::FE_Q_Bubbles(const unsigned int q_degree)
  : FE_Q_Base<dim, spacedim>(TensorProductPolynomialsBubbles<dim>(
                               Polynomials::generate_complete_Lagrange_basis(
                                 QGaussLobatto<1>(q_degree + 1).get_points())),
                             FiniteElementData<dim>(get_dpo_vector(q_degree),
                                                    1,
                                                    q_degree + 1,
                                                    FiniteElementData<dim>::H1),
                             get_riaf_vector(q_degree))
  , n_bubbles((q_degree <= 1) ? 1 : dim)
{
  Assert(q_degree > 0,
         ExcMessage("This element can only be used for polynomial degrees "
                    "greater than zero"));

  this->initialize(QGaussLobatto<1>(q_degree + 1).get_points());

  // adjust unit support point for discontinuous node
  Point<dim> point;
  for (unsigned int d = 0; d < dim; ++d)
    point[d] = 0.5;
  for (unsigned int i = 0; i < n_bubbles; ++i)
    this->unit_support_points.push_back(point);
  AssertDimension(this->n_dofs_per_cell(), this->unit_support_points.size());

  this->reinit_restriction_and_prolongation_matrices();
  if (dim == spacedim)
    {
      internal::FE_Q_Bubbles::compute_embedding_matrices(*this,
                                                         this->prolongation,
                                                         false);
      // Fill restriction matrices with L2-projection
      FETools::compute_projection_matrices(*this, this->restriction);
    }
}



template <int dim, int spacedim>
FE_Q_Bubbles<dim, spacedim>::FE_Q_Bubbles(const Quadrature<1> &points)
  : FE_Q_Base<dim, spacedim>(
      TensorProductPolynomialsBubbles<dim>(
        Polynomials::generate_complete_Lagrange_basis(points.get_points())),
      FiniteElementData<dim>(get_dpo_vector(points.size() - 1),
                             1,
                             points.size(),
                             FiniteElementData<dim>::H1),
      get_riaf_vector(points.size() - 1))
  , n_bubbles((points.size() - 1 <= 1) ? 1 : dim)
{
  Assert(points.size() > 1,
         ExcMessage("This element can only be used for polynomial degrees "
                    "at least one"));

  this->initialize(points.get_points());

  // adjust unit support point for discontinuous node
  Point<dim> point;
  for (unsigned int d = 0; d < dim; ++d)
    point[d] = 0.5;
  for (unsigned int i = 0; i < n_bubbles; ++i)
    this->unit_support_points.push_back(point);
  AssertDimension(this->n_dofs_per_cell(), this->unit_support_points.size());

  this->reinit_restriction_and_prolongation_matrices();
  if (dim == spacedim)
    {
      internal::FE_Q_Bubbles::compute_embedding_matrices(*this,
                                                         this->prolongation,
                                                         false);
      // Fill restriction matrices with L2-projection
      FETools::compute_projection_matrices(*this, this->restriction);
    }
}



template <int dim, int spacedim>
std::string
FE_Q_Bubbles<dim, spacedim>::get_name() const
{
  // note that the FETools::get_fe_by_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch

  std::ostringstream             namebuf;
  bool                           type     = true;
  const unsigned int             n_points = this->degree;
  std::vector<double>            points(n_points);
  const unsigned int             dofs_per_cell = this->n_dofs_per_cell();
  const std::vector<Point<dim>> &unit_support_points =
    this->unit_support_points;
  unsigned int index = 0;

  // Decode the support points in one coordinate direction.
  for (unsigned int j = 0; j < dofs_per_cell; j++)
    {
      if ((dim > 1) ? (unit_support_points[j](1) == 0 &&
                       ((dim > 2) ? unit_support_points[j](2) == 0 : true)) :
                      true)
        {
          if (index == 0)
            points[index] = unit_support_points[j](0);
          else if (index == 1)
            points[n_points - 1] = unit_support_points[j](0);
          else
            points[index - 1] = unit_support_points[j](0);

          index++;
        }
    }
  // Do not consider the discontinuous node for dimension 1
  Assert(index == n_points || (dim == 1 && index == n_points + n_bubbles),
         ExcMessage(
           "Could not decode support points in one coordinate direction."));

  // Check whether the support points are equidistant.
  for (unsigned int j = 0; j < n_points; j++)
    if (std::fabs(points[j] - static_cast<double>(j) / (this->degree - 1)) >
        1e-15)
      {
        type = false;
        break;
      }

  if (type == true)
    {
      if (this->degree > 3)
        namebuf << "FE_Q_Bubbles<" << Utilities::dim_string(dim, spacedim)
                << ">(QIterated(QTrapezoid()," << this->degree - 1 << "))";
      else
        namebuf << "FE_Q_Bubbles<" << Utilities::dim_string(dim, spacedim)
                << ">(" << this->degree - 1 << ")";
    }
  else
    {
      // Check whether the support points come from QGaussLobatto.
      const QGaussLobatto<1> points_gl(n_points);
      type = true;
      for (unsigned int j = 0; j < n_points; j++)
        if (points[j] != points_gl.point(j)(0))
          {
            type = false;
            break;
          }
      if (type == true)
        namebuf << "FE_Q_Bubbles<" << dim << ">(" << this->degree - 1 << ")";
      else
        namebuf << "FE_Q_Bubbles<" << dim << ">(QUnknownNodes(" << this->degree
                << "))";
    }
  return namebuf.str();
}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_Q_Bubbles<dim, spacedim>::clone() const
{
  return std::make_unique<FE_Q_Bubbles<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
void
FE_Q_Bubbles<dim, spacedim>::
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const
{
  Assert(support_point_values.size() == this->unit_support_points.size(),
         ExcDimensionMismatch(support_point_values.size(),
                              this->unit_support_points.size()));
  Assert(nodal_values.size() == this->n_dofs_per_cell(),
         ExcDimensionMismatch(nodal_values.size(), this->n_dofs_per_cell()));
  Assert(support_point_values[0].size() == this->n_components(),
         ExcDimensionMismatch(support_point_values[0].size(),
                              this->n_components()));

  for (unsigned int i = 0; i < this->n_dofs_per_cell() - 1; ++i)
    {
      const std::pair<unsigned int, unsigned int> index =
        this->system_to_component_index(i);
      nodal_values[i] = support_point_values[i](index.first);
    }

  // We don't use the bubble functions for local interpolation
  for (unsigned int i = 0; i < n_bubbles; ++i)
    nodal_values[nodal_values.size() - i - 1] = 0.;
}



template <int dim, int spacedim>
void
FE_Q_Bubbles<dim, spacedim>::get_interpolation_matrix(
  const FiniteElement<dim, spacedim> &x_source_fe,
  FullMatrix<double> &                interpolation_matrix) const
{
  // We don't know how to do this properly, yet.
  // However, for SolutionTransfer to work we need to provide an implementation
  // for the case that the x_source_fe is identical to this FE
  using FEQBUBBLES = FE_Q_Bubbles<dim, spacedim>;

  AssertThrow(
    (x_source_fe.get_name().find("FE_Q_Bubbles<") == 0) ||
      (dynamic_cast<const FEQBUBBLES *>(&x_source_fe) != nullptr),
    (typename FiniteElement<dim, spacedim>::ExcInterpolationNotImplemented()));
  Assert(interpolation_matrix.m() == this->n_dofs_per_cell(),
         ExcDimensionMismatch(interpolation_matrix.m(),
                              this->n_dofs_per_cell()));
  Assert(interpolation_matrix.n() == x_source_fe.n_dofs_per_cell(),
         ExcDimensionMismatch(interpolation_matrix.m(),
                              x_source_fe.n_dofs_per_cell()));

  // Provide a short cut in case we are just inquiring the identity
  auto casted_fe = dynamic_cast<const FEQBUBBLES *>(&x_source_fe);
  if (casted_fe != nullptr && casted_fe->degree == this->degree)
    for (unsigned int i = 0; i < interpolation_matrix.m(); ++i)
      interpolation_matrix.set(i, i, 1.);
  // else we need to do more...
  else
    Assert(
      false,
      (typename FiniteElement<dim,
                              spacedim>::ExcInterpolationNotImplemented()));
}



template <int dim, int spacedim>
std::vector<bool>
FE_Q_Bubbles<dim, spacedim>::get_riaf_vector(const unsigned int q_deg)
{
  const unsigned int n_cont_dofs = Utilities::fixed_power<dim>(q_deg + 1);
  const unsigned int n_bubbles   = (q_deg <= 1 ? 1 : dim);
  return std::vector<bool>(n_cont_dofs + n_bubbles, true);
}



template <int dim, int spacedim>
std::vector<unsigned int>
FE_Q_Bubbles<dim, spacedim>::get_dpo_vector(const unsigned int q_deg)
{
  std::vector<unsigned int> dpo(dim + 1, 1U);
  for (unsigned int i = 1; i < dpo.size(); ++i)
    dpo[i] = dpo[i - 1] * (q_deg - 1);

  // Then add the bubble functions; they are all associated with the
  // cell interior
  dpo[dim] += (q_deg <= 1 ? 1 : dim);
  return dpo;
}



template <int dim, int spacedim>
bool
FE_Q_Bubbles<dim, spacedim>::has_support_on_face(
  const unsigned int shape_index,
  const unsigned int face_index) const
{
  // discontinuous functions have no support on faces
  if (shape_index >= this->n_dofs_per_cell() - n_bubbles)
    return false;
  else
    return FE_Q_Base<dim, spacedim>::has_support_on_face(shape_index,
                                                         face_index);
}



template <int dim, int spacedim>
const FullMatrix<double> &
FE_Q_Bubbles<dim, spacedim>::get_prolongation_matrix(
  const unsigned int         child,
  const RefinementCase<dim> &refinement_case) const
{
  AssertIndexRange(refinement_case,
                   RefinementCase<dim>::isotropic_refinement + 1);
  Assert(refinement_case != RefinementCase<dim>::no_refinement,
         ExcMessage(
           "Prolongation matrices are only available for refined cells!"));
  AssertIndexRange(child, GeometryInfo<dim>::n_children(refinement_case));

  Assert(this->prolongation[refinement_case - 1][child].n() != 0,
         ExcMessage("This prolongation matrix has not been computed yet!"));
  // finally return the matrix
  return this->prolongation[refinement_case - 1][child];
}



template <int dim, int spacedim>
const FullMatrix<double> &
FE_Q_Bubbles<dim, spacedim>::get_restriction_matrix(
  const unsigned int         child,
  const RefinementCase<dim> &refinement_case) const
{
  AssertIndexRange(refinement_case,
                   RefinementCase<dim>::isotropic_refinement + 1);
  Assert(refinement_case != RefinementCase<dim>::no_refinement,
         ExcMessage(
           "Restriction matrices are only available for refined cells!"));
  AssertIndexRange(child, GeometryInfo<dim>::n_children(refinement_case));

  Assert(this->restriction[refinement_case - 1][child].n() != 0,
         ExcMessage("This restriction matrix has not been computed yet!"));

  // finally return the matrix
  return this->restriction[refinement_case - 1][child];
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_Q_Bubbles<dim, spacedim>::compare_for_domination(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  codim) const
{
  Assert(codim <= dim, ExcImpossibleInDim(dim));

  // vertex/line/face domination
  // (if fe_other is derived from FE_DGQ)
  // ------------------------------------
  if (codim > 0)
    if (dynamic_cast<const FE_DGQ<dim, spacedim> *>(&fe_other) != nullptr)
      // there are no requirements between continuous and discontinuous elements
      return FiniteElementDomination::no_requirements;

  // vertex/line/face domination
  // (if fe_other is not derived from FE_DGQ)
  // & cell domination
  // ----------------------------------------
  if (const FE_Q_Bubbles<dim, spacedim> *fe_bubbles_other =
        dynamic_cast<const FE_Q_Bubbles<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_bubbles_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_bubbles_other->degree)
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

  Assert(false, ExcNotImplemented());
  return FiniteElementDomination::neither_element_dominates;
}


// explicit instantiations
#include "fe_q_bubbles.inst"

DEAL_II_NAMESPACE_CLOSE

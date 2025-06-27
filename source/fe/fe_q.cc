// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_hermite.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_wedge_p.h>

#include <deal.II/lac/vector.h>

#include <memory>
#include <sstream>
#include <vector>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace FE_Q
  {
    namespace
    {
      std::vector<Point<1>>
      get_QGaussLobatto_points(const unsigned int degree)
      {
        if (degree > 0)
          return QGaussLobatto<1>(degree + 1).get_points();
        else
          {
            using FEQ = dealii::FE_Q_Base<1, 1>;
            AssertThrow(false, FEQ::ExcFEQCannotHaveDegree0());
          }
        return std::vector<Point<1>>();
      }
    } // namespace
  }   // namespace FE_Q
} // namespace internal



template <int dim, int spacedim>
FE_Q<dim, spacedim>::FE_Q(const unsigned int degree)
  : FE_Q_Base<dim, spacedim>(
      TensorProductPolynomials<dim>(
        Polynomials::generate_complete_Lagrange_basis(
          internal::FE_Q::get_QGaussLobatto_points(degree))),
      FiniteElementData<dim>(this->get_dpo_vector(degree),
                             1,
                             degree,
                             FiniteElementData<dim>::H1),
      std::vector<bool>(1, false))
{
  this->initialize(internal::FE_Q::get_QGaussLobatto_points(degree));
}



template <int dim, int spacedim>
FE_Q<dim, spacedim>::FE_Q(const Quadrature<1> &points)
  : FE_Q_Base<dim, spacedim>(
      TensorProductPolynomials<dim>(
        Polynomials::generate_complete_Lagrange_basis(points.get_points())),
      FiniteElementData<dim>(this->get_dpo_vector(points.size() - 1),
                             1,
                             points.size() - 1,
                             FiniteElementData<dim>::H1),
      std::vector<bool>(1, false))
{
  this->initialize(points.get_points());
}



template <int dim, int spacedim>
std::string
FE_Q<dim, spacedim>::get_name() const
{
  // note that the FETools::get_fe_by_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch

  std::ostringstream  namebuf;
  bool                equidistant = true;
  std::vector<double> points(this->degree + 1);

  // Decode the support points in one coordinate direction.
  TensorProductPolynomials<dim> *poly_space_derived_ptr =
    dynamic_cast<TensorProductPolynomials<dim> *>(this->poly_space.get());
  std::vector<unsigned int> lexicographic =
    poly_space_derived_ptr->get_numbering_inverse();
  for (unsigned int j = 0; j <= this->degree; ++j)
    points[j] = this->unit_support_points[lexicographic[j]][0];

  // Check whether the support points are equidistant.
  for (unsigned int j = 0; j <= this->degree; ++j)
    if (std::fabs(points[j] - static_cast<double>(j) / this->degree) > 1e-15)
      {
        equidistant = false;
        break;
      }

  if (equidistant == true)
    {
      if (this->degree > 2)
        namebuf << "FE_Q<" << Utilities::dim_string(dim, spacedim)
                << ">(QIterated(QTrapezoid()," << this->degree << "))";
      else
        namebuf << "FE_Q<" << Utilities::dim_string(dim, spacedim) << ">("
                << this->degree << ")";
    }
  else
    {
      // Check whether the support points come from QGaussLobatto.
      const QGaussLobatto<1> points_gl(this->degree + 1);
      bool                   gauss_lobatto = true;
      for (unsigned int j = 0; j <= this->degree; ++j)
        if (points[j] != points_gl.point(j)[0])
          {
            gauss_lobatto = false;
            break;
          }
      if (gauss_lobatto == true)
        namebuf << "FE_Q<" << Utilities::dim_string(dim, spacedim) << ">("
                << this->degree << ")";
      else
        namebuf << "FE_Q<" << Utilities::dim_string(dim, spacedim)
                << ">(QUnknownNodes(" << this->degree << "))";
    }
  return namebuf.str();
}



template <int dim, int spacedim>
void
FE_Q<dim, spacedim>::convert_generalized_support_point_values_to_dof_values(
  const std::vector<Vector<double>> &support_point_values,
  std::vector<double>               &nodal_values) const
{
  AssertDimension(support_point_values.size(),
                  this->get_unit_support_points().size());
  AssertDimension(support_point_values.size(), nodal_values.size());
  AssertDimension(this->n_dofs_per_cell(), nodal_values.size());

  for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
    {
      AssertDimension(support_point_values[i].size(), 1);

      nodal_values[i] = support_point_values[i](0);
    }
}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_Q<dim, spacedim>::clone() const
{
  return std::make_unique<FE_Q<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_Q<dim, spacedim>::compare_for_domination(
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
  if (const FE_Q<dim, spacedim> *fe_q_other =
        dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_q_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_q_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_SimplexP<dim, spacedim> *fe_p_other =
             dynamic_cast<const FE_SimplexP<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_p_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_p_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_WedgeP<dim, spacedim> *fe_wp_other =
             dynamic_cast<const FE_WedgeP<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_wp_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_wp_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_PyramidP<dim, spacedim> *fe_pp_other =
             dynamic_cast<const FE_PyramidP<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_pp_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_pp_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_Nothing<dim, spacedim> *fe_nothing =
             dynamic_cast<const FE_Nothing<dim, spacedim> *>(&fe_other))
    {
      if (fe_nothing->is_dominating())
        return FiniteElementDomination::other_element_dominates;
      else
        // the FE_Nothing has no degrees of freedom and it is typically used
        // in a context where we don't require any continuity along the
        // interface
        return FiniteElementDomination::no_requirements;
    }
  else if (const FE_Hermite<dim, spacedim> *fe_hermite_other =
             dynamic_cast<const FE_Hermite<dim, spacedim> *>(&fe_other))
    {
      if (this->degree == 1)
        {
          if (fe_hermite_other->degree > 1)
            return FiniteElementDomination::this_element_dominates;
          else
            return FiniteElementDomination::either_element_can_dominate;
        }
      else if (this->degree >= fe_hermite_other->degree)
        return FiniteElementDomination::other_element_dominates;
      else
        return FiniteElementDomination::neither_element_dominates;
    }

  DEAL_II_NOT_IMPLEMENTED();
  return FiniteElementDomination::neither_element_dominates;
}


// explicit instantiations
#include "fe/fe_q.inst"

DEAL_II_NAMESPACE_CLOSE

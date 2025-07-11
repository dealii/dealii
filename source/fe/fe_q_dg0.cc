// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2025 by the deal.II authors
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
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q_dg0.h>

#include <memory>
#include <sstream>
#include <vector>

DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim>
FE_Q_DG0<dim, spacedim>::FE_Q_DG0(const unsigned int degree)
  : FE_Q_Base<dim, spacedim>(TensorProductPolynomialsConst<dim>(
                               Polynomials::generate_complete_Lagrange_basis(
                                 QGaussLobatto<1>(degree + 1).get_points())),
                             FiniteElementData<dim>(get_dpo_vector(degree),
                                                    1,
                                                    degree,
                                                    FiniteElementData<dim>::L2),
                             get_riaf_vector(degree))
{
  Assert(degree > 0,
         ExcMessage("This element can only be used for polynomial degrees "
                    "greater than zero"));

  this->initialize(QGaussLobatto<1>(degree + 1).get_points());

  // adjust unit support point for discontinuous node
  Point<dim> point;
  for (unsigned int d = 0; d < dim; ++d)
    point[d] = 0.5;
  this->unit_support_points.push_back(point);
  AssertDimension(this->n_dofs_per_cell(), this->unit_support_points.size());
}



template <int dim, int spacedim>
FE_Q_DG0<dim, spacedim>::FE_Q_DG0(const Quadrature<1> &points)
  : FE_Q_Base<dim, spacedim>(
      TensorProductPolynomialsConst<dim>(
        Polynomials::generate_complete_Lagrange_basis(points.get_points())),
      FiniteElementData<dim>(get_dpo_vector(points.size() - 1),
                             1,
                             points.size() - 1,
                             FiniteElementData<dim>::L2),
      get_riaf_vector(points.size() - 1))
{
  const int degree = points.size() - 1;
  Assert(degree > 0,
         ExcMessage("This element can only be used for polynomial degrees "
                    "at least zero"));

  this->initialize(points.get_points());

  // adjust unit support point for discontinuous node
  Point<dim> point;
  for (unsigned int d = 0; d < dim; ++d)
    point[d] = 0.5;
  this->unit_support_points.push_back(point);
  AssertDimension(this->n_dofs_per_cell(), this->unit_support_points.size());
}



template <int dim, int spacedim>
std::string
FE_Q_DG0<dim, spacedim>::get_name() const
{
  // note that the FETools::get_fe_by_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch

  std::ostringstream             namebuf;
  bool                           type     = true;
  const unsigned int             n_points = this->degree + 1;
  std::vector<double>            points(n_points);
  const unsigned int             dofs_per_cell = this->n_dofs_per_cell();
  const std::vector<Point<dim>> &unit_support_points =
    this->unit_support_points;
  unsigned int index = 0;

  // Decode the support points in one coordinate direction.
  for (unsigned int j = 0; j < dofs_per_cell; ++j)
    {
      if ((dim > 1) ? (unit_support_points[j][1] == 0 &&
                       ((dim > 2) ? unit_support_points[j][2] == 0 : true)) :
                      true)
        {
          if (index == 0)
            points[index] = unit_support_points[j][0];
          else if (index == 1)
            points[n_points - 1] = unit_support_points[j][0];
          else
            points[index - 1] = unit_support_points[j][0];

          ++index;
        }
    }
  // Do not consider the discontinuous node for dimension 1
  Assert(index == n_points || (dim == 1 && index == n_points + 1),
         ExcMessage(
           "Could not decode support points in one coordinate direction."));

  // Check whether the support points are equidistant.
  for (unsigned int j = 0; j < n_points; ++j)
    if (std::fabs(points[j] - static_cast<double>(j) / this->degree) > 1e-15)
      {
        type = false;
        break;
      }

  if (type == true)
    {
      if (this->degree > 2)
        namebuf << "FE_Q_DG0<" << Utilities::dim_string(dim, spacedim)
                << ">(QIterated(QTrapezoid()," << this->degree << "))";
      else
        namebuf << "FE_Q_DG0<" << Utilities::dim_string(dim, spacedim) << ">("
                << this->degree << ")";
    }
  else
    {
      // Check whether the support points come from QGaussLobatto.
      const QGaussLobatto<1> points_gl(n_points);
      type = true;
      for (unsigned int j = 0; j < n_points; ++j)
        if (points[j] != points_gl.point(j)[0])
          {
            type = false;
            break;
          }
      if (type == true)
        namebuf << "FE_Q_DG0<" << Utilities::dim_string(dim, spacedim) << ">("
                << this->degree << ")";
      else
        namebuf << "FE_Q_DG0<" << Utilities::dim_string(dim, spacedim)
                << ">(QUnknownNodes(" << this->degree << "))";
    }
  return namebuf.str();
}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_Q_DG0<dim, spacedim>::clone() const
{
  return std::make_unique<FE_Q_DG0<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
void
FE_Q_DG0<dim, spacedim>::convert_generalized_support_point_values_to_dof_values(
  const std::vector<Vector<double>> &support_point_values,
  std::vector<double>               &nodal_dofs) const
{
  Assert(support_point_values.size() == this->unit_support_points.size(),
         ExcDimensionMismatch(support_point_values.size(),
                              this->unit_support_points.size()));
  Assert(nodal_dofs.size() == this->n_dofs_per_cell(),
         ExcDimensionMismatch(nodal_dofs.size(), this->n_dofs_per_cell()));
  Assert(support_point_values[0].size() == this->n_components(),
         ExcDimensionMismatch(support_point_values[0].size(),
                              this->n_components()));

  for (unsigned int i = 0; i < this->n_dofs_per_cell() - 1; ++i)
    {
      const std::pair<unsigned int, unsigned int> index =
        this->system_to_component_index(i);
      nodal_dofs[i] = support_point_values[i](index.first);
    }

  // We don't need the discontinuous function for local interpolation
  nodal_dofs.back() = 0.;
}



template <int dim, int spacedim>
void
FE_Q_DG0<dim, spacedim>::get_interpolation_matrix(
  const FiniteElement<dim, spacedim> &x_source_fe,
  FullMatrix<double>                 &interpolation_matrix) const
{
  // this is only implemented, if the source FE is also a Q_DG0 element
  using FEQDG0 = FE_Q_DG0<dim, spacedim>;

  AssertThrow(
    (x_source_fe.get_name().find("FE_Q_DG0<") == 0) ||
      (dynamic_cast<const FEQDG0 *>(&x_source_fe) != nullptr),
    (typename FiniteElement<dim, spacedim>::ExcInterpolationNotImplemented()));

  Assert(interpolation_matrix.m() == this->n_dofs_per_cell(),
         ExcDimensionMismatch(interpolation_matrix.m(),
                              this->n_dofs_per_cell()));
  Assert(interpolation_matrix.n() == x_source_fe.n_dofs_per_cell(),
         ExcDimensionMismatch(interpolation_matrix.m(),
                              x_source_fe.n_dofs_per_cell()));

  this->FE_Q_Base<dim, spacedim>::get_interpolation_matrix(
    x_source_fe, interpolation_matrix);
}



template <int dim, int spacedim>
std::vector<bool>
FE_Q_DG0<dim, spacedim>::get_riaf_vector(const unsigned int deg)
{
  std::vector<bool> riaf(Utilities::fixed_power<dim>(deg + 1) + 1, false);
  riaf.back() = true;
  return riaf;
}



template <int dim, int spacedim>
std::vector<unsigned int>
FE_Q_DG0<dim, spacedim>::get_dpo_vector(const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim + 1, 1U);
  for (unsigned int i = 1; i < dpo.size(); ++i)
    dpo[i] = dpo[i - 1] * (deg - 1);

  dpo[dim]++; // we need an additional DG0-node for a dim-dimensional object
  return dpo;
}



template <int dim, int spacedim>
bool
FE_Q_DG0<dim, spacedim>::has_support_on_face(
  const unsigned int shape_index,
  const unsigned int face_index) const
{
  // discontinuous function has support on all faces
  if (shape_index == this->n_dofs_per_cell() - 1)
    return true;
  else
    return FE_Q_Base<dim, spacedim>::has_support_on_face(shape_index,
                                                         face_index);
}



template <int dim, int spacedim>
std::pair<Table<2, bool>, std::vector<unsigned int>>
FE_Q_DG0<dim, spacedim>::get_constant_modes() const
{
  Table<2, bool> constant_modes(2, this->n_dofs_per_cell());

  // 1 represented by FE_Q part
  for (unsigned int i = 0; i < this->n_dofs_per_cell() - 1; ++i)
    constant_modes(0, i) = true;

  // 1 represented by DG0 part
  constant_modes(1, this->n_dofs_per_cell() - 1) = true;

  return std::pair<Table<2, bool>, std::vector<unsigned int>>(
    constant_modes, std::vector<unsigned int>(2, 0));
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_Q_DG0<dim, spacedim>::compare_for_domination(
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
  if (const FE_Q_DG0<dim, spacedim> *fe_dg0_other =
        dynamic_cast<const FE_Q_DG0<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_dg0_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_dg0_other->degree)
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


// explicit instantiations
#include "fe/fe_q_dg0.inst"

DEAL_II_NAMESPACE_CLOSE

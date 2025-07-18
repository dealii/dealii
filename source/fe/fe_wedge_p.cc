// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/config.h>

#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/qprojector.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_wedge_p.h>

DEAL_II_NAMESPACE_OPEN

namespace
{
  /**
   * Helper function to set up the dpo vector of FE_WedgeP for a given @p degree.
   */
  internal::GenericDoFsPerObject
  get_dpo_vector_fe_wedge_p(const unsigned int degree)
  {
    internal::GenericDoFsPerObject dpo;

    if (degree == 1)
      {
        dpo.dofs_per_object_exclusive  = {{1}, {0}, {0, 0, 0, 0, 0}, {0}};
        dpo.dofs_per_object_inclusive  = {{1}, {2}, {3, 3, 4, 4, 4}, {6}};
        dpo.object_index               = {{}, {6}, {6}, {6}};
        dpo.first_object_index_on_face = {{}, {3, 3, 4, 4, 4}, {3, 3, 4, 4, 4}};
      }
    else if (degree == 2)
      {
        dpo.dofs_per_object_exclusive  = {{1}, {1}, {0, 0, 1, 1, 1}, {0}};
        dpo.dofs_per_object_inclusive  = {{1}, {3}, {6, 6, 9, 9, 9}, {18}};
        dpo.object_index               = {{}, {6}, {15, 15, 15, 16, 17}, {18}};
        dpo.first_object_index_on_face = {{}, {3, 3, 4, 4, 4}, {6, 6, 8, 8, 8}};
      }
    else
      {
        DEAL_II_NOT_IMPLEMENTED();
      }

    return dpo;
  }

  /**
   * Helper function to set up the dpo vector of FE_WedgeDGP for a given @p degree.
   */
  internal::GenericDoFsPerObject
  get_dpo_vector_fe_wedge_dgp(const unsigned int degree)
  {
    unsigned int n_dofs = 0;

    if (degree == 1)
      n_dofs = 6;
    else if (degree == 2)
      n_dofs = 18;
    else
      DEAL_II_NOT_IMPLEMENTED();

    return internal::expand(3, {{0, 0, 0, n_dofs}}, ReferenceCells::Wedge);
  }
} // namespace



template <int dim, int spacedim>
FE_WedgePoly<dim, spacedim>::FE_WedgePoly(
  const unsigned int                                degree,
  const internal::GenericDoFsPerObject             &dpos,
  const typename FiniteElementData<dim>::Conformity conformity)
  : dealii::FE_Poly<dim, spacedim>(
      ScalarLagrangePolynomialWedge<dim>(degree),
      FiniteElementData<dim>(dpos,
                             ReferenceCells::Wedge,
                             1,
                             degree,
                             conformity),
      std::vector<bool>(
        FiniteElementData<dim>(dpos, ReferenceCells::Wedge, 1, degree)
          .dofs_per_cell,
        true),
      std::vector<ComponentMask>(
        FiniteElementData<dim>(dpos, ReferenceCells::Wedge, 1, degree)
          .dofs_per_cell,
        ComponentMask(std::vector<bool>(1, true))))
{
  AssertDimension(dim, 3);

  Assert(1 <= degree && degree <= 2, ExcNotImplemented());

  FE_SimplexP<2> fe_triangle(degree);
  FE_Q<1>        fe_line(degree);
  FE_Q<2>        fe_quad(degree);

  for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
    {
      const auto pair = this->degree == 1 ? internal::wedge_table_1[i] :
                                            internal::wedge_table_2[i];

      this->unit_support_points.emplace_back(
        fe_triangle.get_unit_support_points()[pair[0]][0],
        fe_triangle.get_unit_support_points()[pair[0]][1],
        fe_line.get_unit_support_points()[pair[1]][0]);
    }

  this->unit_face_support_points.resize(this->reference_cell().n_faces());

  for (const auto f : this->reference_cell().face_indices())
    if (this->reference_cell().face_reference_cell(f) ==
        ReferenceCells::Triangle)
      for (const auto &p : fe_triangle.get_unit_support_points())
        this->unit_face_support_points[f].emplace_back(p[0], p[1]);
    else if (this->reference_cell().face_reference_cell(f) ==
             ReferenceCells::Quadrilateral)
      for (const auto &p : fe_quad.get_unit_support_points())
        this->unit_face_support_points[f].emplace_back(p[0], p[1]);
    else
      DEAL_II_ASSERT_UNREACHABLE();
}



template <int dim, int spacedim>
void
FE_WedgePoly<dim, spacedim>::
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double>               &nodal_values) const
{
  AssertDimension(support_point_values.size(),
                  this->get_unit_support_points().size());
  AssertDimension(support_point_values.size(), nodal_values.size());
  AssertDimension(this->dofs_per_cell, nodal_values.size());

  for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
    {
      AssertDimension(support_point_values[i].size(), 1);

      nodal_values[i] = support_point_values[i](0);
    }
}



template <int dim, int spacedim>
FE_WedgeP<dim, spacedim>::FE_WedgeP(const unsigned int degree)
  : FE_WedgePoly<dim, spacedim>(degree,
                                get_dpo_vector_fe_wedge_p(degree),
                                FiniteElementData<dim>::H1)
{}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_WedgeP<dim, spacedim>::clone() const
{
  return std::make_unique<FE_WedgeP<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
std::string
FE_WedgeP<dim, spacedim>::get_name() const
{
  std::ostringstream namebuf;
  namebuf << "FE_WedgeP<" << Utilities::dim_string(dim, spacedim) << ">("
          << this->degree << ")";

  return namebuf.str();
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_WedgeP<dim, spacedim>::compare_for_domination(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  codim) const
{
  Assert(codim <= dim, ExcImpossibleInDim(dim));

  // vertex/line/face domination
  // (if fe_other is derived from FE_SimplexDGP)
  // ------------------------------------
  if (codim > 0)
    if (dynamic_cast<const FE_SimplexDGP<dim, spacedim> *>(&fe_other) !=
        nullptr)
      // there are no requirements between continuous and discontinuous
      // elements
      return FiniteElementDomination::no_requirements;


  // vertex/line/face domination
  // (if fe_other is not derived from FE_SimplexDGP)
  // & cell domination
  // ----------------------------------------
  if (const FE_WedgeP<dim, spacedim> *fe_wp_other =
        dynamic_cast<const FE_WedgeP<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_wp_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_wp_other->degree)
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
  else if (const FE_Q<dim, spacedim> *fe_q_other =
             dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_q_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_q_other->degree)
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
std::vector<std::pair<unsigned int, unsigned int>>
FE_WedgeP<dim, spacedim>::hp_vertex_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  (void)fe_other;

  Assert((dynamic_cast<const FE_SimplexP<dim, spacedim> *>(&fe_other)) ||
           (dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other)),
         ExcNotImplemented());

  return {{0, 0}};
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_WedgeP<dim, spacedim>::hp_line_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  (void)fe_other;

  Assert((dynamic_cast<const FE_SimplexP<dim, spacedim> *>(&fe_other)) ||
           (dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other)),
         ExcNotImplemented());

  std::vector<std::pair<unsigned int, unsigned int>> result;

  result.reserve(this->degree - 1);
  for (unsigned int i = 0; i < this->degree - 1; ++i)
    result.emplace_back(i, i);

  return result;
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_WedgeP<dim, spacedim>::hp_quad_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  face_no) const
{
  (void)fe_other;

  AssertIndexRange(face_no, 5);

  if (face_no < 2)
    {
      Assert((dynamic_cast<const FE_SimplexP<dim, spacedim> *>(&fe_other)),
             ExcNotImplemented());
    }
  else
    {
      Assert((dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other)),
             ExcNotImplemented());
    }

  std::vector<std::pair<unsigned int, unsigned int>> result;

  result.reserve(this->n_dofs_per_quad(face_no));
  for (unsigned int i = 0; i < this->n_dofs_per_quad(face_no); ++i)
    result.emplace_back(i, i);

  return result;
}



template <int dim, int spacedim>
FE_WedgeDGP<dim, spacedim>::FE_WedgeDGP(const unsigned int degree)
  : FE_WedgePoly<dim, spacedim>(degree,
                                get_dpo_vector_fe_wedge_dgp(degree),
                                FiniteElementData<dim>::L2)
{}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_WedgeDGP<dim, spacedim>::clone() const
{
  return std::make_unique<FE_WedgeDGP<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
std::string
FE_WedgeDGP<dim, spacedim>::get_name() const
{
  std::ostringstream namebuf;
  namebuf << "FE_WedgeDGP<" << Utilities::dim_string(dim, spacedim) << ">("
          << this->degree << ")";

  return namebuf.str();
}

// explicit instantiations
#include "fe/fe_wedge_p.inst"

DEAL_II_NAMESPACE_CLOSE

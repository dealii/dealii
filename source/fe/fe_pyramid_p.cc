// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
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
   * Helper function to set up the dpo vector of FE_PyramidP for a given @p degree.
   */
  unsigned int
  compute_n_dofs(const unsigned int dim, const unsigned int degree)
  {
    AssertDimension(dim, 3);
    return (degree + 1) * (degree + 2) * (2 * degree + 3) / 6;
  }

  /**
   * Helper function to set up the dpo vector of FE_PyramidP for a given @p degree.
   */
  template <int dim>
  std::pair<const internal::GenericDoFsPerObject, const std::vector<Point<dim>>>
  get_support_points_dpo_vector_fe_pyramid_p(const unsigned int degree)
  {
    AssertDimension(dim, 3);

    internal::GenericDoFsPerObject dpo;

    // set dpo for linear case and then add for higher orders
    dpo.dofs_per_object_exclusive = {{1}, {degree - 1}, {0, 0, 0, 0, 0}, {0}};
    dpo.dofs_per_object_inclusive = {{1},
                                     {degree + 1},
                                     {4, 3, 3, 3, 3},
                                     {compute_n_dofs(dim, degree)}};
    dpo.object_index = {{}, {5, 5, 5, 5, 5, 5, 5, 5}, {5, 5, 5, 5, 5}, {5}};
    dpo.first_object_index_on_face = {{}, {4, 3, 3, 3, 3}, {4, 3, 3, 3, 3}};

    std::vector<Point<dim>> support_points;
    std::vector<Point<dim>> support_points_unordered;

    std::vector<unsigned int> n_dofs_total_per_object(4, 0);
    std::vector<unsigned int> n_dofs_per_object(4, 0);

    support_points.resize(compute_n_dofs(dim, degree));
    const double z_equidistance = 1.0 / degree;

    for (unsigned int current_degree = degree; current_degree > 0;
         --current_degree)
      {
        FE_Q<2> fe_q(current_degree);

        // save all support points in an unordered fashion
        for (unsigned int support_point_index = 0;
             support_point_index < fe_q.get_unit_support_points().size();
             ++support_point_index)
          {
            Point<dim> support_point;
            for (unsigned int d = 0; d < dim - 1; ++d)
              support_point[d] =
                (current_degree * z_equidistance) *
                (2 * fe_q.get_unit_support_points()[support_point_index][d] -
                 1);
            support_point[dim - 1] = (degree - current_degree) * z_equidistance;
            support_points_unordered.emplace_back(support_point);
          }
        if (current_degree == degree)
          {
            // the first are the support points at the vertices
            n_dofs_total_per_object[0] = fe_q.reference_cell().n_vertices() + 1;

            // lines are on lines
            n_dofs_total_per_object[1] =
              fe_q.reference_cell().n_lines() * fe_q.n_dofs_per_line();
            // faces are on faces
            n_dofs_total_per_object[2] = fe_q.n_dofs_per_quad();
            // nothing per hex
            n_dofs_total_per_object[3] = 0;

            n_dofs_per_object[0] =
              fe_q.reference_cell().n_lines() * fe_q.n_dofs_per_line();
            n_dofs_per_object[2] = fe_q.n_dofs_per_quad();

            dpo.dofs_per_object_exclusive[2][0] = fe_q.n_dofs_per_quad();
            dpo.dofs_per_object_inclusive[2][0] = fe_q.n_dofs_per_cell();

            dpo.dofs_per_object_inclusive[2][1] += fe_q.n_dofs_per_line();
            dpo.dofs_per_object_inclusive[2][2] += fe_q.n_dofs_per_line();
            dpo.dofs_per_object_inclusive[2][3] += fe_q.n_dofs_per_line();
            dpo.dofs_per_object_inclusive[2][4] += fe_q.n_dofs_per_line();

            dpo.first_object_index_on_face[2][0] += n_dofs_per_object[0];
            dpo.first_object_index_on_face[2][1] += fe_q.n_dofs_per_line();
            dpo.first_object_index_on_face[2][2] += fe_q.n_dofs_per_line();
            dpo.first_object_index_on_face[2][3] += fe_q.n_dofs_per_line();
            dpo.first_object_index_on_face[2][4] += fe_q.n_dofs_per_line();


            dpo.object_index[1][1] += fe_q.n_dofs_per_line();
            dpo.object_index[1][2] += 2 * fe_q.n_dofs_per_line();
            dpo.object_index[1][3] += 3 * fe_q.n_dofs_per_line();
          }
        else
          {
            // on the vertex here is on the line in the pyramid
            n_dofs_total_per_object[1] += fe_q.reference_cell().n_vertices();

            dpo.first_object_index_on_face[2][1] += 2;
            dpo.first_object_index_on_face[2][2] += 2;
            dpo.first_object_index_on_face[2][3] += 2;
            dpo.first_object_index_on_face[2][4] += 2;

            // on the lines here are on faces in the pyramid
            n_dofs_total_per_object[2] +=
              fe_q.reference_cell().n_lines() * fe_q.n_dofs_per_line();

            dpo.dofs_per_object_exclusive[2][1] += fe_q.n_dofs_per_line();
            dpo.dofs_per_object_exclusive[2][2] += fe_q.n_dofs_per_line();
            dpo.dofs_per_object_exclusive[2][3] += fe_q.n_dofs_per_line();
            dpo.dofs_per_object_exclusive[2][4] += fe_q.n_dofs_per_line();
            dpo.dofs_per_object_inclusive[2][1] += fe_q.n_dofs_per_line() + 2;
            dpo.dofs_per_object_inclusive[2][2] += fe_q.n_dofs_per_line() + 2;
            dpo.dofs_per_object_inclusive[2][3] += fe_q.n_dofs_per_line() + 2;
            dpo.dofs_per_object_inclusive[2][4] += fe_q.n_dofs_per_line() + 2;


            // on the faces here is on the hex on the pyramid
            n_dofs_total_per_object[3] += fe_q.n_dofs_per_quad();
            dpo.dofs_per_object_exclusive[3][0] += fe_q.n_dofs_per_quad();


            n_dofs_per_object[1] += 1;
            n_dofs_per_object[3] += fe_q.n_dofs_per_line();
          }
      }

    unsigned int global_counter = 0;

    std::vector<unsigned int> start_lines(4);
    start_lines[0] = n_dofs_total_per_object[0] + n_dofs_per_object[0];
    start_lines[1] = start_lines[0] + n_dofs_per_object[1];
    start_lines[2] = start_lines[1] + n_dofs_per_object[1];
    start_lines[3] = start_lines[2] + n_dofs_per_object[1];

    dpo.object_index[1][4] = start_lines[0];
    dpo.object_index[1][5] = start_lines[1];
    dpo.object_index[1][6] = start_lines[2];
    dpo.object_index[1][7] = start_lines[3];

    std::vector<unsigned int> start_faces(4);
    start_faces[0] = n_dofs_total_per_object[0] + n_dofs_total_per_object[1] +
                     n_dofs_per_object[2];
    start_faces[1] = start_faces[0] + n_dofs_per_object[3];
    start_faces[2] = start_faces[1] + n_dofs_per_object[3];
    start_faces[3] = start_faces[2] + n_dofs_per_object[3];

    dpo.object_index[2][0] =
      n_dofs_total_per_object[0] + n_dofs_total_per_object[1];
    dpo.object_index[2][1] = start_faces[0];
    dpo.object_index[2][2] = start_faces[1];
    dpo.object_index[2][3] = start_faces[2];
    dpo.object_index[2][4] = start_faces[3];


    unsigned int start_hex = n_dofs_total_per_object[0] +
                             n_dofs_total_per_object[1] +
                             n_dofs_total_per_object[2];

    dpo.object_index[3][0] = start_hex;

    for (unsigned int current_degree = degree; current_degree > 0;
         --current_degree)
      {
        FE_Q<2> fe_q(current_degree);

        if (current_degree == degree)
          {
            // this gives all info on the vertices, the first 4 edges and the
            // first face

            // vertices
            for (unsigned int counter = 0;
                 counter < fe_q.reference_cell().n_vertices();
                 ++counter)
              {
                support_points[counter] =
                  support_points_unordered[global_counter++];
              }
            // lines
            for (unsigned int counter = 0;
                 counter <
                 fe_q.reference_cell().n_lines() * fe_q.n_dofs_per_line();
                 ++counter)
              {
                support_points[counter + n_dofs_total_per_object[0]] =
                  support_points_unordered[global_counter++];
              }
            // quad
            for (unsigned int counter = 0; counter < fe_q.n_dofs_per_quad();
                 ++counter)
              {
                support_points[counter + n_dofs_total_per_object[0] +
                               n_dofs_total_per_object[1]] =
                  support_points_unordered[global_counter++];
              }
          }
        else
          {
            // vertices are on lines
            for (unsigned int line = 0;
                 line < fe_q.reference_cell().n_vertices();
                 ++line)
              {
                support_points[start_lines[line]++] =
                  support_points_unordered[global_counter++];
              }
            // lines are on face
            for (unsigned int face = 0; face < fe_q.reference_cell().n_lines();
                 ++face)
              {
                for (unsigned int n_dof = 0; n_dof < fe_q.n_dofs_per_line();
                     ++n_dof)
                  support_points[start_faces[face]++] =
                    support_points_unordered[global_counter++];
              }
            // faces are on hex
            for (unsigned int hex = 0; hex < fe_q.n_dofs_per_quad(); ++hex)
              {
                support_points[start_hex++] =
                  support_points_unordered[global_counter++];
              }
          }
      }
    Point<dim> tip;
    for (unsigned int d = 0; d < dim; ++d)
      {
        if (d == 2)
          tip[d] = 1.0;
        else
          tip[d] = 0.0;
      }
    support_points[4] = tip;


    std::pair<internal::GenericDoFsPerObject, std::vector<Point<dim>>>
      dpos_support_points;
    dpos_support_points.first  = dpo;
    dpos_support_points.second = support_points;

    return dpos_support_points;
  }


  /**
   * Helper function to set up the dpo vector of FE_PyramidP for a given @p degree.
   */
  template <int dim>
  std::pair<const internal::GenericDoFsPerObject, const std::vector<Point<dim>>>
  get_support_points_dpo_vector_fe_pyramid_dgp(const unsigned int degree)
  {
    AssertDimension(dim, 3);
    const auto support_points =
      get_support_points_dpo_vector_fe_pyramid_p<dim>(degree).second;

    std::pair<internal::GenericDoFsPerObject, std::vector<Point<dim>>>
      dpos_support_points;
    dpos_support_points.second = support_points;
    dpos_support_points.first =
      internal::expand(3,
                       {{0, 0, 0, compute_n_dofs(dim, degree)}},
                       ReferenceCells::Pyramid);

    return dpos_support_points;
  }
} // namespace


template <int dim, int spacedim>
FE_PyramidPoly<dim, spacedim>::FE_PyramidPoly(
  const unsigned int degree,
  std::pair<const internal::GenericDoFsPerObject, const std::vector<Point<dim>>>
                                                    dpos_support_points,
  const bool                                        prolongation_is_additive,
  const typename FiniteElementData<dim>::Conformity conformity)
  : dealii::FE_Poly<dim, spacedim>(
      ScalarLagrangePolynomialPyramid<dim>(degree,
                                           compute_n_dofs(dim, degree),
                                           dpos_support_points.second),
      FiniteElementData<dim>(dpos_support_points.first,
                             ReferenceCells::Pyramid,
                             1,
                             degree,
                             conformity),
      std::vector<bool>(FiniteElementData<dim>(dpos_support_points.first,
                                               ReferenceCells::Pyramid,
                                               1,
                                               degree)
                          .dofs_per_cell,
                        prolongation_is_additive),
      std::vector<ComponentMask>(
        FiniteElementData<dim>(dpos_support_points.first,
                               ReferenceCells::Pyramid,
                               1,
                               degree)
          .dofs_per_cell,
        ComponentMask(std::vector<bool>(1, true))))
{
  AssertDimension(dim, 3);
    
  for (auto &support_point : dpos_support_points.second)
    this->unit_support_points.emplace_back(support_point);
  
  if (conformity == FiniteElementData<dim>::H1)
    {
      // face support points
      this->unit_face_support_points.resize(
        this->reference_cell().n_faces());

      for (const auto f : this->reference_cell().face_indices())
        {
          const auto face_reference_cell =
            this->reference_cell().face_reference_cell(f);

          if (face_reference_cell == ReferenceCells::Quadrilateral)
            {
              FE_Q<2> fe_face(degree);
              for (auto &face_support_point :
                    fe_face.get_unit_support_points())
                {
                  Point<dim - 1> p;
                  for (unsigned int d = 0; d < dim - 1; ++d)
                    p[d] = face_support_point[d];
                  this->unit_face_support_points[f].emplace_back(p);
                }
            }
          else if (face_reference_cell == ReferenceCells::Triangle)
            {
              FE_SimplexP<2> fe_face(degree);
              for (auto &face_support_point :
                    fe_face.get_unit_support_points())
                {
                  Point<dim - 1> p;
                  for (unsigned int d = 0; d < dim - 1; ++d)
                    p[d] = face_support_point[d];
                  this->unit_face_support_points[f].emplace_back(p);
                }
            }
        }
    }
}



template <int dim, int spacedim>
void
FE_PyramidPoly<dim, spacedim>::
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
FE_PyramidP<dim, spacedim>::FE_PyramidP(const unsigned int degree)
  : FE_PyramidPoly<dim, spacedim>(
      degree,
      get_support_points_dpo_vector_fe_pyramid_p<dim>(degree),
      false,
      FiniteElementData<dim>::H1)
{}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_PyramidP<dim, spacedim>::clone() const
{
  return std::make_unique<FE_PyramidP<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
std::string
FE_PyramidP<dim, spacedim>::get_name() const
{
  std::ostringstream namebuf;
  namebuf << "FE_PyramidP<" << Utilities::dim_string(dim, spacedim) << ">("
          << this->degree << ")";

  return namebuf.str();
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_PyramidP<dim, spacedim>::compare_for_domination(
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
  if (const FE_PyramidP<dim, spacedim> *fe_pp_other =
        dynamic_cast<const FE_PyramidP<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_pp_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_pp_other->degree)
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

  DEAL_II_NOT_IMPLEMENTED();
  return FiniteElementDomination::neither_element_dominates;
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_PyramidP<dim, spacedim>::hp_vertex_dof_identities(
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
FE_PyramidP<dim, spacedim>::hp_line_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  (void)fe_other;

  Assert((dynamic_cast<const FE_SimplexP<dim, spacedim> *>(&fe_other)) ||
           (dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other)),
         ExcNotImplemented());

  std::vector<std::pair<unsigned int, unsigned int>> result;

  for (unsigned int i = 0; i < this->degree - 1; ++i)
    result.emplace_back(i, i);

  return result;
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_PyramidP<dim, spacedim>::hp_quad_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  face_no) const
{
  (void)fe_other;


  AssertIndexRange(face_no, 5);

  if (face_no == 0)
    {
      Assert((dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other)),
             ExcNotImplemented());
    }
  else
    {
      Assert((dynamic_cast<const FE_SimplexP<dim, spacedim> *>(&fe_other)),
             ExcNotImplemented());
    }

  std::vector<std::pair<unsigned int, unsigned int>> result;
  for (unsigned int i = 0; i < this->n_dofs_per_quad(face_no); ++i)
    result.emplace_back(i, i);

  return result;
}



template <int dim, int spacedim>
FE_PyramidDGP<dim, spacedim>::FE_PyramidDGP(const unsigned int degree)
  : FE_PyramidPoly<dim, spacedim>(
      degree,
      get_support_points_dpo_vector_fe_pyramid_dgp<dim>(degree),
      true,
      FiniteElementData<dim>::L2)
{}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_PyramidDGP<dim, spacedim>::clone() const
{
  return std::make_unique<FE_PyramidDGP<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
std::string
FE_PyramidDGP<dim, spacedim>::get_name() const
{
  std::ostringstream namebuf;
  namebuf << "FE_PyramidDGP<" << Utilities::dim_string(dim, spacedim) << ">("
          << this->degree << ")";

  return namebuf.str();
}

// explicit instantiations
#include "fe_pyramid_p.inst"

DEAL_II_NAMESPACE_CLOSE

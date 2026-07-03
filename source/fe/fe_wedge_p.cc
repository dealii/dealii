// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2021 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/base/config.h>

#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>

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

    return internal::expand<3>({{0, 0, 0, n_dofs}}, ReferenceCells::Wedge);
  }

  /**
   * Helper function that returns a vector that contains the unit cell support
   * points for FE_WedgePoly.
   */
  template <int dim>
  std::vector<Point<dim>>
  unit_support_points_fe_wedge_p(const unsigned int degree)
  {
    Assert(degree > 0, ExcNotImplemented());

    if constexpr (dim == 3)
      {
        std::vector<Point<dim>> unit_points;

        const auto reference_cell = ReferenceCells::Wedge;

        // construct the support points as the tensor product of a triangle and
        // a line
        // FE_SimplexP returns equidistant points while FE_Q gives Gauss-Lobatto
        // points for degree > 2, so the edges of the faces are not compatible
        // at degree > 2 if the "normal" FE_Q is used here, thus, use an
        // equidistant FE_Q instead
        // when FE_SimplexP has electrostatic support points switch to FE_Q
        // until then assert that we are not at a degree larger than 2
        Assert(degree < 3, ExcInternalError());

        const FE_SimplexP<1> fe_line(degree);
        const FE_SimplexP<2> fe_triangle(degree);
        const FE_Q<2>        fe_quad(QIterated<1>(QTrapezoid<1>(), degree));

        // start with the vertices
        for (const unsigned int v : reference_cell.vertex_indices())
          unit_points.push_back(reference_cell.vertex(v));

        // lines
        for (const unsigned int l : reference_cell.line_indices())
          {
            const Point<dim> v0 =
              reference_cell.vertex(reference_cell.line_to_cell_vertices(l, 0));
            const Point<dim> v1 =
              reference_cell.vertex(reference_cell.line_to_cell_vertices(l, 1));

            for (unsigned int i = 0; i < degree - 1; ++i)
              {
                // shift the point on the line such that the support points on
                // the edges are compatible
                const double distance = fe_line.unit_support_point(
                  fe_line.get_first_line_index() + i)[0];
                unit_points.push_back(v0 + distance * (v1 - v0));
              }
          }

        // faces
        for (const unsigned int f : reference_cell.face_indices())
          {
            const auto face_reference_cell =
              reference_cell.face_reference_cell(f);

            const bool is_triangular_face =
              reference_cell.face_reference_cell(f).is_simplex();

            const unsigned int n_dofs_per_quad =
              is_triangular_face ? fe_triangle.n_dofs_per_quad() :
                                   fe_quad.n_dofs_per_quad();

            const std::vector<Point<2>> face_support_points =
              is_triangular_face ? fe_triangle.get_unit_support_points() :
                                   fe_quad.get_unit_support_points();

            const unsigned int first_quad_index =
              is_triangular_face ? fe_triangle.get_first_quad_index() :
                                   fe_quad.get_first_quad_index();

            // go over all DoFs on the face
            for (unsigned int i = 0; i < n_dofs_per_quad; ++i)
              {
                Point<dim> p;
                // linear interpolate the point form the vertices of the face
                // looking at the triangle it is the same as using barycentric
                // coordinates as the linear shape functions are: 1-x-y, x, y
                for (unsigned int v = 0; v < face_reference_cell.n_vertices();
                     ++v)
                  {
                    const auto vertex = reference_cell.vertex(
                      reference_cell.face_to_cell_vertices(
                        f, v, numbers::default_geometric_orientation));

                    p += face_reference_cell.d_linear_shape_function(
                           face_support_points[first_quad_index + i], v) *
                         vertex;
                  }
                unit_points.push_back(p);
              }
          }

        // interior, this is just the tensor product of the interior nodes of
        // the triangle with the line
        for (unsigned int i = 0; i < degree - 1; ++i)
          for (unsigned int j = 0; j < fe_triangle.n_dofs_per_quad(); ++j)
            {
              const double z = fe_line.unit_support_point(
                fe_line.get_first_line_index() + i)[0];

              const Point<2> x_y = fe_triangle.unit_support_point(
                fe_triangle.get_first_quad_index() + j);

              unit_points.push_back(Point<dim>(x_y[0], x_y[1], z));
            }

        return unit_points;
      }
    else
      DEAL_II_ASSERT_UNREACHABLE();
    return {};
  }
} // namespace



template <int dim, int spacedim>
FE_WedgePoly<dim, spacedim>::FE_WedgePoly(
  const unsigned int                                degree,
  const internal::GenericDoFsPerObject             &dpos,
  const bool                                        prolongation_is_additive,
  const typename FiniteElementData<dim>::Conformity conformity)
  : dealii::FE_Poly<dim, spacedim>(
      ScalarLagrangePolynomialWedge<dim>(degree,
                                         unit_support_points_fe_wedge_p<dim>(
                                           degree)),
      FiniteElementData<dim>(dpos,
                             reinterpret_cast<const ReferenceCell<dim> &>(
                               ReferenceCells::Wedge),
                             1,
                             degree,
                             conformity),
      std::vector<bool>(
        FiniteElementData<dim>(dpos,
                               reinterpret_cast<const ReferenceCell<dim> &>(
                                 ReferenceCells::Wedge),
                               1,
                               degree)
          .dofs_per_cell,
        prolongation_is_additive),
      std::vector<ComponentMask>(
        FiniteElementData<dim>(dpos,
                               reinterpret_cast<const ReferenceCell<dim> &>(
                                 ReferenceCells::Wedge),
                               1,
                               degree)
          .dofs_per_cell,
        ComponentMask(std::vector<bool>(1, true))))
{
  AssertDimension(dim, 3);

  Assert(1 <= degree && degree <= 2, ExcNotImplemented());

  const FE_SimplexP<2> fe_triangle(degree);
  const FE_Q<2>        fe_quad(degree);

  this->unit_support_points = unit_support_points_fe_wedge_p<dim>(degree);

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
                                false,
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
  else if (const FE_PyramidP<dim, spacedim> *fe_p_other =
             dynamic_cast<const FE_PyramidP<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_p_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_p_other->degree)
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
           (dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other)) ||
           (dynamic_cast<const FE_PyramidP<dim, spacedim> *>(&fe_other)) ||
           (dynamic_cast<const FE_WedgeP<dim, spacedim> *>(&fe_other)),
         ExcNotImplemented());

  return {{0, 0}};
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_WedgeP<dim, spacedim>::hp_line_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  Assert((dynamic_cast<const FE_SimplexP<dim, spacedim> *>(&fe_other)) ||
           (dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other)) ||
           (dynamic_cast<const FE_PyramidP<dim, spacedim> *>(&fe_other)) ||
           (dynamic_cast<const FE_WedgeP<dim, spacedim> *>(&fe_other)),
         ExcNotImplemented());

  std::vector<std::pair<unsigned int, unsigned int>> identities;
  // check if the support points are the same location on the line
  // to avoid rescaling for pyramids use the support points on the faces
  const auto face_support_points = this->get_unit_face_support_points(0);
  const auto face_support_points_other =
    fe_other.get_unit_face_support_points(0);

  // now just compare the DoFs on the line going from [0,0] to [1,0]
  // for a triangular face that is the first line
  // for a quad face that is the third line
  // adjust the offsets accordingly
  // face number 0 of the wedge is a triangle
  const unsigned int offset =
    this->reference_cell().face_reference_cell(0).n_vertices();

  const unsigned int offset_other =
    fe_other.reference_cell().face_reference_cell(0).is_hyper_cube() ?
      fe_other.reference_cell().face_reference_cell(0).n_vertices() +
        2 * fe_other.n_dofs_per_line() :
      fe_other.reference_cell().face_reference_cell(0).n_vertices();

  // now get the identities
  for (unsigned int i = 0; i < this->n_dofs_per_line(); ++i)
    for (unsigned int j = 0; j < fe_other.n_dofs_per_line(); ++j)
      if (face_support_points[i + offset].distance(
            face_support_points_other[j + offset_other]) < 1e-14)
        identities.emplace_back(i, j);

  return identities;
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_WedgeP<dim, spacedim>::hp_quad_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  face_no) const
{
  AssertIndexRange(face_no, 5);

  std::vector<std::pair<unsigned int, unsigned int>> result;
  unsigned int                                       face_no_neighbor;

  if (face_no < 2)
    {
      // triangular face
      Assert((dynamic_cast<const FE_SimplexP<dim, spacedim> *>(&fe_other)) ||
               (dynamic_cast<const FE_PyramidP<dim, spacedim> *>(&fe_other)) ||
               (dynamic_cast<const FE_WedgeP<dim, spacedim> *>(&fe_other)),
             ExcNotImplemented());
      if ((dynamic_cast<const FE_PyramidP<dim, spacedim> *>(&fe_other)))
        face_no_neighbor = 1;
      else
        face_no_neighbor = 0;
    }
  else
    {
      // quad face
      Assert((dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other)) ||
               (dynamic_cast<const FE_PyramidP<dim, spacedim> *>(&fe_other)) ||
               (dynamic_cast<const FE_WedgeP<dim, spacedim> *>(&fe_other)),
             ExcNotImplemented());
      if ((dynamic_cast<const FE_WedgeP<dim, spacedim> *>(&fe_other)))
        face_no_neighbor = 2;
      else
        face_no_neighbor = 0;
    }

  // compare the face support points
  const auto face_support_points = this->get_unit_face_support_points(face_no);
  const auto face_support_points_other =
    fe_other.get_unit_face_support_points(face_no_neighbor);

  // get the offsets to skip vertices and lines
  const auto face_reference_cell =
    this->reference_cell().face_reference_cell(face_no);
  Assert(face_reference_cell ==
           fe_other.reference_cell().face_reference_cell(face_no_neighbor),
         ExcInternalError());

  const auto offset = face_reference_cell.n_vertices() +
                      face_reference_cell.n_lines() * this->n_dofs_per_line();

  const auto offset_other =
    face_reference_cell.n_vertices() +
    face_reference_cell.n_lines() * fe_other.n_dofs_per_line();

  // now compare the points
  for (unsigned int i = 0; i < this->n_dofs_per_quad(face_no); ++i)
    for (unsigned int j = 0; j < fe_other.n_dofs_per_quad(face_no_neighbor);
         ++j)
      if (face_support_points[i + offset].distance(
            face_support_points_other[j + offset_other]) < 1e-14)
        result.emplace_back(i, j);
  return result;
}



template <int dim, int spacedim>
FE_WedgeDGP<dim, spacedim>::FE_WedgeDGP(const unsigned int degree)
  : FE_WedgePoly<dim, spacedim>(degree,
                                get_dpo_vector_fe_wedge_dgp(degree),
                                true,
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

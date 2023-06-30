// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2023 by the deal.II authors
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

#include <deal.II/base/config.h>

#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/qprojector.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_tools.h>

DEAL_II_NAMESPACE_OPEN

namespace
{
  // TODO: replace this and the equivalent inlined function in
  // FE_SimplexP_Bubbles with QProjector
  template <int dim>
  Point<dim>
  face_midpoint(const unsigned int face_no)
  {
    const auto reference_cell = ReferenceCells::get_simplex<dim>();
    const auto face_reference_cell =
      reference_cell.face_reference_cell(face_no);

    Point<dim> midpoint;
    for (const auto face_vertex_no : face_reference_cell.vertex_indices())
      {
        const auto vertex_no = reference_cell.face_to_cell_vertices(
          face_no,
          face_vertex_no,
          ReferenceCell::default_combined_face_orientation());

        midpoint += reference_cell.template vertex<dim>(vertex_no);
      }

    midpoint /= reference_cell.face_reference_cell(0).n_vertices();
    return midpoint;
  }

  /**
   * Helper function to set up the dpo vector of FE_SimplexP for a given @p dim and
   * @p degree.
   */
  std::vector<unsigned int>
  get_dpo_vector_fe_p(const unsigned int dim, const unsigned int degree)
  {
    switch (dim)
      {
        case 1:
          switch (degree)
            {
              case 1:
                return {1, 0};
              case 2:
                return {1, 1};
              case 3:
                return {1, 2};
              default:
                Assert(false, ExcNotImplemented());
            }
        case 2:
          switch (degree)
            {
              case 1:
                return {1, 0, 0};
              case 2:
                return {1, 1, 0};
              case 3:
                return {1, 2, 1};
              default:
                Assert(false, ExcNotImplemented());
            }
        case 3:
          switch (degree)
            {
              case 1:
                return {1, 0, 0, 0};
              case 2:
                return {1, 1, 0, 0};
              case 3:
                return {1, 2, 1, 0};
              default:
                Assert(false, ExcNotImplemented());
            }
      }

    Assert(false, ExcNotImplemented());
    return {};
  }



  /**
   * Set up a vector that contains the unit (reference) cell support points
   * for FE_SimplexPoly and sufficiently similar elements.
   */
  template <int dim>
  std::vector<Point<dim>>
  unit_support_points_fe_p(const unsigned int degree)
  {
    Assert(dim != 0, ExcInternalError());
    Assert(degree <= 3, ExcNotImplemented());
    std::vector<Point<dim>> unit_points;
    const auto              reference_cell = ReferenceCells::get_simplex<dim>();

    // Piecewise constants are a special case: use a support point at the
    // centroid and only the centroid
    if (degree == 0)
      {
        unit_points.emplace_back(reference_cell.template barycenter<dim>());
        return unit_points;
      }

    // otherwise write everything as linear combinations of vertices
    const auto dpo = get_dpo_vector_fe_p(dim, degree);
    Assert(dpo.size() == dim + 1, ExcInternalError());
    Assert(dpo[0] == 1, ExcNotImplemented());
    // vertices:
    for (unsigned int d : reference_cell.vertex_indices())
      unit_points.push_back(reference_cell.template vertex<dim>(d));
    // lines:
    for (unsigned int l : reference_cell.line_indices())
      {
        const Point<dim> p0 =
          unit_points[reference_cell.line_to_cell_vertices(l, 0)];
        const Point<dim> p1 =
          unit_points[reference_cell.line_to_cell_vertices(l, 1)];
        for (unsigned int p = 0; p < dpo[1]; ++p)
          unit_points.push_back((double(dpo[1] - p) / (dpo[1] + 1)) * p0 +
                                (double(p + 1) / (dpo[1] + 1)) * p1);
      }
    // quads:
    if (dim > 1 && dpo[2] > 0)
      {
        Assert(dpo[2] == 1, ExcNotImplemented());
        if (dim == 2)
          unit_points.push_back(reference_cell.template barycenter<dim>());
        if (dim == 3)
          for (unsigned int f : reference_cell.face_indices())
            unit_points.push_back(face_midpoint<dim>(f));
      }

    return unit_points;
  }

  template <>
  std::vector<Point<0>>
  unit_support_points_fe_p(const unsigned int /*degree*/)
  {
    std::vector<Point<0>> unit_points;
    unit_points.emplace_back();
    return unit_points;
  }

  /**
   * Set up a vector that contains the unit (reference) cell's faces support
   * points for FE_SimplexP and sufficiently similar elements.
   */
  template <int dim>
  std::vector<std::vector<Point<dim - 1>>>
  unit_face_support_points_fe_p(
    const unsigned int                          degree,
    typename FiniteElementData<dim>::Conformity conformity)
  {
    // Discontinuous elements don't have face support points
    if (conformity == FiniteElementData<dim>::Conformity::L2)
      return {};

    // this concept doesn't exist in 1d so just return an empty vector
    if (dim == 1)
      return {};

    std::vector<std::vector<Point<dim - 1>>> unit_face_points;

    // all faces have the same support points
    for (auto face_n : ReferenceCells::get_simplex<dim>().face_indices())
      {
        (void)face_n;
        unit_face_points.emplace_back(
          unit_support_points_fe_p<dim - 1>(degree));
      }

    return unit_face_points;
  }

  /**
   * Specify the constraints which the dofs on the two sides of a cell
   * interface underlie if the line connects two cells of which one is refined
   * once.
   */
  template <int dim>
  FullMatrix<double>
  constraints_fe_p(const unsigned int /*degree*/)
  {
    // no constraints in 1d
    // constraints in 3d not implemented yet
    return FullMatrix<double>();
  }

  template <>
  FullMatrix<double>
  constraints_fe_p<2>(const unsigned int degree)
  {
    const unsigned int dim = 2;

    // the following implements the 2d case
    // (the 3d case is not implemented yet)
    //
    // consult FE_Q_Base::Implementation::initialize_constraints()
    // for more information

    std::vector<Point<dim - 1>> constraint_points;
    // midpoint
    constraint_points.emplace_back(0.5);
    if (degree == 2)
      {
        // midpoint on subface 0
        constraint_points.emplace_back(0.25);
        // midpoint on subface 1
        constraint_points.emplace_back(0.75);
      }

    // Now construct relation between destination (child) and source (mother)
    // dofs.

    const unsigned int n_dofs_constrained = constraint_points.size();
    unsigned int       n_dofs_per_face    = degree + 1;
    FullMatrix<double> interface_constraints(n_dofs_constrained,
                                             n_dofs_per_face);

    const auto poly = BarycentricPolynomials<dim - 1>::get_fe_p_basis(degree);

    for (unsigned int i = 0; i < n_dofs_constrained; ++i)
      for (unsigned int j = 0; j < n_dofs_per_face; ++j)
        {
          interface_constraints(i, j) =
            poly.compute_value(j, constraint_points[i]);

          // if the value is small up to round-off, then simply set it to zero
          // to avoid unwanted fill-in of the constraint matrices (which would
          // then increase the number of other DoFs a constrained DoF would
          // couple to)
          if (std::fabs(interface_constraints(i, j)) < 1e-13)
            interface_constraints(i, j) = 0;
        }
    return interface_constraints;
  }



  /**
   * Helper function to set up the dpo vector of FE_SimplexDGP for a given
   * @p dim and @p degree.
   */
  std::vector<unsigned int>
  get_dpo_vector_fe_dgp(const unsigned int dim, const unsigned int degree)
  {
    // First treat the case of piecewise constant elements:
    if (degree == 0)
      {
        std::vector<unsigned int> dpo(dim + 1, 0U);
        dpo[dim] = 1;
        return dpo;
      }
    else
      {
        // This element has the same degrees of freedom as the continuous one,
        // but they are all counted for the interior of the cell because
        // it is continuous. Rather than hard-code how many DoFs the element
        // has, we just get the numbers from the continuous case and add them
        // up
        const auto continuous_dpo = get_dpo_vector_fe_p(dim, degree);

        switch (dim)
          {
            case 1:
              return {0U,
                      ReferenceCells::Line.n_vertices() * continuous_dpo[0] +
                        continuous_dpo[dim]};

            case 2:
              return {0U,
                      0U,
                      ReferenceCells::Triangle.n_vertices() *
                          continuous_dpo[0] +
                        ReferenceCells::Triangle.n_lines() * continuous_dpo[1] +
                        continuous_dpo[dim]};

            case 3:
              return {
                0U,
                0U,
                0U,
                ReferenceCells::Tetrahedron.n_vertices() * continuous_dpo[0] +
                  ReferenceCells::Tetrahedron.n_lines() * continuous_dpo[1] +
                  ReferenceCells::Tetrahedron.n_faces() * continuous_dpo[2] +
                  continuous_dpo[dim]};
          }

        Assert(false, ExcNotImplemented());
        return {};
      }
  }
} // namespace



template <int dim, int spacedim>
FE_SimplexPoly<dim, spacedim>::FE_SimplexPoly(
  const BarycentricPolynomials<dim>              polynomials,
  const FiniteElementData<dim> &                 fe_data,
  const std::vector<Point<dim>> &                unit_support_points,
  const std::vector<std::vector<Point<dim - 1>>> unit_face_support_points,
  const FullMatrix<double> &                     interface_constraints)
  : dealii::FE_Poly<dim, spacedim>(
      polynomials,
      fe_data,
      std::vector<bool>(fe_data.dofs_per_cell),
      std::vector<ComponentMask>(fe_data.dofs_per_cell,
                                 std::vector<bool>(1, true)))
{
  this->unit_support_points      = unit_support_points;
  this->unit_face_support_points = unit_face_support_points;
  this->interface_constraints    = interface_constraints;
}



template <int dim, int spacedim>
std::pair<Table<2, bool>, std::vector<unsigned int>>
FE_SimplexPoly<dim, spacedim>::get_constant_modes() const
{
  Table<2, bool> constant_modes(1, this->n_dofs_per_cell());
  constant_modes.fill(true);
  return std::pair<Table<2, bool>, std::vector<unsigned int>>(
    constant_modes, std::vector<unsigned int>(1, 0));
}



template <int dim, int spacedim>
const FullMatrix<double> &
FE_SimplexPoly<dim, spacedim>::get_prolongation_matrix(
  const unsigned int         child,
  const RefinementCase<dim> &refinement_case) const
{
  Assert(refinement_case == RefinementCase<dim>::isotropic_refinement,
         ExcNotImplemented());
  AssertDimension(dim, spacedim);

  // initialization upon first request
  if (this->prolongation[refinement_case - 1][child].n() == 0)
    {
      std::lock_guard<std::mutex> lock(this->mutex);

      // if matrix got updated while waiting for the lock
      if (this->prolongation[refinement_case - 1][child].n() ==
          this->n_dofs_per_cell())
        return this->prolongation[refinement_case - 1][child];

      // now do the work. need to get a non-const version of data in order to
      // be able to modify them inside a const function
      auto &this_nonconst = const_cast<FE_SimplexPoly<dim, spacedim> &>(*this);

      std::vector<std::vector<FullMatrix<double>>> isotropic_matrices(
        RefinementCase<dim>::isotropic_refinement);
      isotropic_matrices.back().resize(
        GeometryInfo<dim>::n_children(RefinementCase<dim>(refinement_case)),
        FullMatrix<double>(this->n_dofs_per_cell(), this->n_dofs_per_cell()));

      FETools::compute_embedding_matrices(*this, isotropic_matrices, true);

      this_nonconst.prolongation[refinement_case - 1].swap(
        isotropic_matrices.back());
    }

  // finally return the matrix
  return this->prolongation[refinement_case - 1][child];
}



template <int dim, int spacedim>
const FullMatrix<double> &
FE_SimplexPoly<dim, spacedim>::get_restriction_matrix(
  const unsigned int         child,
  const RefinementCase<dim> &refinement_case) const
{
  Assert(refinement_case == RefinementCase<dim>::isotropic_refinement,
         ExcNotImplemented());
  AssertDimension(dim, spacedim);

  // initialization upon first request
  if (this->restriction[refinement_case - 1][child].n() == 0)
    {
      std::lock_guard<std::mutex> lock(this->mutex);

      // if matrix got updated while waiting for the lock
      if (this->restriction[refinement_case - 1][child].n() ==
          this->n_dofs_per_cell())
        return this->restriction[refinement_case - 1][child];

      // now do the work. need to get a non-const version of data in order to
      // be able to modify them inside a const function
      auto &this_nonconst = const_cast<FE_SimplexPoly<dim, spacedim> &>(*this);

      std::vector<std::vector<FullMatrix<double>>> isotropic_matrices(
        RefinementCase<dim>::isotropic_refinement);
      isotropic_matrices.back().resize(
        GeometryInfo<dim>::n_children(RefinementCase<dim>(refinement_case)),
        FullMatrix<double>(this->n_dofs_per_cell(), this->n_dofs_per_cell()));

      FETools::compute_projection_matrices(*this, isotropic_matrices, true);

      this_nonconst.restriction[refinement_case - 1].swap(
        isotropic_matrices.back());
    }

  // finally return the matrix
  return this->restriction[refinement_case - 1][child];
}



template <int dim, int spacedim>
void
FE_SimplexPoly<dim, spacedim>::get_face_interpolation_matrix(
  const FiniteElement<dim, spacedim> &source_fe,
  FullMatrix<double> &                interpolation_matrix,
  const unsigned int                  face_no) const
{
  Assert(interpolation_matrix.m() == source_fe.n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.m(),
                              source_fe.n_dofs_per_face(face_no)));

  // see if source is a P or Q element
  if ((dynamic_cast<const FE_SimplexPoly<dim, spacedim> *>(&source_fe) !=
       nullptr) ||
      (dynamic_cast<const FE_Q_Base<dim, spacedim> *>(&source_fe) != nullptr))
    {
      const Quadrature<dim - 1> quad_face_support(
        source_fe.get_unit_face_support_points(face_no));

      const double eps = 2e-13 * this->degree * (dim - 1);

      std::vector<Point<dim>> face_quadrature_points(quad_face_support.size());
      QProjector<dim>::project_to_face(this->reference_cell(),
                                       quad_face_support,
                                       face_no,
                                       face_quadrature_points);

      for (unsigned int i = 0; i < source_fe.n_dofs_per_face(face_no); ++i)
        for (unsigned int j = 0; j < this->n_dofs_per_face(face_no); ++j)
          {
            double matrix_entry =
              this->shape_value(this->face_to_cell_index(j, 0),
                                face_quadrature_points[i]);

            // Correct the interpolated value. I.e. if it is close to 1 or
            // 0, make it exactly 1 or 0. Unfortunately, this is required to
            // avoid problems with higher order elements.
            if (std::fabs(matrix_entry - 1.0) < eps)
              matrix_entry = 1.0;
            if (std::fabs(matrix_entry) < eps)
              matrix_entry = 0.0;

            interpolation_matrix(i, j) = matrix_entry;
          }

#ifdef DEBUG
      for (unsigned int j = 0; j < source_fe.n_dofs_per_face(face_no); ++j)
        {
          double sum = 0.;

          for (unsigned int i = 0; i < this->n_dofs_per_face(face_no); ++i)
            sum += interpolation_matrix(j, i);

          Assert(std::fabs(sum - 1) < eps, ExcInternalError());
        }
#endif
    }
  else if (dynamic_cast<const FE_Nothing<dim> *>(&source_fe) != nullptr)
    {
      // nothing to do here, the FE_Nothing has no degrees of freedom anyway
    }
  else
    AssertThrow(
      false,
      (typename FiniteElement<dim,
                              spacedim>::ExcInterpolationNotImplemented()));
}



template <int dim, int spacedim>
void
FE_SimplexPoly<dim, spacedim>::get_subface_interpolation_matrix(
  const FiniteElement<dim, spacedim> &source_fe,
  const unsigned int                  subface,
  FullMatrix<double> &                interpolation_matrix,
  const unsigned int                  face_no) const
{
  Assert(interpolation_matrix.m() == source_fe.n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.m(),
                              source_fe.n_dofs_per_face(face_no)));

  // see if source is a P or Q element
  if ((dynamic_cast<const FE_SimplexPoly<dim, spacedim> *>(&source_fe) !=
       nullptr) ||
      (dynamic_cast<const FE_Q_Base<dim, spacedim> *>(&source_fe) != nullptr))
    {
      const Quadrature<dim - 1> quad_face_support(
        source_fe.get_unit_face_support_points(face_no));

      const double eps = 2e-13 * this->degree * (dim - 1);

      std::vector<Point<dim>> subface_quadrature_points(
        quad_face_support.size());
      QProjector<dim>::project_to_subface(this->reference_cell(),
                                          quad_face_support,
                                          face_no,
                                          subface,
                                          subface_quadrature_points);

      for (unsigned int i = 0; i < source_fe.n_dofs_per_face(face_no); ++i)
        for (unsigned int j = 0; j < this->n_dofs_per_face(face_no); ++j)
          {
            double matrix_entry =
              this->shape_value(this->face_to_cell_index(j, 0),
                                subface_quadrature_points[i]);

            // Correct the interpolated value. I.e. if it is close to 1 or
            // 0, make it exactly 1 or 0. Unfortunately, this is required to
            // avoid problems with higher order elements.
            if (std::fabs(matrix_entry - 1.0) < eps)
              matrix_entry = 1.0;
            if (std::fabs(matrix_entry) < eps)
              matrix_entry = 0.0;

            interpolation_matrix(i, j) = matrix_entry;
          }

#ifdef DEBUG
      for (unsigned int j = 0; j < source_fe.n_dofs_per_face(face_no); ++j)
        {
          double sum = 0.;

          for (unsigned int i = 0; i < this->n_dofs_per_face(face_no); ++i)
            sum += interpolation_matrix(j, i);

          Assert(std::fabs(sum - 1) < eps, ExcInternalError());
        }
#endif
    }
  else if (dynamic_cast<const FE_Nothing<dim> *>(&source_fe) != nullptr)
    {
      // nothing to do here, the FE_Nothing has no degrees of freedom anyway
    }
  else
    AssertThrow(
      false,
      (typename FiniteElement<dim,
                              spacedim>::ExcInterpolationNotImplemented()));
}



template <int dim, int spacedim>
bool
FE_SimplexPoly<dim, spacedim>::hp_constraints_are_implemented() const
{
  return true;
}



template <int dim, int spacedim>
void
FE_SimplexPoly<dim, spacedim>::
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const
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
FE_SimplexP<dim, spacedim>::FE_SimplexP(const unsigned int degree)
  : FE_SimplexPoly<dim, spacedim>(
      BarycentricPolynomials<dim>::get_fe_p_basis(degree),
      FiniteElementData<dim>(get_dpo_vector_fe_p(dim, degree),
                             ReferenceCells::get_simplex<dim>(),
                             1,
                             degree,
                             FiniteElementData<dim>::H1),
      unit_support_points_fe_p<dim>(degree),
      unit_face_support_points_fe_p<dim>(degree, FiniteElementData<dim>::H1),
      constraints_fe_p<dim>(degree))
{
  if (degree > 2)
    for (unsigned int i = 0; i < this->n_dofs_per_line(); ++i)
      this->adjust_line_dof_index_for_line_orientation_table[i] =
        this->n_dofs_per_line() - 1 - i - i;
  // We do not support multiple DoFs per quad yet
  Assert(degree <= 3, ExcNotImplemented());
}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_SimplexP<dim, spacedim>::clone() const
{
  return std::make_unique<FE_SimplexP<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
std::string
FE_SimplexP<dim, spacedim>::get_name() const
{
  std::ostringstream namebuf;
  namebuf << "FE_SimplexP<" << dim << ">(" << this->degree << ")";

  return namebuf.str();
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_SimplexP<dim, spacedim>::compare_for_domination(
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
  if (const FE_SimplexP<dim, spacedim> *fe_p_other =
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

  Assert(false, ExcNotImplemented());
  return FiniteElementDomination::neither_element_dominates;
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_SimplexP<dim, spacedim>::hp_vertex_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  AssertDimension(dim, 2);

  if (dynamic_cast<const FE_SimplexP<dim, spacedim> *>(&fe_other) != nullptr)
    {
      // there should be exactly one single DoF of each FE at a vertex, and
      // they should have identical value
      return {{0U, 0U}};
    }
  else if (dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other) != nullptr)
    {
      // there should be exactly one single DoF of each FE at a vertex, and
      // they should have identical value
      return {{0U, 0U}};
    }
  else if (dynamic_cast<const FE_Nothing<dim> *>(&fe_other) != nullptr)
    {
      // the FE_Nothing has no degrees of freedom, so there are no
      // equivalencies to be recorded
      return {};
    }
  else if (fe_other.n_unique_faces() == 1 && fe_other.n_dofs_per_face(0) == 0)
    {
      // if the other element has no DoFs on faces at all,
      // then it would be impossible to enforce any kind of
      // continuity even if we knew exactly what kind of element
      // we have -- simply because the other element declares
      // that it is discontinuous because it has no DoFs on
      // its faces. in that case, just state that we have no
      // constraints to declare
      return {};
    }
  else
    {
      Assert(false, ExcNotImplemented());
      return {};
    }
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_SimplexP<dim, spacedim>::hp_line_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  AssertDimension(dim, 2);

  if (const FE_SimplexP<dim, spacedim> *fe_p_other =
        dynamic_cast<const FE_SimplexP<dim, spacedim> *>(&fe_other))
    {
      // dofs are located along lines, so two dofs are identical if they are
      // located at identical positions.
      // Therefore, read the points in unit_support_points for the
      // first coordinate direction. For FE_SimplexP, they are currently
      // hard-coded and we iterate over points on the first line which begin
      // after the 3 vertex points in the complete list of unit support points

      std::vector<std::pair<unsigned int, unsigned int>> identities;

      for (unsigned int i = 0; i < this->degree - 1; ++i)
        for (unsigned int j = 0; j < fe_p_other->degree - 1; ++j)
          if (std::fabs(this->unit_support_points[i + 3][0] -
                        fe_p_other->unit_support_points[i + 3][0]) < 1e-14)
            identities.emplace_back(i, j);
          else
            {
              // If nodes are not located in the same place, we have to
              // interpolate. This is then not handled through the
              // current function, but via interpolation matrices that
              // result in constraints, rather than identities. Since
              // that happens in a different function, there is nothing
              // for us to do here.
            }

      return identities;
    }
  else if (const FE_Q<dim, spacedim> *fe_q_other =
             dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other))
    {
      // dofs are located along lines, so two dofs are identical if they are
      // located at identical positions. if we had only equidistant points, we
      // could simply check for similarity like (i+1)*q == (j+1)*p, but we
      // might have other support points (e.g. Gauss-Lobatto
      // points). Therefore, read the points in unit_support_points for the
      // first coordinate direction. For FE_Q, we take the lexicographic
      // ordering of the line support points in the first direction (i.e.,
      // x-direction), which we access between index 1 and p-1 (index 0 and p
      // are vertex dofs). For FE_SimplexP, they are currently hard-coded and we
      // iterate over points on the first line which begin after the 3 vertex
      // points in the complete list of unit support points

      const std::vector<unsigned int> &index_map_inverse_q_other =
        fe_q_other->get_poly_space_numbering_inverse();

      std::vector<std::pair<unsigned int, unsigned int>> identities;

      for (unsigned int i = 0; i < this->degree - 1; ++i)
        for (unsigned int j = 0; j < fe_q_other->degree - 1; ++j)
          if (std::fabs(this->unit_support_points[i + 3][0] -
                        fe_q_other->get_unit_support_points()
                          [index_map_inverse_q_other[j + 1]][0]) < 1e-14)
            identities.emplace_back(i, j);
          else
            {
              // If nodes are not located in the same place, we have to
              // interpolate. This will then also
              // capture the case where the FE_Q has a different polynomial
              // degree than the current element. In either case, the resulting
              // constraints are computed elsewhere, rather than via the
              // identities this function returns: Since
              // that happens in a different function, there is nothing
              // for us to do here.
            }

      return identities;
    }
  else if (dynamic_cast<const FE_Nothing<dim> *>(&fe_other) != nullptr)
    {
      // The FE_Nothing has no degrees of freedom, so there are no
      // equivalencies to be recorded. (If the FE_Nothing is dominating,
      // then this will also leads to constraints, but we are not concerned
      // with this here.)
      return {};
    }
  else if (fe_other.n_unique_faces() == 1 && fe_other.n_dofs_per_face(0) == 0)
    {
      // if the other element has no elements on faces at all,
      // then it would be impossible to enforce any kind of
      // continuity even if we knew exactly what kind of element
      // we have -- simply because the other element declares
      // that it is discontinuous because it has no DoFs on
      // its faces. in that case, just state that we have no
      // constraints to declare
      return {};
    }
  else
    {
      Assert(false, ExcNotImplemented());
      return {};
    }
}



template <int dim, int spacedim>
FE_SimplexDGP<dim, spacedim>::FE_SimplexDGP(const unsigned int degree)
  : FE_SimplexPoly<dim, spacedim>(
      BarycentricPolynomials<dim>::get_fe_p_basis(degree),
      FiniteElementData<dim>(get_dpo_vector_fe_dgp(dim, degree),
                             ReferenceCells::get_simplex<dim>(),
                             1,
                             degree,
                             FiniteElementData<dim>::L2),
      unit_support_points_fe_p<dim>(degree),
      unit_face_support_points_fe_p<dim>(degree, FiniteElementData<dim>::L2),
      constraints_fe_p<dim>(degree))
{}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_SimplexDGP<dim, spacedim>::clone() const
{
  return std::make_unique<FE_SimplexDGP<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
std::string
FE_SimplexDGP<dim, spacedim>::get_name() const
{
  std::ostringstream namebuf;
  namebuf << "FE_SimplexDGP<" << dim << ">(" << this->degree << ")";

  return namebuf.str();
}


template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_SimplexDGP<dim, spacedim>::compare_for_domination(
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
  if (const FE_SimplexDGP<dim, spacedim> *fe_dgp_other =
        dynamic_cast<const FE_SimplexDGP<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_dgp_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_dgp_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_DGQ<dim, spacedim> *fe_dgq_other =
             dynamic_cast<const FE_DGQ<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_dgq_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_dgq_other->degree)
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

  Assert(false, ExcNotImplemented());
  return FiniteElementDomination::neither_element_dominates;
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_SimplexDGP<dim, spacedim>::hp_vertex_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  (void)fe_other;

  return {};
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_SimplexDGP<dim, spacedim>::hp_line_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  (void)fe_other;

  return {};
}

// explicit instantiations
#include "fe_simplex_p.inst"

DEAL_II_NAMESPACE_CLOSE

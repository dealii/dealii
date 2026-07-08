// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2020 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/types.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_orientation.h>

DEAL_II_NAMESPACE_OPEN

namespace
{
  /**
   * Helper function to set up the dpo vector of FE_SimplexP for a given @p dim and
   * @p degree.
   */
  std::vector<unsigned int>
  get_dpo_vector_fe_p(const unsigned int dim, const unsigned int degree)
  {
    Assert(degree != 0, ExcNotImplemented());

    switch (dim)
      {
        case 1:
          return {1, degree - 1};
        case 2:
          // the number of support points on the face is
          // \sum_{i=1}^{degree - 2} i = (degree-2)*(degree-1)/2
          return {1, degree - 1, (degree - 2) * (degree - 1) / 2};
        case 3:
          // the number of support points in the volume are that of a tet
          // with a lower degree (degree-4)
          return {1,
                  degree - 1,
                  (degree - 2) * (degree - 1) / 2,
                  (degree - 3) * (degree - 2) * (degree - 1) / 6};
      }

    DEAL_II_ASSERT_UNREACHABLE();
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
    std::vector<Point<dim>> unit_points;
    const auto              reference_cell = ReferenceCells::get_simplex<dim>();

    // Piecewise constants are a special case: use a support point at the
    // centroid and only the centroid
    if (degree == 0)
      {
        unit_points.emplace_back(reference_cell.barycenter());
        return unit_points;
      }

    // otherwise write everything as linear combinations of vertices
    const auto dpo = get_dpo_vector_fe_p(dim, degree);
    Assert(dpo.size() == dim + 1, ExcInternalError());
    Assert(dpo[0] == 1, ExcNotImplemented());

    // vertices:
    for (const unsigned int d : reference_cell.vertex_indices())
      unit_points.push_back(reference_cell.vertex(d));

    // lines:
    for (const unsigned int l : reference_cell.line_indices())
      {
        const Point<dim> p0 =
          unit_points[reference_cell.line_to_cell_vertices(l, 0)];
        const Point<dim> p1 =
          unit_points[reference_cell.line_to_cell_vertices(l, 1)];
        for (unsigned int p = 0; p < dpo[1]; ++p)
          unit_points.push_back((double(dpo[1] - p) / (dpo[1] + 1)) * p0 +
                                (double(p + 1) / (dpo[1] + 1)) * p1);
      }

    // faces:
    if constexpr (dim == 2)
      {
        unsigned int counter = 0;
        for (unsigned int i = 1; i < degree; ++i)
          for (unsigned int j = 1; j < degree - i; ++j, ++counter)
            {
              const double x = static_cast<double>(j) / degree;
              const double y = static_cast<double>(i) / degree;

              unit_points.push_back(Point<dim>(x, y));
            }
        Assert(counter == dpo[2], ExcInternalError());
      }

    if constexpr (dim == 3)
      for (const unsigned int f : reference_cell.face_indices())
        {
          const Point<dim> p0 =
            unit_points[reference_cell.face_to_cell_vertices(
              f, 0, numbers::default_geometric_orientation)];
          const Point<dim> p1 =
            unit_points[reference_cell.face_to_cell_vertices(
              f, 1, numbers::default_geometric_orientation)];
          const Point<dim> p2 =
            unit_points[reference_cell.face_to_cell_vertices(
              f, 2, numbers::default_geometric_orientation)];

          unsigned int counter = 0;
          for (unsigned int i = 1; i < degree; ++i)
            for (unsigned int j = 1; j < degree - i; ++j, ++counter)
              {
                const double a = static_cast<double>(j) / degree;
                const double b = static_cast<double>(i) / degree;
                const double c = 1.0 - a - b;
                unit_points.push_back(c * p0 + a * p1 + b * p2);
              }
          Assert(counter == dpo[2], ExcInternalError());
        }

    // interior
    if constexpr (dim == 3)
      {
        unsigned int counter = 0;
        for (unsigned int i = 1; i < degree; ++i)
          for (unsigned int j = 1; j < degree - i; ++j)
            for (unsigned int k = 1; k < degree - i - j; ++k, ++counter)
              {
                const double x = static_cast<double>(i) / degree;
                const double y = static_cast<double>(j) / degree;
                const double z = static_cast<double>(k) / degree;

                unit_points.push_back(Point<dim>(x, y, z));
              }
        Assert(counter == dpo[3], ExcInternalError());
      }

    return unit_points;
  }

  template <>
  std::vector<Point<0>>
  unit_support_points_fe_p(const unsigned int /*degree*/)
  {
    return {Point<0>()};
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
    constexpr int dim = 2;

    // the following implements the 2d case
    // (the 3d case is not implemented yet)
    //
    // consult FE_Q_Base::Implementation::initialize_constraints()
    // for more information

    std::vector<Point<dim - 1>> constraint_points;
    // midpoint
    constraint_points.emplace_back(0.5);
    // subface 0
    for (unsigned int i = 1; i < degree; ++i)
      constraint_points.push_back(
        GeometryInfo<dim - 1>::child_to_cell_coordinates(
          Point<dim - 1>(i / double(degree)), 0));
    // subface 1
    for (unsigned int i = 1; i < degree; ++i)
      constraint_points.push_back(
        GeometryInfo<dim - 1>::child_to_cell_coordinates(
          Point<dim - 1>(i / double(degree)), 1));

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

        DEAL_II_NOT_IMPLEMENTED();
        return {};
      }
  }
} // namespace



template <int dim, int spacedim>
FE_SimplexPoly<dim, spacedim>::FE_SimplexPoly(
  const BarycentricPolynomials<dim>              polynomials,
  const FiniteElementData<dim>                  &fe_data,
  const bool                                     prolongation_is_additive,
  const std::vector<Point<dim>>                 &unit_support_points,
  const std::vector<std::vector<Point<dim - 1>>> unit_face_support_points,
  const FullMatrix<double>                      &interface_constraints)
  : dealii::FE_Poly<dim, spacedim>(
      polynomials,
      fe_data,
      std::vector<bool>(fe_data.dofs_per_cell, prolongation_is_additive),
      std::vector<ComponentMask>(fe_data.dofs_per_cell,
                                 ComponentMask(std::vector<bool>(1, true))))
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
  if (dim == 3)
    Assert(RefinementCase<dim>(refinement_case) ==
               RefinementCase<dim>(
                 static_cast<char>(IsotropicRefinementChoice::cut_tet_68)) ||
             RefinementCase<dim>(refinement_case) ==
               RefinementCase<dim>(
                 static_cast<char>(IsotropicRefinementChoice::cut_tet_57)) ||
             RefinementCase<dim>(refinement_case) ==
               RefinementCase<dim>(
                 static_cast<char>(IsotropicRefinementChoice::cut_tet_49)),
           ExcNotImplemented());
  else
    Assert(refinement_case ==
             RefinementCase<dim>(RefinementCase<dim>::isotropic_refinement),
           ExcNotImplemented());
  AssertDimension(dim, spacedim);

  // initialization upon first request
  if (this->prolongation[refinement_case - 1][child].n() == 0)
    {
      std::lock_guard<std::mutex> lock(prolongation_matrix_mutex);

      // if matrix got updated while waiting for the lock
      if (this->prolongation[refinement_case - 1][child].n() ==
          this->n_dofs_per_cell())
        return this->prolongation[refinement_case - 1][child];

      // now do the work. need to get a non-const version of data in order to
      // be able to modify them inside a const function
      auto &this_nonconst = const_cast<FE_SimplexPoly<dim, spacedim> &>(*this);

      if (dim == 2)
        {
          std::vector<std::vector<FullMatrix<double>>> isotropic_matrices(
            RefinementCase<dim>::isotropic_refinement);
          isotropic_matrices.back().resize(
            this->reference_cell().n_children(
              RefinementCase<dim>(refinement_case)),
            FullMatrix<double>(this->n_dofs_per_cell(),
                               this->n_dofs_per_cell()));

          FETools::compute_embedding_matrices(*this, isotropic_matrices, true);

          this_nonconst.prolongation[refinement_case - 1] =
            std::move(isotropic_matrices.back());
        }
      else if (dim == 3)
        {
          std::vector<std::vector<FullMatrix<double>>> matrices(
            static_cast<unsigned int>(IsotropicRefinementChoice::cut_tet_49),
            std::vector<FullMatrix<double>>(
              this->reference_cell().n_children(
                RefinementCase<dim>(refinement_case)),
              FullMatrix<double>(this->n_dofs_per_cell(),
                                 this->n_dofs_per_cell())));
          FETools::compute_embedding_matrices(*this, matrices, true);
          for (unsigned int refinement_direction = static_cast<unsigned int>(
                 IsotropicRefinementChoice::cut_tet_68);
               refinement_direction <=
               static_cast<unsigned int>(IsotropicRefinementChoice::cut_tet_49);
               refinement_direction++)
            this_nonconst.prolongation[refinement_direction - 1] =
              std::move(matrices[refinement_direction - 1]);
        }
      else
        DEAL_II_ASSERT_UNREACHABLE();
    }

  // finally return the matrix
  return this->prolongation[refinement_case - 1][child];
}



template <int dim, int spacedim>
unsigned int
FE_SimplexPoly<dim, spacedim>::face_to_cell_index(
  const unsigned int                 face_dof_index,
  const unsigned int                 face,
  const types::geometric_orientation combined_orientation) const
{
  return FETools::face_to_cell_index(*this,
                                     face_dof_index,
                                     face,
                                     combined_orientation);
}



template <int dim, int spacedim>
const FullMatrix<double> &
FE_SimplexPoly<dim, spacedim>::get_restriction_matrix(
  const unsigned int         child,
  const RefinementCase<dim> &refinement_case) const
{
  if (dim == 3)
    Assert(RefinementCase<dim>(refinement_case) ==
               RefinementCase<dim>(
                 static_cast<char>(IsotropicRefinementChoice::cut_tet_68)) ||
             RefinementCase<dim>(refinement_case) ==
               RefinementCase<dim>(
                 static_cast<char>(IsotropicRefinementChoice::cut_tet_57)) ||
             RefinementCase<dim>(refinement_case) ==
               RefinementCase<dim>(
                 static_cast<char>(IsotropicRefinementChoice::cut_tet_49)),
           ExcNotImplemented());
  else
    Assert(refinement_case == RefinementCase<dim>::isotropic_refinement,
           ExcNotImplemented());
  AssertDimension(dim, spacedim);

  // initialization upon first request
  if (this->restriction[refinement_case - 1][child].n() == 0)
    {
      std::lock_guard<std::mutex> lock(restriction_matrix_mutex);

      // if matrix got updated while waiting for the lock
      if (this->restriction[refinement_case - 1][child].n() ==
          this->n_dofs_per_cell())
        return this->restriction[refinement_case - 1][child];

      // get the restriction matrix
      // Refine a unit cell. As the parent cell is a unit
      // cell, the reference cell of the children equals the parent, i.e. they
      // have the support points at the same locations. So we just have to check
      // if a support point of the parent is one of the interpolation points of
      // the child. If this is not the case we find the interpolation of the
      // point.

      const double       eps = 1e-12;
      FullMatrix<double> restriction_mat(this->n_dofs_per_cell(),
                                         this->n_dofs_per_cell());

      // first get all support points on the reference cell
      const std::vector<Point<dim>> unit_support_points =
        this->get_unit_support_points();

      // now create children on the reference cell
      Triangulation<dim> tria;
      GridGenerator::reference_cell(tria, this->reference_cell());
      tria.begin_active()->set_refine_flag(
        RefinementCase<dim>::isotropic_refinement);
      if (dim == 3)
        tria.begin_active()->set_refine_choice(refinement_case);
      tria.execute_coarsening_and_refinement();

      const auto &child_cell = tria.begin(0)->child(child);

      // iterate over all support points and transform them to the unit cell of
      // the child
      for (unsigned int i = 0; i < unit_support_points.size(); i++)
        {
          std::vector<Point<dim>>            transformed_point(1);
          const std::vector<Point<spacedim>> unit_support_point = {
            dim == 2 ? Point<spacedim>(unit_support_points[i][0],
                                       unit_support_points[i][1]) :
                       Point<spacedim>(unit_support_points[i][0],
                                       unit_support_points[i][1],
                                       unit_support_points[i][2])};
          this->reference_cell()
            .template get_default_linear_mapping<spacedim>()
            .transform_points_real_to_unit_cell(
              child_cell,
              make_array_view(unit_support_point),
              make_array_view(transformed_point));

          // if point is inside the unit cell iterate over all shape functions
          if (this->reference_cell().contains_point(transformed_point[0], eps))
            for (unsigned int j = 0; j < this->n_dofs_per_cell(); j++)
              restriction_mat[i][j] =
                this->shape_value(j, transformed_point[0]);
        }
      if constexpr (running_in_debug_mode())
        {
          for (unsigned int i = 0; i < this->n_dofs_per_cell(); i++)
            {
              double sum = 0.;

              for (unsigned int j = 0; j < this->n_dofs_per_cell(); j++)
                sum += restriction_mat[i][j];

              Assert(std::fabs(sum - 1) < eps || std::fabs(sum) < eps,
                     ExcInternalError(
                       "The entries in a row of the local "
                       "restriction matrix do not add to zero or one. "
                       "This typically indicates that the "
                       "polynomial interpolation is "
                       "ill-conditioned such that round-off "
                       "prevents the sum to be one."));
            }
        }

      // Remove small entries from the matrix
      for (unsigned int i = 0; i < restriction_mat.m(); ++i)
        for (unsigned int j = 0; j < restriction_mat.n(); ++j)
          {
            if (std::fabs(restriction_mat(i, j)) < eps)
              restriction_mat(i, j) = 0.;
            if (std::fabs(restriction_mat(i, j) - 1) < eps)
              restriction_mat(i, j) = 1.;
          }

      const_cast<FullMatrix<double> &>(
        this->restriction[refinement_case - 1][child]) =
        std::move(restriction_mat);
    }

  // finally return the matrix
  return this->restriction[refinement_case - 1][child];
}



template <int dim, int spacedim>
void
FE_SimplexPoly<dim, spacedim>::get_face_interpolation_matrix(
  const FiniteElement<dim, spacedim> &source_fe,
  FullMatrix<double>                 &interpolation_matrix,
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

      const std::vector<Point<dim>> face_quadrature_points =
        QProjector<dim>::project_to_face(this->reference_cell(),
                                         quad_face_support,
                                         face_no,
                                         numbers::default_geometric_orientation)
          .get_points();

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

      if constexpr (running_in_debug_mode())
        {
          for (unsigned int j = 0; j < source_fe.n_dofs_per_face(face_no); ++j)
            {
              double sum = 0.;

              for (unsigned int i = 0; i < this->n_dofs_per_face(face_no); ++i)
                sum += interpolation_matrix(j, i);

              Assert(std::fabs(sum - 1) < eps, ExcInternalError());
            }
        }
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
  FullMatrix<double>                 &interpolation_matrix,
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

      const Quadrature<dim> subface_quadrature =
        QProjector<dim>::project_to_subface(
          this->reference_cell(),
          quad_face_support,
          face_no,
          subface,
          numbers::default_geometric_orientation,
          RefinementCase<dim - 1>::isotropic_refinement);

      for (unsigned int i = 0; i < source_fe.n_dofs_per_face(face_no); ++i)
        for (unsigned int j = 0; j < this->n_dofs_per_face(face_no); ++j)
          {
            double matrix_entry =
              this->shape_value(this->face_to_cell_index(j, 0),
                                subface_quadrature.point(i));

            // Correct the interpolated value. I.e. if it is close to 1 or
            // 0, make it exactly 1 or 0. Unfortunately, this is required to
            // avoid problems with higher order elements.
            if (std::fabs(matrix_entry - 1.0) < eps)
              matrix_entry = 1.0;
            if (std::fabs(matrix_entry) < eps)
              matrix_entry = 0.0;

            interpolation_matrix(i, j) = matrix_entry;
          }

      if constexpr (running_in_debug_mode())
        {
          for (unsigned int j = 0; j < source_fe.n_dofs_per_face(face_no); ++j)
            {
              double sum = 0.;

              for (unsigned int i = 0; i < this->n_dofs_per_face(face_no); ++i)
                sum += interpolation_matrix(j, i);

              Assert(std::fabs(sum - 1) < eps, ExcInternalError());
            }
        }
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
FE_SimplexP<dim, spacedim>::FE_SimplexP(const unsigned int degree)
  : FE_SimplexPoly<dim, spacedim>(
      BarycentricPolynomials<dim>::get_fe_p_basis(degree),
      FiniteElementData<dim>(get_dpo_vector_fe_p(dim, degree),
                             ReferenceCells::get_simplex<dim>(),
                             1,
                             degree,
                             FiniteElementData<dim>::H1),
      false,
      unit_support_points_fe_p<dim>(degree),
      unit_face_support_points_fe_p<dim>(degree, FiniteElementData<dim>::H1),
      constraints_fe_p<dim>(degree))
{
  if (degree > 2)
    for (unsigned int i = 0; i < this->n_dofs_per_line(); ++i)
      this->adjust_line_dof_index_for_line_orientation_table[i] =
        this->n_dofs_per_line() - 1 - i - i;

  // for 1d and 2d or if there are no DoFs on the quads
  // we can skip adjust_quad_dof_index_for_face_orientation_table
  if (dim < 3 || degree < 3)
    return;

  // do some sanity checks
  AssertDimension(this->n_unique_faces(), 1);
  const unsigned int face_no = 0;

  Assert(
    this->adjust_quad_dof_index_for_face_orientation_table[0].n_elements() ==
      this->reference_cell().n_face_orientations(face_no) *
        this->n_dofs_per_quad(face_no),
    ExcInternalError());

  Assert((degree - 2) * (degree - 1) / 2 == this->n_dofs_per_quad(face_no),
         ExcInternalError());

  const auto face_reference_cell =
    this->reference_cell().face_reference_cell(face_no);

  // the interior nodes build a new triangle with lower degree r
  const unsigned int r = degree - 3;

  // now loop over all DoFs on the triangle
  // 0 <= i + j <= r holds on the triangle
  for (unsigned int j = 0, dof_index = 0; j <= r; ++j)
    for (unsigned int i = 0; i <= r - j; ++i, ++dof_index)
      {
        // index in the style of barycentric coordinates
        // the first entry is the remainder as i + j <= r has to hold
        const std::array<unsigned int, 3> local_indices{{r - i - j, i, j}};

        // go over all possible orientations
        for (types::geometric_orientation orientation = 0;
             orientation < this->reference_cell().n_face_orientations(face_no);
             ++orientation)
          {
            // get the correct permutation for the current orientation
            const auto permuted_indices =
              face_reference_cell.permute_by_combined_orientation(
                make_array_view(local_indices),
                face_reference_cell.get_inverse_combined_orientation(
                  orientation));

            // now reconstruct the index of from the permuted i and j
            // take orientation 0 which is the standard orientation
            // then the index k is k =  i + j*(r+1) - (j*(j-1))/2
            const unsigned int k =
              permuted_indices[1] + permuted_indices[2] * (r + 1) -
              (permuted_indices[2] * (permuted_indices[2] - 1)) / 2;

            const int offset =
              static_cast<int>(k) - static_cast<int>(dof_index);
            this->adjust_quad_dof_index_for_face_orientation_table[face_no](
              dof_index, orientation) = offset;
          }
      }
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
  namebuf << "FE_SimplexP<" << Utilities::dim_string(dim, spacedim) << ">("
          << this->degree << ")";

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
  else if (const FE_WedgeP<dim, spacedim> *fe_p_other =
             dynamic_cast<const FE_WedgeP<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_p_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_p_other->degree)
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
FE_SimplexP<dim, spacedim>::hp_vertex_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  if ((dynamic_cast<const FE_SimplexP<dim, spacedim> *>(&fe_other) !=
       nullptr) ||
      (dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other) != nullptr) ||
      (dynamic_cast<const FE_PyramidP<dim, spacedim> *>(&fe_other) !=
       nullptr) ||
      (dynamic_cast<const FE_WedgeP<dim, spacedim> *>(&fe_other) != nullptr))
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
      DEAL_II_NOT_IMPLEMENTED();
      return {};
    }
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_SimplexP<dim, spacedim>::hp_line_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  if ((dynamic_cast<const FE_SimplexP<dim, spacedim> *>(&fe_other)) ||
      (dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other)) ||
      (dynamic_cast<const FE_PyramidP<dim, spacedim> *>(&fe_other)) ||
      (dynamic_cast<const FE_WedgeP<dim, spacedim> *>(&fe_other)))
    {
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
      // face number 0 of the tet is a triangle
      const unsigned int offset =
        this->reference_cell().face_reference_cell(0).n_vertices();

      const unsigned int offset_other =
        fe_other.reference_cell().face_reference_cell(0).is_simplex() ?
          fe_other.reference_cell().face_reference_cell(0).n_vertices() :
          fe_other.reference_cell().face_reference_cell(0).n_vertices() +
            2 * fe_other.n_dofs_per_line();

      // now get the identities
      for (unsigned int i = 0; i < this->degree - 1; ++i)
        for (unsigned int j = 0; j < fe_other.degree - 1; ++j)
          if (face_support_points[i + offset].distance(
                face_support_points_other[j + offset_other]) < 1e-14)
            identities.emplace_back(i, j);

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
      DEAL_II_NOT_IMPLEMENTED();
      return {};
    }
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_SimplexP<dim, spacedim>::hp_quad_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  face_no,
  const unsigned int                  face_no_other) const
{
  AssertDimension(dim, 3);
  AssertIndexRange(face_no, 4);

  if ((dynamic_cast<const FE_SimplexP<dim> *>(&fe_other) != nullptr) ||
      (dynamic_cast<const FE_WedgeP<dim> *>(&fe_other) != nullptr) ||
      (dynamic_cast<const FE_PyramidP<dim> *>(&fe_other) != nullptr))
    {
      return FETools::hp_quad_dof_identities(*this,
                                             fe_other,
                                             face_no,
                                             face_no_other);
    }
  else if (dynamic_cast<const FE_Nothing<dim> *>(&fe_other) != nullptr)
    {
      // the FE_Nothing has no degrees of freedom, so there are no
      // equivalencies to be recorded
      return std::vector<std::pair<unsigned int, unsigned int>>();
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
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
      return std::vector<std::pair<unsigned int, unsigned int>>();
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
      true,
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
  namebuf << "FE_SimplexDGP<" << Utilities::dim_string(dim, spacedim) << ">("
          << this->degree << ")";

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

  DEAL_II_NOT_IMPLEMENTED();
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



template <int dim, int spacedim>
const FullMatrix<double> &
FE_SimplexDGP<dim, spacedim>::get_restriction_matrix(
  const unsigned int         child,
  const RefinementCase<dim> &refinement_case) const
{
  if (dim == 3)
    Assert(RefinementCase<dim>(refinement_case) ==
               RefinementCase<dim>(
                 static_cast<char>(IsotropicRefinementChoice::cut_tet_68)) ||
             RefinementCase<dim>(refinement_case) ==
               RefinementCase<dim>(
                 static_cast<char>(IsotropicRefinementChoice::cut_tet_57)) ||
             RefinementCase<dim>(refinement_case) ==
               RefinementCase<dim>(
                 static_cast<char>(IsotropicRefinementChoice::cut_tet_49)),
           ExcNotImplemented());
  else
    Assert(refinement_case == RefinementCase<dim>::isotropic_refinement,
           ExcNotImplemented());
  AssertDimension(dim, spacedim);

  // initialization upon first request
  if (this->restriction[refinement_case - 1][child].n() == 0)
    {
      std::lock_guard<std::mutex> lock(this->restriction_matrix_mutex);

      // if matrix got updated while waiting for the lock
      if (this->restriction[refinement_case - 1][child].n() ==
          this->n_dofs_per_cell())
        return this->restriction[refinement_case - 1][child];

      // now do the work. need to get a non-const version of data in order to
      // be able to modify them inside a const function
      auto &this_nonconst = const_cast<FE_SimplexDGP<dim, spacedim> &>(*this);

      if (dim == 2)
        {
          std::vector<std::vector<FullMatrix<double>>> isotropic_matrices(
            RefinementCase<dim>::isotropic_refinement);
          isotropic_matrices.back().resize(
            this->reference_cell().n_children(
              RefinementCase<dim>(refinement_case)),
            FullMatrix<double>(this->n_dofs_per_cell(),
                               this->n_dofs_per_cell()));

          FETools::compute_projection_matrices(*this, isotropic_matrices, true);

          this_nonconst.restriction[refinement_case - 1] =
            std::move(isotropic_matrices.back());
        }
      else if (dim == 3)
        {
          std::vector<std::vector<FullMatrix<double>>> matrices(
            static_cast<unsigned int>(IsotropicRefinementChoice::cut_tet_49),
            std::vector<FullMatrix<double>>(
              this->reference_cell().n_children(
                RefinementCase<dim>(refinement_case)),
              FullMatrix<double>(this->n_dofs_per_cell(),
                                 this->n_dofs_per_cell())));
          FETools::compute_projection_matrices(*this, matrices, true);
          for (unsigned int refinement_direction = static_cast<unsigned int>(
                 IsotropicRefinementChoice::cut_tet_68);
               refinement_direction <=
               static_cast<unsigned int>(IsotropicRefinementChoice::cut_tet_49);
               refinement_direction++)
            this_nonconst.restriction[refinement_direction - 1] =
              std::move(matrices[refinement_direction - 1]);
        }
      else
        DEAL_II_ASSERT_UNREACHABLE();
    }

  // finally return the matrix
  return this->restriction[refinement_case - 1][child];
}

// explicit instantiations
#include "fe/fe_simplex_p.inst"

DEAL_II_NAMESPACE_CLOSE

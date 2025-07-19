// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_closest_surface_point
#define dealii_closest_surface_point


#include <deal.II/base/config.h>

#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

DEAL_II_NAMESPACE_OPEN

namespace NonMatching
{


  /**
   * @brief A class for computing the closest points on a surface defined by a level set function.
   *
   * This class implements algorithms to find the closest points on a surface
   * (zero level set) to a given set of query points. The surface is implicitly
   * defined by a level set function discretized on a finite element mesh. The
   * class uses Newton's method to iteratively find the closest surface points
   * by minimizing the distance function subject to the constraint that the
   * point lies on the zero level set.
   *
   * - Level set functions discretized on DoFHandler with cartesian grid.
   * - Configurable Newton iteration parameters (tolerance, maximum iterations)
   * - Supports active cells and MG cells
   *
   * @tparam dim The spatial dimension of the problem
   * @tparam VECTOR The vector type used to store the level set function values
   *
   * @note The level set function should be negative inside the domain, positive outside,
   *       and zero on the surface boundary.
   *
   * @warning The Newton iteration may not converge for points that are too far from the surface
   *          or in regions where the level set function has poor conditioning.
   */
  template <int dim, class VECTOR>
  class ClosestSurfacePoint : public EnableObserverPointer
  {
  public:
    /**
     * Additional data for the closest surface point computation.
     */
    struct AdditionalData
    {
      /**
       * The level of the DoFHandler to use. If set to
       * numbers::invalid_unsigned_int, the finest level is used.
       */
      unsigned int level = numbers::invalid_unsigned_int;
      /**
       * The tolerance for the Newton iteration.
       */
      double tolerance = 1e-10;
      /**
       * The maximum number of iterations for the Newton iteration.
       */
      unsigned int n_iterations = 20;
    };

    using VectorType = VECTOR;

    /**
     * Constructor for the ClosestSurfacePoint class.
     */
    ClosestSurfacePoint(const VectorType      &level_set,
                        const DoFHandler<dim> &dof_handler,
                        const AdditionalData  &data = AdditionalData());



    /**
     * @brief Compute closest points to given set of points.
     *
     * @return A pair of vectors, where the first vector contains the shifted quadrature points in absolute
     * coordinates, and the second vector contains the corresponding shifted
     * points in coordinates of the reference cell.
     *
     * @param search_cell The cell in the triangulation where the search for closest surface points is performed.
     * @param reference_cell The reference cell in which local coordinates unit points are outputed.
     * @param quadrature_points The original quadrature points to be shifted.
     */
    std::pair<std::vector<Point<dim>>, std::vector<Point<dim>>>
    compute_closest_surface_points(
      const typename Triangulation<dim>::cell_iterator &search_cell,
      const typename Triangulation<dim>::cell_iterator &reference_cell,
      const std::vector<Point<dim>> &quadrature_points) const;

  private:
    AdditionalData         data;
    const DoFHandler<dim> &dof_handler;
    const VectorType      &level_set;
    MappingCartesian<dim>  mapping;



    /**
     * Find the closest point on a surface using Newton's method with a
     * monolithic approach.
     *
     * This function implements Newton's method to find the point on the surface
     * defined by the finite element and DOF values that is closest to the given
     * input point. The monolithic approach solves the constrained optimization
     * problem:
     *
     * @f[
     * \min_{x} \frac{1}{2} \|x - x_0 \|^2 \quad \text{subject to} \quad \phi(x)
     * = 0
     * @f]
     *
     * where @f$\phi(x)@f$ is the level set function. Using Lagrange
     * multipliers, this becomes the unconstrained problem:
     *
     * @f[
     * \min_{x, \lambda} L(x, \lambda) = \frac{1}{2} \|x - x_0\|^2 + \lambda
     * \phi(x)
     * @f]
     *
     * Newton's method iteratively solves this system using the Hessian.
     *
     * @param[in] point The reference point for which to find the closest
     * surface point
     * @param[in] fe The finite element used to define the surface geometry
     * @param[in] dof_values The degrees of freedom values that define the
     * surface
     * @param[out] closest_point The computed closest point on the surface
     */
    void
    newton_monolithic(const Point<dim>          &point,
                      const FiniteElement<dim>  &fe,
                      const std::vector<double> &dof_values,
                      Point<dim>                &closest_point) const;
  };



  template <int dim, class VECTOR>
  ClosestSurfacePoint<dim, VECTOR>::ClosestSurfacePoint(
    const VectorType      &level_set,
    const DoFHandler<dim> &dof_handler,
    const AdditionalData  &data)
    : data(data)
    , dof_handler(dof_handler)
    , level_set(level_set)
  {
    if (data.level != numbers::invalid_unsigned_int)
      {
        AssertThrow(data.level <
                      dof_handler.get_triangulation().n_global_levels(),
                    dealii::ExcMessage("Level is larger than number of levels "
                                       "in the Triangulation"));
      }
  }



  template <int dim, class VECTOR>
  std::pair<std::vector<Point<dim>>, std::vector<Point<dim>>>
  ClosestSurfacePoint<dim, VECTOR>::compute_closest_surface_points(
    const typename Triangulation<dim>::cell_iterator &search_cell,
    const typename Triangulation<dim>::cell_iterator &reference_cell,
    const std::vector<Point<dim>>                    &quadrature_points) const
  {
    std::vector<Point<dim>> closest_unit_search_points(quadrature_points);
    for (unsigned int q = 0; q < quadrature_points.size(); ++q)
      closest_unit_search_points[q] =
        mapping.transform_real_to_unit_cell(search_cell, quadrature_points[q]);

    std::vector<double> dof_values_level_set(
      dof_handler.get_fe().dofs_per_cell);
    std::vector<types::global_dof_index> level_set_dof_indices(
      dof_handler.get_fe().dofs_per_cell);

    typename DoFHandler<dim>::cell_iterator dof_cell(
      &search_cell->get_triangulation(),
      search_cell->level(),
      search_cell->index(),
      &dof_handler);

    if (data.level != numbers::invalid_unsigned_int)
      dof_cell->get_mg_dof_indices(level_set_dof_indices);
    else
      dof_cell->get_dof_indices(level_set_dof_indices);

    level_set.extract_subvector_to(level_set_dof_indices.begin(),
                                   level_set_dof_indices.end(),
                                   dof_values_level_set.begin());

    for (size_t i = 0; i < closest_unit_search_points.size(); ++i)
      {
        newton_monolithic(closest_unit_search_points[i],
                          dof_handler.get_fe(),
                          dof_values_level_set,
                          closest_unit_search_points[i]);
      }
    std::vector<Point<dim>> closest_real_points(quadrature_points.size());
    std::vector<Point<dim>> closest_unit_reference_points(
      quadrature_points.size());
    // back to absolute coordinates
    for (unsigned int q = 0; q < quadrature_points.size(); ++q)
      closest_real_points[q] =
        mapping.transform_unit_to_real_cell(search_cell,
                                            closest_unit_search_points[q]);

    for (unsigned int q = 0; q < quadrature_points.size(); ++q)
      closest_unit_reference_points[q] =
        mapping.transform_real_to_unit_cell(reference_cell,
                                            closest_real_points[q]);

    return {closest_real_points, closest_unit_reference_points};
  }



  template <int dim, class VECTOR>
  void
  ClosestSurfacePoint<dim, VECTOR>::newton_monolithic(
    const Point<dim>          &point,
    const FiniteElement<dim>  &fe,
    const std::vector<double> &dof_values,
    Point<dim>                &closest_point) const
  {
    AssertDimension(dof_values.size(), fe.dofs_per_cell);

    Assert(fe.degree > 1,
           dealii::ExcMessage(
             " The Newton iteration to find closest surface points "
             "requires hessians that are not available when the finite element degree is 1."));

    // X, Y, Z, lambda
    Vector<double> current_solution(dim + 1);
    Vector<double> solution_update(dim + 1);

    Vector<double>     residual(dim + 1);
    FullMatrix<double> hessian(dim + 1, dim + 1);

    for (unsigned int i = 0; i < dim; ++i)
      current_solution[i] = closest_point[i];


    for (unsigned int newton_iter = 0; newton_iter < data.n_iterations;
         ++newton_iter)
      {
        hessian             = 0.0;
        residual            = 0.0;
        const double lambda = current_solution[dim];
        for (unsigned int k = 0; k < dof_values.size(); ++k)
          {
            const auto value_k = fe.shape_value(k, closest_point);
            const auto grad_k  = fe.shape_grad(k, closest_point);
            const auto hess_k  = fe.shape_grad_grad(k, closest_point);
            for (unsigned int i = 0; i < dim; ++i)
              {
                for (unsigned int j = 0; j < dim; ++j)
                  hessian(i, j) += lambda * dof_values[k] * hess_k[i][j];

                hessian(i, dim) += dof_values[k] * grad_k[i];
                hessian(dim, i) += dof_values[k] * grad_k[i];
              }

            for (unsigned int i = 0; i < dim; ++i)
              residual[i] -= lambda * dof_values[k] * grad_k[i];

            residual[dim] -= dof_values[k] * value_k;
          }


        for (unsigned int i = 0; i < dim; ++i)
          {
            residual[i] -= current_solution[i] - point[i];
            hessian[i][i] += 1.0;
          }

        if (residual.l2_norm() < data.tolerance)
          break;


        hessian.gauss_jordan();
        hessian.vmult(solution_update, residual);
        current_solution += solution_update;

        for (unsigned int i = 0; i < dim; ++i)
          closest_point[i] = current_solution[i];
      }

    Assert(residual.l2_norm() < 1e-10,
           dealii::ExcMessage("Newton iteration did not converge"));
  }

} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE

#endif

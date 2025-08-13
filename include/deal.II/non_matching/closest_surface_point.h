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

#include <deal.II/base/point.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/lac/read_vector.h>

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
   * @tparam Number The scalar type of vector used to store the level set function values
   *
   * @note The level set function should be negative inside the domain, positive outside,
   *       and zero on the surface boundary.
   *
   * @warning The Newton iteration may not converge for points that are too far from the surface
   *          or in regions where the level set function has poor conditioning.
   */
  template <int dim, class Number>
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
    /**
     * Constructor for the ClosestSurfacePoint class.
     */
    ClosestSurfacePoint(const ReadVector<Number> &level_set,
                        const DoFHandler<dim>    &dof_handler,
                        Mapping<dim>             &mapping,
                        const AdditionalData     &data = AdditionalData());



    /**
     * Compute closest points to given set of points.
     *
     * @return A pair of vectors, where the first vector contains the shifted quadrature points in absolute
     * coordinates, and the second vector contains the corresponding shifted
     * points in coordinates of the reference cell.
     *
     * @param search_cell The cell in the triangulation where the search for closest surface points is performed.
     * @param reference_cell The reference cell in which local coordinates unit points are outputted.
     * @param quadrature_points The original quadrature points to be shifted.
     */
    std::pair<std::vector<Point<dim>>, std::vector<Point<dim>>>
    compute_closest_surface_points(
      const typename Triangulation<dim>::cell_iterator &search_cell,
      const typename Triangulation<dim>::cell_iterator &reference_cell,
      const std::vector<Point<dim>> &quadrature_points) const;

  private:
    AdditionalData                            data;
    ObserverPointer<const DoFHandler<dim>>    dof_handler;
    ObserverPointer<const ReadVector<Number>> level_set;
    ObserverPointer<Mapping<dim>>             mapping;



    /**
     * Find the closest point on a surface using Newton's method with a
     * monolithic approach.
     *
     * This function implements Newton's method to find the point on the surface
     * defined by the finite element and DOF values that is closest to the given
     * input point. The monolithic approach solves the constrained optimization
     * problem:
     *
     * \f[
     * \min_{x} \frac{1}{2} \|x - x_0 \|^2 \quad \text{subject to} \quad \phi(x)
     * = 0
     * \f]
     *
     * where $\phi(x)$ is the level set function. Using Lagrange
     * multipliers, this becomes the unconstrained problem:
     *
     * \f[
     * \min_{x, \lambda} L(x, \lambda) = \frac{1}{2} \|x - x_0\|^2 + \lambda
     * \phi(x)
     * \f]
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
                      const std::vector<Number> &dof_values,
                      Point<dim>                &closest_point) const;
  };



} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE

#endif

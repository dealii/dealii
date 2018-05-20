// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii_non_matching_immersed_surface_quadrature
#define dealii_non_matching_immersed_surface_quadrature

#include <deal.II/base/config.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN
namespace NonMatching
{
  /**
   * Defines a quadrature formula for integration over an oriented surface,
   * $\hat{S}$, immersed in the unit cell. By immersed it is meant that the
   * surface may intersect the unit cell in an arbitrary way. The quadrature
   * formula is described by a set of quadrature points, $\hat{x}_q$, weights,
   * $w_q$, and normalized surface normals, $\hat{n}_q$.
   *
   * We typically want to compute surface integrals in real space.
   * A surface $S$ intersecting a cell $K$ in real space, can be mapped onto a
   * surface $\hat{S}$ intersecting the unit cell $\hat{K}$.
   * Thus a surface integral over $S\cap K$ in real space can be transformed to
   * a surface integral over $\hat{S} \cap \hat{K}$ according to
   * @f[
   * \int_{S\cap K} f(x) dS =
   * \int_{S\cap K} f(x) |d\bar{S}| =
   * \int_{\hat{S}\cap\hat{K}} f(F_{K}(\hat{x})) \det(J) |\left( J^{-1} \right )^T d\hat{S}|,
   * @f]
   * where $F_K$ is the mapping from reference to real space and $J$ is its
   * Jacobian. This transformation is possible since the continuous surface
   * elements are vectors: $d\bar{S}, d\hat{S} \in \mathbb{R}^{dim}$ which are
   * parallel to the normals of $S$ and $\hat{S}$. So in order to compute the
   * integral in real space one needs information about the normal to do the
   * transformation.
   *
   * Thus, in addition to storing points and weights, this quadrature stores
   * also the normalized normal for each quadrature point. This can be viewed
   * as storing a discrete surface element,
   * @f[
   * \Delta \hat{S}_q := w_q \hat{n}_q \approx d\hat{S}(\hat{x}_q),
   * @f]
   * for each quadrature point. The surface integral in real space would then be
   * approximated as
   * @f[
   * \int_{S\cap K} f(x) dS \approx
   * \sum_{q} f \left(F_{K}(\hat{x}_{q}) \right) \det(J_q)
   * |\left( J_q^{-1} \right)^T \hat{n}_q| w_q.
   * @f]
   *
   * @image html immersed_surface_quadrature.svg
   *
   * @author Simon Sticko, 2017
   */
  template <int dim>
  class ImmersedSurfaceQuadrature : public Quadrature<dim>
  {
  public:
    /**
     * Default constructor to initialize the quadrature with no quadrature
     * points.
     */
    ImmersedSurfaceQuadrature() = default;

    /**
     * Construct a quadrature formula from vectors of points, weights and
     * surface normals. The points, weights and normals should be with respect
     * to reference space, and the normals should be normalized.
     */
    ImmersedSurfaceQuadrature(const std::vector<Point<dim>>&     points,
                              const std::vector<double>&         weights,
                              const std::vector<Tensor<1, dim>>& normals);

    /**
     * Extend the given formula by an additional quadrature point.
     * The point, weight and normal should be with respect to reference space,
     * and the normal should be normalized.
     *
     * This function exists since immersed quadrature rules can be rather
     * complicated to construct. Often the construction is done by
     * partitioning the cell into regions and constructing points on each
     * region separately. This can make it cumbersome to create the quadrature
     * from the constructor since all quadrature points have to be known at
     * time of creation of the object.
     *
     * @note This function should only be used during construction of the
     * quadrature formula.
     */
    void
    push_back(const Point<dim>&     point,
              const double          weight,
              const Tensor<1, dim>& normal);

    /**
     * Return a reference to the <tt>i</tt>th surface normal.
     */
    const Tensor<1, dim>&
    normal_vector(const unsigned int i) const;

    /**
     * Return a reference to the whole %vector of normals.
     */
    const std::vector<Tensor<1, dim>>&
    get_normal_vectors() const;

  protected:
    /**
     * %Vector of surface normals at each quadrature point.
     */
    std::vector<Tensor<1, dim>> normals;
  };

} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE

#endif

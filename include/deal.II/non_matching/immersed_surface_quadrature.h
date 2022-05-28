// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2021 by the deal.II authors
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
   * This class defines a quadrature formula to integrate over the intersection
   * between an oriented surface, $\hat{S}$, and a cell or face. The word
   * "immersed" in the class name reflects that the surface may intersect the
   * cell/face in an arbitrary way.
   *
   * The spacedim template parameter of this class is the dimension that the
   * (spacedim-1)-dimensional surface is embedded in:
   * $\hat{S} \subset \mathbb{R}^{\text{spacedim}}$. The dim parameter
   * describes the dimension of the "object" that the surface intersects. That
   * is, dim = spacedim corresponds to the surface intersecting a cell and
   * dim = spacedim - 1 corresponds to the surface intersecting a face. The
   * quadrature formula is described by a set of quadrature points,
   * $\hat{x}_q \in \mathbb{R}^{\text{dim}}$, weights, $w_q$, and normalized
   * surface normals, $\hat{n}_q \in \mathbb{R}^{\text{spacedim}}$.
   *
   * Consider first the case dim = spacedim. We typically want to compute
   * integrals in real space. A surface, $S$, intersecting a cell, $K$, in
   * real space can be mapped onto a surface, $\hat{S}$, intersecting the unit
   * cell, $\hat{K}$. Thus an integral over $S\cap K$ in real space can be
   * transformed to an integral over $\hat{S} \cap \hat{K}$ according to
   * @f[
   * \int_{S\cap K} f dS =
   * \int_{S\cap K} f |d\bar{S}| =
   * \int_{\hat{S}\cap\hat{K}} f \circ F_{K} \det(J) |\left( J^{-1} \right
   * )^T d\hat{S}|,
   * @f]
   * where $F_K$ is the mapping from reference to real space and $J$ is its
   * Jacobian matrix. This transformation is possible since the continuous
   * surface elements are vectors:
   * $d\bar{S}, d\hat{S} \in \mathbb{R}^{spacedim}$, which are parallel to the
   * normals of $S$ and $\hat{S}$. That is, the normal is needed to do the
   * transformation. Thus, in addition to storing points and weights, this
   * quadrature stores also the normalized normal for each quadrature point.
   * This can be viewed as storing a discrete surface element,
   * @f[
   * \Delta \hat{S}_q \dealcoloneq w_q \hat{n}_q \approx d\hat{S}(\hat{x}_q),
   * @f]
   * for each quadrature point. The surface integral in real space would then be
   * approximated as
   * @f[
   * \int_{S\cap K} f dS \approx
   * \sum_{q} f \left(F_{K}(\hat{x}_{q}) \right) \det(J_q)
   * |\left( J_q^{-1} \right)^T \hat{n}_q| w_q.
   * @f]
   *
   * @image html immersed_surface_quadrature.svg
   *
   * When dim = spacedim - 1, this class represents a (spacedim-2)-dimensional
   * integral. That is, if spacedim = 3 we have a line integral immersed in a
   * face. Here, the transformation between the face, $F$, and reference face,
   * $\hat{F}$, reads
   * @f[
   * \int_{S\cap F} f dr =
   * \int_{S\cap F} f |d\bar{r}| =
   * \int_{\hat{S}\cap\hat{F}} f \circ F_{K} | J d\hat{r}|
   * \approx \sum_{q} f \left(F_{K}(\hat{x}_{q}) \right) |J_q \hat{t}_q| w_q,
   * @f]
   * where $\hat{t}_q = \hat{n}_q \times \hat{n}_F$ is the tangent to the curve
   * at $\hat{x}_q$ and $\hat{n}_F$ is the face normal. It would be possible to
   * compute the tangent by only knowing the normal to the curve in the face
   * plane (i.e. the dim-dimensional normal). However, when these quadratures
   * are used, the weak form typically involves the so-called conormal, which
   * can not be computed without knowing the surface normal in
   * $\mathbb{R}^{\text{spacedim}}$. The conormal is the unit vector parallel to
   * the projection of the face normal into the surface plane. This is
   * essentially the same thing as the normalized
   * @ref GlossBoundaryForm "boundary form".
   */
  template <int dim, int spacedim = dim>
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
    ImmersedSurfaceQuadrature(const std::vector<Point<dim>> &         points,
                              const std::vector<double> &             weights,
                              const std::vector<Tensor<1, spacedim>> &normals);

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
    push_back(const Point<dim> &         point,
              const double               weight,
              const Tensor<1, spacedim> &normal);

    /**
     * Return a reference to the <tt>i</tt>th surface normal.
     */
    const Tensor<1, spacedim> &
    normal_vector(const unsigned int i) const;

    /**
     * Return a reference to the whole %vector of normals.
     */
    const std::vector<Tensor<1, spacedim>> &
    get_normal_vectors() const;

  protected:
    /**
     * %Vector of surface normals at each quadrature point.
     */
    std::vector<Tensor<1, spacedim>> normals;
  };

} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE

#endif

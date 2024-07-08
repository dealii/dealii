// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mapping_c1_h
#define dealii_mapping_c1_h


#include <deal.II/base/config.h>

#include <deal.II/fe/mapping_q.h>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup mapping
 * @{
 */

/**
 * Mapping class that uses C1 (continuously differentiable) cubic mappings of
 * the boundary. This class is built atop of MappingQ by simply
 * determining the interpolation points for a cubic mapping of the boundary
 * differently: MappingQ chooses them such that they interpolate the boundary,
 * while this class chooses them such that the discretized boundary is
 * globally continuously differentiable.
 */
template <int dim, int spacedim = dim>
class MappingC1 : public MappingQ<dim, spacedim>
{
public:
  /**
   * Constructor. Pass the fixed degree @p 3 down to the base class, as a
   * cubic mapping suffices to generate a continuous mapping of the boundary.
   */
  MappingC1();

  /**
   * Return a pointer to a copy of the present object. The caller of this copy
   * then assumes ownership of it.
   */
  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const override;

  /**
   * For <tt>dim=2,3</tt>. Append the support points of all shape functions
   * located on bounding lines to the vector @p a. Points located on the
   * line but on vertices are not included.
   *
   * This function chooses the respective points not such that they are
   * interpolating the boundary (as does the base class), but rather such
   * that the resulting cubic mapping is a continuous one.
   */
  virtual void
  add_line_support_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    std::vector<Point<spacedim>> &a) const override;

  /**
   * For <tt>dim=3</tt>. Append the support points of all shape functions
   * located on bounding faces (quads in 3d) to the vector @p a. Points
   * located on the line but on vertices are not included.
   *
   * This function chooses the respective points not such that they are
   * interpolating the boundary (as does the base class), but rather such
   * that the resulting cubic mapping is a continuous one.
   */
  virtual void
  add_quad_support_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    std::vector<Point<spacedim>> &a) const override;
};

/** @} */

/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN

template <>
void
MappingC1<1>::add_line_support_points(const Triangulation<1>::cell_iterator &,
                                      std::vector<Point<1>> &) const;
template <>
void
MappingC1<2>::add_line_support_points(
  const Triangulation<2>::cell_iterator &cell,
  std::vector<Point<2>>                 &a) const;

template <>
void
MappingC1<1>::add_quad_support_points(const Triangulation<1>::cell_iterator &,
                                      std::vector<Point<1>> &) const;
template <>
void
MappingC1<2>::add_quad_support_points(const Triangulation<2>::cell_iterator &,
                                      std::vector<Point<2>> &) const;


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

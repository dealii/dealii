// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2018 by the deal.II authors
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

#ifndef dealii_mapping_c1_h
#define dealii_mapping_c1_h


#include <deal.II/base/config.h>

#include <deal.II/fe/mapping_q.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup mapping */
/*@{*/

/**
 * Mapping class that uses C1 (continuously differentiable) cubic mappings of
 * the boundary. This class is built atop of MappingQ by simply determining
 * the interpolation points for a cubic mapping of the boundary differently:
 * MappingQ chooses them such that they interpolate the boundary, while this
 * class chooses them such that the discretized boundary is globally
 * continuously differentiable.
 *
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

protected:
  /**
   * A class derived from MappingQGeneric that provides the generic mapping
   * with support points on boundary objects so that the corresponding Q3
   * mapping ends up being C1.
   */
  class MappingC1Generic : public MappingQGeneric<dim, spacedim>
  {
  public:
    /**
     * Constructor.
     */
    MappingC1Generic();

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
      const typename Triangulation<dim>::cell_iterator &cell,
      std::vector<Point<dim>> &                         a) const override;

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
      const typename Triangulation<dim>::cell_iterator &cell,
      std::vector<Point<dim>> &                         a) const override;
  };
};

/*@}*/

/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN

template <>
void
MappingC1<1>::MappingC1Generic::add_line_support_points(
  const Triangulation<1>::cell_iterator &,
  std::vector<Point<1>> &) const;
template <>
void
MappingC1<2>::MappingC1Generic::add_line_support_points(
  const Triangulation<2>::cell_iterator &cell,
  std::vector<Point<2>> &                a) const;

template <>
void
MappingC1<1>::MappingC1Generic::add_quad_support_points(
  const Triangulation<1>::cell_iterator &,
  std::vector<Point<1>> &) const;
template <>
void
MappingC1<2>::MappingC1Generic::add_quad_support_points(
  const Triangulation<2>::cell_iterator &,
  std::vector<Point<2>> &) const;


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

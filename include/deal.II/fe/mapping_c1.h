// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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

#ifndef __deal2__mapping_c1_h
#define __deal2__mapping_c1_h


#include <deal.II/base/config.h>
#include <deal.II/fe/mapping_q.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup mapping */
/*@{*/

/**
 * Mapping class that uses C1 (continuously differentiable) cubic
 * mappings of the boundary. This class is built atop of
 * MappingQ by simply determining the interpolation points for a
 * cubic mapping of the boundary differently: MappingQ chooses
 * them such that they interpolate the boundary, while this class
 * chooses them such that the discretized boundary is globally
 * continuously differentiable.
 *
 * To use this class, make sure that the
 * Boundary::@p get_normals_at_vertices function is implemented
 * for the user's boundary object.
 *
 * For more information about the <tt>spacedim</tt> template parameter
 * check the documentation of FiniteElement or the one of
 * Triangulation.
 *
 * @author Wolfgang Bangerth, 2001
 */
template<int dim, int spacedim=dim>
class MappingC1 : public MappingQ<dim,spacedim>
{
public:
  /**
   * Constructor. Pass the fixed degree @p 3 down to the base class, as a
   * cubic mapping suffices to generate a continuous mapping of the boundary.
   */
  MappingC1 ();

  /**
   * Return a pointer to a copy of the present object. The caller of this copy
   * then assumes ownership of it.
   */
  virtual
  Mapping<dim,spacedim> *clone () const;

protected:
  /**
   * For <tt>dim=2,3</tt>. Append the support points of all shape functions
   * located on bounding lines to the vector @p a. Points located on the line
   * but on vertices are not included.
   *
   * Needed by the <tt>compute_support_points_simple(laplace)</tt>
   * functions. For <tt>dim=1</tt> this function is empty.
   *
   * This function chooses the respective points not such that they are
   * interpolating the boundary (as does the base class), but rather such that
   * the resulting cubic mapping is a continuous one.
   */
  virtual void
  add_line_support_points (const typename Triangulation<dim>::cell_iterator &cell,
                           std::vector<Point<dim> > &a) const;

  /**
   * For <tt>dim=3</tt>. Append the support points of all shape functions
   * located on bounding faces (quads in 3d) to the vector @p a. Points
   * located on the line but on vertices are not included.
   *
   * Needed by the @p compute_support_points_laplace function. For
   * <tt>dim=1</tt> and 2 this function is empty.
   *
   * This function chooses the respective points not such that they are
   * interpolating the boundary (as does the base class), but rather such that
   * the resulting cubic mapping is a continuous one.
   */
  virtual void
  add_quad_support_points(const typename Triangulation<dim>::cell_iterator &cell,
                          std::vector<Point<dim> > &a) const;
};

/*@}*/

/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN

template <> void MappingC1<1>::add_line_support_points (
  const Triangulation<1>::cell_iterator &,
  std::vector<Point<1> > &) const;
template <> void MappingC1<2>::add_line_support_points (
  const Triangulation<2>::cell_iterator &cell,
  std::vector<Point<2> > &a) const;

template <> void MappingC1<1>::add_quad_support_points (
  const Triangulation<1>::cell_iterator &,
  std::vector<Point<1> > &) const;
template <> void MappingC1<2>::add_quad_support_points (
  const Triangulation<2>::cell_iterator &,
  std::vector<Point<2> > &) const;


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

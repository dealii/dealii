// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2015 by the deal.II authors
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

#ifndef dealii__mapping_q1_h
#define dealii__mapping_q1_h


#include <deal.II/base/derivative_form.h>
#include <deal.II/base/config.h>
#include <deal.II/base/table.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup mapping */
/*@{*/


/**
 * Mapping of the reference to cell to a general
 * quadrilateral/hexahedra by $d$-linear shape functions.
 *
 * This mapping implemented by this class maps the reference (unit) cell
 * to a general grid cell with
 * straight lines in $d$ dimensions. (Note, however, that in 3D the
 * <i>faces</i> of a general, trilinearly mapped cell may be curved, even if the
 * edges are not). This is the standard mapping used for polyhedral domains. It
 * is also the mapping used throughout deal.II for many functions that two
 * variants, one that allows to pass a mapping argument explicitly and one
 * that simply falls back to the MappingQ1 class declared here.
 *
 * The shape functions for this mapping are the same as for the finite element FE_Q
 * of order 1. Therefore, coupling these two yields an isoparametric element.
 *
 * @author Guido Kanschat, 2000, 2001; Ralf Hartmann, 2000, 2001, 2005
 */
template <int dim, int spacedim=dim>
class MappingQ1 : public MappingQGeneric<dim,spacedim>
{
public:
  /**
   * Default constructor.
   */
  MappingQ1 ();

  // for documentation, see the Mapping base class
  virtual
  Mapping<dim,spacedim> *clone () const;

  /**
   * @name Mapping points between reference and real cells
   * @{
   */

  // for documentation, see the Mapping base class
  virtual
  Point<spacedim>
  transform_unit_to_real_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                               const Point<dim>                                 &p) const;

  // for documentation, see the Mapping base class
  virtual
  Point<dim>
  transform_real_to_unit_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                               const Point<spacedim>                            &p) const;

  /**
   * @}
   */


  /**
   * @name Interface with FEValues
   * @{
   */

public:
  /**
   * Use the InternalData class of the base class without modification
   * and additions.
   */
  typedef typename MappingQGeneric<dim,spacedim>::InternalData InternalData;

protected:


  /**
   * @}
   */

protected:
  /* Trick to templatize transform_real_to_unit_cell<dim, dim+1> */
  template<int dim_>
  Point<dim_>
  transform_real_to_unit_cell_internal_codim1
  (const typename Triangulation<dim_,dim_+1>::cell_iterator &cell,
   const Point<dim_+1> &p,
   const Point<dim_>         &initial_p_unit,
   InternalData        &mdata) const;

  /**
   * Compute an initial guess to pass to the Newton method in
   * transform_real_to_unit_cell.  For the initial guess we proceed in the
   * following way:
   * <ul>
   * <li> find the least square dim-dimensional plane approximating the cell
   * vertices, i.e. we find and affine map A x_hat + b from the reference cell
   * to the real space.
   * <li> Solve the equation A x_hat + b = p for x_hat
   * <li> This x_hat is the initial solution used for the Newton Method.
   * </ul>
   * @note if dim<spacedim we first project p onto the plane. @note if dim==1
   * (for any spacedim) the initial guess is the exact solution and no Newton
   * iteration is needed.   Some details about how we compute the least square
   * plane. We look for a  spacedim x (dim + 1) matrix  X such that  X * M = Y
   * where M is a (dim+1) x n_vertices  matrix and Y a spacedim x n_vertices.
   * And: The i-th column of M is unit_vertex[i] and the last row all 1's. The
   * i-th column of Y is real_vertex[i].  If we split X=[A|b], the least
   * square approx is A x_hat+b  Classically  X = Y * (M^t (M M^t)^{-1})  Let
   * K = M^t * (M M^t)^{-1} = [KA Kb] this can be precomputed, and that is
   * exactly what we do.  Finally A = Y*KA  and  b = Y*Kb.
   */
  Point<dim>
  transform_real_to_unit_cell_initial_guess (const std::vector<Point<spacedim> > &vertex,
                                             const Point<spacedim>                            &p) const;


  /**
   * Transforms a point @p p on the unit cell to the point @p p_real on the
   * real cell @p cell and returns @p p_real.
   *
   * This function is called by @p transform_unit_to_real_cell and multiple
   * times (through the Newton iteration) by @p
   * transform_real_to_unit_cell_internal.
   *
   * Takes a reference to an @p InternalData that must already include the
   * shape values at point @p p and the mapping support points of the cell.
   *
   * This @p InternalData argument avoids multiple computations of the shape
   * values at point @p p and especially multiple computations of the mapping
   * support points.
   */
  Point<spacedim>
  transform_unit_to_real_cell_internal (const InternalData &mdata) const;

  /**
   * Transforms the point @p p on the real cell to the corresponding point on
   * the unit cell @p cell by a Newton iteration.
   *
   * Takes a reference to an @p InternalData that is assumed to be previously
   * created by the @p get_data function with @p UpdateFlags including @p
   * update_transformation_values and @p update_transformation_gradients and a
   * one point Quadrature that includes the given initial guess for the
   * transformation @p initial_p_unit.  Hence this function assumes that @p
   * mdata already includes the transformation shape values and gradients
   * computed at @p initial_p_unit.
   *
   * @p mdata will be changed by this function.
   */
  Point<dim>
  transform_real_to_unit_cell_internal (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                                        const Point<spacedim> &p,
                                        const Point<dim> &initial_p_unit,
                                        InternalData &mdata) const;

  /**
   * Computes the support points of the mapping. For @p MappingQ1 these are
   * the vertices. However, other classes may override this function. In
   * particular, the MappingQ1Eulerian class does exactly this by not
   * computing the support points from the geometry of the current cell but
   * instead evaluating an externally given displacement field in addition to
   * the geometry of the cell.
   */
  virtual void compute_mapping_support_points(
    const typename Triangulation<dim,spacedim>::cell_iterator &cell,
    std::vector<Point<spacedim> > &a) const;
};


#ifndef DOXYGEN
// explicit specializations

template<>
Point<2>
MappingQ1<2,3>::
transform_real_to_unit_cell_internal
(const Triangulation<2,3>::cell_iterator &cell,
 const Point<3> &p,
 const Point<2> &initial_p_unit,
 InternalData    &mdata) const;

template<>
Point<1>
MappingQ1<1,2>::
transform_real_to_unit_cell_internal
(const Triangulation<1,2>::cell_iterator &cell,
 const Point<2> &p,
 const Point<1> &initial_p_unit,
 InternalData    &mdata) const;

template<>
Point<1>
MappingQ1<1,3>::
transform_real_to_unit_cell_internal
(const Triangulation<1,3>::cell_iterator &cell,
 const Point<3> &p,
 const Point<1> &initial_p_unit,
 InternalData    &mdata) const;

#endif

/**
 * In order to avoid creation of static MappingQ1 objects at several places in
 * the library (in particular in backward compatibility functions), we define
 * a static MappingQ1 objects once and for all places where it is needed.
 */
template <int dim, int spacedim=dim>
struct StaticMappingQ1
{
  static MappingQ1<dim, spacedim> mapping;
};


/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif

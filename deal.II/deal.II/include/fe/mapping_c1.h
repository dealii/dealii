//----------------------------  mapping_c1.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mapping_c1.h  ---------------------------
#ifndef __deal2__mapping_c1_h
#define __deal2__mapping_c1_h


#include <base/config.h>
#include <fe/mapping_q.h>


/*!@addtogroup fe */
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
 * for the users boundary object.
 *
 * @author Wolfgang Bangerth, 2001
 */
template <int dim>
class MappingC1 : public MappingQ<dim>
{
  public:
				     /**
				      * Constructor. Pass the fixed
				      * degree @p 3 down to the base
				      * class, as a cubic mapping
				      * suffices to generate a
				      * continuous mapping of the
				      * boundary.
				      */
    MappingC1 ();

  protected:
				     /**
				      * For <tt>dim=2,3</tt>. Append the
				      * support points of all shape
				      * functions located on bounding
				      * lines to the vector
				      * @p a. Points located on the
				      * line but on vertices are not
				      * included.
				      *
				      * Needed by the
				      * <tt>compute_support_points_simple(laplace)</tt>
				      * functions. For <tt>dim=1</tt> this
				      * function is empty.
				      *
				      * This function chooses the
				      * respective points not such
				      * that they are interpolating
				      * the boundary (as does the base
				      * class), but rather such that
				      * the resulting cubic mapping is
				      * a continuous one.
				      */
    virtual void
    add_line_support_points (const typename Triangulation<dim>::cell_iterator &cell,
			     std::vector<Point<dim> > &a) const;

				     /**
				      * For <tt>dim=3</tt>. Append the
				      * support points of all shape
				      * functions located on bounding
				      * faces (quads in 3d) to the
				      * vector @p a. Points located
				      * on the line but on vertices
				      * are not included.
				      *
				      * Needed by the
				      * @p compute_support_points_laplace
				      * function. For <tt>dim=1</tt> and 2
				      * this function is empty.
				      *
				      * This function chooses the
				      * respective points not such
				      * that they are interpolating
				      * the boundary (as does the base
				      * class), but rather such that
				      * the resulting cubic mapping is
				      * a continuous one.
				      */
    virtual void
    add_quad_support_points(const typename Triangulation<dim>::cell_iterator &cell,
			    std::vector<Point<dim> > &a) const;
};

/*@}*/

/* -------------- declaration of explicit specializations ------------- */

/// @if NoDoc

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


/// @endif

#endif

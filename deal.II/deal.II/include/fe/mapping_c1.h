//----------------------------  mapping_c1.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mapping_c1.h  ---------------------------
#ifndef __deal2__mapping_c1_h
#define __deal2__mapping_c1_h


#include <fe/mapping_q.h>


/**
 * Mapping class that uses C1 (continuously differentiable) cubic
 * mappings of the boundary. This class is built atop of
 * @ref{MappingQ} by simply determining the interpolation points for a
 * cubic mapping of the boundary differently: @ref{MappingQ} chooses
 * them such that they interpolate the boundary, while this class
 * chooses them such that the discretized boundary is globally
 * continuous.
 *  
 * @author Wolfgang Bangerth, 2001
 */
template <int dim>
class MappingC1 : public MappingQ<dim>
{
  public:
				     /**
				      * Constructor. Pass the fixed
				      * degree @p{3} down to the base
				      * class, as a cubic mapping
				      * suffices to generate a
				      * continuous mapping of the
				      * boundary.
				      */
    MappingC1 ();

  protected:
				     /**
				      * For @p{dim=2,3}. Append the
				      * support points of all shape
				      * functions located on bounding
				      * lines to the vector
				      * @p{a}. Points located on the
				      * line but on vertices are not
				      * included.
				      *
				      * Needed by the
				      * @p{compute_support_points_simple(laplace)}
				      * functions. For @p{dim=1} this
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
			     typename std::vector<Point<dim> > &a) const;

				     /**
				      * For @p{dim=3}. Append the
				      * support points of all shape
				      * functions located on bounding
				      * faces (quads in 3d) to the
				      * vector @p{a}. Points located
				      * on the line but on vertices
				      * are not included.
				      *
				      * Needed by the
				      * @p{compute_support_points_laplace}
				      * function. For @p{dim=1} and 2
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
			    typename std::vector<Point<dim> > &a) const;
};


#endif

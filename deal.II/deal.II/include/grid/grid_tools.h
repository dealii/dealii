//----------------------------  grid_tools.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  grid_tools.h  ---------------------------
#ifndef __deal2__grid_tools_H
#define __deal2__grid_tools_H


#include <base/config.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>




/**
 * This class is a collection of algorithms working on
 * triangulations. See the descriptions of the individual functions
 * for more information.
 *
 * @author Wolfgang Bangerth, 2001
 */
class GridTools
{
  public:
				     /**
				      * Return the diameter of a
				      * triangulation. The diameter is
				      * computed using only the
				      * vertices, i.e. if the diameter
				      * should be larger than the
				      * maximal distance between
				      * boundary vertices due to a
				      * higher order mapping, then
				      * this function will not catch
				      * this.
				      */
    template <int dim>
    static
    double diameter (const Triangulation<dim> &tria);

				     /**
				      * Same function, but for 1d.
				      */
    static
    double diameter (const Triangulation<1> &tria);

				     /**
				      * Transform the vertices of the
				      * given triangulation by
				      * applying the predicate to all
				      * its vertices. Since the
				      * internal consistency of a
				      * triangulation can only be
				      * guaranteed if the
				      * transformation is applied to
				      * the vertices of only one level
				      * of a hierarchically refined
				      * cells, this function may only
				      * be used on coarse grids,
				      * i.e. before any refinement of
				      * it has taken place.
				      *
				      * The predicate given as
				      * argument is used to transform
				      * each vertex. Its respective
				      * type has to offer a
				      * function-like syntax, i.e. the
				      * predicate is either an object
				      * of a type that has an
				      * @p{operator()}, or it is a
				      * pointer to the function. In
				      * either case, argument and
				      * return value have to be of
				      * type @p{Point<dim>}.
				      */
    template <int dim, typename Predicate>
    static
    void transform (const Predicate    &predicate,
		    Triangulation<dim> &triangulation);

				     /**
				      * Shift each vertex of the
				      * triangulation by the given
				      * shift vector. This function
				      * uses the @ref{transform}
				      * function above, so the
				      * requirements on the
				      * triangulation stated there
				      * hold for this function as
				      * well.
				      */
    template <int dim>
    static
    void shift (const Point<dim>   &shift_vector,
		Triangulation<dim> &triangulation);


				     /**
				      * Rotate all vertices of the
				      * given two-dimensional
				      * triangulation in
				      * counter-clockwise sense around
				      * the origin of the coordinate
				      * system by the given angle
				      * (given in radians, rather than
				      * degrees). This function uses
				      * the @ref{transform} function
				      * above, so the requirements on
				      * the triangulation stated there
				      * hold for this function as
				      * well.
				      */
    static
    void rotate (const double      angle,
		 Triangulation<2> &triangulation);

				     /**
				      * Scale the entire triangulation
				      * by the given factor. To
				      * preserve the orientation of
				      * the triangulation, the factor
				      * must be positive.
				      *
				      * This function uses the
				      * @ref{transform} function
				      * above, so the requirements on
				      * the triangulation stated there
				      * hold for this function as
				      * well.
				      */
    template <int dim>
    static
    void scale (const double        scaling_factor,
		Triangulation<dim> &triangulation);
    
				     /**
				      * Exception
				      */
    DeclException0 (ExcTriangulationHasBeenRefined);

				     /**
				      * Exception
				      */
    DeclException1 (ExcScalingFactorNotPositive,
		    double,
		    << "The scaling factor must be positive, but is " << arg1);
};



/* ----------------- Template function --------------- */

template <int dim, typename Predicate>
void GridTools::transform (const Predicate    &predicate,
			   Triangulation<dim> &triangulation)
{
  Assert (triangulation.n_levels() == 1,
	  ExcTriangulationHasBeenRefined());
  
  std::vector<bool> treated_vertices (triangulation.n_vertices(),
				      false);

				   // loop over all active cells, and
				   // transform those vertices that
				   // have not yet been touched. note
				   // that we get to all vertices in
				   // the triangulation by only
				   // visiting the active cells.
  typename Triangulation<dim>::active_cell_iterator
    cell = triangulation.begin_active (),
    endc = triangulation.end ();
  for (; cell!=endc; ++cell)
    for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
      if (treated_vertices[cell->vertex_index(v)] == false)
	{
					   // transform this vertex
	  cell->vertex(v) = predicate(cell->vertex(v));
					   // and mark it as treated
	  treated_vertices[cell->vertex_index(v)] = true;
	};
};




/*----------------------------   grid_tools.h     ---------------------------*/
/* end of #ifndef __deal2__grid_tools_H */
#endif
/*----------------------------   grid_tools.h     ---------------------------*/

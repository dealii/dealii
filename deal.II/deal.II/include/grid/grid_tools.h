//----------------------------  grid_tools.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004 by the deal.II authors
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
 * This class is a collection of algorithms working on triangulations,
 * such as shifting or rotating triangulations, but also finding a
 * cell that contains a given point. See the descriptions of the
 * individual functions for more information.
 *
 * @author Wolfgang Bangerth, 2001, 2003, 2004
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
                                      * Find and return an iterator to
                                      * the active cell that surrounds
                                      * a given point @p{ref}. The
                                      * type of the first parameter
                                      * may be either
                                      * @ref{Triangulation},
                                      * @ref{DoFHandler}, or
                                      * @ref{MGDoFHandler}, i.e. we
                                      * can find the cell around a
                                      * point for iterators into each
                                      * of these classes.
                                      *
                                      * The algorithm used in this
                                      * function proceeds by first
                                      * looking for the surrounding
                                      * cell on the coarse grid, and
                                      * then recursively checking its
                                      * sibling cells. The complexity
                                      * is thus @p{O(M+log N)} where
                                      * @p{M} is the number of coarse
                                      * grid cells, and @p{N} the
                                      * total number of cells.
                                      *
                                      * There are cases where this
                                      * function will not found a
                                      * given point in space
                                      * dimensions higher than one,
                                      * even though it is inside the
                                      * domain being discretized, or
                                      * will find a point that is
                                      * actually outside the
                                      * domain. The reason for this is
                                      * that we use piecewise d-linear
                                      * mappings of the unit cell to
                                      * real cells. Thus, if a point
                                      * is close to a convex boundary
                                      * or on it, it may not be inside
                                      * any of the cells since they
                                      * have straight boundaries that
                                      * lie entirely inside the
                                      * domain.
                                      *
                                      * Another case for this is that
                                      * a point may not be found even
                                      * though it is actually in one
                                      * of the cells. This may happen,
                                      * if the point is not in one of
                                      * the coarse grid cells, even
                                      * though it is in one of the
                                      * cells on finer levels of the
                                      * triangulation. Note that this
                                      * of course implies that mother
                                      * and child cells do not exactly
                                      * overlap, a case that is
                                      * frequent along curved
                                      * boundaries. In this latter
                                      * case, a different algorithm
                                      * may be used instead that uses
                                      * a linear search over all
                                      * active cells, rather than
                                      * first searchin for a coarse
                                      * grid cell. Note, however, that
                                      * such an algorithm has a
                                      * significantly higher numerical
                                      * cost than the logarithmic
                                      * algorithm used here.
                                      *
                                      * Lastly, if a point lies on the
                                      * boundary of two or more cells,
                                      * then the algorithm may return
                                      * with any of these cells. While
                                      * this is in general not really
                                      * problem, if may be a nuisance
                                      * if the point lies at the
                                      * boundary of cells with
                                      * different refinement levels
                                      * and one would rather like to
                                      * evaluate a solution on the
                                      * cell with more refinement. For
                                      * this, more sophisticated
                                      * algorithms would be necessary,
                                      * though.
                                      */
    template <int dim, typename Container>
    static
    typename Container::active_cell_iterator
    find_active_cell_around_point (const Container  &container,
                                   const Point<dim> &p);

                                     /**
                                      * Use the METIS partitioner to generate
                                      * a partitioning of the active cells
                                      * making up the entire domain. After
                                      * calling this function, the subdomain
                                      * ids of all active cells will have
                                      * values between zero and
                                      * @p{n_partitions-1}. You can access the
                                      * subdomain id of a cell by using
                                      * @p{cell->subdomain_id()}.
                                      *
                                      * This function will generate an error
                                      * if METIS is not installed unless
                                      * @p{n_partitions} is one. I.e., you can
                                      * write a program so that it runs in the
                                      * single-processor single-partition case
                                      * without METIS installed, and only
                                      * requires METIS when multiple
                                      * partitions are required.
                                      */
    template <int dim>
    static
    void partition_triangulation (Triangulation<dim> &triangulation,
                                  const unsigned int  n_partitions);
                                     /**
                                      * Exception
                                      */
    DeclException1 (ExcInvalidNumberOfPartitions,
                    int,
                    << "The number of partitions you gave is " << arg1
                    << ", but must be greater than zero.");
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
				     /**
				      * Exception
				      */
    template <int N>
    DeclException1 (ExcPointNotFoundInCoarseGrid,
		    Point<N>,
		    << "The point <" << arg1
                    << "> could not be found inside any of the "
                    << "coarse grid cells.");
				     /**
				      * Exception
				      */
    template <int N>
    DeclException1 (ExcPointNotFound,
		    Point<N>,
		    << "The point <" << arg1
                    << "> could not be found inside any of the "
                    << "subcells of a coarse grid cell.");
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
}




/*----------------------------   grid_tools.h     ---------------------------*/
/* end of #ifndef __deal2__grid_tools_H */
#endif
/*----------------------------   grid_tools.h     ---------------------------*/

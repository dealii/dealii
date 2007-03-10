//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__grid_tools_H
#define __deal2__grid_tools_H


#include <base/config.h>
#include <fe/mapping.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>

#include <list>

DEAL_II_NAMESPACE_OPEN


template <int> class DoFHandler;
template <int> class Mapping;
namespace hp
{
  template <int> class DoFHandler;
  template <int> class MappingCollection;
}


/**
 * This class is a collection of algorithms working on triangulations,
 * such as shifting or rotating triangulations, but also finding a
 * cell that contains a given point. See the descriptions of the
 * individual functions for more information.
 *
 * @ingroup grid
 * @author Wolfgang Bangerth, 2001, 2003, 2004, Ralf Hartmann, 2005
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
				      * Given a list of vertices (typically
				      * obtained using
				      * Triangulation::get_vertices) as the
				      * first, and a list of vertex indices
				      * that characterize a single cell as the
				      * second argument, return the measure
				      * (area, volume) of this cell. If this
				      * is a real cell, then you can get the
				      * same result using
				      * <code>cell-@>measure()</code>, but
				      * this function also works for cells
				      * that do not exist except that you make
				      * it up by naming its vertices from the
				      * list.
				      */
    template <int dim>
    static
    double cell_measure (const std::vector<Point<dim> > &all_vertices,
			 const unsigned int (&vertex_indices)[GeometryInfo<dim>::vertices_per_cell]);

				     /**
				      * Remove vertices that are not
				      * referenced by any of the
				      * cells. This function is called
				      * by all <tt>GridIn::read_*</tt>
				      * functions to eliminate
				      * vertices that are listed in
				      * the input files but are not
				      * used by the cells in the input
				      * file. While these vertices
				      * should not be in the input
				      * from the beginning, they
				      * sometimes are, most often when
				      * some cells have been removed
				      * by hand without wanting to
				      * update the vertex lists, as
				      * they might be lengthy.
				      *
				      * This function is called by all
				      * <tt>GridIn::read_*</tt>
				      * functions as the triangulation
				      * class requires them to be
				      * called with used vertices
				      * only. This is so, since the
				      * vertices are copied verbatim
				      * by that class, so we have to
				      * eliminate unused vertices
				      * beforehand.
				      */
    template <int dim>
    static
    void delete_unused_vertices (std::vector<Point<dim> >    &vertices,
				 std::vector<CellData<dim> > &cells,
				 SubCellData                 &subcelldata);
    
				     /**
				      * Remove vertices that are duplicated,
				      * due to the input of a structured grid,
				      * for example. If these vertices are not
				      * removed, the faces bounded by these
				      * vertices become part of the boundary,
				      * even if they are in the interior of
				      * the mesh.
				      *
				      * This function is called by some
				      * <tt>GridIn::read_*</tt> functions. Only
				      * the vertices with indices in @p
				      * considered_vertices are tested for
				      * equality. This speeds up the algorithm,
				      * which is quadratic and thus quite slow
				      * to begin with. However, if you wish to
				      * consider all vertices, simply pass an
				      * empty vector.
				      *
				      * Two vertices are considered equal if
				      * their difference in each coordinate
				      * direction is less than @p tol.
				      */
    template <int dim>
    static
    void delete_duplicated_vertices (std::vector<Point<dim> >    &all_vertices,
				     std::vector<CellData<dim> > &cells,
				     SubCellData                 &subcelldata,
				     std::vector<unsigned int>   &considered_vertices,
				     const double                 tol=1e-12);
    
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
				      * <tt>operator()</tt>, or it is a
				      * pointer to the function. In
				      * either case, argument and
				      * return value have to be of
				      * type <tt>Point<dim></tt>.
				      */
    template <int dim, typename Predicate>
    static
    void transform (const Predicate    &predicate,
		    Triangulation<dim> &triangulation);

				     /**
				      * Shift each vertex of the
				      * triangulation by the given
				      * shift vector. This function
				      * uses the transform()
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
				      * the transform() function
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
				      * transform() function
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
                                      * Find and return the number of
                                      * the used vertex in a given
                                      * Container that is located closest
                                      * to a given point @p p. The
                                      * type of the first parameter
                                      * may be either Triangulation,
                                      * DoFHandler, hp::DoFHandler, or
                                      * MGDoFHandler.
                                      *
                                      * @author Ralf B. Schulz, 2006
                                      */
    template <int dim, template <int> class Container>
    static
    unsigned int
    find_closest_vertex (const Container<dim> &container,
                         const Point<dim>     &p);

                                     /**
                                      * Find and return a vector of
                                      * iterators to active cells that
                                      * surround a given vertex @p vertex.
                                      * The type of the first parameter
                                      * may be either Triangulation,
                                      * DoFHandler, hp::DoFHandler, or
                                      * MGDoFHandler.
                                      *
                                      * For locally refined grids, the
                                      * vertex itself might not be a vertex
                                      * of all adjacent cells, but will
                                      * always be located on a face or an
                                      * edge of the adjacent cells returned.
                                      *
                                      * @author Ralf B. Schulz,
                                      * Wolfgang Bangerth, 2006
                                      */
   template<int dim, template <int> class Container>
   static
   std::vector<typename Container<dim>::active_cell_iterator>
   find_cells_adjacent_to_vertex (const Container<dim> &container,
				  const unsigned int    vertex);


                                     /**
                                      * Find and return an iterator to
                                      * the active cell that surrounds
                                      * a given point @p ref. The
                                      * type of the first parameter
                                      * may be either
                                      * Triangulation,
                                      * DoFHandler, or
                                      * MGDoFHandler, i.e. we
                                      * can find the cell around a
                                      * point for iterators into each
                                      * of these classes.
                                      *
                                      * This is solely a wrapper function
                                      * for the @p interpolate function
                                      * given below,
                                      * providing backward compatibility.
                                      * A Q1 mapping is used for the
                                      * boundary, and the iterator to
                                      * the cell in which the point
                                      * resides is returned.
                                      *
                                      * It is recommended to use the
                                      * other version of this function,
                                      * as it simultaneously delivers the
                                      * local coordinate of the given point
                                      * without additional computational cost.
                                      */
    template <int dim, typename Container>
    static
    typename Container::active_cell_iterator
    find_active_cell_around_point (const Container  &container,
                                   const Point<dim> &p);

                                     /**
                                      * Find and return an iterator to
                                      * the active cell that surrounds
                                      * a given point @p p. The
                                      * type of the first parameter
                                      * may be either
                                      * Triangulation,
                                      * DoFHandler, hp::DoFHandler, or
                                      * MGDoFHandler, i.e., we
                                      * can find the cell around a
                                      * point for iterators into each
                                      * of these classes.
                                      *
                                      * The algorithm used in this
                                      * function proceeds by first
                                      * looking for vertex located
                                      * closest to the given point, see
                                      * find_closest_vertex(). Secondly,
                                      * all adjacent cells to this point
                                      * are found in the mesh, see
                                      * find_cells_adjacent_to_vertex().
                                      * Lastly, for each of these cells,
                                      * it is tested whether the point is
                                      * inside. This check is performed
                                      * using arbitrary boundary mappings.
                                      * Still, it is possible that due
                                      * to roundoff errors, the point
                                      * cannot be located exactly inside
                                      * the unit cell. In this case,
                                      * even points at a very small
                                      * distance outside the unit cell
                                      * are allowed. 
                                      *
                                      * If a point lies on the
                                      * boundary of two or more cells,
                                      * then the algorithm tries to identify
                                      * the cell that is of highest
                                      * refinement level.
                                      *
                                      * The function returns an
                                      * iterator to the cell, as well
                                      * as the local position of the
                                      * point inside the unit
                                      * cell. This local position
                                      * might be located slightly
                                      * outside an actual unit cell,
                                      * due to numerical roundoff.
                                      * Therefore, the point returned
                                      * by this function should 
                                      * be projected onto the unit cell,
                                      * using GeometryInfo::project_to_unit_cell.
                                      * This is not automatically performed
                                      * by the algorithm.
                                      */
    template <int dim, template<int> class Container>
    static
    std::pair<typename Container<dim>::active_cell_iterator, Point<dim> >
    find_active_cell_around_point (const Mapping<dim>   &mapping,
                                   const Container<dim> &container,
                                   const Point<dim>     &p);

				     /**
				      * A version of the previous function
				      * where we use that mapping on a given
				      * cell that corresponds to the active
				      * finite element index of that
				      * cell. This is obviously only useful
				      * for hp problems, since the active
				      * finite element index for all other DoF
				      * handlers is always zero.
				      */
    template <int dim>
    static
    std::pair<typename hp::DoFHandler<dim>::active_cell_iterator, Point<dim> >
    find_active_cell_around_point (const hp::MappingCollection<dim>   &mapping,
                                   const hp::DoFHandler<dim> &container,
                                   const Point<dim>     &p);

				     /**
				      * Return a list of all descendents of
				      * the given cell that are active. For
				      * example, if the current cell is once
				      * refined but none of its children are
				      * any further refined, then the returned
				      * list will contain all its children.
				      *
				      * If the current cell is already active,
				      * then the returned list is empty
				      * (because the cell has no children that
				      * may be active).
				      *
				      * Since in C++ the type of the Container
				      * template argument (which can be
				      * Triangulation, DoFHandler,
				      * MGDoFHandler, or hp::DoFHandler) can
				      * not be deduced from a function call,
				      * you will have to specify it after the
				      * function name, as for example in
				      * <code>GridTools::get_active_child_cells@<DoFHandler@<dim@>
				      * @> (cell)</code>.
				      */
    template <class Container>
    static
    std::vector<typename Container::active_cell_iterator>
    get_active_child_cells (const typename Container::cell_iterator &cell);
    
    
                                     /**
                                      * Use the METIS partitioner to generate
                                      * a partitioning of the active cells
                                      * making up the entire domain. After
                                      * calling this function, the subdomain
                                      * ids of all active cells will have
                                      * values between zero and
                                      * @p n_partitions-1. You can access the
                                      * subdomain id of a cell by using
                                      * <tt>cell->subdomain_id()</tt>.
                                      *
                                      * This function will generate an error
                                      * if METIS is not installed unless
                                      * @p n_partitions is one. I.e., you can
                                      * write a program so that it runs in the
                                      * single-processor single-partition case
                                      * without METIS installed, and only
                                      * requires METIS when multiple
                                      * partitions are required.
                                      */
    template <int dim>
    static
    void partition_triangulation (const unsigned int  n_partitions,
                                  Triangulation<dim> &triangulation);
    
                                     /**
                                      * For each active cell, return in the
                                      * output array to which subdomain (as
                                      * given by the <tt>cell->subdomain_id()</tt>
                                      * function) it belongs. The output array
                                      * is supposed to have the right size
                                      * already when calling this function.
				      *
				      * This function returns the association
				      * of each cell with one subdomain. If
				      * you are looking for the association of
				      * each @em DoF with a subdomain, use the
				      * <tt>DoFTools::get_subdomain_association</tt>
				      * function.
                                      */
    template <int dim>
    static void
    get_subdomain_association (const Triangulation<dim>  &triangulation,
                               std::vector<unsigned int> &subdomain);

                                     /**
                                      * Count how many cells are uniquely
                                      * associated with the given @p subdomain
                                      * index.
                                      *
                                      * This function will generate an
                                      * exception if there are no cells with
                                      * the given @p subdomain index.
				      *
				      * This function returns the number of
				      * cells associated with one
				      * subdomain. If you are looking for the
				      * association of @em DoFs with this
				      * subdomain, use the
				      * <tt>DoFTools::count_dofs_with_subdomain_association</tt>
				      * function.
                                      */
    template <int dim>
    static unsigned int
    count_cells_with_subdomain_association (const Triangulation<dim> &triangulation,
                                            const unsigned int        subdomain);

                                     /**
                                      * Given two mesh containers
                                      * (i.e. objects of type
                                      * Triangulation, DoFHandler,
                                      * hp::DoFHandler, or
                                      * MGDoFHandler) that are based
                                      * on the same coarse mesh, this
                                      * function figures out a set of
                                      * cells that are matched between
                                      * the two meshes and where at
                                      * most one of the meshes is more
                                      * refined on this cell. In other
                                      * words, it finds the smallest
                                      * cells that are common to both
                                      * meshes, and that together
                                      * completely cover the domain.
                                      *
                                      * This function is useful, for
                                      * example, in time-dependent or
                                      * nonlinear application, where
                                      * one has to integrate a
                                      * solution defined on one mesh
                                      * (e.g., the one from the
                                      * previous time step or
                                      * nonlinear iteration) against
                                      * the shape functions of another
                                      * mesh (the next time step, the
                                      * next nonlinear iteration). If,
                                      * for example, the new mesh is
                                      * finer, then one has to obtain
                                      * the solution on the coarse
                                      * mesh (mesh_1) and interpolate
                                      * it to the children of the
                                      * corresponding cell of
                                      * mesh_2. Conversely, if the new
                                      * mesh is coarser, one has to
                                      * express the coarse cell shape
                                      * function by a linear
                                      * combination of fine cell shape
                                      * functions. In either case, one
                                      * needs to loop over the finest
                                      * cells that are common to both
                                      * triangulations. This function
                                      * returns a list of pairs of
                                      * matching iterators to cells in
                                      * the two meshes that can be
                                      * used to this end.
                                      *
                                      * Note that the list of these
                                      * iterators is not necessarily
                                      * order, and does also not
                                      * necessarily coincide with the
                                      * order in which cells are
                                      * traversed in one, or both, of
                                      * the meshes given as arguments.
                                      */
    template <typename Container>
    static
    std::list<std::pair<typename Container::cell_iterator,
                        typename Container::cell_iterator> >
    get_finest_common_cells (const Container &mesh_1,
                             const Container &mesh_2);

                                     /**
                                      * Return true if the two
                                      * triangulations are based on
                                      * the same coarse mesh. This is
                                      * determined by checking whether
                                      * they have the same number of
                                      * cells on the coarsest level,
                                      * and then checking that they
                                      * have the same vertices.
                                      *
                                      * The two meshes may have
                                      * different refinement histories
                                      * beyond the coarse mesh.
                                      */
    template <int dim>
    static
    bool
    have_same_coarse_mesh (const Triangulation<dim> &mesh_1,
                           const Triangulation<dim> &mesh_2);

                                     /**
                                      * The same function as above,
                                      * but working on arguments of
                                      * type DoFHandler,
                                      * hp::DoFHandler, or
                                      * MGDoFHandler. This function is
                                      * provided to allow calling
                                      * have_same_coarse_mesh for all
                                      * types of containers
                                      * representing triangulations or
                                      * the classes built on
                                      * triangulations.
                                      */
    template <typename Container>
    static
    bool
    have_same_coarse_mesh (const Container &mesh_1,
                           const Container &mesh_2);

				     /**
				      * Return the diamater of the smallest
				      * active cell of a triangulation. See
				      * @ref step_24 "step-24" for an example
				      * of use of this function.
				      */
    template <int dim>
    static
    double
    minimal_cell_diameter (const Triangulation<dim> &triangulation);

				     /**
				      * Return the diamater of the largest
				      * active cell of a triangulation.
				      */
    template <int dim>
    static
    double
    maximal_cell_diameter (const Triangulation<dim> &triangulation);

				     /**
				      * Given the two triangulations
				      * specified as the first two
				      * arguments, create the
				      * triangulation that contains
				      * the finest cells of both
				      * triangulation and store it in
				      * the third parameter. Previous
				      * content of @p result will be
				      * deleted.
				      *
				      * The function assumes that the
				      * two input triangulations are
				      * derived from the same coarse
				      * grid.
				      */
    template <int dim>
    static
    void
    create_union_triangulation (const Triangulation<dim> &triangulation_1,
				const Triangulation<dim> &triangulation_2,
				Triangulation<dim>       &result);
    
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
    DeclException1 (ExcNonExistentSubdomain,
                    int,
                    << "The subdomain id " << arg1
                    << " has no cells associated with it.");
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

    DeclException1 (ExcVertexNotUsed,
		    unsigned int,
		    << "The given vertex " << arg1
		    << " is not used in the given triangulation");
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



template <class DH>
std::vector<typename DH::active_cell_iterator>
GridTools::get_active_child_cells (const typename DH::cell_iterator& cell)
{
  std::vector<typename DH::active_cell_iterator> child_cells;
  
  if (cell->has_children())
    {
      for (unsigned int child=0;
	   child<GeometryInfo<DH::dimension>::children_per_cell; ++child)
	if (cell->child (child)->has_children())
	  {
	    const std::vector<typename DH::active_cell_iterator>
	      children = get_active_child_cells<DH> (cell->child(child));
	    child_cells.insert (child_cells.end(),
				children.begin(), children.end());
	  }
	else
	  child_cells.push_back (cell->child(child));
    }

  return child_cells;
}






DEAL_II_NAMESPACE_CLOSE

/*----------------------------   grid_tools.h     ---------------------------*/
/* end of #ifndef __deal2__grid_tools_H */
#endif
/*----------------------------   grid_tools.h     ---------------------------*/

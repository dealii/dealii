//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010 by the deal.II authors
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


template <int, int> class DoFHandler;
template <int, int> class Mapping;
namespace hp
{
  template <int, int> class DoFHandler;
  template <int, int> class MappingCollection;
}

class SparsityPattern;


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
    template <int dim, int spacedim>
    static
    double diameter (const Triangulation<dim, spacedim> &tria);

				     /**
				      * Same function, but for 1d.
				      */
    static
    double diameter (const Triangulation<1> &tria);
				     /**
				      * Same function, but for 1d, 2sd.
				      */
    static
    double diameter (const Triangulation<1,2> &tria);

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
				      *
				      * Not implemented for the
				      * codimension one case.
				      */
    template <int dim, int spacedim>
    static
    void delete_unused_vertices (std::vector<Point<spacedim> >    &vertices,
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
    template <int dim, int spacedim>
    static
    void delete_duplicated_vertices (std::vector<Point<spacedim> >    &all_vertices,
				     std::vector<CellData<dim> > &cells,
				     SubCellData                 &subcelldata,
				     std::vector<unsigned int>   &considered_vertices,
				     const double                 tol=1e-12);

				     /**
				      * Transform the vertices of the given
				      * triangulation by applying the
				      * predicate to all its vertices. Since
				      * the internal consistency of a
				      * triangulation can only be guaranteed
				      * if the transformation is applied to
				      * the vertices of only one level of a
				      * hierarchically refined cells, this
				      * function may only be used if all cells
				      * of the triangulation are on the same
				      * refinement level.
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
    template <int dim, typename Predicate, int spacedim>
    static
    void transform (const Predicate    &predicate,
		    Triangulation<dim,spacedim> &triangulation);

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
    template <int dim, int spacedim>
    static
    void shift (const Point<spacedim>   &shift_vector,
		Triangulation<dim,spacedim> &triangulation);


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
    template <int dim, int spacedim>
    static
    void scale (const double        scaling_factor,
		Triangulation<dim, spacedim> &triangulation);

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
    template <int dim, template <int, int> class Container, int spacedim>
    static
    unsigned int
    find_closest_vertex (const Container<dim, spacedim> &container,
                         const Point<spacedim>     &p);

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
   template<int dim, template <int, int> class Container, int spacedim>
   static
   std::vector<typename Container<dim,spacedim>::active_cell_iterator>
   find_cells_adjacent_to_vertex (const Container<dim,spacedim> &container,
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
    template <int dim, template <int,int> class Container, int spacedim>
    static
    typename Container<dim,spacedim>::active_cell_iterator
    find_active_cell_around_point (const Container<dim,spacedim>  &container,
                                   const Point<spacedim> &p);

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
    template <int dim, template<int, int> class Container, int spacedim>
    static
    std::pair<typename Container<dim,spacedim>::active_cell_iterator, Point<spacedim> >
    find_active_cell_around_point (const Mapping<dim,spacedim>   &mapping,
                                   const Container<dim,spacedim> &container,
                                   const Point<spacedim>     &p);

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
    template <int dim, int spacedim>
    static
    std::pair<typename hp::DoFHandler<dim,spacedim>::active_cell_iterator, Point<spacedim> >
    find_active_cell_around_point (const hp::MappingCollection<dim,spacedim>   &mapping,
                                   const hp::DoFHandler<dim,spacedim> &container,
                                   const Point<spacedim>     &p);

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
				      * Extract the active cells around a given
				      * cell @p cell and return them in the
				      * vector @p active_neighbors.
				      */
    template <class Container>
    static void
    get_active_neighbors (const typename Container::active_cell_iterator        &cell,
			  std::vector<typename Container::active_cell_iterator> &active_neighbors);

				     /**
				      * Produce a sparsity pattern in which
				      * nonzero entries indicate that two
				      * cells are connected via a common
				      * face. The diagonal entries of the
				      * sparsity pattern are also set.
				      *
				      * The rows and columns refer to the
				      * cells as they are traversed in their
				      * natural order using cell iterators.
				      */
    template <int dim, int spacedim>
    static void
    get_face_connectivity_of_cells (const Triangulation<dim, spacedim> &triangulation,
				    SparsityPattern                    &connectivity);

                                     /**
                                      * Use the METIS partitioner to generate
                                      * a partitioning of the active cells
                                      * making up the entire domain. After
                                      * calling this function, the subdomain
                                      * ids of all active cells will have
                                      * values between zero and
                                      * @p n_partitions-1. You can access the
                                      * subdomain id of a cell by using
                                      * <tt>cell-@>subdomain_id()</tt>.
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
    template <int dim, int spacedim>
    static
    void
    partition_triangulation (const unsigned int  n_partitions,
			     Triangulation<dim, spacedim> &triangulation);

				     /**
				      * This function does the same as the
				      * previous one, i.e. it partitions a
				      * triangulation using METIS into a
				      * number of subdomains identified by the
				      * <code>cell-@>subdomain_id()</code>
				      * flag.
				      *
				      * The difference to the previous
				      * function is the second argument, a
				      * sparsity pattern that represents the
				      * connectivity pattern between cells.
				      *
				      * While the function above builds it
				      * directly from the triangulation by
				      * considering which cells neighbor each
				      * other, this function can take a more
				      * refined connectivity graph. The
				      * sparsity pattern needs to be of size
				      * $N\times N$, where $N$ is the number
				      * of active cells in the
				      * triangulation. If the sparsity pattern
				      * contains an entry at position $(i,j)$,
				      * then this means that cells $i$ and $j$
				      * (in the order in which they are
				      * traversed by active cell iterators)
				      * are to be considered connected; METIS
				      * will then try to partition the domain
				      * in such a way that (i) the subdomains
				      * are of roughly equal size, and (ii) a
				      * minimal number of connections are
				      * broken.
				      *
				      * This function is mainly useful in
				      * cases where connections between cells
				      * exist that are not present in the
				      * triangulation alone (otherwise the
				      * previous function would be the simpler
				      * one to use). Such connections may
				      * include that certain parts of the
				      * boundary of a domain are coupled
				      * through symmetric boundary conditions
				      * or integrals (e.g. friction contact
				      * between the two sides of a crack in
				      * the domain), or if a numerical scheme
				      * is used that not only connects
				      * immediate neighbors but a larger
				      * neighborhood of cells (e.g. when
				      * solving integral equations).
				      *
				      * In addition, this function may be
				      * useful in cases where the default
				      * sparsity pattern is not entirely
				      * sufficient. This can happen because
				      * the default is to just consider face
				      * neighbors, not neighboring cells that
				      * are connected by edges or
				      * vertices. While the latter couple when
				      * using continuous finite elements, they
				      * are typically still closely connected
				      * in the neighborship graph, and METIS
				      * will not usually cut important
				      * connections in this case. However, if
				      * there are vertices in the mesh where
				      * many cells (many more than the common
				      * 4 or 6 in 2d and 3d, respectively)
				      * come together, then there will be a
				      * significant number of cells that are
				      * connected across a vertex, but several
				      * degrees removed in the connectivity
				      * graph built only using face
				      * neighbors. In a case like this, METIS
				      * may sometimes make bad decisions and
				      * you may want to build your own
				      * connectivity graph.
				      */
    template <int dim, int spacedim>
    static
    void
    partition_triangulation (const unsigned int     n_partitions,
			     const SparsityPattern &cell_connection_graph,
			     Triangulation<dim,spacedim>    &triangulation);

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
    template <int dim, int spacedim>
    static void
    get_subdomain_association (const Triangulation<dim, spacedim>  &triangulation,
                               std::vector<types::subdomain_id_t> &subdomain);

                                     /**
                                      * Count how many cells are uniquely
                                      * associated with the given @p subdomain
                                      * index.
                                      *
                                      * This function may return zero
                                      * if there are no cells with the
                                      * given @p subdomain index. This
                                      * can happen, for example, if
                                      * you try to partition a coarse
                                      * mesh into more partitions (one
                                      * for each processor) than there
                                      * are cells in the mesh.
				      *
				      * This function returns the number of
				      * cells associated with one
				      * subdomain. If you are looking for the
				      * association of @em DoFs with this
				      * subdomain, use the
				      * <tt>DoFTools::count_dofs_with_subdomain_association</tt>
				      * function.
                                      */
    template <int dim, int spacedim>
    static unsigned int
    count_cells_with_subdomain_association (const Triangulation<dim, spacedim> &triangulation,
                                            const types::subdomain_id_t         subdomain);

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
    template <int dim, int spacedim>
    static
    bool
    have_same_coarse_mesh (const Triangulation<dim, spacedim> &mesh_1,
                           const Triangulation<dim, spacedim> &mesh_2);

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
				      * step-24 for an example
				      * of use of this function.
				      */
    template <int dim, int spacedim>
    static
    double
    minimal_cell_diameter (const Triangulation<dim, spacedim> &triangulation);

				     /**
				      * Return the diamater of the largest
				      * active cell of a triangulation.
				      */
    template <int dim, int spacedim>
    static
    double
    maximal_cell_diameter (const Triangulation<dim, spacedim> &triangulation);

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
    template <int dim, int spacedim>
    static
    void
    create_union_triangulation (const Triangulation<dim, spacedim> &triangulation_1,
				const Triangulation<dim, spacedim> &triangulation_2,
				Triangulation<dim, spacedim>       &result);

				     /**
				      * Given a triangulation and a
				      * list of cells whose children
				      * have become distorted as a
				      * result of mesh refinement, try
				      * to fix these cells up by
				      * moving the center node around.
				      *
				      * The function returns a list of
				      * cells with distorted children
				      * that couldn't be fixed up for
				      * whatever reason. The returned
				      * list is therefore a subset of
				      * the input argument.
				      *
				      * For a definition of the
				      * concept of distorted cells,
				      * see the
				      * @ref GlossDistorted "glossary entry".
				      * The first argument passed to the
				      * current function is typically
				      * the exception thrown by the
				      * Triangulation::execute_coarsening_and_refinement
				      * function.
				      */
    template <int dim, int spacedim>
    static
    typename Triangulation<dim,spacedim>::DistortedCellList
    fix_up_distorted_child_cells (const typename Triangulation<dim,spacedim>::DistortedCellList &distorted_cells,
				  Triangulation<dim,spacedim> &triangulation);

				     /**
				      * This function implements a boundary
				      * subgrid extraction.  Given a
				      * <dim,spacedim>-Triangulation (the
				      * "volume mesh") the function extracts a
				      * subset of its boundary (the "surface
				      * mesh").  The boundary to be extracted
				      * is specified by a list of
				      * boundary_ids.  If none is specified
				      * the whole boundary will be
				      * extracted. The function is used in
				      * step-38.
				      *
				      * It also builds a mapping linking the
				      * cells on the surface mesh to the
				      * corresponding faces on the volume
				      * one. This mapping is the return value
				      * of the function.
				      *
				      * @note The function builds the surface
				      * mesh by creating a coarse mesh from
				      * the selected faces of the coarse cells
				      * of the volume mesh. It copies the
				      * boundary indicators of these faces to
				      * the cells of the coarse surface
				      * mesh. The surface mesh is then refined
				      * in the same way as the faces of the
				      * volume mesh are. In order to ensure
				      * that the surface mesh has the same
				      * vertices as the volume mesh, it is
				      * therefore important that you assign
				      * appropriate boundary objects through
				      * Triangulation::set_boundary to the
				      * surface mesh object before calling
				      * this function. If you don't, the
				      * refinement will happen under the
				      * assumption that all faces are straight
				      * (i.e using the StraightBoundary class)
				      * rather than any curved boundary object
				      * you may want to use to determine the
				      * location of new vertices.
				      *
				      * @note Oftentimes, the
				      * <code>Container</code>
				      * template type will be of kind
				      * Triangulation; in that case,
				      * the map that is returned will
				      * be between Triangulation cell
				      * iterators of the surface mesh
				      * and Triangulation face
				      * iterators of the volume
				      * mesh. However, one often needs
				      * to have this mapping between
				      * DoFHandler (or hp::DoFHandler)
				      * iterators. In that case, you
				      * can pass DoFHandler arguments
				      * as first and second parameter;
				      * the function will in that case
				      * re-build the triangulation
				      * underlying the second argument
				      * and return a map between
				      * DoFHandler iterators. However,
				      * the function will not actually
				      * distribute degrees of freedom
				      * on this newly created surface
				      * mesh.
				      *
				      * @note The algorithm outlined
				      * above assumes that all faces
				      * on higher refinement levels
				      * always have exactly the same
				      * boundary indicator as their
				      * parent face. Consequently, we
				      * can start with coarse level
				      * faces and build the surface
				      * mesh based on that. It would
				      * not be very difficult to
				      * extend the function to also
				      * copy boundary indicators from
				      * finer level faces to their
				      * corresponding surface mesh
				      * cells, for example to
				      * accomodate different geometry
				      * descriptions in the case of
				      * curved boundaries.
				      */
    template <template <int,int> class Container, int dim, int spacedim>
    static
    std::map<typename Container<dim-1,spacedim>::cell_iterator,
	     typename Container<dim,spacedim>::face_iterator>
    extract_boundary_mesh (const Container<dim,spacedim> &volume_mesh,
			   Container<dim-1,spacedim>     &surface_mesh,
			   const std::set<unsigned char> &boundary_ids
			   = std::set<unsigned char>());

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

template <int dim, typename Predicate, int spacedim>
void GridTools::transform (const Predicate    &predicate,
			   Triangulation<dim, spacedim> &triangulation)
{
				   // ensure that all the cells of the
				   // triangulation are on the same level
  Assert (triangulation.n_levels() ==
	  static_cast<unsigned int>(triangulation.begin_active()->level()+1),
  	  ExcMessage ("Not all cells of this triangulation are at the same "
		      "refinement level, as is required for this function."));

  std::vector<bool> treated_vertices (triangulation.n_vertices(),
				      false);

				   // loop over all active cells, and
				   // transform those vertices that
				   // have not yet been touched. note
				   // that we get to all vertices in
				   // the triangulation by only
				   // visiting the active cells.
  typename Triangulation<dim, spacedim>::active_cell_iterator
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
	   child<cell->n_children(); ++child)
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



template <class Container>
void
GridTools::get_active_neighbors(const typename Container::active_cell_iterator        &cell,
				std::vector<typename Container::active_cell_iterator> &active_neighbors)
{
  active_neighbors.clear ();
  for (unsigned int n=0; n<GeometryInfo<Container::dimension>::faces_per_cell; ++n)
    if (! cell->at_boundary(n))
      {
	if (Container::dimension == 1)
	  {
					 // check children of neighbor. note
					 // that in 1d children of the neighbor
					 // may be further refined. In 1d the
					 // case is simple since we know what
					 // children bound to the present cell
	    typename Container::cell_iterator
	      neighbor_child = cell->neighbor(n);
	    if (!neighbor_child->active())
	      {
		while (neighbor_child->has_children())
		  neighbor_child = neighbor_child->child (n==0 ? 1 : 0);

		Assert (neighbor_child->neighbor(n==0 ? 1 : 0)==cell,
			ExcInternalError());
	      }
	    active_neighbors.push_back (neighbor_child);
	  }
	else
	  {
	    if (cell->face(n)->has_children())
					   // this neighbor has children. find
					   // out which border to the present
					   // cell
	      for (unsigned int c=0; c<cell->face(n)->number_of_children(); ++c)
		active_neighbors.push_back (cell->neighbor_child_on_subface(n,c));
	    else
	      {
					     // the neighbor must be active
					     // himself
		Assert(cell->neighbor(n)->active(), ExcInternalError());
		active_neighbors.push_back(cell->neighbor(n));
	      }
	  }
      }
}




// declaration of explicit specializations

template <>
double
GridTools::cell_measure<3>(const std::vector<Point<3> > &all_vertices,
			const unsigned int (&vertex_indices) [GeometryInfo<3>::vertices_per_cell]);

template <>
double
GridTools::cell_measure<2>(const std::vector<Point<2> > &all_vertices,
			const unsigned int (&vertex_indices) [GeometryInfo<2>::vertices_per_cell]);

// double
// GridTools::cell_measure<2,3>(const std::vector<Point<3> > &all_vertices,
// 			const unsigned int (&vertex_indices) [GeometryInfo<2>::vertices_per_cell]);







DEAL_II_NAMESPACE_CLOSE

/*----------------------------   grid_tools.h     ---------------------------*/
/* end of #ifndef __deal2__grid_tools_H */
#endif
/*----------------------------   grid_tools.h     ---------------------------*/

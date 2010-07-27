//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__tria_h
#define __deal2__tria_h


#include <base/config.h>
#include <base/point.h>
#include <base/subscriptor.h>
#include <base/smartpointer.h>
#include <base/geometry_info.h>
#include <grid/tria_iterator_selector.h>

#include <vector>
#include <list>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim> class Boundary;
template <int dim, int spacedim> class StraightBoundary;

template <int, int, int> class TriaAccessor;
template <int, int, int> class TriaAccessor;

namespace internal
{
  namespace Triangulation
  {
    template <int dim> class TriaLevel;
    template <int dim> class TriaFaces;

				     /**
				      * Forward declaration of a class into
				      * which we put much of the
				      * implementation of the Triangulation
				      * class. See the .cc file for more
				      * information.
				      */
    struct Implementation;
  }

  namespace TriaAccessor
  {
    struct Implementation;
  }
}

template <int dim, int spacedim> class DoFHandler;
namespace hp
{
  template <int dim, int spacedim> class DoFHandler;
}
template <int dim, int spacedim> class MGDoFHandler;

/*------------------------------------------------------------------------*/

/**
 * Structure which is passed to Triangulation::create_triangulation. It
 * contains all data needed to construct a cell, namely the indices of the
 * vertices and the material indicator.
 *
 * This structure is also used to represent data for faces and edge as part of
 * the SubCellData class. In that case the #vertices array needs to represent
 * the vertices of a face or edge of a cell listed in the argument to
 * Triangulation::create_triangulation that denotes the cells. It can be used
 * to attach boundary indicators to faces.
 *
 * @ingroup grid
 */
template <int structdim>
struct CellData
{
				     /**
				      * Indices of the vertices of
				      * this cell.
				      */
    unsigned int vertices[GeometryInfo<structdim>::vertices_per_cell];

				     /**
				      * Material indicator of this
				      * cell. May be used to denote
				      * different coefficients, etc.
				      *
				      * Note that if this object is part of a
				      * SubCellData object, then it represents
				      * a face or edge of a cell. Since
				      * deal.II only allows material_id tags
				      * on cells, this field is re-interpreted
				      * to indicate the boundary_indicator of
				      * a face or edge. In other words, for
				      * faces, this field is a misnomer.
				      */
    unsigned char material_id;
};



/**
 *  Structure to be passed to Triangulation::create_triangulation function to
 *  describe boundary information.
 *
 *  This structure is the same for all dimensions, since we use an input
 *  function which is the same for all dimensions. The content of objects
 *  of this structure varies with the dimensions, however.
 *
 *  Since in one dimension, there is no boundary information apart
 *  from the two end points of the interval, this structure does not contain
 *  anything and exists only for consistency, to allow a common interface
 *  for all space dimensions. All fields should always be empty.
 *
 *  Boundary data in 2D consists
 *  of a list of lines which belong to a given boundary component. A
 *  boundary component is a list of lines which are given a common
 *  number describing the boundary condition to hold on this part of the
 *  boundary. The triangulation creation function gives lines not in this
 *  list either the boundary indicator zero (if on the boundary) or 255
 *  (if in the interior). Explicitely giving a line the indicator 255
 *  will result in an error, as well as giving an interior line a boundary
 *  indicator.
 *
 * @ingroup grid
 */
struct SubCellData
{
				     /**
				      * Each record of this vector describes
				      * a line on the boundary and its boundary
				      * indicator.
				      */
    std::vector<CellData<1> > boundary_lines;

				     /**
				      * Each record of this vector describes
				      * a quad on the boundary and its boundary
				      * indicator.
				      */
    std::vector<CellData<2> > boundary_quads;

				     /**
				      * This function checks whether the vectors
				      * which may not be used in a given
				      * dimension are really empty. I.e.,
				      * whether the <tt>boundary_*</tt> arrays are
				      * empty when in one space dimension
				      * and whether the @p boundary_quads
				      * array is empty when in two dimensions.
				      *
				      * Since this structure is the same for all
				      * dimensions, the actual dimension has
				      * to be given as a parameter.
				      */
    bool check_consistency (const unsigned int dim) const;
};


/*------------------------------------------------------------------------*/


namespace internal
{
				 /**
				  * A namespace for classes internal to the
				  * triangulation classes and helpers.
				  */
  namespace Triangulation
  {

/**
 * Cache class used to store the number of used and active elements
 * (lines or quads etc) within the levels of a triangulation. This
 * is only the declaration of the template, concrete instantiations
 * are below.
 *
 * In the old days, whenever one wanted to access one of these
 * numbers, one had to perform a loop over all lines, e.g., and count
 * the elements until we hit the end iterator. This is time consuming
 * and since access to the number of lines etc is a rather frequent
 * operation, this was not an optimal solution.
 *
 * @ingroup grid
 * @author Wolfgang Bangerth, 1999
 */
    template <int dim>
    struct NumberCache
    {
    };

/**
 * Cache class used to store the number of used and active elements
 * (lines or quads etc) within the levels of a triangulation. This
 * specialization stores the numbers of lines.
 *
 * In the old days, whenever one wanted to access one of these
 * numbers, one had to perform a loop over all lines, e.g., and count
 * the elements until we hit the end iterator. This is time consuming
 * and since access to the number of lines etc is a rather frequent
 * operation, this was not an optimal solution.
 *
 * @ingroup grid
 * @author Wolfgang Bangerth, 1999
 */
    template <>
    struct NumberCache<1>
    {
					 /**
					  * The number of levels on
					  * which we have used
					  * objects.
					  */
	unsigned int n_levels;

                                         /**
                                          * Number of used lines in the whole
                                          * triangulation.
                                          */
        unsigned int n_lines;

                                         /**
                                          * Array holding the number of used
                                          * lines on each level.
                                          */
        std::vector<unsigned int> n_lines_level;

                                         /**
                                          * Number of active lines in the
                                          * whole triangulation.
                                          */
        unsigned int n_active_lines;

                                         /**
                                          * Array holding the number of active
                                          * lines on each level.
                                          */
        std::vector<unsigned int> n_active_lines_level;

                                         /**
                                          * Constructor. Set values to zero
                                          * by default.
                                          */
        NumberCache ();

                                         /**
                                          * Determine an estimate for the
                                          * memory consumption (in bytes)
                                          * of this object.
                                          */
        unsigned int memory_consumption () const;
    };


/**
 * Cache class used to store the number of used and active elements
 * (lines or quads etc) within the levels of a triangulation. This
 * specialization stores the numbers of quads. Due to the inheritance
 * from the base class NumberCache<1>, the numbers of lines
 * are also within this class.
 *
 * In the old days, whenever one wanted to access one of these
 * numbers, one had to perform a loop over all lines, e.g., and count
 * the elements until we hit the end iterator. This is time consuming
 * and since access to the number of lines etc is a rather frequent
 * operation, this was not an optimal solution.
 *
 * @ingroup grid
 * @author Wolfgang Bangerth, 1999
 */
    template <>
    struct NumberCache<2> : public NumberCache<1>
    {
                                         /**
                                          * Number of used quads in the whole
                                          * triangulation.
                                          */
        unsigned int n_quads;

                                         /**
                                          * Array holding the number of used
                                          * quads on each level.
                                          */
        std::vector<unsigned int> n_quads_level;

                                         /**
                                          * Number of active quads in the
                                          * whole triangulation.
                                          */
        unsigned int n_active_quads;

                                         /**
                                          * Array holding the number of active
                                          * quads on each level.
                                          */
        std::vector<unsigned int> n_active_quads_level;

                                         /**
                                          * Constructor. Set values to zero
                                          * by default.
                                          */
        NumberCache ();

                                         /**
                                          * Determine an estimate for the
                                          * memory consumption (in bytes)
                                          * of this object.
                                          */
        unsigned int memory_consumption () const;
    };

//TODO: Replace boundary[255] by a std::vector so we can use constructor of SmartPointer

/**
 * Cache class used to store the number of used and active elements
 * (lines or quads etc) within the levels of a triangulation. This
 * specialization stores the numbers of hexes. Due to the inheritance
 * from the base class NumberCache<2>, the numbers of lines
 * and quads are also within this class.
 *
 * In the old days, whenever one wanted to access one of these
 * numbers, one had to perform a loop over all lines, e.g., and count
 * the elements until we hit the end . This is time consuming
 * and since access to the number of lines etc is a rather frequent
 * operation, this was not an optimal solution.
 *
 * @ingroup grid
 * @author Wolfgang Bangerth, 1999
 */
    template <>
    struct NumberCache<3> : public NumberCache<2>
    {
                                         /**
                                          * Number of used hexes in the whole
                                          * triangulation.
                                          */
        unsigned int n_hexes;

                                         /**
                                          * Array holding the number of used
                                          * hexes on each level.
                                          */
        std::vector<unsigned int> n_hexes_level;

                                         /**
                                          * Number of active hexes in the
                                          * whole triangulation.
                                          */
        unsigned int n_active_hexes;

                                         /**
                                          * Array holding the number of active
                                          * hexes on each level.
                                          */
        std::vector<unsigned int> n_active_hexes_level;

                                         /**
                                          * Constructor. Set values to zero
                                          * by default.
                                          */
        NumberCache ();

                                         /**
                                          * Determine an estimate for the
                                          * memory consumption (in bytes)
                                          * of this object.
                                          */
        unsigned int memory_consumption () const;
    };
  }
}


/*------------------------------------------------------------------------*/


/**
 *  Triangulations denote a hierarchy of levels of elements which together
 *  form a @p dim -dimensional manifold in @p spacedim spatial dimensions
 *  (if spacedim is not specified it takes the default value @ spacedim=dim).
 *
 *  Thus, for example, an object of type @p Triangulation<1,1> (or simply
 *  @p Triangulation<1> since @p spacedim==dim by default) is used to represent
 *  and handle the usual one-dimensional triangulation used in the finite
 *  element method (so, segments on a straight line). On the other hand,
 *  objects such as @p Triangulation<1,2> or @p Triangulation<2,3> (that
 *  are associated with curves in 2D or surfaces in 3D)
 *  are the ones one wants to use in the boundary element method.
 *
 *  This class is written to be as independent of the dimension as possible
 *  (thus the complex construction of the internal::Triangulation::TriaLevel
 *  classes) to allow code-sharing, to allow reducing the need to mirror
 *  changes in the code for one dimension to the code for other
 *  dimensions. Nonetheless, some of the functions are dependent of the
 *  dimension and there only exist specialized versions for distinct
 *  dimensions.
 *
 *
 *  <h3>Structure and iterators</h3>
 *
 *  The actual data structure of a Triangulation object is rather complex
 *  and quite inconvenient if one attempted to operate on it directly, since
 *  data is spread over quite a lot of arrays and other places. However,
 *  there are ways powerful enough to work on these data structures
 *  without knowing their exact relations. This is done through the
 *  concept of iterators (see the STL documentation and TriaRawIterator).
 *  In order to make things as easy and dimension independent as possible,
 *  use of class local typedefs is made, see below.
 *
 *  The Triangulation class provides iterator which enable looping over all
 *  lines, cells, etc without knowing the exact representation used to
 *  describe them. Their names are typedefs imported from the Iterators
 *  class (thus making them local types to this class) and are as follows:
 *
 *  <ul>
 *  <li> @p raw_line_iterator: loop over all lines, used or not (declared for
 *  all dimensions).
 *
 *  <li> @p line_iterator: loop over all used lines (declared for all dimensions).
 *
 *  <li> @p active_line_iterator: loop over all active lines (declared for all
 *  dimensions).
 *
 *  <li> @p raw_quad_iterator: loop over all quads, used or not (declared only
 *  for <tt>dim>=2</tt>).
 *
 *  <li> @p quad_iterator: loop over all quads (declared only for @p dim>=2).
 *
 *  <li> @p active_quad_iterator: loop over all active quads (declared only for
 *  @p dim>=2).
 *  </ul>
 *
 *  Additionaly, for @p dim==1, the following identities hold:
 *  @verbatim
 *    typedef raw_line_iterator    raw_cell_iterator;
 *    typedef line_iterator        cell_iterator;
 *    typedef active_line_iterator active_cell_iterator;
 *  @endverbatim
 *  while for @p dim==2
 *  @verbatim
 *    typedef quad_line_iterator   raw_cell_iterator;
 *    typedef quad_iterator        cell_iterator;
 *    typedef active_quad_iterator active_cell_iterator;
 *
 *    typedef raw_line_iterator    raw_face_iterator;
 *    typedef line_iterator        face_iterator;
 *    typedef active_line_iterator active_face_iterator;
 *  @endverbatim
 *
 *  By using the cell iterators, you can write code nearly independent of
 *  the spatial dimension. The same applies for substructure iterators,
 *  where a substructure is defined as a face of a cell. The face of a
 *  cell is a vertex in 1D and a line in 2D; however, vertices are
 *  handled in a different way and therefore lines have no faces.
 *
 *  The Triangulation class offers functions like @p begin_active which gives
 *  you an iterator to the first active cell. There are quite a lot of functions
 *  returning iterators. Take a look at the class doc to get an overview.
 *
 *  Usage of these iterators works mostly like with the STL iterators. Some
 *  examples taken from the Triangulation source code follow (notice that in the last
 *  two examples the template parameter @p spacedim has been omitted, so it takes
 *  the default value <code>dim</code>).
 *  <ul>
 *  <li> <em>Counting the number of cells on a specific level</em>
 *    @verbatim
 *     template <int dim, int spacedim>
 *     int Triangulation<dim, spacedim>::n_cells (const int level) const {
 *        cell_iterator cell = begin (level),
 *                      endc = end(level);
 *        int n=0;
 *        for (; cell!=endc; ++cell)
 *          ++n;
 *        return n;
 *      };
 *    @endverbatim
 *    Another way which uses the STL @p distance function would be to write
 *    @verbatim
 *      template <int dim>
 *      int Triangulation<dim>::n_cells (const int level) const {
 *        int n=0;
 *        distance (begin(level),
 *                  (level == levels.size()-1 ?
 *                   cell_iterator(end()) :
 *                   begin (level+1)),
 *                  n);
 *        return n;
 *      };
 *    @endverbatim
 *
 *  <li> <em>Refining all cells of a triangulation</em>
 *    @verbatim
 *      template <int dim>
 *      void Triangulation<dim>::refine_global () {
 *        active_cell_iterator cell = begin_active(),
 *                             endc = end();
 *
 *        for (; cell != endc; ++cell)
 *          cell->set_refine_flag ();
 *        execute_coarsening_and_refinement ();
 *      };
 *    @endverbatim
 *  </ul>
 *
 *
 *  <h3>Usage</h3>
 *
 *  Usage of a Triangulation is mainly done through the use of iterators.
 *  An example probably shows best how to use it:
 *  @verbatim
 *  void main () {
 *    Triangulation<2> tria;
 *
 *    // read in a coarse grid file
 *
 *                                     // we want to log the
 *                                     // refinement history
 *    ofstream history ("mesh.history");
 *
 *                                     // refine first cell
 *    tria.begin_active()->set_refine_flag();
 *    tria.save_refine_flags (history);
 *    tria.execute_coarsening_and_refinement ();
 *
 *                                     // refine first active cell
 *                                     // on coarsest level
 *    tria.begin_active()->set_refine_flag ();
 *    tria.save_refine_flags (history);
 *    tria.execute_coarsening_and_refinement ();
 *
 *    Triangulation<2>::active_cell_iterator cell;
 *    for (int i=0; i<17; ++i)
 *      {
 *                                         // refine the presently
 *                                         // second last cell 17
 *                                         // times
 *        cell = tria.last_active(tria.n_levels()-1);
 *        --cell;
 *        cell->set_refine_flag ();
 *        tria.save_refine_flags (history);
 *        tria.execute_coarsening_and_refinement ();
 *      };
 *                                       // output the grid
 *    ofstream out("grid.1");
 *    GridOut::write_gnuplot (tria, out);
 *  };
 *  @endverbatim
 *
 *
 *  <h3>Creating a triangulation</h3>
 *
 *  There are several possibilities to create a triangulation:
 *  <ul>
 *    <li> The most common domains, such as hypercubes (i.e. lines, squares,
 *       cubes, etc), hyper-balls (circles, balls, ...) and some other, more
 *       weird domains such as the L-shape region and higher dimensional
 *       generalizations and others, are provided by the GridGenerator
 *       class which takes a triangulation and fills it by a division
 *       of the required domain.
 *
 *     <li> Reading in a triangulation: By using an object of the GridIn
 *        class, you can read in fairly general triangulations. See there for
 *        more information. The mentioned class uses the interface described
 *        directly below to transfer the data into the triangulation.
 *
 *     <li> Explicitely creating a triangulation: you can create a triangulation
 *        by providing a list of vertices and a list of cells. Each such cell
 *        consists of a vector storing the indices of the vertices of this cell
 *        in the vertex list. To see how this works, you can take a look at the
 *        GridIn<dim>::read_* functions. The appropriate function to be
 *        called is create_triangulation().
 *
 *        Creating the hierarchical information needed for this
 *        library from cells storing only vertex information can be
 *        quite a complex task.  For example in 2D, we have to create
 *        lines between vertices (but only once, though there are two
 *        cells which link these two vertices) and we have to create
 *        neighborship information. Grids being read in should
 *        therefore not be too large, reading refined grids would be
 *        inefficient (although there is technically no problem in
 *        reading grids with several 10.000 or 100.000 cells; the
 *        library can handle this without much problems). Apart from
 *        the performance aspect, refined grids do not lend too well
 *        to multigrid algorithms, since solving on the coarsest level
 *        is expensive. It is wiser in any case to read in a grid as
 *        coarse as possible and then do the needed refinement steps.
 *
 *        It is your duty to guarantee that cells have the correct orientation.
 *        To guarantee this, in the input vector keeping the cell list, the
 *        vertex indices for each cell have to be in a defined order, see the
 *        documentation of GeometryInfo<dim>. In one dimension, the first vertex
 *        index must refer to that vertex with the lower coordinate value. In 2D
 *        and 3D, the correspondoing conditions are not easy to verify and no
 *        full attempt to do so is made.
 *        If you violate this condition, you may end up with matrix entries
 *        having the wrong sign (clockwise vertex numbering, which results in
 *        a negative area element) of with wrong matrix elements (twisted
 *        quadrilaterals, i.e. two vertices interchanged; this results in
 *        a wrong area element).
 *
 *        There are more subtle conditions which must be imposed upon
 *        the vertex numbering within cells. They do not only hold for
 *        the data read from an UCD or any other input file, but also
 *        for the data passed to the
 *        <tt>Triangulation<dim,spacedim>::create_triangulation ()</tt>
 *        function. See the documentation for the GridIn class
 *        for more details on this, and above all to the
 *        GridReordering class that explains many of the
 *        problems and an algorithm to reorder cells such that they
 *        satisfy the conditions outlined above.
 *
 *     <li> Copying a triangulation: when computing on time dependent meshes
 *        or when using adaptive refinement, you will often want to create a
 *        new triangulation to be the same as another one. This is facilitated
 *        by the @p copy_triangulation function.
 *
 *        It is guaranteed that vertex, line or cell numbers in the two
 *        triangulations are the same and that two iterators walking on the
 *        two triangulations visit matching cells if they are incremented in
 *        parallel. It may be conceivable to implement a clean-up in the copy
 *        operation, which eliminates holes of unused memory, re-joins
 *        scattered data and so on. In principle this would be a useful
 *        operation but guaranteeing some parallelity in the two triangulations
 *        seems more important since usually data will have to be transferred
 *        between the grids.
 *   </ul>
 *
 *   Finally, there is a special function for folks who like bad grids:
 *   <tt>Triangulation<dim,spacedim>::distort_random</tt>. It moves all the vertices in the
 *   grid a bit around by a random value, leaving behind a distorted mesh.
 *   Note that you should apply this function to the final mesh, since
 *   refinement smoothes the mesh a bit.
 *
 *   The function will make sure that vertices on restricted faces (hanging
 *   nodes) will result in the correct place, i.e. in the middle of the two
 *   other vertices of the mother line, and the analogue in higher space
 *   dimensions (vertices on the boundary are not corrected, so don't distort
 *   boundary vertices in more than two space dimension, i.e. in dimensions
 *   where boundary vertices can be hanging nodes). Applying the algorithm
 *   has another drawback related to the
 *   placement of cells, however: the children of a cell will not occupy the
 *   same region of the domain as the mother cell does. While this is the
 *   usual behavior with cells at the boundary, here you may get into trouble
 *   when using multigrid algorithms or when transferring solutions from coarse
 *   to fine grids and back. In general, the use of this function is only safe
 *   if you only use the most refined level of the triangulation for
 *   computations.
 *
 *
 *
 *   <h3>Refinement and coarsening of a triangulation</h3>
 *
 *   Refinement of a triangulation may be done through several ways. The most
 *   low-level way is directly through iterators: let @p i be an iterator to
 *   an active cell (i.e. the cell pointed to has no children), then the
 *   function call <tt>i->set_refine_flag()</tt> marks the respective cell for
 *   refinement. Marking non-active cells results in an error.
 *
 *   After all the cells you wanted to mark for refinement, call the
 *   @p execute_coarsening_and_refinement function to actually perform
 *   the refinement. This function itself first calls the
 *   @p prepare_coarsening_and_refinement function to regularize the resulting
 *   triangulation: since a face between two adjacent cells may only
 *   be subdivided once (i.e. the levels of two adjacent cells may
 *   differ by one at most; it is not possible to have a cell refined
 *   twice while the neighboring one is not refined), some additional
 *   cells are flagged for refinement to smooth the grid. This
 *   enlarges the number of resulting cells but makes the grid more
 *   regular, thus leading to better approximation properties and,
 *   above all, making the handling of data structures and algorithms
 *   much much easier. To be honest, this is mostly an algorithmic
 *   step than one needed by the finite element method.
 *
 *   To coarsen a grid, the same way as above is possible by using
 *   <tt>i->set_coarsen_flag</tt> and calling @p execute_coarsening_and_refinement.
 *
 *   The reason for first coarsening, then refining is that the
 *   refinement usually adds some additional cells to keep the triangulation
 *   regular and thus satisfies all refinement requests, while the coarsening
 *   does not delete cells not requested for; therefore the refinement will
 *   often revert some effects of coarsening while the opposite is not true.
 *   The stated order of coarsening before refinement will thus normally
 *   lead to a result closer to the intended one.
 *
 *   Marking cells for refinement 'by hand' through iterators is one way to
 *   produce a new grid, especially if you know what kind of grid you are
 *   looking for, e.g. if you want to have a grid successively refined
 *   towards the boundary or always at the center (see the example programs,
 *   they do exactly these things). There are more advanced functions,
 *   however, which are more suitable for automatic generation of hierarchical
 *   grids in the context of a-posteriori error estimation and adaptive finite
 *   elements. These functions can be found in the GridRefinement class.
 *
 *
 *   <h3>Smoothing of a triangulation</h3>
 *
 *   Some degradation of approximation properties has been observed
 *   for grids which are too unstructured. Therefore, the
 *   @p prepare_coarsening_and_refinement function which is automatically called
 *   by the @p execute_coarsening_and_refinement function can do some
 *   smoothing of the triangulation. Note that mesh smoothing is only
 *   done for two or more space dimensions, no smoothing is available
 *   at present for one spatial dimension. In the sequel, let
 *   <tt>execute_*</tt> stand for @p execute_coarsening_and_refinement.
 *
 *   For the purpose of smoothing, the
 *   Triangulation constructor takes an argument specifying whether a
 *   smoothing step shall be performed on the grid each time <tt>execute_*</tt>
 *   is called. The default is that such a step not be done, since this results
 *   in additional cells being produced, which may not be necessary in all
 *   cases. If switched on, calling <tt>execute_*</tt> results in
 *   flagging additional cells for refinement to avoid
 *   vertices as the ones mentioned. The algorithms for both regularization
 *   and smoothing of triangulations are described below in the section on
 *   technical issues. The reason why this parameter must be given to the
 *   constructor rather than to <tt>execute_*</tt> is that it would result
 *   in algorithmic problems if you called <tt>execute_*</tt> once without
 *   and once with smoothing, since then in some refinement steps would need
 *   to be refined twice.
 *
 *   The parameter taken by the constructor is an integer which may be composed
 *   bitwise by the constants defined in the <tt>enum MeshSmoothing</tt>. The meaning
 *   of these constants is explained in the following:
 *   <ul>
 *   <li> @p limit_level_difference_at_vertices:
 *     It can be shown, that degradation of approximation occurs if the
 *     triangulation contains vertices which are member of cells with levels
 *     differing by more than one. One such example is the following:
 *
 *     @image html limit_level_difference_at_vertices.png ""
 *
 *     It would seem that in two space dimensions, the maximum jump in levels
 *     between cells sharing a common vertex is two (as in the example
 *     above). However, this is not true if more than four cells meet at a
 *     vertex. It is not uncommon that a coarse (initial) mesh contains
 *     vertices at which six or even eight cells meet, when small features of
 *     the domain have to be resolved even on the coarsest mesh. In that case,
 *     the maximum difference in levels is three or four, respectively. The
 *     problem gets even worse in three space dimensions.
 *
 *     Looking at an interpolation of the second derivative of the finite
 *     element solution (assuming bilinear finite elements), one sees that the
 *     numerical solution is almost totally wrong, compared with the true
 *     second derivative. Indeed, on regular meshes, there exist sharp
 *     estimations that the $H^2$-error is only $O(1)$, so we should not be
 *     surprised; however, the numerical solution may show a value for the
 *     second derivative which may be a factor of ten away from the true
 *     value. These problems are located on the small cell adjacent to the
 *     center vertex, where cells of non-subsequent levels meet, as well as on
 *     the upper and right neighbor of this cell (but with a less degree of
 *     deviation from the true value).
 *
 *     If the smoothing indicator given to the constructor contains the bit for
 *     @p limit_level_difference_at_vertices, situations as the above one are
 *     eliminated by also marking the lower left cell for refinement.
 *
 *     In case of anisotropic refinement, the level of a cell is not linked to
 *     the refinement of a cell as directly as in case of isotropic
 *     refinement. Furthermore, a cell can be strongly refined in one
 *     direction and not or at least much less refined in another. Therefore,
 *     it is very difficult to decide, which cases should be excluded from the
 *     refinement process. As a consequence, when using anisotropic
 *     refinement, the @p limit_level_difference_at_vertices flag must not be
 *     set. On the other hand, the implementation of multigrid methods in
 *     deal.II requires that this bit be set.
 *
 *   <li> @p eliminate_unrefined_islands:
 *     Single cells which are not refined and are surrounded by cells which are
 *     refined usually also lead to a sharp decline in approximation properties
 *     locally. The reason is that the nodes on the faces between unrefined and
 *     refined cells are not real degrees of freedom but carry constraints. The
 *     patch without additional degrees of freedom is thus significantly larger
 *     then the unrefined cell itself. If in the parameter passed to the
 *     constructor the bit for @p eliminate_unrefined_islands is set, all cells
 *     which are not flagged for refinement but which are surrounded by more
 *     refined cells than unrefined cells are flagged for refinement. Cells
 *     which are not yet refined but flagged for that are accounted for the
 *     number of refined neighbors. Cells on the boundary are not accounted for
 *     at all. An unrefined island is, by this definition
 *     also a cell which (in 2D) is surrounded by three refined cells and one
 *     unrefined one, or one surrounded by two refined cells, one unrefined one
 *     and is at the boundary on one side. It is thus not a true island, as the
 *     name of the flag may indicate. However, no better name came to mind to
 *     the author by now.
 *
 *   <li> <tt>eliminate_refined_*_islands</tt>:
 *     This algorithm seeks for isolated cells which are refined or flagged
 *     for refinement. This definition is unlike that for
 *     @p eliminate_unrefined_islands, which would mean that an island is
 *     defined as a cell which
 *     is refined but more of its neighbors are not refined than are refined.
 *     For example, in 2D, a cell's refinement would be reverted if at most
 *     one of its neighbors is also refined (or refined but flagged for
 *     coarsening).
 *
 *     The reason for the change in definition of an island is, that this
 *     option would be a bit dangerous, since if you consider a
 *     chain of refined cells (e.g. along a kink in the solution), the cells
 *     at the two ends would be coarsened, after which the next outermost cells
 *     would need to be coarsened. Therefore, only one loop of flagging cells
 *     like this could be done to avoid eating up the whole chain of refined
 *     cells (`chain reaction'...).
 *
 *     This algorithm also takes into account cells which are not actually
 *     refined but are flagged for refinement. If necessary, it takes away the
 *     refinement flag.
 *
 *     Actually there are two versions of this flag,
 *     @p eliminate_refined_inner_islands and @p eliminate_refined_boundary_islands.
 *     There first eliminates islands defined by the definition above which are
 *     in the interior of the domain, while the second eliminates only those
 *     islands if the cell is at the boundary. The reason for this split of
 *     flags is that one often wants to eliminate such islands in the interior
 *     while those at the boundary may well be wanted, for example if one
 *     refines the mesh according to a criterion associated with a boundary
 *     integral or if one has rough boundary data.
 *
 *   <li> @p do_not_produce_unrefined_islands:
 *     This flag prevents the occurrence of unrefined islands. In more detail:
 *     It prohibits the coarsening of a cell if 'most of the neighbors' will
 *     be refined after the step.
 *
 *   <li> @p patch_level_1:
 *     A triangulation of patch level 1 consists of patches, i.e. of
 *     cells that are refined once. This flag ensures that a mesh of
 *     patch level 1 is still of patch level 1 after coarsening and
 *     refinement. It is, however, the user's responsibility to ensure
 *     that the mesh is of patch level 1 before calling
 *     execute_coarsening_and_refinement the first time. The easiest
 *     way to achieve this is by calling global_refine(1) straight
 *     after creation of the triangulation.  It follows that if at
 *     least one of the children of a cell is or will be refined than
 *     all children need to be refined. If the @p patch_level_1 flag
 *     is set, than the flags @p eliminate_unrefined_islands, @p
 *     eliminate_refined_inner_islands and @p
 *     eliminate_refined_boundary_islands will be ignored as they will
 *     be fulfilled automatically.
 *
 *   <li> @p coarsest_level_1:
 *     Each coarse grid cell is refined at least once, i.e. the
 *     triangulation might have active cells on level 1 but not on
 *     level 0. This flag ensures that a mesh which has
 *     coarsest_level_1 has still coarsest_level_1 after coarsening
 *     and refinement. It is, however, the user's responsibility to
 *     ensure that the mesh has coarsest_level_1 before calling
 *     execute_coarsening_and_refinement the first time. The easiest
 *     way to achieve this is by calling global_refine(1) straight
 *     after creation of the triangulation. It follows that active
 *     cells on level 1 may not be coarsenend.
 *
 *     The main use of this flag is to ensure that each cell has at least one
 *     neighbor in each coordinate direction (i.e. each cell has at least a
 *     left or right, and at least an upper or lower neighbor in 2d). This is
 *     a necessary precondition for some algorihms that compute finite
 *     differences between cells. The DerivativeApproximation class is one of
 *     these algorithms that require that a triangulation is coarsest_level_1
 *     unless all cells already have at least one neighbor in each coordinate
 *     direction on the coarsest level.
 *
 *   <li> @p smoothing_on_refinement:
 *     This flag sums up all smoothing algorithms which may be performed upon
 *     refinement by flagging some more cells for refinement.
 *
 *   <li> @p smoothing_on_coarsening:
 *     This flag sums up all smoothing algorithms which may be performed upon
 *     coarsening by flagging some more cells for coarsening.
 *
 *   <li> @p maximum_smoothing:
 *     This flag includes all the above ones and therefore combines all
 *     smoothing algorithms implemented.
 *
 *   <li> @p allow_anisotropic_smoothing:
 *     This flag is not included in @p maximum_smoothing. The flag is
 *     concerned with the following case: consider the case that an
 *     unrefined and a refined cell share a common face and that one
 *     of the children of the refined cell along the common face is
 *     flagged for further refinement. In that case, the resulting
 *     mesh would have more than one hanging node along one or more of
 *     the edges of the triangulation, a situation that is not
 *     allowed. Consequently, in order to perform the refinement, the
 *     coarser of the two original cells is also going to be refined.
 *
 *     However, in many cases it is sufficient to refine the coarser
 *     of the two original cells in an anisotropic way to avoid the
 *     case of multiple hanging vertices on a single edge. Doing only
 *     the minimal anisotropic refinement can save cells and degrees
 *     of freedom. By specifying this flag, the library can produce
 *     these anisotropic refinements.
 *
 *     The flag is not included by default since it may lead to
 *     anisotropically refined meshes even though no cell has ever
 *     been refined anisotropically explicitly by a user command. This
 *     surprising fact may lead to programs that do the wrong thing
 *     since they are not written for the additional cases that can
 *     happen with anisotropic meshes, see the discussion in the
 *     introduction to step-30.
 *
 *   <li> @p none:
 *     Select no smoothing at all.
 *   </ul>
 *
 *
 *   <h3>Material and boundary information</h3>
 *
 *   Each line, quad, etc stores one byte of information denoting the
 *   material or the part of the boundary that an object
 *   belongs to. The material of a cell may be used
 *   during matrix generation in order to implement different
 *   coefficients in different parts of the domain. It is not used by
 *   functions of the grid and dof handling libraries.
 *
 *   This material_id may be set upon construction of a
 *   triangulation (through the CellData data structure), or later
 *   through use of cell iterators. For a typical use of this
 *   functionality, see the step-28 tutorial
 *   program. The functions of the GridGenerator namespace typically
 *   set the material ID of all cells to zero. When reading a
 *   triangulation, the material id must be specified in the input
 *   file (UCD format) or is otherwise set to zero. Material IDs are
 *   inherited by child cells from their parent upon mesh refinement.
 *
 *   Boundary indicators on lower dimensional objects (these have no
 *   material id) indicate the number of a boundary component. These
 *   are used for two purposes: First, they specify a boundary
 *   curve. When a cell is refined, a function can be used to place
 *   new vertices on this curve. See the section on boundary
 *   approximation below. Furthermore, the weak formulation of the
 *   partial differential equation may have different boundary
 *   conditions on different parts of the boundary. The boundary
 *   indicator can be used in creating the matrix or the right hand
 *   side vector to indicate these different parts of the model (this
 *   use is like the material id of cells).

 *   Material and boundary indicators may be in the range from zero to
 *   254. The value 255 is reserved to denote interior lines (in 2D)
 *   and interior lines and quads (in 3D), which do not have a
 *   boundary or material indicator. This way, a program can easily
 *   determine, whether such an object is at the boundary or not.
 *
 *   Lines in two dimensions and quads in three dimensions inherit their
 *   boundary indicator to their children upon refinement. You should therefore
 *   make sure that if you have different boundary parts, the different parts
 *   are separated by a vertex (in 2D) or a line (in 3D) such that each boundary
 *   line or quad has a unique boundary indicator.
 *
 *   Since in one dimension, no substructures of lower dimension exist to
 *   cells (of course apart from vertices, but these are handled
 *   in another way than the structures and substructures with dimension one
 *   or greater), there is no way to denote boundary indicators to boundary
 *   vertices (the endpoints). This is not a big thing, however, since you
 *   will normally not want to do a loop over all vertices, but rather work
 *   on the cells, which in turn have a possibility to find out whether they
 *   are at one of the endpoints. Only the handling of boundary values gets
 *   a bit more complicated, but this seems to be the price to be paid for
 *   the different handling of vertices from lines and quads.
 *
 *
 *
 *   <h3>History of a triangulation</h3>
 *
 *   It is possible to reconstruct a grid from its refinement history, which
 *   can be stored and loaded through the @p save_refine_flags and
 *   @p load_refine_flags functions. Normally, the code will look like this:
 *   @verbatim
 *                                 // open output file
 *     ofstream history("mesh.history");
 *                                 // do 10 refinement steps
 *     for (int step=0; step<10; ++step) {
 *       ...;
 *       // flag cells according to some criterion
 *       ...;
 *       tria.save_refine_flags (history);
 *       tria.execute_coarsening_and_refinement ();
 *     };
 *   @endverbatim
 *
 *   If you want to re-create the grid from the stored information, you write:
 *   @verbatim
 *                                 // open input file
 *     ifstream history("mesh.history");
 *                                 // do 10 refinement steps
 *     for (int step=0; step<10; ++step) {
 *       tria.load_refine_flags (history);
 *       tria.execute_coarsening_and_refinement ();
 *     };
 *   @endverbatim
 *
 *   The same scheme is employed for coarsening and the coarsening flags.
 *
 *   You may write other information to the output file between different sets
 *   of refinement information, as long as you read it upon re-creation of the
 *   grid. You should make sure that the other information in the new
 *   triangulation which is to be created from the saved flags, matches that of
 *   the old triangulation, for example the smoothing level; if not, the
 *   cells actually created from the flags may be other ones, since smoothing
 *   adds additional cells, but their number may be depending on the smoothing
 *   level.
 *
 *   There actually are two sets of <tt>save_*_flags</tt> and <tt>load_*_flags</tt> functions.
 *   One takes a stream as argument and reads/writes the information from/to the
 *   stream, thus enabling storing flags to files. The other set takes an
 *   argument of type <tt>vector<bool></tt>. This enables the user to temporarily store
 *   some flags, e.g. if another function needs them, and restore them
 *   afterwards.
 *
 *
 *   <h3>User flags and data</h3>
 *
 *   A triangulation offers one bit per line, quad, etc for user flags.
 *   This field can be
 *   accessed as all other data using iterators. Normally, this user flag is
 *   used if an algorithm walks over all cells and needs information whether
 *   another cell, e.g. a neighbor, has already been processed. It can also
 *   be used to flag the lines subject to constraints in 2D, as for example the
 *   functions in the DoFHandler classes do.
 *
 *   There are two functions, @p save_user_flags and @p load_user_flags which
 *   write and read these flags to and from a stream. Unlike
 *   @p save_refine_flags and @p load_refine_flags, these two functions store
 *   and read the flags of all used lines, quads, etc, not only of the
 *   active ones (well, activity is a concept which really only applies to
 *   cells, not for example to lines in 2D, so the abovementioned generalization
 *   to <em>all</em> lines, quads, etc seems plausible).
 *
 *   If you want to store more specific user flags, you can use the functions
 *   @p save_user_flags_line and @p load_user_flags_line and the generalizations
 *   for quads, etc.
 *
 *   As for the refinement and coarsening flags, there exist two versions of these
 *   functions, one which reads/writes from a stream and one which does so from
 *   a <tt>vector<bool></tt>. The latter is used to store flags temporarily, while the
 *   first is used to store them in a file.
 *
 *   It is convention to clear the user flags using the
 *   <tt>Triangulation<>::clear_user_flags()</tt> function before usage, since it is
 *   often necessary to use the flags in more than one function consecutively and
 *   is then error prone to dedicate one of these to clear the flags.
 *
 *   It is recommended that a functions using the flags states so in
 *   its documentation. For example, the
 *   execute_coarsening_and_refinement() function uses some of the user
 *   flags and will therefore destroy any content stored in them.
 *
 *   There is another set of user data, which can be either an
 *   <tt>unsigned int</tt> or a <tt>void *</tt>, for
 *   each line, quad, etc. You can access these through
 *   the functions listed under <tt>User data</tt> in
 *   the accessor classes. These pointers are not used nor changed in
 *   many places of the library, and those classes and functions that
 *   do use them should document this clearly; the most prominent user
 *   of these pointers is the SolutionTransfer class which uses the
 *   cell->user_pointers.
 *
 *   The value of these user indices or pointers is @p NULL by default. Note that
 *   the pointers are not inherited to children upon
 *   refinement. Still, after a remeshing they are available on all
 *   cells, where they were set on the previous mesh (unless, of
 *   course, overwritten by the SolutionTransfer class).
 *
 *   The usual warning about the missing type safety of @p void pointers are
 *   obviously in place here; responsibility for correctness of types etc
 *   lies entirely with the user of the pointer.
 *
 *   @note User pointers and user indices are stored in the same
 *   place. In order to avoid unwanted conversions, Triangulation
 *   checks which one of them is in use and does not allow access to
 *   the other one, until clear_user_data() has been called.
 *
 *   Just like the user flags, this field is not available for vertices,
 *   which does no harm since the vertices have a unique and continuous number
 *   unlike the structured objects lines and quads.
 *
 *
 *   <h3>Boundary approximation</h3>
 *
 *   You can specify a boundary function for each boundary
 *   component. If a new vertex is created on a side or face at the
 *   boundary, this function is used to compute where it will be
 *   placed. The boundary indicator of the face will be used to
 *   determine the proper component. See Boundary for the
 *   details. Usage with the Triangulation object is then like this
 *   (let @p Ball be a class derived from Boundary<tt><2></tt>):
 *
 *   @verbatim
 *     void main () {
 *       Triangulation<2> tria;
 *                                        // set the boundary function
 *                                        // for all boundaries with
 *                                        // boundary indicator 0
 *       Ball ball;
 *       tria.set_boundary (0, ball);
 *
 *       // read some coarse grid
 *
 *
 *       Triangulation<2>::active_cell_iterator cell, endc;
 *       for (int i=0; i<8; ++i)
 *         {
 *           cell = tria.begin_active();
 *           endc = tria.end();
 *
 *                                            // refine all
 *                                            // boundary cells
 *           for (; cell!=endc; ++cell)
 *             if (cell->at_boundary())
 *               cell->set_refine_flag();
 *
 *           tria.execute_coarsening_and_refinement();
 *         };
 *     };
 *   @endverbatim
 *
 *   You should take note of one caveat: if you have concave
 *   boundaries, you must make sure that a new boundary vertex does
 *   not lie too much inside the cell which is to be refined. The
 *   reason is that the center vertex is placed at the point which is
 *   the arithmetic mean of the vertices of the original cell.
 *   Therefore if your new boundary vertex is too near the center of
 *   the old quadrilateral or hexahedron, the distance to the midpoint
 *   vertex will become too small, thus generating distorted
 *   cells. This issue is discussed extensively in
 *   @ref GlossDistorted "distorted cells".
 *
 *
 *   <h3>Technical details</h3>
 *
 *   <h4>Algorithms for mesh regularization and smoothing upon refinement</h4>
 *
 *   We chose an inductive point of view: since upon creation of the
 *   triangulation all cells are on the same level, all regularity assumptions
 *   regarding the maximum difference in level of cells sharing a common face,
 *   edge or vertex hold. Since we use the regularization and smoothing in
 *   each step of the mesh history, when coming to the point of refining it
 *   further the assumptions also hold.
 *
 *   The regularization and smoothing is done in the
 *   @p prepare_coarsening_and_refinement function, which is called by
 *   @p execute_coarsening_and_refinement at the very beginning.  It
 *   decides which additional cells to flag for refinement by looking
 *   at the old grid and the refinement flags for each cell.
 *
 *   <ul>
 *   <li> <em>Regularization:</em> The algorithm walks over all cells checking
 *     whether the present cell is flagged for refinement and a neighbor of the
 *     present cell is refined once less than the present one. If so, flag the
 *     neighbor for refinement. Because of the induction above, there may be no
 *     neighbor with level two less than the present one.
 *
 *     The neighbor thus flagged for refinement may induce more cells which need
 *     to be refined. However, such cells which need additional refinement always
 *     are on one level lower than the present one, so we can get away with only
 *     one sweep over all cells if we do the loop in the reverse way, starting
 *     with those on the highest level. This way, we may flag additional cells
 *     on lower levels, but if these induce more refinement needed, this is
 *     performed later on when we visit them in out backward running loop.
 *
 *   <li> <em>Smoothing:</em>
 *     <ul>
 *     <li> @p limit_level_difference_at_vertices:
 *       First a list is set up which stores for each vertex
 *       the highest level one of the adjacent cells belongs to. Now, since we did
 *       smoothing in the previous refinement steps also, each cell may only have
 *       vertices with levels at most one greater than the level of the present
 *       cell.
 *
 *       However, if we store the level plus one for cells marked for refinement,
 *       we may end up with cells which have vertices of level two greater than
 *       the cells level. We need to refine this cell also, and need thus also
 *       update the levels of its vertices. This itself may lead to cells needing
 *       refinement, but these are on lower levels, as above, which is why we
 *       may do all kinds of additional flagging in one loop only.
 *
 *     <li> @p eliminate_unrefined_islands:
 *       For each cell we count the number of neighbors which are refined or
 *       flagged for refinement. If this exceeds the number of neighbors
 *       which are not refined and not flagged for refinement, then the current
 *       cell is flagged for
 *       refinement. Since this may lead to cells on the same level which also
 *       will need refinement, we will need additional loops of regularization
 *       and smoothing over all cells until nothing changes any more.
 *
 *     <li> <tt>eliminate_refined_*_islands</tt>:
 *       This one does much the same as the above one, but for coarsening. If
 *       a cell is flagged for refinement or if all of its children are active
 *       and if the number of neighbors which are either active and not flagged
 *       for refinement, or not active but all children flagged for coarsening
 *       equals the total number of neighbors, then this cell's children
 *       are flagged for coarsening or (if this cell was flagged for refinement)
 *       the refine flag is cleared.
 *
 *       For a description of the distinction between the two versions of the
 *       flag see above in the section about mesh smoothing in the general part
 *       of this classes description.
 *
 *       The same applies as above: several loops may be necessary.
 *     </ul>
 *   </ul>
 *
 *   Regularization and smoothing are a bit complementary in that we check
 *   whether we need to set additional refinement flags when being on a cell
 *   flagged for refinement (regularization) or on a cell not flagged for
 *   refinement. This makes readable programming easier.
 *
 *   All the described algorithms apply only for more than one space dimension,
 *   since for one dimension no restrictions apply. It may be necessary to
 *   apply some smoothing for multigrid algorithms, but this has to be decided
 *   upon later.
 *
 *
 *   <h3>Warning</h3>
 *
 *   It seems impossible to preserve @p constness of a triangulation through
 *   iterator usage. Thus, if you declare pointers to a @p const triangulation
 *   object, you should be well aware that you might involuntarily alter the
 *   data stored in the triangulation.
 *
 * @ingroup grid aniso
 * @author Wolfgang Bangerth, 1998; Ralf Hartmann, 2005
 */
template <int dim, int spacedim=dim>
class Triangulation : public Subscriptor
{
  private:

				     /**
				      * An internal typedef to make
				      * the definition of the iterator
				      * classes simpler.
				      */
    typedef internal::Triangulation::Iterators<dim, spacedim> IteratorSelector;

  public:
				     /**
				      * Default boundary object. This is used
				      * for those boundaries for which no
				      * boundary object has been explicitly
				      * set using set_boundary().
				      */
    static const StraightBoundary<dim,spacedim> straight_boundary;

				     /**
				      * Default boundary object to be
				      * used in the codimension 1
				      * case.  This is used for those
				      * boundaries for which no
				      * boundary object has been
				      * explicitly set using
				      * set_boundary().
				      */
    static const StraightBoundary<spacedim,spacedim> straight_manifold_description;
				     /**
				      * Declare some symbolic names
				      * for mesh smoothing
				      * algorithms. The meaning of
				      * these flags is documented in
				      * the Triangulation class.
				      */
    enum MeshSmoothing
    {
	  none                               = 0x0,
	  limit_level_difference_at_vertices = 0x1,
	  eliminate_unrefined_islands        = 0x2,
	  patch_level_1                      = 0x4,
	  coarsest_level_1                   = 0x8,

	  allow_anisotropic_smoothing        = 0x10,

	  eliminate_refined_inner_islands    = 0x100,
	  eliminate_refined_boundary_islands = 0x200,
	  do_not_produce_unrefined_islands   = 0x400,

	  smoothing_on_refinement            = (limit_level_difference_at_vertices |
						eliminate_unrefined_islands),
	  smoothing_on_coarsening            = (eliminate_refined_inner_islands |
						eliminate_refined_boundary_islands |
						do_not_produce_unrefined_islands),

	  maximum_smoothing                  = 0xffff ^ allow_anisotropic_smoothing
    };

    typedef typename IteratorSelector::CellAccessor         cell_accessor;
    typedef typename IteratorSelector::FaceAccessor         face_accessor;

    typedef typename IteratorSelector::raw_line_iterator    raw_line_iterator;
    typedef typename IteratorSelector::line_iterator        line_iterator;
    typedef typename IteratorSelector::active_line_iterator active_line_iterator;

    typedef typename IteratorSelector::raw_quad_iterator    raw_quad_iterator;
    typedef typename IteratorSelector::quad_iterator        quad_iterator;
    typedef typename IteratorSelector::active_quad_iterator active_quad_iterator;

    typedef typename IteratorSelector::raw_hex_iterator    raw_hex_iterator;
    typedef typename IteratorSelector::hex_iterator        hex_iterator;
    typedef typename IteratorSelector::active_hex_iterator active_hex_iterator;

    typedef typename IteratorSelector::raw_cell_iterator    raw_cell_iterator;
    typedef typename IteratorSelector::cell_iterator        cell_iterator;
    typedef typename IteratorSelector::active_cell_iterator active_cell_iterator;

    typedef typename IteratorSelector::raw_face_iterator    raw_face_iterator;
    typedef typename IteratorSelector::face_iterator        face_iterator;
    typedef typename IteratorSelector::active_face_iterator active_face_iterator;

				     /**
				      *  Base class for refinement listeners.
				      *  Other classes, which need to be
				      *  informed about refinements of the
				      *  Triangulation,
				      *  can be derived from
				      *  RefinementListener.
				      *  If these classes add theirself to
				      *  the refinement listener list by
				      *  calling add_refinement_listener(),
				      *  they will then be informed about
				      *  the beginning of the refinement
				      *  by a call of pre_refinement_notification()
				      *  and about the end of the refinement
				      *  process by a call of
				      *  post_refinement_notification().
				      *  Both methods are called in the method
				      *  execute_coarsening_and_refinement().
				      *  For an example see hp::DoFHandler().
				      */
    class RefinementListener
    {
      public:
                                         /**
                                          * Destructor. Does nothing, but is
                                          * declared virtual because this
                                          * class also has virtual functions.
                                          */
        virtual ~RefinementListener ();

                                         /**
                                          * Before refinement is actually
                                          * performed, the triangulation class
                                          * calls this method on all objects
                                          * derived from this class and
                                          * registered with the triangulation.
                                          */
        virtual
        void
        pre_refinement_notification (const Triangulation<dim, spacedim> &tria);

                                         /**
                                          * After refinement is actually
                                          * performed, the triangulation class
                                          * calls this method on all objects
                                          * derived from this class and
                                          * registered with the triangulation.
                                          */
        virtual
        void
        post_refinement_notification (const Triangulation<dim, spacedim> &tria);

                                         /**
                                          * At the end of a call to
                                          * copy_triangulation() the
                                          * Triangulation class calls this
                                          * method on all objects derived from
                                          * this class and registered with the
                                          * original Triangulation @p old_tria
                                          * so that they might subscribe to the
                                          * copied one @p new_tria as well, if
                                          * that is desired. By default this
                                          * method does nothing, a different
                                          * behavior has to be implemented in
                                          * derived classes.
                                          */
        virtual
        void
        copy_notification (const Triangulation<dim, spacedim> &old_tria,
			   const Triangulation<dim, spacedim> &new_tria);
    };

				     /**
				      * A structure that is used as an
				      * exception object by the
				      * create_triangulation() function to
				      * indicate which cells among the coarse
				      * mesh cells are inverted or severely
				      * distorted (see the entry on
				      * @ref GlossDistorted "distorted cells"
				      * in the glossary).
				      *
				      * Objects of this kind are
				      * thrown by the
				      * create_triangulation() and
				      * execute_coarsening_and_refinement()
				      * functions, and they can be
				      * caught in user code if this
				      * condition is to be
				      * ignored. Note, however, that
				      * such exceptions are only
				      * produced if the necessity for
				      * this check was indicated when
				      * calling the constructor of the
				      * Triangulation class.
				      *
				      * A cell is called
				      * <i>deformed</i> if the
				      * determinant of the Jacobian of
				      * the mapping from reference
				      * cell to real cell is negative
				      * at least at one vertex. This
				      * computation is done using the
				      * GeometryInfo::jacobian_determinants_at_vertices
				      * function.
				      */
    struct DistortedCellList : public dealii::ExceptionBase
    {
					 /**
					  * Destructor. Empty, but needed
					  * for the sake of exception
					  * specification, since the base
					  * class has this exception
					  * specification and the
					  * automatically generated
					  * destructor would have a
					  * different one due to member
					  * objects.
					  */
	virtual ~DistortedCellList () throw();

					 /**
					  * A list of those cells
					  * among the coarse mesh
					  * cells that are deformed or
					  * whose children are
					  * deformed.
					  */
	std::list<typename Triangulation<dim,spacedim>::cell_iterator>
	distorted_cells;
    };


				     /**
				      * Make the dimension available
				      * in function templates.
				      */
    static const unsigned int dimension = dim;

				     /**
				      * Make the space-dimension available
				      * in function templates.
				      */
    static const unsigned int space_dimension = spacedim;

				     /**
				      *  Create an empty
				      *  triangulation. Do not create
				      *  any cells.
				      *
				      * @param smooth_grid Determines
				      * the level of smoothness of the
				      * mesh size function that should
				      * be enforced upon mesh
				      * refinement.
				      *
				      * @param
				      * check_for_distorted_cells
				      * Determines whether the
				      * triangulation should check
				      * whether any of the cells that
				      * are created by
				      * create_triangulation() or
				      * execute_coarsening_and_refinement()
				      * are distorted (see
				      * @ref GlossDistorted "distorted cells").
				      * If set, these two
				      * functions may throw an
				      * exception if they encounter
				      * distorted cells.
				      */
    Triangulation (const MeshSmoothing smooth_grid = none,
		   const bool check_for_distorted_cells = false);

				     /**
				      *  Copy constructor.
				      *
				      *  You should really use the @p
				      *  copy_triangulation function,
				      *  so we declare this function
				      *  but let it throw an internal
				      *  error. The reason for this is
				      *  that we may want to use
				      *  triangulation objects in
				      *  collections. However, C++
				      *  containers require that the
				      *  objects stored in them are
				      *  copyable, so we need to
				      *  provide a copy
				      *  constructor. On the other
				      *  hand, copying triangulations
				      *  is so expensive that we do
				      *  not want such objects copied
				      *  by accident, for example in
				      *  compiler-generated temporary
				      *  objects. By defining a copy
				      *  constructor but throwing an
				      *  error, we satisfy the formal
				      *  requirements of containers,
				      *  but at the same time disallow
				      *  actual copies. Finally,
				      *  through the exception, one
				      *  easily finds the places where
				      *  code has to be changed to
				      *  avoid copies.
				      */
    Triangulation (const Triangulation<dim, spacedim> &t);

				     /**
				      *  Delete the object and all levels of
				      *  the hierarchy.
				      */
    virtual ~Triangulation ();

				     /**
				      * Reset this triangulation into a
				      * virgin state by deleting all data.
				      *
				      * Note that this operation is only allowed
				      * if no subscriptions to this object exist
				      * any more, such as DoFHandler objects
				      * using it.
				      */
    void clear ();

				     /**
				      * Sets the mesh smoothing to @p
				      * mesh_smoothing. This overrides
				      * the MeshSmoothing given to the
				      * constructor. It is allow to
				      * call this function only if
				      * triangulation is empty.
				      */
    void set_mesh_smoothing(const MeshSmoothing mesh_smoothing);

				     /**
				      * If @p dim==spacedim, assign a boundary
				      * object to a certain part of the
				      * boundary of a the triangulation. If a
				      * face with boundary number @p number is
				      * refined, this object is used to find
				      * the location of new vertices on the
				      * boundary. It is also used for for
				      * non-linear (i.e.: non-Q1)
				      * transformations of cells to the unit
				      * cell in shape function calculations.
				      *
				      * If @p dim!=spacedim the boundary object
				      * is in fact the exact manifold that the
				      * triangulation is approximating (for
				      * example a circle approximated by a
				      * polygon triangulation). As above, the
				      * refinement is made in such a way that
				      * the new points are located on the exact
				      * manifold.
				      *
				      * Numbers of boundary objects correspond
				      * to material numbers of faces at the
				      * boundary, for instance the material id
				      * in a UCD input file. They are not
				      * necessarily consecutive but must be in
				      * the range 0-254.  Material IDs on
				      * boundaries are also called boundary
				      * indicators and are accessed with
				      * accessor functions of that name.
				      *
				      * The @p boundary_object is not copied
				      * and MUST persist until the
				      * triangulation is destroyed. This is
				      * also true for triangulations generated
				      * from this one by @p
				      * copy_triangulation.
				      *
				      * It is possible to remove or replace
				      * the boundary object during the
				      * lifetime of a non-empty
				      * triangulation. Usually, this is done
				      * before the first refinement and is
				      * dangerous afterwards. Removal of a
				      * boundary object is done by
				      * <tt>set_boundary(number)</tt>,
				      * i.e. the function of same name but
				      * only one argument. This operation then
				      * replaces the boundary object given
				      * before by a straight boundary
				      * approximation.
				      */
    void set_boundary (const unsigned int   number,
		       const Boundary<dim,spacedim> &boundary_object);

                                     /**
                                      * Reset those parts of the boundary with
                                      * the given number to use a straight
                                      * boundary approximation. This is the
                                      * default state of a triangulation, and
                                      * undoes assignment of a different
                                      * boundary object by the function of
                                      * same name and two arguments.
                                      */
    void set_boundary (const unsigned int number);

				     /**
				      * Return a constant reference to
				      * a boundary object used for
				      * this triangulation.  Number is
				      * the same as in
				      * @p set_boundary
				      */
    const Boundary<dim,spacedim> & get_boundary (const unsigned int number) const;

				     /**
				      * Returns a vector containing
				      * all boundary indicators
				      * assigned to boundary faces of
				      * this Triangulation
				      * object. Note, that each
				      * boundary indicator is reported
				      * only once. The size of the
				      * return vector will represent
				      * the number of different
				      * indicators (which is greater
				      * or equal one).
				      */
    std::vector<unsigned char> get_boundary_indicators() const;

				     /**
				      *  Copy a triangulation. This
				      *  operation is not cheap, so
				      *  you should be careful with
				      *  using this. We do not
				      *  implement this function as a
				      *  copy constructor, since it
				      *  makes it easier to maintain
				      *  collections of triangulations
				      *  if you can assign them values
				      *  later on.
				      *
				      *  Keep in mind that this
				      *  function also copies the
				      *  pointer to the boundary
				      *  descriptor previously set by
				      *  the @p set_boundary
				      *  function. You must therefore
				      *  also guarantee that the
				      *  boundary objects has a
				      *  lifetime at least as long as
				      *  the copied triangulation.
				      *
				      *  This triangulation must be
				      *  empty beforehand.
				      *
				      *  The function is made
				      *  @p virtual since some
				      *  derived classes might want to
				      *  disable or extend the
				      *  functionality of this
				      *  function.
				      *
				      *  Note, that the list of
				      *  RefinementListeners is not
				      *  copied. However, each
				      *  RefinementListener is notified of the
				      *  copy operation through the
				      *  RefinementListener::copy_notification()
				      *  function, so it might subscribe to the
				      *  new Triangulation as well, if that is
				      *  desired.
				      */
    virtual void copy_triangulation (const Triangulation<dim, spacedim> &old_tria);

				     /**
				      * Create a triangulation from a
				      * list of vertices and a list of
				      * cells, each of the latter
				      * being a list of <tt>1<<dim</tt>
				      * vertex indices. The
				      * triangulation must be empty
				      * upon calling this function and
				      * the cell list should be useful
				      * (connected domain, etc.).
				      *
				      * Material data for the cells is
				      * given within the @p cells
				      * array, while boundary
				      * information is given in the
				      * @p subcelldata field.
				      *
				      * The numbering of vertices
				      * within the @p cells array is
				      * subject to some constraints;
				      * see the general class
				      * documentation for this.
				      *
				      * For conditions when this
				      * function can generate a valid
				      * triangulation, see the
				      * documentation of this class,
				      * and the GridIn and
				      * GridReordering class.
				      *
				      * If the
				      * <code>check_for_distorted_cells</code>
				      * flag was specified upon
				      * creation of this object, at
				      * the very end of its operation,
				      * the current function walks
				      * over all cells and verifies
				      * that none of the cells is
				      * deformed (see the entry on
				      * @ref GlossDistorted "distorted cells"
				      in the glossary), where
				      * we call a cell deformed if the
				      * determinant of the Jacobian of
				      * the mapping from reference
				      * cell to real cell is negative
				      * at least at one of the
				      * vertices (this computation is
				      * done using the
				      * GeometryInfo::jacobian_determinants_at_vertices
				      * function). If there are
				      * deformed cells, this function
				      * throws an exception of kind
				      * DistortedCellList. Since this
				      * happens after all data
				      * structures have been set up,
				      * you can catch and ignore this
				      * exception if you know what you
				      * do -- for example, it may be
				      * that the determinant is zero
				      * (indicating that you have
				      * collapsed edges in a cell) but
				      * that this is ok because you
				      * didn't intend to integrate on
				      * this cell anyway. On the other
				      * hand, deformed cells are often
				      * a sign of a mesh that is too
				      * coarse to resolve the geometry
				      * of the domain, and in this
				      * case ignoring the exception is
				      * probably unwise.
				      *
				      * @note: The check for distorted
				      * cells is only done if
				      * dim==spacedim, as otherwise
				      * cells can legitimately be
				      * twisted if the manifold they
				      * describe is twisted.
				      */
    virtual void create_triangulation (const std::vector<Point<spacedim> >    &vertices,
				       const std::vector<CellData<dim> > &cells,
				       const SubCellData                 &subcelldata);

				     /**
				      * For backward compatibility,
				      * only. This function takes the
				      * cell data in the ordering as
				      * requested by deal.II versions
				      * up to 5.2, converts it to the
				      * new (lexicographic) ordering
				      * and calls
				      * create_triangulation().
				      *
				      * @note This function internally
				      * calls create_triangulation and
				      * therefore can throw the same
				      * exception as the other
				      * function.
				      */
    virtual void create_triangulation_compatibility (
      const std::vector<Point<spacedim> >    &vertices,
      const std::vector<CellData<dim> > &cells,
      const SubCellData                 &subcelldata);

				     /**
				      * Distort the grid by randomly
				      * moving around all the vertices
				      * of the grid.  The direction of
				      * moving is random, while the
				      * length of the shift vector has
				      * a value of @p factor times
				      * the minimal length of the
				      * active lines adjacent to this
				      * vertex. Note that @p factor
				      * should obviously be well below
				      * <tt>0.5</tt>.
				      *
				      * If @p keep_boundary is set to
				      * @p true (which is the
				      * default), then boundary
				      * vertices are not moved.
				      */
    void distort_random (const double factor,
			 const bool   keep_boundary=true);


				     /**
				      *  @name Mesh refinement
				      */
				     /*@{*/
				     /**
				      *  Flag all active cells for
				      *  refinement.  This will refine
				      *  all cells of all levels which
				      *  are not already refined
				      *  (i.e. only cells are refined
				      *  which do not yet have
				      *  children). The cells are only
				      *  flagged, not refined, thus
				      *  you have the chance to save
				      *  the refinement flags.
				      */
    void set_all_refine_flags ();

				     /**
				      * Refine all cells @p times times, by
				      * alternatingly calling
				      * set_all_refine_flags and
				      * execute_coarsening_and_refinement.
				      *
				      * The latter function may throw an
				      * exception if it creates cells that are
				      * distorted (see its documentation for
				      * an explanation). This exception will
				      * be propagated through this function if
				      * that happens, and you may not get the
				      * actual number of refinement steps in
				      * that case.
				      */
    void refine_global (const unsigned int times);

				     /**
				      * Execute both refinement and
				      * coarsening of the
				      * triangulation.
				      *
				      * The function resets all refinement and
				      * coarsening flags to false. It uses the
				      * user flags for internal purposes. They
				      * will therefore be overwritten by
				      * undefined content.
				      *
				      * If the boundary description is
				      * sufficiently irregular, it can
				      * happen that some of the
				      * children produced by mesh
				      * refinement are distorted (see
				      * the extensive discussion on
				      * @ref GlossDistorted "distorted cells").
				      *
				      * To allow user programs to fix
				      * up these cells if that is
				      * desired, this function after
				      * completing all other work may
				      * throw an exception of type
				      * DistortedCellList that
				      * contains a list of those cells
				      * that have been refined and
				      * have at least one child that
				      * is distorted. The function
				      * does not create such an
				      * exception if no cells have
				      * created distorted children.
				      * Note that for the check for
				      * distorted cells to happen, the
				      * <code>check_for_distorted_cells</code>
				      * flag has to be specified upon
				      * creation of a triangulation
				      * object.
				      *
                                      * See the general docs for more
                                      * information.
				      *
				      * Note that this function is
				      * <tt>virtual</tt> to allow
				      * derived classes to insert
				      * hooks, such as saving
				      * refinement flags and the like
				      * (see e.g. the
				      * PersistentTriangulation
				      * class).
				      */
    virtual void execute_coarsening_and_refinement ();

				     /**
				      * Do both preparation for
				      * refinement and coarsening as
				      * well as mesh smoothing.
				      *
				      * Regarding the refinement
				      * process it fixes the closure
				      * of the refinement in
				      * <tt>dim>=2</tt> (make sure that no
				      * two cells are adjacent with a
				      * refinement level differing
				      * with more than one), etc.  It
				      * performs some mesh smoothing
				      * if the according flag was
				      * given to the constructor of
				      * this class.  The function
				      * returns whether additional
				      * cells have been flagged for
				      * refinement.
				      *
				      * See the general doc of this
				      * class for more information on
				      * smoothing upon refinement.
				      *
				      * Regarding the coarsening part,
				      * flagging and deflagging cells
				      * in preparation of the actual
				      * coarsening step are done. This
				      * includes deleting coarsen
				      * flags from cells which may not
				      * be deleted (e.g. because one
				      * neighbor is more refined
				      * than the cell), doing some
				      * smoothing, etc.
				      *
				      * The effect is that only those
				      * cells are flagged for
				      * coarsening which will actually
				      * be coarsened. This includes
				      * the fact that all flagged
				      * cells belong to parent cells
				      * of which all children are
				      * flagged.
				      *
				      * The function returns whether
				      * some cells' flagging has been
				      * changed in the process.
				      *
				      * This function uses the user
				      * flags, so store them if you
				      * still need them afterwards.
				      */
    bool prepare_coarsening_and_refinement ();

				     /**
				      *  Add a
				      *  RefinementListener. Adding
				      *  listeners to the
				      *  Triangulation allows other
				      *  classes to be informed, when
				      *  the Triangulation is refined.
				      */
    void add_refinement_listener (RefinementListener &listener) const;

				     /**
				      *  Removes a
				      *  RefinementListener. When some
				      *  class needs no longer to be
				      *  informed about refinements,
				      *  the listener should be
				      *  removed from the
				      *  Triangulation.
				      */
    void remove_refinement_listener (RefinementListener &listener) const;

				     /*@}*/

				     /**
				      *  @name History of a triangulation
				      */
    				     /*@{*/
				     /**
				      *  Save the addresses of the
				      *  cells which are flagged for
				      *  refinement to @p out.  For
				      *  usage, read the general
				      *  documentation for this class.
				      */
    void save_refine_flags (std::ostream &out) const;

				     /**
				      * Same as above, but store the flags to
				      * a bitvector rather than to a file.
				      */
    void save_refine_flags (std::vector<bool> &v) const;

				     /**
				      *  Read the information stored by
				      *  @p save_refine_flags.
				      */
    void load_refine_flags (std::istream &in);

    				     /**
				      *  Read the information stored by
				      *  @p save_refine_flags.
				      */
    void load_refine_flags (const std::vector<bool> &v);

				     /**
				      * Analogue to @p save_refine_flags.
				      */
    void save_coarsen_flags (std::ostream &out) const;

				     /**
				      * Same as above, but store the flags to
				      * a bitvector rather than to a file.
				      */
    void save_coarsen_flags (std::vector<bool> &v) const;

    				     /**
				      * Analogue to @p load_refine_flags.
				      */
    void load_coarsen_flags (std::istream &out);

				     /**
				      * Analogue to @p load_refine_flags.
				      */
    void load_coarsen_flags (const std::vector<bool> &v);

				     /**
				      * Return whether this triangulation has
				      * ever undergone anisotropic (as opposed
				      * to only isotropic) refinement.
				      */
    bool get_anisotropic_refinement_flag() const;
    				     /*@}*/


				     /**
				      *  @name User data
				      */
				     /*@{*/
    				     /**
				      *  Clear all user flags.
				      */
    void clear_user_flags ();

				     /**
				      *  Save all user flags. See the general
				      *  documentation for this class
				      *  and the documentation for the
				      *  @p save_refine_flags for more
				      *  details.
				      */
    void save_user_flags (std::ostream &out) const;

				     /**
				      * Same as above, but store the flags to
				      * a bitvector rather than to a file.
				      * The output vector is resized if
				      * necessary.
				      */
    void save_user_flags (std::vector<bool> &v) const;

				     /**
				      *  Read the information stored by
				      *  @p save_user_flags.
				      */
    void load_user_flags (std::istream &in);

    				     /**
				      *  Read the information stored by
				      *  @p save_user_flags.
				      */
    void load_user_flags (const std::vector<bool> &v);

    				     /**
				      *  Clear all user flags on lines.
				      */
    void clear_user_flags_line ();

    				     /**
				      * Save the user flags on lines.
				      */
    void save_user_flags_line (std::ostream &out) const;

				     /**
				      * Same as above, but store the flags to
				      * a bitvector rather than to a file.
				      * The output vector is resized if
				      * necessary.
				      */
    void save_user_flags_line (std::vector<bool> &v) const;

				     /**
				      * Load the user flags located on lines.
				      */
    void load_user_flags_line (std::istream &in);

    				     /**
				      * Load the user flags located on lines.
				      */
    void load_user_flags_line (const std::vector<bool> &v);

    				     /**
				      *  Clear all user flags on quads.
				      */
    void clear_user_flags_quad ();

    				     /**
				      * Save the user flags on quads.
				      */
    void save_user_flags_quad (std::ostream &out) const;

				     /**
				      * Same as above, but store the flags to
				      * a bitvector rather than to a file.
				      * The output vector is resized if
				      * necessary.
				      */
    void save_user_flags_quad (std::vector<bool> &v) const;

				     /**
				      * Load the user flags located on quads.
				      */
    void load_user_flags_quad (std::istream &in);

				     /**
				      * Load the user flags located on quads.
				      */
    void load_user_flags_quad (const std::vector<bool> &v);


    				     /**
				      *  Clear all user flags on quads.
				      */
    void clear_user_flags_hex ();

    				     /**
				      * Save the user flags on hexs.
				      */
    void save_user_flags_hex (std::ostream &out) const;

				     /**
				      * Same as above, but store the flags to
				      * a bitvector rather than to a file.
				      * The output vector is resized if
				      * necessary.
				      */
    void save_user_flags_hex (std::vector<bool> &v) const;

				     /**
				      * Load the user flags located on hexs.
				      */
    void load_user_flags_hex (std::istream &in);

				     /**
				      * Load the user flags located on hexs.
				      */
    void load_user_flags_hex (const std::vector<bool> &v);

				     /**
				      * Clear all user pointers and
				      * indices and allow the use of
				      * both for next access.
				      */
    void clear_user_data ();

				     /**
				      * @deprecated User
				      * clear_user_data() instead.
				      *
				      *  Clear all user pointers.
				      */
    void clear_user_pointers ();

				     /**
				      * Save all user indices. The
				      * output vector is resized if
				      * necessary.
				      */
    void save_user_indices (std::vector<unsigned int> &v) const;

    				     /**
				      * Read the information stored by
				      * save_user_indices().
				      */
    void load_user_indices (const std::vector<unsigned int> &v);

				     /**
				      * Save all user pointers. The
				      * output vector is resized if
				      * necessary.
				      */
    void save_user_pointers (std::vector<void *> &v) const;

    				     /**
				      * Read the information stored by
				      * save_user_pointers().
				      */
    void load_user_pointers (const std::vector<void *> &v);

				     /**
				      * Save the user indices on
				      * lines. The output vector is
				      * resized if necessary.
				      */
    void save_user_indices_line (std::vector<unsigned int> &v) const;

    				     /**
				      * Load the user indices located
				      * on lines.
				      */
    void load_user_indices_line (const std::vector<unsigned int> &v);

				     /**
				      * Save the user indices on
				      * quads. The output vector is
				      * resized if necessary.
				      */
    void save_user_indices_quad (std::vector<unsigned int> &v) const;

				     /**
				      * Load the user indices located
				      * on quads.
				      */
    void load_user_indices_quad (const std::vector<unsigned int> &v);

				     /**
				      * Save the user indices on
				      * hexes. The output vector is
				      * resized if necessary.
				      */
    void save_user_indices_hex (std::vector<unsigned int> &v) const;

				     /**
				      * Load the user indices located
				      * on hexs.
				      */
    void load_user_indices_hex (const std::vector<unsigned int> &v);
				     /**
				      * Save the user indices on
				      * lines. The output vector is
				      * resized if necessary.
				      */
    void save_user_pointers_line (std::vector<void *> &v) const;

    				     /**
				      * Load the user pointers located
				      * on lines.
				      */
    void load_user_pointers_line (const std::vector<void *> &v);

				     /**
				      * Save the user pointers on
				      * quads. The output vector is
				      * resized if necessary.
				      */
    void save_user_pointers_quad (std::vector<void *> &v) const;

				     /**
				      * Load the user pointers located
				      * on quads.
				      */
    void load_user_pointers_quad (const std::vector<void *> &v);

				     /**
				      * Save the user pointers on
				      * hexes. The output vector is
				      * resized if necessary.
				      */
    void save_user_pointers_hex (std::vector<void *> &v) const;

				     /**
				      * Load the user pointers located
				      * on hexs.
				      */
    void load_user_pointers_hex (const std::vector<void *> &v);
				     /*@}*/


				     /**
				      *  @name Cell iterator functions
				      */
				     /*@{*/
				     /**
				      *  Iterator to the first cell, used
				      *  or not, on level @p level. If a level
				      *  has no cells, a past-the-end iterator
				      *  is returned.
				      *
				      *  This function calls begin_raw_line()
				      *  in 1D and begin_raw_quad() in 2D.
				      */
    raw_cell_iterator    begin_raw   (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first used cell
				      *  on level @p level.
				      *
				      *  This function calls @p begin_line
				      *  in 1D and @p begin_quad in 2D.
				      */
    cell_iterator        begin       (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first active
				      *  cell on level @p level.
				      *
				      *  This function calls @p begin_active_line
				      *  in 1D and @p begin_active_quad in 2D.
				      */
    active_cell_iterator begin_active(const unsigned int level = 0) const;

				     /**
				      *  Iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      *
				      *  This function calls @p end_line
				      *  in 1D and @p end_quad in 2D.
				      */
    raw_cell_iterator    end () const;

				     /**
				      * Return an iterator which is the first
				      * iterator not on level. If @p level is
				      * the last level, then this returns
				      * <tt>end()</tt>.
				      */
    cell_iterator        end (const unsigned int level) const;

				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If @p level is
				      * the last level, then this returns
				      * <tt>end()</tt>.
				      */
    raw_cell_iterator    end_raw (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If @p level
				      * is the last level, then this returns
				      * <tt>end()</tt>.
				      */
    active_cell_iterator end_active (const unsigned int level) const;


				     /**
				      *  Return an iterator pointing to the
				      *  last cell, used or not.
				      *
				      *  This function calls @p last_raw_line
				      *  in 1D and @p last_raw_quad in 2D.
				      */
    raw_cell_iterator    last_raw () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  cell of the level @p level, used or not.
				      *
				      *  This function calls @p last_raw_line
				      *  in 1D and @p last_raw_quad in 2D.
				      */
    raw_cell_iterator    last_raw (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used cell.
				      *
				      *  This function calls @p last_line
				      *  in 1D and @p last_quad in 2D.
				      */
    cell_iterator        last () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used cell on level @p level.
				      *
				      *  This function calls @p last_line
				      *  in 1D and @p last_quad in 2D.
				      */
    cell_iterator        last (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active cell.
				      *
				      *  This function calls @p last_active_line
				      *  in 1D and @p last_active_quad in 2D.
				      */
    active_cell_iterator last_active () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active cell on level @p level.
				      *
				      *  This function calls @p last_active_line
				      *  in 1D and @p last_active_quad in 2D.
				      */
    active_cell_iterator last_active (const unsigned int level) const;
				     /*@}*/

    				     /*---------------------------------------*/
    				     /*---------------------------------------*/

    				     /**
				      *  @name Face iterator functions
				      */
				     /*@{*/
				     /**
				      *  Iterator to the first face, used
				      *  or not. As faces have no level,
				      *  no argument can be given.
				      *
				      *  This function calls @p begin_raw_line
				      *  in 2D and @p begin_raw_quad in 3D.
				      */
    raw_face_iterator    begin_raw_face   () const;

				     /**
				      *  Iterator to the first used face.
				      *
				      *  This function calls @p begin_line
				      *  in 2D and @p begin_quad in 3D.
				      */
    face_iterator        begin_face       () const;

				     /**
				      *  Iterator to the first active
				      *  face.
				      *
				      *  This function calls @p begin_active_line
				      *  in 2D and @p begin_active_quad in 3D.
				      */
    active_face_iterator begin_active_face() const;

				     /**
				      *  Iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      *
				      *  This function calls @p end_line
				      *  in 2D and @p end_quad in 3D.
				      */
    raw_face_iterator    end_face () const;

				     /**
				      * Return a raw iterator which is past the end.
				      * This is the same as <tt>end()</tt> and is
				      * only for combatibility with older versions.
				      */
    raw_face_iterator    end_raw_face () const;

    				     /**
				      * Return an active iterator which is past the end.
				      * This is the same as <tt>end()</tt> and is
				      * only for combatibility with older versions.
				      */
    active_face_iterator end_active_face () const;

				     /**
				      *  Return an iterator pointing to the
				      *  last face, used or not.
				      *
				      *  This function calls @p last_raw_line
				      *  in 2D and @p last_raw_quad in 3D.
				      */
    raw_face_iterator    last_raw_face () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used face.
				      *
				      *  This function calls @p last_line
				      *  in 2D and @p last_quad in 3D.
				      */
    face_iterator        last_face () const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active face.
				      *
				      *  This function calls @p last_active_line
				      *  in 2D and @p last_active_quad in 3D.
				      */
    active_face_iterator last_active_face () const;

				     /*@}*/


				     /*---------------------------------------*/

				     /**
				      *  @name Line iterator functions
				      */
				     /*@{*/
				     /**
				      *  Iterator to the first line, used
				      *  or not, on level @p level. If a level
				      *  has no lines, a past-the-end iterator
				      *  is returned.
				      *  If lines are no cells, i.e. for @p dim>1
				      *  no @p level argument must be given.
				      *  The same applies for all the other functions
				      *  above, of course.
				      */
    raw_line_iterator
    begin_raw_line   (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first used line
				      *  on level @p level.
				      */
    line_iterator
    begin_line       (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first active
				      *  line on level @p level.
				      */
    active_line_iterator
    begin_active_line(const unsigned int level = 0) const;

				     /**
				      *  Iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      */
    raw_line_iterator
    end_line () const;

				     /**
				      * Return an iterator which is the first
				      * iterator not on level. If @p level is
				      * the last level, then this returns
				      * <tt>end()</tt>.
				      */
    line_iterator        end_line (const unsigned int level) const;

				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If @p level is
				      * the last level, then this returns
				      * <tt>end()</tt>.
				      */
    raw_line_iterator    end_raw_line (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If @p level
				      * is the last level, then this returns
				      * <tt>end()</tt>.
				      */
    active_line_iterator end_active_line (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the
				      *  last line, used or not.
				      */
    raw_line_iterator
    last_raw_line () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  line of the level @p level, used or not.
				      */
    raw_line_iterator
    last_raw_line (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used line.
				      */
    line_iterator
    last_line () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used line on level @p level.
				      */
    line_iterator
    last_line (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active line.
				      */
    active_line_iterator
    last_active_line () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active line on level @p level.
				      */
    active_line_iterator
    last_active_line (const unsigned int level) const;
				     /*@}*/

				     /*---------------------------------------*/

				     /**
				      *  @name Quad iterator functions
				      */
    				     /*@{
				      */
    				     /**
				      *  Iterator to the first quad, used
				      *  or not, on the given level. If a level
				      *  has no quads, a past-the-end iterator
				      *  is returned.
				      *  If quads are no cells, i.e. for $dim>2$
				      *  no level argument must be given.

				      */
    raw_quad_iterator
    begin_raw_quad   (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first used quad
				      *  on level @p level.
				      */
    quad_iterator
    begin_quad       (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first active
				      *  quad on level @p level.
				      */
    active_quad_iterator
    begin_active_quad(const unsigned int level = 0) const;

				     /**
				      *  Iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      */
    raw_quad_iterator    end_quad () const;

				     /**
				      * Return an iterator which is the first
				      * iterator not on level. If @p level is
				      * the last level, then this returns
				      * <tt>end()</tt>.
				      */
    quad_iterator        end_quad (const unsigned int level) const;

				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If @p level is
				      * the last level, then this returns
				      * <tt>end()</tt>.
				      */
    raw_quad_iterator    end_raw_quad (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If @p level
				      * is the last level, then this returns
				      * <tt>end()</tt>.
				      */
    active_quad_iterator end_active_quad (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the
				      *  last quad, used or not.
				      */
    raw_quad_iterator
    last_raw_quad () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  quad of the level @p level, used or not.

				      */
    raw_quad_iterator
    last_raw_quad (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used quad.
				      */
    quad_iterator
    last_quad () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used quad on level @p level.
				      */
    quad_iterator
    last_quad (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active quad.
				      */
    active_quad_iterator
    last_active_quad () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active quad on level @p level.
				      */
    active_quad_iterator
    last_active_quad (const unsigned int level) const;
				     /*@}*/

				     /*---------------------------------------*/

				     /**
				      *  @name Hex iterator functions
				      */
    				     /*@{
				      */
    				     /**
				      *  Iterator to the first hex, used
				      *  or not, on level @p level. If a level
				      *  has no hexs, a past-the-end iterator
				      *  is returned.
				      */
    raw_hex_iterator
    begin_raw_hex   (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first used hex
				      *  on level @p level.
				      */
    hex_iterator
    begin_hex       (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first active
				      *  hex on level @p level.
				      */
    active_hex_iterator
    begin_active_hex(const unsigned int level = 0) const;

				     /**
				      *  Iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      */
    raw_hex_iterator
    end_hex () const;

				     /**
				      * Return an iterator which is the first
				      * iterator not on level. If @p level is
				      * the last level, then this returns
				      * <tt>end()</tt>.
				      */
    hex_iterator        end_hex (const unsigned int level) const;

				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If @p level is
				      * the last level, then this returns
				      * <tt>end()</tt>.
				      */
    raw_hex_iterator    end_raw_hex (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If @p level
				      * is the last level, then this returns
				      * <tt>end()</tt>.
				      */
    active_hex_iterator end_active_hex (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the
				      *  last hex, used or not.
				      */
    raw_hex_iterator
    last_raw_hex () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  hex of the level @p level, used or not.

				      */
    raw_hex_iterator
    last_raw_hex (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used hex.
				      */
    hex_iterator
    last_hex () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used hex on level @p level.
				      */
    hex_iterator
    last_hex (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active hex.
				      */
    active_hex_iterator
    last_active_hex () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active hex on level @p level.
				      */
    active_hex_iterator
    last_active_hex (const unsigned int level) const;
				     /*@}*/

				     /*---------------------------------------*/

				     /**
				      * @name Information about the triangulation
				      */
				     /*@{*/

				     /**
				      * In the following, most
				      * functions are provided in two
				      * versions, with and without an
				      * argument describing the
				      * level. The versions with this
				      * argument are only applicable
				      * for objects descibing the
				      * cells of the present
				      * triangulation. For example: in
				      * 2D <tt>n_lines(level)</tt>
				      * cannot be called, only
				      * <tt>n_lines()</tt>, as lines
				      * are faces in 2D and therefore
				      * have no level.
				      */

				     /**
				      * Total Number of lines, used or
				      * unused.
				      */
    unsigned int n_raw_lines () const;

				     /**
				      * Number of lines, used or
				      * unused, on the given level.
				      */
    unsigned int n_raw_lines (const unsigned int level) const;

				     /**
				      *  Return total number of used lines,
				      *  active or not.
				      */
    unsigned int n_lines () const;

				     /**
				      *  Return total number of used lines,
				      *  active or not on level @p level.
				      */
    unsigned int n_lines (const unsigned int level) const;

				     /**
				      * Return total number of active lines.
				      */
    unsigned int n_active_lines () const;

				     /**
				      *  Return total number of active lines,
				      *  on level @p level.
				      */
    unsigned int n_active_lines (const unsigned int level) const;

				     /**
				      * Total number of quads, used or
				      * unused.
				      */
    unsigned int n_raw_quads () const;

				     /**
				      * Number of quads, used or
				      * unused, on the given level.
				      */
    unsigned int n_raw_quads (const unsigned int level) const;

				     /**
				      *  Return total number of used quads,
				      *  active or not.
				      */
    unsigned int n_quads () const;

				     /**
				      *  Return total number of used quads,
				      *  active or not on level @p level.
				      */
    unsigned int n_quads (const unsigned int level) const;

				     /**
				      *  Return total number of active quads,
				      *  active or not.
				      */
    unsigned int n_active_quads () const;

				     /**
				      *  Return total number of active quads,
				      *  active or not on level @p level.
				      */
    unsigned int n_active_quads (const unsigned int level) const;

				     /**
				      * Total number of hexs, used or
				      * unused.
				      */
    unsigned int n_raw_hexs () const;

				     /**
				      * Number of hexs, used or
				      * unused, on the given level.
				      */
    unsigned int n_raw_hexs (const unsigned int level) const;

				     /**
				      *  Return total number of used hexahedra,
				      *  active or not.
				      */
    unsigned int n_hexs() const;

				     /**
				      *  Return total number of used hexahedra,
				      *  active or not on level @p level.
				      */
    unsigned int n_hexs(const unsigned int level) const;

				     /**
				      *  Return total number of active hexahedra,
				      *  active or not.
				      */
    unsigned int n_active_hexs() const;

				     /**
				      *  Return total number of active hexahedra,
				      *  active or not on level @p level.
				      */
    unsigned int n_active_hexs(const unsigned int level) const;

				     /**
				      * Number of cells, used or
				      * unused, on the given level.
				      */
    unsigned int n_raw_cells (const unsigned int level) const;

    				     /**
				      *  Return total number of used cells,
				      *  active or not.
				      *  Maps to <tt>n_lines()</tt> in one space
				      *  dimension and so on.
				      */
    unsigned int n_cells () const;

    				     /**
				      *  Return total number of used cells,
				      *  active or not, on level @p level.
				      *  Maps to <tt>n_lines(level)</tt> in one space
				      *  dimension and so on.
				      */
    unsigned int n_cells (const unsigned int level) const;

    				     /**
				      *  Return total number of active cells.
				      *  Maps to <tt>n_active_lines()</tt> in one space
				      *  dimension and so on.
				      */
    unsigned int n_active_cells () const;

    				     /**
				      * Return total number of active cells
				      * on level @p level.
				      * Maps to <tt>n_active_lines(level)</tt> in one
				      * space dimension and so on.
				      */
    unsigned int n_active_cells (const unsigned int level) const;

				     /**
				      *  Return total number of faces,
				      *  used or not. In 2d, the result
				      *  equals n_raw_lines(), while in 3d it
				      *  equals n_raw_quads().
				      */
    unsigned int n_raw_faces () const;

				     /**
				      *  Return total number of used faces,
				      *  active or not.  In 2D, the result
				      *  equals n_lines(), while in 3D it
				      *  equals n_quads(). Since there are no
				      *  face objects in 1d, the function
				      *  returns zero in 1d.
				      */
    unsigned int n_faces () const;

    				     /**
				      *  Return total number of active faces,
				      *  active or not.  In 2D, the result
				      *  equals n_active_lines(), while in 3D
				      *  it equals n_active_quads(). Since
				      *  there are no face objects in 1d, the
				      *  function returns zero in 1d.
				      */
    unsigned int n_active_faces () const;

				     /**
				      * Return number of levels in use. This
				      * may be less than the number of levels
				      * existing in the triangulation if by
				      * coarsening the highest level was
				      * completely depopulated. That level is
				      * not removed, since it will most likely
				      * be repopulated soon by the next
				      * refinement process.
				      */
    unsigned int n_levels () const;

				     /**
				      * Return the total number of
				      * vertices.  Some of them may
				      * not be used, which usually
				      * happens upon coarsening of a
				      * triangulation when some
				      * vertices are discarded, but we
				      * do not want to renumber the
				      * remaining ones, leading to
				      * holes in the numbers of used
				      * vertices.  You can get the
				      * number of used vertices using
				      * @p n_used_vertices function.
				      */
    unsigned int n_vertices () const;

				     /**
				      * Return a constant reference to
				      * all the vertices present in
				      * this triangulation. Note that
				      * not necessarily all vertices
				      * in this array are actually
				      * used; for example, if you
				      * coarsen a mesh, then some
				      * vertices are deleted, but
				      * their positions in this array
				      * are unchanged as the indices
				      * of vertices are only allocated
				      * once. You can find out about
				      * which vertices are actually
				      * used by the function
				      * get_used_vertices().
				      */
    const std::vector<Point<spacedim> > &
    get_vertices () const;

				     /**
				      * Return the number of vertices
				      * that are presently in use,
				      * i.e. belong to at least one
				      * used element.
				      */
    unsigned int n_used_vertices () const;

				     /**
				      * Return @p true if the vertex
				      * with this @p index is used.
				      */
    bool vertex_used (const unsigned int index) const;

				     /**
				      * Return a constant reference to
				      * the array of @p bools
				      * indicating whether an entry in
				      * the vertex array is used or
				      * not.
				      */
    const std::vector<bool> &
    get_used_vertices () const;

				     /**
				      * Return the maximum number of
				      * cells meeting at a common
				      * vertex. Since this number is
				      * an invariant under refinement,
				      * only the cells on the coarsest
				      * level are considered. The
				      * operation is thus reasonably
				      * fast. The invariance is only
				      * true for sufficiently many
				      * cells in the coarsest
				      * triangulation (e.g. for a
				      * single cell one would be
				      * returned), so a minimum of
				      * four is returned in two
				      * dimensions, 8 in three
				      * dimensions, etc, which is how
				      * many cells meet if the
				      * triangulation is refined.
				      *
				      * In one space dimension, two is
				      * returned.
				      */
    unsigned int max_adjacent_cells () const;
    				     /*@}*/

				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      *
				      * This function is made virtual,
				      * since a triangulation object
				      * might be accessed through a
				      * pointer to this base class,
				      * even if the actual object is a
				      * derived class.
				      */
    virtual unsigned int memory_consumption () const;


				     /**
				      *  @name Exceptions
				      */
				     /*@{*/
				     /**
				      *  Exception
				      * @ingroup Exceptions
				      */
    DeclException1 (ExcInvalidLevel,
		    int,
		    << "The given level " << arg1
		    << " is not in the valid range!");
				     /**
				      * The function raising this
				      * exception can only operate on
				      * an empty Triangulation, i.e.,
				      * a Triangulation without grid
				      * cells.
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcTriangulationNotEmpty);
				     /**
				      * Trying to re-read a grid, an error occured.
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcGridReadError);
				     /**
				      * Exception
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcFacesHaveNoLevel);
				     /**
				      * The triangulation level you
				      * accessed is empty.
				      *
				      * @ingroup Exceptions
				      */
    DeclException1 (ExcEmptyLevel,
		    int,
		    << "You tried to do something on level " << arg1
		    << ", but this level is empty.");
				     /*@}*/
  protected:

				     /**
				      *  Write a bool vector to the given stream,
				      *  writing a pre- and a postfix magic
				      *  number. The vector is written in an
				      *  almost binary format, i.e. the bool
				      *  flags are packed but the data is written
				      *  as ASCII text.
				      *
				      *  The flags are stored in a
				      *  binary format: for each @p
				      *  true, a @p 1 bit is stored, a
				      *  @p 0 bit otherwise.  The bits
				      *  are stored as <tt>unsigned
				      *  char</tt>, thus avoiding
				      *  endianess. They are written
				      *  to @p out in plain text, thus
				      *  amounting to 3.6 bits in the
				      *  output per bits in the input
				      *  on the average. Other
				      *  information (magic numbers
				      *  and number of elements of the
				      *  input vector) is stored as
				      *  plain text as well. The
				      *  format should therefore be
				      *  interplatform compatible.
				      */
    static void write_bool_vector (const unsigned int       magic_number1,
				   const std::vector<bool> &v,
				   const unsigned int       magic_number2,
				   std::ostream            &out);

				     /**
				      * Re-read a vector of bools previously
				      * written by @p write_bool_vector and
				      * compare with the magic numbers.
				      */
    static void read_bool_vector (const unsigned int       magic_number1,
				  std::vector<bool>       &v,
				  const unsigned int       magic_number2,
				  std::istream            &in);

  private:

				     /**
				      * The (public) clear() will only
				      * work when the triangulation is
				      * not subscriped to by other
				      * users. The
				      * clear_despite_subscriptions()
				      * function now allows the
				      * triangulation being cleared
				      * even when there are
				      * subscriptions.
				      *
				      * Make sure, you know what you
				      * do, when calling this
				      * function, as its use is
				      * reasonable in very rare cases,
				      * only. For example, when the
				      * subscriptions were for the
				      * initially empty Triangulation
				      * and the Triangulation object
				      * wants to release its memory
				      * before throwing an assertion
				      * due to input errors (e.g. in
				      * the create_triangulation()
				      * function).
				      */
    void clear_despite_subscriptions ();

				     /**
				      *  Refine all cells on all levels which
				      *  were previously flagged for refinement.
				      *
				      *  Note, that this function uses
				      *  the <tt>line->user_flags</tt>
				      *  for <tt>dim=2,3</tt> and the
				      *  <tt>quad->user_flags</tt> for
				      *  <tt>dim=3</tt>.
				      *
				      *  The function returns a list
				      *  of cells that have produced
				      *  children that satisfy the
				      *  criteria of
				      *  @ref GlossDistorted "distorted cells"
				      * if the
				      * <code>check_for_distorted_cells</code>
				      * flag was specified upon
				      * creation of this object, at
				      */
    DistortedCellList execute_refinement ();

				     /**
				      * Coarsen all cells which were flagged for
				      * coarsening, or rather: delete all
				      * children of those cells of which all
				      * child cells are flagged for coarsening
				      * and several other constraints hold (see
				      * the general doc of this class).
				      */
    void execute_coarsening ();

				     /**
				      * Make sure that either all or none of
				      * the children of a cell are tagged for
				      * coarsening.
				      */
    void fix_coarsen_flags ();

				     /**
				      *  Array of pointers pointing to the
				      *  objects storing the cell data on the
				      *  different levels.
				      */
    std::vector<internal::Triangulation::TriaLevel<dim>*> levels;

				     /**
				      *  Pointer to the faces of the triangulation. In 1d
				      *  this contains nothing, in 2D it contains data
				      *  concerning lines and in 3D quads and lines.  All of
				      *  these have no level and are therefore treated
				      *  seperately.
				      */
    internal::Triangulation::TriaFaces<dim> * faces;


				     /**
				      *  Array of the vertices of this
				      *  triangulation.
				      */
    std::vector<Point<spacedim> >              vertices;

				     /**
				      *  Array storing a bit-pattern which
				      *  vertices are used.
				      */
    std::vector<bool>                     vertices_used;

				     /**
				      *  Collection of boundary
				      *  objects. We only need 255
				      *  objects rather than 256,
				      *  since the indicator 255 is
				      *  reserved for interior faces
				      *  and can thus never be
				      *  associated with a boundary.
				      */
    SmartPointer<const Boundary<dim,spacedim>,Triangulation<dim,spacedim> > boundary[255];

				     /**
				      *  Collection of ly
				      *  objects. We only need 255
				      *  objects rather than 256,
				      *  since the indicator 255 is
				      *  reserved for interior faces
				      *  and can thus never be
				      *  associated with a boundary.
				      */
    SmartPointer<const Boundary<spacedim,spacedim>,Triangulation<dim,spacedim> > manifold_description[255];

				     /**
				      * Flag indicating whether
				      * anisotropic refinement took
				      * place.
				      */
    bool                             anisotropic_refinement;

				     /**
				      *  Do some smoothing in the process
				      *  of refining the triangulation. See
				      *  the general doc of this class for
				      *  more information about this.
				      */
    MeshSmoothing                    smooth_grid;

				     /**
				      * A flag that determines whether
				      * we are to check for distorted
				      * cells upon creation and
				      * refinement of a mesh.
				      */
    const bool check_for_distorted_cells;

				     /**
				      * Cache to hold the numbers of lines,
				      * quads, hexes, etc. These numbers
				      * are set at the end of the refinement
				      * and coarsening functions and enable
				      * faster access later on. In the old
				      * days, whenever one wanted to access
				      * one of these numbers, one had to
				      * perform a loop over all lines, e.g.,
				      * and count the elements until we hit
				      * the end iterator. This is time
				      * consuming and since access to the
				      * number of lines etc is a rather
				      * frequent operation, this was not
				      * an optimal solution.
				      */
    internal::Triangulation::NumberCache<dim> number_cache;

				     /**
				      *  List of RefinementListeners,
				      *  which want to be informed if the
				      *  Triangulation is refined.
				      */
    mutable std::list<RefinementListener *> refinement_listeners;

				     // make a couple of classes
				     // friends
    template <int,int,int> friend class TriaAccessorBase;
    template <int,int,int> friend class TriaAccessor;
    friend class CellAccessor<dim, spacedim>;

    friend class internal::TriaAccessor::Implementation;

    friend class hp::DoFHandler<dim,spacedim>;

    friend class internal::Triangulation::Implementation;
};


#ifndef DOXYGEN

template <int dim, int spacedim>
inline
bool
Triangulation<dim,spacedim>::vertex_used(const unsigned int index) const
{
  Assert (index < vertices_used.size(),
	  ExcIndexRange(index, 0, vertices_used.size()));
  return vertices_used[index];
}



template <int dim, int spacedim>
inline
unsigned int Triangulation<dim, spacedim>::n_levels () const
{
  return number_cache.n_levels;
}



template <int dim, int spacedim>
inline
unsigned int
Triangulation<dim, spacedim>::n_vertices () const
{
  return vertices.size();
}



template <int dim, int spacedim>
inline
const std::vector<Point<spacedim> > &
Triangulation<dim, spacedim>::get_vertices () const
{
  return vertices;
}



/* -------------- declaration of explicit specializations ------------- */

template <> unsigned int Triangulation<1,1>::n_raw_lines (const unsigned int level) const;
template <> unsigned int Triangulation<1,1>::n_quads () const;
template <> unsigned int Triangulation<1,1>::n_quads (const unsigned int level) const;
template <> unsigned int Triangulation<1,1>::n_raw_quads (const unsigned int level) const;
template <> unsigned int Triangulation<2,2>::n_raw_quads (const unsigned int level) const;
template <> unsigned int Triangulation<1,1>::n_raw_hexs (const unsigned int level) const;
template <> unsigned int Triangulation<1,1>::n_active_quads (const unsigned int level) const;
template <> unsigned int Triangulation<1,1>::n_active_quads () const;
template <> unsigned int Triangulation<1,1>::max_adjacent_cells () const;


// -------------------------------------------------------------------
// -- Explicit specializations for codimension one grids


template <> unsigned int Triangulation<1,2>::n_raw_lines (const unsigned int level) const;
template <> unsigned int Triangulation<1,2>::n_quads () const;
template <> unsigned int Triangulation<1,2>::n_quads (const unsigned int level) const;
template <> unsigned int Triangulation<1,2>::n_raw_quads (const unsigned int level) const;
template <> unsigned int Triangulation<2,3>::n_raw_quads (const unsigned int level) const;
template <> unsigned int Triangulation<1,2>::n_raw_hexs (const unsigned int level) const;
template <> unsigned int Triangulation<1,2>::n_active_quads (const unsigned int level) const;
template <> unsigned int Triangulation<1,2>::n_active_quads () const;
template <> unsigned int Triangulation<1,2>::max_adjacent_cells () const;

// -------------------------------------------------------------------
// -- Explicit specializations for codimension two grids


template <> unsigned int Triangulation<1,3>::n_raw_lines (const unsigned int level) const;
template <> unsigned int Triangulation<1,3>::n_quads () const;
template <> unsigned int Triangulation<1,3>::n_quads (const unsigned int level) const;
template <> unsigned int Triangulation<1,3>::n_raw_quads (const unsigned int level) const;
template <> unsigned int Triangulation<2,3>::n_raw_quads (const unsigned int level) const;
template <> unsigned int Triangulation<1,3>::n_raw_hexs (const unsigned int level) const;
template <> unsigned int Triangulation<1,3>::n_active_quads (const unsigned int level) const;
template <> unsigned int Triangulation<1,3>::n_active_quads () const;
template <> unsigned int Triangulation<1,3>::max_adjacent_cells () const;


// -------------------------------------------------------------------




#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

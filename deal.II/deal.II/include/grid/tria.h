//----------------------------  tria.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tria.h  ---------------------------
#ifndef __deal2__tria_h
#define __deal2__tria_h


#include <vector>
#include <base/point.h>
#include <grid/geometry_info.h>
#include <base/subscriptor.h>

template <int dim> class Boundary;
template <int dim> class StraightBoundary;
template <int dim> class CellAccessor;
template<int celldim, int dim> class TriaObjectAccessor;
template <int dim, typename Accessor> class TriaRawIterator;
template <int dim, typename Accessor> class TriaIterator;
template <int dim, typename Accessor> class TriaActiveIterator;
template <int dim> class TriangulationLevel;
template <int dim> class DoFHandler;
template <int dim> class MGDoFHandler;

/*------------------------------------------------------------------------*/


/**
 *  Structure which is passed to the @ref{Triangulation}@p{<dim>::create_triangulation}
 *  function. It contains all data needed to construct a cell, namely the
 *  indices of the vertices and the material indicator.
 */
template <int dim>
struct CellData
{
#if !((__GNUC__==2) && (__GNUC_MINOR__==95))
    int           vertices[GeometryInfo<dim>::vertices_per_cell];
#else
    int           vertices[1 << dim];
#endif
    unsigned char material_id;
};


/**
 *  Structure to be passed to the @ref{Triangulation}@p{<dim>::create_triangulation}
 *  function to describe boundary information.
 *
 *  This structure is the same for all dimensions, since we use an input
 *  function which is the same for all dimensions. The content of objects
 *  of this structure varies with the dimensions, however.
 *
 *  Since in one space dimension, there is no boundary information apart
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
 */
struct SubCellData {
				     /**
				      * Each record of this vector describes
				      * a line on the boundary and its boundary
				      * indicator.
				      */
    vector<CellData<1> > boundary_lines;

				     /**
				      * Each record of this vector describes
				      * a quad on the boundary and its boundary
				      * indicator.
				      */
    vector<CellData<2> > boundary_quads;

				     /**
				      * This function checks whether the vectors
				      * which may not be used in a given
				      * dimension are really empty. I.e.,
				      * whether the @p{boundary_*} arrays are
				      * empty when in one space dimension
				      * and whether the @p{boundary_quads}
				      * array is empty when in two dimensions.
				      *
				      * Since this structure is the same for all
				      * dimensions, the actual dimension has
				      * to be given as a parameter.
				      */
    bool check_consistency (const unsigned int dim) const;
};


/*------------------------------------------------------------------------*/


/**
 *  This class implements some types which differ between the
 *  dimensions.  Declare it to have a template parameter, but do not
 *  actually declare anything concrete apart from the other classes
 *  which are explicitely instantiated ones with the same name.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class TriaDimensionInfo
{};



/**
 *  This class implements some types which differ between the dimensions.
 *  These are the declararions for the 1D case only.
 *
 *  A @p{line_iterator} is typdef'd to an iterator operating on the
 *  @p{lines} member variable of a @p{Triangulation<1>} object. An
 *  @p{active_line_iterator} only operates on the active lines.
 *  @p{raw_line_iterator} objects operate on all lines, used or not.
 *
 *  Since we are in one dimension, the following identities are declared:
 *  @begin{verbatim}
 *    typedef raw_line_iterator    raw_cell_iterator;
 *    typedef line_iterator        cell_iterator;
 *    typedef active_line_iterator active_cell_iterator;
 *  @end{verbatim}
 *
 *  To enable the declaration of @p{begin_quad} and the like in
 *  @p{Triangulation<1>}, the @p{quad_iterator}s are declared as
 *  @p{void *}. Thus these types exist, but are useless and will
 *  certainly make any involuntary use visible. The same holds
 *  for hexahedron iterators.
 *
 *  The same applies for the @p{face_iterator} types, since lines
 *  have no substructures apart from vertices, which are handled in
 *  a different way, however.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <>
class TriaDimensionInfo<1>
{
  public:
    typedef TriaRawIterator<1,CellAccessor<1> >    raw_line_iterator;
    typedef TriaIterator<1,CellAccessor<1> >       line_iterator;
    typedef TriaActiveIterator<1,CellAccessor<1> > active_line_iterator;

    typedef void * raw_quad_iterator;
    typedef void * quad_iterator;
    typedef void * active_quad_iterator;

    typedef void * raw_hex_iterator;
    typedef void * hex_iterator;
    typedef void * active_hex_iterator;

    typedef raw_line_iterator    raw_cell_iterator;
    typedef line_iterator        cell_iterator;
    typedef active_line_iterator active_cell_iterator;

    typedef void * raw_face_iterator;
    typedef void * face_iterator;
    typedef void * active_face_iterator;
};



/**
 *  This class implements some types which differ between the dimensions.
 *  These are the declararions for the 2D case only.
 *
 *  A @p{line_iterator} is typdef'd to an iterator operating on the
 *  @p{lines} member variable of a @p{Triangulation<2>} object. An
 *  @p{active_line_iterator} only operates on the active lines.
 *  @p{raw_line_iterator} objects operate on all lines, used or not.
 *  Using @p{active_line_iterator}s may not be particularly in 2D useful since it
 *  only operates on unrefined lines. However, also refined lines may bound
 *  unrefined cells if the neighboring cell is refined once more than the
 *  present one.
 *
 *  Similarly to line iterators, @p{quad_iterator}, @p{raw_quad_iterator} and
 *  @p{active_quad_iterator} are declared.
 *  
 *  To enable the declaration of @p{begin_hex} and the like in
 *  @p{Triangulation<[12]>}, the @p{hex_iterator}s are declared as
 *  @p{void *}. Thus these types exist, but are useless and will
 *  certainly make any involuntary use visible.
 *
 *  Since we are in two dimension, the following identities are declared:
 *  @begin{verbatim}
 *    typedef raw_quad_iterator    raw_cell_iterator;
 *    typedef quad_iterator        cell_iterator;
 *    typedef active_quad_iterator active_cell_iterator;
 *
 *    typedef raw_line_iterator    raw_face_iterator;
 *    typedef line_iterator        face_iterator;
 *    typedef active_line_iterator active_face_iterator;    
 *  @end{verbatim}
 *
 * @author Wolfgang Bangerth, 1998
 */
template <>
class TriaDimensionInfo<2>
{
  public:
    typedef TriaRawIterator<2,TriaObjectAccessor<1, 2> >    raw_line_iterator;
    typedef TriaIterator<2,TriaObjectAccessor<1, 2> >       line_iterator;
    typedef TriaActiveIterator<2,TriaObjectAccessor<1, 2> > active_line_iterator;
    
    typedef TriaRawIterator<2,CellAccessor<2> >    raw_quad_iterator;
    typedef TriaIterator<2,CellAccessor<2> >       quad_iterator;
    typedef TriaActiveIterator<2,CellAccessor<2> > active_quad_iterator;

    typedef void * raw_hex_iterator;
    typedef void * hex_iterator;
    typedef void * active_hex_iterator;

    typedef raw_quad_iterator    raw_cell_iterator;
    typedef quad_iterator        cell_iterator;
    typedef active_quad_iterator active_cell_iterator;

    typedef raw_line_iterator    raw_face_iterator;
    typedef line_iterator        face_iterator;
    typedef active_line_iterator active_face_iterator;    
};



/**
 *  This class implements some types which differ between the dimensions.
 *  These are the declararions for the 3D case only.
 *
 *  For the declarations of the data types, more or less the same holds
 *  as for lower dimensions (see @p{TriaDimensionInfo<[12]>}). The
 *  dimension specific data types are here, since we are in three dimensions:
 *  @begin{verbatim}
 *    typedef raw_hex_iterator    raw_cell_iterator;
 *    typedef hex_iterator        cell_iterator;
 *    typedef active_hex_iterator active_cell_iterator;
 *
 *    typedef raw_quad_iterator    raw_face_iterator;
 *    typedef quad_iterator        face_iterator;
 *    typedef active_quad_iterator active_face_iterator;    
 *  @end{verbatim}
 *
 * @author Wolfgang Bangerth, 1998
 */
template <>
class TriaDimensionInfo<3>
{
  public:
    typedef TriaRawIterator<3,TriaObjectAccessor<1, 3> >    raw_line_iterator;
    typedef TriaIterator<3,TriaObjectAccessor<1, 3> >       line_iterator;
    typedef TriaActiveIterator<3,TriaObjectAccessor<1, 3> > active_line_iterator;
    
    typedef TriaRawIterator<3,TriaObjectAccessor<2, 3> >    raw_quad_iterator;
    typedef TriaIterator<3,TriaObjectAccessor<2, 3> >       quad_iterator;
    typedef TriaActiveIterator<3,TriaObjectAccessor<2, 3> > active_quad_iterator;

    typedef TriaRawIterator<3,CellAccessor<3> >    raw_hex_iterator;
    typedef TriaIterator<3,CellAccessor<3> >       hex_iterator;
    typedef TriaActiveIterator<3,CellAccessor<3> > active_hex_iterator;

    typedef raw_hex_iterator    raw_cell_iterator;
    typedef hex_iterator        cell_iterator;
    typedef active_hex_iterator active_cell_iterator;

    typedef raw_quad_iterator    raw_face_iterator;
    typedef quad_iterator        face_iterator;
    typedef active_quad_iterator active_face_iterator;    
};



/*------------------------------------------------------------------------*/

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
 * Note: these class should be made local to the triangulation class
 * once the compiler supports that (gcc2.95 does not at present).
 *
 * @author Wolfgang Bangerth, 1999
 */
template <int dim>
struct TriaNumberCache;

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
 * Note: these class should be made local to the triangulation class
 * once the compiler supports that (gcc2.95 does not at present).
 *
 * @author Wolfgang Bangerth, 1999
 */
template <>
struct TriaNumberCache<1> 
{
				     /**
				      * Number of used lines in the whole
				      * triangulation.
				      */
    unsigned int n_lines;

				     /**
				      * Array holding the number of used
				      * lines on each level.
				      */
    vector<unsigned int> n_lines_level;
    
				     /**
				      * Number of active lines in the
				      * whole triangulation.
				      */
    unsigned int n_active_lines;

				     /**
				      * Array holding the number of active
				      * lines on each level.
				      */
    vector<unsigned int> n_active_lines_level;

				     /**
				      * Constructor. Set values to zero
				      * by default.
				      */
    TriaNumberCache ();
};


/**
 * Cache class used to store the number of used and active elements
 * (lines or quads etc) within the levels of a triangulation. This
 * specialization stores the numbers of quads. Due to the inheritance
 * from the base class @ref{TriaNumberCache<1>}, the numbers of lines
 * are also within this class.
 *
 * In the old days, whenever one wanted to access one of these
 * numbers, one had to perform a loop over all lines, e.g., and count
 * the elements until we hit the end iterator. This is time consuming
 * and since access to the number of lines etc is a rather frequent
 * operation, this was not an optimal solution.
 *
 * Note: these class should be made local to the triangulation class
 * once the compiler supports that (gcc2.95 does not at present).
 *
 * @author Wolfgang Bangerth, 1999
 */
template <>
struct TriaNumberCache<2> : public TriaNumberCache<1>
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
    vector<unsigned int> n_quads_level;
    
				     /**
				      * Number of active quads in the
				      * whole triangulation.
				      */
    unsigned int n_active_quads;

				     /**
				      * Array holding the number of active
				      * quads on each level.
				      */
    vector<unsigned int> n_active_quads_level;

				     /**
				      * Constructor. Set values to zero
				      * by default.
				      */
    TriaNumberCache ();
};


/**
 * Cache class used to store the number of used and active elements
 * (lines or quads etc) within the levels of a triangulation. This
 * specialization stores the numbers of hexes. Due to the inheritance
 * from the base class @ref{TriaNumberCache<2>}, the numbers of lines
 * and quads are also within this class.
 *
 * In the old days, whenever one wanted to access one of these
 * numbers, one had to perform a loop over all lines, e.g., and count
 * the elements until we hit the end iterator. This is time consuming
 * and since access to the number of lines etc is a rather frequent
 * operation, this was not an optimal solution.
 *
 * Note: these class should be made local to the triangulation class
 * once the compiler supports that (gcc2.95 does not at present).
 *
 * @author Wolfgang Bangerth, 1999
 */
template <>
struct TriaNumberCache<3> : public TriaNumberCache<2>
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
    vector<unsigned int> n_hexes_level;
    
				     /**
				      * Number of active hexes in the
				      * whole triangulation.
				      */
    unsigned int n_active_hexes;

				     /**
				      * Array holding the number of active
				      * hexes on each level.
				      */
    vector<unsigned int> n_active_hexes_level;

				     /**
				      * Constructor. Set values to zero
				      * by default.
				      */
    TriaNumberCache ();
};


/*------------------------------------------------------------------------*/


/**
 *  @ref{Triangulation}s denote a hierarchy of levels of elements which together
 *  form a region in @p{dim} spatial dimensions.
 *
 *  This class is written to be as independent of the dimension as possible
 *  (thus the complex construction of the @ref{TriangulationLevel} classes) to
 *  allow code-sharing, to allow reducing the need to mirror changes in the code
 *  for one dimension to the code for other dimensions. Nonetheless, some of
 *  the functions are dependent of the dimension and there only exist
 *  specialized versions for distinct dimensions.
 *
 *
 *  @sect3{Structure and iterators}
 *
 *  The actual data structure of a @ref{Triangulation} object is rather complex
 *  and quite inconvenient if one attempted to operate on it directly, since
 *  data is spread over quite a lot of arrays and other places. However,
 *  there are ways powerful enough to work on these data structures
 *  without knowing their exact relations. This is done through the
 *  concept of iterators (see the STL documentation and @ref{TriaRawIterator}).
 *  In order to make things as easy and dimension independent as possible,
 *  use of class local typedefs is made, see below.
 *  
 *  In the base class @ref{TriaDimensionInfo}, a @p{Cell} is typedef'd to be whatever
 *  is reasonable for a cell in the respective dimension, i.e. a @p{Line} in
 *  one dimension, a @p{Quad} in two dimensions, and so on.
 *
 *  The @ref{Triangulation} class provides iterator which enable looping over all
 *  lines, cells,
 *  etc without knowing the exact representation used to describe them. Their
 *  names are typedefs in the @ref{TriaDimensionInfo} base class (thus making them
 *  local types to this class) and are as follows:
 *
 *  @begin{itemize}
 *  @item @p{raw_line_iterator}: loop over all lines, used or not (declared for
 *  all dimensions).
 *  
 *  @item @p{line_iterator}: loop over all used lines (declared for all dimensions).
 *
 *  @item @p{active_line_iterator}: loop over all active lines (declared for all
 *  dimensions).
 *
 *  @item @p{raw_quad_iterator}: loop over all quads, used or not (declared only
 *  for @p{dim>=2}).
 *  
 *  @item @p{quad_iterator}: loop over all quads (declared only for @p{dim}>=2).
 *
 *  @item @p{active_quad_iterator}: loop over all active quads (declared only for
 *  @p{dim}>=2).
 *  @end{itemize}
 *
 *  Additionaly, for @p{dim}==1, the following identities hold:
 *  @begin{verbatim}
 *    typedef raw_line_iterator    raw_cell_iterator;
 *    typedef line_iterator        cell_iterator;
 *    typedef active_line_iterator active_cell_iterator;
 *  @end{verbatim}
 *  while for @p{dim}==2
 *  @begin{verbatim}
 *    typedef quad_line_iterator   raw_cell_iterator;    
 *    typedef quad_iterator        cell_iterator;
 *    typedef active_quad_iterator active_cell_iterator;
 *
 *    typedef raw_line_iterator    raw_face_iterator;
 *    typedef line_iterator        face_iterator;
 *    typedef active_line_iterator active_face_iterator;    
 *  @end{verbatim}
 *
 *  By using the cell iterators, you can write code nearly independent of
 *  the spatial dimension. The same applies for substructure iterators,
 *  where a substructure is defined as a face of a cell. The face of a
 *  cell is be a vertex in 1D and a line in 2D; however, vertices are
 *  handled in a different way and therefore lines have no faces.
 *
 *  The @ref{Triangulation} class offers functions like @p{begin_active} which gives
 *  you an iterator to the first active cell. There are quite a lot of functions
 *  returning iterators. Take a look at the class doc to get an overview.
 *
 *  Usage of these iterators works mostly like with the STL iterators. Some
 *  examples taken from the @ref{Triangulation} source code follow.
 *  @begin{itemize}
 *  @item @em{Counting the number of cells on a specific level}
 *    @begin{verbatim}
 *     template <int dim>
 *     int Triangulation<dim>::n_cells (const int level) const {
 *        cell_iterator cell = begin (level),
 *                      endc = end(level);
 *        int n=0;
 *        for (; cell!=endc; ++cell)
 *          ++n;
 *        return n;
 *      };
 *    @end{verbatim}
 *    Another way which uses the STL @p{distance} function would be to write
 *    @begin{verbatim}
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
 *    @end{verbatim}
 *    
 *  @item @em{Refining all cells of a triangulation}
 *    @begin{verbatim}
 *      template <int dim>
 *      void Triangulation<dim>::refine_global () {
 *        active_cell_iterator cell = begin_active(),
 *                             endc = end();
 *
 *        for (; cell != endc; ++cell)
 *          cell->set_refine_flag ();
 *        execute_coarsening_and_refinement ();
 *      };
 *    @end{verbatim}
 *  @end{itemize}
 *
 *
 *  @sect3{Usage}
 *
 *  Usage of a @ref{Triangulation} is mainly done through the use of iterators.
 *  An example probably shows best how to use it:
 *  @begin{verbatim}
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
 *  @end{verbatim}
 *
 *  
 *  @sect3{Creating a triangulation}
 *
 *  There are several possibilities to create a triangulation:
 *  @begin{itemize}
 *    @item The most common domains, such as hypercubes (i.e. lines, squares,
 *       cubes, etc), hyper-balls (circles, balls, ...) and some other, more
 *       weird domains such as the L-shape region and higher dimensional
 *       generalizations and others, are provided by the @ref{GridGenerator}
 *       class which takes a triangulation and fills it by a division
 *       of the required domain.
 *   
 *     @item Reading in a triangulation: By using an object of the @ref{GridIn}
 *        class, you can read in fairly general triangulations. See there for
 *        more information. The mentioned class uses the interface described
 *        directly below to transfer the data into the triangulation.
 *
 *     @item Explicitely creating a triangulation: you can create a triangulation
 *        by providing a list of vertices and a list of cells. Each such cell
 *        consists of a vector storing the indices of the vertices of this cell
 *        in the vertex list. To see how this works, you can take a look at the
 *        @ref{GridIn}@p{<dim>::read_*} functions. The appropriate function to be
 *        called is @p{Triangulation<dim>::create_triangulation (2)}.
 *
 *        Creating the hierarchical information needed for this library from
 *        cells storing only vertex information can be quite a complex task.
 *        For example in 2d, we have to create lines between vertices (but only
 *        once, though there are two cells which link these two vertices) and
 *        we have to create neighborship information. Grids being read in
 *        should therefore not be too large, reading refined grids would be
 *        inefficient. Apart from the performance aspect, refined grids do not
 *        lend too well to multigrid algorithms, since solving on the coarsest
 *        level is expensive. It is wiser in any case to read in a grid as coarse
 *        as possible and then do the needed refinement steps.
 *
 *        It is your duty to guarantee that cells have the correct orientation.
 *        To guarantee this, in the input vector keeping the cell list, the
 *        vertex indices for each cell have to be in a defined order. In one
 *        dimension, the first vertex index must refer to that vertex with the
 *        lower coordinate value. In two dimensions, the four vertices must be
 *        given in an order representing a counterclockwise sense. This
 *        condition is not easy to verify and no full attempt to do so is made.
 *        If you violate this condition, you may end up with matrix entries
 *        having the wrong sign (clockwise vertex numbering, which results in
 *        a negative area element) of with wrong matrix elements (twisted
 *        quadrilaterals, i.e. two vertices interchanged; this results in
 *        a wrong area element).
 *
 *        There are more subtle conditions which must be imposed upon the
 *        vertex numbering within cells. See the documentation for the
 *        @ref{GridIn} class for more details on this. They do not only
 *        hold for the data read from an UCD or any other input file, but
 *        also for the data passed to the
 *        @p{Triangulation<dim>::create_triangulation (2)} function.
 *
 *     @item Copying a triangulation: when computing on time dependant meshes
 *        of when using adaptive refinement, you will often want to create a
 *        new triangulation to be the same as another one. This is facilitated
 *        by the @p{copy_triangulation} function.
 *
 *        It is guaranteed that vertex, line or cell numbers in the two
 *        triangulations are the same and that two iterators walking on the
 *        two triangulations visit matching cells if the are incremented in
 *        parallel. It may be conceivable to implement a clean-up in the copy
 *        operation, which eliminates holes of unused memory, re-joins
 *        scattered data and so on. In principle this would be a useful
 *        operation but guaranteeing some parallelity in the two triangulations
 *        seems more important since usually data will have to be transferred
 *        between the grids.
 *   @end{itemize}
 *
 *   The material id for each cell must be specified upon construction of
 *   a triangulation. (There is a special section on material identifier and
 *   boundary indicators. See there for more information.)
 *   The standard region functions (for hypercube, hyper-ball,
 *   etc.) denote all cells the material id zero. You may change that afterwards,
 *   but you should not use the material id 255. When reading a triangulation,
 *   the material id must be specified in the input file (UCD format) or is
 *   otherwise set to zero. When creating explicitely, the material id must
 *   be given to the creation function.
 *
 *   Regarding the boundary indicator for lines in two dimensions and quads
 *   in three (subsumed by the word "faces"), all interior faces are denoted
 *   the value 255. Trying to give an interior face another value results in
 *   an error if in debug mode. Faces at the boundary of the domain are preset
 *   with the boundary indicator zero, but you can give a list of faces with
 *   different boundary indicators to the triangulation creation function.
 *   The standard domain functions assume all faces to have boundary indicator
 *   zero, which you may change manually afterwards. When reading from a file,
 *   you have to give boundary indicators other than zero explicitely, e.g. in
 *   UCD format by giving a list of lines with material id in the input file.
 *
 *   Lines in two dimensions and quads in three dimensions inherit their
 *   boundary indicator to their children upon refinement. You should therefore
 *   make sure that if you have different boundary parts, the different parts
 *   are separated by a vertex (in 2D) or a line (in 3D) such that each boundary
 *   line or quad has a unique boundary indicator.
 *
 *   Likewise, material data is inherited from mother to child cells. Place your
 *   coarse level cells so, that the interface between cells is also the
 *   interface between regions of different materials.
 *
 *   Finally, there is a special function for folks who like bad grids:
 *   @p{Triangulation<dim>::distort_random}. It moves all the vertices in the
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
 *   @sect3{Refinement and coarsening of a triangulation}
 *
 *   Refinement of a triangulation may be done through several ways. The most
 *   low-level way is directly through iterators: let @p{i} be an iterator to
 *   an active cell (i.e. the cell pointed to has no children), then the
 *   function call @p{i->set_refine_flag()} marks the respective cell for
 *   refinement. Marking non-active cells results in an error.
 *
 *   After all the cells you wanted to mark for refinement, call the
 *   @p{execute_coarsening_and_refinement} function to actually perform
 *   the refinement. This function itself first calls the
 *   @p{prepare_coarsening_and_refinement} function to regularize the resulting
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
 *   @p{i->set_coarsen_flag} and calling @p{execute_coarsening_and_refinement}.
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
 *   elements. These functions can be found in the @ref{GridRefinement} class.
 *
 *
 *   @sect3{Smoothing of a triangulation}
 *
 *   Some degradation of approximation properties has been observed
 *   for grids which are too unstructured. Therefore, the
 *   @p{prepare_coarsening_and_refinement} function which is automatically called
 *   by the @p{execute_coarsening_and_refinement} function can do some
 *   smoothing of the triangulation. Note that mesh smoothing is only
 *   done for two or more space dimensions, no smoothing is available
 *   at present for one spatial dimension. In the sequel, let
 *   @p{execute_*} stand for @p{execute_coarsening_and_refinement}.
 *
 *   For the purpose of smoothing, the
 *   @ref{Triangulation} constructor takes an argument specifying whether a
 *   smoothing step shall be performed on the grid each time @p{execute_*}
 *   is called. The default is that such a step not be done, since this results
 *   in additional cells being produced, which may not be necessary in all
 *   cases. If switched on, calling @p{execute_*} results in
 *   flagging additional cells for refinement to avoid
 *   vertices as the ones mentioned. The algorithms for both regularization
 *   and smoothing of triangulations are described below in the section on
 *   technical issues. The reason why this parameter must be given to the
 *   constructor rather than to @p{execute_*} is that it would result
 *   in algorithmic problems if you called @p{execute_*} once without
 *   and once with smoothing, since then in some refinement steps would need
 *   to be refined twice.
 *
 *   The parameter taken by the constructor is an integer which may be composed
 *   bitwise by the constants defined in the @p{enum MeshSmoothing}. The meaning
 *   of these constants is explained in the following:
 *   @begin{itemize}
 *   @item @p{limit_level_difference_at_vertices}:
 *     It can be shown, that degradation of approximation occurs if the
 *     triangulation contains vertices which are member of cells with levels
 *     differing by more than one. One such example is the following:
 *     @begin{verbatim}
 *       |     |     |     |
 *       x-----x-----x--x--x--
 *       |     |     |  |  |
 *       |     |     x--x--x
 *       |     |     |  |  |
 *       x-----x-----x--x--x--
 *       |           |     |
 *       |           |     |
 *       |           |     |
 *       |           x-----x--
 *       |           |     |
 *       |           |     |
 *       |           |     |
 *       x-----------x-----x--
 *     @end{verbatim}
 *     It seems that in two space dimensions, the maximum jump in levels between
 *     cells sharing a common vertex is two (as in the example above). This is
 *     not true if more than four cells meet at a vertex. It is not uncommon
 *     that a coarse (initial) mesh contains vertices at which six or even eight
 *     cells meet, when small features of the domain have to be resolved even on
 *     the coarsest mesh. In that case, the maximum difference in levels is
 *     three or four, respectively. The problem gets even worse in three space
 *     dimensions.
 *
 *     Looking at an interpolation of the second derivative of the finite
 *     element solution (assuming bilinear finite elements), one sees that the
 *     numerical solution is almost totally wrong, compared with the true second
 *     derivative. Indeed, on regular meshes, there exist sharp estimations that
 *     the $H^2$-error is only $O(1)$, so we should not be surprised; however, the
 *     numerical solution may show a value for the second derivative which may
 *     be a factor of ten away from the true value. These problems are located
 *     on the small cell adjacent to the center vertex, where cells of
 *     non-subsequent levels meet, as well as on the upper and right neighbor
 *     of this cell (but with a less degree of deviation from the true value).
 *
 *     If the smoothing indicator given to the constructor contains the bit for
 *     @p{limit_level_difference_at_vertices}, situations as the above one are
 *     eliminated by also marking the lower left cell for refinement.
 *
 *   @item @p{eliminate_unrefined_islands}:
 *     Single cells which are not refined and are surrounded by cells which are
 *     refined usually also lead to a sharp decline in approximation properties
 *     locally. The reason is that the nodes on the faces between unrefined and
 *     refined cells are not real degrees of freedom but carry constraints. The
 *     patch without additional degrees of freedom is thus significantly larger
 *     then the unrefined cell itself. If in the parameter passed to the
 *     constructor the bit for @p{eliminate_unrefined_islands} is set, all cells
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
 *   @item @p{eliminate_refined_*_islands}:
 *     This algorithm seeks for isolated cells which are refined or flagged
 *     for refinement. This definition is unlike that for
 *     @p{eliminate_unrefined_islands}, which would mean that an island is
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
 *     @p{eliminate_refined_inner_islands} and @p{eliminate_refined_boundary_islands}.
 *     There first eliminates islands defined by the definition above which are
 *     in the interior of the domain, while the second eliminates only those
 *     islands if the cell is at the boundary. The reason for this split of
 *     flags is that one often wants to eliminate such islands in the interior
 *     while those at the boundary may well be wanted, for example if one
 *     refines the mesh according to a criterion associated with a boundary
 *     integral or if one has rough boundary data.
 *
 *   @item @p{do_not_produce_unrefined_islands}:
 *     This flag prevents the occurrence of unrefined islands. In more detail:
 *     It prohibits the coarsening of a cell if 'most of the neighbors' will
 *     be refined after the step.
 *
 * @item @p{patch_level_1}:
 *     Ensures patch level 1. As result the triangulation consists of
 *     patches, i.e. of cells that are refined once. It follows that
 *     if at least one of the children of a cell is or will be refined
 *     than all children need to be refined. If the @p{patch_level_1} flag
 *     is set, than the flags @p{eliminate_unrefined_islands},
 *     @p{eliminate_refined_inner_islands} and
 *     @p{eliminate_refined_boundary_islands} will be ignored as they will
 *     be fulfilled automatically.
 *
 *   @item @p{smoothing_on_refinement}:
 *     This flag sums up all smoothing algorithms which may be performed upon
 *     refinement by flagging some more cells for refinement.
 *
 *   @item @p{smoothing_on_coarsening}:
 *     This flag sums up all smoothing algorithms which may be performed upon
 *     coarsening by flagging some more cells for coarsening.
 *
 *   @item @p{maximum_smoothing}:
 *     This flag includes all the above ones and therefore combines all
 *     smoothing algorithms implemented.
 *
 *   @item @p{none}:
 *     Select no smoothing at all.
 *   @end{itemize}
 *
 *
 *   @sect3{Material and boundary information}
 *
 *   Each line, quad, etc stores one byte of information denoting the
 *   material of a cell or the part of the boundary, a lower
 *   dimensional object belongs to. The material of a cell may be used
 *   during matrix generation in order to implement different
 *   coefficients in different parts of the domain. It is not used by
 *   functions of the grid and dof handling libraries.
 *
 *   Boundary indicators on lower dimensional objects (these have no
 *   material id) indicate the number of a boundary component. These
 *   are used for two purposes: First, they specify a boundary
 *   curve. When a cell is refined, a function can be used to place
 *   new vertices on this curve. See the section on boundary
 *   approximation below. Furthermore, the the weak formulation of the
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
 *   @sect3{History of a triangulation}
 *
 *   It is possible to reconstruct a grid from its refinement history, which
 *   can be stored and loaded through the @p{save_refine_flags} and
 *   @p{load_refine_flags} functions. Normally, the code will look like this:
 *   @begin{verbatim}
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
 *   @end{verbatim}
 *
 *   If you want to re-create the grid from the stored information, you write:
 *   @begin{verbatim}
 *                                 // open input file
 *     ifstream history("mesh.history");
 *                                 // do 10 refinement steps
 *     for (int step=0; step<10; ++step) {
 *       tria.load_refine_flags (history);
 *       tria.execute_coarsening_and_refinement ();
 *     };        
 *   @end{verbatim}
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
 *   There actually are two sets of @p{save_*_flags} and @p{load_*_flags} functions.
 *   One takes a stream as argument and reads/writes the information from/to the
 *   stream, thus enabling storing flags to files. The other set takes an
 *   argument of type @p{vector<bool>}. This enables the user to temporarily store
 *   some flags, e.g. if another function needs them, and restore them
 *   afterwards.
 *
 *
 *   @sect3{User flags and data}
 *
 *   A triangulation offers one bit per line, quad, etc for user data.
 *   This field can be
 *   accessed as all other data using iterators. Normally, this user flag is
 *   used if an algorithm walks over all cells and needs information whether
 *   another cell, e.g. a neighbor, has already been processed. It can also
 *   be used to flag the lines subject to constraints in 2D, as for example the
 *   functions in the @ref{DoFHandler} classes do.
 *
 *   There are two functions, @p{save_user_flags} and @p{load_user_flags} which
 *   write and read these flags to and from a stream. Unlike
 *   @p{save_refine_flags} and @p{load_refine_flags}, these two functions store
 *   and read the flags of all used lines, quads, etc, not only of the
 *   active ones (well, activity is a concept which really only applies to
 *   cells, not for example to lines in 2D, so the abovementioned generalization
 *   to @em{all} lines, quads, etc seems plausible).
 *
 *   If you want to store more specific user flags, you can use the functions
 *   @p{save_user_flags_line} and @p{load_user_flags_line} and the generalizations
 *   for quads, etc.
 *
 *   As for the refinement and coarsening flags, there exist two versions of these
 *   functions, one which reads/writes from a stream and one which does so from
 *   a @p{vector<bool>}. The latter is used to store flags temporarily, while the
 *   first is used to store them in a file.
 *
 *   It is convention to clear the user flags using the
 *   @p{Triangulation<>::clear_user_flags()} function before usage, since it is
 *   often necessary to use the flags in more than one function consecutively and
 *   is then error prone to dedicate one of these to clear the flags.
 *
 *   It is recommended that a functions using the flags states so in its
 *   documentation.
 *
 *   There is another set of user data, namely a @p{void *}, for each
 *   line, quad, etc. You can access these user pointers through the
 *   functions @p{user_pointer()}, @p{clear_user_pointer()} and
 *   @p{set_user_pointer(p)} in the accessor classes. These pointers are
 *   neither used nor changed by the library and are @p{NULL} by
 *   default. Thus, their handling is the sole responsibility of the
 *   application program.  Especially, the pointers are not inherited
 *   to children upon refinement. Still, after a remeshing they are
 *   available on all cells, where they were set on the previous mesh.
 *
 *   The usual warning about the missing type safety of @p{void} pointers are
 *   obviously in place here; responsibility for correctness of types etc
 *   lies entirely with the user of the pointer.
 *
 *   Just like the user flags, this field is not available for vertices,
 *   which does no harm since the vertices have a unique and continuous number
 *   unlike the structured objects lines and quads. 
 *
 *   @sect3{Boundary approximation}
 *
 *   You can specify a boundary function for each boundary
 *   component. If a new vertex is created on a side or face at the
 *   boundary, this function is used to compute where it will be
 *   placed. The boundary indicator of the face will be used to
 *   determine the proper component. See @ref{Boundary} for the
 *   details. Usage with the @ref{Triangulation} object is then like this
 *   (let @p{Ball} be a class derived from @ref{Boundary}@p{<2>}):
 * 
 *   @begin{verbatim}
 *     void main () {
 *       Triangulation<2> tria;
 *                                        // set the boundary function
 *                                        // for all boundaries with
 *                                        // boundary indicator 0
 *       Ball ball;
 *       tria.set_boundary (0, &ball);
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
 *   @end{verbatim}
 *
 *   You should take note of one caveat: if you have concave
 *   boundaries, you must make sure that a new boundary vertex does
 *   not lie too much inside the cell which is to be refined. The
 *   reason is that the center vertex is placed at the point which is
 *   the arithmetic mean of the vertices of the original cell.
 *   Therefore if your new boundary vertex is too near the center of
 *   the old quadrilateral or hexahedron, the distance to the midpoint
 *   vertex will become too small, thus generating distorted
 *   cells. Remedy: take care of such situations when defining the
 *   coarse grid.
 *
 *
 *   @sect3{Technical details}
 *
 *   @sect4{Algorithms for mesh regularization and smoothing upon refinement}
 *
 *   We chose an inductive point of view: since upon creation of the
 *   triangulation all cells are on the same level, all regularity assumptions
 *   regarding the maximum difference in level of cells sharing a common face,
 *   edge or vertex hold. Since we use the regularization and smoothing in
 *   each step of the mesh history, when coming to the point of refining it
 *   further the assumptions also hold.
 *
 *   The regularization and smoothing is done in the
 *   @p{prepare_coarsening_and_refinement} function, which is called by
 *   @p{execute_coarsening_and_refinement} at the very beginning.  It
 *   decides which additional cells to flag for refinement by looking
 *   at the old grid and the refinement flags for each cell.
 *
 *   @begin{itemize}
 *   @item @em{Regularization:} The algorithm walks over all cells checking
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
 *   @item @em{Smoothing:}
 *     @begin{itemize}
 *     @item @p{limit_level_difference_at_vertices}:
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
 *     @item @p{eliminate_unrefined_islands}:
 *       For each cell we count the number of neighbors which are refined or
 *       flagged for refinement. If this exceeds the total number of neighbors
 *       (which is the number of faces minus the number of faces of this cell
 *       which are located on the boundary), then this cell is flagged for
 *       refinement. Since this may lead to cells on the same level which also
 *       will need refinement, we will need additional loops of regularization
 *       and smoothing over all cells until nothing changes any more.
 *
 *     @item @p{eliminate_refined_*_islands}:
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
 *     @end{itemize}
 *   @end{itemize}
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
 *   @sect4{Implementation conventions for two spatial dimensions}
 *
 *   There is a convention about the direction of the bounding lines of quads in
 *   2D. The direction of a line is the direction of point 0 towards point 1. We
 *   define, that allowed cells contain of lines of which the direction is
 *   as follows:
 *   @begin{verbatim}
 *          2
 *      3--->---2
 *      |       |
 *     3^       ^1
 *      |       |
 *      0--->---1
 *          0
 *   @end{verbatim}
 *   The number of the vertices and lines is also indicated. This orientation of
 *   lines has to be checked/generated upon construction of a grid and is
 *   preserved upon refinement.
 *
 *   Further we define, that child lines have the same direction as their parent,
 *   i.e. that @p{subline(0).vertex(0)==line.vertex(0)} and
 *   @p{subline(1).vertex(1)==line.vertex(1)}. This also implies, that the
 *   first subline (@p{subline(0)}) is the one at vertex(0) of the old line.
 *
 *   Similarly we define, that the four children of a quad are adjacent to the
 *   vertex with the same number of the old quad.
 *
 *   @sect4{Coordinate systems}
 *
 *   When explicit coordinates are required for points in a cell (e.g for
 *   quadrature formulae or the point of definition of trial functions), we
 *   define the following coordinate system for the unit cell:
 *   @begin{verbatim}
 *    y^   3-------2
 *     |   |       |
 *     |   |       |
 *     |   |       |
 *     |   0-------1
 *     *-------------->x
 *   @end{verbatim}
 *   with vertex 0 being the origin of the coordinate system, vertex 1 having
 *   coordinates @p{(1,0)}, vertex 2 at @p{(1,1)} and vertex 3 at @p{(0,1)}.
 *
 *
 *   @sect3{Implementation conventions for three spatial dimensions}
 *
 *   By convention, we will use the following numbering for vertices, lines and
 *   faces of hexahedra in three space dimensions. Before giving these 
 *   conventions we declare the following sketch to be the standard way of
 *   drawing 3d pictures of hexahedra:
 *   @begin{verbatim}
 *         *-------*        *-------*
 *        /|       |       /       /|
 *       / |       |      /       / |
 *      /  |       |     /       /  |
 *     *   |       |    *-------*   |
 *     |   *-------*    |       |   *
 *     |  /       /     |       |  /
 *     | /       /      |       | /
 *     |/       /       |       |/
 *     *-------*        *-------*
 *   @end{verbatim}
 *   The left part of the picture shows the left, bottom and back face of the
 *   cube, while the right one shall be the top, right and front face. You may
 *   recover the whole cube by moving the two parts together into one.
 *
 *   @sect4{Vertices}
 *
 *   The vertices on the front face are numbered exactly the same way as are
 *   the vertices on a quadrilateral. The vertices on the back face are numbered
 *   similarly by moving the front face to the back (no turning, no twisting, 
 *   just a shift):
 *   @begin{verbatim}
 *         7-------6        7-------6
 *        /|       |       /       /|
 *       / |       |      /       / |
 *      /  |       |     /       /  |
 *     3   |       |    3-------2   |
 *     |   4-------5    |       |   5
 *     |  /       /     |       |  /
 *     | /       /      |       | /
 *     |/       /       |       |/
 *     0-------1        0-------1
 *   @end{verbatim}
 *
 *   @sect4{Lines}
 *
 *   Here, the same holds as for the vertices: the lines of the front face are
 *   numbered as for the quadrilateral, for the back face they are just shifted.
 *   Finally, the four lines connecting front and back face are numbered:
 *   @begin{verbatim}
 *         *---6---*        *---6---*
 *        /|       |       /       /|
 *      11 |       5      11     10 5
 *      /  7       |     /       /  |
 *     *   |       |    *---2---*   |
 *     |   *---4---*    |       |   *
 *     |  /       /     |       1  /
 *     3 8       9      3       | 9
 *     |/       /       |       |/
 *     *---0---*        *---0---*
 *   @end{verbatim}
 *   The directions of the front and back lines is as for the respective faces, while
 *   the connecting lines always point to the back:
 *   @begin{verbatim}
 *         *--->---*        *--->---*
 *        /|       |       /       /|
 *       ^ |       ^      ^       ^ ^
 *      /  ^       |     /       /  |
 *     *   |       |    *--->---*   |
 *     |   *--->---*    |       |   *
 *     |  /       /     |       ^  /
 *     ^ ^       ^      ^       | ^
 *     |/       /       |       |/
 *     *--->---*        *--->---*
 *   @end{verbatim}
 *
 *   @sect4{Faces}
 *
 *   The faces are numbered in the same order as the lines were numbered: front
 *   face, back face, then the four side faces:
 *   @begin{verbatim}
 *         *-------*        *-------*
 *        /|       |       /       /|
 *       / |   1   |      /   4   / |
 *      /  |       |     /       /  |
 *     *   |       |    *-------*   |
 *     | 5 *-------*    |       | 3 *
 *     |  /       /     |       |  /
 *     | /   2   /      |   0   | /
 *     |/       /       |       |/
 *     *-------*        *-------*
 *   @end{verbatim}
 *
 *   The direction of the faces is determined by the numbers the lines have within
 *   a given face. This is like follows:
 *   @begin{itemize}
 *   @item Faces 0 and 1:
 *    @begin{verbatim}
 *          *---2---*        *-------*
 *         /|       |       /       /|
 *        / |       1      /       / |
 *       /  3       |     /       /  |
 *      *   |       |    *---2---*   |
 *      |   *---0---*    |       |   *
 *      |  /       /     |       1  /
 *      | /       /      3       | /
 *      |/       /       |       |/
 *      *-------*        *---0---*
 *    @end{verbatim}
 * 
 *   @item Faces 2 and 4:
 *    @begin{verbatim}
 *          *-------*        *---2---*
 *         /|       |       /       /|
 *        / |       |      3       1 |
 *       /  |       |     /       /  |
 *      *   |       |    *---0---*   |
 *      |   *---2---*    |       |   *
 *      |  /       /     |       |  /
 *      | 3       1      |       | /
 *      |/       /       |       |/
 *      *---0---*        *-------*
 *    @end{verbatim} 
 * 
 *   @item Faces 3 and 5:
 *    @begin{verbatim}
 *          *-------*        *-------*
 *         /|       |       /       /|
 *        2 1       |      /       2 1
 *       /  |       |     /       /  |
 *      *   |       |    *-------*   |
 *      |   *-------*    |       |   *
 *      3  /       /     |       3  /
 *      | 0       /      |       | 0
 *      |/       /       |       |/
 *      *-------*        *-------*
 *    @end{verbatim}
 *   @end{itemize}
 * 
 *   Due to this numbering, the following lines are identical:
 *   @begin{itemize}
 *   @item Line 0 of face 0, and line 0 of face 2;
 *   @item Line 1 of face 0, and line 3 of face 3;
 *   @item Line 2 of face 0, and line 0 of face 4;
 *   @item Line 3 of face 0, and line 3 of face 5;
 *   @item Line 0 of face 1, and line 2 of face 2;
 *   @item Line 1 of face 1, and line 1 of face 3;
 *   @item Line 2 of face 1, and line 2 of face 4;
 *   @item Line 3 of face 1, and line 1 of face 5;
 *   @item Line 3 of face 2, and line 0 of face 5;
 *   @item Line 1 of face 2, and line 0 of face 3;
 *   @item Line 1 of face 4, and line 2 of face 3;
 *   @item Line 3 of face 4, and line 2 of face 5.
 *   @end{itemize}
 *
 *
 *   @sect4{Children}
 *
 *   The eight children of a cell are numbered as follows:
 *   @begin{verbatim}
 *         *-------*        *-------*
 *        /| 7   6 |       / 7   6 /|
 *       /7|       |      /       /6|
 *      /  |       |     / 3   2 /  |
 *     *   | 4   5 |    *-------*2 5|
 *     |3 4*-------*    | 3   2 |   *
 *     |  / 4   5 /     |       |  /
 *     |0/       /      |       |1/
 *     |/0    1 /       | 0   1 |/
 *     *-------*        *-------*
 *   @end{verbatim}
 *
 *   Taking into account the orientation of the faces, the following
 *   children are adjacent to the respective faces:
 *   @begin{itemize}
 *   @item Face 0: children 0, 1, 2, 3;
 *   @item Face 1: children 4, 5, 6, 7;
 *   @item Face 2: children 0, 1, 5, 4;
 *   @item Face 3: children 1, 5, 6, 2;
 *   @item Face 4: children 3, 2, 6, 7;
 *   @item Face 5: children 0, 4, 7, 3.
 *   @end{itemize}
 *   You can get these numbers using the @ref{GeometryInfo<3>}@p{::child_cell_on_face}
 *   function. Each child is adjacent to the vertex with the same number.
 *
 *
 *   @sect4{Coordinate systems}
 *
 *   We define the following coordinate system for the explicit coordinates of
 *   the vertices of the unit cell:
 *   @begin{verbatim}
 *                         7-------6        7-------6
 *                        /|       |       /       /|
 *                       / |       |      /       / |
 *    z                 /  |       |     /       /  |
 *    ^                3   |       |    3-------2   |
 *    |   ^y           |   4-------5    |       |   5
 *    |  /             |  /       /     |       |  /
 *    | /              | /       /      |       | /
 *    |/               |/       /       |       |/
 *    *------>x        0-------1        0-------1
 *   @end{verbatim}
 *   This convention in conjunction with the numbering of the vertices is a bit
 *   unfortunate, since the vertices 0 through 3 have the coordinates @p{(x,0,z)}
 *   with @p{x} and @p{z} being the same as the @p{x} and @p{y} coordinates of a quad
 *   in the plane; more intuitive would have been if they had the coordinates
 *   @p{(x,y,0)}. However, the vertex numbering was historically chosen as shown.
 *
 *   By the convention laid down as above, the vertices have the following
 *   coordinates:
 *   @begin{itemize}
 *      @item Vertex 0: @p{(0,0,0)};
 *      @item Vertex 1: @p{(1,0,0)};
 *      @item Vertex 2: @p{(1,0,1)};
 *      @item Vertex 3: @p{(0,0,1)};
 *      @item Vertex 4: @p{(0,1,0)};
 *      @item Vertex 5: @p{(1,1,0)};
 *      @item Vertex 6: @p{(1,1,1)};
 *      @item Vertex 7: @p{(0,1,1)}.
 *   @end{itemize}
 *
 *
 *   @sect3{Warning}
 *
 *   It seems impossible to preserve @p{const}ness of a triangulation through
 *   iterator usage. Thus, if you declare pointers to a @p{const} triangulation
 *   object, you should be well aware that you might involuntarily alter the
 *   data stored in the triangulation.
 *
 *   @see TriaRawIterator
 *   @author Wolfgang Bangerth, 1998
 */
template <int dim>
class Triangulation : public TriaDimensionInfo<dim>,
		      public Subscriptor
{
  private:
				     /**
				      * Default boundary object. This
				      * declaration is used for the
				      * default argument in
				      * @p{set_boundary}.
				      */
    static const StraightBoundary<dim>& straight_boundary;

  public:
    
				     /**
				      * Declare some symbolic names
				      * for mesh smoothing
				      * algorithms. The meaning of
				      * these flags is documented in
				      * the @ref{Triangulation} class.
				      */
    enum MeshSmoothing
    {
	  none                               = 0x0,
	  limit_level_difference_at_vertices = 0x1,     
	  eliminate_unrefined_islands        = 0x2,
	  patch_level_1                      = 0x4,
	  
	  eliminate_refined_inner_islands    = 0x100,
	  eliminate_refined_boundary_islands = 0x200,
	  do_not_produce_unrefined_islands   = 0x400,
	  
	  smoothing_on_refinement            = (limit_level_difference_at_vertices |
						eliminate_unrefined_islands),
	  smoothing_on_coarsening            = (eliminate_refined_inner_islands |
						eliminate_refined_boundary_islands |
						do_not_produce_unrefined_islands),
	  
	  maximum_smoothing                  = 0xffff
    };
    

    typedef typename TriaDimensionInfo<dim>::raw_line_iterator raw_line_iterator;
    typedef typename TriaDimensionInfo<dim>::line_iterator line_iterator;
    typedef typename TriaDimensionInfo<dim>::active_line_iterator active_line_iterator;

    typedef typename TriaDimensionInfo<dim>::raw_quad_iterator raw_quad_iterator;
    typedef typename TriaDimensionInfo<dim>::quad_iterator quad_iterator;
    typedef typename TriaDimensionInfo<dim>::active_quad_iterator active_quad_iterator;

    typedef typename TriaDimensionInfo<dim>::raw_hex_iterator raw_hex_iterator;
    typedef typename TriaDimensionInfo<dim>::hex_iterator hex_iterator;
    typedef typename TriaDimensionInfo<dim>::active_hex_iterator active_hex_iterator;

    typedef typename TriaDimensionInfo<dim>::raw_cell_iterator raw_cell_iterator;
    typedef typename TriaDimensionInfo<dim>::cell_iterator cell_iterator;
    typedef typename TriaDimensionInfo<dim>::active_cell_iterator active_cell_iterator;

    typedef typename TriaDimensionInfo<dim>::raw_face_iterator raw_face_iterator;
    typedef typename TriaDimensionInfo<dim>::face_iterator face_iterator;
    typedef typename TriaDimensionInfo<dim>::active_face_iterator active_face_iterator;

				     /**
				      *  Create a triangulation and create
				      *  the first level of the hierarchy.
				      *  Do not create any cells.
				      */
    Triangulation (const MeshSmoothing smooth_grid = none);

				     /**
				      *  Copy constructor. You should really
				      *  use the @p{copy_triangulation} function,
				      *  so we declare this function but let
				      *  it throw an internal error. The
				      *  reason for this is that we may use
				      *  triangulation objects in collection,
				      *  but inside them sometimes strange
				      *  operations like copying happen which
				      *  should be avoided for objects as
				      *  large as triangulations. By throwing
				      *  an error, one easily finds these
				      *  places and can find other ways.
				      */
    Triangulation (const Triangulation<dim> &t);
    
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
				      * any more, such as @ref{DoFHandler} objects
				      * using it.
				      */
    void clear ();
    
				     /**					
				      * Assign a boundary object to
				      * the triangulation. If a face
				      * with boundary number @p{number}
				      * is refined, this object is
				      * used to find the location of
				      * new vertices on the
				      * boundary. This will also be
				      * true for non-linear
				      * transformations of cells to
				      * the unit cell in shape
				      * function
				      * calculations. Multiple calls
				      * to this function are allowed
				      * to store different boundary
				      * curves or surfaces.
				      *
				      * Numbers of boundary objects
				      * correspond to material numbers
				      * of faces at the boundary, for
				      * instance the material id in a
				      * UCD input file. They are not
				      * necessarily consecutive but
				      * must be in the range 0-254.
				      * Material IDs on boundaries are
				      * also called boundary indicators
				      * and are accessed with accessor
				      * functions of that name.
				      *
				      * The @p{boundary_object} is not
				      * copied and MUST persist until
				      * the triangulation is
				      * destroyed. Otherwise, the
				      * @ref{Subscriptor} class will issue
				      * @p{ExcObjectInUse}.  This is
				      * also true for triangulations
				      * generated from this one by
				      * @p{copy_triangulation}.
				      *
				      * It is possible to remove or
				      * replace the boundary object
				      * during the lifetime of a
				      * non-empty
				      * triangulation. Usually, this
				      * is done before the first
				      * refinement and is dangerous
				      * afterwards.  Removal of a
				      * boundary object is done by
				      * @p{set_boundary(number)}, which
				      * uses the default argument of this
				      * function and replaces the boundary
				      * approximation by a piecewise
				      * straight line.
				      */
    void set_boundary (unsigned int         number,
		       const Boundary<dim> &boundary_object = straight_boundary);

				     /**
				      * Return a constant reference to
				      * a boundary object used for
				      * this triangulation.  Number is
				      * the same as in
				      * @p{set_boundary}
				      */
    const Boundary<dim> & get_boundary (unsigned int number) const;
    
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
				      *  the @p{set_boundary}
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
				      *  @p{virtual} since some
				      *  derived classes might want to
				      *  disable or extend the
				      *  functionality of this
				      *  function.
				      */
    virtual void copy_triangulation (const Triangulation<dim> &old_tria);

				     /**
				      * Create a triangulation from a
				      * list of vertices and a list of
				      * cells, each of the latter
				      * being a list of @p{1<<dim}
				      * vertex indices. The
				      * triangulation must be empty
				      * upon calling this function and
				      * the cell list should be useful
				      * (connected domain, etc.).
				      *
				      * Material data for the cells is
				      * given within the @p{cells}
				      * array, while boundary
				      * information is given in the
				      * @p{subcelldata} field.
				      *
				      * The numbering of vertices
				      * within the @p{cells} array is
				      * subject to some constraints;
				      * see the general class
				      * documentation for this.
				      *
				      * This function is made
				      * @p{virtual} to allow derived
				      * classes to set up some data
				      * structures as well.
				      */
    virtual void create_triangulation (const vector<Point<dim> >    &vertices,
				       const vector<CellData<dim> > &cells,
				       const SubCellData            &subcelldata);

				     /**
				      * Distort the grid by randomly
				      * moving around all the vertices
				      * of the grid.  The direction of
				      * moving is random, while the
				      * length of the shift vector has
				      * a value of @p{factor} times
				      * the minimal length of the
				      * active lines adjacent to this
				      * vertex. Note that @p{factor}
				      * should obviously be well below
				      * @p{0.5}.
				      *
				      * If @p{keep_boundary} is set to
				      * @p{true} (which is the
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
				      *  Refine all cells @p{times}
				      *  times, by alternatingly
				      *  calling
				      *  @p{set_all_refine_flags()}
				      *  and
				      *  @p{execute_coarsening_and_refinement()}.
				      *  This function actually starts
				      *  the refinement process, so
				      *  you have no way to store the
				      *  refinement flags unless you
				      *  overload the
				      *  @p{execute_coarsening_and_refinement}
				      *  function.
				      */
    void refine_global (const unsigned int times);

				     /**
				      * Execute both refinement and
				      * coarsening of the
				      * triangulation.
				      *
				      * The function resets all
				      * refinement and coarsening
				      * flags to false. It uses the
				      * user flags.
				      *
                                      * See the general docs for more
                                      * information.
				      *
				      * Note that this function is
				      * @p{virtual} to allow derived
				      * classes to insert hooks, such
				      * as saving refinement flags and
				      * the like.
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
				      * @p{dim>=2} (make sure that no
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
				      * This part of the function is
				      * mostly dimension
				      * independent. However, for some
				      * dimension dependent things, it
				      * calls
				      * @p{prepare_refinement_dim_dependent}.
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
    
				     /*@}*/

				     /**
				      *  @name History of a triangulation
				      */
    				     /*@{*/
				     /**
				      *  Save the addresses of the
				      *  cells which are flagged for
				      *  refinement to @p{out}.  For
				      *  usage, read the general
				      *  documentation for this class.
				      */
    void save_refine_flags (ostream &out) const;

				     /**
				      * Same as above, but store the flags to
				      * a bitvector rather than to a file.
				      */
    void save_refine_flags (vector<bool> &v) const;

				     /**
				      *  Read the information stored by
				      *  @p{save_refine_flags}.
				      */
    void load_refine_flags (istream &in);

    				     /**
				      *  Read the information stored by
				      *  @p{save_refine_flags}.
				      */
    void load_refine_flags (const vector<bool> &v);

				     /**
				      * Analogue to @p{save_refine_flags}.
				      */
    void save_coarsen_flags (ostream &out) const;

				     /**
				      * Same as above, but store the flags to
				      * a bitvector rather than to a file.
				      */
    void save_coarsen_flags (vector<bool> &v) const;

    				     /**
				      * Analogue to @p{load_refine_flags}.
				      */
    void load_coarsen_flags (istream &out);

				     /**
				      * Analogue to @p{load_refine_flags}.
				      */
    void load_coarsen_flags (const vector<bool> &v);

    				     /*@}*/


				     /**
				      *  @name User flag handling
				      */
				     /*@{*/
				     /**
				      *  Clear all user pointers.
				      */
    void clear_user_pointers ();

    				     /**
				      *  Clear all user flags.
				      */
    void clear_user_flags ();

				     /**
				      *  Save all user flags. See the general
				      *  documentation for this class
				      *  and the documentation for the
				      *  @p{save_refine_flags} for more
				      *  details.
				      */
    void save_user_flags (ostream &out) const;

				     /**
				      * Same as above, but store the flags to
				      * a bitvector rather than to a file.
				      * The output vector is resized if
				      * necessary.
				      */
    void save_user_flags (vector<bool> &v) const;

				     /**
				      *  Read the information stored by
				      *  @p{save_user_flags}.
				      */
    void load_user_flags (istream &in);

    				     /**
				      *  Read the information stored by
				      *  @p{save_user_flags}.
				      */
    void load_user_flags (const vector<bool> &v);

				     /**
				      * Save the user flags on lines.
				      */
    void save_user_flags_line (ostream &out) const;

				     /**
				      * Same as above, but store the flags to
				      * a bitvector rather than to a file.
				      * The output vector is resized if
				      * necessary.
				      */
    void save_user_flags_line (vector<bool> &v) const;

				     /**
				      * Load the user flags located on lines.
				      */
    void load_user_flags_line (istream &in);

    				     /**
				      * Load the user flags located on lines.
				      */
    void load_user_flags_line (const vector<bool> &v);

    				     /**
				      * Save the user flags on quads.
				      */
    void save_user_flags_quad (ostream &out) const;

				     /**
				      * Same as above, but store the flags to
				      * a bitvector rather than to a file.
				      * The output vector is resized if
				      * necessary.
				      */
    void save_user_flags_quad (vector<bool> &v) const;

				     /**
				      * Load the user flags located on quads.
				      */
    void load_user_flags_quad (istream &in);

				     /**
				      * Load the user flags located on quads.
				      */
    void load_user_flags_quad (const vector<bool> &v);

				     /**
				      * Save the user flags on hexs.
				      */
    void save_user_flags_hex (ostream &out) const;

				     /**
				      * Same as above, but store the flags to
				      * a bitvector rather than to a file.
				      * The output vector is resized if
				      * necessary.
				      */
    void save_user_flags_hex (vector<bool> &v) const;

				     /**
				      * Load the user flags located on hexs.
				      */
    void load_user_flags_hex (istream &in);

				     /**
				      * Load the user flags located on hexs.
				      */
    void load_user_flags_hex (const vector<bool> &v);
				     /*@}*/

				     /* ------------------------------------ */
    
				     /**
				      *  @name Cell iterator functions
				      */
				     /*@{*/
				     /**
				      *  Iterator to the first cell, used
				      *  or not, on level @p{level}. If a level
				      *  has no cells, a past-the-end iterator
				      *  is returned.
				      *
				      *  This function calls @p{begin_raw_line}
				      *  in 1D and @p{begin_raw_quad} in 2D.
				      */
    raw_cell_iterator    begin_raw   (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first used cell
				      *  on level @p{level}.
				      *
				      *  This function calls @p{begin_line}
				      *  in 1D and @p{begin_quad} in 2D.
				      */
    cell_iterator        begin       (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first active
				      *  cell on level @p{level}.
				      *
				      *  This function calls @p{begin_active_line}
				      *  in 1D and @p{begin_active_quad} in 2D.
				      */
    active_cell_iterator begin_active(const unsigned int level = 0) const;

				     /**
				      *  Iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      *
				      *  This function calls @p{end_line}
				      *  in 1D and @p{end_quad} in 2D.
				      */
    raw_cell_iterator    end () const;

				     /**
				      * Return an iterator which is the first
				      * iterator not on level. If @p{level} is
				      * the last level, then this returns
				      * @p{end()}.
				      */
    cell_iterator        end (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If @p{level} is
				      * the last level, then this returns
				      * @p{end()}.
				      */
    raw_cell_iterator    end_raw (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If @p{level}
				      * is the last level, then this returns
				      * @p{end()}.
				      */
    active_cell_iterator end_active (const unsigned int level) const;


				     /**
				      *  Return an iterator pointing to the
				      *  last cell, used or not.
				      *
				      *  This function calls @p{last_raw_line}
				      *  in 1D and @p{last_raw_quad} in 2D.
				      */
    raw_cell_iterator    last_raw () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  cell of the level @p{level}, used or not.
				      *
				      *  This function calls @p{last_raw_line}
				      *  in 1D and @p{last_raw_quad} in 2D.
				      */
    raw_cell_iterator    last_raw (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used cell.
				      *
				      *  This function calls @p{last_line}
				      *  in 1D and @p{last_quad} in 2D.
				      */
    cell_iterator        last () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used cell on level @p{level}.
				      *
				      *  This function calls @p{last_line}
				      *  in 1D and @p{last_quad} in 2D.
				      */
    cell_iterator        last (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active cell.
				      *
				      *  This function calls @p{last_active_line}
				      *  in 1D and @p{last_active_quad} in 2D.
				      */
    active_cell_iterator last_active () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active cell on level @p{level}.
				      *
				      *  This function calls @p{last_active_line}
				      *  in 1D and @p{last_active_quad} in 2D.
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
				      *  or not, on level @p{level}. If a level
				      *  has no faces, a past-the-end iterator
				      *  is returned.
				      *
				      *  This function calls @p{begin_raw_line}
				      *  in 2D and @p{begin_raw_quad} in 3D.
				      */
    raw_face_iterator    begin_raw_face   (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first used face
				      *  on level @p{level}.
				      *
				      *  This function calls @p{begin_line}
				      *  in 2D and @p{begin_quad} in 3D.
				      */
    face_iterator        begin_face       (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first active
				      *  face on level @p{level}.
				      *
				      *  This function calls @p{begin_active_line}
				      *  in 2D and @p{begin_active_quad} in 3D.
				      */
    active_face_iterator begin_active_face(const unsigned int level = 0) const;

				     /**
				      *  Iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      *
				      *  This function calls @p{end_line}
				      *  in 2D and @p{end_quad} in 3D.
				      */
    raw_face_iterator    end_face () const;

				     /**
				      * Return an iterator which is the first
				      * iterator not on level. If @p{level} is
				      * the last level, then this returns
				      * @p{end()}.
				      */
    face_iterator        end_face (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If @p{level} is
				      * the last level, then this returns
				      * @p{end()}.
				      */
    raw_face_iterator    end_raw_face (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If @p{level}
				      * is the last level, then this returns
				      * @p{end()}.
				      */
    active_face_iterator end_active_face (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the
				      *  last face, used or not.
				      *
				      *  This function calls @p{last_raw_line}
				      *  in 2D and @p{last_raw_quad} in 3D.
				      */
    raw_face_iterator    last_raw_face () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  face of the level @p{level}, used or not.
				      *
				      *  This function calls @p{last_raw_line}
				      *  in 2D and @p{last_raw_quad} in 3D.
				      */
    raw_face_iterator    last_raw_face (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used face.
				      *
				      *  This function calls @p{last_line}
				      *  in 2D and @p{last_quad} in 3D.
				      */
    face_iterator        last_face () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used face on level @p{level}.
				      *
				      *  This function calls @p{last_line}
				      *  in 2D and @p{last_quad} in 3D.
				      */
    face_iterator        last_face (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active face.
				      *
				      *  This function calls @p{last_active_line}
				      *  in 2D and @p{last_active_quad} in 3D.
				      */
    active_face_iterator last_active_face () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active face on level @p{level}.
				      *
				      *  This function calls @p{last_active_line}
				      *  in 2D and @p{last_active_quad} in 3D.
				      */
    active_face_iterator last_active_face (const unsigned int level) const;
				     /*@}*/


				     /*---------------------------------------*/

				     /**
				      *  @name Line iterator functions
				      */
				     /*@{*/
				     /**
				      *  Iterator to the first line, used
				      *  or not, on level @p{level}. If a level
				      *  has no lines, a past-the-end iterator
				      *  is returned.
				      */
    raw_line_iterator
    begin_raw_line   (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first used line
				      *  on level @p{level}.
				      */
    line_iterator
    begin_line       (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first active
				      *  line on level @p{level}.
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
				      * iterator not on level. If @p{level} is
				      * the last level, then this returns
				      * @p{end()}.
				      */
    line_iterator        end_line (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If @p{level} is
				      * the last level, then this returns
				      * @p{end()}.
				      */
    raw_line_iterator    end_raw_line (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If @p{level}
				      * is the last level, then this returns
				      * @p{end()}.
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
				      *  line of the level @p{level}, used or not.

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
				      *  used line on level @p{level}.
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
				      *  active line on level @p{level}.
				      */
    active_line_iterator
    last_active_line (const unsigned int level) const;
				     /*@}*/	  

				     /*---------------------------------------*/

				     /**
				      *  @name Quad iterator functions*/
    				     /*@{
				      */
    				     /**
				      *  Iterator to the first quad, used
				      *  or not, on level @p{level}. If a level
				      *  has no quads, a past-the-end iterator
				      *  is returned.
				      */
    raw_quad_iterator
    begin_raw_quad   (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first used quad
				      *  on level @p{level}.
				      */
    quad_iterator
    begin_quad       (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first active
				      *  quad on level @p{level}.
				      */
    active_quad_iterator
    begin_active_quad(const unsigned int level = 0) const;

				     /**
				      *  Iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      */
    raw_quad_iterator
    end_quad () const;

				     /**
				      * Return an iterator which is the first
				      * iterator not on level. If @p{level} is
				      * the last level, then this returns
				      * @p{end()}.
				      */
    quad_iterator        end_quad (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If @p{level} is
				      * the last level, then this returns
				      * @p{end()}.
				      */
    raw_quad_iterator    end_raw_quad (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If @p{level}
				      * is the last level, then this returns
				      * @p{end()}.
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
				      *  quad of the level @p{level}, used or not.

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
				      *  used quad on level @p{level}.
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
				      *  active quad on level @p{level}.
				      */
    active_quad_iterator
    last_active_quad (const unsigned int level) const;
				     /*@}*/

				     /*---------------------------------------*/

				     /**
				      *  @name Hex iterator functions*/
    				     /*@{
				      */
    				     /**
				      *  Iterator to the first hex, used
				      *  or not, on level @p{level}. If a level
				      *  has no hexs, a past-the-end iterator
				      *  is returned.
				      */
    raw_hex_iterator
    begin_raw_hex   (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first used hex
				      *  on level @p{level}.
				      */
    hex_iterator
    begin_hex       (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first active
				      *  hex on level @p{level}.
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
				      * iterator not on level. If @p{level} is
				      * the last level, then this returns
				      * @p{end()}.
				      */
    hex_iterator        end_hex (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If @p{level} is
				      * the last level, then this returns
				      * @p{end()}.
				      */
    raw_hex_iterator    end_raw_hex (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If @p{level}
				      * is the last level, then this returns
				      * @p{end()}.
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
				      *  hex of the level @p{level}, used or not.

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
				      *  used hex on level @p{level}.
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
				      *  active hex on level @p{level}.
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
				      *  Return total number of used lines,
				      *  active or not.
				      */
    unsigned int n_lines () const;
    
				     /**
				      *  Return total number of used lines,
				      *  active or not on level @p{level}.
				      */
    unsigned int n_lines (const unsigned int level) const;
    
				     /**
				      * Return total number of active lines.
				      */
    unsigned int n_active_lines () const;
    
				     /**
				      *  Return total number of active lines,
				      *  on level @p{level}.
				      */
    unsigned int n_active_lines (const unsigned int level) const;
    
				     /**
				      *  Return total number of used quads,
				      *  active or not.
				      */
    unsigned int n_quads () const;

				     /**
				      *  Return total number of used quads,
				      *  active or not on level @p{level}.
				      */
    unsigned int n_quads (const unsigned int level) const;
    
				     /**
				      *  Return total number of active quads,
				      *  active or not.
				      */
    unsigned int n_active_quads () const;
    
				     /**
				      *  Return total number of active quads,
				      *  active or not on level @p{level}.
				      */
    unsigned int n_active_quads (const unsigned int level) const;
    
				     /**
				      *  Return total number of used hexahedra,
				      *  active or not.
				      */
    unsigned int n_hexs() const;

				     /**
				      *  Return total number of used hexahedra,
				      *  active or not on level @p{level}.
				      */
    unsigned int n_hexs(const unsigned int level) const;
    
				     /**
				      *  Return total number of active hexahedra,
				      *  active or not.
				      */
    unsigned int n_active_hexs() const;
    
				     /**
				      *  Return total number of active hexahedra,
				      *  active or not on level @p{level}.
				      */
    unsigned int n_active_hexs(const unsigned int level) const;

				     /**
				      *  Return total number of used cells,
				      *  active or not.
				      *  Maps to @p{n_lines()} in one space
				      *  dimension and so on.
				      */
    unsigned int n_cells () const;

    				     /**
				      *  Return total number of used cells,
				      *  active or not, on level @p{level}.
				      *  Maps to @p{n_lines(level)} in one space
				      *  dimension and so on.
				      */
    unsigned int n_cells (const unsigned int level) const;

    				     /**
				      *  Return total number of active cells.
				      *  Maps to @p{n_active_lines()} in one space
				      *  dimension and so on.
				      */
    unsigned int n_active_cells () const;

    				     /**
				      * Return total number of active cells
				      * on level @p{level}.
				      * Maps to @p{n_active_lines(level)} in one
				      * space dimension and so on.
				      */
    unsigned int n_active_cells (const unsigned int level) const;

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
				      * Return the total number of vertices.
				      * Some of them may not be used, which
				      * usually happens upon coarsening of
				      * a triangulation when some vertices are
				      * discarded, but we do not want to
				      * renumber the remaining one, leading to
				      * holes in the numbers of used vertices.
				      * You can get the number of used vertices
				      * using @p{n_used_vertices} function.
				      */
    unsigned int n_vertices () const;
    
				     /**
				      * Return the number of vertices that are
				      * presently in use, i.e. belong to at least
				      * one used element.
				      */
    unsigned int n_used_vertices () const;
    
				     /**
				      * Return the maximum number of cells
				      * meeting at a common vertex. Since this
				      * number is an invariant under refinement,
				      * only the cells on
				      * the coarsest level are considered. The
				      * operation is thus reasonably fast. The
				      * invariance is only true for sufficiently
				      * many cells in the coarsest triangulation
				      * (e.g. for a single cell one would be
				      * returned),
				      * so a minimum of four is returned in
				      * two dimensions, 8 in three dimensions,
				      * etc, which is how many cells meet if the
				      * triangulation is refined.
				      *
				      * In one space dimension, two is returned.
				      */
    unsigned int max_adjacent_cells () const;
    				     /*@}*/


				     /**
				      *  @name Exceptions
				      */
				     /*@{*/
				     /**
				      *  Exception
				      */
    DeclException1 (ExcInvalidLevel,
		    int,
		    << "The given level " << arg1
		    << " is not in the valid range!");
				     /**
				      *  Exception
				      */
    DeclException0 (ExcCellShouldBeUnused);
				     /**
				      *  Exception
				      */
    DeclException0 (ExcTooFewVerticesAllocated);
				     /**
				      *  Exception
				      */
    DeclException0 (ExcUncaughtState);
				     /**
				      *  Exception
				      */
    DeclException0 (ExcFunctionNotUseful);
				     /**
				      *  Exception
				      */
    DeclException2 (ExcGridsDoNotMatch,
		    int, int,
		    << "The present grid has " << arg1 << " active cells, "
		    << "but the one in the file had " << arg2);
				     /**
				      *  Exception
				      */
    DeclException0 (ExcGridReadError);
				     /**
				      *  Exception
				      */
    DeclException0 (ExcTriangulationNotEmpty);
				     /**
				      *  Exception
				      */
    DeclException1 (ExcGridHasInvalidCell,
		    int,
		    << "Something went wrong when making cell " << arg1
		    << ". Read the docs and the source code "
		    << "for more information.");
				     /**
				      *  Exception
				      */
    DeclException0 (ExcGridHasInvalidVertices);
				     /**
				      *  Exception
				      */
    DeclException1 (ExcInternalErrorOnCell,
		    int,
		    << "Something went wrong upon construction of cell "
		    << arg1);
				     /**
				      * Exception
				      */
    DeclException3 (ExcInvalidVertexIndex,
		    int, int, int,
		    << "Error while creating cell " << arg1
		    << ": the vertex index " << arg2 << " must be between 0 and "
		    << arg3 << ".");
				     /**
				      * Exception
				      */
    DeclException2 (ExcLineInexistant,
		    int, int,
		    << "When trying to give a boundary indicator to a line: "
		    << "the line with end vertices " << arg1 << " and "
		    << arg2 << " does not exist.");
				     /**
				      * Exception
				      */
    DeclException0 (ExcInteriorLineCantBeBoundary);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInternalError);
				     /**
				      * Exception
				      */
    DeclException1 (ExcEmptyLevel,
		    int,
		    << "You tried to do something on level " << arg1
		    << ", but this level is empty.");
				     /**
				      * Exception
				      */
    DeclException0 (ExcIO);
				     /*@}*/
  protected:

				     /**
				      *  Write a bool vector to the given stream,
				      *  writing a pre- and a postfix magica
				      *  number. The vector is written in an
				      *  almost binary format, i.e. the bool
				      *  flags are packed but the data is written
				      *  as ASCII text.
				      *
				      *  The flags are stored in a binary
				      *  format: for each @p{true}, a @p{1}
				      *  bit is stored, a @p{0} bit otherwise.
				      *  The bits are stored as @p{unsigned char},
				      *  thus avoiding endianess. They are
				      *  written to @p{out} in plain text, thus
				      *  amounting to 3.6 bits per active cell
				      *  on the overage. Other information
				      *  (magic numbers ans number of elements)
				      *  is stored as plain text
				      *  as well. The format should therefore be
				      *  interplatform compatible.
				      */
    static void write_bool_vector (const unsigned int  magic_number1,
				   const vector<bool> &v,
				   const unsigned int  magic_number2,
				   ostream            &out);

				     /**
				      * Re-read a vector of bools previously
				      * written by @p{write_bool_vector} and
				      * compare with the magic numbers.
				      */
    static void read_bool_vector (const unsigned int  magic_number1,
				  vector<bool>       &v,
				  const unsigned int  magic_number2,
				  istream            &in);
    
  private:
				     /**
				      *  Refine all cells on all levels which
				      *  were previously flagged for refinement.
				      */ 
    void execute_refinement ();

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
				      * Actually delete a cell, or rather all
				      * it children, which is the
				      * main step for the coarsening process.
				      * This is the dimension dependent part
				      * of @p{execute_coarsening}.
				      */
    void delete_children (cell_iterator &cell);

				     /**
				      * Some dimension dependent stuff for
				      * mesh smoothing.
				      *
				      * At present, this function does nothing
				      * in 1d and 2d, but makes sure no two
				      * cells with a level difference greater
				      * than one share one line in 3d. This
				      * is a requirement needed for the
				      * interpolation of hanging nodes, since
				      * otherwise to steps of interpolation
				      * would be necessary. This would make
				      * the processes implemented in the
				      * @p{ConstraintMatrix} class much more
				      * complex, since these two steps of
				      * interpolation do not commute.
				      */
    void prepare_refinement_dim_dependent ();

				     /**
				      * Make sure that either all or none of
				      * the children of a cell are tagged for
				      * coarsening.
				      */
    void fix_coarsen_flags ();

				     /**
				      * Re-compute the number of lines, quads,
				      * etc. This function is called by
				      * @p{execute_{coarsening,refinement}} and
				      * by @p{create_triangulation} after the
				      * grid was changed.
				      *
				      * This function simply delegates to the
				      * functions below, which count
				      * only a certain class of objects.
				      */
    void update_number_cache ();

				     /**
				      * Recompute the number of lines.
				      */
    void update_number_cache_lines ();

    				     /**
				      * Recompute the number of quads.
				      */
    void update_number_cache_quads ();

    				     /**
				      * Recompute the number of hexes.
				      */
    void update_number_cache_hexes ();


				     /**
				      *  Array of pointers pointing to the
				      *  @p{TriangulationLevel<dim>} objects
				      *  storing the data on the different
				      *  levels.
				      *
				      *  Usage is like @p{levels[3]->quads}.
				      */
    vector<TriangulationLevel<dim>*> levels;

				     /**
				      *  Array of the vertices of this
				      *  triangulation.
				      */
    vector<Point<dim> >              vertices;

				     /**
				      *  Array storing a bit-pattern which
				      *  vertices are used.
				      */
    vector<bool>                     vertices_used;

				     /**
				      *  Collection of boundary objects.
				      */
    const Boundary<dim>* boundary[255];

				     /**
				      *  Do some smoothing in the process
				      *  of refining the triangulation. See
				      *  the general doc of this class for
				      *  more information about this.
				      */
    MeshSmoothing                    smooth_grid;

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
    TriaNumberCache<dim>             number_cache;
    
				     // Friendship includes local classes.
#if (__GNUC__==2) && (__GNUC_MINOR__ < 95)
				     // this seems to be disallowed
				     // by the standard, so gcc2.95
				     // does not accept it
    friend class TriaAccessor<dim>;
    friend class TriaObjectAccessor<1, dim>;
    friend class TriaObjectAccessor<2, dim>;
    friend class TriaObjectAccessor<3, dim>;
#else
				     // this, however, may grant
				     // access to too many classes,
				     // but ...
    template <int dim1> friend class TriaAccessor;

    template <int dim1, int dim2>
    friend class TriaObjectAccessor;
#endif
    
    friend class CellAccessor<dim>;
    
    friend class TriaRawIterator<1,TriaObjectAccessor<1, 1> >;
    friend class TriaRawIterator<1,CellAccessor<1> >;

    friend class TriaRawIterator<2,TriaObjectAccessor<1, 2> >;
    friend class TriaRawIterator<2,TriaObjectAccessor<2, 2> >;
    friend class TriaRawIterator<2,CellAccessor<2> >;

    friend class TriaRawIterator<3,TriaObjectAccessor<1, 3> >;
    friend class TriaRawIterator<3,TriaObjectAccessor<2, 3> >;
    friend class TriaRawIterator<3,TriaObjectAccessor<3, 3> >;
    friend class TriaRawIterator<3,CellAccessor<3> >;

    friend class DoFHandler<dim>;
    friend class MGDoFHandler<dim>;
};


#endif



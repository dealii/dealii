/*----------------------------   tria.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __tria_H
#define __tria_H
/*----------------------------   tria.h     ---------------------------*/

#include <vector>
#include <grid/point.h>
#include <grid/geometry_info.h>



//forward declaration needed
template <int dim> class Boundary;

template <int dim> class TriaAccessor;
template <int dim> class LineAccessor;
template <int dim> class QuadAccessor;
template <int dim> class CellAccessor;

template <int dim> class TriangulationLevel;

template <int dim, class Accessor> class TriaRawIterator;
template <int dim, class Accessor> class TriaIterator;
template <int dim, class Accessor> class TriaActiveIterator;

template <int dim> class DoFHandler;
template <int dim> class MGDoFHandler;

template <int dim> struct CellData;
struct SubCellData;

class dVector;


class istream;
class ostream;



/**
 * Declare some symbolic names for mesh smoothing algorithms. The meaning of
 * these flags is documented in the #Triangulation# class.
 */
enum MeshSmoothing {
      none                               = 0x0,
      limit_level_difference_at_vertices = 0x1,     
      eliminate_unrefined_islands        = 0x2,

      eliminate_refined_inner_islands    = 0x100,
      eliminate_refined_boundary_islands = 0x200,
      
      smoothing_on_refinement            = (limit_level_difference_at_vertices |
					    eliminate_unrefined_islands),
      smoothing_on_coarsening            = (eliminate_refined_inner_islands |
					    eliminate_refined_boundary_islands),
      maximum_smoothing                  = 0xffff
};

				       



    






/*------------------------------------------------------------------------*/


/**
 *  Structure which is passed to the #Triangulation<dim>::create_triangulation#
 *  function. It contains all data needed to construct a cell, namely the
 *  indices of the vertices and the material indicator.
 */
template <int dim>
struct CellData {
    int           vertices[GeometryInfo<dim>::vertices_per_cell];
    unsigned char material_id;
};





/**
 *  Structure to be passed to the #Triangulation<dim>::create_triangulation#
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
				      * whether the #boundary_*# arrays are
				      * empty when in one space dimension
				      * and whether the #boundary_quads#
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
 *  This class implements some types which differ between the dimensions.
 *  Declare it to have a template parameter, but do not actually declare
 *  other types than those explicitely instantiated.
 */
template <int dim>
class TriaDimensionInfo;



/**
 *  This class implements some types which differ between the dimensions.
 *
 *  A #line_iterator# is typdef'd to an iterator operating on the
 *  #lines# member variable of a #Triangulation<1># object. An
 *  #active_line_iterator# only operates on the active lines.
 *  #raw_line_iterator# objects operate on all lines, used or not.
 *
 *  Since we are in one dimension, the following identities are declared:
 *  \begin{verbatim}
 *    typedef raw_line_iterator    raw_cell_iterator;
 *    typedef line_iterator        cell_iterator;
 *    typedef active_line_iterator active_cell_iterator;
 *  \end{verbatim}
 *
 *  To enable the declaration of #begin_quad# and the like in
 *  #Triangulation<1>#, the #quad_iterator#s are declared as
 *  #void *#. Thus these types exist, but are useless and will
 *  certainly make any involuntary use visible.
 *
 *  The same applies for the #face_iterator# types, since lines
 *  have no substructures apart from vertices, which are handled in
 *  a different way, however.
 */
template <>
class TriaDimensionInfo<1> {
  public:
    typedef TriaRawIterator<1,CellAccessor<1> >    raw_line_iterator;
    typedef TriaIterator<1,CellAccessor<1> >       line_iterator;
    typedef TriaActiveIterator<1,CellAccessor<1> > active_line_iterator;

    typedef void * raw_quad_iterator;
    typedef void * quad_iterator;
    typedef void * active_quad_iterator;

    typedef raw_line_iterator    raw_cell_iterator;
    typedef line_iterator        cell_iterator;
    typedef active_line_iterator active_cell_iterator;

    typedef void * raw_face_iterator;
    typedef void * face_iterator;
    typedef void * active_face_iterator;
};


/**
 *  This class implements some types which differ between the dimensions.
 *
 *  A #line_iterator# is typdef'd to an iterator operating on the
 *  #lines# member variable of a #Triangulation<1># object. An
 *  #active_line_iterator# only operates on the active lines.
 *  #raw_line_iterator# objects operate on all lines, used or not.
 *  Using #active_line_iterator#s may not be particularly useful since it
 *  only operates on unrefined lines. However, also refined lines may bound
 *  unrefined cells if the neighboring cell is refined once more than the
 *  present one.
 *
 *  Similarly, #quad_iterator#, #raw_quad_iterator# and
 *  #active_quad_iterator# are declared.
 *  
 *  Since we are in two dimension, the following identities are declared:
 *  \begin{verbatim}
 *    typedef raw_quad_iterator    raw_cell_iterator;
 *    typedef quad_iterator        cell_iterator;
 *    typedef active_quad_iterator active_cell_iterator;
 *
 *    typedef raw_line_iterator    raw_face_iterator;
 *    typedef line_iterator        face_iterator;
 *    typedef active_line_iterator active_face_iterator;    
 *  \end{verbatim}
 */
template <>
class TriaDimensionInfo<2> {
  public:
    typedef TriaRawIterator<2,LineAccessor<2> >    raw_line_iterator;
    typedef TriaIterator<2,LineAccessor<2> >       line_iterator;
    typedef TriaActiveIterator<2,LineAccessor<2> > active_line_iterator;
    
    typedef TriaRawIterator<2,CellAccessor<2> >    raw_quad_iterator;
    typedef TriaIterator<2,CellAccessor<2> >       quad_iterator;
    typedef TriaActiveIterator<2,CellAccessor<2> > active_quad_iterator;

    typedef raw_quad_iterator    raw_cell_iterator;
    typedef quad_iterator        cell_iterator;
    typedef active_quad_iterator active_cell_iterator;

    typedef raw_line_iterator    raw_face_iterator;
    typedef line_iterator        face_iterator;
    typedef active_line_iterator active_face_iterator;    
};







/*------------------------------------------------------------------------*/


/**
 *  #Triangulation#s denote a hierarchy of levels of elements which together
 *  form a region in #dim# spatial dimensions.
 *
 *  This class is written to be as independent of the dimension as possible
 *  (thus the complex construction of the #TriangulationLevel# classes) to
 *  allow code-sharing, to allow reducing the need to mirror changes in the code
 *  for one dimenion to the code for other dimensions. Nonetheless, some of
 *  the functions are dependent of the dimension and there only exist
 *  specialized versions for distinct dimensions.
 *
 *
 *  \subsection{Structure and iterators}
 *
 *  The actual data structure of a #Triangulation# object is rather complex
 *  and quite inconvenient if one attempted to operate on it directly, since
 *  data is spread over quite a lot of arrays and other places. However,
 *  there are ways powerful enough to work on these data structures
 *  without knowing their exact relations. This is done through the
 *  concept of iterators (see the STL documentation and \Ref{TriaRawIterator}).
 *  In order to make things as easy and dimension independent as possible,
 *  use of class local typedefs is made, see below.
 *  
 *  In the base class #TriaDimensionInfo#, a #Cell# is typedef'd to be whatever
 *  is reasonable for a cell in the respective dimension, i.e. a #Line# in
 *  one dimension, a #Quad# in two dimensions, and so on.
 *
 *  The #Triangulation# class provides iterator which enable looping over all
 *  lines, cells,
 *  etc without knowing the exact representation used to describe them. Their
 *  names are typedef's in the #TriaDimensionInfo# base class (thus making them
 *  local types to this class) and are as follows:
 *
 *  #raw_line_iterator#: loop over all lines, used or not (declared for
 *  all dimensions).
 *  
 *  #line_iterator#: loop over all used lines (declared for all dimensions).
 *
 *  #active_line_iterator#: loop over all active lines (declared for all
 *  dimensions).
 *
 *  #raw_quad_iterator#: loop over all quads, used or not (declared only
 *  for #dim>=2#).
 *  
 *  #quad_iterator#: loop over all quads (declared only for #dim#>=2).
 *
 *  #active_quad_iterator#: loop over all active quads (declared only for
 *  #dim#>=2).
 *
 *  Additionaly, for #dim#==1, the following identities hold:
 *  \begin{verbatim}
 *    typedef raw_line_iterator    raw_cell_iterator;
 *    typedef line_iterator        cell_iterator;
 *    typedef active_line_iterator active_cell_iterator;
 *  \end{verbatim}
 *  while for #dim#==2
 *  \begin{verbatim}
 *    typedef quad_line_iterator   raw_cell_iterator;    
 *    typedef quad_iterator        cell_iterator;
 *    typedef active_quad_iterator active_cell_iterator;
 *
 *    typedef raw_line_iterator    raw_face_iterator;
 *    typedef line_iterator        face_iterator;
 *    typedef active_line_iterator active_face_iterator;    
 *  \end{verbatim}
 *
 *  By using the cell iterators, you can write code nearly independent of
 *  the spatial dimension. The same applies for substructure iterators,
 *  where a substructure is defined as a face of a cell. The face of a
 *  cell is be a vertex in 1D and a line in 2D; however, vertices are
 *  handled in a different way and therefore lines have no faces.
 *
 *  The #Triangulation# class offers functions like #begin_active# which gives
 *  you an iterator to the first active cell. There are quite a lot of functions
 *  returning iterators. Take a look at the class doc to get an overview.
 *
 *  Usage of these iterators works mostly like with the STL iterators. Some
 *  examples taken from the #Triangulation# source code follow.
 *  \begin{itemize}
 *  \item {\it Counting the number of cells on a specific level}
 *    \begin{verbatim}
 *     template <int dim>
 *     int Triangulation<dim>::n_cells (const int level) const {
 *        cell_iterator cell = begin (level),
 *                      endc = end(level);
 *        int n=0;
 *        for (; cell!=endc; ++cell)
 *          ++n;
 *        return n;
 *      };
 *    \end{verbatim}
 *    Another way which uses the STL #distance# function would be to write
 *    \begin{verbatim}
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
 *    \end{verbatim}
 *    Unfortunately, #g++# presently (version 2.7.2) fails to find the right
 *    #distance# template instantiation; it seems we have to wait for future
 *    #g++# versions :-(
 *    
 *  \item {\it Refining all cells of a triangulation}
 *    \begin{verbatim}
 *      template <int dim>
 *      void Triangulation<dim>::refine_global () {
 *        active_cell_iterator cell = begin_active(),
 *                             endc = end();
 *
 *        for (; cell != endc; ++cell)
 *          cell->set_refine_flag ();
 *        execute_refinement ();
 *      };
 *    \end{verbatim}
 *  \end{itemize}
 *
 *
 *  \subsection{Usage}
 *
 *  Usage of a #Triangulation# is mainly done through the use of iterators.
 *  An example probably shows best how to use it:
 *  \begin{verbatim}
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
 *    tria.execute_refinement ();
 *
 *                                     // refine first active cell
 *                                     // on coarsest level
 *    tria.begin_active()->set_refine_flag ();
 *    tria.save_refine_flags (history);
 *    tria.execute_refinement ();
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
 *        tria.execute_refinement ();
 *      };
 *                                       // output the grid
 *    ofstream out("grid.1");
 *    tria.print_gnuplot (out);
 *  };  
 *  \end{verbatim}
 *
 *  
 *  \subsection{Creating a triangulation}
 *
 *  There are several possibilities to create a triangulation:
 *  \begin{itemize}
 *    \item Hypercube triangulations: a hypercube triangulation is a
 *       domain which is the tensor product of an interval $[a,b]$ in
 *       the given number of spatial dimensions. If you want to create such
 *       a domain, which is a common test case for model problems, call
 *       #Triangulation<dim>::create_hypercube (a,b)#, which produces a
 *       hypercube domain triangulated with exactly one element. You can
 *       get tensor product meshes by successive refinement of this cell.
 *
 *    \item Other standard regions: you can get the generalized L-shape domain
 *      using the #Triangulation<dim>::create_L_region (a,b)# function, which
 *      is the hypercube with the interval $[a,b]$ without the hypercube
 *      made out of the interval $[(a+b)/2,b]$. Let, for example, be $a=-1$
 *      and $b=1$, then the hpyer-L in two dimensions is the region
 *      $[-1,1]^2 - [0,1]^2$. To create a hyper-L in one dimension results in
 *      an error.
 *
 *      You get the circle or ball (or generalized: hyperball) around origin
 *      #p# and with radius #r# by calling
 *      #Triangulation<dim>::create_hyper_ball (p, r)#. The circle is triangulated
 *       by five cells, the ball by seven cells. The diameter of the center cell is
 *       chosen so that the aspect ratio of the boundary cells after one refinement
 *       is minimized in some way. To create a hyperball in one dimension results in
 *       an error.
 *
 *       Do not forget to attach a suitable
 *       boundary approximation object if you want the triangulation to be refined
 *       at the outer boundaries.
 *   
 *     \item Reading in a triangulation: By using an object of the \Ref{#DataIn#}
 *        class, you can read in fairly general triangulations. See there for
 *        more information. The mentioned class uses the interface described
 *        directly below to transfer the data into the triangulation.
 *
 *     \item Explicitely creating a triangulation: you can create a triangulation
 *        by providing a list of vertices and a list of cells. Each such cell
 *        consists of a vector storing the indices of the vertices of this cell
 *        in the vertex list. To see how this works, you can take a look at the
 *        #DataIn<dim>::read_*# functions. The appropriate function to be
 *        called is #Triangulation<dim>::create_triangulation (2)#.
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
 *        \Ref{DataIn} class for more details on this. They do not only
 *        hold for the data read from an UCD or any other input file, but
 *        also for the data passed to the
 *        #Triangulation<dim>::create_triangulation (2)# function.
 *
 *     \item Copying a triangulation: when computing on time dependant meshes
 *        of when using adaptive refinement, you will often want to create a
 *        new triangulation to be the same as another one. This is facilitated
 *        by the #copy_triangulation# function.
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
 *   \end{itemize}
 *
 *   The material id for each cell must be specified upon construction of
 *   a triangulation. (There is a special section on material ids and
 *   boundary indicators. See there for more information.)
 *   The standard region functions (for hypercube, hyperball,
 *   etc.) denote all cells the material id zero. You may change that afterwards,
 *   but you should not use the material id 255. When reading a triangulation,
 *   the material id must be specified in the input file (UCD format) or is
 *   otherwise set to zero. When creating explicitely, the material id must
 *   be given to the creation function.
 *
 *   Regarding the boundary indicator for lines in two dimensions and quads
 *   in three (subsummed by the word "faces"), all interior faces are denoted
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
 *   #Triangulation<dim>::distort_random#. It moves all the vertices in the
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
 *   usual behaviour with cells at the boundary, here you may get into trouble
 *   when using multigrid algorithms or when transferring solutions from coarse
 *   to fine grids and back. In general, the use of this function is only safe
 *   if you only use the most refined level of the triangulation for
 *   computations.
 *
 *
 *
 *   \subsection{Refinement and coarsening of a triangulation}
 *
 *   Refinement of a triangulation may be done through several ways. The most
 *   low-level way is directly through iterators: let #i# be an iterator to
 *   an active cell (i.e. the cell pointed to has no children), then the
 *   function call #i->set_refine_flag()# marks the respective cell for
 *   refinement. Marking non-active cells results in an error.
 *
 *   After all the cells you wanted to mark for refinement, call the
 *   #execute_refinement# function to actually perform the refinement. This
 *   function itself first calls the #prepare_refinement# function to regularise
 *   the resulting triangulation: since a face between two adjacent cells may
 *   only be subdivided once (i.e. the levels of two adjacent cells may
 *   differ by one at most; it is not possible to have a cell refined twice
 *   while the neighboring one is not refined), some additional cells are
 *   flagged for refinement to smooth the grid. This enlarges the number of
 *   resulting cells but makes the grid more regular, thus leading to better
 *   approximationonal properties and, above all, making the handling of data
 *   structures and algorithms much much easier. To be honest, this is mostly
 *   an algorithmic step than one needed by the finite element method.
 *
 *   To coarsen a grid, the same way as above is possible by using
 *   #i->set_coarsen_flag# and calling #execute_coarsening#. You can use
 *   #execute_coarsening_and_refinement# to get both actions done, first
 *   coarsening and refinement. The reason for this order is that the
 *   refinement usually adds some additional cells to keep the triangulation
 *   regular and thus satifies all refinement requests, while the coarsening
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
 *   elements.
 *
 *   The central function to this is
 *   #refine (const dVector &criterion, const double threshold)#: it takes a
 *   vector of values, one per active cell, which denote the criterion according
 *   to which the triangulation is to be refined. It marks all cells for which
 *   the criterion is greater than the threshold being given as the second
 *   argument. Analogously,
 *   #coarsen (const dVector &criterion, const double threshold)# flags those
 *   cells for coarsening for which the criterion is less than the treshold.
 *
 *   There are two variations of these functions, which rely on #refine# and
 *   coarsen by computing the thresholds from other information:
 *   \begin{itemize}
 *   \item #refine_and_coarsen_fixed_number#: this function takes a vector as
 *     above and two values between zero and one denoting the fractions of cells to
 *     be refined and coarsened. For this purpose, it sorts the criteria per cell
 *     and takes the threshold to be the one belonging to the cell with the
 *     #fraction times n_active_cells# highest criterion. For example, if
 *     the fraction is $0.3$, the threshold is computed to a value such that
 *     30 per cent of cells have a criterion higher than the threshold and are
 *     thus flagged for refinement. The flagging for refinement is done through
 *     the central #refine# function. For coarsening, the same holds.
 *
 *     The sorting of criteria is not done actually, since we only need one
 *     value, in the example above the criterion of the cell which is at
 *     30 per cent in the sorted list of cells. The order of cells with higher
 *     and of those with lower criteria is irrelevant. Getting this value is
 *     accomplished by the #nth_element# function of the #C++# standard
 *     library, which takes only linear time in the number of elements, rather
 *     than #N log N# for sorting all values.
 *
 *     A typical value for the fraction of cells to be refined is 0.3.
 *     However, for singular functions or singular error functionals, you may
 *     want to chose a smaller value to avoid overrefinement in regions which
 *     do not contribute much to the error.
 *
 *   \item #refine_and_coarsen_fixed_fraction#: this function computes the
 *     threshold such that the number of cells getting flagged for refinement
 *     makes up for a certain fraction of the total error. If this fraction is 50
 *     per cent, for example, the threshold is computed such that the cells with
 *     a criterion greater than the threshold together account for half of the
 *     total error. The definition of the fraction is a bit unintuitive, since
 *     the total error is the sum over all cells of the local contribution
 *     squared. We define that the fraction $\alpha$ be such that those
 *     elements with the greatest error are refined for which the condition
 *     $\sum \eta_K^2 \le \alpha\eta^2$ holds. Note that $\alpha$ is not
 *     squared. The sum runs over the mentioned
 *     cells, $\eta_K$ are the local error indicators and $\eta$ is the global
 *     indicator with $\eta^2 = \sum \eta_K^2$, with here the sum running over
 *     all cells.
 *
 *     For the bottom fraction the same holds: the treshold for coarsening is
 *     computed such that the cells with criterion less than the threshold
 *     together make up for the fraction of the total error specified.
 *
 *     This strategy is more suited for singular functions and error
 *     functionals, but may lead to very slow convergence of the grid
 *     if only few cells are refined in each step.
 *
 *     From the implementational point, this time we really need to
 *     sort the array of criteria.
 *     Just like the other strategy described above, this function only
 *     computes the threshold values and then passes over to #refine# and
 *     #coarsen#.
 *
 *     A typical value for the fraction of the total error is 0.5.
 *   \end{itemize}
 *
 *   For a more thorough discussion of advantages and disadvantages of the
 *   different strategies for refinement, see the paper of R. Becker and
 *   R. Rannacher titled "A Feed-Back Approach to Error Control in Finite
 *   Element Methods: Basic Analysis and Examples".
 *
 *   It is assumed that the criterion is a value in a certain norm over each
 *   element, such that the square of the total error is the sum over the
 *   squares of the criteria on the cells. The criteria shall be positive.
 *
 *   You can suppress coarsening or refining by giving zero as the fraction
 *   for one of the operations.
 *
 *
 *   \subsection{Smoothing of a triangulation}
 *
 *   Some degradation of approximation properties has been observed
 *   for grids which are too unstructured. Therefore, the #prepare_refinement#
 *   function which is automatically called by the #execute_refinement# function
 *   can do some smoothing of the triangulation. The same holds for the
 *   #prepare_coarsening# function called by #execute_coarsening#.Note that
 *   mesh smoothing is only done for two or more space dimensions, no smoothing
 *   is available at present for one spatial dimension. In the sequel,
 *   let #execute_*# stand for any of #execute_refinement#, #execute_coarsening#
 *   or #execute_refinement_and_coarsening#.
 *
 *   For the purpose of smoothing, the
 *   #Triangulation# constructor takes an argument specifying whether a
 *   smoothing step shall be performed on the grid each time #execute_*#
 *   is called. The default is that such a step not be done, since this results
 *   in additional cells being produced, which may not be necessary in all
 *   cases. If switched on, calling #execute_*# results in
 *   flagging additional cells for refinement to avoid
 *   vertices as the ones mentioned. The algorithms for both regularisation
 *   and smoothing of triangulations are described below in the section on
 *   technical issues. The reason why this parameter must be given to the
 *   constructor rather than to #execute_*# is that it would result
 *   in algorithmic problems if you called #execute_*# once without
 *   and once with smoothing, since then in some refinement steps would need
 *   to be refined twice.
 *
 *   The parameter taken by the constructor is an integer which may be composed
 *   bitwise by the constants defined in the #enum MeshSmoothing#. The meaning
 *   of these constants is explained in the following:
 *   \begin{itemize}
 *   \item #limit_level_difference_at_vertices#:
 *     It can be shown, that degradation of approximation occurs if the
 *     triangulation contains vertices which are member of cells with levels
 *     differing by more than one. One such example is the following:
 *     \begin{verbatim}
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
 *     \end{verbatim}
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
 *     element solution (asuming bilinear finite elements), one sees that the
 *     numerical solution is almost totally wrong, compared with the true second
 *     derivative. Indeed, on regular meshes, there exist sharp estimations that
 *     the $H^2$-error is only $O(1)$, so we should not be suprised; however, the
 *     numerical solution may show a value for the second derivative which may
 *     be a factor of ten away from the true value. These problems are located
 *     on the small cell adjacent to the center vertex, where cells of
 *     non-subsequent levels meet, as well as on the upper and right neighbor
 *     of this cell (but with a less degree of deviation from the true value).
 *
 *     If the smoothing indicator given to the constructor contains the bit for
 *     #limit_level_difference_at_vertices#, situations as the above one are
 *     eliminated by also marking the lower left cell for refinement.
 *
 *   \item #eliminate_unrefined_islands#:
 *     Single cells which are not refined and are surrounded by cells which are
 *     refined usually also lead to a sharp decline in approximation properties
 *     locally. The reason is that the nodes on the faces between unrefined and
 *     refined cells are not real degrees of freedom but carry constraints. The
 *     patch without additional degrees of freedom is thus significantly larger
 *     then the unrefined cell itself. If in the parameter passed to the
 *     constructor the bit for #eliminate_unrefined_islands# is set, all cells
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
 *   \item #eliminate_refined_*_islands#:
 *     This algorithm seeks for isolated cells which are refined or flagged
 *     for refinement. This definition is unlike that for
 *     #eliminate_unrefined_islands#, which would mean that an island is
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
 *     #eliminate_refined_inner_islands# and #eliminate_refined_boundary_islands#.
 *     There first eliminates islands defined by the definition above which are
 *     in the interior of the domain, while the second eliminates only those
 *     islands if the cell is at the boundary. The reason for this split of
 *     flags is that one often wants to eliminate such islands in the interior
 *     while those at the boundary may well be wanted, for example if one
 *     refines the mesh according to a criterion associated with a boundary
 *     integral or if one has rough boundary data.
 *
 *   \item #smoothing_on_refinement#:
 *     This flag sums up all smoothing algorithms which may be performed upon
 *     refinement by flagging some more cells for refinement.
 *
 *   \item #smoothing_on_coarsening#:
 *     This flag sums up all smoothing algorithms which may be performed upon
 *     coarsening by flagging some more cells for coarsening.
 *
 *   \item #maximum_smoothing#:
 *     This flag includes all the above ones and therefore combines all
 *     smoothing algorithms implemented.
 *
 *   \item #none#:
 *     Select no smoothing at all.
 *   \end{itemize}
 *
 *
 *   \subsection{Material and boundary information}
 *
 *   Each line, quad, etc stores one byte of information denoting the material
 *   a cell is made of (used in the computation of stiffness matrices) and to
 *   which part of the boundary it belongs. Obviously, the material id is what
 *   is needed for a cell, while for all structures with a dimension less than
 *   the dimension of the domain (i.e. lines in 2D, lines and quads in 3D), the
 *   boundary information is what is needed. Since either material or boundary
 *   information is needed, but never both, only one field is used to store this
 *   data, namely the #TriangulationLevel<1>::LinesData.material_id# and
 *   #TriangulationLevel<2>::QuadsData.material_id# vectors. The data can be
 *   read and written using line, quad and cell iterators.
 *
 *   Material and boundary indicators are stored as one byte (an
 *   #unsigned char#). They can therefore be in the range zero to 255, but
 *   only zero to 254 is allowed. The value 255 is reserved to denote
 *   interior lines (in 2D) and interior lines and quads (in 3D), which need
 *   not have a boundary or material indicator. However, using this value, it
 *   is possible to say whether a line in 2D is interior or not, which would
 *   otherwise be impossible because the hierarchical structure of a
 *   triangulation stores neighborship information and the like only with
 *   cells. Finding out whether a line is an interior one would then only be
 *   possible by looking at the cell it belongs to. There would be no way to
 *   loop over all lines and for example do a contour integral, since there
 *   would be no way to find out which of the lines we loop over actually are
 *   on the contour.
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
 *   \subsection{History of a triangulation}
 *
 *   It is possible to reconstruct a grid from its refinement history, which
 *   can be stored and loaded through the #save_refine_flags# and
 *   #load_refine_flags# functions. Normally, the code will look like this:
 *   \begin{verbatim}
 *                                 // open output file
 *     ofstream history("mesh.history");
 *                                 // do 10 refinement steps
 *     for (int step=0; step<10; ++step) {
 *       ...;
 *       // flag cells according to some criterion
 *       ...;
 *       tria.save_refine_flags (history);
 *       tria.execute_refinement ();
 *     };        
 *   \end{verbatim}
 *
 *   If you want to re-create the grid from the stored information, you write:
 *   \begin{verbatim}
 *                                 // open input file
 *     ifstream history("mesh.history");
 *                                 // do 10 refinement steps
 *     for (int step=0; step<10; ++step) {
 *       tria.load_refine_flags (history);
 *       tria.execute_refinement ();
 *     };        
 *   \end{verbatim}
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
 *   There actually are two sets of #save_*_flags# and #load_*_flags# functions.
 *   One takes a stream as argument and reads/writes the information from/to the
 *   stream, thus enabling storing flags to files. The other set takes an
 *   argument of type #vector<bool>#. This enables the user to temporarily store
 *   some flags, e.g. if another function needs them, and restore them
 *   afterwards.
 *
 *
 *   \subsection{User flags}
 *
 *   A triangulation offers one bit per line, quad, etc for user data.
 *   This field can be
 *   accessed as all other data using iterators. Normally, this user flag is
 *   used if an algorithm walks over all cells and needs information whether
 *   another cell, e.g. a neighbor, has already been processed. It can also
 *   be used to flag the lines subject to constraints in 2D, as for example the
 *   functions in the #DoFHandler# classes do.
 *
 *   There are two functions, #save_user_flags# and #load_user_flags# which
 *   write and read these flags to and from a stream. Unlike
 *   #save_refine_flags# and #load_refine_flags#, these two functions store
 *   and read the flags of all used lines, quads, etc, not only of the
 *   active ones (well, activity is a concept which really only applies to
 *   cells, not for example to lines in 2D, so the abovementioned generalisation
 *   to {\it all} lines, quads, etc seems plausible).
 *
 *   If you want to store more specific user flags, you can use the functions
 *   #save_user_flags_line# and #load_user_flags_line# and the generalizations
 *   for quads, etc.
 *
 *   As for the refinement and coarsening flags, there exist two versions of these
 *   functions, one which reads/writes from a stream and one which does so from
 *   a #vector<bool>#. The latter is used to store flags temporarily, while the
 *   first is used to store them in a file.
 *
 *   It is convention to clear the user flags using the
 *   #Triangulation<>::clear_user_flags()# function before usage, since it is
 *   often necessary to use the flags in more than one function consecutively and
 *   is then error prone to dedicate one of these to clear the flags.
 *
 *   It is recommended that a functions using the flags states so in its
 *   documentation.
 *
 *   There is another set of user data, namely a #void *#, for each line, quad,
 *   etc. Just like the user flags, this field is not available for vertices,
 *   which does no harm since the vertices have a unique and continuous number
 *   unlike the structured objects lines and quads. You can access these user
 *   pointers through the functions #user_pointer()#, #clear_user_pointer()#
 *   and #set_user_pointer(p)# in the accessor classes.
 *
 *   You should not rely on any specific value for the user pointer if you
 *   haven't set it yourself before. In special, the pointers are not inherited
 *   upon refinement. In principle, all user pointers should be #NULL# pointers
 *   at the time a line, quad, etc comes into life. #NULL# pointers are always
 *   a good choice, since they raise an exception when being dereferenced.
 *   The usual warning about the missing type safety of #void# pointers are
 *   obviously in place here; responsibility for correctness of types etc
 *   lies entirely with the user of the pointer.
 *
 *
 *   \subsection{Boundary approximation}
 *
 *   You can specify a boundary function: if a new vertex is created on a
 *   side or face at the boundary, this function is used to compute where
 *   it will be placed. See \Ref{Boundary} for the details. Usage with
 *   the #Triangulation# object is then like this (let #Ball# be a class
 *   derived from #Boundary<2>#):
 *   \begin{verbatim}
 *     void main () {
 *       Triangulation<2> tria;
 *                                        // set the boundary function
 *       Ball ball;
 *       tria.set_boundary (&ball);
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
 *           tria.execute_refinement();
 *         };
 *     };            
 *   \end{verbatim}
 *
 *   You should take note of one caveat: if you have concave boundaries, you
 *   must make sure that a new boundary vertex does not lie to much inside the
 *   to be refined cell. The reason is that the center vertex is placed at the
 *   point which is the arithmetic mean of the eight surrounding vertices.
 *   Therefore if your new boundary vertex is too near the center of the old
 *   quadrilateral or hexahedron, the distance to the midpoint vertex will become
 *   too small, thus generating distorted cells. Remedy: you have to take care
 *   of such situations when defining the coarse grid.
 *
 *
 *   \subsection{Technical details}
 *
 *   \subsubsection{Algorithms for mesh regularisation and smoothing upon refinement}
 *
 *   We chose an inductive point of view: since upon creation of the
 *   triangulation all cells are on the same level, all regularity assumptions
 *   regarding the maximum difference in level of cells sharing a common face,
 *   edge or vertex hold. Since we use the regularisation and smoothing in
 *   each step of the mesh history, when coming to the point of refining it
 *   further the assumptions also hold.
 *
 *   The regularisation and smoothing is done in the #prepare_refinement#
 *   function, which is called by #execute_refinement# at the very beginning.
 *   It decides which additional cells to flag for refinement by looking at the
 *   old grid and the refinement flags for each cell.
 *
 *   \begin{itemize}
 *   \item {\it Regularisation:} The algorithm walks over all cells checking
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
 *   \item {\it Smoothing:}
 *     \begin{itemize}
 *     \item #limit_level_difference_at_vertices#:
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
 *     \item #eliminate_unrefined_islands#:
 *       For each cell we count the number of neighbors which are refined or
 *       flagged for refinement. If this exceeds the total number of neighbors
 *       (which is the number of faces minus the number of faces of this cell
 *       which are located on the boundary), then this cell is flagged for
 *       refinement. Since this may lead to cells on the same level which also
 *       will need refinement, we will need additional loops of regularisation
 *       and smoothing over all cells until nothing changes any more.
 *
 *     \item #eliminate_refined_*_islands#:
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
 *       of this class's description.
 *
 *       The same applies as above: several loops may be necessary.
 *     \end{itemize}
 *   \end{itemize}
 *
 *   Regularisation and smoothing are a bit complementary in that we check
 *   whether we need to set additional refinement flags when being on a cell
 *   flagged for refinement (regularisation) or on a cell not flagged for
 *   refinement. This makes readable programming easier.
 *
 *   All the described algorithms apply only for more than one space dimension,
 *   since for one dimension no restrictions apply. It may be necessary to
 *   apply some smoothing for multigrid algorithms, but this has to be decided
 *   upon later.
 *
 *
 *   \subsubsection{Implementational conventions for two spatial dimensions}
 *
 *   There is a convention about the direction of the bounding lines of quads in
 *   2D. The direction of a line is the direction of point 0 towards point 1. We
 *   define, that allowed cells contain of lines of which the direction is
 *   as follows:
 *   \begin{verbatim}
 *          2
 *      3--->---2
 *      |       |
 *     3^       ^1
 *      |       |
 *      0--->---1
 *          0
 *   \end{verbatim}
 *   The number of the vertices and lines is also indicated. This orientation of
 *   lines has to be checked/generated upon construction of a grid and is
 *   preserved upon refinement.
 *
 *   Further we define, that child lines have the same direction as their parent,
 *   i.e. that #subline(0).vertex(0)==line.vertex(0)# and
 *   #subline(1).vertex(1)==line.vertex(1)#. This also implies, that the
 *   first subline (#subline(0)#) is the one at vertex(0) of the old line.
 *
 *   Similarly we define, that the four children of a quad are adjacent to the
 *   vertex with the same number of the old quad.
 *
 *
 *   \subsection{Warning}
 *
 *   It seems impossible to preserve #const#ness of a triangulation through
 *   iterator usage. Thus, if you declare pointers to a #const# triangulation
 *   object, you should be well aware that you might involuntarily alter the
 *   data stored in the triangulation.
 *
 *   @memo Implementation of a multilevel triangulation of a domain
 *   @see TriaRawIterator
 *   @author Wolfgang Bangerth, 1998
 */
template <int dim>
class Triangulation : public TriaDimensionInfo<dim> {
  public:
				     // insert these definitions for gcc2.8,
				     // since it can't inherit typedefs (I
				     // believe it should, but it can't)
    typedef typename TriaDimensionInfo<dim>::raw_line_iterator raw_line_iterator;
    typedef typename TriaDimensionInfo<dim>::line_iterator line_iterator;
    typedef typename TriaDimensionInfo<dim>::active_line_iterator active_line_iterator;

    typedef typename TriaDimensionInfo<dim>::raw_quad_iterator raw_quad_iterator;
    typedef typename TriaDimensionInfo<dim>::quad_iterator quad_iterator;
    typedef typename TriaDimensionInfo<dim>::active_quad_iterator active_quad_iterator;

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
				      *  use the #copy_triangulation# function,
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
    ~Triangulation ();
    
				     /**
				      *  Assign a boundary object to the
				      *  triangulation, which is used to find
				      *  the point where to place new vertices
				      *  on the boundary. Ownership of the
				      *  object remains with the caller of this
				      *  function, i.e. you have to make sure
				      *  that it is not destroyed before
				      *  #Triangulation<>::execute_refinement()#
				      *  is called the last time.
				      *
				      *  If you copy this triangulation to
				      *  another one using the
				      *  #copy_triangulation# function, you must
				      *  make sure that the lifetime of boundary
				      *  object extends to the lifetime of the
				      *  new triangulation as well.
				      */
    void set_boundary (const Boundary<dim> *boundary_object);

				     /**
				      *  Copy a triangulation. This operation is
				      *  not cheap, so you should be careful
				      *  with using this. We do not implement
				      *  this function as a copy constructor,
				      *  since it makes it easier to maintain
				      *  collections of triangulations if you
				      *  can assign them values later on.
				      *
				      *  Keep in mind that this function also
				      *  copies the pointer to the boundary
				      *  descriptor previously set by the
				      *  #set_boundary# function. You must
				      *  therefore also guarantee that the
				      *  boundary objects has a lifetime at
				      *  least as long as the copied
				      *  triangulation.
				      *
				      *  This triangulation must be empty
				      *  beforehand.
				      */
    void copy_triangulation (const Triangulation<dim> &old_tria);
    
				     /**
				      * Create a triangulation from a list
				      * of vertices and a list of cells, each of
				      * the latter being a list of #1<<dim#
				      * vertex indices. The triangulation must
				      * be empty upon calling this function and
				      * the cell list should be useful (connected
				      * domain, etc.).
				      *
				      * Material data for the cells is given
				      * within the #cells# array, while boundary
				      * information is given in the
				      * #subcelldata# field.
				      *
				      * The numbering of vertices within the
				      * #cells# array is subject to some
				      * constraints; see the general class
				      * documentation for this.
				      */
    void create_triangulation (const vector<Point<dim> >    &vertices,
			       const vector<CellData<dim> > &cells,
			       const SubCellData            &subcelldata);
    
				     /**
				      * Initialize the triangulation with a
				      * hypercube (line in 1D, square in 2D, etc)
				      * consisting of exactly one cell. The
				      * hypercube volume is the tensor product
				      * of the intervall $[left,right]$ in the
				      * present number of dimensions, where
				      * the limits are given as arguments. They
				      * default to zero and unity, then producing
				      * the unit hypercube.
				      *
				      * The triangulation needs to be void
				      * upon calling this function.
				      */
    void create_hypercube (const double left = 0.,
			   const double right= 1.);

    				     /**
				      * Initialize the triangulation with a
				      * hyper-L consisting of exactly #2^dim-1#
				      * cells. See the general documentation for a
				      * description of the L-region. The limits
				      * default to minus unity and unity.
				      *
				      * The triangulation needs to be void
				      * upon calling this function.
				      */
    void create_hyper_L (const double left = -1.,
			 const double right= 1.);

				     /**
				      * Initialize the triangulation with a
				      * hyperball, i.e. a circle or a ball.
				      * See the general documentation for a
				      * more concise description. The center of
				      * the hyperball default to the origin,
				      * the radius defaults to unity.
				      *
				      * The triangulation needs to be void
				      * upon calling this function.
				      */
    void create_hyper_ball (const Point<dim> &center = Point<dim>(),
			    const double radius = 1.);

				     /**
				      * Distort the grid by randomly moving
				      * around all the vertices of the grid.
				      * The direction of moving is random,
				      * while the length of the shift vector
				      * has a value of #factor# times the
				      * minimal length of the active lines
				      * adjacent to this vertex. Note that
				      * #factor# should obviously be well
				      * below #0.5#.
				      *
				      * If #keep_boundary# is set to #true#
				      * (which is the default), then boundary
				      * vertices are not moved.
				      */
    void distort_random (const double factor,
			 const bool   keep_boundary=true);
    
				      
				     /**
				      *  @name Mesh refinement
				      */
				     /*@{*/
				     /**
				      *  Flag all active cells for refinement.
				      *  This will refine
				      *  all cells of all levels which are not
				      *  already refined (i.e. only cells are
				      *  refined which do not yet have
				      *  children). The cells are only flagged,
				      *  not refined, thus you have the chance
				      *  to save the refinement flags.
				      */
    void set_all_refine_flags ();

				     /**
				      *  Refine all cells #times# times, by
				      *  alternatingly calling #refine_global()#
				      *  and #execute_refinement()#.
				      *  This function actually starts the
				      *  refinement process, so you have no way
				      *  to store the refinement flags.
				      */
    void refine_global (const unsigned int times);

				     /**
				      * Refine the triangulation according to
				      * the given criteria. The criterion is a
				      * #double# value for each cell which
				      * determines which cells are to be refine
				      * by comparison with the threshold: if the
				      * value for a cell is larger than the
				      * threshold, the cell is flagged for
				      * refinement. It is your duty to guarantee
				      * that the threshold value is in a
				      * resonable range.
				      *
				      * The cells are only flagged for
				      * refinement, they are not actually
				      * refined. To do so, you have to call the
				      * #execute_refinement# function.
				      *
				      * There are more sophisticated strategies
				      * for mesh refinement; refer to the
				      * following functions and to the general
				      * doc for this class for more information.
				      */
    void refine (const dVector &criteria,
		 const double   threshold);

				     /**
				      * Analogue to the #refine# function:
				      * flag all cells for coarsening for
				      * which the criterion is less than the
				      * given threshold.
				      */
    void coarsen (const dVector &criteria,
		  const double   threshold);
    
				     /**
				      * Refine the triangulation by refining
				      * a certain fraction #top_fraction_of_cells#
				      * with the highest error. Likewise coarsen
				      * the fraction #bottom_fraction_of_cells#
				      * with the least error. To actually
				      * perform the refinement, call
				      * #execute_refinement#.
				      *
				      * #fraction_of_cells# shall be a value
				      * between zero and one.
				      *
				      * Refer to the general doc of this class
				      * for more information.
				      */
    void refine_and_coarsen_fixed_number (const dVector &criteria,
					  const double   top_fraction_of_cells,
					  const double   bottom_fraction_of_cells);

				     /**
				      * Refine the triangulation by flagging
				      * those cells which make up a certain
				      * #top_fraction# of the total error.
				      * Likewise, coarsen all cells which
				      * make up only #bottom_fraction#.
				      * To actually perform the refinement, call
				      * #execute_coarsening_and_refinement#.
				      *
				      * #*_fraction# shall be a values
				      * between zero and one.
				      *
				      * Refer to the general doc of this class
				      * for more information.
				      */
    void refine_and_coarsen_fixed_fraction (const dVector &criteria,
					    const double   top_fraction,
					    const double   bottom_fraction);
    
				     /**
				      *  Refine all cells on all levels which
				      *  were previously flagged for refinement.
				      *
				      *  The function resets all refinement
				      *  flags to false.
				      *
                                      *  See the general docs for more
                                      *  information.
                                      *
				      *  This function is dimension specific.
				      */ 
    void execute_refinement ();

				     /**
				      *  Prepare the refinement process by
				      *  fixing the
				      *  closure of the refinement in #dim>=2#
				      *  (make sure that no two cells are
				      *  adjacent with a refinement level
				      *  differing with more than one), etc.
				      *  It performs some mesh smoothing if
				      *  the according flag was given to the
				      *  constructor of this class.
				      *  The function returns whether additional
				      *  cells have been flagged for refinement.
				      *  
				      *  See the general
				      *  doc of this class for more information.
				      *
				      *  This function is mostly dimension
				      *  independent.
				      */
    bool prepare_refinement ();

				     /**
				      * Coarsen all cells which were flagged for
				      * coarsening, or rather: delete all
				      * children of those cells of which all
				      * child cells are flagged for coarsening
				      * and several other constraints hold (see
				      * the general doc of this class).
				      *
				      * The function resets all coarsen
				      * flags to false. It uses the user flags,
                                      * so make sure to save them if you still
                                      * need them after calling this function.
				      *
                                      * See the general docs for more
                                      * information.
				      */
    void execute_coarsening ();

    				     /**
				      * This function does the analogue
				      * of the #prepare_refinement# function:
				      * flagging and deflagging cells in
				      * preparation of the actual coarsening
				      * step. This includes deleting coarsen
                                      * flags from cells which may not be
                                      * deleted (e.g. because one neighbor is
                                      * more refined than the cell), doing
                                      * some smoothing, etc.
				      *
				      * The effect is that only those cells
				      * are flagged for coarsening which
				      * will actually be coarsened. This
				      * includes the fact that all flagged
				      * cells belong to parent cells of which
				      * all children are flagged.
				      *
				      * The function returns whether some
				      * cells' flagging has been changed in
				      * the process.
				      *
				      * This function uses the user flags, so
				      * store them if you still need them
				      * afterwards.
                                      */
    bool prepare_coarsening ();


				     /**
				      * Execute both refinement and coarsening
				      * of the triangulation.
				      */
    void execute_coarsening_and_refinement ();
				     /*@}*/

				     /**
				      *  @name History of a triangulation
				      */
    				     /*@{*/
				     /**
				      *  Save the addresses of the cells which
				      *  are flagged for refinement to #out#.
				      *  For usage, read the general
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
				      *  #save_refine_flags#.
				      */
    void load_refine_flags (istream &in);

    				     /**
				      *  Read the information stored by
				      *  #save_refine_flags#.
				      */
    void load_refine_flags (const vector<bool> &v);

				     /**
				      * Analogue to #save_refine_flags#.
				      */
    void save_coarsen_flags (ostream &out) const;

				     /**
				      * Same as above, but store the flags to
				      * a bitvector rather than to a file.
				      */
    void save_coarsen_flags (vector<bool> &v) const;

    				     /**
				      * Analogue to #load_refine_flags#.
				      */
    void load_coarsen_flags (istream &out);

				     /**
				      * Analogue to #load_refine_flags#.
				      */
    void load_coarsen_flags (const vector<bool> &v);

    				     /*@}*/


				     /**
				      *  @name User flag handling
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
				      *  #save_refine_flags# for more
				      *  details.
				      */
    void save_user_flags (ostream &out) const;

				     /**
				      * Same as above, but store the flags to
				      * a bitvector rather than to a file.
				      */
    void save_user_flags (vector<bool> &v) const;

				     /**
				      *  Read the information stored by
				      *  #save_user_flags#.
				      */
    void load_user_flags (istream &in);

    				     /**
				      *  Read the information stored by
				      *  #save_user_flags#.
				      */
    void load_user_flags (const vector<bool> &v);

				     /**
				      * Save the user flags on lines.
				      */
    void save_user_flags_line (ostream &out) const;

				     /**
				      * Same as above, but store the flags to
				      * a bitvector rather than to a file.
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
				     /*@}*/

				     /* ------------------------------------ */
    
				     /**
				      *  @name Cell iterator functions
				      */
				     /*@{*/
				     /**
				      *  Return iterator to the first cell, used
				      *  or not, on level #level#. If a level
				      *  has no cells, a past-the-end iterator
				      *  is returned.
				      *
				      *  This function calls #begin_raw_line#
				      *  in 1D and #begin_raw_quad# in 2D.
				      */
    raw_cell_iterator    begin_raw   (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first used cell
				      *  on level #level#.
				      *
				      *  This function calls #begin_line#
				      *  in 1D and #begin_quad# in 2D.
				      */
    cell_iterator        begin       (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first active
				      *  cell on level #level#.
				      *
				      *  This function calls #begin_active_line#
				      *  in 1D and #begin_active_quad# in 2D.
				      */
    active_cell_iterator begin_active(const unsigned int level = 0) const;

				     /**
				      *  Return iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      *
				      *  This function calls #end_line#
				      *  in 1D and #end_quad# in 2D.
				      */
    raw_cell_iterator    end () const;

				     /**
				      * Return an iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    cell_iterator        end (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    raw_cell_iterator    end_raw (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If #level#
				      * is the last level, then this returns
				      * #end()#.
				      */
    active_cell_iterator end_active (const unsigned int level) const;

    
				     /**
				      *  Return an iterator pointing to the
				      *  last cell, used or not.
				      *
				      *  This function calls #last_raw_line#
				      *  in 1D and #last_raw_quad# in 2D.
				      */
    raw_cell_iterator    last_raw () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  cell of the level #level#, used or not.
				      *
				      *  This function calls #last_raw_line#
				      *  in 1D and #last_raw_quad# in 2D.
				      */
    raw_cell_iterator    last_raw (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used cell.
				      *
				      *  This function calls #last_line#
				      *  in 1D and #last_quad# in 2D.
				      */
    cell_iterator        last () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used cell on level #level#.
				      *
				      *  This function calls #last_line#
				      *  in 1D and #last_quad# in 2D.
				      */
    cell_iterator        last (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active cell.
				      *
				      *  This function calls #last_active_line#
				      *  in 1D and #last_active_quad# in 2D.
				      */
    active_cell_iterator last_active () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active cell on level #level#.
				      *
				      *  This function calls #last_active_line#
				      *  in 1D and #last_active_quad# in 2D.
				      */
    active_cell_iterator last_active (const unsigned int level) const;
				     //@}

    				     /*---------------------------------------*/
    				     /*---------------------------------------*/

    				     /**
				      *  @name Face iterator functions
				      */
				     /*@{*/
				     /**
				      *  Return iterator to the first face, used
				      *  or not, on level #level#. If a level
				      *  has no faces, a past-the-end iterator
				      *  is returned.
				      *
				      *  This function calls #begin_raw_line#
				      *  in 2D and #begin_raw_quad# in 3D.
				      */
    raw_face_iterator    begin_raw_face   (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first used face
				      *  on level #level#.
				      *
				      *  This function calls #begin_line#
				      *  in 2D and #begin_quad# in 3D.
				      */
    face_iterator        begin_face       (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first active
				      *  face on level #level#.
				      *
				      *  This function calls #begin_active_line#
				      *  in 2D and #begin_active_quad# in 3D.
				      */
    active_face_iterator begin_active_face(const unsigned int level = 0) const;

				     /**
				      *  Return iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      *
				      *  This function calls #end_line#
				      *  in 2D and #end_quad# in 3D.
				      */
    raw_face_iterator    end_face () const;

				     /**
				      * Return an iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    face_iterator        end_face (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    raw_face_iterator    end_raw_face (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If #level#
				      * is the last level, then this returns
				      * #end()#.
				      */
    active_face_iterator end_active_face (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the
				      *  last face, used or not.
				      *
				      *  This function calls #last_raw_line#
				      *  in 2D and #last_raw_quad# in 3D.
				      */
    raw_face_iterator    last_raw_face () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  face of the level #level#, used or not.
				      *
				      *  This function calls #last_raw_line#
				      *  in 2D and #last_raw_quad# in 3D.
				      */
    raw_face_iterator    last_raw_face (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used face.
				      *
				      *  This function calls #last_line#
				      *  in 2D and #last_quad# in 3D.
				      */
    face_iterator        last_face () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used face on level #level#.
				      *
				      *  This function calls #last_line#
				      *  in 2D and #last_quad# in 3D.
				      */
    face_iterator        last_face (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active face.
				      *
				      *  This function calls #last_active_line#
				      *  in 2D and #last_active_quad# in 3D.
				      */
    active_face_iterator last_active_face () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active face on level #level#.
				      *
				      *  This function calls #last_active_line#
				      *  in 2D and #last_active_quad# in 3D.
				      */
    active_face_iterator last_active_face (const unsigned int level) const;
				     //@}

    
				     /*---------------------------------------*/

				     /**
				      *  @name Line iterator functions
				      */
				     /*@{*/
				     /**
				      *  Return iterator to the first line, used
				      *  or not, on level #level#. If a level
				      *  has no lines, a past-the-end iterator
				      *  is returned.
				      */
    raw_line_iterator
    begin_raw_line   (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first used line
				      *  on level #level#.
				      */
    line_iterator
    begin_line       (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first active
				      *  line on level #level#.
				      */
    active_line_iterator
    begin_active_line(const unsigned int level = 0) const;

				     /**
				      *  Return iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      */
    raw_line_iterator
    end_line () const;

				     /**
				      * Return an iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    line_iterator        end_line (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    raw_line_iterator    end_raw_line (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If #level#
				      * is the last level, then this returns
				      * #end()#.
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
				      *  line of the level #level#, used or not.

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
				      *  used line on level #level#.
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
				      *  active line on level #level#.
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
				      *  Return iterator to the first quad, used
				      *  or not, on level #level#. If a level
				      *  has no quads, a past-the-end iterator
				      *  is returned.
				      */
    raw_quad_iterator
    begin_raw_quad   (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first used quad
				      *  on level #level#.
				      */
    quad_iterator
    begin_quad       (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first active
				      *  quad on level #level#.
				      */
    active_quad_iterator
    begin_active_quad(const unsigned int level = 0) const;

				     /**
				      *  Return iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      */
    raw_quad_iterator
    end_quad () const;

				     /**
				      * Return an iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    quad_iterator        end_quad (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    raw_quad_iterator    end_raw_quad (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If #level#
				      * is the last level, then this returns
				      * #end()#.
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
				      *  quad of the level #level#, used or not.

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
				      *  used quad on level #level#.
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
				      *  active quad on level #level#.
				      */
    active_quad_iterator
    last_active_quad (const unsigned int level) const;
				     /*@}*/

				     /*---------------------------------------*/

				     /**
				      *  @name Input/Output functions
				      */
				     /*@{*/
				     /**
				      *  Print the triangulation in GNUPLOT
				      *  format to #out#.
				      */
    void print_gnuplot (ostream &) const;

				     /**
				      *  Print level #level# in GNUPLOT
				      *  format to #out#.
				      */
    void print_gnuplot (ostream &, const unsigned int level) const;

				     /**
				      *  Print cell #cell# in GNUPLOT format
				      *  to #out#.
				      */
    void print_gnuplot (ostream &, const active_cell_iterator &cell) const;
    				     /*@}
				      */

    
				     /**
				      * @name Information about the triangulation
				      */
				     /*@{*/
				     /**
				      *  Return total number of used lines,
				      *  active or not.
				      *
				      * Regarding the computational effort of
				      * this function, the same applies as
				      * for the #n_cells()# function.
				      */
    unsigned int n_lines () const;
    
				     /**
				      *  Return total number of used lines,
				      *  active or not on level #level#.
				      *
				      * Regarding the computational effort of
				      * this function, the same applies as
				      * for the #n_cells()# function.
				      */
    unsigned int n_lines (const unsigned int level) const;
    
				     /**
				      *  Return total number of active lines,
				      *  active or not.
				      *
				      * Regarding the computational effort of
				      * this function, the same applies as
				      * for the #n_cells()# function.
				      */
    unsigned int n_active_lines () const;
    
				     /**
				      *  Return total number of active lines,
				      *  active or not on level #level#.
				      *
				      * Regarding the computational effort of
				      * this function, the same applies as
				      * for the #n_cells()# function.
				      */
    unsigned int n_active_lines (const unsigned int level) const;
    
				     /**
				      *  Return total number of used quads,
				      *  active or not.
				      *
				      * Regarding the computational effort of
				      * this function, the same applies as
				      * for the #n_cells()# function.
				      */
    unsigned int n_quads () const;
    
				     /**
				      *  Return total number of used quads,
				      *  active or not on level #level#.
				      *
				      * Regarding the computational effort of
				      * this function, the same applies as
				      * for the #n_cells()# function.
				      */
    unsigned int n_quads (const unsigned int level) const;
    
				     /**
				      *  Return total number of active quads,
				      *  active or not.
				      *
				      * Regarding the computational effort of
				      * this function, the same applies as
				      * for the #n_cells()# function.
				      */
    unsigned int n_active_quads () const;
    
				     /**
				      *  Return total number of active quads,
				      *  active or not on level #level#.
				      *
				      * Regarding the computational effort of
				      * this function, the same applies as
				      * for the #n_cells()# function.
				      */
    unsigned int n_active_quads (const unsigned int level) const;
    
				     /**
				      *  Return total number of used cells,
				      *  active or not.
				      *  Maps to #n_lines()# in one space
				      *  dimension and so on.
				      *
				      * Be a bit careful with this function,
				      * since it needs quite some computational
				      * ressources: the number of cells is
				      * computed by looping over all cells,
				      * rather than by returning a pre-stored
				      * value. It is therefore often better
				      * to rewrite lines like
				      * #for (i=0; i<tria.n_cells() const; ++i) ...#
				      * by
				      * #const unsigned int n=tria.n_cells();#
				      * #for (i=0; i<n; ++i) ...#.
				      */
    unsigned int n_cells () const;

    				     /**
				      *  Return total number of used cells,
				      *  active or not, on level #level#.
				      *  Maps to #n_lines(level)# in one space
				      *  dimension and so on.
				      *
				      * Regarding the computational effort of
				      * this function, the same applies as
				      * for the #n_cells()# function.
				      */
    unsigned int n_cells (const unsigned int level) const;

    				     /**
				      *  Return total number of active cells.
				      *  Maps to #n_active_lines()# in one space
				      *  dimension and so on.
				      *
				      * Regarding the computational effort of
				      * this function, the same applies as
				      * for the #n_cells()# function.
				      */
    unsigned int n_active_cells () const;

    				     /**
				      * Return total number of active cells
				      * on level #level#.
				      * Maps to #n_active_lines(level)# in one
				      * space dimension and so on.
				      *
				      * Regarding the computational effort of
				      * this function, the same applies as
				      * for the #n_cells()# function.
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
				     //@{
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
    DeclException2 (ExcInvalidVectorSize,
		    int, int,
		    << "The given vector has " << arg1
		    << " elements, but " << arg2 << " were expected.");
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidParameterValue);
				     //@}
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
				      *  format: for each #true#, a #1#
				      *  bit is stored, a #0# bit otherwise.
				      *  The bits are stored as #unsigned char#,
				      *  thus avoiding endianess. They are
				      *  written to #out# in plain text, thus
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
				      * written by #write_bool_vector# and
				      * compare with the magic numbers.
				      */
    static void read_bool_vector (const unsigned int  magic_number1,
				  vector<bool>       &v,
				  const unsigned int  magic_number2,
				  istream            &in);
    
				     /**
				      * Actually delete a cell, which is the
				      * main step for the coarsening process.
				      * This is the dimension dependent part
				      * of #execute_coarsening#.
				      */
    void delete_cell (cell_iterator &cell);
    
				     /**
				      *  Array of pointers pointing to the
				      *  #TriangulationLevel<dim># objects
				      *  storing the data on the different
				      *  levels.
				      *
				      *  Usage is like #levels[3]->quads#.
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
				      *  Pointer to a boundary object.
				      */
    const Boundary<dim>             *boundary;

				     /**
				      *  Do some smoothing in the process
				      *  of refining the triangulation. See
				      *  the general doc of this class for
				      *  more information about this.
				      */
    MeshSmoothing                    smooth_grid;

				     // Friendship includes local classes.
    friend class TriaAccessor<dim>;
    friend class LineAccessor<dim>;
    friend class QuadAccessor<dim>;

    friend class CellAccessor<1>;
    friend class CellAccessor<2>;
    
    friend class TriaRawIterator<1,LineAccessor<1> >;
    friend class TriaRawIterator<1,CellAccessor<1> >;
    friend class TriaRawIterator<2,LineAccessor<2> >;
    friend class TriaRawIterator<2,QuadAccessor<2> >;
    friend class TriaRawIterator<2,CellAccessor<2> >;

    friend class DoFHandler<dim>;
    friend class MGDoFHandler<dim>;
};






/*----------------------------   tria.h     ---------------------------*/
/* end of #ifndef __tria_H */
#endif
/*----------------------------   tria.h     ---------------------------*/





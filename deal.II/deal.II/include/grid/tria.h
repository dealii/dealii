/*----------------------------   tria.h     ---------------------------*/
/*      $Id$                 */
#ifndef __tria_H
#define __tria_H
/*----------------------------   tria.h     ---------------------------*/

#include <vector>
#include <grid/tria_line.h>
#include <grid/tria_quad.h>
#include <grid/point.h>



//forward declaration needed
template <int dim> class Boundary;

template <int dim> class TriaAccessor;
template <int dim> class LineAccessor;
template <int dim> class QuadAccessor;
template <int dim> class CellAccessor;


template <int dim, class Accessor> class TriaRawIterator;
template <int dim, class Accessor> class TriaIterator;
template <int dim, class Accessor> class TriaActiveIterator;

template <int dim> class DoFHandler;

template <int dim> struct CellData;
struct SubCellData;

class dVector;


class istream;
class ostream;



/**
 *  General template for information belonging to one level of a multilevel
 *  hierarchy of a triangulation. This template is only declared to allow
 *  specializations for different dimensions.
 *
 *  @see TriangulationLevel<1>
 *  @see TriangulationLevel<2>
 */
template <int dim>
class TriangulationLevel;


/**
 *  Store all information which belongs to one level of the multilevel hierarchy.
 *
 *  In #TriangulationLevel<0># all data is stored which is not
 *  dependant on the dimension, e.g. a field to store the
 *  refinement flag for the cells (what a cell actually is
 *  is declared elsewhere), etc. Actually, it is only cell-based
 *  data, like neighborship info or refinement flags. There is another
 *  field, which may fit in here, namely the material data (for cells)
 *  or the boundary indicators (for faces), but since we need for a line
 *  or quad either boundary information or material data, we store them
 *  with the lines and quads rather than with the common data. We may,
 *  however lose some memory in three dimensions, when we need the
 *  material data for cell, boundary data for the quads, but nothing
 *  for the lines. Since we only store one byte per line, quad or hex,
 *  this is a minor loss and we can live with that.
 *
 *  @memo Information belonging to one level of the multilevel hierarchy.
 */
class TriangulationLevel<0> {
  public:
				     /**
				      *  Flags for the cells whether they are
				      *  to be refined or not. The meaning
				      *  what a cell is, is dimension specific,
				      *  therefore also the length of this
				      *  vector depends on the dimension: in
				      *  one dimension, the length of this
				      *  vector equals the length of the
				      *  #lines# vector, in two dimensions
				      *  that of the #quads# vector, etc.
				      */
    vector<bool> refine_flags;
    
				     /**
				      *  Levels and indices of the neighbors
				      *  of the cells. Convention is, that the
				      *  neighbors of the cell with index #i#
				      *  are stored in the fields following
				      *  $i*(2*real_space_dimension)$, e.g. in
				      *  one spatial dimension, the neighbors
				      *  of cell 0 are stored in #neighbors[0]#
				      *  and #neighbors[1]#, the neighbors of
				      *  cell 1 are stored in #neighbors[2]#
				      *  and #neighbors[3]#, and so on.
				      *
				      *  In neighbors, #neighbors[i].first# is
				      *  the level, while #neighbors[i].first#
				      *  is the index of the neighbor.
				      *
				      *  If a neighbor does not exist (cell is
				      *  at the boundary), #level=index=-1#
				      *  is set.
				      *
				      *  {\bf Conventions:} The #i#th neighbor
				      *  of a cell is the one which shares
				      *  the #i#th face (#Line# in 2D, #Quad#
				      *  in 3D) of this cell.
				      *
				      *  The neighbor of a cell has at most the
				      *  same level as this cell, i.e. it may
				      *  or may not be refined.
				      *
				      *  In one dimension, a neighbor may have
				      *  any level less or equal the level of
				      *  this cell. If it has the same level,
				      *  it may be refined an arbitrary number
				      *  of times, but the neighbor pointer
				      *  still points to the cell on the same
				      *  level, while the neighbors of the
				      *  childs of the neighbor may point to
				      *  this cell or its children.
				      *
				      *  In two and more dimensions, the
				      *  neighbor is either on the same level
				      *  and refined (in which case its children
				      *  have neighbor pointers to this cell or
				      *  its direct children), unrefined on
				      *  the same level or one level down (in
				      *  which case its neighbor pointer points
				      *  to the mother cell of this cell).
				      */
    vector<pair<int,int> > neighbors;
    
				     /**
				      *  Reserve enough space to accomodate
				      *  #total_cells# cells on this level.
				      *  Since there are no #used# flags on this
				      *  level, you have to give to total number
				      *  of cells, not only the number of newly
				      *  to accomodate ones, like in the
				      *  #TriangulationLevel<N>::reserve_space#
				      *  functions, with #N>0#.
				      *
				      *  Since the
				      *  number of neighbors per cell depends
				      *  on the dimensions, you have to pass
				      *  that additionally.
				      */
    void reserve_space (const unsigned int total_cells,
			const unsigned int dimension);

				     /**
				      *  Check the memory consistency of the
				      *  different containers. Should only be
				      *  called with the prepro flag #DEBUG#
				      *  set. The function should be called from
				      *  the functions of the higher
				      *  #TriangulationLevel# classes.
				      */
    void monitor_memory (const unsigned int true_dimension) const;

				     /**
				      *  Exception
				      */
    DeclException3 (ExcMemoryWasted,
		    char*, int, int,
		    << "The container " << arg1 << " contains "
		    << arg2 << " elements, but it`s capacity is "
		    << arg3 << ".");
				     /**
				      *  Exception
				      */
    DeclException2 (ExcMemoryInexact,
		    int, int,
		    << "The containers have sizes " << arg1 << " and "
		    << arg2 << ", which is not as expected.");
				     /**
				      *  Exception
				      */
    DeclException0 (ExcUnusedMemoryAtEnd);
};



/**
 *  Store all information which belongs to one level of the multilevel hierarchy.
 *  
 *  In one dimension, this is a list of the lines associated with this level,
 *  as well as a list with the indices of the children of these lines.
 *  The #TriangulationsLevel# objects of higher dimensions are derived from
 *  this one.
 *
 *  @memo Information belonging to one level of the multilevel hierarchy.
 */
class TriangulationLevel<1> : public TriangulationLevel<0> {
  private:

				     /**
				      *  This subclass groups together all that
				      *  is needed to describe the lines on one
				      *  level.
				      */
    struct LinesData {
					 /**
					  *  Vector of the lines belonging to
					  *  this level. The index of the line
					  *  on this level equals the index in
					  *  this container, while the global
					  *  index of a line is stored in the
					  *  line itself.
					  */
	vector<Line> lines;
					 /**
					  *  Index of the first child of a line
					  *  in the list on the next level.
					  *  Since when lines are refined, both
					  *  children are created at the same
					  *  time, they are appended to the list
					  *  on the next level after each other.
					  *  We therefore only store the index
					  *  of the first child, the second
					  *  follows immediately afterwards.
					  *
					  *  If a line has no children, -1 is
					  *  stored in this list. A line is
					  *  called active if it has no
					  *  children. The function
					  *  #TriaIterator::active()#
					  *  tests for this.
					  */
	vector<int>  children;
	
					 /**
					  *  Vector storing whether a line is
					  *  used in the #lines# vector.
					  *
					  *  Since it is difficult to delete
					  *  elements in a #vector#, when an
					  *  element is not needed any more
					  *  (e.g. after derefinement), it is
					  *  not deleted from the list, but
					  *  rather the according #used# flag
					  *  is set to #false#.
					  */
	vector<bool> used;

					 /**
					  *  Make available a field for user data,
					  *  one bit per line. This field is usually
					  *  used when an operation runs over all
					  *  cells and needs information whether
					  *  another cell (e.g. a neighbor) has
					  *  already been processed.
					  *
					  *  You can clear all used flags using
					  *  #Triangulation<>::clear_user_flags()#.
					  */
	vector<bool> user_flags;

					 /**
					  * Store boundary and material data. In
					  * one dimension, this field stores
					  * the material id of a line, which is
					  * number between 0 and 254. In more
					  * than one dimension, lines have no
					  * material id, but they may be at the
					  * boundary; then, we store the
					  * boundary indicator in this field,
					  * which denotes to which part of the
					  * boundary this line belongs and which
					  * boundary conditions hold on this
					  * part. The boundary indicator also
					  * is a number between zero and 254;
					  * the id 255 is reserved for lines
					  * in the interior and may be used
					  * to check whether a line is at the
					  * boundary or not, which otherwise
					  * is not possible if you don't know
					  * which cell it belongs to.
					  */
	vector<unsigned char> material_id;
    };
    
  public:
    				     /**
				      *  Data about the lines.
				      */
    LinesData lines;

    				     /**
				      *  Assert that enough space is allocated
				      *  to accomodate #new_lines# new lines.
				      *  This function does not only call
				      *  #vector::reserve()#, but does really
				      *  append the needed elements.
				      *  There are pendants for higher
				      *  dimensions, which you have to call
				      *  explicitely (they can't hand down the
				      *  call because there is no easy relation
				      *  between the number of new quads and
				      *  the number of new lines, etc.). Also
				      *  don't forget to call the
				      *  #TriangulationLevel<0>::reserve_space#
				      *  function.
				      */
    void reserve_space (const unsigned int new_lines);

				     /**
				      *  Check the memory consistency of the
				      *  different containers. Should only be
				      *  called with the prepro flag #DEBUG#
				      *  set. The function should be called from
				      *  the functions of the higher
				      *  #TriangulationLevel# classes.
				      */
    void monitor_memory (const unsigned int true_dimension) const;
};




/**
 *  Store all information which belongs to one level of the multilevel hierarchy.
 *
 *  In 2D this is a vector of the lines and one of the
 *  quads on this levels, as well as a the two associated vectors holding
 *  information about the children of these lines and quads.
 *
 *  The vector of lines and their children is derived from
 *  #TriangulationLevel<1>#.
 *
 *  @memo Information belonging to one level of the multilevel hierarchy.
 */
class TriangulationLevel<2> :  public TriangulationLevel<1>
{
				     /**
				      *  This subclass groups together all that
				      *  is needed to describe the quads on one
				      *  level.
				      *
				      *  It is fully analogous to the
				      *  #LinesData# structure inherited from
				      *  #Triangulation<1>#.
				      */
    struct QuadsData {
					 /**
					  *  Same as for the #lines# array.
					  */
	vector<Quad> quads;
					 /**
					  *  Same as for the #LineData::chilren#
					  *  array, but since there are four
					  *  children, the index points to the
					  *  first while the other three are
					  *  following immediately afterwards.
					  */
	vector<int>  children;

					 /**
					  *  Same as for #LineData::used#.
					  */
	vector<bool> used;

					 /**
					  *  Same as for #LineData::used#.
					  */
	vector<bool> user_flags;

					 /**
					  * Store boundary and material data. In
					  * two dimension, this field stores
					  * the material id of a quad, which is
					  * number between 0 and 254. In more
					  * than two dimensions, quads have no
					  * material id, but they may be at the
					  * boundary; then, we store the
					  * boundary indicator in this field,
					  * which denotes to which part of the
					  * boundary this line belongs and which
					  * boundary conditions hold on this
					  * part. The boundary indicator also
					  * is a number between zero and 254;
					  * the id 255 is reserved for quads
					  * in the interior and may be used
					  * to check whether a quad is at the
					  * boundary or not, which otherwise
					  * is not possible if you don't know
					  * which cell it belongs to.
					  */
	vector<unsigned char> material_id;
    };
    
  public:
				     /**
				      *  Data about the quads.
				      */
    QuadsData quads;

    				     /**
				      *  Assert that enough space is allocated
				      *  to accomodate #new_quads# new quads.
				      */
    void reserve_space (const unsigned int new_quads);

				     /**
				      *  Check the memory consistency of the
				      *  different containers. Should only be
				      *  called with the prepro flag #DEBUG#
				      *  set. The function should be called from
				      *  the functions of the higher
				      *  #TriangulationLevel# classes.
				      */
    void monitor_memory (const unsigned int true_dimension) const;
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
 *                      endc = (level == levels.size()-1 ?
 *                              cell_iterator(end()) :
 *                              begin (level+1));
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
 *
 *   \subsection{Refinement of a triangulation}
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
 *   However, some degradation of approximation properties has been observed
 *   for grids which were refined more than once across a face, as descibed
 *   above, so there is also a practical justification for the above.
 *   It can also be shown, that such degradation occurs if the
 *   triangulation contains vertices which are member of cells with levels
 *   differing by more than one. One such example is the following:
 *   \begin{verbatim}
 *     |     |     |     |
 *     x-----x-----x--x--x--
 *     |     |     |  |  |
 *     |     |     x--x--x
 *     |     |     |  |  |
 *     x-----x-----x--x--x--
 *     |           |     |
 *     |           |     |
 *     |           |     |
 *     |           x-----x--
 *     |           |     |
 *     |           |     |
 *     |           |     |
 *     x-----------x-----x--
 *   \end{verbatim}
 *   It seems that in two space dimensions, the maximum jump in levels between
 *   cells sharing a common vertex is two (as in the example above). This is
 *   not true if more than four cells meet at a vertex. It is not uncommon
 *   that a coarse (initial) mesh contains vertices at which six or even eight
 *   cells meet, when small features of the domain have to be resolved even on
 *   the coarsest mesh. In that case, the maximum difference in levels is
 *   three or four, respectively. The problem gets even worse in three space
 *   dimensions.
 *
 *   Looking at an interpolation of the second derivative of the finite
 *   element solution (asuming bilinear finite elements), one sees that the
 *   numerical solution is almost totally wrong, compared with the true second
 *   derivative. Indeed, on regular meshes, there exist sharp estimations that
 *   the $H^2$-error is only $O(1)$, so we should not be suprised; however, the
 *   numerical solution may show a value for the second derivative which may
 *   be a factor of ten away from the true value. These problems are located
 *   on the small cell adjacent to the center vertex, where cells of
 *   non-subsequent levels meet, as well as on the upper and right neighbor
 *   of this cell (but with a less degree of deviation from the true value).
 *
 *   Due to the approximational problems described above, the
 *   #Triangulation# constructor takes an argument specifying whether a
 *   smoothing step shall be performed on the grid each time #execute_refinement#
 *   is called. The default is that such a step not be done, since this results
 *   in additional cells being produced, which may not be necessary in all
 *   cases. If switched on, calling #execute_refinement# results in
 *   flagging additional cells for refinement to avoid
 *   vertices as the ones mentioned. The algorithms for both regularisation
 *   and smoothing of triangulations are described below in the section on
 *   technical issues. The reason why this parameter must be given to the
 *   constructor rather than to #execute_refinement# is that it would result
 *   in algorithmic problems if you called #execute_refinement# once without
 *   and once with smoothing, since then in some refinement steps would need
 *   to be refined twice.
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
 *   argument.
 *
 *   There are two variations of this function, which rely on #refine# by
 *   computing the threshold from other information:
 *   \begin{itemize}
 *   \item #refine_fixed_number#: this function takes a vector as above and
 *     a value between zero and one denoting the fraction of cells to be
 *     refined. For this purpose, it sorts the criteria per cell and takes
 *     the threshold to be the one belonging to the cell with the
 *     #fraction times n_active_cells# highest criterion. For example, if
 *     the fraction is $0.3$, the threshold is computed to a value such that
 *     30 per cent of cells have a criterion higher than the threshold and are
 *     thus flagged for refinement. The flagging for refinement is done through
 *     the central #refine# function.
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
 *     However, for singular functions or error functionals, you may want to
 *     chose a smaller value to avoid overrefinement in regions which do not
 *     contribute much to the error.
 *
 *   \item #refine_fixed_fraction#: this function computes the threshold such
 *     that the number of cells getting flagged for refinement makes up for a
 *     certain fraction of the total error. If this fraction is 50 per cent,
 *     for example, the threshold is computed such that the cells with a
 *     criterion greater than the threshold together account for half of the
 *     total error.
 *
 *     This strategy is more suited for singular functions and error
 *     functionals, but may lead to very slow convergence of the grid
 *     if only few cells are refined in each step.
 *
 *     From the implementational point, this time we really need to
 *     sort the array of criteria. However, it is not necessary to sort
 *     the whole array, since for example if you chose the fraction at
 *     50 per cent of the total error, it is only necessary to sort at
 *     most the 50 per cent of cells ranking topmost in the list of error
 *     per cell. It is thus reasonable to use an algorithm like
 *     #partial_sort# of the C++ standard library, which only sorts part
 *     of the array and lets the rest unsorted. However, in many cases
 *     much fewer than 50 per cent of the cells account for 50 per cent
 *     of the error, so it may be possible to get away with sorting less
 *     than 50 per cent of the cells. We therefore divide the whole lot
 *     of 50 per cent of cells into, say, 5 parts, first sort for the
 *     10 per cent with highest error; look whether they together make up
 *     for 50 per cent and if so thats ok, we can leave the rest unsorted;
 *     if not, sort the next 10 per cent, and so on. The default is to
 *     devide the maximum number of cells which may get refined (which
 *     equals the fraction of the total error, as explained above) into
 *     five parts, but this value may be given as a parameter to the
 *     #refine_fixed_fraction# function. For highly singular error
 *     functionals, it may be more efficient to chose a greater number
 *     than five. Chosing a value which is too large should not lead to
 *     a large performance drawback; chosing too small a value however
 *     may lead to significantly higher computational costs for sorting
 *     than necessary.
 *
 *     Just like the other strategy described above, this function only
 *     computes the threshold value and then passes over to #refine#.
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
 *   You may write other information to the output file between different sets
 *   of refinement information, as long as you read it upon re-creation of the
 *   grid.
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
 *
 *   It is convention to clear the user flags using the
 *   #Triangulation<>::clear_user_flags()# function before usage, since it is
 *   often necessary to use the flags in more than one function consecutively and
 *   is then error prone to dedicate one of these to clear the flags.
 *
 *   It is recommended that a functions using the flags states so in its
 *   documentation.
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
 *   whether the present cell is flagged for refinement and a neighbor of the
 *   present cell is refined once less than the present one. If so, flag the
 *   neighbor for refinement. Because of the induction above, there may be no
 *   neighbor with level two less than the present one.
 *
 *   The neighbor thus flagged for refinement may induce more cells which need
 *   to be refined. However, such cells which need additional refinement always
 *   are on one level lower than the present one, so we can get away with only
 *   one sweep over all cells if we do the loop in the reverse way, starting
 *   with those on the highest level. This way, we may flag additional cells
 *   on lower levels, but if these induce more refinement needed, this is
 *   performed later on when we visit them in out backward running loop.
 *
 *   \item {\it Smoothing:} First a list is set up which stores for each vertex
 *   the highest level one of the adjacent cells belongs to. Now, since we did
 *   smoothing in the previous refinement steps also, each cell may only have
 *   vertices with levels at most one greater than the level of the present
 *   cell.
 *
 *   However, if we store the level plus one for cells marked for refinement,
 *   we may end up with cells which have vertices of level two greater than
 *   the cells level. We need to refine this cell also, and need thus also
 *   update the levels of its vertices. This itself may lead to cells needing
 *   refinement, but these are on lower levels, as above, which is why we
 *   may do all kinds of additional flagging in one loop only.
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
    Triangulation (const bool smooth_grid = false);

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
				      */
    void set_boundary (const Boundary<dim> *boundary_object);

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
    void create_hyper_ball (const Point<dim> center = Point<dim>(),
			    const double radius = 1.);
    
				      
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
				      * Refine the triangulation by refining
				      * a certain fraction #fraction_of_cells#
				      * with the highest error. To actually
				      * perform the refinement, call
				      * #execute_refinement#.
				      *
				      * #fraction_of_cells# shall be a value
				      * between zero and one.
				      *
				      * Refer to the general doc of this class
				      * for more information.
				      */
    void refine_fixed_number (const dVector &criteria,
			      const double   fraction_of_cells);

				     /**
				      * Refine the triangulation by flagging
				      * those cells which make up a certain
				      * #fraction_of_error# of the total error.
				      * To actually perform the refinement, call
				      * #execute_refinement#.
				      *
				      * #fraction_of_error# shall be a value
				      * between zero and one.
				      * #n_sorting_parts# shall be one or
				      * greater.
				      *
				      * Refer to the general doc of this class
				      * for more information.
				      */
    void refine_fixed_fraction (const dVector      &criteria,
				const double        fraction_of_error,
				const unsigned int  n_sorting_parts = 5);
    
				     /**
				      *  Refine all cells on all levels which
				      *  were previously flagged for refinement.
				      *
				      *  The function resets all refinement
				      *  flags to false.
				      *
				      *  This function is dimension specific.
				      */ 
    void execute_refinement ();
				     /*@}*/

				     /**
				      *  @name History of a triangulation
				      */
    				     /*@{*/
				     /**
				      *  Save the addresses of the cells which
				      *  are flagged for refinement to #out#.
				      *  The addresses are stored in a binary
				      *  format: for each active cell, a #1#
				      *  bit is stored if the flag is flagged
				      *  for refinement, a #0# bit otherwise.
				      *  The bits are stored as #unsigned char#,
				      *  thus avoiding endianess. They are
				      *  written to #out# in plain text, thus
				      *  amounting to 3.6 bits per active cell
				      *  on the overage. Other information
				      *  (three numbers) is stored as plain text
				      *  as well. The format should therefore be
				      *  interplatform compatible.
				      *
				      *  For usage, read the general
				      *  documentation for this class.
				      */
    void save_refine_flags (ostream &out) const;

				     /**
				      *  Read the information stored by
				      *  #save_refine_flags#.
				      */
    void load_refine_flags (istream &in);
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
				      *  Read the information stored by
				      *  #save_user_flags#.
				      */
    void load_user_flags (istream &in);

				     /**
				      * Save the user flags on lines.
				      */
    void save_user_flags_line (ostream &out) const;

				     /**
				      * Load the user flags located on lines.
				      */
    void load_user_flags_line (istream &in);

    				     /**
				      * Save the user flags on quads.
				      */
    void save_user_flags_quad (ostream &out) const;

				     /**
				      * Load the user flags located on quads.
				      */
    void load_user_flags_quad (istream &in);
				     /*@}*/

    
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
				      *  @name Information about the triangulation
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
				      * #for (i=0; i<tria.n_cells(); ++i) ...#
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
				      *  Return total number of active cells
				      *  on level #level#.
				      *  Maps to #n_active_lines(level)# in one space
				      *  dimension and so on.
				      *
				      * Regarding the computational effort of
				      * this function, the same applies as
				      * for the #n_cells()# function.
				      */
    unsigned int n_active_cells (const unsigned int level) const;

				     /**
				      *  Return number of levels in use.
				      */
    unsigned int n_levels () const;

				     /**
				      * Return the maximum number of cells meeting
				      * at a common vertex. Since this number is
				      * an invariant under refinement, only the cells
				      * on the coarsest level are considered. The
				      * operation is thus reasonably fast. The
				      * invariance is only true for sufficiently
				      * many cells in the coarsest triangulation
				      * (e.g. for a single cell one would be returned),
				      * so a minimum of four is returned in
				      * two dimensions, 8 in three dimensions, etc,
				      * which is how many cells meet if the
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
				      *  Prepare the refinement process by
				      *  allocating enough space, fixing the
				      *  closure of the refinement in #dim>=2#
				      *  (make sure that no two cells are
				      *  adjacent with a refinement level
				      *  differing with more than one), etc.
				      *  It performs some mesh smoothing if
				      *  the according flag was given to the
				      *  constructor of this class.
				      *  See the general
				      *  doc of this class for more information.
				      *
				      *  This function is mostly dimension
				      *  independent.
				      */

    void prepare_refinement ();

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
    bool                             smooth_grid;

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
};






/*----------------------------   tria.h     ---------------------------*/
/* end of #ifndef __tria_H */
#endif
/*----------------------------   tria.h     ---------------------------*/





/*----------------------------   tria.h     ---------------------------*/
/*      $Id$                 */
#ifndef __tria_H
#define __tria_H
/*----------------------------   tria.h     ---------------------------*/

#include <vector.h>
#include <grid/tria_line.h>
#include <grid/tria_quad.h>
#include <grid/point.h>



//#ifdef __GNUC__
//#  include <bvector.h>
//class vector<bool> :
//public bit_vector {};
//#endif


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


class istream;
class ostream;



/**
    General template for information belonging to one level of a multilevel
    hierarchy of a triangulation. This template is only declared to allow
    specializations for different dimensions.

    @see TriangulationLevel<1>
    @see TriangulationLevel<2>
    */
template <int dim>
class TriangulationLevel;


/**
    Store all information which belongs to one level of the multilevel hierarchy.

    In #TriangulationLevel<0># all data is stored which is not
    dependant on the dimension, e.g. a field to store the
    refinement flag for the cells (what a cell actually is
    is declared elsewhere), etc.

    @memo Information belonging to one level of the multilevel hierarchy.
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
    Store all information which belongs to one level of the multilevel hierarchy.
    
    In one dimension, this is a list of the lines associated with this level,
    as well as a list with the indices of the children of these lines.
    The #TriangulationsLevel# objects of higher dimensions are derived from
    this one.

    @memo Information belonging to one level of the multilevel hierarchy.
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
    Store all information which belongs to one level of the multilevel hierarchy.

    In 2D this is a vector of the lines and one of the
    quads on this levels, as well as a the two associated vectors holding
    information about the children of these lines and quads.

    The vector of lines and their children is derived from
    #TriangulationLevel<1>#.

    @memo Information belonging to one level of the multilevel hierarchy.
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
    This class implements some types which differ between the dimensions.
    Declare it to have a template parameter, but do not actually declare
    other types than those explicitely instantiated.
    */
template <int dim>
class TriaDimensionInfo;


/**
    This class implements some types which differ between the dimensions.

    A #line_iterator# is typdef'd to an iterator operating on the
    #lines# member variable of a #Triangulation<1># object. An
    #active_line_iterator# only operates on the active lines.
    #raw_line_iterator# objects operate on all lines, used or not.

    Since we are in one dimension, the following identities are declared:
    \begin{verbatim}
      typedef raw_line_iterator    raw_cell_iterator;
      typedef line_iterator        cell_iterator;
      typedef active_line_iterator active_cell_iterator;
    \end{verbatim}

    To enable the declaration of #begin_quad# and the like in
    #Triangulation<1>#, the #quad_iterator#s are declared as
    #void *#. Thus these types exist, but are useless and will
    certainly make any involuntary use visible.
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
};


/**
    This class implements some types which differ between the dimensions.

    A #line_iterator# is typdef'd to an iterator operating on the
    #lines# member variable of a #Triangulation<1># object. An
    #active_line_iterator# only operates on the active lines.
    #raw_line_iterator# objects operate on all lines, used or not.
    Using #active_line_iterator#s may not be particularly useful since it
    only operates on unrefined lines. However, also refined lines may bound
    unrefined cells if the neighboring cell is refined once more than the
    present one.

    Similarly, #quad_iterator#, #raw_quad_iterator# and
    #active_quad_iterator# are declared.
    
    Since we are in two dimension, the following identities are declared:
    \begin{verbatim}
      typedef raw_quad_iterator    raw_cell_iterator;
      typedef quad_iterator        cell_iterator;
      typedef active_quad_iterator active_cell_iterator;
    \end{verbatim}
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
};





/*------------------------------------------------------------------------*/



/**
    #Triangulation#s denote a hierarchy of levels of elements which together
    form a region in #dim# spatial dimensions.

    This class is written to be as independent of the dimension as possible
    (thus the complex construction of the #TriangulationLevel# classes) to
    allow code-sharing, to allow reducing the need to mirror changes in the code
    for one dimenion to the code for other dimensions. Nonetheless, some of
    the functions are dependent of the dimension and there only exist
    specialized versions for distinct dimensions.


    {\bf Structure and iterators}

    The actual data structure of a #Triangulation# object is rather complex
    and quite inconvenient if one attempted to operate on it directly, since
    data is spread over quite a lot of arrays and other places. However,
    there are ways powerful enough to work on these data structures
    without knowing their exact relations. This is done through the
    concept of iterators (see the STL documentation and \Ref{TriaRawIterator}).
    In order to make things as easy and dimension independent as possible,
    use of class local typedefs is made, see below.
    
    In the base class #TriaDimensionInfo#, a #Cell# is typedef'd to be whatever
    is reasonable for a cell in the respective dimension, i.e. a #Line# in
    one dimension, a #Quad# in two dimensions, and so on.

    The #Triangulation# class provides iterator which enable looping over all
    lines, cells,
    etc without knowing the exact representation used to describe them. Their
    names are typedef's in the #TriaDimensionInfo# base class (thus making them
    local types to this class) and are as follows:

    #raw_line_iterator#: loop over all lines, used or not (declared for
    all dimensions).
    
    #line_iterator#: loop over all used lines (declared for all dimensions).

    #active_line_iterator#: loop over all active lines (declared for all
    dimensions).

    #raw_quad_iterator#: loop over all quads, used or not (declared only
    for #dim>=2#).
    
    #quad_iterator#: loop over all quads (declared only for #dim#>=2).

    #active_quad_iterator#: loop over all active quads (declared only for
    #dim#>=2).

    Additionaly, for #dim#==1, the following identities hold:
    \begin{verbatim}
      typedef raw_line_iterator    raw_cell_iterator;
      typedef line_iterator        cell_iterator;
      typedef active_line_iterator active_cell_iterator;
    \end{verbatim}
    while for #dim#==2
    \begin{verbatim}
      typedef quad_line_iterator   raw_cell_iterator;    
      typedef quad_iterator        cell_iterator;
      typedef active_quad_iterator active_cell_iterator;
    \end{verbatim}

    By using the cell iterators, you can write code nearly independent of
    the spatial dimension.

    The #Triangulation# class offers functions like #begin_active# which gives
    you an iterator to the first active cell. There are quite a lot of functions
    returning iterators. Take a look at the class doc to get an overview.

    Usage of these iterators works mostly like with the STL iterators. Some
    examples taken from the #Triangulation# source code follow.
    \begin{itemize}
    \item {\it Counting the number of cells on a specific level}
      \begin{verbatim}
       template <int dim>
       int Triangulation<dim>::n_cells (const int level) const {
          cell_iterator cell = begin (level),
                        endc = (level == levels.size()-1 ?
                                cell_iterator(end()) :
                                begin (level+1));
          int n=0;
          for (; cell!=endc; ++cell)
            ++n;
          return n;
        };
      \end{verbatim}
      Another way which uses the STL #distance# function would be to write
      \begin{verbatim}
        template <int dim>
        int Triangulation<dim>::n_cells (const int level) const {
          int n=0;
          distance (begin(level),
                    (level == levels.size()-1 ?
                     cell_iterator(end()) :
                     begin (level+1)),
                    n);
          return n;
        };  
      \end{verbatim}
      Unfortunately, #g++# presently (version 2.7.2) fails to find the right
      #distance# template instantiation; it seems we have to wait for future
      #g++# versions :-(
      
    \item {\it Refining all cells of a triangulation}
      \begin{verbatim}
        template <int dim>
        void Triangulation<dim>::refine_global () {
          active_cell_iterator cell = begin_active(),
                               endc = end();

          for (; cell != endc; ++cell)
            cell->set_refine_flag ();
          execute_refinement ();
        };
      \end{verbatim}
    \end{itemize}


    {\bf Usage}

    Usage of a #Triangulation# is mainly done through the use of iterators.
    An example probably shows best how to use it:
    \begin{verbatim}
    void main () {
      Triangulation<2> tria;

      // read in a coarse grid file

                                       // we want to log the
                                       // refinement history
      ofstream history ("mesh.history");
      
                                       // refine first cell
      tria.begin_active()->set_refine_flag();
      tria.save_refine_flags (history);
      tria.execute_refinement ();

                                       // refine first active cell
                                       // on coarsest level
      tria.begin_active()->set_refine_flag ();
      tria.save_refine_flags (history);
      tria.execute_refinement ();
      
      Triangulation<2>::active_cell_iterator cell;
      for (int i=0; i<17; ++i) 
        {
                                           // refine the presently
                                           // second last cell 17
                                           // times
          cell = tria.last_active(tria.n_levels()-1);
          --cell;
          cell->set_refine_flag ();
          tria.save_refine_flags (history);
          tria.execute_refinement ();
        };
                                         // output the grid
      ofstream out("grid.1");
      tria.print_gnuplot (out);
    };  
    \end{verbatim}

    
    {\bf Creating a triangulation}

    There are several possibilities to create a triangulation:
    \begin{itemize}
      \item Hypercube triangulations: a hypercube triangulation is a
         domain which is the tensor product of an interval $[a,b]$ in
         the given number of spatial dimensions. If you want to create such
         a domain, which is a common test case for model problems, call
         #Triangulation<dim>::create_hypercube (a,b)#, which produces a
         hypercube domain triangulated with exactly one element. You can
         get tensor product meshes by successive refinement of this cell.

      \item Other standard regions: you can get the generalized L-shape domain
        using the #Triangulation<dim>::create_L_region (a,b)# function, which
	is the hypercube with the interval $[a,b]$ without the hypercube
	made out of the interval $[(a+b)/2,b]$. Let, for example, be $a=-1$
	and $b=1$, then the hpyer-L in two dimensions is the region
	$[-1,1]^2 - [0,1]^2$. To create a hyper-L in one dimension results in
	an error.

	You get the circle or ball (or generalized: hyperball) around origin
	#p# and with radius #r# by calling
	#Triangulation<dim>::create_hyper_ball (p, r)#. The circle is triangulated
	by five cells, the ball by seven cells. The diameter of the center cell is
	chosen so that the aspect ratio of the boundary cells after one refinement
	is minimized in some way. To create a hyperball in one dimension results in
	an error.

	Do not forget to attach a suitable
	boundary approximation object if you want the triangulation to be refined
	at the outer boundaries.
    
      \item Reading in a triangulation: By using an object of the \Ref{#DataIn#}
         class, you can read in fairly general triangulations. See there for
         more information. The mentionned class uses the interface described
         directly below to transfer the data into the triangulation.
    
      \item Explicitely creating a triangulation: you can create a triangulation
         by providing a list of vertices and a list of cells. Each such cell
         consists of a vector storing the indices of the vertices of this cell
         in the vertex list. To see how this works, you can take a look at the
         #DataIn<dim>::read_*# functions. The appropriate function to be
         called is #Triangulation<dim>::create_triangulation (2)#.

         Creating the hierarchical information needed for this library from
         cells storing only vertex information can be quite a complex task.
         For example in 2d, we have to create lines between vertices (but only
         once, though there are two cells which link these two vertices) and
         we have to create neighborship information. Grids being read in
         should therefore not be too large, reading refined grids would be
         inefficient. Apart from the performance aspect, refined grids do not
         lend too well to multigrid algorithms, since solving on the coarsest
         level is expensive. It is wiser in any case to read in a grid as coarse
         as possible and then do the needed refinement steps.

	 It is your duty to guarantee that cells have the correct orientation.
	 To guarantee this, in the input vector keeping the cell list, the
	 vertex indices for each cell have to be in a defined order. In one
	 dimension, the first vertex index must refer to that vertex with the
	 lower coordinate value. In two dimensions, the four vertices must be
	 given in an order representing a counterclockwise sense. This
	 condition is not easy to verify and no full attempt to do so is made.
	 If you violate this condition, you may end up with matrix entries
	 having the wrong sign (clockwise vertex numbering, which results in
	 a negative area element) of with wrong matrix elements (twisted
	 quadrilaterals, i.e. two vertices interchanged; this results in
	 a wrong area element).

	 There are more subtle conditions which must be imposed upon the
	 vertex numbering within cells. See the documentation for the
	 \Ref{DataIn} class for more details on this. They do not only
	 hold for the data read from an UCD or any other input file, but
	 also for the data passed to the
	 #Triangulation<dim>::create_triangulation (2)# function.
    \end{itemize}



    {\bf History of a triangulation}
    
    It is possible to reconstruct a grid from its refinement history, which
    can be stored and loaded through the #save_refine_flags# and
    #load_refine_flags# functions. Normally, the code will look like this:
    \begin{verbatim}
                                  // open output file
      ofstream history("mesh.history");
                                  // do 10 refinement steps
      for (int step=0; step<10; ++step) {
        ...;
        // flag cells according to some criterion
        ...;
        tria.save_refine_flags (history);
        tria.execute_refinement ();
      };        
    \end{verbatim}

    If you want to re-create the grid from the stored information, you write:
    \begin{verbatim}
                                  // open input file
      ifstream history("mesh.history");
                                  // do 10 refinement steps
      for (int step=0; step<10; ++step) {
        tria.load_refine_flags (history);
        tria.execute_refinement ();
      };        
    \end{verbatim}

    You may write other information to the output file between different sets
    of refinement information, as long as you read it upon re-creation of the
    grid.


    {\bf User flags}

    A triangulation offers one bit per line, quad, etc for user data.
    This field can be
    accessed as all other data using iterators. Normally, this user flag is
    used if an algorithm walks over all cells and needs information whether
    another cell, e.g. a neighbor, has already been processed. It can also
    be used to flag the lines subject to constraints in 2D, as for example the
    functions in the #DoFHandler# classes do.

    There are two functions, #save_user_flags# and #load_user_flags# which
    write and read these flags to and from a stream. Unlike
    #save_refine_flags# and #load_refine_flags#, these two functions store
    and read the flags of all used lines, quads, etc, not only of the
    active ones (well, activity is a concept which really only applies to
    cells, not for example to lines in 2D, so the abovementioned generalisation
    to {\it all} lines, quads, etc seems plausible).

    If you want to store more specific user flags, you can use the functions
    #save_user_flags_line# and #load_user_flags_line# and the generalizations
    for quads, etc.

    
    It is convention to clear the user flags using the
    #Triangulation<>::clear_user_flags()# function before usage, since it is
    often necessary to use the flags in more than one function consecutively and
    is then error prone to dedicate one of these to clear the flags.

    It is recommended that a functions using the flags states so in its
    documentation.
    
    
    {\bf Boundary approximation}

    You can specify a boundary function: if a new vertex is created on a
    side or face at the boundary, this function is used to compute where
    it will be placed. See \Ref{Boundary} for the details. Usage with
    the #Triangulation# object is then like this (let #Ball# be a class
    derived from #Boundary<2>#):
    \begin{verbatim}
      void main () {
        Triangulation<2> tria;
                                         // set the boundary function
        Ball ball;
        tria.set_boundary (&ball);

        // read some coarse grid

        
        Triangulation<2>::active_cell_iterator cell, endc;
        for (int i=0; i<8; ++i) 
          {
            cell = tria.begin_active();
            endc = tria.end();

                                             // refine all
                                             // boundary cells
            for (; cell!=endc; ++cell)
              if (cell->at_boundary())
                cell->set_refine_flag();

            tria.execute_refinement();
          };
      };            
    \end{verbatim}

    You should take note of one caveat: if you have concave boundaries, you
    must make sure that a new boundary vertex does not lie to much inside the
    to be refined cell. The reason is that the center vertex is placed at the
    point which is the arithmetic mean of the eight surrounding vertices.
    Therefore if your new boundary vertex is too near the center of the old
    quadrilateral or hexahedron, the distance to the midpoint vertex will become
    too small, thus generating distorted cells. Remedy: you have to take care
    of such situations when defining the coarse grid.


    {\bf Implementational conventions for two spatial dimensions}

    There is a convention about the direction of the bounding lines of quads in
    2D. The direction of a line is the direction of point 0 towards point 1. We
    define, that allowed cells contain of lines of which the direction is
    as follows:

           2
       3--->---2
       |       |
      3^       ^1
       |       |
       0--->---1
           0
    The number of the vertices and lines is also indicated. This orientation of
    lines has to be checked/generated upon construction of a grid and is
    preserved upon refinement.

    Further we define, that child lines have the same direction as their parent,
    i.e. that #subline(0).vertex(0)==line.vertex(0)# and
    #subline(1).vertex(1)==line.vertex(1)#. This also implies, that the
    first subline (#subline(0)#) is the one at vertex(0) of the old line.

    Similarly we define, that the four children of a quad are adjacent to the
    vertex with the same number of the old quad.

    
    {\bf Warning}
    
    It seems impossible to preserve #const#ness of a triangulation through
    iterator usage. Thus, if you declare pointers to a #const# triangulation
    object, you should be well aware that you might involuntarily alter the
    data stored in the triangulation.

    @memo Implementation of a multilevel triangulation of a domain
    @see TriaRawIterator
    @author Wolfgang Bangerth, 1998
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

				     /**
				      *  Create a triangulation and create
				      *  the first level of the hierarchy.
				      *  Do not create any cells.
				      */
    Triangulation ();

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
				      * The numbering of vertices within the
				      * #cells# array is subject to some
				      * constraints; see the general class
				      * documentation for this.
				      */
    void create_triangulation (const vector<Point<dim> >  &vertices,
			       const vector<vector<int> > &cells);
    
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
    DeclException0 (ExcInternalError);
				     //@}
  protected:
				     /**
				      *  Prepare the refinement process by
				      *  allocating enough space, fixing the
				      *  closure of the refinement in #dim>=2#
				      *  (make sure that no two cells are
				      *  adjacent with a refinement level
				      *  differing with more than one), etc.
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





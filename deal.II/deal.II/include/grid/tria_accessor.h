/*----------------------------   tria_accessor.h     ---------------------------*/
/*      $Id$                 */
#ifndef __tria_accessor_H
#define __tria_accessor_H
/*----------------------------   tria_accessor.h     ---------------------------*/


#include <base/exceptions.h>


// forward declaration needed
class Line;
class Quad;
template <int dim> class Point;


template <int dim> class TriaAccessor;
template <int dim> class LineAccessor;
template <int dim> class QuadAccessor;
template <int dim> class CellAccessor;

template <int dim> class DoFCellAccessor;
template <int dim, class BaseClass> class DoFLineAccessor;

template <int dim, class Accessor> class TriaRawIterator;
template <int dim, class Accessor> class TriaIterator;
template <int dim, class Accessor> class TriaActiveIterator;

template <int dim> class Triangulation;







/**
    The three states an iterator can be in: valid, past-the-end and
    invalid.
    */
enum IteratorState { valid, past_the_end, invalid };





/**
    Implements the accessor class descibed in the documentation of
    the iterator classes (see \Ref{TriaRawIterator}.

    This class offers only the basic functionality (stores the necessary
    data members, offers comparison operators and the like), but has no
    functionality to actually dereference data. This is done in the derived
    classes.
 */
template <int dim>
class TriaAccessor {
  protected:
				     /**
				      *  Constructor. Protected, thus
				      *  only callable from friend
				      *  classes.
				      */
    TriaAccessor (Triangulation<dim> *parent     = 0,
		  const int           level      = -1,
		  const int           index      = -1,
		  const void         * = 0) :
		    present_level (level),
		    present_index (index),
		    tria (parent) {};


	
				     /**
				      *  Copy operator. Since this is only
				      *  called from iterators, do not
				      *  return anything, since the
				      *  iterator will return itself.
				      *
				      *  This method is protected, since it
				      *  is only to be called from the
				      *  iterator class.
				      */
    void copy_from (const TriaAccessor &);

				     /**
				      *  Copy operator. This is normally
				      *  used in a context like
				      *  #iterator a,b;  *a=*b;#. Since
				      *  the meaning is to copy the object
				      *  pointed to by #b# to the object
				      *  pointed to by #a# and since
				      *  accessors are not real but
				      *  virtual objects, this operation
				      *  is not useful for iterators on
				      *  triangulations. We declare this
				      *  function here private, thus it may
				      *  not be used from outside.
				      *  Furthermore it is not implemented
				      *  and will give a linker error if
				      *  used anyway.
				      */
    void operator = (const TriaAccessor *);
    
				     /**
				      *  Same as above.
				      */
    TriaAccessor &operator = (const TriaAccessor &);
	
  public:
				     /**
				      *  Compare for equality.            
				      *
				      *  This method should be protected,
				      *  since it is only to be called from
				      *  the iterator class. Due to problems
				      *  with the STL, we have to make it
				      *  public, so don't use it from
				      *  non-friend classes!
				      */
    bool operator == (const TriaAccessor &) const;
	
				     /**
				      *  Compare for inequality.
				      *
				      *  This method should be protected,
				      *  since it is only to be called from
				      *  the iterator class. Due to problems
				      *  with the STL, we have to make it
				      *  public, so don't use it from
				      *  non-friend classes!
				      */
    bool operator != (const TriaAccessor &) const;
    
  public:
    typedef void* LocalData;
				     /**@ name Iterator address and state
				      */
				     /*@{*/
				     /**
				      *  Returns the level the element
				      *  pointed to belongs to.
				      */
    int level () const;
    
				     /**
				      *  Returns the index of the element
				      *  presently pointed to on the
				      *  present level.
				      */
    int index () const;
				     /**
				      *  Return the state of the iterator.
				      *  For the different states an accessor 
				      *  can be in, refer to the
				      *  \Ref{TriaRawIterator} documentation.
				      */
    IteratorState state () const;
				     /*@}*/

				     /**@name Exceptions for derived classes
				      */
				     /*@{*/
				     /**
				      *  Exception
				      */
    DeclException0 (ExcCellNotUsed);
				     /**
				      *  Exception
				      */
    DeclException0 (ExcRefineCellNotUsed);
				     /**
				      *  Exception
				      */
    DeclException0 (ExcRefineCellNotActive);
				     /**
				      *  Exception
				      */
    DeclException1 (ExcInvalidNeighbor,
		    int,
		    << "Neighbor indices must be between 0 and 2^dim-1, but "
		    << "yours was " << arg1);
				     /**
				      *  Exception
				      */
    DeclException0 (ExcUnusedCellAsChild);
				     /**
				      *  Exception
				      */
    DeclException0 (ExcUnusedCellAsNeighbor);
				     /**
				      *  Exception
				      */
    DeclException3 (ExcInvalidIndex,
		    int,
		    int,
		    int,
		    << "Invalid index " << arg1
		    << ", index must be between " << arg2
		    << " and " << arg3 << ".");
				     /**
				      *  Exception
				      */
    DeclException0 (ExcUncaughtCase);
				     /**
				      *  Exception
				      */
    DeclException0 (ExcDereferenceInvalidObject);
				     /*@}*/
	
  protected:
				     /**
				      *  Used to store the level presently
				      *  pointed to.
				      */
    int present_level;
    
				     /**
				      *  Used to store the index of the
				      *  element presently pointed to on
				      *  the level presentl used.
				      */
    int present_index;
    
				     /**
				      *  Pointer to the triangulation which
				      *  we act on.
				      */
    Triangulation<dim> *tria;

    friend class TriaRawIterator<1,LineAccessor<1> >;
    friend class TriaRawIterator<1,CellAccessor<1> >;
    friend class TriaRawIterator<2,LineAccessor<2> >;
    friend class TriaRawIterator<2,QuadAccessor<2> >;
    friend class TriaRawIterator<2,CellAccessor<2> >;

    friend class TriaIterator<1,LineAccessor<1> >;
    friend class TriaIterator<1,CellAccessor<1> >;
    friend class TriaIterator<2,LineAccessor<2> >;
    friend class TriaIterator<2,QuadAccessor<2> >;
    friend class TriaIterator<2,CellAccessor<2> >;

    friend class TriaActiveIterator<1,LineAccessor<1> >;
    friend class TriaActiveIterator<1,CellAccessor<1> >;
    friend class TriaActiveIterator<2,LineAccessor<2> >;
    friend class TriaActiveIterator<2,QuadAccessor<2> >;
    friend class TriaActiveIterator<2,CellAccessor<2> >;

    friend class TriaRawIterator<1,DoFCellAccessor<1> >;
    friend class TriaRawIterator<2,DoFCellAccessor<2> >;
    friend class TriaRawIterator<2,DoFLineAccessor<2,LineAccessor<2> > >;
    friend class TriaIterator<1,DoFCellAccessor<1> >;
    friend class TriaIterator<2,DoFCellAccessor<2> >;
    friend class TriaIterator<2,DoFLineAccessor<2,LineAccessor<2> > >;
};







/**
    Accessor to dereference the data of lines. This accessor is used to
    point to lines in #dim# space dimensions. There is a derived class
    for lines in one space dimension, in which case a line is also a cell
    and thus has much more functionality than in lower dimensions.
    */
template <int dim>
class LineAccessor :  public TriaAccessor<dim> {
  public:
				     /**
				      *  Constructor.
				      */
    LineAccessor (Triangulation<dim> *parent     = 0,
		  const int           level      = -1,
		  const int           index      = -1,
		  const void         *local_data = 0) :
		    TriaAccessor<dim> (parent, level, index, local_data) {};

				     /**
				      *  Copy the data of the given line.
				      */
    void set (const Line &l) const;

				     /**
				      *  Return the index of vertex #i=0,1#
				      *  of a line.
				      */ 
    int vertex_index (const unsigned int i) const;

				     /**
				      *  Return a reference (not an iterator!)
				      *  to the #i#th vertex.
				      */
    Point<dim> & vertex (const unsigned int i) const;

				     /**
				      *  Test for the element being used
				      *  or not.
				      */
    bool used () const;

				     /**
				      *  Set the #used# flag.
				      */
    void set_used_flag () const;

				     /**
				      *  Clear the #used# flag.
				      */
    void clear_used_flag () const;

    				     /**
				      *  Return whether the user flag
				      *  is set or not.
				      */
    bool user_flag_set () const;

				     /**
				      *  Flag the user flag for this cell.
				      */
    void set_user_flag () const;

				     /**
				      *  Clear the user flag.
				      */
    void clear_user_flag () const;

				     /**
				      *  Return a pointer to the #i#th
				      *  child.
				      */
    TriaIterator<dim,LineAccessor<dim> > child (const unsigned int i) const;

				     /**
				      *  Return the index of the #i#th child.
				      *  The level of the child is one higher
				      *  than that of the present cell.
				      *  If the child does not exist, -1
				      *  is returned.
				      */
    int child_index (const unsigned int i) const;

				     /**
				      *  Set the child field. Since we
				      *  only store the index of the first
				      *  child (the others follow directly)
				      *  only one child index is to be
				      *  given. The level of the child is
				      *  one level up of the level of the
				      *  cell to which this iterator points.
				      */
    void set_children (const int index) const;
	
				     /**
				      *  Clear the child field, i.e. set
				      *  it to a value which indicates
				      *  that this cell has no children.
				      */
    void clear_children () const;

				     /**
				      *  Test whether the line has children.
				      */
    bool has_children () const;

  private:
    				     /**
				      *  Copy operator. This is normally
				      *  used in a context like
				      *  #iterator a,b;  *a=*b;#. Since
				      *  the meaning is to copy the object
				      *  pointed to by #b# to the object
				      *  pointed to by #a# and since
				      *  accessors are not real but
				      *  virtual objects, this operation
				      *  is not useful for iterators on
				      *  triangulations. We declare this
				      *  function here private, thus it may
				      *  not be used from outside.
				      *  Furthermore it is not implemented
				      *  and will give a linker error if
				      *  used anyway.
				      */
    void operator = (const LineAccessor &);

  public:
//  protected:
				     // would be better to make this private,
				     // but that would require making all
				     // TriaRawIterator<...>s friend.->wait for
				     // gcc 2.30687
    
				     /**@name Advancement of iterators*/
				     /*@{*/
				     /**
				      *  This operator advances the iterator to
				      *  the next element.
				      *
				      *  The next element is next on this
				      *  level if there are more. If the
				      *  present element is the last on
				      *  this level, the first on the
				      *  next level is accessed.
				      */
    void operator ++ ();

    				     /**
				      *  This operator moves the iterator to
				      *  the previous element.
				      *
				      *  The previous element is previous on
				      *  this level if #index>0#. If the
				      *  present element is the first on
				      *  this level, the last on the
				      *  previous level is accessed.
				      */
    void operator -- ();
				     /*@}*/
};



/**
    Accessor to dereference the data of quads. This accessor is used to
    point to quads in #dim# space dimensions (only #dim>=2# seems reasonable
    to me). There is a derived class
    for quads in two space dimension, in which case a quad is also a cell
    and thus has much more functionality than in lower dimensions.
    */
template <int dim>
class QuadAccessor :  public TriaAccessor<dim> {
  public:
				     /**
				      *  Constructor.
				      */
    QuadAccessor (Triangulation<dim> *parent     = 0,
		  const int           level      = -1,
		  const int           index      = -1,
		  const void         *local_data = 0) :
		    TriaAccessor<dim> (parent, level, index, local_data) {};

				     /**
				      *  Copy the data of the given quad.
				      */
    void set (const Quad &q) const;
    
				     /**
				      *  Return index of vertex #i=0,1,2,3# of
				      *  a quad. The #i#th vertex is the common
				      *  one of line #i# and #(i+3)%4#. See also
				      *  the introduced convention
				      *  (\Ref{Triangulation}).
				      */ 
    int vertex_index (const unsigned int i) const;

    				     /**
				      *  Return a reference (not an iterator!)
				      *  to the #i#th vertex.
				      */
    Point<dim> & vertex (const unsigned int i) const;

				     /**
				      *  Return a pointer to the #i#th line
				      *  bounding this #Quad#.
				      */
    TriaIterator<dim,LineAccessor<dim> > line (const unsigned int i) const;

				     /**
				      * Return the line index of the #i#th
				      * side (a line). The level is naturally
				      * the same as that of the quad.
				      */
    unsigned int line_index (const unsigned int i) const;
    
				     /**
				      *  Test for the element being used
				      *  or not.
				      */
    bool used () const;

				     /**
				      *  Set the #used# flag.
				      */
    void set_used_flag () const;

				     /**
				      *  Clear the #used# flag.
				      */
    void clear_used_flag () const;

    				     /**
				      *  Return whether the user flag
				      *  is set or not.
				      */
    bool user_flag_set () const;

				     /**
				      *  Flag the user flag for this cell.
				      */
    void set_user_flag () const;

				     /**
				      *  Clear the user flag.
				      */
    void clear_user_flag () const;

				     /**
				      *  Return a pointer to the #i#th
				      *  child.
				      */
    TriaIterator<dim,QuadAccessor<dim> > child (const unsigned int i) const;

				     /**
				      *  Return the index of the #i#th child.
				      *  The level of the child is one higher
				      *  than that of the present cell.
				      *  If the child does not exist, -1
				      *  is returned.
				      */
    int child_index (const unsigned int i) const;

				     /**
				      *  Set the child field. Since we
				      *  only store the index of the first
				      *  child (the others follow directly)
				      *  only one child index is to be
				      *  given. The level of the child is
				      *  one level up of the level of the
				      *  cell to which this iterator points.
				      */
    void set_children (const int index) const;
	
				     /**
				      *  Clear the child field, i.e. set
				      *  it to a value which indicates
				      *  that this cell has no children.
				      */
    void clear_children () const;
    
				     /**
				      *  Test whether the quad has children.
				      */
    bool has_children () const;

  private:
    				     /**
				      *  Copy operator. This is normally
				      *  used in a context like
				      *  #iterator a,b;  *a=*b;#. Since
				      *  the meaning is to copy the object
				      *  pointed to by #b# to the object
				      *  pointed to by #a# and since
				      *  accessors are not real but
				      *  virtual objects, this operation
				      *  is not useful for iterators on
				      *  triangulations. We declare this
				      *  function here private, thus it may
				      *  not be used from outside.
				      *  Furthermore it is not implemented
				      *  and will give a linker error if
				      *  used anyway.
				      */
    void operator = (const QuadAccessor &);

  public:
//  protected:
				     // would be better to make this private,
				     // but that would require making all
				     // TriaRawIterator<...>s friend.->wait for
				     // gcc 2.30687

				     /**@name Advancement of iterators*/
				     /*@{*/
				     /**
				      *  This operator advances the iterator to
				      *  the next element.
				      *
				      *  The next element is next on this
				      *  level if there are more. If the
				      *  present element is the last on
				      *  this level, the first on the
				      *  next level is accessed.
				      */
    void operator ++ ();

    				     /**
				      *  This operator moves the iterator to
				      *  the previous element.
				      *
				      *  The previous element is previous on
				      *  this level if #index>0#. If the
				      *  present element is the first on
				      *  this level, the last on the
				      *  previous level is accessed.
				      */
    void operator -- ();
				     /*@}*/
};




/**
  This class allows access to a {\bf cell}, which is a line in 1D and a quad in
  2D. Declare it to have a template parameter, but do not actually declare
  other types than those explicitely instantiated.
   */
template <int dim>
class CellAccessor;





/**
  This class allows access to a cell, i.e. a line on the present dimension.

  The following refers to any space dimension:
  
  This class allows access to a {\bf cell}, which is a line in 1D and a quad in
  2D. Cells have more functionality than lines or quads by themselves, for
  example they can be flagged for refinement, they have neighbors, they have
  the possibility to check whether they are at the boundary etc. This class
  offers access to all this data.
  
  The are specialized versions, \Ref{CellAccessor<1>} and
  \Ref{CellAccessor<2>}, which
  offer access to cells and one and two dimensions respectively. They have
  more or less the same functionality, but return different types in some
  cases.
 */
class CellAccessor<1> :  public LineAccessor<1> {
  public:
				     /**
				      *  Constructor.
				      */
    CellAccessor (Triangulation<1>   *parent     = 0,
		  const int           level      = -1,
		  const int           index      = -1,
		  const void         *local_data = 0) :
		    LineAccessor<1> (parent, level, index, local_data) {};

				     /**
				      *  Return a pointer to the #i#th
				      *  neighbor.
				      *  If the neighbor does not exist, an
				      *  invalid iterator is returned.
				      */
    TriaIterator<1,CellAccessor<1> > neighbor (const unsigned int i) const;

				     /**
				      *  Return the index of the #i#th neighbor.
				      *  If the neighbor does not exist, its
				      *  index is -1.
				      */
    int neighbor_index (const unsigned int i) const;

    				     /**
				      *  Return the level of the #i#th neighbor.
				      *  If the neighbor does not exist, its
				      *  level is -1.
				      */
    int neighbor_level (const unsigned int i) const;

				     /**
				      *  Set the neighbor #i# of this cell
				      *  to the cell pointed to by #pointer#.
				      *  This line must be used.
				      */
    void set_neighbor (const unsigned int i,
		       const TriaIterator<1,CellAccessor<1> > &pointer) const;

				     /**
				      *  Return whether the #i#th vertex is
				      *  part of the boundary. This is true, if
				      *  the #i#th neighbor does not exist.
				      */
    bool at_boundary (const unsigned int i) const;

				     /**
				      *  Return whether the cell is at the
				      *  boundary.
				      */
    bool at_boundary () const;

				     /**
				      *  Return whether the refinement flag
				      *  is set or not.
				      */
    bool refine_flag_set () const;

				     /**
				      *  Flag the cell pointed to
				      *  for refinement.
				      */
    void set_refine_flag () const;

				     /**
				      *  Clear the refinement flag.
				      */
    void clear_refine_flag () const;

				     /**
				      *  Return a pointer to the #i#th
				      *  child. Overloaded version which returns
				      *  a more reasonable iterator class.
				      */
    TriaIterator<1,CellAccessor<1> > child (const unsigned int i) const;

				     /**
				      * Test whether the cell has children
				      * (this is the criterion for activity
				      * of a cell).
				      */
    bool active () const;
    
  private:
    				     /**
				      *  Copy operator. This is normally
				      *  used in a context like
				      *  #iterator a,b;  *a=*b;#. Since
				      *  the meaning is to copy the object
				      *  pointed to by #b# to the object
				      *  pointed to by #a# and since
				      *  accessors are not real but
				      *  virtual objects, this operation
				      *  is not useful for iterators on
				      *  triangulations. We declare this
				      *  function here private, thus it may
				      *  not be used from outside.
				      *  Furthermore it is not implemented
				      *  and will give a linker error if
				      *  used anyway.
				      */
    void operator = (const CellAccessor<1> &);
};






/**
    Specialization of a #CellAccessor# for two space dimensions.
    @see CellAccessor<1>
 */
class CellAccessor<2> :  public QuadAccessor<2> {
  public:
				     /**
				      *  Constructor.
				      */
    CellAccessor (Triangulation<2>   *parent     = 0,
		  const int           level      = -1,
		  const int           index      = -1,
		  const void         *local_data = 0) :
		    QuadAccessor<2> (parent, level, index, local_data) {};

				     /**
				      *  Return a pointer to the #i#th
				      *  neighbor.
				      */
    TriaIterator<2,CellAccessor<2> > neighbor (const unsigned int i) const;

				     /**
				      *  Return the index of the #i#th neighbor.
				      *  If the neighbor does not exist, its
				      *  index is -1.
				      */
    int neighbor_index (const unsigned int i) const;

    				     /**
				      *  Return the level of the #i#th neighbor.
				      *  If the neighbor does not exist, its
				      *  level is -1.
				      */
    int neighbor_level (const unsigned int i) const;

				     /**
				      *  Set the neighbor #i# of this cell
				      *  to the cell pointed to by
				      *  #pointer#.
				      */
    void set_neighbor (const unsigned int i,
		       const TriaIterator<2,CellAccessor<2> > &pointer) const;

				     /**
				      *  Return whether the #i#th side of this
				      *  quad is part of the
				      *  boundary. This is true, if the
				      *  #i#th neighbor does not exist.
				      */
    bool at_boundary (const unsigned int i) const;

				     /**
				      *  Return whether the cell is at the
				      *  boundary.
				      */
    bool at_boundary () const;

				     /**
				      *  Return whether the refinement flag
				      *  is set or not.
				      */
    bool refine_flag_set () const;

				     /**
				      *  Flag the cell pointed to
				      *  for refinement.
				      */
    void set_refine_flag () const;

				     /**
				      *  Clear the refinement flag.
				      */
    void clear_refine_flag () const;
    
				     /**
				      *  Return a pointer to the #i#th
				      *  child. Overloaded version which returns
				      *  a more reasonable iterator class.
				      */
    TriaIterator<2,CellAccessor<2> > child (const unsigned int i) const;
				     /**
				      * Test whether the cell has children
				      * (this is the criterion for activity
				      * of a cell).
				      */
    bool active () const;

  private:
    				     /**
				      *  Copy operator. This is normally
				      *  used in a context like
				      *  #iterator a,b;  *a=*b;#. Since
				      *  the meaning is to copy the object
				      *  pointed to by #b# to the object
				      *  pointed to by #a# and since
				      *  accessors are not real but
				      *  virtual objects, this operation
				      *  is not useful for iterators on
				      *  triangulations. We declare this
				      *  function here private, thus it may
				      *  not be used from outside.
				      *  Furthermore it is not implemented
				      *  and will give a linker error if
				      *  used anyway.
				      */
    void operator = (const CellAccessor<2> &);
};




/*----------------------------   tria_accessor.h     ---------------------------*/
/* end of #ifndef __tria_accessor_H */
#endif
/*----------------------------   tria_accessor.h     ---------------------------*/

/*----------------------------   tria_accessor.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __tria_accessor_H
#define __tria_accessor_H
/*----------------------------   tria_accessor.h     ---------------------------*/


#include <base/exceptions.h>
#include <base/forward-declarations.h>
#include <basic/forward-declarations.h>
#include <grid/tria_iterator_base.h>


// note: in non-debug mode, i.e. with optimizations, the file
// tria_accessor.templates.h is included at the end of this file.
// this includes a lot of templates and thus makes compilation
// slower, but at the same time allows for more aggressive
// inlining and thus faster code.






/**
 *   Implements the accessor class descibed in the documentation of
 *   the iterator classes (see \Ref{TriaRawIterator}.
 *
 *   This class offers only the basic functionality (stores the necessary
 *   data members, offers comparison operators and the like), but has no
 *   functionality to actually dereference data. This is done in the derived
 *   classes.
 *
 *   @author Wolfgang Bangerth, 1998
 */
template <int dim>
class TriaAccessor {
  protected:
				     /**
				      * Declare the data type that this accessor
				      * class expects to get passed from the
				      * iterator classes. Since the pure
				      * triangulation iterators need no
				      * additional data, this data type is
				      * #void#.
				      */
    typedef void AccessorData;
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

  protected:
    
				     /**
				      *  Compare for equality.            
				      */
    bool operator == (const TriaAccessor &) const;
	
				     /**
				      *  Compare for inequality.
				      *
				      * Note that at times, there is a problem
				      * with egcs 1.1 that makes it choose
				      * the global STL operator != (which
				      * does only !(a==b)) over this
				      * one, which then results in an
				      * error because the operator == is
				      * not made public. Strange...
				      */
    bool operator != (const TriaAccessor &) const;
    
  public:
				     /**
				      * Data type to be used for passing
				      * parameters from iterators to the
				      * accessor classes in a unified
				      * way, no matter what the type of
				      * number of these parameters is.
				      */
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

				     /**
				      * Return a pointer to the triangulation
				      * which the object pointed to by this
				      * class belongs to.
				      */
    const Triangulation<dim> & get_triangulation () const;
    
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
    DeclException1 (ExcCantSetChildren,
		    int,
		    << "You can only set the child index if the cell has no "
		    << "children, or clear it. The given "
		    << "index was " << arg1 << " (-1 means: clear children)");
				     /**
				      * Exception
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
				     /**
				      *  Exception
				      */
    DeclException0 (ExcNotUsefulForThisDimension);
				     /**
				      * Exception
				      */
    DeclException0 (ExcCantCompareIterators);
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

    template <int dim, typename Accessor> friend class TriaRawIterator<dim,Accessor>;
    template <int dim, typename Accessor> friend class TriaIterator<dim,Accessor>;
    template <int dim, typename Accessor> friend class TriaActiveIterator<dim,Accessor>;
};







/**
 *   Accessor to dereference the data of lines. This accessor is used to
 *   point to lines in #dim# space dimensions. There is a derived class
 *   for lines in one space dimension, in which case a line is also a cell
 *   and thus has much more functionality than in other dimensions.
 *
 *   @author Wolfgang Bangerth, 1998
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
				      *  set the user flag of this object
				      *  and of all its children and their
				      *  children, etc.
				      */
    void recursively_set_user_flag () const;

    				     /**
				      *  Clear the user flag of this object
				      *  and of all its children and their
				      *  children, etc.
				      */
    void recursively_clear_user_flag () const;

				     /**
				      * Set the user pointer of this line
				      * to #p#.
				      */
    void set_user_pointer (void *p) const;

				     /**
				      * Reset the user pointer of this line
				      * to a #NULL# pointer.
				      */
    void clear_user_pointer () const;

				     /**
				      * Access the value of the user pointer
				      * of this line. It is in the
				      * responsibility of the user to make
				      * sure that the pointer points to
				      * something useful. You should use the
				      * new style cast operator to maintain
				      * a minimum of typesafety, e.g.
				      * #A *a=static_cast<A*>(cell->user_pointer());#.
				      */
    void * user_pointer () const;
    
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

				     /**
				      * Return the number of times that this
				      * cell is refined. Note that not all
				      * its children are refined that often
				      * (which is why we prepend #max_#), 
				      * the returned number is rather the
				      * maximum number of refinement in
				      * any branch of children of this object.
				      */
    unsigned int max_refinement_depth () const;
    
				     /**
				      * Return the boundary indicator of this
				      * line. Since boundary data is only useful
				      * for structures with a dimension less
				      * than the dimension of a cell, this
				      * function issues an error if #dim<2#.
				      *
				      * If the return value is 255, then this
				      * line is in the interior of the domain.
				      */
    unsigned char boundary_indicator () const;

				     /**
				      * Set the boundary indicator of this line.
				      * The same applies as for the
				      * #boundary_indicator()# function.
				      *
				      * You should be careful with this function
				      * and especially never try to set the
				      * boundary indicator to 255, unless
				      * you exactly know what you are doing,
				      * since this value is reserved for
				      * another purpose and algorithms may
				      * not work if boundary cells have
				      * this boundary indicator or if interior
				      * cells have boundary indicators other
				      * than 255.
				      */
    void set_boundary_indicator (unsigned char) const;

				     /**
				      * Return whether this line is at the
				      * boundary. This is checked via the
				      * the boundary indicator field, which
				      * is always 255 if the line is in the
				      * interior of the domain. Obviously,
				      * this is only possible for #dim>1#;
				      * however, for #dim==1#, a line is
				      * a cell and the #CellAccessor# class
				      * offers another possibility to
				      * determine whether a cell is at the
				      * boundary or not.
				      */
    bool at_boundary () const;

    				     /**
				      * Return the length of the line. If the
				      * line describes part of the boundary
				      * (e.g. if it is face to a cell in 2D)
				      * and is not a straight one, ask the
				      * finite element class for the correct
				      * length!
				      */
    double diameter () const;

				     /**
				      * Return the center of the line. This
				      * is the average of the two vertices,
				      * which is the obvious definition for
				      * straight lines. However, if you use
				      * higher order mappings from the unit
				      * cell to the real cell (in more than
				      * one space dimension), the bounding
				      * lines may not necessarily be straight.
				      * In that case ask the finite element
				      * class for the correct place of the
				      * midpoint of the line in real space.
				      */
    Point<dim> center () const;

				     /**
				      * Return the barycenter of the line,
				      * which is the midpoint. The same
				      * applies as for the #center# function
				      * with regard to lines at the boundary.
				      */
    Point<dim> barycenter () const;
    
				     /**
				      * Return the length of the line.
				      * The same
				      * applies as for the #center# function
				      * with regard to lines at the boundary.
				      */
    double measure () const;

				     /**
				      * Compute and return the number of
				      * children of this line. Actually,
				      * this function only counts the number
				      * of active children, i.e. the number
				      * if lines which are not further
				      * refined. Thus, if both of the two
				      * children of a line are further
				      * refined exactly once, the returned
				      * number will be four, not six.
				      *
				      * If the present cell is not refined,
				      * one is returned.
				      */
    unsigned int number_of_children () const;
    
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

    

  protected:
    
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

				     /**
				      * Declare some friends.
				      */
    template <int dim> friend class TriaRawIterator<dim,LineAccessor<dim> >;
    template <>        friend class TriaRawIterator<1,CellAccessor<1> >;
};



/**
 *   Accessor to dereference the data of quads. This accessor is used to
 *   point to quads in #dim# space dimensions (only #dim>=2# seems reasonable
 *   to me). There is a derived class
 *   for quads in two space dimension, in which case a quad is also a cell
 *   and thus has much more functionality than in other dimensions.
 *
 *   @author Wolfgang Bangerth, 1998
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
				      *  Return index of vertex #i=0 through 3# of
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
				      *  Set the #used# flag. You should know
				      *  quite exactly what you are doing of you
				      *  touch this function. It is exclusively
				      *  for internal use in the library.
				      */
    void set_used_flag () const;

				     /**
				      *  Clear the #used# flag. You should know
				      *  quite exactly what you are doing of you
				      *  touch this function. It is exclusively
				      *  for internal use in the library.
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
				      *  set the user flag of this object
				      *  and of all its children and their
				      *  children, etc.
				      */
    void recursively_set_user_flag () const;

    				     /**
				      *  Clear the user flag of this object
				      *  and of all its children and their
				      *  children, etc.
				      */
    void recursively_clear_user_flag () const;

				     /**
				      * Set the user pointer of this quad
				      * to #p#.
				      */
    void set_user_pointer (void *p) const;

				     /**
				      * Reset the user pointer of this quad
				      * to a #NULL# pointer.
				      */
    void clear_user_pointer () const;

				     /**
				      * Access the value of the user pointer
				      * of this quad. It is in the
				      * responsibility of the user to make
				      * sure that the pointer points to
				      * something useful. You should use the
				      * new style cast operator to maintain
				      * a minimum of typesafety, e.g.
				      * #A *a=static_cast<A*>(cell->user_pointer());#.
				      */
    void * user_pointer () const;

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

				     /**
				      * Return the number of times that this
				      * cell is refined. Note that not all
				      * its children are refined that often
				      * (which is why we prepend #max_#), 
				      * the returned number is rather the
				      * maximum number of refinement in
				      * any branch of children of this object.
				      */
    unsigned int max_refinement_depth () const;
    
				     /**
				      * Return the boundary indicator of this
				      * quad. Since boundary data is only useful
				      * for structures with a dimension less
				      * than the dimension of a cell, this
				      * function issues an error if #dim<3#.
				      *
				      * If the return value is 255, then this
				      * quad is in the interior of the domain.
				      */
    unsigned char boundary_indicator () const;

				     /**
				      * Set the boundary indicator of this quad.
				      * The same applies as for the
				      * #boundary_indicator()# function.
				      *
				      * You should be careful with this function
				      * and especially never try to set the
				      * boundary indicator to 255, unless
				      * you exactly know what you are doing,
				      * since this value is reserved for
				      * another purpose and algorithms may
				      * not work if boundary cells have have
				      * this boundary indicator or if interior
				      * cells have boundary indicators other
				      * than 255.
				      */
    void set_boundary_indicator (unsigned char) const;

				     /**
				      * Return whether this line is at the
				      * boundary. This is checked via the
				      * the boundary indicator field, which
				      * is always 255 if the quad is in the
				      * interior of the domain. Obviously,
				      * this function is only useful for #dim>2#;
				      * however, for #dim==2#, a quad is
				      * a cell and the #CellAccessor# class
				      * offers another possibility to
				      * determine whether a cell is at the
				      * boundary or not.
				      */
    bool at_boundary () const;
				     /**
				      * Return the diameter of the quad. If the
				      * quad describes part of the boundary
				      * (e.g. if it is face to a cell in 3D)
				      * and is not a plane one, ask the
				      * finite element class for the correct
				      * length!
				      *
				      * The diameter of a quad is computed to
				      * be the larger of the two diagonals. You
				      * should absolutely be clear about the
				      * fact that this definitely is not the
				      * diameter of all quadrilaterals; however
				      * it may serve as an approximation and
				      * is exact in many cases, especially
				      * if the quadrilateral is not too much
				      * distorted.
				      */
    double diameter () const;

    				     /**
				      * Return the center of the quad. The
				      * center of a quad is defined to be
				      * the average of the four vertices,
				      * which is also the point where the
				      * bilinear mapping places the midpoint
				      * of the unit quad in real space.
				      * However, this may not be the point
				      * of the barycenter of the quad.
				      *
				      * Also note that if you use
				      * higher order mappings from the unit
				      * cell to the real cell (in more than
				      * two space dimension), the bounding
				      * quads may not necessarily be straight.
				      * In that case ask the finite element
				      * class for the correct place of the
				      * midpoint of the quad in real space.
				      */
    Point<dim> center () const;

				     /**
				      * Return the barycenter of the qaud. The
				      * same applies as for the #center# function
				      * with regard to quads at the boundary.
				      */
    Point<dim> barycenter () const;

				     /**
				      * Return the area of the quad. With
				      * area, we mean the area included by
				      * the straight lines connecting the
				      * four vertices, i.e. no information
				      * about the boundary is used; if
				      * you want other than bilinearly
				      * mapped unit quadrilaterals, ask the
				      * appropriate finite element class
				      * for the area of this quad.
				      */
    double measure () const;

				     /**
				      * Compute and return the number of
				      * children of this quad. Actually,
				      * this function only counts the number
				      * of active children, i.e. the number
				      * if quads which are not further
				      * refined. Thus, if all of the four
				      * children of a quad are further
				      * refined exactly once, the returned
				      * number will be 16, not 20.
				      *
				      * If the present cell is not refined,
				      * one is returned.
				      */
    unsigned int number_of_children () const;

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

  protected:
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

				     /**
				      * Declare some friends.
				      */
    template <int dim> friend class TriaRawIterator<dim,QuadAccessor<dim> >;
    template <>        friend class TriaRawIterator<2,CellAccessor<2> >;
};



/**
 *   Accessor to dereference the data of hexahedra. This accessor is used to
 *   point to hexs in #dim# space dimensions (only #dim>=3# seems reasonable
 *   to me). There is a derived class
 *   for hexs in three space dimension, in which case a hex is also a cell
 *   and thus has much more functionality than in other dimensions.
 *
 *   @author Wolfgang Bangerth, 1998
 */
template <int dim>
class HexAccessor :  public TriaAccessor<dim> {
  public:
				     /**
				      *  Constructor.
				      */
    HexAccessor (Triangulation<dim> *parent     = 0,
		 const int           level      = -1,
		 const int           index      = -1,
		 const void         *local_data = 0) :
		    TriaAccessor<dim> (parent, level, index, local_data) {};

				     /**
				      *  Copy the data of the given hex.
				      */
    void set (const Hexahedron &h) const;
    
				     /**
				      *  Return index of vertex #i=0 through 7# of
				      *  a hex. The convention regarding the
				      *  numbering of vertices is laid down
				      *  in the documentation of the
				      *  #Triangulation# class.
				      */ 
    int vertex_index (const unsigned int i) const;

    				     /**
				      *  Return a reference (not an iterator!)
				      *  to the #i#th vertex.
				      */
    Point<dim> & vertex (const unsigned int i) const;

				     /**
				      *  Return a pointer to the #i#th line
				      *  bounding this #Hex#.
				      */
    TriaIterator<dim,LineAccessor<dim> > line (const unsigned int i) const;

				     /**
				      * Return the line index of the #i#th
				      * line. The level is naturally
				      * the same as that of the hex.
				      */
    unsigned int line_index (const unsigned int i) const;
    
    				     /**
				      *  Return a pointer to the #i#th quad
				      *  bounding this #Hex#.
				      */
    TriaIterator<dim,QuadAccessor<dim> > quad (const unsigned int i) const;

				     /**
				      * Return the quad index of the #i#th
				      * side (a quad). The level is naturally
				      * the same as that of the hex.
				      */
    unsigned int quad_index (const unsigned int i) const;

				     /**
				      *  Test for the element being used
				      *  or not.
				      */
    bool used () const;

				     /**
				      *  Set the #used# flag. You should know
				      *  quite exactly what you are doing of you
				      *  touch this function. It is exclusively
				      *  for internal use in the library.
				      */
    void set_used_flag () const;

				     /**
				      *  Clear the #used# flag. You should know
				      *  quite exactly what you are doing of you
				      *  touch this function. It is exclusively
				      *  for internal use in the library.
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
				      *  set the user flag of this object
				      *  and of all its children and their
				      *  children, etc.
				      */
    void recursively_set_user_flag () const;

    				     /**
				      *  Clear the user flag of this object
				      *  and of all its children and their
				      *  children, etc.
				      */
    void recursively_clear_user_flag () const;

				     /**
				      * Set the user pointer of this hex
				      * to #p#.
				      */
    void set_user_pointer (void *p) const;

				     /**
				      * Reset the user pointer of this hex
				      * to a #NULL# pointer.
				      */
    void clear_user_pointer () const;

				     /**
				      * Access the value of the user pointer
				      * of this hex. It is in the
				      * responsibility of the user to make
				      * sure that the pointer points to
				      * something useful. You should use the
				      * new style cast operator to maintain
				      * a minimum of typesafety, e.g.
				      * #A *a=static_cast<A*>(cell->user_pointer());#.
				      */
    void * user_pointer () const;

				     /**
				      *  Return a pointer to the #i#th
				      *  child.
				      */
    TriaIterator<dim,HexAccessor<dim> > child (const unsigned int i) const;

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
				      *  Test whether the hex has children.
				      */
    bool has_children () const;

				     /**
				      * Return the number of times that this
				      * cell is refined. Note that not all
				      * its children are refined that often
				      * (which is why we prepend #max_#), 
				      * the returned number is rather the
				      * maximum number of refinement in
				      * any branch of children of this object.
				      */
    unsigned int max_refinement_depth () const;    
    
				     /**
				      * Return the boundary indicator of this
				      * hex. Since boundary data is only useful
				      * for structures with a dimension less
				      * than the dimension of a cell, this
				      * function issues an error if #dim<4#.
				      *
				      * If the return value is 255, then this
				      * line is in the interior of the domain.
				      */
    unsigned char boundary_indicator () const;

				     /**
				      * Set the boundary indicator of this hex.
				      * The same applies as for the
				      * #boundary_indicator()# function.
				      *
				      * You should be careful with this function
				      * and especially never try to set the
				      * boundary indicator to 255, unless
				      * you exactly know what you are doing,
				      * since this value is reserved for
				      * another purpose and algorithms may
				      * not work if boundary cells have have
				      * this boundary indicator or if interior
				      * cells have boundary indicators other
				      * than 255.
				      */
    void set_boundary_indicator (unsigned char) const;

				     /**
				      * Return whether this hex is at the
				      * boundary. This is checked via
				      * the boundary indicator field, which
				      * is always 255 if the hex is in the
				      * interior of the domain. Obviously,
				      * the use of this function is only
				      * possible for #dim>3#;
				      * however, for #dim==3#, a hex is
				      * a cell and the #CellAccessor# class
				      * offers another possibility to
				      * determine whether a cell is at the
				      * boundary or not.
				      */
    bool at_boundary () const;

				     /**
				      * Return the diameter of the hex.
				      *
				      * The diameter of a hex is computed to
				      * be the largest diagonal. You
				      * should absolutely be clear about the
				      * fact that this definitely is not the
				      * diameter of all hexahedra; however
				      * it may serve as an approximation and
				      * is exact in many cases, especially
				      * if the hexahedron is not too much
				      * distorted.
				      */
    double diameter () const;

    				     /**
				      * Return the center of the hex. The
				      * center of a hex is defined to be
				      * the average of the vertices,
				      * which is also the point where the
				      * trilinear mapping places the midpoint
				      * of the unit hex in real space.
				      * However, this may not be the point
				      * of the barycenter of the hex.
				      */
    Point<dim> center () const;

				     /**
				      * Return the barycenter of the hex.
				      */
    Point<dim> barycenter () const;

				     /**
				      * Return the volume of the hex. With
				      * volume, we mean the area included by
				      * the hexahedron if its faces are
				      * supposed to be derived by a trilinear
				      * mapping from the unit cell, only using
				      * the location of the vertices.
				      * Therefore, no information
				      * about the boundary is used; if
				      * you want other than trilinearly
				      * mapped unit hexahedra, ask the
				      * appropriate finite element class
				      * for the volume.
				      */
    double measure () const;

				     /**
				      * Compute and return the number of
				      * children of this hex. Actually,
				      * this function only counts the number
				      * of active children, i.e. the number
				      * if hexs which are not further
				      * refined. Thus, if all of the eight
				      * children of a hex are further
				      * refined exactly once, the returned
				      * number will be 64, not 80.
				      *
				      * If the present cell is not refined,
				      * one is returned.
				      */
    unsigned int number_of_children () const;

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
    void operator = (const HexAccessor &);

  protected:
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

				     /**
				      * Declare some friends.
				      */
    template <int dim> friend class TriaRawIterator<dim,HexAccessor<dim> >;
    template <>        friend class TriaRawIterator<3,CellAccessor<3> >;
};







/**
 * Intermediate, "typedef"-class, not for public use.
 */
template <int dim>
class TriaSubstructAccessor;



/**
 * Intermediate, "typedef"-class, not for public use.
 *
 * \subsection{Rationale}
 *
 * This class is only a wrapper class used to do kind of a typedef
 * with template parameters. This class and #TriaSubstructAccessor<2>#
 * wrap the following names:
 * \begin{verbatim}
 *   TriaSubstructAccessor<1> := LineAccessor<1>;
 *   TriaSubstructAccessor<2> := QuadAccessor<2>;
 *   TriaSubstructAccessor<3> := HexAccessor<3>;
 * \end{verbatim}
 * We do this rather complex (and needless, provided C++ the needed constructs!)
 * class hierarchy manipulation, since this way we can declare and implement
 * the \Ref{CellAccessor} dimension independent as an inheritance from
 * #TriaSubstructAccessor<dim>#. If we had not declared these
 * types, we would have to write two class declarations, one for
 * #CellAccessor<1>#, derived from #LineAccessor<1>#
 * and one for #CellAccessor<2>#, derived from
 * #QuadAccessor<2>#.
 */
template <>
class TriaSubstructAccessor<1> :  public LineAccessor<1> {
  public:
				     /**
				      * Propagate the AccessorData type
				      * into the present class.
				      */
    typedef typename LineAccessor<1>::AccessorData AccessorData;
    				     /**
				      * Constructor
				      */
    TriaSubstructAccessor (Triangulation<1> *tria,
			   const int         level,
			   const int         index,
			   const void       *local_data) :
		    LineAccessor<1> (tria,level,index,local_data) {};

    				     // do this here, since this way the
				     // CellAccessor has the possibility to know
				     // what a substruct_iterator is. Otherwise
				     // it would have to ask the DoFHandler<dim>
				     // which would need another big include
				     // file and maybe cyclic interdependence
    typedef void * substruct_iterator;
};



/**
 * Intermediate, "typedef"-class, not for public use.
 *
 * @see TriaSubstructAccessor<1>
 */
template <>
class TriaSubstructAccessor<2> : public QuadAccessor<2> {
  public:
				     /**
				      * Propagate the AccessorData type
				      * into the present class.
				      */
    typedef typename QuadAccessor<2>::AccessorData AccessorData;
    				     /**
				      * Constructor
				      */
    TriaSubstructAccessor (Triangulation<2> *tria,
			   const int         level,
			   const int         index,
			   const void       *local_data) :
		    QuadAccessor<2> (tria,level,index,local_data) {};

    				     // do this here, since this way the
				     // CellAccessor has the possibility to know
				     // what a substruct_iterator is. Otherwise
				     // it would have to ask the DoFHandler<dim>
				     // which would need another big include
				     // file and maybe cyclic interdependence
    typedef TriaIterator<2,LineAccessor<2> > substruct_iterator;
};




/**
 * Intermediate, "typedef"-class, not for public use.
 *
 * @see TriaSubstructAccessor<1>
 */
template <>
class TriaSubstructAccessor<3> : public HexAccessor<3> {
  public:
				     /**
				      * Propagate the AccessorData type
				      * into the present class.
				      */
    typedef typename HexAccessor<3>::AccessorData AccessorData;
    
    				     /**
				      * Constructor
				      */
    TriaSubstructAccessor (Triangulation<3> *tria,
			   const int         level,
			   const int         index,
			   const void       *local_data) :
		    HexAccessor<3> (tria,level,index,local_data) {};

    				     // do this here, since this way the
				     // CellAccessor has the possibility to know
				     // what a substruct_iterator is. Otherwise
				     // it would have to ask the DoFHandler<dim>
				     // which would need another big include
				     // file and maybe cyclic interdependence
    typedef TriaIterator<3,QuadAccessor<3> > substruct_iterator;
};






/**
 * This class allows access to a cell: a line in one dimension, a quad
 * in two dimension, etc.
 *
 * The following refers to any space dimension:
 * 
 * This class allows access to a {\bf cell}, which is a line in 1D and a quad in
 * 2D. Cells have more functionality than lines or quads by themselves, for
 * example they can be flagged for refinement, they have neighbors, they have
 * the possibility to check whether they are at the boundary etc. This class
 * offers access to all this data.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class CellAccessor :  public TriaSubstructAccessor<dim> {
  public:
				     /**
				      * Propagate the AccessorData type
				      * into the present class.
				      */
    typedef typename TriaSubstructAccessor<dim>::AccessorData AccessorData;
    
				     /**
				      *  Constructor.
				      */
    CellAccessor (Triangulation<dim> *parent     = 0,
		  const int           level      = -1,
		  const int           index      = -1,
		  const void         *local_data = 0) :
		    TriaSubstructAccessor<dim> (parent, level, index, local_data) {};

				     /**
				      *  Return a pointer to the #i#th
				      *  neighbor.
				      *  If the neighbor does not exist, an
				      *  invalid iterator is returned.
				      */
    TriaIterator<dim,CellAccessor<dim> > neighbor (const unsigned int i) const;

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
		       const TriaIterator<dim,CellAccessor<dim> > &pointer) const;

				     /**
				      *  Return whether the #i#th vertex or
				      *  face (depending on the dimension) is
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
				      *  Flag the cell pointed to fo
				      *  refinement. This function is only
				      *  allowed for active cells.
				      */
    void set_refine_flag () const;

				     /**
				      *  Clear the refinement flag.
				      */
    void clear_refine_flag () const;


    				     /**
				      *  Return whether the coarsen flag
				      *  is set or not.
				      */
    bool coarsen_flag_set () const;

				     /**
				      *  Flag the cell pointed to for
				      *  coarsening. This function is only
				      *  allowed for active cells.
				      */
    void set_coarsen_flag () const;

				     /**
				      *  Clear the coarsen flag.
				      */
    void clear_coarsen_flag () const;

				     /**
				      *  Return a pointer to the #i#th
				      *  child. Overloaded version which returns
				      *  a more reasonable iterator class.
				      */
    TriaIterator<dim,CellAccessor<dim> > child (const unsigned int i) const;

				     /**
				      * Return an iterator to the #i#th face
				      * of this cell.
				      *
				      * This function is not implemented in 1D,
				      * and maps to QuadAccessor::line in 2D.
				      */
    typename TriaSubstructAccessor<dim>::substruct_iterator
    face (const unsigned int i) const;
    
				     /**
				      * Return the material id of this cell.
				      */
    unsigned char material_id () const;

				     /**
				      * Set the material id of this cell.
				      */
    void set_material_id (const unsigned char) const;

				     /**
				      * Test whether the cell has children
				      * (this is the criterion for activity
				      * of a cell).
				      */
    bool active () const;
    
				     /**
				      *  Exception
				      */
    DeclException0 (ExcRefineCellNotActive);
				     /**
				      *  Exception
				      */
    DeclException0 (ExcCellFlaggedForRefinement);
				     /**
				      *  Exception
				      */
    DeclException0 (ExcCellFlaggedForCoarsening);    
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
    void operator = (const CellAccessor<dim> &);
};



// if in optimized mode: include more templates
#ifndef DEBUG
#  include "tria_accessor.templates.h"
#endif


/*----------------------------   tria_accessor.h     ---------------------------*/
/* end of #ifndef __tria_accessor_H */
#endif
/*----------------------------   tria_accessor.h     ---------------------------*/






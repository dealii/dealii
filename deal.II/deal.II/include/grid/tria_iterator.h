//----------------------------  tria_iterator.h  ---------------------------
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
//----------------------------  tria_iterator.h  ---------------------------
#ifndef __deal2__tria_iterator_h
#define __deal2__tria_iterator_h


/*----------------------------   tria-iterator.h     ---------------------------*/


#include <iterator>
#include <iostream>
#include <base/exceptions.h>
#include <base/point.h>
#include <grid/tria_iterator_base.h>

template <int dim> class Triangulation;

// note: in non-debug mode, i.e. with optimizations, the file
// tria_iterator.templates.h is included at the end of this file.
// this includes a lot of templates and thus makes compilation
// slower, but at the same time allows for more aggressive
// inlining and thus faster code.


/**
 *   This class implements an iterator, analogous to those of the standard
 *   template library (STL). It fulfills the requirements of a bidirectional iterator.
 *   See the C++ documentation for further details of iterator specification and
 *   usage. In addition to the STL
 *   iterators an iterator of this class provides a @p{->} operator, i.e. you can
 *   write statements like @p{i->set_refine_flag ();}.
 *   
 *   @em{Note:} Please read the documentation about the prefix and the
 *   postfix @p{++} operators in this and the derived classes!
 *   
 *   @sect3{Purpose}
 *
 *   @p{iterators} are used whenever a loop over all lines, quads, cells etc.
 *   is to be performed. These loops can then be coded like this:
 *   @begin{verbatim}
 *     cell_iterator i   = tria.begin();
 *     cell_iterator end = tria.end();
 *     for (; i!=end; ++i)
 *       if (cell->at_boundary())
 *         cell->set_refine_flag();
 *   @end{verbatim}
 *   Note the usage of @p{++i} instead of @p{i++} since this does not involve
 *   temporaries and copying. You should also really use a fixed value
 *   @p{end} rather than coding @p{for (; i!=tria.end(); ++i)}, since
 *   the creation and copying of these iterators is rather expensive
 *   compared to normal pointers.
 *   
 *   The objects pointed to by iterators are @p{TriangulationLevel<1>::LinesData},
 *   @p{TriangulationLevel<2>::LinesData}
 *   and @p{TriangulationLevel<2>::QuadsData}. To chose which of those, the
 *   template parameter @p{Pointee} is used.
 *
 *   Since the names as is are quite unhandy, the @p{Triangulation<>} class which
 *   uses these iterators declares typedef'd versions. See there for more
 *   information.
 *
 *   The objects pointed to are, as mentioned, @p{LinesData} etc. To be
 *   more exact, when dereferencing an iterator, you do not get a @p{LineData}
 *   object (or the like, but we will assume that you have a @p{line_iterator}
 *   in the following), but a @em{virtual} object (called @em{accessor}) which
 *   behaves as if it stored the data of a line. It does not contain any data
 *   itself, but provides functions to manipulate the data of the line it
 *   stands for.
 *
 *   Since the data of one line is splitted to
 *   several arrays (@p{lines}, @p{children} and @p{used}) for performance reasons
 *   rather than keeping all information in a @p{Line} struct, access through
 *   an accessor is usually much simpler than handling the exact data structure
 *   and also less error prone since the data structure itself can be changed
 *   in an arbitrary way while the only pieces of code which access these
 *   data structures are the accessors.
 *
 *   On the other hand, iterators are not much slower than operating directly
 *   on the data structures, since they perform the loops that you had
 *   to handcode yourself anyway. Most iterator and accessor functions are
 *   inlined. 
 *
 *   The main functionality of iterators, however, resides in the @p{++} and
 *   @p{--} operators. These move the iterator forward or backward just as if
 *   it were a pointer into an array. Here, this operation is not so easy,
 *   since it may include skipping some elements and the transition between
 *   the triangulation levels. This is completely hidden from the user, though
 *   you can still create an iterator pointing to an arbitrary element.
 *   Actually, the operation of moving iterators back and forth is not done in
 *   the iterator classes, but rather in the accessor classes. Since these are
 *   passed as template arguments, you can write your own versions here to add
 *   more functionality.
 *
 *   Furthermore, the iterators decribed here satisfy the requirement of
 *   input and bidirectional iterators as stated by the C++ standard and
 *   the STL documentation. It is therefore possible to use the functions
 *   from the @em{algorithm section} of the C++ standard, e.g. @p{count_if}
 *   (see the documentation for @ref{Triangulation} for an example) and
 *   several others.
 *   
 *
 *   @sect3{Differences between the classes in this inheritance tree}
 *
 *   @p{TriaRawIterator} objects point to lines, cells, etc in
 *   the lists whether they are used or not (in the vectors, also @em{dead}
 *   objects are stored, since deletion in vectors is expensive and we
 *   also do not want to destroy the ordering induced by the numbering
 *   in the vectors). Therefore not all raw iterators point to valid objects.
 *   
 *   There are two derived versions of this class: @ref{TriaIterator}
 *   objects, which only loop over used (valid) cells and
 *   @p{TriaActiveIterator} objects
 *   which only loop over active cells (not refined).
 *
 *   
 *   @sect3{Implementation}
 *
 *   In principle, the Iterator class does not have much functionality. It
 *   only becomes useful when assigned an @p{Accessor} (the second template
 *   parameter), which really does the access to data. An @p{Accessor} has to
 *   fulfil some requirements:
 *   @begin{itemize}
 *     @item It must have two members named @p{present_level} and @p{present_index}
 *       storing the address of the element in the triangulation presently
 *       pointed to. Furthermore, the three @p{Tria{Raw| |Active}Iterator} classes
 *       have to be friends to the accessor or these data members must be public.
 *     @item It must have a constructor which takes 1. a @p{Triangulation<dim>*},
 *       2. and 3. and integer, denoting the initial level and index.
 *     @item For the @p{TriaIterator} and the @p{TriaActiveIterator} class, it must
 *       have a member function @p{bool used()}, for the latter a member function
 *       @p{bool active()}.
 *     @item It must have void operators @p{++} and @p{--}.
 *     @item It must declare a local @p{typedef} @p{AccessorData} which states
 *       the data type the accessor expects to get passed as fourth constructor
 *       argument. By declaring a local data type, the respective iterator
 *       class may type-safely enforce that data type to be one of its own
 *       constructor argument types. If an accessor class does not need
 *       additional data, this type shall be @p{void}.
 *   @end{itemize}
 *   Then the iterator is able to do what it is supposed to. All of the necessary
 *   functions are implemented in the @p{Accessor} base class, but you may write
 *   your own version (non-virtual, since we use templates) to add functionality.
 *
 *   There is a standard implementation, using classes which are derived from
 *   @ref{TriaAccessor}. These classes point to @p{Line}s, @p{Quad}s and the like.
 *   For advanced use of the iterator classes, derive classes from
 *   @p{{Line|Quad|Cell}Accessor} which also dereference data structures in other
 *   objects, e.g. in a finite element context. An iterator with such an accessor
 *   then simultaneously points to (for example) a cell in the triangulation and
 *   the data stored on it in the finite element class.
 *
 *   Derived accessor classes may need additional data (e.g. the @p{DoFAccessor}
 *   needs a pointer to the @p{DoFHandler} to work on). This data can be
 *   set upon construction through the last argument of the constructors.
 *   The data type of this additional data is given by the local data type
 *   @p{AccessorData} explained above. The iterator classes take a pointer to
 *   an object of that data type; by default, this parameter equals the
 *   @p{NULL} pointer.
 *   
 *   
 *   @sect3{Warning}
 *
 *   It seems impossible to preserve @p{const}ness of a triangulation through
 *   iterator usage. Thus, if you declare pointers to a @p{const} triangulation
 *   object, you should be well aware that you might involuntarily alter the
 *   data stored in the triangulation.
 *
 *
 *   @sect3{Internals}
 *   
 *   There is a representation of past-the-end-pointers, denoted by special
 *   values of the member variables @p{present_level} and @p{present_index}:
 *   If @p{present_level>=0} and @p{present_index>=0}, then the object is valid
 *   (there is no check whether the triangulation really has that
 *   many levels or that many cells on the present level when we investigate
 *   the state of an iterator; however, in many places where an iterator is
 *   dereferenced we make this check);
 *   if @p{present_level==-1} and @p{present_index==-1}, then the iterator points
 *   past the end; in all other cases, the iterator is considered invalid.
 *   You can check this by calling the @p{state()} function.
 *
 *   An iterator is also invalid, if the pointer pointing to the @p{Triangulation}
 *   object is invalid or zero.
 *
 *   Finally, an iterator is invalid, if the element pointed to by
 *   @p{present_level} and @p{present_index} is not used, i.e. if the @p{used}
 *   flag is set to false.
 *
 *   The last two checks are not made in @p{state()} since both cases should only
 *   occur upon unitialized construction through @p{memcpy} and the like (the
 *   parent triangulation can only be set upon construction). If
 *   an iterator is constructed empty through the empty constructor,
 *   @p{present_level==-2} and @p{present_index==-2}. Thus, the iterator is
 *   invalid anyway, regardless of the state of the triangulation pointer
 *   and the state of the element pointed to.
 *
 *   Past-the-end iterators may also be used to compare an iterator with the
 *   @em{before-the-start} value, when running backwards. There is no
 *   distiction between the iterators pointing past the two ends of a vector.
 *   
 *   By defining only one value to be past-the-end and making all other values
 *   invalid porvides a second track of security: if we should have forgotten
 *   a check in the library when an iterator is incremented or decremented,
 *   we automatically convert the iterator from the allowed state "past-the-end"
 *   to the disallowed state "invalid" which increases the chance that somehwen
 *   earlier than for past-the-end iterators an exception is raised.
 *
 *   @see Triangulation
 *   @see TriaDimensionInfo
 *   @author Wolfgang Bangerth, 1998
 */
template <int dim, typename Accessor>
class TriaRawIterator : public bidirectional_iterator<Accessor,int> {
  public:
				     /**
				      * Declare the type of the Accessor for
				      * use in the outside world. This way other
				      * functions can use the Accessor's type
				      * without knowledge of how the exact
				      * implementation actually is.
				      */
    typedef Accessor AccessorType;
    
				     /**
				      *  Empty constructor. Such an object
				      *  is not usable!
				      */
    TriaRawIterator ();

				     /**
				      *  Copy constructor.
				      */
    TriaRawIterator (const TriaRawIterator &);

				     /**
				      * Construct an iterator from the given
				      * accessor; the given accessor needs not
				      * be of the same type as the accessor of
				      * this class is, but it needs to be
				      * convertible.
				      *
				      * Through this constructor, it is also
				      * possible to construct object for
				      * derived iterators:
				      * @begin{verbatim}
				      * DoFCellAccessor dof_accessor;
				      * Triangulation::active_cell_iterator cell
				      *   = accessor;
				      * @end{verbatim}
				      */
    TriaRawIterator (const Accessor &a);
    
				     /**
				      *  Proper constructor, initialized
				      *  with the triangulation, the
				      *  level and index of the object
				      *  pointed to. The last parameter is
				      *  of a type declared by the accessor
				      *  class.
				      */
    TriaRawIterator (Triangulation<dim> *parent,
		     const int           level,
		     const int           index,
		     const typename Accessor::AccessorData *local_data = 0);

				     /**
				      * This is a conversion operator
				      * (constructor) which takes another
				      * iterator type and copies the data;
				      * this conversion works, if there is
				      * a conversion path from the
				      * @p{OtherAccessor} class to the @p{Accessor}
				      * class of this object. One such path
				      * would be derived class to base class,
				      * which for example may be used to get
				      * a @p{Triangulation::raw_cell_iterator} from
				      * a @p{DoFHandler::raw_cell_iterator}, since
				      * the @p{DoFAccessor} class is derived from
				      * the @p{TriaAccessor} class.
				      *
				      * Since @p{TriaIterator} and
				      * @p{TriaActiveIterator} are derived classes
				      * of this class, this constructor also
				      * serves to convert these iterators with
				      * other accessor classes.
				      */
    template <typename OtherAccessor>
    TriaRawIterator (const TriaRawIterator<dim,OtherAccessor> &i);
    
				     /**
				      *  @name Dereferencing
				      */
				     /*@{*/
				     /**
				      *  Dereferencing operator, returns a
				      *  reference to an accessor.
				      *  Usage is thus like @p{(*i).index ();}
				      *
				      *  This function has to be specialized
				      *  explicitely for the different
				      *  @p{Pointee}s, to allow an
				      *  @p{iterator<1,TriangulationLevel<1>::LinesData>}
				      *  to point to @p{tria->lines.lines[index]} while
				      *  for one dimension higher it has
				      *  to point to @p{tria->quads.quads[index]}.
				      *
				      *  You must not dereference invalid or
				      *  past the end iterators.
				      */
    const Accessor & operator * () const;
    
				     /**
				      *  Dereferencing operator, non-@p{const}
				      *  version.
				      */
    Accessor & operator * ();
        
				     /**
				      *  Dereferencing operator, returns a
				      *  reference of the cell pointed to.
				      *  Usage is thus like @p{i->index ();}
				      *
				      *  There is a @p{const} and a non-@p{const}
				      *  version.
				      */
    const Accessor * operator -> () const;
        
				     /**
				      *  Dereferencing operator, non-@p{const}
				      *  version.
				      */
    Accessor * operator -> ();
				     /*@}*/
    
				     /**
				      *  Assignment operator.
				      */
    TriaRawIterator & operator = (const TriaRawIterator &);
    
				     /**
				      *  Compare for equality.
				      */
    bool operator == (const TriaRawIterator &) const;
    
				     /**
				      *  Compare for inequality.
				      */
    bool operator != (const TriaRawIterator &) const;

				     /**
				      * Offer a weak ordering of iterators,
				      * which is needed to make @p{map}s with
				      * iterators being keys. An iterator
				      * pointing to an element @p{a} is
				      * less than another iterator pointing
				      * to an element @p{b} if
				      * level(a)<level(b) or
				      * (level(a)==level(b) and index(a)<index(b)).
				      *
				      * Comparing iterators of which one or
				      * both are invalid is an error. The
				      * past-the-end iterator is always
				      * ordered last. Two past-the-end
				      * iterators rank the same, thus false
				      * is returned in that case.
				      */
    bool operator < (const TriaRawIterator &) const;
    
				     /**@name Advancement of iterators*/
				     /*@{*/
				     /**
				      *  Prefix @p{++} operator: @p{++i}. This
				      *  operator advances the iterator to
				      *  the next element and returns
				      *  a reference to @p{*this}.
				      *
				      *  The next element is next on this
				      *  level if there are more. If the
				      *  present element is the last on
				      *  this level, the first on the
				      *  next level is accessed.
				      */
    TriaRawIterator & operator ++ ();
    
				     /**
				      *  Postfix @p{++} operator: @p{i++}. This
				      *  operator advances the iterator to
				      *  the next element, but
				      *  returns an iterator to the element
				      *  priviously pointed to. Since this
				      *  involves a temporary and a copy
				      *  operation and since an
				      *  @p{iterator} is quite a large
				      *  object for a pointer, use the
				      *  prefix operator @p{++i} whenever
				      *  possible, especially in the head
				      *  of for loops
				      *  (@p{for (; i!=end; ++i)}) since there
				      *  you normally never need the
				      *  returned value.
				      */
    TriaRawIterator operator ++ (int);

    				     /**
				      *  Prefix @p{--} operator: @p{--i}. This
				      *  operator advances the iterator to
				      *  the previous element and returns
				      *  a reference to @p{*this}.
				      *
				      *  The previous element is previous on
				      *  this level if @p{index>0}. If the
				      *  present element is the first on
				      *  this level, the last on the
				      *  previous level is accessed.
				      */
    TriaRawIterator & operator -- ();
    
				     /**
				      *  Postfix @p{--} operator: @p{i--}. This
				      *  operator advances the iterator to
				      *  the previous element, but
				      *  returns an iterator to the element
				      *  priviously pointed to. Since this
				      *  involves a temporary and a copy
				      *  operation and since an
				      *  @p{iterator} is quite a large
				      *  object for a pointer, use the
				      *  prefix operator @p{--i} whenever
				      *  possible, especially in the head
				      *  of for loops
				      *  (@p{for (; i!=end; --i)}) since there
				      *  you normally never need the
				      *  returned value.
				      */
    TriaRawIterator operator -- (int);
				     /*@}*/
    
				     /**
				      *  Return the state of the iterator.
				      */
    IteratorState state () const;

				     /**
				      * Print the iterator to @p{out}. The
				      * format is like @p{level.index}.
				      */
    void print (ostream &out) const;


				     /**@name Exceptions*/
				     /*@{*/
				     /**
				      *  Exception
				      */
    DeclException0 (ExcDereferenceInvalidObject);
				     /**
				      *  Exception
				      */
    DeclException0 (ExcAdvanceInvalidObject);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidComparison);
    
				     /*@}*/
  protected:
				     /**
				      * Object holding the real data.
				      */
    Accessor accessor;


				     /**
				      * Make all other iterator class templates
				      * friends of this class. This is
				      * necessary for the implementation of
				      * conversion constructors.
				      *
				      * In fact, we would not need them to
				      * be friends if they were for different
				      * dimensions, but the compiler dislikes
				      * giving a fixed dimension and variable
				      * accessor since then it says that would
				      * be a artial specialization.
				      */
    template <int anydim, typename SomeAccessor> friend class TriaRawIterator;
    template <int anydim, typename SomeAccessor> friend class TriaIterator;
    template <int anydim, typename SomeAccessor> friend class TriaActiveIterator;
};


/**
 *   This specialization of @ref{TriaRawIterator} provides access only to the
 *   @em{used} lines, quads, cells, etc.
 */
template <int dim, typename Accessor>
class TriaIterator : public TriaRawIterator<dim,Accessor> {
  public:
				     /**
				      *  Empty constructor. Such an object
				      *  is not usable!
				      */
    TriaIterator ();

				     /**
				      *  Copy constructor.
				      */
    TriaIterator (const TriaIterator<dim,Accessor> &);

    				     /**
				      *  Cross copy constructor from
				      *  iterators pointing also to
				      *  non-active objects.
				      *
				      *  If the object pointed to is not
				      *  past-the-end and is not
				      *  used, the debug version raises
				      *  an error!
				      */
    TriaIterator (const TriaRawIterator<dim,Accessor> &);

				     /**
				      *  Proper constructor, initialized
				      *  with the triangulation, the
				      *  level and index of the object
				      *  pointed to. The last parameter is
				      *  of a type declared by the accessor
				      *  class.
				      *
				      *  If the object pointed to is not
				      *  past-the-end and is not
				      *  used, the debug version raises
				      *  an error!
				      */
    TriaIterator (Triangulation<dim> *parent,
		  const int           level,
		  const int           index,
		  const typename Accessor::AccessorData *local_data = 0);
    
				     /**
				      * This is a conversion operator
				      * (constructor) which takes another
				      * iterator type and copies the data;
				      * this conversion works, if there is
				      * a conversion path from the
				      * @p{OtherAccessor} class to the @p{Accessor}
				      * class of this object. One such path
				      * would be derived class to base class,
				      * which for example may be used to get
				      * a @p{Triangulation::cell_iterator} from
				      * a @p{DoFHandler::cell_iterator}, since
				      * the @p{DoFAccessor} class is derived from
				      * the @p{TriaAccessor} class.
				      *
				      * Since @p{TriaActiveIterator} is a derived
				      * class of this class, this constructor
				      * also serves to convert these iterators 
				      * with other accessor classes.
				      */
    template <typename OtherAccessor>
    TriaIterator (const TriaIterator<dim,OtherAccessor> &i);

				     /**
				      * Similar conversion operator to the above
				      * one, but does a check whether the
				      * iterator points to a used element,
				      * which is necessary for raw iterators.
				      */
    template <typename OtherAccessor>
    TriaIterator (const TriaRawIterator<dim,OtherAccessor> &i);

    				     /**
				      *  Assignment operator.
				      */
    TriaIterator<dim,Accessor> &
    operator = (const TriaIterator<dim,Accessor> &);

    				     /**
				      *  Cross assignment operator. This
				      *  assignment is only valid if the
				      *  given iterator points to a used
				      *  element.
				      */
    TriaIterator<dim,Accessor> &
    operator = (const TriaRawIterator<dim,Accessor> &);

				     /**@name Advancement of iterators*/
				     /*@{*/
				     /**
				      *  Prefix @p{++} operator: @p{++i}. This
				      *  operator advances the iterator to
				      *  the next used element and returns
				      *  a reference to @p{*this}.
				      */
    TriaIterator<dim,Accessor> & operator ++ ();

				     /**
				      *  Postfix @p{++} operator: @p{i++}. This
				      *  operator advances the iterator to
				      *  the next used element, but
				      *  returns an iterator to the element
				      *  previously pointed to. Since this
				      *  involves a temporary and a copy
				      *  operation and since an
				      *  @p{active_iterator} is quite a large
				      *  object for a pointer, use the
				      *  prefix operator @p{++i} whenever
				      *  possible, especially in the head
				      *  of for loops
				      *  (@p{for (; i!=end; ++i)}) since there
				      *  you normally never need the
				      *  returned value.
				      */
    TriaIterator<dim,Accessor> operator ++ (int);

    				     /**
				      *  Prefix @p{--} operator: @p{--i}. This
				      *  operator advances the iterator to
				      *  the previous used element and returns
				      *  a reference to @p{*this}.
				      */
    TriaIterator<dim,Accessor> & operator -- ();

				     /**
				      *  Postfix @p{--} operator: @p{i--}.
				      */
    TriaIterator<dim,Accessor> operator -- (int);
				     /*@}*/
    
				     /**
				      *  Exception
				      */
    DeclException0 (ExcAssignmentOfUnusedObject);
};


/**
 *   This specialization of @ref{TriaIterator} provides access only to the
 *   @em{active} lines, quads, cells, etc. An active cell is a cell which is not
 *   refined and thus a cell on which calculations on the finest level are done.
 */
template <int dim, typename Accessor>
class TriaActiveIterator : public TriaIterator<dim,Accessor> {
  public:
				     /**
				      *  Empty constructor. Such an object
				      *  is not usable!
				      */
    TriaActiveIterator ();

				     /**
				      *  Copy constructor.
				      */
    TriaActiveIterator (const TriaActiveIterator<dim,Accessor> &);

    				     /**
				      *  Cross copy constructor from
				      *  iterators pointing also to
				      *  non-active objects.
				      *
				      *  If the object pointed to is not
				      *  past-the-end and is not
				      *  active, the debug version raises
				      *  an error!
				      */
    TriaActiveIterator (const TriaRawIterator<dim,Accessor> &);

				     /**
				      *  Cross copy constructor from
				      *  iterators pointing also to
				      *  non-active objects.
				      *
				      *  If the object pointed to is not
				      *  past-the-end and is not
				      *  active, the debug version raises
				      *  an error!
				      */
    TriaActiveIterator (const TriaIterator<dim,Accessor> &);

				     /**
				      *  Proper constructor, initialized
				      *  with the triangulation, the
				      *  level and index of the object
				      *  pointed to. The last parameter is
				      *  of a type declared by the accessor
				      *  class.
				      *
				      *  If the object pointed to is not
				      *  past-the-end and is not
				      *  active, the debug version raises
				      *  an error!
				      */
    TriaActiveIterator (Triangulation<dim> *parent,
			const int           level,
			const int           index,
			const typename Accessor::AccessorData *local_data = 0);

				     /**
				      * This is a conversion operator
				      * (constructor) which takes another
				      * iterator type and copies the data;
				      * this conversion works, if there is
				      * a conversion path from the
				      * @p{OtherAccessor} class to the @p{Accessor}
				      * class of this object. One such path
				      * would be derived class to base class,
				      * which for example may be used to get
				      * a @p{Triangulation::active_cell_iterator} from
				      * a @p{DoFHandler::active_cell_iterator}, since
				      * the @p{DoFAccessor} class is derived from
				      * the @p{TriaAccessor} class.
				      */
    template <typename OtherAccessor>
    TriaActiveIterator (const TriaActiveIterator<dim,OtherAccessor> &i);

				     /**
				      * Similar conversion operator to the above
				      * one, but does a check whether the
				      * iterator points to a used element,
				      * and is active, which is necessary for
				      * raw iterators. Since usual iterators
				      * are also raw iterators, this constructor
				      * works also for parameters of type
				      * @p{TriaIterator<dim,OtherAccessor>}.
				      */
    template <typename OtherAccessor>
    TriaActiveIterator (const TriaRawIterator<dim,OtherAccessor> &i);

    				     /**
				      *  Assignment operator.
				      */
    TriaActiveIterator<dim,Accessor> &
    operator = (const TriaActiveIterator<dim,Accessor> &);

    				     /**
				      *  Cross assignment operator. This
				      *  assignment is only valid if the
				      *  given iterator points to an active
				      *  element or past the end.
				      */
    TriaActiveIterator<dim,Accessor> &
    operator = (const TriaRawIterator<dim,Accessor> &);

				     /**
				      *  Cross assignment operator. This
				      *  assignment is only valid if the
				      *  given iterator points to an active
				      *  element or past the end.
				      */
    TriaActiveIterator<dim,Accessor> &
    operator = (const TriaIterator<dim,Accessor> &);

				     /**
				      *  Prefix @p{++} operator: @p{++i}. This
				      *  operator advances the iterator to
				      *  the next active element and returns
				      *  a reference to @p{*this}.
				      */
    TriaActiveIterator<dim,Accessor> & operator ++ ();

				     /**@name Advancement of iterators*/
				     /*@{*/
				     /**
				      *  Postfix @p{++} operator: @p{i++}. This
				      *  operator advances the iterator to
				      *  the next active element, but
				      *  returns an iterator to the element
				      *  previously pointed to. Since this
				      *  involves a temporary and a copy
				      *  operation and since an
				      *  @p{active_iterator} is quite a large
				      *  object for a pointer, use the
				      *  prefix operator @p{++i} whenever
				      *  possible, especially in the head
				      *  of for loops
				      *  (@p{for (; i!=end; ++i)}) since there
				      *  you normally never need the
				      *  returned value.
				      */
    TriaActiveIterator<dim,Accessor> operator ++ (int);

    				     /**
				      *  Prefix @p{--} operator: @p{--i}. This
				      *  operator advances the iterator to
				      *  the previous active element and
				      *  returns a reference to @p{*this}.
				      */
    TriaActiveIterator<dim,Accessor> & operator -- ();

				     /**
				      *  Postfix @p{--} operator: @p{i--}.
				      */
    TriaActiveIterator<dim,Accessor> operator -- (int);
				     /*@}*/
    
				     /**
				      *  Exception
				      */
    DeclException0 (ExcAssignmentOfInactiveObject);
};


/*----------------------- Inline functions -------------------*/


template <int dim, typename Accessor>
inline
TriaRawIterator<dim,Accessor>::TriaRawIterator (const Accessor &a) :
		accessor (a)
{};



template <int dim, typename Accessor>
template <typename OtherAccessor>
inline
TriaRawIterator<dim,Accessor>::TriaRawIterator (const TriaRawIterator<dim,OtherAccessor> &i)
		: accessor (i.accessor)
{};



template <int dim, typename Accessor>
inline
const Accessor &
TriaRawIterator<dim,Accessor>::operator * () const
{
  Assert (state() == valid, ExcDereferenceInvalidObject());
  return accessor;
};



template <int dim, typename Accessor>
inline
Accessor &
TriaRawIterator<dim,Accessor>::operator * ()
{
  Assert (state() == valid, ExcDereferenceInvalidObject());
  return accessor;
};



template <int dim, typename Accessor>
inline
const Accessor *
TriaRawIterator<dim,Accessor>::operator -> () const
{
  return &(this->operator* ());
};



template <int dim, typename Accessor>
inline
Accessor *
TriaRawIterator<dim,Accessor>::operator -> ()
{
  return &(this->operator* ());
};



template <int dim, typename Accessor>
inline
IteratorState
TriaRawIterator<dim,Accessor>::state () const
{
  return accessor.state ();
};



template <int dim, typename Accessor>
inline
bool
TriaRawIterator<dim,Accessor>::operator < (const TriaRawIterator &i) const
{
  Assert (state() != invalid, ExcDereferenceInvalidObject());
  Assert (i.state() != invalid, ExcDereferenceInvalidObject());
  Assert (&accessor.get_triangulation() == &i.accessor.get_triangulation(),
	  ExcInvalidComparison());
  
  return ((((accessor.level() < i.accessor.level()) ||
	    ((accessor.level() == i.accessor.level()) &&
	     (accessor.index() < i.accessor.index()))        ) &&
	   (state()==valid)                                  &&
	   (i.state()==valid)                                  ) ||
	  ((state()==valid) && (i.state()==past_the_end)));
};


template <int dim, typename Accessor>
inline
TriaRawIterator<dim,Accessor> &
TriaRawIterator<dim,Accessor>::operator ++ ()
{
  Assert (state() == valid, ExcAdvanceInvalidObject());

  ++accessor;
  return *this;
};



template <int dim, typename Accessor>
inline
TriaRawIterator<dim,Accessor> &
TriaRawIterator<dim,Accessor>::operator -- ()
{
  Assert (state() == valid, ExcAdvanceInvalidObject());

  --accessor;
  return *this;
};



template <int dim, typename Accessor>
inline
void
TriaRawIterator<dim,Accessor>::print (ostream &out) const
{
  out << accessor.level() << "." << accessor.index();
};



template <int dim, typename Accessor>
template <typename OtherAccessor>
inline
TriaIterator<dim,Accessor>::TriaIterator (const TriaIterator<dim,OtherAccessor> &i)
		: TriaRawIterator<dim,Accessor> (static_cast<TriaIterator<dim,OtherAccessor> >(i))
{};



template <int dim, typename Accessor>
template <typename OtherAccessor>
inline
TriaIterator<dim,Accessor>::TriaIterator (const TriaRawIterator<dim,OtherAccessor> &i)
		: TriaRawIterator<dim,Accessor> (static_cast<TriaRawIterator<dim,OtherAccessor> >(i))
{
#ifdef DEBUG
				   // do this like this, because:
				   // if we write
				   // "Assert (past_the_end || used)"
				   // used() is called anyway, even if
				   // state==past_the_end, and will then
				   // throw the exception!
  if (state() != past_the_end)
    Assert (accessor.used(),
	    ExcAssignmentOfUnusedObject());
#endif
};



template <int dim, typename Accessor>
template <typename OtherAccessor>
inline
TriaActiveIterator<dim,Accessor>::TriaActiveIterator (const TriaActiveIterator<dim,OtherAccessor> &i)
		: TriaIterator<dim,Accessor> (static_cast<TriaIterator<dim,OtherAccessor> >(i))
{};



template <int dim, typename Accessor>
template <typename OtherAccessor>
inline
TriaActiveIterator<dim,Accessor>::TriaActiveIterator (const TriaRawIterator<dim,OtherAccessor> &i)
		: TriaIterator<dim,Accessor> (i)
{
#ifdef DEBUG
				   // do this like this, because:
				   // if we write
				   // "Assert (past_the_end || used)"
				   // has_children() is called anyway, even if
				   // state==past_the_end, and will then
				   // throw the exception!
  if (state() != past_the_end) 
    Assert (accessor.has_children()==false,
	    ExcAssignmentOfInactiveObject());
#endif
};



template <int dim, typename Accessor>
inline
ostream & operator << (ostream &out, const TriaRawIterator<dim,Accessor> &i)
{
  i.print(out);
  return out;
};



template <int dim, typename Accessor>
inline
ostream & operator << (ostream &out, const TriaIterator<dim,Accessor> &i)
{
  i.print(out);
  return out;
};



template <int dim, typename Accessor>
inline
ostream & operator << (ostream &out, const TriaActiveIterator<dim,Accessor> &i)
{
  i.print(out);
  return out;
};



// if in optimized mode: include more templates
#ifndef DEBUG
#  include "tria_iterator.templates.h"
#endif


/*----------------------------   tria-iterator.h     ---------------------------*/

#endif
/*----------------------------   tria-iterator.h     ---------------------------*/



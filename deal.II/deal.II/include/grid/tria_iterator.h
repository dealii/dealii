/*----------------------------   tria-iterator.h     ---------------------------*/
/*      $Id$                 */
#ifndef __tria_iterator_H
#define __tria_iterator_H
/*----------------------------   tria-iterator.h     ---------------------------*/



#include <iterator>
#include <iostream>
#include <base/exceptions.h>
#include <grid/point.h>
#include <grid/tria_accessor.h>





// forward declaration needed
template <int dim> class Triangulation;








/**
 *   This class implements an iterator, analogous to those of the standard
 *   template library (STL). It fulfills the requirements of a bidirectional iterator.
 *   See the C++ documentation for further details of iterator specification and
 *   usage. In addition to the STL
 *   iterators an iterator of this class provides a #-># operator, i.e. you can
 *   write statements like #i->set_refine_flag ();#.
 *   
 *   {\bf Note:} Please read the documentation about the prefix and the
 *   postfix #++# operators in this and the derived classes!
 *   
 *   \subsection{Purpose}
 *
 *   #iterators# are used whenever a loop over all lines, quads, cells etc.
 *   is to be performed. These loops can then be coded like this:
 *   \begin{verbatim}
 *     cell_iterator i   = tria.begin();
 *     cell_iterator end = tria.end();
 *     for (; i!=end; ++i)
 *       if (cell->at_boundary())
 *         cell->set_refine_flag();
 *   \end{verbatim}
 *   Note the usage of #++i# instead of #i++# since this does not involve
 *   temporaries and copying. You should also really use a fixed value
 *   #end# rather than coding #for (; i!=tria.end(); ++i)#, since
 *   the creation and copying of these iterators is rather expensive
 *   compared to normal pointers.
 *   
 *   The objects pointed to by iterators are #TriangulationLevel<1>::LinesData#,
 *   #TriangulationLevel<2>::LinesData#
 *   and #TriangulationLevel<2>::QuadsData#. To chose which of those, the
 *   template parameter #Pointee# is used.
 *
 *   Since the names as is are quite unhandy, the #Triangulation<># class which
 *   uses these iterators declares typedef'd versions. See there for more
 *   information.
 *
 *   The objects pointed to are, as mentioned, #LinesData# etc. To be
 *   more exact, when dereferencing an iterator, you do not get a #LineData#
 *   object (or the like, but we will assume that you have a #line_iterator#
 *   in the following), but a {\it virtual} object (called {\it accessor}) which
 *   behaves as if it stored the data of a line. It does not contain any data
 *   itself, but provides functions to manipulate the data of the line it
 *   stands for.
 *
 *   Since the data of one line is splitted to
 *   several arrays (#lines#, #children# and #used#) for performance reasons
 *   rather than keeping all information in a #Line# struct, access through
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
 *   The main functionality of iterators, however, resides in the #++# and
 *   #--# operators. These move the iterator forward or backward just as if
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
 *   from the {\it algorithm section} of the C++ standard, e.g. #count_if#
 *   (see the documentation for \Ref{Triangulation} for an example) and
 *   several others. Unfortunately, with some of them (e.g. #distance#),
 *   g++2.7 has some problems and we will have to wait for g++2.8.
 *   
 *
 *   \subsection{Differences between the classes in this inheritance tree}
 *
 *   #TriaRawIterator# objects point to lines, cells, etc in
 *   the lists whether they are used or not (in the vectors, also {\it dead}
 *   objects are stored, since deletion in vectors is expensive and we
 *   also do not want to destroy the ordering induced by the numbering
 *   in the vectors). Therefore not all raw iterators point to valid objects.
 *   
 *   There are two derived versions of this class: \Ref{TriaIterator}
 *   objects, which only loop over used (valid) cells and
 *   #TriaActiveIterator# objects
 *   which only loop over active cells (not refined).
 *
 *   
 *   \subsection{Implementation}
 *
 *   In principle, the Iterator class does not have much functionality. It
 *   only becomes useful when assigned an #Accessor# (the second template
 *   parameter), which really does the access to data. An #Accessor# has to
 *   fulfil some requirements:
 *   \begin{itemize}
 *     \item It must have two members named #present_level# and #present_index#
 *       storing the address of the element in the triangulation presently
 *       pointed to. Furthermore, the three #Tria{Raw| |Active}Iterator# classes
 *       have to be friends to the accessor or these data members must be public.
 *     \item It must have a constructor which takes 1. a #Triangulation<dim>*#,
 *       2. and 3. and integer, denoting the initial level and index.
 *     \item For the #TriaIterator# and the #TriaActiveIterator# class, it must
 *       have a member function #bool used()#, for the latter a member function
 *       #bool active()#.
 *     \item It should not modify the #present_level# and #present_index# fields,
 *       since this is what the iterator classes do, but it should use them to
 *       dereference the data it points to.
 *     \item It must have void operators #++# and #--#.	
 *   \end{itemize}
 *   Then the iterator is able to do what it is supposed to. All of the necessary
 *   functions are implemented in the #Accessor# base class, but you may write
 *   your own version (non-virtual, since we use templates) to add functionality.
 *
 *   There is a standard implementation, using classes which are derived from
 *   \Ref{TriaAccessor}. These classes point to #Line#s, #Quad#s and the like.
 *   For advanced use of the iterator classes, derive classes from
 *   #{Line|Quad|Cell}Accessor# which also dereference data structures in other
 *   objects, e.g. in a finite element context. An iterator with such an accessor
 *   then simultaneously points to (for example) a cell in the triangulation and
 *   the data stored on it in the finite element class.
 *
 *   Derived accessor classes may need additional data (e.g. the #DoFAccessor#
 *   needs a pointer to the #DoFHandler# to work on). This data can be
 *   set upon construction through the last argument of the constructors.
 *   Ideally, its type is a local type to the accessor and must have the name
 *   #Accessor::LocalData#. In the standard implementation, this type is
 *   declared to be a void pointer. The iterator constructors take their
 *   last argument carrying the additional data by default as zero, so unless
 *   #Accessor::LocalData# is a number or a pointer you may not construct
 *   such an iterator without giving the last argument. If you want to use
 *   the additional data, you also have to overload the #TriaAccessor::copy_data#
 *   function.
 *
 *   Unfortunately, the skeched way does not work, since gcc is not able to
 *   recognize the type defined local to the template argument (it does not
 *   suport the #typename# keyword at present), so we can only pass a voie
 *   pointer. You may, however, convert this to any type, normally to another
 *   pointer or to a pointer to a structure pointing to the data to be passed.
 *   The mechanism may be changed if the mentioned features appear in gcc.
 *
 *   Another possibility would be to have a function, say #set_local_data(...)#
 *   in the accessor classes which need additional data. You could then create
 *   an iterator like this:
 *   \begin{verbatim}
 *     TriaIterator<1,MyAccesor> i;
 *     i->set_local_data (1,2,3);
 *   \end{verbatim}
 *   But this will not always work: if the iterator #i# is not a valid one, then
 *   the library will forbid you to dereference it (which normally is a good
 *   idea), thus resulting in an error when you dereference it in the second
 *   line.
 *   
 *   
 *   \subsection{Warning}
 *
 *   It seems impossible to preserve #const#ness of a triangulation through
 *   iterator usage. Thus, if you declare pointers to a #const# triangulation
 *   object, you should be well aware that you might involuntarily alter the
 *   data stored in the triangulation.
 *
 *
 *   \subsection{Internals}
 *   
 *   There is a representation of past-the-end-pointers, denoted by special
 *   values of the member variables #present_level# and #present_index#:
 *   If #present_level>=0# and #present_index>=0#, then the object is valid
 *   (there is no check whether the triangulation really has that
 *   many levels or that many cells on the present level when we investigate
 *   the state of an iterator; however, in many places where an iterator is
 *   dereferenced we make this check);
 *   if #present_level==-1# and #present_index==-1#, then the iterator points
 *   past the end; in all other cases, the iterator is considered invalid.
 *   You can check this by calling the #state()# function.
 *
 *   An iterator is also invalid, if the pointer pointing to the #Triangulation#
 *   object is invalid or zero.
 *
 *   Finally, an iterator is invalid, if the element pointed to by
 *   #present_level# and #present_index# is not used, i.e. if the #used#
 *   flag is set to false.
 *
 *   The last two checks are not made in #state()# since both cases should only
 *   occur upon unitialized construction through #memcpy# and the like (the
 *   parent triangulation can only be set upon construction). If
 *   an iterator is constructed empty through the empty constructor,
 *   #present_level==-2# and #present_index==-2#. Thus, the iterator is
 *   invalid anyway, regardless of the state of the triangulation pointer
 *   and the state of the element pointed to.
 *
 *   Past-the-end iterators may also be used to compare an iterator with the
 *   {\it before-the-start} value, when running backwards. There is no
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
class TriaRawIterator : public bidirectional_iterator<Accessor,int>{
  public:
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
				      *  Proper constructor, initialized
				      *  with the triangulation, the
				      *  level and index of the object
				      *  pointed to.
				      */
    TriaRawIterator (Triangulation<dim> *parent,
		     const int           level,
		     const int           index,
		     const void         *local_data=0);

				     /**
				      *  @name Dereferencing
				      */
				     /*@{*/
				     /**
				      *  Dereferencing operator, returns a
				      *  reference to an accessor.
				      *  Usage is thus like #(*i).index ();#
				      *
				      *  This function has to be specialized
				      *  explicitely for the different
				      *  #Pointee#s, to allow an
				      *  #iterator<1,TriangulationLevel<1>::LinesData>#
				      *  to point to #tria->lines.lines[index]# while
				      *  for one dimension higher it has
				      *  to point to #tria->quads.quads[index]#.
				      *
				      *  You must not dereference invalid or
				      *  past the end iterators.
				      */
    const Accessor & operator * () const;
    
				     /**
				      *  Dereferencing operator, non-#const#
				      *  version.
				      */
    Accessor & operator * ();
        
				     /**
				      *  Dereferencing operator, returns a
				      *  reference of the cell pointed to.
				      *  Usage is thus like #i->index ();#
				      *
				      *  There is a #const# and a non-#const#
				      *  version.
				      */
    const Accessor * operator -> () const;
        
				     /**
				      *  Dereferencing operator, non-#const#
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
				      * which is needed to make #map#s with
				      * iterators being keys. An iterator
				      * pointing to an element #a# is
				      * less than another iterator pointing
				      * to an element #b# if
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
				      *  Prefix #++# operator: #++i#. This
				      *  operator advances the iterator to
				      *  the next element and returns
				      *  a reference to #*this#.
				      *
				      *  The next element is next on this
				      *  level if there are more. If the
				      *  present element is the last on
				      *  this level, the first on the
				      *  next level is accessed.
				      */
    TriaRawIterator & operator ++ ();
    
				     /**
				      *  Postfix #++# operator: #i++#. This
				      *  operator advances the iterator to
				      *  the next element, but
				      *  returns an iterator to the element
				      *  priviously pointed to. Since this
				      *  involves a temporary and a copy
				      *  operation and since an
				      *  #iterator# is quite a large
				      *  object for a pointer, use the
				      *  prefix operator #++i# whenever
				      *  possible, especially in the head
				      *  of for loops
				      *  (#for (; i!=end; ++i)#) since there
				      *  you normally never need the
				      *  returned value.
				      */
    TriaRawIterator operator ++ (int);

    				     /**
				      *  Prefix #--# operator: #--i#. This
				      *  operator advances the iterator to
				      *  the previous element and returns
				      *  a reference to #*this#.
				      *
				      *  The previous element is previous on
				      *  this level if #index>0#. If the
				      *  present element is the first on
				      *  this level, the last on the
				      *  previous level is accessed.
				      */
    TriaRawIterator & operator -- ();
    
				     /**
				      *  Postfix #--# operator: #i--#. This
				      *  operator advances the iterator to
				      *  the previous element, but
				      *  returns an iterator to the element
				      *  priviously pointed to. Since this
				      *  involves a temporary and a copy
				      *  operation and since an
				      *  #iterator# is quite a large
				      *  object for a pointer, use the
				      *  prefix operator #--i# whenever
				      *  possible, especially in the head
				      *  of for loops
				      *  (#for (; i!=end; --i)#) since there
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
				      * Print the iterator in the form
				      * #level.index# to #out#.
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
				     /*@}*/
//  protected:  // don't know why we can't declare this protected (gcc2.7/8 chokes!?)
    Accessor accessor;
};







/**
 *   This specialization of \Ref{TriaRawIterator} provides access only to the
 *   {\it used} lines, quads, cells, etc.
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
				      *  pointed to.
				      *
				      *  If the object pointed to is not
				      *  past-the-end and is not
				      *  used, the debug version raises
				      *  an error!
				      */
    TriaIterator (Triangulation<dim> *parent,
		  const int           level,
		  const int           index,
		  const void         *local_data=0);
    
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
				      *  Prefix #++# operator: #++i#. This
				      *  operator advances the iterator to
				      *  the next used element and returns
				      *  a reference to #*this#.
				      */
    TriaIterator<dim,Accessor> & operator ++ ();

				     /**
				      *  Postfix #++# operator: #i++#. This
				      *  operator advances the iterator to
				      *  the next used element, but
				      *  returns an iterator to the element
				      *  previously pointed to. Since this
				      *  involves a temporary and a copy
				      *  operation and since an
				      *  #active_iterator# is quite a large
				      *  object for a pointer, use the
				      *  prefix operator #++i# whenever
				      *  possible, especially in the head
				      *  of for loops
				      *  (#for (; i!=end; ++i)#) since there
				      *  you normally never need the
				      *  returned value.
				      */
    TriaIterator<dim,Accessor> operator ++ (int);

    				     /**
				      *  Prefix #--# operator: #--i#. This
				      *  operator advances the iterator to
				      *  the previous used element and returns
				      *  a reference to #*this#.
				      */
    TriaIterator<dim,Accessor> & operator -- ();

				     /**
				      *  Postfix #--# operator: #i--#.
				      */
    TriaIterator<dim,Accessor> operator -- (int);
				     /*@}*/
    
				     /**
				      *  Exception
				      */
    DeclException0 (ExcAssignmentOfUnusedObject);
};







/**
 *   This specialization of \Ref{TriaIterator} provides access only to the
 *   {\it active} lines, quads, cells, etc. An active cell is a cell which is not
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
				      *  pointed to.
				      *
				      *  If the object pointed to is not
				      *  past-the-end and is not
				      *  active, the debug version raises
				      *  an error!
				      */
    TriaActiveIterator (Triangulation<dim> *parent,
			const int           level,
			const int           index,
			const void         *local_data=0);

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
				      *  Prefix #++# operator: #++i#. This
				      *  operator advances the iterator to
				      *  the next active element and returns
				      *  a reference to #*this#.
				      */
    TriaActiveIterator<dim,Accessor> & operator ++ ();

				     /**@name Advancement of iterators*/
				     /*@{*/
				     /**
				      *  Postfix #++# operator: #i++#. This
				      *  operator advances the iterator to
				      *  the next active element, but
				      *  returns an iterator to the element
				      *  previously pointed to. Since this
				      *  involves a temporary and a copy
				      *  operation and since an
				      *  #active_iterator# is quite a large
				      *  object for a pointer, use the
				      *  prefix operator #++i# whenever
				      *  possible, especially in the head
				      *  of for loops
				      *  (#for (; i!=end; ++i)#) since there
				      *  you normally never need the
				      *  returned value.
				      */
    TriaActiveIterator<dim,Accessor> operator ++ (int);

    				     /**
				      *  Prefix #--# operator: #--i#. This
				      *  operator advances the iterator to
				      *  the previous active element and
				      *  returns a reference to #*this#.
				      */
    TriaActiveIterator<dim,Accessor> & operator -- ();

				     /**
				      *  Postfix #--# operator: #i--#.
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
const Accessor &
TriaRawIterator<dim,Accessor>::operator * () const {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  return accessor;
};




template <int dim, typename Accessor>
inline
Accessor &
TriaRawIterator<dim,Accessor>::operator * () {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  return accessor;
};



template <int dim, typename Accessor>
inline
const Accessor *
TriaRawIterator<dim,Accessor>::operator -> () const {
  return &(this->operator* ());
};



template <int dim, typename Accessor>
inline
Accessor *
TriaRawIterator<dim,Accessor>::operator -> () {
  return &(this->operator* ());
};




template <int dim, typename Accessor>
inline
IteratorState TriaRawIterator<dim,Accessor>::state () const {
  return accessor.state ();
};



template <int dim, typename Accessor>
inline
bool TriaRawIterator<dim,Accessor>::operator < (const TriaRawIterator &i) const {
  Assert (state() != invalid, ExcDereferenceInvalidObject());
  Assert (i.state() != invalid, ExcDereferenceInvalidObject());
  
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
TriaRawIterator<dim,Accessor>::operator ++ () {
  Assert (state() == valid, ExcAdvanceInvalidObject());

  ++accessor;
  return *this;
};



template <int dim, typename Accessor>
inline
TriaRawIterator<dim,Accessor> &
TriaRawIterator<dim,Accessor>::operator -- () {
  Assert (state() == valid, ExcAdvanceInvalidObject());

  --accessor;
  return *this;
};



template <int dim, typename Accessor>
inline
void TriaRawIterator<dim,Accessor>::print (ostream &out) const {
  out << accessor.level() << "." << accessor.index();
};



template <int dim, typename Accessor>
inline
ostream & operator << (ostream &out, const TriaRawIterator<dim,Accessor> &i) {
  i.print(out);
  return out;
};



template <int dim, typename Accessor>
inline
ostream & operator << (ostream &out, const TriaIterator<dim,Accessor> &i) {
  i.print(out);
  return out;
};



template <int dim, typename Accessor>
inline
ostream & operator << (ostream &out, const TriaActiveIterator<dim,Accessor> &i) {
  i.print(out);
  return out;
};



/*----------------------------   tria-iterator.h     ---------------------------*/
/* end of #ifndef __tria_iterator_H */
#endif
/*----------------------------   tria-iterator.h     ---------------------------*/



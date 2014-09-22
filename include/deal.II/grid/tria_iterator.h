// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__tria_iterator_h
#define __deal2__tria_iterator_h


/*----------------------------   tria-iterator.h     ---------------------------*/


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>

#include <deal.II/base/point.h>
#include <deal.II/grid/tria_iterator_base.h>

#include <iterator>

#include <ostream>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim> class Triangulation;
template <int, int, int> class TriaAccessorBase;

template <typename> class TriaRawIterator;
template <typename> class TriaIterator;
template <typename> class TriaActiveIterator;



// note: in non-debug mode, i.e. with optimizations, the file
// tria_iterator.templates.h is included at the end of this file.
// this includes a lot of templates and thus makes compilation
// slower, but at the same time allows for more aggressive
// inlining and thus faster code.


/**
 * This class implements an iterator, analogous to those of the
 * standard template library (STL). It fulfills the requirements of a
 * bidirectional iterator.  See the C++ documentation for further
 * details of iterator specification and usage.
 *
 *
 * In addition to the STL
 * iterators an iterator of this class provides a <tt>-@></tt>
 * operator, i.e. you can write statements like
 * @code
 * i->set_refine_flag ();
 * @endcode
 *
 * Iterators are used whenever a loop over all lines, quads, cells
 * etc.  is to be performed. These loops can then be coded like this:
 * @code
 *   cell_iterator i   = tria.begin();
 *   cell_iterator end = tria.end();
 *   for (; i!=end; ++i)
 *     if (cell->at_boundary())
 *       cell->set_refine_flag();
 * @endcode
 *
 * Note the usage of <tt>++i</tt> instead of <tt>i++</tt> since this
 * does not involve temporaries and copying. It is recommended to use
 * a fixed value <tt>end</tt> inside the loop instead of
 * <tt>tria.end()</tt>, since the creation and copying of these
 * iterators is rather expensive compared to normal pointers.
 *
 * The objects pointed to are accessors, derived from
 * TriaAccessorBase. Which kind of accessor is determined by the template
 * argument <em>Accessor</em>. These accessors are not so much data
 * structures as they are a collection of functions providing access
 * to the data stored in Tringulation or DoFHandler objects. Using
 * these accessors, the structure of these classes is hidden from the
 * application program.
 *
 * <h3>Which iterator to use when</h3>
 *
 * @attention Application programs will rarely use TriaRawIterator,
 * but rather one of the derived classes TriaIterator or
 * TriaActiveIterator.
 *
 * <ul>
 * <li> TriaRawIterator objects point to lines, cells, etc in the
 * lists whether they are used or not (in the vectors, also
 * <i>dead</i> objects are stored, since deletion in vectors is
 * expensive and we also do not want to destroy the ordering induced
 * by the numbering in the vectors). Therefore not all raw iterators
 * point to valid objects.
 *
 * <li> The derived class TriaIterator selects the valid cells, that
 * is, cells used somewhere in the triangulation hierarchy.
 *
 * <li> TriaActiveIterator objects which only loop over active cells.
 * </ul>
 *
 * <h3>Purpose</h3>
 *
 * Iterators are not much slower than operating directly on the data
 * structures, since they perform the loops that you had to handcode
 * yourself anyway. Most iterator and accessor functions are inlined.
 *
 * The main functionality of iterators, resides in the <tt>++</tt> and
 * <tt>--</tt> operators. These move the iterator forward or backward
 * just as if it were a pointer into an array. Here, this operation is
 * not so easy, since it may include skipping some elements and the
 * transition between the triangulation levels. This is completely
 * hidden from the user, though you can still create an iterator
 * pointing to an arbitrary element.  Actually, the operation of
 * moving iterators back and forth is not done in the iterator
 * classes, but rather in the accessor classes. Since these are passed
 * as template arguments, you can write your own versions here to add
 * more functionality.
 *
 * Furthermore, the iterators described here satisfy the requirement of
 * input and bidirectional iterators as stated by the C++ standard and
 * the STL documentation. It is therefore possible to use the
 * functions from the algorithm section of the C++ standard,
 * e.g. <em>count_if</em> (see the documentation for Triangulation for
 * an example) and several others.
 *
 * <h3>Implementation</h3>
 *
 * The iterator class itself does not have much functionality. It only
 * becomes useful when assigned an Accessor (the second template
 * parameter), which really does the access to data. An Accessor
 * has to fulfil some requirements:
 *
 * <ul>
 * <li> It must have two members named <tt>present_level</tt> and
 * <tt>present_index</tt> storing the address of the element in the
 * triangulation presently pointed to. These data have to be
 * accessible by all triangulation iterators listed above.
 *
 * <li> It must have a constructor which takes a Triangulation* and
 * two unsigned integers, denoting the initial level and index, as
 * well as a data object depending on its type.
 *
 * <li> For the TriaIterator and the TriaActiveIterator class, it must
 * have a member function <tt>bool used()</tt>, for the latter a
 * member function <tt>bool active()</tt>.
 *
 * <li> It must have void operators <tt>++</tt> and <tt>--</tt>.
 *
 * <li> It must declare a local <tt>typedef AccessorData</tt> which
 * states the data type the accessor expects to get passed as fourth
 * constructor argument. By declaring a local data type, the
 * respective iterator class may type-safely enforce that data type to
 * be one of its own constructor argument types. If an accessor class
 * does not need additional data, this type shall be <tt>void</tt>.
 * </ul>
 *
 * Then the iterator is able to do what it is supposed to. All of the
 * necessary functions are implemented in the <tt>Accessor</tt> base
 * class, but you may write your own version (non-virtual, since we
 * use templates) to add functionality.
 *
 * The accessors provided by the library are distributed in three
 * groups, determined by whether they access the data of
 * Triangulation, DoFHandler or MGDoFHandler. They are derived from
 * TriaAccessor, DoFAccessor and MGDoFAccessor,
 * respectively. In each group, there is an accessor to cells, which
 * have more functionality.
 *
 * @attention It seems impossible to preserve constness of a
 * triangulation through iterator usage. Thus, if you declare pointers
 * to a <tt>const</tt> triangulation object, you should be well aware
 * that you might involuntarily alter the data stored in the
 * triangulation.
 *
 * @note More information on valid and invalid iterators can be found
 * in the documentation of TriaAccessorBase, where the iterator states are
 * checked and implemented.
 *
 *
 * <h3>Past-the-end iterators</h3>
 *
 * There is a representation of past-the-end-pointers, denoted by special
 * values of the member variables @p present_level and @p present_index:
 * If <tt>present_level>=0</tt> and <tt>present_index>=0</tt>, then the object is valid
 * (there is no check whether the triangulation really has that
 * many levels or that many cells on the present level when we investigate
 * the state of an iterator; however, in many places where an iterator is
 * dereferenced we make this check);
 * if <tt>present_level==-1</tt> and <tt>present_index==-1</tt>, then the iterator points
 * past the end; in all other cases, the iterator is considered invalid.
 * You can check this by calling the <tt>state()</tt> function.
 *
 * An iterator is also invalid, if the pointer pointing to the Triangulation
 * object is invalid or zero.
 *
 * Finally, an iterator is invalid, if the element pointed to by
 * @p present_level and @p present_index is not used, i.e. if the @p used
 * flag is set to false.
 *
 * The last two checks are not made in <tt>state()</tt> since both cases should only
 * occur upon unitialized construction through @p memcpy and the like (the
 * parent triangulation can only be set upon construction). If
 * an iterator is constructed empty through the empty constructor,
 * <tt>present_level==-2</tt> and <tt>present_index==-2</tt>. Thus, the iterator is
 * invalid anyway, regardless of the state of the triangulation pointer
 * and the state of the element pointed to.
 *
 * Past-the-end iterators may also be used to compare an iterator with the
 * <i>before-the-start</i> value, when running backwards. There is no
 * distinction between the iterators pointing past the two ends of a vector.
 *
 * By defining only one value to be past-the-end and making all other values
 * invalid provides a second track of security: if we should have forgotten
 * a check in the library when an iterator is incremented or decremented,
 * we automatically convert the iterator from the allowed state "past-the-end"
 * to the disallowed state "invalid" which increases the chance that somehwen
 * earlier than for past-the-end iterators an exception is raised.
 *
 * @ref Triangulation
 * @ingroup grid
 * @ingroup Iterators
 * @author Wolfgang Bangerth, 1998
 * @author documentation update Guido Kanschat, 2004
 */
template <typename Accessor>
class TriaRawIterator : public std::iterator<std::bidirectional_iterator_tag,Accessor>
{
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
   * possible to construct objects for
   * derived iterators:
   * @code
   * DoFCellAccessor dof_accessor;
   * Triangulation::active_cell_iterator cell
   *   = accessor;
   * @endcode
   */
  explicit TriaRawIterator (const Accessor &a);

  /**
   * Constructor. Assumes that the
   * other accessor type is
   * convertible to the current
   * one.
   */
  template <typename OtherAccessor>
  explicit TriaRawIterator (const OtherAccessor &a);

  /**
   *  Proper constructor, initialized
   *  with the triangulation, the
   *  level and index of the object
   *  pointed to. The last parameter is
   *  of a type declared by the accessor
   *  class.
   */
  TriaRawIterator (const Triangulation<Accessor::dimension,Accessor::space_dimension> *parent,
                   const int level,
                   const int index,
                   const typename AccessorType::AccessorData *local_data = 0);

  /**
   * This is a conversion operator
   * (constructor) which takes another
   * iterator type and copies the data;
   * this conversion works, if there is
   * a conversion path from the
   * @p OtherAccessor class to the @p Accessor
   * class of this object. One such path
   * would be derived class to base class,
   * which for example may be used to get
   * a Triangulation::raw_cell_iterator from
   * a DoFHandler::raw_cell_iterator, since
   * the DoFAccessor class is derived from
   * the TriaAccessorBase class.
   */
  template <typename OtherAccessor>
  TriaRawIterator (const TriaRawIterator<OtherAccessor> &i);

  /**
   * Another conversion operator,
   * where we use the pointers to
   * the Triangulation from a
   * TriaAccessorBase object, while
   * the additional data is used
   * according to the actual type
   * of Accessor.
   */
  TriaRawIterator (const TriaAccessorBase<Accessor::structure_dimension,Accessor::dimension,Accessor::space_dimension> &tria_accessor,
                   const typename Accessor::AccessorData *local_data);

  /**
   * Conversion constructor. Same
   * as above with the difference
   * that it converts from
   * TriaIterator classes (not
   * TriaRawIterator).
   */
  template <typename OtherAccessor>
  TriaRawIterator (const TriaIterator<OtherAccessor> &i);

  /**
   * Conversion constructor. Same
   * as above with the difference
   * that it converts from
   * TriaActiveIterator classes
   * (not TriaRawIterator).
   */
  template <typename OtherAccessor>
  TriaRawIterator (const TriaActiveIterator<OtherAccessor> &i);

  /**
   *  @name Dereferencing
   */
  /*@{*/
  /**
   *  Dereferencing operator, returns a
   *  reference to an accessor.
   *  Usage is thus like <tt>(*i).index ();</tt>
   *
   *  This function has to be specialized
   *  explicitly for the different
   *  @p Pointees, to allow an
   *  <tt>iterator<1,TriangulationLevel<1>::LinesData></tt>
   *  to point to <tt>tria->lines.cells[index]</tt> while
   *  for one dimension higher it has
   *  to point to <tt>tria->quads.cells[index]</tt>.
   *
   *  You must not dereference invalid or
   *  past the end iterators.
   */
  const Accessor &operator * () const;

  /**
   *  Dereferencing operator, non-@p const
   *  version.
   */
  Accessor &operator * ();

  /**
   *  Dereferencing operator, returns a
   *  reference of the cell pointed to.
   *  Usage is thus like <tt>i->index ();</tt>
   *
   *  There is a @p const and a non-@p const
   *  version.
   */
  const Accessor *operator -> () const;

  /**
   *  Dereferencing operator, non-@p const
   *  version.
   */
  Accessor *operator -> ();


  /**
  * In order be able to assign
  * end-iterators for different
  * accessors to each other, we
  * need an access function which
  * returns the accessor
  * regardless of its state.
  *
  * @warning This function should
  * not be used in application
  * programs. It is only intended
  * for limited purposes inside
  * the library and it makes
  * debugging much harder.
  */
  const Accessor &access_any () const;

  /*@}*/

  /**
   *  Assignment operator.
   */
  TriaRawIterator &operator = (const TriaRawIterator &);

  /**
   *  Assignment operator.
   */
//    template <class OtherAccessor>
//    TriaRawIterator & operator = (const TriaRawIterator<OtherAccessor>&);

  /**
   *  Assignment operator.
   */
//    template <class OtherAccessor>
//    TriaRawIterator & operator = (const TriaIterator<OtherAccessor>&);

  /**
   *  Assignment operator.
   */
//    template <class OtherAccessor>
//    TriaRawIterator & operator = (const TriaActiveIterator<OtherAccessor>&);

  /**
   *  Compare for equality.
   */
  bool operator == (const TriaRawIterator &) const;

  /**
   *  Compare for inequality.
   */
  bool operator != (const TriaRawIterator &) const;

  /**
   * Ordering relation for iterators.
   *
   * This relation attempts a total ordering of cells.
   *
   * The relation is defined as follows:
   *
   * For objects of <tt>Accessor::structure_dimension <
   * Accessor::dimension</tt>, we simply compare the index of such an
   * object. The ordering is lexicographic
   * according to the following hierarchy (in the sense, that the next
   * test is only applied if the previous was inconclusive):
   *
   * <ol>
   * <li> The past-the-end iterator is always ordered last. Two
   * past-the-end iterators rank the same, thus false is returned in
   * that case.</li>
   *
   * <li> The level of the cell.</li>
   * <li> The index of a cell inside the level.</li>
   * </ol>
   *
   * @note The ordering is not consistent between different processor in
   * a parallel::distributed::Triangulation because we rely on index(),
   * which is likely not the same.
   */
  bool operator < (const TriaRawIterator &) const;

  /**@name Advancement of iterators*/
  /*@{*/
  /**
   *  Prefix <tt>++</tt> operator: <tt>++iterator</tt>. This
   *  operator advances the iterator to
   *  the next element and returns
   *  a reference to <tt>*this</tt>.
   */
  TriaRawIterator &operator ++ ();

  /**
   *  Postfix <tt>++</tt> operator: <tt>iterator++</tt>. This
   *  operator advances the iterator to
   *  the next element, but
   *  returns an iterator to the element
   *  previously pointed to.
   *
   *  Since this operation
   *  involves a temporary and a copy
   *  operation and since an
   *  @p iterator is quite a large
   *  object for a pointer, use the
   *  prefix operator <tt>++iterator</tt> whenever
   *  possible, especially in the header
   *  of for loops
   *  (<tt>for (; iterator!=end; ++iterator)</tt>) since there
   *  you normally never need the
   *  returned value.
   */
  TriaRawIterator operator ++ (int);

  /**
   *  Prefix @p -- operator: @p --iterator. This
   *  operator moves the iterator to
   *  the previous element and returns
   *  a reference to <tt>*this</tt>.
   */
  TriaRawIterator &operator -- ();

  /**
   *  Postfix @p -- operator: @p iterator--. This
   *  operator moves the iterator to
   *  the previous element, but
   *  returns an iterator to the element
   *  previously pointed to.
   *
   *  The same applies as for the postfix operator++: If possible,
   *  avoid it by using the prefix operator form to avoid the use
   *  of a temporary variable.
   */
  TriaRawIterator operator -- (int);
  /*@}*/

  /**
   *  Return the state of the iterator.
   */
  IteratorState::IteratorStates state () const;

  /**
   * Print the iterator to a stream
   * <code>out</code>. The
   * format is <tt>level.index</tt>.
   */
  template <class STREAM>
  void print (STREAM &out) const;


  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   */
  std::size_t memory_consumption () const;


  /**@name Exceptions*/
  /*@{*/
  /**
   *  Exception for TriaObjects with
   *  level, i.e. cells.
   */
  DeclException1 (ExcDereferenceInvalidCell,
                  Accessor,
                  << "You tried to dereference a cell iterator for which this "
                  << "is not possible. More information on this iterator: "
                  << "level=" << arg1.level()
                  << ", index=" << arg1.index()
                  << ", state="
                  << (arg1.state() == IteratorState::valid ? "valid" :
                      (arg1.state() == IteratorState::past_the_end ?
                       "past_the_end" : "invalid")));

  /**
   *  Exception for
   *  lower-dimensional TriaObjects
   *  without level, i.e. objects
   *  faces are constructed with.
   */
  DeclException1 (ExcDereferenceInvalidObject,
                  Accessor,
                  << "You tried to dereference an iterator for which this "
                  << "is not possible. More information on this iterator: "
                  << "index=" << arg1.index()
                  << ", state="
                  << (arg1.state() == IteratorState::valid ? "valid" :
                      (arg1.state() == IteratorState::past_the_end ?
                       "past_the_end" : "invalid")));

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
  template <typename SomeAccessor> friend class TriaRawIterator;
  template <typename SomeAccessor> friend class TriaIterator;
  template <typename SomeAccessor> friend class TriaActiveIterator;
};


/**
 *   This specialization of TriaRawIterator provides access only to
 *   the <em>used</em> lines, quads, cells, etc.
 *
 * @ingroup grid
 * @ingroup Iterators
 */
template <typename Accessor>
class TriaIterator : public TriaRawIterator<Accessor>
{
public:
  /**
   *  Empty constructor. Such an object
   *  is not usable!
   */
  TriaIterator ();

  /**
   *  Copy constructor.
   */
  TriaIterator (const TriaIterator<Accessor> &);

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
  TriaIterator (const TriaRawIterator<Accessor> &);

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
  TriaIterator (const Triangulation<Accessor::dimension,Accessor::space_dimension> *parent,
                const int                 level,
                const int                 index,
                const typename Accessor::AccessorData *local_data = 0);

  /**
   * Construct from an accessor of type OtherAccessor that is convertible
   * to the type Accessor.
   */
  template <typename OtherAccessor>
  explicit TriaIterator (const OtherAccessor &a);

  /**
   * This is a conversion operator
   * (constructor) which takes another
   * iterator type and copies the data;
   * this conversion works, if there is
   * a conversion path from the
   * @p OtherAccessor class to the @p Accessor
   * class of this object. One such path
   * would be derived class to base class,
   * which for example may be used to get
   * a Triangulation::cell_iterator from
   * a DoFHandler::cell_iterator, since
   * the DoFAccessor class is derived from
   * the TriaAccessorBase class.
   */
  template <typename OtherAccessor>
  TriaIterator (const TriaIterator<OtherAccessor> &i);

  /**
   * Another conversion operator,
   * where we use the pointers to
   * the Triangulation from a
   * TriaAccessorBase object, while
   * the additional data is used
   * according to the actual type
   * of Accessor.
   */
  TriaIterator (const TriaAccessorBase<Accessor::structure_dimension,Accessor::dimension,Accessor::space_dimension> &tria_accessor,
                const typename Accessor::AccessorData *local_data);

  /**
   * Similar conversion operator to the above
   * one, but does a check whether the
   * iterator points to a used element,
   * which is necessary for raw iterators.
   */
  template <typename OtherAccessor>
  TriaIterator (const TriaRawIterator<OtherAccessor> &i);

  /**
   * Similar conversion operator to
   * the above one, but for
   * conversion from active
   * iterators.
   */
  template <typename OtherAccessor>
  TriaIterator (const TriaActiveIterator<OtherAccessor> &i);

  /**
   *  Assignment operator.
   */
  TriaIterator<Accessor> &
  operator = (const TriaIterator<Accessor> &);

  /**
   *  Cross assignment operator. This
   *  assignment is only valid if the
   *  given iterator points to a used
   *  element.
   */
  TriaIterator<Accessor> &
  operator = (const TriaRawIterator<Accessor> &);

  /**
   *  Assignment
   *  operator. Requires, that
   *  Accessor can be copied from
   *  OtherAccessor.
   */
  template <class OtherAccessor>
  TriaIterator<Accessor> &
  operator = (const TriaIterator<OtherAccessor> &);

  /**
   *  Cross assignment operator. This
   *  assignment is only valid if the
   *  given iterator points to a used
   *  element. Requires, that
   *  Accessor can be copied from
   *  OtherAccessor.
   */
  template <class OtherAccessor>
  TriaIterator<Accessor> &
  operator = (const TriaRawIterator<OtherAccessor> &);

  /**@name Advancement of iterators*/
  /*@{*/
  /**
   *  Prefix <tt>++</tt> operator: <tt>++i</tt>. This
   *  operator advances the iterator to
   *  the next used element and returns
   *  a reference to <tt>*this</tt>.
   */
  TriaIterator<Accessor> &operator ++ ();

  /**
   *  Postfix <tt>++</tt> operator: <tt>i++</tt>. This
   *  operator advances the iterator to
   *  the next used element, but
   *  returns an iterator to the element
   *  previously pointed to. Since this
   *  involves a temporary and a copy
   *  operation and since an
   *  @p active_iterator is quite a large
   *  object for a pointer, use the
   *  prefix operator <tt>++i</tt> whenever
   *  possible, especially in the head
   *  of for loops
   *  (<tt>for (; i!=end; ++i)</tt>) since there
   *  you normally never need the
   *  returned value.
   */
  TriaIterator<Accessor> operator ++ (int);

  /**
   *  Prefix @p -- operator: @p --i. This
   *  operator advances the iterator to
   *  the previous used element and returns
   *  a reference to <tt>*this</tt>.
   */
  TriaIterator<Accessor> &operator -- ();

  /**
   *  Postfix @p -- operator: @p i--.
   */
  TriaIterator<Accessor> operator -- (int);
  /*@}*/

  /**
   *  Exception
   */
  DeclException0 (ExcAssignmentOfUnusedObject);
};


/**
 *   This specialization of TriaIterator provides access only to the
 *   <em>active</em> lines, quads, cells, etc. An active cell is a
 *   cell which is not refined and thus a cell on which calculations
 *   on the finest level are done.
 *
 * @ingroup grid
 * @ingroup Iterators
 */
template <typename Accessor>
class TriaActiveIterator : public TriaIterator<Accessor>
{
public:
  /**
   *  Empty constructor. Such an object
   *  is not usable!
   */
  TriaActiveIterator ();

  /**
   *  Copy constructor.
   */
  TriaActiveIterator (const TriaActiveIterator<Accessor> &);

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
  TriaActiveIterator (const TriaRawIterator<Accessor> &);

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
  TriaActiveIterator (const TriaIterator<Accessor> &);

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
  TriaActiveIterator (const Triangulation<Accessor::dimension,Accessor::space_dimension> *parent,
                      const int level,
                      const int index,
                      const typename Accessor::AccessorData *local_data = 0);

  /**
   * This is a conversion operator
   * (constructor) which takes another
   * iterator type and copies the data;
   * this conversion works, if there is
   * a conversion path from the
   * @p OtherAccessor class to the @p Accessor
   * class of this object. One such path
   * would be derived class to base class,
   * which for example may be used to get
   * a Triangulation::active_cell_iterator from
   * a DoFHandler::active_cell_iterator, since
   * the DoFAccessor class is derived from
   * the TriaAccessorBase class.
   */
  template <typename OtherAccessor>
  TriaActiveIterator (const TriaActiveIterator<OtherAccessor> &i);

  /**
  * Another conversion operator,
  * where we use the pointers to
  * the Triangulation from a
  * TriaAccessorBase object, while
  * the additional data is used
  * according to the actual type
  * of Accessor.
  */
  TriaActiveIterator (const TriaAccessorBase<Accessor::structure_dimension,Accessor::dimension,Accessor::space_dimension> &tria_accessor,
                      const typename Accessor::AccessorData *local_data);

  /**
   * Similar conversion operator to the above
   * one, but does a check whether the
   * iterator points to a used element,
   * and is active, which is necessary for
   * raw iterators. Since usual iterators
   * are also raw iterators, this constructor
   * works also for parameters of type
   * <tt>TriaIterator<OtherAccessor></tt>.
   */
  template <typename OtherAccessor>
  TriaActiveIterator (const TriaRawIterator<OtherAccessor> &i);

  /**
   *  Assignment operator.
   */
  TriaActiveIterator<Accessor> &
  operator = (const TriaActiveIterator<Accessor> &);

  /**
   *  Cross assignment operator. This
   *  assignment is only valid if the
   *  given iterator points to an active
   *  element.
   */
  TriaActiveIterator<Accessor> &
  operator = (const TriaIterator<Accessor> &);

  /**
   *  Cross assignment operator. This
   *  assignment is only valid if the
   *  given iterator points to an active
   *  element or past the end.
   */
  TriaActiveIterator<Accessor> &
  operator = (const TriaRawIterator<Accessor> &);

  /**
   *  Assignment operator. Requires, that
   *  Accessor can be copied from
   *  OtherAccessor.
   */
  template <class OtherAccessor>
  TriaActiveIterator<Accessor> &
  operator = (const TriaActiveIterator<OtherAccessor> &);

  /**
   *  Cross assignment operator. This
   *  assignment is only valid if the
   *  given iterator points to an active
   *  element or past the end. Requires, that
   *  Accessor can be copied from
   *  OtherAccessor.
   */
  template <class OtherAccessor>
  TriaActiveIterator<Accessor> &
  operator = (const TriaRawIterator<OtherAccessor> &);

  /**
   *  Cross assignment operator. This
   *  assignment is only valid if the
   *  given iterator points to an active
   *  element. Requires, that
   *  Accessor can be copied from
   *  OtherAccessor.
   */
  template <class OtherAccessor>
  TriaActiveIterator<Accessor> &
  operator = (const TriaIterator<OtherAccessor> &);

  /**
   *  Prefix <tt>++</tt> operator: <tt>++i</tt>. This
   *  operator advances the iterator to
   *  the next active element and returns
   *  a reference to <tt>*this</tt>.
   */
  TriaActiveIterator<Accessor> &operator ++ ();

  /**@name Advancement of iterators*/
  /*@{*/
  /**
   *  Postfix <tt>++</tt> operator: <tt>i++</tt>. This
   *  operator advances the iterator to
   *  the next active element, but
   *  returns an iterator to the element
   *  previously pointed to. Since this
   *  involves a temporary and a copy
   *  operation and since an
   *  @p active_iterator is quite a large
   *  object for a pointer, use the
   *  prefix operator <tt>++i</tt> whenever
   *  possible, especially in the head
   *  of for loops
   *  (<tt>for (; i!=end; ++i)</tt>) since there
   *  you normally never need the
   *  returned value.
   */
  TriaActiveIterator<Accessor> operator ++ (int);

  /**
   *  Prefix @p -- operator: @p --i. This
   *  operator advances the iterator to
   *  the previous active element and
   *  returns a reference to <tt>*this</tt>.
   */
  TriaActiveIterator<Accessor> &operator -- ();

  /**
   *  Postfix @p -- operator: @p i--.
   */
  TriaActiveIterator<Accessor> operator -- (int);
  /*@}*/

  /**
   *  Exception
   */
  DeclException0 (ExcAssignmentOfInactiveObject);
};


/*----------------------- Inline functions -------------------*/


template <typename Accessor>
inline
TriaRawIterator<Accessor>::
TriaRawIterator (const Accessor &a)
  :
  accessor (a)
{}



template <typename Accessor>
template <typename OtherAccessor>
inline
TriaRawIterator<Accessor>::
TriaRawIterator (const OtherAccessor &a)
  :
  accessor (a)
{}



template <typename Accessor>
template <typename OtherAccessor>
inline
TriaRawIterator<Accessor>::
TriaRawIterator (const TriaRawIterator<OtherAccessor> &i)
  :
  accessor (i.accessor)
{}



template <typename Accessor>
template <typename OtherAccessor>
inline
TriaRawIterator<Accessor>::
TriaRawIterator (const TriaIterator<OtherAccessor> &i)
  :
  accessor (i.accessor)
{}



template <typename Accessor>
template <typename OtherAccessor>
inline
TriaRawIterator<Accessor>::
TriaRawIterator (const TriaActiveIterator<OtherAccessor> &i)
  :
  accessor (i.accessor)
{}



template <typename Accessor>
inline
const Accessor &
TriaRawIterator<Accessor>::operator * () const
{
  Assert (Accessor::structure_dimension!=Accessor::dimension ||
          state() == IteratorState::valid,
          ExcDereferenceInvalidCell(accessor));
  Assert (Accessor::structure_dimension==Accessor::dimension ||
          state() == IteratorState::valid,
          ExcDereferenceInvalidObject(accessor));

  return accessor;
}



template <typename Accessor>
inline
Accessor &
TriaRawIterator<Accessor>::operator * ()
{
  Assert (Accessor::structure_dimension!=Accessor::dimension ||
          state() == IteratorState::valid,
          ExcDereferenceInvalidCell(accessor));
  Assert (Accessor::structure_dimension==Accessor::dimension ||
          state() == IteratorState::valid,
          ExcDereferenceInvalidObject(accessor));

  return accessor;
}



template <typename Accessor>
inline
const Accessor &
TriaRawIterator<Accessor>::access_any () const
{
  return accessor;
}



template <typename Accessor>
inline
const Accessor *
TriaRawIterator<Accessor>::operator -> () const
{
  return &(this->operator* ());
}



template <typename Accessor>
inline
Accessor *
TriaRawIterator<Accessor>::operator -> ()
{
  return &(this->operator* ());
}



template <typename Accessor>
inline
IteratorState::IteratorStates
TriaRawIterator<Accessor>::state () const
{
  return accessor.state ();
}



template <typename Accessor>
inline
bool
TriaRawIterator<Accessor>::operator < (const TriaRawIterator<Accessor> &other) const
{
  Assert (state() != IteratorState::invalid, ExcDereferenceInvalidObject(accessor));
  Assert (other.state() != IteratorState::invalid, ExcDereferenceInvalidObject(other.accessor));

  Assert (&accessor.get_triangulation() == &other.accessor.get_triangulation(),
          ExcInvalidComparison());

  // Deal with iterators past end
  if (state()==IteratorState::past_the_end)
    return false;
  if (other.state()==IteratorState::past_the_end)
    return true;

  return ((**this) < (*other));
}


template <typename Accessor>
inline
TriaRawIterator<Accessor> &
TriaRawIterator<Accessor>::operator ++ ()
{
  Assert (state() == IteratorState::valid, ExcAdvanceInvalidObject());

  ++accessor;
  return *this;
}



template <typename Accessor>
inline
TriaRawIterator<Accessor> &
TriaRawIterator<Accessor>::operator -- ()
{
  Assert (state() == IteratorState::valid, ExcAdvanceInvalidObject());

  --accessor;
  return *this;
}



template <typename Accessor>
template <class STREAM>
inline
void
TriaRawIterator<Accessor>::print (STREAM &out) const
{
  if (Accessor::structure_dimension==Accessor::dimension)
    out << accessor.level() << "." << accessor.index();
  else
    out << accessor.index();
}



template <typename Accessor>
inline
std::size_t
TriaRawIterator<Accessor>::memory_consumption () const
{
  return sizeof(TriaRawIterator<Accessor>);
}



template <typename Accessor>
template <typename OtherAccessor>
inline
TriaIterator<Accessor>::TriaIterator (const TriaIterator<OtherAccessor> &i)
  :
  TriaRawIterator<Accessor> (i.accessor)
{}



template <typename Accessor>
template <typename OtherAccessor>
inline
TriaIterator<Accessor>::TriaIterator (const TriaActiveIterator<OtherAccessor> &i)
  :
  TriaRawIterator<Accessor> (i.accessor)
{}



template <typename Accessor>
template <typename OtherAccessor>
inline
TriaIterator<Accessor>::TriaIterator (const TriaRawIterator<OtherAccessor> &i)
  :
  TriaRawIterator<Accessor> (i.accessor)
{
#ifdef DEBUG
  // do this like this, because:
  // if we write
  // "Assert (IteratorState::past_the_end || used)"
  // used() is called anyway, even if
  // state==IteratorState::past_the_end, and will then
  // throw the exception!
  if (this->state() != IteratorState::past_the_end)
    Assert (this->accessor.used(),
            ExcAssignmentOfUnusedObject());
#endif
}

template <typename Accessor>
template <typename OtherAccessor>
TriaIterator<Accessor>::TriaIterator (const OtherAccessor &a)
  :
  TriaRawIterator<Accessor> (a)
{
#ifdef DEBUG
  // do this like this, because:
  // if we write
  // "Assert (IteratorState::past_the_end || used)"
  // used() is called anyway, even if
  // state==IteratorState::past_the_end, and will then
  // throw the exception!
  if (this->state() != IteratorState::past_the_end)
    Assert (this->accessor.used(),
            ExcAssignmentOfUnusedObject());
#endif
}

template <typename Accessor>
template <typename OtherAccessor>
inline
TriaActiveIterator<Accessor>::TriaActiveIterator (const TriaActiveIterator<OtherAccessor> &i)
  :
  TriaIterator<Accessor> (i.accessor)
{}



template <typename Accessor>
template <typename OtherAccessor>
inline
TriaActiveIterator<Accessor>::TriaActiveIterator (const TriaRawIterator<OtherAccessor> &i)
  :
  TriaIterator<Accessor> (i)
{
#ifdef DEBUG
  // do this like this, because:
  // if we write
  // "Assert (IteratorState::past_the_end || !has_children())"
  // has_children() is called anyway, even if
  // state==IteratorState::past_the_end, and will then
  // throw the exception!
  if (this->state() != IteratorState::past_the_end)
    Assert (this->accessor.has_children()==false,
            ExcAssignmentOfInactiveObject());
#endif
}



/**
 * Print the address to which this iterator points to @p out. The
 * address is given by the pair <tt>(level,index)</tt>, where @p index is
 * an index relative to the level in which the object is that is
 * pointed to.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <typename Accessor>
inline
std::ostream &operator << (std::ostream                        &out,
                           const TriaRawIterator<Accessor> &i)
{
  i.print(out);
  return out;
}



/**
 * Print the address to which this iterator points to @p out. The
 * address is given by the pair <tt>(level,index)</tt>, where @p index is
 * an index relative to the level in which the object is that is
 * pointed to.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <typename Accessor>
inline
std::ostream &operator << (std::ostream                     &out,
                           const TriaIterator<Accessor> &i)
{
  i.print(out);
  return out;
}



/**
 * Print the address to which this iterator points to @p out. The
 * address is given by the pair <tt>(level,index)</tt>, where @p index is
 * an index relative to the level in which the object is that is
 * pointed to.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <typename Accessor>
inline
std::ostream &operator << (std::ostream                           &out,
                           const TriaActiveIterator<Accessor> &i)
{
  i.print(out);
  return out;
}


DEAL_II_NAMESPACE_CLOSE


// if in optimized mode: include more templates
#ifndef DEBUG
#  include "tria_iterator.templates.h"
#endif


/*----------------------------   tria-iterator.h     ---------------------------*/
#endif
/*----------------------------   tria-iterator.h     ---------------------------*/

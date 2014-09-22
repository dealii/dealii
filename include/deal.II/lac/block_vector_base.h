// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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

#ifndef __deal2__block_vector_base_h
#define __deal2__block_vector_base_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/block_indices.h>
#include <deal.II/lac/vector.h>
#include <vector>
#include <iterator>
#include <cmath>
#include <cstddef>

DEAL_II_NAMESPACE_OPEN


/*! @addtogroup Vectors
 *@{
 */

template <typename> class BlockVectorBase;


/**
 * A class that can be used to determine whether a given type is a block
 * vector type or not. For example,
 * @code
 *   IsBlockVector<Vector<double> >::value
 * @endcode
 * has the value false, whereas
 * @code
 *   IsBlockVector<BlockVector<double> >::value
 * @endcode
 * is true. This is sometimes useful in template contexts where we may
 * want to do things differently depending on whether a template type
 * denotes a regular or a block vector type.
 *
 * @author Wolfgang Bangerth, 2010
 */
template <typename VectorType>
struct IsBlockVector
{
private:
  struct yes_type
  {
    char c[1];
  };
  struct no_type
  {
    char c[2];
  };

  /**
   * Overload returning true if the class
   * is derived from BlockVectorBase,
   * which is what block vectors do.
   */
  template <typename T>
  static yes_type check_for_block_vector (const BlockVectorBase<T> *);

  /**
   * Catch all for all other potential
   * vector types that are not block
   * matrices.
   */
  static no_type check_for_block_vector (...);

public:
  /**
   * A statically computable value that
   * indicates whether the template
   * argument to this class is a block
   * vector (in fact whether the type is
   * derived from BlockVectorBase<T>).
   */
  static const bool value = (sizeof(check_for_block_vector
                                    ((VectorType *)0))
                             ==
                             sizeof(yes_type));
};


// instantiation of the static member
template <typename VectorType>
const bool IsBlockVector<VectorType>::value;




namespace internal
{

  /**
   * Namespace in which iterators in block vectors are implemented.
   *
   * @author Wolfgang Bangerth, 2001
   */
  namespace BlockVectorIterators
  {
    /**
     * Declaration of the general
     * template of a structure which is
     * used to determine some types
     * based on the template arguments
     * of other classes.
     */
    template <class BlockVectorType, bool constness>
    struct Types
    {
    };



    /**
     * Declaration of a specialized
     * template of a structure which is
     * used to determine some types
     * based on the template arguments
     * of other classes.
     *
     * This is for the use of non-const
     * iterators.
     */
    template <class BlockVectorType>
    struct Types<BlockVectorType,false>
    {
      /**
       * Type of the vector
       * underlying the block vector
       * used in non-const
       * iterators. There, the
       * vector must not be constant.
       */
      typedef typename BlockVectorType::BlockType Vector;

      /**
       * Type of the block vector
       * used in non-const
       * iterators. There, the block
       * vector must not be constant.
       */
      typedef BlockVectorType BlockVector;

      /**
       * Type of the numbers we point
       * to. Here, they are not
       * constant.
       */
      typedef typename BlockVector::value_type value_type;

      /**
       * Typedef the result of a
       * dereferencing operation for an
       * iterator of the underlying
       * iterator.
       */
      typedef typename Vector::reference dereference_type;
    };



    /**
     * Declaration of a specialized
     * template of a structure which is
     * used to determine some types
     * based on the template arguments
     * of other classes.
     *
     * This is for the use of
     * const_iterator.
     */
    template <class BlockVectorType>
    struct Types<BlockVectorType,true>
    {
      /**
       * Type of the vector
       * underlying the block vector
       * used in
       * const_iterator. There,
       * the vector must be
       * constant.
       */
      typedef const typename BlockVectorType::BlockType Vector;

      /**
       * Type of the block vector
       * used in
       * const_iterator. There,
       * the block vector must be
       * constant.
       */
      typedef const BlockVectorType BlockVector;

      /**
       * Type of the numbers we point
       * to. Here, they are constant
       * since the block vector we
       * use is constant.
       */
      typedef const typename BlockVector::value_type value_type;

      /**
       * Typedef the result of a
       * dereferencing operation for an
       * iterator of the underlying
       * iterator. Since this is for
       * constant iterators, we can only
       * return values, no actual
       * references.
       */
      typedef value_type dereference_type;
    };


    /**
     * General random-access iterator
     * class for block vectors. Since
     * we do not want to have two
     * classes for non-const
     * iterator and
     * const_iterator, we take a
     * second template argument which
     * denotes whether the vector we
     * point into is a constant object
     * or not. The first template
     * argument is always the number
     * type of the block vector in use.
     *
     * This class satisfies all
     * requirements of random access
     * iterators defined in the C++
     * standard. Operations on these
     * iterators are constant in the
     * number of elements in the block
     * vector. However, they are
     * sometimes linear in the number
     * of blocks in the vector, but
     * since that does rarely change
     * dynamically within an
     * application, this is a constant
     * and we again have that the
     * iterator satisfies the
     * requirements of a random access
     * iterator.
     *
     * The implementation of this class
     * has to work around some problems
     * in compilers and standard
     * libraries. One of these requires
     * us to write all comparison
     * operators twice, once comparison
     * with iterators of the same type
     * and once with iterators pointing
     * to numbers of opposite constness
     * specification. The reason is
     * that if we would have written
     * the comparison operators as a
     * template on the constness of the
     * right hand side, then gcc2.95
     * signals an error that these
     * operators ambiguate operators
     * declared somewhere within the
     * standard library. Likewise, we
     * have to work around some
     * problems with granting other
     * iterators friendship. This makes
     * the implementation somewhat
     * non-optimal at places, but at
     * least everything works.
     *
     * @author Wolfgang Bangerth, 2001
     */
    template <class BlockVectorType, bool constness>
    class Iterator :
      public std::iterator<std::random_access_iterator_tag,
      typename Types<BlockVectorType,constness>::value_type>
    {
    private:
      /**
       * Typedef an iterator with
       * opposite constness
       * requirements on the elements
       * it points to.
       */
      typedef Iterator<BlockVectorType,!constness> InverseConstnessIterator;

    public:
      /**
       * Declare the type for container size.
       */
      typedef types::global_dof_index size_type;

      /**
       * Type of the number this
       * iterator points
       * to. Depending on the value
       * of the second template
       * parameter, this is either a
       * constant or non-const
       * number.
       */
      typedef
      typename Types<BlockVectorType,constness>::value_type
      value_type;

      /**
       * Declare some typedefs which
       * are standard for iterators
       * and are used by algorithms
       * to enquire about the
       * specifics of the iterators
       * they work on.
       */
      typedef std::random_access_iterator_tag               iterator_type;
      typedef std::ptrdiff_t                                difference_type;
      typedef typename BlockVectorType::reference           reference;
      typedef value_type                                   *pointer;

      typedef
      typename Types<BlockVectorType,constness>::dereference_type
      dereference_type;

      /**
       * Typedef the type of the
       * block vector (which differs
       * in constness, depending on
       * the second template
       * parameter).
       */
      typedef
      typename Types<BlockVectorType,constness>::BlockVector
      BlockVector;

      /**
       * Construct an iterator from
       * a vector to which we point
       * and the global index of
       * the element pointed to.
       *
       * Depending on the value of
       * the <tt>constness</tt> template
       * argument of this class,
       * the first argument of this
       * constructor is either is a
       * const or non-const
       * reference.
       */
      Iterator (BlockVector     &parent,
                const size_type  global_index);

      /**
       * Copy constructor.
       */
      Iterator (const Iterator<BlockVectorType,constness> &c);

      /**
       * Copy constructor for
       * conversion between iterators
       * with different constness
       * requirements. This
       * constructor throws an error
       * if an attempt is made at
       * converting a constant to a
       * non-constant iterator.
       */
      Iterator (const InverseConstnessIterator &c);

    private:
      /**
       * Constructor used internally
       * in this class. The arguments
       * match exactly the values of
       * the respective member
       * variables.
       */
      Iterator (BlockVector     &parent,
                const size_type  global_index,
                const size_type  current_block,
                const size_type  index_within_block,
                const size_type  next_break_forward,
                const size_type  next_break_backward);

    public:

      /**
       * Copy operator.
       */
      Iterator &operator = (const Iterator &c);

      /**
       * Dereferencing operator. If the
       * template argument
       * <tt>constness</tt> is
       * <tt>true</tt>, then no writing to
       * the result is possible, making
       * this a const_iterator.
       */
      dereference_type operator * () const;

      /**
       * Random access operator,
       * grant access to arbitrary
       * elements relative to the one
       * presently pointed to.
       */
      dereference_type operator [] (const difference_type d) const;

      /**
       * Prefix increment operator. This
       * operator advances the iterator to
       * the next element and returns a
       * reference to <tt>*this</tt>.
       */
      Iterator &operator ++ ();

      /**
       * Postfix increment
       * operator. This operator
       * advances the iterator to
       * the next element and
       * returns a copy of the old
       * value of this iterator.
       */
      Iterator operator ++ (int);

      /**
       * Prefix decrement operator. This
       * operator retracts the iterator to
       * the previous element and returns a
       * reference to <tt>*this</tt>.
       */
      Iterator &operator -- ();

      /**
       * Postfix decrement
       * operator. This operator
       * retracts the iterator to
       * the previous element and
       * returns a copy of the old
       * value of this iterator.
       */
      Iterator operator -- (int);

      /**
       * Compare for equality of
       * iterators. This operator
       * checks whether the vectors
       * pointed to are the same,
       * and if not it throws an
       * exception.
       */
      bool operator == (const Iterator &i) const;

      /**
       * Same, but compare with an
       * iterator of different
       * constness.
       */
      bool operator == (const InverseConstnessIterator &i) const;

      /**
       * Compare for inequality of
       * iterators. This operator
       * checks whether the vectors
       * pointed to are the same,
       * and if not it throws an
       * exception.
       */
      bool operator != (const Iterator &i) const;

      /**
       * Same, but compare with an
       * iterator of different
       * constness.
       */
      bool operator != (const InverseConstnessIterator &i) const;

      /**
       * Check whether this
       * iterators points to an
       * element previous to the
       * one pointed to by the
       * given argument. This
       * operator checks whether
       * the vectors pointed to are
       * the same, and if not it
       * throws an exception.
       */
      bool operator < (const Iterator &i) const;

      /**
       * Same, but compare with an
       * iterator of different
       * constness.
       */
      bool operator < (const InverseConstnessIterator &i) const;

      /**
       * Comparison operator alike
       * to the one above.
       */
      bool operator <= (const Iterator &i) const;

      /**
       * Same, but compare with an
       * iterator of different
       * constness.
       */
      bool operator <= (const InverseConstnessIterator &i) const;

      /**
       * Comparison operator alike
       * to the one above.
       */
      bool operator > (const Iterator &i) const;

      /**
       * Same, but compare with an
       * iterator of different
       * constness.
       */
      bool operator > (const InverseConstnessIterator &i) const;

      /**
       * Comparison operator alike
       * to the one above.
       */
      bool operator >= (const Iterator &i) const;

      /**
       * Same, but compare with an
       * iterator of different
       * constness.
       */
      bool operator >= (const InverseConstnessIterator &i) const;

      /**
       * Return the distance between
       * the two iterators, in
       * elements.
       */
      difference_type operator - (const Iterator &i) const;

      /**
       * Same, but for iterators of
       * opposite constness.
       */
      difference_type operator - (const InverseConstnessIterator &i) const;

      /**
       * Return an iterator which is
       * the given number of elements
       * in front of the present one.
       */
      Iterator operator + (const difference_type &d) const;

      /**
       * Return an iterator which is
       * the given number of elements
       * behind the present one.
       */
      Iterator operator - (const difference_type &d) const;

      /**
       * Move the iterator <tt>d</tt>
       * elements forward at once,
       * and return the result.
       */
      Iterator &operator += (const difference_type &d);

      /**
       * Move the iterator <tt>d</tt>
       * elements backward at once,
       * and return the result.
       */
      Iterator &operator -= (const difference_type &d);

      /** @addtogroup Exceptions
       * @{ */

      /**
       * Exception.
       */
      DeclException0 (ExcPointerToDifferentVectors);
      /**
       * Exception.
       */
      DeclException0 (ExcCastingAwayConstness);
      //@}
    private:
      /**
       * Pointer to the block
       * vector object to which
       * this iterator
       * points. Depending on the
       * value of the <tt>constness</tt>
       * template argument of this
       * class, this is a <tt>const</tt>
       * or non-<tt>const</tt> pointer.
       */
      BlockVector *parent;

      /**
       * Global index of the
       * element to which we
       * presently point.
       */
      size_type     global_index;

      /**
       * Current block and index
       * within this block of the
       * element presently pointed
       * to.
       */
      unsigned int current_block;
      size_type index_within_block;

      /**
       * Indices of the global
       * element address at which
       * we have to move on to
       * another block when moving
       * forward and
       * backward. These indices
       * are kept as a cache since
       * this is much more
       * efficient than always
       * asking the parent object.
       */
      size_type next_break_forward;
      size_type next_break_backward;

      /**
       * Move forward one element.
       */
      void move_forward ();

      /**
       * Move backward one element.
       */
      void move_backward ();


      /**
       * Mark all other instances of
       * this template as friends. In
       * fact, we only need the
       * inverse constness iterator
       * as friend, but this is
       * something that ISO C++ does
       * not allow to specify.
       */
      template <typename N, bool C>
      friend class Iterator;
    };
  }  // namespace BlockVectorIterators
}  // namespace internal


/**
 * A vector composed of several blocks each representing a vector of
 * its own.
 *
 * The BlockVector is a collection of Vectors (e.g. of either deal.II Vector
 * objects or PETScWrappers::Vector object). Each of the vectors inside can
 * have a different size.
 *
 * The functionality of BlockVector includes everything a Vector can do, plus
 * the access to a single Vector inside the BlockVector by
 * <tt>block(i)</tt>. It also has a complete random access iterator, just as
 * the other Vector classes or the standard C++ library template
 * <tt>std::vector</tt>. Therefore, all algorithms working on iterators also
 * work with objects of this class.
 *
 * While this base class implements most of the functionality by dispatching
 * calls to its member functions to the respective functions on each of the
 * individual blocks, this class does not actually allocate some memory or
 * change the size of vectors. For this, the constructors, assignment
 * operators and reinit() functions of derived classes are responsible. This
 * class only handles the common part that is independent of the actual vector
 * type the block vector is built on.
 *
 *
 * <h3>Accessing individual blocks, and resizing vectors</h3>
 *
 * Apart from using this object as a whole, you can use each block
 * separately as a vector, using the block() function.  There
 * is a single caveat: if you have changed the size of one or several
 * blocks, you must call the function collect_sizes() of the block
 * vector to update its internal structures.
 *
 * @attention Warning: If you change the sizes of single blocks
 * without calling collect_sizes(), results may be unpredictable. The
 * debug version does not check consistency here for performance
 * reasons!
 *
 * @see @ref GlossBlockLA "Block (linear algebra)"
 * @author Wolfgang Bangerth, Guido Kanschat, 1999, 2000, 2001, 2002, 2004
 */
template <class VectorType>
class BlockVectorBase : public Subscriptor
{
public:
  /**
   * Typedef the type of the underlying
   * vector.
   */
  typedef VectorType BlockType;

  /*
   * Declare standard types used in
   * all containers. These types
   * parallel those in the
   * <tt>C++</tt> standard
   * libraries
   * <tt>std::vector<...></tt>
   * class. This includes iterator
   * types.
   */
  typedef typename BlockType::value_type  value_type;
  typedef value_type                     *pointer;
  typedef const value_type               *const_pointer;
  typedef dealii::internal::BlockVectorIterators::Iterator<BlockVectorBase,false> iterator;
  typedef dealii::internal::BlockVectorIterators::Iterator<BlockVectorBase,true>  const_iterator;
  typedef typename BlockType::reference       reference;
  typedef typename BlockType::const_reference const_reference;
  typedef types::global_dof_index             size_type;

  /**
   * Declare a type that has holds
   * real-valued numbers with the
   * same precision as the template
   * argument to this class. If the
   * template argument of this
   * class is a real data type,
   * then real_type equals the
   * template argument. If the
   * template argument is a
   * std::complex type then
   * real_type equals the type
   * underlying the complex
   * numbers.
   *
   * This typedef is used to
   * represent the return type of
   * norms.
   */
  typedef typename BlockType::real_type real_type;

  /**
   * A variable that indicates whether this vector
   * supports distributed data storage. If true, then
   * this vector also needs an appropriate compress()
   * function that allows communicating recent set or
   * add operations to individual elements to be communicated
   * to other processors.
   *
   * For the current class, the variable equals the
   * value declared for the type of the individual blocks.
   */
  static const bool supports_distributed_data = BlockType::supports_distributed_data;

  /**
   * Default constructor.
   */
  BlockVectorBase ();

  /**
   * Update internal structures
   * after resizing
   * vectors. Whenever you reinited
   * a block of a block vector, the
   * internal data structures are
   * corrupted. Therefore, you
   * should call this function
   * after al blocks got their new
   * size.
   */
  void collect_sizes ();

  /**
   * Call the compress() function on all
   * the subblocks of the matrix.
  *
  * This functionality only needs to be
  * called if using MPI based vectors and
  * exists in other objects for
  * compatibility.
  *
  * See @ref GlossCompress "Compressing
  * distributed objects" for more
  * information.
   */
  void compress (::dealii::VectorOperation::values operation);

  /**
   * @deprecated: use compress(VectorOperation::values) instead.
   */
  void compress () DEAL_II_DEPRECATED;

  /**
   * Access to a single block.
   */
  BlockType &
  block (const unsigned int i);

  /**
   * Read-only access to a single block.
   */
  const BlockType &
  block (const unsigned int i) const;

  /**
   * Return a reference on the
   * object that describes the
   * mapping between block and
   * global indices. The use of
   * this function is highly
   * deprecated and it should
   * vanish in one of the next
   * versions
   */
  const BlockIndices &
  get_block_indices () const;

  /**
   * Number of blocks.
   */
  unsigned int n_blocks () const;

  /**
   * Return dimension of the vector. This
   * is the sum of the dimensions of all
   * components.
   */
  std::size_t size () const;

  /**
   * Return an index set that describes which elements of this vector
   * are owned by the current processor. Note that this index set does
   * not include elements this vector may store locally as ghost
   * elements but that are in fact owned by another processor.
   * As a consequence, the index sets returned on different
   * processors if this is a distributed vector will form disjoint
   * sets that add up to the complete index set.
   * Obviously, if a vector is created on only one processor, then
   * the result would satisfy
   * @code
   *   vec.locally_owned_elements() == complete_index_set (vec.size())
   * @endcode
   *
   * For block vectors, this function returns the union of the
   * locally owned elements of the individual blocks, shifted by
   * their respective index offsets.
   */
  IndexSet locally_owned_elements () const;

  /**
   * Return an iterator pointing to
   * the first element.
   */
  iterator begin ();

  /**
   * Return an iterator pointing to
   * the first element of a
   * constant block vector.
   */
  const_iterator begin () const;

  /**
   * Return an iterator pointing to
   * the element past the end.
   */
  iterator end ();

  /**
   * Return an iterator pointing to
   * the element past the end of a
   * constant block vector.
   */
  const_iterator end () const;

  /**
   * Access components, returns U(i).
   */
  value_type operator() (const size_type i) const;

  /**
   * Access components, returns U(i)
   * as a writeable reference.
   */
  reference operator() (const size_type i);

  /**
   * Access components, returns U(i).
   *
   * Exactly the same as operator().
   */
  value_type operator[] (const size_type i) const;

  /**
   * Access components, returns U(i)
   * as a writeable reference.
   *
   * Exactly the same as operator().
   */
  reference operator[] (const size_type i);

  /**
   * A collective get operation: instead
   * of getting individual elements of a
   * vector, this function allows to get
   * a whole set of elements at once. The
   * indices of the elements to be read
   * are stated in the first argument,
   * the corresponding values are returned in the
   * second.
   */
  template <typename OtherNumber>
  void extract_subvector_to (const std::vector<size_type> &indices,
                             std::vector<OtherNumber> &values) const;

  /**
   * Just as the above, but with pointers.
   * Useful in minimizing copying of data around.
   */
  template <typename ForwardIterator, typename OutputIterator>
  void extract_subvector_to (ForwardIterator          indices_begin,
                             const ForwardIterator    indices_end,
                             OutputIterator           values_begin) const;

  /**
   * Copy operator: fill all components of
   * the vector with the given scalar
   * value.
   */
  BlockVectorBase &operator = (const value_type s);

  /**
   * Copy operator for arguments of the
   * same type.
   */
  BlockVectorBase &
  operator= (const BlockVectorBase &V);

  /**
   * Copy operator for template arguments
   * of different types.
   */
  template <class VectorType2>
  BlockVectorBase &
  operator= (const BlockVectorBase<VectorType2> &V);

  /**
   * Copy operator from non-block
   * vectors to block vectors.
   */
  BlockVectorBase &
  operator = (const VectorType &v);

  /**
   * Check for equality of two block vector
   * types. This operation is only allowed
   * if the two vectors already have the
   * same block structure.
   */
  template <class VectorType2>
  bool
  operator == (const BlockVectorBase<VectorType2> &v) const;

  /**
   * $U = U * V$: scalar product.
   */
  value_type operator* (const BlockVectorBase &V) const;

  /**
   * Return square of the $l_2$-norm.
   */
  real_type norm_sqr () const;

  /**
   * Return the mean value of the elements
   * of this vector.
   */
  value_type mean_value () const;

  /**
   * Return the $l_1$-norm of the vector,
   * i.e. the sum of the absolute values.
   */
  real_type l1_norm () const;

  /**
   * Return the $l_2$-norm of the vector,
   * i.e. the square root of the sum of
   * the squares of the elements.
   */
  real_type l2_norm () const;

  /**
   * Return the maximum absolute value of
   * the elements of this vector, which is
   * the $l_\infty$-norm of a vector.
   */
  real_type linfty_norm () const;

  /**
   * Returns true if the given global index is
   * in the local range of this processor.
   * Asks the corresponding block.
   */
  bool in_local_range (const size_type global_index) const;

  /**
   * Return whether the vector contains only
   * elements with value zero. This function
   * is mainly for internal consistency
   * check and should seldom be used when
   * not in debug mode since it uses quite
   * some time.
   */
  bool all_zero () const;

  /**
   * Return @p true if the vector has no
   * negative entries, i.e. all entries are
   * zero or positive. This function is
   * used, for example, to check whether
   * refinement indicators are really all
   * positive (or zero).
   */
  bool is_non_negative () const;

  /**
   * Addition operator.  Fast equivalent to
   * <tt>U.add(1, V)</tt>.
   */
  BlockVectorBase &
  operator += (const BlockVectorBase &V);

  /**
   * Subtraction operator.  Fast equivalent
   * to <tt>U.add(-1, V)</tt>.
   */
  BlockVectorBase &
  operator -= (const BlockVectorBase &V);


  /**
   * A collective add operation:
   * This funnction adds a whole
   * set of values stored in @p
   * values to the vector
   * components specified by @p
   * indices.
   */
  template <typename Number>
  void add (const std::vector<size_type> &indices,
            const std::vector<Number>    &values);

  /**
   * This is a second collective
   * add operation. As a
   * difference, this function
   * takes a deal.II vector of
   * values.
   */
  template <typename Number>
  void add (const std::vector<size_type> &indices,
            const Vector<Number>         &values);

  /**
   * Take an address where
   * <tt>n_elements</tt> are stored
   * contiguously and add them into
   * the vector. Handles all cases
   * which are not covered by the
   * other two <tt>add()</tt>
   * functions above.
   */
  template <typename Number>
  void add (const size_type  n_elements,
            const size_type *indices,
            const Number    *values);

  /**
   * $U(0-DIM)+=s$.  Addition of <tt>s</tt>
   * to all components. Note that
   * <tt>s</tt> is a scalar and not a
   * vector.
   */
  void add (const value_type s);

  /**
   * U+=V.
   * Simple vector addition, equal to the
   * <tt>operator +=</tt>.
   */
  void add (const BlockVectorBase &V);

  /**
   * U+=a*V.
   * Simple addition of a scaled vector.
   */
  void add (const value_type a, const BlockVectorBase &V);

  /**
   * U+=a*V+b*W.
   * Multiple addition of scaled vectors.
   */
  void add (const value_type a, const BlockVectorBase &V,
            const value_type b, const BlockVectorBase &W);

  /**
   * U=s*U+V.
   * Scaling and simple vector addition.
   */
  void sadd (const value_type s, const BlockVectorBase &V);

  /**
   * U=s*U+a*V.
   * Scaling and simple addition.
   */
  void sadd (const value_type s, const value_type a, const BlockVectorBase &V);

  /**
   * U=s*U+a*V+b*W.
   * Scaling and multiple addition.
   */
  void sadd (const value_type s, const value_type a,
             const BlockVectorBase &V,
             const value_type b, const BlockVectorBase &W);

  /**
   * U=s*U+a*V+b*W+c*X.
   * Scaling and multiple addition.
   */
  void sadd (const value_type s, const value_type a,
             const BlockVectorBase &V,
             const value_type b, const BlockVectorBase &W,
             const value_type c, const BlockVectorBase &X);

  /**
   * Scale each element of the
   * vector by a constant
   * value.
   */
  BlockVectorBase &operator *= (const value_type factor);

  /**
   * Scale each element of the
   * vector by the inverse of the
   * given value.
   */
  BlockVectorBase &operator /= (const value_type factor);

  /**
   * Multiply each element of this
   * vector by the corresponding
   * element of <tt>v</tt>.
   */
  template <class BlockVector2>
  void scale (const BlockVector2 &v);

  /**
   *  U=a*V. Assignment.
   */
  template <class BlockVector2>
  void equ (const value_type a, const BlockVector2 &V);

  /**
   * U=a*V+b*W.
   * Replacing by sum.
   */
  void equ (const value_type a, const BlockVectorBase &V,
            const value_type b, const BlockVectorBase &W);

  /**
   * This function does nothing but is
   * there for compatibility with the
   * @p PETScWrappers::Vector class.
   *
   * For the PETSc vector wrapper class,
   * this function updates the ghost
   * values of the PETSc vector. This
   * is necessary after any modification
   * before reading ghost values.
   *
   * However, for the implementation of
   * this class, it is immaterial and thus
   * an empty function.
   */
  void update_ghost_values () const;

  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   */
  std::size_t memory_consumption () const;

protected:
  /**
   * Pointer to the array of components.
   */
  std::vector<VectorType> components;

  /**
   * Object managing the
   * transformation between global
   * indices and indices within the
   * different blocks.
   */
  BlockIndices block_indices;

  /**
   * Make the iterator class a
   * friend.
   */
  template <typename N, bool C>
  friend class dealii::internal::BlockVectorIterators::Iterator;

  template <typename> friend class BlockVectorBase;
};


/*@}*/

/*----------------------- Inline functions ----------------------------------*/


#ifndef DOXYGEN
namespace internal
{
  namespace BlockVectorIterators
  {

    template <class BlockVectorType, bool constness>
    inline
    Iterator<BlockVectorType,constness>::
    Iterator (const Iterator<BlockVectorType,constness> &c)
      :
      parent (c.parent),
      global_index (c.global_index),
      current_block (c.current_block),
      index_within_block (c.index_within_block),
      next_break_forward (c.next_break_forward),
      next_break_backward (c.next_break_backward)
    {}



    template <class BlockVectorType, bool constness>
    inline
    Iterator<BlockVectorType,constness>::
    Iterator (const InverseConstnessIterator &c)
      :
      parent (const_cast<BlockVectorType *>(c.parent)),
      global_index (c.global_index),
      current_block (c.current_block),
      index_within_block (c.index_within_block),
      next_break_forward (c.next_break_forward),
      next_break_backward (c.next_break_backward)
    {
      // if constness==false, then the
      // constness of the iterator we
      // got is true and we are trying
      // to cast away the
      // constness. disallow this
      Assert (constness==true, ExcCastingAwayConstness());
    }



    template <class BlockVectorType, bool constness>
    inline
    Iterator<BlockVectorType,constness>::
    Iterator (BlockVector        &parent,
              const size_type  global_index,
              const size_type  current_block,
              const size_type  index_within_block,
              const size_type  next_break_forward,
              const size_type  next_break_backward)
      :
      parent (&parent),
      global_index (global_index),
      current_block (current_block),
      index_within_block (index_within_block),
      next_break_forward (next_break_forward),
      next_break_backward (next_break_backward)
    {
    }



    template <class BlockVectorType, bool constness>
    inline
    Iterator<BlockVectorType,constness> &
    Iterator<BlockVectorType,constness>::
    operator = (const Iterator &c)
    {
      parent              = c.parent;
      global_index        = c.global_index;
      index_within_block  = c.index_within_block;
      current_block       = c.current_block;
      next_break_forward  = c.next_break_forward;
      next_break_backward = c.next_break_backward;

      return *this;
    }



    template <class BlockVectorType, bool constness>
    inline
    typename Iterator<BlockVectorType,constness>::dereference_type
    Iterator<BlockVectorType,constness>::operator * () const
    {
      return parent->block(current_block)(index_within_block);
    }



    template <class BlockVectorType, bool constness>
    inline
    typename Iterator<BlockVectorType,constness>::dereference_type
    Iterator<BlockVectorType,constness>::operator [] (const difference_type d) const
    {
      // if the index pointed to is
      // still within the block we
      // currently point into, then we
      // can save the computation of
      // the block
      if ((global_index+d >= next_break_backward) &&
          (global_index+d <= next_break_forward))
        return parent->block(current_block)(index_within_block + d);

      // if the index is not within the
      // block of the block vector into
      // which we presently point, then
      // there is no way: we have to
      // search for the block. this can
      // be done through the parent
      // class as well.
      return (*parent)(global_index+d);
    }



    template <class BlockVectorType, bool constness>
    inline
    Iterator<BlockVectorType,constness> &
    Iterator<BlockVectorType,constness>::operator ++ ()
    {
      move_forward ();
      return *this;
    }



    template <class BlockVectorType, bool constness>
    inline
    Iterator<BlockVectorType,constness>
    Iterator<BlockVectorType,constness>::operator ++ (int)
    {
      const Iterator old_value = *this;
      move_forward ();
      return old_value;
    }



    template <class BlockVectorType, bool constness>
    inline
    Iterator<BlockVectorType,constness> &
    Iterator<BlockVectorType,constness>::operator -- ()
    {
      move_backward ();
      return *this;
    }



    template <class BlockVectorType, bool constness>
    inline
    Iterator<BlockVectorType,constness>
    Iterator<BlockVectorType,constness>::operator -- (int)
    {
      const Iterator old_value = *this;
      move_backward ();
      return old_value;
    }



    template <class BlockVectorType, bool constness>
    inline
    bool
    Iterator<BlockVectorType,constness>::
    operator == (const Iterator &i) const
    {
      Assert (parent == i.parent, ExcPointerToDifferentVectors());

      return (global_index == i.global_index);
    }



    template <class BlockVectorType, bool constness>
    inline
    bool
    Iterator<BlockVectorType,constness>::
    operator == (const InverseConstnessIterator &i) const
    {
      Assert (parent == i.parent, ExcPointerToDifferentVectors());

      return (global_index == i.global_index);
    }



    template <class BlockVectorType, bool constness>
    inline
    bool
    Iterator<BlockVectorType,constness>::
    operator != (const Iterator &i) const
    {
      Assert (parent == i.parent, ExcPointerToDifferentVectors());

      return (global_index != i.global_index);
    }



    template <class BlockVectorType, bool constness>
    inline
    bool
    Iterator<BlockVectorType,constness>::
    operator != (const InverseConstnessIterator &i) const
    {
      Assert (parent == i.parent, ExcPointerToDifferentVectors());

      return (global_index != i.global_index);
    }



    template <class BlockVectorType, bool constness>
    inline
    bool
    Iterator<BlockVectorType,constness>::
    operator < (const Iterator &i) const
    {
      Assert (parent == i.parent, ExcPointerToDifferentVectors());

      return (global_index < i.global_index);
    }



    template <class BlockVectorType, bool constness>
    inline
    bool
    Iterator<BlockVectorType,constness>::
    operator < (const InverseConstnessIterator &i) const
    {
      Assert (parent == i.parent, ExcPointerToDifferentVectors());

      return (global_index < i.global_index);
    }



    template <class BlockVectorType, bool constness>
    inline
    bool
    Iterator<BlockVectorType,constness>::
    operator <= (const Iterator &i) const
    {
      Assert (parent == i.parent, ExcPointerToDifferentVectors());

      return (global_index <= i.global_index);
    }



    template <class BlockVectorType, bool constness>
    inline
    bool
    Iterator<BlockVectorType,constness>::
    operator <= (const InverseConstnessIterator &i) const
    {
      Assert (parent == i.parent, ExcPointerToDifferentVectors());

      return (global_index <= i.global_index);
    }



    template <class BlockVectorType, bool constness>
    inline
    bool
    Iterator<BlockVectorType,constness>::
    operator > (const Iterator &i) const
    {
      Assert (parent == i.parent, ExcPointerToDifferentVectors());

      return (global_index > i.global_index);
    }



    template <class BlockVectorType, bool constness>
    inline
    bool
    Iterator<BlockVectorType,constness>::
    operator > (const InverseConstnessIterator &i) const
    {
      Assert (parent == i.parent, ExcPointerToDifferentVectors());

      return (global_index > i.global_index);
    }



    template <class BlockVectorType, bool constness>
    inline
    bool
    Iterator<BlockVectorType,constness>::
    operator >= (const Iterator &i) const
    {
      Assert (parent == i.parent, ExcPointerToDifferentVectors());

      return (global_index >= i.global_index);
    }



    template <class BlockVectorType, bool constness>
    inline
    bool
    Iterator<BlockVectorType,constness>::
    operator >= (const InverseConstnessIterator &i) const
    {
      Assert (parent == i.parent, ExcPointerToDifferentVectors());

      return (global_index >= i.global_index);
    }



    template <class BlockVectorType, bool constness>
    inline
    typename Iterator<BlockVectorType,constness>::difference_type
    Iterator<BlockVectorType,constness>::
    operator - (const Iterator &i) const
    {
      Assert (parent == i.parent, ExcPointerToDifferentVectors());

      return (static_cast<signed int>(global_index) -
              static_cast<signed int>(i.global_index));
    }



    template <class BlockVectorType, bool constness>
    inline
    typename Iterator<BlockVectorType,constness>::difference_type
    Iterator<BlockVectorType,constness>::
    operator - (const InverseConstnessIterator &i) const
    {
      Assert (parent == i.parent, ExcPointerToDifferentVectors());

      return (static_cast<signed int>(global_index) -
              static_cast<signed int>(i.global_index));
    }



    template <class BlockVectorType, bool constness>
    inline
    Iterator<BlockVectorType,constness>
    Iterator<BlockVectorType,constness>::
    operator + (const difference_type &d) const
    {
      // if the index pointed to is
      // still within the block we
      // currently point into, then we
      // can save the computation of
      // the block
      if ((global_index+d >= next_break_backward) &&
          (global_index+d <= next_break_forward))
        return Iterator (*parent, global_index+d, current_block,
                         index_within_block+d,
                         next_break_forward, next_break_backward);
      else
        // outside present block, so
        // have to seek new block
        // anyway
        return Iterator (*parent, global_index+d);
    }



    template <class BlockVectorType, bool constness>
    inline
    Iterator<BlockVectorType,constness>
    Iterator<BlockVectorType,constness>::
    operator - (const difference_type &d) const
    {
      // if the index pointed to is
      // still within the block we
      // currently point into, then we
      // can save the computation of
      // the block
      if ((global_index-d >= next_break_backward) &&
          (global_index-d <= next_break_forward))
        return Iterator (*parent, global_index-d, current_block,
                         index_within_block-d,
                         next_break_forward, next_break_backward);
      else
        // outside present block, so
        // have to seek new block
        // anyway
        return Iterator (*parent, global_index-d);
    }



    template <class BlockVectorType, bool constness>
    inline
    Iterator<BlockVectorType,constness> &
    Iterator<BlockVectorType,constness>::
    operator += (const difference_type &d)
    {
      // if the index pointed to is
      // still within the block we
      // currently point into, then we
      // can save the computation of
      // the block
      if ((global_index+d >= next_break_backward) &&
          (global_index+d <= next_break_forward))
        {
          global_index       += d;
          index_within_block += d;
        }
      else
        // outside present block, so
        // have to seek new block
        // anyway
        *this = Iterator (*parent, global_index+d);

      return *this;
    }



    template <class BlockVectorType, bool constness>
    inline
    Iterator<BlockVectorType,constness> &
    Iterator<BlockVectorType,constness>::
    operator -= (const difference_type &d)
    {
      // if the index pointed to is
      // still within the block we
      // currently point into, then we
      // can save the computation of
      // the block
      if ((global_index-d >= next_break_backward) &&
          (global_index-d <= next_break_forward))
        {
          global_index       -= d;
          index_within_block -= d;
        }
      else
        // outside present block, so
        // have to seek new block
        // anyway
        *this = Iterator (*parent, global_index-d);

      return *this;
    }


    template <class BlockVectorType, bool constness>
    Iterator<BlockVectorType,constness>::
    Iterator (BlockVector    &parent,
              const size_type global_index)
      :
      parent (&parent),
      global_index (global_index)
    {
      // find which block we are
      // in. for this, take into
      // account that it happens at
      // times that people want to
      // initialize iterators
      // past-the-end
      if (global_index < parent.size())
        {
          const std::pair<size_type, size_type>
          indices = parent.block_indices.global_to_local(global_index);
          current_block      = indices.first;
          index_within_block = indices.second;

          next_break_backward
            = parent.block_indices.local_to_global (current_block, 0);
          next_break_forward
            = (parent.block_indices.local_to_global (current_block, 0)
               +parent.block_indices.block_size(current_block)-1);
        }
      else
        // past the end. only have one
        // value for this
        {
          this->global_index  = parent.size ();
          current_block       = parent.n_blocks();
          index_within_block  = 0;
          next_break_backward = global_index;
          next_break_forward  = numbers::invalid_size_type;
        };
    }



    template <class BlockVectorType, bool constness>
    void
    Iterator<BlockVectorType,constness>::move_forward ()
    {
      if (global_index != next_break_forward)
        ++index_within_block;
      else
        {
          // ok, we traverse a boundary
          // between blocks:
          index_within_block = 0;
          ++current_block;

          // break backwards is now old
          // break forward
          next_break_backward = next_break_forward+1;

          // compute new break forward
          if (current_block < parent->block_indices.size())
            next_break_forward
            += parent->block_indices.block_size(current_block);
          else
            // if we are beyond the end,
            // then move the next
            // boundary arbitrarily far
            // away
            next_break_forward = numbers::invalid_size_type;
        };

      ++global_index;
    }



    template <class BlockVectorType, bool constness>
    void
    Iterator<BlockVectorType,constness>::move_backward ()
    {
      if (global_index != next_break_backward)
        --index_within_block;
      else if (current_block != 0)
        {
          // ok, we traverse a boundary
          // between blocks:
          --current_block;
          index_within_block = parent->block_indices.block_size(current_block)-1;

          // break forwards is now old
          // break backward
          next_break_forward = next_break_backward-1;

          // compute new break forward
          next_break_backward
          -= parent->block_indices.block_size (current_block);
        }
      else
        // current block was 0, we now
        // get into unspecified terrain
        {
          --current_block;
          index_within_block = numbers::invalid_size_type;
          next_break_forward = 0;
          next_break_backward = 0;
        };

      --global_index;
    }


  } // namespace BlockVectorIterators

} //namespace internal


template <class VectorType>
inline
BlockVectorBase<VectorType>::BlockVectorBase ()
{}



template <class VectorType>
inline
std::size_t
BlockVectorBase<VectorType>::size () const
{
  return block_indices.total_size();
}



template <class VectorType>
inline
IndexSet
BlockVectorBase<VectorType>::locally_owned_elements () const
{
  IndexSet is (size());

  // copy index sets from blocks into the global one, shifted
  // by the appropriate amount for each block
  for (unsigned int b=0; b<n_blocks(); ++b)
    {
      IndexSet x = block(b).locally_owned_elements();
      is.add_indices(x, block_indices.block_start(b));
    }

  is.compress();

  return is;
}



template <class VectorType>
inline
unsigned int
BlockVectorBase<VectorType>::n_blocks () const
{
  return block_indices.size();
}


template <class VectorType>
inline
typename BlockVectorBase<VectorType>::BlockType &
BlockVectorBase<VectorType>::block (const unsigned int i)
{
  Assert(i<n_blocks(), ExcIndexRange(i,0,n_blocks()));

  return components[i];
}



template <class VectorType>
inline
const typename BlockVectorBase<VectorType>::BlockType &
BlockVectorBase<VectorType>::block (const unsigned int i) const
{
  Assert(i<n_blocks(), ExcIndexRange(i,0,n_blocks()));

  return components[i];
}



template <class VectorType>
inline
const BlockIndices &
BlockVectorBase<VectorType>::get_block_indices () const
{
  return block_indices;
}


template <class VectorType>
inline
void
BlockVectorBase<VectorType>::collect_sizes ()
{
  std::vector<size_type> sizes (n_blocks());

  for (size_type i=0; i<n_blocks(); ++i)
    sizes[i] = block(i).size();

  block_indices.reinit(sizes);
}



template <class VectorType>
inline
void
BlockVectorBase<VectorType>::compress (::dealii::VectorOperation::values operation)
{
  for (unsigned int i=0; i<n_blocks(); ++i)
    block(i).compress (operation);
}



template <class VectorType>
inline
void
BlockVectorBase<VectorType>::compress ()
{
  compress(VectorOperation::unknown);
}



template <class VectorType>
inline
typename BlockVectorBase<VectorType>::iterator
BlockVectorBase<VectorType>::begin()
{
  return iterator(*this, 0U);
}



template <class VectorType>
inline
typename BlockVectorBase<VectorType>::const_iterator
BlockVectorBase<VectorType>::begin() const
{
  return const_iterator(*this, 0U);
}


template <class VectorType>
inline
typename BlockVectorBase<VectorType>::iterator
BlockVectorBase<VectorType>::end()
{
  return iterator(*this, size());
}



template <class VectorType>
inline
typename BlockVectorBase<VectorType>::const_iterator
BlockVectorBase<VectorType>::end() const
{
  return const_iterator(*this, size());
}


template <class VectorType>
inline
bool
BlockVectorBase<VectorType>::in_local_range
(const size_type global_index) const
{
  const std::pair<size_type,size_type> local_index
    = block_indices.global_to_local (global_index);

  return components[local_index.first].in_local_range (global_index);
}


template <class VectorType>
bool
BlockVectorBase<VectorType>::all_zero () const
{
  for (size_type i=0; i<n_blocks(); ++i)
    if (components[i].all_zero() == false)
      return false;

  return true;
}



template <class VectorType>
bool
BlockVectorBase<VectorType>::is_non_negative () const
{
  for (size_type i=0; i<n_blocks(); ++i)
    if (components[i].is_non_negative() == false)
      return false;

  return true;
}



template <class VectorType>
typename BlockVectorBase<VectorType>::value_type
BlockVectorBase<VectorType>::
operator * (const BlockVectorBase<VectorType> &v) const
{
  Assert (n_blocks() == v.n_blocks(),
          ExcDimensionMismatch(n_blocks(), v.n_blocks()));

  value_type sum = 0.;
  for (size_type i=0; i<n_blocks(); ++i)
    sum += components[i]*v.components[i];

  return sum;
}


template <class VectorType>
typename BlockVectorBase<VectorType>::real_type
BlockVectorBase<VectorType>::norm_sqr () const
{
  real_type sum = 0.;
  for (size_type i=0; i<n_blocks(); ++i)
    sum += components[i].norm_sqr();

  return sum;
}



template <class VectorType>
typename BlockVectorBase<VectorType>::value_type
BlockVectorBase<VectorType>::mean_value () const
{
  value_type sum = 0.;
  for (size_type i=0; i<n_blocks(); ++i)
    sum += components[i].mean_value() * components[i].size();

  return sum/size();
}



template <class VectorType>
typename BlockVectorBase<VectorType>::real_type
BlockVectorBase<VectorType>::l1_norm () const
{
  real_type sum = 0.;
  for (size_type i=0; i<n_blocks(); ++i)
    sum += components[i].l1_norm();

  return sum;
}



template <class VectorType>
typename BlockVectorBase<VectorType>::real_type
BlockVectorBase<VectorType>::l2_norm () const
{
  return std::sqrt(norm_sqr());
}



template <class VectorType>
typename BlockVectorBase<VectorType>::real_type
BlockVectorBase<VectorType>::linfty_norm () const
{
  real_type sum = 0.;
  for (size_type i=0; i<n_blocks(); ++i)
    {
      value_type newval = components[i].linfty_norm();
      if (sum<newval)
        sum = newval;
    }
  return sum;
}



template <class VectorType>
BlockVectorBase<VectorType> &
BlockVectorBase<VectorType>::operator += (const BlockVectorBase<VectorType> &v)
{
  add (v);
  return *this;
}



template <class VectorType>
BlockVectorBase<VectorType> &
BlockVectorBase<VectorType>::operator -= (const BlockVectorBase<VectorType> &v)
{
  Assert (n_blocks() == v.n_blocks(),
          ExcDimensionMismatch(n_blocks(), v.n_blocks()));

  for (size_type i=0; i<n_blocks(); ++i)
    {
      components[i] -= v.components[i];
    }
  return *this;
}



template <class VectorType>
template <typename Number>
inline
void
BlockVectorBase<VectorType>::add (const std::vector<size_type> &indices,
                                  const std::vector<Number>    &values)
{
  Assert (indices.size() == values.size(),
          ExcDimensionMismatch(indices.size(), values.size()));
  add (indices.size(), &indices[0], &values[0]);
}



template <class VectorType>
template <typename Number>
inline
void
BlockVectorBase<VectorType>::add (const std::vector<size_type> &indices,
                                  const Vector<Number>         &values)
{
  Assert (indices.size() == values.size(),
          ExcDimensionMismatch(indices.size(), values.size()));
  const size_type n_indices = indices.size();
  for (size_type i=0; i<n_indices; ++i)
    (*this)(indices[i]) += values(i);
}



template <class VectorType>
template <typename Number>
inline
void
BlockVectorBase<VectorType>::add (const size_type  n_indices,
                                  const size_type *indices,
                                  const Number    *values)
{
  for (size_type i=0; i<n_indices; ++i)
    (*this)(indices[i]) += values[i];
}



template <class VectorType>
void BlockVectorBase<VectorType>::add (const value_type a)
{
  Assert (numbers::is_finite(a), ExcNumberNotFinite());

  for (size_type i=0; i<n_blocks(); ++i)
    {
      components[i].add(a);
    }
}



template <class VectorType>
void BlockVectorBase<VectorType>::add (const BlockVectorBase<VectorType> &v)
{
  Assert (n_blocks() == v.n_blocks(),
          ExcDimensionMismatch(n_blocks(), v.n_blocks()));

  for (size_type i=0; i<n_blocks(); ++i)
    {
      components[i].add(v.components[i]);
    }
}



template <class VectorType>
void BlockVectorBase<VectorType>::add (const value_type a,
                                       const BlockVectorBase<VectorType> &v)
{

  Assert (numbers::is_finite(a), ExcNumberNotFinite());

  Assert (n_blocks() == v.n_blocks(),
          ExcDimensionMismatch(n_blocks(), v.n_blocks()));

  for (size_type i=0; i<n_blocks(); ++i)
    {
      components[i].add(a, v.components[i]);
    }
}



template <class VectorType>
void BlockVectorBase<VectorType>::add (const value_type a,
                                       const BlockVectorBase<VectorType> &v,
                                       const value_type b,
                                       const BlockVectorBase<VectorType> &w)
{

  Assert (numbers::is_finite(a), ExcNumberNotFinite());
  Assert (numbers::is_finite(b), ExcNumberNotFinite());

  Assert (n_blocks() == v.n_blocks(),
          ExcDimensionMismatch(n_blocks(), v.n_blocks()));
  Assert (n_blocks() == w.n_blocks(),
          ExcDimensionMismatch(n_blocks(), w.n_blocks()));


  for (size_type i=0; i<n_blocks(); ++i)
    {
      components[i].add(a, v.components[i], b, w.components[i]);
    }
}



template <class VectorType>
void BlockVectorBase<VectorType>::sadd (const value_type x,
                                        const BlockVectorBase<VectorType> &v)
{

  Assert (numbers::is_finite(x), ExcNumberNotFinite());

  Assert (n_blocks() == v.n_blocks(),
          ExcDimensionMismatch(n_blocks(), v.n_blocks()));

  for (size_type i=0; i<n_blocks(); ++i)
    {
      components[i].sadd(x, v.components[i]);
    }
}



template <class VectorType>
void BlockVectorBase<VectorType>::sadd (const value_type x, const value_type a,
                                        const BlockVectorBase<VectorType> &v)
{

  Assert (numbers::is_finite(x), ExcNumberNotFinite());
  Assert (numbers::is_finite(a), ExcNumberNotFinite());

  Assert (n_blocks() == v.n_blocks(),
          ExcDimensionMismatch(n_blocks(), v.n_blocks()));

  for (size_type i=0; i<n_blocks(); ++i)
    {
      components[i].sadd(x, a, v.components[i]);
    }
}



template <class VectorType>
void BlockVectorBase<VectorType>::sadd (const value_type x, const value_type a,
                                        const BlockVectorBase<VectorType> &v,
                                        const value_type b,
                                        const BlockVectorBase<VectorType> &w)
{

  Assert (numbers::is_finite(x), ExcNumberNotFinite());
  Assert (numbers::is_finite(a), ExcNumberNotFinite());
  Assert (numbers::is_finite(b), ExcNumberNotFinite());

  Assert (n_blocks() == v.n_blocks(),
          ExcDimensionMismatch(n_blocks(), v.n_blocks()));
  Assert (n_blocks() == w.n_blocks(),
          ExcDimensionMismatch(n_blocks(), w.n_blocks()));

  for (size_type i=0; i<n_blocks(); ++i)
    {
      components[i].sadd(x, a, v.components[i], b, w.components[i]);
    }
}



template <class VectorType>
void BlockVectorBase<VectorType>::sadd (const value_type x, const value_type a,
                                        const BlockVectorBase<VectorType> &v,
                                        const value_type b,
                                        const BlockVectorBase<VectorType> &w,
                                        const value_type c,
                                        const BlockVectorBase<VectorType> &y)
{

  Assert (numbers::is_finite(x), ExcNumberNotFinite());
  Assert (numbers::is_finite(a), ExcNumberNotFinite());
  Assert (numbers::is_finite(b), ExcNumberNotFinite());
  Assert (numbers::is_finite(c), ExcNumberNotFinite());

  Assert (n_blocks() == v.n_blocks(),
          ExcDimensionMismatch(n_blocks(), v.n_blocks()));
  Assert (n_blocks() == w.n_blocks(),
          ExcDimensionMismatch(n_blocks(), w.n_blocks()));
  Assert (n_blocks() == y.n_blocks(),
          ExcDimensionMismatch(n_blocks(), y.n_blocks()));

  for (size_type i=0; i<n_blocks(); ++i)
    {
      components[i].sadd(x, a, v.components[i],
                         b, w.components[i], c, y.components[i]);
    }
}



template <class VectorType>
template <class BlockVector2>
void BlockVectorBase<VectorType>::scale (const BlockVector2 &v)
{
  Assert (n_blocks() == v.n_blocks(),
          ExcDimensionMismatch(n_blocks(), v.n_blocks()));
  for (size_type i=0; i<n_blocks(); ++i)
    components[i].scale(v.block(i));
}



template <class VectorType>
void BlockVectorBase<VectorType>::equ (const value_type a,
                                       const BlockVectorBase<VectorType> &v,
                                       const value_type b,
                                       const BlockVectorBase<VectorType> &w)
{

  Assert (numbers::is_finite(a), ExcNumberNotFinite());
  Assert (numbers::is_finite(b), ExcNumberNotFinite());

  Assert (n_blocks() == v.n_blocks(),
          ExcDimensionMismatch(n_blocks(), v.n_blocks()));
  Assert (n_blocks() == w.n_blocks(),
          ExcDimensionMismatch(n_blocks(), w.n_blocks()));

  for (size_type i=0; i<n_blocks(); ++i)
    {
      components[i].equ( a, v.components[i], b, w.components[i]);
    }
}



template <class VectorType>
std::size_t
BlockVectorBase<VectorType>::memory_consumption () const
{
  std::size_t mem = sizeof(this->n_blocks());
  for (size_type i=0; i<this->components.size(); ++i)
    mem += MemoryConsumption::memory_consumption (this->components[i]);
  mem += MemoryConsumption::memory_consumption (this->block_indices);
  return mem;
}



template <class VectorType>
template <class BlockVector2>
void BlockVectorBase<VectorType>::equ (const value_type    a,
                                       const BlockVector2 &v)
{

  Assert (numbers::is_finite(a), ExcNumberNotFinite());

  Assert (n_blocks() == v.n_blocks(),
          ExcDimensionMismatch(n_blocks(), v.n_blocks()));

  for (size_type i=0; i<n_blocks(); ++i)
    components[i].equ( a, v.components[i]);
}



template <class VectorType>
void BlockVectorBase<VectorType>::update_ghost_values () const
{
  for (size_type i=0; i<n_blocks(); ++i)
    block(i).update_ghost_values ();
}



template <class VectorType>
BlockVectorBase<VectorType> &
BlockVectorBase<VectorType>::operator = (const value_type s)
{

  Assert (numbers::is_finite(s), ExcNumberNotFinite());

  for (size_type i=0; i<n_blocks(); ++i)
    components[i] = s;

  return *this;
}


template <class VectorType>
BlockVectorBase<VectorType> &
BlockVectorBase<VectorType>::operator = (const BlockVectorBase<VectorType> &v)
{
  AssertDimension(n_blocks(), v.n_blocks());

  for (size_type i=0; i<n_blocks(); ++i)
    components[i] = v.components[i];

  return *this;
}


template <class VectorType>
template <class VectorType2>
BlockVectorBase<VectorType> &
BlockVectorBase<VectorType>::operator = (const BlockVectorBase<VectorType2> &v)
{
  AssertDimension(n_blocks(), v.n_blocks());

  for (size_type i=0; i<n_blocks(); ++i)
    components[i] = v.components[i];

  return *this;
}



template <class VectorType>
BlockVectorBase<VectorType> &
BlockVectorBase<VectorType>::operator = (const VectorType &v)
{
  Assert (size() == v.size(),
          ExcDimensionMismatch(size(), v.size()));

  size_type index_v = 0;
  for (size_type b=0; b<n_blocks(); ++b)
    for (size_type i=0; i<block(b).size(); ++i, ++index_v)
      block(b)(i) = v(index_v);

  return *this;
}



template <class VectorType>
template <class VectorType2>
inline
bool
BlockVectorBase<VectorType>::
operator == (const BlockVectorBase<VectorType2> &v) const
{
  Assert (block_indices == v.block_indices, ExcDifferentBlockIndices());

  for (size_type i=0; i<n_blocks(); ++i)
    if ( ! (components[i] == v.components[i]))
      return false;

  return true;
}



template <class VectorType>
inline
BlockVectorBase<VectorType> &
BlockVectorBase<VectorType>::operator *= (const value_type factor)
{

  Assert (numbers::is_finite(factor), ExcNumberNotFinite());

  for (size_type i=0; i<n_blocks(); ++i)
    components[i] *= factor;

  return *this;
}



template <class VectorType>
inline
BlockVectorBase<VectorType> &
BlockVectorBase<VectorType>::operator /= (const value_type factor)
{

  Assert (numbers::is_finite(factor), ExcNumberNotFinite());
  Assert (factor > 0., ExcDivideByZero() );

  for (size_type i=0; i<n_blocks(); ++i)
    components[i] /= factor;

  return *this;
}


template <class VectorType>
inline
typename BlockVectorBase<VectorType>::value_type
BlockVectorBase<VectorType>::operator() (const size_type i) const
{
  const std::pair<size_type,size_type> local_index
    = block_indices.global_to_local (i);
  return components[local_index.first](local_index.second);
}



template <class VectorType>
inline
typename BlockVectorBase<VectorType>::reference
BlockVectorBase<VectorType>::operator() (const size_type i)
{
  const std::pair<size_type,size_type> local_index
    = block_indices.global_to_local (i);
  return components[local_index.first](local_index.second);
}



template <class VectorType>
inline
typename BlockVectorBase<VectorType>::value_type
BlockVectorBase<VectorType>::operator[] (const size_type i) const
{
  return operator()(i);
}



template <class VectorType>
inline
typename BlockVectorBase<VectorType>::reference
BlockVectorBase<VectorType>::operator[] (const size_type i)
{
  return operator()(i);
}



template <typename VectorType>
template <typename OtherNumber>
inline
void BlockVectorBase<VectorType>::extract_subvector_to (const std::vector<size_type> &indices,
                                                        std::vector<OtherNumber> &values) const
{
  for (size_type i = 0; i < indices.size(); ++i)
    values[i] = operator()(indices[i]);
}



template <typename VectorType>
template <typename ForwardIterator, typename OutputIterator>
inline
void BlockVectorBase<VectorType>::extract_subvector_to (ForwardIterator          indices_begin,
                                                        const ForwardIterator    indices_end,
                                                        OutputIterator           values_begin) const
{
  while (indices_begin != indices_end)
    {
      *values_begin = operator()(*indices_begin);
      indices_begin++;
      values_begin++;
    }
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

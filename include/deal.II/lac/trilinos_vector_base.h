// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2014 by the deal.II authors
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

#ifndef __deal2__trilinos_vector_base_h
#define __deal2__trilinos_vector_base_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS

#include <deal.II/base/utilities.h>
#  include <deal.II/base/std_cxx11/shared_ptr.h>
#  include <deal.II/base/subscriptor.h>
#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/vector.h>

#  include <vector>
#  include <utility>
#  include <memory>

#  define TrilinosScalar double
#  include "Epetra_ConfigDefs.h"
#  ifdef DEAL_II_WITH_MPI // only if MPI is installed
#    include "mpi.h"
#    include "Epetra_MpiComm.h"
#  else
#    include "Epetra_SerialComm.h"
#  endif
#  include "Epetra_FEVector.h"

DEAL_II_NAMESPACE_OPEN

// forward declaration
template <typename number> class Vector;


/**
 * @addtogroup TrilinosWrappers
 *@{
 */
namespace TrilinosWrappers
{
  // forward declaration
  class VectorBase;

  /**
   * @cond internal
   */

  /**
   * A namespace for internal implementation details of the
   * TrilinosWrapper members.
   *
   * @ingroup TrilinosWrappers
   */
  namespace internal
  {
    /**
     * Declare type for container size.
     */
    typedef dealii::types::global_dof_index size_type;

    /**
     * This class implements a wrapper for accessing the Trilinos vector in
     * the same way as we access deal.II objects: it is initialized with a
     * vector and an element within it, and has a conversion operator to
     * extract the scalar value of this element. It also has a variety of
     * assignment operator for writing to this one element.
     *
     * @ingroup TrilinosWrappers
     */
    class VectorReference
    {
    private:
      /**
       * Constructor. It is made private so as to only allow the actual vector
       * class to create it.
       */
      VectorReference (VectorBase     &vector,
                       const size_type  index);

    public:

      /**
       * This looks like a copy operator, but does something different than
       * usual. In particular, it does not copy the member variables of this
       * reference. Rather, it handles the situation where we have two vectors
       * @p v and @p w, and assign elements like in <tt>v(i)=w(i)</tt>. Here,
       * both left and right hand side of the assignment have data type
       * VectorReference, but what we really mean is to assign the vector
       * elements represented by the two references. This operator implements
       * this operation. Note also that this allows us to make the assignment
       * operator const.
       */
      const VectorReference &
      operator = (const VectorReference &r) const;

      /**
       * Same as above but for non-const reference objects.
       */
      const VectorReference &
      operator = (const VectorReference &r);

      /**
       * Set the referenced element of the vector to <tt>s</tt>.
       */
      const VectorReference &
      operator = (const TrilinosScalar &s) const;

      /**
       * Add <tt>s</tt> to the referenced element of the vector->
       */
      const VectorReference &
      operator += (const TrilinosScalar &s) const;

      /**
       * Subtract <tt>s</tt> from the referenced element of the vector->
       */
      const VectorReference &
      operator -= (const TrilinosScalar &s) const;

      /**
       * Multiply the referenced element of the vector by <tt>s</tt>.
       */
      const VectorReference &
      operator *= (const TrilinosScalar &s) const;

      /**
       * Divide the referenced element of the vector by <tt>s</tt>.
       */
      const VectorReference &
      operator /= (const TrilinosScalar &s) const;

      /**
       * Convert the reference to an actual value, i.e. return the value of
       * the referenced element of the vector.
       */
      operator TrilinosScalar () const;

      /**
       * Exception
       */
      DeclException1 (ExcTrilinosError,
                      int,
                      << "An error with error number " << arg1
                      << " occurred while calling a Trilinos function");

    private:
      /**
       * Point to the vector we are referencing.
       */
      VectorBase   &vector;

      /**
       * Index of the referenced element of the vector.
       */
      const size_type  index;

      /**
       * Make the vector class a friend, so that it can create objects of the
       * present type.
       */
      friend class ::dealii::TrilinosWrappers::VectorBase;
    };
  }
  /**
   * @endcond
   */


  /**
   * Base class for the two types of Trilinos vectors, the distributed
   * memory vector MPI::Vector and a localized vector Vector. The latter
   * is designed for use in either serial implementations or as a
   * localized copy on each processor.  The implementation of this class
   * is based on the Trilinos vector class Epetra_FEVector, the (parallel)
   * partitioning of which is governed by an Epetra_Map. This means that
   * the vector type is generic and can be done in this base class, while
   * the definition of the partition map (and hence, the constructor and
   * reinit function) will have to be done in the derived classes. The
   * Epetra_FEVector is precisely the kind of vector we deal with all the
   * time - we probably get it from some assembly process, where also
   * entries not locally owned might need to written and hence need to be
   * forwarded to the owner. The only requirement for this class to work
   * is that Trilinos is installed with the same compiler as is used for
   * compilation of deal.II.
   *
   * The interface of this class is modeled after the existing Vector
   * class in deal.II. It has almost the same member functions, and is
   * often exchangable. However, since Trilinos only supports a single
   * scalar type (double), it is not templated, and only works with that
   * type.
   *
   * Note that Trilinos only guarantees that operations do what you expect
   * if the function @p GlobalAssemble has been called after vector
   * assembly in order to distribute the data. Therefore, you need to call
   * Vector::compress() before you actually use the vectors.
   *
   * @ingroup TrilinosWrappers
   * @ingroup Vectors
   * @author Martin Kronbichler, 2008
   */
  class VectorBase : public Subscriptor
  {
  public:
    /**
     * Declare some of the standard types used in all containers. These types
     * parallel those in the <tt>C</tt> standard libraries
     * <tt>vector<...></tt> class.
     */
    typedef TrilinosScalar                  value_type;
    typedef TrilinosScalar                  real_type;
    typedef dealii::types::global_dof_index size_type;
    typedef value_type                     *iterator;
    typedef const value_type               *const_iterator;
    typedef internal::VectorReference       reference;
    typedef const internal::VectorReference const_reference;

    /**
     * @name 1: Basic Object-handling
     */
    //@{

    /**
     * Default constructor that generates an empty (zero size) vector. The
     * function <tt>reinit()</tt> will have to give the vector the correct
     * size and distribution among processes in case of an MPI run.
     */
    VectorBase ();

    /**
     * Copy constructor. Sets the dimension to that of the given vector, and
     * copies all the elements.
     */
    VectorBase (const VectorBase &v);

    /**
     * Destructor
     */
    virtual ~VectorBase ();

    /**
     * Release all memory and return to a state just like after having called
     * the default constructor.
     */
    void clear ();

    /**
     * Reinit functionality, sets the dimension and possibly the parallel
     * partitioning (Epetra_Map) of the calling vector to the settings of the
     * input vector.
     */
    void reinit (const VectorBase &v,
                 const bool        fast = false);

    /**
     * Compress the underlying representation of the Trilinos object,
     * i.e. flush the buffers of the vector object if it has any. This
     * function is necessary after writing into a vector element-by-element
     * and before anything else can be done on it.
     *
     * The (defaulted) argument can be used to specify the compress mode
     * (<code>Add</code> or <code>Insert</code>) in case the vector has not
     * been written to since the last time this function was called. The
     * argument is ignored if the vector has been added or written to since
     * the last time compress() was called.
     *
     * See @ref GlossCompress "Compressing distributed objects"
     * for more information.
     */
    void compress (::dealii::VectorOperation::values operation);

    /**
     * @deprecated: Use the compress(VectorOperation::values) function
     * above instead.
     */
    void compress() DEAL_II_DEPRECATED;

    /**
    * @deprecated Use compress(dealii::VectorOperation::values) instead.
    */
    void compress (const Epetra_CombineMode last_action) DEAL_II_DEPRECATED;

    /**
     * Returns the state of the vector, i.e., whether compress() has already
     * been called after an operation requiring data exchange.
     */
    bool is_compressed () const;

    /**
     * Set all components of the vector to the given number @p s. Simply pass
     * this down to the Trilinos Epetra object, but we still need to declare
     * this function to make the example given in the discussion about making
     * the constructor explicit work.
     *
     * Since the semantics of assigning a scalar to a vector are not
     * immediately clear, this operator should really only be used if you want
     * to set the entire vector to zero. This allows the intuitive notation
     * <tt>v=0</tt>. Assigning other values is deprecated and may be
     * disallowed in the future.
     */
    VectorBase &
    operator = (const TrilinosScalar s);

    /**
     * Copy function. This function takes a VectorBase vector and copies all
     * the elements. The target vector will have the same parallel
     * distribution as the calling vector.
     */
    VectorBase &
    operator = (const VectorBase &v);

    /**
     * Another copy function. This one takes a deal.II vector and copies it
     * into a TrilinosWrapper vector. Note that since we do not provide any
     * Epetra_map that tells about the partitioning of the vector among the
     * MPI processes, the size of the TrilinosWrapper vector has to be the
     * same as the size of the input vector. In order to change the map, use
     * the reinit(const Epetra_Map &input_map) function.
     */
    template <typename Number>
    VectorBase &
    operator = (const ::dealii::Vector<Number> &v);

    /**
     * Test for equality. This function assumes that the present vector and
     * the one to compare with have the same size already, since comparing
     * vectors of different sizes makes not much sense anyway.
     */
    bool operator == (const VectorBase &v) const;

    /**
     * Test for inequality. This function assumes that the present vector and
     * the one to compare with have the same size already, since comparing
     * vectors of different sizes makes not much sense anyway.
     */
    bool operator != (const VectorBase &v) const;

    /**
     * Return the global dimension of the vector.
     */
    size_type size () const;

    /**
     * Return the local dimension of the vector, i.e. the number of elements
     * stored on the present MPI process. For sequential vectors, this number
     * is the same as size(), but for parallel vectors it may be smaller.
     *
     * To figure out which elements exactly are stored locally, use
     * local_range().
     *
     * If the vector contains ghost elements, they are included in this
     * number.
     */
    size_type local_size () const;

    /**
     * Return a pair of indices indicating which elements of this vector are
     * stored locally. The first number is the index of the first element
     * stored, the second the index of the one past the last one that is
     * stored locally. If this is a sequential vector, then the result will be
     * the pair <code>(0,N)</code>, otherwise it will be a pair
     * <code>(i,i+n)</code>, where <code>n=local_size()</code> and
     * <code>i</code> is the first element of the vector stored on this
     * processor, corresponding to the half open interval $[i,i+n)$
     *
     * @note The description above is true most of the time, but
     * not always. In particular, Trilinos vectors need not store
     * contiguous ranges of elements such as $[i,i+n)$. Rather, it
     * can store vectors where the elements are distributed in
     * an arbitrary way across all processors and each processor
     * simply stores a particular subset, not necessarily contiguous.
     * In this case, this function clearly makes no sense since it
     * could, at best, return a range that includes all elements
     * that are stored locally. Thus, the function only succeeds
     * if the locally stored range is indeed contiguous. It will
     * trigger an assertion if the local portion of the vector
     * is not contiguous.
     */
    std::pair<size_type, size_type> local_range () const;

    /**
     * Return whether @p index is in the local range or not, see also
     * local_range().
     *
     * @note The same limitation for the applicability of this
     * function applies as listed in the documentation of local_range().
     */
    bool in_local_range (const size_type index) const;

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
     */
    IndexSet locally_owned_elements () const;

    /**
     * Return if the vector contains ghost elements. This answer is true if
     * there are ghost elements on at least one process.
     */
    bool has_ghost_elements() const;

    /**
     * Return the scalar (inner) product of two vectors. The vectors must have
     * the same size.
     */
    TrilinosScalar operator * (const VectorBase &vec) const;

    /**
     * Return square of the $l_2$-norm.
     */
    real_type norm_sqr () const;

    /**
     * Mean value of the elements of this vector.
     */
    TrilinosScalar mean_value () const;

    /**
     * Compute the minimal value of the elements of this vector.
     */
    TrilinosScalar minimal_value () const;

    /**
     * $l_1$-norm of the vector.  The sum of the absolute values.
     */
    real_type l1_norm () const;

    /**
     * $l_2$-norm of the vector.  The square root of the sum of the squares of
     * the elements.
     */
    real_type l2_norm () const;

    /**
     * $l_p$-norm of the vector. The <i>p</i>th root of the sum of the
     * <i>p</i>th powers of the absolute values of the elements.
     */
    real_type lp_norm (const TrilinosScalar p) const;

    /**
     * Maximum absolute value of the elements.
     */
    real_type linfty_norm () const;

    /**
     * Return whether the vector contains only elements with value
     * zero. This is a collective operation. This function is expensive, because
     * potentially all elements have to be checked.
     */
    bool all_zero () const;

    /**
     * Return @p true if the vector has no negative entries, i.e. all entries
     * are zero or positive. This function is used, for example, to check
     * whether refinement indicators are really all positive (or zero).
     */
    bool is_non_negative () const;
    //@}


    /**
     * @name 2: Data-Access
     */
    //@{

    /**
     * Provide access to a given element, both read and write.
     *
     * When using a vector distributed with MPI, this operation only makes
     * sense for elements that are actually present on the calling
     * processor. Otherwise, an exception is thrown. This is different from
     * the <code>el()</code> function below that always succeeds (but returns
     * zero on non-local elements).
     */
    reference
    operator () (const size_type index);

    /**
     * Provide read-only access to an element.
     *
     * When using a vector distributed with MPI, this operation only makes
     * sense for elements that are actually present on the calling
     * processor. Otherwise, an exception is thrown. This is different from
     * the <code>el()</code> function below that always succeeds (but returns
     * zero on non-local elements).
     */
    TrilinosScalar
    operator () (const size_type index) const;

    /**
     * Provide access to a given element, both read and write.
     *
     * Exactly the same as operator().
     */
    reference
    operator [] (const size_type index);

    /**
     * Provide read-only access to an element.
     *
     * Exactly the same as operator().
     */
    TrilinosScalar
    operator [] (const size_type index) const;

    /**
     * A collective get operation: instead of getting individual elements of a
     * vector, this function allows to get a whole set of elements at
     * once. The indices of the elements to be read are stated in the first
     * argument, the corresponding values are returned in the second.
     */
    void extract_subvector_to (const std::vector<size_type> &indices,
                               std::vector<TrilinosScalar> &values) const;

    /**
     * Just as the above, but with pointers.  Useful in minimizing copying of
     * data around.
     */
    template <typename ForwardIterator, typename OutputIterator>
    void extract_subvector_to (ForwardIterator          indices_begin,
                               const ForwardIterator    indices_end,
                               OutputIterator           values_begin) const;

    /**
     * Return the value of the vector entry <i>i</i>. Note that this function
     * does only work properly when we request a data stored on the local
     * processor. In case the elements sits on another process, this function
     * returns 0 which might or might not be appropriate in a given
     * situation. If you rely on consistent results, use the access functions
     * () or [] that throw an assertion in case a non-local element is used.
     */
    TrilinosScalar el (const size_type index) const;

    /**
     * Make the Vector class a bit like the <tt>vector<></tt> class of the C++
     * standard library by returning iterators to the start and end of the
     * locally owned elements of this vector. The ordering of local elements
     * corresponds to the one given by the global indices in case the vector
     * is constructed from an IndexSet or other methods in deal.II (note that
     * an Epetra_Map can contain elements in arbitrary orders, though).
     *
     * It holds that end() - begin() == local_size().
     */
    iterator begin ();

    /**
     * Return constant iterator to the start of the locally owned elements
     * of the vector.
     */
    const_iterator begin () const;

    /**
     * Return an iterator pointing to the element past the end of the array
     * of locally owned entries.
     */
    iterator end ();

    /**
     * Return a constant iterator pointing to the element past the end of the
     * array of the locally owned entries.
     */
    const_iterator end () const;

    //@}


    /**
     * @name 3: Modification of vectors
     */
    //@{

    /**
     * A collective set operation: instead of setting individual elements of a
     * vector, this function allows to set a whole set of elements at
     * once. The indices of the elements to be set are stated in the first
     * argument, the corresponding values in the second.
     */
    void set (const std::vector<size_type>    &indices,
              const std::vector<TrilinosScalar>  &values);

    /**
     * This is a second collective set operation. As a difference, this
     * function takes a deal.II vector of values.
     */
    void set (const std::vector<size_type>        &indices,
              const ::dealii::Vector<TrilinosScalar> &values);

    /**
     * This collective set operation is of lower level and can handle anything
     * else &mdash; the only thing you have to provide is an address where all
     * the indices are stored and the number of elements to be set.
     */
    void set (const size_type       n_elements,
              const size_type      *indices,
              const TrilinosScalar *values);

    /**
     * A collective add operation: This function adds a whole set of values
     * stored in @p values to the vector components specified by @p indices.
     */
    void add (const std::vector<size_type>      &indices,
              const std::vector<TrilinosScalar> &values);

    /**
     * This is a second collective add operation. As a difference, this
     * function takes a deal.II vector of values.
     */
    void add (const std::vector<size_type>           &indices,
              const ::dealii::Vector<TrilinosScalar> &values);

    /**
     * Take an address where <tt>n_elements</tt> are stored contiguously and
     * add them into the vector. Handles all cases which are not covered by
     * the other two <tt>add()</tt> functions above.
     */
    void add (const size_type       n_elements,
              const size_type      *indices,
              const TrilinosScalar *values);

    /**
     * Multiply the entire vector by a fixed factor.
     */
    VectorBase &operator *= (const TrilinosScalar factor);

    /**
     * Divide the entire vector by a fixed factor.
     */
    VectorBase &operator /= (const TrilinosScalar factor);

    /**
     * Add the given vector to the present one.
     */
    VectorBase &operator += (const VectorBase &V);

    /**
     * Subtract the given vector from the present one.
     */
    VectorBase &operator -= (const VectorBase &V);

    /**
     * Addition of @p s to all components. Note that @p s is a scalar and not
     * a vector.
     */
    void add (const TrilinosScalar s);

    /**
     * Simple vector addition, equal to the <tt>operator +=</tt>.
     *
     * Though, if the second argument <tt>allow_different_maps</tt> is set,
     * then it is possible to add data from a vector that uses a different
     * map, i.e., a vector whose elements are split across processors
     * differently. This may include vectors with ghost elements, for example.
     * In general, however, adding vectors with a different element-to-processor
     * map requires communicating data among processors and, consequently,
     * is a slower operation than when using vectors using the same map.
     */
    void add (const VectorBase &V,
              const bool        allow_different_maps = false);

    /**
     * Simple addition of a multiple of a vector, i.e. <tt>*this += a*V</tt>.
     */
    void add (const TrilinosScalar  a,
              const VectorBase     &V);

    /**
     * Multiple addition of scaled vectors, i.e. <tt>*this += a*V + b*W</tt>.
     */
    void add (const TrilinosScalar  a,
              const VectorBase     &V,
              const TrilinosScalar  b,
              const VectorBase     &W);

    /**
     * Scaling and simple vector addition, i.e.  <tt>*this = s*(*this) +
     * V</tt>.
     */
    void sadd (const TrilinosScalar  s,
               const VectorBase     &V);

    /**
     * Scaling and simple addition, i.e.  <tt>*this = s*(*this) + a*V</tt>.
     */
    void sadd (const TrilinosScalar  s,
               const TrilinosScalar  a,
               const VectorBase     &V);

    /**
     * Scaling and multiple addition.
     */
    void sadd (const TrilinosScalar  s,
               const TrilinosScalar  a,
               const VectorBase     &V,
               const TrilinosScalar  b,
               const VectorBase     &W);

    /**
     * Scaling and multiple addition.  <tt>*this = s*(*this) + a*V + b*W +
     * c*X</tt>.
     */
    void sadd (const TrilinosScalar  s,
               const TrilinosScalar  a,
               const VectorBase     &V,
               const TrilinosScalar  b,
               const VectorBase     &W,
               const TrilinosScalar  c,
               const VectorBase     &X);

    /**
     * Scale each element of this vector by the corresponding element in the
     * argument. This function is mostly meant to simulate multiplication (and
     * immediate re-assignment) by a diagonal scaling matrix.
     */
    void scale (const VectorBase &scaling_factors);

    /**
     * Assignment <tt>*this = a*V</tt>.
     */
    void equ (const TrilinosScalar  a,
              const VectorBase     &V);

    /**
     * Assignment <tt>*this = a*V + b*W</tt>.
     */
    void equ (const TrilinosScalar  a,
              const VectorBase     &V,
              const TrilinosScalar  b,
              const VectorBase     &W);

    /**
     * Compute the elementwise ratio of the two given vectors, that is let
     * <tt>this[i] = a[i]/b[i]</tt>. This is useful for example if you want to
     * compute the cellwise ratio of true to estimated error.
     *
     * This vector is appropriately scaled to hold the result.
     *
     * If any of the <tt>b[i]</tt> is zero, the result is undefined. No
     * attempt is made to catch such situations.
     */
    void ratio (const VectorBase &a,
                const VectorBase &b);
    //@}


    /**
     * @name 4: Mixed stuff
     */
    //@{

    /**
     * Return a const reference to the underlying Trilinos Epetra_MultiVector
     * class.
     */
    const Epetra_MultiVector &trilinos_vector () const;

    /**
     * Return a (modifyable) reference to the underlying Trilinos
     * Epetra_FEVector class.
     */
    Epetra_FEVector &trilinos_vector ();

    /**
     * Return a const reference to the underlying Trilinos Epetra_Map that
     * sets the parallel partitioning of the vector.
     */
    const Epetra_Map &vector_partitioner () const;

    /**
     *  Output of vector in user-defined format in analogy to the
     *  dealii::Vector class.
     */
    void print (const char *format = 0) const;

    /**
     * Print to a stream. @p precision denotes the desired precision with
     * which values shall be printed, @p scientific whether scientific
     * notation shall be used. If @p across is @p true then the vector is
     * printed in a line, while if @p false then the elements are printed on a
     * separate line each.
     */
    void print (std::ostream       &out,
                const unsigned int  precision  = 3,
                const bool          scientific = true,
                const bool          across     = true) const;

    /**
     * Swap the contents of this vector and the other vector @p v. One could
     * do this operation with a temporary variable and copying over the data
     * elements, but this function is significantly more efficient since it
     * only swaps the pointers to the data of the two vectors and therefore
     * does not need to allocate temporary storage and move data around. Note
     * that the vectors need to be of the same size and base on the same map.
     *
     * This function is analog to the the @p swap function of all C standard
     * containers. Also, there is a global function <tt>swap(u,v)</tt> that
     * simply calls <tt>u.swap(v)</tt>, again in analogy to standard
     * functions.
     */
    void swap (VectorBase &v);

    /**
     * Estimate for the memory consumption in bytes.
     */
    std::size_t memory_consumption () const;

    /**
     * Return a reference to the MPI communicator object in use with this
     * object.
     */
    const MPI_Comm &get_mpi_communicator () const;
    //@}

    /**
     * Exception
     */
    DeclException0 (ExcGhostsPresent);

    /**
     * Exception
     */
    DeclException0 (ExcDifferentParallelPartitioning);

    /**
     * Exception
     */
    DeclException1 (ExcTrilinosError,
                    int,
                    << "An error with error number " << arg1
                    << " occurred while calling a Trilinos function");

    /**
     * Exception
     */
    DeclException4 (ExcAccessToNonLocalElement,
                    size_type, size_type, size_type, size_type,
                    << "You tried to access element " << arg1
                    << " of a distributed vector, but this element is not stored "
                    << "on the current processor. Note: There are "
                    << arg2 << " elements stored "
                    << "on the current processor from within the range "
                    << arg3 << " through " << arg4
                    << " but Trilinos vectors need not store contiguous "
                    << "ranges on each processor, and not every element in "
                    << "this range may in fact be stored locally.");


  private:
    /**
     * Trilinos doesn't allow to mix additions to matrix entries and
     * overwriting them (to make synchronisation of parallel computations
     * simpler). The way we do it is to, for each access operation, store
     * whether it is an insertion or an addition. If the previous one was of
     * different type, then we first have to flush the Trilinos buffers;
     * otherwise, we can simply go on.  Luckily, Trilinos has an object for
     * this which does already all the parallel communications in such a case,
     * so we simply use their model, which stores whether the last operation
     * was an addition or an insertion.
     */
    Epetra_CombineMode last_action;

    /**
     * A boolean variable to hold information on whether the vector is
     * compressed or not.
     */
    bool compressed;

    /**
     * Whether this vector has ghost elements. This is true on all processors
     * even if only one of them has any ghost elements.
     */
    bool has_ghosts;

    /**
     * Pointer to the actual Epetra vector object. This may represent a
     * vector that is in fact distributed among multiple processors. The
     * object requires an existing Epetra_Map for
     * storing data when setting it up.
     */
    std_cxx11::shared_ptr<Epetra_FEVector> vector;

    /**
     * A vector object in Trilinos to be used for collecting the non-local
     * elements if the vector was constructed with an additional IndexSet
     * describing ghost elements.
     */
    std_cxx11::shared_ptr<Epetra_MultiVector> nonlocal_vector;

    /**
     * Make the reference class a friend.
     */
    friend class internal::VectorReference;
    friend class Vector;
    friend class MPI::Vector;
  };




// ------------------- inline and template functions --------------

  /**
   * Global function swap which overloads the default implementation of
   * the C standard library which uses a temporary object. The function
   * simply exchanges the data of the two vectors.
   *
   * @relates TrilinosWrappers::VectorBase
   * @author Martin Kronbichler, Wolfgang Bangerth, 2008
   */
  inline
  void swap (VectorBase &u, VectorBase &v)
  {
    u.swap (v);
  }


#ifndef DOXYGEN

  namespace internal
  {
    inline
    VectorReference::VectorReference (VectorBase      &vector,
                                      const size_type  index)
      :
      vector (vector),
      index (index)
    {}


    inline
    const VectorReference &
    VectorReference::operator = (const VectorReference &r) const
    {
      // as explained in the class
      // documentation, this is not the copy
      // operator. so simply pass on to the
      // "correct" assignment operator
      *this = static_cast<TrilinosScalar> (r);

      return *this;
    }



    inline
    const VectorReference &
    VectorReference::operator = (const VectorReference &r)
    {
      // as above
      *this = static_cast<TrilinosScalar> (r);

      return *this;
    }


    inline
    const VectorReference &
    VectorReference::operator = (const TrilinosScalar &value) const
    {
      vector.set (1, &index, &value);
      return *this;
    }



    inline
    const VectorReference &
    VectorReference::operator += (const TrilinosScalar &value) const
    {
      vector.add (1, &index, &value);
      return *this;
    }



    inline
    const VectorReference &
    VectorReference::operator -= (const TrilinosScalar &value) const
    {
      TrilinosScalar new_value = -value;
      vector.add (1, &index, &new_value);
      return *this;
    }



    inline
    const VectorReference &
    VectorReference::operator *= (const TrilinosScalar &value) const
    {
      TrilinosScalar new_value = static_cast<TrilinosScalar>(*this) * value;
      vector.set (1, &index, &new_value);
      return *this;
    }



    inline
    const VectorReference &
    VectorReference::operator /= (const TrilinosScalar &value) const
    {
      TrilinosScalar new_value = static_cast<TrilinosScalar>(*this) / value;
      vector.set (1, &index, &new_value);
      return *this;
    }
  }



  inline
  bool
  VectorBase::is_compressed () const
  {
    return compressed;
  }



  inline
  bool
  VectorBase::in_local_range (const size_type index) const
  {
    std::pair<size_type, size_type> range = local_range();

    return ((index >= range.first) && (index <  range.second));
  }



  inline
  IndexSet
  VectorBase::locally_owned_elements() const
  {
    IndexSet is (size());

    // easy case: local range is contiguous
    if (vector->Map().LinearMap())
      {
        const std::pair<size_type, size_type> x = local_range();
        is.add_range (x.first, x.second);
      }
    else if (vector->Map().NumMyElements() > 0)
      {
        const size_type n_indices = vector->Map().NumMyElements();
#ifndef DEAL_II_WITH_64BIT_INDICES
        unsigned int *vector_indices = (unsigned int *)vector->Map().MyGlobalElements();
#else
        size_type *vector_indices = (size_type *)vector->Map().MyGlobalElements64();
#endif
        is.add_indices(vector_indices, vector_indices+n_indices);
        is.compress();
      }

    return is;
  }



  inline
  bool
  VectorBase::has_ghost_elements() const
  {
    return has_ghosts;
  }



  inline
  internal::VectorReference
  VectorBase::operator () (const size_type index)
  {
    return internal::VectorReference (*this, index);
  }



  inline
  internal::VectorReference
  VectorBase::operator [] (const size_type index)
  {
    return operator() (index);
  }


  inline
  TrilinosScalar
  VectorBase::operator [] (const size_type index) const
  {
    return operator() (index);
  }



  inline
  void VectorBase::extract_subvector_to (const std::vector<size_type> &indices,
                                         std::vector<TrilinosScalar>  &values) const
  {
    for (size_type i = 0; i < indices.size(); ++i)
      values[i] = operator()(indices[i]);
  }



  template <typename ForwardIterator, typename OutputIterator>
  inline
  void VectorBase::extract_subvector_to (ForwardIterator          indices_begin,
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



  inline
  VectorBase::iterator
  VectorBase::begin()
  {
    return (*vector)[0];
  }



  inline
  VectorBase::iterator
  VectorBase::end()
  {
    return (*vector)[0]+local_size();
  }



  inline
  VectorBase::const_iterator
  VectorBase::begin() const
  {
    return (*vector)[0];
  }



  inline
  VectorBase::const_iterator
  VectorBase::end() const
  {
    return (*vector)[0]+local_size();
  }



  inline
  void
  VectorBase::reinit (const VectorBase &v,
                      const bool        fast)
  {
    Assert (vector.get() != 0,
            ExcMessage("Vector has not been constructed properly."));

    if (fast == false ||
        vector_partitioner().SameAs(v.vector_partitioner())==false)
      vector.reset (new Epetra_FEVector(*v.vector));

    if (v.nonlocal_vector.get() != 0)
      nonlocal_vector.reset(new Epetra_MultiVector(v.nonlocal_vector->Map(), 1));
  }



  inline
  void
  VectorBase::compress (const Epetra_CombineMode last_action)
  {
    ::dealii::VectorOperation::values last_action_ =
      ::dealii::VectorOperation::unknown;
    if (last_action == Add)
      last_action_ = ::dealii::VectorOperation::add;
    else if (last_action == Insert)
      last_action_ = ::dealii::VectorOperation::insert;
    else
      AssertThrow(false, ExcNotImplemented());

    compress(last_action_);
  }



  inline
  void
  VectorBase::compress ()
  {
    compress(VectorOperation::unknown);
  }



  inline
  VectorBase &
  VectorBase::operator = (const TrilinosScalar s)
  {
    // if we have ghost values, do not allow
    // writing to this vector at all.
    Assert (!has_ghost_elements(), ExcGhostsPresent());

    Assert (numbers::is_finite(s), ExcNumberNotFinite());

    const int ierr = vector->PutScalar(s);

    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    if (nonlocal_vector.get() != 0)
      nonlocal_vector->PutScalar(0.);

    return *this;
  }



  inline
  void
  VectorBase::set (const std::vector<size_type>      &indices,
                   const std::vector<TrilinosScalar>  &values)
  {
    // if we have ghost values, do not allow
    // writing to this vector at all.
    Assert (!has_ghost_elements(), ExcGhostsPresent());

    Assert (indices.size() == values.size(),
            ExcDimensionMismatch(indices.size(),values.size()));

    set (indices.size(), &indices[0], &values[0]);
  }



  inline
  void
  VectorBase::set (const std::vector<size_type>           &indices,
                   const ::dealii::Vector<TrilinosScalar> &values)
  {
    // if we have ghost values, do not allow
    // writing to this vector at all.
    Assert (!has_ghost_elements(), ExcGhostsPresent());

    Assert (indices.size() == values.size(),
            ExcDimensionMismatch(indices.size(),values.size()));

    set (indices.size(), &indices[0], values.begin());
  }



  inline
  void
  VectorBase::set (const size_type       n_elements,
                   const size_type      *indices,
                   const TrilinosScalar *values)
  {
    // if we have ghost values, do not allow
    // writing to this vector at all.
    Assert (!has_ghost_elements(), ExcGhostsPresent());

    if (last_action == Add)
      vector->GlobalAssemble(Add);

    if (last_action != Insert)
      last_action = Insert;

    for (size_type i=0; i<n_elements; ++i)
      {
        const size_type row = indices[i];
        const TrilinosWrappers::types::int_type local_row = vector->Map().LID(static_cast<TrilinosWrappers::types::int_type>(row));
        if (local_row != -1)
          (*vector)[0][local_row] = values[i];
        else
          {
            const int ierr = vector->ReplaceGlobalValues (1,
                                                          (const TrilinosWrappers::types::int_type *)(&row),
                                                          &values[i]);
            AssertThrow (ierr == 0, ExcTrilinosError(ierr));
            compressed = false;
          }
        // in set operation, do not use the pre-allocated vector for nonlocal
        // entries even if it exists. This is to ensure that we really only
        // set the elements touched by the set() method and not all contained
        // in the nonlocal entries vector (there is no way to distinguish them
        // on the receiving processor)
      }
  }



  inline
  void
  VectorBase::add (const std::vector<size_type>      &indices,
                   const std::vector<TrilinosScalar>  &values)
  {
    // if we have ghost values, do not allow
    // writing to this vector at all.
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (indices.size() == values.size(),
            ExcDimensionMismatch(indices.size(),values.size()));

    add (indices.size(), &indices[0], &values[0]);
  }



  inline
  void
  VectorBase::add (const std::vector<size_type>           &indices,
                   const ::dealii::Vector<TrilinosScalar> &values)
  {
    // if we have ghost values, do not allow
    // writing to this vector at all.
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (indices.size() == values.size(),
            ExcDimensionMismatch(indices.size(),values.size()));

    add (indices.size(), &indices[0], values.begin());
  }



  inline
  void
  VectorBase::add (const size_type       n_elements,
                   const size_type      *indices,
                   const TrilinosScalar *values)
  {
    // if we have ghost values, do not allow
    // writing to this vector at all.
    Assert (!has_ghost_elements(), ExcGhostsPresent());

    if (last_action != Add)
      {
        if (last_action == Insert)
          vector->GlobalAssemble(Insert);
        last_action = Add;
      }

    for (size_type i=0; i<n_elements; ++i)
      {
        const size_type row = indices[i];
        const TrilinosWrappers::types::int_type local_row = vector->Map().LID(static_cast<TrilinosWrappers::types::int_type>(row));
        if (local_row != -1)
          (*vector)[0][local_row] += values[i];
        else if (nonlocal_vector.get() == 0)
          {
            const int ierr = vector->SumIntoGlobalValues (1,
                                                          (const TrilinosWrappers::types::int_type *)(&row),
                                                          &values[i]);
            AssertThrow (ierr == 0, ExcTrilinosError(ierr));
            compressed = false;
          }
        else
          {
            // use pre-allocated vector for non-local entries if it exists for
            // addition operation
            const TrilinosWrappers::types::int_type my_row = nonlocal_vector->Map().LID(static_cast<TrilinosWrappers::types::int_type>(row));
            Assert(my_row != -1,
                   ExcMessage("Attempted to write into off-processor vector entry "
                              "that has not be specified as being writable upon "
                              "initialization"));
            (*nonlocal_vector)[0][my_row] += values[i];
            compressed = false;
          }
      }
  }



  inline
  VectorBase::size_type
  VectorBase::size () const
  {
#ifndef DEAL_II_WITH_64BIT_INDICES
    return (size_type) (vector->Map().MaxAllGID() + 1 - vector->Map().MinAllGID());
#else
    return (size_type) (vector->Map().MaxAllGID64() + 1 - vector->Map().MinAllGID64());
#endif
  }



  inline
  VectorBase::size_type
  VectorBase::local_size () const
  {
    return (size_type) vector->Map().NumMyElements();
  }



  inline
  std::pair<VectorBase::size_type, VectorBase::size_type>
  VectorBase::local_range () const
  {
#ifndef DEAL_II_WITH_64BIT_INDICES
    const TrilinosWrappers::types::int_type begin = vector->Map().MinMyGID();
    const TrilinosWrappers::types::int_type end = vector->Map().MaxMyGID()+1;
#else
    const TrilinosWrappers::types::int_type begin = vector->Map().MinMyGID64();
    const TrilinosWrappers::types::int_type end = vector->Map().MaxMyGID64()+1;
#endif

    Assert (end-begin == vector->Map().NumMyElements(),
            ExcMessage ("This function only makes sense if the elements that this "
                        "vector stores on the current processor form a contiguous range. "
                        "This does not appear to be the case for the current vector."));

    return std::make_pair (begin, end);
  }



  inline
  TrilinosScalar
  VectorBase::operator * (const VectorBase &vec) const
  {
    Assert (vector->Map().SameAs(vec.vector->Map()),
            ExcDifferentParallelPartitioning());
    Assert (!has_ghost_elements(), ExcGhostsPresent());

    TrilinosScalar result;

    const int ierr = vector->Dot(*(vec.vector), &result);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return result;
  }



  inline
  VectorBase::real_type
  VectorBase::norm_sqr () const
  {
    const TrilinosScalar d = l2_norm();
    return d*d;
  }



  inline
  TrilinosScalar
  VectorBase::mean_value () const
  {
    Assert (!has_ghost_elements(), ExcGhostsPresent());

    TrilinosScalar mean;
    const int ierr = vector->MeanValue (&mean);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return mean;
  }



  inline
  TrilinosScalar
  VectorBase::minimal_value () const
  {
    TrilinosScalar min_value;
    const int ierr = vector->MinValue (&min_value);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return min_value;
  }



  inline
  VectorBase::real_type
  VectorBase::l1_norm () const
  {
    Assert (!has_ghost_elements(), ExcGhostsPresent());

    TrilinosScalar d;
    const int ierr = vector->Norm1 (&d);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return d;
  }



  inline
  VectorBase::real_type
  VectorBase::l2_norm () const
  {
    Assert (!has_ghost_elements(), ExcGhostsPresent());

    TrilinosScalar d;
    const int ierr = vector->Norm2 (&d);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return d;
  }



  inline
  VectorBase::real_type
  VectorBase::lp_norm (const TrilinosScalar p) const
  {
    Assert (!has_ghost_elements(), ExcGhostsPresent());

    TrilinosScalar norm = 0;
    TrilinosScalar sum=0;
    const size_type n_local = local_size();

    // loop over all the elements because
    // Trilinos does not support lp norms
    for (size_type i=0; i<n_local; i++)
      sum += std::pow(std::fabs((*vector)[0][i]), p);

    norm = std::pow(sum, static_cast<TrilinosScalar>(1./p));

    return norm;
  }



  inline
  VectorBase::real_type
  VectorBase::linfty_norm () const
  {
    // while we disallow the other
    // norm operations on ghosted
    // vectors, this particular norm
    // is safe to run even in the
    // presence of ghost elements
    TrilinosScalar d;
    const int ierr = vector->NormInf (&d);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return d;
  }



  // inline also scalar products, vector
  // additions etc. since they are all
  // representable by a single Trilinos
  // call. This reduces the overhead of the
  // wrapper class.
  inline
  VectorBase &
  VectorBase::operator *= (const TrilinosScalar a)
  {
    Assert (numbers::is_finite(a), ExcNumberNotFinite());

    const int ierr = vector->Scale(a);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return *this;
  }



  inline
  VectorBase &
  VectorBase::operator /= (const TrilinosScalar a)
  {
    Assert (numbers::is_finite(a), ExcNumberNotFinite());

    const TrilinosScalar factor = 1./a;

    Assert (numbers::is_finite(factor), ExcNumberNotFinite());

    const int ierr = vector->Scale(factor);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return *this;
  }



  inline
  VectorBase &
  VectorBase::operator += (const VectorBase &v)
  {
    Assert (size() == v.size(),
            ExcDimensionMismatch(size(), v.size()));
    Assert (vector->Map().SameAs(v.vector->Map()),
            ExcDifferentParallelPartitioning());

    const int ierr = vector->Update (1.0, *(v.vector), 1.0);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return *this;
  }



  inline
  VectorBase &
  VectorBase::operator -= (const VectorBase &v)
  {
    Assert (size() == v.size(),
            ExcDimensionMismatch(size(), v.size()));
    Assert (vector->Map().SameAs(v.vector->Map()),
            ExcDifferentParallelPartitioning());

    const int ierr = vector->Update (-1.0, *(v.vector), 1.0);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return *this;
  }



  inline
  void
  VectorBase::add (const TrilinosScalar s)
  {
    // if we have ghost values, do not allow
    // writing to this vector at all.
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (numbers::is_finite(s), ExcNumberNotFinite());

    size_type n_local = local_size();
    for (size_type i=0; i<n_local; i++)
      (*vector)[0][i] += s;
  }



  inline
  void
  VectorBase::add (const TrilinosScalar  a,
                   const VectorBase     &v)
  {
    // if we have ghost values, do not allow
    // writing to this vector at all.
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (local_size() == v.local_size(),
            ExcDimensionMismatch(local_size(), v.local_size()));

    Assert (numbers::is_finite(a), ExcNumberNotFinite());

    const int ierr = vector->Update(a, *(v.vector), 1.);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  inline
  void
  VectorBase::add (const TrilinosScalar  a,
                   const VectorBase     &v,
                   const TrilinosScalar  b,
                   const VectorBase     &w)
  {
    // if we have ghost values, do not allow
    // writing to this vector at all.
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (local_size() == v.local_size(),
            ExcDimensionMismatch(local_size(), v.local_size()));
    Assert (local_size() == w.local_size(),
            ExcDimensionMismatch(local_size(), w.local_size()));

    Assert (numbers::is_finite(a), ExcNumberNotFinite());
    Assert (numbers::is_finite(b), ExcNumberNotFinite());

    const int ierr = vector->Update(a, *(v.vector), b, *(w.vector), 1.);

    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  inline
  void
  VectorBase::sadd (const TrilinosScalar  s,
                    const VectorBase     &v)
  {
    // if we have ghost values, do not allow
    // writing to this vector at all.
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (size() == v.size(),
            ExcDimensionMismatch (size(), v.size()));

    Assert (numbers::is_finite(s), ExcNumberNotFinite());

    if (local_size() == v.local_size())
      {
        const int ierr = vector->Update(1., *(v.vector), s);
        AssertThrow (ierr == 0, ExcTrilinosError(ierr));
      }
    else
      {
        VectorBase tmp = v;
        tmp *= s;
        this->add(tmp, true);
      }
  }



  inline
  void
  VectorBase::sadd (const TrilinosScalar  s,
                    const TrilinosScalar  a,
                    const VectorBase     &v)
  {
    // if we have ghost values, do not allow
    // writing to this vector at all.
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (size() == v.size(),
            ExcDimensionMismatch (size(), v.size()));
    Assert (numbers::is_finite(s), ExcNumberNotFinite());
    Assert (numbers::is_finite(a), ExcNumberNotFinite());

    if (local_size() == v.local_size())
      {
        const int ierr = vector->Update(a, *(v.vector), s);
        AssertThrow (ierr == 0, ExcTrilinosError(ierr));
      }
    else
      {
        (*this)*=s;
        VectorBase tmp = v;
        tmp *= a;
        this->add(tmp, true);
      }
  }



  inline
  void
  VectorBase::sadd (const TrilinosScalar  s,
                    const TrilinosScalar  a,
                    const VectorBase     &v,
                    const TrilinosScalar  b,
                    const VectorBase     &w)
  {
    // if we have ghost values, do not allow
    // writing to this vector at all.
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (size() == v.size(),
            ExcDimensionMismatch (size(), v.size()));
    Assert (size() == w.size(),
            ExcDimensionMismatch (size(), w.size()));
    Assert (numbers::is_finite(s), ExcNumberNotFinite());
    Assert (numbers::is_finite(a), ExcNumberNotFinite());
    Assert (numbers::is_finite(b), ExcNumberNotFinite());

    if (local_size() == v.local_size() && local_size() == w.local_size())
      {
        const int ierr = vector->Update(a, *(v.vector), b, *(w.vector), s);
        AssertThrow (ierr == 0, ExcTrilinosError(ierr));
      }
    else
      {
        (*this)*=s;
        {
          VectorBase tmp = v;
          tmp *= a;
          this->add(tmp, true);
        }
        {
          VectorBase tmp = w;
          tmp *= b;
          this->add(tmp, true);
        }
      }
  }



  inline
  void
  VectorBase::sadd (const TrilinosScalar  s,
                    const TrilinosScalar  a,
                    const VectorBase     &v,
                    const TrilinosScalar  b,
                    const VectorBase     &w,
                    const TrilinosScalar  c,
                    const VectorBase     &x)
  {
    // if we have ghost values, do not allow
    // writing to this vector at all.
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (size() == v.size(),
            ExcDimensionMismatch (size(), v.size()));
    Assert (size() == w.size(),
            ExcDimensionMismatch (size(), w.size()));
    Assert (size() == x.size(),
            ExcDimensionMismatch (size(), x.size()));
    Assert (numbers::is_finite(s), ExcNumberNotFinite());
    Assert (numbers::is_finite(a), ExcNumberNotFinite());
    Assert (numbers::is_finite(b), ExcNumberNotFinite());
    Assert (numbers::is_finite(c), ExcNumberNotFinite());

    if (local_size() == v.local_size()
        && local_size() == w.local_size()
        && local_size() == x.local_size())
      {
        // Update member can only
        // input two other vectors so
        // do it in two steps
        const int ierr = vector->Update(a, *(v.vector), b, *(w.vector), s);
        AssertThrow (ierr == 0, ExcTrilinosError(ierr));

        const int jerr = vector->Update(c, *(x.vector), 1.);
        Assert (jerr == 0, ExcTrilinosError(jerr));
        (void)jerr; // removes -Wunused-parameter warning in optimized mode
      }
    else
      {
        (*this)*=s;
        {
          VectorBase tmp = v;
          tmp *= a;
          this->add(tmp, true);
        }
        {
          VectorBase tmp = w;
          tmp *= b;
          this->add(tmp, true);
        }
        {
          VectorBase tmp = x;
          tmp *= c;
          this->add(tmp, true);
        }
      }
  }



  inline
  void
  VectorBase::scale (const VectorBase &factors)
  {
    // if we have ghost values, do not allow
    // writing to this vector at all.
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (local_size() == factors.local_size(),
            ExcDimensionMismatch(local_size(), factors.local_size()));

    const int ierr = vector->Multiply (1.0, *(factors.vector), *vector, 0.0);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  inline
  void
  VectorBase::equ (const TrilinosScalar  a,
                   const VectorBase     &v)
  {
    // if we have ghost values, do not allow
    // writing to this vector at all.
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (numbers::is_finite(a), ExcNumberNotFinite());

    // If we don't have the same map, copy.
    if (vector->Map().SameAs(v.vector->Map())==false)
      {
        this->sadd(0., a, v);
      }
    else
      {
        // Otherwise, just update
        int ierr = vector->Update(a, *v.vector, 0.0);
        AssertThrow (ierr == 0, ExcTrilinosError(ierr));

        last_action = Zero;
      }

  }



  inline
  void
  VectorBase::equ (const TrilinosScalar  a,
                   const VectorBase     &v,
                   const TrilinosScalar  b,
                   const VectorBase     &w)
  {
    // if we have ghost values, do not allow
    // writing to this vector at all.
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (v.local_size() == w.local_size(),
            ExcDimensionMismatch (v.local_size(), w.local_size()));

    Assert (numbers::is_finite(a), ExcNumberNotFinite());
    Assert (numbers::is_finite(b), ExcNumberNotFinite());

    // If we don't have the same map, copy.
    if (vector->Map().SameAs(v.vector->Map())==false)
      {
        sadd(0., a, v, b, w);
      }
    else
      {
        // Otherwise, just update. verify
        // that *this does not only have
        // the same map as v (the
        // if-condition above) but also as
        // w
        Assert (vector->Map().SameAs(w.vector->Map()),
                ExcDifferentParallelPartitioning());
        int ierr = vector->Update(a, *v.vector, b, *w.vector, 0.0);
        AssertThrow (ierr == 0, ExcTrilinosError(ierr));

        last_action = Zero;
      }
  }



  inline
  void
  VectorBase::ratio (const VectorBase &v,
                     const VectorBase &w)
  {
    Assert (v.local_size() == w.local_size(),
            ExcDimensionMismatch (v.local_size(), w.local_size()));
    Assert (local_size() == w.local_size(),
            ExcDimensionMismatch (local_size(), w.local_size()));

    const int ierr = vector->ReciprocalMultiply(1.0, *(w.vector),
                                                *(v.vector), 0.0);

    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  inline
  const Epetra_MultiVector &
  VectorBase::trilinos_vector () const
  {
    return static_cast<const Epetra_MultiVector &>(*vector);
  }



  inline
  Epetra_FEVector &
  VectorBase::trilinos_vector ()
  {
    return *vector;
  }



  inline
  const Epetra_Map &
  VectorBase::vector_partitioner () const
  {
    return static_cast<const Epetra_Map &>(vector->Map());
  }



  inline
  const MPI_Comm &
  VectorBase::get_mpi_communicator () const
  {
    static MPI_Comm comm;

#ifdef DEAL_II_WITH_MPI

    const Epetra_MpiComm *mpi_comm
      = dynamic_cast<const Epetra_MpiComm *>(&vector->Map().Comm());
    comm = mpi_comm->Comm();

#else

    comm = MPI_COMM_SELF;

#endif

    return comm;
  }



#endif // DOXYGEN

}

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS

/*----------------------------   trilinos_vector_base.h     ---------------------------*/

#endif
/*----------------------------   trilinos_vector_base.h     ---------------------------*/

// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_petsc_vector_base_h
#  define dealii_petsc_vector_base_h


#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_PETSC

#    include <deal.II/base/index_set.h>
#    include <deal.II/base/subscriptor.h>

#    include <deal.II/lac/exceptions.h>
#    include <deal.II/lac/vector.h>
#    include <deal.II/lac/vector_operation.h>

#    include <petscvec.h>

#    include <utility>
#    include <vector>

DEAL_II_NAMESPACE_OPEN

// forward declaration
#    ifndef DOXYGEN
template <typename number>
class Vector;

namespace PETScWrappers
{
  class VectorBase;
}
#    endif

/**
 * A namespace in which wrapper classes for PETSc objects reside.
 *
 * @ingroup PETScWrappers
 * @ingroup Vectors
 */
namespace PETScWrappers
{
  /**
   * @cond internal
   */

  /**
   * A namespace for internal implementation details of the PETScWrapper
   * members.
   * @ingroup PETScWrappers
   */
  namespace internal
  {
    /**
     * Since access to PETSc vectors only goes through functions, rather than
     * by obtaining a reference to a vector element, we need a wrapper class
     * that acts as if it was a reference, and basically redirects all
     * accesses (read and write) to member functions of this class.
     *
     * This class implements such a wrapper: it is initialized with a vector
     * and an element within it, and has a conversion operator to extract the
     * scalar value of this element. It also has a variety of assignment
     * operator for writing to this one element.
     * @ingroup PETScWrappers
     */
    class VectorReference
    {
    public:
      /**
       * Declare type for container size.
       */
      using size_type = types::global_dof_index;

    private:
      /**
       * Constructor. It is made private so as to only allow the actual vector
       * class to create it.
       */
      VectorReference(const VectorBase &vector, const size_type index);

    public:
      /*
       * Copy constrcutor.
       */
      VectorReference(const VectorReference &vector) = default;

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
      operator=(const VectorReference &r) const;

      /**
       * The same function as above, but for non-const reference objects. The
       * function is needed since the compiler might otherwise automatically
       * generate a copy operator for non-const objects.
       */
      VectorReference &
      operator=(const VectorReference &r);

      /**
       * Set the referenced element of the vector to <tt>s</tt>.
       */
      const VectorReference &
      operator=(const PetscScalar &s) const;

      /**
       * Add <tt>s</tt> to the referenced element of the vector.
       */
      const VectorReference &
      operator+=(const PetscScalar &s) const;

      /**
       * Subtract <tt>s</tt> from the referenced element of the vector.
       */
      const VectorReference &
      operator-=(const PetscScalar &s) const;

      /**
       * Multiply the referenced element of the vector by <tt>s</tt>.
       */
      const VectorReference &
      operator*=(const PetscScalar &s) const;

      /**
       * Divide the referenced element of the vector by <tt>s</tt>.
       */
      const VectorReference &
      operator/=(const PetscScalar &s) const;

      /**
       * Return the real part of the value of the referenced element.
       */
      PetscReal
      real() const;

      /**
       * Return the imaginary part of the value of the referenced element.
       *
       * @note This operation is not defined for real numbers and an exception
       * is thrown.
       */
      PetscReal
      imag() const;

      /**
       * Convert the reference to an actual value, i.e. return the value of
       * the referenced element of the vector.
       */
      operator PetscScalar() const;
      /**
       * Exception
       */
      DeclException3(
        ExcAccessToNonlocalElement,
        int,
        int,
        int,
        << "You tried to access element " << arg1
        << " of a distributed vector, but only elements in range [" << arg2
        << "," << arg3 << "] are stored locally and can be accessed."
        << "\n\n"
        << "A common source for this kind of problem is that you "
        << "are passing a 'fully distributed' vector into a function "
        << "that needs read access to vector elements that correspond "
        << "to degrees of freedom on ghost cells (or at least to "
        << "'locally active' degrees of freedom that are not also "
        << "'locally owned'). You need to pass a vector that has these "
        << "elements as ghost entries.");
      /**
       * Exception.
       */
      DeclException2(ExcWrongMode,
                     int,
                     int,
                     << "You tried to do a "
                     << (arg1 == 1 ? "'set'" : (arg1 == 2 ? "'add'" : "???"))
                     << " operation but the vector is currently in "
                     << (arg2 == 1 ? "'set'" : (arg2 == 2 ? "'add'" : "???"))
                     << " mode. You first have to call 'compress()'.");

    private:
      /**
       * Point to the vector we are referencing.
       */
      const VectorBase &vector;

      /**
       * Index of the referenced element of the vector.
       */
      const size_type index;

      // Make the vector class a friend, so that it can create objects of the
      // present type.
      friend class ::dealii::PETScWrappers::VectorBase;
    };
  } // namespace internal
  /**
   * @endcond
   */


  /**
   * Base class for all vector classes that are implemented on top of the
   * PETSc vector types. Since in PETSc all vector types (i.e. sequential and
   * parallel ones) are built by filling the contents of an abstract object
   * that is only referenced through a pointer of a type that is independent
   * of the actual vector type, we can implement almost all functionality of
   * vectors in this base class. As such, this class can also be used as a
   * deal.II-compatible wrapper for a PETSc <code>Vec</code> object of any
   * type. Derived classes will then only have to provide the functionality to
   * create one or the other kind of vector.
   *
   * The interface of this class is modeled after the existing Vector class in
   * deal.II. It has almost the same member functions, and is often
   * exchangeable. However, since PETSc only supports a single scalar type
   * (either double, float, or a complex data type), it is not templated, and
   * only works with whatever your PETSc installation has defined the data
   * type @p PetscScalar to.
   *
   * Note that PETSc only guarantees that operations do what you expect if the
   * functions @p VecAssemblyBegin and @p VecAssemblyEnd have been called
   * after vector assembly. Therefore, you need to call Vector::compress()
   * before you actually use the vector.
   *
   * @ingroup PETScWrappers
   */
  class VectorBase : public Subscriptor
  {
  public:
    /**
     * Declare some of the standard types used in all containers. These types
     * parallel those in the <tt>C++</tt> standard libraries
     * <tt>vector<...></tt> class.
     */
    using value_type      = PetscScalar;
    using real_type       = PetscReal;
    using size_type       = types::global_dof_index;
    using reference       = internal::VectorReference;
    using const_reference = const internal::VectorReference;

    /**
     * Default constructor. It doesn't do anything, derived classes will have
     * to initialize the data.
     */
    VectorBase();

    /**
     * Copy constructor. Sets the dimension to that of the given vector, and
     * copies all elements.
     */
    VectorBase(const VectorBase &v);

    /**
     * Initialize a Vector from a PETSc Vec object. Note that we do not copy
     * the vector and we do not obtain ownership, so we do not destroy the
     * PETSc object in the destructor.
     */
    explicit VectorBase(const Vec &v);

    /**
     * The copy assignment operator is deleted to avoid accidental usage with
     * unexpected behavior.
     */
    VectorBase &
    operator=(const VectorBase &) = delete;

    /**
     * Destructor.
     */
    virtual ~VectorBase() override;

    /**
     * Release all memory and return to a state just like after having called
     * the default constructor.
     */
    virtual void
    clear();

    /**
     * Compress the underlying representation of the PETSc object, i.e. flush
     * the buffers of the vector object if it has any. This function is
     * necessary after writing into a vector element-by-element and before
     * anything else can be done on it.
     *
     * See
     * @ref GlossCompress "Compressing distributed objects"
     * for more information.
     */
    void
    compress(const VectorOperation::values operation);

    /**
     * Set all components of the vector to the given number @p s. Simply pass
     * this down to the individual block objects, but we still need to declare
     * this function to make the example given in the discussion about making
     * the constructor explicit work.
     *
     *
     * Since the semantics of assigning a scalar to a vector are not
     * immediately clear, this operator should really only be used if you want
     * to set the entire vector to zero. This allows the intuitive notation
     * <tt>v=0</tt>. Assigning other values is deprecated and may be
     * disallowed in the future.
     */
    VectorBase &
    operator=(const PetscScalar s);

    /**
     * Test for equality. This function assumes that the present vector and
     * the one to compare with have the same size already, since comparing
     * vectors of different sizes makes not much sense anyway.
     */
    bool
    operator==(const VectorBase &v) const;

    /**
     * Test for inequality. This function assumes that the present vector and
     * the one to compare with have the same size already, since comparing
     * vectors of different sizes makes not much sense anyway.
     */
    bool
    operator!=(const VectorBase &v) const;

    /**
     * Return the global dimension of the vector.
     */
    size_type
    size() const;

    /**
     * Return the local dimension of the vector, i.e. the number of elements
     * stored on the present MPI process. For sequential vectors, this number
     * is the same as size(), but for parallel vectors it may be smaller.
     *
     * To figure out which elements exactly are stored locally, use
     * local_range() or locally_owned_elements().
     *
     * @deprecated use locally_owned_size() instead.
     */
    DEAL_II_DEPRECATED_EARLY
    size_type
    local_size() const;

    /**
     * Return the local dimension of the vector, i.e. the number of elements
     * stored on the present MPI process. For sequential vectors, this number
     * is the same as size(), but for parallel vectors it may be smaller.
     *
     * To figure out which elements exactly are stored locally, use
     * local_range() or locally_owned_elements().
     */
    size_type
    locally_owned_size() const;

    /**
     * Return a pair of indices indicating which elements of this vector are
     * stored locally. The first number is the index of the first element
     * stored, the second the index of the one past the last one that is
     * stored locally. If this is a sequential vector, then the result will be
     * the pair (0,N), otherwise it will be a pair (i,i+n), where
     * <tt>n=locally_owned_size()</tt>.
     */
    std::pair<size_type, size_type>
    local_range() const;

    /**
     * Return whether @p index is in the local range or not, see also
     * local_range().
     */
    bool
    in_local_range(const size_type index) const;

    /**
     * Return an index set that describes which elements of this vector are
     * owned by the current processor. Note that this index set does not
     * include elements this vector may store locally as ghost elements but
     * that are in fact owned by another processor. As a consequence, the
     * index sets returned on different processors if this is a distributed
     * vector will form disjoint sets that add up to the complete index set.
     * Obviously, if a vector is created on only one processor, then the
     * result would satisfy
     * @code
     *   vec.locally_owned_elements() == complete_index_set (vec.size())
     * @endcode
     */
    IndexSet
    locally_owned_elements() const;

    /**
     * Return if the vector contains ghost elements.
     *
     * @see
     * @ref GlossGhostedVector "vectors with ghost elements"
     */
    bool
    has_ghost_elements() const;

    /**
     * This function only exists for compatibility with the @p
     * LinearAlgebra::distributed::Vector class and does nothing: this class
     * implements ghost value updates in a different way that is a better fit
     * with the underlying PETSc vector object.
     */
    void
    update_ghost_values() const;

    /**
     * Provide access to a given element, both read and write.
     */
    reference
    operator()(const size_type index);

    /**
     * Provide read-only access to an element.
     */
    PetscScalar
    operator()(const size_type index) const;

    /**
     * Provide access to a given element, both read and write.
     *
     * Exactly the same as operator().
     */
    reference operator[](const size_type index);

    /**
     * Provide read-only access to an element.
     *
     * Exactly the same as operator().
     */
    PetscScalar operator[](const size_type index) const;

    /**
     * A collective set operation: instead of setting individual elements of a
     * vector, this function allows to set a whole set of elements at once.
     * The indices of the elements to be set are stated in the first argument,
     * the corresponding values in the second.
     *
     * @deprecated Use import() instead.
     */
    DEAL_II_DEPRECATED
    void
    set(const std::vector<size_type> &  indices,
        const std::vector<PetscScalar> &values);

    /**
     * Instead of getting individual elements of a vector via operator(),
     * this function allows getting a whole set of elements at once. The
     * indices of the elements to be read are stated in the first argument, the
     * corresponding values are returned in the second.
     *
     * If the current vector is called @p v, then this function is the equivalent
     * to the code
     * @code
     *   for (unsigned int i=0; i<indices.size(); ++i)
     *     values[i] = v[indices[i]];
     * @endcode
     *
     * @pre The sizes of the @p indices and @p values arrays must be identical.
     */
    void
    extract_subvector_to(const std::vector<size_type> &indices,
                         std::vector<PetscScalar> &    values) const;

    /**
     * Instead of getting individual elements of a vector via operator(),
     * this function allows getting a whole set of elements at once. In
     * contrast to the previous function, this function obtains the
     * indices of the elements by dereferencing all elements of the iterator
     * range provided by the first two arguments, and puts the vector
     * values into memory locations obtained by dereferencing a range
     * of iterators starting at the location pointed to by the third
     * argument.
     *
     * If the current vector is called @p v, then this function is the equivalent
     * to the code
     * @code
     *   ForwardIterator indices_p = indices_begin;
     *   OutputIterator  values_p  = values_begin;
     *   while (indices_p != indices_end)
     *   {
     *     *values_p = v[*indices_p];
     *     ++indices_p;
     *     ++values_p;
     *   }
     * @endcode
     *
     * @pre It must be possible to write into as many memory locations
     *   starting at @p values_begin as there are iterators between
     *   @p indices_begin and @p indices_end.
     */
    template <typename ForwardIterator, typename OutputIterator>
    void
    extract_subvector_to(const ForwardIterator indices_begin,
                         const ForwardIterator indices_end,
                         OutputIterator        values_begin) const;

    /**
     * A collective add operation: This function adds a whole set of values
     * stored in @p values to the vector components specified by @p indices.
     */
    void
    add(const std::vector<size_type> &  indices,
        const std::vector<PetscScalar> &values);

    /**
     * This is a second collective add operation. As a difference, this
     * function takes a deal.II vector of values.
     */
    void
    add(const std::vector<size_type> &       indices,
        const ::dealii::Vector<PetscScalar> &values);

    /**
     * Take an address where <tt>n_elements</tt> are stored contiguously and
     * add them into the vector. Handles all cases which are not covered by
     * the other two <tt>add()</tt> functions above.
     */
    void
    add(const size_type    n_elements,
        const size_type *  indices,
        const PetscScalar *values);

    /**
     * Return the scalar product of two vectors. The vectors must have the
     * same size.
     *
     * For complex valued vector, this gives$\left(v^\ast,vec\right)$.
     */
    PetscScalar operator*(const VectorBase &vec) const;

    /**
     * Return the square of the $l_2$-norm.
     */
    real_type
    norm_sqr() const;

    /**
     * Return the mean value of the elements of this vector.
     */
    PetscScalar
    mean_value() const;

    /**
     * $l_1$-norm of the vector. The sum of the absolute values.
     *
     * @note In complex-valued PETSc priori to 3.7.0 this norm is implemented
     * as the sum of absolute values of real and imaginary parts of elements
     * of a complex vector.
     */
    real_type
    l1_norm() const;

    /**
     * $l_2$-norm of the vector.  The square root of the sum of the squares of
     * the elements.
     */
    real_type
    l2_norm() const;

    /**
     * $l_p$-norm of the vector. The pth root of the sum of the pth powers of
     * the absolute values of the elements.
     */
    real_type
    lp_norm(const real_type p) const;

    /**
     * $l_\infty$-norm of the vector. Return the value of the vector element
     * with the maximum absolute value.
     */
    real_type
    linfty_norm() const;

    /**
     * Performs a combined operation of a vector addition and a subsequent
     * inner product, returning the value of the inner product. In other
     * words, the result of this function is the same as if the user called
     * @code
     * this->add(a, V);
     * return_value = *this * W;
     * @endcode
     *
     * The reason this function exists is for compatibility with deal.II's own
     * vector classes which can implement this functionality with less memory
     * transfer. However, for PETSc vectors such a combined operation is not
     * natively supported and thus the cost is completely equivalent as
     * calling the two methods separately.
     *
     * For complex-valued vectors, the scalar product in the second step is
     * implemented as
     * $\left<v,w\right>=\sum_i v_i \bar{w_i}$.
     */
    PetscScalar
    add_and_dot(const PetscScalar a, const VectorBase &V, const VectorBase &W);

    /**
     * Return the value of the vector element with the largest negative value.
     *
     * @deprecated This function has been deprecated to improve compatibility
     * with other classes inheriting from VectorSpaceVector. If you need to
     * use this functionality then use the PETSc function VecMin instead.
     */
    DEAL_II_DEPRECATED
    real_type
    min() const;

    /**
     * Return the value of the vector element with the largest positive value.
     *
     * @deprecated This function has been deprecated to improve compatibility
     * with other classes inheriting from VectorSpaceVector. If you need to
     * use this functionality then use the PETSc function VecMax instead.
     */
    DEAL_II_DEPRECATED
    real_type
    max() const;

    /**
     * Return whether the vector contains only elements with value zero. This
     * is a collective operation. This function is expensive, because
     * potentially all elements have to be checked.
     */
    bool
    all_zero() const;

    /**
     * Return @p true if the vector has no negative entries, i.e. all entries
     * are zero or positive. This function is used, for example, to check
     * whether refinement indicators are really all positive (or zero).
     *
     * @deprecated This function has been deprecated to improve compatibility
     * with other classes inheriting from VectorSpaceVector.
     */
    DEAL_II_DEPRECATED
    bool
    is_non_negative() const;

    /**
     * Multiply the entire vector by a fixed factor.
     */
    VectorBase &
    operator*=(const PetscScalar factor);

    /**
     * Divide the entire vector by a fixed factor.
     */
    VectorBase &
    operator/=(const PetscScalar factor);

    /**
     * Add the given vector to the present one.
     */
    VectorBase &
    operator+=(const VectorBase &V);

    /**
     * Subtract the given vector from the present one.
     */
    VectorBase &
    operator-=(const VectorBase &V);

    /**
     * Addition of @p s to all components. Note that @p s is a scalar and not
     * a vector.
     */
    void
    add(const PetscScalar s);

    /**
     * Simple addition of a multiple of a vector, i.e. <tt>*this += a*V</tt>.
     */
    void
    add(const PetscScalar a, const VectorBase &V);

    /**
     * Multiple addition of scaled vectors, i.e. <tt>*this += a*V+b*W</tt>.
     */
    void
    add(const PetscScalar a,
        const VectorBase &V,
        const PetscScalar b,
        const VectorBase &W);

    /**
     * Scaling and simple vector addition, i.e. <tt>*this = s*(*this)+V</tt>.
     */
    void
    sadd(const PetscScalar s, const VectorBase &V);

    /**
     * Scaling and simple addition, i.e. <tt>*this = s*(*this)+a*V</tt>.
     */
    void
    sadd(const PetscScalar s, const PetscScalar a, const VectorBase &V);

    /**
     * Scale each element of this vector by the corresponding element in the
     * argument. This function is mostly meant to simulate multiplication (and
     * immediate re-assignment) by a diagonal scaling matrix.
     */
    void
    scale(const VectorBase &scaling_factors);

    /**
     * Assignment <tt>*this = a*V</tt>.
     */
    void
    equ(const PetscScalar a, const VectorBase &V);

    /**
     * Prints the PETSc vector object values using PETSc internal vector
     * viewer function <tt>VecView</tt>. The default format prints the
     * vector's contents, including indices of vector elements. For other
     * valid view formats, consult
     * http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/VecView.html
     */
    void
    write_ascii(const PetscViewerFormat format = PETSC_VIEWER_DEFAULT);

    /**
     * Print to a stream. @p precision denotes the desired precision with
     * which values shall be printed, @p scientific whether scientific
     * notation shall be used. If @p across is @p true then the vector is
     * printed in a line, while if @p false then the elements are printed on a
     * separate line each.
     */
    void
    print(std::ostream &     out,
          const unsigned int precision  = 3,
          const bool         scientific = true,
          const bool         across     = true) const;

    /**
     * Swap the contents of this vector and the other vector @p v. One could
     * do this operation with a temporary variable and copying over the data
     * elements, but this function is significantly more efficient since it
     * only swaps the pointers to the data of the two vectors and therefore
     * does not need to allocate temporary storage and move data around.
     *
     * This function is analogous to the @p swap function of all C++
     * standard containers. Also, there is a global function
     * <tt>swap(u,v)</tt> that simply calls <tt>u.swap(v)</tt>, again in
     * analogy to standard functions.
     */
    void
    swap(VectorBase &v);

    /**
     * Conversion operator to gain access to the underlying PETSc type. If you
     * do this, you cut this class off some information it may need, so this
     * conversion operator should only be used if you know what you do. In
     * particular, it should only be used for read-only operations into the
     * vector.
     */
    operator const Vec &() const;

    /**
     * Estimate for the memory consumption (not implemented for this class).
     */
    std::size_t
    memory_consumption() const;

    /**
     * Return a reference to the MPI communicator object in use with this
     * object.
     */
    virtual const MPI_Comm &
    get_mpi_communicator() const;

  protected:
    /**
     * A generic vector object in PETSc. The actual type, a sequential vector,
     * is set in the constructor.
     */
    Vec vector;

    /**
     * Denotes if this vector has ghost indices associated with it. This means
     * that at least one of the processes in a parallel program has at least
     * one ghost index.
     */
    bool ghosted;

    /**
     * This vector contains the global indices of the ghost values. The
     * location in this vector denotes the local numbering, which is used in
     * PETSc.
     */
    IndexSet ghost_indices;

    /**
     * Store whether the last action was a write or add operation. This
     * variable is @p mutable so that the accessor classes can write to it,
     * even though the vector object they refer to is constant.
     */
    mutable VectorOperation::values last_action;

    // Make the reference class a friend.
    friend class internal::VectorReference;

    /**
     * Specifies if the vector is the owner of the PETSc Vec. This is true if
     * it got created by this class and determines if it gets destroyed in
     * the destructor.
     */
    bool obtained_ownership;

    /**
     * Collective set or add operation: This function is invoked by the
     * collective @p set and @p add with the @p add_values flag set to the
     * corresponding value.
     */
    void
    do_set_add_operation(const size_type    n_elements,
                         const size_type *  indices,
                         const PetscScalar *values,
                         const bool         add_values);
  };



  // ------------------- inline and template functions --------------

  /**
   * Global function @p swap which overloads the default implementation of the
   * C++ standard library which uses a temporary object. The function simply
   * exchanges the data of the two vectors.
   *
   * @relatesalso PETScWrappers::VectorBase
   */
  inline void
  swap(VectorBase &u, VectorBase &v)
  {
    u.swap(v);
  }

#    ifndef DOXYGEN
  namespace internal
  {
    inline VectorReference::VectorReference(const VectorBase &vector,
                                            const size_type   index)
      : vector(vector)
      , index(index)
    {}


    inline const VectorReference &
    VectorReference::operator=(const VectorReference &r) const
    {
      // as explained in the class
      // documentation, this is not the copy
      // operator. so simply pass on to the
      // "correct" assignment operator
      *this = static_cast<PetscScalar>(r);

      return *this;
    }



    inline VectorReference &
    VectorReference::operator=(const VectorReference &r)
    {
      // as explained in the class
      // documentation, this is not the copy
      // operator. so simply pass on to the
      // "correct" assignment operator
      *this = static_cast<PetscScalar>(r);

      return *this;
    }



    inline const VectorReference &
    VectorReference::operator=(const PetscScalar &value) const
    {
      Assert((vector.last_action == VectorOperation::insert) ||
               (vector.last_action == VectorOperation::unknown),
             ExcWrongMode(VectorOperation::insert, vector.last_action));

      Assert(!vector.has_ghost_elements(), ExcGhostsPresent());

      const PetscInt petsc_i = index;

      const PetscErrorCode ierr =
        VecSetValues(vector, 1, &petsc_i, &value, INSERT_VALUES);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      vector.last_action = VectorOperation::insert;

      return *this;
    }



    inline const VectorReference &
    VectorReference::operator+=(const PetscScalar &value) const
    {
      Assert((vector.last_action == VectorOperation::add) ||
               (vector.last_action == VectorOperation::unknown),
             ExcWrongMode(VectorOperation::add, vector.last_action));

      Assert(!vector.has_ghost_elements(), ExcGhostsPresent());

      vector.last_action = VectorOperation::add;

      // we have to do above actions in any
      // case to be consistent with the MPI
      // communication model (see the
      // comments in the documentation of
      // PETScWrappers::MPI::Vector), but we
      // can save some work if the addend is
      // zero
      if (value == PetscScalar())
        return *this;

      // use the PETSc function to add something
      const PetscInt       petsc_i = index;
      const PetscErrorCode ierr =
        VecSetValues(vector, 1, &petsc_i, &value, ADD_VALUES);
      AssertThrow(ierr == 0, ExcPETScError(ierr));


      return *this;
    }



    inline const VectorReference &
    VectorReference::operator-=(const PetscScalar &value) const
    {
      Assert((vector.last_action == VectorOperation::add) ||
               (vector.last_action == VectorOperation::unknown),
             ExcWrongMode(VectorOperation::add, vector.last_action));

      Assert(!vector.has_ghost_elements(), ExcGhostsPresent());

      vector.last_action = VectorOperation::add;

      // we have to do above actions in any
      // case to be consistent with the MPI
      // communication model (see the
      // comments in the documentation of
      // PETScWrappers::MPI::Vector), but we
      // can save some work if the addend is
      // zero
      if (value == PetscScalar())
        return *this;

      // use the PETSc function to
      // add something
      const PetscInt       petsc_i     = index;
      const PetscScalar    subtractand = -value;
      const PetscErrorCode ierr =
        VecSetValues(vector, 1, &petsc_i, &subtractand, ADD_VALUES);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      return *this;
    }



    inline const VectorReference &
    VectorReference::operator*=(const PetscScalar &value) const
    {
      Assert((vector.last_action == VectorOperation::insert) ||
               (vector.last_action == VectorOperation::unknown),
             ExcWrongMode(VectorOperation::insert, vector.last_action));

      Assert(!vector.has_ghost_elements(), ExcGhostsPresent());

      vector.last_action = VectorOperation::insert;

      // we have to do above actions in any
      // case to be consistent with the MPI
      // communication model (see the
      // comments in the documentation of
      // PETScWrappers::MPI::Vector), but we
      // can save some work if the factor is
      // one
      if (value == 1.)
        return *this;

      const PetscInt    petsc_i   = index;
      const PetscScalar new_value = static_cast<PetscScalar>(*this) * value;

      const PetscErrorCode ierr =
        VecSetValues(vector, 1, &petsc_i, &new_value, INSERT_VALUES);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      return *this;
    }



    inline const VectorReference &
    VectorReference::operator/=(const PetscScalar &value) const
    {
      Assert((vector.last_action == VectorOperation::insert) ||
               (vector.last_action == VectorOperation::unknown),
             ExcWrongMode(VectorOperation::insert, vector.last_action));

      Assert(!vector.has_ghost_elements(), ExcGhostsPresent());

      vector.last_action = VectorOperation::insert;

      // we have to do above actions in any
      // case to be consistent with the MPI
      // communication model (see the
      // comments in the documentation of
      // PETScWrappers::MPI::Vector), but we
      // can save some work if the factor is
      // one
      if (value == 1.)
        return *this;

      const PetscInt    petsc_i   = index;
      const PetscScalar new_value = static_cast<PetscScalar>(*this) / value;

      const PetscErrorCode ierr =
        VecSetValues(vector, 1, &petsc_i, &new_value, INSERT_VALUES);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      return *this;
    }



    inline PetscReal
    VectorReference::real() const
    {
#      ifndef PETSC_USE_COMPLEX
      return static_cast<PetscScalar>(*this);
#      else
      return PetscRealPart(static_cast<PetscScalar>(*this));
#      endif
    }



    inline PetscReal
    VectorReference::imag() const
    {
#      ifndef PETSC_USE_COMPLEX
      return PetscReal(0);
#      else
      return PetscImaginaryPart(static_cast<PetscScalar>(*this));
#      endif
    }

  } // namespace internal

  inline bool
  VectorBase::in_local_range(const size_type index) const
  {
    PetscInt             begin, end;
    const PetscErrorCode ierr =
      VecGetOwnershipRange(static_cast<const Vec &>(vector), &begin, &end);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return ((index >= static_cast<size_type>(begin)) &&
            (index < static_cast<size_type>(end)));
  }


  inline IndexSet
  VectorBase::locally_owned_elements() const
  {
    IndexSet is(size());

    // PETSc only allows for contiguous local ranges, so this is simple
    const std::pair<size_type, size_type> x = local_range();
    is.add_range(x.first, x.second);
    return is;
  }



  inline bool
  VectorBase::has_ghost_elements() const
  {
    return ghosted;
  }



  inline void
  VectorBase::update_ghost_values() const
  {}



  inline internal::VectorReference
  VectorBase::operator()(const size_type index)
  {
    return internal::VectorReference(*this, index);
  }



  inline PetscScalar
  VectorBase::operator()(const size_type index) const
  {
    return static_cast<PetscScalar>(internal::VectorReference(*this, index));
  }



  inline internal::VectorReference VectorBase::operator[](const size_type index)
  {
    return operator()(index);
  }



  inline PetscScalar VectorBase::operator[](const size_type index) const
  {
    return operator()(index);
  }

  inline const MPI_Comm &
  VectorBase::get_mpi_communicator() const
  {
    static MPI_Comm comm;
    PetscObjectGetComm(reinterpret_cast<PetscObject>(vector), &comm);
    return comm;
  }

  inline void
  VectorBase::extract_subvector_to(const std::vector<size_type> &indices,
                                   std::vector<PetscScalar> &    values) const
  {
    Assert(indices.size() <= values.size(),
           ExcDimensionMismatch(indices.size(), values.size()));
    extract_subvector_to(indices.begin(), indices.end(), values.begin());
  }

  template <typename ForwardIterator, typename OutputIterator>
  inline void
  VectorBase::extract_subvector_to(const ForwardIterator indices_begin,
                                   const ForwardIterator indices_end,
                                   OutputIterator        values_begin) const
  {
    const PetscInt n_idx = static_cast<PetscInt>(indices_end - indices_begin);
    if (n_idx == 0)
      return;

    // if we are dealing
    // with a parallel vector
    if (ghosted)
      {
        // there is the possibility
        // that the vector has
        // ghost elements. in that
        // case, we first need to
        // figure out which
        // elements we own locally,
        // then get a pointer to
        // the elements that are
        // stored here (both the
        // ones we own as well as
        // the ghost elements). in
        // this array, the locally
        // owned elements come
        // first followed by the
        // ghost elements whose
        // position we can get from
        // an index set
        PetscInt       begin, end;
        PetscErrorCode ierr = VecGetOwnershipRange(vector, &begin, &end);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        Vec locally_stored_elements = nullptr;
        ierr = VecGhostGetLocalForm(vector, &locally_stored_elements);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        PetscInt lsize;
        ierr = VecGetSize(locally_stored_elements, &lsize);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        PetscScalar *ptr;
        ierr = VecGetArray(locally_stored_elements, &ptr);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        for (PetscInt i = 0; i < n_idx; ++i)
          {
            const unsigned int index = *(indices_begin + i);
            if (index >= static_cast<unsigned int>(begin) &&
                index < static_cast<unsigned int>(end))
              {
                // local entry
                *(values_begin + i) = *(ptr + index - begin);
              }
            else
              {
                // ghost entry
                const unsigned int ghostidx =
                  ghost_indices.index_within_set(index);

                AssertIndexRange(ghostidx + end - begin, lsize);
                *(values_begin + i) = *(ptr + ghostidx + end - begin);
              }
          }

        ierr = VecRestoreArray(locally_stored_elements, &ptr);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        ierr = VecGhostRestoreLocalForm(vector, &locally_stored_elements);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
      }
    // if the vector is local or the
    // caller, then simply access the
    // element we are interested in
    else
      {
        PetscInt       begin, end;
        PetscErrorCode ierr = VecGetOwnershipRange(vector, &begin, &end);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        PetscScalar *ptr;
        ierr = VecGetArray(vector, &ptr);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        for (PetscInt i = 0; i < n_idx; ++i)
          {
            const unsigned int index = *(indices_begin + i);

            Assert(index >= static_cast<unsigned int>(begin) &&
                     index < static_cast<unsigned int>(end),
                   ExcInternalError());

            *(values_begin + i) = *(ptr + index - begin);
          }

        ierr = VecRestoreArray(vector, &ptr);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
      }
  }

#    endif // DOXYGEN
} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_WITH_PETSC

#endif
/*---------------------------- petsc_vector_base.h --------------------------*/

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_psblas_vector_h
#define dealii_psblas_vector_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/types.h>

#include <deal.II/lac/read_vector.h>
#include <deal.II/lac/vector_operation.h>
#include <deal.II/lac/vector_type_traits.h>

#include <memory.h>

#include <cstddef>

#ifdef DEAL_II_WITH_PSBLAS

#  include <psb_base_cbind.h>
#  include <psb_c_base.h>
#  include <psb_c_dbase.h>

DEAL_II_NAMESPACE_OPEN

namespace PSCToolkitWrappers
{

  namespace internal
  {
    /*
     * Custom deleter for PSBLAS descriptor.
     */
    struct PSBLASDescriptorDeleter
    {
      void
      operator()(psb_c_descriptor *p) const
      {
        if (p)
          psb_c_cdfree(p);
      }
    };

    /**
     * Enum to indicate the state of the vector (building or assembled).
     * TODO[MF]: use the same also when I'll introduce the matrix class. Move to
     * a common .h file?
     */

    enum State
    {
      /**
       * State entered after the default constructor, before any allocation.
       * In this state, no operations are possible.
       */
      Default,
      /**
       * State entered after the first allocation, and before the first
       * assembly; in this state it is possible to add communication
       * requirements among different processes.
       */
      Build,
      /*
       * State entered after the assembly; computations such as matrix-vector
       * products, are only possible in this state.
       */
      Assembled
    };

  } // namespace internal

  class Vector : public ReadVector<double>
  {
  private:
    /**
     * This class provides a wrappers for accessing psblas vector elements.
     */
    class VectorReference
    {
    private:
      using size_type = types::global_dof_index;

      using value_type = double;

      /**
       * Constructor.
       */
      VectorReference(Vector &vector, const size_type index)
        : vector(vector)
        , index(index)
      {}

    public:
      /**
       * Set the referenced element of the vector to <tt>s</tt>.
       */
      const VectorReference &
      operator=(const value_type &s) const
      {
        Assert(!vector.has_ghost_elements(), ExcGhostsPresent());
        Assert(
          vector.owned_elements.is_element(index),
          ExcMessage(
            "You are trying to write to an element of the vector that is not "
            "locally owned. This is not allowed for the current interface to"
            " PSBLAS vectors."));

        // Make sure the operation is consistent with the last one
        Assert(vector.last_action == VectorOperation::insert ||
                 vector.last_action == VectorOperation::unknown,
               ExcWrongMode(VectorOperation::insert, vector.last_action));

        std::vector<size_type>  idx{index};
        std::vector<value_type> value{s};
        vector.set(idx, value);
        vector.last_action = VectorOperation::insert;
        return *this;
      }

      /**
       * Add <tt>s</tt> to the referenced element of the vector.
       */
      const VectorReference &
      operator+=(const value_type &s) const
      {
        Assert(!vector.has_ghost_elements(), ExcGhostsPresent());
        Assert(vector.last_action == VectorOperation::add ||
                 vector.last_action == VectorOperation::unknown,
               ExcWrongMode(VectorOperation::add, vector.last_action));

        vector.last_action = VectorOperation::add;

        // First check for early return
        if (s == 0.)
          return *this;

        std::vector<size_type>  idx{index};
        std::vector<value_type> value{s};
        vector.add(idx, value);
        return *this;
      }

      /**
       * Subtract <tt>s</tt> to the referenced element of the vector.
       */
      const VectorReference &
      operator-=(const value_type &s) const
      {
        Assert(!vector.has_ghost_elements(), ExcGhostsPresent());
        Assert(vector.last_action == VectorOperation::add ||
                 vector.last_action == VectorOperation::unknown,
               ExcWrongMode(VectorOperation::add, vector.last_action));

        vector.last_action = VectorOperation::add;

        // First check for early return
        if (s == 0.)
          return *this;

        std::vector<size_type>  idx{index};
        std::vector<value_type> value{-s};
        vector.add(idx, value);
        return *this;
      }

      /**
       * Multiply <tt>s</tt> to the referenced element of the vector.
       */
      const VectorReference &
      operator*=(const value_type &s) const
      {
        Assert(!vector.has_ghost_elements(), ExcGhostsPresent());
        Assert((vector.last_action == VectorOperation::insert) ||
                 (vector.last_action == VectorOperation::unknown),
               ExcWrongMode(VectorOperation::insert, vector.last_action));

        vector.last_action = VectorOperation::insert;
        if (s == 1.)
          return *this;

        std::vector<size_type>  idx{index};
        value_type              new_value = static_cast<value_type>(*this) * s;
        std::vector<value_type> value{new_value};
        vector.set(idx, value);

        return *this;
      }

      /**
       * Divide <tt>s</tt> to the referenced element of the vector.
       */
      const VectorReference &
      operator/=(const value_type &s) const
      {
        Assert(!vector.has_ghost_elements(), ExcGhostsPresent());

        Assert((vector.last_action == VectorOperation::insert) ||
                 (vector.last_action == VectorOperation::unknown),
               ExcWrongMode(VectorOperation::insert, vector.last_action));
        std::vector<size_type>  idx{index};
        value_type              new_value = static_cast<value_type>(*this) / s;
        std::vector<value_type> value{new_value};
        vector.set(idx, value);
        vector.last_action = VectorOperation::insert;
        return *this;
      }

      /*
       * Convert the reference to an actual value, i.e. return the value of
       * the referenced element of the vector.
       */
      operator value_type() const
      {
        AssertIndexRange(index, vector.size());
        if (vector.ghosted)
          {
            AssertThrow(vector.ghost_indices.is_element(index) ||
                          vector.owned_elements.is_element(index),
                        ExcMessage(
                          "You are trying to access an element of a vector "
                          "that is neither a locally owned element nor a "
                          "ghost element of the vector."));
          }
        else
          {
            AssertThrow(
              vector.owned_elements.is_element(index),
              ExcAccessToNonlocalElement(index,
                                         *vector.owned_elements.begin(),
                                         (*vector.owned_elements.begin() +
                                          vector.locally_owned_size())));
          }
        return psb_c_dgetelem(vector.psblas_vector,
                              index,
                              vector.psblas_descriptor.get());
      };

    private:
      Vector &vector;

      const size_type index;

      friend class Vector;
    };

  public:
    using size_type = dealii::types::global_dof_index;

    using value_type = double;

    /**
     * Exception
     */
    DeclException1(ExcInitializePSBLASVector,
                   int,
                   << "An error with error number " << arg1
                   << " occurred while initializing a PSBLAS vector.");

    /**
     * Exception
     */
    DeclException1(ExcFreePSBLASVector,
                   int,
                   << "An error with error number " << arg1
                   << " occurred while freeing a PSBLAS vector.");

    /**
     * Exception
     */
    DeclException1(ExcInitializePSBLASDescriptor,
                   int,
                   << "An error with error number " << arg1
                   << " occurred while initializing a PSBLAS descriptor.");

    /**
     * Exception
     */
    DeclException1(ExcAssemblePSBLASVector,
                   int,
                   << "An error with error number " << arg1
                   << " occurred while assembling a PSBLAS vector.");

    /**
     * Exception
     */
    DeclException1(ExcAssemblePSBLASDescriptor,
                   int,
                   << "An error with error number " << arg1
                   << " occurred while assembling a PSBLAS descriptor.");


    /**
     * Exception
     */
    DeclException1(ExcInvalidState,
                   int,
                   << "Vector's state is invalid. It is in "
                   << (arg1 == 0 ? "Default" :
                       arg1 == 1 ? "Build" :
                                   "Assembled")
                   << " state. Did you forget to call reinit() or compress()?");

    /**
     * Exception
     */
    DeclException1(ExcInvalidStateBuild,
                   int,
                   << "Vector's state is invalid. It should be in "
                   << "Build state but it is in "
                   << (arg1 == 0 ? "Default" : "Assembled")
                   << " state. Did you forget to call reinit()?");

    /**
     * Exception
     */
    DeclException1(ExcInvalidStateAssembled,
                   int,
                   << "Vector's state is invalid. It should be in "
                   << "Assembled state but it is in "
                   << (arg1 == 0 ? "Default" : "Build")
                   << " state. Did you forget to call compress()?");

    /**
     * Exception
     */
    DeclExceptionMsg(
      ExcInvalidDefault,
      "Vector's state is invalid. It should be in "
      "Build or Assembled state but it is in Default state. Did you forget"
      " to reinit it()?");

    /**
     * Exception
     */
    DeclException2(ExcCallingPSBLASFunction,
                   int,
                   std::string,
                   << "An error with error number " << arg1
                   << " occurred while calling a PSBLAS function: " << arg2
                   << std::endl);

    /**
     * Exception
     */
    DeclException1(ExcAXPBY,
                   int,
                   << "An error with error number " << arg1
                   << " occurred while performing the AXPBY operation.");

    /**
     * Exception
     */
    DeclException1(ExcInsertionInPSBLASVector,
                   int,
                   << "An error with error number " << arg1
                   << " occurred while inserting values into a PSBLAS vector.");

    /**
     * Exception
     */
    DeclException2(ExcWrongMode,
                   int,
                   int,
                   << "You tried to do a "
                   << (arg1 == 1 ? "'set'" : (arg1 == 2 ? "'add'" : "???"))
                   << " operation but the vector is currently in "
                   << (arg2 == 1 ? "'set'" : (arg2 == 2 ? "'add'" : "???"))
                   << " mode. You first have to call 'compress()'.");

    /**
     * Exception
     */
    DeclException3(
      ExcAccessToNonlocalElement,
      int,
      int,
      int,
      << "You tried to access element " << arg1
      << " of a distributed vector, but only elements in range [" << arg2 << ','
      << arg3 << "] are stored locally and can be accessed."
      << "\n\n"
      << "A common source for this kind of problem is that you "
      << "are passing a 'fully distributed' vector into a function "
      << "that needs read access to vector elements that correspond "
      << "to degrees of freedom on ghost cells (or at least to "
      << "'locally active' degrees of freedom that are not also "
      << "'locally owned'). You need to pass a vector that has these "
      << "elements as ghost entries.");

    /**
     *Default constructor. Generates an empty (zero-size) vector.
     */
    Vector();

    /**
     * Copy constructor. Sets the dimension to that of the given vector, and
     * copies all elements.
     */
    Vector(const Vector &);

    /**
     * Construct a new parallel PSBLAS vector without ghost elements from an
     * IndexSet.
     */
    explicit Vector(const IndexSet &local_partitioning,
                    const MPI_Comm  communicator);

    /**
     * Construct a new parallel ghosted PSBLAS vector from IndexSets.
     *
     * @note This operation always creates a ghosted vector, which is considered
     * read-only.
     *
     * @see
     * @ref GlossGhostedVector "vectors with ghost elements"
     */
    Vector(const IndexSet &local_partitioning,
           const IndexSet &ghost_indices,
           const MPI_Comm  communicator);

    /**
     * Destructor. Internally, its frees the PSBLAS vector and descriptor.
     */
    ~Vector();

    /**
     * Reinit as a vector without ghost elements. If omit_zeroing_entries is
     * false, the vector is filled by zeros. Otherwise, the elements are left an
     * unspecified state.
     *
     * @see
     * @ref GlossGhostedVector "vectors with ghost elements"
     */
    void
    reinit(const IndexSet &local_partitioning,
           const MPI_Comm  communicator,
           const bool      omit_zeroing_entries = false);

    /**
     * Change the dimension to that of the vector @p v, and also take over
     * the partitioning into local sizes as well as the MPI communicator.
     * The same applies as for the other @p reinit function.
     *
     * The elements of @p v are not copied, i.e. this function is the same
     * as calling <tt>reinit(v.size(), v.locally_owned_size(),
     * omit_zeroing_entries)</tt>.
     */
    void
    reinit(const Vector &v, const bool omit_zeroing_entries = false);

    /**
     * Construct a new parallel ghosted PSBLAS vector from IndexSets.
     *
     * Note that the @p ghost IndexSet may be empty and that any indices
     * already contained in @p local are ignored during construction. The
     * global indices in ghost are supplied as ghost indices so that they can be
     * read locally.
     *
     * @note This operation always creates a ghosted vector, which is considered
     * read-only.
     */
    void
    reinit(const IndexSet &local_partitioning,
           const IndexSet &ghost_indices,
           const MPI_Comm  communicator);


    /**
     * Copy assignment. The behavior of this class is the same as that of
     * PETScWrappers::MPI::Vector class.
     *
     * Resize the present vector if necessary. Also take
     * over the MPI communicator of v. The semantics of this operator are
     * complex. If the two vectors have the same size, and if either the left or
     * right hand side vector of the assignment (i.e., either the input vector
     * on the right hand side, or the calling vector to the left of the
     * assignment operator) currently has ghost elements, then the left hand
     * side vector will also have ghost values and will consequently be a
     * read-only vector (see also the glossary entry on the issue). Otherwise,
     * the left hand vector will be a writable vector after this operation.
     * These semantics facilitate having a vector with ghost elements on the
     * left hand side of the assignment, and a vector without ghost elements on
     * the right hand side, with the resulting left hand side vector having the
     * correct values in both its locally owned and its ghost elements. On the
     * other hand, if the left hand side vector does not have the correct size
     * yet, or is perhaps an entirely uninitialized vector, then the assignment
     * is simply a copy operation in the usual sense: In that case, if the right
     * hand side has no ghost elements (i.e., is a completely distributed
     * vector), then the left hand side will have no ghost elements either. And
     * if the right hand side has ghost elements (and is consequently
     * read-only), then the left hand side will have these same properties after
     * the operation.
     */
    Vector &
    operator=(const Vector &v);

    /**
     * Return the global size of the vector, i.e. the sum of the local sizes
     * over all MPI processes.
     */
    size_type
    size() const override;

    /**
     * Extract a range of elements all at once.
     */
    virtual void
    extract_subvector_to(
      const ArrayView<const types::global_dof_index> &indices,
      const ArrayView<value_type>                    &elements) const override;

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
                         std::vector<value_type>      &values) const;

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
    extract_subvector_to(ForwardIterator indices_begin,
                         ForwardIterator indices_end,
                         OutputIterator  values_begin) const;


    /**
     * Return the local dimension of the vector, i.e. the number of elements
     * stored on the present MPI process. For sequential vectors, this number
     * is the same as size(), but for parallel vectors it may be smaller.
     *
     * To figure out which elements exactly are stored locally, use
     * locally_owned_elements().
     */
    size_type
    locally_owned_size() const;

    /**
     * A collective set operation: instead of setting individual elements of a
     * vector, this function allows to set a whole set of elements at once.
     * The indices of the elements to be set are stated in the first argument,
     * the corresponding values in the second.
     */
    void
    set(const std::vector<size_type>  &indices,
        const std::vector<value_type> &values);

    /**
     * A collective add operation: This function adds a whole set of values
     * stored in @p values to the vector components specified by @p indices.
     */
    void
    add(const std::vector<size_type>  &indices,
        const std::vector<value_type> &values);

    /**
     * Addition of a multiple of a vector, i.e. <tt>*this += s*V</tt>.
     */
    void
    add(const value_type s, const Vector &V);

    /**
     * Addition of <tt>s</tt> to all components.
     */
    void
    add(const value_type s);

    /**
     * Scale each element of this vector by the corresponding element in the
     * argument. This function is mostly meant to simulate multiplication (and
     * immediate re-assignment) by a diagonal scaling matrix.
     */
    void
    scale(const Vector &v);

    /**
     * Performs a combined operation of a vector addition and a subsequent
     * inner product, returning the value of the inner product. In other
     * words, the result of this function is the same as if the user called
     * @code
     * this->add(a, V);
     * return_value = *this * W;
     * @endcode
     *
     */
    value_type
    add_and_dot(const value_type a, const Vector &v, const Vector &W);

    /*
     * Scaling and vector addition, i.e.  <tt>*this = s*(*this)+V</tt>.
     */
    void
    sadd(const value_type s, const Vector &V);

    /*
     * Scaling and vector addition, i.e.  <tt>*this = s*(*this)+a*V</tt>.
     */
    void
    sadd(const value_type s, const value_type a, const Vector &V);

    /*
     * Assignment *this = a*V.
     */
    void
    equ(const value_type a, const Vector &v);

    /**
     * Provide read-only access to an element.
     */
    value_type
    operator()(const size_type index) const;

    /**
     * Provide read-write access to an element of the vector.
     */
    VectorReference
    operator()(const size_type index);

    /**
     * Provide read-only access to an element.
     */
    value_type
    operator[](const size_type index) const;

    /**
     * Provide read-write access to an element of the vector.
     */
    VectorReference
    operator[](const size_type index);

    /**
     * Dot product of the vector with another vector.
     */
    value_type
    operator*(const Vector &v) const;

    /**
     * Subtract the given vector from the present one.
     */
    Vector &
    operator-=(const Vector &v);

    /**
     * Add the given vector from the present one.
     */
    Vector &
    operator+=(const Vector &v);

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
     * <tt>v=0</tt>.
     */
    Vector &
    operator=(const value_type s);

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
    const IndexSet &
    locally_owned_elements() const;

    /**
     * Return the IndexSet of ghost elements.
     */
    const IndexSet &
    ghost_elements() const;

    /**
     * Return if the vector contains ghost elements.
     *
     * @see
     * @ref GlossGhostedVector "vectors with ghost elements"
     */
    bool
    has_ghost_elements() const;

    /**
     * Update ghosted elements.
     */
    void
    update_ghost_values() const;

    /**
     * Swap the contents of this vector and the other vector @p v. This function only
     * swaps the pointers to the data of the two vectors and therefore does not
     * need to allocate temporary storage and move data around.
     *
     * This function is analogous to the @p swap function of all C++
     * standard containers.
     */
    void
    swap(Vector &v);

    /**
     * Compress the underlying representation of the PSBLAS object, i.e. flush
     * the buffers of the vector object if it has any. This function is
     * necessary after writing into a vector element-by-element and before
     * anything else can be done on it.
     *
     */
    void
    compress(const VectorOperation::values operation);

    /**
     * Return an iterator to the start of the locally owned elements of
     * the vector.
     */
    value_type *
    begin();

    /**
     * Return a constant iterator to the start of the locally owned elements of
     * the vector.
     */
    const value_type *
    begin() const;

    /**
     * Return an iterator to the element past the end of the locally owned
     * elements of the vector.
     */
    value_type *
    end();

    /**
     * Return a constant iterator to the element past the end of the locally
     * owned elements of the vector.
     */
    const value_type *
    end() const;

    MPI_Comm
    get_mpi_communicator() const;

    /**
     * Get the underlying PSBLAS descriptor. Use it only when you know what you
     * are doing.
     */
    psb_c_descriptor *
    get_psblas_descriptor() const;

    /**
     * Get a pointer to the underlying PSBLAS vector. Use it only when you know
     * what you are doing.
     */
    psb_c_dvector *
    get_psblas_vector() const;

    /**
     * Release all memory and return to a state just like after having
     * called the default constructor.
     */
    void
    clear();

    /**
     * $l_\infty$-norm of the vector. Return the value of the vector element
     * with the maximum absolute value.
     */
    Vector::value_type
    linfty_norm() const;

    /**
     * $l_1$-norm of the vector. The sum of the absolute values.
     */
    Vector::value_type
    l1_norm() const;

    /**
     * $l_2$-norm of the vector.  The square root of the sum of the squares of
     * the elements.
     */
    Vector::value_type
    l2_norm() const;

    /**
     * Compute the mean value of the vector.
     */
    Vector::value_type
    mean_value() const;

    /**
     * Return whether the vector contains only elements with value zero. This
     * is a @ref GlossCollectiveOperation "collective operation". This function is expensive, because
     * potentially all elements have to be checked.
     */
    bool
    all_zero() const;

    /**
     * Estimate for the memory consumption (not implemented for this class).
     */
    std::size_t
    memory_consumption() const;

  private:
    /*
     * Pointer to the underlying PSBLAS vector.
     */
    psb_c_dvector *psblas_vector;

    /*
     * Pointer to PSBLAS context.
     */
    psb_c_ctxt *psblas_context;

    /*
     * Pointer to PSBLAS descriptor.
     */
    std::shared_ptr<psb_c_descriptor> psblas_descriptor;

    /**
     * The MPI communicator over which the vector is distributed.
     */
    MPI_Comm communicator;

    /**
     * This IndexSet contains the global indices of the locally owned values.
     */
    IndexSet owned_elements;

    /**
     * This IndexSet contains the global indices of the ghost values.
     */
    IndexSet ghost_indices;

    /**
     * Denotes if this vector has ghost indices associated with it. This means
     * that at least one of the processes in a parallel program has at least
     * one ghost index.
     */
    bool ghosted;

    /**
     * State of the descriptor associated with the vector. Its state can be
     * either default, building or assembled).
     */
    internal::State state;

    /**
     * Member variable storing if the last operation done on the vector was a
     * write or an add operation.
     */
    VectorOperation::values last_action;

    // TODO[MF]: uncomment when the matrix class will be introduced
    // friend class SparseMatrix;
  };


  /* ----------------------------- Inline functions ---------------- */


  inline PSCToolkitWrappers::Vector::value_type *
  PSCToolkitWrappers::Vector::begin()
  {
    return psb_c_dvect_f_get_pnt(psblas_vector);
  }



  inline const PSCToolkitWrappers::Vector::value_type *
  PSCToolkitWrappers::Vector::begin() const
  {
    return psb_c_dvect_f_get_pnt(psblas_vector);
  }



  inline PSCToolkitWrappers::Vector::value_type *
  PSCToolkitWrappers::Vector::end()
  {
    return psb_c_dvect_f_get_pnt(psblas_vector) + locally_owned_size();
  }



  inline const PSCToolkitWrappers::Vector::value_type *
  PSCToolkitWrappers::Vector::end() const
  {
    return psb_c_dvect_f_get_pnt(psblas_vector) + locally_owned_size();
  }



  inline Vector::size_type
  Vector::size() const
  {
    return owned_elements.size();
  }



  inline const IndexSet &
  Vector::locally_owned_elements() const
  {
    return owned_elements;
  }



  inline const IndexSet &
  Vector::ghost_elements() const
  {
    return ghost_indices;
  }



  inline bool
  Vector::has_ghost_elements() const
  {
    return ghosted;
  }



  inline Vector::value_type
  Vector::operator()(const Vector::size_type index) const
  {
    return psb_c_dgetelem(psblas_vector, index, psblas_descriptor.get());
  }



  inline Vector::VectorReference
  Vector::operator()(const size_type index)
  {
    return VectorReference(*this, index);
  }



  inline Vector::value_type
  Vector::operator[](const Vector::size_type index) const
  {
    return operator()(index);
  }



  inline Vector::VectorReference
  Vector::operator[](const size_type index)
  {
    return operator()(index);
  }



  inline void
  Vector::extract_subvector_to(
    const ArrayView<const types::global_dof_index> &indices,
    const ArrayView<double>                        &elements) const
  {
    AssertDimension(indices.size(), elements.size());
    extract_subvector_to(indices.begin(), indices.end(), elements.begin());
  }



  inline void
  Vector::extract_subvector_to(const std::vector<size_type>    &indices,
                               std::vector<Vector::value_type> &values) const
  {
    AssertDimension(indices.size(), values.size());
    extract_subvector_to(indices.begin(), indices.end(), values.begin());
  }



  template <typename ForwardIterator, typename OutputIterator>
  inline void
  Vector::extract_subvector_to(ForwardIterator indices_begin,
                               ForwardIterator indices_end,
                               OutputIterator  output) const
  {
    if (indices_begin == indices_end)
      return;

    if (ghosted)
      {
        types::global_dof_index begin = *owned_elements.begin();
        types::global_dof_index end   = begin + owned_elements.n_elements();

        auto input = indices_begin;
        while (input != indices_end)
          {
            const auto index = static_cast<psb_l_t>(*input);
            AssertThrow((index >= begin && index < end) ||
                          ghost_indices.is_element(index),
                        ExcMessage(
                          "You are trying to access an element of a vector "
                          "that is neither a locally owned element nor a "
                          "ghost element of the vector."));
            *output =
              psb_c_dgetelem(psblas_vector, index, psblas_descriptor.get());


            ++input;
            ++output;
          }
      }
    else
      {
        // no ghost elements, so we can
        // just access the local
        // elements directly
        while (indices_begin != indices_end)
          {
            const size_type index = *indices_begin;
            Assert(owned_elements.is_element(index),
                   ExcMessage("You are accessing elements of a vector without "
                              "ghost elements that are not actually owned by "
                              "this vector. A typical case where this may "
                              "happen is if you are passing a non-ghosted "
                              "(completely distributed) vector to a function "
                              "that expects a vector that stores ghost "
                              "elements for all locally relevant or locally "
                              "active vector entries."));

            *output =
              psb_c_dgetelem(psblas_vector, index, psblas_descriptor.get());

            ++indices_begin;
            ++output;
          }
      }
  }


} // namespace PSCToolkitWrappers

/**
 * Declare PSCToolkitWrappers::Vector as distributed vector.
 */
template <>
struct is_serial_vector<PSCToolkitWrappers::Vector> : std::false_type
{};
DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PSBLAS
#endif

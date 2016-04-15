// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2016 by the deal.II authors
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

#ifndef dealii__petsc_parallel_vector_h
#define dealii__petsc_parallel_vector_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/base/subscriptor.h>
#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/vector.h>
#  include <deal.II/lac/petsc_vector_base.h>
#  include <deal.II/base/index_set.h>

DEAL_II_NAMESPACE_OPEN


/*! @addtogroup PETScWrappers
 *@{
 */
namespace PETScWrappers
{
  /**
   * Namespace for PETSc classes that work in parallel over MPI, such as
   * distributed vectors and matrices.
   *
   * @ingroup PETScWrappers
   * @author Wolfgang Bangerth, 2004
   */
  namespace MPI
  {

    /**
     * Implementation of a parallel vector class based on PETSC and using MPI
     * communication to synchronise distributed operations. All the
     * functionality is actually in the base class, except for the calls to
     * generate a parallel vector. This is possible since PETSc only works on
     * an abstract vector type and internally distributes to functions that do
     * the actual work depending on the actual vector type (much like using
     * virtual functions). Only the functions creating a vector of specific
     * type differ, and are implemented in this particular class.
     *
     *
     * <h3>Parallel communication model</h3>
     *
     * The parallel functionality of PETSc is built on top of the Message
     * Passing Interface (MPI). MPI's communication model is built on
     * collective communications: if one process wants something from another,
     * that other process has to be willing to accept this communication. A
     * process cannot query data from another process by calling a remote
     * function, without that other process expecting such a transaction. The
     * consequence is that most of the operations in the base class of this
     * class have to be called collectively. For example, if you want to
     * compute the l2 norm of a parallel vector, @em all processes across
     * which this vector is shared have to call the @p l2_norm function. If
     * you don't do this, but instead only call the @p l2_norm function on one
     * process, then the following happens: This one process will call one of
     * the collective MPI functions and wait for all the other processes to
     * join in on this. Since the other processes don't call this function,
     * you will either get a time-out on the first process, or, worse, by the
     * time the next a call to a PETSc function generates an MPI message on
     * the other processes, you will get a cryptic message that only a subset
     * of processes attempted a communication. These bugs can be very hard to
     * figure out, unless you are well-acquainted with the communication model
     * of MPI, and know which functions may generate MPI messages.
     *
     * One particular case, where an MPI message may be generated unexpectedly
     * is discussed below.
     *
     *
     * <h3>Accessing individual elements of a vector</h3>
     *
     * PETSc does allow read access to individual elements of a vector, but in
     * the distributed case only to elements that are stored locally. We
     * implement this through calls like <tt>d=vec(i)</tt>. However, if you
     * access an element outside the locally stored range, an exception is
     * generated.
     *
     * In contrast to read access, PETSc (and the respective deal.II wrapper
     * classes) allow to write (or add) to individual elements of vectors,
     * even if they are stored on a different process. You can do this
     * writing, for example, <tt>vec(i)=d</tt> or <tt>vec(i)+=d</tt>, or
     * similar operations. There is one catch, however, that may lead to very
     * confusing error messages: PETSc requires application programs to call
     * the compress() function when they switch from adding, to elements to
     * writing to elements. The reasoning is that all processes might
     * accumulate addition operations to elements, even if multiple processes
     * write to the same elements. By the time we call compress() the next
     * time, all these additions are executed. However, if one process adds to
     * an element, and another overwrites to it, the order of execution would
     * yield non-deterministic behavior if we don't make sure that a
     * synchronisation with compress() happens in between.
     *
     * In order to make sure these calls to compress() happen at the
     * appropriate time, the deal.II wrappers keep a state variable that store
     * which is the presently allowed operation: additions or writes. If it
     * encounters an operation of the opposite kind, it calls compress() and
     * flips the state. This can sometimes lead to very confusing behavior, in
     * code that may for example look like this:
     * @code
     *   PETScWrappers::MPI::Vector vector;
     *   ...
     *                   // do some write operations on the vector
     *   for (unsigned int i=0; i<vector.size(); ++i)
     *     vector(i) = i;
     *
     *                   // do some additions to vector elements, but
     *                   // only for some elements
     *   for (unsigned int i=0; i<vector.size(); ++i)
     *     if (some_condition(i) == true)
     *       vector(i) += 1;
     *
     *                   // do another collective operation
     *   const double norm = vector.l2_norm();
     * @endcode
     *
     * This code can run into trouble: by the time we see the first addition
     * operation, we need to flush the overwrite buffers for the vector, and
     * the deal.II library will do so by calling compress(). However, it will
     * only do so for all processes that actually do an addition -- if the
     * condition is never true for one of the processes, then this one will
     * not get to the actual compress() call, whereas all the other ones do.
     * This gets us into trouble, since all the other processes hang in the
     * call to flush the write buffers, while the one other process advances
     * to the call to compute the l2 norm. At this time, you will get an error
     * that some operation was attempted by only a subset of processes. This
     * behavior may seem surprising, unless you know that write/addition
     * operations on single elements may trigger this behavior.
     *
     * The problem described here may be avoided by placing additional calls
     * to compress(), or making sure that all processes do the same type of
     * operations at the same time, for example by placing zero additions if
     * necessary.
     *
     * @see
     * @ref GlossGhostedVector "vectors with ghost elements"
     *
     * @ingroup PETScWrappers
     * @ingroup Vectors
     * @author Wolfgang Bangerth, 2004
     */
    class Vector : public VectorBase
    {
    public:
      /**
       * Declare type for container size.
       */
      typedef types::global_dof_index size_type;

      /**
       * A variable that indicates whether this vector supports distributed
       * data storage. If true, then this vector also needs an appropriate
       * compress() function that allows communicating recent set or add
       * operations to individual elements to be communicated to other
       * processors.
       *
       * For the current class, the variable equals true, since it does
       * support parallel data storage.
       */
      static const bool supports_distributed_data = true;

      /**
       * Default constructor. Initialize the vector as empty.
       */
      Vector ();

      /**
       * Constructor. Set dimension to @p n and initialize all elements with
       * zero.
       *
       * @arg local_size denotes the size of the chunk that shall be stored on
       * the present process.
       *
       * @arg communicator denotes the MPI communicator over which the
       * different parts of the vector shall communicate
       *
       * The constructor is made explicit to avoid accidents like this:
       * <tt>v=0;</tt>. Presumably, the user wants to set every element of the
       * vector to zero, but instead, what happens is this call:
       * <tt>v=Vector@<number@>(0);</tt>, i.e. the vector is replaced by one
       * of length zero.
       */
      explicit Vector (const MPI_Comm  &communicator,
                       const size_type  n,
                       const size_type  local_size);


      /**
       * Copy-constructor from deal.II vectors. Sets the dimension to that of
       * the given vector, and copies all elements.
       *
       * @arg local_size denotes the size of the chunk that shall be stored on
       * the present process.
       *
       * @arg communicator denotes the MPI communicator over which the
       * different parts of the vector shall communicate
       */
      template <typename Number>
      explicit Vector (const MPI_Comm               &communicator,
                       const dealii::Vector<Number> &v,
                       const size_type               local_size);


      /**
       * Copy-constructor the values from a PETSc wrapper vector class.
       *
       * @arg local_size denotes the size of the chunk that shall be stored on
       * the present process.
       *
       * @arg communicator denotes the MPI communicator over which the
       * different parts of the vector shall communicate
       */
      explicit Vector (const MPI_Comm     &communicator,
                       const VectorBase   &v,
                       const size_type     local_size);

      /**
       * Constructs a new parallel ghosted PETSc vector from an IndexSet. Note
       * that @p local must be contiguous and the global size of the vector is
       * determined by local.size(). The global indices in @p ghost are
       * supplied as ghost indices that can also be read locally.
       *
       * Note that the @p ghost IndexSet may be empty and that any indices
       * already contained in @p local are ignored during construction. That
       * way, the ghost parameter can equal the set of locally relevant
       * degrees of freedom, see step-32.
       *
       * @note This operation always creates a ghosted vector.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      Vector (const IndexSet &local,
              const IndexSet &ghost,
              const MPI_Comm &communicator);

      /**
       * Constructs a new parallel PETSc vector from an IndexSet. This creates
       * a non ghosted vector.
       */
      explicit Vector (const IndexSet &local,
                       const MPI_Comm &communicator);

      /**
       * Release all memory and return to a state just like after having
       * called the default constructor.
       */
      void clear ();

      /**
       * Copy the given vector. Resize the present vector if necessary. Also
       * take over the MPI communicator of @p v.
       */
      Vector &operator= (const Vector &v);

      /**
       * Copy the given sequential (non-distributed) vector into the present
       * parallel vector. It is assumed that they have the same size, and this
       * operation does not change the partitioning of the parallel vector by
       * which its elements are distributed across several MPI processes. What
       * this operation therefore does is to copy that chunk of the given
       * vector @p v that corresponds to elements of the target vector that
       * are stored locally, and copies them. Elements that are not stored
       * locally are not touched.
       *
       * This being a parallel vector, you must make sure that @em all
       * processes call this function at the same time. It is not possible to
       * change the local part of a parallel vector on only one process,
       * independent of what other processes do, with this function.
       */
      Vector &operator= (const PETScWrappers::Vector &v);

      /**
       * Set all components of the vector to the given number @p s. Simply
       * pass this down to the base class, but we still need to declare this
       * function to make the example given in the discussion about making the
       * constructor explicit work.
       */
      Vector &operator= (const PetscScalar s);

      /**
       * Copy the values of a deal.II vector (as opposed to those of the PETSc
       * vector wrapper class) into this object.
       *
       * Contrary to the case of sequential vectors, this operators requires
       * that the present vector already has the correct size, since we need
       * to have a partition and a communicator present which we otherwise
       * can't get from the source vector.
       */
      template <typename number>
      Vector &operator= (const dealii::Vector<number> &v);

      /**
       * Change the dimension of the vector to @p N. It is unspecified how
       * resizing the vector affects the memory allocation of this object;
       * i.e., it is not guaranteed that resizing it to a smaller size
       * actually also reduces memory consumption, or if for efficiency the
       * same amount of memory is used
       *
       * @p local_size denotes how many of the @p N values shall be stored
       * locally on the present process. for less data.
       *
       * @p communicator denotes the MPI communicator henceforth to be used
       * for this vector.
       *
       * If @p omit_zeroing_entries is false, the vector is filled by zeros.
       * Otherwise, the elements are left an unspecified state.
       */
      void reinit (const MPI_Comm  &communicator,
                   const size_type  N,
                   const size_type  local_size,
                   const bool       omit_zeroing_entries = false);

      /**
       * Change the dimension to that of the vector @p v, and also take over
       * the partitioning into local sizes as well as the MPI communicator.
       * The same applies as for the other @p reinit function.
       *
       * The elements of @p v are not copied, i.e. this function is the same
       * as calling <tt>reinit(v.size(), v.local_size(),
       * omit_zeroing_entries)</tt>.
       */
      void reinit (const Vector &v,
                   const bool    omit_zeroing_entries = false);

      /**
       * Reinit as a vector without ghost elements. See the constructor with
       * same signature for more details.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      void reinit (const IndexSet &local,
                   const IndexSet &ghost,
                   const MPI_Comm &communicator);

      /**
       * Reinit as a vector without ghost elements. See constructor with same
       * signature for more details.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      void reinit (const IndexSet &local,
                   const MPI_Comm &communicator);

      /**
       * Return a reference to the MPI communicator object in use with this
       * vector.
       */
      const MPI_Comm &get_mpi_communicator () const;

      /**
       * Print to a stream. @p precision denotes the desired precision with
       * which values shall be printed, @p scientific whether scientific
       * notation shall be used. If @p across is @p true then the vector is
       * printed in a line, while if @p false then the elements are printed on
       * a separate line each.
       *
       * @note This function overloads the one in the base class to ensure
       * that the right thing happens for parallel vectors that are
       * distributed across processors.
       */
      void print (std::ostream       &out,
                  const unsigned int  precision  = 3,
                  const bool          scientific = true,
                  const bool          across     = true) const;

      /**
       * @copydoc PETScWrappers::VectorBase::all_zero()
       *
       * @note This function overloads the one in the base class to make this
       * a collective operation.
       */
      bool all_zero () const;

    protected:
      /**
       * Create a vector of length @p n. For this class, we create a parallel
       * vector. @p n denotes the total size of the vector to be created. @p
       * local_size denotes how many of these elements shall be stored
       * locally.
       */
      virtual void create_vector (const size_type n,
                                  const size_type local_size);



      /**
       * Create a vector of global length @p n, local size @p local_size and
       * with the specified ghost indices. Note that you need to call
       * update_ghost_values() before accessing those.
       */
      virtual void create_vector (const size_type n,
                                  const size_type local_size,
                                  const IndexSet &ghostnodes);


    private:
      /**
       * Copy of the communicator object to be used for this parallel vector.
       */
      MPI_Comm communicator;
    };


// ------------------ template and inline functions -------------


    /**
     * Global function @p swap which overloads the default implementation of
     * the C++ standard library which uses a temporary object. The function
     * simply exchanges the data of the two vectors.
     *
     * @relates PETScWrappers::MPI::Vector
     * @author Wolfgang Bangerth, 2004
     */
    inline
    void swap (Vector &u, Vector &v)
    {
      u.swap (v);
    }


#ifndef DOXYGEN

    template <typename number>
    Vector::Vector (const MPI_Comm         &communicator,
                    const dealii::Vector<number> &v,
                    const size_type         local_size)
      :
      communicator (communicator)
    {
      Vector::create_vector (v.size(), local_size);

      *this = v;
    }



    inline
    Vector &
    Vector::operator= (const PetscScalar s)
    {
      VectorBase::operator= (s);

      return *this;
    }



    inline
    Vector &
    Vector::operator= (const Vector &v)
    {
      // make sure left- and right-hand side of the assignment are compress()'ed:
      Assert(v.last_action == VectorOperation::unknown,
             internal::VectorReference::ExcWrongMode (VectorOperation::unknown,
                                                      v.last_action));
      Assert(last_action == VectorOperation::unknown,
             internal::VectorReference::ExcWrongMode (VectorOperation::unknown,
                                                      last_action));


      if (v.size()==0)
        {
          // this happens if v has not been initialized to something useful:
          // Vector x,v;x=v;
          // we skip the code below and create a simple serial vector of
          // length 0

          int ierr;
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
          ierr = VecDestroy (vector);
#else
          ierr = VecDestroy (&vector);
#endif
          AssertThrow (ierr == 0, ExcPETScError(ierr));

          const int n = 0;
          ierr = VecCreateSeq (PETSC_COMM_SELF, n, &vector);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
          ghosted = false;
          ghost_indices.clear();
          return *this;
        }

      // if the vectors have different sizes,
      // then first resize the present one
      if (size() != v.size())
        {
          if (v.has_ghost_elements())
            reinit( v.locally_owned_elements(), v.ghost_indices, v.communicator);
          else
            reinit (v.communicator, v.size(), v.local_size(), true);
        }

      const int ierr = VecCopy (v.vector, vector);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      if (has_ghost_elements())
        {
          int ierr;

          ierr = VecGhostUpdateBegin(vector, INSERT_VALUES, SCATTER_FORWARD);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
          ierr = VecGhostUpdateEnd(vector, INSERT_VALUES, SCATTER_FORWARD);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
        }
      return *this;
    }



    template <typename number>
    inline
    Vector &
    Vector::operator= (const dealii::Vector<number> &v)
    {
      Assert (size() == v.size(),
              ExcDimensionMismatch (size(), v.size()));

      // FIXME: the following isn't necessarily fast, but this is due to
      // the fact that PETSc doesn't offer an inlined access operator.
      //
      // if someone wants to contribute some code: to make this code
      // faster, one could either first convert all values to PetscScalar,
      // and then set them all at once using VecSetValues. This has the
      // drawback that it could take quite some memory, if the vector is
      // large, and it would in addition allocate memory on the heap, which
      // is expensive. an alternative would be to split the vector into
      // chunks of, say, 128 elements, convert a chunk at a time and set it
      // in the output vector using VecSetValues. since 128 elements is
      // small enough, this could easily be allocated on the stack (as a
      // local variable) which would make the whole thing much more
      // efficient.
      //
      // a second way to make things faster is for the special case that
      // number==PetscScalar. we could then declare a specialization of
      // this template, and omit the conversion. the problem with this is
      // that the best we can do is to use VecSetValues, but this isn't
      // very efficient either: it wants to see an array of indices, which
      // in this case a) again takes up a whole lot of memory on the heap,
      // and b) is totally dumb since its content would simply be the
      // sequence 0,1,2,3,...,n. the best of all worlds would probably be a
      // function in Petsc that would take a pointer to an array of
      // PetscScalar values and simply copy n elements verbatim into the
      // vector...
      for (size_type i=0; i<v.size(); ++i)
        (*this)(i) = v(i);

      compress (::dealii::VectorOperation::insert);

      return *this;
    }



    inline
    const MPI_Comm &
    Vector::get_mpi_communicator () const
    {
      return communicator;
    }

#endif // DOXYGEN
  }
}

namespace internal
{
  namespace LinearOperator
  {
    template <typename> class ReinitHelper;

    /**
     * A helper class used internally in linear_operator.h. Specialization for
     * PETScWrappers::MPI::Vector.
     */
    template<>
    class ReinitHelper<PETScWrappers::MPI::Vector>
    {
    public:
      template <typename Matrix>
      static
      void reinit_range_vector (const Matrix &matrix,
                                PETScWrappers::MPI::Vector &v,
                                bool omit_zeroing_entries)
      {
        v.reinit(matrix.locally_owned_range_indices(), matrix.get_mpi_communicator());
      }

      template <typename Matrix>
      static
      void reinit_domain_vector(const Matrix &matrix,
                                PETScWrappers::MPI::Vector &v,
                                bool omit_zeroing_entries)
      {
        v.reinit(matrix.locally_owned_domain_indices(), matrix.get_mpi_communicator());
      }
    };

  } /* namespace LinearOperator */
} /* namespace internal */

/**@}*/

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC

/*----------------------------   petsc_parallel_vector.h     ---------------------------*/

#endif
/*----------------------------   petsc_parallel_vector.h     ---------------------------*/

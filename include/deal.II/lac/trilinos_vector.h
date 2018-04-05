// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2018 by the deal.II authors
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

#ifndef dealii_trilinos_vector_h
#define dealii_trilinos_vector_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS
#  include <deal.II/base/subscriptor.h>
#  include <deal.II/base/index_set.h>
#  include <deal.II/base/utilities.h>
#  include <deal.II/base/mpi.h>
#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/vector.h>
#  include <deal.II/lac/vector_operation.h>
#  include <deal.II/lac/vector_type_traits.h>

#  include <vector>
#  include <utility>
#  include <memory>

#  include <Epetra_ConfigDefs.h>
#  ifdef DEAL_II_WITH_MPI // only if MPI is installed
#    include <mpi.h>
#    include <Epetra_MpiComm.h>
#  else
#    include <Epetra_SerialComm.h>
#  endif
#  include <Epetra_FEVector.h>
#  include <Epetra_Map.h>
#  include <Epetra_LocalMap.h>

DEAL_II_NAMESPACE_OPEN


/**
 * @addtogroup TrilinosWrappers
 * @{
 */

/**
 * A namespace in which wrapper classes for Trilinos objects reside.
 *
 * @ingroup TrilinosWrappers
 */
namespace TrilinosWrappers
{
  class SparseMatrix;

  /**
   * @cond internal
   */

  /**
   * A namespace for internal implementation details of the TrilinosWrapper
   * members.
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
      VectorReference (MPI::Vector     &vector,
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
      operator= (const VectorReference &r) const;

      /**
       * Same as above but for non-const reference objects.
       */
      VectorReference &
      operator= (const VectorReference &r);

      /**
       * Set the referenced element of the vector to <tt>s</tt>.
       */
      const VectorReference &
      operator= (const TrilinosScalar &s) const;

      /**
       * Add <tt>s</tt> to the referenced element of the vector->
       */
      const VectorReference &
      operator+= (const TrilinosScalar &s) const;

      /**
       * Subtract <tt>s</tt> from the referenced element of the vector->
       */
      const VectorReference &
      operator-= (const TrilinosScalar &s) const;

      /**
       * Multiply the referenced element of the vector by <tt>s</tt>.
       */
      const VectorReference &
      operator*= (const TrilinosScalar &s) const;

      /**
       * Divide the referenced element of the vector by <tt>s</tt>.
       */
      const VectorReference &
      operator/= (const TrilinosScalar &s) const;

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
      MPI::Vector   &vector;

      /**
       * Index of the referenced element of the vector.
       */
      const size_type  index;

      /**
       * Make the vector class a friend, so that it can create objects of the
       * present type.
       */
      friend class ::dealii::TrilinosWrappers::MPI::Vector;
    };
  }
  /**
   * @endcond
   */

  namespace
  {
#ifndef DEAL_II_WITH_64BIT_INDICES
    // define a helper function that queries the global ID of local ID of
    // an Epetra_BlockMap object  by calling either the 32- or 64-bit
    // function necessary.
    inline
    int gid(const Epetra_BlockMap &map, int i)
    {
      return map.GID(i);
    }
#else
    // define a helper function that queries the global ID of local ID of
    // an Epetra_BlockMap object  by calling either the 32- or 64-bit
    // function necessary.
    inline
    long long int gid(const Epetra_BlockMap &map, int i)
    {
      return map.GID64(i);
    }
#endif
  }

  /**
   * Namespace for Trilinos vector classes that work in parallel over MPI.
   *
   * @ingroup TrilinosWrappers
   * @author Martin Kronbichler, Wolfgang Bangerth, Daniel Arndt, 2008, 2017
   */
  namespace MPI
  {
    class BlockVector;

    /**
     * This class implements a wrapper to use the Trilinos distributed vector
     * class Epetra_FEVector, the (parallel) partitioning of which
     * is governed by an Epetra_Map.
     * The Epetra_FEVector is precisely the kind of vector
     * we deal with all the time - we probably get it from some assembly
     * process, where also entries not locally owned might need to written and
     * hence need to be forwarded to the owner.
     *
     * The interface of this class is modeled after the existing Vector class in
     * deal.II. It has almost the same member functions, and is often
     * exchangeable. However, since Trilinos only supports a single scalar type
     * (double), it is not templated, and only works with that type.
     *
     * Note that Trilinos only guarantees that operations do what you expect
     * if the function @p GlobalAssemble has been called after vector assembly
     * in order to distribute the data. This is necessary since some processes
     * might have accumulated data of elements that are not owned by
     * themselves, but must be sent to the owning process. In order to avoid
     * using the wrong data, you need to call Vector::compress() before you
     * actually use the vectors.
     *
     * <h3>Parallel communication model</h3>
     *
     * The parallel functionality of Trilinos is built on top of the Message
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
     * time the next a call to a Trilinos function generates an MPI message on
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
     * Trilinos does of course allow read access to individual elements of a
     * vector, but in the distributed case only to elements that are stored
     * locally. We implement this through calls like <tt>d=vec(i)</tt>.
     * However, if you access an element outside the locally stored range, an
     * exception is generated.
     *
     * In contrast to read access, Trilinos (and the respective deal.II
     * wrapper classes) allow to write (or add) to individual elements of
     * vectors, even if they are stored on a different process. You can do
     * this by writing into or adding to elements using the syntax
     * <tt>vec(i)=d</tt> or <tt>vec(i)+=d</tt>, or similar operations. There
     * is one catch, however, that may lead to very confusing error messages:
     * Trilinos requires application programs to call the compress() function
     * when they switch from performing a set of operations that add to
     * elements, to performing a set of operations that write to elements. The
     * reasoning is that all processes might accumulate addition operations to
     * elements, even if multiple processes write to the same elements. By the
     * time we call compress() the next time, all these additions are
     * executed. However, if one process adds to an element, and another
     * overwrites to it, the order of execution would yield non-deterministic
     * behavior if we don't make sure that a synchronization with compress()
     * happens in between.
     *
     * In order to make sure these calls to compress() happen at the
     * appropriate time, the deal.II wrappers keep a state variable that store
     * which is the presently allowed operation: additions or writes. If it
     * encounters an operation of the opposite kind, it calls compress() and
     * flips the state. This can sometimes lead to very confusing behavior, in
     * code that may for example look like this:
     *
     * @code
     * TrilinosWrappers::MPI::Vector vector;
     * // do some write operations on the vector
     * for (size_type i=0; i<vector->size(); ++i)
     *   vector(i) = i;
     *
     *                   // do some additions to vector elements, but
     *                   // only for some elements
     *   for (size_type i=0; i<vector->size(); ++i)
     *     if (some_condition(i) == true)
     *       vector(i) += 1;
     *
     *                   // do another collective operation
     *   const double norm = vector->l2_norm();
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
     *
     * <h3>Ghost elements of vectors</h3>
     *
     * Parallel vectors come in two kinds: without and with ghost elements.
     * Vectors without ghost elements uniquely partition the vector elements
     * between processors: each vector entry has exactly one processor that
     * owns it. For such vectors, you can read those elements that the
     * processor you are currently on owns, and you can write into any element
     * whether you own it or not: if you don't own it, the value written or
     * added to a vector element will be shipped to the processor that owns
     * this vector element the next time you call compress(), as described
     * above.
     *
     * What we call a 'ghosted' vector (see
     * @ref GlossGhostedVector "vectors with ghost elements"
     * ) is simply a view of the parallel vector where the element
     * distributions overlap. The 'ghosted' Trilinos vector in itself has no
     * idea of which entries are ghosted and which are locally owned. In fact,
     * a ghosted vector may not even store all of the elements a non- ghosted
     * vector would store on the current processor.  Consequently, for
     * Trilinos vectors, there is no notion of an 'owner' of vector elements
     * in the way we have it in the non-ghost case view.
     *
     * This explains why we do not allow writing into ghosted vectors on the
     * Trilinos side: Who would be responsible for taking care of the
     * duplicated entries, given that there is not such information as locally
     * owned indices? In other words, since a processor doesn't know which
     * other processors own an element, who would it send a value to if one
     * were to write to it? The only possibility would be to send this
     * information to <i>all</i> other processors, but that is clearly not
     * practical. Thus, we only allow reading from ghosted vectors, which
     * however we do very often.
     *
     * So how do you fill a ghosted vector if you can't write to it? This only
     * happens through the assignment with a non-ghosted vector. It can go
     * both ways (non-ghosted is assigned to a ghosted vector, and a ghosted
     * vector is assigned to a non-ghosted one; the latter one typically only
     * requires taking out the locally owned part as most often ghosted
     * vectors store a superset of elements of non-ghosted ones). In general,
     * you send data around with that operation and it all depends on the
     * different views of the two vectors. Trilinos also allows you to get
     * subvectors out of a big vector that way.
     *
     *
     * <h3>Thread safety of Trilinos vectors</h3>
     *
     * When writing into Trilinos vectors from several threads in shared
     * memory, several things must be kept in mind as there is no built-in
     * locks in this class to prevent data races. Simultaneous access to the
     * same vector entry at the same time results in data races and must be
     * explicitly avoided by the user. However, it is possible to access
     * <b>different</b> entries of the vector from several threads
     * simultaneously when only one MPI process is present or the vector has
     * been constructed with an additional index set for ghost entries in
     * write mode.
     *
     * @ingroup TrilinosWrappers
     * @ingroup Vectors
     * @author Martin Kronbichler, Wolfgang Bangerth, Daniel Arndt,
     *         2008, 2009, 2017
     */
    class Vector : public Subscriptor
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
      Vector ();

      /**
       * Copy constructor using the given vector.
       */
      Vector (const Vector &v);

      /**
       * This constructor takes an IndexSet that defines how to distribute the
       * individual components among the MPI processors. Since it also
       * includes information about the size of the vector, this is all we
       * need to generate a %parallel vector.
       *
       * Depending on whether the @p parallel_partitioning argument uniquely
       * subdivides elements among processors or not, the resulting vector may
       * or may not have ghost elements. See the general documentation of this
       * class for more information.
       *
       * In case the provided IndexSet forms an overlapping partitioning,
       * it is not clear which elements are owned by which process and
       * locally_owned_elements() will return an IndexSet of size zero.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      explicit Vector (const IndexSet &parallel_partitioning,
                       const MPI_Comm &communicator = MPI_COMM_WORLD);

      /**
       * Creates a ghosted parallel vector.
       *
       * Depending on whether the @p ghost argument uniquely subdivides
       * elements among processors or not, the resulting vector may or may not
       * have ghost elements. See the general documentation of this class for
       * more information.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      Vector (const IndexSet &local,
              const IndexSet &ghost,
              const MPI_Comm &communicator = MPI_COMM_WORLD);

      /**
       * Copy constructor from the TrilinosWrappers vector class. Since a
       * vector of this class does not necessarily need to be distributed
       * among processes, the user needs to supply us with an IndexSet and an
       * MPI communicator that set the partitioning details.
       *
       * Depending on whether the @p parallel_partitioning argument uniquely
       * subdivides elements among processors or not, the resulting vector may
       * or may not have ghost elements. See the general documentation of this
       * class for more information.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      Vector (const IndexSet   &parallel_partitioning,
              const Vector &v,
              const MPI_Comm   &communicator = MPI_COMM_WORLD);

      /**
       * Copy-constructor from deal.II vectors. Sets the dimension to that of
       * the given vector, and copies all the elements.
       *
       * Depending on whether the @p parallel_partitioning argument uniquely
       * subdivides elements among processors or not, the resulting vector may
       * or may not have ghost elements. See the general documentation of this
       * class for more information.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      template <typename Number>
      Vector (const IndexSet               &parallel_partitioning,
              const dealii::Vector<Number> &v,
              const MPI_Comm               &communicator = MPI_COMM_WORLD);

      /**
       * Move constructor. Creates a new vector by stealing the internal data
       * of the vector @p v.
       */
      Vector (Vector &&v) noexcept;

      /**
       * Destructor.
       */
      ~Vector () = default;

      /**
       * Release all memory and return to a state just like after having called
       * the default constructor.
       */
      void clear ();

      /**
       * Reinit functionality. This function sets the calling vector to the
       * dimension and the parallel distribution of the input vector, but does
       * not copy the elements in <tt>v</tt>. If <tt>omit_zeroing_entries</tt>
       * is not <tt>true</tt>, the elements in the vector are initialized with
       * zero. If it is set to <tt>true</tt>, the vector entries are in an
       * unspecified state and the user has to set all elements. In the
       * current implementation, this method does not touch the vector entries
       * in case the vector layout is unchanged from before, otherwise entries
       * are set to zero.  Note that this behavior might change between
       * releases without notification.
       *
       * This function has a third argument, <tt>allow_different_maps</tt>,
       * that allows for an exchange of data between two equal-sized vectors
       * (but being distributed differently among the processors). A trivial
       * application of this function is to generate a replication of a whole
       * vector on each machine, when the calling vector is built with a map
       * consisting of all indices on each process, and <tt>v</tt>
       * is a distributed vector. In this case, the variable
       * <tt>omit_zeroing_entries</tt> needs to be set to <tt>false</tt>,
       * since it does not make sense to exchange data between differently
       * parallelized vectors without touching the elements.
       */
      void reinit (const Vector &v,
                   const bool        omit_zeroing_entries = false,
                   const bool        allow_different_maps = false);

      /**
       * Reinit functionality. This function destroys the old vector content
       * and generates a new one based on the input partitioning.  The flag
       * <tt>omit_zeroing_entries</tt> determines whether the vector should be
       * filled with zero (false). If the flag is set to <tt>true</tt>, the
       * vector entries are in an unspecified state and the user has to set
       * all elements. In the current implementation, this method still sets
       * the entries to zero, but this might change between releases without
       * notification.
       *
       * Depending on whether the @p parallel_partitioning argument uniquely
       * subdivides elements among processors or not, the resulting vector may
       * or may not have ghost elements. See the general documentation of this
       * class for more information.
       *
       * In case @p parallel_partitioning is overlapping, it is not clear which
       * process should own which elements. Hence, locally_owned_elements()
       * returns an empty IndexSet in this case.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      void reinit (const IndexSet &parallel_partitioning,
                   const MPI_Comm &communicator = MPI_COMM_WORLD,
                   const bool      omit_zeroing_entries = false);

      /**
       * Reinit functionality. This function destroys the old vector content
       * and generates a new one based on the input partitioning. In addition
       * to just specifying one index set as in all the other methods above,
       * this method allows to supply an additional set of ghost entries.
       * There are two different versions of a vector that can be created. If
       * the flag @p vector_writable is set to @p false, the vector only
       * allows read access to the joint set of @p parallel_partitioning and
       * @p ghost_entries. The effect of the reinit method is then equivalent
       * to calling the other reinit method with an index set containing both
       * the locally owned entries and the ghost entries.
       *
       * If the flag @p vector_writable is set to true, this creates an
       * alternative storage scheme for ghost elements that allows multiple
       * threads to write into the vector (for the other reinit methods, only
       * one thread is allowed to write into the ghost entries at a time).
       *
       * Depending on whether the @p ghost_entries argument uniquely
       * subdivides elements among processors or not, the resulting vector may
       * or may not have ghost elements. See the general documentation of this
       * class for more information.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      void reinit (const IndexSet &locally_owned_entries,
                   const IndexSet &ghost_entries,
                   const MPI_Comm &communicator = MPI_COMM_WORLD,
                   const bool      vector_writable = false);

      /**
       * Create vector by merging components from a block vector.
       */
      void reinit (const BlockVector &v,
                   const bool         import_data = false);

      /**
       * Compress the underlying representation of the Trilinos object, i.e.
       * flush the buffers of the vector object if it has any. This function is
       * necessary after writing into a vector element-by-element and before
       * anything else can be done on it.
       *
       * The (defaulted) argument can be used to specify the compress mode
       * (<code>Add</code> or <code>Insert</code>) in case the vector has not
       * been written to since the last time this function was called. The
       * argument is ignored if the vector has been added or written to since
       * the last time compress() was called.
       *
       * See
       * @ref GlossCompress "Compressing distributed objects"
       * for more information.
       */
      void compress (::dealii::VectorOperation::values operation);

      /**
       * Set all components of the vector to the given number @p s. Simply
       * pass this down to the base class, but we still need to declare this
       * function to make the example given in the discussion about making the
       * constructor explicit work.
       * the constructor explicit work.
       *
       * Since the semantics of assigning a scalar to a vector are not
       * immediately clear, this operator can only be used if you want
       * to set the entire vector to zero. This allows the intuitive notation
       * <tt>v=0</tt>.
       */
      Vector &operator= (const TrilinosScalar s);

      /**
       * Copy the given vector. Resize the present vector if necessary. In
       * this case, also the Epetra_Map that designs the parallel partitioning
       * is taken from the input vector.
       */
      Vector &operator= (const Vector &v);

      /**
       * Move the given vector. This operator replaces the present vector with
       * @p v by efficiently swapping the internal data structures.
       */
      Vector &operator= (Vector &&v) noexcept;

      /**
       * Another copy function. This one takes a deal.II vector and copies it
       * into a TrilinosWrapper vector. Note that since we do not provide any
       * Epetra_map that tells about the partitioning of the vector among the
       * MPI processes, the size of the TrilinosWrapper vector has to be the
       * same as the size of the input vector.
       */
      template <typename Number>
      Vector &operator= (const ::dealii::Vector<Number> &v);

      /**
       * This reinit function is meant to be used for parallel calculations
       * where some non-local data has to be used. The typical situation where
       * one needs this function is the call of the
       * FEValues<dim>::get_function_values function (or of some derivatives)
       * in parallel. Since it is usually faster to retrieve the data in
       * advance, this function can be called before the assembly forks out to
       * the different processors. What this function does is the following:
       * It takes the information in the columns of the given matrix and looks
       * which data couples between the different processors. That data is
       * then queried from the input vector. Note that you should not write to
       * the resulting vector any more, since the some data can be stored
       * several times on different processors, leading to unpredictable
       * results. In particular, such a vector cannot be used for matrix-
       * vector products as for example done during the solution of linear
       * systems.
       */
      void import_nonlocal_data_for_fe
      (const dealii::TrilinosWrappers::SparseMatrix &matrix,
       const Vector                                 &vector);

      /**
       * Test for equality. This function assumes that the present vector and
       * the one to compare with have the same size already, since comparing
       * vectors of different sizes makes not much sense anyway.
       */
      bool operator== (const Vector &v) const;

      /**
       * Test for inequality. This function assumes that the present vector and
       * the one to compare with have the same size already, since comparing
       * vectors of different sizes makes not much sense anyway.
       */
      bool operator!= (const Vector &v) const;

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
       * @note The description above is true most of the time, but not always.
       * In particular, Trilinos vectors need not store contiguous ranges of
       * elements such as $[i,i+n)$. Rather, it can store vectors where the
       * elements are distributed in an arbitrary way across all processors and
       * each processor simply stores a particular subset, not necessarily
       * contiguous. In this case, this function clearly makes no sense since it
       * could, at best, return a range that includes all elements that are
       * stored locally. Thus, the function only succeeds if the locally stored
       * range is indeed contiguous. It will trigger an assertion if the local
       * portion of the vector is not contiguous.
       */
      std::pair<size_type, size_type> local_range () const;

      /**
       * Return whether @p index is in the local range or not, see also
       * local_range().
       *
       * @note The same limitation for the applicability of this function
       * applies as listed in the documentation of local_range().
       */
      bool in_local_range (const size_type index) const;

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
      IndexSet locally_owned_elements () const;

      /**
       * Return if the vector contains ghost elements. This answer is true if
       * there are ghost elements on at least one process.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      bool has_ghost_elements() const;

      /**
       * This function only exists for compatibility with the @p
       * LinearAlgebra::distributed::Vector class and does nothing: this class
       * implements ghost value updates in a different way that is a better fit
       * with the underlying Trilinos vector object.
       */
      void update_ghost_values () const;

      /**
       * Return the scalar (inner) product of two vectors. The vectors must have
       * the same size.
       */
      TrilinosScalar operator* (const Vector &vec) const;

      /**
       * Return the square of the $l_2$-norm.
       */
      real_type norm_sqr () const;

      /**
       * Mean value of the elements of this vector.
       */
      TrilinosScalar mean_value () const;

      /**
       * Compute the minimal value of the elements of this vector.
       */
      TrilinosScalar min () const;

      /**
       * Compute the maximal value of the elements of this vector.
       */
      TrilinosScalar max () const;

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
       * transfer. However, for Trilinos vectors such a combined operation is
       * not natively supported and thus the cost is completely equivalent as
       * calling the two methods separately.
       *
       * For complex-valued vectors, the scalar product in the second step is implemented as
       * $\left<v,w\right>=\sum_i v_i \bar{w_i}$.
       */
      TrilinosScalar add_and_dot (const TrilinosScalar a,
                                  const Vector    &V,
                                  const Vector    &W);

      /**
       * Return whether the vector contains only elements with value zero. This
       * is a collective operation. This function is expensive, because
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
       * sense for elements that are actually present on the calling processor.
       * Otherwise, an exception is thrown.
       */
      reference
      operator() (const size_type index);

      /**
       * Provide read-only access to an element.
       *
       * When using a vector distributed with MPI, this operation only makes
       * sense for elements that are actually present on the calling processor.
       * Otherwise, an exception is thrown.
       */
      TrilinosScalar
      operator() (const size_type index) const;

      /**
       * Provide access to a given element, both read and write.
       *
       * Exactly the same as operator().
       */
      reference
      operator[] (const size_type index);

      /**
       * Provide read-only access to an element.
       *
       * Exactly the same as operator().
       */
      TrilinosScalar
      operator[] (const size_type index) const;

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
      void extract_subvector_to (const std::vector<size_type> &indices,
                                 std::vector<TrilinosScalar> &values) const;

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
      void extract_subvector_to (ForwardIterator          indices_begin,
                                 const ForwardIterator    indices_end,
                                 OutputIterator           values_begin) const;

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
       * Return constant iterator to the start of the locally owned elements of
       * the vector.
       */
      const_iterator begin () const;

      /**
       * Return an iterator pointing to the element past the end of the array of
       * locally owned entries.
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
       * vector, this function allows to set a whole set of elements at once.
       * The indices of the elements to be set are stated in the first argument,
       * the corresponding values in the second.
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
      Vector &operator*= (const TrilinosScalar factor);

      /**
       * Divide the entire vector by a fixed factor.
       */
      Vector &operator/= (const TrilinosScalar factor);

      /**
       * Add the given vector to the present one.
       */
      Vector &operator+= (const Vector &V);

      /**
       * Subtract the given vector from the present one.
       */
      Vector &operator-= (const Vector &V);

      /**
       * Addition of @p s to all components. Note that @p s is a scalar and not
       * a vector.
       */
      void add (const TrilinosScalar s);

      /**
       * Simple vector addition, equal to the <tt>operator+=</tt>.
       *
       * Though, if the second argument <tt>allow_different_maps</tt> is set,
       * then it is possible to add data from a vector that uses a different
       * map, i.e., a vector whose elements are split across processors
       * differently. This may include vectors with ghost elements, for example.
       * In general, however, adding vectors with a different element-to-
       * processor map requires communicating data among processors and,
       * consequently, is a slower operation than when using vectors using the
       * same map.
       */
      void add (const Vector &V,
                const bool        allow_different_maps = false);

      /**
       * Simple addition of a multiple of a vector, i.e. <tt>*this += a*V</tt>.
       */
      void add (const TrilinosScalar  a,
                const Vector     &V);

      /**
       * Multiple addition of scaled vectors, i.e. <tt>*this += a*V + b*W</tt>.
       */
      void add (const TrilinosScalar  a,
                const Vector     &V,
                const TrilinosScalar  b,
                const Vector     &W);

      /**
       * Scaling and simple vector addition, i.e.  <tt>*this = s*(*this) +
       * V</tt>.
       */
      void sadd (const TrilinosScalar  s,
                 const Vector     &V);

      /**
       * Scaling and simple addition, i.e.  <tt>*this = s*(*this) + a*V</tt>.
       */
      void sadd (const TrilinosScalar  s,
                 const TrilinosScalar  a,
                 const Vector     &V);

      /**
       * Scale each element of this vector by the corresponding element in the
       * argument. This function is mostly meant to simulate multiplication (and
       * immediate re-assignment) by a diagonal scaling matrix.
       */
      void scale (const Vector &scaling_factors);

      /**
       * Assignment <tt>*this = a*V</tt>.
       */
      void equ (const TrilinosScalar  a,
                const Vector     &V);
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
       * This function is analogous to the @p swap function of all C++
       * standard containers. Also, there is a global function
       * <tt>swap(u,v)</tt> that simply calls <tt>u.swap(v)</tt>, again in
       * analogy to standard functions.
       */
      void swap (Vector &v);

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
       * overwriting them (to make synchronization of parallel computations
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
       * Pointer to the actual Epetra vector object. This may represent a vector
       * that is in fact distributed among multiple processors. The object
       * requires an existing Epetra_Map for storing data when setting it up.
       */
      std::unique_ptr<Epetra_FEVector> vector;

      /**
       * A vector object in Trilinos to be used for collecting the non-local
       * elements if the vector was constructed with an additional IndexSet
       * describing ghost elements.
       */
      std::unique_ptr<Epetra_MultiVector> nonlocal_vector;

      /**
       * An IndexSet storing the indices this vector owns exclusively.
       */
      IndexSet owned_elements;

      /**
       * Make the reference class a friend.
       */
      friend class internal::VectorReference;
    };




// ------------------- inline and template functions --------------


    /**
     * Global function @p swap which overloads the default implementation of
     * the C++ standard library which uses a temporary object. The function
     * simply exchanges the data of the two vectors.
     *
     * @relatesalso TrilinosWrappers::MPI::Vector
     * @author Martin Kronbichler, Wolfgang Bangerth, 2008
     */
    inline
    void swap (Vector &u, Vector &v)
    {
      u.swap (v);
    }
  }

#ifndef DOXYGEN

  namespace internal
  {
    inline
    VectorReference::VectorReference (MPI::Vector      &vector,
                                      const size_type  index)
      :
      vector (vector),
      index (index)
    {}


    inline
    const VectorReference &
    VectorReference::operator= (const VectorReference &r) const
    {
      // as explained in the class
      // documentation, this is not the copy
      // operator. so simply pass on to the
      // "correct" assignment operator
      *this = static_cast<TrilinosScalar> (r);

      return *this;
    }



    inline
    VectorReference &
    VectorReference::operator= (const VectorReference &r)
    {
      // as above
      *this = static_cast<TrilinosScalar> (r);

      return *this;
    }


    inline
    const VectorReference &
    VectorReference::operator= (const TrilinosScalar &value) const
    {
      vector.set (1, &index, &value);
      return *this;
    }



    inline
    const VectorReference &
    VectorReference::operator+= (const TrilinosScalar &value) const
    {
      vector.add (1, &index, &value);
      return *this;
    }



    inline
    const VectorReference &
    VectorReference::operator-= (const TrilinosScalar &value) const
    {
      TrilinosScalar new_value = -value;
      vector.add (1, &index, &new_value);
      return *this;
    }



    inline
    const VectorReference &
    VectorReference::operator*= (const TrilinosScalar &value) const
    {
      TrilinosScalar new_value = static_cast<TrilinosScalar>(*this) * value;
      vector.set (1, &index, &new_value);
      return *this;
    }



    inline
    const VectorReference &
    VectorReference::operator/= (const TrilinosScalar &value) const
    {
      TrilinosScalar new_value = static_cast<TrilinosScalar>(*this) / value;
      vector.set (1, &index, &new_value);
      return *this;
    }
  }

  namespace MPI
  {
    inline
    bool
    Vector::in_local_range (const size_type index) const
    {
      std::pair<size_type, size_type> range = local_range();

      return ((index >= range.first) && (index <  range.second));
    }



    inline
    IndexSet
    Vector::locally_owned_elements() const
    {
      Assert(owned_elements.size()==size(),
             ExcMessage("The locally owned elements have not been properly initialized!"
                        " This happens for example if this object has been initialized"
                        " with exactly one overlapping IndexSet."));
      return owned_elements;
    }



    inline
    bool
    Vector::has_ghost_elements() const
    {
      return has_ghosts;
    }



    inline
    void
    Vector::update_ghost_values () const
    {}



    inline
    internal::VectorReference
    Vector::operator() (const size_type index)
    {
      return internal::VectorReference (*this, index);
    }



    inline
    internal::VectorReference
    Vector::operator[] (const size_type index)
    {
      return operator() (index);
    }



    inline
    TrilinosScalar
    Vector::operator[] (const size_type index) const
    {
      return operator() (index);
    }



    inline
    void Vector::extract_subvector_to (const std::vector<size_type> &indices,
                                       std::vector<TrilinosScalar>  &values) const
    {
      for (size_type i = 0; i < indices.size(); ++i)
        values[i] = operator()(indices[i]);
    }



    template <typename ForwardIterator, typename OutputIterator>
    inline
    void Vector::extract_subvector_to (ForwardIterator          indices_begin,
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
    Vector::iterator
    Vector::begin()
    {
      return (*vector)[0];
    }



    inline
    Vector::iterator
    Vector::end()
    {
      return (*vector)[0]+local_size();
    }



    inline
    Vector::const_iterator
    Vector::begin() const
    {
      return (*vector)[0];
    }



    inline
    Vector::const_iterator
    Vector::end() const
    {
      return (*vector)[0]+local_size();
    }



    inline
    void
    Vector::set (const std::vector<size_type>      &indices,
                 const std::vector<TrilinosScalar>  &values)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert (!has_ghost_elements(), ExcGhostsPresent());

      Assert (indices.size() == values.size(),
              ExcDimensionMismatch(indices.size(),values.size()));

      set (indices.size(), indices.data(), values.data());
    }



    inline
    void
    Vector::set (const std::vector<size_type>           &indices,
                 const ::dealii::Vector<TrilinosScalar> &values)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert (!has_ghost_elements(), ExcGhostsPresent());

      Assert (indices.size() == values.size(),
              ExcDimensionMismatch(indices.size(),values.size()));

      set (indices.size(), indices.data(), values.begin());
    }



    inline
    void
    Vector::set (const size_type       n_elements,
                 const size_type      *indices,
                 const TrilinosScalar *values)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert (!has_ghost_elements(), ExcGhostsPresent());

      if (last_action == Add)
        {
          const int ierr = vector->GlobalAssemble(Add);
          AssertThrow (ierr == 0, ExcTrilinosError(ierr));
        }

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
    Vector::add (const std::vector<size_type>      &indices,
                 const std::vector<TrilinosScalar>  &values)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert (!has_ghost_elements(), ExcGhostsPresent());
      Assert (indices.size() == values.size(),
              ExcDimensionMismatch(indices.size(),values.size()));

      add (indices.size(), indices.data(), values.data());
    }



    inline
    void
    Vector::add (const std::vector<size_type>           &indices,
                 const ::dealii::Vector<TrilinosScalar> &values)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert (!has_ghost_elements(), ExcGhostsPresent());
      Assert (indices.size() == values.size(),
              ExcDimensionMismatch(indices.size(),values.size()));

      add (indices.size(), indices.data(), values.begin());
    }



    inline
    void
    Vector::add (const size_type       n_elements,
                 const size_type      *indices,
                 const TrilinosScalar *values)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert (!has_ghost_elements(), ExcGhostsPresent());

      if (last_action != Add)
        {
          if (last_action == Insert)
            {
              const int ierr = vector->GlobalAssemble(Insert);
              AssertThrow (ierr == 0, ExcTrilinosError(ierr));
            }
          last_action = Add;
        }

      for (size_type i=0; i<n_elements; ++i)
        {
          const size_type row = indices[i];
          const TrilinosWrappers::types::int_type local_row = vector->Map().LID(static_cast<TrilinosWrappers::types::int_type>(row));
          if (local_row != -1)
            (*vector)[0][local_row] += values[i];
          else if (nonlocal_vector.get() == nullptr)
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
    Vector::size_type
    Vector::size () const
    {
#ifndef DEAL_II_WITH_64BIT_INDICES
      return (size_type) (vector->Map().MaxAllGID() + 1 - vector->Map().MinAllGID());
#else
      return (size_type) (vector->Map().MaxAllGID64() + 1 - vector->Map().MinAllGID64());
#endif
    }



    inline
    Vector::size_type
    Vector::local_size () const
    {
      return (size_type) vector->Map().NumMyElements();
    }



    inline
    std::pair<Vector::size_type, Vector::size_type>
    Vector::local_range () const
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
    Vector::operator* (const Vector &vec) const
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
    Vector::real_type
    Vector::norm_sqr () const
    {
      const TrilinosScalar d = l2_norm();
      return d*d;
    }



    inline
    TrilinosScalar
    Vector::mean_value () const
    {
      Assert (!has_ghost_elements(), ExcGhostsPresent());

      TrilinosScalar mean;
      const int ierr = vector->MeanValue (&mean);
      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

      return mean;
    }



    inline
    TrilinosScalar
    Vector::min () const
    {
      TrilinosScalar min_value;
      const int ierr = vector->MinValue (&min_value);
      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

      return min_value;
    }



    inline
    TrilinosScalar
    Vector::max () const
    {
      TrilinosScalar max_value;
      const int ierr = vector->MaxValue (&max_value);
      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

      return max_value;
    }



    inline
    Vector::real_type
    Vector::l1_norm () const
    {
      Assert (!has_ghost_elements(), ExcGhostsPresent());

      TrilinosScalar d;
      const int ierr = vector->Norm1 (&d);
      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

      return d;
    }



    inline
    Vector::real_type
    Vector::l2_norm () const
    {
      Assert (!has_ghost_elements(), ExcGhostsPresent());

      TrilinosScalar d;
      const int ierr = vector->Norm2 (&d);
      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

      return d;
    }



    inline
    Vector::real_type
    Vector::lp_norm (const TrilinosScalar p) const
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
    Vector::real_type
    Vector::linfty_norm () const
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



    inline
    TrilinosScalar
    Vector::add_and_dot (const TrilinosScalar a,
                         const Vector &V,
                         const Vector &W)
    {
      this->add(a, V);
      return *this * W;
    }



    // inline also scalar products, vector
    // additions etc. since they are all
    // representable by a single Trilinos
    // call. This reduces the overhead of the
    // wrapper class.
    inline
    Vector &
    Vector::operator*= (const TrilinosScalar a)
    {
      AssertIsFinite(a);

      const int ierr = vector->Scale(a);
      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

      return *this;
    }



    inline
    Vector &
    Vector::operator/= (const TrilinosScalar a)
    {
      AssertIsFinite(a);

      const TrilinosScalar factor = 1./a;

      AssertIsFinite(factor);

      const int ierr = vector->Scale(factor);
      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

      return *this;
    }



    inline
    Vector &
    Vector::operator+= (const Vector &v)
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
    Vector &
    Vector::operator-= (const Vector &v)
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
    Vector::add (const TrilinosScalar s)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert (!has_ghost_elements(), ExcGhostsPresent());
      AssertIsFinite(s);

      size_type n_local = local_size();
      for (size_type i=0; i<n_local; i++)
        (*vector)[0][i] += s;
    }



    inline
    void
    Vector::add (const TrilinosScalar  a,
                 const Vector     &v)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert (!has_ghost_elements(), ExcGhostsPresent());
      Assert (local_size() == v.local_size(),
              ExcDimensionMismatch(local_size(), v.local_size()));

      AssertIsFinite(a);

      const int ierr = vector->Update(a, *(v.vector), 1.);
      AssertThrow (ierr == 0, ExcTrilinosError(ierr));
    }



    inline
    void
    Vector::add (const TrilinosScalar  a,
                 const Vector     &v,
                 const TrilinosScalar  b,
                 const Vector     &w)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert (!has_ghost_elements(), ExcGhostsPresent());
      Assert (local_size() == v.local_size(),
              ExcDimensionMismatch(local_size(), v.local_size()));
      Assert (local_size() == w.local_size(),
              ExcDimensionMismatch(local_size(), w.local_size()));

      AssertIsFinite(a);
      AssertIsFinite(b);

      const int ierr = vector->Update(a, *(v.vector), b, *(w.vector), 1.);

      AssertThrow (ierr == 0, ExcTrilinosError(ierr));
    }



    inline
    void
    Vector::sadd (const TrilinosScalar  s,
                  const Vector     &v)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert (!has_ghost_elements(), ExcGhostsPresent());
      Assert (size() == v.size(),
              ExcDimensionMismatch (size(), v.size()));

      AssertIsFinite(s);

      // We assume that the vectors have the same Map
      // if the local size is the same and if the vectors are not ghosted
      if (local_size() == v.local_size() && !v.has_ghost_elements())
        {
          Assert (this->vector->Map().SameAs(v.vector->Map())==true,
                  ExcDifferentParallelPartitioning());
          const int ierr = vector->Update(1., *(v.vector), s);
          AssertThrow (ierr == 0, ExcTrilinosError(ierr));
        }
      else
        {
          (*this) *= s;
          this->add(v, true);
        }
    }



    inline
    void
    Vector::sadd (const TrilinosScalar  s,
                  const TrilinosScalar  a,
                  const Vector     &v)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert (!has_ghost_elements(), ExcGhostsPresent());
      Assert (size() == v.size(),
              ExcDimensionMismatch (size(), v.size()));
      AssertIsFinite(s);
      AssertIsFinite(a);

      // We assume that the vectors have the same Map
      // if the local size is the same and if the vectors are not ghosted
      if (local_size() == v.local_size() && !v.has_ghost_elements())
        {
          Assert (this->vector->Map().SameAs(v.vector->Map())==true,
                  ExcDifferentParallelPartitioning());
          const int ierr = vector->Update(a, *(v.vector), s);
          AssertThrow (ierr == 0, ExcTrilinosError(ierr));
        }
      else
        {
          (*this) *= s;
          Vector tmp = v;
          tmp *= a;
          this->add(tmp, true);
        }
    }



    inline
    void
    Vector::scale (const Vector &factors)
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
    Vector::equ (const TrilinosScalar  a,
                 const Vector     &v)
    {
      // if we have ghost values, do not allow
      // writing to this vector at all.
      Assert (!has_ghost_elements(), ExcGhostsPresent());
      AssertIsFinite(a);

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
    const Epetra_MultiVector &
    Vector::trilinos_vector () const
    {
      return static_cast<const Epetra_MultiVector &>(*vector);
    }



    inline
    Epetra_FEVector &
    Vector::trilinos_vector ()
    {
      return *vector;
    }



    inline
    const Epetra_Map &
    Vector::vector_partitioner () const
    {
      return static_cast<const Epetra_Map &>(vector->Map());
    }



    inline
    const MPI_Comm &
    Vector::get_mpi_communicator () const
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

    template <typename number>
    Vector::Vector (const IndexSet               &parallel_partitioner,
                    const dealii::Vector<number> &v,
                    const MPI_Comm               &communicator)
    {
      *this = Vector(parallel_partitioner.make_trilinos_map (communicator, true),
                     v);
      owned_elements = parallel_partitioner;
    }



    inline
    Vector &
    Vector::operator= (const TrilinosScalar s)
    {
      AssertIsFinite(s);

      int ierr = vector->PutScalar(s);
      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

      if (nonlocal_vector.get() != nullptr)
        {
          ierr = nonlocal_vector->PutScalar(0.);
          AssertThrow (ierr == 0, ExcTrilinosError(ierr));
        }

      return *this;
    }
  } /* end of namespace MPI */

#endif /* DOXYGEN */

} /* end of namespace TrilinosWrappers */

/*@}*/


namespace internal
{
  namespace LinearOperatorImplementation
  {
    template <typename> class ReinitHelper;

    /**
     * A helper class used internally in linear_operator.h. Specialization for
     * TrilinosWrappers::MPI::Vector.
     */
    template <>
    class ReinitHelper<TrilinosWrappers::MPI::Vector>
    {
    public:
      template <typename Matrix>
      static
      void reinit_range_vector (const Matrix &matrix,
                                TrilinosWrappers::MPI::Vector &v,
                                bool omit_zeroing_entries)
      {
        v.reinit(matrix.locally_owned_range_indices(), matrix.get_mpi_communicator(), omit_zeroing_entries);
      }

      template <typename Matrix>
      static
      void reinit_domain_vector(const Matrix &matrix,
                                TrilinosWrappers::MPI::Vector &v,
                                bool omit_zeroing_entries)
      {
        v.reinit(matrix.locally_owned_domain_indices(), matrix.get_mpi_communicator(), omit_zeroing_entries);
      }
    };

  } /* namespace LinearOperator */
} /* namespace internal */



/**
 * Declare dealii::TrilinosWrappers::MPI::Vector as distributed vector.
 *
 * @author Uwe Koecher, 2017
 */
template <>
struct is_serial_vector< TrilinosWrappers::MPI::Vector > : std::false_type
{
};


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS

/*----------------------------   trilinos_vector.h     ---------------------------*/

#endif // dealii_trilinos_vector_h

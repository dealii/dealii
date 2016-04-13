// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2016 by the deal.II authors
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

#ifndef dealii__trilinos_vector_h
#define dealii__trilinos_vector_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/std_cxx11/shared_ptr.h>
#  include <deal.II/base/subscriptor.h>
#  include <deal.II/base/index_set.h>
#  include <deal.II/base/utilities.h>
#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/vector.h>
#  include <deal.II/lac/trilinos_vector_base.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include "Epetra_Map.h"
#  include "Epetra_LocalMap.h"
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

DEAL_II_NAMESPACE_OPEN


// forward declaration
template <typename> class Vector;

/**
 * @addtogroup TrilinosWrappers
 * @{
 */
namespace TrilinosWrappers
{
  class SparseMatrix;

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
   * This namespace is restricted to vectors only, whereas matrices are always
   * MPI based when run on more than one processor.
   *
   * @ingroup TrilinosWrappers
   * @author Martin Kronbichler, Wolfgang Bangerth, 2008
   */
  namespace MPI
  {
    class BlockVector;

    /**
     * This class implements a wrapper to use the Trilinos distributed vector
     * class Epetra_FEVector. This class is derived from the
     * TrilinosWrappers::VectorBase class and provides all functionality
     * included there.
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
     * TrilinosWrappers::Vector vector;
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
     * in the way we have it in the the non-ghost case view.
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
     * @author Martin Kronbichler, Wolfgang Bangerth, 2008, 2009
     */
    class Vector : public VectorBase
    {
    public:
      /**
       * Declare type for container size.
       */
      typedef dealii::types::global_dof_index size_type;

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
       * @name Basic constructors and initialization.
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

#ifdef DEAL_II_WITH_CXX11
      /**
       * Move constructor. Creates a new vector by stealing the internal data
       * of the vector @p v.
       *
       * @note This constructor is only available if deal.II is configured
       * with C++11 support.
       */
      Vector (Vector &&v);
#endif

      /**
       * Destructor.
       */
      ~Vector ();

      /**
       * Reinit functionality. This function sets the calling vector to the
       * dimension and the parallel distribution of the input vector, but does
       * not copy the elements in <tt>v</tt>. If <tt>omit_zeroing_entries</tt>
       * is not <tt>true</tt>, the elements in the vector are initialized with
       * zero, otherwise the content will be left unchanged and the user has
       * to set all elements.
       *
       * This function has a third argument, <tt>allow_different_maps</tt>,
       * that allows for an exchange of data between two equal-sized vectors
       * (but being distributed differently among the processors). A trivial
       * application of this function is to generate a replication of a whole
       * vector on each machine, when the calling vector is built according to
       * the localized vector class TrilinosWrappers::Vector, and <tt>v</tt>
       * is a distributed vector. In this case, the variable
       * <tt>omit_zeroing_entries</tt> needs to be set to <tt>false</tt>,
       * since it does not make sense to exchange data between differently
       * parallelized vectors without touching the elements.
       */
      void reinit (const VectorBase &v,
                   const bool        omit_zeroing_entries = false,
                   const bool        allow_different_maps = false);

      /**
       * Create vector by merging components from a block vector.
       */
      void reinit (const BlockVector &v,
                   const bool         import_data = false);

      /**
       * Set all components of the vector to the given number @p s. Simply
       * pass this down to the base class, but we still need to declare this
       * function to make the example given in the discussion about making the
       * constructor explicit work.
       */
      Vector &operator= (const TrilinosScalar s);

      /**
       * Copy the given vector. Resize the present vector if necessary. In
       * this case, also the Epetra_Map that designs the parallel partitioning
       * is taken from the input vector.
       */
      Vector &operator= (const Vector &v);

#ifdef DEAL_II_WITH_CXX11
      /**
       * Move the given vector. This operator replaces the present vector with
       * @p v by efficiently swapping the internal data structures.
       *
       * @note This operator is only available if deal.II is configured with
       * C++11 support.
       */
      Vector &operator= (Vector &&v);
#endif

      /**
       * Copy operator from a given localized vector (present on all
       * processes) in TrilinosWrappers format to the current distributed
       * vector. This function assumes that the calling vector (left hand
       * object) already is of the same size as the right hand side vector.
       * Otherwise, an exception will be thrown.
       */
      Vector &operator= (const ::dealii::TrilinosWrappers::Vector &v);

      /**
       * Another copy function. This one takes a deal.II vector and copies it
       * into a TrilinosWrapper vector. Note that since we do not provide any
       * Epetra_map that tells about the partitioning of the vector among the
       * MPI processes, the size of the TrilinosWrapper vector has to be the
       * same as the size of the input vector. In order to change the map, use
       * the reinit(const Epetra_Map &input_map) function.
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
//@}
      /**
       * @name Initialization with an Epetra_Map
       */
//@{
      /**
       * This constructor takes an Epetra_Map that already knows how to
       * distribute the individual components among the MPI processors. Since
       * it also includes information about the size of the vector, this is
       * all we need to generate a parallel vector.
       *
       * Depending on whether the @p parallel_partitioning argument uniquely
       * subdivides elements among processors or not, the resulting vector may
       * or may not have ghost elements. See the general documentation of this
       * class for more information.
       *
       * This function is deprecated.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      explicit Vector (const Epetra_Map &parallel_partitioning)  DEAL_II_DEPRECATED;

      /**
       * Copy constructor from the TrilinosWrappers vector class. Since a
       * vector of this class does not necessarily need to be distributed
       * among processes, the user needs to supply us with an Epetra_Map that
       * sets the partitioning details.
       *
       * Depending on whether the @p parallel_partitioning argument uniquely
       * subdivides elements among processors or not, the resulting vector may
       * or may not have ghost elements. See the general documentation of this
       * class for more information.
       *
       * This function is deprecated.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      Vector (const Epetra_Map &parallel_partitioning,
              const VectorBase &v) DEAL_II_DEPRECATED;

      /**
       * Reinitialize from a deal.II vector. The Epetra_Map specifies the
       * %parallel partitioning.
       *
       * Depending on whether the @p parallel_partitioning argument uniquely
       * subdivides elements among processors or not, the resulting vector may
       * or may not have ghost elements. See the general documentation of this
       * class for more information.
       *
       * This function is deprecated.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      template <typename number>
      void reinit (const Epetra_Map             &parallel_partitioner,
                   const dealii::Vector<number> &v) DEAL_II_DEPRECATED;

      /**
       * Reinit functionality. This function destroys the old vector content
       * and generates a new one based on the input map.
       *
       * Depending on whether the @p parallel_partitioning argument uniquely
       * subdivides elements among processors or not, the resulting vector may
       * or may not have ghost elements. See the general documentation of this
       * class for more information.
       *
       * This function is deprecated.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      void reinit (const Epetra_Map &parallel_partitioning,
                   const bool        omit_zeroing_entries = false) DEAL_II_DEPRECATED;

      /**
       * Copy-constructor from deal.II vectors. Sets the dimension to that of
       * the given vector, and copies all elements.
       *
       * Depending on whether the @p parallel_partitioning argument uniquely
       * subdivides elements among processors or not, the resulting vector may
       * or may not have ghost elements. See the general documentation of this
       * class for more information.
       *
       * This function is deprecated.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      template <typename Number>
      Vector (const Epetra_Map             &parallel_partitioning,
              const dealii::Vector<Number> &v) DEAL_II_DEPRECATED;
//@}
      /**
       * @name Initialization with an IndexSet
       */
//@{
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
              const VectorBase &v,
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
       * Reinit functionality. This function destroys the old vector content
       * and generates a new one based on the input partitioning.  The flag
       * <tt>omit_zeroing_entries</tt> determines whether the vector should be
       * filled with zero (false) or left untouched (true).
       *
       *
       * Depending on whether the @p parallel_partitioning argument uniquely
       * subdivides elements among processors or not, the resulting vector may
       * or may not have ghost elements. See the general documentation of this
       * class for more information.
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
//@}
    };




// ------------------- inline and template functions --------------


    /**
     * Global function @p swap which overloads the default implementation of
     * the C++ standard library which uses a temporary object. The function
     * simply exchanges the data of the two vectors.
     *
     * @relates TrilinosWrappers::MPI::Vector
     * @author Martin Kronbichler, Wolfgang Bangerth, 2008
     */
    inline
    void swap (Vector &u, Vector &v)
    {
      u.swap (v);
    }


#ifndef DOXYGEN

    template <typename number>
    Vector::Vector (const Epetra_Map             &input_map,
                    const dealii::Vector<number> &v)
    {
      reinit (input_map, v);
    }



    template <typename number>
    Vector::Vector (const IndexSet               &parallel_partitioner,
                    const dealii::Vector<number> &v,
                    const MPI_Comm               &communicator)
    {
      *this = Vector(parallel_partitioner.make_trilinos_map (communicator, true),
                     v);
    }




    template <typename number>
    void Vector::reinit (const Epetra_Map             &parallel_partitioner,
                         const dealii::Vector<number> &v)
    {
      if (vector.get() == 0 || vector->Map().SameAs(parallel_partitioner) == false)
        vector.reset (new Epetra_FEVector(parallel_partitioner));

      has_ghosts = vector->Map().UniqueGIDs()==false;

      const int size = parallel_partitioner.NumMyElements();

      // Need to copy out values, since the deal.II might not use doubles, so
      // that a direct access is not possible.
      for (int i=0; i<size; ++i)
        (*vector)[0][i] = v(gid(parallel_partitioner,i));
    }


    inline
    Vector &
    Vector::operator= (const TrilinosScalar s)
    {
      VectorBase::operator= (s);

      return *this;
    }


    template <typename Number>
    Vector &
    Vector::operator= (const ::dealii::Vector<Number> &v)
    {
      if (size() != v.size())
        {
          vector.reset (new Epetra_FEVector(Epetra_Map
                                            (static_cast<TrilinosWrappers::types::int_type>(v.size()), 0,
#ifdef DEAL_II_WITH_MPI
                                             Epetra_MpiComm(MPI_COMM_SELF)
#else
                                             Epetra_SerialComm()
#endif
                                            )));
        }

      reinit (vector_partitioner(), v);
      return *this;
    }


#endif /* DOXYGEN */

  } /* end of namespace MPI */



  /**
   * This class is a specialization of a Trilinos vector to a localized
   * version. The purpose of this class is to provide a copy interface from
   * the possibly parallel Vector class to a local vector on each processor,
   * in order to be able to access all elements in the vector or to apply
   * certain deal.II functions.
   *
   * This class is deprecated, use TrilinosWrappers::MPI::Vector instead.
   *
   * @ingroup TrilinosWrappers
   * @ingroup Vectors
   * @author Martin Kronbichler, 2008
   */
  class Vector : public VectorBase
  {
  public:
    /**
     * Declare type for container size.
     */
    typedef dealii::types::global_dof_index size_type;

    /**
     * A variable that indicates whether this vector supports distributed data
     * storage. If true, then this vector also needs an appropriate compress()
     * function that allows communicating recent set or add operations to
     * individual elements to be communicated to other processors.
     *
     * For the current class, the variable equals false, since it does not
     * support parallel data storage.  If you do need parallel data storage,
     * use TrilinosWrappers::MPI::Vector.
     */
    static const bool supports_distributed_data = false;

    /**
     * Default constructor that generates an empty (zero size) vector. The
     * function <tt>reinit()</tt> will have to give the vector the correct
     * size.
     */
    Vector () DEAL_II_DEPRECATED;

    /**
     * This constructor takes as input the number of elements in the vector.
     */
    explicit Vector (const size_type n) DEAL_II_DEPRECATED;

    /**
     * This constructor takes as input the number of elements in the vector.
     * If the map is not localized, i.e., if there are some elements that are
     * not present on all processes, only the global size of the map will be
     * taken and a localized map will be generated internally. In other words,
     * which element of the @p partitioning argument are set is in fact
     * ignored, the only thing that matters is the size of the index space
     * described by this argument.
     */
    explicit Vector (const Epetra_Map &partitioning) DEAL_II_DEPRECATED;

    /**
     * This constructor takes as input the number of elements in the vector.
     * If the index set is not localized, i.e., if there are some elements
     * that are not present on all processes, only the global size of the
     * index set will be taken and a localized version will be generated
     * internally. In other words, which element of the @p partitioning
     * argument are set is in fact ignored, the only thing that matters is the
     * size of the index space described by this argument.
     */
    explicit Vector (const IndexSet &partitioning,
                     const MPI_Comm &communicator = MPI_COMM_WORLD) DEAL_II_DEPRECATED;

    /**
     * This constructor takes a (possibly parallel) Trilinos Vector and
     * generates a localized version of the whole content on each processor.
     */
    explicit Vector (const VectorBase &V) DEAL_II_DEPRECATED;

    /**
     * Copy-constructor from deal.II vectors. Sets the dimension to that of
     * the given vector, and copies all elements.
     */
    template <typename Number>
    explicit Vector (const dealii::Vector<Number> &v) DEAL_II_DEPRECATED;

    /**
     * Reinit function that resizes the vector to the size specified by
     * <tt>n</tt>.
     */
    void reinit (const size_type n,
                 const bool      omit_zeroing_entries = false);

    /**
     * Initialization with an Epetra_Map. Similar to the call in the other
     * class MPI::Vector, with the difference that now a copy on all processes
     * is generated. This initialization function is appropriate when the data
     * in the localized vector should be imported from a distributed vector
     * that has been initialized with the same communicator. The variable
     * <tt>omit_zeroing_entries</tt> determines whether the vector should be
     * filled with zero or left untouched.
     *
     * Which element of the @p input_map argument are set is in fact ignored,
     * the only thing that matters is the size of the index space described by
     * this argument.
     */
    void reinit (const Epetra_Map &input_map,
                 const bool        omit_zeroing_entries = false);

    /**
     * Initialization with an IndexSet. Similar to the call in the other class
     * MPI::Vector, with the difference that now a copy on all processes is
     * generated. This initialization function is appropriate in case the data
     * in the localized vector should be imported from a distributed vector
     * that has been initialized with the same communicator. The variable
     * <tt>omit_zeroing_entries</tt> determines whether the vector should be
     * filled with zero (false) or left untouched (true).
     *
     * Which element of the @p input_map argument are set is in fact ignored,
     * the only thing that matters is the size of the index space described by
     * this argument.
     */
    void reinit (const IndexSet   &input_map,
                 const MPI_Comm   &communicator = MPI_COMM_WORLD,
                 const bool        omit_zeroing_entries = false);

    /**
     * Reinit function. Takes the information of a Vector and copies
     * everything to the calling vector, now also allowing different maps.
     */
    void reinit (const VectorBase &V,
                 const bool        omit_zeroing_entries = false,
                 const bool        allow_different_maps = false);

    /**
     * Set all components of the vector to the given number @p s. Simply pass
     * this down to the base class, but we still need to declare this function
     * to make the example given in the discussion about making the
     * constructor explicit work.
     */
    Vector &operator= (const TrilinosScalar s);

    /**
     * Sets the left hand argument to the (parallel) Trilinos Vector.
     * Equivalent to the @p reinit function.
     */
    Vector &operator= (const MPI::Vector &v);

    /**
     * Sets the left hand argument to the deal.II vector.
     */
    template <typename Number>
    Vector &operator= (const ::dealii::Vector<Number> &v);

    /**
     * Copy operator. Copies both the dimension and the content in the right
     * hand argument.
     */
    Vector &operator= (const Vector &v);

    /**
     * This function does nothing but is there for compatibility with the @p
     * PETScWrappers::Vector class.
     *
     * For the PETSc vector wrapper class, this function updates the ghost
     * values of the PETSc vector. This is necessary after any modification
     * before reading ghost values.
     *
     * However, for the implementation of this class, it is immaterial and
     * thus an empty function.
     */
    void update_ghost_values () const;
  };



// ------------------- inline and template functions --------------


  /**
   * Global function @p swap which overloads the default implementation of the
   * C++ standard library which uses a temporary object. The function simply
   * exchanges the data of the two vectors.
   *
   * @relates TrilinosWrappers::Vector
   * @author Martin Kronbichler, Wolfgang Bangerth, 2008
   */
  inline
  void swap (Vector &u, Vector &v)
  {
    u.swap (v);
  }


#ifndef DOXYGEN

  template <typename number>
  Vector::Vector (const dealii::Vector<number> &v)
  {
    Epetra_LocalMap map ((TrilinosWrappers::types::int_type)v.size(), 0, Utilities::Trilinos::comm_self());
    vector.reset (new Epetra_FEVector(map));
    *this = v;
  }



  inline
  Vector &
  Vector::operator= (const TrilinosScalar s)
  {
    VectorBase::operator= (s);

    return *this;
  }



  template <typename Number>
  Vector &
  Vector::operator= (const ::dealii::Vector<Number> &v)
  {
    if (size() != v.size())
      {
        vector.reset();

        Epetra_LocalMap map ((TrilinosWrappers::types::int_type)v.size(), 0,
                             Utilities::Trilinos::comm_self());
        vector.reset (new Epetra_FEVector(map));
      }

    const Epetra_Map &map = vector_partitioner();
    const TrilinosWrappers::types::int_type size = map.NumMyElements();

    Assert (map.MaxLID() == size-1,
            ExcDimensionMismatch(map.MaxLID(), size-1));

    // Need to copy out values, since the
    // deal.II might not use doubles, so
    // that a direct access is not possible.
    for (TrilinosWrappers::types::int_type i=0; i<size; ++i)
      (*vector)[0][i] = v(i);

    return *this;
  }



  inline
  void
  Vector::update_ghost_values () const
  {}


#endif

} /* namespace TrilinosWrappers */

/*@}*/


namespace internal
{
  namespace LinearOperator
  {
    template <typename> class ReinitHelper;

    /**
     * A helper class used internally in linear_operator.h. Specialization for
     * TrilinosWrappers::MPI::Vector.
     */
    template<>
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

    /**
     * A helper class used internally in linear_operator.h. Specialization for
     * TrilinosWrappers::Vector.
     */
    template<>
    class ReinitHelper<TrilinosWrappers::Vector>
    {
    public:
      template <typename Matrix>
      static
      void reinit_range_vector (const Matrix &matrix,
                                TrilinosWrappers::Vector &v,
                                bool omit_zeroing_entries)
      {
        v.reinit(matrix.locally_owned_range_indices(),
                 matrix.get_mpi_communicator(),
                 omit_zeroing_entries);
      }

      template <typename Matrix>
      static
      void reinit_domain_vector(const Matrix &matrix,
                                TrilinosWrappers::Vector &v,
                                bool omit_zeroing_entries)
      {
        v.reinit(matrix.locally_owned_domain_indices(),
                 matrix.get_mpi_communicator(),
                 omit_zeroing_entries);
      }
    };

  } /* namespace LinearOperator */
} /* namespace internal */


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS

/*----------------------------   trilinos_vector.h     ---------------------------*/

#endif
/*----------------------------   trilinos_vector.h     ---------------------------*/

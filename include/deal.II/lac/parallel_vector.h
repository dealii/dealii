// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2014 by the deal.II authors
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

#ifndef __deal2__parallel_vector_h
#define __deal2__parallel_vector_h

#include <deal.II/base/config.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/types.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/lac/vector_view.h>

#include <cstring>
#include <iomanip>


DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_PETSC
namespace PETScWrappers
{
  namespace MPI
  {
    class Vector;
  }
}
#endif

#ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers
{
  namespace MPI
  {
    class Vector;
  }
}
#endif


namespace parallel
{
  namespace distributed
  {
    template <typename Number> class BlockVector;

    /*! @addtogroup Vectors
     *@{
     */


    /**
     * Implementation of a parallel vector class. The design of this class is
     * similar to the standard ::dealii::Vector class in deal.II, with the
     * exception that storage is distributed with MPI.
     *
     * The vector is designed for the following scheme of parallel partitioning:
     * - The indices held by individual processes (locally owned part) in the
     *   MPI parallelization form a contiguous range
     *   <code>[my_first_index,my_last_index)</code>.
     * - Ghost indices residing on arbitrary positions of other processors are
     *   allowed. It is in general more efficient if ghost indices are
     *   clustered, since they are stored as a set of intervals. The
     *   communication pattern of the ghost indices is determined when calling
     *   the function <code>reinit (locally_owned, ghost_indices,
     *   communicator)</code>, and retained until the partitioning is changed
     *   again. This allows for efficient parallel communication of indices. In
     *   particular, it stores the communication pattern, rather than having
     *   to compute it again for every communication.
     * - Besides the usual global access operator () it is also possible to
     *   access vector entries in the local index space with the function @p
     *   local_element(). Locally owned indices are placed first, [0,
     *   local_size()), and then all ghost indices follow after them
     *   contiguously, [local_size(), local_size()+n_ghost_entries()).
     *
     * Functions related to parallel functionality:
     * - The function <code>compress()</code> goes through the data associated
     *   with ghost indices and communicates it to the owner process, which can
     *   then add it to the correct position. This can be used e.g. after
     *   having run an assembly routine involving ghosts that fill this vector.
     *   Note that the @p insert mode of @p compress() does not set the
     *   elements included in ghost entries but simply discards them, assuming
     *   that the owning processor has set them to the desired value already.
     * - The <code>update_ghost_values()</code> function imports the data from
     *   the owning processor to the ghost indices in order to provide read
     *   access to the data associated with ghosts.
     * - It is possible to split the above functions into two phases, where
     *   the first initiates the communication and the second one finishes
     *   it. These functions can be used to overlap communication with
     *   computations in other parts of the code.
     * - Of course, reduction operations (like norms) make use of collective
     *   all-to-all MPI communications.
     *
     * This vector can take two different states with respect to ghost
     * elements:
     * - After creation and whenever zero_out_ghosts() is called (or
     *   <code>operator = (0.)</code>), the vector does only allow writing
     *   into ghost elements but not reading from ghost elements.
     * - After a call to update_ghost_values(), the vector does not allow
     *   writing into ghost elements but only reading from them. This is in
     *   order to avoid undesired ghost data artifacts when calling compress()
     *   after modifying some vector entries.
     * The current status of the ghost entries (read mode or write mode) can
     * be queried by the method has_ghost_elements(), which returns
     * <code>true</code> exactly when ghost elements have been updated and
     * <code>false</code> otherwise, irrespective of the actual number of
     * ghost entries in the vector layout (for that information, use
     * n_ghost_entries() instead).
     *
     * This vector uses the facilities of the class dealii::Vector<Number> for
     * implementing the operations on the local range of the vector. In
     * particular, it also inherits thread parallelism that splits most
     * vector-vector operations into smaller chunks if the program uses
     * multiple threads. This may or may not be desired when working also with
     * MPI.
     *
     * @author Katharina Kormann, Martin Kronbichler, 2010, 2011
     */
    template <typename Number>
    class Vector : public Subscriptor
    {
    public:
      /**
       * Declare standard types used in all containers. These types parallel
       * those in the <tt>C++</tt> standard libraries <tt>vector<...></tt>
       * class.
       */
      typedef Number                                            value_type;
      typedef value_type                                       *pointer;
      typedef const value_type                                 *const_pointer;
      typedef value_type                                       *iterator;
      typedef const value_type                                 *const_iterator;
      typedef value_type                                       &reference;
      typedef const value_type                                 &const_reference;
      typedef types::global_dof_index                           size_type;
      typedef typename numbers::NumberTraits<Number>::real_type real_type;

      /**
       * A variable that indicates whether this vector
       * supports distributed data storage. If true, then
       * this vector also needs an appropriate compress()
       * function that allows communicating recent set or
       * add operations to individual elements to be communicated
       * to other processors.
       *
       * For the current class, the variable equals
       * true, since it does support parallel data storage.
       */
      static const bool supports_distributed_data = true;

      /**
       * @name 1: Basic Object-handling
       */
      //@{
      /**
       * Empty constructor.
       */
      Vector ();

      /**
       * Copy constructor. Uses the parallel partitioning of @p in_vector.
       */
      Vector (const Vector<Number> &in_vector);

      /**
       * Constructs a parallel vector of the given global size without any
       * actual parallel distribution.
       */
      Vector (const size_type size);

      /**
       * Constructs a parallel vector. The local range is specified by @p
       * locally_owned_set (note that this must be a contiguous interval,
       * multiple intervals are not possible). The IndexSet @p ghost_indices
       * specifies ghost indices, i.e., indices which one might need to read
       * data from or accumulate data from. It is allowed that the set of
       * ghost indices also contains the local range, but it does not need to.
       *
       * This function involves global communication, so it should only be
       * called once for a given layout. Use the constructor with
       * Vector<Number> argument to create additional vectors with the same
       * parallel layout.
       */
      Vector (const IndexSet &local_range,
              const IndexSet &ghost_indices,
              const MPI_Comm  communicator);

      /**
       * Same constructor as above but without any ghost indices.
       */
      Vector (const IndexSet &local_range,
              const MPI_Comm  communicator);

      /**
       * Create the vector based on the parallel partitioning described in @p
       * partitioner. The input argument is a shared pointer, which store the
       * partitioner data only once and share it between several vectors with
       * the same layout.
       */
      Vector (const std_cxx11::shared_ptr<const Utilities::MPI::Partitioner> &partitioner);

      /**
       * Destructor.
       */
      ~Vector ();

      /**
       * Sets the global size of the vector to @p size without any actual
       * parallel distribution.
       */
      void reinit (const size_type size,
                   const bool      fast = false);

      /**
       * Uses the parallel layout of the input vector @p in_vector and
       * allocates memory for this vector. Recommended initialization function
       * when several vectors with the same layout should be created.
       *
       * If the flag @p fast is set to false, the memory will be initialized
       * with zero, otherwise the memory will be untouched (and the user must
       * make sure to fill it with reasonable data before using it).
       */
      template <typename Number2>
      void reinit(const Vector<Number2> &in_vector,
                  const bool             fast = false);

      /**
       * Initialize the vector. The local range is specified by @p
       * locally_owned_set (note that this must be a contiguous interval,
       * multiple intervals are not possible). The IndexSet @p ghost_indices
       * specifies ghost indices, i.e., indices which one might need to read
       * data from or accumulate data from. It is allowed that the set of
       * ghost indices also contains the local range, but it does not need to.
       *
       * This function involves global communication, so it should only be
       * called once for a given layout. Use the @p reinit function with
       * Vector<Number> argument to create additional vectors with the same
       * parallel layout.
       */
      void reinit (const IndexSet &local_range,
                   const IndexSet &ghost_indices,
                   const MPI_Comm  communicator);

      /**
       * Same as above, but without ghost entries.
       */
      void reinit (const IndexSet &local_range,
                   const MPI_Comm  communicator);

      /**
       * Initialize the vector given to the parallel partitioning described in
       * @p partitioner. The input argument is a shared pointer, which store
       * the partitioner data only once and share it between several vectors
       * with the same layout.
       */
      void reinit (const std_cxx11::shared_ptr<const Utilities::MPI::Partitioner> &partitioner);

      /**
       * Swap the contents of this vector and the other vector @p v. One could
       * do this operation with a temporary variable and copying over the data
       * elements, but this function is significantly more efficient since it
       * only swaps the pointers to the data of the two vectors and therefore
       * does not need to allocate temporary storage and move data around.
       *
       * This function is analog to the the @p swap function of all C++
       * standard containers. Also, there is a global function
       * <tt>swap(u,v)</tt> that simply calls <tt>u.swap(v)</tt>, again in
       * analogy to standard functions.
       *
       * This function is virtual in order to allow for derived classes to
       * handle memory separately.
       */
      void swap (Vector<Number> &v);

      /**
       * Assigns the vector to the parallel partitioning of the input vector
       * @p in_vector, and copies all the data.
       */
      Vector<Number> &
      operator = (const Vector<Number>  &in_vector);

      /**
       * Assigns the vector to the parallel partitioning of the input vector
       * @p in_vector, and copies all the data.
       */
      template <typename Number2>
      Vector<Number> &
      operator = (const Vector<Number2> &in_vector);

#ifdef DEAL_II_WITH_PETSC
      /**
       * Copy the content of a PETSc vector into the calling vector. This
       * function assumes that the vectors layouts have already been
       * initialized to match.
       *
       * This operator is only available if deal.II was configured with PETSc.
       */
      Vector<Number> &
      operator = (const PETScWrappers::MPI::Vector &petsc_vec);
#endif

#ifdef DEAL_II_WITH_TRILINOS
      /**
       * Copy the content of a Trilinos vector into the calling vector. This
       * function assumes that the vectors layouts have already been
       * initialized to match.
       *
       * This operator is only available if deal.II was configured with
       * Trilinos.
       */
      Vector<Number> &
      operator = (const TrilinosWrappers::MPI::Vector &trilinos_vec);
#endif

      /**
       * This method copies the local range from another vector with the same
       * local range, but possibly different layout of ghost indices.
       */
      void copy_from (const Vector<Number> &in_vector,
                      const bool            call_update_ghost_values = false);

      /**
       * Sets all elements of the vector to the scalar @p s. If the scalar is
       * zero, also ghost elements are set to zero, otherwise they remain
       * unchanged.
       */
      Vector<Number> &operator = (const Number s);

      /**
       * This function copies the data that has accumulated in the data buffer
       * for ghost indices to the owning processor. For the meaning of the
       * argument @p operation, see the entry on @ref GlossCompress
       * "Compressing distributed vectors and matrices" in the glossary.
       *
       * There are two variants for this function. If called with argument @p
       * VectorOperation::add adds all the data accumulated in ghost elements
       * to the respective elements on the owning processor and clears the
       * ghost array afterwards. If called with argument @p
       * VectorOperation::insert, a set operation is performed. Since setting
       * elements in a vector with ghost elements is ambiguous (as one can set
       * both the element on the ghost site as well as the owning site), this
       * operation makes the assumption that all data is set correctly on the
       * owning processor. Upon call of compress(VectorOperation::insert), all
       * ghost entries are thus simply zeroed out (using
       * zero_ghost_values()). In debug mode, a check is performed for whether
       * the data set is actually consistent between processors,
       * i.e., whenever a non-zero ghost element is found, it is compared to
       * the value on the owning processor and an exception is thrown if these
       * elements do not agree.
       *
       */
      void compress (::dealii::VectorOperation::values operation);

      /**
       * @deprecated: use compress(VectorOperation::values) instead.
       */
      void compress () DEAL_II_DEPRECATED;

      /**
       * Fills the data field for ghost indices with the values stored in the
       * respective positions of the owning processor. This function is needed
       * before reading from ghosts. The function is @p const even though
       * ghost data is changed. This is needed to allow functions with a @p
       * const vector to perform the data exchange without creating
       * temporaries.
       *
       * After calling this method, write access to ghost elements of the
       * vector is forbidden and an exception is thrown. Only read access to
       * ghost elements is allowed in this state. Note that all subsequent
       * operations on this vector, like global vector addition, etc., will
       * also update the ghost values by a call to this method after the
       * operation. However, global reduction operations like norms or the
       * inner product will always ignore ghost elements in order to avoid
       * counting the ghost data more than once. To allow writing to ghost
       * elements again, call zero_out_ghosts().
       */
      void update_ghost_values () const;

      /**
       * Initiates communication for the @p compress() function with
       * non-blocking communication. This function does not wait for the
       * transfer to finish, in order to allow for other computations during
       * the time it takes until all data arrives.
       *
       * Before the data is actually exchanged, the function must be followed
       * by a call to @p compress_finish().
       *
       * In case this function is called for more than one vector before @p
       * compress_finish() is invoked, it is mandatory to specify a unique
       * communication channel to each such call, in order to avoid several
       * messages with the same ID that will corrupt this operation.
       */
      void compress_start (const unsigned int communication_channel = 0,
                           ::dealii::VectorOperation::values operation = VectorOperation::add);

      /**
       * For all requests that have been initiated in compress_start, wait for
       * the communication to finish. Once it is finished, add or set the data
       * (depending on the flag operation) to the respective positions in the
       * owning processor, and clear the contents in the ghost data
       * fields. The meaning of this argument is the same as in
       * compress().
       *
       * This function should be called exactly once per vector after calling
       * compress_start, otherwise the result is undefined. In particular, it
       * is not well-defined to call compress_start on the same vector again
       * before compress_finished has been called. However, there is no
       * warning to prevent this situation.
       *
       * Must follow a call to the @p compress_start function.
       */
      void compress_finish (::dealii::VectorOperation::values operation);

      /**
       * @deprecated: use compress_finish(VectorOperation::values) instead.
       */
      void compress_finish (const bool add_ghost_data = true) DEAL_II_DEPRECATED;

      /**
       * Initiates communication for the @p update_ghost_values() function
       * with non-blocking communication. This function does not wait for the
       * transfer to finish, in order to allow for other computations during
       * the time it takes until all data arrives.
       *
       * Before the data is actually exchanged, the function must be followed
       * by a call to @p update_ghost_values_finish().
       *
       * In case this function is called for more than one vector before @p
       * update_ghost_values_finish() is invoked, it is mandatory to specify a
       * unique communication channel to each such call, in order to avoid
       * several messages with the same ID that will corrupt this operation.
       */
      void update_ghost_values_start (const unsigned int communication_channel = 0) const;


      /**
       * For all requests that have been started in update_ghost_values_start,
       * wait for the communication to finish.
       *
       * Must follow a call to the @p update_ghost_values_start function
       * before reading data from ghost indices.
       */
      void update_ghost_values_finish () const;

      /**
       * This method zeros the entries on ghost dofs, but does not touch
       * locally owned DoFs.
       *
       * After calling this method, read access to ghost elements of the
       * vector is forbidden and an exception is thrown. Only write access to
       * ghost elements is allowed in this state.
       */
      void zero_out_ghosts ();

      /**
       * Returns whether the vector currently is in a state where ghost values
       * can be read or not. This is the same functionality as other parallel
       * vectors have. If this method returns false, this only means that
       * read-access to ghost elements is prohibited whereas write access is
       * still possible (to those entries specified as ghosts during
       * initialization), not that there are no ghost elements at all.
       */
      bool has_ghost_elements() const;

      /**
       * Return whether the vector contains only elements with value
       * zero. This is a collective operation. This function is expensive, because
       * potentially all elements have to be checked.
       */
      bool all_zero () const;

      /**
       * Return @p true if the vector has no negative entries, i.e. all
       * entries are zero or positive. This function is used, for example, to
       * check whether refinement indicators are really all positive (or
       * zero).
       *
       * The function obviously only makes sense if the template argument of
       * this class is a real type. If it is a complex type, then an exception
       * is thrown.
       */
      bool is_non_negative () const;

      /**
       * Checks for equality of the two vectors.
       */
      template <typename Number2>
      bool operator == (const Vector<Number2> &v) const;

      /**
       * Checks for inequality of the two vectors.
       */
      template <typename Number2>
      bool operator != (const Vector<Number2> &v) const;

      /**
       * Perform the inner product of two vectors.
       */
      template <typename Number2>
      Number operator * (const Vector<Number2> &V) const;

      /**
       * Computes the square of the l<sub>2</sub> norm of the vector (i.e.,
       * the sum of the squares of all entries among all processors).
       */
      real_type norm_sqr () const;

      /**
       * Computes the mean value of all the entries in the vector.
       */
      Number mean_value () const;

      /**
       * Returns the l<sub>1</sub> norm of the vector (i.e., the sum of the
       * absolute values of all entries among all processors).
       */
      real_type l1_norm () const;

      /**
       * Returns the l<sub>2</sub> norm of the vector (i.e., square root of
       * the sum of the square of all entries among all processors).
       */
      real_type l2_norm () const;

      /**
       * Returns the l<sub>p</sub> norm with real @p p of the vector (i.e.,
       * the pth root of sum of the pth power of all entries among all
       * processors).
       */
      real_type lp_norm (const real_type p) const;

      /**
       * Returns the maximum norm of the vector (i.e., maximum absolute value
       * among all entries among all processors).
       */
      real_type linfty_norm () const;

      /**
       * Returns the global size of the vector, equal to the sum of the number
       * of locally owned indices among all the processors.
       */
      size_type size () const;

      /**
       * Returns the local size of the vector, i.e., the number of indices
       * owned locally.
       */
      size_type local_size() const;

      /**
       * Returns the half-open interval that specifies the locally owned range
       * of the vector. Note that <code>local_size() == local_range().second -
       * local_range().first</code>.
       */
      std::pair<size_type, size_type> local_range () const;

      /**
       * Returns true if the given global index is in the local range of this
       * processor.
       */
      bool in_local_range (const size_type global_index) const;

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
       * Returns the number of ghost elements present on the vector.
       */
      size_type n_ghost_entries () const;

      /**
       * Return an index set that describes which elements of this vector are
       * not owned by the current processor but can be written into or read
       * from locally (ghost elements).
       */
      const IndexSet &ghost_elements() const;

      /**
       * Returns whether the given global index is a ghost index on the
       * present processor. Returns false for indices that are owned locally
       * and for indices not present at all.
       */
      bool is_ghost_entry (const types::global_dof_index global_index) const;

      /**
       * Make the @p Vector class a bit like the <tt>vector<></tt> class of
       * the C++ standard library by returning iterators to the start and end
       * of the <i>locally owned</i> elements of this vector.
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
       * Return a constant iterator pointing to the element past the end of
       * the array of the locally owned entries.
       */
      const_iterator end () const;
      //@}


      /**
       * @name 2: Data-Access
       */
      //@{

      /**
       * Read access to the data in the position corresponding to @p
       * global_index. The index must be either in the local range of the
       * vector or be specified as a ghost index at construction.
       *
       * Performance: <tt>O(1)</tt> for locally owned elements that represent
       * a contiguous range and <tt>O(log(n<sub>ranges</sub>))</tt> for ghost
       * elements (quite fast, but slower than local_element()).
       */
      Number operator () (const size_type global_index) const;

      /**
       * Read and write access to the data in the position corresponding to @p
       * global_index. The index must be either in the local range of the
       * vector or be specified as a ghost index at construction.
       *
       * Performance: <tt>O(1)</tt> for locally owned elements that represent
       * a contiguous range and <tt>O(log(n<sub>ranges</sub>))</tt> for ghost
       * elements (quite fast, but slower than local_element()).
       */
      Number &operator () (const size_type global_index);

      /**
       * Read access to the data in the position corresponding to @p
       * global_index. The index must be either in the local range of the
       * vector or be specified as a ghost index at construction.
       *
       * This function does the same thing as operator().
       */
      Number operator [] (const size_type global_index) const;

      /**
       * Read and write access to the data in the position corresponding to @p
       * global_index. The index must be either in the local range of the
       * vector or be specified as a ghost index at construction.
       *
       * This function does the same thing as operator().
       */
      Number &operator [] (const size_type global_index);

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
       * Read access to the data field specified by @p local_index. Locally
       * owned indices can be accessed with indices
       * <code>[0,local_size)</code>, and ghost indices with indices
       * <code>[local_size,local_size+ n_ghost_entries]</code>.
       *
       * Performance: Direct array access (fast).
       */
      Number local_element (const size_type local_index) const;

      /**
       * Read and write access to the data field specified by @p
       * local_index. Locally owned indices can be accessed with indices
       * <code>[0,local_size)</code>, and ghost indices with indices
       * <code>[local_size,local_size+n_ghosts]</code>.
       *
       * Performance: Direct array access (fast).
       */
      Number &local_element (const size_type local_index);
      //@}


      /**
       * @name 3: Modification of vectors
       */
      //@{

      /**
       * Add the given vector to the present one.
       */
      Vector<Number> &operator += (const Vector<Number> &V);

      /**
       * Subtract the given vector from the present one.
       */
      Vector<Number> &operator -= (const Vector<Number> &V);

      /**
       * A collective add operation: This funnction adds a whole set of values
       * stored in @p values to the vector components specified by @p indices.
       */
      template <typename OtherNumber>
      void add (const std::vector<size_type>   &indices,
                const std::vector<OtherNumber>  &values);

      /**
       * This is a second collective add operation. As a difference, this
       * function takes a deal.II vector of values.
       */
      template <typename OtherNumber>
      void add (const std::vector<size_type>        &indices,
                const ::dealii::Vector<OtherNumber> &values);

      /**
       * Take an address where <tt>n_elements</tt> are stored contiguously and
       * add them into the vector. Handles all cases which are not covered by
       * the other two <tt>add()</tt> functions above.
       */
      template <typename OtherNumber>
      void add (const size_type    n_elements,
                const size_type   *indices,
                const OtherNumber  *values);

      /**
       * Addition of @p s to all components. Note that @p s is a scalar and
       * not a vector.
       */
      void add (const Number s);

      /**
       * Simple vector addition, equal to the <tt>operator +=</tt>.
       */
      void add (const Vector<Number> &V);

      /**
       * Simple addition of a multiple of a vector, i.e. <tt>*this +=
       * a*V</tt>.
       */
      void add (const Number a, const Vector<Number> &V);

      /**
       * Multiple addition of scaled vectors, i.e. <tt>*this += a*V+b*W</tt>.
       */
      void add (const Number a, const Vector<Number> &V,
                const Number b, const Vector<Number> &W);

      /**
       * Scaling and simple vector addition, i.e.  <tt>*this =
       * s*(*this)+V</tt>.
       */
      void sadd (const Number          s,
                 const Vector<Number> &V);

      /**
       * Scaling and simple addition, i.e.  <tt>*this = s*(*this)+a*V</tt>.
       */
      void sadd (const Number          s,
                 const Number          a,
                 const Vector<Number> &V);

      /**
       * Scaling and multiple addition.
       */
      void sadd (const Number          s,
                 const Number          a,
                 const Vector<Number> &V,
                 const Number          b,
                 const Vector<Number> &W);

      /**
       * Scaling and multiple addition.  <tt>*this = s*(*this)+a*V + b*W +
       * c*X</tt>.
       */
      void sadd (const Number          s,
                 const Number          a,
                 const Vector<Number> &V,
                 const Number          b,
                 const Vector<Number> &W,
                 const Number          c,
                 const Vector<Number> &X);

      /**
       * Scale each element of the vector by the given factor.
       *
       * @deprecated This function is deprecated and will be removed in a
       * future version. Use <tt>operator *=</tt> and <tt>operator /=</tt>
       * instead.
       */
      void scale (const Number factor) DEAL_II_DEPRECATED;

      /**
       * Scale each element of the vector by a constant value.
       */
      Vector<Number> &operator *= (const Number factor);

      /**
       * Scale each element of the vector by the inverse of the given value.
       */
      Vector<Number> &operator /= (const Number factor);

      /**
       * Scale each element of this vector by the corresponding element in the
       * argument. This function is mostly meant to simulate multiplication
       * (and immediate re-assignment) by a diagonal scaling matrix.
       */
      void scale (const Vector<Number> &scaling_factors);

      /**
       * Scale each element of this vector by the corresponding element in the
       * argument. This function is mostly meant to simulate multiplication
       * (and immediate re-assignment) by a diagonal scaling matrix.
       */
      template <typename Number2>
      void scale (const Vector<Number2> &scaling_factors);

      /**
       * Assignment <tt>*this = a*u</tt>.
       */
      void equ (const Number a, const Vector<Number> &u);

      /**
       * Assignment <tt>*this = a*u</tt>.
       */
      template <typename Number2>
      void equ (const Number a, const Vector<Number2> &u);

      /**
       * Assignment <tt>*this = a*u + b*v</tt>.
       */
      void equ (const Number a, const Vector<Number> &u,
                const Number b, const Vector<Number> &v);

      /**
       * Assignment <tt>*this = a*u + b*v + b*w</tt>.
       */
      void equ (const Number a, const Vector<Number> &u,
                const Number b, const Vector<Number> &v,
                const Number c, const Vector<Number> &w);

      /**
       * Compute the elementwise ratio of the two given vectors, that is let
       * <tt>this[i] = a[i]/b[i]</tt>. This is useful for example if you want
       * to compute the cellwise ratio of true to estimated error.
       *
       * This vector is appropriately scaled to hold the result.
       *
       * If any of the <tt>b[i]</tt> is zero, the result is undefined. No
       * attempt is made to catch such situations.
       */
      void ratio (const Vector<Number> &a,
                  const Vector<Number> &b);
      //@}


      /**
       * @name 4: Mixed stuff
       */
      //@{
      /**
       * Return a reference to the MPI communicator object in use with this
       * vector.
       */
      const MPI_Comm &get_mpi_communicator () const;

      /**
       * Checks whether the given partitioner is compatible with the
       * partitioner used for this vector. Two partitioners are compatible if
       * the have the same local size and the same ghost indices. They do not
       * necessarily need to be the same data field. This is a local operation
       * only, i.e., if only some processors decide that the partitioning is
       * not compatible, only these processors will return @p false, whereas
       * the other processors will return @p true.
       */
      bool
      partitioners_are_compatible (const Utilities::MPI::Partitioner &part) const;


      /**
       * Prints the vector to the output stream @p out.
       */
      void print (std::ostream       &out,
                  const unsigned int  precision  = 3,
                  const bool          scientific = true,
                  const bool          across     = true) const;

      /**
       * Returns the memory consumption of this class in bytes.
       */
      std::size_t memory_consumption () const;
      //@}

      /**
       * Exception
       */
      DeclException3 (ExcNonMatchingElements,
                      double, double, unsigned int,
                      << "Called compress(VectorOperation::insert), but"
                      << " the element received from a remote processor, value "
                      << std::setprecision(16) << arg1
                      << ", does not match with the value "
                      << std::setprecision(16) << arg2
                      << " on the owner processor " << arg3);

    private:
      /**
       * Local part of all_zero().
       */
      bool all_zero_local () const;

      /**
       * Local part of is_non_negative().
       */
      bool is_non_negative_local () const;

      /**
       * Local part of operator==.
       */
      template <typename Number2>
      bool vectors_equal_local (const Vector<Number2> &v) const;

      /**
       * Local part of the inner product of two vectors.
       */
      template <typename Number2>
      Number inner_product_local (const Vector<Number2> &V) const;

      /**
       * Local part of norm_sqr().
       */
      real_type norm_sqr_local () const;

      /**
       * Local part of mean_value().
       */
      Number mean_value_local () const;

      /**
       * Local part of l1_norm().
       */
      real_type l1_norm_local () const;

      /**
       * Local part of lp_norm().
       */
      real_type lp_norm_local (const real_type p) const;

      /**
       * Local part of linfty_norm().
       */
      real_type linfty_norm_local () const;

      /**
       * Shared pointer to store the parallel partitioning information. This
       * information can be shared between several vectors that have the same
       * partitioning.
       */
      std_cxx11::shared_ptr<const Utilities::MPI::Partitioner> partitioner;

      /**
       * The size that is currently allocated in the val array.
       */
      size_type allocated_size;

      /**
       * Pointer to the array of local elements of this vector.
       */
      Number         *val;

      /**
       * Temporary storage that holds the data that is sent to this processor
       * in @p compress() or sent from this processor in @p
       * update_ghost_values.
       */
      mutable Number *import_data;

      /**
       * Stores whether the vector currently allows for reading ghost elements
       * or not. Note that this is to ensure consistent ghost data and does
       * not indicate whether the vector actually can store ghost elements. In
       * particular, when assembling a vector we do not allow reading
       * elements, only writing them.
       */
      mutable bool vector_is_ghosted;

      /**
       * Provide this class with all functionality of ::dealii::Vector by
       * creating a VectorView object.
       */
      VectorView<Number> vector_view;

#ifdef DEAL_II_WITH_MPI
      /**
       * A vector that collects all requests from @p compress()
       * operations. This class uses persistent MPI communicators, i.e., the
       * communication channels are stored during successive calls to a given
       * function. This reduces the overhead involved with setting up the MPI
       * machinery, but it does not remove the need for a receive operation to
       * be posted before the data can actually be sent.
       */
      std::vector<MPI_Request>   compress_requests;

      /**
       * A vector that collects all requests from @p update_ghost_values()
       * operations. This class uses persistent MPI communicators.
       */
      mutable std::vector<MPI_Request>   update_ghost_values_requests;
#endif

      /**
       * A lock that makes sure that the @p compress and @p
       * update_ghost_values functions give reasonable results also when used
       * with several threads.
       */
      mutable Threads::Mutex mutex;

      /**
       * A helper function that clears the compress_requests and
       * update_ghost_values_requests field. Used in reinit functions.
       */
      void clear_mpi_requests ();

      /**
       * A helper function that is used to resize the val array.
       */
      void resize_val (const size_type new_allocated_size);

      /*
       * Make all other vector types friends.
       */
      template <typename Number2> friend class Vector;

      /**
       * Make BlockVector type friends.
       */
      template <typename Number2> friend class BlockVector;
    };

    /*@}*/


    /*----------------------- Inline functions ----------------------------------*/

#ifndef DOXYGEN

    template <typename Number>
    inline
    Vector<Number>::Vector ()
      :
      partitioner (new Utilities::MPI::Partitioner()),
      allocated_size (0),
      val (0),
      import_data (0),
      vector_is_ghosted (false),
      vector_view (0, static_cast<Number *>(0))
    {}



    template <typename Number>
    inline
    Vector<Number>::Vector (const Vector<Number> &v)
      :
      Subscriptor(),
      allocated_size (0),
      val (0),
      import_data (0),
      vector_is_ghosted (false),
      vector_view (0, static_cast<Number *>(0))
    {
      reinit (v, true);
      vector_view = v.vector_view;
      zero_out_ghosts();
    }



    template <typename Number>
    inline
    Vector<Number>::Vector (const IndexSet &local_range,
                            const IndexSet &ghost_indices,
                            const MPI_Comm  communicator)
      :
      allocated_size (0),
      val (0),
      import_data (0),
      vector_is_ghosted (false),
      vector_view (0, static_cast<Number *>(0))
    {
      reinit (local_range, ghost_indices, communicator);
    }



    template <typename Number>
    inline
    Vector<Number>::Vector (const IndexSet &local_range,
                            const MPI_Comm  communicator)
      :
      allocated_size (0),
      val (0),
      import_data (0),
      vector_is_ghosted (false),
      vector_view (0, static_cast<Number *>(0))
    {
      IndexSet ghost_indices(local_range.size());
      reinit (local_range, ghost_indices, communicator);
    }



    template <typename Number>
    inline
    Vector<Number>::Vector (const size_type size)
      :
      allocated_size (0),
      val (0),
      import_data (0),
      vector_is_ghosted (false),
      vector_view (0, static_cast<Number *>(0))
    {
      reinit (size, false);
    }



    template <typename Number>
    inline
    Vector<Number>::
    Vector (const std_cxx11::shared_ptr<const Utilities::MPI::Partitioner> &partitioner)
      :
      allocated_size (0),
      val (0),
      import_data (0),
      vector_is_ghosted (false),
      vector_view (0, static_cast<Number *>(0))
    {
      reinit (partitioner);
    }



    template <typename Number>
    inline
    Vector<Number>::~Vector ()
    {
      if (val != 0)
        delete[] val;
      val = 0;

      if (import_data != 0)
        delete[] import_data;
      import_data = 0;

      clear_mpi_requests();
    }



    template <typename Number>
    inline
    Vector<Number> &
    Vector<Number>::operator = (const Vector<Number> &c)
    {
      Assert (c.partitioner.get() != 0, ExcNotInitialized());

      // we update ghost values whenever one of the input or output vector
      // already held ghost values or when we import data from a vector with
      // the same local range but different ghost layout
      bool must_update_ghost_values = true;

      // check whether the two vectors use the same parallel partitioner. if
      // not, check if all local ranges are the same (that way, we can
      // exchange data between different parallel layouts)
      if (partitioner.get() == 0)
        reinit (c, true);
      else if (partitioner.get() != c.partitioner.get())
        {
          size_type local_ranges_different_loc = (local_range() !=
                                                  c.local_range());
          if ((partitioner->n_mpi_processes() > 1 &&
               Utilities::MPI::max(local_ranges_different_loc,
                                   partitioner->get_communicator()) != 0)
              ||
              local_ranges_different_loc)
            reinit (c, true);
        }
      else
        must_update_ghost_values = vector_is_ghosted || c.vector_is_ghosted;

      vector_view = c.vector_view;
      if (must_update_ghost_values)
        update_ghost_values();
      return *this;
    }



    template <typename Number>
    template <typename Number2>
    inline
    Vector<Number> &
    Vector<Number>::operator = (const Vector<Number2> &c)
    {
      Assert (c.partitioner.get() != 0, ExcNotInitialized());

      // check whether the two vectors use the same parallel partitioner. if
      // not, check if all local ranges are the same (that way, we can
      // exchange data between different parallel layouts)
      if (partitioner.get() == 0)
        reinit (c, true);
      else if (partitioner.get() != c.partitioner.get())
        {
          size_type local_ranges_different_loc = (local_range() !=
                                                  c.local_range());
          if ((partitioner->n_mpi_processes() > 1 &&
               Utilities::MPI::max(local_ranges_different_loc,
                                   partitioner->get_communicator()) != 0)
              ||
              local_ranges_different_loc)
            reinit (c, true);
        }
      vector_view.reinit (partitioner->local_size(), val);

      if (partitioner->local_size())
        vector_view.equ (1., c.vector_view);

      if (vector_is_ghosted || c.vector_is_ghosted)
        update_ghost_values();
      return *this;
    }



    template <typename Number>
    inline
    void
    Vector<Number>::compress (::dealii::VectorOperation::values operation)
    {
      compress_start (0, operation);
      compress_finish(operation);
    }



    template <typename Number>
    inline
    void
    Vector<Number>::compress ()
    {
      compress(VectorOperation::unknown);
    }



    template <typename Number>
    inline
    void
    Vector<Number>::compress_finish (const bool add_value)
    {
      if (add_value)
        compress_finish(VectorOperation::add);
      else
        compress_finish(VectorOperation::insert);
    }



    template <typename Number>
    inline
    void
    Vector<Number>::update_ghost_values () const
    {
      update_ghost_values_start ();
      update_ghost_values_finish ();
    }



    template <typename Number>
    inline
    void
    Vector<Number>::zero_out_ghosts ()
    {
      std::fill_n (&val[partitioner->local_size()],
                   partitioner->n_ghost_indices(),
                   Number());
      vector_is_ghosted = false;
    }



    template <typename Number>
    inline
    bool
    Vector<Number>::has_ghost_elements () const
    {
      return vector_is_ghosted;
    }



    template <typename Number>
    inline
    bool
    Vector<Number>::all_zero_local () const
    {
      return partitioner->local_size()>0 ? vector_view.all_zero () : true;
    }



    template <typename Number>
    inline
    bool
    Vector<Number>::all_zero () const
    {
      // use int instead of bool. in order to make global reduction operations
      // work also when MPI_Init was not called, only call MPI_Allreduce
      // commands when there is more than one processor (note that reinit()
      // functions handle this case correctly through the job_supports_mpi()
      // query). this is the same in all the functions below
      int local_result = -static_cast<int>(all_zero_local());
      if (partitioner->n_mpi_processes() > 1)
        return -Utilities::MPI::max(local_result,
                                    partitioner->get_communicator());
      else
        return -local_result;
    }



    template <typename Number>
    inline
    bool
    Vector<Number>::is_non_negative_local () const
    {
      return partitioner->local_size()>0 ? vector_view.is_non_negative () : true;
    }



    template <typename Number>
    inline
    bool
    Vector<Number>::is_non_negative () const
    {
      int local_result = -static_cast<int>(is_non_negative_local());
      if (partitioner->n_mpi_processes() > 1)
        return -Utilities::MPI::max(local_result,
                                    partitioner->get_communicator());
      else
        return -local_result;
    }



    template <typename Number>
    template <typename Number2>
    inline
    bool
    Vector<Number>::vectors_equal_local (const Vector<Number2> &v) const
    {
      return partitioner->local_size()>0 ?
             vector_view.template operator == <Number2>(v.vector_view)
             : true;
    }



    template <typename Number>
    template <typename Number2>
    inline
    bool
    Vector<Number>::operator == (const Vector<Number2> &v) const
    {
      // MPI does not support bools, so use unsigned int instead. Two vectors
      // are equal if the check for non-equal fails on all processors
      unsigned int local_result = static_cast<int>(!vectors_equal_local(v));
      unsigned int result =
        partitioner->n_mpi_processes() > 1
        ?
        Utilities::MPI::max(local_result, partitioner->get_communicator())
        :
        local_result;
      return result==0;
    }



    template <typename Number>
    template <typename Number2>
    inline
    bool
    Vector<Number>::operator != (const Vector<Number2> &v) const
    {
      return !(operator == (v));
    }



    template <typename Number>
    template <typename Number2>
    inline
    Number
    Vector<Number>::inner_product_local(const Vector<Number2> &V) const
    {
      // on some processors, the size might be zero, which is not allowed by
      // the dealii::Vector class. Therefore, insert a check here
      return (partitioner->local_size()>0 ?
              vector_view.operator* (V.vector_view)
              : Number());
    }



    template <typename Number>
    template <typename Number2>
    inline
    Number
    Vector<Number>::operator * (const Vector<Number2> &V) const
    {
      Number local_result = inner_product_local(V);
      if (partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum (local_result,
                                    partitioner->get_communicator());
      else
        return local_result;
    }



    template <typename Number>
    inline
    typename Vector<Number>::real_type
    Vector<Number>::norm_sqr_local () const
    {
      return partitioner->local_size()>0 ? vector_view.norm_sqr() : real_type();
    }



    template <typename Number>
    inline
    typename Vector<Number>::real_type
    Vector<Number>::norm_sqr () const
    {
      real_type local_result = norm_sqr_local();
      if (partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum(local_result,
                                   partitioner->get_communicator());
      else
        return local_result;
    }



    template <typename Number>
    inline
    Number
    Vector<Number>::mean_value_local () const
    {
      Assert (partitioner->size()!=0, ExcEmptyObject());
      return (partitioner->local_size() ?
              vector_view.mean_value()
              : Number());
    }



    template <typename Number>
    inline
    Number
    Vector<Number>::mean_value () const
    {
      Number local_result = mean_value_local();
      if (partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum (local_result *
                                    (real_type)partitioner->local_size(),
                                    partitioner->get_communicator())
               /(real_type)partitioner->size();
      else
        return local_result;
    }



    template <typename Number>
    inline
    typename Vector<Number>::real_type
    Vector<Number>::l1_norm_local () const
    {
      return partitioner->local_size() ? vector_view.l1_norm() : real_type();
    }



    template <typename Number>
    inline
    typename Vector<Number>::real_type
    Vector<Number>::l1_norm () const
    {
      real_type local_result = l1_norm_local();
      if (partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::sum(local_result,
                                   partitioner->get_communicator());
      else
        return local_result;
    }



    template <typename Number>
    inline
    typename Vector<Number>::real_type
    Vector<Number>::l2_norm () const
    {
      return std::sqrt(norm_sqr());
    }



    template <typename Number>
    inline
    typename Vector<Number>::real_type
    Vector<Number>::lp_norm_local (const real_type p) const
    {
      return partitioner->local_size() ? vector_view.lp_norm(p) : real_type();
    }



    template <typename Number>
    inline
    typename Vector<Number>::real_type
    Vector<Number>::lp_norm (const real_type p) const
    {
      const real_type local_result = lp_norm_local(p);
      if (partitioner->n_mpi_processes() > 1)
        return std::pow (Utilities::MPI::sum(std::pow(local_result,p),
                                             partitioner->get_communicator()),
                         static_cast<real_type>(1.0/p));
      else
        return local_result;
    }



    template <typename Number>
    inline
    typename Vector<Number>::real_type
    Vector<Number>::linfty_norm_local () const
    {
      return partitioner->local_size() ? vector_view.linfty_norm() : real_type();
    }



    template <typename Number>
    inline
    typename Vector<Number>::real_type
    Vector<Number>::linfty_norm () const
    {
      const real_type local_result = linfty_norm_local();
      if (partitioner->n_mpi_processes() > 1)
        return Utilities::MPI::max (local_result,
                                    partitioner->get_communicator());
      else
        return local_result;
    }



    template <typename Number>
    inline
    typename Vector<Number>::size_type
    Vector<Number>::size () const
    {
      return partitioner->size();
    }



    template <typename Number>
    inline
    typename Vector<Number>::size_type
    Vector<Number>::local_size () const
    {
      return partitioner->local_size();
    }



    template <typename Number>
    inline
    std::pair<typename Vector<Number>::size_type,
        typename Vector<Number>::size_type>
        Vector<Number>::local_range () const
    {
      return partitioner->local_range();
    }



    template <typename Number>
    inline
    bool
    Vector<Number>::in_local_range
    (const size_type global_index) const
    {
      return partitioner->in_local_range (global_index);
    }



    template <typename Number>
    inline
    IndexSet
    Vector<Number>::locally_owned_elements() const
    {
      IndexSet is (size());

      is.add_range (local_range().first, local_range().second);

      return is;
    }



    template <typename Number>
    inline
    typename Vector<Number>::size_type
    Vector<Number>::n_ghost_entries () const
    {
      return partitioner->n_ghost_indices();
    }



    template <typename Number>
    inline
    const IndexSet &
    Vector<Number>::ghost_elements() const
    {
      return partitioner->ghost_indices();
    }



    template <typename Number>
    inline
    bool
    Vector<Number>::is_ghost_entry (const size_type global_index) const
    {
      return partitioner->is_ghost_entry (global_index);
    }



    template <typename Number>
    inline
    typename Vector<Number>::iterator
    Vector<Number>::begin ()
    {
      return vector_view.begin();
    }



    template <typename Number>
    inline
    typename Vector<Number>::const_iterator
    Vector<Number>::begin () const
    {
      return vector_view.begin();
    }



    template <typename Number>
    inline
    typename Vector<Number>::iterator
    Vector<Number>::end ()
    {
      return vector_view.end();
    }



    template <typename Number>
    inline
    typename Vector<Number>::const_iterator
    Vector<Number>::end () const
    {
      return vector_view.end();
    }



    template <typename Number>
    inline
    Number
    Vector<Number>::operator() (const size_type global_index) const
    {
      // do not allow reading a vector which is not in ghost mode
      Assert (in_local_range (global_index) || vector_is_ghosted == true,
              ExcMessage("You tried to read a ghost element of this vector, "
                         "but it has not imported its ghost values."));
      return val[partitioner->global_to_local(global_index)];
    }



    template <typename Number>
    inline
    Number &
    Vector<Number>::operator() (const size_type global_index)
    {
      // we would like to prevent reading ghosts from a vector that does not
      // have them imported, but this is not possible because we might be in a
      // part of the code where the vector has enabled ghosts but is non-const
      // (then, the compiler picks this method according to the C++ rule book
      // even if a human would pick the const method when this subsequent use
      // is just a read)
      return val[partitioner->global_to_local (global_index)];
    }



    template <typename Number>
    inline
    Number
    Vector<Number>::operator[] (const size_type global_index) const
    {
      return operator()(global_index);
    }



    template <typename Number>
    inline
    Number &
    Vector<Number>::operator[] (const size_type global_index)
    {
      return operator()(global_index);
    }



    template <typename Number>
    template <typename OtherNumber>
    inline
    void Vector<Number>::extract_subvector_to (const std::vector<size_type> &indices,
                                               std::vector<OtherNumber> &values) const
    {
      for (size_type i = 0; i < indices.size(); ++i)
        values[i] = operator()(indices[i]);
    }



    template <typename Number>
    template <typename ForwardIterator, typename OutputIterator>
    inline
    void Vector<Number>::extract_subvector_to (ForwardIterator          indices_begin,
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



    template <typename Number>
    inline
    Number
    Vector<Number>::local_element (const size_type local_index) const
    {
      AssertIndexRange (local_index,
                        partitioner->local_size()+
                        partitioner->n_ghost_indices());
      // do not allow reading a vector which is not in ghost mode
      Assert (local_index < local_size() || vector_is_ghosted == true,
              ExcMessage("You tried to read a ghost element of this vector, "
                         "but it has not imported its ghost values."));
      return val[local_index];
    }



    template <typename Number>
    inline
    Number &
    Vector<Number>::local_element (const size_type local_index)
    {
      AssertIndexRange (local_index,
                        partitioner->local_size()+
                        partitioner->n_ghost_indices());
      return val[local_index];
    }



    template <typename Number>
    inline
    Vector<Number> &
    Vector<Number>::operator = (const Number s)
    {
      // if we call Vector::operator=0, we want to zero out all the entries
      // plus ghosts.
      if (partitioner->local_size() > 0)
        vector_view.dealii::template Vector<Number>::operator= (s);
      if (s==Number())
        zero_out_ghosts();

      return *this;
    }



    template <typename Number>
    inline
    Vector<Number> &
    Vector<Number>::operator += (const Vector<Number> &v)
    {
      AssertDimension (local_size(), v.local_size());

      // dealii::Vector does not allow empty fields but this might happen on
      // some processors for parallel implementation
      if (local_size()>0)
        vector_view += v.vector_view;

      if (vector_is_ghosted)
        update_ghost_values();

      return *this;
    }



    template <typename Number>
    inline
    Vector<Number> &
    Vector<Number>::operator -= (const Vector<Number> &v)
    {
      AssertDimension (local_size(), v.local_size());

      // dealii::Vector does not allow empty fields but this might happen on
      // some processors for parallel implementation
      if (local_size()>0)
        vector_view -= v.vector_view;

      if (vector_is_ghosted)
        update_ghost_values();

      return *this;
    }



    template <typename Number>
    template <typename OtherNumber>
    inline
    void
    Vector<Number>::add (const std::vector<size_type> &indices,
                         const std::vector<OtherNumber>  &values)
    {
      AssertDimension (indices.size(), values.size());
      add (indices.size(), &indices[0], &values[0]);
    }



    template <typename Number>
    template <typename OtherNumber>
    inline
    void
    Vector<Number>::add (const std::vector<size_type>    &indices,
                         const ::dealii::Vector<OtherNumber> &values)
    {
      AssertDimension (indices.size(), values.size());
      add (indices.size(), &indices[0], values.begin());
    }



    template <typename Number>
    template <typename OtherNumber>
    inline
    void
    Vector<Number>::add (const size_type    n_indices,
                         const size_type   *indices,
                         const OtherNumber *values)
    {
      for (size_type i=0; i<n_indices; ++i)
        {
          Assert (numbers::is_finite(values[i]),
                  ExcMessage("The given value is not finite but either infinite or Not A Number (NaN)"));
          this->operator()(indices[i]) += values[i];
        }
    }



    template <typename Number>
    inline
    void
    Vector<Number>::add (const Number a)
    {
      // dealii::Vector does not allow empty fields but this might happen on
      // some processors for parallel implementation
      if (local_size())
        vector_view.add (a);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    inline
    void
    Vector<Number>::add (const Vector<Number> &v)
    {
      // dealii::Vector does not allow empty fields but this might happen on
      // some processors for parallel implementation
      if (local_size())
        vector_view.add (v.vector_view);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    inline
    void
    Vector<Number>::add (const Number a,
                         const Vector<Number> &v)
    {
      // dealii::Vector does not allow empty fields but this might happen on
      // some processors for parallel implementation
      if (local_size())
        vector_view.add (a, v.vector_view);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    inline
    void
    Vector<Number>::add (const Number a,
                         const Vector<Number> &v,
                         const Number b,
                         const Vector<Number> &w)
    {
      // dealii::Vector does not allow empty fields but this might happen on
      // some processors for parallel implementation
      if (local_size())
        vector_view.add (a, v.vector_view, b, w.vector_view);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    inline
    void
    Vector<Number>::sadd (const Number x,
                          const Vector<Number> &v)
    {
      // dealii::Vector does not allow empty fields but this might happen on
      // some processors for parallel implementation
      if (local_size())
        vector_view.sadd (x, v.vector_view);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    inline
    void
    Vector<Number>::sadd (const Number x,
                          const Number a,
                          const Vector<Number> &v)
    {
      // dealii::Vector does not allow empty fields but this might happen on
      // some processors for parallel implementation
      if (local_size())
        vector_view.sadd (x, a, v.vector_view);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    inline
    void
    Vector<Number>::sadd (const Number x,
                          const Number a,
                          const Vector<Number> &v,
                          const Number b,
                          const Vector<Number> &w)
    {
      // dealii::Vector does not allow empty fields but this might happen on
      // some processors for parallel implementation
      if (local_size())
        vector_view.sadd (x, a, v.vector_view, b, w.vector_view);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    inline
    void
    Vector<Number>::sadd (const Number s,
                          const Number a,
                          const Vector<Number> &v,
                          const Number b,
                          const Vector<Number> &w,
                          const Number c,
                          const Vector<Number> &x)
    {
      // dealii::Vector does not allow empty fields but this might happen on
      // some processors for parallel implementation
      if (local_size())
        vector_view.sadd (s, a, v.vector_view, b, w.vector_view,
                          c, x.vector_view);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    inline
    void
    Vector<Number>::scale (const Number factor)
    {
      operator *=(factor);
    }



    template <typename Number>
    inline
    Vector<Number> &
    Vector<Number>::operator *= (const Number factor)
    {
      // dealii::Vector does not allow empty fields but this might happen on
      // some processors for parallel implementation
      if (local_size())
        vector_view *= factor;

      if (vector_is_ghosted)
        update_ghost_values();

      return *this;
    }



    template <typename Number>
    inline
    Vector<Number> &
    Vector<Number>::operator /= (const Number factor)
    {
      operator *= (1./factor);
      return *this;
    }



    template <typename Number>
    inline
    void
    Vector<Number>::scale (const Vector<Number> &scaling_factors)
    {
      // dealii::Vector does not allow empty fields but this might happen on
      // some processors for parallel implementation
      if (local_size())
        vector_view.scale (scaling_factors.vector_view);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    template <typename Number2>
    inline
    void
    Vector<Number>::scale (const Vector<Number2> &scaling_factors)
    {
      if (local_size())
        vector_view.template scale<Number2> (scaling_factors.vector_view);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    inline
    void
    Vector<Number>::equ (const Number a,
                         const Vector<Number> &v)
    {
      // dealii::Vector does not allow empty fields but this might happen on
      // some processors for parallel implementation
      if (local_size())
        vector_view.equ (a, v.vector_view);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    template <typename Number2>
    inline
    void
    Vector<Number>::equ (const Number a,
                         const Vector<Number2> &v)
    {
      // dealii::Vector does not allow empty fields but this might happen on
      // some processors for parallel implementation
      if (local_size())
        vector_view.equ (a, v.vector_view);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    inline
    void
    Vector<Number>::equ (const Number a,
                         const Vector<Number> &v,
                         const Number b,
                         const Vector<Number> &w)
    {
      // dealii::Vector does not allow empty fields but this might happen on
      // some processors for parallel implementation
      if (local_size())
        vector_view.equ (a, v.vector_view, b, w.vector_view);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    inline
    void
    Vector<Number>::equ (const Number a,
                         const Vector<Number> &v,
                         const Number b,
                         const Vector<Number> &w,
                         const Number c,
                         const Vector<Number> &x)
    {
      // dealii::Vector does not allow empty fields but this might happen on
      // some processors for parallel implementation
      if (local_size())
        vector_view.equ (a, v.vector_view, b, w.vector_view,
                         c, x.vector_view);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    inline
    void
    Vector<Number>::ratio (const Vector<Number> &a,
                           const Vector<Number> &b)
    {
      // dealii::Vector does not allow empty fields but this might happen on
      // some processors for parallel implementation
      if (local_size())
        vector_view.ratio (a.vector_view, b.vector_view);

      if (vector_is_ghosted)
        update_ghost_values();
    }



    template <typename Number>
    inline
    const MPI_Comm &
    Vector<Number>::get_mpi_communicator() const
    {
      return partitioner->get_communicator();
    }



    template <typename Number>
    inline
    bool
    Vector<Number>::partitioners_are_compatible
    (const Utilities::MPI::Partitioner &part) const
    {
      return partitioner->is_compatible (part);
    }

#endif  // ifndef DOXYGEN

  } // end of namespace distributed

} // end of namespace parallel



/**
 * Global function @p swap which overloads the default implementation
 * of the C++ standard library which uses a temporary object. The
 * function simply exchanges the data of the two vectors.
 *
 * @relates Vector
 * @author Katharina Kormann, Martin Kronbichler, 2011
 */
template <typename Number>
inline
void swap (parallel::distributed::Vector<Number> &u,
           parallel::distributed::Vector<Number> &v)
{
  u.swap (v);
}


DEAL_II_NAMESPACE_CLOSE

#endif

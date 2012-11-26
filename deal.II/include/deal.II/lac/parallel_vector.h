//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2011, 2012 by deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
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


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace distributed
  {

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
     *   contiguously, [local_size(), local_size()+n_ghost_indices()).
     *
     * Functions related to parallel functionality:
     * - The function <code>compress()</code> goes through the data associated
     *   with ghost indices and communicates it to the owner process, which can
     *   then add/set it to the correct position. This can be used e.g. after
     *   having run an assembly routine involving ghosts that fill this vector.
     * - The <code>update_ghost_values()</code> function imports the data from the owning
     *   processor to the ghost indices in order to provide read access to the
     *   data associated with ghosts.
     * - It is possible to split the above functions into two phases, where the first
     *   initiates the communication and the second one finishes it. These
     *   functions can be used to overlap communication with computations in
     *   other parts of the code.
     * - Of course, reduction operations (like norms) make use of collective
     *   all-to-all MPI communications.
     *
     * @author Katharina Kormann, Martin Kronbichler, 2010, 2011
     */
    template <typename Number>
    class Vector : public Subscriptor
    {
    public:
      /**
       * Declare standard types used in all
       * containers. These types parallel those in
       * the <tt>C++</tt> standard libraries
       * <tt>vector<...></tt> class.
       */
      typedef Number                                            value_type;
      typedef value_type                                       *pointer;
      typedef const value_type                                 *const_pointer;
      typedef value_type                                       *iterator;
      typedef const value_type                                 *const_iterator;
      typedef value_type                                       &reference;
      typedef const value_type                                 &const_reference;
      typedef size_t                                            size_type;
      typedef typename numbers::NumberTraits<Number>::real_type real_type;

      /**
       * @name 1: Basic Object-handling
       */
      //@{
      /**
       * Empty constructor.
       */
      Vector ();

      /**
       * Copy constructor. Uses the parallel
       * partitioning of @p in_vector.
       */
      Vector (const Vector<Number> &in_vector);

      /**
       * Constructs a parallel vector of the given
       * global size without any actual parallel
       * distribution.
       */
      Vector (const unsigned int size);

      /**
       * Constructs a parallel vector. The local
       * range is specified by @p locally_owned_set
       * (note that this must be a contiguous
       * interval, multiple intervals are not
       * possible). The IndexSet @p ghost_indices
       * specifies ghost indices, i.e., indices
       * which one might need to read data from or
       * accumulate data from. It is allowed that
       * the set of ghost indices also contains the
       * local range, but it does not need to.
       *
       * This function involves global
       * communication, so it should only be called
       * once for a given layout. Use the
       * constructor with Vector<Number> argument to
       * create additional vectors with the same
       * parallel layout.
       */
      Vector (const IndexSet &local_range,
              const IndexSet &ghost_indices,
              const MPI_Comm  communicator);

      /**
       * Create the vector based on the parallel
       * partitioning described in @p
       * partitioner. The input argument is a shared
       * pointer, which store the partitioner data
       * only once and share it between several
       * vectors with the same layout.
       */
      Vector (const std_cxx1x::shared_ptr<const Utilities::MPI::Partitioner> &partitioner);

      /**
       * Destructor.
       */
      ~Vector ();

      /**
       * Sets the global size of the vector to @p
       * size without any actual parallel
       * distribution.
       */
      void reinit (const unsigned int size,
                   const bool         fast = false);

      /**
       * Uses the parallel layout of the input
       * vector @p in_vector and allocates memory
       * for this vector. Recommended initialization
       * function when several vectors with the same
       * layout should be created.
       *
       * If the flag @p fast is set to false, the
       * memory will be initialized with zero,
       * otherwise the memory will be untouched (and
       * the user must make sure to fill it with
       * reasonable data before using it).
       */
      template <typename Number2>
      void reinit(const Vector<Number2> &in_vector,
                  const bool             fast = false);

      /**
       * Initialize the vector. The local range is
       * specified by @p locally_owned_set (note
       * that this must be a contiguous interval,
       * multiple intervals are not possible). The
       * IndexSet @p ghost_indices specifies ghost
       * indices, i.e., indices which one might need
       * to read data from or accumulate data
       * from. It is allowed that the set of ghost
       * indices also contains the local range, but
       * it does not need to.
       *
       * This function involves global
       * communication, so it should only be called
       * once for a given layout. Use the @p reinit
       * function with Vector<Number> argument to
       * create additional vectors with the same
       * parallel layout.
       */
      void reinit (const IndexSet &local_range,
                   const IndexSet &ghost_indices,
                   const MPI_Comm  communicator);

      /**
       * Initialize the vector given to the parallel
       * partitioning described in @p
       * partitioner. The input argument is a shared
       * pointer, which store the partitioner data
       * only once and share it between several
       * vectors with the same layout.
       */
      void reinit (const std_cxx1x::shared_ptr<const Utilities::MPI::Partitioner> &partitioner);

      /**
       * Swap the contents of this
       * vector and the other vector
       * @p v. One could do this
       * operation with a temporary
       * variable and copying over the
       * data elements, but this
       * function is significantly more
       * efficient since it only swaps
       * the pointers to the data of
       * the two vectors and therefore
       * does not need to allocate
       * temporary storage and move
       * data around.
       *
       * This function is analog to the
       * the @p swap function of all C++
       * standard containers. Also,
       * there is a global function
       * <tt>swap(u,v)</tt> that simply calls
       * <tt>u.swap(v)</tt>, again in analogy
       * to standard functions.
       *
       * This function is virtual in
       * order to allow for derived
       * classes to handle memory
       * separately.
       */
      void swap (Vector<Number> &v);

      /**
       * Assigns the vector to the parallel
       * partitioning of the input vector @p
       * in_vector, and copies all the data.
       */
      Vector<Number> &
      operator = (const Vector<Number>  &in_vector);

      /**
       * Assigns the vector to the parallel
       * partitioning of the input vector @p
       * in_vector, and copies all the data.
       */
      template <typename Number2>
      Vector<Number> &
      operator = (const Vector<Number2> &in_vector);

      /**
       * This method copies the local range from
       * another vector with the same local range,
       * but possibly different layout of ghost
       * indices.
       */
      void copy_from (const Vector<Number> &in_vector,
                      const bool            call_update_ghost_values = false);

      /**
       * Sets all elements of the vector to the
       * scalar @p s. If the scalar is zero, also
       * ghost elements are set to zero, otherwise
       * they remain unchanged.
       */
      Vector<Number> &operator = (const Number s);

      /**
       * This function copies the data that has
       * accumulated in the data buffer for ghost
       * indices to the owning processor.
       *
       * For the meaning of this argument,
       * see the entry on @ref
       * GlossCompress "Compressing
       * distributed vectors and matrices"
       * in the glossary.
       */
      void compress (::dealii::VectorOperation::values operation
                     =::dealii::VectorOperation::unknown);


      /**
       * Fills the data field for ghost indices with
       * the values stored in the respective
       * positions of the owning processor. This
       * function is needed before reading from
       * ghosts. The function is @p const even
       * though ghost data is changed. This is
       * needed to allow functions with a @p const
       * vector to perform the data exchange without
       * creating temporaries.
       */
      void update_ghost_values () const;

      /**
       * Initiates communication for the @p
       * compress() function with non-blocking
       * communication. This function does not wait
       * for the transfer to finish, in order to
       * allow for other computations during the
       * time it takes until all data arrives.
       *
       * Before the data is actually exchanged, the
       * function must be followed by a call to @p
       * compress_finish().
       *
       * In case this function is called for more
       * than one vector before @p
       * compress_finish() is invoked, it is
       * mandatory to specify a unique
       * communication channel to each such call, in
       * order to avoid several messages with the
       * same ID that will corrupt this operation.
       */
      void compress_start (const unsigned int communication_channel = 0);

      /**
       * For all requests that have been initiated
       * in compress_start, wait for the
       * communication to finish. Once it is
       * finished, add or set the data (depending on
       * whether @p add_ghost_data is @p true or @p
       * false) to the respective positions in the
       * owning processor, and clear the contents in
       * the ghost data fields. The meaning of
       * this argument is the same as in compress().
       *
       * Must follow a call to the @p compress_start
       * function.
       */
      void compress_finish (const bool add_ghost_data = true);


      /**
       * Initiates communication for the @p
       * update_ghost_values() function with non-blocking
       * communication. This function does not wait
       * for the transfer to finish, in order to
       * allow for other computations during the
       * time it takes until all data arrives.
       *
       * Before the data is actually exchanged, the
       * function must be followed by a call to @p
       * update_ghost_values_finish().
       *
       * In case this function is called for more
       * than one vector before @p
       * update_ghost_values_finish() is invoked, it is
       * mandatory to specify a unique communication
       * channel to each such call, in order to
       * avoid several messages with the same ID
       * that will corrupt this operation.
       */
      void update_ghost_values_start (const unsigned int communication_channel = 0) const;


      /**
       * For all requests that have been started in
       * update_ghost_values_start, wait for the communication
       * to finish.
       *
       * Must follow a call to the @p
       * update_ghost_values_start function before reading
       * data from ghost indices.
       */
      void update_ghost_values_finish () const;

      /**
       * This method zeros the entries on ghost
       * dofs, but does not touch locally owned
       * DoFs.
       */
      void zero_out_ghosts ();

      /**
       * Return whether the vector contains only
       * elements with value zero. This function
       * is mainly for internal consistency
       * checks and should seldom be used when
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
       *
       * The function obviously only makes
       * sense if the template argument of this
       * class is a real type. If it is a
       * complex type, then an exception is
       * thrown.
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
       * Computes the square of the l<sub>2</sub>
       * norm of the vector (i.e., the sum of the
       * squares of all entries among all
       * processors).
       */
      real_type norm_sqr () const;

      /**
       * Computes the mean value of all the entries
       * in the vector.
       */
      Number mean_value () const;

      /**
       * Returns the l<sub>1</sub> norm of the
       * vector (i.e., the sum of the absolute
       * values of all entries among all
       * processors).
       */
      real_type l1_norm () const;

      /**
       * Returns the l<sub>2</sub> norm of the
       * vector (i.e., square root of the sum of the
       * square of all entries among all
       * processors).
       */
      real_type l2_norm () const;

      /**
       * Returns the l<sub>p</sub> norm with real @p
       * p of the vector (i.e., the pth root of sum
       * of the pth power of all entries among all
       * processors).
       */
      real_type lp_norm (const real_type p) const;

      /**
       * Returns the maximum norm of the vector
       * (i.e., maximum absolute value among all
       * entries among all processors).
       */
      real_type linfty_norm () const;

      /**
       * Returns the global size of the vector,
       * equal to the sum of the number of locally
       * owned indices among all the processors.
       */
      types::global_dof_index size () const;

      /**
       * Returns the local size of the vector, i.e.,
       * the number of indices owned locally.
       */
      unsigned int local_size() const;

      /**
       * Returns the half-open interval that
       * specifies the locally owned range of the
       * vector. Note that <code>local_size() ==
       * local_range().second -
       * local_range().first</code>.
       */
      std::pair<types::global_dof_index, types::global_dof_index> local_range () const;

      /**
       * Returns true if the given global index is
       * in the local range of this processor.
       */
      bool in_local_range (const types::global_dof_index global_index) const;

      /**
       * Returns the number of ghost elements
       * present on the vector.
       */
      unsigned int n_ghost_entries () const;

      /**
       * Returns whether the given global index is a
       * ghost index on the present
       * processor. Returns false for indices that
       * are owned locally and for indices not
       * present at all.
       */
      bool is_ghost_entry (const types::global_dof_index global_index) const;

      /**
       * Make the @p Vector class a bit like
       * the <tt>vector<></tt> class of the C++
       * standard library by returning
       * iterators to the start and end of the
       * locally owned elements of this vector.
       */
      iterator begin ();

      /**
       * Return constant iterator to the start of
       * the vector.
       */
      const_iterator begin () const;

      /**
       * Return an iterator pointing to the
       * element past the end of the array of
       * locally owned entries.
       */
      iterator end ();

      /**
       * Return a constant iterator pointing to
       * the element past the end of the array
       * of the locally owned entries.
       */
      const_iterator end () const;
      //@}


      /**
       * @name 2: Data-Access
       */
      //@{

      /**
       * Read access to the data in the
       * position corresponding to @p
       * global_index. The index must be
       * either in the local range of the
       * vector or be specified as a ghost
       * index at construction.
       */
      Number operator () (const types::global_dof_index global_index) const;

      /**
       * Read and write access to the data
       * in the position corresponding to
       * @p global_index. The index must be
       * either in the local range of the
       * vector or be specified as a ghost
       * index at construction.
       */
      Number &operator () (const types::global_dof_index global_index);

      /**
       * Read access to the data in the
       * position corresponding to @p
       * global_index. The index must be
       * either in the local range of the
       * vector or be specified as a ghost
       * index at construction.
       *
       * This function does the same thing
       * as operator().
       */
      Number operator [] (const types::global_dof_index global_index) const;

      /**
       * Read and write access to the data
       * in the position corresponding to
       * @p global_index. The index must be
       * either in the local range of the
       * vector or be specified as a ghost
       * index at construction.
       *
       * This function does the same thing
       * as operator().
       */
      Number &operator [] (const types::global_dof_index global_index);

      /**
       * Read access to the data field specified by
       * @p local_index. Locally owned indices can
       * be accessed with indices
       * <code>[0,local_size)</code>, and ghost
       * indices with indices
       * <code>[local_size,local_size+
       * n_ghost_entries]</code>.
       */
      Number local_element (const unsigned int local_index) const;

      /**
       * Read and write access to the data field
       * specified by @p local_index. Locally owned
       * indices can be accessed with indices
       * <code>[0,local_size)</code>, and ghost
       * indices with indices
       * <code>[local_size,local_size+n_ghosts]</code>.
       */
      Number &local_element (const unsigned int local_index);
      //@}


      /**
       * @name 3: Modification of vectors
       */
      //@{

      /**
       * Add the given vector to the present
       * one.
       */
      Vector<Number> &operator += (const Vector<Number> &V);

      /**
       * Subtract the given vector from the
       * present one.
       */
      Vector<Number> &operator -= (const Vector<Number> &V);

      /**
       * A collective add operation:
       * This funnction adds a whole
       * set of values stored in @p
       * values to the vector
       * components specified by @p
       * indices.
       */
      template <typename OtherNumber>
      void add (const std::vector<unsigned int> &indices,
                const std::vector<OtherNumber>  &values);

      /**
       * This is a second collective
       * add operation. As a
       * difference, this function
       * takes a deal.II vector of
       * values.
       */
      template <typename OtherNumber>
      void add (const std::vector<unsigned int>     &indices,
                const ::dealii::Vector<OtherNumber> &values);

      /**
       * Take an address where
       * <tt>n_elements</tt> are stored
       * contiguously and add them into
       * the vector. Handles all cases
       * which are not covered by the
       * other two <tt>add()</tt>
       * functions above.
       */
      template <typename OtherNumber>
      void add (const unsigned int  n_elements,
                const unsigned int *indices,
                const OtherNumber  *values);

      /**
       * Addition of @p s to all
       * components. Note that @p s is a
       * scalar and not a vector.
       */
      void add (const Number s);

      /**
       * Simple vector addition, equal to the
       * <tt>operator +=</tt>.
       */
      void add (const Vector<Number> &V);

      /**
       * Simple addition of a multiple of a
       * vector, i.e. <tt>*this += a*V</tt>.
       */
      void add (const Number a, const Vector<Number> &V);

      /**
       * Multiple addition of scaled vectors,
       * i.e. <tt>*this += a*V+b*W</tt>.
       */
      void add (const Number a, const Vector<Number> &V,
                const Number b, const Vector<Number> &W);

      /**
       * Scaling and simple vector addition,
       * i.e.
       * <tt>*this = s*(*this)+V</tt>.
       */
      void sadd (const Number          s,
                 const Vector<Number> &V);

      /**
       * Scaling and simple addition, i.e.
       * <tt>*this = s*(*this)+a*V</tt>.
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
       * Scaling and multiple addition.
       * <tt>*this = s*(*this)+a*V + b*W + c*X</tt>.
       */
      void sadd (const Number          s,
                 const Number          a,
                 const Vector<Number> &V,
                 const Number          b,
                 const Vector<Number> &W,
                 const Number          c,
                 const Vector<Number> &X);

      /**
       * Scale each element of the
       * vector by the given factor.
       *
       * This function is deprecated
       * and will be removed in a
       * future version. Use
       * <tt>operator *=</tt> and
       * <tt>operator /=</tt> instead.
       */
      void scale (const Number factor) DEAL_II_DEPRECATED;


      /**
       * Scale each element of the
       * vector by a constant
       * value.
       */
      Vector<Number> &operator *= (const Number factor);

      /**
       * Scale each element of the
       * vector by the inverse of the
       * given value.
       */
      Vector<Number> &operator /= (const Number factor);

      /**
       * Scale each element of this
       * vector by the corresponding
       * element in the argument. This
       * function is mostly meant to
       * simulate multiplication (and
       * immediate re-assignment) by a
       * diagonal scaling matrix.
       */
      void scale (const Vector<Number> &scaling_factors);

      /**
       * Scale each element of this
       * vector by the corresponding
       * element in the argument. This
       * function is mostly meant to
       * simulate multiplication (and
       * immediate re-assignment) by a
       * diagonal scaling matrix.
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
       * Compute the elementwise ratio of the
       * two given vectors, that is let
       * <tt>this[i] = a[i]/b[i]</tt>. This is
       * useful for example if you want to
       * compute the cellwise ratio of true to
       * estimated error.
       *
       * This vector is appropriately
       * scaled to hold the result.
       *
       * If any of the <tt>b[i]</tt> is
       * zero, the result is
       * undefined. No attempt is made
       * to catch such situations.
       */
      void ratio (const Vector<Number> &a,
                  const Vector<Number> &b);
      //@}


      /**
       * @name 4: Mixed stuff
       */
      //@{
      /**
       * Checks whether the given
       * partitioner is compatible with the
       * partitioner used for this
       * vector. Two partitioners are
       * compatible if the have the same
       * local size and the same ghost
       * indices. They do not necessarily
       * need to be the same data
       * field. This is a local operation
       * only, i.e., if only some
       * processors decide that the
       * partitioning is not compatible,
       * only these processors will return
       * @p false, whereas the other
       * processors will return @p true.
       */
      bool
      partitioners_are_compatible (const Utilities::MPI::Partitioner &part) const;


      /**
       * Prints the vector to the output stream @p
       * out.
       */
      void print (std::ostream       &out,
                  const unsigned int  precision  = 3,
                  const bool          scientific = true,
                  const bool          across     = true) const;

      /**
       * Returns the memory consumption of this
       * class in bytes.
       */
      std::size_t memory_consumption () const;
      //@}

    private:
      /**
       * Shared pointer to store the parallel
       * partitioning information. This information
       * can be shared between several vectors that
       * have the same partitioning.
       */
      std_cxx1x::shared_ptr<const Utilities::MPI::Partitioner> partitioner;

      /**
       * The size that is currently allocated in the
       * val array.
       */
      unsigned int    allocated_size;

      /**
       * Pointer to the array of
       * local elements of this vector.
       */
      Number         *val;

      /**
       * Temporary storage that holds the data that
       * is sent to this processor in @p compress()
       * or sent from this processor in @p
       * update_ghost_values.
       */
      mutable Number *import_data;

      /**
       * Provide this class with all functionality
       * of ::dealii::Vector by creating a
       * VectorView object.
       */
      VectorView<Number> vector_view;

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
      /**
       * A vector that collects all requests from @p
       * compress() operations. This class uses
       * persistent MPI communicators, i.e., the
       * communication channels are stored during
       * successive calls to a given function. This
       * reduces the overhead involved with setting
       * up the MPI machinery, but it does not
       * remove the need for a receive operation to
       * be posted before the data can actually be
       * sent.
       */
      std::vector<MPI_Request>   compress_requests;

      /**
       * A vector that collects all requests from @p
       * update_ghost_values() operations. This class uses
       * persistent MPI communicators.
       */
      mutable std::vector<MPI_Request>   update_ghost_values_requests;
#endif

      /**
       * A lock that makes sure that
       * the @p compress and @p
       * update_ghost_values functions
       * give reasonable results also
       * when used with several
       * threads.
       */
      mutable Threads::ThreadMutex mutex;

      /**
       * A helper function that clears the
       * compress_requests and update_ghost_values_requests
       * field. Used in reinit functions.
       */
      void clear_mpi_requests ();

      /**
       * A helper function that is used to resize
       * the val array.
       */
      void resize_val (const unsigned int new_allocated_size);

      /*
       * Make all other vector types
       * friends.
       */
      template <typename Number2> friend class Vector;
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
      vector_view (0, static_cast<Number *>(0))
    {
      reinit (v, true);
      vector_view = v.vector_view;
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
      vector_view (0, static_cast<Number *>(0))
    {
      reinit (local_range, ghost_indices, communicator);
    }



    template <typename Number>
    inline
    Vector<Number>::Vector (const unsigned int size)
      :
      allocated_size (0),
      val (0),
      import_data (0),
      vector_view (0, static_cast<Number *>(0))
    {
      reinit (size, false);
    }



    template <typename Number>
    inline
    Vector<Number>::
    Vector (const std_cxx1x::shared_ptr<const Utilities::MPI::Partitioner> &partitioner)
      :
      allocated_size (0),
      val (0),
      import_data (0),
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

      // check whether the two vectors use the same
      // parallel partitioner. if not, check if all
      // local ranges are the same (that way, we can
      // exchange data between different parallel
      // layouts)
      if (partitioner.get() == 0)
        reinit (c, true);
      else if (partitioner.get() != c.partitioner.get())
        {
          unsigned int local_ranges_different_loc = (local_range() !=
                                                     c.local_range());
          if (Utilities::MPI::max(local_ranges_different_loc,
                                  partitioner->get_communicator()) != 0)
            reinit (c, true);
        }
      vector_view = c.vector_view;
      return *this;
    }



    template <typename Number>
    template <typename Number2>
    inline
    Vector<Number> &
    Vector<Number>::operator = (const Vector<Number2> &c)
    {
      Assert (c.partitioner.get() != 0, ExcNotInitialized());

      // check whether the two vectors use the same
      // parallel partitioner. if not, check if all
      // local ranges are the same (that way, we can
      // exchange data between different parallel
      // layouts)
      if (partitioner.get() == 0)
        reinit (c, true);
      else if (partitioner.get() != c.partitioner.get())
        {
          unsigned int local_ranges_different_loc = (local_range() !=
                                                     c.local_range());
          if (Utilities::MPI::max(local_ranges_different_loc,
                                  partitioner->get_communicator()) != 0)
            reinit (c, true);
        }
      vector_view.reinit (partitioner->local_size(), val);
      if (partitioner->local_size() > 0)
        vector_view.equ (1., c.vector_view);
      return *this;
    }



    template <typename Number>
    inline
    void
    Vector<Number>::compress (::dealii::VectorOperation::values operation)
    {
      compress_start ();
      if (operation == ::dealii::VectorOperation::insert)
        compress_finish (false);
      else
        compress_finish (true);
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
    }



    template <typename Number>
    inline
    bool
    Vector<Number>::all_zero () const
    {
      // use int instead of bool
      int local_result = (partitioner->local_size()>0 ?
                          -vector_view.all_zero () : -1);
      return - Utilities::MPI::max(local_result,partitioner->get_communicator());
    }



    template <typename Number>
    inline
    bool
    Vector<Number>::is_non_negative () const
    {
      // use int instead of bool
      int local_result = (partitioner->local_size()>0 ?
                          -vector_view.is_non_negative () : -1);
      return  - Utilities::MPI::max(local_result,partitioner->get_communicator());
    }



    template <typename Number>
    template <typename Number2>
    inline
    bool
    Vector<Number>::operator == (const Vector<Number2> &v) const
    {
      AssertDimension (local_size(), v.local_size());

      // MPI does not support bools, so use unsigned
      // int instead. Two vectors are equal if the
      // check for non-equal fails on all processors
      unsigned int local_result = (partitioner->local_size()>0 ?
                                   vector_view.template operator !=
                                   <Number2>(v.vector_view)
                                   : 0 );
      unsigned int result = Utilities::MPI::max(local_result,
                                                partitioner->get_communicator());
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
    Vector<Number>::operator * (const Vector<Number2> &V) const
    {
      Number local_result = (partitioner->local_size()>0 ?
                             vector_view.operator* (V.vector_view)
                             : 0);
      return Utilities::MPI::sum (local_result, partitioner->get_communicator());
    }



    template <typename Number>
    inline
    typename Vector<Number>::real_type
    Vector<Number>::norm_sqr () const
    {
      // on some processors, the size might be zero,
      // which is not allowed by the base
      // class. Therefore, insert a check here
      Number local_result = (partitioner->local_size()>0 ?
                             vector_view.norm_sqr()
                             : 0);
      return Utilities::MPI::sum(local_result,partitioner->get_communicator());
    }



    template <typename Number>
    inline
    Number
    Vector<Number>::mean_value () const
    {
      Number local_result =
        (partitioner->local_size()>0 ?
         vector_view.mean_value()
         : 0)
        *((real_type)partitioner->local_size()/(real_type)partitioner->size());
      return Utilities::MPI::sum (local_result, partitioner->get_communicator());
    }



    template <typename Number>
    inline
    typename Vector<Number>::real_type
    Vector<Number>::l1_norm () const
    {
      Number local_result = (partitioner->local_size()>0 ?
                             vector_view.l1_norm()
                             : 0);
      return Utilities::MPI::sum(local_result, partitioner->get_communicator());
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
    Vector<Number>::lp_norm (const real_type p) const
    {
      const Number local_result = (partitioner->local_size()>0 ?
                                   std::pow(vector_view.lp_norm(p),p)
                                   : 0);
      return std::pow (Utilities::MPI::sum(local_result,
                                           partitioner->get_communicator()),
                       static_cast<Number>(1.0/p));
    }



    template <typename Number>
    inline
    typename Vector<Number>::real_type
    Vector<Number>::linfty_norm () const
    {
      const Number local_result = (partitioner->local_size()>0 ?
                                   vector_view.linfty_norm()
                                   : 0);
      return Utilities::MPI::max (local_result, partitioner->get_communicator());
    }



    template <typename Number>
    inline
    types::global_dof_index Vector<Number>::size () const
    {
      return partitioner->size();
    }



    template <typename Number>
    inline
    unsigned int Vector<Number>::local_size () const
    {
      return partitioner->local_size();
    }



    template <typename Number>
    inline
    std::pair<types::global_dof_index, types::global_dof_index>
    Vector<Number>::local_range () const
    {
      return partitioner->local_range();
    }



    template <typename Number>
    inline
    bool
    Vector<Number>::in_local_range
    (const types::global_dof_index global_index) const
    {
      return partitioner->in_local_range (global_index);
    }



    template <typename Number>
    inline
    unsigned int
    Vector<Number>::n_ghost_entries () const
    {
      return partitioner->n_ghost_indices();
    }



    template <typename Number>
    inline
    bool
    Vector<Number>::is_ghost_entry (const types::global_dof_index global_index) const
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
    Vector<Number>::operator() (const types::global_dof_index global_index) const
    {
      return val[partitioner->global_to_local(global_index)];
    }



    template <typename Number>
    inline
    Number &
    Vector<Number>::operator() (const types::global_dof_index global_index)
    {
      return val[partitioner->global_to_local (global_index)];
    }



    template <typename Number>
    inline
    Number
    Vector<Number>::operator[] (const types::global_dof_index global_index) const
    {
      return operator()(global_index);
    }



    template <typename Number>
    inline
    Number &
    Vector<Number>::operator[] (const types::global_dof_index global_index)
    {
      return operator()(global_index);
    }



    template <typename Number>
    inline
    Number
    Vector<Number>::local_element (const unsigned int local_index) const
    {
      AssertIndexRange (local_index,
                        partitioner->local_size()+
                        partitioner->n_ghost_indices());
      return val[local_index];
    }



    template <typename Number>
    inline
    Number &
    Vector<Number>::local_element (const unsigned int local_index)
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
      // if we call Vector::operator=0, we want to
      // zero out all the entries plus ghosts.
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
      // dealii::Vector does not allow empty fields
      // but this might happen on some processors
      // for parallel implementation
      if (local_size()>0)
        vector_view += v.vector_view;
      return *this;
    }



    template <typename Number>
    inline
    Vector<Number> &
    Vector<Number>::operator -= (const Vector<Number> &v)
    {
      AssertDimension (local_size(), v.local_size());
      // dealii::Vector does not allow empty fields
      // but this might happen on some processors
      // for parallel implementation
      if (local_size()>0)
        vector_view -= v.vector_view;
      return *this;
    }



    template <typename Number>
    template <typename OtherNumber>
    inline
    void
    Vector<Number>::add (const std::vector<unsigned int> &indices,
                         const std::vector<OtherNumber>  &values)
    {
      AssertDimension (indices.size(), values.size());
      add (indices.size(), &indices[0], &values[0]);
    }



    template <typename Number>
    template <typename OtherNumber>
    inline
    void
    Vector<Number>::add (const std::vector<unsigned int>    &indices,
                         const ::dealii::Vector<OtherNumber> &values)
    {
      AssertDimension (indices.size(), values.size());
      add (indices.size(), &indices[0], values.begin());
    }



    template <typename Number>
    template <typename OtherNumber>
    inline
    void
    Vector<Number>::add (const unsigned int  n_indices,
                         const unsigned int *indices,
                         const OtherNumber  *values)
    {
      for (unsigned int i=0; i<n_indices; ++i)
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
      // dealii::Vector does not allow empty fields
      // but this might happen on some processors
      // for parallel implementation
      if (local_size()>0)
        vector_view.add (a);
    }



    template <typename Number>
    inline
    void
    Vector<Number>::add (const Vector<Number> &v)
    {
      AssertDimension (local_size(), v.local_size());
      // dealii::Vector does not allow empty fields
      // but this might happen on some processors
      // for parallel implementation
      if (local_size()>0)
        vector_view.add (v.vector_view);
    }



    template <typename Number>
    inline
    void
    Vector<Number>::add (const Number a,
                         const Vector<Number> &v)
    {
      AssertDimension (local_size(), v.local_size());
      // dealii::Vector does not allow empty fields
      // but this might happen on some processors
      // for parallel implementation
      if (local_size()>0)
        vector_view.add (a, v.vector_view);
    }



    template <typename Number>
    inline
    void
    Vector<Number>::add (const Number a,
                         const Vector<Number> &v,
                         const Number b,
                         const Vector<Number> &w)
    {
      AssertDimension (local_size(), v.local_size());
      AssertDimension (local_size(), w.local_size());
      // dealii::Vector does not allow empty fields
      // but this might happen on some processors
      // for parallel implementation
      if (local_size()>0)
        vector_view.add (a, v.vector_view, b, w.vector_view);
    }



    template <typename Number>
    inline
    void
    Vector<Number>::sadd (const Number x,
                          const Vector<Number> &v)
    {
      AssertDimension (local_size(), v.local_size());
      // dealii::Vector does not allow empty fields
      // but this might happen on some processors
      // for parallel implementation
      if (local_size()>0)
        vector_view.sadd (x, v.vector_view);
    }



    template <typename Number>
    inline
    void
    Vector<Number>::sadd (const Number x,
                          const Number a,
                          const Vector<Number> &v)
    {
      AssertDimension (local_size(), v.local_size());
      // dealii::Vector does not allow empty fields
      // but this might happen on some processors
      // for parallel implementation
      if (local_size()>0)
        vector_view.sadd (x, a, v.vector_view);
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
      AssertDimension (local_size(), v.local_size());
      AssertDimension (local_size(), w.local_size());
      // dealii::Vector does not allow empty fields
      // but this might happen on some processors
      // for parallel implementation
      if (local_size()>0)
        vector_view.sadd (x, a, v.vector_view, b, w.vector_view);
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
      AssertDimension (local_size(), v.local_size());
      AssertDimension (local_size(), w.local_size());
      AssertDimension (local_size(), x.local_size());
      // dealii::Vector does not allow empty fields
      // but this might happen on some processors
      // for parallel implementation
      if (local_size()>0)
        vector_view.sadd (s, a, v.vector_view, b, w.vector_view,
                          c, x.vector_view);
    }



    template <typename Number>
    inline
    void
    Vector<Number>::scale (const Number factor)
    {
      // dealii::Vector does not allow empty fields
      // but this might happen on some processors
      // for parallel implementation
      if (local_size()>0)
        vector_view.scale (factor);
    }



    template <typename Number>
    inline
    Vector<Number> &
    Vector<Number>::operator *= (const Number factor)
    {
      // dealii::Vector does not allow empty fields
      // but this might happen on some processors
      // for parallel implementation
      if (local_size()>0)
        vector_view.operator *= (factor);
      return *this;
    }



    template <typename Number>
    inline
    Vector<Number> &
    Vector<Number>::operator /= (const Number factor)
    {
      // dealii::Vector does not allow empty fields
      // but this might happen on some processors
      // for parallel implementation
      if (local_size()>0)
        vector_view.operator /= (factor);
      return *this;
    }



    template <typename Number>
    inline
    void
    Vector<Number>::scale (const Vector<Number> &scaling_factors)
    {
      // dealii::Vector does not allow empty fields
      // but this might happen on some processors
      // for parallel implementation
      if (local_size()>0)
        vector_view.scale (scaling_factors.vector_view);
    }



    template <typename Number>
    template <typename Number2>
    inline
    void
    Vector<Number>::scale (const Vector<Number2> &scaling_factors)
    {
      vector_view.template scale<Number2> (scaling_factors.vector_view);
    }



    template <typename Number>
    inline
    void
    Vector<Number>::equ (const Number a,
                         const Vector<Number> &v)
    {
      AssertDimension (local_size(), v.local_size());
      // dealii::Vector does not allow empty fields
      // but this might happen on some processors
      // for parallel implementation
      if (local_size()>0)
        vector_view.equ (a, v.vector_view);
    }



    template <typename Number>
    template <typename Number2>
    inline
    void
    Vector<Number>::equ (const Number a,
                         const Vector<Number2> &v)
    {
      AssertDimension (local_size(), v.local_size());
      // dealii::Vector does not allow empty fields
      // but this might happen on some processors
      // for parallel implementation
      if (local_size()>0)
        vector_view.equ (a, v.vector_view);
    }



    template <typename Number>
    inline
    void
    Vector<Number>::equ (const Number a,
                         const Vector<Number> &v,
                         const Number b,
                         const Vector<Number> &w)
    {
      AssertDimension (local_size(), v.local_size());
      AssertDimension (local_size(), w.local_size());
      // dealii::Vector does not allow empty fields
      // but this might happen on some processors
      // for parallel implementation
      if (local_size()>0)
        vector_view.equ (a, v.vector_view, b, w.vector_view);
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
      AssertDimension (local_size(), v.local_size());
      AssertDimension (local_size(), w.local_size());
      AssertDimension (local_size(), w.local_size());
      // dealii::Vector does not allow empty fields
      // but this might happen on some processors
      // for parallel implementation
      if (local_size()>0)
        vector_view.equ (a, v.vector_view, b, w.vector_view,
                         c, x.vector_view);
    }



    template <typename Number>
    inline
    void
    Vector<Number>::ratio (const Vector<Number> &a,
                           const Vector<Number> &b)
    {
      AssertDimension (local_size(), a.local_size());
      AssertDimension (local_size(), b.local_size());
      // dealii::Vector does not allow empty fields
      // but this might happen on some processors
      // for parallel implementation
      if (local_size()>0)
        vector_view.ratio (a.vector_view, b.vector_view);
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

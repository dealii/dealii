//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__trilinos_vector_h
#define __deal2__trilinos_vector_h


#include <base/config.h>
#include <base/subscriptor.h>
#include <lac/exceptions.h>
#include <lac/vector.h>

#include <vector>
#include <utility>
#include <memory>

#ifdef DEAL_II_USE_TRILINOS

#define TrilinosScalar double
#  include "Epetra_ConfigDefs.h"
#  ifdef DEAL_II_COMPILER_SUPPORTS_MPI // only if MPI is installed
#    include "mpi.h"
#    include "Epetra_MpiComm.h"
#  else
#  include "Epetra_SerialComm.h"
#  endif
#  include "Epetra_FEVector.h"
#  include "Epetra_Map.h"

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
  class Vector;

				   /**
				    * @cond internal
				    */

/**
 * A namespace for internal implementation details of the TrilinosWrapper
 * members.
 * @ingroup TrilinosWrappers
 */
  namespace internal
  {
                                     /**
                                      * This class implements a wrapper for
				      * accessing the Trilinos vector
				      * in the same way as we access 
				      * deal.II objects:
                                      * it is initialized with a vector and an
                                      * element within it, and has a
                                      * conversion operator to extract the
                                      * scalar value of this element. It also
                                      * has a variety of assignment operator
                                      * for writing to this one element.
				      * @ingroup TrilinosWrappers
                                      */
    class VectorReference
    {
      private:
                                         /**
                                          * Constructor. It is made private so
                                          * as to only allow the actual vector
                                          * class to create it.
                                          */
        VectorReference (Vector  &vector,
                         const unsigned int index);

      public:
                                         /**
                                          * This looks like a copy operator,
                                          * but does something different than
                                          * usual. In particular, it does not
                                          * copy the member variables of this
                                          * reference. Rather, it handles the
                                          * situation where we have two
                                          * vectors @p v and @p w, and assign
                                          * elements like in
                                          * <tt>v(i)=w(i)</tt>. Here, both
                                          * left and right hand side of the
                                          * assignment have data type
                                          * VectorReference, but what we
                                          * really mean is to assign the
                                          * vector elements represented by the
                                          * two references. This operator
                                          * implements this operation. Note
                                          * also that this allows us to make
                                          * the assignment operator const.
                                          */
        const VectorReference & operator = (const VectorReference &r) const;

                                         /**
                                          * Set the referenced element of the
                                          * vector to <tt>s</tt>.
                                          */
        const VectorReference & operator = (const TrilinosScalar &s) const;

                                         /**
                                          * Add <tt>s</tt> to the referenced
                                          * element of the vector->
                                          */
        const VectorReference & operator += (const TrilinosScalar &s) const;

                                         /**
                                          * Subtract <tt>s</tt> from the
                                          * referenced element of the vector->
                                          */
        const VectorReference & operator -= (const TrilinosScalar &s) const;

                                         /**
                                          * Multiply the referenced element of
                                          * the vector by <tt>s</tt>.
                                          */
        const VectorReference & operator *= (const TrilinosScalar &s) const;

                                         /**
                                          * Divide the referenced element of
                                          * the vector by <tt>s</tt>.
                                          */
        const VectorReference & operator /= (const TrilinosScalar &s) const;

                                         /**
                                          * Convert the reference to an actual
                                          * value, i.e. return the value of
                                          * the referenced element of the
                                          * vector.
                                          */
        operator TrilinosScalar () const;

                                         /**
                                          * Exception
                                          */
        DeclException1 (ExcTrilinosError,
                        int,
                        << "An error with error number " << arg1
                        << " occured while calling a Trilinos function");

                                         /**
                                          * Exception
                                          */
        DeclException3 (ExcAccessToNonLocalElement,
                        int, int, int,
                        << "You tried to access element " << arg1
                        << " of a distributed vector, but only elements "
                        << arg2 << " through " << arg3
                        << " are stored locally and can be accessed.");

      private:
                                         /**
                                          * Point to the vector we are
                                          * referencing.
                                          */
        Vector   &vector;

                                         /**
                                          * Index of the referenced element
                                          * of the vector.
                                          */
        const unsigned int  index;

                                         /**
                                          * Make the vector class a friend, so
                                          * that it can create objects of the
                                          * present type.
                                          */
        friend class ::dealii::TrilinosWrappers::Vector;
    };
  }
                                   /**
                                    * @endcond
                                    */


/** 
 * This class implements a wrapper to use the Trilinos distributed vector
 * class Epetra_FEVector. This is precisely the kind of vector we deal with
 * all the time - we probably get it from some assembly process, where also
 * entries not locally owned might need to written and hence need to be
 * forwarded to the owner. This class is designed to be used in a distributed
 * memory architecture with an MPI compiler on the bottom, but works equally
 * well also for serial processes. The only requirement for this class to
 * work is that Trilinos is installed with the respective compiler as a
 * basis.
 *
 * The interface of this class is modeled after the existing Vector
 * class in deal.II. It has almost the same member functions, and is often
 * exchangable. However, since Trilinos only supports a single scalar type
 * (double), it is not templated, and only works with that type.
 *
 * Note that Trilinos only guarantees that operations do what you expect if the
 * function @p GlobalAssemble has been called after vector assembly in order to
 * distribute the data. Therefore, you need to call Vector::compress()
 * before you actually use the vectors.
 *
  * <h3>Parallel communication model</h3>
 *
 * The parallel functionality of Trilinos is built on top of the Message Passing
 * Interface (MPI). MPI's communication model is built on collective
 * communications: if one process wants something from another, that other
 * process has to be willing to accept this communication. A process cannot
 * query data from another process by calling a remote function, without that
 * other process expecting such a transaction. The consequence is that most of
 * the operations in the base class of this class have to be called
 * collectively. For example, if you want to compute the l2 norm of a parallel
 * vector, @em all processes across which this vector is shared have to call
 * the @p l2_norm function. If you don't do this, but instead only call the @p
 * l2_norm function on one process, then the following happens: This one
 * process will call one of the collective MPI functions and wait for all the
 * other processes to join in on this. Since the other processes don't call
 * this function, you will either get a time-out on the first process, or,
 * worse, by the time the next a callto a Trilinos function generates an MPI
 * message on the other processes , you will get a cryptic message that only a
 * subset of processes attempted a communication. These bugs can be very hard
 * to figure out, unless you are well-acquainted with the communication model
 * of MPI, and know which functions may generate MPI messages.
 *
 * One particular case, where an MPI message may be generated unexpectedly is
 * discussed below.
 *
 * <h3>Accessing individual elements of a vector</h3>
 *
 * Trilinos does allow read access to individual elements of a vector, but in the
 * distributed case only to elements that are stored locally. We implement
 * this through calls like <tt>d=vec(i)</tt>. However, if you access an
 * element outside the locally stored range, an exception is generated.
 *
 * In contrast to read access, Trilinos (and the respective deal.II wrapper
 * classes) allow to write (or add) to individual elements of vectors, even if
 * they are stored on a different process. You can do this writing, for
 * example, <tt>vec(i)=d</tt> or <tt>vec(i)+=d</tt>, or similar
 * operations. There is one catch, however, that may lead to very confusing
 * error messages: Trilinos requires application programs to call the compress()
 * function when they switch from adding, to elements to writing to
 * elements. The reasoning is that all processes might accumulate addition
 * operations to elements, even if multiple processes write to the same
 * elements. By the time we call compress() the next time, all these additions
 * are executed. However, if one process adds to an element, and another
 * overwrites to it, the order of execution would yield non-deterministic
 * behavior if we don't make sure that a synchronisation with compress()
 * happens in between.
 *
 * In order to make sure these calls to compress() happen at the appropriate
 * time, the deal.II wrappers keep a state variable that store which is the
 * presently allowed operation: additions or writes. If it encounters an
 * operation of the opposite kind, it calls compress() and flips the
 * state. This can sometimes lead to very confusing behavior, in code that may
 * for example look like this:
 * @verbatim
 *   TrilinosWrappers::Vector vector;
 *   ...
 *                   // do some write operations on the vector
 *   for (unsigned int i=0; i<vector->size(); ++i)
 *     vector(i) = i;
 *
 *                   // do some additions to vector elements, but
 *                   // only for some elements
 *   for (unsigned int i=0; i<vector->size(); ++i)
 *     if (some_condition(i) == true)
 *       vector(i) += 1;
 *
 *                   // do another collective operation
 *   const double norm = vector->l2_norm();
 * @endverbatim
 *
 * This code can run into trouble: by the time we see the first addition
 * operation, we need to flush the overwrite buffers for the vector, and the
 * deal.II library will do so by calling compress(). However, it will only do
 * so for all processes that actually do an addition -- if the condition is
 * never true for one of the processes, then this one will not get to the
 * actual compress() call, whereas all the other ones do. This gets us into
 * trouble, since all the other processes hang in the call to flush the write
 * buffers, while the one other process advances to the call to compute the l2
 * norm. At this time, you will get an error that some operation was attempted
 * by only a subset of processes. This behavior may seem surprising, unless
 * you know that write/addition operations on single elements may trigger this
 * behavior.
 *
 * The problem described here may be avoided by placing additional calls to
 * compress(), or making sure that all processes do the same type of
 * operations at the same time, for example by placing zero additions if
 * necessary.
 *
 * @ingroup TrilinosWrappers
 * @ingroup Vectors
 * @see @ref SoftwareTrilinos
 * @author Martin Kronbichler, Wolfgang Bangerth, 2008
 */
  class Vector
  {
    public:
                                       /**
                                        * Declare some of the standard types
                                        * used in all containers. These types
                                        * parallel those in the <tt>C</tt>
                                        * standard libraries <tt>vector<...></tt>
                                        * class.
                                        */
      typedef TrilinosScalar            value_type;
      typedef TrilinosScalar            real_type;
      typedef size_t                    size_type;
      typedef internal::VectorReference reference;
      typedef const internal::VectorReference const_reference;

                                       /**
                                        * Default constructor. It doesn't do
                                        * anything, derived classes will have
                                        * to initialize the data.
                                        */
      Vector ();
                                       /**
                                        * One of the constructors that
				        * actually builds a vector. This
				        * one requires prior knowledge
                                        * of the size of the vector and
				        * a communicator from 
                                        * Epetra_CommSerial or Epetra_CommMpi,
				        * depending on whether we use a 
				        * serial or parallel MPI-based program.
				        * This command distributes the
				        * vector linearly among the processes,
				        * from the beginning to the end, 
				        * so you might want to use some
				        * more advanced mapping and the
				        * constructor with argument
				        * Epetra_Map.
                                        */
      Vector (unsigned int GlobalSize, Epetra_Comm &Comm);

                                       /**
				        * This constructor takes an
				        * Epetra_Map that already knows how
				        * to distribute the individual 
				        * components among the MPI processors,
				        * including the size of the vector.
                                        */
      Vector (const Epetra_Map &InputMap);

                                       /**
                                        * Copy constructor. Sets the dimension
                                        * to that of the given vector and uses
				        * the map of that vector, but
                                        * does not copy any element. Instead,
				        * the memory will remain untouched
				        * in case <tt>fast</tt> is false and
				        * initialized with zero otherwise.
                                        */
      Vector (const Vector &v, 
	      const bool    fast = false);

                                       /**
                                        * Destructor
                                        */
      virtual ~Vector ();

				       /** 
				        * Reinit functionality. This function
				        * destroys the old vector content 
				        * and generates a new one based on
				        * the input map.
				        */
      void reinit (const Epetra_Map &input_map);

				       /** 
				        * Reinit functionality. This function
				        * copies the vector v to the current
				        * one.
				        */
      void reinit (const Vector &v,
		   const bool    fast = false);

                                       /**
                                        * Release all memory and return
                                        * to a state just like after
                                        * having called the default
                                        * constructor.
                                        */
      void clear ();

                                       /**
                                        * Compress the underlying
                                        * representation of the Trilinos object,
                                        * i.e. flush the buffers of the vector
                                        * object if it has any. This function
                                        * is necessary after writing into a
                                        * vector element-by-element and before
                                        * anything else can be done on it.
                                        */
      void compress ();

                                       /**
                                        * Set all components of the vector to
                                        * the given number @p s. Simply pass
                                        * this down to the Trilinos Epetra
                                        * object, but we still need to declare
                                        * this function to make the example
                                        * given in the discussion about making
                                        * the constructor explicit work.
                                        *
                                        * Since the semantics of assigning a
                                        * scalar to a vector are not
                                        * immediately clear, this operator
                                        * should really only be used if you
                                        * want to set the entire vector to
                                        * zero. This allows the intuitive
                                        * notation <tt>v=0</tt>. Assigning
                                        * other values is deprecated and may
                                        * be disallowed in the future.
                                        */
      Vector & operator = (const TrilinosScalar s);

                                       /**
                                        * Copy the given vector. Resize the
                                        * present vector if necessary.
                                        */
      Vector & operator = (const Vector &v);

                                       /**
                                        * Test for equality. This function
                                        * assumes that the present vector and
                                        * the one to compare with have the same
                                        * size already, since comparing vectors
                                        * of different sizes makes not much
                                        * sense anyway.
                                        */
      bool operator == (const Vector &v) const;

                                       /**
                                        * Test for inequality. This function
                                        * assumes that the present vector and
                                        * the one to compare with have the same
                                        * size already, since comparing vectors
                                        * of different sizes makes not much
                                        * sense anyway.
                                        */
      bool operator != (const Vector &v) const;

                                       /**
                                        * Return the global dimension of the
                                        * vector.
                                        */
      unsigned int size () const;

                                       /**
                                        * Return the local dimension of the
                                        * vector, i.e. the number of elements
                                        * stored on the present MPI
                                        * process. For sequential vectors,
                                        * this number is the same as size(),
                                        * but for parallel vectors it may be
                                        * smaller.
					*
					* To figure out which elements
					* exactly are stored locally,
					* use local_range().
                                        */
      unsigned int local_size () const;

                                       /**
					* Return a pair of indices
					* indicating which elements of
					* this vector are stored
					* locally. The first number is
					* the index of the first
					* element stored, the second
					* the index of the one past
					* the last one that is stored
					* locally. If this is a
					* sequential vector, then the
					* result will be the pair
					* (0,N), otherwise it will be
					* a pair (i,i+n), where
					* <tt>n=local_size()</tt>.
					*/
      std::pair<unsigned int, unsigned int> local_range () const;

				       /**
					* Return whether @p index is
					* in the local range or not,
					* see also local_range().
					*/
      bool in_local_range (const unsigned int index) const;

                                       /**
                                        * Provide access to a given element,
                                        * both read and write.
                                        */
      reference
      operator () (const unsigned int index);

                                       /**
                                        * Provide read-only access to an
                                        * element. This is equivalent to
				        * the <code>el()</code> command.
                                        */
      TrilinosScalar
      operator () (const unsigned int index) const;

                                       /**
                                        * Return the value of the vector entry
                                        * <i>i</i>. Note that this function
				        * does only work properly when 
				        * we request a data stored on the
				        * local processor. The function will
				        * throw an exception in case the 
				        * elements sits on another process.
                                        */
      TrilinosScalar el (const unsigned int index) const;

                                       /**
                                        * A collective set operation: instead
                                        * of setting individual elements of a
                                        * vector, this function allows to set
                                        * a whole set of elements at once. The
                                        * indices of the elements to be set
                                        * are stated in the first argument,
                                        * the corresponding values in the
                                        * second.
                                        */
      void set (const std::vector<unsigned int>    &indices,
		const std::vector<TrilinosScalar>  &values);

				       /**
				        * This is a second collective set 
				        * operation. As a difference, this
				        * function takes a deal.II vector
				        * of values.
				        */
      void set (const std::vector<unsigned int>        &indices,
		const ::dealii::Vector<TrilinosScalar> &values);

                                       /**
				        * This collective set operation is
				        * of lower level and can handle 
				        * anything else &ndash; the only
				        * thing you have to provide is 
				        * an address where all the indices
				        * are stored and the number of 
				        * elements to be set.
				        */
      void set (const unsigned int    n_elements,
		const unsigned int   *indices,
		const TrilinosScalar *values);

				       /**
                                        * A collective add operation: This
					* function adds a whole set of values
					* stored in @p values to the vector
					* components specified by @p indices.
                                        */
      void add (const std::vector<unsigned int>   &indices,
		const std::vector<TrilinosScalar> &values);

				       /**
				        * This is a second collective add 
				        * operation. As a difference, this
				        * function takes a deal.II vector
				        * of values.
				        */
      void add (const std::vector<unsigned int>        &indices,
		const ::dealii::Vector<TrilinosScalar> &values);

				      /**
				       * Take an address where n_elements
				       * are stored contiguously and add
				       * them into the vector.
				       */
      void add (const unsigned int    n_elements,
		const unsigned int   *indices,
		const TrilinosScalar *values);

                                       /**
                                        * Return the scalar (inner) product of
                                        * two vectors. The vectors must have the
                                        * same size.
                                        */
      TrilinosScalar operator * (const Vector &vec) const;

                                       /**
                                        * Return square of the $l_2$-norm.
                                        */
      real_type norm_sqr () const;

                                       /**
                                        * Mean value of the elements of
                                        * this vector.
                                        */
      TrilinosScalar mean_value () const;

                                       /**
                                        * $l_1$-norm of the vector.
                                        * The sum of the absolute values.
                                        */
      real_type l1_norm () const;

                                       /**
                                        * $l_2$-norm of the vector.  The
                                        * square root of the sum of the
                                        * squares of the elements.
                                        */
      real_type l2_norm () const;

                                       /**
                                        * $l_p$-norm of the vector. The
                                        * <i>p</i>th root of the sum of the 
				        * <i>p</i>th
                                        * powers of the absolute values
                                        * of the elements.
                                        */
      real_type lp_norm (const TrilinosScalar p) const;

                                       /**
                                        * Maximum absolute value of the
                                        * elements.
                                        */
      real_type linfty_norm () const;

                                       /**
                                        * Return whether the vector contains
                                        * only elements with value zero. This
                                        * function is mainly for internal
                                        * consistency checks and should
                                        * seldomly be used when not in debug
                                        * mode since it uses quite some time.
                                        */
      bool all_zero () const;

                                       /**
                                        * Return @p true if the vector has no
                                        * negative entries, i.e. all entries
                                        * are zero or positive. This function
                                        * is used, for example, to check
                                        * whether refinement indicators are
                                        * really all positive (or zero).
                                        */
      bool is_non_negative () const;

                                       /**
                                        * Multiply the entire vector by a
                                        * fixed factor.
                                        */
      Vector & operator *= (const TrilinosScalar factor);

                                       /**
                                        * Divide the entire vector by a
                                        * fixed factor.
                                        */
      Vector & operator /= (const TrilinosScalar factor);

                                       /**
                                        * Add the given vector to the present
                                        * one.
                                        */
      Vector & operator += (const Vector &V);

                                       /**
                                        * Subtract the given vector from the
                                        * present one.
                                        */
      Vector & operator -= (const Vector &V);

                                       /**
                                        * Addition of @p s to all
                                        * components. Note that @p s is a
                                        * scalar and not a vector.
                                        */
      void add (const TrilinosScalar s);

                                       /**
                                        * Simple vector addition, equal to the
                                        * <tt>operator +=</tt>.
                                        */
      void add (const Vector &V);

                                       /**
                                        * Simple addition of a multiple of a
                                        * vector, i.e. <tt>*this = a*V</tt>.
                                        */
      void add (const TrilinosScalar a, const Vector &V);

                                       /**
                                        * Multiple addition of scaled vectors,
                                        * i.e. <tt>*this = a*V + b*W</tt>.
                                        */
      void add (const TrilinosScalar a, const Vector &V,
                const TrilinosScalar b, const Vector &W);

                                       /**
                                        * Scaling and simple vector addition,
                                        * i.e.
                                        * <tt>*this = s*(*this) + V</tt>.
                                        */
      void sadd (const TrilinosScalar s,
                 const Vector        &V);

                                       /**
                                        * Scaling and simple addition, i.e.
                                        * <tt>*this = s*(*this) + a*V</tt>.
                                        */
      void sadd (const TrilinosScalar s,
                 const TrilinosScalar a,
                 const Vector        &V);

                                       /**
                                        * Scaling and multiple addition.
                                        */
      void sadd (const TrilinosScalar s,
                 const TrilinosScalar a,
                 const Vector        &V,
                 const TrilinosScalar b,
                 const Vector        &W);

                                       /**
                                        * Scaling and multiple addition.
                                        * <tt>*this = s*(*this) + a*V + b*W + c*X</tt>.
                                        */
      void sadd (const TrilinosScalar s,
                 const TrilinosScalar a,
                 const Vector        &V,
                 const TrilinosScalar b,
                 const Vector        &W,
                 const TrilinosScalar c,
                 const Vector        &X);

                                       /**
                                        * Scale each element of this
                                        * vector by the corresponding
                                        * element in the argument. This
                                        * function is mostly meant to
                                        * simulate multiplication (and
                                        * immediate re-assignment) by a
                                        * diagonal scaling matrix.
                                        */
      void scale (const Vector &scaling_factors);

                                       /**
                                        * Assignment <tt>*this = a*V</tt>.
                                        */
      void equ (const TrilinosScalar a, const Vector &V);

                                       /**
                                        * Assignment <tt>*this = a*V + b*W</tt>.
                                        */
      void equ (const TrilinosScalar a, const Vector &V,
                const TrilinosScalar b, const Vector &W);

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
      void ratio (const Vector &a,
                  const Vector &b);

				     /**
				      *  Output of vector in user-defined
				      *  format in analogy to the 
				      *  dealii::Vector<number> class.
				      */
    void print (const char* format = 0) const;

                                       /**
                                        * Print to a
                                        * stream. @p precision denotes
                                        * the desired precision with
                                        * which values shall be printed,
                                        * @p scientific whether
                                        * scientific notation shall be
                                        * used. If @p across is
                                        * @p true then the vector is
                                        * printed in a line, while if
                                        * @p false then the elements
                                        * are printed on a separate line
                                        * each.
                                        */
      void print (std::ostream       &out,
                  const unsigned int  precision  = 3,
                  const bool          scientific = true,
                  const bool          across     = true) const;

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
                                        * data around. Note that the 
				        * vectors need to be of the same
				        * size and base on the same 
				        * map.
                                        *
                                        * This function is analog to the
                                        * the @p swap function of all C
                                        * standard containers. Also,
                                        * there is a global function
                                        * <tt>swap(u,v)</tt> that simply calls
                                        * <tt>u.swap(v)</tt>, again in analogy
                                        * to standard functions.
                                        */
      void swap (Vector &v);

				       /**
					* Estimate for the memory
					* consumption (not implemented
					* for this class).
					*/
      unsigned int memory_consumption () const;

				       /**
					* Exception
					*/
      DeclException1 (ExcTrilinosError,
		      int,
		      << "An error with error number " << arg1
		      << " occured while calling a Trilinos function");

                                       /**
                                        * Exception
                                        */
      DeclException3 (ExcAccessToNonlocalElement,
		      int, int, int,
		      << "You tried to access element " << arg1
		      << " of a distributed vector, but only entries "
		      << arg2 << " through " << arg3
		      << " are stored locally and can be accessed.");


                                       /**
                                        * An Epetra map used to map vector data
                                        * accross multiple processes. This is
                                        * the communicator and data distribution
				        * object common to all
                                        * Trilinos objects used by deal.II.
                                        */
      Epetra_Map map;

                                       /**
                                        * An Epetra distibuted vector type.
                                        * Requires an existing Epetra_Map for
                                        * storing data.
                                        */
      std::auto_ptr<Epetra_FEVector> vector;


    private:
                                       /**
                                        * Trilinos doesn't allow to mix additions
                                        * to matrix entries and overwriting
                                        * them (to make synchronisation of
                                        * parallel computations
                                        * simpler). The way we do it is to,
                                        * for each access operation, store
                                        * whether it is an insertion or an
                                        * addition. If the previous one was of
                                        * different type, then we first have
                                        * to flush the Trilinos buffers;
                                        * otherwise, we can simply go on.
				        * Luckily, Trilinos has an object
				        * for this which does already all
				        * the parallel communications in
				        * such a case, so we simply use their
				        * model, which stores whether the
				        * last operation was an addition
				        * or an insertion.
                                        */
      Epetra_CombineMode last_action;

                                       /**
                                        * Make the reference class a friend.
                                        */
      friend class internal::VectorReference;

  };



// ------------------- inline and template functions --------------

/**
 * Global function swap which overloads the default implementation
 * of the C standard library which uses a temporary object. The
 * function simply exchanges the data of the two vectors.
 *
 * @ingroup TrilinosWrappers
 * @relates TrilinosWrappers::Vector
 * @author Martin Kronbichler, Wolfgang Bangerth, 2008
 */
  inline
  void swap (Vector &u, Vector &v)
  {
    u.swap (v);
  }

#ifndef DOXYGEN
  namespace internal
  {
    inline
    VectorReference::VectorReference (Vector            &vector,
                                      const unsigned int index)
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
  Vector::in_local_range (const unsigned int index) const
  {
    std::pair<unsigned int, unsigned int> range = local_range();

    return ((index >= range.first) && (index <  range.second));
  }



  inline
  Vector &
  Vector::operator = (const Vector &v)
  {
				    // if the vectors have different sizes,
				    // then first resize the present one
    if (!map.SameAs(v.map))
      {
	map = v.map;
	  vector = std::auto_ptr<Epetra_FEVector> 
			(new Epetra_FEVector(*v.vector));
      }
    else
      *vector = *v.vector;
	  
    last_action = Insert;
    
    return *this;
  }



  inline
  internal::VectorReference
  Vector::operator () (const unsigned int index)
  {
    return internal::VectorReference (*this, index);
  }



  inline
  TrilinosScalar
  Vector::operator () (const unsigned int index) const
  {
    TrilinosScalar value = el(index);

    return value;
  }


#endif // DOXYGEN
}

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS

/*----------------------------   trilinos_vector.h     ---------------------------*/

#endif
/*----------------------------   trilinos_vector.h     ---------------------------*/

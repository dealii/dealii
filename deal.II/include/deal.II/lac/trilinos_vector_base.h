//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__trilinos_vector_base_h
#define __deal2__trilinos_vector_base_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_USE_TRILINOS

#include <deal.II/base/utilities.h>
#  include <deal.II/base/std_cxx1x/shared_ptr.h>
#  include <deal.II/base/subscriptor.h>
#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/vector.h>

#  include <vector>
#  include <utility>
#  include <memory>

#  define TrilinosScalar double
#  include "Epetra_ConfigDefs.h"
#  ifdef DEAL_II_COMPILER_SUPPORTS_MPI // only if MPI is installed
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
					* This class implements a
					* wrapper for accessing the
					* Trilinos vector in the same
					* way as we access deal.II
					* objects: it is initialized
					* with a vector and an element
					* within it, and has a
					* conversion operator to
					* extract the scalar value of
					* this element. It also has a
					* variety of assignment
					* operator for writing to this
					* one element.  @ingroup
					* TrilinosWrappers
					*/
    class VectorReference
    {
      private:
                                       /**
					* Constructor. It is made
					* private so as to only allow
					* the actual vector class to
					* create it.
					*/
        VectorReference (VectorBase        &vector,
			 const unsigned int index);

      public:
                                       /**
					* This looks like a copy
					* operator, but does something
					* different than usual. In
					* particular, it does not copy
					* the member variables of this
					* reference. Rather, it
					* handles the situation where
					* we have two vectors @p v and
					* @p w, and assign elements
					* like in
					* <tt>v(i)=w(i)</tt>. Here,
					* both left and right hand
					* side of the assignment have
					* data type VectorReference,
					* but what we really mean is
					* to assign the vector
					* elements represented by the
					* two references. This
					* operator implements this
					* operation. Note also that
					* this allows us to make the
					* assignment operator const.
					*/
	const VectorReference &
	  operator = (const VectorReference &r) const;

                                       /**
					* Set the referenced element of the
					* vector to <tt>s</tt>.
					*/
	const VectorReference &
	  operator = (const TrilinosScalar &s) const;

                                       /**
					* Add <tt>s</tt> to the
					* referenced element of the
					* vector->
					*/
	const VectorReference &
	  operator += (const TrilinosScalar &s) const;

                                       /**
					* Subtract <tt>s</tt> from the
					* referenced element of the
					* vector->
					*/
	const VectorReference &
	  operator -= (const TrilinosScalar &s) const;

	                               /**
					* Multiply the referenced
					* element of the vector by
					* <tt>s</tt>.
					*/
	const VectorReference &
	  operator *= (const TrilinosScalar &s) const;

                                       /**
					* Divide the referenced
					* element of the vector by
					* <tt>s</tt>.
					*/
	const VectorReference &
	  operator /= (const TrilinosScalar &s) const;

                                       /**
					* Convert the reference to an
					* actual value, i.e. return
					* the value of the referenced
					* element of the vector.
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
	VectorBase   &vector;

                                       /**
					* Index of the referenced element
					* of the vector.
					*/
	const unsigned int  index;

                                       /**
					* Make the vector class a
					* friend, so that it can
					* create objects of the
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
                                        * Declare some of the standard
                                        * types used in all
                                        * containers. These types
                                        * parallel those in the
                                        * <tt>C</tt> standard libraries
                                        * <tt>vector<...></tt> class.
                                        */
      typedef TrilinosScalar            value_type;
      typedef TrilinosScalar            real_type;
      typedef std::size_t               size_type;
      typedef internal::VectorReference reference;
      typedef const internal::VectorReference const_reference;

				       /**
					* @name 1: Basic Object-handling
					*/
				       //@{

                                       /**
                                        * Default constructor that
                                        * generates an empty (zero size)
                                        * vector. The function
                                        * <tt>reinit()</tt> will have to
                                        * give the vector the correct
                                        * size and distribution among
                                        * processes in case of an MPI
                                        * run.
                                        */
      VectorBase ();

                                       /**
                                        * Copy constructor. Sets the
                                        * dimension to that of the given
                                        * vector, and copies all the
                                        * elements.
                                        */
      VectorBase (const VectorBase &v);

                                       /**
                                        * Destructor
                                        */
      virtual ~VectorBase ();

                                       /**
                                        * Release all memory and return
                                        * to a state just like after
                                        * having called the default
                                        * constructor.
                                        */
      void clear ();

				       /**
					* Reinit functionality, sets the
					* dimension and possibly the
					* parallel partitioning (Epetra_Map)
					* of the calling vector to the
					* settings of the input vector.
					*/
      void reinit (const VectorBase &v,
		   const bool        fast = false);

                                       /**
                                        * Compress the underlying
                                        * representation of the Trilinos
                                        * object, i.e. flush the buffers
                                        * of the vector object if it has
                                        * any. This function is
                                        * necessary after writing into a
                                        * vector element-by-element and
                                        * before anything else can be
                                        * done on it.
					*
					* The (defaulted) argument can
					* be used to specify the
					* compress mode
					* (<code>Add</code> or
					* <code>Insert</code>) in case
					* the vector has not been
					* written to since the last
					* time this function was
					* called. The argument is
					* ignored if the vector has
					* been added or written to
					* since the last time
					* compress() was called.
					*
					* See @ref GlossCompress "Compressing distributed objects"
					* for more information.
					* more information.
                                        */
      void compress (const Epetra_CombineMode last_action = Zero);

				       /**
					* Returns the state of the
					* vector, i.e., whether
					* compress() has already been
					* called after an operation
					* requiring data exchange.
					*/
      bool is_compressed () const;

                                       /**
                                        * Set all components of the
                                        * vector to the given number @p
                                        * s. Simply pass this down to
                                        * the Trilinos Epetra object,
                                        * but we still need to declare
                                        * this function to make the
                                        * example given in the
                                        * discussion about making the
                                        * constructor explicit work.
                                        *
                                        * Since the semantics of
                                        * assigning a scalar to a vector
                                        * are not immediately clear,
                                        * this operator should really
                                        * only be used if you want to
                                        * set the entire vector to
                                        * zero. This allows the
                                        * intuitive notation
                                        * <tt>v=0</tt>. Assigning other
                                        * values is deprecated and may
                                        * be disallowed in the future.
                                        */
      VectorBase &
	operator = (const TrilinosScalar s);

                                       /**
					* Copy function. This function takes
					* a VectorBase vector and copies all
					* the elements. The target vector
					* will have the same parallel
					* distribution as the calling
					* vector.
					*/
      VectorBase &
	operator = (const VectorBase &v);

                                       /**
					* Another copy function. This
					* one takes a deal.II vector and
					* copies it into a
					* TrilinosWrapper vector. Note
					* that since we do not provide
					* any Epetra_map that tells
					* about the partitioning of the
					* vector among the MPI
					* processes, the size of the
					* TrilinosWrapper vector has to
					* be the same as the size of the
					* input vector. In order to
					* change the map, use the
					* reinit(const Epetra_Map
					* &input_map) function.
					*/
      template <typename Number>
      VectorBase &
	operator = (const ::dealii::Vector<Number> &v);

                                       /**
                                        * Test for equality. This
                                        * function assumes that the
                                        * present vector and the one to
                                        * compare with have the same
                                        * size already, since comparing
                                        * vectors of different sizes
                                        * makes not much sense anyway.
                                        */
      bool operator == (const VectorBase &v) const;

                                       /**
                                        * Test for inequality. This
                                        * function assumes that the
                                        * present vector and the one to
                                        * compare with have the same
                                        * size already, since comparing
                                        * vectors of different sizes
                                        * makes not much sense anyway.
                                        */
      bool operator != (const VectorBase &v) const;

                                       /**
                                        * Return the global dimension of
                                        * the vector.
                                        */
      unsigned int size () const;

                                       /**
                                        * Return the local dimension of
                                        * the vector, i.e. the number of
                                        * elements stored on the present
                                        * MPI process. For sequential
                                        * vectors, this number is the
                                        * same as size(), but for
                                        * parallel vectors it may be
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
					* the index of the first element
					* stored, the second the index
					* of the one past the last one
					* that is stored locally. If
					* this is a sequential vector,
					* then the result will be the
					* pair (0,N), otherwise it will
					* be a pair (i,i+n), where
					* <tt>n=local_size()</tt>.
					*/
      std::pair<unsigned int, unsigned int> local_range () const;

				       /**
					* Return whether @p index is in
					* the local range or not, see
					* also local_range().
					*/
      bool in_local_range (const unsigned int index) const;

				       /**
					* Return if the vector contains ghost
					* elements.
					*/
      bool has_ghost_elements() const;

                                       /**
                                        * Return the scalar (inner)
                                        * product of two vectors. The
                                        * vectors must have the same
                                        * size.
                                        */
      TrilinosScalar operator * (const VectorBase &vec) const;

                                       /**
                                        * Return square of the
                                        * $l_2$-norm.
                                        */
      real_type norm_sqr () const;

                                       /**
                                        * Mean value of the elements of
                                        * this vector.
                                        */
      TrilinosScalar mean_value () const;

                                       /**
                                        * Compute the minimal value of
                                        * the elements of this vector.
                                        */
      TrilinosScalar minimal_value () const;

                                       /**
                                        * $l_1$-norm of the vector.  The
                                        * sum of the absolute values.
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
                                        * <i>p</i>th root of the sum of
                                        * the <i>p</i>th powers of the
                                        * absolute values of the
                                        * elements.
                                        */
      real_type lp_norm (const TrilinosScalar p) const;

                                       /**
                                        * Maximum absolute value of the
                                        * elements.
                                        */
      real_type linfty_norm () const;

                                       /**
                                        * Return whether the vector
                                        * contains only elements with
                                        * value zero. This function is
                                        * mainly for internal
                                        * consistency checks and should
                                        * seldomly be used when not in
                                        * debug mode since it uses quite
                                        * some time.
                                        */
      bool all_zero () const;

                                       /**
                                        * Return @p true if the vector
                                        * has no negative entries,
                                        * i.e. all entries are zero or
                                        * positive. This function is
                                        * used, for example, to check
                                        * whether refinement indicators
                                        * are really all positive (or
                                        * zero).
                                        */
      bool is_non_negative () const;
				       //@}


                                       /**
					* @name 2: Data-Access
					*/
				       //@{

                                       /**
                                        * Provide access to a given
                                        * element, both read and write.
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
                                        * Return the value of the vector
                                        * entry <i>i</i>. Note that this
                                        * function does only work
                                        * properly when we request a
                                        * data stored on the local
                                        * processor. The function will
                                        * throw an exception in case the
                                        * elements sits on another
                                        * process.
                                        */
      TrilinosScalar el (const unsigned int index) const;

                                       /**
                                        * A collective set operation:
                                        * instead of setting individual
                                        * elements of a vector, this
                                        * function allows to set a whole
                                        * set of elements at once. The
                                        * indices of the elements to be
                                        * set are stated in the first
                                        * argument, the corresponding
                                        * values in the second.
                                        */
      void set (const std::vector<unsigned int>    &indices,
		const std::vector<TrilinosScalar>  &values);

				       /**
				        * This is a second collective
				        * set operation. As a
				        * difference, this function
				        * takes a deal.II vector of
				        * values.
				        */
      void set (const std::vector<unsigned int>        &indices,
		const ::dealii::Vector<TrilinosScalar> &values);
				       //@}


				       /**
					* @name 3: Modification of vectors
					*/
				       //@{

                                       /**
				        * This collective set operation
				        * is of lower level and can
				        * handle anything else &mdash;
				        * the only thing you have to
				        * provide is an address where
				        * all the indices are stored and
				        * the number of elements to be
				        * set.
				        */
      void set (const unsigned int    n_elements,
		const unsigned int   *indices,
		const TrilinosScalar *values);

				       /**
                                        * A collective add operation:
                                        * This funnction adds a whole
                                        * set of values stored in @p
                                        * values to the vector
                                        * components specified by @p
                                        * indices.
                                        */
      void add (const std::vector<unsigned int>   &indices,
		const std::vector<TrilinosScalar> &values);

				       /**
				        * This is a second collective
				        * add operation. As a
				        * difference, this function
				        * takes a deal.II vector of
				        * values.
				        */
      void add (const std::vector<unsigned int>        &indices,
		const ::dealii::Vector<TrilinosScalar> &values);

				      /**
				       * Take an address where
				       * <tt>n_elements</tt> are stored
				       * contiguously and add them into
				       * the vector. Handles all cases
				       * which are not covered by the
				       * other two <tt>add()</tt>
				       * functions above.
				       */
      void add (const unsigned int    n_elements,
		const unsigned int   *indices,
		const TrilinosScalar *values);

                                       /**
                                        * Multiply the entire vector by
                                        * a fixed factor.
                                        */
      VectorBase & operator *= (const TrilinosScalar factor);

                                       /**
                                        * Divide the entire vector by a
                                        * fixed factor.
                                        */
      VectorBase & operator /= (const TrilinosScalar factor);

                                       /**
                                        * Add the given vector to the
                                        * present one.
                                        */
      VectorBase & operator += (const VectorBase &V);

                                       /**
                                        * Subtract the given vector from
                                        * the present one.
                                        */
      VectorBase & operator -= (const VectorBase &V);

                                       /**
                                        * Addition of @p s to all
                                        * components. Note that @p s is
                                        * a scalar and not a vector.
                                        */
      void add (const TrilinosScalar s);

                                       /**
                                        * Simple vector addition, equal
                                        * to the <tt>operator
                                        * +=</tt>.
					*
					* Though, if the second argument
                                        * <tt>allow_different_maps</tt>
                                        * is set, then it is possible to
                                        * add data from a different map.
                                        */
      void add (const VectorBase &V,
		const bool        allow_different_maps = false);

                                       /**
                                        * Simple addition of a multiple
                                        * of a vector, i.e. <tt>*this =
                                        * a*V</tt>.
                                        */
      void add (const TrilinosScalar  a,
		const VectorBase     &V);

                                       /**
                                        * Multiple addition of scaled
                                        * vectors, i.e. <tt>*this = a*V
                                        * + b*W</tt>.
                                        */
      void add (const TrilinosScalar  a,
		const VectorBase     &V,
		const TrilinosScalar  b,
		const VectorBase     &W);

                                       /**
                                        * Scaling and simple vector
                                        * addition, i.e.  <tt>*this =
                                        * s*(*this) + V</tt>.
                                        */
      void sadd (const TrilinosScalar  s,
		 const VectorBase     &V);

                                       /**
                                        * Scaling and simple addition,
                                        * i.e.  <tt>*this = s*(*this) +
                                        * a*V</tt>.
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
                                        * Scaling and multiple addition.
                                        * <tt>*this = s*(*this) + a*V +
                                        * b*W + c*X</tt>.
                                        */
      void sadd (const TrilinosScalar  s,
		 const TrilinosScalar  a,
		 const VectorBase     &V,
		 const TrilinosScalar  b,
		 const VectorBase     &W,
		 const TrilinosScalar  c,
		 const VectorBase     &X);

                                       /**
                                        * Scale each element of this
                                        * vector by the corresponding
                                        * element in the argument. This
                                        * function is mostly meant to
                                        * simulate multiplication (and
                                        * immediate re-assignment) by a
                                        * diagonal scaling matrix.
                                        */
      void scale (const VectorBase &scaling_factors);

                                       /**
                                        * Assignment <tt>*this =
                                        * a*V</tt>.
                                        */
      void equ (const TrilinosScalar  a,
		const VectorBase     &V);

                                       /**
                                        * Assignment <tt>*this = a*V +
                                        * b*W</tt>.
                                        */
      void equ (const TrilinosScalar  a,
		const VectorBase     &V,
		const TrilinosScalar  b,
		const VectorBase     &W);

                                       /**
                                        * Compute the elementwise ratio
                                        * of the two given vectors, that
                                        * is let <tt>this[i] =
                                        * a[i]/b[i]</tt>. This is useful
                                        * for example if you want to
                                        * compute the cellwise ratio of
                                        * true to estimated error.
                                        *
                                        * This vector is appropriately
                                        * scaled to hold the result.
                                        *
                                        * If any of the <tt>b[i]</tt> is
                                        * zero, the result is
                                        * undefined. No attempt is made
                                        * to catch such situations.
                                        */
      void ratio (const VectorBase &a,
		  const VectorBase &b);
				       //@}


				       /**
					* @name 4: Mixed stuff
					*/
				       //@{

                                       /**
                                        * Return a const reference to the
                                        * underlying Trilinos
                                        * Epetra_MultiVector class.
                                        */
      const Epetra_MultiVector & trilinos_vector () const;

                                       /**
                                        * Return a (modifyable) reference to
                                        * the underlying Trilinos
                                        * Epetra_FEVector class.
                                        */
      Epetra_FEVector & trilinos_vector ();

                                       /**
                                        * Return a const reference to the
                                        * underlying Trilinos Epetra_Map
                                        * that sets the parallel
                                        * partitioning of the vector.
                                        */
      const Epetra_Map & vector_partitioner () const;

				       /**
					*  Output of vector in
					*  user-defined format in analogy
					*  to the dealii::Vector<number>
					*  class.
					*/
      void print (const char* format = 0) const;

                                       /**
                                        * Print to a stream. @p
                                        * precision denotes the desired
                                        * precision with which values
                                        * shall be printed, @p
                                        * scientific whether scientific
                                        * notation shall be used. If @p
                                        * across is @p true then the
                                        * vector is printed in a line,
                                        * while if @p false then the
                                        * elements are printed on a
                                        * separate line each.
                                        */
      void print (std::ostream       &out,
		  const unsigned int  precision  = 3,
		  const bool          scientific = true,
		  const bool          across     = true) const;

                                       /**
                                        * Swap the contents of this
                                        * vector and the other vector @p
                                        * v. One could do this operation
                                        * with a temporary variable and
                                        * copying over the data
                                        * elements, but this function is
                                        * significantly more efficient
                                        * since it only swaps the
                                        * pointers to the data of the
                                        * two vectors and therefore does
                                        * not need to allocate temporary
                                        * storage and move data
                                        * around. Note that the vectors
                                        * need to be of the same size
                                        * and base on the same map.
                                        *
                                        * This function is analog to the
                                        * the @p swap function of all C
                                        * standard containers. Also,
                                        * there is a global function
                                        * <tt>swap(u,v)</tt> that simply
                                        * calls <tt>u.swap(v)</tt>,
                                        * again in analogy to standard
                                        * functions.
                                        */
      void swap (VectorBase &v);

				       /**
					* Estimate for the memory
					* consumption in bytes.
					*/
      std::size_t memory_consumption () const;
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


    private:
				       /**
                                        * Trilinos doesn't allow to
                                        * mix additions to matrix
                                        * entries and overwriting them
                                        * (to make synchronisation of
                                        * parallel computations
                                        * simpler). The way we do it
                                        * is to, for each access
                                        * operation, store whether it
                                        * is an insertion or an
                                        * addition. If the previous
                                        * one was of different type,
                                        * then we first have to flush
                                        * the Trilinos buffers;
                                        * otherwise, we can simply go
                                        * on.  Luckily, Trilinos has
                                        * an object for this which
                                        * does already all the
                                        * parallel communications in
                                        * such a case, so we simply
                                        * use their model, which
                                        * stores whether the last
                                        * operation was an addition or
                                        * an insertion.
                                        */
      Epetra_CombineMode last_action;

				       /**
					* A boolean variable to hold
					* information on whether the
					* vector is compressed or not.
					*/
      bool compressed;

                                       /**
                                        * An Epetra distibuted vector
                                        * type. Requires an existing
                                        * Epetra_Map for storing data.
                                        */
      std_cxx1x::shared_ptr<Epetra_FEVector> vector;


                                       /**
                                        * Make the reference class a
                                        * friend.
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
    VectorReference::VectorReference (VectorBase        &vector,
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
  VectorBase::is_compressed () const
  {
    return compressed;
  }



  inline
  bool
  VectorBase::in_local_range (const unsigned int index) const
  {
    std::pair<unsigned int, unsigned int> range = local_range();

    return ((index >= range.first) && (index <  range.second));
  }



  inline
  bool
  VectorBase::has_ghost_elements() const
  {
    return vector->Map().UniqueGIDs()==false;
  }



  inline
  internal::VectorReference
  VectorBase::operator () (const unsigned int index)
  {
    return internal::VectorReference (*this, index);
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
  }



  inline
  void
  VectorBase::compress (const Epetra_CombineMode given_last_action)
  {
    Epetra_CombineMode mode = (last_action != Zero) ?
			      last_action : given_last_action;
#ifdef DEBUG
#  ifdef DEAL_II_COMPILER_SUPPORTS_MPI
				     // check that every process has decided
				     // to use the same mode. This will
				     // otherwise result in undefined
				     // behaviour in the call to
				     // GlobalAssemble().
    double double_mode = mode;
    Utilities::MPI::MinMaxAvg result
      = Utilities::MPI::min_max_avg (double_mode,
				     dynamic_cast<const Epetra_MpiComm*>
				     (&vector_partitioner().Comm())->GetMpiComm());
    Assert(result.max-result.min<1e-5,
	   ExcMessage ("Not all processors agree whether the last operation on "
		       "this vector was an addition or a set operation. This will "
		       "prevent the compress() operation from succeeding."));

#  endif
#endif

				 // Now pass over the information about
				 // what we did last to the vector.
    const int ierr = vector->GlobalAssemble(mode);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
    last_action = Zero;

    compressed = true;
  }



  inline
  VectorBase &
  VectorBase::operator = (const TrilinosScalar s)
  {

    Assert (numbers::is_finite(s), ExcNumberNotFinite());

    const int ierr = vector->PutScalar(s);

    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return *this;
  }



  inline
  void
  VectorBase::set (const std::vector<unsigned int>    &indices,
		   const std::vector<TrilinosScalar>  &values)
  {
    Assert (indices.size() == values.size(),
	    ExcDimensionMismatch(indices.size(),values.size()));

    set (indices.size(), &indices[0], &values[0]);
  }



  inline
  void
  VectorBase::set (const std::vector<unsigned int>        &indices,
		   const ::dealii::Vector<TrilinosScalar> &values)
  {
    Assert (indices.size() == values.size(),
	    ExcDimensionMismatch(indices.size(),values.size()));

    set (indices.size(), &indices[0], values.begin());
  }



  inline
  void
  VectorBase::set (const unsigned int    n_elements,
		   const unsigned int   *indices,
		   const TrilinosScalar *values)
  {
				     // if we have ghost values, do not allow
				     // writing to this vector at all.
    Assert (!has_ghost_elements(), ExcGhostsPresent());

    if (last_action == Add)
      vector->GlobalAssemble(Add);

    if (last_action != Insert)
      last_action = Insert;

    for (unsigned int i=0; i<n_elements; ++i)
      {
	const unsigned int row = indices[i];
	const int local_row = vector->Map().LID(indices[i]);
	if (local_row == -1)
	  {
	    const int ierr = vector->ReplaceGlobalValues (1,
							  (const int*)(&row),
							  &values[i]);
	    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
	    compressed = false;
	  }
	else
	  (*vector)[0][local_row] = values[i];
      }
  }



  inline
  void
  VectorBase::add (const std::vector<unsigned int>    &indices,
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
  VectorBase::add (const std::vector<unsigned int>        &indices,
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
  VectorBase::add (const unsigned int    n_elements,
		   const unsigned int   *indices,
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

    for (unsigned int i=0; i<n_elements; ++i)
      {
	const unsigned int row = indices[i];
	const int local_row = vector->Map().LID(row);
	if (local_row == -1)
	  {
	    const int ierr = vector->SumIntoGlobalValues (1,
							  (const int*)(&row),
							  &values[i]);
	    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
	    compressed = false;
	  }
	else
	  (*vector)[0][local_row] += values[i];
      }
  }



  inline
  unsigned int
  VectorBase::size () const
  {
    return (unsigned int) (vector->Map().MaxAllGID() + 1 -
			   vector->Map().MinAllGID());
  }



  inline
  unsigned int
  VectorBase::local_size () const
  {
    return (unsigned int) vector->Map().NumMyElements();
  }



  inline
  std::pair<unsigned int, unsigned int>
  VectorBase::local_range () const
  {
    int begin, end;
    begin = vector->Map().MinMyGID();
    end = vector->Map().MaxMyGID()+1;
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
    const unsigned int n_local = local_size();

                                        // loop over all the elements because
                                        // Trilinos does not support lp norms
    for (unsigned int i=0; i<n_local; i++)
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
    Assert (numbers::is_finite(s), ExcNumberNotFinite());

    unsigned int n_local = local_size();
    for (unsigned int i=0; i<n_local; i++)
      (*vector)[0][i] += s;
  }



  inline
  void
  VectorBase::add (const TrilinosScalar  a,
		   const VectorBase     &v)
  {
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
    Assert (local_size() == v.local_size(),
	    ExcDimensionMismatch(local_size(), v.local_size()));

    Assert (numbers::is_finite(s), ExcNumberNotFinite());

    const int ierr = vector->Update(1., *(v.vector), s);

    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  inline
  void
  VectorBase::sadd (const TrilinosScalar  s,
		    const TrilinosScalar  a,
		    const VectorBase     &v)
  {
    Assert (local_size() == v.local_size(),
	    ExcDimensionMismatch(local_size(), v.local_size()));

    Assert (numbers::is_finite(s), ExcNumberNotFinite());
    Assert (numbers::is_finite(a), ExcNumberNotFinite());

    const int ierr = vector->Update(a, *(v.vector), s);

    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  inline
  void
  VectorBase::sadd (const TrilinosScalar  s,
		    const TrilinosScalar  a,
		    const VectorBase     &v,
		    const TrilinosScalar  b,
		    const VectorBase     &w)
  {
    Assert (local_size() == v.local_size(),
	    ExcDimensionMismatch(local_size(), v.local_size()));
    Assert (local_size() == w.local_size(),
	    ExcDimensionMismatch(local_size(), w.local_size()));

    Assert (numbers::is_finite(s), ExcNumberNotFinite());
    Assert (numbers::is_finite(a), ExcNumberNotFinite());
    Assert (numbers::is_finite(b), ExcNumberNotFinite());

    const int ierr = vector->Update(a, *(v.vector), b, *(w.vector), s);

    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
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
    Assert (local_size() == v.local_size(),
	    ExcDimensionMismatch(local_size(), v.local_size()));
    Assert (local_size() == w.local_size(),
	    ExcDimensionMismatch(local_size(), w.local_size()));
    Assert (local_size() == x.local_size(),
	    ExcDimensionMismatch(local_size(), x.local_size()));

    Assert (numbers::is_finite(s), ExcNumberNotFinite());
    Assert (numbers::is_finite(a), ExcNumberNotFinite());
    Assert (numbers::is_finite(b), ExcNumberNotFinite());
    Assert (numbers::is_finite(c), ExcNumberNotFinite());

                                        // Update member can only
				        // input two other vectors so
				        // do it in two steps
    const int ierr = vector->Update(a, *(v.vector), b, *(w.vector), s);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    const int jerr = vector->Update(c, *(x.vector), 1.);
    Assert (jerr == 0, ExcTrilinosError(jerr));
  }



  inline
  void
  VectorBase::scale (const VectorBase &factors)
  {
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
    Assert (numbers::is_finite(a), ExcNumberNotFinite());

				   // If we don't have the same map, copy.
    if (vector->Map().SameAs(v.vector->Map())==false)
      {
	*vector = *v.vector;
	*this *= a;
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
    Assert (v.local_size() == w.local_size(),
	    ExcDimensionMismatch (v.local_size(), w.local_size()));

    Assert (numbers::is_finite(a), ExcNumberNotFinite());
    Assert (numbers::is_finite(b), ExcNumberNotFinite());

				   // If we don't have the same map, copy.
     if (vector->Map().SameAs(v.vector->Map())==false)
      {
	*vector = *v.vector;
	sadd(a, b, w);
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
    return static_cast<const Epetra_MultiVector&>(*vector);
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
    return static_cast<const Epetra_Map&>(vector->Map());
  }


#endif // DOXYGEN

}

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS

/*----------------------------   trilinos_vector_base.h     ---------------------------*/

#endif
/*----------------------------   trilinos_vector_base.h     ---------------------------*/

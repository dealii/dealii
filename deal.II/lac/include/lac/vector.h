//----------------------------  vector.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  vector.h  ---------------------------
#ifndef __deal2__vector_h
#define __deal2__vector_h


#include <base/config.h>
#include <base/exceptions.h>

#include <cstdio>


#ifdef DEAL_II_USE_PETSC
namespace PETScWrappers
{
  class Vector;
  namespace MPI
  {
    class Vector;
  }
}
#endif



/*! @addtogroup Vectors
 *@{
 */


/**
 * Numerical vector of data.  For this class there are different types
 * of functions available. The first type of function mesures the norm
 * of the vector in order to mesure its length in a suitable norm. The
 * second type support the abgebraic operation for vectors. The third
 * und last type helps us to manipulate the components of the vector.
 * As opposed to the array of the C++ standard library called
 * @p vector (with a lowercase "v"), this class implements an element
 * of a vector space suitable for numerical computations.
 *
 *
 * @sect2{On template instantiations}
 *
 * Member functions of this class are either implemented in this file
 * or in a file of the same name with suffix ``.templates.h''. For the
 * most common combinations of the template parameters, instantiations
 * of this class are provided in a file with suffix ``.cc'' in the
 * ``source'' directory. If you need an instantiation that is not
 * listed there, you have to include this file along with the
 * corresponding ``.templates.h'' file and instantiate the respective
 * class yourself.
 *
 * @author Guido Kanschat, Franz-Theo Suttmeier, Wolfgang Bangerth
 */
template <typename Number>
class Vector
{
  public:
				     /**
				      * Declare standard types used in all
				      * containers. These types parallel
				      * those in the <tt>C++</tt> standard libraries
				      * <tt>vector<...></tt> class.
				      */
    typedef Number            value_type;
    typedef value_type       *pointer;
    typedef const value_type *const_pointer;
    typedef value_type       *iterator;
    typedef const value_type *const_iterator;
    typedef value_type       &reference;
    typedef const value_type &const_reference;
    typedef size_t            size_type;


				     /**
				      * @name 1: Basic Object-handling 
				      */
				     //@{
				     /**
				      *  Constructor. Create a vector of
				      *  dimension zero.
				      */
    Vector ();
    
				     /**
				      * Copy-constructor. Sets the dimension
                                      * to that of the given vector, and
                                      * copies all elements.
				      *
				      * We would like to make this
				      * constructor explicit, but STL
				      * insists on using it implicitly.
				      */
    Vector (const Vector<Number> &v);


// note: I disabled this function for the time being, since gcc2.95
// does not respect the "explicit" keyword for template constructors.
// this leads to unwanted conversions and in some places to automatically
// generated temporaries, where this is not a good idea. [WB]
// 				     /**
// 				      * Copy constructor taking a vector of
// 				      * another data type. This will fail if
// 				      * there is no conversion path from
// 				      * @p OtherNumber to @p Number. Note that
// 				      * you may lose accuracy when copying
// 				      * to a vector with data elements with
// 				      * less accuracy.
// 				      */
//     template <typename OtherNumber>
//     explicit
//     Vector (const Vector<OtherNumber> &v);

#ifdef DEAL_II_USE_PETSC
                                     /**
                                      * Another copy constructor: copy the
                                      * values from a sequential PETSc wrapper
                                      * vector class. This copy constructor is
                                      * only available if PETSc was detected
                                      * during configuration time.
                                      */
    Vector (const PETScWrappers::Vector &v);

                                     /**
                                      * Another copy constructor: copy the
                                      * values from a parallel PETSc wrapper
                                      * vector class. This copy constructor is
                                      * only available if PETSc was detected
                                      * during configuration time.
                                      *
                                      * Note that due to the communication
                                      * model used in MPI, this operation can
                                      * only succeed if all processes do it at
                                      * the same time. I.e., it is not
                                      * possible for only one process to
                                      * obtain a copy of a parallel vector
                                      * while the other jobs do something
                                      * else.
                                      */
    Vector (const PETScWrappers::MPI::Vector &v);
#endif
    
				     /**
				      * Constructor. Set dimension to
				      * @p n and initialize all
				      * elements with zero.
				      *
				      * The constructor is made
				      * explicit to avoid accidents
				      * like this:
				      * <tt>v=0;</tt>. Presumably, the user
				      * wants to set every element of
				      * the vector to zero, but
				      * instead, what happens is this
				      * call:
				      * <tt>v=Vector@<number@>(0);</tt>,
				      * i.e. the vector is replaced by
				      * one of length zero.
				      */
    explicit Vector (const unsigned int n);

				     /**
				      * Initialize the vector with a
				      * given range of values pointed
				      * to by the iterators. This
				      * function is there in analogy
				      * to the @p std::vector class.
				      */
    template <typename InputIterator>
    Vector (const InputIterator first,
            const InputIterator last);
    
				     /**
				      * Destructor, deallocates
				      * memory. Made virtual to allow
				      * for derived classes to behave
				      * properly.
				      */
    virtual ~Vector ();

                                     /**
                                      * This function does nothing but is
                                      * there for compatibility with the
                                      * @p PETScWrappers::Vector class.
                                      * 
                                      * For the PETSc vector wrapper class,
                                      * thios function compresses the
                                      * underlying representation of the PETSc
                                      * object, i.e. flushes the buffers of
                                      * the vector object if it has any. This
                                      * function is necessary after writing
                                      * into a vector element-by-element and
                                      * before anything else can be done on
                                      * it.
                                      *
                                      * However, for the implementation of
                                      * this class, it is immaterial and thus
                                      * an empty function.
                                      */
    void compress () const;

				     /**
				      * Change the dimension of the vector to
				      * @p N. The reserved memory for this
				      * vector remains unchanged if possible,
				      * to make things faster; this may waste
				      * some memory, so keep this in mind.
				      * However, if <tt>N==0</tt> all memory is
				      * freed, i.e. if you want to resize the
				      * vector and release the memory not
				      * needed, you have to first call
				      * <tt>reinit(0)</tt> and then
				      * <tt>reinit(N)</tt>. This cited behaviour is
				      * analogous to that of the STL
				      * containers.
				      *
				      * If @p fast is false, the vector is
				      * filled by zeros. Otherwise, the
				      * elements are left an unspecified
				      * state.
				      */ 
    void reinit (const unsigned int N,
		 const bool         fast=false);
    
				     /**
				      * Change the dimension to that of the
				      * vector @p V. The same applies as for
				      * the other @p reinit function.
				      *
				      * The elements of @p V are not copied,
				      * i.e.  this function is the same as
				      * calling <tt>reinit (V.size(), fast)</tt>.
				      */
    template <typename Number2>
    void reinit (const Vector<Number2> &V,
		 const bool            fast=false);

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
				      */
    void swap (Vector<Number> &v);
    
				     /**
				      * Set all components of the vector to
				      * the given number @p s.
				      */
    Vector<Number> & operator= (const Number s);
    
				     /**
				      * Copy the given vector. Resize the
				      * present vector if necessary.
				      */
    Vector<Number> & operator= (const Vector<Number> &c);
    
				     /**
				      * Copy the given vector. Resize the
				      * present vector if necessary.
				      */
    template <typename Number2>
    Vector<Number> & operator= (const Vector<Number2> &v);

#ifdef DEAL_II_USE_PETSC
                                     /**
                                      * Another copy operator: copy the values
                                      * from a sequential PETSc wrapper vector
                                      * class. This operator is only available
                                      * if PETSc was detected during
                                      * configuration time.
                                      */
    Vector<Number> &
    operator = (const PETScWrappers::Vector &v);

                                     /**
                                      * Another copy operator: copy the values
                                      * from a parallel PETSc wrapper vector
                                      * class. This operator is only available
                                      * if PETSc was detected during
                                      * configuration time.
                                      *
                                      * Note that due to the communication
                                      * model used in MPI, this operation can
                                      * only succeed if all processes do it at
                                      * the same time. I.e., it is not
                                      * possible for only one process to
                                      * obtain a copy of a parallel vector
                                      * while the other jobs do something
                                      * else.
                                      */
    Vector<Number> &
    operator = (const PETScWrappers::MPI::Vector &v);
#endif

                                     /**
                                      * Test for equality. This function
                                      * assumes that the present vector and
                                      * the one to compare with have the same
                                      * size already, since comparing vectors
                                      * of different sizes makes not much
                                      * sense anyway.
                                      */
    template <typename Number2>
    bool operator == (const Vector<Number2> &v) const;
    
                                     /**
                                      * Test for inequality. This function
                                      * assumes that the present vector and
                                      * the one to compare with have the same
                                      * size already, since comparing vectors
                                      * of different sizes makes not much
                                      * sense anyway.
                                      */
    template <typename Number2>
    bool operator != (const Vector<Number2> &v) const;
    
				     /**
				      * Return the scalar product of
				      * two vectors.  The return type
				      * is the underlying type of
				      * @p this vector, so the return
				      * type and the accuracy with
				      * which it the result is
				      * computed depend on the order
				      * of the arguments of this
				      * vector.
				      */
    template <typename Number2>
    Number operator* (const Vector<Number2> &V) const;

				     /**
				      * Return square of the $l_2$-norm.
				      */
    Number norm_sqr () const;

				     /**
				      * Mean value of the elements of
				      * this vector.
				      */
    Number mean_value () const;

				     /**
				      * $l_1$-norm of the vector.
				      * The sum of the absolute values.
				      */
    Number l1_norm () const;

				     /**
				      * $l_2$-norm of the vector.  The
				      * square root of the sum of the
				      * squares of the elements.
				      */
    Number l2_norm () const;

				     /**
				      * $l_p$-norm of the vector. The
				      * pth root of the sum of the pth
				      * powers of the absolute values
				      * of the elements.
				      */
    Number lp_norm (const Number p) const;

				     /**
				      * Maximum absolute value of the
				      * elements.
				      */
    Number linfty_norm () const;

				     /**
				      * Return dimension of the vector.
				      */
    unsigned int size () const;

				     /**
				      * Return whether the vector contains only
				      * elements with value zero. This function
				      * is mainly for internal consistency
				      * checks and should seldomly be used when
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
                                      */
    bool is_non_negative () const;
    
				     /**
				      * Make the @p Vector class a bit like the
				      * <tt>vector<></tt> class of the C++ standard
				      * library by returning iterators to
				      * the start and end of the elements of this
				      * vector.
				      */
    iterator begin ();

				     /**
				      * Return constant iterator to the start of
				      * the vectors.
				      */
    const_iterator begin () const;

				     /**
				      * Return an iterator pointing to the
				      * element past the end of the array.
				      */
    iterator end ();

    				     /**
				      * Return a constant iterator pointing to
				      * the element past the end of the array.
				      */
    const_iterator end () const;  
				     //@}


				     /**
				      * @name 2: Data-Access
				      */
				     //@{
				     /**
				      * Access the value of the @p ith
				      * component.
				      */
    Number operator() (const unsigned int i) const;
    
				     /**
				      * Access the @p ith component
				      * as a writeable reference.
				      */
    Number& operator() (const unsigned int i);
				     //@}


				     /**
				      * @name 3: Modification of vectors
				      */
				     //@{
    
				     /**
				      * Add the given vector to the present
				      * one.
				      */
    Vector<Number> & operator += (const Vector<Number> &V);

    				     /**
				      * Subtract the given vector from the
				      * present one.
				      */
    Vector<Number> & operator -= (const Vector<Number> &V);

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
    void scale (const Number factor);

    
				     /**
				      * Scale each element of the
				      * vector by a constant
				      * value.
				      */
    Vector<Number> & operator *= (const Number factor);

				     /**
				      * Scale each element of the
				      * vector by the inverse of the
				      * given value.
				      */
    Vector<Number> & operator /= (const Number factor);
    
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
				      * Assignment <tt>*this = a*V</tt>.
				      */
    template <typename Number2>
    void equ (const Number a, const Vector<Number2>& V);
    
				     /**
				      * Assignment <tt>*this = a*V + b*W</tt>.
				      */
    void equ (const Number a, const Vector<Number>& V,
	      const Number b, const Vector<Number>& W);

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
				      * @name 5: Mixed stuff
				      */
				     //@{
				     /**
				      *  Output of vector in user-defined format.
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
				      * Write the vector en bloc to a
				      * file. This is done in a binary
				      * mode, so the output is neither
				      * readable by humans nor
				      * (probably) by other computers
				      * using a different operating
				      * system or number format.
				      */
    void block_write (std::ostream &out) const;

				     /**
				      * Read a vector en block from a
				      * file. This is done using the
				      * inverse operations to the
				      * above function, so it is
				      * reasonably fast because the
				      * bitstream is not interpreted.
				      *
				      * The vector is resized if
				      * necessary.
				      *
				      * A primitive form of error
				      * checking is performed which
				      * will recognize the bluntest
				      * attempts to interpret some
				      * data as a vector stored
				      * bitwise to a file, but not
				      * more.
				      */
    void block_read (std::istream &in);

				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    unsigned int memory_consumption () const;
				     //@}

				     /**
				      * Exception
				      */
    DeclException0 (ExcEmptyVector);

  protected:

				     /**
				      * Dimension. Actual number of
				      * components contained in the
				      * vector.  Get this number by
				      * calling <tt>size()</tt>.
				      */
    unsigned int dim;

				     /**
				      * Amount of memory actually
				      * reserved for this vector. This
				      * number may be greater than
				      * @p dim if a @p reinit was
				      * called with less memory
				      * requirements than the vector
				      * needed last time. At present
				      * @p reinit does not free
				      * memory when the number of
				      * needed elements is reduced.
				      */
    unsigned int maxdim;

				     /**
				      * Pointer to the array of
				      * elements of this vector.
				      */
    Number *val;

                                     /**
                                      * Make all other vector types
                                      * friends.
                                      */
    template <typename Number2> friend class Vector;    
};

/*@}*/
/*----------------------- Inline functions ----------------------------------*/



template <typename Number>
inline
Vector<Number>::Vector ()
                :
		dim(0),
		maxdim(0),
		val(0)
{}



template <typename Number>
template <typename InputIterator>
Vector<Number>::Vector (const InputIterator first, const InputIterator last)
		:
		dim (0),
		maxdim (0),
		val (0)
{
				   // allocate memory. do not
				   // initialize it, as we will copy
				   // over to it in a second
  reinit (std::distance (first, last), true);
  std::copy (first, last, begin());
}



template <typename Number>
inline
Vector<Number>::Vector (const unsigned int n)
                :
		dim(0),
		maxdim(0),
		val(0)
{
  reinit (n, false);
}



template <typename Number>
inline
Vector<Number>::~Vector ()
{
  if (val)
    {
      delete[] val;
      val=0;
    }
}



template <typename Number>
inline
void Vector<Number>::reinit (const unsigned int n, const bool fast)
{
  if (n==0) 
    {
      if (val) delete[] val;
      val = 0;
      maxdim = dim = 0;
      return;
    };
  
  if (n>maxdim)
    {
      if (val) delete[] val;
      val = new value_type[n];
      Assert (val != 0, ExcOutOfMemory());
      maxdim = n;
    };
  dim = n;
  if (fast == false)
    *this = 0;
}



template <typename Number>
inline
unsigned int Vector<Number>::size () const
{
  return dim;
}



template <typename Number>
inline
typename Vector<Number>::iterator 
Vector<Number>::begin () 
{
  return &val[0];
}



template <typename Number>
inline
typename Vector<Number>::const_iterator 
Vector<Number>::begin () const 
{
  return &val[0];
}



template <typename Number>
inline
typename Vector<Number>::iterator
Vector<Number>::end () 
{
  return &val[dim];
}



template <typename Number>
inline
typename Vector<Number>::const_iterator
Vector<Number>::end () const
{
  return &val[dim];
}



template <typename Number>
inline
Number Vector<Number>::operator() (const unsigned int i) const
{
  Assert (i<dim, ExcIndexRange(i,0,dim));
  return val[i];
}



template <typename Number>
inline
Number& Vector<Number>::operator() (const unsigned int i)
{
  Assert (i<dim, ExcIndexRange(i,0,dim));
  return val[i];
}



template <typename Number>
inline
Vector<Number> & Vector<Number>::operator *= (const Number factor)
{
  scale (factor);
  return *this;
}



template <typename Number>
inline
Vector<Number> & Vector<Number>::operator /= (const Number factor)
{
  scale (1./factor);
  return *this;
}




template <typename Number>
inline
Vector<Number> & Vector<Number>::operator = (const Number s)
{
  Assert (dim!=0, ExcEmptyVector());
  std::fill (begin(), end(), s);
  return *this;
}



template <typename Number>
template <typename Number2>
inline
bool
Vector<Number>::operator != (const Vector<Number2>& v) const
{
  return ! (*this == v);
}



template <typename Number>
inline
void
Vector<Number>::compress () const
{}



/*! @addtogroup Vectors
 *@{
 */


/**
 * Global function @p swap which overloads the default implementation
 * of the C++ standard library which uses a temporary object. The
 * function simply exchanges the data of the two vectors.
 *
 * @relates Vector
 * @author Wolfgang Bangerth, 2000
 */
template <typename Number>
inline
void swap (Vector<Number> &u, Vector<Number> &v)
{
  u.swap (v);
}

/*@}*/

#endif

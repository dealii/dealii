//----------------------------  vector.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  vector.h  ---------------------------
#ifndef __deal2__vector_h
#define __deal2__vector_h


#include <base/exceptions.h>
#include <cstdio>


/**
 * Numerical vector of data.  For this class there are different types
 * of functions available. The first type of function mesures the norm
 * of the vector in order to mesure its length in a suitable norm. The
 * second type support the abgebraic operation for vectors. The third
 * und last type helps us to manipulate the components of the vector.
 * As opposed to the array of the C++ standard library called
 * @p{vector} (with a lowercase "v"), this class implements an element
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
				      * those in the @p{C++} standard libraries
				      * @p{vector<...>} class.
				      */
    typedef Number value_type;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef value_type* iterator;
    typedef const value_type* const_iterator;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef size_t size_type;


				     /**
				      * @name 1: Basic Object-handling 
				      */
				     //@{
				     /**
				      *  Dummy-Constructor. Dimension=0
				      */
    Vector ();
    
				     /**
				      * Copy-Constructor. Dimension set to
				      * that of V, all components are copied
				      * from V
				      */
    Vector (const Vector<Number>& V);


// note: I disabled this function for the time being, since egcs1.1.2
// does not respect the "explicit" keyword for template constructors.
// this leads to unwanted conversions and in some places to automatically
// generated temporaries, where this is not a good idea. [WB]
// 				     /**
// 				      * Copy constructor taking a vector of
// 				      * another data type. This will fail if
// 				      * there is no conversion path from
// 				      * @p{OtherNumber} to @p{Number}. Note that
// 				      * you may lose accuracy when copying
// 				      * to a vector with data elements with
// 				      * less accuracy.
// 				      */
//     template <typename OtherNumber>
//     explicit
//     Vector (const Vector<OtherNumber> &v);
    
				     /**
				      * Constructor. Set dimension to @p{n} and
				      * initialize all elements with zero.
				      */
    Vector (const unsigned int n);
    
				     /**
				      * Destructor, deallocates
				      * memory. Made virtual to allow
				      * for derived classes to behave
				      * properly.
				      */
    virtual ~Vector ();

				     /**
				      * Set all entries to zero. Equivalent to
				      * @p{v = 0}, but more obvious and faster.
				      * Note that this function does not change
				      * the size of the vector, unlike the
				      * STL's @p{vector<>::clear} function.
				      */
    void clear ();    

				     /**
				      * Change the dimension of the vector to
				      * @p{N}. The reserved memory for this vector
				      * remains unchanged if possible, to make
				      * things faster, but this may waste some
				      * memory, so take this in the back of your
				      * head.
				      * However, if @p{N==0} all memory is freed,
				      * i.e. if you want to resize the vector
				      * and release the memory not needed, you
				      * have to first call @p{reinit(0)} and then
				      * @p{reinit(N)}. This cited behaviour is
				      * analogous to that of the STL containers.
				      *
				      * On @p{fast==false}, the vector is filled by
				      * zeros.
				      */ 
    void reinit (const unsigned int N,
		 const bool         fast=false);
    
				     /**
				      * Change the dimension to that of the
				      * vector @p{V}. The same applies as for
				      * the other @p{reinit} function.
				      *
				      * The elements of @p{V} are not copied, i.e.
				      * this function is the same as calling
				      * @p{reinit (V.size(), fast)}.
				      */
    void reinit (const Vector<Number> &V,
		 const bool            fast=false);

				     /**
				      * Swap the contents of this
				      * vector and the other vector
				      * @p{v}. One could do this
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
				      * the @p{swap} function of all C++
				      * standard containers. Also,
				      * there is a global function
				      * @p{swap(u,v)} that simply calls
				      * @p{u.swap(v)}, again in analogy
				      * to standard functions.
				      */
    void swap (Vector<Number> &v);
    
				     /**
				      * $U(0-N) = s$: fill all components.
				      */
    Vector<Number> & operator= (const Number s);
    
				     /**
				      *  $U = V$: copy all components.
				      */
    Vector<Number> & operator= (const Vector<Number>& V);

				     /**
				      * $U = V$ for different types.
				      */
    template<typename Number2>
    Vector<Number> & operator= (const Vector<Number2>& V);
    
				     /**
				      * $U = U * V$: scalar product.
				      */
    Number operator* (const Vector<Number>& V) const;

				     /**
				      * Return square of the $l_2$-norm.
				      */
    Number norm_sqr () const;

				     /**
				      * Return the mean value of the elements of
				      * this vector.
				      */
    Number mean_value () const;

				     /**
				      * Return the $l_1$-norm of the vector, i.e.
				      * the sum of the absolute values.
				      */
    Number l1_norm () const;

				     /**
				      * Return the $l_2$-norm of the vector, i.e.
				      * the square root of the sum of the
				      * squares of the elements.
				      */
    Number l2_norm () const;

				     /**
				      * Return the maximum absolute value of the
				      * elements of this vector, which is the
				      * $l_\infty$-norm of a vector.
				      */
    Number linfty_norm () const;

				     /**
				      * Return dimension of the vector. This
				      * function was formerly called @p{n()}, but
				      * was renamed to get the @p{Vector} class
				      * closer to the C++ standard library's
				      * @p{vector} container.
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
				      * Make the @p{Vector} class a bit like the
				      * @p{vector<>} class of the C++ standard
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
				      * Access components, returns U(i).
				      */
    Number operator() (const unsigned int i) const;
    
				     /**
				      * Access components, returns U(i)
				      * as a writeable reference.
				      */
    Number& operator() (const unsigned int i);
				     //@}


				     /**
				      * @name 3: Modification of vectors
				      */
				     //@{
				     /**
				      * Addition operator.
				      * Fast equivalent to @p{U.add(1, V)}.
				      */
    Vector<Number> & operator += (const Vector<Number> &V);

    				     /**
				      * Subtraction operator.
				      * Fast equivalent to @p{U.add(-1, V)}.
				      */
    Vector<Number> & operator -= (const Vector<Number> &V);

				     /**
				      * $U(0-DIM)+=s$.
				      * Addition of @p{s} to all components. Note
				      * that @p{s} is a scalar and not a vector.
				      */
    void add (const Number s);
    
				     /**
				      * U+=V.
				      * Simple vector addition, equal to the
				      * @p{operator +=}.
				      */
    void add (const Vector<Number>& V);
    
				     /**
				      * U+=a*V.
				      * Simple addition of a scaled vector.
				      */
    void add (const Number a, const Vector<Number>& V);
    
				     /**
				      * U+=a*V+b*W.
				      * Multiple addition of scaled vectors.
				      */
    void add (const Number a, const Vector<Number>& V,
	      const Number b, const Vector<Number>& W);
    
				     /**
				      * U=s*U+V.
				      * Scaling and simple vector addition.
				      */
    void sadd (const Number s, const Vector<Number>& V);
    
				     /**
				      * U=s*U+a*V.
				      * Scaling and simple addition.
				      */
    void sadd (const Number s, const Number a, const Vector<Number>& V);
    
				     /**
				      * U=s*U+a*V+b*W.
				      * Scaling and multiple addition.
				      */
    void sadd (const Number s, const Number a,
	       const Vector<Number>& V, const Number b, const Vector<Number>& W);
    
				     /**
				      * U=s*U+a*V+b*W+c*X.
				      * Scaling and multiple addition.
				      */
    void sadd (const Number s, const Number a,
	       const Vector<Number>& V, const Number b, const Vector<Number>& W, 
	       const Number c, const Vector<Number>& X);
    
				     /**
				      * Scale each element of the vector by the
				      * given factor.
				      */
//TODO:[?] Why not have an operator *= instead of/in addition to `scale'?    
    void scale (const Number factor);
    
				     /**
				      *  U=a*V. Replacing.
				      */
    void equ (const Number a, const Vector<Number>& V);
    
				     /**
				      * U=a*V+b*W.
				      * Replacing by sum.
				      */
    void equ (const Number a, const Vector<Number>& V,
	      const Number b, const Vector<Number>& W);

				     /**
				      * Compute the elementwise ratio of the
				      * two given vectors, that is let
				      * @p{this[i] = a[i]/b[i]}. This is useful
				      * for example if you want to compute
				      * the cellwise ratio of true to estimated
				      * error.
				      *
				      * This vector is appropriately scaled to
				      * hold the result.
				      *
				      * If any of the @p{b[i]} is zero, the result
				      * is undefined. No attempt is made to
				      * catch such situations.
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
    void print (FILE* fp, const char* format = 0) const;
    
				     /**
				      *  Output of vector in user-defined format.
				      */
    void print (const char* format = 0) const;

				     /**
				      * Print to a
				      * stream. @p{precision} denotes
				      * the desired precision with
				      * which values shall be printed,
				      * @p{scientific} whether
				      * scientific notation shall be
				      * used. If @p{across} is
				      * @p{true} then the vector is
				      * printed in a line, while if
				      * @p{false} then the elements
				      * are printed on a separate line
				      * each.
				      */
    void print (std::ostream       &out,
		const unsigned int  precision  = 3,
		const bool          scientific = true,
		const bool          across     = true) const;

				     /**
				      * Write the vector en bloc to a file. This
				      * is done in a binary mode, so the output
				      * is neither readable by humans nor 
				      * (probably) by other computers using
				      * a different operating system of number
				      * format.
				      */
    void block_write (std::ostream &out) const;

				     /**
				      * Read a vector en block from a file. This
				      * is done using the inverse operations to
				      * the above function, so it is reasonably
				      * fast because the bitstream is not
				      * interpreted.
				      *
				      * The vector is resized if necessary.
				      *
				      * A primitive form of error checking is
				      * performed which will recognize the
				      * bluntest attempts to interpret some
				      * data as a vector stored bitwise to a
				      * file, but not more.
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
				      * Dimension. Actual number of components
				      * contained in the vector.
				      * Get this number by calling @p{size()}.
				      */
    unsigned int dim;

				     /**
				      * Amount of memory actually reserved for
				      * this vector. This number may be greater
				      * than @p{dim} if a @p{reinit} was called with
				      * less memory requirements than the vector
				      * needed last time. At present @p{reinit}
				      * does not free memory when the number of
				      * needed elements is reduced.
				      */
    unsigned int maxdim;

				     /**
				      * Pointer to the array of components.
				      */
    Number *val;
};


/*----------------------- Inline functions ----------------------------------*/


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
};


template <typename Number>
inline
typename Vector<Number>::const_iterator 
Vector<Number>::begin () const 
{
  return &val[0];
};


template <typename Number>
inline
typename Vector<Number>::iterator
Vector<Number>::end () 
{
  return &val[dim];
};


template <typename Number>
inline
typename Vector<Number>::const_iterator
Vector<Number>::end () const 
{
  return &val[dim];
};


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



/**
 * Global function @p{swap} which overloads the default implementation
 * of the C++ standard library which uses a temporary object. The
 * function simply exchanges the data of the two vectors.
 *
 * @author Wolfgang Bangerth, 2000
 */
template <typename Number>
inline
void swap (Vector<Number> &u, Vector<Number> &v)
{
  u.swap (v);
};



#endif

//----------------------------  block_vector.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  block_vector.h  ---------------------------
#ifndef __deal2__block_vector_h
#define __deal2__block_vector_h


#include <lac/vector.h>
#include <lac/block_indices.h>
#include <base/exceptions.h>
#include <cstdio>
#include <vector>


/**
 * Several vectors of data. 
 *
 * The BlockVector is a collection of normal LAC-#Vector#s. Each of
 * the vectors inside can have a different size.
 *
 * The functionality of #BlockVector# includes everything a #Vector#
 * can do, plus the access to a single #Vector# inside the
 * #BlockVector# by #block(i)#.
 *
 *
 * \section{On template instantiations}
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
 * @author Guido Kanschat, 1999; Wolfgang Bangerth, 2000
 */
template <int n_blocks, typename Number>
class BlockVector
{
  public:
				     /**
				      *  Dummy-Constructor. Dimension=0
				      */
    BlockVector ();
    
				     /**
				      * Copy-Constructor. Dimension set to
				      * that of V, all components are copied
				      * from V
				      */
    BlockVector (const BlockVector<n_blocks,Number>& V);


// note: I disabled this function for the time being, since egcs1.1.2
// does not respect the "explicit" keyword for template constructors.
// this leads to unwanted conversions and in some places to automatically
// generated temporaries, where this is not a good idea. [WB]
// 				     /**
// 				      * Copy constructor taking a BlockVector of
// 				      * another data type. This will fail if
// 				      * there is no conversion path from
// 				      * #OtherNumber# to #Number#. Note that
// 				      * you may lose accuracy when copying
// 				      * to a BlockVector with data elements with
// 				      * less accuracy.
// 				      */
//     template <typename OtherNumber>
//     explicit
//     BlockVector (const BlockVector<OtherNumber> &v);
    
				     /**
				      * Constructor. Set dimension to #n# and
				      * initialize all elements with zero.
				      */
    BlockVector (const vector<unsigned int> &n);

                                     /**
				      * Destructor. Clears memory
				      */
    ~BlockVector ();

				     /**
				      * Change the dimension of the vector to
				      * #N#. The reserved memory for this vector
				      * remains unchanged if possible, to make
				      * things faster, but this may waste some
				      * memory, so take this in the back of your
				      * head.
				      * However, if #N==0# all memory is freed,
				      * i.e. if you want to resize the vector
				      * and release the memory not needed, you
				      * have to first call #reinit(0)# and then
				      * #reinit(N)#. This cited behaviour is
				      * analogous to that of the STL containers.
				      *
				      * On #fast==false#, the vector is filled by
				      * zeros.
				      */ 
    void reinit (const vector<unsigned int> &N,
		 const bool                  fast=false);
    
				     /**
				      * Change the dimension to that of the
				      * vector #V#. The same applies as for
				      * the other #reinit# function.
				      *
				      * The elements of #V# are not copied, i.e.
				      * this function is the same as calling
				      * #reinit (V.size(), fast)#.
				      */
    void reinit (const BlockVector<n_blocks,Number> &V,
		 const bool                          fast=false);
    
				     /**
				      * Set all entries to zero. Equivalent to
				      * #v = 0#, but more obvious and faster.
				      * Note that this function does not change
				      * the size of the vector, unlike the
				      * STL's #vector<>::clear# function.
				      */
    void clear ();
    
				     /**
				      * Swap the contents of this
				      * vector and the other vector
				      * #v#. One could do this
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
				      * the #swap# function of all C++
				      * standard containers. Also,
				      * there is a global function
				      * #swap(u,v)# that simply calls
				      * #u.swap(v)#, again in analogy
				      * to standard functions.
				      */
    void swap (BlockVector<n_blocks,Number> &v);
    
				     /**
				      * Access to a single block.
				      */
    Vector<Number>& block (const unsigned int i);
    
				     /**
				      * Read-only access to a single block.
				      */
    const Vector<Number>& block (const unsigned int i) const;
    
				     /**
				      * $U(0-N) = s$: fill all components.
				      */
    BlockVector<n_blocks,Number> & operator= (const Number s);
    
				     /**
				      *  $U = V$: copy all components.
				      */
    BlockVector<n_blocks,Number> &
    operator= (const BlockVector<n_blocks,Number>& V);

				     /**
				      * $U = V$ for different types.
				      */
    template<typename Number2>
    BlockVector<n_blocks,Number> &
    operator= (const BlockVector<n_blocks, Number2>& V);
    
				     /**
				      * $U = U * V$: scalar product.
				      */
    Number operator* (const BlockVector<n_blocks,Number>& V) const;

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
  				      * Return dimension of the vector. This is the
				      * sum of the dimensions of all components.
  				      */
    unsigned int size () const;

				     /**
				      * Return whether the vector contains only
				      * elements with value zero. This function
				      * is mainly for internal consistency
				      * check and should seldomly be used when
				      * not in debug mode since it uses quite
				      * some time.
				      */
    bool all_zero () const;

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
				      * Fast equivalent to #U.add(1, V)#.
				      */
    BlockVector<n_blocks,Number> &
    operator += (const BlockVector<n_blocks,Number> &V);

    				     /**
				      * Subtraction operator.
				      * Fast equivalent to #U.add(-1, V)#.
				      */
    BlockVector<n_blocks,Number> &
    operator -= (const BlockVector<n_blocks,Number> &V);

				     /**
				      * $U(0-DIM)+=s$.
				      * Addition of #s# to all components. Note
				      * that #s# is a scalar and not a vector.
				      */
    void add (const Number s);
    
				     /**
				      * U+=V.
				      * Simple vector addition, equal to the
				      * #operator +=#.
				      */
    void add (const BlockVector<n_blocks,Number>& V);
    
				     /**
				      * U+=a*V.
				      * Simple addition of a scaled vector.
				      */
    void add (const Number a, const BlockVector<n_blocks,Number>& V);
    
				     /**
				      * U+=a*V+b*W.
				      * Multiple addition of scaled vectors.
				      */
    void add (const Number a, const BlockVector<n_blocks,Number>& V,
	      const Number b, const BlockVector<n_blocks,Number>& W);
    
				     /**
				      * U=s*U+V.
				      * Scaling and simple vector addition.
				      */
    void sadd (const Number s, const BlockVector<n_blocks,Number>& V);
    
				     /**
				      * U=s*U+a*V.
				      * Scaling and simple addition.
				      */
    void sadd (const Number s, const Number a, const BlockVector<n_blocks,Number>& V);
    
				     /**
				      * U=s*U+a*V+b*W.
				      * Scaling and multiple addition.
				      */
    void sadd (const Number s, const Number a,
	       const BlockVector<n_blocks,Number>& V,
	       const Number b, const BlockVector<n_blocks,Number>& W);
    
				     /**
				      * U=s*U+a*V+b*W+c*X.
				      * Scaling and multiple addition.
				      */
    void sadd (const Number s, const Number a,
	       const BlockVector<n_blocks,Number>& V,
	       const Number b, const BlockVector<n_blocks,Number>& W, 
	       const Number c, const BlockVector<n_blocks,Number>& X);
    
				     /**
				      * Scale each element of the vector by the
				      * given factor. This function was
				      * previously called #equ(Number)#, which
				      * in my eyes is an extremely unintuitive
				      * naming and was thus replaced.
				      */
    void scale (const Number factor);
    
				     /**
				      *  U=a*V. Replacing.
				      */
    void equ (const Number a, const BlockVector<n_blocks,Number>& V);
    
				     /**
				      * U=a*V+b*W.
				      * Replacing by sum.
				      */
    void equ (const Number a, const BlockVector<n_blocks,Number>& V,
	      const Number b, const BlockVector<n_blocks,Number>& W);

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
				      * Print to a stream.
				      * 
				      */
    void print (ostream &, unsigned int precision = 3,
		bool scientific = true,
		bool across = true) const;

				     /**
				      * Write the vector en bloc to a file. This
				      * is done in a binary mode, so the output
				      * is neither readable by humans nor 
				      * (probably) by other computers using
				      * a different operating system of number
				      * format.
				      */
    void block_write (ostream &out) const;

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
    void block_read (istream &in);
				     //@}

				     /**
				      * Exception
				      */
    DeclException2 (ExcDimensionsDontMatch,
		    int, int,
		    << "The dimensions " << arg1 << " and " << arg2
		    << " do not match here.");
				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidIndex,
		    int, int,
		    << "The given index " << arg1
		    << " should be less than " << arg2 << ".");
				     /**
				      * Exception
				      */
    DeclException1 (ExcInvalidNumber,
		    int,
		    << "The provided number is invalid here: " << arg1);
				     /**
				      * Exception
				      */
    DeclException0 (ExcOutOfMemory);
				     /**
				      * Exception
				      */
    DeclException0 (ExcEmptyVector);
				     /**
				      * Exception
				      */
    DeclException0 (ExcIO);

  protected:
				     /**
				      * Pointer to the array of components.
				      */
    Vector<Number> components[n_blocks];

				     /**
				      * Object managing the
				      * transformation between global
				      * indices and indices within the
				      * different blocks.
				      */
    BlockIndices<n_blocks> block_indices;
};


/*----------------------- Inline functions ----------------------------------*/


template <int n_blocks, typename Number>
inline
unsigned int BlockVector<n_blocks,Number>::size () const
{
  return block_indices.total_size();
}


template <int n_blocks, typename Number>
inline
Number BlockVector<n_blocks,Number>::operator() (const unsigned int i) const
{
  const pair<unsigned int,unsigned int> local_index
    = block_indices.global_to_local (i);
  return components[local_index.first](local_index.second);
}


template <int n_blocks, typename Number>
inline
Number& BlockVector<n_blocks,Number>::operator() (const unsigned int i)
{
  const pair<unsigned int,unsigned int> local_index
    = block_indices.global_to_local (i);
  return components[local_index.first](local_index.second);
}



template <int n_blocks, typename Number>
inline
Vector<Number>&
BlockVector<n_blocks,Number>::block(unsigned int i)
{
  Assert(i<n_blocks, ExcIndexRange(i,0,n_blocks));

  return components[i];
}

template <int n_blocks, typename Number>
inline
const Vector<Number>&
BlockVector<n_blocks,Number>::block(unsigned int i) const
{
  Assert(i<n_blocks, ExcIndexRange(i,0,n_blocks));

  return components[i];
}


/**
 * Global function #swap# which overloads the default implementation
 * of the C++ standard library which uses a temporary object. The
 * function simply exchanges the data of the two vectors.
 *
 * @author Wolfgang Bangerth, 2000
 */
template <int n_blocks, typename Number>
inline
void swap (BlockVector<n_blocks,Number> &u,
	   BlockVector<n_blocks,Number> &v)
{
  u.swap (v);
};




#endif

//----------------------------  block_vector.h  ---------------------------
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
//----------------------------  block_vector.h  ---------------------------
#ifndef __deal2__block_vector_h
#define __deal2__block_vector_h


#include <base/config.h>
#include <base/exceptions.h>
#include <lac/vector.h>
#include <lac/block_indices.h>
#include <lac/block_vector_base.h>

#include <cstdio>
#include <vector>



/*! @addtogroup Vectors
 *@{
 */


/**
 * An implementation of block vectors based on deal.II vectors. While the base
 * class provides for most of the interface, this class handles the actual
 * allocation of vectors and provides functions that are specific to the
 * underlying vector type.
 *
 * @ref Instantiations: some (<tt>@<float@> @<double@></tt>)
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999, 2000, 2001, 2002, 2004
 */
template <typename Number>
class BlockVector : public BlockVectorBase<Vector<Number> >
{
  public:
                                     /**
                                      * Typedef the base class for simpler
                                      * access to its own typedefs.
                                      */
    typedef BlockVectorBase<Vector<Number> > BaseClass;
    
                                     /**
                                      * Typedef the type of the underlying
                                      * vector.
                                      */
    typedef typename BaseClass::BlockType  BlockType;

				     /**
				      * Import the typedefs from the base
				      * class.
				      */
    typedef typename BaseClass::value_type      value_type;
    typedef typename BaseClass::pointer         pointer;
    typedef typename BaseClass::const_pointer   const_pointer;
    typedef typename BaseClass::reference       reference;
    typedef typename BaseClass::const_reference const_reference;
    typedef typename BaseClass::size_type       size_type;
    typedef typename BaseClass::iterator        iterator;
    typedef typename BaseClass::const_iterator  const_iterator;

				     /**
				      *  Constructor. There are three
				      *  ways to use this
				      *  constructor. First, without
				      *  any arguments, it generates
				      *  an object with no
				      *  blocks. Given one argument,
				      *  it initializes <tt>num_blocks</tt>
				      *  blocks, but these blocks have
				      *  size zero. The third variant
				      *  finally initializes all
				      *  blocks to the same size
				      *  <tt>block_size</tt>.
				      *
				      *  Confer the other constructor
				      *  further down if you intend to
				      *  use blocks of different
				      *  sizes.
				      */
    explicit BlockVector (const unsigned int num_blocks = 0,
			  const unsigned int block_size = 0);
    
				     /**
				      * Copy-Constructor. Dimension set to
				      * that of V, all components are copied
				      * from V
				      */
    BlockVector (const BlockVector<Number>& V);


#ifndef DEAL_II_EXPLICIT_CONSTRUCTOR_BUG    
				     /**
				      * Copy constructor taking a BlockVector of
				      * another data type. This will fail if
				      * there is no conversion path from
				      * <tt>OtherNumber</tt> to <tt>Number</tt>. Note that
				      * you may lose accuracy when copying
				      * to a BlockVector with data elements with
				      * less accuracy.
				      *
				      * Older versions of gcc did not honor
				      * the @p explicit keyword on template
				      * constructors. In such cases, it is
				      * easy to accidentally write code that
				      * can be very inefficient, since the
				      * compiler starts performing hidden
				      * conversions. To avoid this, this
				      * function is disabled if we have
				      * detected a broken compiler during
				      * configuration.
				      */
    template <typename OtherNumber>
    explicit
    BlockVector (const BlockVector<OtherNumber> &v);
#endif    

                                     /**
                                      * Constructor. Set the number of
                                      * blocks to
                                      * <tt>block_sizes.size()</tt> and
                                      * initialize each block with
                                      * <tt>block_sizes[i]</tt> zero
                                      * elements.
                                      */
    BlockVector (const std::vector<unsigned int> &block_sizes);

				     /**
				      * Constructor. Set the number of
				      * blocks to
				      * <tt>n.size()</tt>. Initialize the
				      * vector with the elements
				      * pointed to by the range of
				      * iterators given as second and
				      * third argument. Apart from the
				      * first argument, this
				      * constructor is in complete
				      * analogy to the respective
				      * constructor of the
				      * <tt>std::vector</tt> class, but the
				      * first argument is needed in
				      * order to know how to subdivide
				      * the block vector into
				      * different blocks.
				      */
    template <typename InputIterator>
    BlockVector (const std::vector<unsigned int> &n,
		 const InputIterator              first,
		 const InputIterator              end);
    
                                     /**
				      * Destructor. Clears memory
				      */
    ~BlockVector ();

				     /**
				      * Copy operator: fill all components of
				      * the vector with the given scalar
				      * value.
				      */
    BlockVector & operator = (const value_type s);

				     /**
				      * Copy operator for arguments of the
				      * same type.
				      */
    BlockVector &
    operator= (const BlockVector &V);

				     /**
				      * Copy operator for template arguments
				      * of different types.
				      */
    template <class Number2>
    BlockVector &
    operator= (const BlockVector<Number2> &V);

    				     /**
				      * Reinitialize the BlockVector to
				      * contain <tt>num_blocks</tt> blocks of
				      * size <tt>block_size</tt> each.
				      *
				      * If <tt>fast==false</tt>, the vector
				      * is filled with zeros.
				      */
    void reinit (const unsigned int num_blocks,
		 const unsigned int block_size,
		 const bool fast = false);
  
				     /**
				      * Reinitialize the BlockVector such that
				      * it contains
				      * <tt>block_sizes.size()</tt>
				      * blocks. Each block is reinitialized to
				      * dimension <tt>block_sizes[i]</tt>.
				      *
				      * If the number of blocks is the
				      * same as before this function
				      * was called, all vectors remain
				      * the same and reinit() is
				      * called for each vector.
				      *
				      * If <tt>fast==false</tt>, the vector
				      * is filled with zeros.
				      *
				      * Note that you must call this
				      * (or the other reinit()
				      * functions) function, rather
				      * than calling the reinit()
				      * functions of an individual
				      * block, to allow the block
				      * vector to update its caches of
				      * vector sizes. If you call
				      * reinit() on one of the
				      * blocks, then subsequent
				      * actions on this object may
				      * yield unpredictable results
				      * since they may be routed to
				      * the wrong block.
				      */ 
    void reinit (const std::vector<unsigned int> &N,
		 const bool                       fast=false);
    
				     /**
				      * Change the dimension to that
				      * of the vector <tt>V</tt>. The same
				      * applies as for the other
				      * reinit() function.
				      *
				      * The elements of <tt>V</tt> are not
				      * copied, i.e.  this function is
				      * the same as calling <tt>reinit
				      * (V.size(), fast)</tt>.
				      *
				      * Note that you must call this
				      * (or the other reinit()
				      * functions) function, rather
				      * than calling the reinit()
				      * functions of an individual
				      * block, to allow the block
				      * vector to update its caches of
				      * vector sizes. If you call
				      * reinit() of one of the
				      * blocks, then subsequent
				      * actions of this object may
				      * yield unpredictable results
				      * since they may be routed to
				      * the wrong block.
				      */
    template <typename Number2>
    void reinit (const BlockVector<Number2> &V,
		 const bool                 fast=false);

				     /**
				      * Scale each element of the
				      * vector by the given factor.
				      *
				      * This function is deprecated
				      * and will be removed in a
				      * future version. Use
				      * <tt>operator *=</tt> and
				      * <tt>operator /=</tt> instead.
				      *
				      * @deprecated Use <tt>operator*=</tt>
				      * instead.
				      */
    void scale (const value_type factor);
    
				     /**
				      * Multiply each element of this
				      * vector by the corresponding
				      * element of <tt>v</tt>.
				      */
    template <class BlockVector2>
    void scale (const BlockVector2 &v);
    
				     /**
				      * Swap the contents of this
				      * vector and the other vector
				      * <tt>v</tt>. One could do this
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
				      * Limitation: right now this
				      * function only works if both
				      * vectors have the same number
				      * of blocks. If needed, the
				      * numbers of blocks should be
				      * exchanged, too.
				      *
				      * This function is analog to the
				      * the swap() function of all C++
				      * standard containers. Also,
				      * there is a global function
				      * swap(u,v) that simply calls
				      * <tt>u.swap(v)</tt>, again in analogy
				      * to standard functions.
				      */
    void swap (BlockVector<Number> &v);    

				     /**
				      *  Output of vector in user-defined
				      *  format.
				      */
    void print (const char* format = 0) const;

				     /**
				      * Print to a stream.
				      */
    void print (std::ostream       &out,
		const unsigned int  precision = 3,
		const bool          scientific = true,
		const bool          across = true) const;

				     /**
				      * Write the vector en bloc to a
				      * stream. This is done in a binary mode,
				      * so the output is neither readable by
				      * humans nor (probably) by other
				      * computers using a different operating
				      * system or number format.
				      */
    void block_write (std::ostream &out) const;

				     /**
				      * Read a vector en block from a
				      * file. This is done using the inverse
				      * operations to the above function, so
				      * it is reasonably fast because the
				      * bitstream is not interpreted.
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

      				     /** @addtogroup Exceptions
				      * @{ */
    
				     /**
				      * Exception
				      */
    DeclException0 (ExcIteratorRangeDoesNotMatchVectorSize);
				     //@}
};

/*@}*/

/*----------------------- Inline functions ----------------------------------*/


template <typename Number>
template <typename InputIterator>
BlockVector<Number>::BlockVector (const std::vector<unsigned int> &n,
				  const InputIterator              first,
				  const InputIterator              end)
{
				   // first set sizes of blocks, but
				   // don't initialize them as we will
				   // copy elements soon
  reinit (n, true);
  InputIterator start = first;
  for (unsigned int b=0; b<n.size(); ++b)
    {
      InputIterator end = start;
      std::advance (end, static_cast<signed int>(n[b]));
      std::copy (start, end, this->block(b).begin());
      start = end;
    };
  Assert (start == end, ExcIteratorRangeDoesNotMatchVectorSize());
}



template <typename Number>
inline
BlockVector<Number> &
BlockVector<Number>::operator = (const value_type s)
{
  BaseClass::operator = (s);
  return *this;
}



template <typename Number>
inline
BlockVector<Number> &
BlockVector<Number>::operator = (const BlockVector &v)
{
  BaseClass::operator = (v);
  return *this;
}



template <typename Number>
template <typename Number2>
inline
BlockVector<Number> &
BlockVector<Number>::operator = (const BlockVector<Number2> &v)
{
  BaseClass::operator = (v);
  return *this;
}


template <typename Number>
void BlockVector<Number>::scale (const value_type factor)
{
  for (unsigned int i=0; i<this->n_blocks();++i)
    this->components[i].scale(factor);
}



template <typename Number>
template <class BlockVector2>
void BlockVector<Number>::scale (const BlockVector2 &v)
{
  BaseClass::scale (v);
}




/*! @addtogroup Vectors
 *@{
 */

/**
 * Global function which overloads the default implementation
 * of the C++ standard library which uses a temporary object. The
 * function simply exchanges the data of the two vectors.
 *
 * @relates BlockVector
 * @author Wolfgang Bangerth, 2000
 */
template <typename Number>
inline
void swap (BlockVector<Number> &u,
	   BlockVector<Number> &v)
{
  u.swap (v);
}


/*@}*/

#endif

//----------------------------  block_vector.h  ---------------------------
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
 * The BlockVector is a collection of normal LAC-@p{Vector}s. Each of
 * the vectors inside can have a different size. The special case of a
 * block vector with constant block size is supported by constructor
 * and reinit functions.
 *
 * The functionality of @p{BlockVector} includes everything a @p{Vector}
 * can do, plus the access to a single @p{Vector} inside the
 * @p{BlockVector} by @p{block(i)}.
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
 * @author Wolfgang Bangerth, Guido Kanschat, 1999, 2000
 */
template <typename Number>
class BlockVector
{
    template <typename Pointee>
    class Iterator
    {
      public:
					 /**
					  * Construct an iterator from
					  * a vector to which we point
					  * and the global index of
					  * the element pointed to.
					  */
	Iterator (BlockVector<Number> &parent,
		  const unsigned       global_index);
	
					 /**
					  * Copy constructor.
					  */
	Iterator (const Iterator &c);

					 /**
					  * Copy operator.
					  */
	Iterator & operator = (const Iterator &c);

					 /**
					  * Dereferencing operator. If
					  * the template argument
					  * @p{Pointee} has a
					  * @p{const} specification,
					  * then no writing to the
					  * result is possible, making
					  * this a @p{const_iterator}.
					  */
	Pointee & operator * () const;

					 /**
					  * Prefix @p{++} operator:
					  * @p{++i}. This operator
					  * advances the iterator to
					  * the next element and
					  * returns a reference to
					  * @p{*this}.
					  */
	Iterator & operator ++ ();

					 /**
					  * Postfix @p{++} operator:
					  * @p{i++}. This operator
					  * advances the iterator to
					  * the next element and
					  * returns a copy of the old
					  * value of this iterator.
					  */
	Iterator operator ++ (int);

					 /**
					  * Prefix @p{--} operator:
					  * @p{--i}. This operator
					  * retracts the iterator to
					  * the previous element and
					  * returns a reference to
					  * @p{*this}.
					  */
	Iterator & operator -- ();

					 /**
					  * Postfix @p{--} operator:
					  * @p{i--}. This operator
					  * retracts the iterator to
					  * the previous element and
					  * returns a copy of the old
					  * value of this iterator.
					  */
	Iterator operator -- (int);

					 /**
					  * Compare for equality of
					  * iterators. This operator
					  * checks whether the vectors
					  * pointed to are the same,
					  * and if not it throws an
					  * exception.					 
					  */
	bool operator == (const Iterator &i) const;

					 /**
					  * Compare for inequality of
					  * iterators. This operator
					  * checks whether the vectors
					  * pointed to are the same,
					  * and if not it throws an
					  * exception.					 
					  */
	bool operator != (const Iterator &i) const;

					 /**
					  * Check whether this
					  * iterators points to an
					  * element previous to the
					  * one pointed to by the
					  * given argument. This
					  * operator checks whether
					  * the vectors pointed to are
					  * the same, and if not it
					  * throws an exception.
					  */
	bool operator < (const Iterator &i) const;

					 /**
					  * Comparison operator alike
					  * to the one above.
					  */
	bool operator <= (const Iterator &i) const;

					 /**
					  * Comparison operator alike
					  * to the one above.
					  */
	bool operator > (const Iterator &i) const;

					 /**
					  * Comparison operator alike
					  * to the one above.
					  */
	bool operator >= (const Iterator &i) const;
	
					 /**
					  * Exception.
					  */
	DeclException0 (ExcPointerToDifferentVectors);

	
      private:
					 /**
					  * Pointer to the block
					  * vector object to which
					  * this iterator points.
					  */
	BlockVector<Number> *parent;

					 /**
					  * Global index of the
					  * element to which we
					  * presently point.
					  */
	unsigned int         global_index;

					 /**
					  * Current block and index
					  * within this block of the
					  * element presently pointed
					  * to.
					  */
	unsigned int current_block;
	unsigned int index_within_block;

					 /**
					  * Indices of the global
					  * element address at which
					  * we have to move on to
					  * another block when moving
					  * forward and
					  * backward. These indices
					  * are kept as a cache since
					  * this is much more
					  * efficient than always
					  * asking the parent object.
					  */
	unsigned int next_break_forward;
	unsigned int next_break_backward;

					 /**
					  * Move forward one element.
					  */
	void move_forward ();

					 /**
					  * Move backward one element.
					  */
	void move_backward ();
    };
    
  public:
				     /**
				      * Declare standard types used in
				      * all containers. These types
				      * parallel those in the @p{C++}
				      * standard libraries
				      * @p{vector<...>} class. The
				      * @p{iterator} types are not
				      * declared at present, since
				      * there are no iterators
				      * implemented that cycle through
				      * the individual sub-vectors.
				      */
    typedef Number                  value_type;
    typedef value_type             *pointer;
    typedef const value_type       *const_pointer;
    typedef Iterator<Number>        iterator;
    typedef Iterator<const Number>  const_iterator;
    typedef value_type             &reference;
    typedef const value_type       &const_reference;
    typedef size_t                  size_type;

				     /**
				      *  Constructor. There are three
				      *  ways to use this
				      *  constructor. First, without
				      *  any arguments, it generates
				      *  an objetct with no
				      *  blocks. Given one argument,
				      *  it initializes @p{num_blocks}
				      *  blocks, but these blocks have
				      *  size zero. The third variant
				      *  finally initializes all
				      *  blocks to the same size
				      *  @p{block_size}.
				      *
				      *  Confer the other constructor
				      *  further down if you intend to
				      *  use blocks of different
				      *  sizes.
				      */
    BlockVector (unsigned int num_blocks = 0,
		 unsigned int block_size = 0);
    
				     /**
				      * Copy-Constructor. Dimension set to
				      * that of V, all components are copied
				      * from V
				      */
    BlockVector (const BlockVector<Number>& V);


// note: I disabled this function for the time being, since egcs1.1.2
// does not respect the "explicit" keyword for template constructors.
// this leads to unwanted conversions and in some places to automatically
// generated temporaries, where this is not a good idea. [WB]
// 				     /**
// 				      * Copy constructor taking a BlockVector of
// 				      * another data type. This will fail if
// 				      * there is no conversion path from
// 				      * @p{OtherNumber} to @p{Number}. Note that
// 				      * you may lose accuracy when copying
// 				      * to a BlockVector with data elements with
// 				      * less accuracy.
// 				      */
//     template <typename OtherNumber>
//     explicit
//     BlockVector (const BlockVector<OtherNumber> &v);
    
				     /**
				      * Constructor. Set the number of
				      * blocks to @p{n.size()} and
				      * initialize each block with
				      * @p{n[i]} zero elements.
				      */
    BlockVector (const std::vector<unsigned int> &n);

				     /**
				      * Constructor. Set the number of
				      * blocks to
				      * @p{n.size()}. Initialize the
				      * vector with the elements
				      * pointed to by the range of
				      * iterators given as second and
				      * third argument. Apart from the
				      * first argument, this
				      * constructor is in complete
				      * analogy to the respective
				      * constructor of the
				      * @p{std::vector} class, but the
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
				      * Reinitialize the BlockVector to
				      * contain @p{num_blocks} blocks of
				      * size @p{block_size} each.
				      *
				      * If @p{fast==false}, the vector
				      * is filled with zeros.
				      */
    void reinit (const unsigned int num_blocks,
		 const unsigned int block_size,
		 const bool fast = false);
  
				     /**
				      * Reinitialize the BlockVector
				      * such that it contains
				      * @p{N.size()} blocks. Each
				      * Block is reinitialized to
				      * dimension @p{N[i]}.
				      *
				      * If the number of blocks is the
				      * same as before this function
				      * was called, all vectors remain
				      * the same and @p{reinit} is
				      * called for each vector. While
				      * reinitailizing a usual vector
				      * can consume a lot of time,
				      * this function here definitely
				      * has a potential to slow down a
				      * program considerably.
				      *
				      * If @p{fast==false}, the vector
				      * is filled with zeros.
				      */ 
    void reinit (const std::vector<unsigned int> &N,
		 const bool                       fast=false);
    
				     /**
				      * Change the dimension to that of the
				      * vector @p{V}. The same applies as for
				      * the other @p{reinit} function.
				      *
				      * The elements of @p{V} are not copied, i.e.
				      * this function is the same as calling
				      * @p{reinit (V.size(), fast)}.
				      */
    void reinit (const BlockVector<Number> &V,
		 const bool                 fast=false);
    
				     /**
				      * Set all entries to zero. Equivalent to
				      * @p{v = 0}, but more obvious and faster.
				      * Note that this function does not change
				      * the size of the vector, unlike the
				      * STL's @p{vector<>::clear} function.
				      */
    void clear ();
    
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
				      * Limitation: right now this
				      * function only works if both
				      * vectors have the same number
				      * of blocks. If needed, the
				      * numbers of blocks should be
				      * exchanged, too.
				      *
				      * This function is analog to the
				      * the @p{swap} function of all C++
				      * standard containers. Also,
				      * there is a global function
				      * @p{swap(u,v)} that simply calls
				      * @p{u.swap(v)}, again in analogy
				      * to standard functions.
				      */
    void swap (BlockVector<Number> &v);
    
				     /**
				      * Access to a single block.
				      */
    Vector<Number>& block (const unsigned int i);
    
				     /**
				      * Read-only access to a single block.
				      */
    const Vector<Number> &
    block (const unsigned int i) const;

				     /**
				      * Return a reference on the
				      * object that describes the
				      * mapping between block and
				      * global indices. The use of
				      * this function is highly
				      * deprecated and it should
				      * vanish in one of the next
				      * versions
				      */
    const BlockIndices &
    get_block_indices () const;
    
				     /**
				      * $U(0-N) = s$: fill all components.
				      */
    BlockVector<Number> & operator= (const Number s);
    
				     /**
				      *  $U = V$: copy all components.
				      */
    BlockVector<Number> &
    operator= (const BlockVector<Number>& V);

				     /**
				      * $U = V$ for different types.
				      */
    template<typename Number2>
    BlockVector<Number> &
    operator= (const BlockVector< Number2>& V);
    
				     /**
				      * $U = U * V$: scalar product.
				      */
    Number operator* (const BlockVector<Number>& V) const;

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
				      * Number of blocks.
				      */
    unsigned int n_blocks () const;
  
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

				     /**
				      * Return an iterator pointing to
				      * the first element.
				      */
    iterator begin ();

				     /**
				      * Return an iterator pointing to
				      * the first element of a
				      * constant block vector.
				      */
    const_iterator begin () const;
    
				     //@}


				     /**
				      * @name 3: Modification of vectors
				      */
				     //@{
				     /**
				      * Addition operator.
				      * Fast equivalent to @p{U.add(1, V)}.
				      */
    BlockVector<Number> &
    operator += (const BlockVector<Number> &V);

    				     /**
				      * Subtraction operator.
				      * Fast equivalent to @p{U.add(-1, V)}.
				      */
    BlockVector<Number> &
    operator -= (const BlockVector<Number> &V);

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
    void add (const BlockVector<Number>& V);
    
				     /**
				      * U+=a*V.
				      * Simple addition of a scaled vector.
				      */
    void add (const Number a, const BlockVector<Number>& V);
    
				     /**
				      * U+=a*V+b*W.
				      * Multiple addition of scaled vectors.
				      */
    void add (const Number a, const BlockVector<Number>& V,
	      const Number b, const BlockVector<Number>& W);
    
				     /**
				      * U=s*U+V.
				      * Scaling and simple vector addition.
				      */
    void sadd (const Number s, const BlockVector<Number>& V);
    
				     /**
				      * U=s*U+a*V.
				      * Scaling and simple addition.
				      */
    void sadd (const Number s, const Number a, const BlockVector<Number>& V);
    
				     /**
				      * U=s*U+a*V+b*W.
				      * Scaling and multiple addition.
				      */
    void sadd (const Number s, const Number a,
	       const BlockVector<Number>& V,
	       const Number b, const BlockVector<Number>& W);
    
				     /**
				      * U=s*U+a*V+b*W+c*X.
				      * Scaling and multiple addition.
				      */
    void sadd (const Number s, const Number a,
	       const BlockVector<Number>& V,
	       const Number b, const BlockVector<Number>& W, 
	       const Number c, const BlockVector<Number>& X);
    
				     /**
				      * Scale each element of the vector by the
				      * given factor. This function was
				      * previously called @p{equ(Number)}, which
				      * in my eyes is an extremely unintuitive
				      * naming and was thus replaced.
				      */
    void scale (const Number factor);
    
				     /**
				      * Scale each element of the
				      * vector by a constant
				      * value. This operator is an
				      * alias to the @ref{scale}
				      * function, except that it
				      * returns a reference to itself.
				      */
    BlockVector<Number> & operator *= (const Number factor);

				     /**
				      *  U=a*V. Replacing.
				      */
    void equ (const Number a, const BlockVector<Number>& V);
    
				     /**
				      * U=a*V+b*W.
				      * Replacing by sum.
				      */
    void equ (const Number a, const BlockVector<Number>& V,
	      const Number b, const BlockVector<Number>& W);

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
				      * Print to a stream.
				      */
    void print (std::ostream       &out,
		const unsigned int  precision = 3,
		const bool          scientific = true,
		const bool          across = true) const;

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
    DeclException0 (ExcIteratorRangeDoesNotMatchVectorSize);

  protected:
				     /**
				      * Pointer to the array of components.
				      */
    typename std::vector<Vector<Number> > components;

				     /**
				      * Object managing the
				      * transformation between global
				      * indices and indices within the
				      * different blocks.
				      */
    BlockIndices block_indices;
    
  private:
				     /**
				      * The number of blocks. This
				      * number is redundant to
				      * @p{components.size()} and stored
				      * here for convenience.
				      */
    unsigned int num_blocks;

				     /**
				      * Make the iterator class a
				      * friend.
				      */
    template <typename Pointee> friend class Iterator;
};


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
      std::copy (start, end, block(b).begin());
      start = end;
    };
  Assert (start == end, ExcIteratorRangeDoesNotMatchVectorSize());
};




template <typename Number>
inline
unsigned int BlockVector<Number>::size () const
{
  return block_indices.total_size();
}


template <typename Number>
inline
unsigned int BlockVector<Number>::n_blocks () const
{
  return num_blocks;
}



template <typename Number>
inline
Number BlockVector<Number>::operator() (const unsigned int i) const
{
  const std::pair<unsigned int,unsigned int> local_index
    = block_indices.global_to_local (i);
  return components[local_index.first](local_index.second);
}



template <typename Number>
inline
Number& BlockVector<Number>::operator() (const unsigned int i)
{
  const std::pair<unsigned int,unsigned int> local_index
    = block_indices.global_to_local (i);
  return components[local_index.first](local_index.second);
}



template <typename Number>
inline
BlockVector<Number> & BlockVector<Number>::operator *= (const Number factor) 
{
  scale (factor);
  return *this;
};



template <typename Number>
inline
Vector<Number> &
BlockVector<Number>::block (const unsigned int i)
{
  Assert(i<num_blocks, ExcIndexRange(i,0,num_blocks));

  return components[i];
}



template <typename Number>
inline
const Vector<Number> &
BlockVector<Number>::block (const unsigned int i) const
{
  Assert(i<num_blocks, ExcIndexRange(i,0,num_blocks));

  return components[i];
}



template <typename Number>
inline
const BlockIndices&
BlockVector<Number>::get_block_indices () const
{
  return block_indices;
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
void swap (BlockVector<Number> &u,
	   BlockVector<Number> &v)
{
  u.swap (v);
};




#endif

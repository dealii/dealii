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


#include <base/config.h>
#include <lac/vector.h>
#include <lac/block_indices.h>
#include <base/exceptions.h>
#include <cstdio>
#include <vector>
#include <iterator>


// forward declaration
template <typename Number>
class BlockVector;


/**
 * Namespace in which iterators in block vectors are implemented.
 *
 * @author Wolfgang Bangerth, 2001
 */
namespace BlockVectorIterators
{
				   /**
				    * Declaration of the general
				    * template of a structure which is
				    * used to determine some types
				    * based on the template arguments
				    * of other classes.
				    */
  template <typename number, bool constness>
  struct Types
  {
  };


  
				   /**
				    * Declaration of a specialized
				    * template of a structure which is
				    * used to determine some types
				    * based on the template arguments
				    * of other classes.
				    *
				    * This is for the use of non-const
				    * iterators.
				    */
  template <typename number>
  struct Types<number,false>
  {
				       /**
					* Type of the vector
					* underlying the block vector
					* used in non-const
					* iterators. There, the
					* vector must not be constant.
					*/
      typedef Vector<number>      Vector;

				       /**
					* Type of the block vector
					* used in non-const
					* iterators. There, the block
					* vector must not be constant.
					*/
      typedef BlockVector<number> BlockVector;

				       /**
					* Type of the numbers we point
					* to. Here, they are not
					* constant.
					*/
      typedef number              NumberType;
  };


  
				   /**
				    * Declaration of a specialized
				    * template of a structure which is
				    * used to determine some types
				    * based on the template arguments
				    * of other classes.
				    *
				    * This is for the use of
				    * @p{const_iterator}s.
				    */
  template <typename number>
  struct Types<number,true>
  {
				       /**
					* Type of the vector
					* underlying the block vector
					* used in
					* @p{const_iterators}. There,
					* the vector must be
					* constant.
					*/
      typedef const Vector<number>      Vector;

				       /**
					* Type of the block vector
					* used in
					* @p{const_iterator}s. There,
					* the block vector must be
					* constant.
					*/
      typedef const BlockVector<number> BlockVector;

				       /**
					* Type of the numbers we point
					* to. Here, they are constant
					* since the block vector we
					* use is constant.
					*/
      typedef const number              NumberType;
  };


				   /**
				    * General random-access iterator
				    * class for block vectors. Since
				    * we do not want to have two
				    * classes for non-const
				    * @p{iterator}s and
				    * @p{const_iterator}s, we take a
				    * second template argument which
				    * denotes whether the vector we
				    * point into is a constant object
				    * or not. The first template
				    * argument is always the number
				    * type of the block vector in use.
				    *
				    * This class satisfies all
				    * requirements of random access
				    * iterators defined in the C++
				    * standard. Operations on these
				    * iterators are constant in the
				    * number of elements in the block
				    * vector. However, they are
				    * sometimes linear in the number
				    * of blocks in the vector, but
				    * since that does rarely change
				    * dynamically within an
				    * application, this is a constant
				    * and we again have that the
				    * iterator satisfies the
				    * requirements of a random access
				    * iterator.
				    *
				    * The implementation of this class
				    * has to work around some problems
				    * in compilers and standard
				    * libraries. One of these requires
				    * us to write all comparison
				    * operators twice, once comparison
				    * with iterators of the same type
				    * and once with iterators pointing
				    * to numbers of opposite constness
				    * specification. The reason is
				    * that if we would have written
				    * the comparison operators as a
				    * template on the constness of the
				    * right hand side, then gcc2.95
				    * signals an error that these
				    * operators ambiguate operators
				    * declared somewhere within the
				    * standard library. Likewise, we
				    * have to work around some
				    * problems with granting other
				    * iterators friendship. This makes
				    * the implementation somewhat
				    * non-optimal at places, but at
				    * least everything works.
				    *
				    * @author Wolfgang Bangerth, 2001
				    */
  template <typename number, bool constness>
  class Iterator :
#ifdef HAVE_STD_ITERATOR_CLASS  
               public std::iterator<std::random_access_iterator_tag,
				    typename Types<number,constness>::NumberType>
#else
               random_access_iterator<typename Types<number,constness>::NumberType,int>
#endif      
  {
    private:
				       /**
					* Typedef an iterator with
					* opposite constness
					* requirements on the elements
					* it points to.
					*/
      typedef Iterator<number,!constness> InverseConstnessIterator;

    public:
				       /**
					* Declare some typedefs which
					* are standard for iterators
					* and are used by algorithms
					* to enquire about the
					* specifics of the iterators
					* they work on.
					*/
      typedef std::random_access_iterator_tag               iterator_type;
      typedef typename Types<number,constness>::NumberType  value_type;
      typedef ptrdiff_t                                     difference_type;
      typedef value_type                                   &reference;
      typedef value_type                                   *pointer;
      
				       /**
					* Typedef the type of the
					* block vector (which differs
					* in constness, depending on
					* the second template
					* parameter).
					*/
      typedef typename Types<number,constness>::BlockVector BlockVectorType;

				       /**
					* Type of the number this
					* iterator points
					* to. Depending on the value
					* of the second template
					* parameter, this is either a
					* constant or non-const
					* number.
					*/
      typedef typename Types<number,constness>::NumberType NumberType;
	
				       /**
					* Construct an iterator from
					* a vector to which we point
					* and the global index of
					* the element pointed to.
					*
					* Depending on the value of
					* the @p{constness} template
					* argument of this class,
					* the first argument of this
					* constructor is either is a
					* @p{const} or non-@p{const}
					* reference.
					*/
      Iterator (BlockVectorType    &parent,
		const unsigned int  global_index);
	
				       /**
					* Copy constructor.
					*/
      Iterator (const Iterator<number,constness> &c);

				       /**
					* Copy constructor for
					* conversion between iterators
					* with different constness
					* requirements. This
					* constructor throws an error
					* if an attempt is made at
					* converting a constant to a
					* non-constant iterator.
					*/
      Iterator (const InverseConstnessIterator &c);

    private:
				       /**
					* Constructor used internally
					* in this class. The arguments
					* match exactly the values of
					* the respective member
					* variables.
					*/
      Iterator (BlockVectorType    &parent,
		const unsigned int  global_index,
		const unsigned int  current_block,
		const unsigned int  index_within_block,
		const unsigned int  next_break_forward,
		const unsigned int  next_break_backward);
      
    public:
      
				       /**
					* Copy operator.
					*/
      Iterator & operator = (const Iterator &c);

				       /**
					* Dereferencing operator. If
					* the template argument
					* @p{constness} is @p{true},
					* then no writing to the
					* result is possible, making
					* this a @p{const_iterator}.
					*/
      reference operator * () const;

				       /**
					* Dereferencing operator. If
					* the template argument
					* @p{constness} is @p{true},
					* then no writing to the
					* result is possible, making
					* this a @p{const_iterator}.
					*/ 
      pointer operator -> () const;

				       /**
					* Random access operator,
					* grant access to arbitrary
					* elements relative to the one
					* presently pointed to.
					*/
      reference operator [] (const difference_type d) const;
      
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
					* Same, but compare with an
					* iterator of different
					* constness.
					*/
      bool operator == (const InverseConstnessIterator &i) const;

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
					* Same, but compare with an
					* iterator of different
					* constness.
					*/
      bool operator != (const InverseConstnessIterator &i) const;      

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
					* Same, but compare with an
					* iterator of different
					* constness.
					*/
      bool operator < (const InverseConstnessIterator &i) const;      

				       /**
					* Comparison operator alike
					* to the one above.
					*/
      bool operator <= (const Iterator &i) const;

				       /**
					* Same, but compare with an
					* iterator of different
					* constness.
					*/
      bool operator <= (const InverseConstnessIterator &i) const;      

				       /**
					* Comparison operator alike
					* to the one above.
					*/
      bool operator > (const Iterator &i) const;

				       /**
					* Same, but compare with an
					* iterator of different
					* constness.
					*/
      bool operator > (const InverseConstnessIterator &i) const;      

				       /**
					* Comparison operator alike
					* to the one above.
					*/
      bool operator >= (const Iterator &i) const;

				       /**
					* Same, but compare with an
					* iterator of different
					* constness.
					*/
      bool operator >= (const InverseConstnessIterator &i) const;

				       /**
					* Return the distance between
					* the two iterators, in
					* elements.
					*/
      difference_type operator - (const Iterator &i) const;

				       /**
					* Same, but for iterators of
					* opposite constness.
					*/
      difference_type operator - (const InverseConstnessIterator &i) const;

				       /**
					* Return an iterator which is
					* the given number of elements
					* in front of the present one.
					*/
      Iterator operator + (const difference_type &d) const;
      
				       /**
					* Return an iterator which is
					* the given number of elements
					* behind the present one.
					*/
      Iterator operator - (const difference_type &d) const;

				       /**
					* Move the iterator @p{d}
					* elements forward at once,
					* and return the result.
					*/
      Iterator & operator += (const difference_type &d);

				       /**
					* Move the iterator @p{d}
					* elements backward at once,
					* and return the result.
					*/
      Iterator & operator -= (const difference_type &d);
      
				       /**
					* Exception.
					*/
      DeclException0 (ExcPointerToDifferentVectors);
				       /**
					* Exception.
					*/
      DeclException0 (ExcCastingAwayConstness);
	
    private:
				       /**
					* Pointer to the block
					* vector object to which
					* this iterator
					* points. Depending on the
					* value of the @p{constness}
					* template argument of this
					* class, this is a @p{const}
					* or non-@p{const} pointer.
					*/
      BlockVectorType *parent;

				       /**
					* Global index of the
					* element to which we
					* presently point.
					*/
      unsigned int     global_index;

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

      friend class InverseConstnessIterator;
  };
};



/**
 * A vector composed of several blocks each representing a vector of
 * its own.
 *
 * The @p{BlockVector} is a collection of normal LAC-@ref{Vector}s. Each of
 * the vectors inside can have a different size. The special case of a
 * block vector with constant block size is supported by constructor
 * and @p{reinit} functions.
 *
 * The functionality of @p{BlockVector} includes everything a
 * @p{Vector} can do, plus the access to a single @p{Vector} inside
 * the @p{BlockVector} by @p{block(i)}. It also has a complete random
 * access iterator, just as the LAC-@ref{Vector} class or the standard
 * C++ library template @p{std::vector}. Therefore, all algorithms
 * working on iterators also work with objects of this class.
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
 * @author Wolfgang Bangerth, Guido Kanschat, 1999, 2000, 2001
 */
template <typename Number>
class BlockVector
{
  public:
				     /**
				      * Declare standard types used in
				      * all containers. These types
				      * parallel those in the @p{C++}
				      * standard libraries
				      * @p{vector<...>} class. This
				      * includes iterator types.
				      */
    typedef Number                  value_type;
    typedef value_type             *pointer;
    typedef const value_type       *const_pointer;
    typedef BlockVectorIterators::Iterator<Number,false> iterator;
    typedef BlockVectorIterators::Iterator<Number,true>  const_iterator;
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
    Vector<Number> & block (const unsigned int i);
    
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
    
				     /**
				      * Return an iterator pointing to
				      * the element past the end.
				      */
    iterator end ();

				     /**
				      * Return an iterator pointing to
				      * the element past the end of a
				      * constant block vector.
				      */
    const_iterator end () const;    
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
    friend class iterator;
    friend class const_iterator;
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


namespace BlockVectorIterators
{

  template <typename number, bool constness>
  inline
  Iterator<number,constness>::
  Iterator (const Iterator<number,constness> &c)
		  :
		  parent (c.parent),
		  global_index (c.global_index),
		  current_block (c.current_block),
		  index_within_block (c.index_within_block),
		  next_break_forward (c.next_break_forward),
		  next_break_backward (c.next_break_backward)
  {};



  template <typename number, bool constness>
  inline
  Iterator<number,constness>::
  Iterator (const InverseConstnessIterator &c)
		  :
		  parent (const_cast<BlockVectorType*>(c.parent)),
		  global_index (c.global_index),
		  current_block (c.current_block),
		  index_within_block (c.index_within_block),
		  next_break_forward (c.next_break_forward),
		  next_break_backward (c.next_break_backward)
  {
				     // if constness==false, then the
				     // constness of the iterator we
				     // got is true and we are trying
				     // to cast away the
				     // constness. disallow this
    Assert (constness==true, ExcCastingAwayConstness());
  };
  


  template <typename number, bool constness>
  inline
  Iterator<number,constness>::
  Iterator (BlockVectorType    &parent,
	    const unsigned int  global_index,
	    const unsigned int  current_block,
	    const unsigned int  index_within_block,
	    const unsigned int  next_break_forward,
	    const unsigned int  next_break_backward)
		  :
		  parent (&parent),
		  global_index (global_index),
		  current_block (current_block),
		  index_within_block (index_within_block),
		  next_break_forward (next_break_forward),
		  next_break_backward (next_break_backward)
  {
  };
  


  template <typename number, bool constness>
  inline
  Iterator<number,constness> &
  Iterator<number,constness>::
  operator = (const Iterator &c)
  {
    parent              = c.parent;
    global_index        = c.global_index;
    index_within_block  = c.index_within_block;
    current_block       = c.current_block;
    next_break_forward  = c.next_break_forward;
    next_break_backward = c.next_break_backward;

    return *this;
  };



  template <typename number, bool constness>
  inline
  typename Iterator<number,constness>::reference
  Iterator<number,constness>::operator * () const
  {
				     // we might want to return the
				     // following directly, but if we do
				     // that we run into a gcc bug where
				     // it complains that we return a
				     // reference to a temporary
    const reference p = parent->block(current_block)(index_within_block);
    return p;
  };



  template <typename number, bool constness>
  inline
  typename Iterator<number,constness>::pointer
  Iterator<number,constness>::operator -> () const
  {
				     // we might want to return the
				     // following directly, but if we do
				     // that we run into a gcc bug where
				     // it complains that we return a
				     // reference to a temporary
    const reference p = parent->block(current_block)(index_within_block);
    return &p;
  };
  


  template <typename number, bool constness>
  inline
  typename Iterator<number,constness>::reference
  Iterator<number,constness>::operator [] (const difference_type d) const
  {
				     // if the index pointed to is
				     // still within the block we
				     // currently point into, then we
				     // can save the computation of
				     // the block
    if ((global_index+d >= next_break_backward) &&
	(global_index+d <= next_break_forward))
      {
	const reference p = parent->block(current_block)(index_within_block + d);
	return p;
      };
    
    
				     // if the index is not within the
				     // block of the block vector into
				     // which we presently point, then
				     // there is no way: we have to
				     // search for the block. this can
				     // be done through the parent
				     // class as well.
    const reference p = (*parent)(global_index+d);
    return p;
  };



  template <typename number, bool constness>
  inline
  Iterator<number,constness> &
  Iterator<number,constness>::operator ++ ()
  {
    move_forward ();
    return *this;
  };



  template <typename number, bool constness>
  inline
  Iterator<number,constness>
  Iterator<number,constness>::operator ++ (int)
  {
    const Iterator old_value = *this;
    move_forward ();
    return old_value;
  };



  template <typename number, bool constness>
  inline
  Iterator<number,constness> &
  Iterator<number,constness>::operator -- ()
  {
    move_backward ();
    return *this;
  };



  template <typename number, bool constness>
  inline
  Iterator<number,constness>
  Iterator<number,constness>::operator -- (int)
  {
    const Iterator old_value = *this;
    move_backward ();
    return old_value;
  };



  template <typename number, bool constness>
  inline
  bool
  Iterator<number,constness>::operator == (const Iterator &i) const
  {
    Assert (parent == i.parent, ExcPointerToDifferentVectors());

    return (global_index == i.global_index);
  };



  template <typename number, bool constness>
  inline
  bool
  Iterator<number,constness>::operator == (const InverseConstnessIterator &i) const
  {
    Assert (parent == i.parent, ExcPointerToDifferentVectors());

    return (global_index == i.global_index);
  };
  


  template <typename number, bool constness>
  inline
  bool
  Iterator<number,constness>::operator != (const Iterator &i) const
  {
    Assert (parent == i.parent, ExcPointerToDifferentVectors());

    return (global_index != i.global_index);
  };



  template <typename number, bool constness>
  inline
  bool
  Iterator<number,constness>::operator != (const InverseConstnessIterator &i) const
  {
    Assert (parent == i.parent, ExcPointerToDifferentVectors());

    return (global_index != i.global_index);
  };
  


  template <typename number, bool constness>
  inline
  bool
  Iterator<number,constness>::operator < (const Iterator &i) const
  {
    Assert (parent == i.parent, ExcPointerToDifferentVectors());

    return (global_index < i.global_index);
  };



  template <typename number, bool constness>
  inline
  bool
  Iterator<number,constness>::operator < (const InverseConstnessIterator &i) const
  {
    Assert (parent == i.parent, ExcPointerToDifferentVectors());

    return (global_index < i.global_index);
  };



  template <typename number, bool constness>
  inline
  bool
  Iterator<number,constness>::operator <= (const Iterator &i) const
  {
    Assert (parent == i.parent, ExcPointerToDifferentVectors());

    return (global_index <= i.global_index);
  };



  template <typename number, bool constness>
  inline
  bool
  Iterator<number,constness>::operator <= (const InverseConstnessIterator &i) const
  {
    Assert (parent == i.parent, ExcPointerToDifferentVectors());

    return (global_index <= i.global_index);
  };



  template <typename number, bool constness>
  inline
  bool
  Iterator<number,constness>::operator > (const Iterator &i) const
  {
    Assert (parent == i.parent, ExcPointerToDifferentVectors());

    return (global_index > i.global_index);
  };



  template <typename number, bool constness>
  inline
  bool
  Iterator<number,constness>::operator > (const InverseConstnessIterator &i) const
  {
    Assert (parent == i.parent, ExcPointerToDifferentVectors());

    return (global_index > i.global_index);
  };



  template <typename number, bool constness>
  inline
  bool
  Iterator<number,constness>::operator >= (const Iterator &i) const
  {
    Assert (parent == i.parent, ExcPointerToDifferentVectors());

    return (global_index >= i.global_index);
  };  



  template <typename number, bool constness>
  inline
  bool
  Iterator<number,constness>::operator >= (const InverseConstnessIterator &i) const
  {
    Assert (parent == i.parent, ExcPointerToDifferentVectors());

    return (global_index >= i.global_index);
  };  



  template <typename number, bool constness>
  inline
  typename Iterator<number,constness>::difference_type
  Iterator<number,constness>::operator - (const Iterator &i) const
  {
    Assert (parent == i.parent, ExcPointerToDifferentVectors());

    return (static_cast<signed int>(global_index) -
	    static_cast<signed int>(i.global_index));
  };  



  template <typename number, bool constness>
  inline
  typename Iterator<number,constness>::difference_type
  Iterator<number,constness>::operator - (const InverseConstnessIterator &i) const
  {
    Assert (parent == i.parent, ExcPointerToDifferentVectors());

    return (static_cast<signed int>(global_index) -
	    static_cast<signed int>(i.global_index));
  };



  template <typename number, bool constness>
  inline
  Iterator<number,constness>
  Iterator<number,constness>::operator + (const difference_type &d) const
  {
				     // if the index pointed to is
				     // still within the block we
				     // currently point into, then we
				     // can save the computation of
				     // the block
    if ((global_index+d >= next_break_backward) &&
	(global_index+d <= next_break_forward))
      return Iterator (*parent, global_index+d, current_block,
		       index_within_block+d,
		       next_break_forward, next_break_backward);
    else
				       // outside present block, so
				       // have to seek new block
				       // anyway
      return Iterator (*parent, global_index+d);
  };



  template <typename number, bool constness>
  inline
  Iterator<number,constness>
  Iterator<number,constness>::operator - (const difference_type &d) const
  {
				     // if the index pointed to is
				     // still within the block we
				     // currently point into, then we
				     // can save the computation of
				     // the block
    if ((global_index-d >= next_break_backward) &&
	(global_index-d <= next_break_forward))
      return Iterator (*parent, global_index-d, current_block,
		       index_within_block-d,
		       next_break_forward, next_break_backward);
    else
				       // outside present block, so
				       // have to seek new block
				       // anyway
      return Iterator (*parent, global_index-d);
  };



  template <typename number, bool constness>
  inline
  Iterator<number,constness> &
  Iterator<number,constness>::operator += (const difference_type &d)
  {
				     // if the index pointed to is
				     // still within the block we
				     // currently point into, then we
				     // can save the computation of
				     // the block
    if ((global_index+d >= next_break_backward) &&
	(global_index+d <= next_break_forward))
      {
	global_index       += d;
	index_within_block += d;
      }
    else
				       // outside present block, so
				       // have to seek new block
				       // anyway
      *this = Iterator (*parent, global_index+d);

    return *this;
  };



  template <typename number, bool constness>
  inline
  Iterator<number,constness> &
  Iterator<number,constness>::operator -= (const difference_type &d)
  {
				     // if the index pointed to is
				     // still within the block we
				     // currently point into, then we
				     // can save the computation of
				     // the block
    if ((global_index-d >= next_break_backward) &&
	(global_index-d <= next_break_forward))
      {
	global_index       -= d;
	index_within_block -= d;
      }
    else
				       // outside present block, so
				       // have to seek new block
				       // anyway
      *this = Iterator (*parent, global_index-d);

    return *this;
  };
  
};


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

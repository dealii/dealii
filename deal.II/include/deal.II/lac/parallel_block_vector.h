//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__parallel_block_vector_h
#define __deal2__parallel_block_vector_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/lac/block_indices.h>
#include <deal.II/lac/block_vector_base.h>

#include <cstdio>
#include <vector>

DEAL_II_NAMESPACE_OPEN


// TODO: global reduction operations (operator *, {l1,l2,lp,linfty}norm, mean
// value) should use one MPI communication with several Number values, not use
// the parallel::distributed::Vector operation directly.

namespace parallel
{
  namespace distributed
  {

    /*! @addtogroup Vectors
     *@{
     */


    /**
     * An implementation of block vectors based on distribued deal.II
     * vectors. While the base class provides for most of the interface, this
     * class handles the actual allocation of vectors and provides functions that
     * are specific to the underlying vector type.
     *
     * @note Instantiations for this template are provided for <tt>@<float@> and
     * @<double@></tt>; others can be generated in application programs (see the
     * section on @ref Instantiations in the manual).
     *
     * @see @ref GlossBlockLA "Block (linear algebra)"
     * @author Katharina Kormann, Martin Kronbichler, 2011
     */
    template <typename Number>
    class BlockVector : public BlockVectorBase<Vector<Number> >
    {
    public:
      /**
       * Declare the type for container size.
       */
      typedef std::size_t size_type;

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
      typedef typename BaseClass::real_type       real_type;
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
      explicit BlockVector (const size_type num_blocks = 0,
                            const size_type block_size = 0);

      /**
       * Copy-Constructor. Dimension set to
       * that of V, all components are copied
       * from V
       */
      BlockVector (const BlockVector<Number> &V);


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
      BlockVector (const std::vector<size_type> &block_sizes);

      /**
       * Destructor. Clears memory
       */
      ~BlockVector ();

      /**
       * Copy operator: fill all components of
       * the vector with the given scalar
       * value.
       */
      BlockVector &operator = (const value_type s);

      /**
       * Copy operator for arguments of the
       * same type. Resize the
       * present vector if necessary.
       */
      BlockVector &
      operator= (const BlockVector &V);

      /**
       * Copy operator for template arguments
       * of different types. Resize the
       * present vector if necessary.
       */
      template <class Number2>
      BlockVector &
      operator= (const BlockVector<Number2> &V);

      /**
       * Copy a regular vector into a
       * block vector.
       */
      BlockVector &
      operator= (const Vector<Number> &V);

      /**
       * Reinitialize the BlockVector to
       * contain <tt>num_blocks</tt> blocks of
       * size <tt>block_size</tt> each.
       *
       * If the second argument is left
       * at its default value, then the
       * block vector allocates the
       * specified number of blocks but
       * leaves them at zero size. You
       * then need to later
       * reinitialize the individual
       * blocks, and call
       * collect_sizes() to update the
       * block system's knowledge of
       * its individual block's sizes.
       *
       * If <tt>fast==false</tt>, the vector
       * is filled with zeros.
       */
      void reinit (const size_type num_blocks,
                   const size_type block_size = 0,
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
      void reinit (const std::vector<size_type> &N,
                   const bool                    fast=false);

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
      void scale (const value_type factor) DEAL_II_DEPRECATED;

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

      /** @addtogroup Exceptions
       * @{ */

      /**
       * Exception
       */
      DeclException0 (ExcIteratorRangeDoesNotMatchVectorSize);
      //@}
    };

    /*@}*/

#ifndef DOXYGEN
    /*----------------------- Inline functions ----------------------------------*/


    template <typename Number>
    inline
    BlockVector<Number>::BlockVector (const size_type n_blocks,
                                      const size_type block_size)
    {
      reinit (n_blocks, block_size);
    }



    template <typename Number>
    inline
    BlockVector<Number>::BlockVector (const std::vector<size_type> &n)
    {
      reinit (n, false);
    }



    template <typename Number>
    inline
    BlockVector<Number>::BlockVector (const BlockVector<Number> &v)
      :
      BlockVectorBase<Vector<Number> > ()
    {
      this->components.resize (v.n_blocks());
      this->block_indices = v.block_indices;

      for (size_type i=0; i<this->n_blocks(); ++i)
        this->components[i] = v.components[i];
    }


#ifndef DEAL_II_EXPLICIT_CONSTRUCTOR_BUG

    template <typename Number>
    template <typename OtherNumber>
    inline
    BlockVector<Number>::BlockVector (const BlockVector<OtherNumber> &v)
    {
      reinit (v, true);
      *this = v;
    }

#endif



    template <typename Number>
    inline
    void BlockVector<Number>::reinit (const size_type n_bl,
                                      const size_type bl_sz,
                                      const bool         fast)
    {
      std::vector<size_type> n(n_bl, bl_sz);
      reinit(n, fast);
    }


    template <typename Number>
    inline
    void BlockVector<Number>::reinit (const std::vector<size_type> &n,
                                      const bool                    fast)
    {
      this->block_indices.reinit (n);
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());

      for (size_type i=0; i<this->n_blocks(); ++i)
        this->components[i].reinit(n[i], fast);
    }



    template <typename Number>
    template <typename Number2>
    inline
    void BlockVector<Number>::reinit (const BlockVector<Number2> &v,
                                      const bool fast)
    {
      this->block_indices = v.get_block_indices();
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());

      for (size_type i=0; i<this->n_blocks(); ++i)
        this->block(i).reinit(v.block(i), fast);
    }



    template <typename Number>
    inline
    BlockVector<Number>::~BlockVector ()
    {}



    template <typename Number>
    inline
    BlockVector<Number> &
    BlockVector<Number>::operator = (const value_type s)
    {

      Assert (numbers::is_finite(s), ExcNumberNotFinite());

      BaseClass::operator = (s);
      return *this;
    }



    template <typename Number>
    inline
    BlockVector<Number> &
    BlockVector<Number>::operator = (const BlockVector &v)
    {
      reinit (v, true);
      BaseClass::operator = (v);
      return *this;
    }



    template <typename Number>
    inline
    BlockVector<Number> &
    BlockVector<Number>::operator = (const Vector<Number> &v)
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
      reinit (v, true);
      BaseClass::operator = (v);
      return *this;
    }



    template <typename Number>
    inline
    void BlockVector<Number>::swap (BlockVector<Number> &v)
    {
      Assert (this->n_blocks() == v.n_blocks(),
              ExcDimensionMismatch(this->n_blocks(), v.n_blocks()));

      for (size_type i=0; i<this->n_blocks(); ++i)
        dealii::swap (this->components[i], v.components[i]);
      dealii::swap (this->block_indices, v.block_indices);
    }



    template <typename Number>
    void BlockVector<Number>::scale (const value_type factor)
    {

      Assert (numbers::is_finite(factor), ExcNumberNotFinite());

      for (size_type i=0; i<this->n_blocks(); ++i)
        this->components[i].scale(factor);
    }



    template <typename Number>
    template <class BlockVector2>
    void BlockVector<Number>::scale (const BlockVector2 &v)
    {
      BaseClass::scale (v);
    }

#endif // DOXYGEN

  } // end of namespace distributed

} // end of namespace parallel

/**
 * Global function which overloads the default implementation
 * of the C++ standard library which uses a temporary object. The
 * function simply exchanges the data of the two vectors.
 *
 * @relates BlockVector
 * @author Katharina Kormann, Martin Kronbichler, 2011
 */
template <typename Number>
inline
void swap (parallel::distributed::BlockVector<Number> &u,
           parallel::distributed::BlockVector<Number> &v)
{
  u.swap (v);
}

DEAL_II_NAMESPACE_CLOSE

#endif

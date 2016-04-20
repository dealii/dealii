// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii__petsc_block_vector_h
#define dealii__petsc_block_vector_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/petsc_vector.h>
#  include <deal.II/lac/petsc_parallel_block_vector.h>
#  include <deal.II/lac/block_indices.h>
#  include <deal.II/lac/block_vector_base.h>
#  include <deal.II/lac/exceptions.h>

DEAL_II_NAMESPACE_OPEN



namespace PETScWrappers
{
  /*! @addtogroup PETScWrappers
   *@{
   */

  /**
   * An implementation of block vectors based on the vector class implemented
   * in PETScWrappers. While the base class provides for most of the
   * interface, this class handles the actual allocation of vectors and
   * provides functions that are specific to the underlying vector type.
   *
   * This class is deprecated, use PETScWrappers::MPI::BlockVector.
   *
   * @ingroup Vectors
   *
   * @see
   * @ref GlossBlockLA "Block (linear algebra)"
   * @author Wolfgang Bangerth, 2004
   */
  class BlockVector : public BlockVectorBase<Vector>
  {
  public:
    /**
     * Typedef the base class for simpler access to its own typedefs.
     */
    typedef BlockVectorBase<Vector> BaseClass;

    /**
     * Typedef the type of the underlying vector.
     */
    typedef BaseClass::BlockType  BlockType;

    /**
     * Import the typedefs from the base class.
     */
    typedef BaseClass::value_type      value_type;
    typedef BaseClass::pointer         pointer;
    typedef BaseClass::const_pointer   const_pointer;
    typedef BaseClass::reference       reference;
    typedef BaseClass::const_reference const_reference;
    typedef BaseClass::size_type       size_type;
    typedef BaseClass::iterator        iterator;
    typedef BaseClass::const_iterator  const_iterator;

    /**
     * Constructor. There are three ways to use this constructor. First,
     * without any arguments, it generates an object with no blocks. Given one
     * argument, it initializes <tt>num_blocks</tt> blocks, but these blocks
     * have size zero. The third variant finally initializes all blocks to the
     * same size <tt>block_size</tt>.
     *
     * Confer the other constructor further down if you intend to use blocks
     * of different sizes.
     */
    explicit BlockVector (const unsigned int num_blocks = 0,
                          const size_type    block_size = 0);

    /**
     * Copy-Constructor. Dimension set to that of V, all components are copied
     * from V
     */
    BlockVector (const BlockVector  &V);

    /**
     * Copy-constructor: copy the values from a PETSc wrapper parallel block
     * vector class.
     *
     *
     * Note that due to the communication model of MPI, @em all processes have
     * to actually perform this operation, even if they do not use the result.
     * It is not sufficient if only one processor tries to copy the elements
     * from the other processors over to its own process space.
     */
    explicit BlockVector (const MPI::BlockVector &v);

    /**
     * Constructor. Set the number of blocks to <tt>n.size()</tt> and
     * initialize each block with <tt>n[i]</tt> zero elements.
     */
    BlockVector (const std::vector<size_type> &n);

    /**
     * Constructor. Set the number of blocks to <tt>n.size()</tt>. Initialize
     * the vector with the elements pointed to by the range of iterators given
     * as second and third argument. Apart from the first argument, this
     * constructor is in complete analogy to the respective constructor of the
     * <tt>std::vector</tt> class, but the first argument is needed in order
     * to know how to subdivide the block vector into different blocks.
     */
    template <typename InputIterator>
    BlockVector (const std::vector<size_type> &n,
                 const InputIterator           first,
                 const InputIterator           end);

    /**
     * Destructor. Clears memory
     */
    ~BlockVector ();

    /**
     * Copy operator: fill all components of the vector with the given scalar
     * value.
     */
    BlockVector &operator = (const value_type s);

    /**
     * Copy operator for arguments of the same type.
     */
    BlockVector &
    operator= (const BlockVector &V);

    /**
     * Copy all the elements of the parallel block vector @p v into this local
     * vector. Note that due to the communication model of MPI, @em all
     * processes have to actually perform this operation, even if they do not
     * use the result. It is not sufficient if only one processor tries to
     * copy the elements from the other processors over to its own process
     * space.
     */
    BlockVector &
    operator = (const MPI::BlockVector &v);

    /**
     * Reinitialize the BlockVector to contain <tt>num_blocks</tt> blocks of
     * size <tt>block_size</tt> each.
     *
     * If <tt>omit_zeroing_entries==false</tt>, the vector is filled with
     * zeros.
     */
    void reinit (const unsigned int num_blocks,
                 const size_type    block_size,
                 const bool omit_zeroing_entries = false);

    /**
     * Reinitialize the BlockVector such that it contains
     * <tt>block_sizes.size()</tt> blocks. Each block is reinitialized to
     * dimension <tt>block_sizes[i]</tt>.
     *
     * If the number of blocks is the same as before this function was called,
     * all vectors remain the same and reinit() is called for each vector.
     *
     * If <tt>omit_zeroing_entries==false</tt>, the vector is filled with
     * zeros.
     *
     * Note that you must call this (or the other reinit() functions)
     * function, rather than calling the reinit() functions of an individual
     * block, to allow the block vector to update its caches of vector sizes.
     * If you call reinit() on one of the blocks, then subsequent actions on
     * this object may yield unpredictable results since they may be routed to
     * the wrong block.
     */
    void reinit (const std::vector<size_type> &N,
                 const bool                   omit_zeroing_entries=false);

    /**
     * Change the dimension to that of the vector <tt>V</tt>. The same applies
     * as for the other reinit() function.
     *
     * The elements of <tt>V</tt> are not copied, i.e.  this function is the
     * same as calling <tt>reinit (V.size(), omit_zeroing_entries)</tt>.
     *
     * Note that you must call this (or the other reinit() functions)
     * function, rather than calling the reinit() functions of an individual
     * block, to allow the block vector to update its caches of vector sizes.
     * If you call reinit() of one of the blocks, then subsequent actions of
     * this object may yield unpredictable results since they may be routed to
     * the wrong block.
     */
    void reinit (const BlockVector &V,
                 const bool         omit_zeroing_entries=false);

    /**
     * Change the number of blocks to <tt>num_blocks</tt>. The individual
     * blocks will get initialized with zero size, so it is assumed that the
     * user resizes the individual blocks by herself in an appropriate way,
     * and calls <tt>collect_sizes</tt> afterwards.
     */
    void reinit (const unsigned int num_blocks);

    /**
     * Swap the contents of this vector and the other vector <tt>v</tt>. One
     * could do this operation with a temporary variable and copying over the
     * data elements, but this function is significantly more efficient since
     * it only swaps the pointers to the data of the two vectors and therefore
     * does not need to allocate temporary storage and move data around.
     *
     * Limitation: right now this function only works if both vectors have the
     * same number of blocks. If needed, the numbers of blocks should be
     * exchanged, too.
     *
     * This function is analog to the the swap() function of all C++ standard
     * containers. Also, there is a global function swap(u,v) that simply
     * calls <tt>u.swap(v)</tt>, again in analogy to standard functions.
     */
    void swap (BlockVector &v);

    /**
     * Print to a stream.
     */
    void print (std::ostream       &out,
                const unsigned int  precision = 3,
                const bool          scientific = true,
                const bool          across = true) const;

    /**
     * @addtogroup Exceptions
     * @{
     */

    /**
     * Exception
     */
    DeclException0 (ExcIteratorRangeDoesNotMatchVectorSize);
    ///@}
  } DEAL_II_DEPRECATED;

  /*@}*/

  /*----------------------- Inline functions ----------------------------------*/



  inline
  BlockVector::BlockVector (const unsigned int n_blocks,
                            const size_type    block_size)
  {
    reinit (n_blocks, block_size);
  }



  inline
  BlockVector::BlockVector (const std::vector<size_type> &n)
  {
    reinit (n, false);
  }


  inline
  BlockVector::BlockVector (const BlockVector &v)
    :
    BlockVectorBase<Vector > ()
  {
    this->components.resize (v.n_blocks());
    block_indices = v.block_indices;

    for (unsigned int i=0; i<this->n_blocks(); ++i)
      this->components[i] = v.components[i];
  }



  inline
  BlockVector::BlockVector (const MPI::BlockVector &v)
    :
    BlockVectorBase<Vector > ()
  {
    this->components.resize (v.get_block_indices().size());
    block_indices = v.get_block_indices();

    for (unsigned int i=0; i<this->n_blocks(); ++i)
      this->components[i] = v.block(i);
  }



  template <typename InputIterator>
  BlockVector::BlockVector (const std::vector<size_type> &n,
                            const InputIterator           first,
                            const InputIterator           end)
  {
    // first set sizes of blocks, but
    // don't initialize them as we will
    // copy elements soon
    (void)end;
    reinit (n, true);
    InputIterator start = first;
    for (unsigned int b=0; b<n.size(); ++b)
      {
        InputIterator end = start;
        std::advance (end, static_cast<signed int>(n[b]));

        for (size_type i=0; i<n[b]; ++i, ++start)
          this->block(b)(i) = *start;
      }
    Assert (start == end, ExcIteratorRangeDoesNotMatchVectorSize());
  }



  inline
  BlockVector &
  BlockVector::operator = (const value_type s)
  {
    BaseClass::operator = (s);
    return *this;
  }



  inline
  BlockVector &
  BlockVector::operator = (const BlockVector &v)
  {
    BaseClass::operator = (v);
    return *this;
  }



  inline
  BlockVector &
  BlockVector::operator = (const MPI::BlockVector &v)
  {
    BaseClass::operator = (v);
    return *this;
  }



  inline
  BlockVector::~BlockVector ()
  {}


  inline
  void
  BlockVector::reinit (const unsigned int n_bl,
                       const size_type    bl_sz,
                       const bool         omit_zeroing_entries)
  {
    std::vector<size_type> n(n_bl, bl_sz);
    reinit(n, omit_zeroing_entries);
  }



  inline
  void
  BlockVector::reinit (const std::vector<size_type> &n,
                       const bool                    omit_zeroing_entries)
  {
    block_indices.reinit (n);
    if (this->components.size() != this->n_blocks())
      this->components.resize(this->n_blocks());

    for (unsigned int i=0; i<this->n_blocks(); ++i)
      this->components[i].reinit(n[i], omit_zeroing_entries);
  }


  inline
  void
  BlockVector::reinit (const BlockVector &v,
                       const bool omit_zeroing_entries)
  {
    block_indices = v.get_block_indices();
    if (this->components.size() != this->n_blocks())
      this->components.resize(this->n_blocks());

    for (unsigned int i=0; i<this->n_blocks(); ++i)
      block(i).reinit(v.block(i), omit_zeroing_entries);
  }



  inline
  void
  BlockVector::reinit (const unsigned int num_blocks)
  {
    reinit (num_blocks, 0, true);
  }



  inline
  void
  BlockVector::swap (BlockVector &v)
  {
    Assert (this->n_blocks() == v.n_blocks(),
            ExcDimensionMismatch(this->n_blocks(), v.n_blocks()));

    for (unsigned int i=0; i<this->n_blocks(); ++i)
      this->components[i].swap (v.components[i]);
    ::dealii::swap (this->block_indices, v.block_indices);
  }



  inline
  void
  BlockVector::print (std::ostream       &out,
                      const unsigned int  precision,
                      const bool          scientific,
                      const bool          across) const
  {
    for (unsigned int i=0; i<this->n_blocks(); ++i)
      {
        if (across)
          out << 'C' << i << ':';
        else
          out << "Component " << i << std::endl;
        this->components[i].print(out, precision, scientific, across);
      }
  }




  /**
   * Global function which overloads the default implementation of the C++
   * standard library which uses a temporary object. The function simply
   * exchanges the data of the two vectors.
   *
   * @relates PETScWrappers::BlockVector
   * @author Wolfgang Bangerth, 2000
   */
  inline
  void swap (BlockVector &u,
             BlockVector &v)
  {
    u.swap (v);
  }

}


namespace internal
{
  namespace LinearOperator
  {
    template <typename> class ReinitHelper;

    /**
     * A helper class used internally in linear_operator.h. Specialization for
     * PETScWrappers::BlockVector.
     */
    template<>
    class ReinitHelper<PETScWrappers::BlockVector>
    {
    public:
      template <typename Matrix>
      static
      void reinit_range_vector (const Matrix &matrix,
                                PETScWrappers::BlockVector &v,
                                bool omit_zeroing_entries)
      {
        v.reinit(matrix.locally_owned_range_sizes(), omit_zeroing_entries);
      }

      template <typename Matrix>
      static
      void reinit_domain_vector(const Matrix &matrix,
                                PETScWrappers::BlockVector &v,
                                bool omit_zeroing_entries)
      {
        v.reinit(matrix.locally_owned_domain_sizes(), omit_zeroing_entries);
      }
    };

  } /* namespace LinearOperator */
} /* namespace internal */

DEAL_II_NAMESPACE_CLOSE

#endif  // DEAL_II_WITH_PETSC

#endif

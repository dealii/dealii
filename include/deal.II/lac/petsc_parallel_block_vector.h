// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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

#ifndef __deal2__petsc_parallel_block_vector_h
#define __deal2__petsc_parallel_block_vector_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/petsc_parallel_vector.h>
#  include <deal.II/lac/block_indices.h>
#  include <deal.II/lac/block_vector_base.h>
#  include <deal.II/lac/exceptions.h>

DEAL_II_NAMESPACE_OPEN


namespace PETScWrappers
{
  // forward declaration
  class BlockVector;

  namespace MPI
  {

    /*! @addtogroup PETScWrappers
     *@{
     */

    /**
     * An implementation of block vectors based on the parallel vector class
     * implemented in PETScWrappers. While the base class provides for most of the
     * interface, this class handles the actual allocation of vectors and provides
     * functions that are specific to the underlying vector type.
     *
     * The model of distribution of data is such that each of the blocks is
     * distributed across all MPI processes named in the MPI communicator. I.e. we
     * don't just distribute the whole vector, but each component. In the
     * constructors and reinit() functions, one therefore not only has to specify
     * the sizes of the individual blocks, but also the number of elements of each
     * of these blocks to be stored on the local process.
     *
     * @ingroup Vectors
     * @see @ref GlossBlockLA "Block (linear algebra)"
     * @author Wolfgang Bangerth, 2004
     */
    class BlockVector : public BlockVectorBase<Vector>
    {
    public:
      /**
       * Typedef the base class for simpler
       * access to its own typedefs.
       */
      typedef BlockVectorBase<Vector> BaseClass;

      /**
       * Typedef the type of the underlying
       * vector.
       */
      typedef BaseClass::BlockType  BlockType;

      /**
       * Import the typedefs from the base
       * class.
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
       * Default constructor. Generate an
       * empty vector without any blocks.
       */
      BlockVector ();

      /**
       *  Constructor. Generate a block
       *  vector with @p n_blocks blocks,
       *  each of which is a parallel
       *  vector across @p communicator
       *  with @p block_size elements of
       *  which @p local_size elements are
       *  stored on the present process.
       */
      explicit BlockVector (const unsigned int  n_blocks,
                            const MPI_Comm     &communicator,
                            const size_type     block_size,
                            const size_type     local_size);

      /**
       * Copy-Constructor. Set all the
       * properties of the parallel vector
       * to those of the given argument and
       * copy the elements.
       */
      BlockVector (const BlockVector  &V);

      /**
       * Constructor. Set the number of
       * blocks to
       * <tt>block_sizes.size()</tt> and
       * initialize each block with
       * <tt>block_sizes[i]</tt> zero
       * elements. The individual blocks
       * are distributed across the given
       * communicator, and each store
       * <tt>local_elements[i]</tt>
       * elements on the present process.
       */
      BlockVector (const std::vector<size_type> &block_sizes,
                   const MPI_Comm               &communicator,
                   const std::vector<size_type> &local_elements);

      /**
       * Create a BlockVector with parallel_partitioning.size() blocks,
       * each initialized with the given IndexSet.
       */
      explicit BlockVector (const std::vector<IndexSet> &parallel_partitioning,
                            const MPI_Comm &communicator = MPI_COMM_WORLD);

      /**
       * Same as above, but include ghost elements
       */
      BlockVector (const std::vector<IndexSet> &parallel_partitioning,
                   const std::vector<IndexSet> &ghost_indices,
                   const MPI_Comm &communicator);



      /**
       * Destructor. Clears memory
       */
      ~BlockVector ();

      /**
       * Copy operator: fill all components
       * of the vector that are locally
       * stored with the given scalar value.
       */
      BlockVector &operator = (const value_type s);

      /**
       * Copy operator for arguments of the
       * same type.
       */
      BlockVector &
      operator= (const BlockVector &V);

      /**
       * Copy the given sequential
       * (non-distributed) block vector
       * into the present parallel block
       * vector. It is assumed that they
       * have the same size, and this
       * operation does not change the
       * partitioning of the parallel
       * vectors by which its elements are
       * distributed across several MPI
       * processes. What this operation
       * therefore does is to copy that
       * chunk of the given vector @p v
       * that corresponds to elements of
       * the target vector that are stored
       * locally, and copies them, for each
       * of the individual blocks of this
       * object. Elements that are not
       * stored locally are not touched.
       *
       * This being a parallel vector, you
       * must make sure that @em all
       * processes call this function at
       * the same time. It is not possible
       * to change the local part of a
       * parallel vector on only one
       * process, independent of what other
       * processes do, with this function.
       */
      BlockVector &
      operator = (const PETScWrappers::BlockVector &v);

      /**
       * Reinitialize the BlockVector to
       * contain @p n_blocks of size @p
       * block_size, each of which stores
       * @p local_size elements
       * locally. The @p communicator
       * argument denotes which MPI channel
       * each of these blocks shall
       * communicate.
       *
       * If <tt>fast==false</tt>, the vector
       * is filled with zeros.
       */
      void reinit (const unsigned int  n_blocks,
                   const MPI_Comm     &communicator,
                   const size_type     block_size,
                   const size_type     local_size,
                   const bool fast = false);

      /**
       * Reinitialize the BlockVector such
       * that it contains
       * <tt>block_sizes.size()</tt>
       * blocks. Each block is
       * reinitialized to dimension
       * <tt>block_sizes[i]</tt>. Each of
       * them stores
       * <tt>local_sizes[i]</tt> elements
       * on the present process.
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
       * reinit() of one of the
       * blocks, then subsequent
       * actions on this object may
       * yield unpredictable results
       * since they may be routed to
       * the wrong block.
       */
      void reinit (const std::vector<size_type> &block_sizes,
                   const MPI_Comm               &communicator,
                   const std::vector<size_type> &local_sizes,
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
       * reinit() on one of the
       * blocks, then subsequent
       * actions on this object may
       * yield unpredictable results
       * since they may be routed to
       * the wrong block.
       */
      void reinit (const BlockVector &V,
                   const bool         fast=false);

      /**
       * Reinitialize the BlockVector using IndexSets. See the constructor
       * with the same arguments for details.
       */
      void reinit (const std::vector<IndexSet> &parallel_partitioning,
                   const MPI_Comm              &communicator);

      /**
       * Same as above but include ghost entries.
       */
      void reinit (const std::vector<IndexSet> &parallel_partitioning,
                   const std::vector<IndexSet> &ghost_entries,
                   const MPI_Comm              &communicator);

      /**
       * Change the number of blocks to
       * <tt>num_blocks</tt>. The individual
       * blocks will get initialized with
       * zero size, so it is assumed that
       * the user resizes the
       * individual blocks by herself
       * in an appropriate way, and
       * calls <tt>collect_sizes</tt>
       * afterwards.
       */
      void reinit (const unsigned int num_blocks);

      /**
       * Returns if this vector is a ghosted vector (and thus read-only).
       */
      bool has_ghost_elements() const;

      /**
       * Return a reference to the MPI
       * communicator object in use with
       * this vector.
       */
      const MPI_Comm &get_mpi_communicator () const;

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
      void swap (BlockVector &v);

      /**
       * Print to a stream.
       */
      void print (std::ostream       &out,
                  const unsigned int  precision = 3,
                  const bool          scientific = true,
                  const bool          across = true) const;

      /**
       * Exception
       */
      DeclException0 (ExcIteratorRangeDoesNotMatchVectorSize);
      /**
       * Exception
       */
      DeclException0 (ExcNonMatchingBlockVectors);
    };

    /*@}*/

    /*----------------------- Inline functions ----------------------------------*/


    inline
    BlockVector::BlockVector ()
    {}



    inline
    BlockVector::BlockVector (const unsigned int  n_blocks,
                              const MPI_Comm     &communicator,
                              const size_type     block_size,
                              const size_type     local_size)
    {
      reinit (n_blocks, communicator, block_size, local_size);
    }



    inline
    BlockVector::BlockVector (const std::vector<size_type> &block_sizes,
                              const MPI_Comm     &communicator,
                              const std::vector<size_type> &local_elements)
    {
      reinit (block_sizes, communicator, local_elements, false);
    }


    inline
    BlockVector::BlockVector (const BlockVector &v)
      :
      BlockVectorBase<Vector > ()
    {
      this->components.resize (v.n_blocks());
      this->block_indices = v.block_indices;

      for (unsigned int i=0; i<this->n_blocks(); ++i)
        this->components[i] = v.components[i];
    }

    inline
    BlockVector::BlockVector (const std::vector<IndexSet> &parallel_partitioning,
                              const MPI_Comm              &communicator)
    {
      reinit(parallel_partitioning, communicator);
    }

    inline
    BlockVector::BlockVector (const std::vector<IndexSet> &parallel_partitioning,
                              const std::vector<IndexSet> &ghost_indices,
                              const MPI_Comm &communicator)
    {
      reinit(parallel_partitioning, ghost_indices, communicator);
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
      // we only allow assignment to vectors with the same number of blocks
      // or to an empty BlockVector
      Assert (n_blocks() == 0 || n_blocks() == v.n_blocks(),
              ExcDimensionMismatch(n_blocks(), v.n_blocks()));

      if (this->n_blocks() != v.n_blocks())
        reinit(v.n_blocks());

      for (size_type i=0; i<this->n_blocks(); ++i)
        this->components[i] = v.block(i);

      collect_sizes();

      return *this;
    }

    inline
    BlockVector::~BlockVector ()
    {}


    inline
    void
    BlockVector::reinit (const unsigned int  n_blocks,
                         const MPI_Comm     &communicator,
                         const size_type     block_size,
                         const size_type     local_size,
                         const bool fast)
    {
      reinit(std::vector<size_type>(n_blocks, block_size),
             communicator,
             std::vector<size_type>(n_blocks, local_size),
             fast);
    }



    inline
    void
    BlockVector::reinit (const std::vector<size_type> &block_sizes,
                         const MPI_Comm               &communicator,
                         const std::vector<size_type> &local_sizes,
                         const bool                    fast)
    {
      this->block_indices.reinit (block_sizes);
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());

      for (unsigned int i=0; i<this->n_blocks(); ++i)
        this->components[i].reinit(communicator, block_sizes[i],
                                   local_sizes[i], fast);
    }


    inline
    void
    BlockVector::reinit (const BlockVector &v,
                         const bool fast)
    {
      this->block_indices = v.get_block_indices();
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());

      for (unsigned int i=0; i<this->n_blocks(); ++i)
        block(i).reinit(v.block(i), fast);
    }

    inline
    void
    BlockVector::reinit (const std::vector<IndexSet> &parallel_partitioning,
                         const MPI_Comm              &communicator)
    {
      std::vector<size_type> sizes(parallel_partitioning.size());
      for (unsigned int i=0; i<parallel_partitioning.size(); ++i)
        sizes[i] = parallel_partitioning[i].size();

      this->block_indices.reinit(sizes);
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());

      for (unsigned int i=0; i<this->n_blocks(); ++i)
        block(i).reinit(parallel_partitioning[i], communicator);
    }

    inline
    void
    BlockVector::reinit (const std::vector<IndexSet> &parallel_partitioning,
                         const std::vector<IndexSet> &ghost_entries,
                         const MPI_Comm              &communicator)
    {
      std::vector<types::global_dof_index> sizes(parallel_partitioning.size());
      for (unsigned int i=0; i<parallel_partitioning.size(); ++i)
        sizes[i] = parallel_partitioning[i].size();

      this->block_indices.reinit(sizes);
      if (this->components.size() != this->n_blocks())
        this->components.resize(this->n_blocks());

      for (unsigned int i=0; i<this->n_blocks(); ++i)
        block(i).reinit(parallel_partitioning[i], ghost_entries[i], communicator);
    }



    inline
    const MPI_Comm &
    BlockVector::get_mpi_communicator () const
    {
      return block(0).get_mpi_communicator();
    }

    inline
    bool
    BlockVector::has_ghost_elements() const
    {
      bool ghosted=block(0).has_ghost_elements();
#ifdef DEBUG
      for (unsigned int i=0; i<this->n_blocks(); ++i)
        Assert(block(i).has_ghost_elements()==ghosted, ExcInternalError());
#endif
      return ghosted;
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
     * Global function which overloads the default implementation
     * of the C++ standard library which uses a temporary object. The
     * function simply exchanges the data of the two vectors.
     *
     * @relates PETScWrappers::MPI::BlockVector
     * @author Wolfgang Bangerth, 2000
     */
    inline
    void swap (BlockVector &u,
               BlockVector &v)
    {
      u.swap (v);
    }

  }

}

DEAL_II_NAMESPACE_CLOSE

#endif  // DEAL_II_WITH_PETSC

#endif

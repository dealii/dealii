// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2016 by the deal.II authors
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

#ifndef dealii__trilinos_parallel_block_vector_h
#define dealii__trilinos_parallel_block_vector_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/lac/trilinos_vector.h>
#  include <deal.II/lac/block_indices.h>
#  include <deal.II/lac/block_vector_base.h>
#  include <deal.II/lac/exceptions.h>

#  include <functional>

DEAL_II_NAMESPACE_OPEN

// forward declaration
template <typename Number> class BlockVector;

/*! @addtogroup TrilinosWrappers
 *@{
 */

namespace TrilinosWrappers
{
  // forward declaration
  namespace MPI
  {
    class BlockVector;
  }
  class BlockVector;
  class BlockSparseMatrix;


  namespace MPI
  {
    /**
     * An implementation of block vectors based on the vector class
     * implemented in TrilinosWrappers. While the base class provides for most
     * of the interface, this class handles the actual allocation of vectors
     * and provides functions that are specific to the underlying vector type.
     *
     * The model of distribution of data is such that each of the blocks is
     * distributed across all MPI processes named in the MPI communicator.
     * I.e. we don't just distribute the whole vector, but each component. In
     * the constructors and reinit() functions, one therefore not only has to
     * specify the sizes of the individual blocks, but also the number of
     * elements of each of these blocks to be stored on the local process.
     *
     * @ingroup Vectors
     * @ingroup TrilinosWrappers @see
     * @ref GlossBlockLA "Block (linear algebra)"
     * @author Martin Kronbichler, Wolfgang Bangerth, 2008, 2009
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
       * Default constructor. Generate an empty vector without any blocks.
       */
      BlockVector ();

      /**
       * Constructor. Generate a block vector with as many blocks as there are
       * entries in @p partitioning. Each Epetra_Map contains the layout of
       * the distribution of data among the MPI processes.
       *
       * This function is deprecated.
       */
      explicit BlockVector (const std::vector<Epetra_Map> &parallel_partitioning) DEAL_II_DEPRECATED;

      /**
       * Constructor. Generate a block vector with as many blocks as there are
       * entries in @p partitioning.  Each IndexSet together with the MPI
       * communicator contains the layout of the distribution of data among
       * the MPI processes.
       */
      explicit BlockVector (const std::vector<IndexSet> &parallel_partitioning,
                            const MPI_Comm              &communicator = MPI_COMM_WORLD);

      /**
       * Creates a BlockVector with ghost elements. See the respective
       * reinit() method for more details. @p ghost_values may contain any
       * elements in @p parallel_partitioning, they will be ignored.
       */
      BlockVector (const std::vector<IndexSet> &parallel_partitioning,
                   const std::vector<IndexSet> &ghost_values,
                   const MPI_Comm              &communicator,
                   const bool                   vector_writable = false);

      /**
       * Copy-Constructor. Set all the properties of the parallel vector to
       * those of the given argument and copy the elements.
       */
      BlockVector (const BlockVector  &v);

#ifdef DEAL_II_WITH_CXX11
      /**
       * Move constructor. Creates a new vector by stealing the internal data
       * of the vector @p v.
       *
       * @note This constructor is only available if deal.II is configured
       * with C++11 support.
       */
      BlockVector (BlockVector &&v);
#endif

      /**
       * Creates a block vector consisting of <tt>num_blocks</tt> components,
       * but there is no content in the individual components and the user has
       * to fill appropriate data using a reinit of the blocks.
       */
      explicit BlockVector (const size_type num_blocks);

      /**
       * Destructor. Clears memory
       */
      ~BlockVector ();

      /**
       * Copy operator: fill all components of the vector that are locally
       * stored with the given scalar value.
       */
      BlockVector &operator= (const value_type s);

      /**
       * Copy operator for arguments of the same type.
       */
      BlockVector &operator= (const BlockVector &v);

#ifdef DEAL_II_WITH_CXX11
      /**
       * Move the given vector. This operator replaces the present vector with
       * @p v by efficiently swapping the internal data structures.
       *
       * @note This operator is only available if deal.II is configured with
       * C++11 support.
       */
      BlockVector &operator= (BlockVector &&v);
#endif

      /**
       * Copy operator for arguments of the localized Trilinos vector type.
       */
      BlockVector &
      operator= (const ::dealii::TrilinosWrappers::BlockVector &v);

      /**
       * Another copy function. This one takes a deal.II block vector and
       * copies it into a TrilinosWrappers block vector. Note that the number
       * of blocks has to be the same in the vector as in the input vector.
       * Use the reinit() command for resizing the BlockVector or for changing
       * the internal structure of the block components.
       *
       * Since Trilinos only works on doubles, this function is limited to
       * accept only one possible number type in the deal.II vector.
       */
      template <typename Number>
      BlockVector &operator= (const ::dealii::BlockVector<Number> &v);

      /**
       * Reinitialize the BlockVector to contain as many blocks as there are
       * Epetra_Maps given in the input argument, according to the parallel
       * distribution of the individual components described in the maps.
       *
       * If <tt>omit_zeroing_entries==false</tt>, the vector is filled with
       * zeros.
       *
       * This function is deprecated.
       */
      void reinit (const std::vector<Epetra_Map> &parallel_partitioning,
                   const bool                     omit_zeroing_entries = false) DEAL_II_DEPRECATED;

      /**
       * Reinitialize the BlockVector to contain as many blocks as there are
       * index sets given in the input argument, according to the parallel
       * distribution of the individual components described in the maps.
       *
       * If <tt>omit_zeroing_entries==false</tt>, the vector is filled with
       * zeros.
       */
      void reinit (const std::vector<IndexSet> &parallel_partitioning,
                   const MPI_Comm              &communicator = MPI_COMM_WORLD,
                   const bool                   omit_zeroing_entries = false);

      /**
       * Reinit functionality. This function destroys the old vector content
       * and generates a new one based on the input partitioning. In addition
       * to just specifying one index set as in all the other methods above,
       * this method allows to supply an additional set of ghost entries.
       * There are two different versions of a vector that can be created. If
       * the flag @p vector_writable is set to @p false, the vector only
       * allows read access to the joint set of @p parallel_partitioning and
       * @p ghost_entries. The effect of the reinit method is then equivalent
       * to calling the other reinit method with an index set containing both
       * the locally owned entries and the ghost entries.
       *
       * If the flag @p vector_writable is set to true, this creates an
       * alternative storage scheme for ghost elements that allows multiple
       * threads to write into the vector (for the other reinit methods, only
       * one thread is allowed to write into the ghost entries at a time).
       */
      void reinit (const std::vector<IndexSet> &partitioning,
                   const std::vector<IndexSet> &ghost_values,
                   const MPI_Comm              &communicator = MPI_COMM_WORLD,
                   const bool      vector_writable = false);


      /**
       * Change the dimension to that of the vector <tt>V</tt>. The same
       * applies as for the other reinit() function.
       *
       * The elements of <tt>V</tt> are not copied, i.e.  this function is the
       * same as calling <tt>reinit (V.size(), omit_zeroing_entries)</tt>.
       *
       * Note that you must call this (or the other reinit() functions)
       * function, rather than calling the reinit() functions of an individual
       * block, to allow the block vector to update its caches of vector
       * sizes. If you call reinit() on one of the blocks, then subsequent
       * actions on this object may yield unpredictable results since they may
       * be routed to the wrong block.
       */
      void reinit (const BlockVector &V,
                   const bool omit_zeroing_entries = false);

      /**
       * Change the number of blocks to <tt>num_blocks</tt>. The individual
       * blocks will get initialized with zero size, so it is assumed that the
       * user resizes the individual blocks by herself in an appropriate way,
       * and calls <tt>collect_sizes</tt> afterwards.
       */
      void reinit (const size_type num_blocks);

      /**
       * This reinit function is meant to be used for parallel calculations
       * where some non-local data has to be used. The typical situation where
       * one needs this function is the call of the
       * FEValues<dim>::get_function_values function (or of some derivatives)
       * in parallel. Since it is usually faster to retrieve the data in
       * advance, this function can be called before the assembly forks out to
       * the different processors. What this function does is the following:
       * It takes the information in the columns of the given matrix and looks
       * which data couples between the different processors. That data is
       * then queried from the input vector. Note that you should not write to
       * the resulting vector any more, since the some data can be stored
       * several times on different processors, leading to unpredictable
       * results. In particular, such a vector cannot be used for matrix-
       * vector products as for example done during the solution of linear
       * systems.
       */
      void import_nonlocal_data_for_fe (const TrilinosWrappers::BlockSparseMatrix &m,
                                        const BlockVector                         &v);

      /**
       * Returns the state of the vector, i.e., whether compress() needs to be
       * called after an operation requiring data exchange. Does only return
       * non-true values when used in <tt>debug</tt> mode, since it is quite
       * expensive to keep track of all operations that lead to the need for
       * compress().
       *
       * This function is deprecated.
       */
      bool is_compressed () const DEAL_II_DEPRECATED;

      /**
       * Returns if this Vector contains ghost elements.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      bool has_ghost_elements() const;

      /**
       * Swap the contents of this vector and the other vector <tt>v</tt>. One
       * could do this operation with a temporary variable and copying over
       * the data elements, but this function is significantly more efficient
       * since it only swaps the pointers to the data of the two vectors and
       * therefore does not need to allocate temporary storage and move data
       * around.
       *
       * Limitation: right now this function only works if both vectors have
       * the same number of blocks. If needed, the numbers of blocks should be
       * exchanged, too.
       *
       * This function is analog to the the swap() function of all C++
       * standard containers. Also, there is a global function swap(u,v) that
       * simply calls <tt>u.swap(v)</tt>, again in analogy to standard
       * functions.
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



    /*----------------------- Inline functions ----------------------------------*/


    inline
    BlockVector::BlockVector ()
    {}



    inline
    BlockVector::BlockVector (const std::vector<IndexSet> &parallel_partitioning,
                              const MPI_Comm              &communicator)
    {
      reinit (parallel_partitioning, communicator, false);
    }



    inline
    BlockVector::BlockVector (const std::vector<IndexSet> &parallel_partitioning,
                              const std::vector<IndexSet> &ghost_values,
                              const MPI_Comm              &communicator,
                              const bool                   vector_writable)
    {
      reinit(parallel_partitioning, ghost_values, communicator,
             vector_writable);
    }



    inline
    BlockVector::BlockVector (const size_type num_blocks)
    {
      reinit (num_blocks);
    }



    inline
    BlockVector::BlockVector (const BlockVector &v)
      :
      BlockVectorBase<Vector > ()
    {
      this->components.resize (v.n_blocks());
      this->block_indices = v.block_indices;

      for (size_type i=0; i<this->n_blocks(); ++i)
        this->components[i] = v.components[i];
    }



#ifdef DEAL_II_WITH_CXX11
    inline
    BlockVector::BlockVector (BlockVector &&v)
    {
      // initialize a minimal, valid object and swap
      reinit (0);
      swap(v);
    }
#endif



    template <typename Number>
    BlockVector &
    BlockVector::operator= (const ::dealii::BlockVector<Number> &v)
    {
      if (n_blocks() != v.n_blocks())
        {
          std::vector<size_type> block_sizes (v.n_blocks(), 0);
          block_indices.reinit (block_sizes);
          if (components.size() != n_blocks())
            components.resize(n_blocks());
        }

      for (size_type i=0; i<this->n_blocks(); ++i)
        this->components[i] = v.block(i);

      collect_sizes();

      return *this;
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
      std::swap(this->components, v.components);

      dealii::swap(this->block_indices, v.block_indices);
    }



    /**
     * Global function which overloads the default implementation of the C++
     * standard library which uses a temporary object. The function simply
     * exchanges the data of the two vectors.
     *
     * @relates TrilinosWrappers::MPI::BlockVector
     * @author Martin Kronbichler, Wolfgang Bangerth, 2008
     */
    inline
    void swap (BlockVector &u,
               BlockVector &v)
    {
      u.swap (v);
    }

  } /* namespace MPI */

} /* namespace TrilinosWrappers */

/*@}*/


namespace internal
{
  namespace LinearOperator
  {
    template <typename> class ReinitHelper;

    /**
     * A helper class used internally in linear_operator.h. Specialization for
     * TrilinosWrappers::MPI::BlockVector.
     */
    template<>
    class ReinitHelper<TrilinosWrappers::MPI::BlockVector>
    {
    public:
      template <typename Matrix>
      static
      void reinit_range_vector (const Matrix &matrix,
                                TrilinosWrappers::MPI::BlockVector &v,
                                bool omit_zeroing_entries)
      {
        v.reinit(matrix.range_partitioner(), omit_zeroing_entries);
      }

      template <typename Matrix>
      static
      void reinit_domain_vector(const Matrix &matrix,
                                TrilinosWrappers::MPI::BlockVector &v,
                                bool omit_zeroing_entries)
      {
        v.reinit(matrix.domain_partitioner(), omit_zeroing_entries);
      }
    };

  } /* namespace LinearOperator */
} /* namespace internal */


DEAL_II_NAMESPACE_CLOSE

#endif  // DEAL_II_WITH_TRILINOS

#endif

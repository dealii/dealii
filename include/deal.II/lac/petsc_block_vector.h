// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_petsc_block_vector_h
#define dealii_petsc_block_vector_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/block_indices.h>
#  include <deal.II/lac/block_vector_base.h>
#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/petsc_vector.h>
#  include <deal.II/lac/vector_type_traits.h>

#  include <cstddef>

DEAL_II_NAMESPACE_OPEN


namespace PETScWrappers
{
  // forward declaration
  class BlockVector;

  namespace MPI
  {
    /**
     * @addtogroup PETScWrappers
     * @{
     */

    /**
     * An implementation of block vectors based on the parallel vector class
     * implemented in PETScWrappers. While the base class provides for most of
     * the interface, this class handles the actual allocation of vectors and
     * provides functions that are specific to the underlying vector type.
     *
     * The model of distribution of data is such that each of the blocks is
     * distributed across all MPI processes named in the MPI communicator.
     * I.e. we don't just distribute the whole vector, but each component. In
     * the constructors and reinit() functions, one therefore not only has to
     * specify the sizes of the individual blocks, but also the number of
     * elements of each of these blocks to be stored on the local process.
     *
     * @ingroup Vectors
     * @see @ref GlossBlockLA "Block (linear algebra)"
     */
    class BlockVector : public BlockVectorBase<Vector>
    {
    public:
      /**
       * Typedef the base class for simpler access to its own alias.
       */
      using BaseClass = BlockVectorBase<Vector>;

      /**
       * Typedef the type of the underlying vector.
       */
      using BlockType = BaseClass::BlockType;

      /**
       * Import the alias from the base class.
       */
      using value_type      = BaseClass::value_type;
      using pointer         = BaseClass::pointer;
      using const_pointer   = BaseClass::const_pointer;
      using reference       = BaseClass::reference;
      using const_reference = BaseClass::const_reference;
      using size_type       = BaseClass::size_type;
      using iterator        = BaseClass::iterator;
      using const_iterator  = BaseClass::const_iterator;

      /**
       * Default constructor. Generate an empty vector without any blocks.
       */
      BlockVector();

      /**
       * Constructor. Generate a block vector with @p n_blocks blocks, each of
       * which is a parallel vector across @p communicator with @p block_size
       * elements of which @p locally_owned_size elements are stored on the
       * present process.
       */
      explicit BlockVector(const unsigned int n_blocks,
                           const MPI_Comm     communicator,
                           const size_type    block_size,
                           const size_type    locally_owned_size);

      /**
       * Copy constructor. Set all the properties of the parallel vector to
       * those of the given argument and copy the elements.
       */
      BlockVector(const BlockVector &V);

      /**
       * Constructor. Set the number of blocks to <tt>block_sizes.size()</tt>
       * and initialize each block with <tt>block_sizes[i]</tt> zero elements.
       * The individual blocks are distributed across the given communicator,
       * and each store <tt>local_elements[i]</tt> elements on the present
       * process.
       */
      BlockVector(const std::vector<size_type> &block_sizes,
                  const MPI_Comm                communicator,
                  const std::vector<size_type> &local_elements);

      /**
       * Create a BlockVector with parallel_partitioning.size() blocks, each
       * initialized with the given IndexSet.
       */
      explicit BlockVector(const std::vector<IndexSet> &parallel_partitioning,
                           const MPI_Comm communicator = MPI_COMM_WORLD);

      /**
       * Same as above, but include ghost elements
       */
      BlockVector(const std::vector<IndexSet> &parallel_partitioning,
                  const std::vector<IndexSet> &ghost_indices,
                  const MPI_Comm               communicator);

      /**
       * Create a BlockVector with a PETSc Vec
       * It infers the number of blocks from the Vec if it is of type VECNEST,
       * otherwise the block vector will only have a single block.
       * Internally, we always store a VECNEST vector.
       */
      explicit BlockVector(Vec v);

      /**
       * Create a BlockVector with an array of PETSc vectors.
       */
      template <std::size_t num_blocks>
      explicit BlockVector(const std::array<Vec, num_blocks> &);

      /**
       * Destructor. Clears memory
       */
      ~BlockVector() override;

      /**
       * Copy operator: fill all components of the vector that are locally
       * stored with the given scalar value.
       */
      BlockVector &
      operator=(const value_type s);

      /**
       * Copy operator for arguments of the same type. Resize the present vector
       * if necessary to the correct number of blocks, then copy the individual
       * blocks from `v` using the copy-assignment operator of the class that
       * represents the individual blocks.
       *
       * Copying the vectors that make up individual blocks can have complex
       * semantics in parallel vector classes. See the information provided
       * by the class used to represent the individual blocks.
       */
      BlockVector &
      operator=(const BlockVector &V);

      /**
       * This method associates the PETSc Vec to the instance of the class.
       * Infers the number of blocks from v if it is of type VECNEST, otherwise
       * the block vector will only have a single block.
       */
      void
      reinit(Vec v);

      /**
       * Reinitialize the BlockVector to contain @p n_blocks of size @p
       * block_size, each of which stores @p locally_owned_size elements
       * locally. The @p communicator argument denotes which MPI channel each
       * of these blocks shall communicate.
       *
       * If <tt>omit_zeroing_entries==false</tt>, the vector is filled with
       * zeros.
       */
      void
      reinit(const unsigned int n_blocks,
             const MPI_Comm     communicator,
             const size_type    block_size,
             const size_type    locally_owned_size,
             const bool         omit_zeroing_entries = false);

      /**
       * Reinitialize the BlockVector such that it contains
       * <tt>block_sizes.size()</tt> blocks. Each block is reinitialized to
       * dimension <tt>block_sizes[i]</tt>. Each of them stores
       * <tt>locally_owned_sizes[i]</tt> elements on the present process.
       *
       * If the number of blocks is the same as before this function was
       * called, all vectors remain the same and reinit() is called for each
       * vector.
       *
       * If <tt>omit_zeroing_entries==false</tt>, the vector is filled with
       * zeros.
       *
       * Note that you must call this (or the other reinit() functions)
       * function, rather than calling the reinit() functions of an individual
       * block, to allow the block vector to update its caches of vector
       * sizes. If you call reinit() of one of the blocks, then subsequent
       * actions on this object may yield unpredictable results since they may
       * be routed to the wrong block.
       */
      void
      reinit(const std::vector<size_type> &block_sizes,
             const MPI_Comm                communicator,
             const std::vector<size_type> &locally_owned_sizes,
             const bool                    omit_zeroing_entries = false);

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
      void
      reinit(const BlockVector &V, const bool omit_zeroing_entries = false);

      /**
       * Reinitialize the BlockVector using IndexSets. See the constructor
       * with the same arguments for details.
       */
      void
      reinit(const std::vector<IndexSet> &parallel_partitioning,
             const MPI_Comm               communicator);

      /**
       * Same as above but include ghost entries.
       */
      void
      reinit(const std::vector<IndexSet> &parallel_partitioning,
             const std::vector<IndexSet> &ghost_entries,
             const MPI_Comm               communicator);

      /**
       * Initialize each block given to each parallel partitioning described in
       * @p partitioners.
       *
       * You can decide whether your vector will contain ghost elements with
       * @p make_ghosted.
       */
      void
      reinit(
        const std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
                  &partitioners,
        const bool make_ghosted = true);

      /**
       * This function collects the sizes of the sub-objects and stores them
       * in internal arrays, in order to be able to relay global indices into
       * the vector to indices into the subobjects. You *must* call this
       * function each time after you have changed the size of the sub-objects.
       */
      void
      collect_sizes();

      /**
       * Call the compress() function on all the subblocks of the vector
       * and update the internal state of the nested PETSc vector.
       *
       * See
       * @ref GlossCompress "Compressing distributed objects"
       * for more information.
       */
      void
      compress(VectorOperation::values operation);

      /**
       * Change the number of blocks to <tt>num_blocks</tt>. The individual
       * blocks will get initialized with zero size, so it is assumed that the
       * user resizes the individual blocks by herself in an appropriate way,
       * and calls <tt>collect_sizes</tt> afterwards.
       */
      void
      reinit(const unsigned int num_blocks);

      /**
       * Return if this vector is a ghosted vector (and thus read-only).
       */
      bool
      has_ghost_elements() const;

      /**
       * Return the underlying MPI communicator.
       */
      MPI_Comm
      get_mpi_communicator() const;

      /**
       * Conversion operator to gain access to the underlying PETSc type. If you
       * do this, you cut this class off some information it may need, so this
       * conversion operator should only be used if you know what you do. In
       * particular, it should only be used for read-only operations into the
       * vector.
       */
      operator const Vec &() const;

      /**
       * Return a reference to the underlying PETSc type. It can be used to
       * modify the underlying data, so use it only when you know what you
       * are doing.
       */
      Vec &
      petsc_vector();

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
       * This function is analogous to the swap() function of all C++
       * standard containers. Also, there is a global function swap(u,v) that
       * simply calls <tt>u.swap(v)</tt>, again in analogy to standard
       * functions.
       */
      void
      swap(BlockVector &v) noexcept;

      /**
       * Print to a stream.
       */
      void
      print(std::ostream      &out,
            const unsigned int precision  = 3,
            const bool         scientific = true,
            const bool         across     = true) const;

      /**
       * Exception
       */
      DeclException0(ExcIteratorRangeDoesNotMatchVectorSize);
      /**
       * Exception
       */
      DeclException0(ExcNonMatchingBlockVectors);

    private:
      /**
       * A PETSc Vec object that describes the entire block vector.
       * Internally, this is done by creating
       * a "nested" vector using PETSc's VECNEST object whose individual
       * blocks are the blocks of this vector.
       */
      Vec petsc_nest_vector;

      /**
       * Utility to set up the VECNEST object
       */
      void
      setup_nest_vec();
    };

    /** @} */

    /*--------------------- Inline functions --------------------------------*/

    inline BlockVector::BlockVector()
      : BlockVectorBase<Vector>()
      , petsc_nest_vector(nullptr)
    {}



    inline BlockVector::BlockVector(const unsigned int n_blocks,
                                    const MPI_Comm     communicator,
                                    const size_type    block_size,
                                    const size_type    locally_owned_size)
      : BlockVector()
    {
      reinit(n_blocks, communicator, block_size, locally_owned_size);
    }



    inline BlockVector::BlockVector(
      const std::vector<size_type> &block_sizes,
      const MPI_Comm                communicator,
      const std::vector<size_type> &local_elements)
      : BlockVector()
    {
      reinit(block_sizes, communicator, local_elements, false);
    }



    inline BlockVector::BlockVector(const BlockVector &v)
      : BlockVector()
    {
      this->block_indices = v.block_indices;

      this->components.resize(this->n_blocks());
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        this->components[i] = v.components[i];

      this->collect_sizes();
    }



    inline BlockVector::BlockVector(
      const std::vector<IndexSet> &parallel_partitioning,
      const MPI_Comm               communicator)
      : BlockVector()
    {
      reinit(parallel_partitioning, communicator);
    }



    inline BlockVector::BlockVector(
      const std::vector<IndexSet> &parallel_partitioning,
      const std::vector<IndexSet> &ghost_indices,
      const MPI_Comm               communicator)
      : BlockVector()
    {
      reinit(parallel_partitioning, ghost_indices, communicator);
    }



    inline BlockVector::BlockVector(Vec v)
      : BlockVector()
    {
      this->reinit(v);
    }



    template <std::size_t num_blocks>
    inline BlockVector::BlockVector(const std::array<Vec, num_blocks> &arrayV)
      : BlockVector()
    {
      this->block_indices.reinit(num_blocks, 0);

      this->components.resize(num_blocks);
      for (auto i = 0; i < num_blocks; ++i)
        this->components[i].reinit(arrayV[i]);
      this->collect_sizes();
    }



    inline BlockVector &
    BlockVector::operator=(const value_type s)
    {
      BaseClass::operator=(s);
      return *this;
    }



    inline BlockVector &
    BlockVector::operator=(const BlockVector &v)
    {
      // we only allow assignment to vectors with the same number of blocks
      // or to an empty BlockVector
      Assert(this->n_blocks() == 0 || this->n_blocks() == v.n_blocks(),
             ExcDimensionMismatch(this->n_blocks(), v.n_blocks()));

      if (this->n_blocks() != v.n_blocks())
        this->block_indices = v.block_indices;

      this->components.resize(this->n_blocks());
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        this->components[i] = v.components[i];

      this->collect_sizes();

      return *this;
    }



    inline void
    BlockVector::reinit(const unsigned int n_blocks,
                        const MPI_Comm     communicator,
                        const size_type    block_size,
                        const size_type    locally_owned_size,
                        const bool         omit_zeroing_entries)
    {
      reinit(std::vector<size_type>(n_blocks, block_size),
             communicator,
             std::vector<size_type>(n_blocks, locally_owned_size),
             omit_zeroing_entries);
    }



    inline void
    BlockVector::reinit(const std::vector<size_type> &block_sizes,
                        const MPI_Comm                communicator,
                        const std::vector<size_type> &locally_owned_sizes,
                        const bool                    omit_zeroing_entries)
    {
      this->block_indices.reinit(block_sizes);

      this->components.resize(this->n_blocks());
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        this->components[i].reinit(communicator,
                                   block_sizes[i],
                                   locally_owned_sizes[i],
                                   omit_zeroing_entries);

      this->collect_sizes();
    }


    inline void
    BlockVector::reinit(const BlockVector &v, const bool omit_zeroing_entries)
    {
      if (this->n_blocks() != v.n_blocks())
        this->block_indices = v.get_block_indices();

      this->components.resize(this->n_blocks());
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        this->components[i].reinit(v.components[i], omit_zeroing_entries);

      this->collect_sizes();
    }



    inline void
    BlockVector::reinit(const std::vector<IndexSet> &parallel_partitioning,
                        const MPI_Comm               communicator)
    {
      // update the number of blocks
      this->block_indices.reinit(parallel_partitioning.size(), 0);

      // initialize each block
      this->components.resize(this->n_blocks());
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        this->components[i].reinit(parallel_partitioning[i], communicator);

      // update block_indices content
      this->collect_sizes();
    }



    inline void
    BlockVector::reinit(const std::vector<IndexSet> &parallel_partitioning,
                        const std::vector<IndexSet> &ghost_entries,
                        const MPI_Comm               communicator)
    {
      AssertDimension(parallel_partitioning.size(), ghost_entries.size());

      // update the number of blocks
      this->block_indices.reinit(parallel_partitioning.size(), 0);

      // initialize each block
      this->components.resize(this->n_blocks());
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        this->components[i].reinit(parallel_partitioning[i],
                                   ghost_entries[i],
                                   communicator);

      // update block_indices content
      this->collect_sizes();
    }



    inline void
    BlockVector::reinit(
      const std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
                &partitioners,
      const bool make_ghosted)
    {
      // update the number of blocks
      this->block_indices.reinit(partitioners.size(), 0);

      // initialize each block
      this->components.resize(this->n_blocks());
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
        this->components[i].reinit(partitioners[i], make_ghosted);

      // update block_indices content
      this->collect_sizes();
    }



    inline MPI_Comm
    BlockVector::get_mpi_communicator() const
    {
      return PetscObjectComm(reinterpret_cast<PetscObject>(petsc_nest_vector));
    }



    inline bool
    BlockVector::has_ghost_elements() const
    {
      bool ghosted = block(0).has_ghost_elements();
      if constexpr (running_in_debug_mode())
        {
          for (unsigned int i = 0; i < this->n_blocks(); ++i)
            Assert(block(i).has_ghost_elements() == ghosted,
                   ExcInternalError());
        }
      return ghosted;
    }



    inline void
    BlockVector::swap(BlockVector &v) noexcept
    {
      std::swap(this->components, v.components);
      std::swap(this->petsc_nest_vector, v.petsc_nest_vector);

      ::dealii::swap(this->block_indices, v.block_indices);
    }



    inline void
    BlockVector::print(std::ostream      &out,
                       const unsigned int precision,
                       const bool         scientific,
                       const bool         across) const
    {
      for (unsigned int i = 0; i < this->n_blocks(); ++i)
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
     * @relatesalso PETScWrappers::MPI::BlockVector
     */
    inline void
    swap(BlockVector &u, BlockVector &v) noexcept
    {
      u.swap(v);
    }
  } // namespace MPI

} // namespace PETScWrappers

namespace internal
{
  namespace LinearOperatorImplementation
  {
    template <typename>
    class ReinitHelper;

    /**
     * A helper class used internally in linear_operator.h. Specialization for
     * PETScWrappers::MPI::BlockVector.
     */
    template <>
    class ReinitHelper<PETScWrappers::MPI::BlockVector>
    {
    public:
      template <typename Matrix>
      static void
      reinit_range_vector(const Matrix                    &matrix,
                          PETScWrappers::MPI::BlockVector &v,
                          bool /*omit_zeroing_entries*/)
      {
        v.reinit(matrix.locally_owned_range_indices(),
                 matrix.get_mpi_communicator());
      }

      template <typename Matrix>
      static void
      reinit_domain_vector(const Matrix                    &matrix,
                           PETScWrappers::MPI::BlockVector &v,
                           bool /*omit_zeroing_entries*/)
      {
        v.reinit(matrix.locally_owned_domain_indices(),
                 matrix.get_mpi_communicator());
      }
    };

  } // namespace LinearOperatorImplementation
} /* namespace internal */


/**
 * Declare dealii::PETScWrappers::MPI::BlockVector as distributed vector.
 */
template <>
struct is_serial_vector<PETScWrappers::MPI::BlockVector> : std::false_type
{};


DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC

#endif

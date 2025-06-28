// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_trilinos_tpetra_block_vector_h
#define dealii_trilinos_tpetra_block_vector_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  include <deal.II/lac/block_indices.h>
#  include <deal.II/lac/block_vector_base.h>
#  include <deal.II/lac/trilinos_tpetra_vector.h>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  /**
   * @addtogroup TpetraWrappers
   * @{
   */
  namespace TpetraWrappers
  {
    /**
     * An implementation of block vectors based on the vector class
     * implemented in LinearAlgebra::TpetraWrappers.
     * While the base class provides for most of the interface, this class
     * handles the actual allocation of vectors and provides functions that
     * are specific to the underlying vector type.
     *
     * The model of distribution of data is such that each of the blocks is
     * distributed across all MPI processes named in the MPI communicator.
     * I.e. we don't just distribute the whole vector, but each component. In
     * the constructors and reinit() functions, one therefore not only has to
     * specify the sizes of the individual blocks, but also the number of
     * elements of each of these blocks to be stored on the local process.
     *
     * @ingroup Vectors
     * @ingroup TrilinosWrappers
     * @see @ref GlossBlockLA "Block (linear algebra)"
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class BlockVector : public dealii::BlockVectorBase<
                          TpetraWrappers::Vector<Number, MemorySpace>>
    {
    public:
      /**
       * Typedef the base class for simpler access to its own alias.
       */
      using BaseClass =
        dealii::BlockVectorBase<TpetraWrappers::Vector<Number, MemorySpace>>;

      /**
       * Typedef the type of the underlying vector.
       */
      using BlockType = typename BaseClass::BlockType;

      /**
       * Import the alias from the base class.
       */
      using value_type      = typename BaseClass::value_type;
      using pointer         = typename BaseClass::pointer;
      using const_pointer   = typename BaseClass::const_pointer;
      using reference       = typename BaseClass::reference;
      using const_reference = typename BaseClass::const_reference;
      using size_type       = typename BaseClass::size_type;
      using iterator        = typename BaseClass::iterator;
      using const_iterator  = typename BaseClass::const_iterator;

      /**
       * Default constructor. Generate an empty vector without any blocks.
       */
      BlockVector() = default;

      /**
       * Constructor. Generate a block vector with as many blocks as there are
       * entries in @p partitioning.  Each IndexSet together with the MPI
       * communicator contains the layout of the distribution of data among
       * the MPI processes.
       */
      BlockVector(const std::vector<IndexSet> &parallel_partitioning,
                  const MPI_Comm               communicator = MPI_COMM_WORLD);

      /**
       * Creates a BlockVector with ghost elements. See the respective
       * reinit() method for more details. @p ghost_values may contain any
       * elements in @p parallel_partitioning, they will be ignored.
       */
      BlockVector(const std::vector<IndexSet> &parallel_partitioning,
                  const std::vector<IndexSet> &ghost_values,
                  const MPI_Comm               communicator,
                  const bool                   vector_writable = false);

      /**
       * Copy-Constructor. Set all the properties of the parallel vector to
       * those of the given argument and copy the elements.
       */
      BlockVector(const BlockVector<Number, MemorySpace> &v);

      /**
       * Creates a block vector consisting of <tt>num_blocks</tt> components,
       * but there is no content in the individual components and the user has
       * to fill appropriate data using a reinit of the blocks.
       */
      explicit BlockVector(const size_type num_blocks);

      /**
       * Destructor. Clears memory
       */
      ~BlockVector() override = default;

      /**
       * Copy operator: fill all components of the vector that are locally
       * stored with the given scalar value.
       */
      BlockVector<Number, MemorySpace> &
      operator=(const Number s);

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
      BlockVector<Number, MemorySpace> &
      operator=(const BlockVector<Number, MemorySpace> &v);

      /**
       * Reinitialize the BlockVector to contain as many blocks as there are
       * index sets given in the input argument, according to the parallel
       * distribution of the individual components described in the maps.
       *
       * If <tt>omit_zeroing_entries==false</tt>, the vector is filled with
       * zeros.
       */
      void
      reinit(const std::vector<IndexSet> &parallel_partitioning,
             const MPI_Comm               communicator         = MPI_COMM_WORLD,
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
      void
      reinit(const std::vector<IndexSet> &partitioning,
             const std::vector<IndexSet> &ghost_values,
             const MPI_Comm               communicator    = MPI_COMM_WORLD,
             const bool                   vector_writable = false);

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
      reinit(const BlockVector<Number, MemorySpace> &V,
             const bool omit_zeroing_entries = false);

      /**
       * Change the number of blocks to <tt>num_blocks</tt>. The individual
       * blocks will get initialized with zero size, so it is assumed that the
       * user resizes the individual blocks by herself in an appropriate way,
       * and calls <tt>collect_sizes</tt> afterwards.
       */
      void
      reinit(const size_type num_blocks);

      /**
       * Return if this Vector contains ghost elements.
       *
       * @see
       * @ref GlossGhostedVector "vectors with ghost elements"
       */
      bool
      has_ghost_elements() const;

      /**
       * Print to a stream.
       */
      void
      print(std::ostream      &out,
            const unsigned int precision  = 3,
            const bool         scientific = true,
            const bool         across     = true) const;
    };
  } // namespace TpetraWrappers

  /** @} */

} // namespace LinearAlgebra

/**
 * Declare dealii::LinearAlgebra::TpetraWrappers::BlockVector as distributed
 * vector.
 */
template <typename Number, typename MemorySpace>
struct is_serial_vector<
  LinearAlgebra::TpetraWrappers::BlockVector<Number, MemorySpace>>
  : std::false_type
{};

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_TPETRA

#endif // dealii_trilinos_tpetra_block_vector_h

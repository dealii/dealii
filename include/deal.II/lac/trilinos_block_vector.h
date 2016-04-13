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

#ifndef dealii__trilinos_block_vector_h
#define dealii__trilinos_block_vector_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/lac/trilinos_vector.h>
#  include <deal.II/lac/trilinos_parallel_block_vector.h>
#  include <deal.II/lac/block_indices.h>
#  include <deal.II/lac/block_vector_base.h>
#  include <deal.II/lac/exceptions.h>

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


  /**
   * An implementation of block vectors based on the vector class implemented
   * in TrilinosWrappers. While the base class provides for most of the
   * interface, this class handles the actual allocation of vectors and
   * provides functions that are specific to the underlying vector type.
   *
   * In contrast to the class MPI::BlockVector, this class is based on a
   * localized version of the vectors, which means that the whole vector is
   * stored on each processor. Note that matrix vector products with this
   * block vector class do only work in case the program is run on only one
   * processor, since the Trilinos matrices are inherently parallel.
   *
   * This class is deprecated, use TrilinosWrappers::MPI::BlockVector instead.
   *
   * @ingroup Vectors
   * @ingroup TrilinosWrappers @see
   * @ref GlossBlockLA "Block (linear algebra)"
   * @author Martin Kronbichler, 2008
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
     * entries in Input_Maps.  For this non-distributed vector, the %parallel
     * partitioning is not used, just the global size of the partitioner.
     */
    explicit BlockVector (const std::vector<Epetra_Map> &partitioner) DEAL_II_DEPRECATED;

    /**
     * Constructor. Generate a block vector with as many blocks as there are
     * entries in Input_Maps.  For this non-distributed vector, the %parallel
     * partitioning is not used, just the global size of the partitioner.
     */
    explicit BlockVector (const std::vector<IndexSet> &partitioner,
                          const MPI_Comm              &communicator = MPI_COMM_WORLD) DEAL_II_DEPRECATED;

    /**
     * Copy-Constructor. Set all the properties of the non-%parallel vector to
     * those of the given %parallel vector and import the elements.
     */
    BlockVector (const MPI::BlockVector &V) DEAL_II_DEPRECATED;

    /**
     * Copy-Constructor. Set all the properties of the vector to those of the
     * given input vector and copy the elements.
     */
    BlockVector (const BlockVector  &V) DEAL_II_DEPRECATED;

    /**
     * Creates a block vector consisting of <tt>num_blocks</tt> components,
     * but there is no content in the individual components and the user has
     * to fill appropriate data using a reinit of the blocks.
     */
    explicit BlockVector (const size_type num_blocks) DEAL_II_DEPRECATED;

    /**
     * Constructor. Set the number of blocks to <tt>n.size()</tt> and
     * initialize each block with <tt>n[i]</tt> zero elements.
     *
     * References BlockVector.reinit().
     */
    explicit BlockVector (const std::vector<size_type> &N) DEAL_II_DEPRECATED;

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
                 const InputIterator           end) DEAL_II_DEPRECATED;

    /**
     * Destructor. Clears memory
     */
    ~BlockVector ();

    /**
     * Copy operator: fill all components of the vector that are locally
     * stored with the given scalar value.
     */
    BlockVector &
    operator = (const value_type s);

    /**
     * Copy operator for a distributed Trilinos vector to a localized one.
     */
    BlockVector &
    operator = (const MPI::BlockVector &V);

    /**
     * Copy operator for arguments of the same type.
     */
    BlockVector &
    operator = (const BlockVector &V);

    /**
     * Another copy function. This one takes a deal.II block vector and copies
     * it into a TrilinosWrappers block vector. Note that the number of blocks
     * has to be the same in the vector as in the input vector. Use the
     * reinit() command for resizing the BlockVector or for changing the
     * internal structure of the block components.
     *
     * Since Trilinos only works on doubles, this function is limited to
     * accept only one possible number type in the deal.II vector.
     */
    template <typename Number>
    BlockVector &
    operator = (const ::dealii::BlockVector<Number> &V);

    /**
     * Reinitialize the BlockVector to contain as many blocks as there are
     * Epetra_Maps given in the input argument, according to the global size
     * of the individual components described in the maps. Note that the
     * resulting vector will be stored completely on each process. The
     * Epetra_Map is useful when data exchange with a distributed vector based
     * on the same Epetra_map is intended. In that case, the same communicator
     * is used for data exchange.
     *
     * If <tt>omit_zeroing_entries==false</tt>, the vector is filled with
     * zeros.
     */
    void reinit (const std::vector<Epetra_Map> &partitioning,
                 const bool                     omit_zeroing_entries = false);

    /**
     * Reinitialize the BlockVector to contain as many blocks as there are
     * index sets given in the input argument, according to the global size of
     * the individual components described in the index set, and using a given
     * MPI communicator. The MPI communicator is useful when data exchange
     * with a distributed vector based on the same initialization is intended.
     * In that case, the same communicator is used for data exchange.
     *
     * If <tt>omit_zeroing_entries==false</tt>, the vector is filled with
     * zeros.
     */
    void reinit (const std::vector<IndexSet> &partitioning,
                 const MPI_Comm              &communicator = MPI_COMM_WORLD,
                 const bool                   omit_zeroing_entries = false);

    /**
     * Reinitialize the BlockVector to contain as many blocks as there are
     * elements in the first argument, and with the respective sizes. Since no
     * distribution map is given, all vectors are local vectors.
     *
     * If <tt>omit_zeroing_entries==false</tt>, the vector is filled with
     * zeros.
     */
    void reinit (const std::vector<size_type> &N,
                 const bool                    omit_zeroing_entries=false);

    /**
     * Reinitialize the vector in the same way as the given to a distributed
     * block vector. The elements will be copied in this process.
     */
    void reinit (const MPI::BlockVector &V);

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
     * If you call reinit() on one of the blocks, then subsequent actions on
     * this object may yield unpredictable results since they may be routed to
     * the wrong block.
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
     * Exception
     */
    DeclException0 (ExcIteratorRangeDoesNotMatchVectorSize);

    /**
     * Exception
     */
    DeclException0 (ExcNonMatchingBlockVectors);

    /**
     * Exception
     */
    DeclException2 (ExcNonLocalizedMap,
                    int, int,
                    << "For the generation of a localized vector the map has "
                    << "to assign all elements to all vectors! "
                    << "local_size = global_size is a necessary condition, but"
                    << arg1 << " != " << arg2 << " was given!");

  };



  /*----------------------- Inline functions ----------------------------------*/



  inline
  BlockVector::BlockVector ()
  {}



  inline
  BlockVector::BlockVector (const std::vector<Epetra_Map> &partitioning)
  {
    reinit (partitioning);
  }



  inline
  BlockVector::BlockVector (const std::vector<IndexSet> &partitioning,
                            const MPI_Comm              &communicator)
  {
    reinit (partitioning, communicator);
  }



  inline
  BlockVector::BlockVector (const std::vector<size_type> &N)
  {
    reinit (N);
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
    for (size_type b=0; b<n.size(); ++b)
      {
        InputIterator end = start;
        std::advance (end, static_cast<size_type>(n[b]));

        for (size_type i=0; i<n[b]; ++i, ++start)
          this->block(b)(i) = *start;
      }
    Assert (start == end, ExcIteratorRangeDoesNotMatchVectorSize());
  }



  inline
  BlockVector::BlockVector (const size_type num_blocks)
  {
    reinit (num_blocks);
  }



  inline
  BlockVector::~BlockVector()
  {}



  inline
  BlockVector::BlockVector (const MPI::BlockVector &v)
  {
    reinit (v);
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


  inline
  void
  BlockVector::swap (BlockVector &v)
  {
    Assert (n_blocks() == v.n_blocks(),
            ExcDimensionMismatch(n_blocks(),v.n_blocks()));

    for (unsigned int row=0; row<n_blocks(); ++row)
      block(row).swap (v.block(row));
  }


  template <typename Number>
  BlockVector &
  BlockVector::operator = (const ::dealii::BlockVector<Number> &v)
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


  /**
   * Global function which overloads the default implementation of the C++
   * standard library which uses a temporary object. The function simply
   * exchanges the data of the two vectors.
   *
   * @relates TrilinosWrappers::BlockVector
   * @author Martin Kronbichler, 2008
   */
  inline
  void swap (BlockVector &u,
             BlockVector &v)
  {
    u.swap (v);
  }

} /* namespace TrilinosWrappers */

/*@}*/


namespace internal
{
  namespace LinearOperator
  {
    template <typename> class ReinitHelper;

    /**
     * A helper class used internally in linear_operator.h. Specialization for
     * TrilinosWrappers::BlockVector.
     */
    template<>
    class ReinitHelper<TrilinosWrappers::BlockVector>
    {
    public:
      template <typename Matrix>
      static
      void reinit_range_vector (const Matrix &matrix,
                                TrilinosWrappers::BlockVector &v,
                                bool omit_zeroing_entries)
      {
        v.reinit(matrix.range_partitioner(), omit_zeroing_entries);
      }

      template <typename Matrix>
      static
      void reinit_domain_vector(const Matrix &matrix,
                                TrilinosWrappers::BlockVector &v,
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

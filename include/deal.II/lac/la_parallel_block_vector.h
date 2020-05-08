// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_la_parallel_block_vector_h
#define dealii_la_parallel_block_vector_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/lac/block_indices.h>
#include <deal.II/lac/block_vector_base.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector_operation.h>
#include <deal.II/lac/vector_type_traits.h>

#include <cstdio>
#include <vector>

DEAL_II_NAMESPACE_OPEN


// Forward declarations
#ifndef DOXYGEN
#  ifdef DEAL_II_WITH_PETSC
namespace PETScWrappers
{
  namespace MPI
  {
    class BlockVector;
  }
} // namespace PETScWrappers
#  endif

#  ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers
{
  namespace MPI
  {
    class BlockVector;
  }
} // namespace TrilinosWrappers
#  endif
#endif

namespace LinearAlgebra
{
  namespace distributed
  {
    /*! @addtogroup Vectors
     *@{
     */


    /**
     * An implementation of block vectors based on distributed deal.II
     * vectors. While the base class provides for most of the interface, this
     * class handles the actual allocation of vectors and provides functions
     * that are specific to the underlying vector type.
     *
     * @note Instantiations for this template are provided for <tt>@<float@>
     * and @<double@></tt>; others can be generated in application programs
     * (see the section on
     * @ref Instantiations
     * in the manual).
     *
     * @see
     * @ref GlossBlockLA "Block (linear algebra)"
     * @author Katharina Kormann, Martin Kronbichler, 2011
     */
    template <typename Number>
    class BlockVector : public BlockVectorBase<Vector<Number>>,
                        public VectorSpaceVector<Number>
    {
    public:
      /**
       * The chunks size to split communication in update_ghost_values()
       * and compress() calls.
       *
       * Most common MPI implementations will get slow when too many
       * messages/requests are outstanding. Even when messages are small,
       * say 1 kB only, we should collect enough data with @p communication_block_size
       * to cover typical infiniband latencies which are around a few
       * microseconds. Sending 20 kB at a throughput of 5 GB/s takes 4
       * microseconds, so we should arrive at the bandwidth dominated regime
       * then which is good enough.
       */
      static constexpr unsigned int communication_block_size = 20;

      /**
       * Typedef the base class for simpler access to its own alias.
       */
      using BaseClass = BlockVectorBase<Vector<Number>>;

      /**
       * Typedef the type of the underlying vector.
       */
      using BlockType = typename BaseClass::BlockType;

      /**
       * Import the alias from the base class.
       */
      using value_type      = typename BaseClass::value_type;
      using real_type       = typename BaseClass::real_type;
      using pointer         = typename BaseClass::pointer;
      using const_pointer   = typename BaseClass::const_pointer;
      using reference       = typename BaseClass::reference;
      using const_reference = typename BaseClass::const_reference;
      using size_type       = typename BaseClass::size_type;
      using iterator        = typename BaseClass::iterator;
      using const_iterator  = typename BaseClass::const_iterator;

      /**
       * @name 1: Basic operations
       */
      //@{

      /**
       * Constructor. There are three ways to use this constructor. First,
       * without any arguments, it generates an object with no blocks. Given
       * one argument, it initializes <tt>num_blocks</tt> blocks, but these
       * blocks have size zero. The third variant finally initializes all
       * blocks to the same size <tt>block_size</tt>.
       *
       * Confer the other constructor further down if you intend to use blocks
       * of different sizes.
       */
      explicit BlockVector(const size_type num_blocks = 0,
                           const size_type block_size = 0);

      /**
       * Copy-Constructor. Dimension set to that of V, all components are
       * copied from V
       */
      BlockVector(const BlockVector<Number> &V);

      /**
       * Copy constructor taking a BlockVector of another data type. This will
       * fail if there is no conversion path from <tt>OtherNumber</tt> to
       * <tt>Number</tt>. Note that you may lose accuracy when copying to a
       * BlockVector with data elements with less accuracy.
       *
       * Older versions of gcc did not honor the @p explicit keyword on
       * template constructors. In such cases, it is easy to accidentally
       * write code that can be very inefficient, since the compiler starts
       * performing hidden conversions. To avoid this, this function is
       * disabled if we have detected a broken compiler during configuration.
       */
      template <typename OtherNumber>
      explicit BlockVector(const BlockVector<OtherNumber> &v);

      /**
       * Constructor. Set the number of blocks to <tt>block_sizes.size()</tt>
       * and initialize each block with <tt>block_sizes[i]</tt> zero elements.
       */
      BlockVector(const std::vector<size_type> &block_sizes);

      /**
       * Construct a block vector with an IndexSet for the local range and
       * ghost entries for each block.
       */
      BlockVector(const std::vector<IndexSet> &local_ranges,
                  const std::vector<IndexSet> &ghost_indices,
                  const MPI_Comm               communicator);

      /**
       * Same as above but the ghost indices are assumed to be empty.
       */
      BlockVector(const std::vector<IndexSet> &local_ranges,
                  const MPI_Comm               communicator);

      /**
       * Destructor.
       *
       * @note We need to explicitly provide a destructor, otherwise the
       *   linker may think it is unused and discards it, although required
       *   in a different section. The Intel compiler is prone to this
       *   behavior.
       */
      virtual ~BlockVector() override = default;

      /**
       * Copy operator: fill all components of the vector with the given
       * scalar value.
       */
      virtual BlockVector &
      operator=(const value_type s) override;

      /**
       * Copy operator for arguments of the same type. Resize the present
       * vector if necessary.
       */
      BlockVector &
      operator=(const BlockVector &V);

      /**
       * Copy operator for template arguments of different types. Resize the
       * present vector if necessary.
       */
      template <class Number2>
      BlockVector &
      operator=(const BlockVector<Number2> &V);

      /**
       * Copy a regular vector into a block vector.
       */
      BlockVector &
      operator=(const Vector<Number> &V);

#ifdef DEAL_II_WITH_PETSC
      /**
       * Copy the content of a PETSc vector into the calling vector. This
       * function assumes that the vectors layouts have already been
       * initialized to match.
       *
       * This operator is only available if deal.II was configured with PETSc.
       */
      BlockVector<Number> &
      operator=(const PETScWrappers::MPI::BlockVector &petsc_vec);
#endif

#ifdef DEAL_II_WITH_TRILINOS
      /**
       * Copy the content of a Trilinos vector into the calling vector. This
       * function assumes that the vectors layouts have already been
       * initialized to match.
       *
       * This operator is only available if deal.II was configured with
       * Trilinos.
       */
      BlockVector<Number> &
      operator=(const TrilinosWrappers::MPI::BlockVector &trilinos_vec);
#endif

      /**
       * Reinitialize the BlockVector to contain <tt>num_blocks</tt> blocks of
       * size <tt>block_size</tt> each.
       *
       * If the second argument is left at its default value, then the block
       * vector allocates the specified number of blocks but leaves them at
       * zero size. You then need to later reinitialize the individual blocks,
       * and call collect_sizes() to update the block system's knowledge of
       * its individual block's sizes.
       *
       * If <tt>omit_zeroing_entries==false</tt>, the vector is filled with
       * zeros.
       */
      void
      reinit(const size_type num_blocks,
             const size_type block_size           = 0,
             const bool      omit_zeroing_entries = false);

      /**
       * Reinitialize the BlockVector such that it contains
       * <tt>block_sizes.size()</tt> blocks. Each block is reinitialized to
       * dimension <tt>block_sizes[i]</tt>.
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
       * sizes. If you call reinit() on one of the blocks, then subsequent
       * actions on this object may yield unpredictable results since they may
       * be routed to the wrong block.
       */
      void
      reinit(const std::vector<size_type> &N,
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
       * sizes. If you call reinit() of one of the blocks, then subsequent
       * actions of this object may yield unpredictable results since they may
       * be routed to the wrong block.
       */
      template <typename Number2>
      void
      reinit(const BlockVector<Number2> &V,
             const bool                  omit_zeroing_entries = false);

      /**
       * This function copies the data that has accumulated in the data buffer
       * for ghost indices to the owning processor. For the meaning of the
       * argument @p operation, see the entry on
       * @ref GlossCompress "Compressing distributed vectors and matrices"
       * in the glossary.
       *
       * There are two variants for this function. If called with argument @p
       * VectorOperation::add adds all the data accumulated in ghost elements
       * to the respective elements on the owning processor and clears the
       * ghost array afterwards. If called with argument @p
       * VectorOperation::insert, a set operation is performed. Since setting
       * elements in a vector with ghost elements is ambiguous (as one can set
       * both the element on the ghost site as well as the owning site), this
       * operation makes the assumption that all data is set correctly on the
       * owning processor. Upon call of compress(VectorOperation::insert), all
       * ghost entries are therefore simply zeroed out (using
       * zero_ghost_values()). In debug mode, a check is performed that makes
       * sure that the data set is actually consistent between processors,
       * i.e., whenever a non-zero ghost element is found, it is compared to
       * the value on the owning processor and an exception is thrown if these
       * elements do not agree.
       *
       */
      virtual void
      compress(::dealii::VectorOperation::values operation) override;

      /**
       * Fills the data field for ghost indices with the values stored in the
       * respective positions of the owning processor. This function is needed
       * before reading from ghosts. The function is @p const even though
       * ghost data is changed. This is needed to allow functions with a @p
       * const vector to perform the data exchange without creating
       * temporaries.
       */
      void
      update_ghost_values() const;

      /**
       * This method zeros the entries on ghost dofs, but does not touch
       * locally owned DoFs.
       *
       * After calling this method, read access to ghost elements of the
       * vector is forbidden and an exception is thrown. Only write access to
       * ghost elements is allowed in this state.
       */
      void
      zero_out_ghosts() const;

      /**
       * Return if this Vector contains ghost elements.
       */
      bool
      has_ghost_elements() const;

      /**
       * This is a collective add operation that adds a whole set of values
       * stored in @p values to the vector components specified by @p indices.
       */
      template <typename OtherNumber>
      void
      add(const std::vector<size_type> &       indices,
          const ::dealii::Vector<OtherNumber> &values);

      /**
       * Scaling and simple vector addition, i.e.  <tt>*this =
       * s*(*this)+V</tt>.
       */
      void
      sadd(const Number s, const BlockVector<Number> &V);

      /**
       * Return whether the vector contains only elements with value zero.
       * This function is mainly for internal consistency checks and should
       * seldom be used when not in debug mode since it uses quite some time.
       */
      virtual bool
      all_zero() const override;

      /**
       * Compute the mean value of all the entries in the vector.
       */
      virtual Number
      mean_value() const override;

      /**
       * $l_p$-norm of the vector. The pth root of the sum of the pth powers
       * of the absolute values of the elements.
       */
      real_type
      lp_norm(const real_type p) const;

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
      swap(BlockVector<Number> &v);
      //@}

      /**
       * @name 2: Implementation of VectorSpaceVector
       */
      //@{

      /**
       * Change the dimension to that of the vector V. The elements of V are not
       * copied.
       */
      virtual void
      reinit(const VectorSpaceVector<Number> &V,
             const bool omit_zeroing_entries = false) override;

      /**
       * Multiply the entire vector by a fixed factor.
       */
      virtual BlockVector<Number> &
      operator*=(const Number factor) override;

      /**
       * Divide the entire vector by a fixed factor.
       */
      virtual BlockVector<Number> &
      operator/=(const Number factor) override;

      /**
       * Add the vector @p V to the present one.
       */
      virtual BlockVector<Number> &
      operator+=(const VectorSpaceVector<Number> &V) override;

      /**
       * Subtract the vector @p V from the present one.
       */
      virtual BlockVector<Number> &
      operator-=(const VectorSpaceVector<Number> &V) override;

      /**
       * Import all the elements present in the vector's IndexSet from the input
       * vector @p V. VectorOperation::values @p operation is used to decide if
       * the elements in @p V should be added to the current vector or replace the
       * current elements. The last parameter can be used if the same
       * communication pattern is used multiple times. This can be used to
       * improve performance.
       */
      virtual void
      import(
        const LinearAlgebra::ReadWriteVector<Number> &  V,
        VectorOperation::values                         operation,
        std::shared_ptr<const CommunicationPatternBase> communication_pattern =
          std::shared_ptr<const CommunicationPatternBase>()) override;

      /**
       * Return the scalar product of two vectors.
       */
      virtual Number
      operator*(const VectorSpaceVector<Number> &V) const override;

      /**
       * Calculate the scalar product between each block of this vector and @p V
       * and store the result in a full matrix @p matrix. This function
       * computes the result by forming $A_{ij}=U_i \cdot V_j$ where $U_i$
       * and $V_j$ indicate the $i$th block (not element!) of $U$ and the
       * $j$th block of $V$, respectively. If @p symmetric is
       * <code>true</code>, it is assumed that inner product results in a
       * square symmetric matrix and almost half of the scalar products can be
       * avoided.
       *
       * Obviously, this function can only be used if all blocks of both vectors
       * are of the same size.
       *
       * @note Internally, a single global reduction will be called to
       * accumulate scalar product between locally owned degrees of freedom.
       */
      template <typename FullMatrixType>
      void
      multivector_inner_product(FullMatrixType &           matrix,
                                const BlockVector<Number> &V,
                                const bool symmetric = false) const;

      /**
       * Calculate the scalar product between each block of this vector and @p V
       * using a metric tensor @p matrix. This function
       * computes the result of $ \sum_{ij} A^{ij} U_i \cdot V_j$ where $U_i$
       * and $V_j$ indicate the $i$th block (not element) of $U$ and the
       * $j$th block of $V$, respectively. If @p symmetric is
       * <code>true</code>, it is assumed that $U_i \cdot V_j$ and $A^{ij}$ are
       * symmetric matrices and almost half of the scalar products can be
       * avoided.
       *
       * Obviously, this function can only be used if all blocks of both vectors
       * are of the same size.
       *
       * @note Internally, a single global reduction will be called to
       * accumulate the scalar product between locally owned degrees of freedom.
       */
      template <typename FullMatrixType>
      Number
      multivector_inner_product_with_metric(const FullMatrixType &     matrix,
                                            const BlockVector<Number> &V,
                                            const bool symmetric = false) const;

      /**
       * Set each block of this vector as follows:
       * $V^i = s V^i + b \sum_{j} U_j A^{ji}$ where $V^i$
       * and $U_j$ indicate the $i$th block (not element) of $V$ and the
       * $j$th block of $U$, respectively.
       *
       * Obviously, this function can only be used if all blocks of both vectors
       * are of the same size.
       */
      template <typename FullMatrixType>
      void
      mmult(BlockVector<Number> & V,
            const FullMatrixType &matrix,
            const Number          s = Number(0.),
            const Number          b = Number(1.)) const;

      /**
       * Add @p a to all components. Note that @p a is a scalar not a vector.
       */
      virtual void
      add(const Number a) override;

      /**
       * Simple addition of a multiple of a vector, i.e. <tt>*this += a*V</tt>.
       */
      virtual void
      add(const Number a, const VectorSpaceVector<Number> &V) override;

      /**
       * Multiple addition of scaled vectors, i.e. <tt>*this += a*V+b*W</tt>.
       */
      virtual void
      add(const Number                     a,
          const VectorSpaceVector<Number> &V,
          const Number                     b,
          const VectorSpaceVector<Number> &W) override;

      /**
       * A collective add operation: This function adds a whole set of values
       * stored in @p values to the vector components specified by @p indices.
       */
      virtual void
      add(const std::vector<size_type> &indices,
          const std::vector<Number> &   values);

      /**
       * Scaling and simple addition of a multiple of a vector, i.e. <tt>*this =
       * s*(*this)+a*V</tt>.
       */
      virtual void
      sadd(const Number                     s,
           const Number                     a,
           const VectorSpaceVector<Number> &V) override;

      /**
       * Scale each element of this vector by the corresponding element in the
       * argument. This function is mostly meant to simulate multiplication (and
       * immediate re-assignment) by a diagonal scaling matrix.
       */
      virtual void
      scale(const VectorSpaceVector<Number> &scaling_factors) override;

      /**
       * Assignment <tt>*this = a*V</tt>.
       */
      virtual void
      equ(const Number a, const VectorSpaceVector<Number> &V) override;

      /**
       * Return the l<sub>1</sub> norm of the vector (i.e., the sum of the
       * absolute values of all entries among all processors).
       */
      virtual real_type
      l1_norm() const override;

      /**
       * Return the $l_2$ norm of the vector (i.e., the square root of
       * the sum of the square of all entries among all processors).
       */
      virtual real_type
      l2_norm() const override;

      /**
       * Return the square of the $l_2$ norm of the vector.
       */
      real_type
      norm_sqr() const;

      /**
       * Return the maximum norm of the vector (i.e., the maximum absolute value
       * among all entries and among all processors).
       */
      virtual real_type
      linfty_norm() const override;

      /**
       * Perform a combined operation of a vector addition and a subsequent
       * inner product, returning the value of the inner product. In other
       * words, the result of this function is the same as if the user called
       * @code
       * this->add(a, V);
       * return_value = *this * W;
       * @endcode
       *
       * The reason this function exists is that this operation involves less
       * memory transfer than calling the two functions separately. This method
       * only needs to load three vectors, @p this, @p V, @p W, whereas calling
       * separate methods means to load the calling vector @p this twice. Since
       * most vector operations are memory transfer limited, this reduces the
       * time by 25\% (or 50\% if @p W equals @p this).
       *
       * For complex-valued vectors, the scalar product in the second step is
       * implemented as
       * $\left<v,w\right>=\sum_i v_i \bar{w_i}$.
       */
      virtual Number
      add_and_dot(const Number                     a,
                  const VectorSpaceVector<Number> &V,
                  const VectorSpaceVector<Number> &W) override;

      /**
       * Return the global size of the vector, equal to the sum of the number of
       * locally owned indices among all processors.
       */
      virtual size_type
      size() const override;

      /**
       * Return an index set that describes which elements of this vector are
       * owned by the current processor. As a consequence, the index sets
       * returned on different processors if this is a distributed vector will
       * form disjoint sets that add up to the complete index set. Obviously, if
       * a vector is created on only one processor, then the result would
       * satisfy
       * @code
       *  vec.locally_owned_elements() == complete_index_set(vec.size())
       * @endcode
       */
      virtual dealii::IndexSet
      locally_owned_elements() const override;

      /**
       * Print the vector to the output stream @p out.
       */
      virtual void
      print(std::ostream &     out,
            const unsigned int precision  = 3,
            const bool         scientific = true,
            const bool         across     = true) const override;

      /**
       * Return the memory consumption of this class in bytes.
       */
      virtual std::size_t
      memory_consumption() const override;
      //@}

      /**
       * @addtogroup Exceptions
       * @{
       */

      /**
       * Attempt to perform an operation between two incompatible vector types.
       *
       * @ingroup Exceptions
       */
      DeclException0(ExcVectorTypeNotCompatible);

      /**
       * Exception
       */
      DeclException0(ExcIteratorRangeDoesNotMatchVectorSize);
      //@}
    };

    /*@}*/

  } // end of namespace distributed

} // end of namespace LinearAlgebra


/**
 * Global function which overloads the default implementation of the C++
 * standard library which uses a temporary object. The function simply
 * exchanges the data of the two vectors.
 *
 * @relatesalso BlockVector
 * @author Katharina Kormann, Martin Kronbichler, 2011
 */
template <typename Number>
inline void
swap(LinearAlgebra::distributed::BlockVector<Number> &u,
     LinearAlgebra::distributed::BlockVector<Number> &v)
{
  u.swap(v);
}


/**
 * Declare dealii::LinearAlgebra::distributed::BlockVector as distributed
 * vector.
 *
 * @author Uwe Koecher, 2017
 */
template <typename Number>
struct is_serial_vector<LinearAlgebra::distributed::BlockVector<Number>>
  : std::false_type
{};


DEAL_II_NAMESPACE_CLOSE

#ifdef DEAL_II_MSVC
#  include <deal.II/lac/la_parallel_block_vector.templates.h>
#endif

#endif

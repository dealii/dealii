// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

#ifndef __deal2__precondition_block_h
#define __deal2__precondition_block_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/lac/precondition_block_base.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

template<class MATRIX, typename inverse_type>
class PreconditionBlockJacobi;

/*! @addtogroup Preconditioners
 *@{
 */


/**
 * Base class for actual block preconditioners. This class assumes the
 * <tt>MATRIX</tt> consisting of invertible blocks of @p blocksize on
 * the diagonal and provides the inversion of the diagonal blocks of
 * the matrix. It is not necessary for this class that the matrix be
 * block diagonal; rather, it applies to
 * matrices of arbitrary structure with the minimal property of having
 * invertible blocks on the diagonal. Still the matrix must have
 * access to single matrix entries. Therefore, BlockMatrixArray and similar
 * classes are not a possible matrix class template arguments.
 *
 * The block matrix structure used by this class is given, e.g., for
 * the DG method for the transport equation. For a downstream
 * numbering the matrices even have got a block lower left matrix
 * structure, i.e. the matrices are empty above the diagonal blocks.
 *
 * @note This class is intended to be used for matrices whose structure
 * is given by local contributions from disjoint cells, such as for DG
 * methods. It is not intended for problems where the block structure
 * results from different physical variables such as in the Stokes
 * equations considered in step-22.
 *
 * For all matrices that are empty above and below the diagonal blocks
 * (i.e. for all block diagonal matrices) the @p BlockJacobi
 * preconditioner is a direct solver. For all matrices that are empty
 * only above the diagonal blocks (e.g. the matrices one gets by the
 * DG method with downstream numbering) @p BlockSOR is a direct
 * solver.
 *
 * This first implementation of the @p PreconditionBlock assumes the
 * matrix has blocks each of the same block size. Varying block sizes
 * within the matrix must still be implemented if needed.
 *
 * The first template parameter denotes the type of number
 * representation in the sparse matrix, the second denotes the type of
 * number representation in which the inverted diagonal block matrices
 * are stored within this class by <tt>invert_diagblocks()</tt>. If you
 * don't want to use the block inversion as an exact solver, but
 * rather as a preconditioner, you may probably want to store the
 * inverted blocks with less accuracy than the original matrix; for
 * example, <tt>number==double, inverse_type=float</tt> might be a viable
 * choice.
 *
 * @see @ref GlossBlockLA "Block (linear algebra)"
 * @author Ralf Hartmann, Guido Kanschat
 * @date 1999, 2000, 2010
 */
template<class MATRIX, typename inverse_type = typename MATRIX::value_type>
class PreconditionBlock
  : public virtual Subscriptor,
    protected PreconditionBlockBase<inverse_type>
{
private:
  /**
   * Define number type of matrix.
   */
  typedef typename MATRIX::value_type number;

  /**
   * Value type for inverse matrices.
   */
  typedef inverse_type value_type;

public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Parameters for block preconditioners.
   */
  class AdditionalData
  {
  public:
    /**
     * Constructor. Block size
     * must be given since there
     * is no reasonable default
     * parameter.
     */
    AdditionalData (const size_type block_size,
                    const double relaxation = 1.,
                    const bool invert_diagonal = true,
                    const bool same_diagonal = false);

    /**
     * Relaxation parameter.
     */
    double relaxation;

    /**
     * Block size.
     */
    size_type block_size;

    /**
     * Invert diagonal during initialization.
     */
    bool invert_diagonal;

    /**
     * Assume all diagonal blocks
     * are equal to save memory.
     */
    bool same_diagonal;
    /**
     * Choose the inversion
     * method for the blocks.
     */
    typename PreconditionBlockBase<inverse_type>::Inversion inversion;

    /**
     * The if #inversion is SVD,
     * the threshold below which
     * a singular value will be
     * considered zero and thus
     * not inverted. This
     * parameter is used in the
     * call to LAPACKFullMatrix::compute_inverse_svd().
     */
    double threshold;
  };


  /**
   * Constructor.
   */
  PreconditionBlock(bool store_diagonals = false);

  /**
   * Destructor.
   */
  ~PreconditionBlock();

  /**
   * Initialize matrix and block
   * size.  We store the matrix and
   * the block size in the
   * preconditioner object. In a
   * second step, the inverses of
   * the diagonal blocks may be
   * computed.
   *
   * Additionally, a relaxation
   * parameter for derived classes
   * may be provided.
   */
  void initialize (const MATRIX &A,
                   const AdditionalData parameters);
protected:
  /**
   * Initialize matrix and block
   * size for permuted
   * preconditioning. Additionally
   * to the parameters of the other
   * initalize() function, we hand
   * over two index vectors with
   * the permutation and its
   * inverse. For the meaning of
   * these vectors see
   * PreconditionBlockSOR.
   *
   * In a second step, the inverses
   * of the diagonal blocks may be
   * computed. Make sure you use
   * invert_permuted_diagblocks()
   * to yield consistent data.
   *
   * Additionally, a relaxation
   * parameter for derived classes
   * may be provided.
   */
  void initialize (const MATRIX &A,
                   const std::vector<size_type> &permutation,
                   const std::vector<size_type> &inverse_permutation,
                   const AdditionalData parameters);

  /**
   * Set either the permutation of
   * rows or the permutation of
   * blocks, depending on the size
   * of the vector.
   *
   * If the size of the permutation
   * vectors is equal to the
   * dimension of the linear
   * system, it is assumed that
   * rows are permuted
   * individually. In this case,
   * set_permutation() must be
   * called before initialize(),
   * since the diagonal blocks are
   * built from the permuted
   * entries of the matrix.
   *
   * If the size of the permutation
   * vector is not equal to the
   * dimension of the system, the
   * diagonal blocks are computed
   * from the unpermuted
   * entries. Instead, the
   * relaxation methods step() and
   * Tstep() are executed applying
   * the blocks in the order given
   * by the permutation
   * vector. They will throw an
   * exception if length of this
   * vector is not equal to the
   * number of blocks.
   *
   * @note Permutation of blocks
   * can only be applied to the
   * relaxation operators step()
   * and Tstep(), not to the
   * preconditioning operators
   * vmult() and Tvmult().
   *
   * @note It is safe to call
   * set_permutation() before
   * initialize(), while the other
   * order is only admissible for
   * block permutation.
   */
  void set_permutation(const std::vector<size_type> &permutation,
                       const std::vector<size_type> &inverse_permutation);

  /**
   * Replacement of
   * invert_diagblocks() for
   * permuted preconditioning.
   */
  void invert_permuted_diagblocks(
    const std::vector<size_type> &permutation,
    const std::vector<size_type> &inverse_permutation);
public:
  /**
   * Deletes the inverse diagonal
   * block matrices if existent,
   * sets the blocksize to 0, hence
   * leaves the class in the state
   * that it had directly after
   * calling the constructor.
   */
  void clear();

  /**
   * Checks whether the object is empty.
   */
  bool empty () const;

  /**
   * Read-only access to entries.
   * This function is only possible
   * if the inverse diagonal blocks
   * are stored.
   */
  value_type el(size_type i,
                size_type j) const;

  /**
   * Stores the inverse of the
   * diagonal blocks in
   * @p inverse. This costs some
   * additional memory - for DG
   * methods about 1/3 (for double
   * inverses) or 1/6 (for float
   * inverses) of that used for the
   * matrix - but it makes the
   * preconditioning much faster.
   *
   * It is not allowed to call this
   * function twice (will produce
   * an error) before a call of
   * <tt>clear(...)</tt>  because at the
   * second time there already
   * exist the inverse matrices.
   *
   * After this function is called,
   * the lock on the matrix given
   * through the @p use_matrix
   * function is released, i.e. you
   * may overwrite of delete it.
   * You may want to do this in
   * case you use this matrix to
   * precondition another matrix.
   */
  void invert_diagblocks();

  /**
   * Perform one block relaxation
   * step in forward numbering.
   *
   * Depending on the arguments @p
   * dst and @p pref, this performs
   * an SOR step (both reference
   * the same vector) of a Jacobi
   * step (botha different
   * vectors). For the Jacobi step,
   * the calling function must copy
   * @p dst to @p pref after this.
   *
   * @note If a permutation is set,
   * it is automatically honored by
   * this function.
   */
  template <typename number2>
  void forward_step (
    Vector<number2>       &dst,
    const Vector<number2> &prev,
    const Vector<number2> &src,
    const bool transpose_diagonal) const;

  /**
   * Perform one block relaxation
   * step in backward numbering.
   *
   * Depending on the arguments @p
   * dst and @p pref, this performs
   * an SOR step (both reference
   * the same vector) of a Jacobi
   * step (botha different
   * vectors). For the Jacobi step,
   * the calling function must copy
   * @p dst to @p pref after this.
   *
   * @note If a permutation is set,
   * it is automatically honored by
   * this function.
   */
  template <typename number2>
  void backward_step (
    Vector<number2>       &dst,
    const Vector<number2> &prev,
    const Vector<number2> &src,
    const bool transpose_diagonal) const;


  /**
   * Return the size of the blocks.
   */
  size_type block_size () const;

  /**
   * @deprecated Use size()
   * instead.
   *
   * The number of blocks of the
   * matrix.
   */
  unsigned int n_blocks() const DEAL_II_DEPRECATED;

  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   */
  std::size_t memory_consumption () const;

  /** @addtogroup Exceptions
   * @{ */

  /**
   * For non-overlapping block
   * preconditioners, the block
   * size must divide the matrix
   * size. If not, this exception
   * gets thrown.
   */
  DeclException2 (ExcWrongBlockSize,
                  int, int,
                  << "The blocksize " << arg1
                  << " and the size of the matrix " << arg2
                  << " do not match.");

  /**
   * Exception
   */
  DeclException0 (ExcInverseMatricesAlreadyExist);

  //@}

protected:
  /**
   * Size of the blocks. Each
   * diagonal block is assumed to
   * be of the same size.
   */
  size_type blocksize;

  /**
   * Pointer to the matrix. Make
   * sure that the matrix exists as
   * long as this class needs it,
   * i.e. until calling
   * @p invert_diagblocks, or (if
   * the inverse matrices should
   * not be stored) until the last
   * call of the preconditoining
   * @p vmult function of the
   * derived classes.
   */
  SmartPointer<const MATRIX,PreconditionBlock<MATRIX,inverse_type> > A;
  /**
   * Relaxation parameter to be
   * used by derived classes.
   */
  double relaxation;

  /**
   * The permutation vector
   */
  std::vector<size_type> permutation;

  /**
   * The inverse permutation vector
   */
  std::vector<size_type> inverse_permutation;

  /**
   * Flag for diagonal compression.
   * @ref set_same_diagonal()
   */
};



/**
 * Block Jacobi preconditioning.
 * See PreconditionBlock for requirements on the matrix.
 *
 * @note Instantiations for this template are provided for <tt>@<float@> and
 * @<double@></tt>; others can be generated in application programs (see the
 * section on @ref Instantiations in the manual).
 *
 * @author Ralf Hartmann, Guido Kanschat, 1999, 2000, 2003
 */
template<class MATRIX, typename inverse_type = typename MATRIX::value_type>
class PreconditionBlockJacobi : public virtual Subscriptor,
  private PreconditionBlock<MATRIX, inverse_type>
{
private:
  /**
   * Define number type of matrix.
   */
  typedef typename MATRIX::value_type number;

public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * STL conforming iterator.
   */
  class const_iterator
  {
  private:
    /**
     * Accessor class for iterators
     */
    class Accessor
    {
    public:
      /**
       * Constructor. Since we use
       * accessors only for read
       * access, a const matrix
       * pointer is sufficient.
       */
      Accessor (const PreconditionBlockJacobi<MATRIX, inverse_type> *matrix,
                const size_type row);

      /**
       * Row number of the element
       * represented by this
       * object.
       */
      size_type row() const;

      /**
       * Column number of the
       * element represented by
       * this object.
       */
      size_type column() const;

      /**
       * Value of this matrix entry.
       */
      inverse_type value() const;

    protected:
      /**
       * The matrix accessed.
       */
      const PreconditionBlockJacobi<MATRIX, inverse_type> *matrix;

      /**
       * Save block size here
       * for further reference.
       */
      size_type bs;

      /**
       * Current block number.
       */
      size_type a_block;

      /**
       * Iterator inside block.
       */
      typename FullMatrix<inverse_type>::const_iterator b_iterator;

      /**
       * End of current block.
       */
      typename FullMatrix<inverse_type>::const_iterator b_end;

      /**
       * Make enclosing class a
       * friend.
       */
      friend class const_iterator;
    };

  public:
    /**
     * Constructor.
     */
    const_iterator(const PreconditionBlockJacobi<MATRIX, inverse_type> *matrix,
                   const size_type row);

    /**
     * Prefix increment.
     */
    const_iterator &operator++ ();

    /**
     * Postfix increment.
     */
    const_iterator &operator++ (int);

    /**
     * Dereferencing operator.
     */
    const Accessor &operator* () const;

    /**
     * Dereferencing operator.
     */
    const Accessor *operator-> () const;

    /**
     * Comparison. True, if
     * both iterators point to
     * the same matrix
     * position.
     */
    bool operator == (const const_iterator &) const;
    /**
     * Inverse of <tt>==</tt>.
     */
    bool operator != (const const_iterator &) const;

    /**
     * Comparison
     * operator. Result is true
     * if either the first row
     * number is smaller or if
     * the row numbers are
     * equal and the first
     * index is smaller.
     */
    bool operator < (const const_iterator &) const;

  private:
    /**
     * Store an object of the
     * accessor class.
     */
    Accessor accessor;
  };

  /**
   * import functions from private base class
   */
  using PreconditionBlock<MATRIX, inverse_type>::initialize;
  using PreconditionBlock<MATRIX, inverse_type>::clear;
  using PreconditionBlock<MATRIX, inverse_type>::empty;
  using PreconditionBlock<MATRIX, inverse_type>::el;
  using PreconditionBlock<MATRIX, inverse_type>::set_same_diagonal;
  using PreconditionBlock<MATRIX, inverse_type>::invert_diagblocks;
  using PreconditionBlock<MATRIX, inverse_type>::block_size;
  using PreconditionBlockBase<inverse_type>::size;
  using PreconditionBlockBase<inverse_type>::inverse;
  using PreconditionBlockBase<inverse_type>::inverse_householder;
  using PreconditionBlockBase<inverse_type>::inverse_svd;
  using PreconditionBlockBase<inverse_type>::log_statistics;
  /**
   * @deprecated Use size() instead
   */
  using PreconditionBlock<MATRIX, inverse_type>::n_blocks;
  using PreconditionBlock<MATRIX, inverse_type>::set_permutation;

  /**
   * Execute block Jacobi
   * preconditioning.
   *
   * This function will
   * automatically use the inverse
   * matrices if they exist, if not
   * then BlockJacobi will need
   * much time inverting the
   * diagonal block matrices in
   * each preconditioning step.
   */
  template <typename number2>
  void vmult (Vector<number2> &, const Vector<number2> &) const;

  /**
   * Same as @p vmult, since Jacobi is symmetric.
   */
  template <typename number2>
  void Tvmult (Vector<number2> &, const Vector<number2> &) const;
  /**
   * Execute block Jacobi
   * preconditioning, adding to @p dst.
   *
   * This function will
   * automatically use the inverse
   * matrices if they exist, if not
   * then BlockJacobi will need
   * much time inverting the
   * diagonal block matrices in
   * each preconditioning step.
   */
  template <typename number2>
  void vmult_add (Vector<number2> &, const Vector<number2> &) const;

  /**
   * Same as @p vmult_add, since
   * Jacobi is symmetric.
   */
  template <typename number2>
  void Tvmult_add (Vector<number2> &, const Vector<number2> &) const;

  /**
   * Perform one step of the Jacobi
   * iteration.
   */
  template <typename number2>
  void step (Vector<number2> &dst, const Vector<number2> &rhs) const;

  /**
   * Perform one step of the Jacobi
   * iteration.
   */
  template <typename number2>
  void Tstep (Vector<number2> &dst, const Vector<number2> &rhs) const;

  /**
   * STL-like iterator with the
   * first entry.
   */
  const_iterator begin () const;

  /**
   * Final iterator.
   */
  const_iterator end () const;

  /**
   * STL-like iterator with the
   * first entry of row @p r.
   */
  const_iterator begin (const size_type r) const;

  /**
   * Final iterator of row @p r.
   */
  const_iterator end (const size_type r) const;


private:
  /**
   * Actual implementation of the
   * preconditioner.
   *
   * Depending on @p adding, the
   * result of preconditioning is
   * added to the destination vector.
   */
  template <typename number2>
  void do_vmult (Vector<number2> &,
                 const Vector<number2> &,
                 bool adding) const;

  friend class Accessor;
  friend class const_iterator;
};



/**
 * Block SOR preconditioning.
 *
 * The functions @p vmult and @p Tvmult execute a (transposed)
 * block-SOR step, based on the blocks in PreconditionBlock. The
 * elements outside the diagonal blocks may be distributed
 * arbitrarily.
 *
 * See PreconditionBlock for requirements on the matrix. The blocks
 * used in this class must be contiguous and non-overlapping. An
 * overlapping Schwarz relaxation method can be found in
 * RelaxationBlockSOR; that class does not offer preconditioning,
 * though.
 *
 * <h3>Permutations</h3>
 *
 * Optionally, the entries of the source vector can be treated in the
 * order of the indices in the permutation vector set by
 * #set_permutation (or the opposite order for Tvmult()). The inverse
 * permutation is used for storing elements back into this
 * vector. This functionality is automatically enabled after a call to
 * set_permutation() with vectors of nonzero size.
 *
 * @note The diagonal blocks, like the matrix, are not permuted!
 * Therefore, the permutation vector can only swap whole blocks. It
 * may not change the order inside blocks or swap single indices
 * between blocks.
 *
 * <h3>Instantiations</h3>
 *
 * @note Instantiations for this template are provided for <tt>@<float@> and
 * @<double@></tt>; others can be generated in application programs (see the
 * section on @ref Instantiations in the manual).
 *
 * @author Ralf Hartmann, Guido Kanschat, 1999, 2000, 2001, 2002, 2003
 */
template<class MATRIX, typename inverse_type = typename MATRIX::value_type>
class PreconditionBlockSOR : public virtual Subscriptor,
  protected PreconditionBlock<MATRIX, inverse_type>
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Default constructor.
   */
  PreconditionBlockSOR();

  /**
   * Define number type of matrix.
   */
  typedef typename MATRIX::value_type number;

  /**
   * import types and functions from protected base class.
   */
  using typename PreconditionBlock<MATRIX,inverse_type>::AdditionalData;
  using PreconditionBlock<MATRIX, inverse_type>::initialize;
  using PreconditionBlock<MATRIX, inverse_type>::clear;
  using PreconditionBlock<MATRIX, inverse_type>::empty;
  using PreconditionBlockBase<inverse_type>::size;
  using PreconditionBlockBase<inverse_type>::inverse;
  using PreconditionBlockBase<inverse_type>::inverse_householder;
  using PreconditionBlockBase<inverse_type>::inverse_svd;
  using PreconditionBlock<MATRIX, inverse_type>::el;
  using PreconditionBlock<MATRIX, inverse_type>::set_same_diagonal;
  using PreconditionBlock<MATRIX, inverse_type>::invert_diagblocks;
  using PreconditionBlock<MATRIX, inverse_type>::set_permutation;
  using PreconditionBlockBase<inverse_type>::log_statistics;

  /**
   * Execute block SOR
   * preconditioning.
   *
   * This function will
   * automatically use the inverse
   * matrices if they exist, if not
   * then BlockSOR will waste much
   * time inverting the diagonal
   * block matrices in each
   * preconditioning step.
   *
   * For matrices which are empty
   * above the diagonal blocks
   * BlockSOR is a direct solver.
   */
  template <typename number2>
  void vmult (Vector<number2> &, const Vector<number2> &) const;

  /**
   * Execute block SOR
   * preconditioning.
   *
   * Warning: this function
   * performs normal @p vmult
   * without adding. The reason for
   * its existence is that
   * BlockMatrixArray
   * requires the adding version by
   * default. On the other hand,
   * adding requires an additional
   * auxiliary vector, which is not
   * desirable.
   *
   * @see vmult
   */
  template <typename number2>
  void vmult_add (Vector<number2> &, const Vector<number2> &) const;

  /**
   * Backward application of vmult().
   *
   * In the current implementation,
   * this is not the transpose of
   * vmult(). It is a
   * transposed Gauss-Seidel
   * algorithm applied to the whole
   * matrix, but the diagonal
   * blocks being inverted are not
   * transposed. Therefore, it is
   * the transposed, if the
   * diagonal blocks are symmetric.
   */
  template <typename number2>
  void Tvmult (Vector<number2> &, const Vector<number2> &) const;

  /**
   * Execute backward block SOR
   * preconditioning.
   *
   * Warning: this function
   * performs normal @p vmult
   * without adding. The reason for
   * its existence is that
   * BlockMatrixArray
   * requires the adding version by
   * default. On the other hand,
   * adding requires an additional
   * auxiliary vector, which is not
   * desirable.
   *
   * @see vmult
   */
  template <typename number2>
  void Tvmult_add (Vector<number2> &, const Vector<number2> &) const;

  /**
   * Perform one step of the SOR
   * iteration.
   */
  template <typename number2>
  void step (Vector<number2> &dst, const Vector<number2> &rhs) const;

  /**
   * Perform one step of the
   * transposed SOR iteration.
   */
  template <typename number2>
  void Tstep (Vector<number2> &dst, const Vector<number2> &rhs) const;

protected:
  /**
   * Constructor to be used by
   * PreconditionBlockSSOR.
   */
  PreconditionBlockSOR(bool store);

  /**
   * Implementation of the forward
   * substitution loop called by
   * vmult() and vmult_add().
   *
   * If a #permutation is set by
   * set_permutation(), it will
   * automatically be honored by
   * this function.
   *
   * The parameter @p adding does
   * not have any function, yet.
   */
  template <typename number2>
  void forward (Vector<number2> &,
                const Vector<number2> &,
                const bool transpose_diagonal,
                const bool adding) const;

  /**
   * Implementation of the backward
   * substitution loop called by
   * Tvmult() and Tvmult_add().
   *
   * If a #permutation is set by
   * set_permutation(), it will
   * automatically be honored by
   * this function.
   *
   * The parameter @p adding does
   * not have any function, yet.
   */
  template <typename number2>
  void backward (Vector<number2> &,
                 const Vector<number2> &,
                 const bool transpose_diagonal,
                 const bool adding) const;
};


/**
 * Block SSOR preconditioning.
 *
 * The functions @p vmult and @p Tvmult execute a block-SSOR step,
 * based on the implementation in PreconditionBlockSOR.  This
 * class requires storage of the diagonal blocks and their inverses.
 *
 * See PreconditionBlock for requirements on the matrix. The blocks
 * used in this class must be contiguous and non-overlapping. An
 * overlapping Schwarz relaxation method can be found in
 * RelaxationBlockSSOR; that class does not offer preconditioning,
 * though.
 *
 * @note Instantiations for this template are provided for <tt>@<float@> and
 * @<double@></tt>; others can be generated in application programs (see the
 * section on @ref Instantiations in the manual).
 *
 * @author Ralf Hartmann, Guido Kanschat, 1999, 2000
 */
template<class MATRIX, typename inverse_type = typename MATRIX::value_type>
class PreconditionBlockSSOR : public virtual Subscriptor,
  private PreconditionBlockSOR<MATRIX, inverse_type>
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Define number type of matrix.
   */
  typedef typename MATRIX::value_type number;

  /**
   * Constructor.
   */
  PreconditionBlockSSOR ();

  // Keep AdditionalData accessible
  using typename PreconditionBlockSOR<MATRIX,inverse_type>::AdditionalData;

  // The following are the
  // functions of the base classes
  // which we want to keep
  // accessible.
  /**
   * Make initialization function
   * publicly available.
   */
  using PreconditionBlockSOR<MATRIX,inverse_type>::initialize;
  using PreconditionBlockSOR<MATRIX,inverse_type>::clear;
  using PreconditionBlockBase<inverse_type>::size;
  using PreconditionBlockBase<inverse_type>::inverse;
  using PreconditionBlockBase<inverse_type>::inverse_householder;
  using PreconditionBlockBase<inverse_type>::inverse_svd;
  using PreconditionBlockBase<inverse_type>::log_statistics;
  using PreconditionBlockSOR<MATRIX,inverse_type>::set_permutation;
  using PreconditionBlockSOR<MATRIX, inverse_type>::empty;
  using PreconditionBlockSOR<MATRIX, inverse_type>::el;
  using PreconditionBlockSOR<MATRIX,inverse_type>::set_same_diagonal;
  using PreconditionBlockSOR<MATRIX,inverse_type>::invert_diagblocks;

  /**
   * Execute block SSOR
   * preconditioning.
   *
   * This function will
   * automatically use the inverse
   * matrices if they exist, if not
   * then BlockSOR will waste much
   * time inverting the diagonal
   * block matrices in each
   * preconditioning step.
   */
  template <typename number2>
  void vmult (Vector<number2> &, const Vector<number2> &) const;

  /**
   * Same as vmult()
   */
  template <typename number2>
  void Tvmult (Vector<number2> &, const Vector<number2> &) const;

  /**
   * Perform one step of the SOR
   * iteration.
   */
  template <typename number2>
  void step (Vector<number2> &dst, const Vector<number2> &rhs) const;

  /**
   * Perform one step of the
   * transposed SOR iteration.
   */
  template <typename number2>
  void Tstep (Vector<number2> &dst, const Vector<number2> &rhs) const;
};

/*@}*/
//---------------------------------------------------------------------------

#ifndef DOXYGEN

template<class MATRIX, typename inverse_type>
inline bool
PreconditionBlock<MATRIX, inverse_type>::empty () const
{
  if (A == 0)
    return true;
  return A->empty();
}


template<class MATRIX, typename inverse_type>
inline unsigned int
PreconditionBlock<MATRIX, inverse_type>::n_blocks () const
{
  return this->size();
}


template<class MATRIX, typename inverse_type>
inline inverse_type
PreconditionBlock<MATRIX, inverse_type>::el (
  size_type i,
  size_type j) const
{
  const size_type bs = blocksize;
  const unsigned int nb = i/bs;

  const FullMatrix<inverse_type> &B = this->inverse(nb);

  const size_type ib = i % bs;
  const size_type jb = j % bs;

  if (jb + nb*bs != j)
    {
      return 0.;
    }

  return B(ib, jb);
}

//---------------------------------------------------------------------------

template<class MATRIX, typename inverse_type>
inline
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::Accessor::
Accessor (const PreconditionBlockJacobi<MATRIX, inverse_type> *matrix,
          const size_type row)
  :
  matrix(matrix),
  b_iterator(&matrix->inverse(0), 0, 0),
  b_end(&matrix->inverse(0), 0, 0)
{
  bs = matrix->block_size();
  a_block = row / bs;

  // This is the end accessor, which
  // does not hava a valid block.
  if (a_block == matrix->size())
    return;

  const size_type r = row % bs;

  b_iterator = matrix->inverse(a_block).begin(r);
  b_end = matrix->inverse(a_block).end();

  Assert (a_block < matrix->size(),
          ExcIndexRange(a_block, 0, matrix->size()));
}


template<class MATRIX, typename inverse_type>
inline
typename PreconditionBlockJacobi<MATRIX, inverse_type>::size_type
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::Accessor::row() const
{
  Assert (a_block < matrix->size(),
          ExcIteratorPastEnd());

  return bs * a_block + b_iterator->row();
}


template<class MATRIX, typename inverse_type>
inline
typename PreconditionBlockJacobi<MATRIX, inverse_type>::size_type
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::Accessor::column() const
{
  Assert (a_block < matrix->size(),
          ExcIteratorPastEnd());

  return bs * a_block + b_iterator->column();
}


template<class MATRIX, typename inverse_type>
inline
inverse_type
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::Accessor::value() const
{
  Assert (a_block < matrix->size(),
          ExcIteratorPastEnd());

  return b_iterator->value();
}


template<class MATRIX, typename inverse_type>
inline
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::
const_iterator(const PreconditionBlockJacobi<MATRIX, inverse_type> *matrix,
               const size_type row)
  :
  accessor(matrix, row)
{}


template<class MATRIX, typename inverse_type>
inline
typename PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator &
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::operator++ ()
{
  Assert (*this != accessor.matrix->end(), ExcIteratorPastEnd());

  ++accessor.b_iterator;
  if (accessor.b_iterator == accessor.b_end)
    {
      ++accessor.a_block;

      if (accessor.a_block < accessor.matrix->size())
        {
          accessor.b_iterator = accessor.matrix->inverse(accessor.a_block).begin();
          accessor.b_end = accessor.matrix->inverse(accessor.a_block).end();
        }
    }
  return *this;
}


template<class MATRIX, typename inverse_type>
inline
const typename PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::Accessor &
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::operator* () const
{
  return accessor;
}


template<class MATRIX, typename inverse_type>
inline
const typename PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::Accessor *
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::operator-> () const
{
  return &accessor;
}


template<class MATRIX, typename inverse_type>
inline
bool
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::
operator == (const const_iterator &other) const
{
  if (accessor.a_block == accessor.matrix->size() &&
      accessor.a_block == other.accessor.a_block)
    return true;

  if (accessor.a_block != other.accessor.a_block)
    return false;

  return (accessor.row() == other.accessor.row() &&
          accessor.column() == other.accessor.column());
}


template<class MATRIX, typename inverse_type>
inline
bool
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::
operator != (const const_iterator &other) const
{
  return ! (*this == other);
}


template<class MATRIX, typename inverse_type>
inline
bool
PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator::
operator < (const const_iterator &other) const
{
  return (accessor.row() < other.accessor.row() ||
          (accessor.row() == other.accessor.row() &&
           accessor.column() < other.accessor.column()));
}


template<class MATRIX, typename inverse_type>
inline
typename PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator
PreconditionBlockJacobi<MATRIX, inverse_type>::begin () const
{
  return const_iterator(this, 0);
}


template<class MATRIX, typename inverse_type>
inline
typename PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator
PreconditionBlockJacobi<MATRIX, inverse_type>::end () const
{
  return const_iterator(this, this->size() * this->block_size());
}


template<class MATRIX, typename inverse_type>
inline
typename PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator
PreconditionBlockJacobi<MATRIX, inverse_type>::begin (
  const size_type r) const
{
  Assert (r < this->A->m(), ExcIndexRange(r, 0, this->A->m()));
  return const_iterator(this, r);
}



template<class MATRIX, typename inverse_type>
inline
typename PreconditionBlockJacobi<MATRIX, inverse_type>::const_iterator
PreconditionBlockJacobi<MATRIX, inverse_type>::end (
  const size_type r) const
{
  Assert (r < this->A->m(), ExcIndexRange(r, 0, this->A->m()));
  return const_iterator(this, r+1);
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

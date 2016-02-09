// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2016 by the deal.II authors
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

#ifndef dealii__relaxation_block_h
#define dealii__relaxation_block_h

#include <deal.II/base/subscriptor.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/precondition_block_base.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <vector>
#include <set>

DEAL_II_NAMESPACE_OPEN

/**
 * Base class for the implementation of overlapping, multiplicative Schwarz
 * relaxation methods and smoothers.
 *
 * This class uses the infrastructure provided by PreconditionBlockBase. It
 * adds functions to initialize with a block list and to do the relaxation
 * step. The actual relaxation method with the interface expected by
 * SolverRelaxation and MGSmootherRelaxation is in the derived classes.
 *
 * This class allows for more general relaxation methods than
 * PreconditionBlock, since the index sets may be arbitrary and overlapping,
 * while there only contiguous, disjoint sets of equal size are allowed. As a
 * drawback, this class cannot be used as a preconditioner, since its
 * implementation relies on a straight forward implementation of the Gauss-
 * Seidel process.
 *
 * @ingroup Preconditioners
 * @author Guido Kanschat
 * @date 2010
 */
template <typename MatrixType, typename inverse_type=typename MatrixType::value_type>
class RelaxationBlock :
  protected PreconditionBlockBase<inverse_type>
{
private:
  /**
   * Define number type of matrix.
   */
  typedef typename MatrixType::value_type number;

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
   * Parameters for block relaxation methods. In addition to typical control
   * parameters like #relaxation, this object also contains the block
   * structure in #block_list and an optional ordering of the blocks in
   * #order.
   */
  class AdditionalData : public Subscriptor
  {
  public:
    /**
     * Constructor.
     */
    AdditionalData (const double relaxation = 1.,
                    const bool invert_diagonal = true,
                    const bool same_diagonal = false);

    /**
     * The mapping from indices to blocks. Each row of this pattern enumerates
     * the indices constituting a diagonal block to be inverted.
     */
    SparsityPattern block_list;

    /**
     * Relaxation parameter.
     */
    double relaxation;

    /**
     * Invert diagonal during initialization. Alternatively, diagonal blocks
     * are inverted on the fly, whenever they are used. While inverting blocks
     * in advance requires more memory, it usually saves a lot of computation.
     * See #same_diagonal on how you can avoid memory overhead.
     */
    bool invert_diagonal;

    /**
     * Assume all diagonal blocks are equal to save memory. If this flag is
     * true, then only the first diagonal block of the matrix is inverted and
     * stored. It is then used for all other blocks.
     *
     * \note Avoid setting this true if your blocks are not equal, in
     * particular if their sizes differ.
     */
    bool same_diagonal;
    /**
     * Choose the inversion method for the blocks.
     */
    typename PreconditionBlockBase<inverse_type>::Inversion inversion;

    /**
     * If #inversion is SVD, we can compute the Penrose-Moore inverse of the
     * blocks. In order to do so, we can specify here the threshold below
     * which a singular value will be considered zero and thus not inverted.
     * This parameter is used in the call to
     * LAPACKFullMatrix::compute_inverse_svd().
     */
    double threshold;

    /**
     * The order in which blocks should be traversed. This vector can initiate
     * several modes of execution:
     *
     * <ol>
     *
     * <li>If the length of the vector is zero, then the relaxation method
     * will be executed from first to last block.</li>
     *
     * <li> If the length is one, then the inner vector must have the same
     * size as the number of blocks. The relaxation method is applied in the
     * order given in this vector.</li>
     *
     * <li> If the outer vector has length greater one, then the relaxation
     * method is applied several times, each time in the order given by the
     * inner vector of the corresponding index. This mode can for instance be
     * used for ADI methods and similar direction sweeps.</li>
     *
     * </ol>
     */
    std::vector<std::vector<unsigned int> > order;
    /**
     * Return the memory allocated in this object.
     */
    std::size_t memory_consumption() const;
  };

  /**
   * Initialize matrix and additional information. In a second step, the
   * inverses of the diagonal blocks may be computed.
   *
   * Note that AdditionalData, different from other preconditioners, defines
   * quite large objects, and that therefore the object is not copied, but
   * rather a pointer is stored. Thus, the lifetime of
   * <code>additional_data</code> hast to exceed the lifetime of this object.
   */
  void initialize (const MatrixType     &A,
                   const AdditionalData &parameters);

  /**
   * Deletes the inverse diagonal block matrices if existent, sets the
   * blocksize to 0, hence leaves the class in the state that it had directly
   * after calling the constructor.
   */
  void clear();

  /**
   * Checks whether the object is empty.
   */
  bool empty () const;

  /**
   * Read-only access to entries.  This function is only possible if the
   * inverse diagonal blocks are stored.
   */
  value_type el(size_type i,
                size_type j) const;

  /**
   * Stores the inverse of the diagonal blocks in @p inverse. This costs some
   * additional memory - for DG methods about 1/3 (for double inverses) or 1/6
   * (for float inverses) of that used for the matrix - but it makes the
   * preconditioning much faster.
   *
   * It is not allowed to call this function twice (will produce an error)
   * before a call of <tt>clear(...)</tt> because at the second time there
   * already exist the inverse matrices.
   *
   * After this function is called, the lock on the matrix given through the
   * @p use_matrix function is released, i.e. you may overwrite of delete it.
   * You may want to do this in case you use this matrix to precondition
   * another matrix.
   */
  void invert_diagblocks();

protected:
  /**
   * Perform one block relaxation step.
   *
   * Depending on the arguments @p dst and @p pref, this performs an SOR step
   * (both reference the same vector) of a Jacobi step (both are different
   * vectors). For the Jacobi step, the calling function must copy @p dst to
   * @p pref after this.
   */
  template <typename number2>
  void do_step (
    Vector<number2>       &dst,
    const Vector<number2> &prev,
    const Vector<number2> &src,
    const bool backward) const;
  /**
   * Pointer to the matrix. Make sure that the matrix exists as long as this
   * class needs it, i.e. until calling @p invert_diagblocks, or (if the
   * inverse matrices should not be stored) until the last call of the
   * preconditioning @p vmult function of the derived classes.
   */
  SmartPointer<const MatrixType,RelaxationBlock<MatrixType,inverse_type> > A;

  /**
   * Control information.
   */
  SmartPointer<const AdditionalData, RelaxationBlock<MatrixType,inverse_type> > additional_data;
};


/**
 * Block Jacobi (additive Schwarz) method with possibly overlapping blocks.
 *
 * This class implements the step() and Tstep() functions expected by the
 * @ref ConceptRelaxationType "relaxation concept".
 * They perform an additive Schwarz method on the blocks provided in the block
 * list of AdditionalData. Differing from PreconditionBlockJacobi, these
 * blocks may be of varying size, non- contiguous, and overlapping. On the
 * other hand, this class does not implement the preconditioner interface
 * expected by Solver objects.
 *
 * @ingroup Preconditioners
 * @author Guido Kanschat
 * @date 2010
 */
template<typename MatrixType, typename inverse_type = typename MatrixType::value_type>
class RelaxationBlockJacobi : public virtual Subscriptor,
  protected RelaxationBlock<MatrixType, inverse_type>
{
public:
  /**
   * Default constructor.
   */
//    RelaxationBlockJacobi();

  /**
   * Define number type of matrix.
   */
  typedef typename MatrixType::value_type number;

  /**
   * Make type publicly available.
   */
  using typename RelaxationBlock<MatrixType,inverse_type>::AdditionalData;

  /**
   * Make initialization function publicly available.
   */
  using RelaxationBlock<MatrixType, inverse_type>::initialize;

  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MatrixType, inverse_type>::clear;

  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MatrixType, inverse_type>::empty;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MatrixType, inverse_type>::size;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MatrixType, inverse_type>::inverse;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MatrixType, inverse_type>::inverse_householder;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MatrixType, inverse_type>::inverse_svd;
  using PreconditionBlockBase<inverse_type>::log_statistics;
  /**
   * Perform one step of the Jacobi iteration.
   */
  template <typename number2>
  void step (Vector<number2> &dst, const Vector<number2> &rhs) const;

  /**
   * Perform one step of the Jacobi iteration.
   */
  template <typename number2>
  void Tstep (Vector<number2> &dst, const Vector<number2> &rhs) const;

  /**
   * Return the memory allocated in this object.
   */
  std::size_t memory_consumption() const;
};


/**
 * Block Gauss-Seidel method with possibly overlapping blocks.
 *
 * This class implements the step() and Tstep() functions expected by the
 * @ref ConceptRelaxationType "relaxation concept".
 * They perform a multiplicative Schwarz method on the blocks provided in the
 * block list of AdditionalData.  Differing from PreconditionBlockSOR, these
 * blocks may be of varying size, non-contiguous, and overlapping. On the
 * other hand, this class does not implement the preconditioner interface
 * expected by Solver objects.
 *
 * @ingroup Preconditioners
 * @author Guido Kanschat
 * @date 2010
 */
template<typename MatrixType, typename inverse_type = typename MatrixType::value_type>
class RelaxationBlockSOR : public virtual Subscriptor,
  protected RelaxationBlock<MatrixType, inverse_type>
{
public:
  /**
   * Default constructor.
   */
//    RelaxationBlockSOR();

  /**
   * Define number type of matrix.
   */
  typedef typename MatrixType::value_type number;

  /**
   * Make type publicly available.
   */
  using typename RelaxationBlock<MatrixType,inverse_type>::AdditionalData;

  /**
   * Make initialization function publicly available.
   */
  using RelaxationBlock<MatrixType, inverse_type>::initialize;

  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MatrixType, inverse_type>::clear;

  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MatrixType, inverse_type>::empty;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MatrixType, inverse_type>::size;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MatrixType, inverse_type>::inverse;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MatrixType, inverse_type>::inverse_householder;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MatrixType, inverse_type>::inverse_svd;
  using PreconditionBlockBase<inverse_type>::log_statistics;
  /**
   * Perform one step of the SOR iteration.
   */
  template <typename number2>
  void step (Vector<number2> &dst, const Vector<number2> &rhs) const;

  /**
   * Perform one step of the transposed SOR iteration.
   */
  template <typename number2>
  void Tstep (Vector<number2> &dst, const Vector<number2> &rhs) const;
};


/**
 * Symmetric block Gauss-Seidel method with possibly overlapping blocks.
 *
 * This class implements the step() and Tstep() functions expected by the
 * @ref ConceptRelaxationType "relaxation concept".
 * They perform a multiplicative Schwarz method on the blocks provided in the
 * block list of AdditionalData in symmetric fashion. Differing from
 * PreconditionBlockSSOR, these blocks may be of varying size, non-contiguous,
 * and overlapping. On the other hand, this class does not implement the
 * preconditioner interface expected by Solver objects.
 *
 * @ingroup Preconditioners
 * @author Guido Kanschat
 * @date 2010
 */
template<typename MatrixType, typename inverse_type = typename MatrixType::value_type>
class RelaxationBlockSSOR : public virtual Subscriptor,
  protected RelaxationBlock<MatrixType, inverse_type>
{
public:
  /**
   * Define number type of matrix.
   */
  typedef typename MatrixType::value_type number;

  /**
   * Make type publicly available.
   */
  using typename RelaxationBlock<MatrixType,inverse_type>::AdditionalData;

  /**
   * Make initialization function publicly available.
   */
  using RelaxationBlock<MatrixType, inverse_type>::initialize;

  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MatrixType, inverse_type>::clear;

  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MatrixType, inverse_type>::empty;

  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MatrixType, inverse_type>::size;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MatrixType, inverse_type>::inverse;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MatrixType, inverse_type>::inverse_householder;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MatrixType, inverse_type>::inverse_svd;
  using PreconditionBlockBase<inverse_type>::log_statistics;
  /**
   * Perform one step of the SOR iteration.
   */
  template <typename number2>
  void step (Vector<number2> &dst, const Vector<number2> &rhs) const;

  /**
   * Perform one step of the transposed SOR iteration.
   */
  template <typename number2>
  void Tstep (Vector<number2> &dst, const Vector<number2> &rhs) const;
};


DEAL_II_NAMESPACE_CLOSE

#endif

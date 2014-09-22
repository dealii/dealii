// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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

#ifndef __deal2__relaxation_block_h
#define __deal2__relaxation_block_h

#include <deal.II/base/subscriptor.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/precondition_block_base.h>
#include <deal.II/lac/block_list.h>

#include <vector>
#include <set>

DEAL_II_NAMESPACE_OPEN

/**
 * @warning The part of the interface based on BlockList may change in
 * a future release.
 *
 * Base class for the implementation of overlapping, multiplicative
 * Schwarz relaxation methods and smoothers.
 *
 * This class uses the infrastructure provided by
 * PreconditionBlockBase. It adds functions to initialize with a block
 * list and to do the relaxation step. The actual relaxation method
 * with the interface expected by SolverRelaxation and
 * MGSmootherRelaxation is in the derived classes.
 *
 * This class allows for more general relaxation methods than
 * PreconditionBlock, since the index sets may be arbitrary and
 * overlapping, while there only contiguous, disjoint sets of equal
 * size are allowed. As a drawback, this class cannot be used as a
 * preconditioner, since its implementation relies on a straight
 * forward implementation of the Gauss-Seidel process.
 *
 * @ingroup Preconditioners
 * @author Guido Kanschat
 * @date 2010
 */
template <class MATRIX, typename inverse_type=typename MATRIX::value_type>
class RelaxationBlock :
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
   * Parameters for block relaxation methods.
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
     * The mapping from indices to blocks.
     */
    SparsityPattern block_list;

    /**
     * Relaxation parameter.
     */
    double relaxation;

    /**
     * Invert diagonal during initialization.
     */
    bool invert_diagonal;

    /**
     * Assume all diagonal blocks are equal to save memory.
     */
    bool same_diagonal;
    /**
     * Choose the inversion method for the blocks.
     */
    typename PreconditionBlockBase<inverse_type>::Inversion inversion;

    /**
     * The if #inversion is SVD, the threshold below which a singular
     * value will be considered zero and thus not inverted. This
     * parameter is used in the call to
     * LAPACKFullMatrix::compute_inverse_svd().
     */
    double threshold;

    /**
     * The order in which blocks should be traversed. This vector can
     * initiate several modes of execution:
     *
     * <ol>
     *
     * <li>If the length of the vector is zero, then the relaxation
     * method will be executed from first to last block.</li>
     *
     * <li> If the length is one, then the inner vector must have the
     * same size as the number of blocks. The relaxation method is
     * applied in the order given in this vector.</li>
     *
     * <li> If the outer vector has length greater one, then the
     * relaxation method is applied several times, each time in the
     * order given by the inner vector of the corresponding
     * index. This mode can for instance be used for ADI methods and
     * similar direction sweeps.</li>
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
   * Initialize matrix and block size.  We store the matrix and the
   * block size in the preconditioner object. In a second step, the
   * inverses of the diagonal blocks may be computed.
   *
   * Additionally, a relaxation parameter for derived classes may be
   * provided.
   */
  void initialize (const MATRIX &A,
                   const AdditionalData &parameters);

  /**
   * Deletes the inverse diagonal block matrices if existent, sets the
   * blocksize to 0, hence leaves the class in the state that it had
   * directly after calling the constructor.
   */
  void clear();

  /**
   * Checks whether the object is empty.
   */
  bool empty () const;

  /**
   * Read-only access to entries.  This function is only possible if
   * the inverse diagonal blocks are stored.
   */
  value_type el(size_type i,
                size_type j) const;

  /**
   * Stores the inverse of the diagonal blocks in @p inverse. This
   * costs some additional memory - for DG methods about 1/3 (for
   * double inverses) or 1/6 (for float inverses) of that used for the
   * matrix - but it makes the preconditioning much faster.
   *
   * It is not allowed to call this function twice (will produce an
   * error) before a call of <tt>clear(...)</tt> because at the second
   * time there already exist the inverse matrices.
   *
   * After this function is called, the lock on the matrix given
   * through the @p use_matrix function is released, i.e. you may
   * overwrite of delete it.  You may want to do this in case you use
   * this matrix to precondition another matrix.
   */
  void invert_diagblocks();

protected:
  /**
   * Perform one block relaxation step.
   *
   * Depending on the arguments @p dst and @p pref, this performs an
   * SOR step (both reference the same vector) of a Jacobi step (both
   * are different vectors). For the Jacobi step, the calling function
   * must copy @p dst to @p pref after this.
   */
  template <typename number2>
  void do_step (
    Vector<number2>       &dst,
    const Vector<number2> &prev,
    const Vector<number2> &src,
    const bool backward) const;
  /**
   * Pointer to the matrix. Make sure that the matrix exists as long
   * as this class needs it, i.e. until calling @p invert_diagblocks,
   * or (if the inverse matrices should not be stored) until the last
   * call of the preconditoining @p vmult function of the derived
   * classes.
   */
  SmartPointer<const MATRIX,RelaxationBlock<MATRIX,inverse_type> > A;

  /**
   * Control information.
   */
  SmartPointer<const AdditionalData, RelaxationBlock<MATRIX,inverse_type> > additional_data;
};


/**
 * Block Jacobi (additive Schwarz) method with possibly overlapping
 * blocks.
 *
 * This class implements the step() and Tstep() functions expected by
 * SolverRelaxation and MGSmootherRelaxation. They perform an additive
 * Schwarz method on the blocks provided in the BlockList of
 * AdditionalData. Differing from PreconditionBlockJacobi, these
 * blocks may be of varying size, non-contiguous, and overlapping. On
 * the other hand, this class does not implement the preconditioner
 * interface expected by Solver objects.
 *
 * @ingroup Preconditioners
 * @author Guido Kanschat
 * @date 2010
 */
template<class MATRIX, typename inverse_type = typename MATRIX::value_type>
class RelaxationBlockJacobi : public virtual Subscriptor,
  protected RelaxationBlock<MATRIX, inverse_type>
{
public:
  /**
   * Default constructor.
   */
//    RelaxationBlockJacobi();

  /**
   * Define number type of matrix.
   */
  typedef typename MATRIX::value_type number;

  /**
   * Make type publicly available.
   */
  using typename RelaxationBlock<MATRIX,inverse_type>::AdditionalData;

  /**
   * Make initialization function
   * publicly available.
   */
  using RelaxationBlock<MATRIX, inverse_type>::initialize;

  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MATRIX, inverse_type>::clear;

  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MATRIX, inverse_type>::empty;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MATRIX, inverse_type>::size;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MATRIX, inverse_type>::inverse;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MATRIX, inverse_type>::inverse_householder;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MATRIX, inverse_type>::inverse_svd;
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
 * This class implements the step() and Tstep() functions expected by
 * SolverRelaxation and MGSmootherRelaxation. They perform a
 * multiplicative Schwarz method on the blocks provided in the
 * BlockList of AdditionalData. Differing from PreconditionBlockSOR,
 * these blocks may be of varying size, non-contiguous, and
 * overlapping. On the other hand, this class does not implement the
 * preconditioner interface expected by Solver objects.
 *
 * @ingroup Preconditioners
 * @author Guido Kanschat
 * @date 2010
 */
template<class MATRIX, typename inverse_type = typename MATRIX::value_type>
class RelaxationBlockSOR : public virtual Subscriptor,
  protected RelaxationBlock<MATRIX, inverse_type>
{
public:
  /**
   * Default constructor.
   */
//    RelaxationBlockSOR();

  /**
   * Define number type of matrix.
   */
  typedef typename MATRIX::value_type number;

  /**
   * Make type publicly available.
   */
  using typename RelaxationBlock<MATRIX,inverse_type>::AdditionalData;

  /**
   * Make initialization function
   * publicly available.
   */
  using RelaxationBlock<MATRIX, inverse_type>::initialize;

  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MATRIX, inverse_type>::clear;

  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MATRIX, inverse_type>::empty;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MATRIX, inverse_type>::size;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MATRIX, inverse_type>::inverse;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MATRIX, inverse_type>::inverse_householder;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MATRIX, inverse_type>::inverse_svd;
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
 * This class implements the step() and Tstep() functions expected by
 * SolverRelaxation and MGSmootherRelaxation. They perform a
 * multiplicative Schwarz method on the blocks provided in the
 * BlockList of AdditionalData in symmetric fashion. Differing from
 * PreconditionBlockSSOR, these blocks may be of varying size,
 * non-contiguous, and overlapping. On the other hand, this class does
 * not implement the preconditioner interface expected by Solver
 * objects.
 *
 * @ingroup Preconditioners
 * @author Guido Kanschat
 * @date 2010
 */
template<class MATRIX, typename inverse_type = typename MATRIX::value_type>
class RelaxationBlockSSOR : public virtual Subscriptor,
  protected RelaxationBlock<MATRIX, inverse_type>
{
public:
  /**
   * Define number type of matrix.
   */
  typedef typename MATRIX::value_type number;

  /**
   * Make type publicly available.
   */
  using typename RelaxationBlock<MATRIX,inverse_type>::AdditionalData;

  /**
   * Make initialization function publicly available.
   */
  using RelaxationBlock<MATRIX, inverse_type>::initialize;

  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MATRIX, inverse_type>::clear;

  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MATRIX, inverse_type>::empty;

  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MATRIX, inverse_type>::size;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MATRIX, inverse_type>::inverse;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MATRIX, inverse_type>::inverse_householder;
  /**
   * Make function of base class public again.
   */
  using RelaxationBlock<MATRIX, inverse_type>::inverse_svd;
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

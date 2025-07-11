// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_sparse_decomposition_h
#define dealii_sparse_decomposition_h

#include <deal.II/base/config.h>

#include <deal.II/lac/sparse_matrix.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup Preconditioners
 * @{
 */

/**
 * Abstract base class for incomplete decompositions of a sparse matrix into
 * sparse factors. This class can't be used by itself, but only as the base
 * class of derived classes that actually implement particular decompositions
 * such as SparseILU or SparseMIC.
 *
 * The decomposition is stored as a sparse matrix which is why this class is
 * derived from the SparseMatrix. Since it is not a matrix in the usual sense
 * (the stored entries are not those of a matrix, but of the two factors of
 * the original matrix), the derivation is <tt>protected</tt> rather than
 * <tt>public</tt>.
 *
 *
 * <h3>Fill-in</h3>
 *
 * Sparse decompositions are frequently used with additional fill-in, i.e.,
 * the sparsity structure of the decomposition is denser than that of the
 * matrix to be decomposed. The initialize() function of this class allows
 * this fill-in via the AdditionalData object as long as all entries present
 * in the original matrix are present in the decomposition also, i.e. the
 * sparsity pattern of the decomposition is a superset of the sparsity pattern
 * in the original matrix.
 *
 * Such fill-in can be accomplished by various ways, one of which is the
 * copy-constructor of the SparsityPattern class that allows the addition of
 * side-diagonals to a given sparsity structure.
 *
 *
 * <h3>Unified use of preconditioners</h3>
 *
 * While objects of this class can not be used directly (this class is only a
 * base class for others implementing actual decompositions), derived classes
 * such as SparseILU and SparseMIC can be used in the usual form as
 * preconditioners. For example, this works:
 * @code
 * SparseILU<double> ilu;
 * ilu.initialize(matrix, SparseILU<double>::AdditionalData(...));
 *
 * somesolver.solve (A, x, f, ilu);
 * @endcode
 *
 * Through the AdditionalData object it is possible to specify additional
 * parameters of the LU decomposition.
 *
 * 1/ The matrix diagonal can be strengthened by adding
 * <code>strengthen_diagonal</code> times the sum of the absolute row entries
 * of each row to the respective diagonal entries. By default no strengthening
 * is performed.
 *
 * 2/ By default, each initialize() function call creates its own sparsity.
 * For that, it copies the sparsity of <code>matrix</code> and adds a specific
 * number of extra off diagonal entries specified by
 * <code>extra_off_diagonals</code>.
 *
 * 3/ By setting <code>use_previous_sparsity=true</code> the sparsity is not
 * recreated but the sparsity of the previous initialize() call is reused
 * (recycled). This might be useful when several linear problems on the same
 * sparsity need to solved, as for example several Newton iteration steps on
 * the same triangulation. The default is <code>false</code>.
 *
 * 4/ It is possible to give a user defined sparsity to
 * <code>use_this_sparsity</code>. Then, no sparsity is created but
 * <code>*use_this_sparsity</code> is used to store the decomposed matrix. For
 * restrictions on the sparsity see section `Fill-in' above).
 *
 *
 * <h3>Particular implementations</h3>
 *
 * It is enough to override the initialize() and vmult() methods to implement
 * particular LU decompositions, like the true LU, or the Cholesky
 * decomposition. Additionally, if that decomposition needs fine tuned
 * diagonal strengthening on a per row basis, it may override the
 * get_strengthen_diagonal() method.
 */
template <typename number>
class SparseLUDecomposition : protected SparseMatrix<number>,
                              public virtual EnableObserverPointer
{
protected:
  /**
   * Constructor. Does nothing.
   *
   * Call the initialize() function before using this object as preconditioner
   * (vmult()).
   */
  SparseLUDecomposition();

public:
  /**
   * Declare type for container size.
   */
  using size_type = typename SparseMatrix<number>::size_type;

  /**
   * Destruction. Mark the destructor pure to ensure that this class isn't
   * used directly, but only its derived classes.
   */
  virtual ~SparseLUDecomposition() override = 0;

  /**
   * Deletes all member variables. Leaves the class in the state that it had
   * directly after calling the constructor
   */
  virtual void
  clear() override;

  /**
   * Parameters for SparseDecomposition.
   */
  class AdditionalData
  {
  public:
    /**
     * Constructor. For the parameters' description, see below.
     */
    explicit AdditionalData(const double       strengthen_diagonal   = 0.,
                            const unsigned int extra_off_diagonals   = 0,
                            const bool         use_previous_sparsity = false,
                            const SparsityPattern *use_this_sparsity = nullptr);

    /**
     * <code>strengthen_diag</code> times the sum of absolute row entries is
     * added to the diagonal entries.
     *
     * Per default, this value is zero, i.e. the diagonal is not strengthened.
     */
    double strengthen_diagonal;

    /**
     * By default, the <code>initialize(matrix, data)</code> function creates
     * its own sparsity. This sparsity has the same SparsityPattern as
     * <code>matrix</code> with some extra off diagonals the number of which
     * is specified by <code>extra_off_diagonals</code>.
     *
     * The user can give a SparsityPattern to <code>use_this_sparsity</code>.
     * Then this sparsity is used and the <code>extra_off_diagonals</code>
     * argument is ignored.
     */
    unsigned int extra_off_diagonals;

    /**
     * If this flag is true the initialize() function uses the same sparsity
     * that was used during the previous initialize() call.
     *
     * This might be useful when several linear problems on the same sparsity
     * need to solved, as for example several Newton iteration steps on the
     * same triangulation.
     */
    bool use_previous_sparsity;

    /**
     * When a SparsityPattern is given to this argument, the initialize()
     * function calls <code>reinit(*use_this_sparsity)</code> causing this
     * sparsity to be used.
     *
     * Note that the sparsity structures of <code>*use_this_sparsity</code>
     * and the matrix passed to the initialize function need not be equal.
     * Fill-in is allowed, as well as filtering out some elements in the
     * matrix.
     */
    const SparsityPattern *use_this_sparsity;
  };

  /**
   * This function needs to be called before an object of this class is used
   * as preconditioner.
   *
   * For more detail about possible parameters, see the class documentation
   * and the documentation of the SparseLUDecomposition::AdditionalData class.
   *
   * According to the <code>parameters</code>, this function creates a new
   * SparsityPattern or keeps the previous sparsity or takes the sparsity
   * given by the user to <code>data</code>. Then, this function performs the
   * LU decomposition.
   *
   * After this function is called the preconditioner is ready to be used
   * (using the <code>vmult</code> function of derived classes).
   */
  template <typename somenumber>
  void
  initialize(const SparseMatrix<somenumber> &matrix,
             const AdditionalData            parameters);

  /**
   * Return whether the object is empty. It calls the inherited
   * SparseMatrix::empty() function.
   */
  bool
  empty() const;

  /**
   * Return the dimension of the codomain (or range) space. It calls the
   * inherited SparseMatrix::m() function. Note that the matrix is of
   * dimension $m \times n$.
   */
  size_type
  m() const;

  /**
   * Return the dimension of the domain space. It calls the  inherited
   * SparseMatrix::n() function. Note that the matrix is of dimension $m
   * \times n$.
   */
  size_type
  n() const;

  /**
   * Adding Matrix-vector multiplication. Add <i>M*src</i> on <i>dst</i> with
   * <i>M</i> being this matrix.
   *
   * Source and destination must not be the same vector.
   */
  template <class OutVector, class InVector>
  void
  vmult_add(OutVector &dst, const InVector &src) const;

  /**
   * Adding Matrix-vector multiplication. Add <i>M<sup>T</sup>*src</i> to
   * <i>dst</i> with <i>M</i> being this matrix. This function does the same
   * as vmult_add() but takes the transposed matrix.
   *
   * Source and destination must not be the same vector.
   */
  template <class OutVector, class InVector>
  void
  Tvmult_add(OutVector &dst, const InVector &src) const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  virtual std::size_t
  memory_consumption() const;

  /**
   * @addtogroup Exceptions
   * @{
   */

  /**
   * Exception
   */
  DeclException1(ExcInvalidStrengthening,
                 double,
                 << "The strengthening parameter " << arg1
                 << " is not greater or equal than zero!");
  /** @} */
protected:
  /**
   * Copies the passed SparseMatrix onto this object. This object's sparsity
   * pattern remains unchanged.
   */
  template <typename somenumber>
  void
  copy_from(const SparseMatrix<somenumber> &matrix);

  /**
   * Performs the strengthening loop. For each row calculates the sum of
   * absolute values of its elements, determines the strengthening factor
   * (through get_strengthen_diagonal()) sf and multiplies the diagonal entry
   * with <code>sf+1</code>.
   */
  virtual void
  strengthen_diagonal_impl();

  /**
   * In the decomposition phase, computes a strengthening factor for the
   * diagonal entry in row <code>row</code> with sum of absolute values of its
   * elements <code>rowsum</code>.
   *
   * @note The default implementation in SparseLUDecomposition returns
   * <code>strengthen_diagonal</code>'s value. This variable is set to
   * a nonzero value in several of the derived classes.
   */
  virtual number
  get_strengthen_diagonal(const number rowsum, const size_type row) const;

  /**
   * The default strengthening value, returned by get_strengthen_diagonal().
   */
  double strengthen_diagonal;

  /**
   * For every row in the underlying SparsityPattern, this array contains a
   * pointer to the row's first afterdiagonal entry. Becomes available after
   * invocation of prebuild_lower_bound().
   */
  std::vector<const size_type *> prebuilt_lower_bound;

  /**
   * Fills the #prebuilt_lower_bound array.
   */
  void
  prebuild_lower_bound();

private:
  /**
   * In general this pointer is zero except for the case that no
   * SparsityPattern is given to this class. Then, a SparsityPattern is
   * created and is passed down to the SparseMatrix base class.
   *
   * Nevertheless, the SparseLUDecomposition needs to keep ownership of this
   * sparsity. It keeps this pointer to it enabling it to delete this sparsity
   * at destruction time.
   */
  SparsityPattern *own_sparsity;
};

/** @} */
//---------------------------------------------------------------------------

#ifndef DOXYGEN

template <typename number>
inline number
SparseLUDecomposition<number>::get_strengthen_diagonal(
  const number /*rowsum*/,
  const size_type /*row*/) const
{
  return strengthen_diagonal;
}



template <typename number>
inline bool
SparseLUDecomposition<number>::empty() const
{
  return SparseMatrix<number>::empty();
}


template <typename number>
inline typename SparseLUDecomposition<number>::size_type
SparseLUDecomposition<number>::m() const
{
  return SparseMatrix<number>::m();
}


template <typename number>
inline typename SparseLUDecomposition<number>::size_type
SparseLUDecomposition<number>::n() const
{
  return SparseMatrix<number>::n();
}

// Note: This function is required for full compatibility with
// the LinearOperator class. ::MatrixInterfaceWithVmultAdd
// picks up the vmult_add function in the protected SparseMatrix
// base class.
template <typename number>
template <class OutVector, class InVector>
inline void
SparseLUDecomposition<number>::vmult_add(OutVector      &dst,
                                         const InVector &src) const
{
  OutVector tmp;
  tmp.reinit(dst);
  this->vmult(tmp, src);
  dst += tmp;
}

// Note: This function is required for full compatibility with
// the LinearOperator class. ::MatrixInterfaceWithVmultAdd
// picks up the vmult_add function in the protected SparseMatrix
// base class.
template <typename number>
template <class OutVector, class InVector>
inline void
SparseLUDecomposition<number>::Tvmult_add(OutVector      &dst,
                                          const InVector &src) const
{
  OutVector tmp;
  tmp.reinit(dst);
  this->Tvmult(tmp, src);
  dst += tmp;
}

//---------------------------------------------------------------------------


template <typename number>
SparseLUDecomposition<number>::AdditionalData::AdditionalData(
  const double           strengthen_diag,
  const unsigned int     extra_off_diag,
  const bool             use_prev_sparsity,
  const SparsityPattern *use_this_spars)
  : strengthen_diagonal(strengthen_diag)
  , extra_off_diagonals(extra_off_diag)
  , use_previous_sparsity(use_prev_sparsity)
  , use_this_sparsity(use_this_spars)
{}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_sparse_decomposition_h

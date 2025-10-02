// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_matrix_scaling_h
#define dealii_matrix_scaling_h

#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#ifdef DEAL_II_WITH_MPI
#  include <deal.II/base/partitioner.h>
#endif

#ifdef DEAL_II_WITH_TRILINOS
#  include <deal.II/lac/trilinos_sparse_matrix.h>
#  include <deal.II/lac/trilinos_vector.h>
#endif

#ifdef DEAL_II_WITH_PETSC
#  include <deal.II/lac/petsc_sparse_matrix.h>
#  include <deal.II/lac/petsc_vector.h>
#endif

DEAL_II_NAMESPACE_OPEN

/**
 * This class provides access to various matrix scaling algorithms. The
 * scaling algorithms scale the matrix in place by the action of two diagonal
 * matrices: one acting on the rows and one acting on the columns. The diagonals
 * of the two scaling matrices applied on the last scaled matrix are stored in
 * the class as vectors. The algorithms implemented here are the Sinkhorn-Knopp
 * algorithm described in this <a
 * href="https://doi.org/10.1137/060659624">article</a> and the scaling
 * algorithm that preserves symmetry described in this <a
 * href="https://doi.org/10.1137/110825753">article</a>.
 *
 * The parallel implementation with MPI of the algorithms follows the
 * description in this <a
 * href="https://doi.org/10.1007/978-3-540-92859-1_27">article</a>. The
 * diagonals of the scaling matrices are distributed between the MPI processes
 * following the row distribution of the matrix. The diagonal of the column
 * scaling is distributed in a balanced way if the scaled matrix is not square.
 *
 *  <h4>Instantiations</h4>
 *
 * There are instantiations of this class for SparseMatrix<double>,
 * SparseMatrix<float>, FullMatrix<double>, FullMatrix<float>,
 * BlockSparseMatrix<float>, BlockSparseMatrix<double>, SparseMatrixEZ<double>
 * and SparseMatrixEZ<float>.
 *
 * The distributed matrices supported are TrilinosWrappers::SparseMatrix and
 * PETScWrappers::MPI::SparseMatrix.
 */
class MatrixScaling : public EnableObserverPointer
{
public:
  /**
   * AdditionalData allows to choose the scaling algorithm and contains
   * parameters for the scaling algorithms.
   */
  struct AdditionalData
  {
    /**
     * Supported scaling algorithms within <tt>MatrixScaling</tt>.
     */
    enum class ScalingAlgorithm
    {
      sinkhorn_knopp,
      l1_linf_symmetry_preserving
    };

    /**
     * Struct that contains parameters for the Sinkhorn-Knopp algorithm.
     */
    struct SKParameters
    {
      /**
       * Scaling norm available to use in the Sinkhorn-Knopp algorithm.
       */
      enum class NormType
      {
        l1,
        l_infinity
      };

      /**
       * Construct a new SKParameters object with the given parameters.
       */
      explicit SKParameters(const NormType     norm_type      = NormType::l1,
                            const unsigned int max_iterations = 20)
        : max_iterations(max_iterations)
        , norm_type(norm_type)
      {}

      /**
       * Maximum number of iterations to perform in the SK algorithm.
       */
      unsigned int max_iterations;

      /**
       * Norm type to use in the SK algorithm.
       */
      NormType norm_type;
    };

    /**
     * Struct that contains parameters for the symmetry preserving scaling
     * algorithm. Recommended usage is to alternate infinity norm steps with l1
     * norm steps.
     */
    struct l1linfParameters
    {
      /**
       * Construct a new l1linfParameters object with the given parameters.
       */
      explicit l1linfParameters(const unsigned int start_inf_norm_steps = 1,
                                const unsigned int l1_norm_steps        = 3,
                                const unsigned int end_inf_norm_steps   = 1)
        : start_inf_norm_steps(start_inf_norm_steps)
        , l1_norm_steps(l1_norm_steps)
        , end_inf_norm_steps(end_inf_norm_steps)
      {}

      /**
       * Number of initial linfty norm scaling steps.
       */
      unsigned int start_inf_norm_steps;

      /**
       * Number of intermediate l1 norm scaling steps.
       */
      unsigned int l1_norm_steps;

      /**
       * Number of final linfty norm scaling steps.
       */
      unsigned int end_inf_norm_steps;
    };

    /**
     * Construct a new AdditionalData object with the given parameters.
     */
    explicit AdditionalData(
      const double           scaling_tolerance = 1e-5,
      const ScalingAlgorithm alg =
        ScalingAlgorithm::l1_linf_symmetry_preserving,
      const SKParameters     sk_params     = SKParameters(),
      const l1linfParameters l1linf_params = l1linfParameters())
      : scaling_tolerance(scaling_tolerance)
      , algorithm(alg)
      , sinkhorn_knopp_parameters(sk_params)
      , l1linf_parameters(l1linf_params)
    {}

    /**
     * Tolerance for the scaling algorithms. Convergence is achieved if the
     * difference in norm between a vector of all ones and the row and column
     * norms vectors is below this tolerance.
     */
    double scaling_tolerance;

    /**
     * Scaling algorithm to use.
     */
    ScalingAlgorithm algorithm;

    /**
     * Parameters for the Sinkhorn-Knopp algorithm.
     */
    SKParameters sinkhorn_knopp_parameters;

    /**
     * Parameters for the symmetry preserving scaling algorithm.
     */
    l1linfParameters l1linf_parameters;
  };

  /**
   * Constructor. The constructor takes a <tt>AdditionalData</tt> object that
   * contains the parameters for the scaling algorithms and the algorithm
   * selected.
   */
  explicit MatrixScaling(const AdditionalData &control = AdditionalData());

  /**
   * Destructor.
   */
  ~MatrixScaling() = default;

  /**
   * Scale the input <tt>matrix</tt> in place according to the selected
   * algorithm in <tt>AdditionalData</tt>.
   */
  template <class Matrix>
  void
  scale_matrix(Matrix &matrix);

  /**
   * Scale the system <tt>matrix</tt> in place according to the selected
   * algorithm in <tt>AdditionalData</tt> and scale in place the right hand side
   * vector <tt>rhs</tt> according to the row scaling. The solution of the
   * scaled system can be scaled back to the original system by calling
   * <tt>scale_system_solution()</tt>.
   */
  template <class Matrix, class VectorType>
  void
  scale_linear_system(Matrix &matrix, VectorType &rhs);

  /**
   * Scale back the linear system solution vector according to the column
   * scaling. This function should be called after solving the scaled system
   * obtained by calling <tt>scale_linear_system()</tt>.
   */
  template <class VectorType>
  void
  scale_system_solution(VectorType &sol) const;

  /**
   * Return a const reference to the (locally owned in distributed setting) row
   * scaling (i.e. the diagonal of the row scaling matrix).
   */
  const Vector<double> &
  get_row_scaling() const;

  /**
   * Return a const reference to the (locally owned in distributed setting)
   * column scaling (i.e. the diagonal of the column scaling matrix).
   */
  const Vector<double> &
  get_column_scaling() const;

  /**
   * Return whether the scaling algorithm has converged or not. Convergence is
   * not necessary to have a good scaling but it can be helpful to know if it
   * was achieved.
   */
  bool
  is_converged() const;

private:
  /**
   * Struct that contains the parameters for the scaling algorithms.
   */
  AdditionalData control;

  /**
   * Vector that contains the (locally owned in distributed setting) row scaling
   * (i.e. the diagonal of the row scaling matrix).
   */
  Vector<double> row_scaling;

  /**
   * Vector that contains the (locally owned in distributed setting) column
   * scaling (i.e. the diagonal of the column scaling matrix).
   */
  Vector<double> column_scaling;

  /**
   * Flag that indicates whether the last scaling algorithm run converged or
   * not.
   */
  bool converged;

  /**
   * IndexSet storing the locally owned rows of the row scaling.
   */
  IndexSet locally_owned_rows;

  /**
   * IndexSet storing the locally owned column indexes of the column scaling.
   */
  IndexSet locally_owned_cols;

  /**
   * IndexSet storing the ghost column indexes of the column scaling. Here are
   * stored the columns indices that will require the scaling from other MPI
   * ranks.
   */
  IndexSet ghost_columns;

  /**
   * Vector storing the owner MPI rank of each ghost column index.
   */
  std::vector<types::global_dof_index> ghost_column_owners;

#ifdef DEAL_II_WITH_MPI
  /**
   * Partitioner for communications.
   */
  Utilities::MPI::Partitioner partitioner;
#endif

  /**
   * Implementation of the l1 scaling iterations in the symmetry preserving
   * algorithm.
   */
  template <class Matrix>
  void
  l1_scaling(Matrix &matrix, const unsigned int nsteps);

  /**
   * Implementation of the linfty scaling iterations in the symmetry preserving
   * algorithm.
   */
  template <class Matrix>
  void
  linfty_scaling(Matrix &matrix, const unsigned int nsteps);

  /**
   * Implementation of the Sinkhorn-Knopp scaling iterations.
   */
  template <class Matrix>
  void
  sk_scaling(Matrix &matrix, const unsigned int nsteps);
};

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_matrix_scaling_h

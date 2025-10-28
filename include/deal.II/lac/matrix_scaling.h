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
#include <deal.II/base/partitioner.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

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
 * matrices: one acting on the rows and one acting on the columns. The
 * algorithms implemented here are the Sinkhorn-Knopp algorithm described in
 * @cite Knight2008 and the scaling algorithm that preserves symmetry described
 * in @cite Knight2014.
 *
 * The parallel implementation with MPI of the algorithms follows the
 * description in @cite Amestoy2008. The diagonals of the scaling matrices are
 * distributed between the MPI processes following the row distribution of the
 * matrix. The diagonal of the column scaling is distributed in a balanced way
 * if the scaled matrix is not square.
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
   * Vectorial norms used to check convergence in the scaling algorithms.
   */
  enum class ConvergenceNormType
  {
    //! l_1 vector norm
    l1,

    //! l_infinity vector norm
    l_infty
  };

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
      //! Sinkhorn-Knopp scaling algorithm
      /**
       * The Sinkhorn-Knopp algorithm (also known as RAS algorithm)
       * alternatingly scales rows and columns by the chosen <tt>NormType</tt>
       * and attempts to make the matrix doubly stochastic (if all the entries
       * are non-negative). The algorithm does not preserve symmetry.
       */
      sinkhorn_knopp,

      //! Symmetry preserving scaling algorithm
      /**
       * The symmetry preserving algorithm simultaneously scales rows and
       * columns preserving symmetry. The symmetry preserving algorithm
       * alternates l_infinity norm scaling with l_1 norm scaling as it has been
       * shown to yield the best results scaling linear systems. The number of
       * scaling steps in each norm is determined by the
       * <tt>l1linfParameters</tt> struct
       */
      l1_linf_symmetry_preserving
    };

    /**
     * Struct that contains parameters for the Sinkhorn-Knopp algorithm.
     */
    struct SKParameters
    {
      /**
       * Scaling norm available to use in the Sinkhorn-Knopp algorithm. By
       * setting this parameter the SK algorithm alternatingly scales rows and
       * column using the selected vector norm type.
       */
      enum class NormType
      {
        //! l_1 vector norm
        l1,

        //! l_infinity vector norm
        l_infty
      };

      /**
       * Construct a new SKParameters object with the given parameters.
       */
      explicit SKParameters(const NormType     norm_type      = NormType::l1,
                            const unsigned int max_iterations = 20);

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
                                const unsigned int end_inf_norm_steps   = 1);

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
      const l1linfParameters l1linf_params = l1linfParameters());


    /**
     * Tolerance for the scaling algorithms. Convergence is achieved if the
     * difference in norm between a vector of all ones and the row and column
     * norms vectors is below this tolerance.
     * Convergence is not needed to obtain a good scaling and reaching it might
     * require a lot of iterations and thus time. Convergence is not guaranteed
     * for all kind of matrices. For more details check out the articles cited
     * in the class description.
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
   *
   * The function returns whether the scaling algorithm has converged or not.
   * Convergence is not necessary to have a good scaling but it can be helpful
   * to know if it was achieved.
   */
  template <class Matrix>
  bool
  find_scaling_and_scale_matrix(Matrix &matrix);

  /**
   * Scale the system <tt>matrix</tt> in place according to the selected
   * algorithm in <tt>AdditionalData</tt> and scale in place the right hand side
   * vector <tt>rhs</tt> according to the row scaling. The solution of the
   * scaled system can be scaled back to the original system by calling
   * <tt>scale_system_solution()</tt>.
   *
   * The function returns whether the scaling algorithm has converged or not.
   * Convergence is not necessary to have a good scaling but it can be helpful
   * to know if it was achieved.
   */
  template <class Matrix, class VectorType>
  bool
  find_scaling_and_scale_linear_system(Matrix &matrix, VectorType &rhs);

  /**
   * Scale back the linear system solution vector according to the column
   * scaling. This function should be called after solving the scaled system
   * obtained by calling <tt>find_scaling_and_scale_linear_system()</tt>.
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
  std::vector<unsigned int> ghost_column_owners;

  /**
   * Partitioner for communications.
   */
  Utilities::MPI::Partitioner partitioner;

  /**
   * Implementation of the l1 scaling iterations in the symmetry preserving
   * algorithm. The function returns whether the scaling algorithm has converged
   * or not.
   */
  template <class Matrix>
  bool
  do_l1_scaling(Matrix &matrix, const unsigned int nsteps);

  /**
   * Implementation of the linfty scaling iterations in the symmetry preserving
   * algorithm. The function returns whether the scaling algorithm has converged
   * or not.
   */
  template <class Matrix>
  bool
  do_linfty_scaling(Matrix &matrix, const unsigned int nsteps);

  /**
   * Implementation of the Sinkhorn-Knopp scaling iterations. The function
   * returns whether the scaling algorithm has converged or not.
   */
  template <class Matrix>
  bool
  do_sk_scaling(Matrix &matrix, const unsigned int nsteps);

  /*
   * Check convergence of the sequential scaling algorithms. Only on either row
   * or columns for SK scaling.
   */
  template <typename Number>
  bool
  check_convergence(const Vector<Number>      &row_col_norm,
                    const ConvergenceNormType &norm_type) const;

  /*
   * Check convergence of the sequential scaling algorithms. On both row
   * and columns for symmetry preserving scaling.
   */
  template <typename Number>
  bool
  check_convergence(const Vector<Number>      &row_norm,
                    const Vector<Number>      &col_norm,
                    const ConvergenceNormType &norm_type) const;

  /*
   * Check convergence of the parallel scaling algorithm. Only on either row
   * or columns for SK scaling.
   */
  bool
  check_convergence(const Vector<double>      &local_row_col_norm,
                    const ConvergenceNormType &norm_type,
                    const MPI_Comm             mpi_communicator) const;

  /*
   * Check convergence of the parallel scaling algorithm. Only on both row
   * and columns for symmetry preserving scaling.
   */
  bool
  check_convergence(const Vector<double>      &local_row_norm,
                    const Vector<double>      &local_col_norm,
                    const ConvergenceNormType &norm_type,
                    const MPI_Comm             mpi_communicator) const;

  /*
   * Fill the send_data map with the partial column norms that need to be sent
   * to other MPI processes. Prepare the local_col_norms vector. This is a
   * helper function for the distributed implementation of the scaling
   * algorithms.
   */
  void
  send_prepare_col_norms(
    const std::map<types::global_dof_index, double> &partial_column_norms,
    std::map<unsigned int,
             std::vector<std::pair<types::global_dof_index, double>>>
                   &send_data,
    Vector<double> &local_col_norms);

  /*
   * Fill the send_column_norms map with the updated column norms that need to
   * be sent to other MPI processes. This is a helper function for the
   * distributed implementation of the scaling algorithms.
   */
  void
  send_prepare_updated_col_norms(
    const std::map<unsigned int,
                   std::vector<std::pair<types::global_dof_index, double>>>
                         &received_data,
    const Vector<double> &local_col_norms,
    std::map<unsigned int,
             std::vector<std::pair<types::global_dof_index, double>>>
      &send_column_norms);
};

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_matrix_scaling_h

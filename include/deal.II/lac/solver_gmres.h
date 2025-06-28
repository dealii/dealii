// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_solver_gmres_h
#define dealii_solver_gmres_h



#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/lac/block_vector_base.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/orthogonalization.h>
#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector.h>

#include <boost/signals2/signal.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// forward declarations
#ifndef DOXYGEN
namespace LinearAlgebra
{
  namespace distributed
  {
    template <typename, typename>
    class Vector;
  } // namespace distributed
} // namespace LinearAlgebra
#endif

/**
 * @addtogroup Solvers
 * @{
 */

namespace internal
{
  /**
   * A namespace for helper classes and functions of the GMRES solver.
   */
  namespace SolverGMRESImplementation
  {
    /**
     * Class to hold temporary vectors.  This class automatically allocates a
     * new vector, once it is needed.
     *
     * A future version should also be able to shift through vectors
     * automatically, avoiding restart.
     */
    template <typename VectorType>
    class TmpVectors
    {
    public:
      /**
       * Constructor. Prepares an array of @p VectorType of length @p
       * max_size.
       */
      TmpVectors(const unsigned int max_size, VectorMemory<VectorType> &vmem);

      /**
       * Destructor. Delete all allocated vectors.
       */
      ~TmpVectors() = default;

      /**
       * Get vector number @p i. If this vector was unused before, an error
       * occurs.
       */
      VectorType &
      operator[](const unsigned int i) const;

      /**
       * Get vector number @p i. Allocate it if necessary.
       *
       * If a vector must be allocated, @p temp is used to reinit it to the
       * proper dimensions.
       */
      VectorType &
      operator()(const unsigned int i, const VectorType &temp);

      /**
       * Return size of data vector. It is used in the solver to store
       * the Arnoldi vectors.
       */
      unsigned int
      size() const;


    private:
      /**
       * Pool where vectors are obtained from.
       */
      VectorMemory<VectorType> &mem;

      /**
       * Field for storing the vectors.
       */
      std::vector<typename VectorMemory<VectorType>::Pointer> data;
    };



    /**
     * Class that performs the Arnoldi orthogonalization process within the
     * SolverGMRES, SolverFGMRES, and SolverMPGMRES classes. It uses one of
     * the algorithms in LinearAlgebra::OrthogonalizationStrategy for the
     * work on the global vectors, transforms the resulting Hessenberg
     * matrix into an upper triangular matrix by Givens rotations, and
     * eventually solves the minimization problem in the projected Krylov
     * space.
     */
    template <typename Number>
    class ArnoldiProcess
    {
    public:
      /**
       * Initialize the data structures in this class with the given
       * parameters for the solution process.
       */
      void
      initialize(const LinearAlgebra::OrthogonalizationStrategy
                                    orthogonalization_strategy,
                 const unsigned int max_basis_size,
                 const bool         force_reorthogonalization);

      /**
       * Orthonormalize the vector at the position @p n within the array
       * @p orthogonal_vectors against the @p n orthonormal vectors with
       * indices <tt>0, ..., n - 1</tt> using the modified or classical
       * Gram-Schmidt algorithm. The class internally stores the factors used
       * for orthogonalization in an upper Hessenberg matrix. For the
       * classical Gram-Schmidt and modified Gram-Schmidt algorithms, loss of
       * orthogonality is checked every fifth step (in case it is not yet
       * already set via the initialize() function). In case this is detected,
       * all subsequent iterations use re-orthogonalization as stored
       * internally in this class, and a call to the optional signal is made.
       *
       * Note that the projected Hessenberg matrix and its factorization are
       * only consistent if @p n is incremented by one for each successive
       * call, or if @p n is zero when starting to build a new orthogonal
       * basis in restarted GMRES; an assertion will be raised if this
       * assumption is not fulfilled.
       *
       * Within this function, the factors for the QR factorization are
       * computed alongside the Hessenberg matrix, and an estimate of the
       * residual in the subspace is returned from this function.
       */
      template <typename VectorType>
      double
      orthonormalize_nth_vector(
        const unsigned int                        n,
        TmpVectors<VectorType>                   &orthogonal_vectors,
        const unsigned int                        accumulated_iterations = 0,
        const boost::signals2::signal<void(int)> &reorthogonalize_signal =
          boost::signals2::signal<void(int)>());

      /**
       * Using the matrix and right hand side computed during the
       * factorization, solve the underlying minimization problem for the
       * residual in the Krylov space, returning the resulting solution as a
       * const reference. Note that the dimension of the vector is set to the
       * size of the Krylov space.
       */
      const Vector<double> &
      solve_projected_system(const bool orthogonalization_finished);

      /**
       * Return the upper Hessenberg matrix resulting from the
       * Gram-Schmidt orthogonalization process.
       */
      const FullMatrix<double> &
      get_hessenberg_matrix() const;

      /**
       * Temporary vector to implement work for deal.II vector types
       */
      std::vector<const Number *> vector_ptrs;

    private:
      /**
       * Projected system matrix in upper Hessenberg form.
       */
      FullMatrix<double> hessenberg_matrix;

      /**
       * Upper triangular matrix that results from performing the QR
       * factorization with Givens rotations on the upper Hessenberg matrix; the
       * matrix Q is contained in the array givens_rotations.
       */
      FullMatrix<double> triangular_matrix;

      /**
       * Representation of the factor Q in the QR factorization of the
       * Hessenberg matrix.
       */
      std::vector<std::pair<double, double>> givens_rotations;

      /**
       * Right-hand side vector for orthogonalization.
       */
      Vector<double> projected_rhs;

      /**
       * Solution vector when computing the minimization in the projected
       * Krylov space.
       */
      Vector<double> projected_solution;

      /**
       * Auxiliary vector for orthogonalization.
       */
      Vector<double> h;

      /**
       * Flag to keep track reorthogonalization, which is checked every fifth
       * iteration by default for
       * LinearAlgebra::OrthogonalizationStrategy::classical_gram_schmidt and
       * LinearAlgebra::OrthogonalizationStrategy::modified_gram_schmidt; for
       * LinearAlgebra::OrthogonalizationStrategy::delayed_classical_gram_schmidt,
       * no check is made.
       */
      bool do_reorthogonalization;

      /**
       * Selected orthogonalization algorithm.
       */
      LinearAlgebra::OrthogonalizationStrategy orthogonalization_strategy;

      /**
       * This is a helper function to perform the incremental computation of
       * the QR factorization of the Hessenberg matrix involved in the Arnoldi
       * process. More precisely, it transforms the member variable
       * @p hessenberg_matrix into an upper triangular matrix R labeled
       * @p matrix, an orthogonal matrix Q represented by a vector of Givens
       * rotations, and the associated right hand side to minimize the norm of
       * the solution in the Krylov subspace.
       *
       * More precisely, this function is called once a new column is added to
       * the Hessenberg matrix and performs all necessary steps for that
       * column. First, all evaluations with the Givens rotations resulting
       * from the previous elimination steps are performed. Then, the single
       * additional entry below the diagonal in the Hessenberg matrix is
       * eliminated by a Givens rotation, appending a new pair of Givens
       * factors, and the right-hand side vector in the projected system is
       * updated. The column number @p col for which the Gram-Schmidt should
       * run needs to be given, because the algorithmic variant with delayed
       * orthogonalization might lag by one step compared to the other sizes
       * in the problem, and needs to perform additional computations.
       *
       * In most cases, the matrices and vectors passed to this function are
       * the member variables of the present class, but also other scenarios
       * are supported. The function returns the modulus of the last entry in
       * the transformed right-hand side, which is the obtained residual of
       * the global vector x after minimization within the Krylov space.
       */
      double
      do_givens_rotation(const bool          delayed_reorthogonalization,
                         const int           col,
                         FullMatrix<double> &matrix,
                         std::vector<std::pair<double, double>> &rotations,
                         Vector<double>                         &rhs);
    };
  } // namespace SolverGMRESImplementation
} // namespace internal



/**
 * Implementation of the Restarted Preconditioned Direct Generalized Minimal
 * Residual Method (GMRES). The stopping criterion is the norm of the residual.
 *
 * The AdditionalData structure allows to control the size of the Arnoldi
 * basis used for orthogonalization (default: 30 vectors). It is related to
 * the number of temporary vectors used, which is the basis size plus
 * two. Additionally, it allows you to choose between right or left
 * preconditioning (default: left preconditioning). Furthermore, it
 * includes a flag indicating whether or not the default residual is used as
 * stopping criterion and an option for the orthogonalization algorithm, see
 * LinearAlgebra::OrthogonalizationStrategy for available options.
 *
 *
 * <h3>Left versus right preconditioning</h3>
 *
 * @p AdditionalData allows you to choose between left and right
 * preconditioning. Left preconditioning, conceptually, corresponds to
 * replacing the linear system $Ax=b$ by $P^{-1}Ax=P^{-1}b$ where
 * $P^{-1}$ is the preconditioner (i.e., an approximation of
 * $A^{-1}$). In contrast, right preconditioning should be understood
 * as replacing $Ax=b$ by $AP^{-1}y=b$, solving for $y$, and then
 * computing the solution of the original problem as $x=P^{-1}y$. Note
 * that in either case, $P^{-1}$ is simply an operator that can be
 * applied to a vector; that is, it is not the inverse of some
 * operator that also separately has to be available. In practice,
 * $P^{-1}$ should be an operator that approximates multiplying by
 * $A^{-1}$.
 *
 * The choice between left and right preconditioning also affects
 * which kind of residual is used to measure convergence. With left
 * preconditioning, this is the <b>preconditioned</b> residual
 * $r_k=P^{-1}b-P^{-1}Ax_k$ given the approximate solution $x_k$ in
 * the $k$th iteration, while with right preconditioning, it is the
 * residual $r_k=b-Ax_k$ of the unpreconditioned system.
 *
 * Optionally, this behavior can be overridden by using the flag
 * AdditionalData::use_default_residual. A <tt>true</tt> value refers to the
 * behavior described in the previous paragraph, while <tt>false</tt> reverts
 * it. Be aware though that additional residuals have to be computed in this
 * case, impeding the overall performance of the solver.
 *
 *
 * <h3> Preconditioners need to be linear operators </h3>
 *
 * GMRES expects the preconditioner to be a *linear* operator, i.e.,
 * the operator $P^{-1}$ used as preconditioner needs to satisfy
 * $P^{-1}(x+y) = P^{-1}x + P^{-1}y$ and $P^{-1}(\alpha x) = \alpha
 * P^{-1}x$. For many preconditioners, this is true. For example, if
 * you used Jacobi preconditioning, then $P^{-1}$ is a diagonal matrix
 * whose diagonal entries equal $\frac{1}{A_{ii}}$. In this case, the
 * operator $P^{-1}$ is clearly linear since it is simply the
 * multiplication of a given vector by a fixed matrix.
 *
 * On the other hand, if $P^{-1}$ involves more complicated
 * operations, it is sometimes *not* linear. The typical case to
 * illustrate this is where $A$ is a block matrix and $P^{-1}$
 * involves multiplication with blocks (as done, for example, in
 * step-20, step-22, and several other preconditioners) where one
 * block involves a linear solve. For example, in a Stokes problem,
 * the preconditioner may involve a linear solve with the upper left
 * $A_{uu}$ block. If this linear solve is done exactly (e.g., via a
 * direct solver, or an iterative solver with a very tight tolerance),
 * then the linear solve corresponds to multiplying by $A^{-1}_{uu}$,
 * which is a linear operation. On the other hand, if one uses an
 * iterative solver with a loose tolerance (e.g.,
 * `1e-3*right_hand_side.l2_norm()`), then many solvers like CG will
 * find the solution in a Krylov subspace of fairly low dimension;
 * crucially, this subspace is built iteratively starting with the
 * initial residual -- in other words, the *subspace depends on the
 * right hand side*, and consequently the solution returned by such a
 * solver *is not a linear operation on the given right hand side* of
 * the linear system being solved.
 *
 * In cases such as these, the preconditioner with its inner, inexact
 * linear solve is not a linear operator. This violates the
 * assumptions of GMRES, and often leads to unnecessarily many GMRES
 * iterations. The solution is to use the SolverFGMRES class instead,
 * which does not rely on the assumption that the preconditioner is a
 * linear operator, and instead explicitly does the extra work
 * necessary to satisfy the assumptions that lead GMRES to implicitly
 * require a linear operator as preconditioner.
 *
 *
 * <h3>The size of the Arnoldi basis</h3>
 *
 * The maximal basis size is controlled by AdditionalData::max_basis_size. If
 * the number of iteration steps exceeds this number, all basis vectors are
 * discarded and the iteration starts anew from the approximation obtained so
 * far. This algorithm strategy is typically called restarted GMRES method.
 *
 * Note that the minimizing property of GMRES only pertains to the Krylov
 * space spanned by the Arnoldi basis. Therefore, restarted GMRES is
 * <b>not</b> minimizing anymore. The choice of the basis length is a
 * trade-off between memory consumption and convergence speed, since a longer
 * basis means minimization over a larger space.
 *
 * For the requirements on matrices and vectors in order to work with this
 * class, see the documentation of the Solver base class.
 *
 *
 * <h3>Observing the progress of linear solver iterations</h3>
 *
 * The solve() function of this class uses the mechanism described in the
 * Solver base class to determine convergence. This mechanism can also be used
 * to observe the progress of the iteration.
 *
 *
 * <h3>Eigenvalue and condition number estimates</h3>
 *
 * This class can estimate eigenvalues and condition number during the
 * solution process. This is done by creating the Hessenberg matrix during the
 * inner iterations. The eigenvalues are estimated as the eigenvalues of the
 * Hessenberg matrix and the condition number is estimated as the ratio of the
 * largest and smallest singular value of the Hessenberg matrix. The estimates
 * can be obtained by connecting a function as a slot using @p
 * connect_condition_number_slot and @p connect_eigenvalues_slot. These slots
 * will then be called from the solver with the estimates as argument.
 */
template <typename VectorType = Vector<double>>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
class SolverGMRES : public SolverBase<VectorType>
{
public:
  /**
   * Standardized data struct to pipe additional data to the solver.
   */
  struct AdditionalData
  {
    /**
     * Constructor. By default, set the size of the Arnoldi basis to 30. Also
     * set preconditioning from left, the residual of the stopping criterion
     * to the default residual, and re-orthogonalization only if
     * necessary. Also, the batched mode with reduced functionality to track
     * information is disabled by default. Finally, the default
     * orthogonalization algorithm is the classical Gram-Schmidt method with
     * delayed reorthogonalization, which combines stability with fast
     * execution, especially in parallel.
     */
    explicit AdditionalData(const unsigned int max_basis_size        = 30,
                            const bool         right_preconditioning = false,
                            const bool         use_default_residual  = true,
                            const bool force_re_orthogonalization    = false,
                            const bool batched_mode                  = false,
                            const LinearAlgebra::OrthogonalizationStrategy
                              orthogonalization_strategy =
                                LinearAlgebra::OrthogonalizationStrategy::
                                  delayed_classical_gram_schmidt);

    /**
     * Maximum number of temporary vectors. Together with max_basis_size, this
     * parameter controls the size of the Arnoldi basis, which corresponds to
     * max_n_tmp_vectors-2 as used in previous versions of the deal.II
     * library. SolverGMRES assumes that there are at least three temporary
     * vectors, so this value must be greater than or equal to three. If both
     * this variable and max_basis_size are set to a non-zero value, the
     * choice in max_basis_size takes precedence.
     *
     * @deprecated Use max_basis_size instead.
     */
    unsigned int max_n_tmp_vectors;

    /**
     * Maximum size of the Arnoldi basis. SolverGMRES assumes that there is at
     * least one vector in the Arnoldi basis, so this value must be greater
     * than or equal to one. Note that whenever this variable is set to a
     * non-zero value, including the value set by the default constructor,
     * this variable takes precedence over max_n_tmp_vectors.
     */
    unsigned int max_basis_size;

    /**
     * Flag for right preconditioning.
     *
     * @note Change between left and right preconditioning will also change
     * the way residuals are evaluated. See the corresponding section in the
     * SolverGMRES.
     */
    bool right_preconditioning;

    /**
     * Flag for the default residual that is used to measure convergence.
     */
    bool use_default_residual;

    /**
     * Flag to force re-orthogonalization of orthonormal basis in every step.
     * If set to false, the solver automatically checks for loss of
     * orthogonality every 5 iterations and enables re-orthogonalization only
     * if necessary.
     */
    bool force_re_orthogonalization;

    /**
     * Flag to control whether a reduced mode of the solver should be
     * run. This is necessary when running (several) SolverGMRES instances
     * involving very small and cheap linear systems where the feedback from
     * all signals, eigenvalue computations, and log stream are disabled.
     */
    bool batched_mode;

    /**
     * Strategy to orthogonalize vectors.
     */
    LinearAlgebra::OrthogonalizationStrategy orthogonalization_strategy;
  };

  /**
   * Constructor.
   */
  SolverGMRES(SolverControl            &cn,
              VectorMemory<VectorType> &mem,
              const AdditionalData     &data = AdditionalData());

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  SolverGMRES(SolverControl &cn, const AdditionalData &data = AdditionalData());

  /**
   * The copy constructor is deleted.
   */
  SolverGMRES(const SolverGMRES<VectorType> &) = delete;

  /**
   * Solve the linear system $Ax=b$ for x.
   */
  template <typename MatrixType, typename PreconditionerType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_linear_operator_on<MatrixType, VectorType> &&
     concepts::is_linear_operator_on<PreconditionerType, VectorType>))
  void solve(const MatrixType         &A,
             VectorType               &x,
             const VectorType         &b,
             const PreconditionerType &preconditioner);

  /**
   * Connect a slot to retrieve the estimated condition number. Called on each
   * outer iteration if every_iteration=true, otherwise called once when
   * iterations are ended (i.e., either because convergence has been achieved,
   * or because divergence has been detected).
   */
  boost::signals2::connection
  connect_condition_number_slot(const std::function<void(double)> &slot,
                                const bool every_iteration = false);

  /**
   * Connect a slot to retrieve the estimated eigenvalues. Called on each
   * outer iteration if every_iteration=true, otherwise called once when
   * iterations are ended (i.e., either because convergence has been achieved,
   * or because divergence has been detected).
   */
  boost::signals2::connection
  connect_eigenvalues_slot(
    const std::function<void(const std::vector<std::complex<double>> &)> &slot,
    const bool every_iteration = false);

  /**
   * Connect a slot to retrieve the Hessenberg matrix obtained by the
   * projection of the initial matrix on the Krylov basis. Called on each
   * outer iteration if every_iteration=true, otherwise called once when
   * iterations are ended (i.e., either because convergence has been achieved,
   * or because divergence has been detected).
   */
  boost::signals2::connection
  connect_hessenberg_slot(
    const std::function<void(const FullMatrix<double> &)> &slot,
    const bool every_iteration = true);

  /**
   * Connect a slot to retrieve the basis vectors of the Krylov space
   * generated by the Arnoldi algorithm. Called at once when iterations
   * are completed (i.e., either because convergence has been achieved,
   * or because divergence has been detected).
   */
  boost::signals2::connection
  connect_krylov_space_slot(
    const std::function<
      void(const internal::SolverGMRESImplementation::TmpVectors<VectorType> &)>
      &slot);


  /**
   * Connect a slot to retrieve a notification when the vectors are
   * re-orthogonalized.
   */
  boost::signals2::connection
  connect_re_orthogonalization_slot(const std::function<void(int)> &slot);


  DeclException1(ExcTooFewTmpVectors,
                 int,
                 << "The number of temporary vectors you gave (" << arg1
                 << ") is too small. It should be at least 10 for "
                 << "any results, and much more for reasonable ones.");

protected:
  /**
   * Includes the maximum number of tmp vectors.
   */
  AdditionalData additional_data;

  /**
   * Signal used to retrieve the estimated condition number. Called once when
   * all iterations are ended.
   */
  boost::signals2::signal<void(double)> condition_number_signal;

  /**
   * Signal used to retrieve the estimated condition numbers. Called on each
   * outer iteration.
   */
  boost::signals2::signal<void(double)> all_condition_numbers_signal;

  /**
   * Signal used to retrieve the estimated eigenvalues. Called once when all
   * iterations are ended.
   */
  boost::signals2::signal<void(const std::vector<std::complex<double>> &)>
    eigenvalues_signal;

  /**
   * Signal used to retrieve the estimated eigenvalues. Called on each outer
   * iteration.
   */
  boost::signals2::signal<void(const std::vector<std::complex<double>> &)>
    all_eigenvalues_signal;

  /**
   * Signal used to retrieve the Hessenberg matrix. Called once when
   * all iterations are ended.
   */
  boost::signals2::signal<void(const FullMatrix<double> &)> hessenberg_signal;

  /**
   * Signal used to retrieve the Hessenberg matrix. Called on each outer
   * iteration.
   */
  boost::signals2::signal<void(const FullMatrix<double> &)>
    all_hessenberg_signal;

  /**
   * Signal used to retrieve the Krylov space basis vectors. Called once
   * when all iterations are ended.
   */
  boost::signals2::signal<void(
    const internal::SolverGMRESImplementation::TmpVectors<VectorType> &)>
    krylov_space_signal;

  /**
   * Signal used to retrieve a notification
   * when the vectors are re-orthogonalized.
   */
  boost::signals2::signal<void(int)> re_orthogonalize_signal;

  /**
   * A reference to the underlying SolverControl object. In the regular case,
   * this is not needed, as the signal from the base class is used, but the
   * batched variant cannot use those mechanisms due to the high costs.
   */
  SolverControl &solver_control;

  /**
   * Implementation of the computation of the norm of the residual.
   */
  virtual double
  criterion();

  /**
   * Estimates the eigenvalues from the Hessenberg matrix, H_orig, generated
   * during the inner iterations for @p n vectors in total. Uses these
   * estimates to compute the condition number. Calls the signals
   * eigenvalues_signal and cond_signal with these estimates as arguments.
   */
  static void
  compute_eigs_and_cond(
    const FullMatrix<double> &H_orig,
    const unsigned int        n,
    const boost::signals2::signal<
      void(const std::vector<std::complex<double>> &)> &eigenvalues_signal,
    const boost::signals2::signal<void(const FullMatrix<double> &)>
                                                &hessenberg_signal,
    const boost::signals2::signal<void(double)> &cond_signal);

  /**
   * Class that performs the actual orthogonalization process and solves the
   * projected linear system.
   */
  internal::SolverGMRESImplementation::ArnoldiProcess<
    typename VectorType::value_type>
    arnoldi_process;
};


/**
 * Implementation of the multiple preconditioned generalized minimal
 * residual method (MPGMRES).
 *
 * This method is a variant of the flexible GMRES, utilizing $N$
 * preconditioners to search for a solution within a multi-Krylov space.
 * These spaces are characterized by by all possible $N$-variate,
 * non-commuting polynomials of the preconditioners and system matrix
 * applied to a residual up to some fixed degree. In contrast, the flexible
 * GMRES method implemented in SolverFGMRES constructs only one "Krylov"
 * subspace, which is formed by univariate polynomials in one
 * preconditioner that may change at each iteration.
 *
 * We implement two strategies, a "full" and a "truncated" MPGMRES version
 * that differ in how they construct Krylov subspaces. The full MPGMRES
 * version constructs a Krylov subspace for every possible combination of
 * preconditioner application to the initial residual. For two
 * preconditioners $P_1$, $P_2$ this looks as follows:
 * @f{align*}{
 *   r, P_1r, P_2r, P_1AP_1r, P_2AP_1r, P_1AP_2r, P_2AP_2r, P_1AP_1AP_1r,
 *   P_2AP_1AP_1r, P_1AP_2AP_1r, P_2AP_2AP_1r, P_1AP_1AP_2r, P_2AP_1AP_2r,
 *   P_1AP_2AP_2r, P_2AP_2AP_2r, \ldots
 * @f}
 * The truncated version constructs independent Krylov subspaces by
 * dropping all "mixing" terms in the series expansion. For the example
 * with two preconditioners $P_1$, $P_2$ this looks as follows:
 * @f{align*}{
 *   r, P_1r, P_2r, P_1AP_1r, P_2AP_2r, P_1AP_1AP_1r, P_2AP_2AP_2r, \ldots
 * @f}
 * For reference, FGMRES uses the following construction:
 * @f{align*}{
 *   r, P_1r, P_2AP_1r, P_1AP_2AP_1r, \ldots
 * @f}
 * By default the truncated variant is used. You can switch to the full
 * version by setting the
 * AdditionalData::use_truncated_mpgmres_strategy option in the
 * AdditionalData object to false.
 *
 * For more details see @cite Greif2017.
 *
 * @note The present implementation of MPGMRES differs from the one
 * outlined in @cite Greif2017 in how one iteration is defined. Our
 * implementation constructs the search space one vector at a time,
 * producing a new iterate with each addition. In contrast, the routine
 * described in @cite Greif2017 constructs an iteration step by combining
 * all possible preconditioner applications corresponding to the total
 * polynomial degree of each multi-Krylov space. For the full MPGMRES
 * strategy this results in an exponential increase of possible
 * preconditioner applications that have to be computed for reach iteration
 * cycle; see Section 2.4 in @cite Greif2017.
 *
 * @note This method always uses right preconditioning, as opposed to
 * SolverGMRES, which allows the user to choose between left and right
 * preconditioning.
 *
 * @note The MPGMRES implementation needs two vectors in each iteration
 * steps yielding a total of
 * <tt>2*SolverMPGMRES::%AdditionalData::%max_basis_size+1</tt> auxiliary
 * vectors. Otherwise, FGMRES requires roughly the same number of
 * operations per iteration compared to GMRES, except one application of
 * the preconditioner less at each restart and at the end of solve().
 */

template <typename VectorType = Vector<double>>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
class SolverMPGMRES : public SolverBase<VectorType>
{
public:
  /**
   * Standardized data struct to pipe additional data to the solver.
   */
  struct AdditionalData
  {
    /**
     * Constructor. By default, set the maximum basis size to 30.
     */
    explicit AdditionalData(const unsigned int max_basis_size = 30,
                            const LinearAlgebra::OrthogonalizationStrategy
                              orthogonalization_strategy =
                                LinearAlgebra::OrthogonalizationStrategy::
                                  delayed_classical_gram_schmidt,
                            const bool use_truncated_mpgmres_strategy = true)
      : max_basis_size(max_basis_size)
      , orthogonalization_strategy(orthogonalization_strategy)
      , use_truncated_mpgmres_strategy(use_truncated_mpgmres_strategy)
    {}

    /**
     * Maximum basis size.
     */
    unsigned int max_basis_size;

    /**
     * Strategy to orthogonalize vectors.
     */
    LinearAlgebra::OrthogonalizationStrategy orthogonalization_strategy;

    /**
     * If set to true (the default) a "truncated" search space is
     * constructed consisting of the span of independent Krylov space
     * associated with each preconditioner. If set to false, the full
     * MPGMRES strategy for constructing the search space is used. This
     * space consists of all possible combinations of iterative
     * preconditioner applications; see the documentation of SolverMPGMRES
     * for details.
     */
    bool use_truncated_mpgmres_strategy;
  };

  /**
   * Constructor.
   */
  SolverMPGMRES(SolverControl            &cn,
                VectorMemory<VectorType> &mem,
                const AdditionalData     &data = AdditionalData());

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  SolverMPGMRES(SolverControl        &cn,
                const AdditionalData &data = AdditionalData());

  /**
   * Solve the linear system $Ax=b$ for x.
   */
  template <typename MatrixType, typename... PreconditionerTypes>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_linear_operator_on<MatrixType, VectorType> &&
     (concepts::is_linear_operator_on<PreconditionerTypes, VectorType> && ...)))
  void solve(const MatrixType &A,
             VectorType       &x,
             const VectorType &b,
             const PreconditionerTypes &...preconditioners);

protected:
  /**
   * Indexing strategy to construct the search space.
   *
   * This enum class is internally used in the implementation of the
   * MPGMRES algorithm to switch between the strategies.
   */
  enum class IndexingStrategy
  {
    fgmres,
    truncated_mpgmres,
    full_mpgmres,
  };

  /**
   * Solve the linear system $Ax=b$ for x.
   */
  template <typename MatrixType, typename... PreconditionerTypes>
  void
  solve_internal(const MatrixType       &A,
                 VectorType             &x,
                 const VectorType       &b,
                 const IndexingStrategy &indexing_strategy,
                 const PreconditionerTypes &...preconditioners);

private:
  /**
   * Additional flags.
   */
  AdditionalData additional_data;

  /**
   * Class that performs the actual orthogonalization process and solves the
   * projected linear system.
   */
  internal::SolverGMRESImplementation::ArnoldiProcess<
    typename VectorType::value_type>
    arnoldi_process;
};



/**
 * Implementation of the generalized minimal residual method with flexible
 * preconditioning (flexible GMRES or FGMRES).
 *
 * This flexible version of the GMRES method allows for the use of a
 * different preconditioner in each iteration step; in particular,
 * this also allows for the use of preconditioners that are not linear
 * operators. Therefore, it is also more robust with respect to
 * inaccurate evaluation of the preconditioner.  An important
 * application is the use of a Krylov space method inside the
 * preconditioner with low solver tolerance. See the documentation of
 * the SolverGMRES class for an elaboration of the issues involved.
 *
 * For more details see @cite Saad1991.
 *
 * @note This method always uses right preconditioning, as opposed to
 * SolverGMRES, which allows the user to choose between left and right
 * preconditioning.
 *
 * @note FGMRES needs two vectors in each iteration steps yielding a total
 * of <tt>2*SolverFGMRES::%AdditionalData::%max_basis_size+1</tt> auxiliary
 * vectors. Otherwise, FGMRES requires roughly the same number of
 * operations per iteration compared to GMRES, except one application of
 * the preconditioner less at each restart and at the end of solve().
 */
template <typename VectorType = Vector<double>>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
class SolverFGMRES : public SolverMPGMRES<VectorType>
{
public:
  /**
   * Standardized data struct to pipe additional data to the solver.
   */
  struct AdditionalData
  {
    /**
     * Constructor. By default, set the maximum basis size to 30.
     */
    explicit AdditionalData( //
      const unsigned int max_basis_size = 30,
      const LinearAlgebra::OrthogonalizationStrategy
        orthogonalization_strategy = LinearAlgebra::OrthogonalizationStrategy::
          delayed_classical_gram_schmidt)
      : max_basis_size(max_basis_size)
      , orthogonalization_strategy(orthogonalization_strategy)
    {}

    /**
     * Maximum basis size.
     */
    unsigned int max_basis_size;

    /**
     * Strategy to orthogonalize vectors.
     */
    LinearAlgebra::OrthogonalizationStrategy orthogonalization_strategy;
  };


  /**
   * Constructor.
   */
  SolverFGMRES(SolverControl            &cn,
               VectorMemory<VectorType> &mem,
               const AdditionalData     &data = AdditionalData());

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  SolverFGMRES(SolverControl        &cn,
               const AdditionalData &data = AdditionalData());

  /**
   * Solve the linear system $Ax=b$ for x.
   */
  template <typename MatrixType, typename... PreconditionerTypes>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_linear_operator_on<MatrixType, VectorType> &&
     (concepts::is_linear_operator_on<PreconditionerTypes, VectorType> && ...)))
  void solve(const MatrixType &A,
             VectorType       &x,
             const VectorType &b,
             const PreconditionerTypes &...preconditioners);
};

/** @} */


/* --------------------- Inline and template functions ------------------- */

#ifndef DOXYGEN

template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
inline SolverGMRES<VectorType>::AdditionalData::AdditionalData(
  const unsigned int                             max_basis_size,
  const bool                                     right_preconditioning,
  const bool                                     use_default_residual,
  const bool                                     force_re_orthogonalization,
  const bool                                     batched_mode,
  const LinearAlgebra::OrthogonalizationStrategy orthogonalization_strategy)
  : max_n_tmp_vectors(0)
  , max_basis_size(max_basis_size)
  , right_preconditioning(right_preconditioning)
  , use_default_residual(use_default_residual)
  , force_re_orthogonalization(force_re_orthogonalization)
  , batched_mode(batched_mode)
  , orthogonalization_strategy(orthogonalization_strategy)
{
  Assert(max_basis_size >= 1,
         ExcMessage("SolverGMRES needs at least one vector in the "
                    "Arnoldi basis."));
}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
SolverGMRES<VectorType>::SolverGMRES(SolverControl            &cn,
                                     VectorMemory<VectorType> &mem,
                                     const AdditionalData     &data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
  , solver_control(cn)
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
SolverGMRES<VectorType>::SolverGMRES(SolverControl        &cn,
                                     const AdditionalData &data)
  : SolverBase<VectorType>(cn)
  , additional_data(data)
  , solver_control(cn)
{}



namespace internal
{
  namespace SolverGMRESImplementation
  {
    template <typename VectorType>
    inline TmpVectors<VectorType>::TmpVectors(const unsigned int max_size,
                                              VectorMemory<VectorType> &vmem)
      : mem(vmem)
      , data(max_size)
    {}



    template <typename VectorType>
    inline VectorType &
    TmpVectors<VectorType>::operator[](const unsigned int i) const
    {
      AssertIndexRange(i, data.size());

      Assert(data[i] != nullptr, ExcNotInitialized());
      return *data[i];
    }



    template <typename VectorType>
    inline VectorType &
    TmpVectors<VectorType>::operator()(const unsigned int i,
                                       const VectorType  &temp)
    {
      AssertIndexRange(i, data.size());
      if (data[i] == nullptr)
        {
          data[i] = std::move(typename VectorMemory<VectorType>::Pointer(mem));
          data[i]->reinit(temp, true);
        }
      return *data[i];
    }



    template <typename VectorType>
    unsigned int
    TmpVectors<VectorType>::size() const
    {
      return (data.size() > 0 ? data.size() - 1 : 0);
    }



    template <typename VectorType, typename Enable = void>
    struct is_dealii_compatible_vector;

    template <typename VectorType>
    struct is_dealii_compatible_vector<
      VectorType,
      std::enable_if_t<!internal::is_block_vector<VectorType>>>
    {
      static constexpr bool value =
        std::is_same_v<
          VectorType,
          LinearAlgebra::distributed::Vector<typename VectorType::value_type,
                                             MemorySpace::Host>> ||
        std::is_same_v<VectorType, Vector<typename VectorType::value_type>>;
    };



    template <typename VectorType>
    struct is_dealii_compatible_vector<
      VectorType,
      std::enable_if_t<internal::is_block_vector<VectorType>>>
    {
      static constexpr bool value =
        std::is_same_v<
          typename VectorType::BlockType,
          LinearAlgebra::distributed::Vector<typename VectorType::value_type,
                                             MemorySpace::Host>> ||
        std::is_same_v<VectorType, Vector<typename VectorType::value_type>>;
    };



    template <typename VectorType,
              std::enable_if_t<!IsBlockVector<VectorType>::value, VectorType>
                * = nullptr>
    unsigned int
    n_blocks(const VectorType &)
    {
      return 1;
    }



    template <typename VectorType,
              std::enable_if_t<IsBlockVector<VectorType>::value, VectorType> * =
                nullptr>
    unsigned int
    n_blocks(const VectorType &vector)
    {
      return vector.n_blocks();
    }



    template <typename VectorType,
              std::enable_if_t<!IsBlockVector<VectorType>::value, VectorType>
                * = nullptr>
    VectorType &
    block(VectorType &vector, const unsigned int b)
    {
      AssertDimension(b, 0);
      return vector;
    }



    template <typename VectorType,
              std::enable_if_t<!IsBlockVector<VectorType>::value, VectorType>
                * = nullptr>
    const VectorType &
    block(const VectorType &vector, const unsigned int b)
    {
      AssertDimension(b, 0);
      return vector;
    }



    template <typename VectorType,
              std::enable_if_t<IsBlockVector<VectorType>::value, VectorType> * =
                nullptr>
    typename VectorType::BlockType &
    block(VectorType &vector, const unsigned int b)
    {
      return vector.block(b);
    }



    template <typename VectorType,
              std::enable_if_t<IsBlockVector<VectorType>::value, VectorType> * =
                nullptr>
    const typename VectorType::BlockType &
    block(const VectorType &vector, const unsigned int b)
    {
      return vector.block(b);
    }



    template <bool delayed_reorthogonalization,
              typename VectorType,
              std::enable_if_t<!is_dealii_compatible_vector<VectorType>::value,
                               VectorType> * = nullptr>
    void
    Tvmult_add(const unsigned int            n,
               const VectorType             &vv,
               const TmpVectors<VectorType> &orthogonal_vectors,
               Vector<double>               &h,
               std::vector<const typename VectorType::value_type *> &)
    {
      for (unsigned int i = 0; i < n; ++i)
        {
          h(i) += vv * orthogonal_vectors[i];
          if (delayed_reorthogonalization)
            h(n + i) += orthogonal_vectors[i] * orthogonal_vectors[n - 1];
        }
      if (delayed_reorthogonalization)
        h(n + n) += vv * vv;
    }



    // worker method for deal.II's vector types implemented in .cc file
    template <bool delayed_reorthogonalization, typename Number>
    void
    do_Tvmult_add(const unsigned int                 n_vectors,
                  const std::size_t                  locally_owned_size,
                  const Number                      *current_vector,
                  const std::vector<const Number *> &orthogonal_vectors,
                  Vector<double>                    &b);



    template <bool delayed_reorthogonalization,
              typename VectorType,
              std::enable_if_t<is_dealii_compatible_vector<VectorType>::value,
                               VectorType> * = nullptr>
    void
    Tvmult_add(
      const unsigned int                                    n,
      const VectorType                                     &vv,
      const TmpVectors<VectorType>                         &orthogonal_vectors,
      Vector<double>                                       &h,
      std::vector<const typename VectorType::value_type *> &vector_ptrs)
    {
      for (unsigned int b = 0; b < n_blocks(vv); ++b)
        {
          vector_ptrs.resize(n);
          for (unsigned int i = 0; i < n; ++i)
            vector_ptrs[i] = block(orthogonal_vectors[i], b).begin();

          do_Tvmult_add<delayed_reorthogonalization>(n,
                                                     block(vv, b).end() -
                                                       block(vv, b).begin(),
                                                     block(vv, b).begin(),
                                                     vector_ptrs,
                                                     h);
        }

      Utilities::MPI::sum(h, block(vv, 0).get_mpi_communicator(), h);
    }



    template <bool delayed_reorthogonalization,
              typename VectorType,
              std::enable_if_t<!is_dealii_compatible_vector<VectorType>::value,
                               VectorType> * = nullptr>
    double
    subtract_and_norm(const unsigned int            n,
                      const TmpVectors<VectorType> &orthogonal_vectors,
                      const Vector<double>         &h,
                      VectorType                   &vv,
                      std::vector<const typename VectorType::value_type *> &)
    {
      Assert(n > 0, ExcInternalError());

      VectorType &last_vector =
        const_cast<VectorType &>(orthogonal_vectors[n - 1]);
      for (unsigned int i = 0; i < n - 1; ++i)
        {
          if (delayed_reorthogonalization && i + 2 < n)
            last_vector.add(-h(n + i), orthogonal_vectors[i]);
          vv.add(-h(i), orthogonal_vectors[i]);
        }

      if (delayed_reorthogonalization)
        {
          if (n > 1)
            last_vector.sadd(1. / h(n + n - 1),
                             -h(n + n - 2) / h(n + n - 1),
                             orthogonal_vectors[n - 2]);

          // h(n + n) is lucky breakdown
          const double scaling_factor_vv = h(n + n) > 0.0 ?
                                             1. / (h(n + n - 1) * h(n + n)) :
                                             1. / (h(n + n - 1) * h(n + n - 1));
          vv.sadd(scaling_factor_vv,
                  -h(n - 1) * scaling_factor_vv,
                  last_vector);

          // the delayed reorthogonalization computes the norm from other
          // quantities
          return std::numeric_limits<double>::signaling_NaN();
        }
      else
        return std::sqrt(
          vv.add_and_dot(-h(n - 1), orthogonal_vectors[n - 1], vv));
    }



    // worker method for deal.II's vector types implemented in .cc file
    template <bool delayed_reorthogonalization, typename Number>
    double
    do_subtract_and_norm(const unsigned int                 n_vectors,
                         const std::size_t                  locally_owned_size,
                         const std::vector<const Number *> &orthogonal_vectors,
                         const Vector<double>              &h,
                         Number                            *current_vector);



    template <bool delayed_reorthogonalization,
              typename VectorType,
              std::enable_if_t<is_dealii_compatible_vector<VectorType>::value,
                               VectorType> * = nullptr>
    double
    subtract_and_norm(
      const unsigned int                                    n,
      const TmpVectors<VectorType>                         &orthogonal_vectors,
      const Vector<double>                                 &h,
      VectorType                                           &vv,
      std::vector<const typename VectorType::value_type *> &vector_ptrs)
    {
      double norm_vv_temp = 0.0;

      for (unsigned int b = 0; b < n_blocks(vv); ++b)
        {
          vector_ptrs.resize(n);
          for (unsigned int i = 0; i < n; ++i)
            vector_ptrs[i] = block(orthogonal_vectors[i], b).begin();

          norm_vv_temp += do_subtract_and_norm<delayed_reorthogonalization>(
            n,
            block(vv, b).end() - block(vv, b).begin(),
            vector_ptrs,
            h,
            block(vv, b).begin());
        }

      return std::sqrt(
        Utilities::MPI::sum(norm_vv_temp, block(vv, 0).get_mpi_communicator()));
    }



    template <typename VectorType,
              std::enable_if_t<!is_dealii_compatible_vector<VectorType>::value,
                               VectorType> * = nullptr>
    void
    add(VectorType                   &p,
        const unsigned int            n,
        const Vector<double>         &h,
        const TmpVectors<VectorType> &tmp_vectors,
        const bool                    zero_out,
        std::vector<const typename VectorType::value_type *> &)
    {
      if (zero_out)
        p.equ(h(0), tmp_vectors[0]);
      else
        p.add(h(0), tmp_vectors[0]);

      for (unsigned int i = 1; i < n; ++i)
        p.add(h(i), tmp_vectors[i]);
    }



    // worker method for deal.II's vector types implemented in .cc file
    template <typename Number>
    void
    do_add(const unsigned int                 n_vectors,
           const std::size_t                  locally_owned_size,
           const std::vector<const Number *> &tmp_vectors,
           const Vector<double>              &h,
           const bool                         zero_out,
           Number                            *output);



    template <typename VectorType,
              std::enable_if_t<is_dealii_compatible_vector<VectorType>::value,
                               VectorType> * = nullptr>
    void
    add(VectorType                                           &p,
        const unsigned int                                    n,
        const Vector<double>                                 &h,
        const TmpVectors<VectorType>                         &tmp_vectors,
        const bool                                            zero_out,
        std::vector<const typename VectorType::value_type *> &vector_ptrs)
    {
      for (unsigned int b = 0; b < n_blocks(p); ++b)
        {
          vector_ptrs.resize(n);
          for (unsigned int i = 0; i < n; ++i)
            vector_ptrs[i] = block(tmp_vectors[i], b).begin();
          do_add(n,
                 block(p, b).end() - block(p, b).begin(),
                 vector_ptrs,
                 h,
                 zero_out,
                 block(p, b).begin());
        }
    }



    template <typename Number>
    inline void
    ArnoldiProcess<Number>::initialize(
      const LinearAlgebra::OrthogonalizationStrategy orthogonalization_strategy,
      const unsigned int                             basis_size,
      const bool                                     force_reorthogonalization)
    {
      this->orthogonalization_strategy = orthogonalization_strategy;
      this->do_reorthogonalization     = force_reorthogonalization;

      hessenberg_matrix.reinit(basis_size + 1, basis_size);
      triangular_matrix.reinit(basis_size + 1, basis_size, true);

      // some additional vectors, also used in the orthogonalization
      projected_rhs.reinit(basis_size + 1, true);
      givens_rotations.reserve(basis_size);

      if (orthogonalization_strategy ==
          LinearAlgebra::OrthogonalizationStrategy::
            delayed_classical_gram_schmidt)
        h.reinit(2 * basis_size + 3);
      else
        h.reinit(basis_size + 1);
    }



    template <typename Number>
    template <typename VectorType>
    inline double
    ArnoldiProcess<Number>::orthonormalize_nth_vector(
      const unsigned int                        n,
      TmpVectors<VectorType>                   &orthogonal_vectors,
      const unsigned int                        accumulated_iterations,
      const boost::signals2::signal<void(int)> &reorthogonalize_signal)
    {
      AssertIndexRange(n, hessenberg_matrix.m());
      AssertIndexRange(n, orthogonal_vectors.size() + 1);

      VectorType &vv = orthogonal_vectors[n];

      double residual_estimate = std::numeric_limits<double>::signaling_NaN();
      if (n == 0)
        {
          givens_rotations.clear();
          residual_estimate = vv.l2_norm();
          if (residual_estimate != 0.)
            vv /= residual_estimate;
          projected_rhs(0) = residual_estimate;
        }
      else if (orthogonalization_strategy ==
               LinearAlgebra::OrthogonalizationStrategy::
                 delayed_classical_gram_schmidt)
        {
          // The algorithm implemented in the following few lines is algorithm
          // 4 of Bielich et al. (2022).

          // To avoid un-scaled numbers as appearing with the original
          // algorithm by Bielich et al., we use a preliminary scaling of the
          // last vector. This will be corrected in the delayed step.
          const double previous_scaling = n > 0 ? h(n + n - 2) : 1.;

          // Reset h to zero
          h.reinit(n + n + 1);

          // global reduction
          Tvmult_add<true>(n, vv, orthogonal_vectors, h, vector_ptrs);

          // delayed correction terms
          double tmp = 0;
          for (unsigned int i = 0; i < n - 1; ++i)
            tmp += h(n + i) * h(n + i);
          const double alpha_j = h(n + n - 1) > tmp ?
                                   std::sqrt(h(n + n - 1) - tmp) :
                                   std::sqrt(h(n + n - 1));
          h(n + n - 1)         = alpha_j;

          tmp = 0;
          for (unsigned int i = 0; i < n - 1; ++i)
            tmp += h(i) * h(n + i);
          h(n - 1) = (h(n - 1) - tmp) / alpha_j;

          // representation of H(j-1)
          if (n > 1)
            {
              for (unsigned int i = 0; i < n - 1; ++i)
                hessenberg_matrix(i, n - 2) += h(n + i) * previous_scaling;
              hessenberg_matrix(n - 1, n - 2) = alpha_j * previous_scaling;
            }
          for (unsigned int i = 0; i < n; ++i)
            {
              double sum = 0;
              for (unsigned int j = (i == 0 ? 0 : i - 1); j < n - 1; ++j)
                sum += hessenberg_matrix(i, j) * h(n + j);
              hessenberg_matrix(i, n - 1) = (h(i) - sum) / alpha_j;
            }

          // compute norm estimate for approximate convergence criterion
          // (value of norm to be corrected in next iteration)
          double sum = 0;
          for (unsigned int i = 0; i < n - 1; ++i)
            sum += h(i) * h(i);
          sum += (2. - 1.) * h(n - 1) * h(n - 1);
          hessenberg_matrix(n, n - 1) =
            std::sqrt(std::abs(h(n + n) - sum)) / alpha_j;

          // projection and delayed reorthogonalization. We scale the vector
          // vv here by the preliminary norm to avoid working with too large
          // values and correct the actual norm in the Hessenberg matrix in
          // high precision in the next iteration.
          h(n + n) = hessenberg_matrix(n, n - 1);
          subtract_and_norm<true>(n, orthogonal_vectors, h, vv, vector_ptrs);

          // transform new column of upper Hessenberg matrix into upper
          // triangular form by computing the respective factor
          residual_estimate = do_givens_rotation(
            true, n - 2, triangular_matrix, givens_rotations, projected_rhs);
        }
      else
        {
          // need initial norm for detection of re-orthogonalization, see below
          double     norm_vv       = 0.0;
          double     norm_vv_start = 0;
          const bool consider_reorthogonalize =
            (do_reorthogonalization == false) && (n % 5 == 0);
          if (consider_reorthogonalize)
            norm_vv_start = vv.l2_norm();

          // Reset h to zero
          h.reinit(n);

          // run two loops with index 0: orthogonalize, 1: reorthogonalize
          for (unsigned int c = 0; c < 2; ++c)
            {
              // Orthogonalization
              if (orthogonalization_strategy ==
                  LinearAlgebra::OrthogonalizationStrategy::
                    modified_gram_schmidt)
                {
                  double htmp = vv * orthogonal_vectors[0];
                  h(0) += htmp;
                  for (unsigned int i = 1; i < n; ++i)
                    {
                      htmp = vv.add_and_dot(-htmp,
                                            orthogonal_vectors[i - 1],
                                            orthogonal_vectors[i]);
                      h(i) += htmp;
                    }

                  norm_vv = std::sqrt(
                    vv.add_and_dot(-htmp, orthogonal_vectors[n - 1], vv));
                }
              else if (orthogonalization_strategy ==
                       LinearAlgebra::OrthogonalizationStrategy::
                         classical_gram_schmidt)
                {
                  Tvmult_add<false>(n, vv, orthogonal_vectors, h, vector_ptrs);
                  norm_vv = subtract_and_norm<false>(
                    n, orthogonal_vectors, h, vv, vector_ptrs);
                }
              else
                {
                  AssertThrow(false, ExcNotImplemented());
                }

              if (c == 1)
                break; // reorthogonalization already performed -> finished

              // Re-orthogonalization if loss of orthogonality detected. For the
              // test, use a strategy discussed in C. T. Kelley, Iterative
              // Methods for Linear and Nonlinear Equations, SIAM, Philadelphia,
              // 1995: Compare the norm of vv after orthogonalization with its
              // norm when starting the orthogonalization. If vv became very
              // small (here: less than the square root of the machine precision
              // times 10), it is almost in the span of the previous vectors,
              // which indicates loss of precision.
              if (consider_reorthogonalize)
                {
                  if (norm_vv >
                      10. * norm_vv_start *
                        std::sqrt(std::numeric_limits<
                                  typename VectorType::value_type>::epsilon()))
                    break;

                  else
                    {
                      do_reorthogonalization = true;
                      if (!reorthogonalize_signal.empty())
                        reorthogonalize_signal(accumulated_iterations);
                    }
                }

              if (do_reorthogonalization == false)
                break; // no reorthogonalization needed -> finished
            }

          for (unsigned int i = 0; i < n; ++i)
            hessenberg_matrix(i, n - 1) = h(i);
          hessenberg_matrix(n, n - 1) = norm_vv;

          // norm_vv is a lucky breakdown, the solver will reach convergence,
          // but we must not divide by zero here.
          if (norm_vv != 0)
            vv /= norm_vv;

          residual_estimate = do_givens_rotation(
            false, n - 1, triangular_matrix, givens_rotations, projected_rhs);
        }

      return residual_estimate;
    }



    template <typename Number>
    inline double
    ArnoldiProcess<Number>::do_givens_rotation(
      const bool                              delayed_reorthogonalization,
      const int                               col,
      FullMatrix<double>                     &matrix,
      std::vector<std::pair<double, double>> &rotations,
      Vector<double>                         &rhs)
    {
      // for the delayed orthogonalization, we can only compute the column of
      // the previous iteration (as there will be correction terms added to the
      // present column for stability reasons), but we still want to compute
      // the residual estimate from the accumulated work; we therefore perform
      // givens rotations on two columns simultaneously
      if (delayed_reorthogonalization)
        {
          if (col >= 0)
            {
              AssertDimension(rotations.size(), static_cast<std::size_t>(col));
              matrix(0, col) = hessenberg_matrix(0, col);
            }
          double H_next = hessenberg_matrix(0, col + 1);
          for (int i = 0; i < col; ++i)
            {
              const double c   = rotations[i].first;
              const double s   = rotations[i].second;
              const double Hi  = matrix(i, col);
              const double Hi1 = hessenberg_matrix(i + 1, col);
              H_next = -s * H_next + c * hessenberg_matrix(i + 1, col + 1);
              matrix(i, col)     = c * Hi + s * Hi1;
              matrix(i + 1, col) = -s * Hi + c * Hi1;
            }

          if (col >= 0)
            {
              const double H_col1 = hessenberg_matrix(col + 1, col);
              const double H_col  = matrix(col, col);
              const double r = 1. / std::sqrt(H_col * H_col + H_col1 * H_col1);
              rotations.emplace_back(H_col * r, H_col1 * r);
              matrix(col, col) =
                rotations[col].first * H_col + rotations[col].second * H_col1;

              rhs(col + 1) = -rotations[col].second * rhs(col);
              rhs(col) *= rotations[col].first;

              H_next =
                -rotations[col].second * H_next +
                rotations[col].first * hessenberg_matrix(col + 1, col + 1);
            }

          const double H_last = hessenberg_matrix(col + 2, col + 1);
          const double r = 1. / std::sqrt(H_next * H_next + H_last * H_last);
          return std::abs(H_last * r * rhs(col + 1));
        }
      else
        {
          AssertDimension(rotations.size(), static_cast<std::size_t>(col));

          matrix(0, col) = hessenberg_matrix(0, col);
          for (int i = 0; i < col; ++i)
            {
              const double c     = rotations[i].first;
              const double s     = rotations[i].second;
              const double Hi    = matrix(i, col);
              const double Hi1   = hessenberg_matrix(i + 1, col);
              matrix(i, col)     = c * Hi + s * Hi1;
              matrix(i + 1, col) = -s * Hi + c * Hi1;
            }

          const double Hi  = matrix(col, col);
          const double Hi1 = hessenberg_matrix(col + 1, col);
          const double r   = 1. / std::sqrt(Hi * Hi + Hi1 * Hi1);
          rotations.emplace_back(Hi * r, Hi1 * r);
          matrix(col, col) =
            rotations[col].first * Hi + rotations[col].second * Hi1;

          rhs(col + 1) = -rotations[col].second * rhs(col);
          rhs(col) *= rotations[col].first;

          return std::abs(rhs(col + 1));
        }
    }



    template <typename Number>
    inline const Vector<double> &
    ArnoldiProcess<Number>::solve_projected_system(
      const bool orthogonalization_finished)
    {
      FullMatrix<double>  tmp_triangular_matrix;
      Vector<double>      tmp_rhs;
      FullMatrix<double> *matrix = &triangular_matrix;
      Vector<double>     *rhs    = &projected_rhs;
      unsigned int        n      = givens_rotations.size();

      // If we solve with the delayed orthogonalization, we still need to
      // perform the elimination of the last column before we can solve the
      // projected system. We distinguish two cases, one where the
      // orthogonalization has finished (i.e., end of inner iteration in
      // GMRES) and we can safely overwrite the content of the tridiagonal
      // matrix and right hand side, and the case during the inner iterations,
      // where we need to create copies of the matrices in the QR
      // decomposition as well as the right hand side.
      if (orthogonalization_strategy ==
          LinearAlgebra::OrthogonalizationStrategy::
            delayed_classical_gram_schmidt)
        {
          n += 1;
          if (!orthogonalization_finished)
            {
              tmp_triangular_matrix = triangular_matrix;
              tmp_rhs               = projected_rhs;
              std::vector<std::pair<double, double>> tmp_givens_rotations(
                givens_rotations);
              do_givens_rotation(false,
                                 givens_rotations.size(),
                                 tmp_triangular_matrix,
                                 tmp_givens_rotations,
                                 tmp_rhs);
              matrix = &tmp_triangular_matrix;
              rhs    = &tmp_rhs;
            }
          else
            do_givens_rotation(false,
                               givens_rotations.size(),
                               triangular_matrix,
                               givens_rotations,
                               projected_rhs);
        }

      // Now solve the triangular system by backward substitution
      projected_solution.reinit(n);
      for (int i = n - 1; i >= 0; --i)
        {
          double s = (*rhs)(i);
          for (unsigned int j = i + 1; j < n; ++j)
            s -= projected_solution(j) * (*matrix)(i, j);
          projected_solution(i) = s / (*matrix)(i, i);
          AssertIsFinite(projected_solution(i));
        }

      return projected_solution;
    }



    template <typename Number>
    inline const FullMatrix<double> &
    ArnoldiProcess<Number>::get_hessenberg_matrix() const
    {
      return hessenberg_matrix;
    }



    // A comparator for better printing eigenvalues
    inline bool
    complex_less_pred(const std::complex<double> &x,
                      const std::complex<double> &y)
    {
      return x.real() < y.real() ||
             (x.real() == y.real() && x.imag() < y.imag());
    }
  } // namespace SolverGMRESImplementation
} // namespace internal



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
inline void SolverGMRES<VectorType>::compute_eigs_and_cond(
  const FullMatrix<double> &H_orig,
  const unsigned int        n,
  const boost::signals2::signal<void(const std::vector<std::complex<double>> &)>
    &eigenvalues_signal,
  const boost::signals2::signal<void(const FullMatrix<double> &)>
                                              &hessenberg_signal,
  const boost::signals2::signal<void(double)> &cond_signal)
{
  // Avoid copying the Hessenberg matrix if it isn't needed.
  if ((!eigenvalues_signal.empty() || !hessenberg_signal.empty() ||
       !cond_signal.empty()) &&
      n > 0)
    {
      LAPACKFullMatrix<double> mat(n, n);
      for (unsigned int i = 0; i < n; ++i)
        for (unsigned int j = 0; j < n; ++j)
          mat(i, j) = H_orig(i, j);
      hessenberg_signal(H_orig);
      // Avoid computing eigenvalues if they are not needed.
      if (!eigenvalues_signal.empty())
        {
          // Copy mat so that we can compute svd below. Necessary since
          // compute_eigenvalues will leave mat in state
          // LAPACKSupport::unusable.
          LAPACKFullMatrix<double> mat_eig(mat);
          mat_eig.compute_eigenvalues();
          std::vector<std::complex<double>> eigenvalues(n);
          for (unsigned int i = 0; i < mat_eig.n(); ++i)
            eigenvalues[i] = mat_eig.eigenvalue(i);
          // Sort eigenvalues for nicer output.
          std::sort(eigenvalues.begin(),
                    eigenvalues.end(),
                    internal::SolverGMRESImplementation::complex_less_pred);
          eigenvalues_signal(eigenvalues);
        }
      // Calculate condition number, avoid calculating the svd if a slot
      // isn't connected. Need at least a 2-by-2 matrix to do the estimate.
      if (!cond_signal.empty() && (mat.n() > 1))
        {
          mat.compute_svd();
          double condition_number =
            mat.singular_value(0) / mat.singular_value(mat.n() - 1);
          cond_signal(condition_number);
        }
    }
}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
template <typename MatrixType, typename PreconditionerType>
DEAL_II_CXX20_REQUIRES(
  (concepts::is_linear_operator_on<MatrixType, VectorType> &&
   concepts::is_linear_operator_on<PreconditionerType, VectorType>))
void SolverGMRES<VectorType>::solve(const MatrixType         &A,
                                    VectorType               &x,
                                    const VectorType         &b,
                                    const PreconditionerType &preconditioner)
{
  std::unique_ptr<LogStream::Prefix> prefix;
  if (!additional_data.batched_mode)
    prefix = std::make_unique<LogStream::Prefix>("GMRES");

  // extra call to std::max to placate static analyzers: coverity rightfully
  // complains that data.max_n_tmp_vectors - 2 may overflow
  const unsigned int basis_size =
    (additional_data.max_basis_size > 0 ?
       additional_data.max_basis_size :
       std::max(additional_data.max_n_tmp_vectors, 3u) - 2);

  // Generate an object where basis vectors are stored.
  internal::SolverGMRESImplementation::TmpVectors<VectorType> basis_vectors(
    basis_size + 2, this->memory);

  // number of the present iteration; this number is not reset to zero upon a
  // restart
  unsigned int accumulated_iterations = 0;

  const bool do_eigenvalues =
    !additional_data.batched_mode &&
    (!condition_number_signal.empty() ||
     !all_condition_numbers_signal.empty() || !eigenvalues_signal.empty() ||
     !all_eigenvalues_signal.empty() || !hessenberg_signal.empty() ||
     !all_hessenberg_signal.empty());

  SolverControl::State iteration_state = SolverControl::iterate;
  double               res             = std::numeric_limits<double>::lowest();

  // switch to determine whether we want a left or a right preconditioner. at
  // present, left is default, but both ways are implemented
  const bool left_precondition = !additional_data.right_preconditioning;

  // Per default the left preconditioned GMRES uses the preconditioned
  // residual and the right preconditioned GMRES uses the unpreconditioned
  // residual as stopping criterion.
  const bool use_default_residual = additional_data.use_default_residual;

  // define an alias
  VectorType &p = basis_vectors(basis_size + 1, x);

  // Following vectors are needed when we are not using the default residuals
  // as stopping criterion
  typename VectorMemory<VectorType>::Pointer r;
  typename VectorMemory<VectorType>::Pointer x_;
  if (!use_default_residual)
    {
      r  = std::move(typename VectorMemory<VectorType>::Pointer(this->memory));
      x_ = std::move(typename VectorMemory<VectorType>::Pointer(this->memory));
      r->reinit(x);
      x_->reinit(x);
    }

  arnoldi_process.initialize(additional_data.orthogonalization_strategy,
                             basis_size,
                             additional_data.force_re_orthogonalization);

  ///////////////////////////////////////////////////////////////////////////
  // outer iteration: loop until we either reach convergence or the maximum
  // number of iterations is exceeded. each cycle of this loop amounts to one
  // restart
  do
    {
      VectorType &v = basis_vectors(0, x);

      // Compute the preconditioned/unpreconditioned residual for left/right
      // preconditioning. If 'x' is the zero vector, then we can bypass the
      // full computation. But 'x' is only likely to be the zero vector if
      // that's what the user provided as the starting guess, so it's only
      // worth checking for this in the first iteration. (Calling all_zero()
      // costs as much in memory transfer and communication as computing the
      // norm of a vector.)
      if (left_precondition)
        {
          if (accumulated_iterations == 0 && x.all_zero())
            preconditioner.vmult(v, b);
          else
            {
              A.vmult(p, x);
              p.sadd(-1., 1., b);
              preconditioner.vmult(v, p);
            }
        }
      else
        {
          if (accumulated_iterations == 0 && x.all_zero())
            v = b;
          else
            {
              A.vmult(v, x);
              v.sadd(-1., 1., b);
            }
        }

      const double norm_v = arnoldi_process.orthonormalize_nth_vector(
        0, basis_vectors, accumulated_iterations, re_orthogonalize_signal);

      // check the residual here as well since it may be that we got the exact
      // (or an almost exact) solution vector at the outset. if we wouldn't
      // check here, the next scaling operation would produce garbage
      if (use_default_residual)
        {
          res = norm_v;
          if (additional_data.batched_mode)
            iteration_state = solver_control.check(accumulated_iterations, res);
          else
            iteration_state =
              this->iteration_status(accumulated_iterations, res, x);

          if (iteration_state != SolverControl::iterate)
            break;
        }
      else
        {
          deallog << "default_res=" << norm_v << std::endl;

          if (left_precondition)
            {
              A.vmult(*r, x);
              r->sadd(-1., 1., b);
            }
          else
            preconditioner.vmult(*r, v);

          res = r->l2_norm();
          if (additional_data.batched_mode)
            iteration_state = solver_control.check(accumulated_iterations, res);
          else
            iteration_state =
              this->iteration_status(accumulated_iterations, res, x);

          if (iteration_state != SolverControl::iterate)
            break;
        }

      // inner iteration doing at most as many steps as the size of the
      // Arnoldi basis
      unsigned int inner_iteration = 0;
      for (; (inner_iteration < basis_size &&
              iteration_state == SolverControl::iterate);
           ++inner_iteration)
        {
          ++accumulated_iterations;
          // yet another alias
          VectorType &vv = basis_vectors(inner_iteration + 1, x);

          if (left_precondition)
            {
              A.vmult(p, basis_vectors[inner_iteration]);
              preconditioner.vmult(vv, p);
            }
          else
            {
              preconditioner.vmult(p, basis_vectors[inner_iteration]);
              A.vmult(vv, p);
            }

          res =
            arnoldi_process.orthonormalize_nth_vector(inner_iteration + 1,
                                                      basis_vectors,
                                                      accumulated_iterations,
                                                      re_orthogonalize_signal);

          if (use_default_residual)
            {
              if (additional_data.batched_mode)
                iteration_state =
                  solver_control.check(accumulated_iterations, res);
              else
                iteration_state =
                  this->iteration_status(accumulated_iterations, res, x);
            }
          else
            {
              if (!additional_data.batched_mode)
                deallog << "default_res=" << res << std::endl;

              *x_ = x;
              const Vector<double> &projected_solution =
                arnoldi_process.solve_projected_system(false);

              if (left_precondition)
                for (unsigned int i = 0; i < inner_iteration + 1; ++i)
                  x_->add(projected_solution(i), basis_vectors[i]);
              else
                {
                  p = 0.;
                  for (unsigned int i = 0; i < inner_iteration + 1; ++i)
                    p.add(projected_solution(i), basis_vectors[i]);
                  preconditioner.vmult(*r, p);
                  x_->add(1., *r);
                };
              A.vmult(*r, *x_);
              r->sadd(-1., 1., b);

              // Now *r contains the unpreconditioned residual!!
              if (left_precondition)
                {
                  res = r->l2_norm();
                  iteration_state =
                    this->iteration_status(accumulated_iterations, res, x);
                }
              else
                {
                  preconditioner.vmult(*x_, *r);
                  res = x_->l2_norm();

                  if (additional_data.batched_mode)
                    iteration_state =
                      solver_control.check(accumulated_iterations, res);
                  else
                    iteration_state =
                      this->iteration_status(accumulated_iterations, res, x);
                }
            }
        }

      // end of inner iteration; now update the global solution vector x with
      // the solution of the projected system (least-squares solution)
      const Vector<double> &projected_solution =
        arnoldi_process.solve_projected_system(true);

      if (do_eigenvalues)
        compute_eigs_and_cond(arnoldi_process.get_hessenberg_matrix(),
                              inner_iteration,
                              all_eigenvalues_signal,
                              all_hessenberg_signal,
                              condition_number_signal);

      if (left_precondition)
        dealii::internal::SolverGMRESImplementation::add(
          x,
          inner_iteration,
          projected_solution,
          basis_vectors,
          false,
          arnoldi_process.vector_ptrs);
      else
        {
          dealii::internal::SolverGMRESImplementation::add(
            p,
            inner_iteration,
            projected_solution,
            basis_vectors,
            true,
            arnoldi_process.vector_ptrs);
          preconditioner.vmult(v, p);
          x.add(1., v);
        }

      // in the last round, print the eigenvalues from the last Arnoldi step
      if (iteration_state != SolverControl::iterate)
        {
          if (do_eigenvalues)
            compute_eigs_and_cond(arnoldi_process.get_hessenberg_matrix(),
                                  inner_iteration,
                                  eigenvalues_signal,
                                  hessenberg_signal,
                                  condition_number_signal);

          if (!additional_data.batched_mode && !krylov_space_signal.empty())
            krylov_space_signal(basis_vectors);

          // end of outer iteration. restart if no convergence and the number of
          // iterations is not exceeded
        }
    }
  while (iteration_state == SolverControl::iterate);

  // in case of failure: throw exception
  AssertThrow(iteration_state == SolverControl::success,
              SolverControl::NoConvergence(accumulated_iterations, res));
}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
boost::signals2::connection
  SolverGMRES<VectorType>::connect_condition_number_slot(
    const std::function<void(double)> &slot,
    const bool                         every_iteration)
{
  if (every_iteration)
    {
      return all_condition_numbers_signal.connect(slot);
    }
  else
    {
      return condition_number_signal.connect(slot);
    }
}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
boost::signals2::connection SolverGMRES<VectorType>::connect_eigenvalues_slot(
  const std::function<void(const std::vector<std::complex<double>> &)> &slot,
  const bool every_iteration)
{
  if (every_iteration)
    {
      return all_eigenvalues_signal.connect(slot);
    }
  else
    {
      return eigenvalues_signal.connect(slot);
    }
}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
boost::signals2::connection SolverGMRES<VectorType>::connect_hessenberg_slot(
  const std::function<void(const FullMatrix<double> &)> &slot,
  const bool                                             every_iteration)
{
  if (every_iteration)
    {
      return all_hessenberg_signal.connect(slot);
    }
  else
    {
      return hessenberg_signal.connect(slot);
    }
}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
boost::signals2::connection SolverGMRES<VectorType>::connect_krylov_space_slot(
  const std::function<void(
    const internal::SolverGMRESImplementation::TmpVectors<VectorType> &)> &slot)
{
  return krylov_space_signal.connect(slot);
}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
boost::signals2::connection
  SolverGMRES<VectorType>::connect_re_orthogonalization_slot(
    const std::function<void(int)> &slot)
{
  return re_orthogonalize_signal.connect(slot);
}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
double SolverGMRES<VectorType>::criterion()
{
  // dummy implementation. this function is not needed for the present
  // implementation of gmres
  DEAL_II_ASSERT_UNREACHABLE();
  return 0;
}



//----------------------------------------------------------------------//



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
SolverMPGMRES<VectorType>::SolverMPGMRES(SolverControl            &cn,
                                         VectorMemory<VectorType> &mem,
                                         const AdditionalData     &data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
SolverMPGMRES<VectorType>::SolverMPGMRES(SolverControl        &cn,
                                         const AdditionalData &data)
  : SolverBase<VectorType>(cn)
  , additional_data(data)
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
template <typename MatrixType, typename... PreconditionerTypes>
DEAL_II_CXX20_REQUIRES(
  (concepts::is_linear_operator_on<MatrixType, VectorType> &&
   (concepts::is_linear_operator_on<PreconditionerTypes, VectorType> && ...)))
void SolverMPGMRES<VectorType>::solve(
  const MatrixType &A,
  VectorType       &x,
  const VectorType &b,
  const PreconditionerTypes &...preconditioners)
{
  LogStream::Prefix prefix("MPGMRES");

  if (additional_data.use_truncated_mpgmres_strategy)
    SolverMPGMRES<VectorType>::solve_internal(
      A, x, b, IndexingStrategy::truncated_mpgmres, preconditioners...);
  else
    SolverMPGMRES<VectorType>::solve_internal(
      A, x, b, IndexingStrategy::full_mpgmres, preconditioners...);
}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
template <typename MatrixType, typename... PreconditionerTypes>
void SolverMPGMRES<VectorType>::solve_internal(
  const MatrixType       &A,
  VectorType             &x,
  const VectorType       &b,
  const IndexingStrategy &indexing_strategy,
  const PreconditionerTypes &...preconditioners)
{
  constexpr std::size_t n_preconditioners = sizeof...(PreconditionerTypes);

  // A lambda for applying the nth preconditioner to a vector src storing
  // the result in dst:

  const auto apply_nth_preconditioner = [&](unsigned int n,
                                            auto        &dst,
                                            const auto  &src) {
    // We cycle through all preconditioners and call the nth one:
    std::size_t i = 0;

    [[maybe_unused]] bool preconditioner_called = false;

    const auto call_matching_preconditioner = [&](const auto &preconditioner) {
      if (i++ == n)
        {
          Assert(!preconditioner_called, dealii::ExcInternalError());
          preconditioner_called = true;
          preconditioner.vmult(dst, src);
        }
    };

    // https://en.cppreference.com/w/cpp/language/fold
    (call_matching_preconditioner(preconditioners), ...);
    Assert(preconditioner_called, dealii::ExcInternalError());
  };

  std::size_t current_index = 0;

  // A lambda that cycles through all preconditioners in sequence while
  // applying exactly one preconditioner with each function invocation to
  // the vector src and storing the result in dst:

  const auto preconditioner_vmult = [&](auto &dst, const auto &src) {
    // We have no preconditioner that we could apply
    if (n_preconditioners == 0)
      dst = src;
    else
      {
        apply_nth_preconditioner(current_index, dst, src);
        current_index = (current_index + 1) % n_preconditioners;
      }
  };

  // Return the correct index for constructing the next vector in the
  // Krylov space sequence according to the chosen indexing strategy

  const auto previous_vector_index =
    [n_preconditioners, indexing_strategy](unsigned int i) -> unsigned int {
    // In the special case of no preconditioners we simply fall back to the
    // FGMRES indexing strategy.
    if (n_preconditioners == 0)
      {
        return i;
      }

    switch (indexing_strategy)
      {
        case IndexingStrategy::fgmres:
          // 0, 1, 2, 3, ...
          return i;
        case IndexingStrategy::full_mpgmres:
          // 0, 0, ..., 1, 1, ..., 2, 2, ..., 3, 3, ...
          return i / n_preconditioners;
        case IndexingStrategy::truncated_mpgmres:
          // 0, 0, ..., 1, 2, 3, ...
          return (1 + i >= n_preconditioners) ? (1 + i - n_preconditioners) : 0;
        default:
          DEAL_II_ASSERT_UNREACHABLE();
          return 0;
      }
  };

  SolverControl::State iteration_state = SolverControl::iterate;

  const unsigned int basis_size = additional_data.max_basis_size;

  // Generate an object where basis vectors are stored.
  typename internal::SolverGMRESImplementation::TmpVectors<VectorType> v(
    basis_size + 1, this->memory);
  typename internal::SolverGMRESImplementation::TmpVectors<VectorType> z(
    basis_size, this->memory);

  // number of the present iteration; this number is not reset to zero upon a
  // restart
  unsigned int accumulated_iterations = 0;

  // matrix used for the orthogonalization process later
  arnoldi_process.initialize(additional_data.orthogonalization_strategy,
                             basis_size,
                             false);

  // Iteration starts here
  double res = std::numeric_limits<double>::lowest();

  do
    {
      // Compute the residual. If 'x' is the zero vector, then we can bypass
      // the full computation. But 'x' is only likely to be the zero vector if
      // that's what the user provided as the starting guess, so it's only
      // worth checking for this in the first iteration. (Calling all_zero()
      // costs as much in memory transfer and communication as computing the
      // norm of a vector.)
      if (accumulated_iterations == 0 && x.all_zero())
        v(0, x) = b;
      else
        {
          A.vmult(v(0, x), x);
          v[0].sadd(-1., 1., b);
        }

      res             = arnoldi_process.orthonormalize_nth_vector(0, v);
      iteration_state = this->iteration_status(accumulated_iterations, res, x);
      if (iteration_state == SolverControl::success)
        break;

      unsigned int inner_iteration = 0;
      for (; (inner_iteration < basis_size &&
              iteration_state == SolverControl::iterate);
           ++inner_iteration)
        {
          preconditioner_vmult(z(inner_iteration, x),
                               v[previous_vector_index(inner_iteration)]);
          A.vmult(v(inner_iteration + 1, x), z[inner_iteration]);

          res =
            arnoldi_process.orthonormalize_nth_vector(inner_iteration + 1, v);

          // check convergence. note that the vector 'x' we pass to the
          // criterion is not the final solution we compute if we
          // decide to jump out of the iteration (we update 'x' again
          // right after the current loop)
          iteration_state =
            this->iteration_status(++accumulated_iterations, res, x);
        }

      // Solve triangular system with projected quantities and update solution
      // vector
      const Vector<double> &projected_solution =
        arnoldi_process.solve_projected_system(true);
      dealii::internal::SolverGMRESImplementation::add(
        x,
        inner_iteration,
        projected_solution,
        z,
        false,
        arnoldi_process.vector_ptrs);
    }
  while (iteration_state == SolverControl::iterate);

  // in case of failure: throw exception
  if (iteration_state != SolverControl::success)
    AssertThrow(false,
                SolverControl::NoConvergence(accumulated_iterations, res));
}



//----------------------------------------------------------------------//



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
SolverFGMRES<VectorType>::SolverFGMRES(SolverControl            &cn,
                                       VectorMemory<VectorType> &mem,
                                       const AdditionalData     &data)
  : SolverMPGMRES<VectorType>(
      cn,
      mem,
      typename SolverMPGMRES<VectorType>::AdditionalData{
        data.max_basis_size,
        data.orthogonalization_strategy,
        true})
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
SolverFGMRES<VectorType>::SolverFGMRES(SolverControl        &cn,
                                       const AdditionalData &data)
  : SolverMPGMRES<VectorType>(
      cn,
      typename SolverMPGMRES<VectorType>::AdditionalData{
        data.max_basis_size,
        data.orthogonalization_strategy,
        true})
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
template <typename MatrixType, typename... PreconditionerTypes>
DEAL_II_CXX20_REQUIRES(
  (concepts::is_linear_operator_on<MatrixType, VectorType> &&
   (concepts::is_linear_operator_on<PreconditionerTypes, VectorType> && ...)))
void SolverFGMRES<VectorType>::solve(
  const MatrixType &A,
  VectorType       &x,
  const VectorType &b,
  const PreconditionerTypes &...preconditioners)
{
  LogStream::Prefix prefix("FGMRES");
  SolverMPGMRES<VectorType>::solve_internal(
    A, x, b, SolverFGMRES::IndexingStrategy::fgmres, preconditioners...);
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

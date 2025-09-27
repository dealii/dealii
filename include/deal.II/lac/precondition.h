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

#ifndef dealii_precondition_h
#define dealii_precondition_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_space.h>
#include <deal.II/base/mutex.h>
#include <deal.II/base/observer_pointer.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_vector_base.h>
#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/identity_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/vector_memory.h>

#include <Kokkos_Core.hpp>

#include <limits>

DEAL_II_NAMESPACE_OPEN

// forward declarations
#ifndef DOXYGEN
template <typename number>
class Vector;
template <typename number>
class SparseMatrix;
namespace LinearAlgebra
{
  namespace distributed
  {
    template <typename, typename>
    class Vector;
    template <typename, typename>
    class BlockVector;
  } // namespace distributed
} // namespace LinearAlgebra
#endif


/**
 * @addtogroup Preconditioners
 * @{
 */

namespace internal
{
  /**
   * An enum to define the available types of eigenvalue estimation
   * algorithms.
   */
  enum class EigenvalueAlgorithm
  {
    /**
     * This option runs the conjugate gradient solver and computes an
     * eigenvalue estimation from the underlying Lanczos space. This only
     * works for symmetric positive definite matrices.
     */
    lanczos,
    /**
     * This option runs a power iteration to estimate the largest
     * eigenvalue. This algorithm also works for non-symmetric matrices,
     * but typically gives less accurate estimates than the option 'lanczos'
     * for the same number of iterations because it does not take the
     * relation between vectors in the iterations into account
     * (roughly speaking the off-diagonal entries in the
     * tri-diagonal matrix of the Lanczos iteration).
     */
    power_iteration
  };

  /**
   * A struct that contains information about the eigenvalue estimation
   * performed by the PreconditionChebyshev class.
   */
  struct EigenvalueInformation
  {
    /**
     * Estimate for the smallest eigenvalue.
     */
    double min_eigenvalue_estimate;
    /**
     * Estimate for the largest eigenvalue.
     */
    double max_eigenvalue_estimate;
    /**
     * Number of CG iterations performed or 0.
     */
    unsigned int cg_iterations;
    /**
     * The degree of the Chebyshev polynomial (either as set using
     * AdditionalData::degree or estimated as described there).
     */
    unsigned int degree;
    /**
     * Constructor initializing with invalid values.
     */
    EigenvalueInformation()
      : min_eigenvalue_estimate{std::numeric_limits<double>::max()}
      , max_eigenvalue_estimate{std::numeric_limits<double>::lowest()}
      , cg_iterations{0}
      , degree{0}
    {}
  };

  /**
   * Base class for PreconditionRelaxation::AdditionalData and
   * PreconditionChebyshev::AdditionalData with relevant information
   * for estimating the eigenvalues.
   */
  template <typename PreconditionerType>
  struct EigenvalueAlgorithmAdditionalData
  {
    /**
     * Constructor.
     */
    EigenvalueAlgorithmAdditionalData(
      const double              smoothing_range,
      const unsigned int        eig_cg_n_iterations,
      const double              eig_cg_residual,
      const double              max_eigenvalue,
      const EigenvalueAlgorithm eigenvalue_algorithm);

    /**
     * Copy assignment operator.
     */
    EigenvalueAlgorithmAdditionalData<PreconditionerType> &
    operator=(
      const EigenvalueAlgorithmAdditionalData<PreconditionerType> &other_data);

    /**
     * This sets the range between the largest eigenvalue in the matrix and
     * the smallest eigenvalue to be treated. If the parameter is set to a
     * number less than 1, an estimate for the largest and for the smallest
     * eigenvalue will be calculated internally. For a smoothing range larger
     * than one, the preconditioner will act in the interval
     * $[\lambda_\mathrm{max}/ \tt{smoothing\_range}, \lambda_\mathrm{max}]$,
     * where $\lambda_\mathrm{max}$ is an estimate of the maximum eigenvalue
     * of the matrix. A choice of <tt>smoothing_range</tt> between 5 and 20 is
     * useful in case the preconditioner is used as a smoother in multigrid.
     */
    double smoothing_range;

    /**
     * Maximum number of CG iterations performed for finding the maximum
     * eigenvalue. If set to zero, no computations are performed. Instead, the
     * user must supply a largest eigenvalue via the variable
     * PreconditionRelaxation::AdditionalData::max_eigenvalue.
     */
    unsigned int eig_cg_n_iterations;

    /**
     * Tolerance for iterations performed for finding the maximum
     * eigenvalue by the eigenvalue algorithm (Lanczos or power iteration).
     */
    double eig_cg_residual;

    /**
     * Maximum eigenvalue to work with. Only in effect if @p
     * eig_cg_n_iterations is set to zero, otherwise this parameter is
     * ignored.
     */
    double max_eigenvalue;

    /**
     * Constraints to be used for the operator given. This variable is used to
     * zero out the correct entries when creating an initial guess.
     */
    dealii::AffineConstraints<double> constraints;

    /**
     * Stores the preconditioner object that the Chebyshev is wrapped around.
     */
    EigenvalueAlgorithm eigenvalue_algorithm;

    /**
     * Preconditioner.
     */
    std::shared_ptr<PreconditionerType> preconditioner;
  };
} // namespace internal


/**
 * This preconditioner implements the identity operator.
 * This class can be used you if you want to use a linear
 * solver without preconditioning, because all linear solvers
 * require the usage of a preconditioner.
 *
 * It can be used in the following way:
 *
 * @code
 * SolverControl           solver_control (1000, 1e-12);
 * SolverCG<>              cg (solver_control);
 * cg.solve (system_matrix, solution, system_rhs, PreconditionIdentity());
 * @endcode
 *
 * See the step-3 tutorial program for an example and additional explanations.
 *
 * Alternatively, the IdentityMatrix class can be used to precondition in this
 * way.
 */
class PreconditionIdentity : public EnableObserverPointer
{
public:
  /**
   * Declare type for container size.
   */
  using size_type = types::global_dof_index;

  /**
   * This function is only present to provide the interface of a
   * preconditioner to be handed to a smoother.  This does nothing.
   */
  struct AdditionalData
  {
    /**
     * Constructor.
     */
    AdditionalData() = default;
  };

  /**
   * Constructor, sets the domain and range sizes to their defaults.
   */
  PreconditionIdentity();

  /**
   * The matrix argument is ignored and here just for compatibility with more
   * complex preconditioners.
   */
  template <typename MatrixType>
  void
  initialize(const MatrixType     &matrix,
             const AdditionalData &additional_data = AdditionalData());

  /**
   * Apply preconditioner.
   */
  template <typename VectorType>
  void
  vmult(VectorType &, const VectorType &) const;

  /**
   * Apply transpose preconditioner. Since this is the identity, this function
   * is the same as vmult().
   */
  template <typename VectorType>
  void
  Tvmult(VectorType &, const VectorType &) const;

  /**
   * Apply preconditioner, adding to the previous value.
   */
  template <typename VectorType>
  void
  vmult_add(VectorType &, const VectorType &) const;

  /**
   * Apply transpose preconditioner, adding. Since this is the identity, this
   * function is the same as vmult_add().
   */
  template <typename VectorType>
  void
  Tvmult_add(VectorType &, const VectorType &) const;

  /**
   * This function is only present to provide the interface of a
   * preconditioner to be handed to a smoother.  This does nothing.
   */
  void
  clear();

  /**
   * Return the dimension of the codomain (or range) space. Note that the
   * matrix is of dimension $m \times n$.
   *
   * @note This function should only be called if the preconditioner has been
   * initialized.
   */
  size_type
  m() const;

  /**
   * Return the dimension of the domain space. Note that the matrix is of
   * dimension $m \times n$.
   *
   * @note This function should only be called if the preconditioner has been
   * initialized.
   */
  size_type
  n() const;

private:
  /**
   * The dimension of the range space.
   */
  size_type n_rows;

  /**
   * The dimension of the domain space.
   */
  size_type n_columns;
};



/**
 * Preconditioning with Richardson's method. This preconditioner just scales
 * the vector with a constant relaxation factor provided by the AdditionalData
 * object.
 *
 * In Krylov-space methods, this preconditioner should not have any effect.
 * Using SolverRichardson, the two relaxation parameters will be just
 * multiplied. Still, this class is useful in multigrid smoother objects
 * (MGSmootherRelaxation).
 */
class PreconditionRichardson : public EnableObserverPointer
{
public:
  /**
   * Declare type for container size.
   */
  using size_type = types::global_dof_index;

  /**
   * Parameters for Richardson preconditioner.
   */
  class AdditionalData
  {
  public:
    /**
     * Constructor. Block size must be given since there is no reasonable
     * default parameter.
     */
    AdditionalData(const double relaxation = 1.);

    /**
     * Relaxation parameter.
     */
    double relaxation;
  };

  /**
   * Constructor, sets the relaxation parameter, domain and range sizes to
   * their default.
   */
  PreconditionRichardson();

  /**
   * Change the relaxation parameter.
   */
  void
  initialize(const AdditionalData &parameters);

  /**
   * Change the relaxation parameter in a way consistent with other
   * preconditioners. The matrix argument is ignored and here just for
   * compatibility with more complex preconditioners.
   */
  template <typename MatrixType>
  void
  initialize(const MatrixType &matrix, const AdditionalData &parameters);

  /**
   * Apply preconditioner.
   */
  template <typename VectorType>
  void
  vmult(VectorType &, const VectorType &) const;

  /**
   * Apply transpose preconditioner. Since this is the identity, this function
   * is the same as vmult().
   */
  template <typename VectorType>
  void
  Tvmult(VectorType &, const VectorType &) const;
  /**
   * Apply preconditioner, adding to the previous value.
   */
  template <typename VectorType>
  void
  vmult_add(VectorType &, const VectorType &) const;

  /**
   * Apply transpose preconditioner, adding. Since this is the identity, this
   * function is the same as vmult_add().
   */
  template <typename VectorType>
  void
  Tvmult_add(VectorType &, const VectorType &) const;

  /**
   * This function is only present to provide the interface of a
   * preconditioner to be handed to a smoother.  This does nothing.
   */
  void
  clear()
  {}

  /**
   * Return the dimension of the codomain (or range) space. Note that the
   * matrix is of dimension $m \times n$.
   *
   * @note This function should only be called if the preconditioner has been
   * initialized.
   */
  size_type
  m() const;

  /**
   * Return the dimension of the domain space. Note that the matrix is of
   * dimension $m \times n$.
   *
   * @note This function should only be called if the preconditioner has been
   * initialized.
   */
  size_type
  n() const;

private:
  /**
   * The relaxation parameter multiplied with the vectors.
   */
  double relaxation;

  /**
   * The dimension of the range space.
   */
  size_type n_rows;

  /**
   * The dimension of the domain space.
   */
  size_type n_columns;
};



/**
 * Preconditioner using a matrix-builtin function.  This class forms a
 * preconditioner suitable for the LAC solver classes. Since many
 * preconditioning methods are based on matrix entries, these have to be
 * implemented as member functions of the underlying matrix implementation.
 * This class now is intended to allow easy access to these member functions
 * from LAC solver classes.
 *
 * It seems that all builtin preconditioners have a relaxation parameter, so
 * please use PreconditionRelaxation for these.
 *
 * You will usually not want to create a named object of this type, although
 * possible. The most common use is like this:
 * @code
 * SolverGMRES<SparseMatrix<double> Vector<double>> gmres(control,memory,500);
 *
 * gmres.solve(
 *   matrix, solution, right_hand_side,
 *   PreconditionUseMatrix<SparseMatrix<double>,Vector<double> >(
 *     matrix, &SparseMatrix<double>::template precondition_Jacobi<double>));
 * @endcode
 * This creates an unnamed object to be passed as the fourth parameter to the
 * solver function of the SolverGMRES class. It assumes that the SparseMatrix
 * class has a function <tt>precondition_Jacobi</tt> taking two vectors
 * (source and destination) as parameters (Actually, there is no function like
 * that, the existing function takes a third parameter, denoting the
 * relaxation parameter; this example is therefore only meant to illustrate
 * the general idea).
 *
 * Note that due to the default template parameters, the above example could
 * be written shorter as follows:
 * @code
 * ...
 * gmres.solve(
 *   matrix, solution, right_hand_side,
 *   PreconditionUseMatrix<>(
 *     matrix,&SparseMatrix<double>::template precondition_Jacobi<double>));
 * @endcode
 */
template <typename MatrixType = SparseMatrix<double>,
          typename VectorType = Vector<double>>
class PreconditionUseMatrix : public EnableObserverPointer
{
public:
  /**
   * Type of the preconditioning function of the matrix.
   */
  using function_ptr = void (MatrixType::*)(VectorType &,
                                            const VectorType &) const;

  /**
   * Constructor.  This constructor stores a reference to the matrix object
   * for later use and selects a preconditioning method, which must be a
   * member function of that matrix.
   */
  PreconditionUseMatrix(const MatrixType &M, const function_ptr method);

  /**
   * Execute preconditioning. Calls the function passed to the constructor of
   * this object with the two arguments given here.
   */
  void
  vmult(VectorType &dst, const VectorType &src) const;

private:
  /**
   * Pointer to the matrix in use.
   */
  const MatrixType &matrix;

  /**
   * Pointer to the preconditioning function.
   */
  const function_ptr precondition;
};



/**
 * Base class for other preconditioners. Here, only some common features
 * Jacobi, SOR and SSOR preconditioners are implemented. For preconditioning,
 * refer to derived classes.
 *
 * One iteration is given as:
 * @f[
 *  x^{n+1} = x^{n} + \alpha P^{-1} (b-Ax^n).
 * @f]
 *
 * The relaxation parameter $\alpha$ has to be in the range:
 * @f[
 *  0 < \alpha < \frac{2}{\lambda_{\max}(P^{-1}A)}.
 * @f]
 * Its theoretically optimal value is given by:
 * @f[
 *    \alpha := \frac{2}{\lambda_{\min}(P^{-1}A)+\lambda_{\max}(P^{-1}A)}.
 * @f]
 *
 * For details on the algorithm, see @cite Varga2009.
 *
 * The relaxation parameter can be set manually or can be automatically
 * determined, for which we use the theoretically optimal value.
 * For this purpose, the user needs to set the relaxation parameter
 * to zero. Internally, the minimum and maximum eigenvalues of the
 * preconditioned system are estimated by an eigenvalue algorithm, and the
 * resulting estimate is multiplied by the 1.2 for safety reasons. For more
 * details on the underlying algorithms, see PreconditionChebyshev.
 */
template <typename MatrixType         = SparseMatrix<double>,
          typename PreconditionerType = IdentityMatrix>
class PreconditionRelaxation : public EnableObserverPointer
{
public:
  /**
   * Declare type for container size.
   */
  using size_type = types::global_dof_index;

  /**
   * Class for parameters.
   */
  class AdditionalData
    : public internal::EigenvalueAlgorithmAdditionalData<PreconditionerType>
  {
  public:
    using EigenvalueAlgorithm = internal::EigenvalueAlgorithm;

    /**
     * Constructor.
     */
    AdditionalData(const double              relaxation          = 1.,
                   const unsigned int        n_iterations        = 1,
                   const double              smoothing_range     = 0.,
                   const unsigned int        eig_cg_n_iterations = 8,
                   const double              eig_cg_residual     = 1e-2,
                   const double              max_eigenvalue      = 1,
                   const EigenvalueAlgorithm eigenvalue_algorithm =
                     EigenvalueAlgorithm::lanczos);

    /**
     * Relaxation parameter.
     */
    double relaxation;

    /**
     * Number of smoothing steps to be performed in an
     * invocation of vmult() or step().
     */
    unsigned int n_iterations;
  };

  /**
   * Initialize matrix and relaxation parameter. The matrix is just stored in
   * the preconditioner object. The relaxation parameter should be larger than
   * zero and smaller than 2 for numerical reasons. It defaults to 1.
   */
  void
  initialize(const MatrixType     &A,
             const AdditionalData &parameters = AdditionalData());

  /**
   * Release the matrix and reset its pointer.
   */
  void
  clear();

  /**
   * Return the dimension of the codomain (or range) space. Note that the
   * matrix is of dimension $m \times n$.
   */
  size_type
  m() const;

  /**
   * Return the dimension of the domain space. Note that the matrix is of
   * dimension $m \times n$.
   */
  size_type
  n() const;

  /**
   * Apply preconditioner.
   */
  template <typename VectorType>
  void
  vmult(VectorType &, const VectorType &) const;

  /**
   * Apply transpose preconditioner. Since this is a symmetric preconditioner,
   * this function is the same as vmult().
   */
  template <typename VectorType>
  void
  Tvmult(VectorType &, const VectorType &) const;

  /**
   * Perform one step of the preconditioned Richardson iteration
   */
  template <typename VectorType>
  void
  step(VectorType &x, const VectorType &rhs) const;

  /**
   * Perform one transposed step of the preconditioned Richardson iteration.
   */
  template <typename VectorType>
  void
  Tstep(VectorType &x, const VectorType &rhs) const;

  using EigenvalueInformation = internal::EigenvalueInformation;

  /**
   * Estimate eigenvalues and set relaxation parameter.
   */
  template <typename VectorType>
  EigenvalueInformation
  estimate_eigenvalues(const VectorType &src) const;

  /**
   * Return the relaxation parameter. This function also allows to return
   * the parameter in case it was internally determined by running an eigenvalue
   * algorithm.
   */
  double
  get_relaxation() const;

protected:
  /**
   * Pointer to the matrix object.
   */
  ObserverPointer<const MatrixType, PreconditionRelaxation<MatrixType>> A;

  /**
   * Stores the additional data passed to the initialize function, obtained
   * through a copy operation.
   */
  AdditionalData data;

  /**
   * Preconditioner.
   */
  std::shared_ptr<PreconditionerType> preconditioner;

  /**
   * Stores whether the preconditioner has been set up and eigenvalues have
   * been computed.
   */
  bool eigenvalues_are_initialized;
};



#ifndef DOXYGEN

namespace internal
{
  // a helper type-trait that leverage SFINAE to figure out if MatrixType has
  // ... MatrixType::vmult(VectorType &, const VectorType&,
  // std::function<...>, std::function<...>) const
  template <typename MatrixType, typename VectorType>
  using vmult_functions_t = decltype(std::declval<const MatrixType>().vmult(
    std::declval<VectorType &>(),
    std::declval<const VectorType &>(),
    std::declval<
      const std::function<void(const unsigned int, const unsigned int)> &>(),
    std::declval<
      const std::function<void(const unsigned int, const unsigned int)> &>()));

  template <typename MatrixType,
            typename VectorType,
            typename PreconditionerType>
  constexpr bool has_vmult_with_std_functions =
    is_supported_operation<vmult_functions_t, MatrixType, VectorType> &&
    std::is_same_v<PreconditionerType, dealii::DiagonalMatrix<VectorType>> &&
    (std::is_same_v<VectorType,
                    dealii::Vector<typename VectorType::value_type>> ||
     std::is_same_v<
       VectorType,
       LinearAlgebra::distributed::Vector<typename VectorType::value_type,
                                          MemorySpace::Host>>);


  template <typename MatrixType, typename VectorType>
  constexpr bool has_vmult_with_std_functions_for_precondition =
    is_supported_operation<vmult_functions_t, MatrixType, VectorType>;

  namespace PreconditionRelaxation
  {
    template <typename T, typename VectorType>
    using Tvmult_t = decltype(std::declval<const T>().Tvmult(
      std::declval<VectorType &>(),
      std::declval<const VectorType &>()));

    template <typename T, typename VectorType>
    constexpr bool has_Tvmult = is_supported_operation<Tvmult_t, T, VectorType>;

    template <typename T, typename VectorType>
    using step_t = decltype(std::declval<const T>().step(
      std::declval<VectorType &>(),
      std::declval<const VectorType &>()));

    template <typename T, typename VectorType>
    constexpr bool has_step = is_supported_operation<step_t, T, VectorType>;

    template <typename T, typename VectorType>
    using step_omega_t =
      decltype(std::declval<const T>().step(std::declval<VectorType &>(),
                                            std::declval<const VectorType &>(),
                                            std::declval<const double>()));

    template <typename T, typename VectorType>
    constexpr bool has_step_omega =
      is_supported_operation<step_omega_t, T, VectorType>;

    template <typename T, typename VectorType>
    using Tstep_t = decltype(std::declval<const T>().Tstep(
      std::declval<VectorType &>(),
      std::declval<const VectorType &>()));

    template <typename T, typename VectorType>
    constexpr bool has_Tstep = is_supported_operation<Tstep_t, T, VectorType>;

    template <typename T, typename VectorType>
    using Tstep_omega_t =
      decltype(std::declval<const T>().Tstep(std::declval<VectorType &>(),
                                             std::declval<const VectorType &>(),
                                             std::declval<const double>()));

    template <typename T, typename VectorType>
    constexpr bool has_Tstep_omega =
      is_supported_operation<Tstep_omega_t, T, VectorType>;

    template <typename T, typename VectorType>
    using jacobi_step_t = decltype(std::declval<const T>().Jacobi_step(
      std::declval<VectorType &>(),
      std::declval<const VectorType &>(),
      std::declval<const double>()));

    template <typename T, typename VectorType>
    constexpr bool has_jacobi_step =
      is_supported_operation<jacobi_step_t, T, VectorType>;

    template <typename T, typename VectorType>
    using SOR_step_t = decltype(std::declval<const T>().SOR_step(
      std::declval<VectorType &>(),
      std::declval<const VectorType &>(),
      std::declval<const double>()));

    template <typename T, typename VectorType>
    constexpr bool has_SOR_step =
      is_supported_operation<SOR_step_t, T, VectorType>;

    template <typename T, typename VectorType>
    using SSOR_step_t = decltype(std::declval<const T>().SSOR_step(
      std::declval<VectorType &>(),
      std::declval<const VectorType &>(),
      std::declval<const double>()));

    template <typename T, typename VectorType>
    constexpr bool has_SSOR_step =
      is_supported_operation<SSOR_step_t, T, VectorType>;

    template <typename MatrixType>
    class PreconditionJacobiImpl
    {
    public:
      PreconditionJacobiImpl(const MatrixType &A, const double relaxation)
        : A(&A)
        , relaxation(relaxation)
      {}

      template <typename VectorType>
      void
      vmult(VectorType &dst, const VectorType &src) const
      {
        this->A->precondition_Jacobi(dst, src, this->relaxation);
      }

      template <typename VectorType>
      void
      Tvmult(VectorType &dst, const VectorType &src) const
      {
        // call vmult, since preconditioner is symmetrical
        this->vmult(dst, src);
      }

      template <typename VectorType,
                std::enable_if_t<has_jacobi_step<MatrixType, VectorType>,
                                 MatrixType> * = nullptr>
      void
      step(VectorType &dst, const VectorType &src) const
      {
        this->A->Jacobi_step(dst, src, this->relaxation);
      }

      template <typename VectorType,
                std::enable_if_t<!has_jacobi_step<MatrixType, VectorType>,
                                 MatrixType> * = nullptr>
      void
      step(VectorType &, const VectorType &) const
      {
        AssertThrow(false,
                    ExcMessage(
                      "Matrix A does not provide a Jacobi_step() function!"));
      }

      template <typename VectorType>
      void
      Tstep(VectorType &dst, const VectorType &src) const
      {
        // call step, since preconditioner is symmetrical
        this->step(dst, src);
      }

    private:
      const ObserverPointer<const MatrixType> A;
      const double                            relaxation;
    };

    template <typename MatrixType>
    class PreconditionSORImpl
    {
    public:
      PreconditionSORImpl(const MatrixType &A, const double relaxation)
        : A(&A)
        , relaxation(relaxation)
      {}

      template <typename VectorType>
      void
      vmult(VectorType &dst, const VectorType &src) const
      {
        this->A->precondition_SOR(dst, src, this->relaxation);
      }

      template <typename VectorType>
      void
      Tvmult(VectorType &dst, const VectorType &src) const
      {
        this->A->precondition_TSOR(dst, src, this->relaxation);
      }

      template <typename VectorType,
                std::enable_if_t<has_SOR_step<MatrixType, VectorType>,
                                 MatrixType> * = nullptr>
      void
      step(VectorType &dst, const VectorType &src) const
      {
        this->A->SOR_step(dst, src, this->relaxation);
      }

      template <typename VectorType,
                std::enable_if_t<!has_SOR_step<MatrixType, VectorType>,
                                 MatrixType> * = nullptr>
      void
      step(VectorType &, const VectorType &) const
      {
        AssertThrow(false,
                    ExcMessage(
                      "Matrix A does not provide a SOR_step() function!"));
      }

      template <typename VectorType,
                std::enable_if_t<has_SOR_step<MatrixType, VectorType>,
                                 MatrixType> * = nullptr>
      void
      Tstep(VectorType &dst, const VectorType &src) const
      {
        this->A->TSOR_step(dst, src, this->relaxation);
      }

      template <typename VectorType,
                std::enable_if_t<!has_SOR_step<MatrixType, VectorType>,
                                 MatrixType> * = nullptr>
      void
      Tstep(VectorType &, const VectorType &) const
      {
        AssertThrow(false,
                    ExcMessage(
                      "Matrix A does not provide a TSOR_step() function!"));
      }

    private:
      const ObserverPointer<const MatrixType> A;
      const double                            relaxation;
    };

    template <typename MatrixType>
    class PreconditionSSORImpl
    {
    public:
      using size_type = typename MatrixType::size_type;

      PreconditionSSORImpl(const MatrixType &A, const double relaxation)
        : A(&A)
        , relaxation(relaxation)
      {
        // in case we have a SparseMatrix class, we can extract information
        // about the diagonal.
        const SparseMatrix<typename MatrixType::value_type> *mat =
          dynamic_cast<const SparseMatrix<typename MatrixType::value_type> *>(
            &*this->A);

        // calculate the positions first after the diagonal.
        if (mat != nullptr)
          {
            const size_type n = this->A->n();
            pos_right_of_diagonal.resize(n, static_cast<std::size_t>(-1));
            for (size_type row = 0; row < n; ++row)
              {
                // find the first element in this line which is on the right of
                // the diagonal.  we need to precondition with the elements on
                // the left only. note: the first entry in each line denotes the
                // diagonal element, which we need not check.
                typename SparseMatrix<
                  typename MatrixType::value_type>::const_iterator it =
                  mat->begin(row) + 1;
                for (; it < mat->end(row); ++it)
                  if (it->column() > row)
                    break;
                pos_right_of_diagonal[row] = it - mat->begin();
              }
          }
      }

      template <typename VectorType>
      void
      vmult(VectorType &dst, const VectorType &src) const
      {
        this->A->precondition_SSOR(dst,
                                   src,
                                   this->relaxation,
                                   pos_right_of_diagonal);
      }

      template <typename VectorType>
      void
      Tvmult(VectorType &dst, const VectorType &src) const
      {
        this->A->precondition_SSOR(dst,
                                   src,
                                   this->relaxation,
                                   pos_right_of_diagonal);
      }

      template <typename VectorType,
                std::enable_if_t<has_SSOR_step<MatrixType, VectorType>,
                                 MatrixType> * = nullptr>
      void
      step(VectorType &dst, const VectorType &src) const
      {
        this->A->SSOR_step(dst, src, this->relaxation);
      }

      template <typename VectorType,
                std::enable_if_t<!has_SSOR_step<MatrixType, VectorType>,
                                 MatrixType> * = nullptr>
      void
      step(VectorType &, const VectorType &) const
      {
        AssertThrow(false,
                    ExcMessage(
                      "Matrix A does not provide a SSOR_step() function!"));
      }

      template <typename VectorType>
      void
      Tstep(VectorType &dst, const VectorType &src) const
      {
        // call step, since preconditioner is symmetrical
        this->step(dst, src);
      }

    private:
      const ObserverPointer<const MatrixType> A;
      const double                            relaxation;

      /**
       * An array that stores for each matrix row where the first position after
       * the diagonal is located.
       */
      std::vector<std::size_t> pos_right_of_diagonal;
    };

    template <typename MatrixType>
    class PreconditionPSORImpl
    {
    public:
      using size_type = typename MatrixType::size_type;

      PreconditionPSORImpl(const MatrixType             &A,
                           const double                  relaxation,
                           const std::vector<size_type> &permutation,
                           const std::vector<size_type> &inverse_permutation)
        : A(&A)
        , relaxation(relaxation)
        , permutation(permutation)
        , inverse_permutation(inverse_permutation)
      {}

      template <typename VectorType>
      void
      vmult(VectorType &dst, const VectorType &src) const
      {
        dst = src;
        this->A->PSOR(dst, permutation, inverse_permutation, this->relaxation);
      }

      template <typename VectorType>
      void
      Tvmult(VectorType &dst, const VectorType &src) const
      {
        dst = src;
        this->A->TPSOR(dst, permutation, inverse_permutation, this->relaxation);
      }

    private:
      const ObserverPointer<const MatrixType> A;
      const double                            relaxation;

      const std::vector<size_type> &permutation;
      const std::vector<size_type> &inverse_permutation;
    };

    template <typename MatrixType,
              typename PreconditionerType,
              typename VectorType,
              std::enable_if_t<has_step_omega<PreconditionerType, VectorType>,
                               PreconditionerType> * = nullptr>
    void
    step(const MatrixType &,
         const PreconditionerType &preconditioner,
         VectorType               &dst,
         const VectorType         &src,
         const double              relaxation,
         VectorType &,
         VectorType &)
    {
      preconditioner.step(dst, src, relaxation);
    }

    template <
      typename MatrixType,
      typename PreconditionerType,
      typename VectorType,
      std::enable_if_t<!has_step_omega<PreconditionerType, VectorType> &&
                         has_step<PreconditionerType, VectorType>,
                       PreconditionerType> * = nullptr>
    void
    step(const MatrixType &,
         const PreconditionerType &preconditioner,
         VectorType               &dst,
         const VectorType         &src,
         const double              relaxation,
         VectorType &,
         VectorType &)
    {
      Assert(relaxation == 1.0, ExcInternalError());

      (void)relaxation;

      preconditioner.step(dst, src);
    }

    template <
      typename MatrixType,
      typename PreconditionerType,
      typename VectorType,
      std::enable_if_t<!has_step_omega<PreconditionerType, VectorType> &&
                         !has_step<PreconditionerType, VectorType>,
                       PreconditionerType> * = nullptr>
    void
    step(const MatrixType         &A,
         const PreconditionerType &preconditioner,
         VectorType               &dst,
         const VectorType         &src,
         const double              relaxation,
         VectorType               &residual,
         VectorType               &tmp)
    {
      residual.reinit(dst, true);
      tmp.reinit(dst, true);

      A.vmult(residual, dst);
      residual.sadd(-1.0, 1.0, src);

      preconditioner.vmult(tmp, residual);
      dst.add(relaxation, tmp);
    }

    template <typename MatrixType,
              typename PreconditionerType,
              typename VectorType,
              std::enable_if_t<has_Tstep_omega<PreconditionerType, VectorType>,
                               PreconditionerType> * = nullptr>
    void
    Tstep(const MatrixType &,
          const PreconditionerType &preconditioner,
          VectorType               &dst,
          const VectorType         &src,
          const double              relaxation,
          VectorType &,
          VectorType &)
    {
      preconditioner.Tstep(dst, src, relaxation);
    }

    template <
      typename MatrixType,
      typename PreconditionerType,
      typename VectorType,
      std::enable_if_t<!has_Tstep_omega<PreconditionerType, VectorType> &&
                         has_Tstep<PreconditionerType, VectorType>,
                       PreconditionerType> * = nullptr>
    void
    Tstep(const MatrixType &,
          const PreconditionerType &preconditioner,
          VectorType               &dst,
          const VectorType         &src,
          const double              relaxation,
          VectorType &,
          VectorType &)
    {
      Assert(relaxation == 1.0, ExcInternalError());

      (void)relaxation;

      preconditioner.Tstep(dst, src);
    }

    template <typename MatrixType,
              typename VectorType,
              std::enable_if_t<has_Tvmult<MatrixType, VectorType>, MatrixType>
                * = nullptr>
    void
    Tvmult(const MatrixType &A, VectorType &dst, const VectorType &src)
    {
      A.Tvmult(dst, src);
    }

    template <typename MatrixType,
              typename VectorType,
              std::enable_if_t<!has_Tvmult<MatrixType, VectorType>, MatrixType>
                * = nullptr>
    void
    Tvmult(const MatrixType &, VectorType &, const VectorType &)
    {
      AssertThrow(false,
                  ExcMessage("Matrix A does not provide a Tvmult() function!"));
    }

    template <
      typename MatrixType,
      typename PreconditionerType,
      typename VectorType,
      std::enable_if_t<!has_Tstep_omega<PreconditionerType, VectorType> &&
                         !has_Tstep<PreconditionerType, VectorType>,
                       PreconditionerType> * = nullptr>
    void
    Tstep(const MatrixType         &A,
          const PreconditionerType &preconditioner,
          VectorType               &dst,
          const VectorType         &src,
          const double              relaxation,
          VectorType               &residual,
          VectorType               &tmp)
    {
      residual.reinit(dst, true);
      tmp.reinit(dst, true);

      Tvmult(A, residual, dst);
      residual.sadd(-1.0, 1.0, src);

      Tvmult(preconditioner, tmp, residual);
      dst.add(relaxation, tmp);
    }

    // 0) general implementation
    template <typename MatrixType,
              typename PreconditionerType,
              typename VectorType,
              std::enable_if_t<!has_vmult_with_std_functions_for_precondition<
                                 PreconditionerType,
                                 VectorType>,
                               int> * = nullptr>
    void
    step_operations(const MatrixType         &A,
                    const PreconditionerType &preconditioner,
                    VectorType               &dst,
                    const VectorType         &src,
                    const double              relaxation,
                    VectorType               &tmp1,
                    VectorType               &tmp2,
                    const unsigned int        i,
                    const bool                transposed)
    {
      if (i == 0)
        {
          if (transposed)
            Tvmult(preconditioner, dst, src);
          else
            preconditioner.vmult(dst, src);

          if (relaxation != 1.0)
            dst *= relaxation;
        }
      else
        {
          if (transposed)
            Tstep(A, preconditioner, dst, src, relaxation, tmp1, tmp2);
          else
            step(A, preconditioner, dst, src, relaxation, tmp1, tmp2);
        }
    }

    // 1) specialized implementation with a preconditioner that accepts
    // ranges
    template <
      typename MatrixType,
      typename PreconditionerType,
      typename VectorType,
      std::enable_if_t<
        has_vmult_with_std_functions_for_precondition<PreconditionerType,
                                                      VectorType> &&
          !has_vmult_with_std_functions_for_precondition<MatrixType,
                                                         VectorType>,
        int> * = nullptr>
    void
    step_operations(const MatrixType         &A,
                    const PreconditionerType &preconditioner,
                    VectorType               &dst,
                    const VectorType         &src,
                    const double              relaxation,
                    VectorType               &tmp,
                    VectorType &,
                    const unsigned int i,
                    const bool         transposed)
    {
      (void)transposed;
      using Number          = typename VectorType::value_type;
      Number       *dst_ptr = dst.begin();
      const Number *src_ptr = src.begin();

      if (i == 0)
        {
          preconditioner.vmult(
            dst,
            src,
            [&](const unsigned int start_range, const unsigned int end_range) {
              // zero 'dst' before running the vmult operation
              if (end_range > start_range)
                std::memset(dst_ptr + start_range,
                            0,
                            sizeof(Number) * (end_range - start_range));
            },
            [&](const unsigned int start_range, const unsigned int end_range) {
              if (relaxation == 1.0)
                return; // nothing to do

              DEAL_II_OPENMP_SIMD_PRAGMA
              for (std::size_t i = start_range; i < end_range; ++i)
                dst_ptr[i] *= relaxation;
            });
        }
      else
        {
          tmp.reinit(src, true);

          Assert(transposed == false, ExcNotImplemented());

          A.vmult(tmp, dst);

          preconditioner.vmult(
            dst,
            tmp,
            [&](const unsigned int start_range, const unsigned int end_range) {
              const auto tmp_ptr = tmp.begin();

              if (relaxation == 1.0)
                {
                  DEAL_II_OPENMP_SIMD_PRAGMA
                  for (std::size_t i = start_range; i < end_range; ++i)
                    tmp_ptr[i] = src_ptr[i] - tmp_ptr[i];
                }
              else
                {
                  // note: we scale the residual here to be able to add into
                  // the dst vector, which contains the solution from the last
                  // iteration
                  DEAL_II_OPENMP_SIMD_PRAGMA
                  for (std::size_t i = start_range; i < end_range; ++i)
                    tmp_ptr[i] = relaxation * (src_ptr[i] - tmp_ptr[i]);
                }
            },
            [&](const unsigned int, const unsigned int) {
              // nothing to do, since scaling by the relaxation factor
              // has been done in the pre operation
            });
        }
    }

    // 2) specialized implementation with a preconditioner and a matrix that
    // accepts ranges
    template <
      typename MatrixType,
      typename PreconditionerType,
      typename VectorType,
      std::enable_if_t<
        has_vmult_with_std_functions_for_precondition<PreconditionerType,
                                                      VectorType> &&
          has_vmult_with_std_functions_for_precondition<MatrixType, VectorType>,
        int> * = nullptr>
    void
    step_operations(const MatrixType         &A,
                    const PreconditionerType &preconditioner,
                    VectorType               &dst,
                    const VectorType         &src,
                    const double              relaxation,
                    VectorType               &tmp,
                    VectorType &,
                    const unsigned int i,
                    const bool         transposed)
    {
      (void)transposed;
      using Number = typename VectorType::value_type;

      Number       *dst_ptr = dst.begin();
      const Number *src_ptr = src.begin();

      if (i == 0)
        {
          preconditioner.vmult(
            dst,
            src,
            [&](const unsigned int start_range, const unsigned int end_range) {
              // zero 'dst' before running the vmult operation
              if (end_range > start_range)
                std::memset(dst_ptr + start_range,
                            0,
                            sizeof(Number) * (end_range - start_range));
            },
            [&](const unsigned int start_range, const unsigned int end_range) {
              if (relaxation == 1.0)
                return; // nothing to do

              DEAL_II_OPENMP_SIMD_PRAGMA
              for (std::size_t i = start_range; i < end_range; ++i)
                dst_ptr[i] *= relaxation;
            });
        }
      else
        {
          tmp.reinit(src, true);
          const auto tmp_ptr = tmp.begin();

          Assert(transposed == false, ExcNotImplemented());

          A.vmult(
            tmp,
            dst,
            [&](const unsigned int start_range, const unsigned int end_range) {
              // zero 'tmp' before running the vmult
              // operation
              if (end_range > start_range)
                std::memset(tmp_ptr + start_range,
                            0,
                            sizeof(Number) * (end_range - start_range));
            },
            [&](const unsigned int start_range, const unsigned int end_range) {
              if (relaxation == 1.0)
                {
                  DEAL_II_OPENMP_SIMD_PRAGMA
                  for (std::size_t i = start_range; i < end_range; ++i)
                    tmp_ptr[i] = src_ptr[i] - tmp_ptr[i];
                }
              else
                {
                  // note: we scale the residual here to be able to add into
                  // the dst vector, which contains the solution from the last
                  // iteration
                  DEAL_II_OPENMP_SIMD_PRAGMA
                  for (std::size_t i = start_range; i < end_range; ++i)
                    tmp_ptr[i] = relaxation * (src_ptr[i] - tmp_ptr[i]);
                }
            });

          preconditioner.vmult(dst, tmp, [](const auto, const auto) {
            // note: `dst` vector does not have to be zeroed out
            // since we add the result into it
          });
        }
    }

    // 3) specialized implementation for inverse-diagonal preconditioner
    template <
      typename MatrixType,
      typename VectorType,
      std::enable_if_t<
        !IsBlockVector<VectorType>::value &&
          !std::is_same_v<
            VectorType,
            LinearAlgebra::distributed::Vector<typename VectorType::value_type,
                                               MemorySpace::Default>> &&
          !has_vmult_with_std_functions<MatrixType,
                                        VectorType,
                                        dealii::DiagonalMatrix<VectorType>>,
        VectorType> * = nullptr>
    void
    step_operations(const MatrixType                         &A,
                    const dealii::DiagonalMatrix<VectorType> &preconditioner,
                    VectorType                               &dst,
                    const VectorType                         &src,
                    const double                              relaxation,
                    VectorType                               &tmp,
                    VectorType &,
                    const unsigned int i,
                    const bool         transposed)
    {
      using Number = typename VectorType::value_type;

      if (i == 0)
        {
          Number       *dst_ptr  = dst.begin();
          const Number *src_ptr  = src.begin();
          const Number *diag_ptr = preconditioner.get_vector().begin();

          if (relaxation == 1.0)
            {
              DEAL_II_OPENMP_SIMD_PRAGMA
              for (unsigned int i = 0; i < dst.locally_owned_size(); ++i)
                dst_ptr[i] = src_ptr[i] * diag_ptr[i];
            }
          else
            {
              DEAL_II_OPENMP_SIMD_PRAGMA
              for (unsigned int i = 0; i < dst.locally_owned_size(); ++i)
                dst_ptr[i] = relaxation * src_ptr[i] * diag_ptr[i];
            }
        }
      else
        {
          tmp.reinit(src, true);

          if (transposed)
            Tvmult(A, tmp, dst);
          else
            A.vmult(tmp, dst);

          Number       *dst_ptr  = dst.begin();
          const Number *src_ptr  = src.begin();
          const Number *tmp_ptr  = tmp.begin();
          const Number *diag_ptr = preconditioner.get_vector().begin();

          if (relaxation == 1.0)
            {
              DEAL_II_OPENMP_SIMD_PRAGMA
              for (unsigned int i = 0; i < dst.locally_owned_size(); ++i)
                dst_ptr[i] += (src_ptr[i] - tmp_ptr[i]) * diag_ptr[i];
            }
          else
            {
              DEAL_II_OPENMP_SIMD_PRAGMA
              for (unsigned int i = 0; i < dst.locally_owned_size(); ++i)
                dst_ptr[i] +=
                  relaxation * (src_ptr[i] - tmp_ptr[i]) * diag_ptr[i];
            }
        }
    }

    // 4) specialized implementation for inverse-diagonal preconditioner and
    // matrix that accepts ranges
    template <typename MatrixType,
              typename VectorType,
              std::enable_if_t<!IsBlockVector<VectorType>::value &&
                                 has_vmult_with_std_functions<
                                   MatrixType,
                                   VectorType,
                                   dealii::DiagonalMatrix<VectorType>>,
                               VectorType> * = nullptr>
    void
    step_operations(const MatrixType                         &A,
                    const dealii::DiagonalMatrix<VectorType> &preconditioner,
                    VectorType                               &dst,
                    const VectorType                         &src,
                    const double                              relaxation,
                    VectorType                               &tmp,
                    VectorType &,
                    const unsigned int i,
                    const bool         transposed)
    {
      (void)transposed;
      using Number = typename VectorType::value_type;

      if (i == 0)
        {
          Number       *dst_ptr  = dst.begin();
          const Number *src_ptr  = src.begin();
          const Number *diag_ptr = preconditioner.get_vector().begin();

          if (relaxation == 1.0)
            {
              DEAL_II_OPENMP_SIMD_PRAGMA
              for (unsigned int i = 0; i < dst.locally_owned_size(); ++i)
                dst_ptr[i] = src_ptr[i] * diag_ptr[i];
            }
          else
            {
              DEAL_II_OPENMP_SIMD_PRAGMA
              for (unsigned int i = 0; i < dst.locally_owned_size(); ++i)
                dst_ptr[i] = relaxation * src_ptr[i] * diag_ptr[i];
            }
        }
      else
        {
          tmp.reinit(src, true);

          Assert(transposed == false, ExcNotImplemented());

          A.vmult(
            tmp,
            dst,
            [&](const unsigned int start_range, const unsigned int end_range) {
              // zero 'tmp' before running the vmult operation
              if (end_range > start_range)
                std::memset(tmp.begin() + start_range,
                            0,
                            sizeof(Number) * (end_range - start_range));
            },
            [&](const unsigned int begin, const unsigned int end) {
              const Number *dst_ptr  = dst.begin();
              const Number *src_ptr  = src.begin();
              Number       *tmp_ptr  = tmp.begin();
              const Number *diag_ptr = preconditioner.get_vector().begin();

              // for efficiency reason, write back to temp_vector that is
              // already read (avoid read-for-ownership)
              if (relaxation == 1.0)
                {
                  DEAL_II_OPENMP_SIMD_PRAGMA
                  for (std::size_t i = begin; i < end; ++i)
                    tmp_ptr[i] =
                      dst_ptr[i] + (src_ptr[i] - tmp_ptr[i]) * diag_ptr[i];
                }
              else
                {
                  DEAL_II_OPENMP_SIMD_PRAGMA
                  for (std::size_t i = begin; i < end; ++i)
                    tmp_ptr[i] = dst_ptr[i] + relaxation *
                                                (src_ptr[i] - tmp_ptr[i]) *
                                                diag_ptr[i];
                }
            });

          tmp.swap(dst);
        }
    }

  } // namespace PreconditionRelaxation
} // namespace internal

#endif



/**
 * Jacobi preconditioner using matrix built-in function.  The
 * <tt>MatrixType</tt> class used is required to have a function
 * <tt>precondition_Jacobi(VectorType&, const VectorType&, double</tt>). This
 * class satisfies the
 * @ref ConceptRelaxationType "relaxation concept".
 *
 * @code
 * // Declare related objects
 *
 * SparseMatrix<double> A;
 * Vector<double> x;
 * Vector<double> b;
 * SolverCG<> solver(...);
 *
 * //...initialize and build A
 *
 * // Define and initialize preconditioner:
 *
 * PreconditionJacobi<SparseMatrix<double> > precondition;
 * precondition.initialize(
 *   A, PreconditionJacobi<SparseMatrix<double>>::AdditionalData(.6));
 *
 * solver.solve (A, x, b, precondition);
 * @endcode
 */
template <typename MatrixType = SparseMatrix<double>>
class PreconditionJacobi
  : public PreconditionRelaxation<
      MatrixType,
      internal::PreconditionRelaxation::PreconditionJacobiImpl<MatrixType>>
{
  using PreconditionerType =
    internal::PreconditionRelaxation::PreconditionJacobiImpl<MatrixType>;
  using BaseClass = PreconditionRelaxation<MatrixType, PreconditionerType>;

public:
  /**
   * An alias to the base class AdditionalData.
   */
  using AdditionalData = typename BaseClass::AdditionalData;

  /**
   * @copydoc PreconditionRelaxation::initialize()
   */
  void
  initialize(const MatrixType     &A,
             const AdditionalData &parameters = AdditionalData());
};


/**
 * SOR preconditioner using matrix built-in function.
 *
 * Assuming the matrix <i>A = D + L + U</i> is split into its diagonal
 * <i>D</i> as well as the strict lower and upper triangles <i>L</i> and
 * <i>U</i>, then the SOR preconditioner with relaxation parameter <i>r</i> is
 * @f[
 *  P^{-1} = r (D+rL)^{-1}.
 * @f]
 * It is this operator <i>P<sup>-1</sup></i>, which is implemented by vmult()
 * through forward substitution. Analogously, Tvmult() implements the
 * operation of <i>r(D+rU)<sup>-1</sup></i>.
 *
 * The SOR iteration itself can be directly written as
 * @f[
 *  x^{k+1} = x^k - r D^{-1} \bigl(L x^{k+1} + U x^k - b\bigr).
 * @f]
 * Using the right hand side <i>b</i> and the previous iterate <i>x</i>, this
 * is the operation implemented by step().
 *
 * The MatrixType class used is required to have functions
 * <tt>precondition_SOR(VectorType&, const VectorType&, double)</tt> and
 * <tt>precondition_TSOR(VectorType&, const VectorType&, double)</tt>. This
 * class satisfies the
 * @ref ConceptRelaxationType "relaxation concept".
 *
 * @code
 * // Declare related objects
 *
 * SparseMatrix<double> A;
 * Vector<double> x;
 * Vector<double> b;
 * SolverCG<> solver(...);
 *
 * //...initialize and build A
 *
 * // Define and initialize preconditioner
 *
 * PreconditionSOR<SparseMatrix<double> > precondition;
 * precondition.initialize(
 *   A, PreconditionSOR<SparseMatrix<double>>::AdditionalData(.6));
 *
 * solver.solve (A, x, b, precondition);
 * @endcode
 */
template <typename MatrixType = SparseMatrix<double>>
class PreconditionSOR
  : public PreconditionRelaxation<
      MatrixType,
      internal::PreconditionRelaxation::PreconditionSORImpl<MatrixType>>
{
  using PreconditionerType =
    internal::PreconditionRelaxation::PreconditionSORImpl<MatrixType>;
  using BaseClass = PreconditionRelaxation<MatrixType, PreconditionerType>;

public:
  /**
   * An alias to the base class AdditionalData.
   */
  using AdditionalData = typename BaseClass::AdditionalData;

  /**
   * @copydoc PreconditionRelaxation::initialize()
   */
  void
  initialize(const MatrixType     &A,
             const AdditionalData &parameters = AdditionalData());
};



/**
 * SSOR preconditioner using matrix built-in function.  The
 * <tt>MatrixType</tt> class used is required to have a function
 * <tt>precondition_SSOR(VectorType&, const VectorType&, double)</tt>. This
 * class satisfies the
 * @ref ConceptRelaxationType "relaxation concept".
 *
 * @code
 * // Declare related objects
 *
 * SparseMatrix<double> A;
 * Vector<double> x;
 * Vector<double> b;
 * SolverCG<> solver(...);
 *
 * //...initialize and build A
 *
 * // Define and initialize preconditioner
 *
 * PreconditionSSOR<SparseMatrix<double> > precondition;
 * precondition.initialize(
 *   A, PreconditionSSOR<SparseMatrix<double>>::AdditionalData(.6));
 *
 * solver.solve (A, x, b, precondition);
 * @endcode
 */
template <typename MatrixType = SparseMatrix<double>>
class PreconditionSSOR
  : public PreconditionRelaxation<
      MatrixType,
      internal::PreconditionRelaxation::PreconditionSSORImpl<MatrixType>>
{
  using PreconditionerType =
    internal::PreconditionRelaxation::PreconditionSSORImpl<MatrixType>;
  using BaseClass = PreconditionRelaxation<MatrixType, PreconditionerType>;

public:
  /**
   * An alias to the base class AdditionalData.
   */
  using AdditionalData = typename BaseClass::AdditionalData;

  /**
   * Initialize matrix and relaxation parameter. The matrix is just stored in
   * the preconditioner object. The relaxation parameter should be larger than
   * zero and smaller than 2 for numerical reasons. It defaults to 1.
   */
  void
  initialize(const MatrixType     &A,
             const AdditionalData &parameters = AdditionalData());
};


/**
 * Permuted SOR preconditioner using matrix built-in function.  The
 * <tt>MatrixType</tt> class used is required to have functions
 * <tt>PSOR(VectorType&, const VectorType&, double)</tt> and
 * <tt>TPSOR(VectorType&, const VectorType&, double)</tt>.
 *
 * @code
 * // Declare related objects
 *
 * SparseMatrix<double> A;
 * Vector<double> x;
 * Vector<double> b;
 * SolverCG<> solver(...);
 *
 * //...initialize and build A
 *
 * std::vector<unsigned int> permutation(x.size());
 * std::vector<unsigned int> inverse_permutation(x.size());
 *
 * //...fill permutation and its inverse with reasonable values
 *
 * // Define and initialize preconditioner
 *
 * PreconditionPSOR<SparseMatrix<double> > precondition;
 * precondition.initialize (A, permutation, inverse_permutation, .6);
 *
 * solver.solve (A, x, b, precondition);
 * @endcode
 */
template <typename MatrixType = SparseMatrix<double>>
class PreconditionPSOR
  : public PreconditionRelaxation<
      MatrixType,
      internal::PreconditionRelaxation::PreconditionPSORImpl<MatrixType>>
{
  using PreconditionerType =
    internal::PreconditionRelaxation::PreconditionPSORImpl<MatrixType>;
  using BaseClass = PreconditionRelaxation<MatrixType, PreconditionerType>;

public:
  /**
   * Declare type for container size.
   */
  using size_type = typename BaseClass::size_type;

  /**
   * Parameters for PreconditionPSOR.
   */
  class AdditionalData
  {
  public:
    /**
     * Constructor. For the parameters' description, see below.
     *
     * The permutation vectors are stored as a reference. Therefore, it has to
     * be assured that the lifetime of the vector exceeds the lifetime of the
     * preconditioner.
     *
     * The relaxation parameter should be larger than zero and smaller than 2
     * for numerical reasons. It defaults to 1.
     */
    AdditionalData(const std::vector<size_type> &permutation,
                   const std::vector<size_type> &inverse_permutation,
                   const typename BaseClass::AdditionalData &parameters =
                     typename BaseClass::AdditionalData());

    /**
     * Storage for the permutation vector.
     */
    const std::vector<size_type> &permutation;
    /**
     * Storage for the inverse permutation vector.
     */
    const std::vector<size_type> &inverse_permutation;
    /**
     * Relaxation parameters
     */
    typename BaseClass::AdditionalData parameters;
  };

  /**
   * Initialize matrix and relaxation parameter. The matrix is just stored in
   * the preconditioner object.
   *
   * The permutation vector is stored as a pointer. Therefore, it has to be
   * assured that the lifetime of the vector exceeds the lifetime of the
   * preconditioner.
   *
   * The relaxation parameter should be larger than zero and smaller than 2
   * for numerical reasons. It defaults to 1.
   */
  void
  initialize(const MatrixType                         &A,
             const std::vector<size_type>             &permutation,
             const std::vector<size_type>             &inverse_permutation,
             const typename BaseClass::AdditionalData &parameters =
               typename BaseClass::AdditionalData());

  /**
   * Initialize matrix and relaxation parameter. The matrix is just stored in
   * the preconditioner object.
   *
   * For more detail about possible parameters, see the class documentation
   * and the documentation of the PreconditionPSOR::AdditionalData class.
   *
   * After this function is called the preconditioner is ready to be used
   * (using the <code>vmult</code> function of derived classes).
   */
  void
  initialize(const MatrixType &A, const AdditionalData &additional_data);
};



/**
 * Preconditioning with a Chebyshev polynomial for symmetric positive definite
 * matrices. This preconditioner is based on an iteration of an inner
 * preconditioner of type @p PreconditionerType with coefficients that are
 * adapted to optimally cover an eigenvalue range between the largest
 * eigenvalue $\lambda_{\max{}}$ down to a given lower eigenvalue
 * $\lambda_{\min{}}$ specified by the optional parameter
 * @p smoothing_range. The algorithm is based on the following three-term
 * recurrence:
 * @f[
 *  x^{n+1} = x^{n} + \alpha^n_0 (x^{n} - x^{n-1}) + \alpha^n_1 P^{-1} (b-Ax^n)
 *    \quad\text{with}\quad
 *    \alpha^0_0 := 0,\;
 *    \alpha^0_1 := \frac{2\rho_0}{\lambda_{\max}-\lambda_{\min}}\;
 *    \alpha^n_0 := \rho_n \rho_{n-1},\;\text{and}\;
 *    \alpha^n_1 := \frac{4\rho_n}{\lambda_{\max}-\lambda_{\min}},
 * @f]
 * where the parameter $\rho_0$ is set to $\rho_0 =
 * \frac{\lambda_{\max{}}-\lambda_{\min{}}}{\lambda_{\max{}}+\lambda_{\min{}}}$
 * for the maximal eigenvalue $\lambda_{\max{}}$ and updated via $\rho_n =
 * \left(2\frac{\lambda_{\max{}}+\lambda_{\min{}}}
 * {\lambda_{\max{}}-\lambda_{\min{}}} - \rho_{n-1}\right)^{-1}$. The
 * Chebyshev polynomial is constructed to strongly damp the eigenvalue range
 * between $\lambda_{\min{}}$ and $\lambda_{\max{}}$ and is visualized e.g. in
 * Utilities::LinearAlgebra::chebyshev_filter().
 *
 * The typical use case for the preconditioner is a Jacobi preconditioner
 * specified through DiagonalMatrix, which is also the default value for the
 * preconditioner. Note that if the degree variable is set to one, the
 * Chebyshev iteration corresponds to a Jacobi preconditioner (or the
 * underlying preconditioner type) with relaxation parameter according to the
 * specified smoothing range.
 *
 * Besides the default choice of a pointwise Jacobi preconditioner, this class
 * also allows for more advanced types of preconditioners, for example
 * iterating block-Jacobi preconditioners in DG methods.
 *
 * Apart from the inner preconditioner object, this iteration does not need
 * access to matrix entries, which makes it an ideal ingredient for
 * matrix-free computations. In that context, this class can be used as a
 * multigrid smoother that is trivially %parallel (assuming that matrix-vector
 * products are %parallel and the inner preconditioner is %parallel). Its use
 * is demonstrated in the step-37 and step-59 tutorial programs.
 *
 * <h4>Estimation of the eigenvalues</h4>
 *
 * The Chebyshev method relies on an estimate of the eigenvalues of the matrix
 * which are computed during the first invocation of vmult(). This class
 * offers several algorithms to this end, see
 * PreconditionChebyshev::AdditionalData::EigenvalueAlgorithm. The default
 * algorithm invokes the Lanczos method via the SolverCG class, which requires
 * symmetry and positive definiteness of the (preconditioned) matrix system
 * are required. Also note that deal.II needs to be configured with LAPACK
 * support to use this option. The eigenvalue algorithm can be controlled by
 * PreconditionChebyshev::AdditionalData::eig_cg_n_iterations specifying how
 * many iterations should be performed. For all algorithms, the iterative
 * process is started from an initial vector that depends on the vector
 * type. For the classes dealii::Vector or
 * dealii::LinearAlgebra::distributed::Vector, which have fast element access,
 * it is a vector with entries `(-5.5, -4.5, -3.5, -2.5, ..., 3.5, 4.5, 5.5)`
 * with appropriate epilogue and adjusted such that its mean is always zero,
 * which works well for the Laplacian. This setup is stable in parallel in the
 * sense that for a different number of processors but the same ordering of
 * unknowns, the same initial vector and thus eigenvalue distribution will be
 * computed, apart from roundoff errors. For other vector types, the initial
 * vector contains all ones, scaled by the length of the vector, except for
 * the very first entry that is zero, triggering high-frequency content again.
 *
 * The computation of eigenvalues happens the first time one of the vmult(),
 * Tvmult(), step() or Tstep() functions is called or when
 * estimate_eigenvalues() is called directly. In the latter case, it is
 * necessary to provide a temporary vector of the same layout as the source
 * and destination vectors used during application of the preconditioner.
 *
 * The estimates for minimum and maximum eigenvalue are taken from the
 * underlying solver or eigenvalue algorithm in the given number of
 * iterations, even if the solver did not converge in the requested number of
 * iterations. Finally, the maximum eigenvalue is multiplied by a safety
 * factor of 1.2.
 *
 * Due to the cost of the eigenvalue estimate, this class is most appropriate
 * if it is applied repeatedly, e.g. in a smoother for a geometric multigrid
 * solver, that can in turn be used to solve several linear systems.
 *
 * <h4>Bypassing the eigenvalue computation</h4>
 *
 * In some contexts, the automatic eigenvalue computation of this class may
 * result in a bad quality, e.g. when the polynomial basis or numbering of
 * unknowns is such that the initial vector described above is a bad
 * choice. It is possible to bypass the automatic eigenvalue computation by
 * setting AdditionalData::eig_cg_n_iterations to zero, and provide the
 * variable AdditionalData::max_eigenvalue instead. The minimal eigenvalue is
 * implicitly specified via `max_eigenvalue/smoothing_range`.
 *
 * <h4>Using the PreconditionChebyshev as a solver</h4>
 *
 * If the range <tt>[max_eigenvalue/smoothing_range, max_eigenvalue]</tt>
 * contains all eigenvalues of the preconditioned matrix system and the degree
 * (i.e., number of iterations) is high enough, this class can also be used as
 * a direct solver. For an error estimation of the Chebyshev iteration that
 * can be used to determine the number of iteration, see @cite Varga2009.
 *
 * In order to use Chebyshev as a solver, set the degree to
 * numbers::invalid_unsigned_int to force the automatic computation of the
 * number of iterations needed to reach a given target tolerance. Note that
 * this currently only works for symmetric positive definite matrices with the
 * eigenvalue algorithm set to the conjugate gradient algorithm. In this case,
 * the target tolerance is read from the variable
 * PreconditionChebyshev::AdditionalData::smoothing_range (it needs to be a
 * number less than one to force any iterations obviously).
 *
 * For details on the algorithm, see section 5.1 of @cite Varga2009.
 *
 * <h4>Requirements on the templated classes</h4>
 *
 * The class `MatrixType` must be derived from
 * EnableObserverPointer because a ObserverPointer to `MatrixType`
 * is held in the class. In particular, this means that the matrix object needs
 * to persist during the lifetime of PreconditionChebyshev. The preconditioner
 * is held in a shared_ptr that is copied into the AdditionalData member
 * variable of the class, so the variable used for initialization can safely be
 * discarded after calling initialize(). Both the matrix and the preconditioner
 * need to provide
 * @p vmult() functions for the matrix-vector product and @p m() functions for
 * accessing the number of rows in the (square) matrix. Furthermore, the
 * matrix must provide <tt>el(i,i)</tt> methods for accessing the matrix
 * diagonal in case the preconditioner type is DiagonalMatrix. Even though
 * it is highly recommended to pass the inverse diagonal entries inside a
 * separate preconditioner object for implementing the Jacobi method (which is
 * the only possible way to operate this class when computing in %parallel
 * with MPI because there is no knowledge about the locally stored range of
 * entries that would be needed from the matrix alone), there is a backward
 * compatibility function that can extract the diagonal in case of a serial
 * computation.
 *
 * <h4>Optimized operations with specific `MatrixType` argument</h4>
 *
 * This class enables to embed the vector updates into the matrix-vector
 * product in case the `MatrixType` supports this. To this end, the
 * `VectorType` needs to be of type LinearAlgebra::distributed::Vector, the
 * `PreconditionerType` needs to be DiagonalMatrix, and the class `MatrixType`
 * needs to provide a function with the signature
 * @code
 * void MatrixType::vmult(
 *    VectorType &,
 *    const VectorType &,
 *    const std::function<void(const unsigned int, const unsigned int)> &,
 *    const std::function<void(const unsigned int, const unsigned int)> &) const
 * @endcode
 * where the two given functions run before and after the matrix-vector
 * product, respectively. They take as arguments a sub-range among the locally
 * owned elements of the vector, defined as half-open intervals. The intervals
 * are designed to be scheduled close to the time the matrix-vector product
 * touches upon the entries in the `src` and `dst` vectors, respectively, with
 * the requirement that
 * <ul>
 * <li> the matrix-vector product may only access an entry in `src` or `dst`
 * once the `operation_before_matrix_vector_product` has been run on that
 * vector entry; </li>
 * <li> `operation_after_matrix_vector_product` may run on a range of entries
 * `[i,j)` once the matrix-vector product does not access the entries `[i,j)`
 * in `src` and `dst` any more. </li>
 * </ul>
 * The motivation for this function is to increase data locality and hence
 * cache usage. For the example of a class similar to the one in the step-37
 * tutorial program, the implementation is
 * @code
 * void
 * vmult(LinearAlgebra::distributed::Vector<number> &      dst,
 *       const LinearAlgebra::distributed::Vector<number> &src,
 *       const std::function<void(const unsigned int, const unsigned int)>
 *         &operation_before_matrix_vector_product,
 *       const std::function<void(const unsigned int, const unsigned int)>
 *         &operation_after_matrix_vector_product) const
 * {
 *   data.cell_loop(&LaplaceOperator::local_apply,
 *                  this,
 *                  dst,
 *                  src,
 *                  operation_before_matrix_vector_product,
 *                  operation_after_matrix_vector_product);
 * }
 * @endcode
 * In terms of the Chebyshev iteration, the operation before the loop will
 * set `dst` to zero, whereas the operation after the loop performs the
 * iteration leading to $x^{n+1}$ described above, modifying the `dst` and
 * `src` vectors.
 */
template <typename MatrixType         = SparseMatrix<double>,
          typename VectorType         = Vector<double>,
          typename PreconditionerType = DiagonalMatrix<VectorType>>
class PreconditionChebyshev : public EnableObserverPointer
{
public:
  /**
   * Declare type for container size.
   */
  using size_type = types::global_dof_index;

  /**
   * Standardized data struct to pipe additional parameters to the
   * preconditioner.
   */
  struct AdditionalData
    : public internal::EigenvalueAlgorithmAdditionalData<PreconditionerType>
  {
    using EigenvalueAlgorithm = internal::EigenvalueAlgorithm;

    /**
     * An enum to define the available types of polynomial types.
     */
    enum class PolynomialType
    {
      /**
       * First-kind Chebyshev polynomials.
       */
      first_kind,
      /**
       * Fourth-kind Chebyshev polynomials according to @cite lottes2022optimal
       * and @cite phillips2022optimal.
       */
      fourth_kind
    };

    /**
     * Constructor.
     */
    AdditionalData(
      const unsigned int        degree              = 1,
      const double              smoothing_range     = 0.,
      const unsigned int        eig_cg_n_iterations = 8,
      const double              eig_cg_residual     = 1e-2,
      const double              max_eigenvalue      = 1,
      const EigenvalueAlgorithm eigenvalue_algorithm =
        EigenvalueAlgorithm::lanczos,
      const PolynomialType polynomial_type = PolynomialType::first_kind);

    /**
     * This determines the degree of the Chebyshev polynomial. The degree of
     * the polynomial gives the number of matrix-vector products to be
     * performed for one application of the step() operation. During vmult(),
     * the method performs `(degree-1)` matrix-vector products. Degree one
     * corresponds to a damped Jacobi method.
     *
     * If the degree is set to numbers::invalid_unsigned_int, the algorithm
     * will automatically determine the number of necessary iterations based
     * on the usual Chebyshev error formula as mentioned in the discussion of
     * the main class.
     */
    unsigned int degree;

    /**
     * Specifies the polynomial type to be used.
     */
    PolynomialType polynomial_type;
  };


  /**
   * Constructor.
   */
  PreconditionChebyshev();

  /**
   * Initialize function. Takes the matrix which is used to form the
   * preconditioner, and additional flags if there are any. This function
   * works only if the input matrix has an operator <tt>el(i,i)</tt> for
   * accessing all the elements in the diagonal. Alternatively, the diagonal
   * can be supplied with the help of the AdditionalData field.
   *
   * This function calculates an estimate of the eigenvalue range of the
   * matrix weighted by its diagonal using a modified CG iteration in case the
   * given number of iterations is positive.
   */
  void
  initialize(const MatrixType     &matrix,
             const AdditionalData &additional_data = AdditionalData());

  /**
   * Compute the action of the preconditioner on <tt>src</tt>, storing the
   * result in <tt>dst</tt>.
   */
  void
  vmult(VectorType &dst, const VectorType &src) const;

  /**
   * Compute the action of the transposed preconditioner on <tt>src</tt>,
   * storing the result in <tt>dst</tt>.
   */
  void
  Tvmult(VectorType &dst, const VectorType &src) const;

  /**
   * Perform one step of the preconditioned Richardson iteration.
   */
  void
  step(VectorType &dst, const VectorType &src) const;

  /**
   * Perform one transposed step of the preconditioned Richardson iteration.
   */
  void
  Tstep(VectorType &dst, const VectorType &src) const;

  /**
   * Resets the preconditioner.
   */
  void
  clear();

  /**
   * Return the dimension of the codomain (or range) space. Note that the
   * matrix is of dimension $m \times n$.
   */
  size_type
  m() const;

  /**
   * Return the dimension of the domain space. Note that the matrix is of
   * dimension $m \times n$.
   */
  size_type
  n() const;

  using EigenvalueInformation = internal::EigenvalueInformation;

  /**
   * Compute eigenvalue estimates required for the preconditioner.
   *
   * This function is called automatically on first use of the preconditioner
   * if it is not called by the user. The layout of the vector @p src is used
   * to create internal temporary vectors and its content does not matter.
   *
   * Initializes the factors theta and delta based on an eigenvalue
   * computation. If the user set provided values for the largest eigenvalue
   * in AdditionalData, no computation is performed and the information given
   * by the user is used.
   */
  EigenvalueInformation
  estimate_eigenvalues(const VectorType &src) const;

private:
  /**
   * A pointer to the underlying matrix.
   */
  ObserverPointer<
    const MatrixType,
    PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>>
    matrix_ptr;

  /**
   * Internal vector used for the <tt>vmult</tt> operation.
   */
  mutable VectorType solution_old;

  /**
   * Internal vector used for the <tt>vmult</tt> operation.
   */
  mutable VectorType temp_vector1;

  /**
   * Internal vector used for the <tt>vmult</tt> operation.
   */
  mutable VectorType temp_vector2;

  /**
   * Stores the additional data passed to the initialize function, obtained
   * through a copy operation.
   */
  AdditionalData data;

  /**
   * Average of the largest and smallest eigenvalue under consideration.
   */
  double theta;

  /**
   * Half the interval length between the largest and smallest eigenvalue
   * under consideration.
   */
  double delta;

  /**
   * Stores whether the preconditioner has been set up and eigenvalues have
   * been computed.
   */
  bool eigenvalues_are_initialized;

  /**
   * A mutex to avoid that multiple vmult() invocations by different threads
   * overwrite the temporary vectors.
   */
  mutable Threads::Mutex mutex;
};



/** @} */
/* ---------------------------------- Inline functions ------------------- */

#ifndef DOXYGEN


namespace internal
{
  template <typename VectorType>
  void
  set_initial_guess(VectorType &vector)
  {
    vector = 1. / std::sqrt(static_cast<double>(vector.size()));
    if (vector.locally_owned_elements().is_element(0))
      vector(0) = 0.;
  }

  template <typename Number>
  void
  set_initial_guess(::dealii::Vector<Number> &vector)
  {
    // Choose a high-frequency mode consisting of numbers between 0 and 1
    // that is cheap to compute (cheaper than random numbers) but avoids
    // obviously re-occurring numbers in multi-component systems by choosing
    // a period of 11
    for (unsigned int i = 0; i < vector.size(); ++i)
      vector(i) = i % 11;

    const Number mean_value = vector.mean_value();
    vector.add(-mean_value);
  }

  template <typename Number, typename MemorySpace>
  void
  set_initial_guess(
    ::dealii::LinearAlgebra::distributed::BlockVector<Number, MemorySpace>
      &vector)
  {
    for (unsigned int block = 0; block < vector.n_blocks(); ++block)
      set_initial_guess(vector.block(block));
  }

  template <typename Number, typename MemorySpace>
  void
  set_initial_guess(
    ::dealii::LinearAlgebra::distributed::Vector<Number, MemorySpace> &vector)
  {
    // Choose a high-frequency mode consisting of numbers between 0 and 1
    // that is cheap to compute (cheaper than random numbers) but avoids
    // obviously re-occurring numbers in multi-component systems by choosing
    // a period of 11.
    // Make initial guess robust with respect to number of processors
    // by operating on the global index.
    types::global_dof_index first_local_range = 0;
    if (!vector.locally_owned_elements().is_empty())
      first_local_range = vector.locally_owned_elements().nth_index_in_set(0);

    const auto n_local_elements = vector.locally_owned_size();
    Number    *values_ptr       = vector.get_values();
    Kokkos::RangePolicy<typename MemorySpace::kokkos_space::execution_space,
                        Kokkos::IndexType<types::global_dof_index>>
      policy(0, n_local_elements);
    Kokkos::parallel_for(
      "dealii::PreconditionChebyshev::set_initial_guess",
      policy,
      KOKKOS_LAMBDA(types::global_dof_index i) {
        values_ptr[i] = (i + first_local_range) % 11;
      });
    const Number mean_value = vector.mean_value();
    vector.add(-mean_value);
  }

  struct EigenvalueTracker
  {
  public:
    void
    slot(const std::vector<double> &eigenvalues)
    {
      values = eigenvalues;
    }

    std::vector<double> values;
  };



  template <typename MatrixType,
            typename VectorType,
            typename PreconditionerType>
  double
  power_iteration(const MatrixType         &matrix,
                  VectorType               &eigenvector,
                  const PreconditionerType &preconditioner,
                  const unsigned int        n_iterations)
  {
    typename VectorType::value_type eigenvalue_estimate = 0.;
    eigenvector /= eigenvector.l2_norm();
    VectorType vector1, vector2;
    vector1.reinit(eigenvector, true);
    if (!std::is_same_v<PreconditionerType, PreconditionIdentity>)
      vector2.reinit(eigenvector, true);
    for (unsigned int i = 0; i < n_iterations; ++i)
      {
        if (!std::is_same_v<PreconditionerType, PreconditionIdentity>)
          {
            matrix.vmult(vector2, eigenvector);
            preconditioner.vmult(vector1, vector2);
          }
        else
          matrix.vmult(vector1, eigenvector);

        eigenvalue_estimate = eigenvector * vector1;

        vector1 /= vector1.l2_norm();
        eigenvector.swap(vector1);
      }
    return std::abs(eigenvalue_estimate);
  }



  template <typename MatrixType,
            typename VectorType,
            typename PreconditionerType>
  EigenvalueInformation
  estimate_eigenvalues(
    const EigenvalueAlgorithmAdditionalData<PreconditionerType> &data,
    const MatrixType                                            *matrix_ptr,
    VectorType                                                  &solution_old,
    VectorType                                                  &temp_vector1,
    const unsigned int                                           degree)
  {
    Assert(data.preconditioner.get() != nullptr, ExcNotInitialized());

    EigenvalueInformation info{};

    if (data.eig_cg_n_iterations > 0)
      {
        Assert(data.eig_cg_n_iterations > 2,
               ExcMessage(
                 "Need to set at least two iterations to find eigenvalues."));

        internal::EigenvalueTracker eigenvalue_tracker;

        // set an initial guess that contains some high-frequency parts (to the
        // extent possible without knowing the discretization and the numbering)
        // to trigger high eigenvalues according to the external function
        internal::set_initial_guess(temp_vector1);
        data.constraints.set_zero(temp_vector1);

        if (data.eigenvalue_algorithm == internal::EigenvalueAlgorithm::lanczos)
          {
            // set a very strict tolerance to force at least two iterations
            IterationNumberControl control(data.eig_cg_n_iterations,
                                           1e-10,
                                           false,
                                           false);

            dealii::SolverCG<VectorType> solver(control);
            solver.connect_eigenvalues_slot(
              [&eigenvalue_tracker](const std::vector<double> &eigenvalues) {
                eigenvalue_tracker.slot(eigenvalues);
              });

            solver.solve(*matrix_ptr,
                         solution_old,
                         temp_vector1,
                         *data.preconditioner);

            info.cg_iterations = control.last_step();
          }
        else if (data.eigenvalue_algorithm ==
                 internal::EigenvalueAlgorithm::power_iteration)
          {
            (void)degree;

            Assert(degree != numbers::invalid_unsigned_int,
                   ExcMessage("Cannot estimate the minimal eigenvalue with the "
                              "power iteration"));

            eigenvalue_tracker.values.push_back(
              internal::power_iteration(*matrix_ptr,
                                        temp_vector1,
                                        *data.preconditioner,
                                        data.eig_cg_n_iterations));
          }
        else
          DEAL_II_NOT_IMPLEMENTED();

        // read the eigenvalues from the attached eigenvalue tracker
        if (eigenvalue_tracker.values.empty())
          info.min_eigenvalue_estimate = info.max_eigenvalue_estimate = 1.;
        else
          {
            info.min_eigenvalue_estimate = eigenvalue_tracker.values.front();

            // include a safety factor since the CG method will in general not
            // be converged
            info.max_eigenvalue_estimate =
              1.2 * eigenvalue_tracker.values.back();
          }
      }
    else
      {
        info.max_eigenvalue_estimate = data.max_eigenvalue;
        info.min_eigenvalue_estimate =
          data.max_eigenvalue / data.smoothing_range;
      }

    return info;
  }
} // namespace internal



inline PreconditionIdentity::PreconditionIdentity()
  : n_rows(0)
  , n_columns(0)
{}

template <typename MatrixType>
inline void
PreconditionIdentity::initialize(const MatrixType &matrix,
                                 const PreconditionIdentity::AdditionalData &)
{
  n_rows    = matrix.m();
  n_columns = matrix.n();
}


template <typename VectorType>
inline void
PreconditionIdentity::vmult(VectorType &dst, const VectorType &src) const
{
  dst = src;
}



template <typename VectorType>
inline void
PreconditionIdentity::Tvmult(VectorType &dst, const VectorType &src) const
{
  dst = src;
}

template <typename VectorType>
inline void
PreconditionIdentity::vmult_add(VectorType &dst, const VectorType &src) const
{
  dst += src;
}



template <typename VectorType>
inline void
PreconditionIdentity::Tvmult_add(VectorType &dst, const VectorType &src) const
{
  dst += src;
}



inline void
PreconditionIdentity::clear()
{}



inline PreconditionIdentity::size_type
PreconditionIdentity::m() const
{
  Assert(n_rows != 0, ExcNotInitialized());
  return n_rows;
}

inline PreconditionIdentity::size_type
PreconditionIdentity::n() const
{
  Assert(n_columns != 0, ExcNotInitialized());
  return n_columns;
}

//---------------------------------------------------------------------------

inline PreconditionRichardson::AdditionalData::AdditionalData(
  const double relaxation)
  : relaxation(relaxation)
{}


inline PreconditionRichardson::PreconditionRichardson()
  : relaxation(0)
  , n_rows(0)
  , n_columns(0)
{
  AdditionalData add_data;
  relaxation = add_data.relaxation;
}



inline void
PreconditionRichardson::initialize(
  const PreconditionRichardson::AdditionalData &parameters)
{
  relaxation = parameters.relaxation;
}



template <typename MatrixType>
inline void
PreconditionRichardson::initialize(
  const MatrixType                             &matrix,
  const PreconditionRichardson::AdditionalData &parameters)
{
  relaxation = parameters.relaxation;
  n_rows     = matrix.m();
  n_columns  = matrix.n();
}



template <typename VectorType>
inline void
PreconditionRichardson::vmult(VectorType &dst, const VectorType &src) const
{
  static_assert(
    std::is_same_v<size_type, typename VectorType::size_type>,
    "PreconditionRichardson and VectorType must have the same size_type.");

  dst.equ(relaxation, src);
}



template <typename VectorType>
inline void
PreconditionRichardson::Tvmult(VectorType &dst, const VectorType &src) const
{
  static_assert(
    std::is_same_v<size_type, typename VectorType::size_type>,
    "PreconditionRichardson and VectorType must have the same size_type.");

  dst.equ(relaxation, src);
}

template <typename VectorType>
inline void
PreconditionRichardson::vmult_add(VectorType &dst, const VectorType &src) const
{
  static_assert(
    std::is_same_v<size_type, typename VectorType::size_type>,
    "PreconditionRichardson and VectorType must have the same size_type.");

  dst.add(relaxation, src);
}



template <typename VectorType>
inline void
PreconditionRichardson::Tvmult_add(VectorType &dst, const VectorType &src) const
{
  static_assert(
    std::is_same_v<size_type, typename VectorType::size_type>,
    "PreconditionRichardson and VectorType must have the same size_type.");

  dst.add(relaxation, src);
}

inline PreconditionRichardson::size_type
PreconditionRichardson::m() const
{
  Assert(n_rows != 0, ExcNotInitialized());
  return n_rows;
}

inline PreconditionRichardson::size_type
PreconditionRichardson::n() const
{
  Assert(n_columns != 0, ExcNotInitialized());
  return n_columns;
}

//---------------------------------------------------------------------------

template <typename MatrixType, typename PreconditionerType>
inline void
PreconditionRelaxation<MatrixType, PreconditionerType>::initialize(
  const MatrixType     &rA,
  const AdditionalData &parameters)
{
  A                           = &rA;
  eigenvalues_are_initialized = false;

  Assert(parameters.preconditioner, ExcNotInitialized());

  this->data = parameters;
}


template <typename MatrixType, typename PreconditionerType>
inline void
PreconditionRelaxation<MatrixType, PreconditionerType>::clear()
{
  eigenvalues_are_initialized = false;
  A                           = nullptr;
  data.relaxation             = 1.0;
  data.preconditioner         = nullptr;
}

template <typename MatrixType, typename PreconditionerType>
inline
  typename PreconditionRelaxation<MatrixType, PreconditionerType>::size_type
  PreconditionRelaxation<MatrixType, PreconditionerType>::m() const
{
  Assert(A != nullptr, ExcNotInitialized());
  return A->m();
}

template <typename MatrixType, typename PreconditionerType>
inline
  typename PreconditionRelaxation<MatrixType, PreconditionerType>::size_type
  PreconditionRelaxation<MatrixType, PreconditionerType>::n() const
{
  Assert(A != nullptr, ExcNotInitialized());
  return A->n();
}

template <typename MatrixType, typename PreconditionerType>
template <typename VectorType>
inline void
PreconditionRelaxation<MatrixType, PreconditionerType>::vmult(
  VectorType       &dst,
  const VectorType &src) const
{
  Assert(this->A != nullptr, ExcNotInitialized());
  Assert(data.preconditioner != nullptr, ExcNotInitialized());

  if (eigenvalues_are_initialized == false)
    estimate_eigenvalues(src);

  VectorType tmp1, tmp2;

  for (unsigned int i = 0; i < data.n_iterations; ++i)
    internal::PreconditionRelaxation::step_operations(*A,
                                                      *data.preconditioner,
                                                      dst,
                                                      src,
                                                      data.relaxation,
                                                      tmp1,
                                                      tmp2,
                                                      i,
                                                      false);
}

template <typename MatrixType, typename PreconditionerType>
template <typename VectorType>
inline void
PreconditionRelaxation<MatrixType, PreconditionerType>::Tvmult(
  VectorType       &dst,
  const VectorType &src) const
{
  Assert(this->A != nullptr, ExcNotInitialized());
  Assert(data.preconditioner != nullptr, ExcNotInitialized());

  if (eigenvalues_are_initialized == false)
    estimate_eigenvalues(src);

  VectorType tmp1, tmp2;

  for (unsigned int i = 0; i < data.n_iterations; ++i)
    internal::PreconditionRelaxation::step_operations(
      *A, *data.preconditioner, dst, src, data.relaxation, tmp1, tmp2, i, true);
}

template <typename MatrixType, typename PreconditionerType>
template <typename VectorType>
inline void
PreconditionRelaxation<MatrixType, PreconditionerType>::step(
  VectorType       &dst,
  const VectorType &src) const
{
  Assert(this->A != nullptr, ExcNotInitialized());
  Assert(data.preconditioner != nullptr, ExcNotInitialized());

  if (eigenvalues_are_initialized == false)
    estimate_eigenvalues(src);

  VectorType tmp1, tmp2;

  for (unsigned int i = 1; i <= data.n_iterations; ++i)
    internal::PreconditionRelaxation::step_operations(*A,
                                                      *data.preconditioner,
                                                      dst,
                                                      src,
                                                      data.relaxation,
                                                      tmp1,
                                                      tmp2,
                                                      i,
                                                      false);
}

template <typename MatrixType, typename PreconditionerType>
template <typename VectorType>
inline void
PreconditionRelaxation<MatrixType, PreconditionerType>::Tstep(
  VectorType       &dst,
  const VectorType &src) const
{
  Assert(this->A != nullptr, ExcNotInitialized());
  Assert(data.preconditioner != nullptr, ExcNotInitialized());

  if (eigenvalues_are_initialized == false)
    estimate_eigenvalues(src);

  VectorType tmp1, tmp2;

  for (unsigned int i = 1; i <= data.n_iterations; ++i)
    internal::PreconditionRelaxation::step_operations(
      *A, *data.preconditioner, dst, src, data.relaxation, tmp1, tmp2, i, true);
}

template <typename MatrixType, typename PreconditionerType>
template <typename VectorType>
inline internal::EigenvalueInformation
PreconditionRelaxation<MatrixType, PreconditionerType>::estimate_eigenvalues(
  const VectorType &src) const
{
  Assert(eigenvalues_are_initialized == false, ExcInternalError());

  EigenvalueInformation info;

  if (data.relaxation == 0.0)
    {
      VectorType solution_old, temp_vector1;

      solution_old.reinit(src);
      temp_vector1.reinit(src, true);

      info = internal::estimate_eigenvalues<MatrixType>(
        data, A, solution_old, temp_vector1, data.n_iterations);

      const double alpha =
        (data.smoothing_range > 1. ?
           info.max_eigenvalue_estimate / data.smoothing_range :
           std::min(0.9 * info.max_eigenvalue_estimate,
                    info.min_eigenvalue_estimate));

      const_cast<PreconditionRelaxation<MatrixType, PreconditionerType> *>(this)
        ->data.relaxation = 2.0 / (alpha + info.max_eigenvalue_estimate);
    }

  const_cast<PreconditionRelaxation<MatrixType, PreconditionerType> *>(this)
    ->eigenvalues_are_initialized = true;

  return info;
}

template <typename MatrixType, typename PreconditionerType>
double
PreconditionRelaxation<MatrixType, PreconditionerType>::get_relaxation() const
{
  return data.relaxation;
}


//---------------------------------------------------------------------------

template <typename MatrixType>
inline void
PreconditionJacobi<MatrixType>::initialize(const MatrixType     &A,
                                           const AdditionalData &parameters_in)
{
  Assert(parameters_in.preconditioner == nullptr, ExcInternalError());
  Assert(
    parameters_in.relaxation != 0.0,
    ExcMessage(
      "Relaxation cannot automatically be determined by PreconditionJacobi."));

  AdditionalData parameters;
  parameters.relaxation   = 1.0;
  parameters.n_iterations = parameters_in.n_iterations;
  parameters.preconditioner =
    std::make_shared<PreconditionerType>(A, parameters_in.relaxation);

  this->BaseClass::initialize(A, parameters);
}

//---------------------------------------------------------------------------

template <typename MatrixType>
inline void
PreconditionSOR<MatrixType>::initialize(const MatrixType     &A,
                                        const AdditionalData &parameters_in)
{
  Assert(parameters_in.preconditioner == nullptr, ExcInternalError());
  Assert(
    parameters_in.relaxation != 0.0,
    ExcMessage(
      "Relaxation cannot automatically be determined by PreconditionSOR."));

  AdditionalData parameters;
  parameters.relaxation   = 1.0;
  parameters.n_iterations = parameters_in.n_iterations;
  parameters.preconditioner =
    std::make_shared<PreconditionerType>(A, parameters_in.relaxation);

  this->BaseClass::initialize(A, parameters);
}

//---------------------------------------------------------------------------

template <typename MatrixType>
inline void
PreconditionSSOR<MatrixType>::initialize(const MatrixType     &A,
                                         const AdditionalData &parameters_in)
{
  Assert(parameters_in.preconditioner == nullptr, ExcInternalError());
  Assert(
    parameters_in.relaxation != 0.0,
    ExcMessage(
      "Relaxation cannot automatically be determined by PreconditionSSOR."));

  AdditionalData parameters;
  parameters.relaxation   = 1.0;
  parameters.n_iterations = parameters_in.n_iterations;
  parameters.preconditioner =
    std::make_shared<PreconditionerType>(A, parameters_in.relaxation);

  this->BaseClass::initialize(A, parameters);
}



//---------------------------------------------------------------------------

template <typename MatrixType>
inline void
PreconditionPSOR<MatrixType>::initialize(
  const MatrixType                         &A,
  const std::vector<size_type>             &p,
  const std::vector<size_type>             &ip,
  const typename BaseClass::AdditionalData &parameters_in)
{
  Assert(parameters_in.preconditioner == nullptr, ExcInternalError());
  Assert(
    parameters_in.relaxation != 0.0,
    ExcMessage(
      "Relaxation cannot automatically be determined by PreconditionPSOR."));

  typename BaseClass::AdditionalData parameters;
  parameters.relaxation   = 1.0;
  parameters.n_iterations = parameters_in.n_iterations;
  parameters.preconditioner =
    std::make_shared<PreconditionerType>(A, parameters_in.relaxation, p, ip);

  this->BaseClass::initialize(A, parameters);
}


template <typename MatrixType>
inline void
PreconditionPSOR<MatrixType>::initialize(const MatrixType     &A,
                                         const AdditionalData &additional_data)
{
  initialize(A,
             additional_data.permutation,
             additional_data.inverse_permutation,
             additional_data.parameters);
}

template <typename MatrixType>
PreconditionPSOR<MatrixType>::AdditionalData::AdditionalData(
  const std::vector<size_type> &permutation,
  const std::vector<size_type> &inverse_permutation,
  const typename PreconditionPSOR<MatrixType>::BaseClass::AdditionalData
    &parameters)
  : permutation(permutation)
  , inverse_permutation(inverse_permutation)
  , parameters(parameters)
{}


//---------------------------------------------------------------------------


template <typename MatrixType, typename VectorType>
PreconditionUseMatrix<MatrixType, VectorType>::PreconditionUseMatrix(
  const MatrixType  &M,
  const function_ptr method)
  : matrix(M)
  , precondition(method)
{}



template <typename MatrixType, typename VectorType>
void
PreconditionUseMatrix<MatrixType, VectorType>::vmult(
  VectorType       &dst,
  const VectorType &src) const
{
  (matrix.*precondition)(dst, src);
}

//---------------------------------------------------------------------------

namespace internal
{

  template <typename PreconditionerType>
  inline EigenvalueAlgorithmAdditionalData<PreconditionerType>::
    EigenvalueAlgorithmAdditionalData(
      const double              smoothing_range,
      const unsigned int        eig_cg_n_iterations,
      const double              eig_cg_residual,
      const double              max_eigenvalue,
      const EigenvalueAlgorithm eigenvalue_algorithm)
    : smoothing_range(smoothing_range)
    , eig_cg_n_iterations(eig_cg_n_iterations)
    , eig_cg_residual(eig_cg_residual)
    , max_eigenvalue(max_eigenvalue)
    , eigenvalue_algorithm(eigenvalue_algorithm)
  {}



  template <typename PreconditionerType>
  inline EigenvalueAlgorithmAdditionalData<PreconditionerType> &
  EigenvalueAlgorithmAdditionalData<PreconditionerType>::operator=(
    const EigenvalueAlgorithmAdditionalData &other_data)
  {
    smoothing_range      = other_data.smoothing_range;
    eig_cg_n_iterations  = other_data.eig_cg_n_iterations;
    eig_cg_residual      = other_data.eig_cg_residual;
    max_eigenvalue       = other_data.max_eigenvalue;
    preconditioner       = other_data.preconditioner;
    eigenvalue_algorithm = other_data.eigenvalue_algorithm;
    constraints.copy_from(other_data.constraints);

    return *this;
  }
} // namespace internal

template <typename MatrixType, typename PreconditionerType>
inline PreconditionRelaxation<MatrixType, PreconditionerType>::AdditionalData::
  AdditionalData(const double              relaxation,
                 const unsigned int        n_iterations,
                 const double              smoothing_range,
                 const unsigned int        eig_cg_n_iterations,
                 const double              eig_cg_residual,
                 const double              max_eigenvalue,
                 const EigenvalueAlgorithm eigenvalue_algorithm)
  : internal::EigenvalueAlgorithmAdditionalData<PreconditionerType>(
      smoothing_range,
      eig_cg_n_iterations,
      eig_cg_residual,
      max_eigenvalue,
      eigenvalue_algorithm)
  , relaxation(relaxation)
  , n_iterations(n_iterations)
{}



//---------------------------------------------------------------------------

namespace internal
{
  namespace PreconditionChebyshevImplementation
  {
    // for deal.II vectors, perform updates for Chebyshev preconditioner all
    // at once to reduce memory transfer. Here, we select between general
    // vectors and deal.II vectors where we expand the loop over the (local)
    // size of the vector

    // generic part for non-deal.II vectors
    template <typename VectorType, typename PreconditionerType>
    inline void
    vector_updates(const VectorType         &rhs,
                   const PreconditionerType &preconditioner,
                   const unsigned int        iteration_index,
                   const double              factor1,
                   const double              factor2,
                   VectorType               &solution_old,
                   VectorType               &temp_vector1,
                   VectorType               &temp_vector2,
                   VectorType               &solution)
    {
      if (iteration_index == 0)
        {
          solution.equ(factor2, rhs);
          preconditioner.vmult(solution_old, solution);
        }
      else if (iteration_index == 1)
        {
          // compute t = P^{-1} * (b-A*x^{n})
          temp_vector1.sadd(-1.0, 1.0, rhs);
          preconditioner.vmult(solution_old, temp_vector1);

          // compute x^{n+1} = x^{n} + f_1 * x^{n} + f_2 * t
          solution_old.sadd(factor2, 1 + factor1, solution);
        }
      else
        {
          // compute t = P^{-1} * (b-A*x^{n})
          temp_vector1.sadd(-1.0, 1.0, rhs);
          preconditioner.vmult(temp_vector2, temp_vector1);

          // compute x^{n+1} = x^{n} + f_1 * (x^{n}-x^{n-1}) + f_2 * t
          solution_old.sadd(-factor1, factor2, temp_vector2);
          solution_old.add(1 + factor1, solution);
        }

      solution.swap(solution_old);
    }

    // generic part for deal.II vectors
    template <
      typename Number,
      typename PreconditionerType,
      std::enable_if_t<
        !has_vmult_with_std_functions_for_precondition<
          PreconditionerType,
          LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>>,
        int> * = nullptr>
    inline void
    vector_updates(
      const LinearAlgebra::distributed::Vector<Number, MemorySpace::Host> &rhs,
      const PreconditionerType &preconditioner,
      const unsigned int        iteration_index,
      const double              factor1_,
      const double              factor2_,
      LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>
        &solution_old,
      LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>
        &temp_vector1,
      LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>
        &temp_vector2,
      LinearAlgebra::distributed::Vector<Number, MemorySpace::Host> &solution)
    {
      const Number factor1        = factor1_;
      const Number factor1_plus_1 = 1. + factor1_;
      const Number factor2        = factor2_;

      if (iteration_index == 0)
        {
          // compute t = P^{-1} * (b)
          preconditioner.vmult(solution_old, rhs);

          // compute x^{n+1} = f_2 * t
          const auto solution_old_ptr = solution_old.begin();
          DEAL_II_OPENMP_SIMD_PRAGMA
          for (unsigned int i = 0; i < solution_old.locally_owned_size(); ++i)
            solution_old_ptr[i] = solution_old_ptr[i] * factor2;
        }
      else if (iteration_index == 1)
        {
          // compute t = P^{-1} * (b-A*x^{n})
          temp_vector1.sadd(-1.0, 1.0, rhs);

          preconditioner.vmult(solution_old, temp_vector1);

          // compute x^{n+1} = x^{n} + f_1 * x^{n} + f_2 * t
          const auto solution_ptr     = solution.begin();
          const auto solution_old_ptr = solution_old.begin();
          DEAL_II_OPENMP_SIMD_PRAGMA
          for (unsigned int i = 0; i < solution_old.locally_owned_size(); ++i)
            solution_old_ptr[i] =
              factor1_plus_1 * solution_ptr[i] + solution_old_ptr[i] * factor2;
        }
      else
        {
          // compute t = P^{-1} * (b-A*x^{n})
          temp_vector1.sadd(-1.0, 1.0, rhs);

          preconditioner.vmult(temp_vector2, temp_vector1);

          // compute x^{n+1} = x^{n} + f_1 * (x^{n}-x^{n-1}) + f_2 * t
          const auto solution_ptr     = solution.begin();
          const auto solution_old_ptr = solution_old.begin();
          const auto temp_vector2_ptr = temp_vector2.begin();
          DEAL_II_OPENMP_SIMD_PRAGMA
          for (unsigned int i = 0; i < solution_old.locally_owned_size(); ++i)
            solution_old_ptr[i] = factor1_plus_1 * solution_ptr[i] -
                                  factor1 * solution_old_ptr[i] +
                                  temp_vector2_ptr[i] * factor2;
        }

      solution.swap(solution_old);
    }

    template <
      typename Number,
      typename PreconditionerType,
      std::enable_if_t<
        has_vmult_with_std_functions_for_precondition<
          PreconditionerType,
          LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>>,
        int> * = nullptr>
    inline void
    vector_updates(
      const LinearAlgebra::distributed::Vector<Number, MemorySpace::Host> &rhs,
      const PreconditionerType &preconditioner,
      const unsigned int        iteration_index,
      const double              factor1_,
      const double              factor2_,
      LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>
        &solution_old,
      LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>
        &temp_vector1,
      LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>
        &temp_vector2,
      LinearAlgebra::distributed::Vector<Number, MemorySpace::Host> &solution)
    {
      const Number factor1        = factor1_;
      const Number factor1_plus_1 = 1. + factor1_;
      const Number factor2        = factor2_;

      const auto rhs_ptr          = rhs.begin();
      const auto temp_vector1_ptr = temp_vector1.begin();
      const auto temp_vector2_ptr = temp_vector2.begin();
      const auto solution_ptr     = solution.begin();
      const auto solution_old_ptr = solution_old.begin();

      if (iteration_index == 0)
        {
          preconditioner.vmult(
            solution,
            rhs,
            [&](const auto start_range, const auto end_range) {
              if (end_range > start_range)
                std::memset(solution_ptr + start_range,
                            0,
                            sizeof(Number) * (end_range - start_range));
            },
            [&](const auto begin, const auto end) {
              DEAL_II_OPENMP_SIMD_PRAGMA
              for (std::size_t i = begin; i < end; ++i)
                solution_ptr[i] *= factor2;
            });
        }
      else
        {
          preconditioner.vmult(
            temp_vector2,
            temp_vector1,
            [&](const auto begin, const auto end) {
              if (end > begin)
                std::memset(temp_vector2_ptr + begin,
                            0,
                            sizeof(Number) * (end - begin));

              DEAL_II_OPENMP_SIMD_PRAGMA
              for (std::size_t i = begin; i < end; ++i)
                temp_vector1_ptr[i] = rhs_ptr[i] - temp_vector1_ptr[i];
            },
            [&](const auto begin, const auto end) {
              if (iteration_index == 1)
                {
                  DEAL_II_OPENMP_SIMD_PRAGMA
                  for (std::size_t i = begin; i < end; ++i)
                    temp_vector2_ptr[i] = factor1_plus_1 * solution_ptr[i] +
                                          factor2 * temp_vector2_ptr[i];
                }
              else
                {
                  DEAL_II_OPENMP_SIMD_PRAGMA
                  for (std::size_t i = begin; i < end; ++i)
                    temp_vector2_ptr[i] = factor1_plus_1 * solution_ptr[i] -
                                          factor1 * solution_old_ptr[i] +
                                          factor2 * temp_vector2_ptr[i];
                }
            });
        }

      if (iteration_index > 0)
        {
          solution_old.swap(temp_vector2);
          solution_old.swap(solution);
        }
    }

    // worker routine for deal.II vectors. Because of vectorization, we need
    // to put the loop into an extra structure because the virtual function of
    // VectorUpdatesRange prevents the compiler from applying vectorization.
    template <typename Number>
    struct VectorUpdater
    {
      VectorUpdater(const Number      *rhs,
                    const Number      *matrix_diagonal_inverse,
                    const unsigned int iteration_index,
                    const Number       factor1,
                    const Number       factor2,
                    Number            *solution_old,
                    Number            *tmp_vector,
                    Number            *solution)
        : rhs(rhs)
        , matrix_diagonal_inverse(matrix_diagonal_inverse)
        , iteration_index(iteration_index)
        , factor1(factor1)
        , factor2(factor2)
        , solution_old(solution_old)
        , tmp_vector(tmp_vector)
        , solution(solution)
      {}

      void
      apply_to_subrange(const std::size_t begin, const std::size_t end) const
      {
        // To circumvent a bug in gcc
        // (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=63945), we create
        // copies of the variables factor1 and factor2 and do not check based on
        // factor1.
        const Number factor1        = this->factor1;
        const Number factor1_plus_1 = 1. + this->factor1;
        const Number factor2        = this->factor2;
        if (iteration_index == 0)
          {
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (std::size_t i = begin; i < end; ++i)
              solution[i] = factor2 * matrix_diagonal_inverse[i] * rhs[i];
          }
        else if (iteration_index == 1)
          {
            // x^{n+1} = x^{n} + f_1 * x^{n} + f_2 * P^{-1} * (b-A*x^{n})
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (std::size_t i = begin; i < end; ++i)
              // for efficiency reason, write back to temp_vector that is
              // already read (avoid read-for-ownership)
              tmp_vector[i] =
                factor1_plus_1 * solution[i] +
                factor2 * matrix_diagonal_inverse[i] * (rhs[i] - tmp_vector[i]);
          }
        else
          {
            // x^{n+1} = x^{n} + f_1 * (x^{n}-x^{n-1})
            //           + f_2 * P^{-1} * (b-A*x^{n})
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (std::size_t i = begin; i < end; ++i)
              // for efficiency reason, write back to temp_vector, which is
              // already modified during vmult (in best case, the modified
              // values are not written back to main memory yet so that
              // we do not have to pay additional costs for writing and
              // read-for-ownership)
              tmp_vector[i] =
                factor1_plus_1 * solution[i] - factor1 * solution_old[i] +
                factor2 * matrix_diagonal_inverse[i] * (rhs[i] - tmp_vector[i]);
          }
      }

      const Number      *rhs;
      const Number      *matrix_diagonal_inverse;
      const unsigned int iteration_index;
      const Number       factor1;
      const Number       factor2;
      mutable Number    *solution_old;
      mutable Number    *tmp_vector;
      mutable Number    *solution;
    };

    template <typename Number>
    struct VectorUpdatesRange : public ::dealii::parallel::ParallelForInteger
    {
      VectorUpdatesRange(const VectorUpdater<Number> &updater,
                         const std::size_t            size)
        : updater(updater)
      {
        if (size < internal::VectorImplementation::minimum_parallel_grain_size)
          VectorUpdatesRange::apply_to_subrange(0, size);
        else
          apply_parallel(
            0,
            size,
            internal::VectorImplementation::minimum_parallel_grain_size);
      }

      ~VectorUpdatesRange() override = default;

      virtual void
      apply_to_subrange(const std::size_t begin,
                        const std::size_t end) const override
      {
        updater.apply_to_subrange(begin, end);
      }

      const VectorUpdater<Number> &updater;
    };

    // selection for diagonal matrix around deal.II vector
    template <typename Number>
    inline void
    vector_updates(
      const ::dealii::Vector<Number>                         &rhs,
      const dealii::DiagonalMatrix<::dealii::Vector<Number>> &jacobi,
      const unsigned int                                      iteration_index,
      const double                                            factor1,
      const double                                            factor2,
      ::dealii::Vector<Number>                               &solution_old,
      ::dealii::Vector<Number>                               &temp_vector1,
      ::dealii::Vector<Number> &,
      ::dealii::Vector<Number> &solution)
    {
      VectorUpdater<Number> upd(rhs.begin(),
                                jacobi.get_vector().begin(),
                                iteration_index,
                                factor1,
                                factor2,
                                solution_old.begin(),
                                temp_vector1.begin(),
                                solution.begin());
      VectorUpdatesRange<Number>(upd, rhs.size());

      // swap vectors x^{n+1}->x^{n}, given the updates in the function above
      if (iteration_index == 0)
        {
          // nothing to do here because we can immediately write into the
          // solution vector without remembering any of the other vectors
        }
      else
        {
          solution.swap(temp_vector1);
          solution_old.swap(temp_vector1);
        }
    }

    // selection for diagonal matrix around parallel deal.II vector
    template <typename Number>
    inline void
    vector_updates(
      const LinearAlgebra::distributed::Vector<Number, MemorySpace::Host> &rhs,
      const dealii::DiagonalMatrix<
        LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>> &jacobi,
      const unsigned int iteration_index,
      const double       factor1,
      const double       factor2,
      LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>
        &solution_old,
      LinearAlgebra::distributed::Vector<Number, MemorySpace::Host>
        &temp_vector1,
      LinearAlgebra::distributed::Vector<Number, MemorySpace::Host> &,
      LinearAlgebra::distributed::Vector<Number, MemorySpace::Host> &solution)
    {
      VectorUpdater<Number> upd(rhs.begin(),
                                jacobi.get_vector().begin(),
                                iteration_index,
                                factor1,
                                factor2,
                                solution_old.begin(),
                                temp_vector1.begin(),
                                solution.begin());
      VectorUpdatesRange<Number>(upd, rhs.locally_owned_size());

      // swap vectors x^{n+1}->x^{n}, given the updates in the function above
      if (iteration_index == 0)
        {
          // nothing to do here because we can immediately write into the
          // solution vector without remembering any of the other vectors
        }
      else
        {
          solution.swap(temp_vector1);
          solution_old.swap(temp_vector1);
        }
    }

    // We need to have a separate declaration for static const members

    // general case and the case that the preconditioner can work on
    // ranges (covered by vector_updates())
    template <
      typename MatrixType,
      typename VectorType,
      typename PreconditionerType,
      std::enable_if_t<
        !has_vmult_with_std_functions<MatrixType,
                                      VectorType,
                                      PreconditionerType> &&
          !(has_vmult_with_std_functions_for_precondition<PreconditionerType,
                                                          VectorType> &&
            has_vmult_with_std_functions_for_precondition<MatrixType,
                                                          VectorType>),
        int> * = nullptr>
    inline void
    vmult_and_update(const MatrixType         &matrix,
                     const PreconditionerType &preconditioner,
                     const VectorType         &rhs,
                     const unsigned int        iteration_index,
                     const double              factor1,
                     const double              factor2,
                     VectorType               &solution,
                     VectorType               &solution_old,
                     VectorType               &temp_vector1,
                     VectorType               &temp_vector2)
    {
      if (iteration_index > 0)
        matrix.vmult(temp_vector1, solution);
      vector_updates(rhs,
                     preconditioner,
                     iteration_index,
                     factor1,
                     factor2,
                     solution_old,
                     temp_vector1,
                     temp_vector2,
                     solution);
    }

    // case that both the operator and the preconditioner can work on
    // subranges
    template <
      typename MatrixType,
      typename VectorType,
      typename PreconditionerType,
      std::enable_if_t<
        !has_vmult_with_std_functions<MatrixType,
                                      VectorType,
                                      PreconditionerType> &&
          (has_vmult_with_std_functions_for_precondition<PreconditionerType,
                                                         VectorType> &&
           has_vmult_with_std_functions_for_precondition<MatrixType,
                                                         VectorType>),
        int> * = nullptr>
    inline void
    vmult_and_update(const MatrixType         &matrix,
                     const PreconditionerType &preconditioner,
                     const VectorType         &rhs,
                     const unsigned int        iteration_index,
                     const double              factor1_,
                     const double              factor2_,
                     VectorType               &solution,
                     VectorType               &solution_old,
                     VectorType               &temp_vector1,
                     VectorType               &temp_vector2)
    {
      using Number = typename VectorType::value_type;

      const Number factor1        = factor1_;
      const Number factor1_plus_1 = 1. + factor1_;
      const Number factor2        = factor2_;

      if (iteration_index == 0)
        {
          preconditioner.vmult(
            solution,
            rhs,
            [&](const unsigned int start_range, const unsigned int end_range) {
              // zero 'solution' before running the vmult operation
              if (end_range > start_range)
                std::memset(solution.begin() + start_range,
                            0,
                            sizeof(Number) * (end_range - start_range));
            },
            [&](const unsigned int start_range, const unsigned int end_range) {
              const auto solution_ptr = solution.begin();

              DEAL_II_OPENMP_SIMD_PRAGMA
              for (std::size_t i = start_range; i < end_range; ++i)
                solution_ptr[i] *= factor2;
            });
        }
      else
        {
          temp_vector1.reinit(rhs, true);
          temp_vector2.reinit(rhs, true);

          // 1) compute residual (including operator application)
          matrix.vmult(
            temp_vector1,
            solution,
            [&](const unsigned int start_range, const unsigned int end_range) {
              // zero 'temp_vector1' before running the vmult
              // operation
              if (end_range > start_range)
                std::memset(temp_vector1.begin() + start_range,
                            0,
                            sizeof(Number) * (end_range - start_range));
            },
            [&](const unsigned int start_range, const unsigned int end_range) {
              const auto rhs_ptr = rhs.begin();
              const auto tmp_ptr = temp_vector1.begin();

              DEAL_II_OPENMP_SIMD_PRAGMA
              for (std::size_t i = start_range; i < end_range; ++i)
                tmp_ptr[i] = rhs_ptr[i] - tmp_ptr[i];
            });

          // 2) perform vector updates (including preconditioner application)
          preconditioner.vmult(
            temp_vector2,
            temp_vector1,
            [&](const unsigned int start_range, const unsigned int end_range) {
              // zero 'temp_vector2' before running the vmult
              // operation
              if (end_range > start_range)
                std::memset(temp_vector2.begin() + start_range,
                            0,
                            sizeof(Number) * (end_range - start_range));
            },
            [&](const unsigned int start_range, const unsigned int end_range) {
              const auto solution_ptr     = solution.begin();
              const auto solution_old_ptr = solution_old.begin();
              const auto tmp_ptr          = temp_vector2.begin();

              if (iteration_index == 1)
                {
                  DEAL_II_OPENMP_SIMD_PRAGMA
                  for (std::size_t i = start_range; i < end_range; ++i)
                    tmp_ptr[i] =
                      factor1_plus_1 * solution_ptr[i] + factor2 * tmp_ptr[i];
                }
              else
                {
                  DEAL_II_OPENMP_SIMD_PRAGMA
                  for (std::size_t i = start_range; i < end_range; ++i)
                    tmp_ptr[i] = factor1_plus_1 * solution_ptr[i] -
                                 factor1 * solution_old_ptr[i] +
                                 factor2 * tmp_ptr[i];
                }
            });

          solution.swap(temp_vector2);
          solution_old.swap(temp_vector2);
        }
    }

    // case that the operator can work on subranges and the preconditioner
    // is a diagonal
    template <typename MatrixType,
              typename VectorType,
              typename PreconditionerType,
              std::enable_if_t<has_vmult_with_std_functions<MatrixType,
                                                            VectorType,
                                                            PreconditionerType>,
                               int> * = nullptr>
    inline void
    vmult_and_update(const MatrixType         &matrix,
                     const PreconditionerType &preconditioner,
                     const VectorType         &rhs,
                     const unsigned int        iteration_index,
                     const double              factor1,
                     const double              factor2,
                     VectorType               &solution,
                     VectorType               &solution_old,
                     VectorType               &temp_vector1,
                     VectorType &)
    {
      using Number = typename VectorType::value_type;
      VectorUpdater<Number> updater(rhs.begin(),
                                    preconditioner.get_vector().begin(),
                                    iteration_index,
                                    factor1,
                                    factor2,
                                    solution_old.begin(),
                                    temp_vector1.begin(),
                                    solution.begin());
      if (iteration_index > 0)
        matrix.vmult(
          temp_vector1,
          solution,
          [&](const unsigned int start_range, const unsigned int end_range) {
            // zero 'temp_vector1' before running the vmult
            // operation
            if (end_range > start_range)
              std::memset(temp_vector1.begin() + start_range,
                          0,
                          sizeof(Number) * (end_range - start_range));
          },
          [&](const unsigned int start_range, const unsigned int end_range) {
            if (end_range > start_range)
              updater.apply_to_subrange(start_range, end_range);
          });
      else
        updater.apply_to_subrange(0U, solution.locally_owned_size());

      // swap vectors x^{n+1}->x^{n}, given the updates in the function above
      if (iteration_index == 0)
        {
          // nothing to do here because we can immediately write into the
          // solution vector without remembering any of the other vectors
        }
      else
        {
          solution.swap(temp_vector1);
          solution_old.swap(temp_vector1);
        }
    }

    template <typename MatrixType, typename PreconditionerType>
    inline void
    initialize_preconditioner(
      const MatrixType & /*matrix*/,
      std::shared_ptr<PreconditionerType> &preconditioner)
    {
      (void)preconditioner;
      AssertThrow(preconditioner.get() != nullptr, ExcNotInitialized());
    }

    template <typename MatrixType, typename VectorType>
    inline void
    initialize_preconditioner(
      const MatrixType                                    &matrix,
      std::shared_ptr<dealii::DiagonalMatrix<VectorType>> &preconditioner)
    {
      if (preconditioner.get() == nullptr || preconditioner->m() != matrix.m())
        {
          if (preconditioner.get() == nullptr)
            preconditioner =
              std::make_shared<dealii::DiagonalMatrix<VectorType>>();

          Assert(
            preconditioner->m() == 0,
            ExcMessage(
              "Preconditioner appears to be initialized but not sized correctly"));

          // This part only works in serial
          if (preconditioner->m() != matrix.m())
            {
              preconditioner->get_vector().reinit(matrix.m());
              for (typename VectorType::size_type i = 0; i < matrix.m(); ++i)
                preconditioner->get_vector()(i) = 1. / matrix.el(i, i);
            }
        }
    }
  } // namespace PreconditionChebyshevImplementation
} // namespace internal



template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::
  AdditionalData::AdditionalData(const unsigned int        degree,
                                 const double              smoothing_range,
                                 const unsigned int        eig_cg_n_iterations,
                                 const double              eig_cg_residual,
                                 const double              max_eigenvalue,
                                 const EigenvalueAlgorithm eigenvalue_algorithm,
                                 const PolynomialType      polynomial_type)
  : internal::EigenvalueAlgorithmAdditionalData<PreconditionerType>(
      smoothing_range,
      eig_cg_n_iterations,
      eig_cg_residual,
      max_eigenvalue,
      eigenvalue_algorithm)
  , degree(degree)
  , polynomial_type(polynomial_type)
{}



template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::
  PreconditionChebyshev()
  : theta(1.)
  , delta(1.)
  , eigenvalues_are_initialized(false)
{
  static_assert(
    std::is_same_v<size_type, typename VectorType::size_type>,
    "PreconditionChebyshev and VectorType must have the same size_type.");
}



template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline void
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::initialize(
  const MatrixType     &matrix,
  const AdditionalData &additional_data)
{
  matrix_ptr = &matrix;
  data       = additional_data;
  Assert(data.degree > 0,
         ExcMessage("The degree of the Chebyshev method must be positive."));
  internal::PreconditionChebyshevImplementation::initialize_preconditioner(
    matrix, data.preconditioner);
  eigenvalues_are_initialized = false;
}



template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline void
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::clear()
{
  eigenvalues_are_initialized = false;
  theta = delta = 1.0;
  matrix_ptr    = nullptr;
  {
    VectorType empty_vector;
    solution_old.reinit(empty_vector);
    temp_vector1.reinit(empty_vector);
    temp_vector2.reinit(empty_vector);
  }
  data.preconditioner.reset();
}



template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline typename internal::EigenvalueInformation
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::
  estimate_eigenvalues(const VectorType &src) const
{
  Assert(eigenvalues_are_initialized == false, ExcInternalError());

  solution_old.reinit(src);
  temp_vector1.reinit(src, true);

  auto info = internal::estimate_eigenvalues<MatrixType>(
    data, matrix_ptr, solution_old, temp_vector1, data.degree);

  const double alpha = (data.smoothing_range > 1. ?
                          info.max_eigenvalue_estimate / data.smoothing_range :
                          std::min(0.9 * info.max_eigenvalue_estimate,
                                   info.min_eigenvalue_estimate));

  // in case the user set the degree to invalid unsigned int, we have to
  // determine the number of necessary iterations from the Chebyshev error
  // estimate, given the target tolerance specified by smoothing_range. This
  // estimate is based on the error formula given in section 5.1 of
  // R. S. Varga, Matrix iterative analysis, 2nd ed., Springer, 2009
  if (data.degree == numbers::invalid_unsigned_int)
    {
      const double actual_range = info.max_eigenvalue_estimate / alpha;
      const double sigma        = (1. - std::sqrt(1. / actual_range)) /
                           (1. + std::sqrt(1. / actual_range));
      const double eps = data.smoothing_range;
      const_cast<
        PreconditionChebyshev<MatrixType, VectorType, PreconditionerType> *>(
        this)
        ->data.degree =
        1 + static_cast<unsigned int>(
              std::log(1. / eps + std::sqrt(1. / eps / eps - 1.)) /
              std::log(1. / sigma));
    }

  info.degree = data.degree;

  const_cast<
    PreconditionChebyshev<MatrixType, VectorType, PreconditionerType> *>(this)
    ->delta =
    (data.polynomial_type == AdditionalData::PolynomialType::fourth_kind) ?
      (info.max_eigenvalue_estimate) :
      ((info.max_eigenvalue_estimate - alpha) * 0.5);
  const_cast<
    PreconditionChebyshev<MatrixType, VectorType, PreconditionerType> *>(this)
    ->theta = (info.max_eigenvalue_estimate + alpha) * 0.5;

  // We do not need the second temporary vector in case we have a
  // DiagonalMatrix as preconditioner and use deal.II's own vectors
  using NumberType = typename VectorType::value_type;
  if (std::is_same_v<PreconditionerType, dealii::DiagonalMatrix<VectorType>> ==
        false ||
      (std::is_same_v<VectorType, dealii::Vector<NumberType>> == false &&
       ((std::is_same_v<
           VectorType,
           LinearAlgebra::distributed::Vector<NumberType, MemorySpace::Host>> ==
         false) ||
        (std::is_same_v<VectorType,
                        LinearAlgebra::distributed::
                          Vector<NumberType, MemorySpace::Default>> == false))))
    temp_vector2.reinit(src, true);
  else
    {
      VectorType empty_vector;
      temp_vector2.reinit(empty_vector);
    }

  const_cast<
    PreconditionChebyshev<MatrixType, VectorType, PreconditionerType> *>(this)
    ->eigenvalues_are_initialized = true;

  return info;
}



template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline void
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::vmult(
  VectorType       &solution,
  const VectorType &rhs) const
{
  std::lock_guard<std::mutex> lock(mutex);
  if (eigenvalues_are_initialized == false)
    estimate_eigenvalues(rhs);

  internal::PreconditionChebyshevImplementation::vmult_and_update(
    *matrix_ptr,
    *data.preconditioner,
    rhs,
    0,
    0.,
    (data.polynomial_type == AdditionalData::PolynomialType::fourth_kind) ?
      (4. / (3. * delta)) :
      (1. / theta),
    solution,
    solution_old,
    temp_vector1,
    temp_vector2);

  // if delta is zero, we do not need to iterate because the updates will be
  // zero
  if (data.degree < 2 || std::abs(delta) < 1e-40)
    return;

  double rhok = delta / theta, sigma = theta / delta;
  for (unsigned int k = 0; k < data.degree - 1; ++k)
    {
      double factor1 = 0.0;
      double factor2 = 0.0;

      if (data.polynomial_type == AdditionalData::PolynomialType::fourth_kind)
        {
          factor1 = (2 * k + 1.) / (2 * k + 5.);
          factor2 = (8 * k + 12.) / (delta * (2 * k + 5.));
        }
      else
        {
          const double rhokp = 1. / (2. * sigma - rhok);
          factor1            = rhokp * rhok;
          factor2            = 2. * rhokp / delta;
          rhok               = rhokp;
        }

      internal::PreconditionChebyshevImplementation::vmult_and_update(
        *matrix_ptr,
        *data.preconditioner,
        rhs,
        k + 1,
        factor1,
        factor2,
        solution,
        solution_old,
        temp_vector1,
        temp_vector2);
    }
}



template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline void
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::Tvmult(
  VectorType       &solution,
  const VectorType &rhs) const
{
  std::lock_guard<std::mutex> lock(mutex);
  if (eigenvalues_are_initialized == false)
    estimate_eigenvalues(rhs);

  internal::PreconditionChebyshevImplementation::vector_updates(
    rhs,
    *data.preconditioner,
    0,
    0.,
    (data.polynomial_type == AdditionalData::PolynomialType::fourth_kind) ?
      (4. / (3. * delta)) :
      (1. / theta),
    solution_old,
    temp_vector1,
    temp_vector2,
    solution);

  if (data.degree < 2 || std::abs(delta) < 1e-40)
    return;

  double rhok = delta / theta, sigma = theta / delta;
  for (unsigned int k = 0; k < data.degree - 1; ++k)
    {
      double factor1 = 0.0;
      double factor2 = 0.0;

      if (data.polynomial_type == AdditionalData::PolynomialType::fourth_kind)
        {
          factor1 = (2 * k + 1.) / (2 * k + 5.);
          factor2 = (8 * k + 12.) / (delta * (2 * k + 5.));
        }
      else
        {
          const double rhokp = 1. / (2. * sigma - rhok);
          factor1            = rhokp * rhok;
          factor2            = 2. * rhokp / delta;
          rhok               = rhokp;
        }

      matrix_ptr->Tvmult(temp_vector1, solution);
      internal::PreconditionChebyshevImplementation::vector_updates(
        rhs,
        *data.preconditioner,
        k + 1,
        factor1,
        factor2,
        solution_old,
        temp_vector1,
        temp_vector2,
        solution);
    }
}



template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline void
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::step(
  VectorType       &solution,
  const VectorType &rhs) const
{
  std::lock_guard<std::mutex> lock(mutex);
  if (eigenvalues_are_initialized == false)
    estimate_eigenvalues(rhs);

  internal::PreconditionChebyshevImplementation::vmult_and_update(
    *matrix_ptr,
    *data.preconditioner,
    rhs,
    1,
    0.,
    (data.polynomial_type == AdditionalData::PolynomialType::fourth_kind) ?
      (4. / (3. * delta)) :
      (1. / theta),
    solution,
    solution_old,
    temp_vector1,
    temp_vector2);

  if (data.degree < 2 || std::abs(delta) < 1e-40)
    return;

  double rhok = delta / theta, sigma = theta / delta;
  for (unsigned int k = 0; k < data.degree - 1; ++k)
    {
      double factor1 = 0.0;
      double factor2 = 0.0;

      if (data.polynomial_type == AdditionalData::PolynomialType::fourth_kind)
        {
          factor1 = (2 * k + 1.) / (2 * k + 5.);
          factor2 = (8 * k + 12.) / (delta * (2 * k + 5.));
        }
      else
        {
          const double rhokp = 1. / (2. * sigma - rhok);
          factor1            = rhokp * rhok;
          factor2            = 2. * rhokp / delta;
          rhok               = rhokp;
        }

      internal::PreconditionChebyshevImplementation::vmult_and_update(
        *matrix_ptr,
        *data.preconditioner,
        rhs,
        k + 2,
        factor1,
        factor2,
        solution,
        solution_old,
        temp_vector1,
        temp_vector2);
    }
}



template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline void
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::Tstep(
  VectorType       &solution,
  const VectorType &rhs) const
{
  std::lock_guard<std::mutex> lock(mutex);
  if (eigenvalues_are_initialized == false)
    estimate_eigenvalues(rhs);

  matrix_ptr->Tvmult(temp_vector1, solution);
  internal::PreconditionChebyshevImplementation::vector_updates(
    rhs,
    *data.preconditioner,
    1,
    0.,
    (data.polynomial_type == AdditionalData::PolynomialType::fourth_kind) ?
      (4. / (3. * delta)) :
      (1. / theta),
    solution_old,
    temp_vector1,
    temp_vector2,
    solution);

  if (data.degree < 2 || std::abs(delta) < 1e-40)
    return;

  double rhok = delta / theta, sigma = theta / delta;
  for (unsigned int k = 0; k < data.degree - 1; ++k)
    {
      double factor1 = 0.0;
      double factor2 = 0.0;

      if (data.polynomial_type == AdditionalData::PolynomialType::fourth_kind)
        {
          factor1 = (2 * k + 1.) / (2 * k + 5.);
          factor2 = (8 * k + 12.) / (delta * (2 * k + 5.));
        }
      else
        {
          const double rhokp = 1. / (2. * sigma - rhok);
          factor1            = rhokp * rhok;
          factor2            = 2. * rhokp / delta;
          rhok               = rhokp;
        }

      matrix_ptr->Tvmult(temp_vector1, solution);
      internal::PreconditionChebyshevImplementation::vector_updates(
        rhs,
        *data.preconditioner,
        k + 2,
        factor1,
        factor2,
        solution_old,
        temp_vector1,
        temp_vector2,
        solution);
    }
}



template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline typename PreconditionChebyshev<MatrixType,
                                      VectorType,
                                      PreconditionerType>::size_type
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::m() const
{
  Assert(matrix_ptr != nullptr, ExcNotInitialized());
  return matrix_ptr->m();
}



template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline typename PreconditionChebyshev<MatrixType,
                                      VectorType,
                                      PreconditionerType>::size_type
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::n() const
{
  Assert(matrix_ptr != nullptr, ExcNotInitialized());
  return matrix_ptr->n();
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

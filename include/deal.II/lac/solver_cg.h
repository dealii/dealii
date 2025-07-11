// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_solver_cg_h
#define dealii_solver_cg_h


#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/tridiagonal_matrix.h>

#include <boost/signals2.hpp>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

// forward declaration
#ifndef DOXYGEN
class PreconditionIdentity;
namespace LinearAlgebra
{
  namespace distributed
  {
    template <typename, typename>
    class Vector;
  }
} // namespace LinearAlgebra
#endif


/** @addtogroup Solvers */
/** @{ */

/**
 * This class implements the preconditioned Conjugate Gradients (CG)
 * method that can be used to solve linear systems with a symmetric positive
 * definite matrix. This
 * class is used first in step-3 and step-4, but is used in many other
 * tutorial programs as well. Like all other solver classes, it can work on
 * any kind of vector and matrix as long as they satisfy certain requirements
 * (for the requirements on matrices and vectors in order to work with this
 * class, see the documentation of the Solver base class). The type of the
 * solution vector must be passed as template argument, and defaults to
 * dealii::Vector<double>. The AdditionalData structure allows to control the
 * type of residual for the stopping condition.
 *
 * @note The CG method requires a symmetric preconditioner (i.e., for example,
 * SOR is not a possible choice). There is a variant of the solver,
 * SolverFlexibleCG, that allows to use a variable preconditioner or a
 * preconditioner with some slight non-symmetry (like weighted Schwarz
 * methods), by using a different formula for the step length in the
 * computation of the next search direction.
 *
 *
 * <h3>Eigenvalue computation</h3>
 *
 * The cg-method performs an orthogonal projection of the original
 * preconditioned linear system to another system of smaller dimension.
 * Furthermore, the projected matrix @p T is tri-diagonal. Since the
 * projection is orthogonal, the eigenvalues of @p T approximate those of the
 * original preconditioned matrix @p PA. In fact, after @p n steps, where @p n
 * is the dimension of the original system, the eigenvalues of both matrices
 * are equal. But, even for small numbers of iteration steps, the condition
 * number of @p T is a good estimate for the one of @p PA.
 *
 * After @p m steps the matrix T_m can be written in terms of the coefficients
 * @p alpha and @p beta as the tri-diagonal matrix with diagonal elements
 * <tt>1/alpha_0</tt>, <tt>1/alpha_1 + beta_0/alpha_0</tt>, ...,
 * <tt>1/alpha_{m-1}+beta_{m-2}/alpha_{m-2}</tt> and off-diagonal elements
 * <tt>sqrt(beta_0)/alpha_0</tt>, ..., <tt>sqrt(beta_{m-2})/alpha_{m-2}</tt>.
 * The eigenvalues of this matrix can be computed by postprocessing.
 *
 * @see Y. Saad: "Iterative methods for Sparse Linear Systems", section 6.7.3
 * for details.
 *
 * The coefficients, eigenvalues and condition number (computed as the ratio
 * of the largest over smallest eigenvalue) can be obtained by connecting a
 * function as a slot to the solver using one of the functions @p
 * connect_coefficients_slot, @p connect_eigenvalues_slot and @p
 * connect_condition_number_slot. These slots will then be called from the
 * solver with the estimates as argument.
 *
 * <h3>Observing the progress of linear solver iterations</h3>
 *
 * The solve() function of this class uses the mechanism described in the
 * Solver base class to determine convergence. This mechanism can also be used
 * to observe the progress of the iteration.
 *
 * <h4>Optimized operations with specific `MatrixType` argument</h4>
 *
 * This class enables to embed the vector updates into the matrix-vector
 * product in case the `MatrixType` and `PreconditionerType` support such a
 * mode of operation. To this end, the `VectorType` needs to be
 * LinearAlgebra::distributed::Vector, the class `MatrixType` needs to provide
 * a function with the signature
 * @code
 * void MatrixType::vmult(
 *    VectorType &,
 *    const VectorType &,
 *    const std::function<void(const unsigned int, const unsigned int)> &,
 *    const std::function<void(const unsigned int, const unsigned int)> &) const
 * @endcode
 * where the two given functions run before and after the matrix-vector
 * product, respectively, and the `PreconditionerType` needs to provide a
 * function either the signature
 * @code
 * Number PreconditionerType::apply(unsigned int index, const Number src) const
 * @endcode
 * to apply the action of the preconditioner on a single element (effectively
 * being a diagonal preconditioner), or the signature
 * @code
 * void PreconditionerType::apply_to_subrange(unsigned int start_range,
 *                                            unsigned int end_range,
 *                                            const Number* src_ptr_to_subrange,
 *                                            Number* dst_ptr_to_subrange)
 * @endcode
 * where the pointers `src_ptr_to_subrange` and `dst_ptr_to_subrange` point to
 * the location in the vector where the operation should be applied to. If both
 * functions are given, the more optimized `apply` path is selected. The
 * functions passed to `MatrixType::vmult` take as arguments a sub-range among
 * the locally owned elements of the vector, defined as half-open
 * intervals. The intervals are designed to be scheduled close to the time the
 * matrix-vector product touches those entries in the `src` and `dst` vectors,
 * respectively, with the requirement that
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
 *
 * In terms of the SolverCG implementation, the operation before the loop will
 * run the updates on the vectors according to a variant presented in
 * Algorithm 2.2 of @cite Chronopoulos1989 (but for a preconditioner), whereas
 * the operation after the loop performs a total of 7 reductions in parallel.
 *
 * <h3>Preconditioned residual</h3>
 *
 * @p AdditionalData allows you to choose between using the explicit
 * or implicit residual as a stopping condition for the iterative
 * solver. This behavior can be overridden by using the flag
 * AdditionalData::use_default_residual. A <tt>true</tt> value refers to the
 * implicit residual, while <tt>false</tt> reverts
 * it. The former uses the result of the matrix-vector
 * product already computed in other algorithm steps to derive the residual by
 * a mere vector update, whereas the latter explicitly calculates the system
 * residual with an additional matrix-vector product.
 * More information on explicit and implicit residual stopping
 * criteria can be found
 * <a
 * href="https://en.wikipedia.org/wiki/Conjugate_gradient_method#Explicit_residual_calculation">link
 * here</a>.
 */
template <typename VectorType = Vector<double>>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
class SolverCG : public SolverBase<VectorType>
{
public:
  /**
   * Declare type for container size.
   */
  using size_type = types::global_dof_index;


  /**
   * Standardized data struct to pipe additional data to the solver.
   */
  struct AdditionalData
  {
    /**
     * Constructor. By default, set the residual of the stopping criterion
     * to the implicit residual. A <tt>true</tt> value of
     * AdditionalData::use_default_residual refers to the
     * implicit residual, while <tt>false</tt> reverts
     * it. The former uses the result of the matrix-vector
     * product already computed in other algorithm steps to derive the residual
     * by a mere vector update, whereas the latter explicitly calculates the
     * system residual with an additional matrix-vector product. More
     * information on explicit and implicit residual stopping criteria can be
     * found <a
     * href="https://en.wikipedia.org/wiki/Conjugate_gradient_method#Explicit_residual_calculation">link
     * here</a>.
     */
    explicit AdditionalData(const bool use_default_residual = true);

    /**
     * Flag for the default residual that is used to measure convergence.
     */
    bool use_default_residual;
  };


  /**
   * Constructor.
   */
  SolverCG(SolverControl            &cn,
           VectorMemory<VectorType> &mem,
           const AdditionalData     &data = AdditionalData());

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  SolverCG(SolverControl &cn, const AdditionalData &data = AdditionalData());

  /**
   * Virtual destructor.
   */
  virtual ~SolverCG() override = default;

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
   * Connect a slot to retrieve the CG coefficients. The slot will be called
   * with alpha as the first argument and with beta as the second argument,
   * where alpha and beta follow the notation in Y. Saad: "Iterative methods
   * for Sparse Linear Systems", section 6.7. Called once per iteration
   */
  boost::signals2::connection
  connect_coefficients_slot(
    const std::function<void(typename VectorType::value_type,
                             typename VectorType::value_type)> &slot);

  /**
   * Connect a slot to retrieve the estimated condition number. Called on each
   * iteration if every_iteration=true, otherwise called once when iterations
   * are ended (i.e., either because convergence has been achieved, or because
   * divergence has been detected).
   */
  boost::signals2::connection
  connect_condition_number_slot(const std::function<void(double)> &slot,
                                const bool every_iteration = false);

  /**
   * Connect a slot to retrieve the estimated eigenvalues. Called on each
   * iteration if every_iteration=true, otherwise called once when iterations
   * are ended (i.e., either because convergence has been achieved, or because
   * divergence has been detected).
   */
  boost::signals2::connection
  connect_eigenvalues_slot(
    const std::function<void(const std::vector<double> &)> &slot,
    const bool every_iteration = false);

protected:
  /**
   * Interface for derived class. This function gets the current iteration
   * vector, the residual and the update vector in each step. It can be used
   * for graphical output of the convergence history.
   */
  virtual void
  print_vectors(const unsigned int step,
                const VectorType  &x,
                const VectorType  &r,
                const VectorType  &d) const;

  /**
   * Estimates the eigenvalues from diagonal and offdiagonal. Uses these
   * estimate to compute the condition number. Calls the signals
   * eigenvalues_signal and cond_signal with these estimates as arguments.
   */
  static void
  compute_eigs_and_cond(
    const std::vector<typename VectorType::value_type> &diagonal,
    const std::vector<typename VectorType::value_type> &offdiagonal,
    const boost::signals2::signal<void(const std::vector<double> &)>
                                                &eigenvalues_signal,
    const boost::signals2::signal<void(double)> &cond_signal);

  /**
   * Additional parameters.
   */
  AdditionalData additional_data;

  /**
   * Signal used to retrieve the CG coefficients. Called on each iteration.
   */
  boost::signals2::signal<void(typename VectorType::value_type,
                               typename VectorType::value_type)>
    coefficients_signal;

  /**
   * Signal used to retrieve the estimated condition number. Called once when
   * all iterations are ended.
   */
  boost::signals2::signal<void(double)> condition_number_signal;

  /**
   * Signal used to retrieve the estimated condition numbers. Called on each
   * iteration.
   */
  boost::signals2::signal<void(double)> all_condition_numbers_signal;

  /**
   * Signal used to retrieve the estimated eigenvalues. Called once when all
   * iterations are ended.
   */
  boost::signals2::signal<void(const std::vector<double> &)> eigenvalues_signal;

  /**
   * Signal used to retrieve the estimated eigenvalues. Called on each
   * iteration.
   */
  boost::signals2::signal<void(const std::vector<double> &)>
    all_eigenvalues_signal;

  /**
   * Flag to indicate whether the classical Fletcher--Reeves update formula
   * for the parameter $\beta_k$ (standard CG algorithm, minimal storage
   * needs) or the flexible conjugate gradient method with Polak-Ribiere
   * formula for $\beta_k$ should be used. This base class implementation of
   * SolverCG will always use the former method, whereas the derived class
   * SolverFlexibleCG will use the latter.
   */
  bool determine_beta_by_flexible_formula;
};



/**
 * This class implements a flexible variant of the conjugate gradient method,
 * which is based on a different formula to compute $\beta_k$ in the process
 * of constructing a new search direction that is A-orthogonal against the
 * previous one. Rather than using the Fletcher--Reeves update formula with
 * $\beta_k = \frac{\mathbf{r}^T_{k+1} \mathbf{z}_{k+1}}{\mathbf{r}^T_{k}
 * \mathbf{z}_{k}}$ for computing the new search direction (here
 * $\mathbf{r}_{k+1}$ is the residual in step $k+1$ and $\mathbf{z}_{k+1} =
 * P^{-1} \mathbf{r}_{k+1}$) as in the classical conjugate gradient algorithm,
 * this class selects the Polak-Ribiere formula $\beta_k =
 * \frac{\mathbf{r}^T_{k+1} \left(\mathbf{z}_{k+1} -
 * \mathbf{z}_{k}\right)}{\mathbf{r}^T_{k} \mathbf{z}_{k}}$. The
 * additional term $\mathbf{r}^T_{k+1} \mathbf{z}_{k}$ is zero for linear
 * symmetric-positive definite preconditioners due to the construction of the
 * search directions, so the behavior of SolverFlexibleCG is equivalent for
 * those kinds of situations and merely increases costs by requiring an
 * additional stored vector and associated vector operations. While there are
 * no theoretical guarantees for convergence as in the classical CG algorithm,
 * the current class has been documented to be much more robust for variable
 * preconditioners (e.g., involving some iterative inverse that is not fully
 * converged) or a preconditioner with some slight non-symmetry (like weighted
 * Schwarz methods), which results from the local optimality of the search
 * direction with at least as good progress as the locally optimal steepest
 * descent method.
 */
template <typename VectorType = Vector<double>>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
class SolverFlexibleCG : public SolverCG<VectorType>
{
public:
  /**
   * Declare type for container size.
   */
  using size_type = types::global_dof_index;


  /**
   * Standardized data struct to pipe additional data to the solver.
   */
  struct AdditionalData
  {
    /**
     * Constructor. By default, set the residual of the stopping criterion
     * to the implicit residual. A <tt>true</tt> value of
     * AdditionalData::use_default_residual refers to the
     * implicit residual, while <tt>false</tt> reverts
     * it. The former uses the result of the matrix-vector
     * product already computed in other algorithm steps to derive the residual
     * by a mere vector update, whereas the latter explicitly calculates the
     * system residual with an additional matrix-vector product. More
     * information on explicit and implicit residual stopping criteria can be
     * found <a
     * href="https://en.wikipedia.org/wiki/Conjugate_gradient_method#Explicit_residual_calculation">link
     * here</a>.
     */
    explicit AdditionalData(const bool use_default_residual = true);

    /**
     * Flag for the default residual that is used to measure convergence.
     */
    bool use_default_residual;
  };


  /**
   * Constructor.
   */
  SolverFlexibleCG(SolverControl            &cn,
                   VectorMemory<VectorType> &mem,
                   const AdditionalData     &data = AdditionalData());

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  SolverFlexibleCG(SolverControl        &cn,
                   const AdditionalData &data = AdditionalData());

protected:
  /**
   * Additional parameters.
   */
  AdditionalData additional_data;
};


/** @} */

/*------------------------- Implementation ----------------------------*/

#ifndef DOXYGEN


template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
inline SolverCG<VectorType>::AdditionalData::AdditionalData(
  const bool use_default_residual)
  : use_default_residual(use_default_residual)
{}

template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
SolverCG<VectorType>::SolverCG(SolverControl            &cn,
                               VectorMemory<VectorType> &mem,
                               const AdditionalData     &data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
  , determine_beta_by_flexible_formula(false)
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
SolverCG<VectorType>::SolverCG(SolverControl &cn, const AdditionalData &data)
  : SolverBase<VectorType>(cn)
  , additional_data(data)
  , determine_beta_by_flexible_formula(false)
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
void SolverCG<VectorType>::print_vectors(const unsigned int,
                                         const VectorType &,
                                         const VectorType &,
                                         const VectorType &) const
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
inline void SolverCG<VectorType>::compute_eigs_and_cond(
  const std::vector<typename VectorType::value_type> &diagonal,
  const std::vector<typename VectorType::value_type> &offdiagonal,
  const boost::signals2::signal<void(const std::vector<double> &)>
                                              &eigenvalues_signal,
  const boost::signals2::signal<void(double)> &cond_signal)
{
  // Avoid computing eigenvalues unless they are needed.
  if (!cond_signal.empty() || !eigenvalues_signal.empty())
    {
      TridiagonalMatrix<typename VectorType::value_type> T(diagonal.size(),
                                                           true);
      for (size_type i = 0; i < diagonal.size(); ++i)
        {
          T(i, i) = diagonal[i];
          if (i < diagonal.size() - 1)
            T(i, i + 1) = offdiagonal[i];
        }
      T.compute_eigenvalues();
      // Need two eigenvalues to estimate the condition number.
      if (diagonal.size() > 1)
        {
          auto condition_number = T.eigenvalue(T.n() - 1) / T.eigenvalue(0);
          // Condition number is real valued and nonnegative; simply take
          // the absolute value:
          cond_signal(std::abs(condition_number));
        }
      // Avoid copying the eigenvalues of T to a vector unless a signal is
      // connected.
      if (!eigenvalues_signal.empty())
        {
          std::vector<double> eigenvalues(T.n());
          for (unsigned int j = 0; j < T.n(); ++j)
            {
              // for a hermitian matrix, all eigenvalues are real-valued
              // and non-negative, simply return the absolute value:
              eigenvalues[j] = std::abs(T.eigenvalue(j));
            }
          eigenvalues_signal(eigenvalues);
        }
    }
}



namespace internal
{

  namespace SolverCG
  {
    // This base class is used to select different variants of the conjugate
    // gradient solver. The default variant is used for standard matrix and
    // preconditioner arguments, as provided by the derived class
    // IterationWork below, but there is also a specialized variant further
    // down that uses SFINAE to identify whether matrices and preconditioners
    // support special operations on sub-ranges of the vectors.
    template <typename VectorType,
              typename MatrixType,
              typename PreconditionerType>
    struct IterationWorkerBase
    {
      using Number = typename VectorType::value_type;

      const MatrixType         &A;
      const PreconditionerType &preconditioner;
      const bool                flexible;
      VectorType               &x;
      const VectorType         &b;

      typename VectorMemory<VectorType>::Pointer r_pointer;
      typename VectorMemory<VectorType>::Pointer p_pointer;
      typename VectorMemory<VectorType>::Pointer v_pointer;
      typename VectorMemory<VectorType>::Pointer z_pointer;
      typename VectorMemory<VectorType>::Pointer explicit_r_pointer;

      // Define some aliases for simpler access, using the variables 'r' for
      // the residual b - A*x, 'p' for the search direction, and 'v' for the
      // auxiliary vector. This naming convention is used e.g. by the
      // description on
      // https://en.wikipedia.org/wiki/Conjugate_gradient_method. The variable
      // 'z' gets only used for the flexible variant of the CG method.
      VectorType &r;
      VectorType &p;
      VectorType &v;
      VectorType &z;
      VectorType &explicit_r;

      Number     r_dot_preconditioner_dot_r;
      Number     alpha;
      Number     beta;
      double     residual_norm;
      Number     previous_alpha;
      const bool use_default_residual;


      IterationWorkerBase(const MatrixType         &A,
                          const PreconditionerType &preconditioner,
                          const bool                flexible,
                          VectorMemory<VectorType> &memory,
                          VectorType               &x,
                          const VectorType         &b,
                          const bool               &use_default_residual)
        : A(A)
        , preconditioner(preconditioner)
        , flexible(flexible)
        , x(x)
        , b(b)
        , r_pointer(memory)
        , p_pointer(memory)
        , v_pointer(memory)
        , z_pointer(memory)
        , explicit_r_pointer(memory)
        , r(*r_pointer)
        , p(*p_pointer)
        , v(*v_pointer)
        , z(*z_pointer)
        , explicit_r(*explicit_r_pointer)
        , r_dot_preconditioner_dot_r(Number())
        , alpha(Number())
        , beta(Number())
        , residual_norm(0.0)
        , previous_alpha(Number())
        , use_default_residual(use_default_residual)
      {}

      void
      startup()
      {
        // Initialize without setting the vector entries, as those would soon
        // be overwritten anyway
        r.reinit(x, true);
        p.reinit(x, true);
        v.reinit(x, true);
        if (flexible)
          z.reinit(x, true);
        if (!use_default_residual)
          explicit_r.reinit(x, true);

        // compute residual. if vector is zero, then short-circuit the full
        // computation
        if (!x.all_zero())
          {
            A.vmult(r, x);
            r.sadd(-1., 1., b);
          }
        else
          r.equ(1., b);

        residual_norm = r.l2_norm();
      }
    };



    // Implementation of a conjugate gradient operation with matrices and
    // preconditioners without special capabilities
    template <typename VectorType,
              typename MatrixType,
              typename PreconditionerType,
              typename = int>
    struct IterationWorker
      : public IterationWorkerBase<VectorType, MatrixType, PreconditionerType>
    {
      using BaseClass =
        IterationWorkerBase<VectorType, MatrixType, PreconditionerType>;


      IterationWorker(const MatrixType         &A,
                      const PreconditionerType &preconditioner,
                      const bool                flexible,
                      VectorMemory<VectorType> &memory,
                      VectorType               &x,
                      const VectorType         &b,
                      const bool               &use_default_residual)
        : BaseClass(A,
                    preconditioner,
                    flexible,
                    memory,
                    x,
                    b,
                    use_default_residual)
      {}

      using BaseClass::A;
      using BaseClass::alpha;
      using BaseClass::b;
      using BaseClass::beta;
      using BaseClass::explicit_r;
      using BaseClass::p;
      using BaseClass::preconditioner;
      using BaseClass::r;
      using BaseClass::r_dot_preconditioner_dot_r;
      using BaseClass::residual_norm;
      using BaseClass::use_default_residual;
      using BaseClass::v;
      using BaseClass::x;
      using BaseClass::z;

      void
      do_iteration(const unsigned int iteration_index)
      {
        using Number = typename VectorType::value_type;

        const Number previous_r_dot_preconditioner_dot_r =
          r_dot_preconditioner_dot_r;

        if (std::is_same_v<PreconditionerType, PreconditionIdentity> == false)
          {
            preconditioner.vmult(v, r);
            r_dot_preconditioner_dot_r = r * v;
          }
        else
          r_dot_preconditioner_dot_r = residual_norm * residual_norm;

        const VectorType &direction =
          std::is_same_v<PreconditionerType, PreconditionIdentity> ? r : v;

        if (iteration_index > 1)
          {
            Assert(std::abs(previous_r_dot_preconditioner_dot_r) != 0.,
                   ExcDivideByZero());
            beta =
              r_dot_preconditioner_dot_r / previous_r_dot_preconditioner_dot_r;
            if (this->flexible)
              beta -= (r * z) / previous_r_dot_preconditioner_dot_r;
            p.sadd(beta, 1., direction);
          }
        else
          p.equ(1., direction);

        if (this->flexible)
          z.swap(v);

        A.vmult(v, p);

        const Number p_dot_A_dot_p = p * v;
        Assert(std::abs(p_dot_A_dot_p) != 0., ExcDivideByZero());

        this->previous_alpha = alpha;
        alpha                = r_dot_preconditioner_dot_r / p_dot_A_dot_p;

        x.add(alpha, p);

        // compute the residual norm with implicit residual
        if (use_default_residual)
          {
            residual_norm = std::sqrt(std::abs(r.add_and_dot(-alpha, v, r)));
          }
        // compute the residual norm with the explicit residual, i.e.
        // compute l2 norm of Ax - b.
        else
          {
            // compute the residual conjugate gradient update
            r.add(-alpha, v);
            // compute explicit residual
            A.vmult(explicit_r, x);
            explicit_r.add(-1, b);
            residual_norm = explicit_r.l2_norm();
          }
      }

      void
      finalize_after_convergence(const unsigned int)
      {}
    };


    // In the following, we provide a specialization of the above
    // IterationWorker class that picks up particular features in the matrix
    // and preconditioners.

    // a helper type-trait that leverage SFINAE to figure out if MatrixType has
    // ... MatrixType::vmult(VectorType &, const VectorType&,
    // std::function<...>, std::function<...>) const
    template <typename MatrixType, typename VectorType>
    using vmult_functions_t = decltype(std::declval<const MatrixType>().vmult(
      std::declval<VectorType &>(),
      std::declval<const VectorType &>(),
      std::declval<
        const std::function<void(const unsigned int, const unsigned int)> &>(),
      std::declval<const std::function<void(const unsigned int,
                                            const unsigned int)> &>()));

    template <typename MatrixType, typename VectorType>
    constexpr bool has_vmult_functions =
      is_supported_operation<vmult_functions_t, MatrixType, VectorType>;

    // a helper type-trait that leverage SFINAE to figure out if
    // PreconditionerType has ... PreconditionerType::apply_to_subrange(const
    // unsigned int, const unsigned int, const Number*, Number*) const
    template <typename PreconditionerType>
    using apply_to_subrange_t =
      decltype(std::declval<const PreconditionerType>()
                 .apply_to_subrange(0U, 0U, nullptr, nullptr));

    template <typename PreconditionerType>
    constexpr bool has_apply_to_subrange =
      is_supported_operation<apply_to_subrange_t, PreconditionerType>;

    // a helper type-trait that leverage SFINAE to figure out if
    // PreconditionerType has ... PreconditionerType::apply(const
    // unsigned int, const Number) const
    template <typename PreconditionerType>
    using apply_t =
      decltype(std::declval<const PreconditionerType>().apply(0U, 0.0));

    template <typename PreconditionerType>
    constexpr bool has_apply =
      is_supported_operation<apply_t, PreconditionerType>;


    // Internal function to run one iteration of the conjugate gradient solver
    // for matrices and preconditioners that support interleaving the vector
    // updates with the matrix-vector product.
    template <typename VectorType,
              typename MatrixType,
              typename PreconditionerType>
    struct IterationWorker<
      VectorType,
      MatrixType,
      PreconditionerType,
      std::enable_if_t<has_vmult_functions<MatrixType, VectorType> &&
                         (has_apply_to_subrange<PreconditionerType> ||
                          has_apply<PreconditionerType>)&&std::
                           is_same_v<VectorType,
                                     LinearAlgebra::distributed::Vector<
                                       typename VectorType::value_type,
                                       MemorySpace::Host>>,
                       int>>
      : public IterationWorkerBase<VectorType, MatrixType, PreconditionerType>
    {
      using Number = typename VectorType::value_type;

      Number next_r_dot_preconditioner_dot_r;
      Number previous_beta;

      IterationWorker(const MatrixType         &A,
                      const PreconditionerType &preconditioner,
                      const bool                flexible,
                      VectorMemory<VectorType> &memory,
                      VectorType               &x,
                      const VectorType         &b,
                      const bool               &use_default_residual)
        : IterationWorkerBase<VectorType, MatrixType, PreconditionerType>(
            A,
            preconditioner,
            flexible,
            memory,
            x,
            b,
            use_default_residual)
        , next_r_dot_preconditioner_dot_r(0.)
        , previous_beta(0.)
      {}

      // This is the main iteration function, that will use some of the
      // specialized functions below
      void
      do_iteration(const unsigned int iteration_index)
      {
        if (iteration_index > 1)
          {
            previous_beta = this->beta;
            this->beta    = next_r_dot_preconditioner_dot_r /
                         this->r_dot_preconditioner_dot_r;
          }

        std::array<VectorizedArray<Number>, 7> vectorized_sums = {};

        this->A.vmult(
          this->v,
          this->p,
          [&](const unsigned int begin, const unsigned int end) {
            operation_before_loop(iteration_index, begin, end);
          },
          [&](const unsigned int begin, const unsigned int end) {
            operation_after_loop(begin, end, vectorized_sums);
          });

        std::array<Number, 7> scalar_sums;
        for (unsigned int i = 0; i < 7; ++i)
          scalar_sums[i] = vectorized_sums[i][0];
        for (unsigned int l = 1; l < VectorizedArray<Number>::size(); ++l)
          for (unsigned int i = 0; i < 7; ++i)
            scalar_sums[i] += vectorized_sums[i][l];

        Utilities::MPI::sum(dealii::ArrayView<const Number>(scalar_sums.data(),
                                                            7),
                            this->r.get_mpi_communicator(),
                            dealii::ArrayView<Number>(scalar_sums.data(), 7));

        this->r_dot_preconditioner_dot_r = scalar_sums[6];

        const Number p_dot_A_dot_p = scalar_sums[0];
        Assert(std::abs(p_dot_A_dot_p) != 0., ExcDivideByZero());

        this->previous_alpha = this->alpha;
        this->alpha          = this->r_dot_preconditioner_dot_r / p_dot_A_dot_p;

        // Round-off errors near zero might yield negative values, so take
        // the absolute value in the next two formulas
        this->residual_norm = std::sqrt(std::abs(
          scalar_sums[3] +
          this->alpha * (-2. * scalar_sums[2] + this->alpha * scalar_sums[1])));

        next_r_dot_preconditioner_dot_r = std::abs(
          this->r_dot_preconditioner_dot_r +
          this->alpha * (-2. * scalar_sums[4] + this->alpha * scalar_sums[5]));
      }

      // Function that we use if the PreconditionerType implements an apply()
      // function
      template <typename U = void>
      std::enable_if_t<has_apply<PreconditionerType>, U>
      operation_before_loop(const unsigned int iteration_index,
                            const unsigned int start_range,
                            const unsigned int end_range) const
      {
        Number                *x       = this->x.begin();
        Number                *r       = this->r.begin();
        Number                *p       = this->p.begin();
        Number                *v       = this->v.begin();
        const Number           alpha   = this->alpha;
        const Number           beta    = this->beta;
        constexpr unsigned int n_lanes = VectorizedArray<Number>::size();
        const unsigned int     end_regular =
          start_range + (end_range - start_range) / n_lanes * n_lanes;
        if (iteration_index == 1)
          {
            // Vectorize by hand since compilers are often pretty bad at
            // doing these steps automatically even with
            // DEAL_II_OPENMP_SIMD_PRAGMA
            for (unsigned int j = start_range; j < end_regular; j += n_lanes)
              {
                VectorizedArray<Number> rj, pj;
                rj.load(r + j);
                DEAL_II_OPENMP_SIMD_PRAGMA
                for (unsigned int l = 0; l < n_lanes; ++l)
                  pj[l] = this->preconditioner.apply(j + l, rj[l]);
                pj.store(p + j);
                rj = VectorizedArray<Number>();
                rj.store(v + j);
              }
            for (unsigned int j = end_regular; j < end_range; ++j)
              {
                p[j] = this->preconditioner.apply(j, r[j]);
                v[j] = Number();
              }
          }
        else if (iteration_index % 2 == 0 && beta != Number())
          {
            for (unsigned int j = start_range; j < end_regular; j += n_lanes)
              {
                VectorizedArray<Number> rj, vj, pj, prec_rj;
                rj.load(r + j);
                vj.load(v + j);
                rj -= alpha * vj;
                rj.store(r + j);
                DEAL_II_OPENMP_SIMD_PRAGMA
                for (unsigned int l = 0; l < n_lanes; ++l)
                  prec_rj[l] = this->preconditioner.apply(j + l, rj[l]);
                pj.load(p + j);
                pj = beta * pj + prec_rj;
                pj.store(p + j);
                rj = VectorizedArray<Number>();
                rj.store(v + j);
              }
            for (unsigned int j = end_regular; j < end_range; ++j)
              {
                r[j] -= alpha * v[j];
                p[j] = beta * p[j] + this->preconditioner.apply(j, r[j]);
                v[j] = Number();
              }
          }
        else if (iteration_index % 2 == 0 && beta == Number())
          {
            // Case where beta is zero: we cannot reconstruct p_{j-1} in the
            // next iteration, and must hence compute x here. This can happen
            // before termination
            for (unsigned int j = start_range; j < end_regular; j += n_lanes)
              {
                VectorizedArray<Number> rj, xj, vj, pj, prec_rj;
                rj.load(r + j);
                vj.load(v + j);
                xj.load(x + j);
                pj.load(p + j);
                xj += alpha * pj;
                xj.store(x + j);
                rj -= alpha * vj;
                rj.store(r + j);
                DEAL_II_OPENMP_SIMD_PRAGMA
                for (unsigned int l = 0; l < n_lanes; ++l)
                  prec_rj[l] = this->preconditioner.apply(j + l, rj[l]);
                prec_rj.store(p + j);
                rj = VectorizedArray<Number>();
                rj.store(v + j);
              }
            for (unsigned int j = end_regular; j < end_range; ++j)
              {
                r[j] -= alpha * v[j];
                x[j] += alpha * p[j];
                p[j] = this->preconditioner.apply(j, r[j]);
                v[j] = Number();
              }
          }
        else
          {
            const Number alpha_plus_previous_alpha_over_beta =
              alpha + this->previous_alpha / this->previous_beta;
            const Number previous_alpha_over_beta =
              this->previous_alpha / this->previous_beta;
            for (unsigned int j = start_range; j < end_regular; j += n_lanes)
              {
                VectorizedArray<Number> rj, vj, pj, xj, prec_rj, prec_vj;
                xj.load(x + j);
                pj.load(p + j);
                xj += alpha_plus_previous_alpha_over_beta * pj;
                rj.load(r + j);
                vj.load(v + j);
                DEAL_II_OPENMP_SIMD_PRAGMA
                for (unsigned int l = 0; l < n_lanes; ++l)
                  {
                    prec_rj[l] = this->preconditioner.apply(j + l, rj[l]);
                    prec_vj[l] = this->preconditioner.apply(j + l, vj[l]);
                  }
                xj -= previous_alpha_over_beta * prec_rj;
                xj.store(x + j);
                rj -= alpha * vj;
                rj.store(r + j);
                prec_rj -= alpha * prec_vj;
                pj = beta * pj + prec_rj;
                pj.store(p + j);
                rj = VectorizedArray<Number>();
                rj.store(v + j);
              }
            for (unsigned int j = end_regular; j < end_range; ++j)
              {
                x[j] += alpha_plus_previous_alpha_over_beta * p[j];
                x[j] -= previous_alpha_over_beta *
                        this->preconditioner.apply(j, r[j]);
                r[j] -= alpha * v[j];
                p[j] = beta * p[j] + this->preconditioner.apply(j, r[j]);
                v[j] = Number();
              }
          }
      }

      // Function that we use if the PreconditionerType implements an apply()
      // function
      template <typename U = void>
      std::enable_if_t<has_apply<PreconditionerType>, U>
      operation_after_loop(
        const unsigned int                      start_range,
        const unsigned int                      end_range,
        std::array<VectorizedArray<Number>, 7> &vectorized_sums) const
      {
        const Number                          *r       = this->r.begin();
        const Number                          *p       = this->p.begin();
        const Number                          *v       = this->v.begin();
        std::array<VectorizedArray<Number>, 7> my_sums = {};
        constexpr unsigned int n_lanes = VectorizedArray<Number>::size();
        const unsigned int     end_regular =
          start_range + (end_range - start_range) / n_lanes * n_lanes;
        for (unsigned int j = start_range; j < end_regular; j += n_lanes)
          {
            VectorizedArray<Number> pj, vj, rj, prec_vj, prec_rj;
            pj.load(p + j);
            vj.load(v + j);
            rj.load(r + j);
            DEAL_II_OPENMP_SIMD_PRAGMA
            for (unsigned int l = 0; l < n_lanes; ++l)
              {
                prec_vj[l] = this->preconditioner.apply(j + l, vj[l]);
                prec_rj[l] = this->preconditioner.apply(j + l, rj[l]);
              }
            my_sums[0] += pj * vj;
            my_sums[1] += vj * vj;
            my_sums[2] += rj * vj;
            my_sums[3] += rj * rj;
            my_sums[4] += rj * prec_vj;
            my_sums[5] += vj * prec_vj;
            my_sums[6] += rj * prec_rj;
          }
        for (unsigned int j = end_regular; j < end_range; ++j)
          {
            const Number prec_v = this->preconditioner.apply(j, v[j]);
            const Number prec_r = this->preconditioner.apply(j, r[j]);
            my_sums[0][0] += p[j] * v[j];
            my_sums[1][0] += v[j] * v[j];
            my_sums[2][0] += r[j] * v[j];
            my_sums[3][0] += r[j] * r[j];
            my_sums[4][0] += r[j] * prec_v;
            my_sums[5][0] += v[j] * prec_v;
            my_sums[6][0] += r[j] * prec_r;
          }
        for (unsigned int i = 0; i < vectorized_sums.size(); ++i)
          vectorized_sums[i] += my_sums[i];
      }

      // Function that we use if the PreconditionerType implements an apply()
      // function
      template <typename U = void>
      std::enable_if_t<has_apply<PreconditionerType>, U>
      finalize_after_convergence(const unsigned int iteration_index)
      {
        if (iteration_index % 2 == 1 || this->beta == Number())
          this->x.add(this->alpha, this->p);
        else
          {
            using Number                 = typename VectorType::value_type;
            const unsigned int end_range = this->x.locally_owned_size();

            Number *x = this->x.begin();
            Number *r = this->r.begin();
            Number *p = this->p.begin();

            // Note that we use 'beta' here rather than 'previous_beta' in the
            // formula above, which is because the shift in beta ->
            // previous_beta has not been applied at this stage, allowing us
            // to recover the previous search direction
            const Number alpha_plus_previous_alpha_over_beta =
              this->alpha + this->previous_alpha / this->beta;
            const Number previous_alpha_over_beta =
              this->previous_alpha / this->beta;

            DEAL_II_OPENMP_SIMD_PRAGMA
            for (unsigned int j = 0; j < end_range; ++j)
              {
                x[j] += alpha_plus_previous_alpha_over_beta * p[j] -
                        previous_alpha_over_beta *
                          this->preconditioner.apply(j, r[j]);
              }
          }
      }

      // Function that we use if the PreconditionerType does not implement an
      // apply() function, where we instead need to choose the
      // apply_to_subrange function
      template <typename U = void>
      std::enable_if_t<!has_apply<PreconditionerType>, U>
      operation_before_loop(const unsigned int iteration_index,
                            const unsigned int start_range,
                            const unsigned int end_range) const
      {
        Number                        *x     = this->x.begin() + start_range;
        Number                        *r     = this->r.begin() + start_range;
        Number                        *p     = this->p.begin() + start_range;
        Number                        *v     = this->v.begin() + start_range;
        const Number                   alpha = this->alpha;
        const Number                   beta  = this->beta;
        constexpr unsigned int         grain_size = 128;
        std::array<Number, grain_size> prec_r;
        if (iteration_index == 1)
          {
            for (unsigned int j = start_range; j < end_range; j += grain_size)
              {
                const unsigned int length = std::min(grain_size, end_range - j);
                this->preconditioner.apply_to_subrange(j,
                                                       j + length,
                                                       r,
                                                       prec_r.data());
                DEAL_II_OPENMP_SIMD_PRAGMA
                for (unsigned int i = 0; i < length; ++i)
                  {
                    p[i] = prec_r[i];
                    v[i] = Number();
                  }
                p += length;
                r += length;
                v += length;
              }
          }
        else if (iteration_index % 2 == 0 && beta != Number())
          {
            for (unsigned int j = start_range; j < end_range; j += grain_size)
              {
                const unsigned int length = std::min(grain_size, end_range - j);
                DEAL_II_OPENMP_SIMD_PRAGMA
                for (unsigned int i = 0; i < length; ++i)
                  r[i] -= this->alpha * v[i];
                this->preconditioner.apply_to_subrange(j,
                                                       j + length,
                                                       r,
                                                       prec_r.data());
                DEAL_II_OPENMP_SIMD_PRAGMA
                for (unsigned int i = 0; i < length; ++i)
                  {
                    p[i] = this->beta * p[i] + prec_r[i];
                    v[i] = Number();
                  }
                p += length;
                r += length;
                v += length;
              }
          }
        else if (iteration_index % 2 == 0 && beta == Number())
          {
            // Case where beta is zero: We cannot reconstruct p_{j-1} in the
            // next iteration, and must hence compute x here. This can happen
            // before termination
            for (unsigned int j = start_range; j < end_range; j += grain_size)
              {
                const unsigned int length = std::min(grain_size, end_range - j);
                DEAL_II_OPENMP_SIMD_PRAGMA
                for (unsigned int i = 0; i < length; ++i)
                  r[i] -= this->alpha * v[i];
                this->preconditioner.apply_to_subrange(j,
                                                       j + length,
                                                       r,
                                                       prec_r.data());
                DEAL_II_OPENMP_SIMD_PRAGMA
                for (unsigned int i = 0; i < length; ++i)
                  {
                    x[i] += this->alpha * p[i];
                    p[i] = prec_r[i];
                    v[i] = Number();
                  }
                p += length;
                r += length;
                v += length;
              }
          }
        else
          {
            const Number alpha_plus_previous_alpha_over_beta =
              this->alpha + this->previous_alpha / this->previous_beta;
            const Number previous_alpha_over_beta =
              this->previous_alpha / this->previous_beta;
            for (unsigned int j = start_range; j < end_range; j += grain_size)
              {
                const unsigned int length = std::min(grain_size, end_range - j);
                this->preconditioner.apply_to_subrange(j,
                                                       j + length,
                                                       r,
                                                       prec_r.data());
                DEAL_II_OPENMP_SIMD_PRAGMA
                for (unsigned int i = 0; i < length; ++i)
                  {
                    x[i] += alpha_plus_previous_alpha_over_beta * p[i] -
                            previous_alpha_over_beta * prec_r[i];
                    r[i] -= this->alpha * v[i];
                  }
                this->preconditioner.apply_to_subrange(j,
                                                       j + length,
                                                       r,
                                                       prec_r.data());
                DEAL_II_OPENMP_SIMD_PRAGMA
                for (unsigned int i = 0; i < length; ++i)
                  {
                    p[i] = this->beta * p[i] + prec_r[i];
                    v[i] = Number();
                  }
                p += length;
                r += length;
                v += length;
                x += length;
              }
          }
      }

      // Function that we use if the PreconditionerType does not implement an
      // apply() function and where we instead need to use the
      // apply_to_subrange function
      template <typename U = void>
      std::enable_if_t<!has_apply<PreconditionerType>, U>
      operation_after_loop(
        const unsigned int                      start_range,
        const unsigned int                      end_range,
        std::array<VectorizedArray<Number>, 7> &vectorized_sums) const
      {
        const Number                          *r          = this->r.begin();
        const Number                          *p          = this->p.begin();
        const Number                          *v          = this->v.begin();
        std::array<VectorizedArray<Number>, 7> my_sums    = {};
        constexpr unsigned int                 grain_size = 128;
        Assert(grain_size % VectorizedArray<Number>::size() == 0,
               ExcNotImplemented());
        const unsigned int end_regular =
          start_range + (end_range - start_range) / grain_size * grain_size;
        std::array<Number, grain_size> prec_r;
        std::array<Number, grain_size> prec_v;
        for (unsigned int j = start_range; j < end_regular; j += grain_size)
          {
            this->preconditioner.apply_to_subrange(j,
                                                   j + grain_size,
                                                   r + j,
                                                   prec_r.data());
            this->preconditioner.apply_to_subrange(j,
                                                   j + grain_size,
                                                   v + j,
                                                   prec_v.data());
            VectorizedArray<Number> pj, vj, rj, prec_vj, prec_rj;
            for (unsigned int i = 0; i < grain_size;
                 i += VectorizedArray<Number>::size())
              {
                pj.load(p + j + i);
                vj.load(v + j + i);
                rj.load(r + j + i);
                prec_rj.load(prec_r.data() + i);
                prec_vj.load(prec_v.data() + i);

                my_sums[0] += pj * vj;
                my_sums[1] += vj * vj;
                my_sums[2] += rj * vj;
                my_sums[3] += rj * rj;
                my_sums[4] += rj * prec_vj;
                my_sums[5] += vj * prec_vj;
                my_sums[6] += rj * prec_rj;
              }
          }
        const unsigned int length = end_range - end_regular;
        AssertIndexRange(length, grain_size);
        this->preconditioner.apply_to_subrange(end_regular,
                                               end_regular + length,
                                               r + end_regular,
                                               prec_r.data());
        this->preconditioner.apply_to_subrange(end_regular,
                                               end_regular + length,
                                               v + end_regular,
                                               prec_v.data());
        for (unsigned int j = end_regular; j < end_range; ++j)
          {
            my_sums[0][0] += p[j] * v[j];
            my_sums[1][0] += v[j] * v[j];
            my_sums[2][0] += r[j] * v[j];
            my_sums[3][0] += r[j] * r[j];
            my_sums[4][0] += r[j] * prec_v[j - end_regular];
            my_sums[5][0] += v[j] * prec_v[j - end_regular];
            my_sums[6][0] += r[j] * prec_r[j - end_regular];
          }
        for (unsigned int i = 0; i < vectorized_sums.size(); ++i)
          vectorized_sums[i] += my_sums[i];
      }

      // Function that we use if the PreconditionerType does not implement an
      // apply() function, where we instead need to choose the
      // apply_to_subrange function
      template <typename U = void>
      std::enable_if_t<!has_apply<PreconditionerType>, U>
      finalize_after_convergence(const unsigned int iteration_index)
      {
        if (iteration_index % 2 == 1 || this->beta == Number())
          this->x.add(this->alpha, this->p);
        else
          {
            const unsigned int end_range = this->x.locally_owned_size();

            Number      *x = this->x.begin();
            Number      *r = this->r.begin();
            Number      *p = this->p.begin();
            const Number alpha_plus_previous_alpha_over_beta =
              this->alpha + this->previous_alpha / this->beta;
            const Number previous_alpha_over_beta =
              this->previous_alpha / this->beta;

            constexpr unsigned int         grain_size = 128;
            std::array<Number, grain_size> prec_r;
            for (unsigned int j = 0; j < end_range; j += grain_size)
              {
                const unsigned int length = std::min(grain_size, end_range - j);
                this->preconditioner.apply_to_subrange(j,
                                                       j + length,
                                                       r,
                                                       prec_r.data());
                DEAL_II_OPENMP_SIMD_PRAGMA
                for (unsigned int i = 0; i < length; ++i)
                  x[i] += alpha_plus_previous_alpha_over_beta * p[i] -
                          previous_alpha_over_beta * prec_r[i];

                x += length;
                r += length;
                p += length;
              }
          }
      }
    };
  } // namespace SolverCG
} // namespace internal



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
template <typename MatrixType, typename PreconditionerType>
DEAL_II_CXX20_REQUIRES(
  (concepts::is_linear_operator_on<MatrixType, VectorType> &&
   concepts::is_linear_operator_on<PreconditionerType, VectorType>))
void SolverCG<VectorType>::solve(const MatrixType         &A,
                                 VectorType               &x,
                                 const VectorType         &b,
                                 const PreconditionerType &preconditioner)
{
  using number = typename VectorType::value_type;

  SolverControl::State solver_state = SolverControl::iterate;

  LogStream::Prefix prefix("cg");

  // Should we build the matrix for eigenvalue computations?
  const bool do_eigenvalues =
    !condition_number_signal.empty() || !all_condition_numbers_signal.empty() ||
    !eigenvalues_signal.empty() || !all_eigenvalues_signal.empty();

  // vectors used for eigenvalue computations
  std::vector<typename VectorType::value_type> diagonal;
  std::vector<typename VectorType::value_type> offdiagonal;

  typename VectorType::value_type eigen_beta_alpha = 0;

  int it = 0;

  internal::SolverCG::
    IterationWorker<VectorType, MatrixType, PreconditionerType>
      worker(A,
             preconditioner,
             determine_beta_by_flexible_formula,
             this->memory,
             x,
             b,
             additional_data.use_default_residual);

  worker.startup();

  solver_state = this->iteration_status(0, worker.residual_norm, x);
  if (solver_state != SolverControl::iterate)
    return;

  while (solver_state == SolverControl::iterate)
    {
      ++it;

      worker.do_iteration(it);

      print_vectors(it, x, worker.r, worker.p);

      if (it > 1)
        {
          this->coefficients_signal(worker.previous_alpha, worker.beta);
          // set up the vectors containing the diagonal and the off diagonal
          // of the projected matrix.
          if (do_eigenvalues)
            {
              diagonal.push_back(number(1.) / worker.previous_alpha +
                                 eigen_beta_alpha);
              eigen_beta_alpha = worker.beta / worker.previous_alpha;
              offdiagonal.push_back(std::sqrt(worker.beta) /
                                    worker.previous_alpha);
            }
          compute_eigs_and_cond(diagonal,
                                offdiagonal,
                                all_eigenvalues_signal,
                                all_condition_numbers_signal);
        }

      solver_state = this->iteration_status(it, worker.residual_norm, x);
    }

  worker.finalize_after_convergence(it);

  compute_eigs_and_cond(diagonal,
                        offdiagonal,
                        eigenvalues_signal,
                        condition_number_signal);

  AssertThrow(solver_state == SolverControl::success,
              SolverControl::NoConvergence(it, worker.residual_norm));
}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
boost::signals2::connection SolverCG<VectorType>::connect_coefficients_slot(
  const std::function<void(typename VectorType::value_type,
                           typename VectorType::value_type)> &slot)
{
  return coefficients_signal.connect(slot);
}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
boost::signals2::connection SolverCG<VectorType>::connect_condition_number_slot(
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
boost::signals2::connection SolverCG<VectorType>::connect_eigenvalues_slot(
  const std::function<void(const std::vector<double> &)> &slot,
  const bool                                              every_iteration)
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
SolverFlexibleCG<VectorType>::AdditionalData::AdditionalData(
  const bool use_default_residual)
  : use_default_residual(use_default_residual)
{}


template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
SolverFlexibleCG<VectorType>::SolverFlexibleCG(SolverControl            &cn,
                                               VectorMemory<VectorType> &mem,
                                               const AdditionalData     &data)
  : SolverCG<VectorType>(cn, mem)
{
  this->determine_beta_by_flexible_formula = true;
  this->additional_data                    = data;
}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
SolverFlexibleCG<VectorType>::SolverFlexibleCG(SolverControl        &cn,
                                               const AdditionalData &data)
  : SolverCG<VectorType>(cn)
{
  this->determine_beta_by_flexible_formula = true;
  this->additional_data                    = data;
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

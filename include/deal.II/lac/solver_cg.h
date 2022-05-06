// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_solver_cg_h
#define dealii_solver_cg_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/tridiagonal_matrix.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

// forward declaration
#ifndef DOXYGEN
class PreconditionIdentity;
#endif


/*!@addtogroup Solvers */
/*@{*/

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
 * dealii::Vector<double>.
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
 */
template <typename VectorType = Vector<double>>
class SolverCG : public SolverBase<VectorType>
{
public:
  /**
   * Declare type for container size.
   */
  using size_type = types::global_dof_index;

  /**
   * Standardized data struct to pipe additional data to the solver.
   * Here, it doesn't store anything but just exists for consistency
   * with the other solver classes.
   */
  struct AdditionalData
  {};

  /**
   * Constructor.
   */
  SolverCG(SolverControl &           cn,
           VectorMemory<VectorType> &mem,
           const AdditionalData &    data = AdditionalData());

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
  void
  solve(const MatrixType &        A,
        VectorType &              x,
        const VectorType &        b,
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
                const VectorType & x,
                const VectorType & r,
                const VectorType & d) const;

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
      &                                          eigenvalues_signal,
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
 * $\beta_k = \frac{\mathbf{r}^T_{k+1} \mathbf{z}^T_{k+1}}{\mathbf{r}^T_{k}
 * \mathbf{z}^T_{k}}$ for computing the new search direction (here
 * $\mathbf{r}_{k+1}$ is the residual in step $k+1$ and $\mathbf{z}_{k+1} =
 * P^{-1} \mathbf{r}_{k+1}$) as in the classical conjugate gradient algorithm,
 * this class selects the Polak-Ribiere formula $\beta_k =
 * \frac{\mathbf{r}^T_{k+1} \left(\mathbf{z}^T_{k+1} -
 * \mathbf{z}^T_{k}\right)}{\mathbf{r}^T_{k} \mathbf{z}^T_{k}}$. The
 * additional term $\mathbf{r}^T_{k+1} \mathbf{z}^T_{k}$ is zero for linear
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
class SolverFlexibleCG : public SolverCG<VectorType>
{
public:
  /**
   * Declare type for container size.
   */
  using size_type = types::global_dof_index;

  /**
   * Standardized data struct to pipe additional data to the solver.
   * Here, it doesn't store anything but just exists for consistency
   * with the other solver classes.
   */
  struct AdditionalData
  {};

  /**
   * Constructor.
   */
  SolverFlexibleCG(SolverControl &           cn,
                   VectorMemory<VectorType> &mem,
                   const AdditionalData &    data = AdditionalData());

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  SolverFlexibleCG(SolverControl &       cn,
                   const AdditionalData &data = AdditionalData());
};


/*@}*/

/*------------------------- Implementation ----------------------------*/

#ifndef DOXYGEN



template <typename VectorType>
SolverCG<VectorType>::SolverCG(SolverControl &           cn,
                               VectorMemory<VectorType> &mem,
                               const AdditionalData &    data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
  , determine_beta_by_flexible_formula(false)
{}



template <typename VectorType>
SolverCG<VectorType>::SolverCG(SolverControl &cn, const AdditionalData &data)
  : SolverBase<VectorType>(cn)
  , additional_data(data)
  , determine_beta_by_flexible_formula(false)
{}



template <typename VectorType>
void
SolverCG<VectorType>::print_vectors(const unsigned int,
                                    const VectorType &,
                                    const VectorType &,
                                    const VectorType &) const
{}



template <typename VectorType>
inline void
SolverCG<VectorType>::compute_eigs_and_cond(
  const std::vector<typename VectorType::value_type> &diagonal,
  const std::vector<typename VectorType::value_type> &offdiagonal,
  const boost::signals2::signal<void(const std::vector<double> &)>
    &                                          eigenvalues_signal,
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
    // Detector class to find out whether the MatrixType in SolverCG has a
    // vmult function that takes two additional std::function objects, which
    // we use to run the operation on slices of the vector during the
    // matrix-vector product, and whether PreconditionerType can run
    // operations on an individual element
    template <typename MatrixType,
              typename VectorType,
              typename PreconditionerType>
    struct supports_vmult_with_std_functions
    {
    private:
      // this will work always
      static bool
      detect_matrix(...);

      // this detector will work only if we have
      // "... MatrixType::vmult(VectorType, const VectorType,
      // const std::function<...>&, const std::function<...>&) const"
      template <typename MatrixType2>
      static decltype(std::declval<MatrixType2 const>().vmult(
        std::declval<VectorType &>(),
        std::declval<const VectorType &>(),
        std::declval<const std::function<void(const unsigned int,
                                              const unsigned int)> &>(),
        std::declval<const std::function<void(const unsigned int,
                                              const unsigned int)> &>()))
      detect_matrix(const MatrixType2 &);

      // this will work always
      static bool
      detect_preconditioner(...);

      // this detector will work only if we have
      // "... PreconditionerType::vmult(std::size_t, std::size_t, Number, const
      // VectorType, const std::function<...>&, const std::function<...>&)
      // const"
      template <typename PreconditionerType2>
      static decltype(
        std::declval<PreconditionerType2 const>().apply_to_subrange(
          std::size_t(),
          std::size_t(),
          std::declval<const typename PreconditionerType2::value_type *>(),
          std::declval<typename PreconditionerType2::value_type *>()))
      detect_preconditioner(const PreconditionerType2 &);

    public:
      // finally here we check if both our detectors have void return
      // type. This will happen if the compiler can use the second detector,
      // otherwise SFINAE let's it work with the more general first one that
      // is bool
      static const bool value =
        !std::is_same<decltype(detect_matrix(std::declval<MatrixType>())),
                      bool>::value &&
        !std::is_same<decltype(detect_preconditioner(
                        std::declval<PreconditionerType>())),
                      bool>::value;
    };

    // We need to have a separate declaration for static const members
    template <typename T, typename U, typename V>
    const bool supports_vmult_with_std_functions<T, U, V>::value;

    /**
     * Internal function to run one iteration of the conjugate gradient solver
     * for standard matrix and preconditioner arguments.
     */
    template <typename Number,
              typename VectorType,
              typename MatrixType,
              typename PreconditionerType,
              typename std::enable_if<
                !supports_vmult_with_std_functions<MatrixType,
                                                   VectorType,
                                                   PreconditionerType>::value,
                int>::type * = nullptr>
    void
    do_cg_iteration(const MatrixType &        A,
                    const PreconditionerType &preconditioner,
                    const typename dealii::SolverCG<VectorType>::AdditionalData
                      &                additional_data,
                    const unsigned int iteration,
                    VectorType &       x,
                    VectorType &       r,
                    VectorType &       p,
                    VectorType &       v,
                    VectorType &       z,
                    Number &           r_dot_preconditioner_dot_r,
                    Number &           alpha,
                    Number &           beta,
                    double &           residual_norm,
                    const Number /*old_alpha*/,
                    const Number /*old_beta*/)
    {
      const Number old_r_dot_preconditioner_dot_r = r_dot_preconditioner_dot_r;

      if (std::is_same<PreconditionerType, PreconditionIdentity>::value ==
          false)
        {
          preconditioner.vmult(v, r);
          r_dot_preconditioner_dot_r = r * v;
        }
      else
        r_dot_preconditioner_dot_r = residual_norm * residual_norm;

      const VectorType &direction =
        std::is_same<PreconditionerType, PreconditionIdentity>::value ? r : v;

      if (iteration > 1)
        {
          Assert(std::abs(old_r_dot_preconditioner_dot_r) != 0.,
                 ExcDivideByZero());
          beta = r_dot_preconditioner_dot_r / old_r_dot_preconditioner_dot_r;
          if (additional_data.use_flexible_variant)
            beta -= (r * z) / old_r_dot_preconditioner_dot_r;
          p.sadd(beta, 1., direction);
        }
      else
        p.equ(1., direction);

      if (additional_data.use_flexible_variant)
        z.swap(v);

      A.vmult(v, p);

      const Number p_dot_A_dot_p = p * v;
      Assert(std::abs(p_dot_A_dot_p) != 0., ExcDivideByZero());
      alpha = r_dot_preconditioner_dot_r / p_dot_A_dot_p;

      x.add(alpha, p);
      residual_norm = std::sqrt(std::abs(r.add_and_dot(-alpha, v, r)));
    }

    /**
     * Internal function to run one iteration of the conjugate gradient solver
     * for matrices and preconditioners that support interleaving the vector
     * updates with the matrix-vector product.
     */
    template <typename Number,
              typename VectorType,
              typename MatrixType,
              typename PreconditionerType,
              typename std::enable_if<
                supports_vmult_with_std_functions<MatrixType,
                                                  VectorType,
                                                  PreconditionerType>::value,
                int>::type * = nullptr>
    void
    do_cg_iteration(
      const MatrixType &        A,
      const PreconditionerType &preconditioner,
      const typename dealii::SolverCG<VectorType>::AdditionalData &,
      const unsigned int iteration,
      VectorType &       x_vector,
      VectorType &       r_vector,
      VectorType &       p_vector,
      VectorType &       v_vector,
      VectorType &,
      Number &     r_dot_preconditioner_dot_r,
      Number &     alpha,
      Number &     beta,
      double &     residual_norm,
      const Number old_alpha,
      const Number old_beta)
    {
      const auto operation_before_loop = [&](const unsigned int start_range,
                                             const unsigned int end_range) {
        Number *                       x = x_vector.begin() + start_range;
        Number *                       r = r_vector.begin() + start_range;
        Number *                       p = p_vector.begin() + start_range;
        Number *                       v = v_vector.begin() + start_range;
        constexpr unsigned int         grain_size = 32;
        std::array<Number, grain_size> prec_r;
        if (iteration == 1)
          {
            for (unsigned int j = start_range; j < end_range; j += grain_size)
              {
                const unsigned int length = std::min(grain_size, end_range - j);
                preconditioner.apply_on_subrange(j, length, r, prec_r.data());
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
        else if (iteration % 2 == 0)
          {
            for (unsigned int j = start_range; j < end_range; j += grain_size)
              {
                const unsigned int length = std::min(grain_size, end_range - j);
                DEAL_II_OPENMP_SIMD_PRAGMA
                for (unsigned int i = 0; i < length; ++i)
                  r[i] -= alpha * v[i];
                preconditioner.apply_on_subrange(j, length, r, prec_r.data());
                DEAL_II_OPENMP_SIMD_PRAGMA
                for (unsigned int i = 0; i < length; ++i)
                  {
                    p[i] = beta * p[i] + prec_r[i];
                    v[i] = Number();
                  }
                p += length;
                r += length;
                v += length;
              }
          }
        else
          {
            const Number alpha_plus_alpha_old = alpha + old_alpha / old_beta;
            const Number alpha_old_beta_old   = old_alpha / old_beta;
            for (unsigned int j = start_range; j < end_range; j += grain_size)
              {
                const unsigned int length = std::min(grain_size, end_range - j);
                preconditioner.apply_on_subrange(j, length, r, prec_r.data());
                DEAL_II_OPENMP_SIMD_PRAGMA
                for (unsigned int i = 0; i < length; ++i)
                  {
                    x[i] += alpha_plus_alpha_old * p[i] +
                            alpha_old_beta_old * prec_r[i];
                    r[i] -= alpha * v[i];
                  }
                preconditioner.apply_on_subrange(j, length, r, prec_r.data());
                DEAL_II_OPENMP_SIMD_PRAGMA
                for (unsigned int i = 0; i < length; ++i)
                  {
                    p[i] = beta * p[i] + prec_r[i];
                    v[i] = 0.;
                  }
                p += length;
                r += length;
                v += length;
                x += length;
              }
          }
      };

      std::array<Number, 7> local_sums = {};
      const auto operation_after_loop  = [&](const unsigned int start_range,
                                            const unsigned int end_range) {
        Number *                       x = x_vector.begin() + start_range;
        Number *                       r = r_vector.begin() + start_range;
        Number *                       p = p_vector.begin() + start_range;
        Number *                       v = v_vector.begin() + start_range;
        constexpr unsigned int         grain_size = 32;
        std::array<Number, grain_size> prec_r;
        std::array<Number, grain_size> prec_v;
        for (unsigned int j = start_range; j < end_range; j += grain_size)
          {
            const unsigned int length = std::min(grain_size, end_range - j);
            preconditioner.apply_on_subrange(j, length, r, prec_r.data());
            preconditioner.apply_on_subrange(j, length, v, prec_v.data());
            for (unsigned int i = 0; i < length; ++i)
              {
                local_sums[0] += p[i] * v[i];
                local_sums[1] += v[i] * v[i];
                local_sums[2] += r[i] * v[i];
                local_sums[3] += r[i] * r[i];
                local_sums[4] += r[i] * prec_v[i];
                local_sums[5] += v[i] * prec_v[i];
                local_sums[6] += r[i] * prec_r[i];
              }
          }
      };

      A.vmult(v_vector, p_vector, operation_before_loop, operation_after_loop);

      Utilities::MPI::sum(dealii::ArrayView<const double>(local_sums.begin(),
                                                          7),
                          r_vector.get_mpi_communicator(),
                          dealii::ArrayView<double>(local_sums.begin_raw(), 7));

      const Number p_dot_A_dot_p = local_sums[0];
      Assert(std::abs(p_dot_A_dot_p) != 0., ExcDivideByZero());

      const Number old_r_dot_preconditioner_dot_r = local_sums[6];
      alpha         = old_r_dot_preconditioner_dot_r / p_dot_A_dot_p;
      residual_norm = std::sqrt(local_sums[3] + 2 * local_sums[2] +
                                alpha * alpha * local_sums[1]);

      r_dot_preconditioner_dot_r =
        old_r_dot_preconditioner_dot_r -
        alpha * (2. * local_sums[4] - alpha * local_sums[5]);
      beta = r_dot_preconditioner_dot_r / old_r_dot_preconditioner_dot_r;
    }
  } // namespace SolverCG
} // namespace internal



template <typename VectorType>
template <typename MatrixType, typename PreconditionerType>
void
SolverCG<VectorType>::solve(const MatrixType &        A,
                            VectorType &              x,
                            const VectorType &        b,
                            const PreconditionerType &preconditioner)
{
  using number = typename VectorType::value_type;

  SolverControl::State solver_state = SolverControl::iterate;

  LogStream::Prefix prefix("cg");

  // Memory allocation
  typename VectorMemory<VectorType>::Pointer r_pointer(this->memory);
  typename VectorMemory<VectorType>::Pointer p_pointer(this->memory);
  typename VectorMemory<VectorType>::Pointer v_pointer(this->memory);
  typename VectorMemory<VectorType>::Pointer z_pointer(this->memory);

  // Define some aliases for simpler access, using the variables 'r' for the
  // residual b - A*x, 'p' for the search direction, and 'v' for the auxiliary
  // vector. This naming convention is used e.g. by the description on
  // https://en.wikipedia.org/wiki/Conjugate_gradient_method. The variable 'z'
  // gets only used for the flexible variant of the CG method.
  VectorType &r = *r_pointer;
  VectorType &p = *p_pointer;
  VectorType &v = *v_pointer;
  VectorType &z = *z_pointer;

  // Should we build the matrix for eigenvalue computations?
  const bool do_eigenvalues =
    !condition_number_signal.empty() || !all_condition_numbers_signal.empty() ||
    !eigenvalues_signal.empty() || !all_eigenvalues_signal.empty();

  // vectors used for eigenvalue computations
  std::vector<typename VectorType::value_type> diagonal;
  std::vector<typename VectorType::value_type> offdiagonal;

  typename VectorType::value_type eigen_beta_alpha = 0;

  // resize the vectors, but do not set the values since they'd be overwritten
  // soon anyway.
  r.reinit(x, true);
  p.reinit(x, true);
  v.reinit(x, true);
  if (determine_beta_by_flexible_formula)
    z.reinit(x, true);

  int    it                         = 0;
  number r_dot_preconditioner_dot_r = number();
  number beta                       = number();
  number alpha                      = number();
  number old_beta                   = number();
  number old_alpha                  = number();

  // compute residual. if vector is zero, then short-circuit the full
  // computation
  if (!x.all_zero())
    {
      A.vmult(r, x);
      r.sadd(-1., 1., b);
    }
  else
    r.equ(1., b);

  double residual_norm = r.l2_norm();
  solver_state         = this->iteration_status(0, residual_norm, x);
  if (solver_state != SolverControl::iterate)
    return;

  while (solver_state == SolverControl::iterate)
    {
      it++;

      internal::SolverCG::do_cg_iteration(A,
                                          preconditioner,
                                          additional_data,
                                          it,
                                          x,
                                          r,
                                          p,
                                          v,
                                          z,
                                          r_dot_preconditioner_dot_r,
                                          alpha,
                                          beta,
                                          residual_norm,
                                          old_alpha,
                                          old_beta);

      old_alpha = alpha;
      old_beta = beta;

      print_vectors(it, x, r, p);

      if (it > 1)
        {
          this->coefficients_signal(old_alpha, beta);
          // set up the vectors containing the diagonal and the off diagonal of
          // the projected matrix.
          if (do_eigenvalues)
            {
              diagonal.push_back(number(1.) / old_alpha + eigen_beta_alpha);
              eigen_beta_alpha = beta / old_alpha;
              offdiagonal.push_back(std::sqrt(beta) / old_alpha);
            }
          compute_eigs_and_cond(diagonal,
                                offdiagonal,
                                all_eigenvalues_signal,
                                all_condition_numbers_signal);
        }

      solver_state = this->iteration_status(it, residual_norm, x);
    }

  compute_eigs_and_cond(diagonal,
                        offdiagonal,
                        eigenvalues_signal,
                        condition_number_signal);

  AssertThrow(solver_state == SolverControl::success,
              SolverControl::NoConvergence(it, residual_norm));
}



template <typename VectorType>
boost::signals2::connection
SolverCG<VectorType>::connect_coefficients_slot(
  const std::function<void(typename VectorType::value_type,
                           typename VectorType::value_type)> &slot)
{
  return coefficients_signal.connect(slot);
}



template <typename VectorType>
boost::signals2::connection
SolverCG<VectorType>::connect_condition_number_slot(
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
boost::signals2::connection
SolverCG<VectorType>::connect_eigenvalues_slot(
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
SolverFlexibleCG<VectorType>::SolverFlexibleCG(SolverControl &           cn,
                                               VectorMemory<VectorType> &mem,
                                               const AdditionalData &)
  : SolverCG<VectorType>(cn, mem)
{
  this->determine_beta_by_flexible_formula = true;
}



template <typename VectorType>
SolverFlexibleCG<VectorType>::SolverFlexibleCG(SolverControl &cn,
                                               const AdditionalData &)
  : SolverCG<VectorType>(cn)
{
  this->determine_beta_by_flexible_formula = true;
}



#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

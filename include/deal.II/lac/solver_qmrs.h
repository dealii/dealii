// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_solver_qmrs_h
#define dealii_solver_qmrs_h

#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup Solvers
 * @{
 */

/**
 * <h3>Quasi-minimal method for symmetric matrices (SQMR)</h3>
 *
 * The SQMR (symmetric quasi-minimal residual) method is supposed to solve
 * symmetric indefinite linear systems with symmetric, not necessarily definite
 * preconditioners. It is a variant of the original quasi-minimal residual
 * method (QMR) and produces the same iterative solution. This version of SQMR
 * is adapted from the respective symmetric QMR-from-BiCG algorithm given by
 * both Freund/Nachtigal: A new Krylov-subspace method for symmetric indefinite
 * linear systems, NASA STI/Recon Technical Report N, 95 (1994) and
 * Freund/Nachtigal: Software for simplified Lanczos and QMR algorithms, Appl.
 * Num. Math. 19 (1995), pp. 319-341 and provides both right and left (but not
 * split) preconditioning.
 *
 *
 * <h3>Trade off of stability to simplicity</h3>
 *
 * Note, that the QMR implementation that the given algorithm is based on is
 * derived from classical BiCG. It can be shown (Freund/Szeto: A transpose-free
 * quasi-minimal residual squared algorithm for non-Hermitian linear systems,
 * Advances in Computer Methods for Partial Differential Equations VII
 * (IMACS, New Brunswick, NJ, 1992) pp. 258-264) that the QMR iterates can
 * be generated from the BiCG iteration through one additional vector and
 * some scalar updates. Possible breakdowns (or precisely, divisions by
 * zero) of BiCG therefore obviously transfer to this simple no-look-ahead
 * algorithm.
 *
 * In return the algorithm is cheap compared to classical QMR or BiCGStab,
 * using only one matrix-vector product with the system matrix and
 * one application of the preconditioner per iteration respectively.
 *
 * The residual used for measuring convergence is only approximately calculated
 * by an upper bound. If this value comes below a threshold prescribed within
 * the AdditionalData struct, then the exact residual of the current QMR iterate
 * will be calculated using another multiplication with the system matrix. By
 * experience (according to Freund and Nachtigal) this technique is useful for a
 * threshold that is ten times the solving tolerance, and in that case will be
 * only used in the last one or two steps of the complete iteration.
 *
 * For the requirements on matrices and vectors in order to work with this
 * class, see the documentation of the Solver base class.
 *
 * Like all other solver classes, this class has a local structure called @p
 * AdditionalData which is used to pass additional parameters to the solver,
 * like damping parameters or the number of temporary vectors. We use this
 * additional structure instead of passing these values directly to the
 * constructor because this makes the use of the @p SolverSelector and other
 * classes much easier and guarantees that these will continue to work even if
 * number or type of the additional parameters for a certain solver changes.
 *
 *
 * <h3>Observing the progress of linear solver iterations</h3>
 *
 * The solve() function of this class uses the mechanism described in the
 * Solver base class to determine convergence. This mechanism can also be used
 * to observe the progress of the iteration.
 */
template <typename VectorType = Vector<double>>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
class SolverQMRS : public SolverBase<VectorType>
{
public:
  /**
   * Standardized data struct to pipe additional data to the solver.
   *
   * The user is able to switch between right and left preconditioning, that
   * means solving the systems <i>P<sup>-1</sup>A</i> and <i>AP<sup>-1</sup></i>
   * respectively, using the corresponding parameter. Note that left
   * preconditioning means to employ the preconditioned (BiCG-)residual and
   * otherwise the unpreconditioned one. The default is the application from the
   * right side.
   *
   * The @p solver_tolerance threshold is used to define the said bound below which the residual
   * is computed exactly. See the class documentation for more information. The
   * default value is 1e-9, that is the default solving precision multiplied by
   * ten.
   *
   * SQMR is susceptible to breakdowns (divisions by zero), so we need a
   * parameter telling us which numbers are considered zero. The proper
   * breakdown criterion is very unclear, so experiments may be necessary here.
   * It is even possible to achieve convergence despite of dividing through by
   * small numbers. There are even cases in which it is advantageous to accept
   * such divisions because the cheap iteration cost makes the algorithm the
   * fastest of all available indefinite iterative solvers. Nonetheless, the
   * default breakdown threshold value is 1e-16.
   */
  struct AdditionalData
  {
    /**
     * Constructor.
     *
     * The default is right preconditioning, with the @p solver_tolerance chosen to be 1e-9 and
     * the @p breakdown_threshold set at 1e-16.
     */
    explicit AdditionalData(const bool   left_preconditioning = false,
                            const double solver_tolerance     = 1.e-9,
                            const bool   breakdown_testing    = true,
                            const double breakdown_threshold  = 1.e-16)
      : left_preconditioning(left_preconditioning)
      , solver_tolerance(solver_tolerance)
      , breakdown_testing(breakdown_testing)
      , breakdown_threshold(breakdown_threshold)
    {}

    /**
     * Flag for using a left-preconditioned version.
     */
    bool left_preconditioning;

    /**
     * The threshold below which the current residual is computed exactly.
     */
    double solver_tolerance;

    /**
     * Flag for breakdown testing.
     */
    bool breakdown_testing;

    /**
     * Breakdown threshold. Scalars measured to this bound are used for
     * divisions.
     */
    double breakdown_threshold;
  };

  /**
   * Constructor.
   */
  SolverQMRS(SolverControl            &cn,
             VectorMemory<VectorType> &mem,
             const AdditionalData     &data = AdditionalData());

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  SolverQMRS(SolverControl &cn, const AdditionalData &data = AdditionalData());

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
   * Interface for derived class. This function gets the current iteration
   * vector, the residual and the update vector in each step. It can be used
   * for a graphical output of the convergence history.
   */
  virtual void
  print_vectors(const unsigned int step,
                const VectorType  &x,
                const VectorType  &r,
                const VectorType  &d) const;

protected:
  /**
   * Additional parameters.
   */
  AdditionalData additional_data;

private:
  /**
   * A structure returned by the iterate() function representing what it found
   * is happening during the iteration.
   */
  struct IterationResult
  {
    SolverControl::State state;
    double               last_residual;

    IterationResult(const SolverControl::State state,
                    const double               last_residual);
  };

  /**
   * The iteration loop itself. The function returns a structure indicating
   * what happened in this function.
   */
  template <typename MatrixType, typename PreconditionerType>
  IterationResult
  iterate(const MatrixType         &A,
          VectorType               &x,
          const VectorType         &b,
          const PreconditionerType &preconditioner,
          VectorType               &r,
          VectorType               &u,
          VectorType               &q,
          VectorType               &t,
          VectorType               &d);

  /**
   * Number of the current iteration (accumulated over restarts)
   */
  unsigned int step;
};

/** @} */
/*------------------------- Implementation ----------------------------*/

#ifndef DOXYGEN


template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
SolverQMRS<VectorType>::IterationResult::IterationResult(
  const SolverControl::State state,
  const double               last_residual)
  : state(state)
  , last_residual(last_residual)
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
SolverQMRS<VectorType>::SolverQMRS(SolverControl            &cn,
                                   VectorMemory<VectorType> &mem,
                                   const AdditionalData     &data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
  , step(0)
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
SolverQMRS<VectorType>::SolverQMRS(SolverControl        &cn,
                                   const AdditionalData &data)
  : SolverBase<VectorType>(cn)
  , additional_data(data)
  , step(0)
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
void SolverQMRS<VectorType>::print_vectors(const unsigned int,
                                           const VectorType &,
                                           const VectorType &,
                                           const VectorType &) const
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
template <typename MatrixType, typename PreconditionerType>
DEAL_II_CXX20_REQUIRES(
  (concepts::is_linear_operator_on<MatrixType, VectorType> &&
   concepts::is_linear_operator_on<PreconditionerType, VectorType>))
void SolverQMRS<VectorType>::solve(const MatrixType         &A,
                                   VectorType               &x,
                                   const VectorType         &b,
                                   const PreconditionerType &preconditioner)
{
  LogStream::Prefix prefix("SQMR");


  // temporary vectors, allocated through the @p VectorMemory object at the
  // start of the actual solution process and deallocated at the end.
  typename VectorMemory<VectorType>::Pointer Vr(this->memory);
  typename VectorMemory<VectorType>::Pointer Vu(this->memory);
  typename VectorMemory<VectorType>::Pointer Vq(this->memory);
  typename VectorMemory<VectorType>::Pointer Vt(this->memory);
  typename VectorMemory<VectorType>::Pointer Vd(this->memory);


  // resize the vectors, but do not set
  // the values since they'd be overwritten
  // soon anyway.
  Vr->reinit(x, true);
  Vu->reinit(x, true);
  Vq->reinit(x, true);
  Vt->reinit(x, true);
  Vd->reinit(x, true);

  step = 0;

  IterationResult state(SolverControl::failure, 0);

  do
    {
      if (step > 0)
        deallog << "Restart step " << step << std::endl;
      state = iterate(A, x, b, preconditioner, *Vr, *Vu, *Vq, *Vt, *Vd);
    }
  while (state.state == SolverControl::iterate);


  // in case of failure: throw exception
  AssertThrow(state.state == SolverControl::success,
              SolverControl::NoConvergence(step, state.last_residual));
  // otherwise exit as normal
}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
template <typename MatrixType, typename PreconditionerType>
typename SolverQMRS<VectorType>::IterationResult
  SolverQMRS<VectorType>::iterate(const MatrixType         &A,
                                  VectorType               &x,
                                  const VectorType         &b,
                                  const PreconditionerType &preconditioner,
                                  VectorType               &r,
                                  VectorType               &u,
                                  VectorType               &q,
                                  VectorType               &t,
                                  VectorType               &d)
{
  SolverControl::State state = SolverControl::iterate;

  int it = 0;

  double tau, rho, theta = 0;
  double res;

  // Compute the start residual
  A.vmult(r, x);
  r.sadd(-1., 1., b);

  // Doing the initial preconditioning
  if (additional_data.left_preconditioning)
    {
      // Left preconditioning
      preconditioner.vmult(t, r);
      q = t;
    }
  else
    {
      // Right preconditioning
      t = r;
      preconditioner.vmult(q, t);
    }

  tau = t.norm_sqr();
  res = std::sqrt(tau);

  if (this->iteration_status(step, res, x) == SolverControl::success)
    return IterationResult(SolverControl::success, res);

  rho = q * r;

  while (state == SolverControl::iterate)
    {
      ++step;
      ++it;
      //--------------------------------------------------------------
      // Step 1: apply the system matrix and compute one inner product
      //--------------------------------------------------------------
      A.vmult(t, q);
      const double sigma = q * t;

      // Check the breakdown criterion
      if (additional_data.breakdown_testing == true &&
          std::fabs(sigma) < additional_data.breakdown_threshold)
        return IterationResult(SolverControl::iterate, res);
      // Update the residual
      const double alpha = rho / sigma;
      r.add(-alpha, t);

      //--------------------------------------------------------------
      // Step 2: update the solution vector
      //--------------------------------------------------------------
      const double theta_old = theta;

      // Apply the preconditioner
      if (additional_data.left_preconditioning)
        {
          // Left Preconditioning
          preconditioner.vmult(t, r);
        }
      else
        {
          // Right Preconditioning
          t = r;
        }

      // Double updates
      theta            = t * t / tau;
      const double psi = 1. / (1. + theta);
      tau *= theta * psi;

      // Actual update of the solution vector
      d.sadd(psi * theta_old, psi * alpha, q);
      x += d;

      print_vectors(step, x, r, d);

      // Check for convergence
      // Compute a simple and cheap upper bound of the norm of the residual
      // vector b-Ax
      res = std::sqrt((it + 1) * tau);
      // If res lies close enough, within the desired tolerance, calculate the
      // exact residual
      if (res < additional_data.solver_tolerance)
        {
          A.vmult(u, x);
          u.sadd(-1., 1., b);
          res = u.l2_norm();
        }
      state = this->iteration_status(step, res, x);
      if ((state == SolverControl::success) ||
          (state == SolverControl::failure))
        return IterationResult(state, res);

      //--------------------------------------------------------------
      // Step 3: check breakdown criterion and update the vectors
      //--------------------------------------------------------------
      if (additional_data.breakdown_testing == true &&
          std::fabs(sigma) < additional_data.breakdown_threshold)
        return IterationResult(SolverControl::iterate, res);

      const double rho_old = rho;

      // Applying the preconditioner
      if (additional_data.left_preconditioning)
        {
          // Left preconditioning
          u = t;
        }
      else
        {
          // Right preconditioning
          preconditioner.vmult(u, t);
        }

      // Double and vector updates
      rho               = u * r;
      const double beta = rho / rho_old;
      q.sadd(beta, 1., u);
    }
  return IterationResult(SolverControl::success, res);
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

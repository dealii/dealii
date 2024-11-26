// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_solver_bicgstab_h
#define dealii_solver_bicgstab_h


#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>

#include <cmath>
#include <limits>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup Solvers
 * @{
 */

/**
 * Bicgstab algorithm by van der Vorst.
 *
 * For the requirements on matrices and vectors in order to work with this
 * class, see the documentation of the SolverBase base class.
 *
 * Like all other solver classes, this class has a local structure called @p
 * AdditionalData which is used to pass additional parameters to the solver,
 * like damping parameters or the number of temporary vectors. We use this
 * additional structure instead of passing these values directly to the
 * constructor because this makes the use of the @p SolverSelector and other
 * classes much easier and guarantees that these will continue to work even if
 * number or type of the additional parameters for a certain solver changes.
 *
 * The Bicgstab method has two additional parameters found in the
 * SolverBicgstab::AdditionalData struct: the first, @p exact_residual is a boolean,
 * deciding whether to compute the actual residual in each step (@p true) or
 * to use the length of the computed orthogonal residual (@p false). Note that
 * computing the residual causes a third matrix-vector multiplication, though
 * no additional preconditioning, in each step. The reason for doing this is,
 * that the size of the orthogonalized residual computed during the iteration
 * may be larger by orders of magnitude than the true residual. This is due to
 * numerical instabilities related to badly conditioned matrices. Since this
 * instability results in a bad stopping criterion, the default for this
 * parameter is @p true. Whenever the user knows that the estimated residual
 * works reasonably as well, the flag should be set to @p false in order to
 * increase the performance of the solver.
 *
 * The second parameter @p breakdown is the size of a breakdown criterion. It is difficult
 * to find a general good criterion, so if things do not work for you, try to
 * change this value.
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
class SolverBicgstab : public SolverBase<VectorType>
{
public:
  /**
   * There are two possibilities to compute the residual: one is an estimate
   * using the computed value @p tau. The other is exact computation using
   * another matrix vector multiplication. This increases the costs of the
   * algorithm, so it is should be set to false whenever the problem allows
   * it.
   *
   * Bicgstab is susceptible to breakdowns, so we need a parameter telling us,
   * which numbers are considered zero.
   */
  struct AdditionalData
  {
    /**
     * Constructor.
     *
     * The default is to perform an exact residual computation and breakdown
     * parameter is the minimum finite value representable by the value_type of
     * VectorType.
     */
    explicit AdditionalData(
      const bool   exact_residual = true,
      const double breakdown =
        std::numeric_limits<typename VectorType::value_type>::min())
      : exact_residual(exact_residual)
      , breakdown(breakdown)
    {}
    /**
     * Flag for exact computation of residual.
     */
    bool exact_residual;
    /**
     * Breakdown threshold.
     */
    double breakdown;
  };

  /**
   * Constructor.
   */
  SolverBicgstab(SolverControl            &cn,
                 VectorMemory<VectorType> &mem,
                 const AdditionalData     &data = AdditionalData());

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  SolverBicgstab(SolverControl        &cn,
                 const AdditionalData &data = AdditionalData());

  /**
   * Virtual destructor.
   */
  virtual ~SolverBicgstab() override = default;

  /**
   * Solve primal problem only.
   */
  template <typename MatrixType, typename PreconditionerType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_linear_operator_on<MatrixType, VectorType> &&
     concepts::is_linear_operator_on<PreconditionerType, VectorType>))
  void solve(const MatrixType         &A,
             VectorType               &x,
             const VectorType         &b,
             const PreconditionerType &preconditioner);

protected:
  /**
   * Computation of the stopping criterion.
   */
  template <typename MatrixType>
  double
  criterion(const MatrixType &A,
            const VectorType &x,
            const VectorType &b,
            VectorType       &t);

  /**
   * Interface for derived class.  This function gets the current iteration
   * vector, the residual and the update vector in each step. It can be used
   * for graphical output of the convergence history.
   */
  virtual void
  print_vectors(const unsigned int step,
                const VectorType  &x,
                const VectorType  &r,
                const VectorType  &d) const;

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
    bool                 breakdown;
    SolverControl::State state;
    unsigned int         last_step;
    double               last_residual;

    IterationResult(const bool                 breakdown,
                    const SolverControl::State state,
                    const unsigned int         last_step,
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
          const unsigned int        step);
};


/** @} */
/*-------------------------Inline functions -------------------------------*/

#ifndef DOXYGEN


template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
SolverBicgstab<VectorType>::IterationResult::IterationResult(
  const bool                 breakdown,
  const SolverControl::State state,
  const unsigned int         last_step,
  const double               last_residual)
  : breakdown(breakdown)
  , state(state)
  , last_step(last_step)
  , last_residual(last_residual)
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
SolverBicgstab<VectorType>::SolverBicgstab(SolverControl            &cn,
                                           VectorMemory<VectorType> &mem,
                                           const AdditionalData     &data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
SolverBicgstab<VectorType>::SolverBicgstab(SolverControl        &cn,
                                           const AdditionalData &data)
  : SolverBase<VectorType>(cn)
  , additional_data(data)
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
template <typename MatrixType>
double SolverBicgstab<VectorType>::criterion(const MatrixType &A,
                                             const VectorType &x,
                                             const VectorType &b,
                                             VectorType       &t)
{
  A.vmult(t, x);
  return std::sqrt(t.add_and_dot(-1.0, b, t));
}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
void SolverBicgstab<VectorType>::print_vectors(const unsigned int,
                                               const VectorType &,
                                               const VectorType &,
                                               const VectorType &) const
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
template <typename MatrixType, typename PreconditionerType>
typename SolverBicgstab<VectorType>::IterationResult
  SolverBicgstab<VectorType>::iterate(const MatrixType         &A,
                                      VectorType               &x,
                                      const VectorType         &b,
                                      const PreconditionerType &preconditioner,
                                      const unsigned int        last_step)
{
  // Allocate temporary memory.
  typename VectorMemory<VectorType>::Pointer Vr(this->memory);
  typename VectorMemory<VectorType>::Pointer Vrbar(this->memory);
  typename VectorMemory<VectorType>::Pointer Vp(this->memory);
  typename VectorMemory<VectorType>::Pointer Vy(this->memory);
  typename VectorMemory<VectorType>::Pointer Vz(this->memory);
  typename VectorMemory<VectorType>::Pointer Vt(this->memory);
  typename VectorMemory<VectorType>::Pointer Vv(this->memory);

  // Define a few aliases for simpler use of the vectors
  VectorType &r    = *Vr;
  VectorType &rbar = *Vrbar;
  VectorType &p    = *Vp;
  VectorType &y    = *Vy;
  VectorType &z    = *Vz;
  VectorType &t    = *Vt;
  VectorType &v    = *Vv;

  r.reinit(x, true);
  rbar.reinit(x, true);
  p.reinit(x, true);
  y.reinit(x, true);
  z.reinit(x, true);
  t.reinit(x, true);
  v.reinit(x, true);

  using value_type = typename VectorType::value_type;
  using real_type  = typename numbers::NumberTraits<value_type>::real_type;

  A.vmult(r, x);
  r.sadd(-1., 1., b);
  value_type res = r.l2_norm();

  unsigned int step = last_step;

  SolverControl::State state = this->iteration_status(step, res, x);
  if (state == SolverControl::State::success)
    return IterationResult(false, state, step, res);

  rbar = r;

  value_type alpha = 1.;
  value_type rho   = 1.;
  value_type omega = 1.;

  do
    {
      ++step;

      const value_type rhobar = (step == 1 + last_step) ? res * res : r * rbar;

      if (std::fabs(rhobar) < additional_data.breakdown)
        {
          return IterationResult(true, state, step, res);
        }

      const value_type beta = rhobar * alpha / (rho * omega);
      rho                   = rhobar;
      if (step == last_step + 1)
        {
          p = r;
        }
      else
        {
          p.sadd(beta, 1., r);
          p.add(-beta * omega, v);
        }

      preconditioner.vmult(y, p);
      A.vmult(v, y);
      const value_type rbar_dot_v = rbar * v;
      if (std::fabs(rbar_dot_v) < additional_data.breakdown)
        {
          return IterationResult(true, state, step, res);
        }

      alpha = rho / rbar_dot_v;

      res = std::sqrt(real_type(r.add_and_dot(-alpha, v, r)));

      // check for early success, see the lac/bicgstab_early testcase as to
      // why this is necessary
      //
      // note: the vector *Vx we pass to the iteration_status signal here is
      // only the current approximation, not the one we will return with, which
      // will be x=*Vx + alpha*y
      if (this->iteration_status(step, res, x) == SolverControl::success)
        {
          x.add(alpha, y);
          print_vectors(step, x, r, y);
          return IterationResult(false, SolverControl::success, step, res);
        }

      preconditioner.vmult(z, r);
      A.vmult(t, z);
      const value_type t_dot_r   = t * r;
      const real_type  t_squared = t * t;
      if (t_squared < additional_data.breakdown)
        {
          return IterationResult(true, state, step, res);
        }
      omega = t_dot_r / t_squared;
      x.add(alpha, y, omega, z);

      if (additional_data.exact_residual)
        {
          r.add(-omega, t);
          res = criterion(A, x, b, t);
        }
      else
        res = std::sqrt(real_type(r.add_and_dot(-omega, t, r)));

      state = this->iteration_status(step, res, x);
      print_vectors(step, x, r, y);
    }
  while (state == SolverControl::iterate);

  return IterationResult(false, state, step, res);
}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
template <typename MatrixType, typename PreconditionerType>
DEAL_II_CXX20_REQUIRES(
  (concepts::is_linear_operator_on<MatrixType, VectorType> &&
   concepts::is_linear_operator_on<PreconditionerType, VectorType>))
void SolverBicgstab<VectorType>::solve(const MatrixType         &A,
                                       VectorType               &x,
                                       const VectorType         &b,
                                       const PreconditionerType &preconditioner)
{
  LogStream::Prefix prefix("Bicgstab");

  IterationResult state(false, SolverControl::failure, 0, 0);
  do
    {
      state = iterate(A, x, b, preconditioner, state.last_step);
    }
  while (state.state == SolverControl::iterate);

  // In case of failure: throw exception
  AssertThrow(state.state == SolverControl::success,
              SolverControl::NoConvergence(state.last_step,
                                           state.last_residual));
  // Otherwise exit as normal
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2000 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_solver_minres_h
#define dealii_solver_minres_h


#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/signaling_nan.h>
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
 * Minimal residual method for symmetric matrices.
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
 * However, since the MinRes method does not need additional data, the
 * respective structure is empty and does not offer any functionality. The
 * constructor has a default argument, so you may call it without the
 * additional parameter.
 *
 * The preconditioner has to be positive definite and symmetric
 *
 * The algorithm is taken from the Master thesis of Astrid Battermann
 * @cite Battermann1996 with some changes.
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
class SolverMinRes : public SolverBase<VectorType>
{
public:
  /**
   * Standardized data struct to pipe additional data to the solver. This
   * solver does not need additional data yet.
   */
  struct AdditionalData
  {};

  /**
   * Constructor.
   */
  SolverMinRes(SolverControl            &cn,
               VectorMemory<VectorType> &mem,
               const AdditionalData     &data = AdditionalData());

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  SolverMinRes(SolverControl        &cn,
               const AdditionalData &data = AdditionalData());

  /**
   * Virtual destructor.
   */
  virtual ~SolverMinRes() override = default;

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
   * @addtogroup Exceptions
   * @{
   */

  /**
   * Exception
   */
  DeclExceptionMsg(ExcPreconditionerNotDefinite,
                   "The preconditioner for MinRes must be a symmetric and "
                   "definite operator, even though MinRes can solve linear "
                   "systems with symmetric and *indefinite* operators. "
                   "During iterations, MinRes has detected that the "
                   "preconditioner is apparently not definite.");
  /** @} */

protected:
  /**
   * Implementation of the computation of the norm of the residual.
   */
  virtual double
  criterion();

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
   * Within the iteration loop, the square of the residual vector is stored in
   * this variable. The function @p criterion uses this variable to compute
   * the convergence value, which in this class is the norm of the residual
   * vector and thus the square root of the @p res2 value.
   */
  double res2;
};

/** @} */
/*------------------------- Implementation ----------------------------*/

#ifndef DOXYGEN

template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
SolverMinRes<VectorType>::SolverMinRes(SolverControl            &cn,
                                       VectorMemory<VectorType> &mem,
                                       const AdditionalData &)
  : SolverBase<VectorType>(cn, mem)
  , res2(numbers::signaling_nan<double>())
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
SolverMinRes<VectorType>::SolverMinRes(SolverControl &cn,
                                       const AdditionalData &)
  : SolverBase<VectorType>(cn)
  , res2(numbers::signaling_nan<double>())
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
double SolverMinRes<VectorType>::criterion()
{
  return res2;
}


template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
void SolverMinRes<VectorType>::print_vectors(const unsigned int,
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
void SolverMinRes<VectorType>::solve(const MatrixType         &A,
                                     VectorType               &x,
                                     const VectorType         &b,
                                     const PreconditionerType &preconditioner)
{
  LogStream::Prefix prefix("minres");

  // Memory allocation
  typename VectorMemory<VectorType>::Pointer Vu0(this->memory);
  typename VectorMemory<VectorType>::Pointer Vu1(this->memory);
  typename VectorMemory<VectorType>::Pointer Vu2(this->memory);

  typename VectorMemory<VectorType>::Pointer Vm0(this->memory);
  typename VectorMemory<VectorType>::Pointer Vm1(this->memory);
  typename VectorMemory<VectorType>::Pointer Vm2(this->memory);

  typename VectorMemory<VectorType>::Pointer Vv(this->memory);

  // define some aliases for simpler access
  using vecptr     = VectorType *;
  vecptr      u[3] = {Vu0.get(), Vu1.get(), Vu2.get()};
  vecptr      m[3] = {Vm0.get(), Vm1.get(), Vm2.get()};
  VectorType &v    = *Vv;

  // resize the vectors, but do not set the values since they'd be overwritten
  // soon anyway.
  u[0]->reinit(b, true);
  u[1]->reinit(b, true);
  u[2]->reinit(b, true);
  m[0]->reinit(b, true);
  m[1]->reinit(b, true);
  m[2]->reinit(b, true);
  v.reinit(b, true);

  // some values needed
  double delta[3] = {0, 0, 0};
  double f[2]     = {0, 0};
  double e[2]     = {0, 0};

  double r_l2 = 0;
  double r0   = 0;
  double tau  = 0;
  double c    = 0;
  double s    = 0;
  double d_   = 0;

  // The iteration step.
  unsigned int j = 1;


  // Start of the solution process
  A.vmult(*m[0], x);
  *u[1] = b;
  *u[1] -= *m[0];
  // Precondition is applied.
  // The preconditioner has to be
  // positive definite and symmetric

  // M v = u[1]
  preconditioner.vmult(v, *u[1]);

  delta[1] = v * (*u[1]);
  // Preconditioner positive
  Assert(delta[1] >= 0, ExcPreconditionerNotDefinite());

  r0   = std::sqrt(delta[1]);
  r_l2 = r0;


  u[0]->reinit(b);
  delta[0] = 1.;
  m[0]->reinit(b);
  m[1]->reinit(b);
  m[2]->reinit(b);

  SolverControl::State conv = this->iteration_status(0, r_l2, x);
  while (conv == SolverControl::iterate)
    {
      if (delta[1] != 0)
        v *= 1. / std::sqrt(delta[1]);
      else
        v.reinit(b);

      A.vmult(*u[2], v);
      u[2]->add(-std::sqrt(delta[1] / delta[0]), *u[0]);

      const double gamma = *u[2] * v;
      u[2]->add(-gamma / std::sqrt(delta[1]), *u[1]);
      *m[0] = v;

      // precondition: solve M v = u[2]
      // Preconditioner has to be positive
      // definite and symmetric.
      preconditioner.vmult(v, *u[2]);

      delta[2] = v * (*u[2]);

      Assert(delta[2] >= 0, ExcPreconditionerNotDefinite());

      if (j == 1)
        {
          d_   = gamma;
          e[1] = std::sqrt(delta[2]);
        }
      if (j > 1)
        {
          d_   = s * e[0] - c * gamma;
          e[0] = c * e[0] + s * gamma;
          f[1] = s * std::sqrt(delta[2]);
          e[1] = -c * std::sqrt(delta[2]);
        }

      const double d = std::sqrt(d_ * d_ + delta[2]);

      if (j > 1)
        tau *= s / c;
      c = d_ / d;
      tau *= c;

      s = std::sqrt(delta[2]) / d;

      if (j == 1)
        tau = r0 * c;

      m[0]->add(-e[0], *m[1]);
      if (j > 1)
        m[0]->add(-f[0], *m[2]);
      *m[0] *= 1. / d;
      x.add(tau, *m[0]);
      r_l2 *= std::fabs(s);

      conv = this->iteration_status(j, r_l2, x);

      // next iteration step
      ++j;
      // All vectors have to be shifted
      // one iteration step.
      // This should be changed one time.
      swap(*m[2], *m[1]);
      swap(*m[1], *m[0]);

      // likewise, but reverse direction:
      //   u[0] = u[1];
      //   u[1] = u[2];
      swap(*u[0], *u[1]);
      swap(*u[1], *u[2]);

      // these are scalars, so need
      // to bother
      f[0]     = f[1];
      e[0]     = e[1];
      delta[0] = delta[1];
      delta[1] = delta[2];
    }

  // in case of failure: throw exception
  AssertThrow(conv == SolverControl::success,
              SolverControl::NoConvergence(j, r_l2));

  // otherwise exit as normal
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

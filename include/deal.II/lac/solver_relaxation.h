// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_solver_relaxation_h
#define dealii_solver_relaxation_h


#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Implementation of an iterative solver based on relaxation methods. The
 * stopping criterion is the norm of the residual.
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
 * AdditionalData of this class currently does not contain any data.
 *
 *
 * <h3>Observing the progress of linear solver iterations</h3>
 *
 * The solve() function of this class uses the mechanism described in the
 * Solver base class to determine convergence. This mechanism can also be used
 * to observe the progress of the iteration.
 *
 *
 * @ingroup Solvers
 */
template <typename VectorType = Vector<double>>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
class SolverRelaxation : public SolverBase<VectorType>
{
public:
  /**
   * Standardized data struct to pipe additional data to the solver. There is
   * no data in here for relaxation methods.
   */
  struct AdditionalData
  {};

  /**
   * Constructor.
   */
  SolverRelaxation(SolverControl        &cn,
                   const AdditionalData &data = AdditionalData());

  /**
   * Solve the system $Ax = b$ using the relaxation method $x_{k+1} =
   * R(x_k,b)$. The matrix <i>A</i> itself is only used to compute the
   * residual.
   */
  template <typename MatrixType, typename RelaxationType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_linear_operator_on<MatrixType, VectorType> &&
     requires(const RelaxationType &R, VectorType &a, VectorType &b) {
       R.step(a, b);
     }))
  void solve(const MatrixType     &A,
             VectorType           &x,
             const VectorType     &b,
             const RelaxationType &R);
};

//----------------------------------------------------------------------//

template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
SolverRelaxation<VectorType>::SolverRelaxation(SolverControl &cn,
                                               const AdditionalData &)
  : SolverBase<VectorType>(cn)
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
template <typename MatrixType, typename RelaxationType>
DEAL_II_CXX20_REQUIRES(
  (concepts::is_linear_operator_on<MatrixType, VectorType> &&
   requires(const RelaxationType &R, VectorType &a, VectorType &b) {
     R.step(a, b);
   }))
void SolverRelaxation<VectorType>::solve(const MatrixType     &A,
                                         VectorType           &x,
                                         const VectorType     &b,
                                         const RelaxationType &R)
{
  GrowingVectorMemory<VectorType> mem;
  SolverControl::State            conv = SolverControl::iterate;

  // Memory allocation
  typename VectorMemory<VectorType>::Pointer Vr(mem);
  VectorType                                &r = *Vr;
  r.reinit(x);
  typename VectorMemory<VectorType>::Pointer Vd(mem);
  VectorType                                &d = *Vd;
  d.reinit(x);

  LogStream::Prefix prefix("Relaxation");

  int iter = 0;
  // Main loop
  for (; conv == SolverControl::iterate; ++iter)
    {
      // Compute residual
      A.vmult(r, x);
      r.sadd(-1., 1., b);

      // The required norm of the
      // (preconditioned)
      // residual is computed in
      // criterion() and stored
      // in res.
      conv = this->iteration_status(iter, r.l2_norm(), x);
      if (conv != SolverControl::iterate)
        break;
      R.step(x, b);
    }

  // in case of failure: throw exception
  AssertThrow(conv == SolverControl::success,
              SolverControl::NoConvergence(iter, r.l2_norm()));
  // otherwise exit as normal
}


DEAL_II_NAMESPACE_CLOSE

#endif

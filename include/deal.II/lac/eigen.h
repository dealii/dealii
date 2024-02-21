// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2000 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_eigen_h
#define dealii_eigen_h


#include <deal.II/base/config.h>

#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/vector_memory.h>

#include <cmath>
#include <limits>

DEAL_II_NAMESPACE_OPEN


/**
 * @addtogroup Solvers
 * @{
 */

/**
 * Power method (von Mises) for eigenvalue computations.
 *
 * This method determines the largest eigenvalue of a matrix by applying
 * increasing powers of this matrix to a vector. If there is an eigenvalue $l$
 * with dominant absolute value, the iteration vectors will become aligned to
 * its eigenspace and $Ax = lx$.
 *
 * A shift parameter allows to shift the spectrum, so it is possible to
 * compute the smallest eigenvalue, too.
 *
 * Convergence of this method is known to be slow.
 */
template <typename VectorType = Vector<double>>
class EigenPower : private SolverBase<VectorType>
{
public:
  /**
   * Declare type of container size.
   */
  using size_type = types::global_dof_index;

  /**
   * Standardized data struct to pipe additional data to the solver.
   */
  struct AdditionalData
  {
    /**
     * Shift parameter. This parameter allows to shift the spectrum to compute
     * a different eigenvalue.
     */
    double shift;
    /**
     * Constructor. Set the shift parameter.
     */
    AdditionalData(const double shift = 0.)
      : shift(shift)
    {}
  };

  /**
   * Constructor.
   */
  EigenPower(SolverControl            &cn,
             VectorMemory<VectorType> &mem,
             const AdditionalData     &data = AdditionalData());


  /**
   * Power method. @p x is the (not necessarily normalized, but nonzero) start
   * vector for the power method. After the iteration, @p value is the
   * approximated eigenvalue and @p x is the corresponding eigenvector,
   * normalized with respect to the l2-norm.
   */
  template <typename MatrixType>
  void
  solve(double &value, const MatrixType &A, VectorType &x);

protected:
  /**
   * Shift parameter.
   */
  AdditionalData additional_data;
};

/**
 * Inverse iteration (Wieland) for eigenvalue computations.
 *
 * This class implements an adaptive version of the inverse iteration by
 * Wieland.
 *
 * There are two choices for the stopping criterion: by default, the norm of
 * the residual $A x - l x$ is computed. Since this might not converge to zero
 * for non-symmetric matrices with non-trivial Jordan blocks, it can be
 * replaced by checking the difference of successive eigenvalues. Use
 * AdditionalData::use_residual for switching this option.
 *
 * Usually, the initial guess entering this method is updated after each step,
 * replacing it with the new approximation of the eigenvalue. Using a
 * parameter AdditionalData::relaxation between 0 and 1, this update can be
 * damped. With relaxation parameter 0, no update is performed. This damping
 * allows for slower adaption of the shift value to make sure that the method
 * converges to the eigenvalue closest to the initial guess. This can be aided
 * by the parameter AdditionalData::start_adaption, which indicates the first
 * iteration step in which the shift value should be adapted.
 */
template <typename VectorType = Vector<double>>
class EigenInverse : private SolverBase<VectorType>
{
public:
  /**
   * Declare type of container size.
   */
  using size_type = types::global_dof_index;

  /**
   * Standardized data struct to pipe additional data to the solver.
   */
  struct AdditionalData
  {
    /**
     * Damping of the updated shift value.
     */
    double relaxation;

    /**
     * Start step of adaptive shift parameter.
     */
    unsigned int start_adaption;
    /**
     * Flag for the stopping criterion.
     */
    bool use_residual;
    /**
     * Constructor.
     */
    AdditionalData(double       relaxation     = 1.,
                   unsigned int start_adaption = 6,
                   bool         use_residual   = true)
      : relaxation(relaxation)
      , start_adaption(start_adaption)
      , use_residual(use_residual)
    {}
  };

  /**
   * Constructor.
   */
  EigenInverse(SolverControl            &cn,
               VectorMemory<VectorType> &mem,
               const AdditionalData     &data = AdditionalData());

  /**
   * Inverse method. @p value is the start guess for the eigenvalue and @p x
   * is the (not necessarily normalized, but nonzero) start vector for the
   * power method. After the iteration, @p value is the approximated
   * eigenvalue and @p x is the corresponding eigenvector, normalized with
   * respect to the l2-norm.
   */
  template <typename MatrixType>
  void
  solve(double &value, const MatrixType &A, VectorType &x);

protected:
  /**
   * Flags for execution.
   */
  AdditionalData additional_data;
};

/** @} */
//---------------------------------------------------------------------------


template <typename VectorType>
EigenPower<VectorType>::EigenPower(SolverControl            &cn,
                                   VectorMemory<VectorType> &mem,
                                   const AdditionalData     &data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
{}



template <typename VectorType>
template <typename MatrixType>
void
EigenPower<VectorType>::solve(double &value, const MatrixType &A, VectorType &x)
{
  SolverControl::State conv = SolverControl::iterate;

  LogStream::Prefix prefix("Power method");

  typename VectorMemory<VectorType>::Pointer Vy(this->memory);
  VectorType                                &y = *Vy;
  y.reinit(x);
  typename VectorMemory<VectorType>::Pointer Vr(this->memory);
  VectorType                                &r = *Vr;
  r.reinit(x);

  double length     = x.l2_norm();
  double old_length = 0.;
  x *= 1. / length;

  A.vmult(y, x);

  // Main loop
  int iter = 0;
  for (; conv == SolverControl::iterate; ++iter)
    {
      y.add(additional_data.shift, x);

      // Compute absolute value of eigenvalue
      old_length = length;
      length     = y.l2_norm();

      // do a little trick to compute the sign
      // with not too much effect of round-off errors.
      double    entry  = 0.;
      size_type i      = 0;
      double    thresh = length / x.size();
      do
        {
          Assert(i < x.size(), ExcInternalError());
          entry = y(i++);
        }
      while (std::fabs(entry) < thresh);

      --i;

      // Compute unshifted eigenvalue
      value = (entry * x(i) < 0.) ? -length : length;
      value -= additional_data.shift;

      // Update normalized eigenvector
      x.equ(1 / length, y);

      // Compute residual
      A.vmult(y, x);

      // Check the change of the eigenvalue
      // Brrr, this is not really a good criterion
      conv = this->iteration_status(iter,
                                    std::fabs(1. / length - 1. / old_length),
                                    x);
    }

  // in case of failure: throw exception
  AssertThrow(conv == SolverControl::success,
              SolverControl::NoConvergence(
                iter, std::fabs(1. / length - 1. / old_length)));

  // otherwise exit as normal
}

//---------------------------------------------------------------------------

template <typename VectorType>
EigenInverse<VectorType>::EigenInverse(SolverControl            &cn,
                                       VectorMemory<VectorType> &mem,
                                       const AdditionalData     &data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
{}



template <typename VectorType>
template <typename MatrixType>
void
EigenInverse<VectorType>::solve(double           &value,
                                const MatrixType &A,
                                VectorType       &x)
{
  LogStream::Prefix prefix("Wielandt");

  SolverControl::State conv = SolverControl::iterate;

  // Prepare matrix for solver
  auto   A_op          = linear_operator(A);
  double current_shift = -value;
  auto   A_s           = A_op + current_shift * identity_operator(A_op);

  // Define solver
  ReductionControl        inner_control(5000, 1.e-16, 1.e-5, false, false);
  PreconditionIdentity    prec;
  SolverGMRES<VectorType> solver(inner_control, this->memory);

  // Next step for recomputing the shift
  unsigned int goal = additional_data.start_adaption;

  // Auxiliary vector
  typename VectorMemory<VectorType>::Pointer Vy(this->memory);
  VectorType                                &y = *Vy;
  y.reinit(x);
  typename VectorMemory<VectorType>::Pointer Vr(this->memory);
  VectorType                                &r = *Vr;
  r.reinit(x);

  double length    = x.l2_norm();
  double old_value = value;

  x *= 1. / length;

  // Main loop
  double    res  = std::numeric_limits<double>::lowest();
  size_type iter = 0;
  for (; conv == SolverControl::iterate; ++iter)
    {
      solver.solve(A_s, y, x, prec);

      // Compute absolute value of eigenvalue
      length = y.l2_norm();

      // do a little trick to compute the sign
      // with not too much effect of round-off errors.
      double    entry  = 0.;
      size_type i      = 0;
      double    thresh = length / x.size();
      do
        {
          Assert(i < x.size(), ExcInternalError());
          entry = y(i++);
        }
      while (std::fabs(entry) < thresh);

      --i;

      // Compute unshifted eigenvalue
      value = (entry * x(i) < 0. ? -1. : 1.) / length - current_shift;

      if (iter == goal)
        {
          const auto  &relaxation = additional_data.relaxation;
          const double new_shift =
            relaxation * (-value) + (1. - relaxation) * current_shift;

          A_s           = A_op + new_shift * identity_operator(A_op);
          current_shift = new_shift;

          ++goal;
        }

      // Update normalized eigenvector
      x.equ(1. / length, y);
      // Compute residual
      if (additional_data.use_residual)
        {
          y.equ(value, x);
          A.vmult(r, x);
          r.sadd(-1., value, x);
          res = r.l2_norm();
          // Check the residual
          conv = this->iteration_status(iter, res, x);
        }
      else
        {
          res  = std::fabs(1. / value - 1. / old_value);
          conv = this->iteration_status(iter, res, x);
        }
      old_value = value;
    }

  // in case of failure: throw
  // exception
  AssertThrow(conv == SolverControl::success,
              SolverControl::NoConvergence(iter, res));
  // otherwise exit as normal
}

DEAL_II_NAMESPACE_CLOSE

#endif

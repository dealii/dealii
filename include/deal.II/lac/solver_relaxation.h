// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2021 by the deal.II authors
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

#ifndef dealii_solver_relaxation_h
#define dealii_solver_relaxation_h


#include <deal.II/base/config.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/subscriptor.h>

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
class SolverRelaxation : public SolverBase<VectorType>
{
public:
  /**
   * Standardized data struct to pipe additional data to the solver. There is
   * no data in here for relaxation methods.
   */
  struct AdditionalData
  {
    /**
     * Constructor.
     */
    AdditionalData(const double relaxation = 1.0);

    /**
     * Relaxation parameter.
     */
    double relaxation;
  };

  /**
   * Constructor.
   */
  SolverRelaxation(SolverControl &       cn,
                   const AdditionalData &data = AdditionalData());

  /**
   * Solve the system $Ax = b$ using the relaxation method $x_{k+1} =
   * preconditioner(x_k,b)$. The matrix @p A itself is only used to
   * compute the residual.
   */
  template <typename MatrixType, class PreconditionerType>
  void
  solve(const MatrixType &        A,
        VectorType &              x,
        const VectorType &        b,
        const PreconditionerType &preconditioner);

private:
  /**
   * Relaxation parameter.
   */
  const double relaxation;
};

//----------------------------------------------------------------------//

template <class VectorType>
SolverRelaxation<VectorType>::SolverRelaxation(SolverControl &       cn,
                                               const AdditionalData &data)
  : SolverBase<VectorType>(cn)
  , relaxation(data.relaxation)
{}



namespace internal
{
  namespace SolverRelaxation
  {
    template <typename T, typename VectorType>
    using step_t = decltype(
      std::declval<T const>().step(std::declval<VectorType &>(),
                                   std::declval<const VectorType &>()));

    template <typename T, typename VectorType>
    constexpr bool has_step = is_supported_operation<step_t, T, VectorType>;

    template <class VectorType,
              typename MatrixType,
              class PreconditionerType,
              std::enable_if_t<has_step<PreconditionerType, VectorType>,
                               PreconditionerType> * = nullptr>
    void
    step(const MatrixType &,
         const PreconditionerType &preconditioner,
         VectorType &              x,
         const VectorType &        b,
         const double              relaxation,
         VectorType &,
         VectorType &)
    {
      AssertThrow(relaxation == 1.0, ExcNotImplemented());

      preconditioner.step(x, b);
    }

    template <class VectorType,
              typename MatrixType,
              class PreconditionerType,
              std::enable_if_t<!has_step<PreconditionerType, VectorType>,
                               PreconditionerType> * = nullptr>
    void
    step(const MatrixType &        A,
         const PreconditionerType &preconditioner,
         VectorType &              x,
         const VectorType &        b,
         const double              relaxation,
         VectorType &              r,
         VectorType &              d)
    {
      A.vmult(r, x);
      r.sadd(-1.0, 1.0, b);

      preconditioner.vmult(d, r);
      x.add(relaxation, d);
    }

  } // namespace SolverRelaxation
} // namespace internal



template <class VectorType>
SolverRelaxation<VectorType>::AdditionalData::AdditionalData(
  const double relaxation)
  : relaxation(relaxation)
{}



template <class VectorType>
template <typename MatrixType, class PreconditionerType>
void
SolverRelaxation<VectorType>::solve(const MatrixType &        A,
                                    VectorType &              x,
                                    const VectorType &        b,
                                    const PreconditionerType &preconditioner)
{
  GrowingVectorMemory<VectorType> mem;
  SolverControl::State            conv = SolverControl::iterate;

  // Memory allocation
  typename VectorMemory<VectorType>::Pointer Vr(mem);
  VectorType &                               r = *Vr;
  r.reinit(x);
  typename VectorMemory<VectorType>::Pointer Vd(mem);
  VectorType &                               d = *Vd;
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

      internal::SolverRelaxation::step(
        A, preconditioner, x, b, relaxation, r, d);
    }

  // in case of failure: throw exception
  AssertThrow(conv == SolverControl::success,
              SolverControl::NoConvergence(iter, r.l2_norm()));
  // otherwise exit as normal
}


DEAL_II_NAMESPACE_CLOSE

#endif

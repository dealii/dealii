// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2017 by the deal.II authors
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

#ifndef dealii_solver_richardson_h
#define dealii_solver_richardson_h


#include <deal.II/base/config.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/signaling_nan.h>

#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup Solvers */
/*@{*/

/**
 * Implementation of the preconditioned Richardson iteration method. The
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
 *
 * For the Richardson method, the additional data is the damping parameter,
 * which is the only content of the @p AdditionalData structure. By default,
 * the constructor of the structure sets it to one.
 *
 *
 * <h3>Observing the progress of linear solver iterations</h3>
 *
 * The solve() function of this class uses the mechanism described in the
 * Solver base class to determine convergence. This mechanism can also be used
 * to observe the progress of the iteration.
 *
 *
 * @author Ralf Hartmann
 */
template <class VectorType = Vector<double>>
class SolverRichardson : public Solver<VectorType>
{
public:
  /**
   * Standardized data struct to pipe additional data to the solver.
   */
  struct AdditionalData
  {
    /**
     * Constructor. By default, set the damping parameter to one.
     */
    explicit AdditionalData(const double omega                       = 1,
                            const bool   use_preconditioned_residual = false);

    /**
     * Relaxation parameter.
     */
    double omega;

    /**
     * Parameter for stopping criterion.
     */
    bool use_preconditioned_residual;
  };

  /**
   * Constructor.
   */
  SolverRichardson(SolverControl &           cn,
                   VectorMemory<VectorType> &mem,
                   const AdditionalData &    data = AdditionalData());

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  SolverRichardson(SolverControl &       cn,
                   const AdditionalData &data = AdditionalData());

  /**
   * Virtual destructor.
   */
  virtual ~SolverRichardson() override = default;

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
   * Solve $A^Tx=b$ for $x$.
   */
  template <typename MatrixType, typename PreconditionerType>
  void
  Tsolve(const MatrixType &        A,
         VectorType &              x,
         const VectorType &        b,
         const PreconditionerType &preconditioner);

  /**
   * Set the damping-coefficient. Default is 1., i.e. no damping.
   */
  void
  set_omega(const double om = 1.);

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

protected:
  /**
   * Implementation of the computation of the norm of the residual.
   * Depending on the flags given to the solver, the default
   * implementation of this function uses either the actual
   * residual, @p r, or the preconditioned residual, @p d.
   */
  virtual typename VectorType::value_type
  criterion(const VectorType &r, const VectorType &d) const;

  /**
   * Control parameters.
   */
  AdditionalData additional_data;
};

/*@}*/
/*----------------- Implementation of the Richardson Method ------------------*/

#ifndef DOXYGEN

template <class VectorType>
inline SolverRichardson<VectorType>::AdditionalData::AdditionalData(
  const double omega,
  const bool   use_preconditioned_residual)
  : omega(omega)
  , use_preconditioned_residual(use_preconditioned_residual)
{}


template <class VectorType>
SolverRichardson<VectorType>::SolverRichardson(SolverControl &           cn,
                                               VectorMemory<VectorType> &mem,
                                               const AdditionalData &    data)
  : Solver<VectorType>(cn, mem)
  , additional_data(data)
{}



template <class VectorType>
SolverRichardson<VectorType>::SolverRichardson(SolverControl &       cn,
                                               const AdditionalData &data)
  : Solver<VectorType>(cn)
  , additional_data(data)
{}



template <class VectorType>
template <typename MatrixType, typename PreconditionerType>
void
SolverRichardson<VectorType>::solve(const MatrixType &        A,
                                    VectorType &              x,
                                    const VectorType &        b,
                                    const PreconditionerType &preconditioner)
{
  SolverControl::State conv = SolverControl::iterate;

  double last_criterion = -std::numeric_limits<double>::max();

  unsigned int iter = 0;

  // Memory allocation.
  // 'Vr' holds the residual, 'Vd' the preconditioned residual
  typename VectorMemory<VectorType>::Pointer Vr(this->memory);
  typename VectorMemory<VectorType>::Pointer Vd(this->memory);

  VectorType &r = *Vr;
  r.reinit(x);

  VectorType &d = *Vd;
  d.reinit(x);

  LogStream::Prefix prefix("Richardson");

  // Main loop
  while (conv == SolverControl::iterate)
    {
      // Do not use residual,
      // but do it in 2 steps
      A.vmult(r, x);
      r.sadd(-1., 1., b);
      preconditioner.vmult(d, r);

      // get the required norm of the (possibly preconditioned)
      // residual
      last_criterion = criterion(r, d);
      conv           = this->iteration_status(iter, last_criterion, x);
      if (conv != SolverControl::iterate)
        break;

      x.add(additional_data.omega, d);
      print_vectors(iter, x, r, d);

      ++iter;
    }

  // in case of failure: throw exception
  if (conv != SolverControl::success)
    AssertThrow(false, SolverControl::NoConvergence(iter, last_criterion));
  // otherwise exit as normal
}



template <class VectorType>
template <typename MatrixType, typename PreconditionerType>
void
SolverRichardson<VectorType>::Tsolve(const MatrixType &        A,
                                     VectorType &              x,
                                     const VectorType &        b,
                                     const PreconditionerType &preconditioner)
{
  SolverControl::State conv           = SolverControl::iterate;
  double               last_criterion = -std::numeric_limits<double>::max();

  unsigned int iter = 0;

  // Memory allocation.
  // 'Vr' holds the residual, 'Vd' the preconditioned residual
  typename VectorMemory<VectorType>::Pointer Vr(this->memory);
  typename VectorMemory<VectorType>::Pointer Vd(this->memory);

  VectorType &r = *Vr;
  r.reinit(x);

  VectorType &d = *Vd;
  d.reinit(x);

  LogStream::Prefix prefix("RichardsonT");

  // Main loop
  while (conv == SolverControl::iterate)
    {
      // Do not use Tresidual,
      // but do it in 2 steps
      A.Tvmult(r, x);
      r.sadd(-1., 1., b);
      preconditioner.Tvmult(d, r);

      last_criterion = criterion(r, d);
      conv           = this->iteration_status(iter, last_criterion, x);
      if (conv != SolverControl::iterate)
        break;

      x.add(additional_data.omega, d);
      print_vectors(iter, x, r, d);

      ++iter;
    }

  // in case of failure: throw exception
  if (conv != SolverControl::success)
    AssertThrow(false, SolverControl::NoConvergence(iter, last_criterion));

  // otherwise exit as normal
}


template <class VectorType>
void
SolverRichardson<VectorType>::print_vectors(const unsigned int,
                                            const VectorType &,
                                            const VectorType &,
                                            const VectorType &) const
{}



template <class VectorType>
inline typename VectorType::value_type
SolverRichardson<VectorType>::criterion(const VectorType &r,
                                        const VectorType &d) const
{
  if (!additional_data.use_preconditioned_residual)
    return r.l2_norm();
  else
    return d.l2_norm();
}


template <class VectorType>
inline void
SolverRichardson<VectorType>::set_omega(const double om)
{
  additional_data.omega = om;
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

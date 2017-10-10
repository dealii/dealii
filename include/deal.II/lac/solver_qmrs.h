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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii_solver_qmrs_h
#define dealii_solver_qmrs_h

#include <deal.II/base/config.h>
#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/base/logstream.h>
#include <cmath>
#include <deal.II/base/subscriptor.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup Solvers */
/*@{*/

/**
 * Quasi-minimal residual method for symmetric matrices.
 *
 * The QMRS method is supposed to solve symmetric indefinite linear systems
 * with symmetric, not necessarily definite preconditioners. This version of
 * QMRS is adapted from Freund/Nachtigal: Software for simplified Lanczos and
 * QMR algorithms, Appl. Num. Math. 19 (1995), pp. 319-341
 *
 * This version is for right preconditioning only, since then only the
 * preconditioner is used: left preconditioning seems to require the inverse.
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
 * However, since the QMRS method does not need additional data, the
 * respective structure is empty and does not offer any functionality. The
 * constructor has a default argument, so you may call it without the
 * additional parameter.
 *
 *
 * <h3>Observing the progress of linear solver iterations</h3>
 *
 * The solve() function of this class uses the mechanism described in the
 * Solver base class to determine convergence. This mechanism can also be used
 * to observe the progress of the iteration.
 *
 *
 * @author Guido Kanschat, 1999
 */
template <typename VectorType = Vector<double> >
class SolverQMRS : public Solver<VectorType>
{
public:
  /**
   * Standardized data struct to pipe additional data to the solver.
   *
   * There are two possibilities to compute the residual: one is an estimate
   * using the computed value @p tau. The other is exact computation using
   * another matrix vector multiplication.
   *
   * QMRS, is susceptible to breakdowns, so we need a parameter telling us,
   * which numbers are considered zero. The proper breakdown criterion is very
   * unclear, so experiments may be necessary here.
   */
  struct AdditionalData
  {
    /**
     * Constructor.
     *
     * The default is no exact residual computation and breakdown parameter
     * 1e-16.
     */
    explicit
    AdditionalData(const bool   exact_residual = false,
                   const double breakdown      = 1.e-16) :
      exact_residual(exact_residual),
      breakdown(breakdown)
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
  SolverQMRS (SolverControl            &cn,
              VectorMemory<VectorType> &mem,
              const AdditionalData     &data=AdditionalData());

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  SolverQMRS (SolverControl        &cn,
              const AdditionalData &data=AdditionalData());

  /**
   * Solve the linear system $Ax=b$ for x.
   */
  template <typename MatrixType, typename PreconditionerType>
  void
  solve (const MatrixType         &A,
         VectorType               &x,
         const VectorType         &b,
         const PreconditionerType &preconditioner);

  /**
   * Interface for derived class. This function gets the current iteration
   * vector, the residual and the update vector in each step. It can be used
   * for graphical output of the convergence history.
   */
  virtual void print_vectors (const unsigned int step,
                              const VectorType   &x,
                              const VectorType   &r,
                              const VectorType   &d) const;
protected:
  /**
   * Implementation of the computation of the norm of the residual.
   */
  virtual double criterion();

  /**
   * Within the iteration loop, the square of the residual vector is stored in
   * this variable. The function @p criterion uses this variable to compute
   * the convergence value, which in this class is the norm of the residual
   * vector and thus the square root of the @p res2 value.
   */
  double res2;

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

    IterationResult (const SolverControl::State state,
                     const double               last_residual);
  };


  /**
   * The iteration loop itself. The function returns a structure indicating
   * what happened in this function.
   */
  template <typename MatrixType, typename PreconditionerType>
  IterationResult
  iterate (const MatrixType         &A,
           VectorType &x,
           const VectorType &b,
           const PreconditionerType &preconditioner,
           VectorType &v,
           VectorType &p,
           VectorType &q,
           VectorType &t,
           VectorType &d);

  /**
   * Number of the current iteration (accumulated over restarts)
   */
  unsigned int step;
};

/*@}*/
/*------------------------- Implementation ----------------------------*/

#ifndef DOXYGEN

template <class VectorType>
SolverQMRS<VectorType>::IterationResult::IterationResult (const SolverControl::State state,
                                                          const double               last_residual)
  :
  state (state),
  last_residual (last_residual)
{}


template <class VectorType>
SolverQMRS<VectorType>::SolverQMRS (SolverControl            &cn,
                                    VectorMemory<VectorType> &mem,
                                    const AdditionalData     &data)
  :
  Solver<VectorType>(cn,mem),
  additional_data(data)
{}



template <class VectorType>
SolverQMRS<VectorType>::SolverQMRS(SolverControl        &cn,
                                   const AdditionalData &data)
  :
  Solver<VectorType>(cn),
  additional_data(data)
{}



template <class VectorType>
double
SolverQMRS<VectorType>::criterion()
{
  return std::sqrt(res2);
}



template <class VectorType>
void
SolverQMRS<VectorType>::print_vectors(const unsigned int,
                                      const VectorType &,
                                      const VectorType &,
                                      const VectorType &) const
{}



template <class VectorType>
template <typename MatrixType, typename PreconditionerType>
void
SolverQMRS<VectorType>::solve (const MatrixType         &A,
                               VectorType               &x,
                               const VectorType         &b,
                               const PreconditionerType &preconditioner)
{
  LogStream::Prefix prefix("QMRS");

  // temporary vectors, allocated through the @p VectorMemory object at the
  // start of the actual solution process and deallocated at the end.
  typename VectorMemory<VectorType>::Pointer Vv(this->memory);
  typename VectorMemory<VectorType>::Pointer Vp(this->memory);
  typename VectorMemory<VectorType>::Pointer Vq(this->memory);
  typename VectorMemory<VectorType>::Pointer Vt(this->memory);
  typename VectorMemory<VectorType>::Pointer Vd(this->memory);


  // resize the vectors, but do not set
  // the values since they'd be overwritten
  // soon anyway.
  Vv->reinit(x, true);
  Vp->reinit(x, true);
  Vq->reinit(x, true);
  Vt->reinit(x, true);

  step = 0;

  IterationResult state (SolverControl::failure,0);

  do
    {
      if (step > 0)
        deallog << "Restart step " << step << std::endl;
      state = iterate(A, x, b, preconditioner, *Vv, *Vp, *Vq, *Vt, *Vd);
    }
  while (state.state == SolverControl::iterate);

  // in case of failure: throw exception
  AssertThrow(state.state == SolverControl::success,
              SolverControl::NoConvergence (step,
                                            state.last_residual));
  // otherwise exit as normal
}



template <class VectorType>
template <typename MatrixType, typename PreconditionerType>
typename SolverQMRS<VectorType>::IterationResult
SolverQMRS<VectorType>::iterate(const MatrixType         &A,
                                VectorType &x,
                                const VectorType &b,
                                const PreconditionerType &preconditioner,
                                VectorType &v,
                                VectorType &p,
                                VectorType &q,
                                VectorType &t,
                                VectorType &d)
{
  /* Remark: the matrix A in the article is the preconditioned matrix.
   * Therefore, we have to precondition x before we compute the first residual.
   * In step 1 we replace p by q to avoid one preconditioning step.
   * There are still two steps left, making this algorithm expensive.
   */

  SolverControl::State state = SolverControl::iterate;

  int  it=0;

  double tau, rho, theta=0;
  double res;

  d.reinit(x);

  // Apply right preconditioning to x
  preconditioner.vmult(q,x);
  // Preconditioned residual
  A.vmult(v,q);
  v.sadd(-1.,1.,b);
  res = v.l2_norm();

  if (this->iteration_status(step, res, x) == SolverControl::success)
    return IterationResult(SolverControl::success, res);

  p = v;

  preconditioner.vmult(q,p);

  tau = v.norm_sqr();
  rho = q*v;

  while (state == SolverControl::iterate)
    {
      step++;
      it++;
      // Step 1
      A.vmult(t,q);
      // Step 2
      const double sigma = q*t;

//TODO:[?] Find a really good breakdown criterion. The absolute one detects breakdown instead of convergence
      if (std::fabs(sigma/rho) < additional_data.breakdown)
        return IterationResult(SolverControl::iterate, std::fabs(sigma/rho));
      // Step 3
      const double alpha = rho/sigma;

      v.add(-alpha,t);
      // Step 4
      const double theta_old = theta;
      theta = v*v/tau;
      const double psi = 1./(1.+theta);
      tau *= theta*psi;

      d.sadd(psi*theta_old, psi*alpha, p);
      x += d;

      print_vectors(step,x,v,d);
      // Step 5
      if (additional_data.exact_residual)
        {
          A.vmult(q,x);
          q.sadd(-1.,1.,b);
          res = q.l2_norm();
        }
      else
        res = std::sqrt((it+1)*tau);
      state = this->iteration_status(step,res,x);
      if ((state == SolverControl::success)
          || (state == SolverControl::failure))
        return IterationResult(state, res);
      // Step 6
      if (std::fabs(rho) < additional_data.breakdown)
        return IterationResult(SolverControl::iterate, std::fabs(rho));
      // Step 7
      const double rho_old = rho;
      preconditioner.vmult(q,v);
      rho = q*v;

      const double beta = rho/rho_old;
      p.sadd(beta,v);
      preconditioner.vmult(q,p);
    }
  return IterationResult(SolverControl::success, std::fabs(rho));
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

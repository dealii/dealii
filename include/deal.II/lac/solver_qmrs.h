// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

#ifndef __deal2__solver_qmrs_h
#define __deal2__solver_qmrs_h

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
 * The QMRS method is supposed to solve symmetric indefinite linear
 * systems with symmetric, not necessarily definite preconditioners.
 * This version of QMRS is adapted from
 * Freund/Nachtigal: Software for simplified Lanczos and QMR
 * algorithms, Appl. Num. Math. 19 (1995), pp. 319-341
 *
 * This version is for right preconditioning only, since then only the
 * preconditioner is used: left preconditioning seems to require the
 * inverse.
 *
 * For the requirements on matrices and vectors in order to work with
 * this class, see the documentation of the Solver base class.
 *
 * Like all other solver classes, this class has a local structure called
 * @p AdditionalData which is used to pass additional parameters to the
 * solver, like damping parameters or the number of temporary vectors. We
 * use this additional structure instead of passing these values directly
 * to the constructor because this makes the use of the @p SolverSelector and
 * other classes much easier and guarantees that these will continue to
 * work even if number or type of the additional parameters for a certain
 * solver changes.
 *
 * However, since the QMRS method does not need additional data, the respective
 * structure is empty and does not offer any functionality. The constructor
 * has a default argument, so you may call it without the additional
 * parameter.
 *
 * @author Guido Kanschat, 1999
 */
template <class VECTOR = Vector<double> >
class SolverQMRS : public Solver<VECTOR>
{
public:
  /**
   * Standardized data struct to
   * pipe additional data to the
   * solver.
   *
   * There are two possibilities to compute
   * the residual: one is an estimate using
   * the computed value @p tau. The other
   * is exact computation using another matrix
   * vector multiplication.
   *
   * QMRS, is susceptible to
   * breakdowns, so we need a
   * parameter telling us, which
   * numbers are considered
   * zero. The proper breakdown
   * criterion is very unclear, so
   * experiments may be necessary
   * here.
   */
  struct AdditionalData
  {
    /**
     * Constructor.
     *
     * The default is no exact residual
     * computation and breakdown
     * parameter 1e-16.
     */
    AdditionalData(bool exact_residual = false,
                   double breakdown=1.e-16) :
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
  SolverQMRS (SolverControl &cn,
              VectorMemory<VECTOR> &mem,
              const AdditionalData &data=AdditionalData());

  /**
   * Constructor. Use an object of
   * type GrowingVectorMemory as
   * a default to allocate memory.
   */
  SolverQMRS (SolverControl        &cn,
              const AdditionalData &data=AdditionalData());

  /**
   * Solve the linear system $Ax=b$
   * for x.
   */
  template<class MATRIX, class PRECONDITIONER>
  void
  solve (const MATRIX         &A,
         VECTOR               &x,
         const VECTOR         &b,
         const PRECONDITIONER &precondition);

  /**
   * Interface for derived class.
   * This function gets the current
   * iteration vector, the residual
   * and the update vector in each
   * step. It can be used for a
   * graphical output of the
   * convergence history.
   */
  virtual void print_vectors(const unsigned int step,
                             const VECTOR &x,
                             const VECTOR &r,
                             const VECTOR &d) const;
protected:
  /**
   * Implementation of the computation of
   * the norm of the residual.
   */
  virtual double criterion();

  /**
   * Temporary vectors, allocated through
   * the @p VectorMemory object at the start
   * of the actual solution process and
   * deallocated at the end.
   */
  VECTOR *Vv;
  VECTOR *Vp;
  VECTOR *Vq;
  VECTOR *Vt;
  VECTOR *Vd;
  /**
   * Iteration vector.
   */
  VECTOR *Vx;
  /**
   * RHS vector.
   */
  const VECTOR *Vb;

  /**
   * Within the iteration loop, the
   * square of the residual vector is
   * stored in this variable. The
   * function @p criterion uses this
   * variable to compute the convergence
   * value, which in this class is the
   * norm of the residual vector and thus
   * the square root of the @p res2 value.
   */
  double res2;
  /**
   * Additional parameters..
   */
  AdditionalData additional_data;
private:
  /**
   * The iteration loop itself.
   */
  template<class MATRIX, class PRECONDITIONER>
  bool
  iterate(const MATRIX &A, const PRECONDITIONER &precondition);
  /**
   * The current iteration step.
   */
  unsigned int step;
};

/*@}*/
/*------------------------- Implementation ----------------------------*/

#ifndef DOXYGEN

template<class VECTOR>
SolverQMRS<VECTOR>::SolverQMRS(SolverControl &cn,
                               VectorMemory<VECTOR> &mem,
                               const AdditionalData &data)
  :
  Solver<VECTOR>(cn,mem),
  additional_data(data)
{}



template<class VECTOR>
SolverQMRS<VECTOR>::SolverQMRS(SolverControl &cn,
                               const AdditionalData &data)
  :
  Solver<VECTOR>(cn),
  additional_data(data)
{}



template<class VECTOR>
double
SolverQMRS<VECTOR>::criterion()
{
  return std::sqrt(res2);
}



template<class VECTOR>
void
SolverQMRS<VECTOR>::print_vectors(const unsigned int,
                                  const VECTOR &,
                                  const VECTOR &,
                                  const VECTOR &) const
{}



template<class VECTOR>
template<class MATRIX, class PRECONDITIONER>
void
SolverQMRS<VECTOR>::solve (const MATRIX         &A,
                           VECTOR               &x,
                           const VECTOR         &b,
                           const PRECONDITIONER &precondition)
{
  deallog.push("QMRS");

  // Memory allocation
  Vv  = this->memory.alloc();
  Vp  = this->memory.alloc();
  Vq  = this->memory.alloc();
  Vt  = this->memory.alloc();
  Vd  = this->memory.alloc();

  Vx = &x;
  Vb = &b;
  // resize the vectors, but do not set
  // the values since they'd be overwritten
  // soon anyway.
  Vv->reinit(x, true);
  Vp->reinit(x, true);
  Vq->reinit(x, true);
  Vt->reinit(x, true);

  step = 0;

  bool state;

  do
    {
      if (step)
        deallog << "Restart step " << step << std::endl;
      state = iterate(A, precondition);
    }
  while (state);

  // Deallocate Memory
  this->memory.free(Vv);
  this->memory.free(Vp);
  this->memory.free(Vq);
  this->memory.free(Vt);
  this->memory.free(Vd);

  // Output
  deallog.pop();

  // in case of failure: throw exception
  if (this->control().last_check() != SolverControl::success)
    AssertThrow(false, SolverControl::NoConvergence (this->control().last_step(),
                                                     this->control().last_value()));
  // otherwise exit as normal
}



template<class VECTOR>
template<class MATRIX, class PRECONDITIONER>
bool
SolverQMRS<VECTOR>::iterate(const MATRIX         &A,
                            const PRECONDITIONER &precondition)
{
  /* Remark: the matrix A in the article is the preconditioned matrix.
   * Therefore, we have to precondition x before we compute the first residual.
   * In step 1 we replace p by q to avoid one preconditioning step.
   * There are still two steps left, making this algorithm expensive.
   */

  SolverControl::State state = SolverControl::iterate;

  // define some aliases for simpler access
  VECTOR &v  = *Vv;
  VECTOR &p  = *Vp;
  VECTOR &q  = *Vq;
  VECTOR &t  = *Vt;
  VECTOR &d  = *Vd;
  VECTOR &x  = *Vx;
  const VECTOR &b = *Vb;

  int  it=0;

  double tau, rho, theta=0, sigma, alpha, psi, theta_old, rho_old, beta;
  double res;

  d.reinit(x);

  // Apply right preconditioning to x
  precondition.vmult(q,x);
  // Preconditioned residual
  A.vmult(v,q);
  v.sadd(-1.,1.,b);
  res = v.l2_norm();

  if (this->control().check(step, res) == SolverControl::success)
    return false;

  p = v;

  precondition.vmult(q,p);

  tau = v.norm_sqr();
  rho = q*v;

  while (state == SolverControl::iterate)
    {
      step++;
      it++;
      // Step 1
      A.vmult(t,q);
      // Step 2
      sigma = q*t;

//TODO:[?] Find a really good breakdown criterion. The absolute one detects breakdown instead of convergence
      if (std::fabs(sigma/rho) < additional_data.breakdown)
        return true;
      // Step 3
      alpha = rho/sigma;

      v.add(-alpha,t);
      // Step 4
      theta_old = theta;
      theta = v*v/tau;
      psi = 1./(1.+theta);
      tau *= theta*psi;

      d.sadd(psi*theta_old, psi*alpha, p);
      x.add(d);

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
      state = this->control().check(step,res);
      if ((state == SolverControl::success)
          || (state == SolverControl::failure))
        return false;
      // Step 6
      if (std::fabs(rho) < additional_data.breakdown)
        return true;
      // Step 7
      rho_old = rho;
      precondition.vmult(q,v);
      rho = q*v;

      beta = rho/rho_old;
      p.sadd(beta,v);
      precondition.vmult(q,p);
    }
  return false;
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

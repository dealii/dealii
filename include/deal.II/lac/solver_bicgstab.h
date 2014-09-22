// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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

#ifndef __deal2__solver_bicgstab_h
#define __deal2__solver_bicgstab_h


#include <deal.II/base/config.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>
#include <cmath>
#include <deal.II/base/subscriptor.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup Solvers */
/*@{*/

/**
 * Bicgstab algorithm by van der Vorst.
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
 * The Bicgstab-method has two additional parameters: the first is a
 * boolean, deciding whether to compute the actual residual in each step (@p
 * true) or to use the length of the computed orthogonal residual (@p
 * false). Note that computing the residual causes a third
 * matrix-vector-multiplication, though no additional preconditioning, in
 * each step. The reason for doing this is, that the size of the
 * orthogonalized residual computed during the iteration may be larger by
 * orders of magnitude than the true residual. This is due to numerical
 * instabilities related to badly conditioned matrices. Since this
 * instability results in a bad stopping criterion, the default for this
 * parameter is @p true. Whenever the user knows that the estimated residual
 * works reasonably as well, the flag should be set to @p false in order
 * to increase the performance of the solver.
 *
 * The second parameter is the size of a breakdown criterion. It is
 * difficult to find a general good criterion, so if things do not
 * work for you, try to change this value.
 */
template <class VECTOR = Vector<double> >
class SolverBicgstab : public Solver<VECTOR>
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
     * parameter 1e-10.
     */
    AdditionalData(const bool   exact_residual = true,
                   const double breakdown      = 1.e-10) :
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
  SolverBicgstab (SolverControl        &cn,
                  VectorMemory<VECTOR> &mem,
                  const AdditionalData &data=AdditionalData());

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  SolverBicgstab (SolverControl        &cn,
                  const AdditionalData &data=AdditionalData());

  /**
   * Virtual destructor.
   */
  virtual ~SolverBicgstab ();

  /**
   * Solve primal problem only.
   */
  template<class MATRIX, class PRECONDITIONER>
  void
  solve (const MATRIX &A,
         VECTOR       &x,
         const VECTOR &b,
         const PRECONDITIONER &precondition);

protected:
  /**
   * Computation of the stopping criterion.
   */
  template <class MATRIX>
  double criterion (const MATRIX &A, const VECTOR &x, const VECTOR &b);

  /**
   * Interface for derived class.  This function gets the current iteration
   * vector, the residual and the update vector in each step. It can be used
   * for a graphical output of the convergence history.
   */
  virtual void print_vectors(const unsigned int step,
                             const VECTOR &x,
                             const VECTOR &r,
                             const VECTOR &d) const;

  /**
   * Auxiliary vector.
   */
  VECTOR *Vx;
  /**
   * Auxiliary vector.
   */
  VECTOR *Vr;
  /**
   * Auxiliary vector.
   */
  VECTOR *Vrbar;
  /**
   * Auxiliary vector.
   */
  VECTOR *Vp;
  /**
   * Auxiliary vector.
   */
  VECTOR *Vy;
  /**
   * Auxiliary vector.
   */
  VECTOR *Vz;
  /**
   * Auxiliary vector.
   */
  VECTOR *Vt;
  /**
   * Auxiliary vector.
   */
  VECTOR *Vv;
  /**
   * Right hand side vector.
   */
  const VECTOR *Vb;

  /**
   * Auxiliary value.
   */
  double alpha;
  /**
   * Auxiliary value.
   */
  double beta;
  /**
   * Auxiliary value.
   */
  double omega;
  /**
   * Auxiliary value.
   */
  double rho;
  /**
   * Auxiliary value.
   */
  double rhobar;

  /**
   * Current iteration step.
   */
  unsigned int step;

  /**
   * Residual.
   */
  double res;

  /**
   * Additional parameters.
   */
  AdditionalData additional_data;

private:
  /**
   * Everything before the iteration loop.
   */
  template <class MATRIX>
  SolverControl::State start(const MATRIX &A);

  /**
   * The iteration loop itself.
   */
  template<class MATRIX, class PRECONDITIONER>
  bool
  iterate(const MATRIX &A, const PRECONDITIONER &precondition);

};

/*@}*/
/*-------------------------Inline functions -------------------------------*/

#ifndef DOXYGEN

template<class VECTOR>
SolverBicgstab<VECTOR>::SolverBicgstab (SolverControl &cn,
                                        VectorMemory<VECTOR> &mem,
                                        const AdditionalData &data)
  :
  Solver<VECTOR>(cn,mem),
  additional_data(data)
{}



template<class VECTOR>
SolverBicgstab<VECTOR>::SolverBicgstab (SolverControl &cn,
                                        const AdditionalData &data)
  :
  Solver<VECTOR>(cn),
  additional_data(data)
{}



template<class VECTOR>
SolverBicgstab<VECTOR>::~SolverBicgstab ()
{}



template <class VECTOR>
template <class MATRIX>
double
SolverBicgstab<VECTOR>::criterion (const MATRIX &A, const VECTOR &x, const VECTOR &b)
{
  A.vmult(*Vt, x);
  Vt->add(-1.,b);
  res = Vt->l2_norm();

  return res;
}



template <class VECTOR >
template <class MATRIX>
SolverControl::State
SolverBicgstab<VECTOR>::start(const MATRIX &A)
{
  A.vmult(*Vr, *Vx);
  Vr->sadd(-1.,1.,*Vb);
  res = Vr->l2_norm();

  return this->control().check(step, res);
}



template<class VECTOR>
void
SolverBicgstab<VECTOR>::print_vectors(const unsigned int,
                                      const VECTOR &,
                                      const VECTOR &,
                                      const VECTOR &) const
{}



template<class VECTOR>
template<class MATRIX, class PRECONDITIONER>
bool
SolverBicgstab<VECTOR>::iterate(const MATRIX &A,
                                const PRECONDITIONER &precondition)
{
//TODO:[GK] Implement "use the length of the computed orthogonal residual" in the BiCGStab method.
  SolverControl::State state = SolverControl::iterate;
  alpha = omega = rho = 1.;

  VECTOR &r = *Vr;
  VECTOR &rbar = *Vrbar;
  VECTOR &p = *Vp;
  VECTOR &y = *Vy;
  VECTOR &z = *Vz;
  VECTOR &t = *Vt;
  VECTOR &v = *Vv;

  rbar = r;
  bool startup = true;

  do
    {
      ++step;

      rhobar = r*rbar;
      beta   = rhobar * alpha / (rho * omega);
      rho    = rhobar;
      if (startup == true)
        {
          p = r;
          startup = false;
        }
      else
        p.sadd(beta, 1., r, -beta*omega, v);

      precondition.vmult(y,p);
      A.vmult(v,y);
      rhobar = rbar * v;

      alpha = rho/rhobar;

//TODO:[?] Find better breakdown criterion

      if (std::fabs(alpha) > 1.e10)
        return true;

      r.add(-alpha, v);

      // check for early success, see the lac/bicgstab_early testcase as to
      // why this is necessary
      if (this->control().check(step, r.l2_norm()) == SolverControl::success)
        {
          Vx->add(alpha, y);
          print_vectors(step, *Vx, r, y);
          return false;
        }

      precondition.vmult(z,r);
      A.vmult(t,z);
      rhobar = t*r;
      omega = rhobar/(t*t);
      Vx->add(alpha, y, omega, z);
      r.add(-omega, t);

      if (additional_data.exact_residual)
        res = criterion(A, *Vx, *Vb);
      else
        res = r.l2_norm();

      state = this->control().check(step, res);
      print_vectors(step, *Vx, r, y);
    }
  while (state == SolverControl::iterate);
  return false;
}


template<class VECTOR>
template<class MATRIX, class PRECONDITIONER>
void
SolverBicgstab<VECTOR>::solve(const MATRIX &A,
                              VECTOR       &x,
                              const VECTOR &b,
                              const PRECONDITIONER &precondition)
{
  deallog.push("Bicgstab");
  Vr    = this->memory.alloc();
  Vr->reinit(x, true);
  Vrbar = this->memory.alloc();
  Vrbar->reinit(x, true);
  Vp    = this->memory.alloc();
  Vp->reinit(x, true);
  Vy    = this->memory.alloc();
  Vy->reinit(x, true);
  Vz    = this->memory.alloc();
  Vz->reinit(x, true);
  Vt    = this->memory.alloc();
  Vt->reinit(x, true);
  Vv    = this->memory.alloc();
  Vv->reinit(x, true);

  Vx = &x;
  Vb = &b;

  step = 0;

  bool state;

  do
    {
      if (step != 0)
        deallog << "Restart step " << step << std::endl;
      if (start(A) == SolverControl::success)
        break;
      state = iterate(A, precondition);
    }
  while (state);

  this->memory.free(Vr);
  this->memory.free(Vrbar);
  this->memory.free(Vp);
  this->memory.free(Vy);
  this->memory.free(Vz);
  this->memory.free(Vt);
  this->memory.free(Vv);

  deallog.pop();

  // in case of failure: throw exception
  if (this->control().last_check() != SolverControl::success)
    AssertThrow(false, SolverControl::NoConvergence (this->control().last_step(),
                                                     this->control().last_value()));
  // otherwise exit as normal
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

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

#ifndef __deal2__solver_richardson_h
#define __deal2__solver_richardson_h


#include <deal.II/base/config.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/base/subscriptor.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup Solvers */
/*@{*/

/**
 * Implementation of the preconditioned Richardson iteration method. The stopping criterion
 * is the norm of the residual.
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
 * For the Richardson method, the additional data is the damping parameter,
 * which is the only content of the @p AdditionalData structure. By default,
 * the constructor of the structure sets it to one.
 *
 * @author Ralf Hartmann
 */
template <class VECTOR = Vector<double> >
class SolverRichardson : public Solver<VECTOR>
{
public:
  /**
   * Standardized data struct to
   * pipe additional data to the
   * solver.
   */
  struct AdditionalData
  {
    /**
     * Constructor. By default,
     * set the damping parameter
     * to one.
     */
    AdditionalData (const double omega                       = 1,
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
  SolverRichardson (SolverControl        &cn,
                    VectorMemory<VECTOR> &mem,
                    const AdditionalData &data=AdditionalData());

  /**
   * Constructor. Use an object of
   * type GrowingVectorMemory as
   * a default to allocate memory.
   */
  SolverRichardson (SolverControl        &cn,
                    const AdditionalData &data=AdditionalData());

  /**
   * Virtual destructor.
   */
  virtual ~SolverRichardson ();

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
   * Solve $A^Tx=b$ for $x$.
   */
  template<class MATRIX, class PRECONDITIONER>
  void
  Tsolve (const MATRIX         &A,
          VECTOR               &x,
          const VECTOR         &b,
          const PRECONDITIONER &precondition);

  /**
   * Set the damping-coefficient.
   * Default is 1., i.e. no damping.
   */
  void set_omega (const double om=1.);

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
  virtual typename VECTOR::value_type criterion();

  /**
   * Residual. Temporary vector allocated through
   * the VectorMemory object at the start
   * of the actual solution process and
   * deallocated at the end.
   */
  VECTOR *Vr;
  /**
   * Preconditioned
   * residual. Temporary vector
   * allocated through the
   * VectorMemory object at the
   * start of the actual solution
   * process and deallocated at the
   * end.
   */
  VECTOR *Vd;

  /**
   * Control parameters.
   */
  AdditionalData additional_data;

  /**
   * Within the iteration loop, the
   * norm of the residual is
   * stored in this variable. The
   * function @p criterion uses this
   * variable to compute the convergence
   * value, which in this class is the
   * norm of the residual vector and thus
   * the square root of the @p res2 value.
   */
  typename VECTOR::value_type res;
};

/*@}*/
/*----------------- Implementation of the Richardson Method ------------------*/

#ifndef DOXYGEN

template <class VECTOR>
inline
SolverRichardson<VECTOR>::AdditionalData::
AdditionalData (const double omega,
                const bool   use_preconditioned_residual)
  :
  omega(omega),
  use_preconditioned_residual(use_preconditioned_residual)
{}


template <class VECTOR>
SolverRichardson<VECTOR>::SolverRichardson(SolverControl &cn,
                                           VectorMemory<VECTOR> &mem,
                                           const AdditionalData &data)
  :
  Solver<VECTOR> (cn,mem),
  additional_data(data)
{}



template <class VECTOR>
SolverRichardson<VECTOR>::SolverRichardson(SolverControl &cn,
                                           const AdditionalData &data)
  :
  Solver<VECTOR> (cn),
  additional_data(data)
{}



template <class VECTOR>
SolverRichardson<VECTOR>::~SolverRichardson()
{}


template <class VECTOR>
template <class MATRIX, class PRECONDITIONER>
void
SolverRichardson<VECTOR>::solve (const MATRIX         &A,
                                 VECTOR               &x,
                                 const VECTOR         &b,
                                 const PRECONDITIONER &precondition)
{
  SolverControl::State conv=SolverControl::iterate;

  // Memory allocation
  Vr  = this->memory.alloc();
  VECTOR &r  = *Vr;
  r.reinit(x);
  Vd  = this->memory.alloc();
  VECTOR &d  = *Vd;
  d.reinit(x);

  deallog.push("Richardson");

  try
    {
      // Main loop
      for (int iter=0; conv==SolverControl::iterate; iter++)
        {
          // Do not use residual,
          // but do it in 2 steps
          A.vmult(r,x);
          r.sadd(-1.,1.,b);
          precondition.vmult(d,r);

          // The required norm of the
          // (preconditioned)
          // residual is computed in
          // criterion() and stored
          // in res.
          conv = this->control().check (iter, criterion());
//        conv = this->control().check (iter, std::sqrt(A.matrix_norm_square(r)));
          if (conv != SolverControl::iterate)
            break;

          x.add(additional_data.omega,d);
          print_vectors(iter,x,r,d);
        }
    }
  catch (...)
    {
      this->memory.free(Vr);
      this->memory.free(Vd);
      deallog.pop();
      throw;
    }
  // Deallocate Memory
  this->memory.free(Vr);
  this->memory.free(Vd);
  deallog.pop();

  // in case of failure: throw exception
  if (this->control().last_check() != SolverControl::success)
    AssertThrow(false, SolverControl::NoConvergence (this->control().last_step(),
                                                     this->control().last_value()));
  // otherwise exit as normal
}


template <class VECTOR>
template <class MATRIX, class PRECONDITIONER>
void
SolverRichardson<VECTOR>::Tsolve (const MATRIX         &A,
                                  VECTOR               &x,
                                  const VECTOR         &b,
                                  const PRECONDITIONER &precondition)
{
  SolverControl::State conv=SolverControl::iterate;

  // Memory allocation
  Vr  = this->memory.alloc();
  VECTOR &r  = *Vr;
  r.reinit(x);
  Vd  =this-> memory.alloc();
  VECTOR &d  = *Vd;
  d.reinit(x);

  deallog.push("RichardsonT");

  try
    {
      // Main loop
      for (int iter=0; conv==SolverControl::iterate; iter++)
        {
          // Do not use Tresidual,
          // but do it in 2 steps
          A.Tvmult(r,x);
          r.sadd(-1.,1.,b);
          precondition.Tvmult(d,r);

          conv = this->control().check (iter, criterion());
          if (conv != SolverControl::iterate)
            break;

          x.add(additional_data.omega,d);
          print_vectors(iter,x,r,d);
        }
    }
  catch (...)
    {
      this->memory.free(Vr);
      this->memory.free(Vd);
      deallog.pop();
      throw;
    }

  // Deallocate Memory
  this->memory.free(Vr);
  this->memory.free(Vd);
  deallog.pop();
  // in case of failure: throw exception
  if (this->control().last_check() != SolverControl::success)
    AssertThrow(false, SolverControl::NoConvergence (this->control().last_step(),
                                                     this->control().last_value()));
  // otherwise exit as normal
}


template <class VECTOR>
void
SolverRichardson<VECTOR>::print_vectors(const unsigned int,
                                        const VECTOR &,
                                        const VECTOR &,
                                        const VECTOR &) const
{}



template <class VECTOR>
inline typename VECTOR::value_type
SolverRichardson<VECTOR>::criterion()
{
  if (!additional_data.use_preconditioned_residual)
    res = Vr->l2_norm();
  else
    res = Vd->l2_norm();
  return res;
}


template <class VECTOR>
inline void
SolverRichardson<VECTOR>::set_omega (const double om)
{
  additional_data.omega=om;
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

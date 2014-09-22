// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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

#ifndef __deal2__eigen_h
#define __deal2__eigen_h


#include <deal.II/base/config.h>
#include <deal.II/lac/shifted_matrix.h>
#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/precondition.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN


/*!@addtogroup Solvers */
/*@{*/

/**
 * Power method (von Mises) for eigenvalue computations.
 *
 * This method determines the largest eigenvalue of a matrix by
 * applying increasing powers of this matrix to a vector. If there is
 * an eigenvalue $l$ with dominant absolute value, the iteration vectors
 * will become aligned to its eigenspace and $Ax = lx$.
 *
 * A shift parameter allows to shift the spectrum, so it is possible
 * to compute the smallest eigenvalue, too.
 *
 * Convergence of this method is known to be slow.
 *
 * @author Guido Kanschat, 2000
 */
template <class VECTOR = Vector<double> >
class EigenPower : private Solver<VECTOR>
{
public:
  /**
   * Declare type of container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Standardized data struct to
   * pipe additional data to the
   * solver.
   */
  struct AdditionalData
  {
    /**
     * Shift parameter. This
     * parameter allows to shift
     * the spectrum to compute a
     * different eigenvalue.
     */
    double shift;
    /**
     * Constructor. Set the shift parameter.
     */
    AdditionalData (const double shift = 0.)
      :
      shift(shift)
    {}

  };

  /**
   * Constructor.
   */
  EigenPower (SolverControl &cn,
              VectorMemory<VECTOR> &mem,
              const AdditionalData &data=AdditionalData());

  /**
   * Virtual destructor.
   */
  virtual ~EigenPower ();

  /**
   * Power method. @p x is the
   * (not necessarily normalized,
   * but nonzero) start vector for
   * the power method. After the
   * iteration, @p value is the
   * approximated eigenvalue and
   * @p x is the corresponding
   * eigenvector, normalized with
   * respect to the l2-norm.
   */
  template <class MATRIX>
  void
  solve (double       &value,
         const MATRIX &A,
         VECTOR       &x);

protected:
  /**
   * Shift parameter.
   */
  AdditionalData additional_data;
};

/**
 * Inverse iteration (Wieland) for eigenvalue computations.
 *
 * This class implements an adaptive version of the inverse iteration by Wieland.
 *
 * There are two choices for the stopping criterion: by default, the
 * norm of the residual $A x - l x$ is computed. Since this might not
 * converge to zero for non-symmetric matrices with non-trivial Jordan
 * blocks, it can be replaced by checking the difference of successive
 * eigenvalues. Use AdditionalData::use_residual for switching
 * this option.
 *
 * Usually, the initial guess entering this method is updated after
 * each step, replacing it with the new approximation of the
 * eigenvalue. Using a parameter AdditionalData::relaxation
 * between 0 and 1, this update can be damped. With relaxation
 * parameter 0, no update is performed. This damping allows for slower
 * adaption of the shift value to make sure that the method converges
 * to the eigenvalue closest to the initial guess. This can be aided
 * by the parameter AdditionalData::start_adaption, which
 * indicates the first iteration step in which the shift value should
 * be adapted.
 *
 * @author Guido Kanschat, 2000, 2003
 */
template <class VECTOR = Vector<double> >
class EigenInverse : private Solver<VECTOR>
{
public:
  /**
   * Declare type of container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Standardized data struct to
   * pipe additional data to the
   * solver.
   */
  struct AdditionalData
  {
    /**
     * Damping of the updated shift value.
     */
    double relaxation;

    /**
     * Start step of adaptive
     * shift parameter.
     */
    unsigned int start_adaption;
    /**
     * Flag for the stopping criterion.
     */
    bool use_residual;
    /**
     * Constructor.
     */
    AdditionalData (double relaxation = 1.,
                    unsigned int start_adaption = 6,
                    bool use_residual = true)
      :
      relaxation(relaxation),
      start_adaption(start_adaption),
      use_residual(use_residual)
    {}

  };

  /**
   * Constructor.
   */
  EigenInverse (SolverControl &cn,
                VectorMemory<VECTOR> &mem,
                const AdditionalData &data=AdditionalData());


  /**
   * Virtual destructor.
   */
  virtual ~EigenInverse ();

  /**
   * Inverse method. @p value is
   * the start guess for the
   * eigenvalue and @p x is the
   * (not necessarily normalized,
   * but nonzero) start vector for
   * the power method. After the
   * iteration, @p value is the
   * approximated eigenvalue and
   * @p x is the corresponding
   * eigenvector, normalized with
   * respect to the l2-norm.
   */
  template <class MATRIX>
  void
  solve (double       &value,
         const MATRIX &A,
         VECTOR       &x);

protected:
  /**
   * Flags for execution.
   */
  AdditionalData additional_data;
};

/*@}*/
//---------------------------------------------------------------------------


template <class VECTOR>
EigenPower<VECTOR>::EigenPower (SolverControl &cn,
                                VectorMemory<VECTOR> &mem,
                                const AdditionalData &data)
  :
  Solver<VECTOR>(cn, mem),
  additional_data(data)
{}



template <class VECTOR>
EigenPower<VECTOR>::~EigenPower ()
{}



template <class VECTOR>
template <class MATRIX>
void
EigenPower<VECTOR>::solve (double       &value,
                           const MATRIX &A,
                           VECTOR       &x)
{
  SolverControl::State conv=SolverControl::iterate;

  deallog.push("Power method");

  VECTOR *Vy = this->memory.alloc ();
  VECTOR &y = *Vy;
  y.reinit (x);
  VECTOR *Vr = this->memory.alloc ();
  VECTOR &r = *Vr;
  r.reinit (x);

  double length = x.l2_norm ();
  double old_length = 0.;
  x.scale(1./length);

  A.vmult (y,x);

  // Main loop
  for (int iter=0; conv==SolverControl::iterate; iter++)
    {
      y.add(additional_data.shift, x);

      // Compute absolute value of eigenvalue
      old_length = length;
      length = y.l2_norm ();

      // do a little trick to compute the sign
      // with not too much effect of round-off errors.
      double entry = 0.;
      size_type i = 0;
      double thresh = length/x.size();
      do
        {
          Assert (i<x.size(), ExcInternalError());
          entry = y (i++);
        }
      while (std::fabs(entry) < thresh);

      --i;

      // Compute unshifted eigenvalue
      value = (entry * x (i) < 0.) ? -length : length;
      value -= additional_data.shift;

      // Update normalized eigenvector
      x.equ (1/length, y);

      // Compute residual
      A.vmult (y,x);

      // Check the change of the eigenvalue
      // Brrr, this is not really a good criterion
      conv = this->control().check (iter, std::fabs(1./length-1./old_length));
    }

  this->memory.free(Vy);
  this->memory.free(Vr);

  deallog.pop();

  // in case of failure: throw exception
  if (this->control().last_check() != SolverControl::success)
    AssertThrow(false, SolverControl::NoConvergence (this->control().last_step(),
                                                     this->control().last_value()));
  // otherwise exit as normal
}

//---------------------------------------------------------------------------

template <class VECTOR>
EigenInverse<VECTOR>::EigenInverse (SolverControl &cn,
                                    VectorMemory<VECTOR> &mem,
                                    const AdditionalData &data)
  :
  Solver<VECTOR>(cn, mem),
  additional_data(data)
{}



template <class VECTOR>
EigenInverse<VECTOR>::~EigenInverse ()
{}



template <class VECTOR>
template <class MATRIX>
void
EigenInverse<VECTOR>::solve (double       &value,
                             const MATRIX &A,
                             VECTOR       &x)
{
  deallog.push("Wielandt");

  SolverControl::State conv=SolverControl::iterate;

  // Prepare matrix for solver
  ShiftedMatrix <MATRIX> A_s(A, -value);

  // Define solver
  ReductionControl inner_control (5000, 1.e-16, 1.e-5, false, false);
  PreconditionIdentity prec;
  SolverGMRES<VECTOR>
  solver(inner_control, this->memory);

  // Next step for recomputing the shift
  unsigned int goal = additional_data.start_adaption;

  // Auxiliary vector
  VECTOR *Vy = this->memory.alloc ();
  VECTOR &y = *Vy;
  y.reinit (x);
  VECTOR *Vr = this->memory.alloc ();
  VECTOR &r = *Vr;
  r.reinit (x);

  double length = x.l2_norm ();
  double old_value = value;

  x.scale(1./length);

  // Main loop
  for (size_type iter=0; conv==SolverControl::iterate; iter++)
    {
      solver.solve (A_s, y, x, prec);

      // Compute absolute value of eigenvalue
      length = y.l2_norm ();

      // do a little trick to compute the sign
      // with not too much effect of round-off errors.
      double entry = 0.;
      size_type i = 0;
      double thresh = length/x.size();
      do
        {
          Assert (i<x.size(), ExcInternalError());
          entry = y (i++);
        }
      while (std::fabs(entry) < thresh);

      --i;

      // Compute unshifted eigenvalue
      value = (entry * x (i) < 0.) ? -length : length;
      value = 1./value;
      value -= A_s.shift ();

      if (iter==goal)
        {
          const double new_shift = - additional_data.relaxation * value
                                   + (1.-additional_data.relaxation) * A_s.shift();
          A_s.shift(new_shift);
          ++goal;
        }

      // Update normalized eigenvector
      x.equ (1./length, y);
      // Compute residual
      if (additional_data.use_residual)
        {
          y.equ (value, x);
          A.vmult(r,x);
          r.sadd(-1., value, x);
          double res = r.l2_norm();
          // Check the residual
          conv = this->control().check (iter, res);
        }
      else
        {
          conv = this->control().check (iter, std::fabs(1./value-1./old_value));
        }
      old_value = value;
    }

  this->memory.free(Vy);
  this->memory.free(Vr);

  deallog.pop();

  // in case of failure: throw
  // exception
  if (this->control().last_check() != SolverControl::success)
    throw SolverControl::NoConvergence (this->control().last_step(),
                                        this->control().last_value());
  // otherwise exit as normal
}

DEAL_II_NAMESPACE_CLOSE

#endif

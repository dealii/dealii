//----------------------------  solver_bicgstab.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  solver_bicgstab.h  ---------------------------
#ifndef __deal2__solver_bicgstab_h
#define __deal2__solver_bicgstab_h


#include <base/logstream.h>
#include <lac/solver.h>
#include <lac/solver_control.h>
#include <cmath>

/**
 * Bicgstab algorithm by van der Vorst.
 *
 * Like all other solver classes, this class has a local structure called
 * #AdditionalData# which is used to pass additional parameters to the
 * solver, like damping parameters or the number of temporary vectors. We
 * use this additional structure instead of passing these values directly
 * to the constructor because this makes the use of the #SolverSelector# and
 * other classes much easier and guarantees that these will continue to
 * work even if number or type of the additional parameters for a certain
 * solver changes.
 *
 * However, since the BiCGStab method does not need additional data, the respective
 * structure is empty and does not offer any functionality. The constructor
 * has a default argument, so you may call it without the additional
 * parameter.
 */
template <class VECTOR = Vector<double> >
class SolverBicgstab : private Solver<VECTOR>
{
  public:
    				     /**
				      * Standardized data struct to
				      * pipe additional data to the
				      * solver. This solver does not
				      * need additional data yet.
				      */
    struct AdditionalData {};

				     /**
				      * Constructor.
				      */
    SolverBicgstab (SolverControl &cn,
		    VectorMemory<VECTOR> &mem,
		    const AdditionalData &data=AdditionalData());

				     /**
				      * Virtual destructor.
				      */
    virtual ~SolverBicgstab ();

				     /**
				      * Solve primal problem only.
				      */
    template<class MATRIX, class Preconditioner>
    typename Solver<VECTOR>::ReturnState
    solve (const MATRIX &A,
	   VECTOR       &x,
	   const VECTOR &b,
	   const Preconditioner& precondition);

  protected:
				     /**
				      * Computation of the stopping criterion.
				      */
    template <class MATRIX>
    double criterion (const MATRIX& A, const VECTOR& x, const VECTOR& b);

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
			       const VECTOR& x,
			       const VECTOR& r,
			       const VECTOR& d) const;

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
    VECTOR *Vs;
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
  
  private:
				     /**
				      * Everything before the iteration loop.
				      */
    template <class MATRIX>
    SolverControl::State start(const MATRIX& A);

				     /**
				      * The iteration loop itself.
				      */
    template<class MATRIX, class PRECONDITIONER>
    typename Solver<VECTOR>::ReturnState 
    iterate(const MATRIX& A, const PRECONDITIONER& precondition);
  
};

/*-------------------------Inline functions -------------------------------*/

template<class VECTOR>
SolverBicgstab<VECTOR>::SolverBicgstab (SolverControl &cn,
					VectorMemory<VECTOR> &mem,
					const AdditionalData &)
		:
		Solver<VECTOR>(cn,mem)
{}



template<class VECTOR>
SolverBicgstab<VECTOR>::~SolverBicgstab ()
{}



template <class VECTOR>
template <class MATRIX>
double
SolverBicgstab<VECTOR>::criterion (const MATRIX& A, const VECTOR& x, const VECTOR& b)
{
  res = A.residual(*Vt, x, b);
  return res;
}



template <class VECTOR >
template <class MATRIX>
SolverControl::State
SolverBicgstab<VECTOR>::start(const MATRIX& A)
{
  res = A.residual(*Vr, *Vx, *Vb);
  Vp->reinit(*Vx);
  Vv->reinit(*Vx);
  Vrbar->equ(1.,*Vr);
  SolverControl::State state = control().check(step, res);
  return state;
}



template<class VECTOR>
void
SolverBicgstab<VECTOR>::print_vectors(const unsigned int,
				      const VECTOR&,
				      const VECTOR&,
				      const VECTOR&) const
{}



template<class VECTOR>
template<class MATRIX, class PRECONDITIONER>
typename Solver<VECTOR>::ReturnState
SolverBicgstab<VECTOR>::iterate(const MATRIX& A,
				const PRECONDITIONER& precondition)
{
  SolverControl::State state = SolverControl::iterate;
  alpha = omega = rho = 1.;

  VECTOR& r = *Vr;
  VECTOR& rbar = *Vrbar;
  VECTOR& p = *Vp;
  VECTOR& y = *Vy;
  VECTOR& z = *Vz;
  VECTOR& s = *Vs;
  VECTOR& t = *Vt;
  VECTOR& v = *Vv;
  
  do
    {
      rhobar = r*rbar;
      beta   = rhobar * alpha / (rho * omega);
      rho    = rhobar;
      p.sadd(beta, 1., r, -beta*omega, v);
      precondition(y,p);
      A.vmult(v,y);
      rhobar = rbar * v;

      alpha = rho/rhobar;

//TODO: Find better breakdown criterion (G)

      if (fabs(alpha) > 1.e10)
	return typename Solver<VECTOR>::ReturnState(breakdown);
    
      s.equ(1., r, -alpha, v);
      precondition(z,s);
      A.vmult(t,z);
      rhobar = t*s;
      omega = rhobar/(t*t);
      Vx->add(alpha, y, omega, z);
      r.equ(1., s, -omega, t);

      state = control().check(++step, criterion(A, *Vx, *Vb));
      print_vectors(step, *Vx, r, y);
    }
  while (state == SolverControl::iterate);
  if (state == SolverControl::success) return success;
  return exceeded;
}


template<class VECTOR>
template<class MATRIX, class PRECONDITIONER>
typename Solver<VECTOR>::ReturnState
SolverBicgstab<VECTOR>::solve(const MATRIX &A,
			      VECTOR       &x,
			      const VECTOR &b,
			      const PRECONDITIONER& precondition)
{
  deallog.push("Bicgstab");
  Vr    = memory.alloc(); Vr->reinit(x);
  Vrbar = memory.alloc(); Vrbar->reinit(x);
  Vp    = memory.alloc();
  Vy    = memory.alloc(); Vy->reinit(x);
  Vz    = memory.alloc(); Vz->reinit(x);
  Vs    = memory.alloc(); Vs->reinit(x);
  Vt    = memory.alloc(); Vt->reinit(x);
  Vv    = memory.alloc();

  Vx = &x;
  Vb = &b;

  step = 0;

  typename Solver<VECTOR>::ReturnState state = breakdown;
  
  do 
    {
      if (step)
	deallog << "Restart step " << step << endl;
      if (start(A) == SolverControl::success) break;  
      state = iterate(A, precondition);
    }
  while (state == breakdown);

  memory.free(Vr);
  memory.free(Vrbar);
  memory.free(Vp);
  memory.free(Vy);
  memory.free(Vz);
  memory.free(Vs);
  memory.free(Vt);
  memory.free(Vv);
  
  deallog.pop();
  return state;
}


#endif

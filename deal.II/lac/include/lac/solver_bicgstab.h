/*----------------------------   solver_bicgstab.h     ---------------------------*/
/*      $Id$                 */
#ifndef __solver_bicgstab_H
#define __solver_bicgstab_H
/*----------------------------   solver_bicgstab.h     ---------------------------*/

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
template <class Matrix = SparseMatrix<double>,
          class Vector = Vector<double> >
class SolverBicgstab : public Solver<Matrix,Vector>
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
    SolverBicgstab(SolverControl &cn,
		   VectorMemory<Vector> &mem,
		   const AdditionalData &data=AdditionalData());

				     /**
				      * Solve primal problem only.
				      */
    template<class Preconditioner>
    typename Solver<Matrix,Vector>::ReturnState
    solve (const Matrix &A,
	   Vector       &x,
	   const Vector &b,
	   const Preconditioner& precondition);

  protected:
				     /**
				      * Computation of the stopping criterion.
				      */
    virtual double criterion();

				     /**
				      * Auxiliary vector.
				      */
    Vector *Vx;
				     /**
				      * Auxiliary vector.
				      */
    Vector *Vr;
				     /**
				      * Auxiliary vector.
				      */
    Vector *Vrbar;
				     /**
				      * Auxiliary vector.
				      */
    Vector *Vp;
				     /**
				      * Auxiliary vector.
				      */
    Vector *Vy;
				     /**
				      * Auxiliary vector.
				      */
    Vector *Vz;
				     /**
				      * Auxiliary vector.
				      */
    Vector *Vs;
				     /**
				      * Auxiliary vector.
				      */
    Vector *Vt;
				     /**
				      * Auxiliary vector.
				      */
    Vector *Vv;
				     /**
				      * Right hand side vector.
				      */
    const Vector *Vb;
  
				     /**
				      * Pointer to the system matrix.
				      */
    const Matrix *MA;
  
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
    SolverControl::State start();

				     /**
				      * The iteration loop itself.
				      */
    template<class Preconditioner>
    typename Solver<Matrix,Vector>::ReturnState 
    iterate(const Preconditioner& precondition);
  
};

/*-------------------------Inline functions -------------------------------*/

template<class Matrix, class Vector>
SolverBicgstab<Matrix, Vector>::SolverBicgstab(SolverControl &cn,
					       VectorMemory<Vector> &mem,
					       const AdditionalData &) :
		Solver<Matrix,Vector>(cn,mem)  {}


template<class Matrix, class Vector>
double
SolverBicgstab<Matrix, Vector>::criterion()
{
  res = MA->residual(*Vt, *Vx, *Vb);
  return res;
}


template < class Matrix, class Vector >
SolverControl::State
SolverBicgstab<Matrix, Vector>::start()
{
  res = MA->residual(*Vr, *Vx, *Vb);
  Vp->reinit(*Vx);
  Vv->reinit(*Vx);
  Vrbar->equ(1.,*Vr);
  SolverControl::State state = control().check(step, res);
  return state;
}



template<class Matrix, class Vector>
template<class Preconditioner>
typename Solver<Matrix,Vector>::ReturnState
SolverBicgstab<Matrix, Vector>::iterate(const Preconditioner& precondition)
{
  SolverControl::State state = SolverControl::iterate;
  alpha = omega = rho = 1.;

  Vector& r = *Vr;
  Vector& rbar = *Vrbar;
  Vector& p = *Vp;
  Vector& y = *Vy;
  Vector& z = *Vz;
  Vector& s = *Vs;
  Vector& t = *Vt;
  Vector& v = *Vv;
  
  do
    {
      rhobar = r*rbar;
      beta   = rhobar * alpha / (rho * omega);
      rho    = rhobar;
      p.sadd(beta, 1., r, -beta*omega, v);
      precondition(y,p);
      MA->vmult(v,y);
      rhobar = rbar * v;

      alpha = rho/rhobar;

//TODO: Find better breakdown criterion (G)

      if (fabs(alpha) > 1.e10) return ReturnState(breakdown);
    
      s.equ(1., r, -alpha, v);
      precondition(z,s);
      MA->vmult(t,z);
      rhobar = t*s;
      omega = rhobar/(t*t);
      Vx->add(alpha, y, omega, z);
      r.equ(1., s, -omega, t);

      state = control().check(++step, criterion());
    }
  while (state == SolverControl::iterate);
  if (state == SolverControl::success) return success;
  return exceeded;
}


template<class Matrix, class Vector>
template<class Preconditioner>
typename Solver<Matrix,Vector>::ReturnState
SolverBicgstab<Matrix, Vector>::solve(const Matrix &A,
				      Vector       &x,
				      const Vector &b,
				      const Preconditioner& precondition)
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

  MA = &A;
  Vx = &x;
  Vb = &b;

  step = 0;

  ReturnState state = breakdown;
  
  do 
    {
      if (step)
	deallog << "Restart step " << step << endl;
      if (start() == SolverControl::success) break;  
      state = iterate(precondition);
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


/*----------------------------   solver_bicgstab.h     ---------------------------*/
/* end of #ifndef __solver_bicgstab_H */
#endif
/*----------------------------   solver_bicgstab.h     ---------------------------*/

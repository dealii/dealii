/*----------------------------   solver_qmrs.h     ---------------------------*/
/*      $Id$                 */
#ifndef __lac__solver_qmrs_H
#define __lac__solver_qmrs_H
/*----------------------------   solver_qmrs.h     ---------------------------*/



#include <lac/solver.h>
#include <lac/solver_control.h>
#include <base/logstream.h>
#include <cmath>


/**
 * QMRS method.
 *
 * The QMRS method is supposed to solve symmetric indefinite linear
 * systems with symmetric, not necessarily definite preconditioners.
 * This version of QMRS is adapted from
 * Freund/Nachtigal: Software for simplified Lanczos and QMR
 * algorithms, Appl. Num. Math. 19 (1995), pp. 319-341
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
 * However, since the QMRS method does not need additional data, the respective
 * structure is empty and does not offer any functionality. The constructor
 * has a default argument, so you may call it without the additional
 * parameter.
 *
 * @author Guido Kanschat, 1999
 */
template<class Matrix, class Vector>
class SolverQMRS : public Solver<Matrix,Vector>
{
  public:
    				     /**
				      * Standardized data struct to
				      * pipe additional data to the
				      * solver. This solver does not
				      * need additional data.
				      */
    struct AdditionalData {};

				     /**
				      * Constructor.
				      */
    SolverQMRS (SolverControl &cn,
	      VectorMemory<Vector> &mem,
	      const AdditionalData &data=AdditionalData());

				     /**
				      * Solver method.
				      */
    template<class Preconditioner>
    typename Solver<Matrix,Vector>::ReturnState
    solve (const Matrix &A,
		       Vector       &x,
		       const Vector &b,
		       const Preconditioner& precondition);

  protected:
				     /**
				      * Implementation of the computation of
				      * the norm of the residual.
				      */
    virtual long double criterion();
    
				     /**
				      * Temporary vectors, allocated through
				      * the #VectorMemory# object at the start
				      * of the actual solution process and
				      * deallocated at the end.
				      */
    Vector *Vv;
    Vector *Vp;
    Vector *Vq;
    Vector *Vt;
    Vector *Vd;
    
				     /**
				      * Within the iteration loop, the
				      * square of the residual vector is
				      * stored in this variable. The
				      * function #criterion# uses this
				      * variable to compute the convergence
				      * value, which in this class is the
				      * norm of the residual vector and thus
				      * the square root of the #res2# value.
				      */
    long double res2;
};




/*------------------------- Implementation ----------------------------*/
 

template<class Matrix, class Vector>
SolverQMRS<Matrix,Vector>::SolverQMRS(SolverControl &cn,
				  VectorMemory<Vector> &mem,
				  const AdditionalData &) :
		Solver<Matrix,Vector>(cn,mem) {};


template<class Matrix, class Vector>
long double
SolverQMRS<Matrix,Vector>::criterion()
{
  return sqrt(res2);
};



template<class Matrix, class Vector>
template<class Preconditioner>
typename Solver<Matrix,Vector>::ReturnState 
SolverQMRS<Matrix,Vector>::solve (const Matrix &A,
				Vector       &x,
				const Vector &b,
				const Preconditioner& precondition)
{
  SolverControl::State conv=SolverControl::iterate;

  deallog.push("QMRS");
  
				   // Memory allocation
  Vv  = memory.alloc();
  Vp  = memory.alloc();
  Vq  = memory.alloc();
  Vt  = memory.alloc();
  Vd  = memory.alloc();
				   // define some aliases for simpler access
  Vector& v  = *Vv; 
  Vector& p  = *Vp; 
  Vector& q  = *Vq; 
  Vector& t  = *Vt;
  Vector& d  = *Vd;
				   // resize the vectors, but do not set
				   // the values since they'd be overwritten
				   // soon anyway.
  v.reinit(x.size(), true);
  p.reinit(x.size(), true);
  q.reinit(x.size(), true);
  t.reinit(x.size(), true);
				   // This vector wants to be zero.
  d.reinit(x.size());

  int  it=0;
 
  double tau, rho, theta=0, sigma, alpha, psi, theta_old, rho_old, beta;
  
  double res = A.residual(v,x,b);
  conv = control().check(0,res);
  if (conv) 
    {
      memory.free(Vv);
      memory.free(Vp);
      memory.free(Vq);
      memory.free(Vt);
      memory.free(Vd);
      deallog.pop();
      return success;
    };

				   // Step 0
  p.equ(1.,v);
  precondition(q,p);
 
  tau = v.norm_sqr();
  rho = q*v;

  while (conv == SolverControl::iterate)
    {
      it++;
				       // Step 1
      A.vmult(t,p);
				       // Step 2
      sigma = q*t;
//      if (fabs(sigma) < ??)
//TODO: Breakdown criteria here and below

				       // Step 3
      alpha = rho/sigma;
      v.add(-alpha,t);
				       // Step 4
      theta_old = theta;
      theta = v*v/tau;
      psi = 1.*(1.+theta);
      tau *= theta*psi;

      d.sadd(psi*theta_old, psi*alpha, p);
      x.add(d);
				       // Step 5
      res = sqrt((it+1)*tau);
      conv = control().check(it,res);
      if (conv) break;
				       // Step 6
//      if (fabs(rho) < ??)
				       // Step 7
      rho_old = rho;
      precondition(q,v);
      rho = q*v;

      beta = rho/rho_old;
      p.sadd(beta,1.,v);
      precondition(q,p);
    };

    
  // Deallocate Memory
 
  memory.free(Vv);
  memory.free(Vp);
  memory.free(Vq);
  memory.free(Vt);
  memory.free(Vd);
 
  // Output

  deallog.pop();
 
  if (conv == SolverControl::failure)
    return exceeded;
  else
    return success;
};




/*----------------------------   solver_qmrs.h     ---------------------------*/
/* end of #ifndef __solver_qmrs_H */
#endif
/*----------------------------   solver_qmrs.h     ---------------------------*/

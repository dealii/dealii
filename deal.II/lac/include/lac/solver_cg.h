/*----------------------------   solver_pcg.h     ---------------------------*/
/*      $Id$                 */
#ifndef __solver_pcg_H
#define __solver_pcg_H
/*----------------------------   solver_pcg.h     ---------------------------*/



#include <lac/solver.h>
#include <lac/solver_control.h>




/**
 * Preconditioned cg method.
 * @author Original implementation by G. Kanschat, R. Becker and F.-T. Suttmeier, reworking and  documentation by Wolfgang Bangerth
 */
template<class Matrix, class Vector>
class SolverPCG : public Solver<Matrix,Vector> {
  public:

				     /**
				      * Constructor.
				      */
    SolverPCG (SolverControl &cn, VectorMemory<Vector> &mem) :
		    Solver<Matrix,Vector>(cn,mem) {};

				     /**
				      * Solver method.
				      */
    virtual ReturnState solve (const Matrix &A,
			       Vector       &x,
			       const Vector &b);

  protected:
				     /**
				      * Implementation of the computation of
				      * the norm of the residual.
				      */
    virtual double criterion();
    
				     /**
				      * Temporary vectors, allocated through
				      * the #VectorMemory# object at the start
				      * of the actual solution process and
				      * deallocated at the end.
				      */
    Vector *Vr, *Vp, *Vz, *VAp;
    
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
    double res2;
};




/*------------------------- Implementation ----------------------------*/
 

template<class Matrix, class Vector>
double SolverPCG<Matrix,Vector>::criterion()
{
  return sqrt(res2);
};



template<class Matrix, class Vector>
Solver<Matrix,Vector>::ReturnState 
SolverPCG<Matrix,Vector>::solve (const Matrix &A,
				 Vector       &x,
				 const Vector &b) {
  SolverControl::State conv=SolverControl::iterate;
  
				   // Memory allocation
  Vr  = memory.alloc();
  Vp  = memory.alloc();
  Vz  = memory.alloc();
  VAp = memory.alloc();
				   // define some aliases for simpler access
  Vector& g  = *Vr; 
  Vector& h  = *Vp; 
  Vector& d  = *Vz; 
  Vector& Ad = *VAp;
				   // resize the vectors, but do not set
				   // the values since they'd be overwritten
				   // soon anyway.
  g.reinit(x.size(), true);
  h.reinit(x.size(), true);
  d.reinit(x.size(), true);
  Ad.reinit(x.size(), true);
				   // Implementation taken from the DEAL
				   // library
  int  it=0;
  double res,gh,alpha,beta;
 
  res = A.residual(g,x,b);
  conv = control().check(0,res);
  if (conv) 
    {
      memory.free(Vr);
      memory.free(Vp);
      memory.free(Vz);
      memory.free(VAp);
      
      return success;
    };
  
  g.scale(-1.);
  A.precondition(h,g);
 
  d.equ(-1.,h);
 
  gh = g*h;
 
  while (conv == SolverControl::iterate)
    {
      A.vmult(Ad,d);
      
      alpha = d*Ad;
      alpha = gh/alpha;
      
      g.add(alpha,Ad);
      x.add(alpha,d );
      res = sqrt(g*g);

      conv = control().check(it,res);
      if (conv) break;
      
      A.precondition(h,g);
      
      beta = gh;
      gh   = g*h;
      beta = gh/beta;
      
      d.sadd(beta,-1.,h);
      it++;
    };

    
  // Deallocate Memory
 
  memory.free(Vr);
  memory.free(Vp);
  memory.free(Vz);
  memory.free(VAp);
 
  // Output
 
  if (conv == SolverControl::failure)
    return exceeded;
  else
    return success;
};




/*----------------------------   solver_pcg.h     ---------------------------*/
/* end of #ifndef __solver_pcg_H */
#endif
/*----------------------------   solver_pcg.h     ---------------------------*/

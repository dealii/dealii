//----------------------------  solver_qmrs.h  ---------------------------
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
//----------------------------  solver_qmrs.h  ---------------------------
#ifndef __deal2__solver_qmrs_h
#define __deal2__solver_qmrs_h


//TODO: Check for programming errors!!!

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
 * This version is for right preconditioning only, since then only the
 * preconditioner is used: left preconditioning seems to require the
 * inverse.
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
template <class Matrix = SparseMatrix<double>,
          class Vector = Vector<double> >
class SolverQMRS : public Solver<Matrix,Vector>
{
  public:
    				     /**
				      * Standardized data struct to
				      * pipe additional data to the
				      * solver. This solver does not
				      * need additional data.
				      *
				      * There are two possibilities to compute
				      * the residual: one is an estimate using
				      * the computed value #tau#. The other
				      * is exact computation using another matrix
				      * vector multiplication.
				      *
				      * QMRS, is susceptible to breakdowns, so
				      * we need a parameter telling us, which
				      * numbers are considered zero.
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
				      * Iteration vector.
				      */
    Vector *Vx;
				     /**
				      * RHS vector.
				      */
    const Vector *Vb;
    
				     /**
				      * Pointer to the matrix to be inverted.
				      */
    const Matrix* MA;
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
				     /**
				      * Breakdown threshold.
				      */
    AdditionalData additional_data;
  private:
				     /**
				      * The iteration loop itself.
				      */
    template<class Preconditioner>
    typename Solver<Matrix,Vector>::ReturnState 
    iterate(const Preconditioner& precondition);
				     /**
				      * The current iteration step.
				      */
    unsigned int step;
};


/*------------------------- Implementation ----------------------------*/


template<class Matrix, class Vector>
SolverQMRS<Matrix,Vector>::SolverQMRS(SolverControl &cn,
				  VectorMemory<Vector> &mem,
				  const AdditionalData &data) :
		Solver<Matrix,Vector>(cn,mem),
		additional_data(data)
{};


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
  deallog.push("QMRS");
  
				   // Memory allocation
  Vv  = memory.alloc();
  Vp  = memory.alloc();
  Vq  = memory.alloc();
  Vt  = memory.alloc();
  Vd  = memory.alloc();

  MA = &A;
  Vx = &x;
  Vb = &b;
				   // resize the vectors, but do not set
				   // the values since they'd be overwritten
				   // soon anyway.
  Vv->reinit(x.size(), true);
  Vp->reinit(x.size(), true);
  Vq->reinit(x.size(), true);
  Vt->reinit(x.size(), true);

  step = 0;
  
  ReturnState state = breakdown;

  do 
    {
      if (step)
	deallog << "Restart step " << step << endl;
      state = iterate(precondition);
    }
  while (state == breakdown);
    
  // Deallocate Memory
 
  memory.free(Vv);
  memory.free(Vp);
  memory.free(Vq);
  memory.free(Vt);
  memory.free(Vd);
 
  // Output

  deallog.pop();
 
  return state;
};

template<class Matrix, class Vector>
template<class Preconditioner>
typename Solver<Matrix,Vector>::ReturnState 
SolverQMRS<Matrix,Vector>::iterate(const Preconditioner& precondition)
{
/* Remark: the matrix A in the article is the preconditioned matrix.
 * Therefore, we have to precondition x before we compute the first residual.
 * In step 1 we replace p by q to avoid one preconditioning step.
 * There are still two steps left, making this algorithm expensive.
 */

  SolverControl::State state = SolverControl::iterate;

				   // define some aliases for simpler access
  Vector& v  = *Vv; 
  Vector& p  = *Vp; 
  Vector& q  = *Vq; 
  Vector& t  = *Vt;
  Vector& d  = *Vd;
  Vector& x  = *Vx;
  const Vector& b = *Vb;
  
  const Matrix& A  = *MA;
  
  int  it=0;
 
  double tau, rho, theta=0, sigma, alpha, psi, theta_old, rho_old, beta;
  double res;

  d.reinit(x);
  
				   // Apply right preconditioning to x
  precondition(q,x);  
				   // Preconditioned residual
  res = A.residual(v, q, b);

  if (control().check(step, res) == SolverControl::success)
    return ReturnState(success);

  p = v;
  
  precondition(q,p);
 
  tau = v.norm_sqr();
  //deallog << "tau:" << tau << endl;
  rho = q*v;
  //deallog << "rho:" << rho << endl;


while (state == SolverControl::iterate)
    {
      step++; it++;
				       // Step 1
      A.vmult(t,q);
				       // Step 2
      sigma = q*t;

//TODO: Find a really good breakdown criterion
// The absolute one detects breakdown instead of convergence
      if (fabs(sigma) < additional_data.breakdown)
	return ReturnState(breakdown);
				       // Step 3
      alpha = rho/sigma;
      //deallog << "alpha:" << alpha << endl;

      v.add(-alpha,t);
				       // Step 4
      theta_old = theta;
      theta = v*v/tau;
      psi = 1./(1.+theta);
      tau *= theta*psi;

      //deallog << "psi:" << psi << endl;
      //deallog << "theta:" << theta << endl;
      //deallog << "tau:" << tau << endl;

      d.sadd(psi*theta_old, psi*alpha, p);
      x.add(d);
				       // Step 5
      if (additional_data.exact_residual)
	res = A.residual(q,x,b);
      else
	res = sqrt((it+1)*tau);
      state = control().check(step,res);
      if (state == SolverControl::success)
	return ReturnState(success);
      else if (state == SolverControl::failure)
	return ReturnState(exceeded);
				       // Step 6
      if (fabs(rho) < additional_data.breakdown)
	return ReturnState(breakdown);
				       // Step 7
      rho_old = rho;
      precondition(q,v);
      rho = q*v;

      beta = rho/rho_old;
      p.sadd(beta,v);
      precondition(q,p);
    }
  return ReturnState(exceeded);
}


#endif

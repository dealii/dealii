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

#include <lac/solver.h>
#include <lac/solver_control.h>
#include <base/logstream.h>
#include <cmath>
#include <base/subscriptor.h>

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
 * @p{AdditionalData} which is used to pass additional parameters to the
 * solver, like damping parameters or the number of temporary vectors. We
 * use this additional structure instead of passing these values directly
 * to the constructor because this makes the use of the @p{SolverSelector} and
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
class SolverQMRS : public Subscriptor, private Solver<VECTOR>
{
  public:
    				     /**
				      * Standardized data struct to
				      * pipe additional data to the
				      * solver.
				      *
				      * There are two possibilities to compute
				      * the residual: one is an estimate using
				      * the computed value @p{tau}. The other
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
	  {};
	
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
				      * Solver method.
				      */
    template<class MATRIX, class PRECONDITIONER>
    typename Solver<VECTOR>::ReturnState
    solve (const MATRIX &A,
		       VECTOR       &x,
		       const VECTOR &b,
		       const PRECONDITIONER& precondition);

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
   protected:
				     /**
				      * Implementation of the computation of
				      * the norm of the residual.
				      */
    virtual double criterion();
    
				     /**
				      * Temporary vectors, allocated through
				      * the @p{VectorMemory} object at the start
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
				      * function @p{criterion} uses this
				      * variable to compute the convergence
				      * value, which in this class is the
				      * norm of the residual vector and thus
				      * the square root of the @p{res2} value.
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
    typename Solver<VECTOR>::ReturnState 
    iterate(const MATRIX& A, const PRECONDITIONER& precondition);
				     /**
				      * The current iteration step.
				      */
    unsigned int step;
};


/*------------------------- Implementation ----------------------------*/


template<class VECTOR>
SolverQMRS<VECTOR>::SolverQMRS(SolverControl &cn,
			       VectorMemory<VECTOR> &mem,
			       const AdditionalData &data) :
		Solver<VECTOR>(cn,mem),
		additional_data(data)
{};


template<class VECTOR>
double
SolverQMRS<VECTOR>::criterion()
{
  return sqrt(res2);
}



template<class VECTOR>
void
SolverQMRS<VECTOR>::print_vectors(const unsigned int,
				  const VECTOR&,
				  const VECTOR&,
				  const VECTOR&) const
{}



template<class VECTOR>
template<class MATRIX, class PRECONDITIONER>
typename Solver<VECTOR>::ReturnState 
SolverQMRS<VECTOR>::solve (const MATRIX &A,
			   VECTOR       &x,
			   const VECTOR &b,
			   const PRECONDITIONER& precondition)
{
  deallog.push("QMRS");
  
				   // Memory allocation
  Vv  = memory.alloc();
  Vp  = memory.alloc();
  Vq  = memory.alloc();
  Vt  = memory.alloc();
  Vd  = memory.alloc();

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
  
  typename Solver<VECTOR>::ReturnState state = breakdown;

  do 
    {
      if (step)
	deallog << "Restart step " << step << endl;
      state = iterate(A, precondition);
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



template<class VECTOR>
template<class MATRIX, class PRECONDITIONER>
typename Solver<VECTOR>::ReturnState 
SolverQMRS<VECTOR>::iterate(const MATRIX& A,
			    const PRECONDITIONER& precondition)
{
/* Remark: the matrix A in the article is the preconditioned matrix.
 * Therefore, we have to precondition x before we compute the first residual.
 * In step 1 we replace p by q to avoid one preconditioning step.
 * There are still two steps left, making this algorithm expensive.
 */

  SolverControl::State state = SolverControl::iterate;

				   // define some aliases for simpler access
  VECTOR& v  = *Vv; 
  VECTOR& p  = *Vp; 
  VECTOR& q  = *Vq; 
  VECTOR& t  = *Vt;
  VECTOR& d  = *Vd;
  VECTOR& x  = *Vx;
  const VECTOR& b = *Vb;
  
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

  if (control().check(step, res) == SolverControl::success)
    return success;

  p = v;
  
  precondition.vmult(q,p);
 
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
      if (fabs(sigma/rho) < additional_data.breakdown)
	return breakdown;
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

      print_vectors(step,x,v,d);
				       // Step 5
      if (additional_data.exact_residual)
	{
	  A.vmult(q,x);
	  q.sadd(-1.,1.,b);
	  res = q.l2_norm();
	}
      else
	res = sqrt((it+1)*tau);
      state = control().check(step,res);
      if (state == SolverControl::success)
	return success;
      else if (state == SolverControl::failure)
	return exceeded;
				       // Step 6
      if (fabs(rho) < additional_data.breakdown)
	return breakdown;
				       // Step 7
      rho_old = rho;
      precondition.vmult(q,v);
      rho = q*v;

      beta = rho/rho_old;
      p.sadd(beta,v);
      precondition.vmult(q,p);
    }
  return exceeded;
}


#endif

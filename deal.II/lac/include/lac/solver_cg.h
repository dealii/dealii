//----------------------------  solver_cg.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  solver_cg.h  ---------------------------
#ifndef __deal2__solver_cg_h
#define __deal2__solver_cg_h


#include <lac/solver.h>
#include <lac/solver_control.h>
#include <base/exceptions.h>
#include <base/logstream.h>
#include <base/subscriptor.h>
#include <cmath>


/**
 * Preconditioned cg method.
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
 * However, since the CG method does not need additional data, the respective
 * structure is empty and does not offer any functionality. The constructor
 * has a default argument, so you may call it without the additional
 * parameter.
 *
 * This version of CG is taken from Braess: "Finite Elements" and is analogous
 * to the one in the SIAM Templates Book. It requires a symmetric preconditioner,
 * i.e. SOR is not sufficient.
 *
 * @author Original implementation by G. Kanschat, R. Becker and F.-T. Suttmeier, reworking and  documentation by Wolfgang Bangerth
 */
template <class VECTOR = Vector<double> >
class SolverCG : public Subscriptor, private Solver<VECTOR>
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
					* Write coefficients alpha and
					* beta to the log file for later use.
					*/
      bool log_coefficients;

				       /**
					* Constructor. Initialize data fields.
					* Confer the description of those.
					*/
      AdditionalData (bool log_coefficients = false)
	:
	log_coefficients (log_coefficients)
	{}
    };

				     /**
				      * Constructor.
				      */
    SolverCG (SolverControl &cn,
	      VectorMemory<VECTOR> &mem,
	      const AdditionalData &data=AdditionalData());

				     /**
				      * Virtual destructor.
				      */
    virtual ~SolverCG ();

				     /**
				      * Solver method.
				      */
    template<class MATRIX, class PRECONDITIONER>
    void
    solve (const MATRIX &A,
	   VECTOR       &x,
	   const VECTOR &b,
	   const PRECONDITIONER& precondition);

  protected:
				     /**
				      * Implementation of the computation of
				      * the norm of the residual. This can be
				      * replaced by a more problem oriented
				      * functional in a derived class.
				      */
    virtual double criterion();

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
				      * Temporary vectors, allocated through
				      * the @p{VectorMemory} object at the start
				      * of the actual solution process and
				      * deallocated at the end.
				      */
    VECTOR *Vr;
    VECTOR *Vp;
    VECTOR *Vz;
    VECTOR *VAp;
    
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
				      * Additional parameters.
				      */
    AdditionalData additional_data;
};


/*------------------------- Implementation ----------------------------*/


template<class VECTOR>
SolverCG<VECTOR>::SolverCG(SolverControl &cn,
			   VectorMemory<VECTOR> &mem,
			   const AdditionalData &data)
		:
  Solver<VECTOR>(cn,mem),
  additional_data(data)
{}



template<class VECTOR>
SolverCG<VECTOR>::~SolverCG ()
{}



template<class VECTOR>
double
SolverCG<VECTOR>::criterion()
{
  return std::sqrt(res2);
}



template<class VECTOR>
void
SolverCG<VECTOR>::print_vectors(const unsigned int,
				const VECTOR&,
				const VECTOR&,
				const VECTOR&) const
{}



template<class VECTOR>
template<class MATRIX, class PRECONDITIONER>
void
SolverCG<VECTOR>::solve (const MATRIX &A,
			 VECTOR       &x,
			 const VECTOR &b,
			 const PRECONDITIONER& precondition)
{
  SolverControl::State conv=SolverControl::iterate;

  deallog.push("cg");
  
				   // Memory allocation
  Vr  = memory.alloc();
  Vp  = memory.alloc();
  Vz  = memory.alloc();
  VAp = memory.alloc();
				   // define some aliases for simpler access
  VECTOR& g  = *Vr; 
  VECTOR& h  = *Vp; 
  VECTOR& d  = *Vz; 
  VECTOR& Ad = *VAp;
				   // resize the vectors, but do not set
				   // the values since they'd be overwritten
				   // soon anyway.
  g.reinit(x, true);
  h.reinit(x, true);
  d.reinit(x, true);
  Ad.reinit(x, true);
				   // Implementation taken from the DEAL
				   // library
  int  it=0;
  double res,gh,alpha,beta;

				   // compute residual. if vector is
				   // zero, then short-circuit the
				   // full computation
  if (!x.all_zero())
    {
      A.vmult(g,x);
      g.sadd(-1.,1.,b);
    }
  else
    g = b;
  res = g.l2_norm();
  
  conv = control().check(0,res);
  if (conv) 
    {
      memory.free(Vr);
      memory.free(Vp);
      memory.free(Vz);
      memory.free(VAp);
      deallog.pop();
      return;
    };
  
  g.scale(-1.);
  precondition.vmult(h,g);
 
  d.equ(-1.,h);
 
  gh = g*h;
 
  while (conv == SolverControl::iterate)
    {
      it++;
      A.vmult(Ad,d);
      
      alpha = d*Ad;
      alpha = gh/alpha;
      
      g.add(alpha,Ad);
      x.add(alpha,d );
      res = g.l2_norm();

      print_vectors(it, x, g, d);
      
      conv = control().check(it,res);
      if (conv)
	break;
      
      precondition.vmult(h,g);
      
      beta = gh;
      gh   = g*h;
      beta = gh/beta;
      
      if (additional_data.log_coefficients)
	deallog << "alpha-beta:" << alpha << '\t' << beta << std::endl;
      
      d.sadd(beta,-1.,h);
    };


				   // Deallocate Memory
  memory.free(Vr);
  memory.free(Vp);
  memory.free(Vz);
  memory.free(VAp);
 
				   // Output
  deallog.pop();
 
				   // in case of failure: throw
				   // exception
  if (control().last_check() != SolverControl::success)
    throw SolverControl::NoConvergence (control().last_step(),
					control().last_value());
				   // otherwise exit as normal
};


#endif

//----------------------------  solver_gmres.h  ---------------------------
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
//----------------------------  solver_gmres.h  ---------------------------
#ifndef __deal2__solver_gmres_h
#define __deal2__solver_gmres_h


/*----------------------------   solver_pgmres.h     ---------------------------*/

#include <base/config.h>
#include <base/subscriptor.h>
#include <base/logstream.h>
#include <lac/solver.h>
#include <lac/solver_control.h>
#include <lac/full_matrix.h>
#include <vector>


/**
 * Implementation of the Restarted Preconditioned Direct Generalized
 * Minimal Residual Method. The stopping criterion is the norm of the
 * residual.
 *
 * You have to give the maximum number of temporary vectors to the
 * constructor which are to be used to do the orthogonalization. If
 * the number of iterations needed to solve the problem to the given
 * criterion, an intermediate solution is computed and a restart is
 * performed. If you don't want to use the restarted method, you can
 * limit the number of iterations (stated in the @p{SolverControl}
 * object given to the constructor) to be below the number of
 * temporary vectors minus three. Note the subtraction, which is due
 * to the fact that three vectors are used for other purposes, so the
 * number of iterations before a restart occurs is less by three than
 * the total number of temporary vectors. If the size of the matrix is
 * smaller than the maximum number of temporary vectors, then fewer
 * vectors are allocated since GMRES can then be used as a direct
 * solver.
 *
 * Note that restarts don't compensate temporary vectors very well,
 * i.e.  giving too few temporary vectors will increase the necessary
 * iteration steps greatly; it is not uncommon that you will need more
 * iterations than there are degrees of freedom in your solution
 * vector, if the number of temporary vectors is lower than the size
 * of the vector, even though GMRES is an exact solver if a sufficient
 * number of temporary vectors is given. You should therefore always
 * give as many temporary vectors as you can, unless you are limited
 * by the available memory; only then should you start to trade
 * computational speed against memory. One of the few other
 * possibilities is to use a good preconditioner.
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
 * For the GMRes method, the @p{AdditionalData} structure contains the number
 * of temporary vectors as commented upon above. By default, the number
 * of these vectors is set to 30.
 *
 * For the requirements on matrices and vectors in order to work with
 * this class, see the documentation of the @ref{Solver} base class.
 *
 * @author Wolfgang Bangerth
 */
template <class VECTOR = Vector<double> >
class SolverGMRES : public Solver<VECTOR>
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
					  * set the number of temporary
					  * vectors to 30.
					  */
	AdditionalData(const unsigned int max_n_tmp_vectors = 30)
			:
			max_n_tmp_vectors(max_n_tmp_vectors)
	  {};
	
					 /**
					  * Maximum number of
					  * tmp vectors.
					  */
	unsigned int    max_n_tmp_vectors;
    };
    
				     /**
				      * Constructor.
				      */
    SolverGMRES (SolverControl        &cn,
		 VectorMemory<VECTOR> &mem,
		 const AdditionalData &data=AdditionalData());

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

    DeclException1 (ExcTooFewTmpVectors,
		    int,
		    << "The number of temporary vectors you gave ("
		    << arg1 << ") is too small. It should be at least 10 for "
		    << "any results, and much more for reasonable ones.");
    
  protected:
				     /**
				      * Includes the maximum number of
				      * tmp vectors.
				      */
    AdditionalData additional_data;
    
				     /**
				      * Implementation of the computation of
				      * the norm of the residual.
				      */
    virtual double criterion();
    
				     /**
				      * Transformation of an upper
				      * Hessenberg matrix into
				      * tridiagonal structure by givens
				      * rotation of the last column
				      */
    void givens_rotation (Vector<double> &h,  Vector<double> &b,
			  Vector<double> &ci, Vector<double> &si,
			  int col) const;
				   /**
				    * Projected system matrix
				    */
  FullMatrix<double> H;
				   /**
				    * Auxiliary matrix for inverting @p{H}
				    */
  FullMatrix<double> H1;

  private:
				   /**
				    * No copy constructor.
				    */
  SolverGMRES (const SolverGMRES<VECTOR>&);
};


/* --------------------- Inline and template functions ------------------- */


template <class VECTOR>
SolverGMRES<VECTOR>::SolverGMRES (SolverControl        &cn,
					 VectorMemory<VECTOR> &mem,
					 const AdditionalData &data) :
		Solver<VECTOR> (cn,mem),
		additional_data(data)
{
  Assert (additional_data.max_n_tmp_vectors >= 10, 
	  ExcTooFewTmpVectors (additional_data.max_n_tmp_vectors));
};


template <class VECTOR>
inline
void
SolverGMRES<VECTOR>::givens_rotation (Vector<double> &h,
				      Vector<double> &b,
				      Vector<double> &ci,
				      Vector<double> &si, 
				      int     col) const
{
  for (int i=0 ; i<col ; i++)
    {
      const double s = si(i);
      const double c = ci(i);
      const double dummy = h(i);
      h(i)   =  c*dummy + s*h(i+1);
      h(i+1) = -s*dummy + c*h(i+1);
    };
  
  const double r = 1./sqrt(h(col)*h(col) + h(col+1)*h(col+1));
  si(col) = h(col+1) *r;
  ci(col) = h(col)   *r;
  h(col)  =  ci(col)*h(col) + si(col)*h(col+1);
  b(col+1)= -si(col)*b(col);
  b(col) *=  ci(col);
}


template<class VECTOR>
template<class MATRIX, class PRECONDITIONER>
void
SolverGMRES<VECTOR>::solve (const MATRIX         &A,
			    VECTOR               &x,
			    const VECTOR         &b,
			    const PRECONDITIONER &precondition)
{
				   // this code was written a very
				   // long time ago by people not
				   // associated with deal.II. we
				   // don't make any guarantees to its
				   // optimality or that it even works
				   // as expected...
  
//TODO:[?] Check, why there are two different start residuals.
//TODO:[?] Allocate vectors only when needed.

  deallog.push("GMRES");
  const unsigned int n_tmp_vectors = additional_data.max_n_tmp_vectors;

				   // allocate an array of n_tmp_vectors
				   // temporary vectors from the VectorMemory
				   // object; resize them but do not set their
				   // values since they will be overwritten soon
				   // anyway.

				   // This is really bad since vectors
				   // should only be allocated if
				   // really needed. (GK)
  std::vector<VECTOR*> tmp_vectors (n_tmp_vectors, 0);
  for (unsigned int tmp=0; tmp<n_tmp_vectors; ++tmp)
    {
      tmp_vectors[tmp] = memory.alloc();
      tmp_vectors[tmp]->reinit (x, true);
    };

				   // number of the present iteration; this
				   // number is not reset to zero upon a
				   // restart
  unsigned int accumulated_iterations = 0;
  
				   // matrix used for the orthogonalization
				   // process later
  H.reinit(n_tmp_vectors, n_tmp_vectors-1);

				   // some additional vectors, also used
				   // in the orthogonalization
  ::Vector<double> gamma(n_tmp_vectors),
                   ci   (n_tmp_vectors-1),
                   si   (n_tmp_vectors-1),
                   h    (n_tmp_vectors-1);


unsigned int dim = 0;

  SolverControl::State iteration_state = SolverControl::iterate;
  
				   // switch to determine whether we want a
				   // left or a right preconditioner. at
				   // present, left is default, but both
				   // ways are implemented
  const bool left_precondition = true;

				   // define two aliases
  VECTOR &v = *tmp_vectors[0];
  VECTOR &p = *tmp_vectors[n_tmp_vectors-1];


///////////////////////////////////
				   // outer iteration: loop until we
				   // either reach convergence or the
				   // maximum number of iterations is
				   // exceeded. each cycle of this
				   // loop amounts to one restart
  do
    {
				       // reset this vector to the
				       // right size
      h.reinit (n_tmp_vectors-1);
      
      if (left_precondition)
	{
	  A.vmult(p,x);
	  p.sadd(-1.,1.,b);
	  precondition.vmult(v,p);
	} else {
	  A.vmult(v,x);
	  v.sadd(-1.,1.,b);
	};

      double rho = v.l2_norm();

				       // check the residual here as
				       // well since it may be that
				       // we got the exact (or an almost
				       // exact) solution vector at
				       // the outset. if we wouldn't
				       // check here, the next scaling
				       // operation would produce
				       // garbage
      iteration_state = control().check (accumulated_iterations, rho);
      if (iteration_state != SolverControl::iterate)
	break;
      
      gamma(0) = rho;
      
      v.scale (1./rho);

				       // inner iteration doing at most as
				       // many steps as there are temporary
				       // vectors. the number of steps actually
				       // been done is propagated outside
				       // through the @p{dim} variable
      for (unsigned int inner_iteration=0;
	   ((inner_iteration < n_tmp_vectors-2)
	    &&
	    (iteration_state==SolverControl::iterate));
	   ++inner_iteration)
	{
	  ++accumulated_iterations;
	  // yet another alias
	  VECTOR& vv = *tmp_vectors[inner_iteration+1];
	  
	  if (left_precondition)
	    {
	      A.vmult(p, *tmp_vectors[inner_iteration]);
	      precondition.vmult(vv,p);
	    } else {
	      precondition.vmult(p,*tmp_vectors[inner_iteration]);
	      A.vmult(vv,p);
	    };
	  
	  dim = inner_iteration+1;
	  
					   /* Orthogonalization */
	  for (unsigned int i=0 ; i<dim ; ++i)
	    {
	      h(i) = vv * *tmp_vectors[i];
	      vv.add(-h(i),*tmp_vectors[i]);
	    };
      
	  double s = vv.l2_norm();
	  h(inner_iteration+1) = s;
    
					   /* Re-orthogonalization */
	  for (unsigned i=0; i<dim; ++i)
	    {
	      double htmp = vv * *tmp_vectors[i];
	      h(i) += htmp;
	      vv.add(-htmp,*tmp_vectors[i]);
	    };
	  
	  s = vv.l2_norm();
	  h(inner_iteration+1) = s;
	  
	  vv.scale(1./s);
	  
					   /*  Transformation into
					       triagonal structure  */
	  givens_rotation(h,gamma,ci,si,inner_iteration);
	  
					   /*  append vector on matrix  */
	  for (unsigned int i=0; i<dim; ++i)
	    H(i,inner_iteration) = h(i);
	  
					   /*  residual  */
	  rho = fabs(gamma(dim));
    
	  iteration_state = control().check (accumulated_iterations, rho);
	};

				       // end of inner iteration. now
				       // calculate the solution from the
				       // temporary vectors
      h.reinit(dim);
      H1.reinit(dim+1,dim);

      for (unsigned int i=0; i<dim+1; ++i)
	for (unsigned int j=0; j<dim; ++j)
	  H1(i,j) = H(i,j);

      H1.backward(h,gamma);

      if (left_precondition)
	for (unsigned int i=0 ; i<dim; ++i)
	  x.add(h(i), *tmp_vectors[i]);
      else
	{
	  p = 0.;
	  for (unsigned int i=0; i<dim; ++i)
	    p.add(h(i), *tmp_vectors[i]);
	  precondition.vmult(v,p);
	  x.add(1.,v);
	};
				       // end of outer iteration. restart if
				       // no convergence and the number of
				       // iterations is not exceeded
    }
  while (iteration_state == SolverControl::iterate);

				   // free the allocated memory before
				   // leaving
  for (unsigned int tmp=0; tmp<n_tmp_vectors; ++tmp)
    memory.free (tmp_vectors[tmp]);

  deallog.pop();

				   // in case of failure: throw
				   // exception
  if (control().last_check() != SolverControl::success)
    throw SolverControl::NoConvergence (control().last_step(),
					control().last_value());
				   // otherwise exit as normal
};


template<class VECTOR>
double
SolverGMRES<VECTOR>::criterion () 
{
				   // dummy implementation. this function is
				   // not needed for the present implementation
				   // of gmres
  Assert (false, ExcInternalError());
  return 0;
};


/*----------------------------   solver_pgmres.h     ---------------------------*/

#endif
/*----------------------------   solver_pgmres.h     ---------------------------*/

/*----------------------------   solver_pgmres.h     ---------------------------*/
/*      $Id$                 */
#ifndef __solver_pgmres_H
#define __solver_pgmres_H
/*----------------------------   solver_pgmres.h     ---------------------------*/

#include <lac/solver.h>
#include <lac/solver_control.h>
#include <lac/fullmatrix.h>
#include <vector>



/**
 * Implementation of the Restarted Preconditioned Direct Generalized
 * Minimal Residual Method. The stopping criterion is the norm of the
 * residual.
 *
 * @author Original implementation by the DEAL authors, adapted by Wolfgang Bangerth
 */
template<class Matrix, class Vector>
class SolverPGMRES : public Solver<Matrix, Vector> {
  public:
				     /**
				      * Constructor.
				      */
    SolverPGMRES (SolverControl        &cn,
		  VectorMemory<Vector> &mem,
		  const unsigned int    n_tmp_vectors) :
		    Solver<Matrix,Vector> (cn,mem),
		    n_tmp_vectors (n_tmp_vectors)
      {};
    
				     /**
				      * Solver method.
				      */
    virtual ReturnState solve (const Matrix &A,
			       Vector       &x,
			       const Vector &b);

  protected:
    const unsigned int n_tmp_vectors;
    
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
    void givens_rotation (Vector &h,  Vector &b,
			  Vector &ci, Vector &si,
			  int col) const;
};




/* ------------------------- Inline functions ----------------------------- */

template <class Matrix, class Vector>
inline
void
SolverPGMRES<Matrix,Vector>::givens_rotation (Vector& h, Vector& b,
					      Vector& ci, Vector& si, 
					      int col) const
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


// // restarted method
// template<class Matrix, class Vector>
// inline int
// SolverPGMRES<Matrix,Vector>::solve (Matrix& A, Vector& x, Vector& b)
// {
//   int reached, kmax = mem.n()-1;
//   for(int j=0;;j++)
//   {
//     info.usediter() = j*(kmax-1);
//     reached = dgmres(A,x,b,mem,info);
//     if(reached) break;
//   }
//   if (reached<0) return 1;
//   return 0;
// }




template<class Matrix, class Vector>
inline
Solver<Matrix,Vector>::ReturnState
SolverPGMRES<Matrix,Vector>::solve (const Matrix& A,
				    Vector& x,
				    const Vector& b)
{
				   // this code was written by the fathers of
				   // DEAL. I take absolutely no guarantees
				   // for any failures or airplane-explosions
				   // or nuclear wars or whatever resulting
				   // from this code. I tried to clean a bit,
				   // but whoever wrote this code in the first
				   // place should get stoned, IMHO! (WB)

  const unsigned int kmax = n_tmp_vectors-1;
				   // allocate an array of n_tmp_vectors
				   // temporary vectors from the VectorMemory
				   // object
  vector<Vector*> tmp_vectors (n_tmp_vectors, 0);
  for (unsigned int tmp=0; tmp<n_tmp_vectors; ++tmp)
    {
      tmp_vectors[tmp] = memory.alloc();
      tmp_vectors[tmp]->reinit (x.size());
    };

// WB  
//  int k0   = info.usediter();
  int k0 = 0;
  
				   // matrix used for the orthogonalization
				   // process later
  FullMatrix<double> H(kmax+1, kmax);

				   // some additional vectors, also used
				   // in the orthogonalization
  ::Vector<double> gamma(kmax+1), ci(kmax), si(kmax), h(kmax);


  unsigned int dim;

  SolverControl::State reached = SolverControl::iterate;
  
				   // switch to determine whether we want a
				   // left or a right preconditioner. at
				   // present, left is default, but both
				   // ways are implemented
  const bool left_precondition = true;

				   // define two aliases
  Vector &v = *tmp_vectors[0];
  Vector &p = *tmp_vectors[kmax];

  if (left_precondition)
    {
      A.residual(p,x,b);
      A.precondition(v,p);
    } else {
      A.residual(v,x,b);
    };
 
  double rho = v.l2_norm();
  gamma(0) = rho;

  v.scale (1./rho);

				   // inner iteration doing at most as
				   // many steps as there are temporary
				   // vectors. the number of steps actually
				   // been done is propagated outside
				   // through the #dim# variable
  for (unsigned int inner_iteration=0;
       inner_iteration<kmax-1 && (reached==SolverControl::iterate);
       inner_iteration++)
    {
				       // yet another alias
      Vector& vv = *tmp_vectors[inner_iteration+1];
      
      if (left_precondition)
	{
	  A.vmult(p, *tmp_vectors[inner_iteration]);
	  A.precondition(vv,p);
	} else {
	  A.precondition(p,*tmp_vectors[inner_iteration]);
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
    
      reached = control().check (k0+inner_iteration, rho);
    };
  
				   /*  Calculate solution  */
  h.reinit(dim);
  FullMatrix<double> H1(dim+1,dim);

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
      A.precondition(v,p);
      x.add(1.,v);
    };

				   // free the allocated memory before
				   // leaving
  for (unsigned int tmp=0; tmp<n_tmp_vectors; ++tmp)
    memory.free (tmp_vectors[tmp]);

  if (reached)
    return success;
  else
    return exceeded;
};


  


template<class Matrix, class Vector>
double
SolverPGMRES<Matrix,Vector>::criterion () 
{
				   // dummy implementation. this function is
				   // not needed for the present implementation
				   // of gmres
  Assert (false, ExcInternalError());
  return 0;
};



/*----------------------------   solver_pgmres.h     ---------------------------*/
/* end of #ifndef __solver_pgmres_H */
#endif
/*----------------------------   solver_pgmres.h     ---------------------------*/

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



/*
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
				   // but whoever wrote this code should get
				   // stone, IMHO! (WB)

  int kmax = n_tmp_vectors;
  FullMatrix<double> H(kmax+1, kmax), H1(kmax+1, kmax);
  
  ::Vector<double> y(kmax), b0(kmax+1);
  int i,k;

  SolverControl::State conv=SolverControl::iterate;
  
  double rho,beta;

				   // allocate an array of n_tmp_vectors
				   // temporary vectors from the VectorMemory
				   // object
  vector<Vector*> tmp_vectors (n_tmp_vectors, 0);
  for (unsigned int tmp=0; tmp<n_tmp_vectors; ++tmp)
    {
      tmp_vectors[tmp] = memory.alloc();
      tmp_vectors[tmp]->reinit (x.size());
    };

  
  A.residual(*tmp_vectors[0],x,b);
  
  rho = tmp_vectors[0]->l2_norm();
  beta = rho;
  
  tmp_vectors[0]->scale (1./rho);
  
  for (k=0 ; k<kmax-1 && (conv==SolverControl::iterate) ; k++)
  {
    A.vmult(*tmp_vectors[k+1], *tmp_vectors[k]);
    
    H.reinit(k+2,k+1);
    if (k>0)  H.fill(H1);
    
    for (i=0 ; i<=k ; i++)
    {
      H(i,k) = *tmp_vectors[k+1] * *tmp_vectors[i];
      tmp_vectors[k+1]->add(-H(i,k),*tmp_vectors[i]);
    }
    
    double s = tmp_vectors[k+1]->l2_norm();
    H(k+1,k) = s;
    
    // Re-orthogonalization
    
    //printf("\n");
    for (i=0 ; i<=k ; i++)
    {
      double htmp = *tmp_vectors[k+1] * *tmp_vectors[i];
      //printf(" %e ",htmp);
      H(i,k) += htmp;
      tmp_vectors[k+1]->add(-htmp,*tmp_vectors[i]);
    }
    //printf("\n");
    
    s = tmp_vectors[k+1]->l2_norm();
    H(k+1,k) = s;
    
    tmp_vectors[k+1]->scale(1./s);
    
    // Least - Squares
    
    y.reinit(k+1);
    b0.reinit(k+2);
    b0(0) = beta;
    H1 = H;
    rho = H.least_squares(y,b0);
    conv = control().check(k,rho);
  }

				   // this will miserably fail if the
				   // loop above was left before k=kmax-1!
  for (i=0 ; i<kmax ; i++)
//  for (i=0 ; i<k ; i++)
    x.add(y(i), *tmp_vectors[i]);


				   // free the allocated memory before
				   // leaving
  for (unsigned int tmp=0; tmp<n_tmp_vectors; ++tmp)
    memory.free (tmp_vectors[tmp]);


  if (conv == SolverControl::failure)
    return exceeded;
  else
    return success;
}
*/




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
				   // but whoever wrote this code should get
				   // stone, IMHO! (WB)

  int kmax = n_tmp_vectors-1;
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
  

  FullMatrix<double> H(kmax+1, kmax);
  ::Vector<double> gamma(kmax+1), ci(kmax), si(kmax), h(kmax);
  int i,k,reached=0,dim;
  int left_precondition = 1;

  Vector& v = *tmp_vectors[0];
  Vector& p = *tmp_vectors[kmax];

  if (left_precondition)
  {
    A.residual(p,x,b);
    A.precondition(v,p);
  }
  else
  {
    A.residual(v,x,b);
  }
 
  double rho = sqrt(v*v);
  gamma(0) = rho;
  
  v.equ(1./rho,v);

  for (k=0 ; k<kmax-1 && (!reached) ; k++)
  {
    Vector& vv = *tmp_vectors[k+1];
    
    if (left_precondition)
    {
      A.vmult(p, *tmp_vectors[k]);
      A.precondition(vv,p);
    }
    else
    {
      A.precondition(p,*tmp_vectors[k]);
      A.vmult(vv,p);
    }

// WB why is this here?    
//    double s0 = sqrt(vv*vv);
    dim = k+1;

    /* Orthogonalization */
    
    for (i=0 ; i<dim ; i++)
    {
      h(i) = vv * *tmp_vectors[i];
      vv.add(-h(i),*tmp_vectors[i]);
    }
    double s = sqrt(vv*vv);
    h(k+1) = s;
    
    /* Re-orthogonalization */
    
    for (i=0 ; i<dim ; i++)
    {
      double htmp = vv * *tmp_vectors[i];
      h(i) += htmp;
      vv.add(-htmp,*tmp_vectors[i]);
    }
    s = sqrt(vv*vv);
    h(k+1) = s;
    
    vv.equ(1./s, vv);
    
    /*  Transformation into triagonal structure  */
    
    givens_rotation(h,gamma,ci,si,k);

    /*  append vector on matrix  */

    for (i=0 ; i<dim ; i++)
      H(i,k) = h(i);

    /*  residual  */

    rho = fabs(gamma(dim));
    
// WB    
//    reached = info.check_residual(k0+k,rho);
    reached = control().check (k0+k, rho);
  }
  
  /*  Calculate solution  */  

  h.reinit(dim);
  FullMatrix<double> H1(dim+1,dim);

  for (i=0 ; i<dim+1 ; i++) {
    for (int j=0 ; j<dim ; j++) {
      H1(i,j) = H(i,j);
    }
  }

  H1.backward(h,gamma);

  if (left_precondition)
  {
    for (i=0 ; i<dim ; i++) {
      x.add(h(i), *tmp_vectors[i]);
    }
  }
  else
  {
    p = 0.;
    for (i=0 ; i<dim ; i++)
      p.add(h(i), *tmp_vectors[i]);
    A.precondition(v,p);
    x.add(1.,v);
  }

				   // free the allocated memory before
				   // leaving
  for (unsigned int tmp=0; tmp<n_tmp_vectors; ++tmp)
    memory.free (tmp_vectors[tmp]);

// WB  
//  return reached;
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

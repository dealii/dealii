/*----------------------   precondition_selector.h     -------------------------*/
/*      $Id$                 */
/*                Ralf Hartmann, University of Heidelberg                       */
#ifndef __precondition_selector_H
#define __precondition_selector_H
/*----------------------   precondition_selector.h     -------------------------*/


#include <base/smartpointer.h>
#include <string>

/**
 * Selects the preconditioner. The constructor of this class takes 
 * the name of the preconditioning as string, the matrix that is used
 * by the matrix-builtin precondition functions and a damping parameter 
 * #omega# of the preconditioning. Each time, the #operator()# function
 * is called, this preselected preconditioner, this matrix and
 * this #omega# is used
 * for the preconditioning. This class is designed for being used as
 * argument of the #solve# function of a #Solver# and it covers the
 * selection of all matrix-builtin precondition functions. The selection
 * of other preconditioners, like BlockSOR or ILU should be handled in
 * derived classes by the user.
 *
 * \subsection{Usage}
 * The simplest use of this class is the following:
 * \begin{verbatim}
 *                                  // generate a #SolverControl# and
 *                                  // a #VectorMemory#
 * SolverControl control;
 * VectorMemory<Vector<double> > memory;
 *                                  // generate a solver
 * SolverCG<SparseMatrix<double>, Vector<double> > solver(control, memory);
 *                                  // generate a #PreconditionSelector#
 * PreconditionSelector<SparseMatrix<double>, Vector<double> >
 *   preconditioning("jacobi", 1.);
 *                                  // give a matrix whose diagonal entries
 *                                  // are to be used for the preconditioning.
 *                                  // Generally the matrix of the linear
 *                                  // equation system Ax=b.
 * preconditioning.use_matrix(A);
 *                                  // call the #solve# function with this
 *                                  // preconditioning as last argument
 * solver.solve(A,x,b,preconditioning);
 * \end{verbatim}
 * The same example where also the #SolverSelector# class is used reads
 * \begin{verbatim}
 *                                  // generate a #SolverControl# and
 *                                  // a #VectorMemory#
 * SolverControl control;
 * VectorMemory<Vector<double> > memory;
 *                                  // generate a #SolverSelector# that
 *                                  // calls the #SolverCG#
 * SolverSelector<SparseMatrix<double>, Vector<double> > 
 *   solver_selector("cg", control, memory);
 *                                  // generate a #PreconditionSelector#
 * PreconditionSelector<SparseMatrix<double>, Vector<double> >
 *   preconditioning("jacobi", 1.);
 *
 * preconditioning.use_matrix(A);
 *
 * solver_selector.solve(A,x,b,preconditioning);
 * \end{verbatim}
 *
 * @author Ralf Hartmann, 1999
 */
template <class Matrix, class Vector>
class PreconditionSelector
{
  public:
    
				     /**
				      * Constructor. #omega# denotes
				      * the damping parameter of
				      * the preconditioning.
				      */
    PreconditionSelector(const string preconditioning,
			 const typename Vector::value_type &omega=1.);
    
				     /**
				      * Destructor.
				      */
    virtual ~PreconditionSelector();

				     /**
				      * Takes the matrix that is needed
				      * for preconditionings that involves a
				      * matrix. e.g. for precondition_jacobi,
				      * ~_sor, ~_ssor.
				      */
    void use_matrix(const Matrix &M);

				     /**
				      * Precondition procedure. Calls the
				      * preconditioning that was specified in
				      * the constructor.
				      */    
    virtual void operator() (Vector &dst, const Vector &src) const;

				     /**
				      * Get the names of all implemented
				      * preconditionings.
				      */
    static string get_precondition_names();

				     /**
				      * Exception.
				      */
    DeclException0 (ExcNoMatrixGivenToUse);


  protected:
    
				     /**
				      * Stores the name of the
				      * preconditioning.
				      */
    string preconditioning;

  private:
				     /**
				      * Matrix that is used for the
				      * matrix-builtin preconditioning function.
				      * cf. also #PreconditionUseMatrix#.
				      */
    SmartPointer<const Matrix> A;
    
				     /**
				      * Stores the damping parameter
				      * of the preconditioner.
				      */
    const typename Vector::value_type omega;
};


/* --------------------- Inline and template functions ------------------- */


template <class Matrix, class Vector>
PreconditionSelector<Matrix, Vector>
::PreconditionSelector(string preconditioning,
		       const typename Vector::value_type &omega) :
		preconditioning(preconditioning),
		omega(omega)  {}


template <class Matrix, class Vector>
PreconditionSelector<Matrix, Vector>::~PreconditionSelector()
{
				   // release the matrix A
  A=0;
}


template <class Matrix, class Vector>
void PreconditionSelector<Matrix, Vector>::use_matrix(const Matrix &M)
{
  A=&M;
}

template <class Matrix, class Vector>
void PreconditionSelector<Matrix, Vector>::operator() (Vector &dst,
						       const Vector &src) const
{
  if (preconditioning=="none")
    {
      dst=src;
    }
  else 
    {
      Assert(A!=0, ExcNoMatrixGivenToUse());
      
      if (preconditioning=="jacobi")
	{
	  A->precondition_Jacobi(dst,src,omega);
	} 
      else if (preconditioning=="sor")
	{
	  A->precondition_SOR(dst,src,omega);      
	}
      else if (preconditioning=="ssor")
	{
	  A->precondition_SSOR(dst,src,omega);      
	}
      else 
	Assert(false,ExcNotImplemented());
    }
}


template <class Matrix, class Vector>
string PreconditionSelector<Matrix, Vector>::get_precondition_names()
{
  return "none|jacobi|sor|ssor";
}



/*----------------------------   precondition_selector.h     ---------------------------*/
/* end of #ifndef __precondition_selector_H */
#endif
/*----------------------------   precondition_selector.h     ---------------------------*/

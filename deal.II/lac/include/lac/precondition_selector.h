//----------------------------  precondition_selector.h  ---------------------------
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
//----------------------------  precondition_selector.h  ---------------------------
#ifndef __deal2__precondition_selector_h
#define __deal2__precondition_selector_h


#include <base/smartpointer.h>
#include <string>

/**
 * Selects the preconditioner. The constructor of this class takes 
 * the name of the preconditioning and the damping parameter 
 * @p{omega} of the preconditioning and the @p{use_matrix} function takes
 * the matrix that is used
 * by the matrix-builtin precondition functions. Each time, the @p{operator()} function
 * is called, this preselected preconditioner, this matrix and
 * this @p{omega} is used
 * for the preconditioning. This class is designed for being used as
 * argument of the @p{solve} function of a @p{Solver} and it covers the
 * selection of all matrix-builtin precondition functions. The selection
 * of other preconditioners, like BlockSOR or ILU should be handled in
 * derived classes by the user.
 *
 * @sect3{Usage}
 * The simplest use of this class is the following:
 * @begin{verbatim}
 *                                  // generate a @p{SolverControl} and
 *                                  // a @p{VectorMemory}
 * SolverControl control;
 * VectorMemory<Vector<double> > memory;
 *                                  // generate a solver
 * SolverCG<SparseMatrix<double>, Vector<double> > solver(control, memory);
 *                                  // generate a @p{PreconditionSelector}
 * PreconditionSelector<SparseMatrix<double>, Vector<double> >
 *   preconditioning("jacobi", 1.);
 *                                  // give a matrix whose diagonal entries
 *                                  // are to be used for the preconditioning.
 *                                  // Generally the matrix of the linear
 *                                  // equation system Ax=b.
 * preconditioning.use_matrix(A);
 *                                  // call the @p{solve} function with this
 *                                  // preconditioning as last argument
 * solver.solve(A,x,b,preconditioning);
 * @end{verbatim}
 * The same example where also the @p{SolverSelector} class is used reads
 * @begin{verbatim}
 *                                  // generate a @p{SolverControl} and
 *                                  // a @p{VectorMemory}
 * SolverControl control;
 * VectorMemory<Vector<double> > memory;
 *                                  // generate a @p{SolverSelector} that
 *                                  // calls the @p{SolverCG}
 * SolverSelector<SparseMatrix<double>, Vector<double> > 
 *   solver_selector("cg", control, memory);
 *                                  // generate a @p{PreconditionSelector}
 * PreconditionSelector<SparseMatrix<double>, Vector<double> >
 *   preconditioning("jacobi", 1.);
 *
 * preconditioning.use_matrix(A);
 *
 * solver_selector.solve(A,x,b,preconditioning);
 * @end{verbatim}
 * Now the use of the @p{SolverSelector} in combination with the @p{PreconditionSelector}
 * allows the user to select both, the solver and the preconditioner, at the
 * beginning of his program and each time the
 * solver is started (that is several times e.g. in a nonlinear iteration) this
 * preselected solver and preconditioner is called.
 *
 * @author Ralf Hartmann, 1999
 */
template <class Matrix = SparseMatrix<double>,
          class Vector = Vector<double> >
class PreconditionSelector : public Subscriptor
{
  public:
    
				     /**
				      * Constructor. @p{omega} denotes
				      * the damping parameter of
				      * the preconditioning.
				      */
    PreconditionSelector(const std::string                 &preconditioning,
			 const typename Vector::value_type &omega=1.);
    
				     /**
				      * Destructor.
				      */
    virtual ~PreconditionSelector();

				     /**
				      * Takes the matrix that is needed
				      * for preconditionings that involves a
				      * matrix. e.g. for @p{precondition_jacobi},
				      * @p{~_sor}, @p{~_ssor}.
				      */
    void use_matrix(const Matrix &M);

				     /**
				      * Precondition procedure. Calls the
				      * preconditioning that was specified in
				      * the constructor.
				      */    
    virtual void vmult (Vector &dst, const Vector&src) const;

				     /**
				      * Get the names of all implemented
				      * preconditionings.
				      */
    static std::string get_precondition_names();

				     /**
				      * Exception.
				      */
    DeclException0 (ExcNoMatrixGivenToUse);


  protected:
    
				     /**
				      * Stores the name of the
				      * preconditioning.
				      */
    std::string preconditioning;

  private:
				     /**
				      * Matrix that is used for the
				      * matrix-builtin preconditioning function.
				      * cf. also @p{PreconditionUseMatrix}.
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
PreconditionSelector<Matrix,Vector>
::PreconditionSelector(const std::string                 &preconditioning,
		       const typename Vector::value_type &omega) :
		preconditioning(preconditioning),
		omega(omega)  {}


template <class Matrix, class Vector>
PreconditionSelector<Matrix,Vector>::~PreconditionSelector()
{
				   // release the matrix A
  A=0;
}


template <class Matrix, class Vector>
void PreconditionSelector<Matrix,Vector>::use_matrix(const Matrix &M)
{
  A=&M;
}

template <class Matrix, class Vector>
void PreconditionSelector<Matrix,Vector>::vmult (Vector &dst,
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
std::string PreconditionSelector<Matrix,Vector>::get_precondition_names()
{
  return "none|jacobi|sor|ssor";
}


#endif

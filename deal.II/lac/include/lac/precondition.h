//----------------------------  precondition.h  ---------------------------
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
//----------------------------  precondition.h  ---------------------------
#ifndef __deal2__precondition_h
#define __deal2__precondition_h

// This file contains simple preconditioners.

#include <lac/vector_memory.h>
#include <base/smartpointer.h>

template <typename number> class Vector;
template <typename number> class SparseMatrix;


/**
 * No preconditioning.  This class helps you, if you want to use a
 * linear solver without preconditioning. All solvers in LAC require a
 * preconditioner. Therefore, you must use the identity provided here
 * to avoid preconditioning.
 *
 * @author Guido Kanschat, 1999
 */
class PreconditionIdentity : public Subscriptor
{
  public:
				     /**
				      * Apply preconditioner.
				      */
    template<class VECTOR>
    void vmult (VECTOR&, const VECTOR&) const;

				     /**
				      * Apply transpose
				      * preconditioner. Since this is
				      * the identity, this function is
				      * the same as
				      * @ref{vmult}.
				      */
    template<class VECTOR>
    void Tvmult (VECTOR&, const VECTOR&) const;
};



/**
 * Preconditioner using a matrix-builtin function.
 * This class forms a preconditioner suitable for the LAC solver
 * classes. Since many preconditioning methods are based on matrix
 * entries, these have to be implemented as member functions of the
 * underlying matrix implementation. This class now is intended to
 * allow easy access to these member functions from LAC solver
 * classes.
 *
 * It seems that all builtin preconditioners have a relaxation
 * parameter, so please use @p{PreconditionRelaxation} for these.
 *
 * @sect3{Use}
 * You will usually not want to create a named object of this type,
 * although possible. The most common use is like this:
 * @begin{verbatim}
 *    SolverGMRES<SparseMatrix<double>,
 *                Vector<double> >      gmres(control,memory,500);
 *
 *    gmres.solve (matrix, solution, right_hand_side,
 *		   PreconditionUseMatrix<SparseMatrix<double>,Vector<double> >
 *		   (matrix,&SparseMatrix<double>::template precondition_Jacobi));
 * @end{verbatim}
 * This creates an unnamed object to be passed as the fourth parameter to
 * the solver function of the @p{SolverGMRES} class. It assumes that the
 * @p{SparseMatrix} class has a function @p{precondition_Jacobi} taking two
 * vectors (source and destination) as parameters. (Actually, there is no
 * function like that, the existing function takes a third parameter,
 * denoting the relaxation parameter; this example is therefore only meant to
 * illustrate the general idea.)
 *
 * Note that due to the default template parameters, the above example
 * could be written shorter as follows:
 * @begin{verbatim}
 *    ...
 *    gmres.solve (matrix, solution, right_hand_side,
 *		   PreconditionUseMatrix<>
 *		     (matrix,&SparseMatrix<double>::template precondition_Jacobi));
 * @end{verbatim}
 *
 * @author Guido Kanschat, Wolfgang Bangerth, 1999
 */
template<class MATRIX = SparseMatrix<double>, class VECTOR = Vector<double> >
class PreconditionUseMatrix : public Subscriptor
{
  public:
				     /**
				      * Type of the preconditioning
				      * function of the matrix.
				      */
    typedef void ( MATRIX::* function_ptr)(VECTOR&, const VECTOR&) const;
    
				     /**
				      * Constructor.
				      * This constructor stores a
				      * reference to the matrix object
				      * for later use and selects a
				      * preconditioning method, which
				      * must be a member function of
				      * that matrix.
				      */
    PreconditionUseMatrix(const MATRIX      &M,
			  const function_ptr method);
    
				     /**
				      * Execute preconditioning. Calls the
				      * function passed to the constructor
				      * of this object with the two
				      * arguments given here.
				      */
    void vmult (VECTOR       &dst,
		const VECTOR &src) const;

  private:
				     /**
				      * Pointer to the matrix in use.
				      */
    const MATRIX& matrix;
    
				     /**
				      * Pointer to the preconditioning
				      * function.
				      */
    const function_ptr precondition;
};



/**
 * Base class for other preconditioners.
 * Here, only some common features Jacobi, SOR and SSOR preconditioners
 * are implemented. For preconditioning, refer to derived classes.
 *
 * @author Guido Kanschat, 2000
 */
template<class MATRIX = SparseMatrix<double> >
class PreconditionRelaxation : public Subscriptor
{
  public:
				     /**
				      * Initialize matrix and
				      * relaxation parameter. The
				      * matrix is just stored in the
				      * preconditioner object. The
				      * relaxation parameter should be
				      * larger than zero and smaller
				      * than 2 for numerical
				      * reasons. It defaults to 1.
				      */
    void initialize (const MATRIX& A, const double omega = 1.);
    
  protected:
				     /**
				      * Pointer to the matrix object.
				      */
    SmartPointer<const MATRIX> A;

				     /**
				      * Relaxation parameter.
				      */
    double omega;
};



/**
 * Jacobi preconditioner using matrix built-in function.  The MATRIX
 * class used is required to have a function
 * @p{precondition_Jacobi(VECTOR&, const VECTOR&, double}
 *
 * @sect2{Usage example}
 * @begin{itemize}
 *     // Declare related objects
 *
 * SparseMatrix<double> A;
 * Vector<double> x;
 * Vector<double> b;
 * SolverCG<> solver(...);
 *
 * //...initialize and build A
 *
 *     // Define and initialize preconditioner
 *
 * PreconditionJacobi<SparseMatrix<double> > precondition;
 * precondition.initialize (A, .6);
 *
 * solver.solve (A, x, b, precondition);
 * @end{itemize}
 *
 * @author Guido Kanschat, 2000
 */
template <class MATRIX = SparseMatrix<double> >
class PreconditionJacobi : public PreconditionRelaxation<MATRIX>
{
  public:
				     /**
				      * Apply preconditioner.
				      */
    template<class VECTOR>
    void vmult (VECTOR&, const VECTOR&) const;
				     /**
				      * Apply transpose
				      * preconditioner. Since this is
				      * a symmetric preconditioner,
				      * this function is the same as
				      * @ref{vmult}.
				      */
    template<class VECTOR>
    void Tvmult (VECTOR&, const VECTOR&) const;
};


/**
 * SOR preconditioner using matrix built-in function.  The MATRIX
 * class used is required to have functions
 * @p{precondition_SOR(VECTOR&, const VECTOR&, double)} and
 * @p{precondition_TSOR(VECTOR&, const VECTOR&, double)}.
 *
 *
 * @sect2{Usage example}
 * @begin{itemize}
 *     // Declare related objects
 *
 * SparseMatrix<double> A;
 * Vector<double> x;
 * Vector<double> b;
 * SolverCG<> solver(...);
 *
 * //...initialize and build A
 *
 *     // Define and initialize preconditioner
 *
 * PreconditionSOR<SparseMatrix<double> > precondition;
 * precondition.initialize (A, .6);
 *
 * solver.solve (A, x, b, precondition);
 * @end{itemize}
 *
 * @author Guido Kanschat, 2000
 */
template <class MATRIX = SparseMatrix<double> >
class PreconditionSOR : public PreconditionRelaxation<MATRIX>
{
  public:
				     /**
				      * Apply preconditioner.
				      */
    template<class VECTOR>
    void vmult (VECTOR&, const VECTOR&) const;

				     /**
				      * Apply transpose
				      * preconditioner.
				      */
    template<class VECTOR>
    void Tvmult (VECTOR&, const VECTOR&) const;
};



/**
 * SSOR preconditioner using matrix built-in function.  The MATRIX
 * class used is required to have a function
 * @p{precondition_SSOR(VECTOR&, const VECTOR&, double}
 *
 *
 * @sect2{Usage example}
 * @begin{itemize}
 *     // Declare related objects
 *
 * SparseMatrix<double> A;
 * Vector<double> x;
 * Vector<double> b;
 * SolverCG<> solver(...);
 *
 * //...initialize and build A
 *
 *     // Define and initialize preconditioner
 *
 * PreconditionSSOR<SparseMatrix<double> > precondition;
 * precondition.initialize (A, .6);
 *
 * solver.solve (A, x, b, precondition);
 * @end{itemize}
 *
 * @author Guido Kanschat, 2000
 */
template <class MATRIX = SparseMatrix<double> >
class PreconditionSSOR : public PreconditionRelaxation<MATRIX>
{
  public:
				     /**
				      * Apply preconditioner.
				      */
    template<class VECTOR>
    void vmult (VECTOR&, const VECTOR&) const;

				     /**
				      * Apply transpose
				      * preconditioner. Since this is
				      * a symmetric preconditioner,
				      * this function is the same as
				      * @ref{vmult}.
				      */
    template<class VECTOR>
    void Tvmult (VECTOR&, const VECTOR&) const;
};


/**
 * Preconditioner using an iterative solver.  This preconditioner uses
 * a fully initialized LAC iterative solver for the approximate
 * inverse of the matrix. Naturally, this solver needs another
 * preconditionig method.
 *
 * Usually, the use of @p{ReductionControl} is preferred over the use of
 * the basic @p{SolverControl} in defining this solver.
 *
 * @sect2{Usage example}
 *
 * Krylov space methods like @ref{SolverCG} or @ref{SolverBicgstab}
 * become inefficient if soution down to machine accuracy is
 * needed. This is due to the fact, that round-off errors spoil the
 * orthogonality of the vector sequences. Therefore, a nested
 * iteration of two methods is proposed: The outer method is
 * @ref{SolverRichardson}, since it is robust with respect to round-of
 * errors. The inner loop is an appropriate Krylov space method, since
 * it is fast.
 *
 * @begin{itemize}
 *     // Declare related objects
 *
 * SparseMatrix<double> A;
 * Vector<double> x;
 * Vector<double> b;
 * GrowingVectorMemory<Vector<double> > mem;

 * ReductionControl inner_control (10, 1.e-30, 1.e-2)
 * SolverCG<Vector<double> > inner_iteration (inner_control, mem);
 * PreconditionSSOR <SparseMatrix<double> > inner_precondition;
 * inner_precondition.initialize (A, 1.2);
 *
 * PreconditionLACSolver precondition;
 * precondition.initialize (inner_iteration, A, inner_precondition);
 *
 * SolverControl outer_control(100, 1.e-16);
 * SolverRichardson<Vector<double> > outer_iteration;
 *
 * outer_iteration.solve (A, x, b, precondition);
 * @end{itemize}
 *
 * Each time we call the inner loop, reduction of the residual by a
 * factor @p{1.e-2} is attempted. Since the right hand side vector of
 * the inner iteration is the residual of the outer loop, the relative
 * errors are far from machine accuracy, even if the errors of the
 * outer loop are in the range of machine accuracy.
 *
 * @author Guido Kanschat, 1999
 */
template<class SOLVER, class MATRIX = SparseMatrix<double>, class PRECONDITION = PreconditionIdentity>
class PreconditionLACSolver : public Subscriptor
{
  public:
				     /**
				      * Constructor. All work is done
				      * in initialize.
				      */
    PreconditionLACSolver ();

				     /**
				      * Initialization
				      * function. Provide a solver
				      * object, a matrix, and another
				      * preconditioner for this.
				      */
    void initialize (SOLVER&,
		     const MATRIX&,
		     const PRECONDITION&);
    
				     /**
				      * Execute preconditioning.
				      */
    template<class VECTOR>
    void vmult (VECTOR&, const VECTOR&) const;

  private:
				     /**
				      * The solver object to use.
				      */
    SmartPointer<SOLVER> solver;

				     /**
				      * The matrix in use.
				      */
    SmartPointer<const MATRIX> matrix;
    
				     /**
				      * The preconditioner to use.
				      */
    SmartPointer<const PRECONDITION> precondition;
};



/**
 * Matrix with preconditioner.
 * Given a matrix $A$ and a preconditioner $P$, this class implements a new matrix
 * with the matrix-vector product $PA$. It needs an auxiliary vector for that.
 *
 * By this time, this is considered a temporary object to be plugged
 * into eigenvalue solvers. Therefore, no @p{SmartPointer} is used for
 * @p{A} and @p{P}.
 *
 * @author Guido Kanschat, 2000
 */
template<class MATRIX, class PRECOND, class VECTOR>
class PreconditionedMatrix : public Subscriptor
{
  public:
				     /**
				      * Constructor. Provide matrix,
				      * preconditioner and a memory
				      * pool to obtain the auxiliary
				      * vector.
				      */
    PreconditionedMatrix (const MATRIX&          A,
			  const PRECOND&         P,
			  VectorMemory<VECTOR>&  mem);

				     /**
				      * Preconditioned
				      * matrix-vector-product.
				      */
    void vmult (VECTOR& dst, const VECTOR& src) const;

				     /**
				      * Transposed preconditioned
				      * matrix-vector-product.
				      */
    void Tvmult (VECTOR& dst, const VECTOR& src) const;

				     /**
				      * Residual $b-PAx$.
				      */
    double residual (VECTOR& dst, const VECTOR& src, const VECTOR& rhs) const;

  private:
				     /**
				      * Storage for the matrix.
				      */
    const MATRIX& A;
				     /**
				      * Storage for preconditioner.
				      */
    const PRECOND& P;
				     /**
				      * Memory pool for vectors.
				      */
    VectorMemory<VECTOR>& mem;
};


    

/* ---------------------------------- Inline functions ------------------- */

template<class VECTOR>
inline void
PreconditionIdentity::vmult (VECTOR& dst, const VECTOR& src) const
{
  dst = src;
}



template<class VECTOR>
inline void
PreconditionIdentity::Tvmult (VECTOR& dst, const VECTOR& src) const
{
  dst = src;
}

//----------------------------------------------------------------------//

template <class MATRIX>
inline void
PreconditionRelaxation<MATRIX>::initialize (const MATRIX &rA,
					    const double  o)
{
  A = &rA;
  omega = o;
}

//----------------------------------------------------------------------//

template <class MATRIX>
template<class VECTOR>
inline void
PreconditionJacobi<MATRIX>::vmult (VECTOR& dst, const VECTOR& src) const
{
  Assert (A!=0, ExcNotInitialized());
  A->precondition_Jacobi (dst, src, omega);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionJacobi<MATRIX>::Tvmult (VECTOR& dst, const VECTOR& src) const
{
  Assert (A!=0, ExcNotInitialized());
  A->precondition_Jacobi (dst, src, omega);
}


//----------------------------------------------------------------------//

template <class MATRIX>
template<class VECTOR>
inline void
PreconditionSOR<MATRIX>::vmult (VECTOR& dst, const VECTOR& src) const
{
  Assert (A!=0, ExcNotInitialized());
  A->precondition_SOR (dst, src, omega);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionSOR<MATRIX>::Tvmult (VECTOR& dst, const VECTOR& src) const
{
  Assert (A!=0, ExcNotInitialized());
  A->precondition_TSOR (dst, src, omega);
}


//----------------------------------------------------------------------//

template <class MATRIX>
template<class VECTOR>
inline void
PreconditionSSOR<MATRIX>::vmult (VECTOR& dst, const VECTOR& src) const
{
  Assert (A!=0, ExcNotInitialized());
  A->precondition_SSOR (dst, src, omega);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionSSOR<MATRIX>::Tvmult (VECTOR& dst, const VECTOR& src) const
{
  Assert (A!=0, ExcNotInitialized());
  A->precondition_SSOR (dst, src, omega);
}


//----------------------------------------------------------------------//


template<class MATRIX, class VECTOR>
PreconditionUseMatrix<MATRIX,VECTOR>::PreconditionUseMatrix(const MATRIX       &M,
							    const function_ptr  method)
		:
		matrix(M), precondition(method)
{}



template<class MATRIX, class VECTOR>
void
PreconditionUseMatrix<MATRIX,VECTOR>::vmult (VECTOR& dst,
					     const VECTOR& src) const
{
  (matrix.*precondition)(dst, src);
}

//////////////////////////////////////////////////////////////////////

template<class SOLVER, class MATRIX, class PRECONDITION>
PreconditionLACSolver<SOLVER,MATRIX,PRECONDITION>
::PreconditionLACSolver ()
		:
		solver(0), matrix(0), precondition(0)
{}


template<class SOLVER, class MATRIX, class PRECONDITION>
void
PreconditionLACSolver<SOLVER,MATRIX,PRECONDITION>
::initialize (SOLVER& s,
	      const MATRIX& m,
	      const PRECONDITION& p)
{
  solver = s;
  matrix = m;
  precondition = p;
}


template<class SOLVER, class MATRIX, class PRECONDITION>
template<class VECTOR>
void
PreconditionLACSolver<SOLVER,MATRIX,PRECONDITION>::vmult (VECTOR& dst,
							  const VECTOR& src) const
{
  Assert (solver !=0 && matrix != 0 && precondition != 0,
	  ExcNotInitialized());
  
  solver->solve(*matrix, dst, src, *precondition);
}

//////////////////////////////////////////////////////////////////////


template<class MATRIX, class PRECOND, class VECTOR>
inline
PreconditionedMatrix<MATRIX, PRECOND, VECTOR>
::PreconditionedMatrix (const MATRIX&  A,
			const PRECOND& P,
			VectorMemory<VECTOR>&  mem):
		A(A), P(P), mem(mem)
{}


template<class MATRIX, class PRECOND, class VECTOR>
inline void
PreconditionedMatrix<MATRIX, PRECOND, VECTOR>
::vmult (VECTOR& dst,
	 const VECTOR& src) const
{
  VECTOR* h = mem.alloc();
  h->reinit(src);
  A.vmult(*h, src);
  P.vmult(dst, *h);
  mem.free(h);
}



template<class MATRIX, class PRECOND, class VECTOR>
inline void
PreconditionedMatrix<MATRIX, PRECOND, VECTOR>
::Tvmult (VECTOR& dst,
	 const VECTOR& src) const
{
  VECTOR* h = mem.alloc();
  h->reinit(src);
  A.Tvmult(*h, src);
  P.Tvmult(dst, *h);
  mem.free(h);
}



template<class MATRIX, class PRECOND, class VECTOR>
inline double
PreconditionedMatrix<MATRIX, PRECOND, VECTOR>
::residual (VECTOR& dst,
	    const VECTOR& src,
	    const VECTOR& rhs) const
{
  VECTOR* h = mem.alloc();
  h->reinit(src);
  A.vmult(*h, src);
  P.vmult(dst, *h);
  mem.free(h);
  dst.sadd(-1.,1.,rhs);
  return dst.l2_norm ();
}

#endif

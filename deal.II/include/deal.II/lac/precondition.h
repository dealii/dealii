//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__precondition_h
#define __deal2__precondition_h

// This file contains simple preconditioners.

#include <deal.II/base/config.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/lac/tridiagonal_matrix.h>
#include <deal.II/lac/vector_memory.h>

DEAL_II_NAMESPACE_OPEN

template <typename number> class Vector;
template <typename number> class SparseMatrix;

/*! @addtogroup Preconditioners
 *@{
 */


/**
 * No preconditioning.  This class helps you, if you want to use a
 * linear solver without preconditioning. All solvers in LAC require a
 * preconditioner. Therefore, you must use the identity provided here
 * to avoid preconditioning. It can be used in the following way:
 *
 @verbatim
  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              cg (solver_control);
  cg.solve (system_matrix, solution, system_rhs,
	    PreconditionIdentity());
 @endverbatim
 *
 * See the step-3 tutorial program for an example and
 * additional explanations.
 *
 * Alternatively, the IdentityMatrix class can be used to precondition
 * in this way.
 *
 * @author Guido Kanschat, 1999
 */
class PreconditionIdentity : public Subscriptor
{
  public:

				     /**
				      * This function is only
                                      * present to
                                      * provide the interface of
                                      * a preconditioner to be
                                      * handed to a smoother.
                                      * This does nothing.
				      */
  struct AdditionalData
  {
                                       /**
					* Constructor.
					*/
    AdditionalData (){}
  };

				     /**
				      * The matrix
				      * argument is ignored and here
				      * just for compatibility with
				      * more complex preconditioners.
				      */
  template <class MATRIX>
  void initialize (const MATRIX         &matrix,
		   const AdditionalData &additional_data = AdditionalData());

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
				      * vmult().
				      */
    template<class VECTOR>
    void Tvmult (VECTOR&, const VECTOR&) const;

				     /**
				      * Apply preconditioner, adding to the previous value.
				      */
    template<class VECTOR>
    void vmult_add (VECTOR&, const VECTOR&) const;

				     /**
				      * Apply transpose
				      * preconditioner, adding. Since this is
				      * the identity, this function is
				      * the same as
				      * vmult_add().
				      */
    template<class VECTOR>
    void Tvmult_add (VECTOR&, const VECTOR&) const;

				     /**
				      * This function is only
                                      * present to
                                      * provide the interface of
                                      * a preconditioner to be
                                      * handed to a smoother.
                                      * This does nothing.
				      */
    void clear (){}
};



/**
 * Preconditioning with Richardson's method. This preconditioner just
 * scales the vector with a constant relaxation factor provided by the
 * AdditionalData object.
 *
 * In Krylov-space methods, this preconditioner should not have any
 * effect. Using SolverRichardson, the two relaxation parameters will
 * be just multiplied. Still, this class is useful in multigrid
 * smoother objects (MGSmootherRelaxation).
 *
 * @author Guido Kanschat, 2005
 */
class PreconditionRichardson : public Subscriptor
{
  public:
				     /**
				      * Parameters for Richardson
				      * preconditioner.
				      */
    class AdditionalData
    {
      public:
					 /**
					  * Constructor. Block size
					  * must be given since there
					  * is no reasonable default
					  * parameter.
					  */
	AdditionalData (const double relaxation = 1.);

					 /**
					  * Relaxation parameter.
					  */
	double relaxation;
    };

				     /**
				      * Constructor, sets the
				      * relaxation parameter to its
				      * default.
				      */
    PreconditionRichardson();

				     /**
				      * Change the relaxaton parameter.
				      */
    void initialize (const AdditionalData &parameters);

				     /**
				      * Change the relaxaton parameter
				      * in a way consistent with other
				      * preconditioners. The matrix
				      * argument is ignored and here
				      * just for compatibility with
				      * more complex preconditioners.
				      */
    template <class MATRIX>
    void initialize (const MATRIX&,
		     const AdditionalData &parameters);

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
				      * vmult().
				      */
    template<class VECTOR>
    void Tvmult (VECTOR&, const VECTOR&) const;
				     /**
				      * Apply preconditioner, adding to the previous value.
				      */
    template<class VECTOR>
    void vmult_add (VECTOR&, const VECTOR&) const;

				     /**
				      * Apply transpose
				      * preconditioner, adding. Since this is
				      * the identity, this function is
				      * the same as
				      * vmult_add().
				      */
    template<class VECTOR>
    void Tvmult_add (VECTOR&, const VECTOR&) const;

				     /**
				      * This function is only
				      * present to
				      * provide the interface of
				      * a preconditioner to be
				      * handed to a smoother.
				      * This does nothing.
				      */
    void clear (){}

  private:
				     /**
				      * The relaxation parameter
				      * multiplied with the vectors.
				      */
    double relaxation;
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
 * parameter, so please use PreconditionRelaxation for these.
 *
 * You will usually not want to create a named object of this type,
 * although possible. The most common use is like this:
 * @code
 *    SolverGMRES<SparseMatrix<double>,
 *                Vector<double> >      gmres(control,memory,500);
 *
 *    gmres.solve (matrix, solution, right_hand_side,
 *		   PreconditionUseMatrix<SparseMatrix<double>,Vector<double> >
 *		   (matrix,&SparseMatrix<double>::template precondition_Jacobi<double>));
 * @endcode
 * This creates an unnamed object to be passed as the fourth parameter to
 * the solver function of the SolverGMRES class. It assumes that the
 * SparseMatrix class has a function <tt>precondition_Jacobi</tt> taking two
 * vectors (source and destination) as parameters (Actually, there is no
 * function like that, the existing function takes a third parameter,
 * denoting the relaxation parameter; this example is therefore only meant to
 * illustrate the general idea).
 *
 * Note that due to the default template parameters, the above example
 * could be written shorter as follows:
 * @code
 *    ...
 *    gmres.solve (matrix, solution, right_hand_side,
 *		   PreconditionUseMatrix<>
 *		     (matrix,&SparseMatrix<double>::template precondition_Jacobi<double>));
 * @endcode
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
    const MATRIX &matrix;

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
				      * Class for parameters.
				      */
    class AdditionalData
    {
      public:
					 /**
					  * Constructor.
					  */
	AdditionalData (const double relaxation = 1.);

					 /**
					  * Relaxation parameter.
					  */
	double relaxation;
    };

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
    void initialize (const MATRIX &A,
		     const AdditionalData & parameters = AdditionalData());

				     /**
				      * Release the matrix and reset
				      * its pointer.
				      */
    void clear();

  protected:
				     /**
				      * Pointer to the matrix object.
				      */
    SmartPointer<const MATRIX, PreconditionRelaxation<MATRIX> > A;

				     /**
				      * Relaxation parameter.
				      */
    double relaxation;
};



/**
 * Jacobi preconditioner using matrix built-in function.  The
 * <tt>MATRIX</tt> class used is required to have a function
 * <tt>precondition_Jacobi(VECTOR&, const VECTOR&, double</tt>)
 *
 * @code
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
 * @endcode
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
				      * vmult().
				      */
    template<class VECTOR>
    void Tvmult (VECTOR&, const VECTOR&) const;

				     /**
				      * Perform one step of the
				      * preconditioned Richardson
				      * iteration.
				      */
    template<class VECTOR>
    void step (VECTOR& x, const VECTOR& rhs) const;

				     /**
				      * Perform one transposed step of
				      * the preconditioned Richardson
				      * iteration.
				      */
    template<class VECTOR>
    void Tstep (VECTOR& x, const VECTOR& rhs) const;
};


/**
 * SOR preconditioner using matrix built-in function.
 *
 * Assuming the matrix <i>A = D + L + U</i> is split into its diagonal
 * <i>D</i> as well as the strict lower and upper triangles <i>L</i>
 * and <i>U</i>, then the SOR preconditioner with relaxation parameter
 * <i>r</i> is
 * @f[
 *  P^{-1} = r (D+rL)^{-1}.
 * @f]
 * It is this operator <i>P<sup>-1</sup></i>, which is implemented by
 * vmult() through forward substitution. Analogously, Tvmult()
 * implements the operation of <i>r(D+rU)<sup>-1</sup></i>.
 *
 * The SOR iteration itself can be directly written as
 * @f[
 *  x^{k+1} = x^k - r D^{-1} \bigl(L x^{k+1} + U x^k - b\bigr).
 * @f]
 * Using the right hand side <i>b</i> and the previous iterate
 * <i>x</i>, this is the operation implemented by step().
 *
 * The MATRIX
 * class used is required to have functions
 * <tt>precondition_SOR(VECTOR&, const VECTOR&, double)</tt> and
 * <tt>precondition_TSOR(VECTOR&, const VECTOR&, double)</tt>.
 *
 * @code
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
 * @endcode
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

				     /**
				      * Perform one step of the
				      * preconditioned Richardson
				      * iteration.
				      */
    template<class VECTOR>
    void step (VECTOR& x, const VECTOR& rhs) const;

				     /**
				      * Perform one transposed step of
				      * the preconditioned Richardson
				      * iteration.
				      */
    template<class VECTOR>
    void Tstep (VECTOR& x, const VECTOR& rhs) const;
};



/**
 * SSOR preconditioner using matrix built-in function.  The
 * <tt>MATRIX</tt> class used is required to have a function
 * <tt>precondition_SSOR(VECTOR&, const VECTOR&, double)</tt>
 *
 * @code
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
 * @endcode
 *
 * @author Guido Kanschat, 2000
 */
template <class MATRIX = SparseMatrix<double> >
class PreconditionSSOR : public PreconditionRelaxation<MATRIX>
{
  public:

				     /**
				      * A typedef to the base class.
				      */
    typedef PreconditionRelaxation<MATRIX> BaseClass;


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
    void initialize (const MATRIX &A,
		     const typename BaseClass::AdditionalData &parameters = typename BaseClass::AdditionalData());

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
				      * vmult().
				      */
    template<class VECTOR>
    void Tvmult (VECTOR&, const VECTOR&) const;


				     /**
				      * Perform one step of the
				      * preconditioned Richardson
				      * iteration
				      */
    template<class VECTOR>
    void step (VECTOR& x, const VECTOR& rhs) const;

				     /**
				      * Perform one transposed step of
				      * the preconditioned Richardson
				      * iteration.
				      */
    template<class VECTOR>
    void Tstep (VECTOR& x, const VECTOR& rhs) const;

  private:
				     /**
				      * An array that stores for each matrix
				      * row where the first position after
				      * the diagonal is located.
				      */
    std::vector<unsigned int> pos_right_of_diagonal;
};


/**
 * Permuted SOR preconditioner using matrix built-in function.  The
 * <tt>MATRIX</tt> class used is required to have functions
 * <tt>PSOR(VECTOR&, const VECTOR&, double)</tt> and
 * <tt>TPSOR(VECTOR&, const VECTOR&, double)</tt>.
 *
 * @code
 *     // Declare related objects
 *
 * SparseMatrix<double> A;
 * Vector<double> x;
 * Vector<double> b;
 * SolverCG<> solver(...);
 *
 * //...initialize and build A
 *
 * std::vector<unsigned int> permutation(x.size());
 * std::vector<unsigned int> inverse_permutation(x.size());
 *
 * //...fill permutation and its inverse with reasonable values
 *
 *     // Define and initialize preconditioner
 *
 * PreconditionPSOR<SparseMatrix<double> > precondition;
 * precondition.initialize (A, permutation, inverse_permutation, .6);
 *
 * solver.solve (A, x, b, precondition);
 * @endcode
 *
 * @author Guido Kanschat, 2003
 */
template <class MATRIX = SparseMatrix<double> >
class PreconditionPSOR : public PreconditionRelaxation<MATRIX>
{
  public:
				     /**
				      * Initialize matrix and
				      * relaxation parameter. The
				      * matrix is just stored in the
				      * preconditioner object.
				      *
				      * The permutation vector is
				      * stored as a
				      * pointer. Therefore, it has to
				      * be assured that the lifetime
				      * of the vector exceeds the
				      * lifetime of the
				      * preconditioner.
				      *
				      * The relaxation parameter
				      * should be larger than zero and
				      * smaller than 2 for numerical
				      * reasons. It defaults to 1.
				      */
    void initialize (const MATRIX &A,
		     const std::vector<unsigned int> &permutation,
		     const std::vector<unsigned int> &inverse_permutation,
		     const typename PreconditionRelaxation<MATRIX>::AdditionalData &
		     parameters = typename PreconditionRelaxation<MATRIX>::AdditionalData());

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
  private:
				     /**
				      * Storage for the permutation vector.
				      */
    const std::vector<unsigned int>* permutation;
				     /**
				      * Storage for the inverse
				      * permutation vector.
				      */
    const std::vector<unsigned int>* inverse_permutation;
};



/**
 * @deprecated This class has been superseded by IterativeInverse,
 * which is more flexible and easier to use.
 *
 * Preconditioner using an iterative solver.  This preconditioner uses
 * a fully initialized LAC iterative solver for the approximate
 * inverse of the matrix. Naturally, this solver needs another
 * preconditionig method.
 *
 * Usually, the use of ReductionControl is preferred over the use of
 * the basic SolverControl in defining this solver.
 *
 * Krylov space methods like SolverCG or SolverBicgstab
 * become inefficient if soution down to machine accuracy is
 * needed. This is due to the fact, that round-off errors spoil the
 * orthogonality of the vector sequences. Therefore, a nested
 * iteration of two methods is proposed: The outer method is
 * SolverRichardson, since it is robust with respect to round-of
 * errors. The inner loop is an appropriate Krylov space method, since
 * it is fast.
 *
 * @code
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
 * @endcode
 *
 * Each time we call the inner loop, reduction of the residual by a
 * factor <tt>1.e-2</tt> is attempted. Since the right hand side vector of
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
    SmartPointer<SOLVER,PreconditionLACSolver<SOLVER,MATRIX,PRECONDITION> > solver;

				     /**
				      * The matrix in use.
				      */
    SmartPointer<const MATRIX,PreconditionLACSolver<SOLVER,MATRIX,PRECONDITION> > matrix;

				     /**
				      * The preconditioner to use.
				      */
    SmartPointer<const PRECONDITION,PreconditionLACSolver<SOLVER,MATRIX,PRECONDITION> > precondition;
};



/**
 * @deprecated Use ProductMatrix instead.
 *
 * Matrix with preconditioner.
 * Given a matrix $A$ and a preconditioner $P$, this class implements a new matrix
 * with the matrix-vector product $PA$. It needs an auxiliary vector for that.
 *
 * By this time, this is considered a temporary object to be plugged
 * into eigenvalue solvers. Therefore, no SmartPointer is used for
 * <tt>A</tt> and <tt>P</tt>.
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
    PreconditionedMatrix (const MATRIX &         A,
			  const PRECOND &        P,
			  VectorMemory<VECTOR> & mem);

				     /**
				      * Preconditioned
				      * matrix-vector-product.
				      */
    void vmult (VECTOR &dst, const VECTOR &src) const;

				     /**
				      * Transposed preconditioned
				      * matrix-vector-product.
				      */
    void Tvmult (VECTOR &dst, const VECTOR &src) const;

				     /**
				      * Residual $b-PAx$.
				      */
    double residual (VECTOR &dst, const VECTOR &src, const VECTOR &rhs) const;

  private:
				     /**
				      * Storage for the matrix.
				      */
    const MATRIX &A;
				     /**
				      * Storage for preconditioner.
				      */
    const PRECOND &P;
				     /**
				      * Memory pool for vectors.
				      */
    VectorMemory<VECTOR> &mem;
};



/**
 * Preconditioning with a Chebyshev polynomial for symmetric positive
 * definite matrices. This preconditioner is similar to a Jacobi
 * preconditioner if the degree variable is set to one, otherwise some
 * higher order polynomial corrections are used. This preconditioner needs
 * access to the diagonal of the matrix its acts on and needs a respective
 * <tt>vmult</tt> implemention. However, it does not need to explicitly know
 * the matrix entries.
 *
 * This class is useful e.g. in multigrid smoother objects, since it is
 * trivially %parallel (assuming that matrix-vector products are %parallel).
 *
 * @author Martin Kronbichler, 2009
 */
template <class MATRIX=SparseMatrix<double>, class VECTOR=Vector<double> >
class PreconditionChebyshev : public Subscriptor
{
  public:
				     /**
				      * Standardized data struct to
				      * pipe additional parameters
				      * to the preconditioner.
				      */
    struct AdditionalData
    {
					 /**
					  * Constructor.
					  */
	AdditionalData (const unsigned int degree              = 0,
			const double       smoothing_range     = 0.,
			const bool         nonzero_starting    = false,
			const unsigned int eig_cg_n_iterations = 8,
			const double       eig_cg_residual     = 1e-2);

					 /**
					  * This determines the degree of the
					  * Chebyshev polynomial. The degree
					  * of the polynomial gives the number
					  * of matrix-vector products to be
					  * performed for one application of
					  * the vmult() operation. Degree zero
					  * corresponds to a damped Jacobi
					  * method.
					  */
	unsigned int degree;

					 /**
					  * This sets the range between the
					  * largest eigenvalue in the matrix
					  * and the smallest eigenvalue to be
					  * treated. If the parameter is zero,
					  * an estimate for the largest and
					  * for the smallest eigenvalue will
					  * be calculated
					  * internally. Otherwise, the
					  * Chebyshev polynomial will act in
					  * the interval
					  * $[\lambda_\mathrm{max}/
					  * \tt{smoothing\_range},
					  * \lambda_\mathrm{max}]$, where
					  * $\lambda_\mathrm{max}$ is an
					  * estimate of the maximum eigenvalue
					  * of the matrix. A choice of
					  * <tt>smoothing_range</tt> between 5
					  * and 20 is useful in case the
					  * preconditioner is used as a
					  * smoother in multigrid.
					  */
	double smoothing_range;

					 /**
					  * When this flag is set to
					  * <tt>true</tt>, it enables the
					  * method <tt>vmult(dst, src)</tt> to
					  * use non-zero data in the vector
					  * <tt>dst</tt>, appending to it the
					  * Chebyshev corrections. This can be
					  * useful in some situations
					  * (e.g. when used for high-frequency
					  * error smoothing in a multigrid
					  * algorithm), but not the way the
					  * solver classes expect a
					  * preconditioner to work (where one
					  * ignores the content in
					  * <tt>dst</tt> for the
					  * preconditioner application).
					  */
	bool nonzero_starting;

					 /**
					  * Maximum number of CG iterations
					  * performed for finding the maximum
					  * eigenvalue.
					  */
	unsigned int eig_cg_n_iterations;

					 /**
					  * Tolerance for CG iterations
					  * performed for finding the maximum
					  * eigenvalue.
					  */
	double eig_cg_residual;

					 /**
					  * Stores the inverse of the diagonal
					  * of the underlying matrix.
					  */
	VECTOR matrix_diagonal_inverse;
    };

    PreconditionChebyshev ();

				     /**
				      * Initialize function. Takes the
				      * matrix which is used to form the
				      * preconditioner, and additional
				      * flags if there are any. This
				      * function works only if the input
				      * matrix has an operator
				      * <tt>el(i,i)</tt> for accessing all
				      * the elements in the
				      * diagonal. Alternatively, the
				      * diagonal can be supplied with the
				      * help of the AdditionalData field.
				      *
				      * This function calculates an
				      * estimate of the eigenvalue range
				      * of the matrix weighted by its
				      * diagonal using a modified CG
				      * iteration.
				      */
    void initialize (const MATRIX         &matrix,
		     const AdditionalData &additional_data = AdditionalData());

				     /**
				      * Computes the action of the
				      * preconditioner on <tt>src</tt>,
				      * storing the result in
				      * <tt>dst</tt>.
				      */
    void vmult (VECTOR       &dst,
		const VECTOR &src) const;

				     /**
				      * Computes the action of the
				      * transposed preconditioner on
				      * <tt>src</tt>, storing the result
				      * in <tt>dst</tt>.
				      */
    void Tvmult (VECTOR       &dst,
		 const VECTOR &src) const;

				     /**
				      * Resets the preconditioner.
				      */
    void clear ();

  private:

				     /**
				      * A pointer to the underlying
				      * matrix.
				      */
    SmartPointer<const MATRIX,PreconditionChebyshev<MATRIX,VECTOR> > matrix_ptr;

				     /**
				      * Internal vector used for the
				      * <tt>vmult</tt> operation.
				      */
    mutable VECTOR update1;

				     /**
				      * Internal vector used for the
				      * <tt>vmult</tt> operation.
				      */
    mutable VECTOR update2;

				     /**
				      * Stores the additional data
				      * provided to the initialize
				      * function.
				      */
    AdditionalData data;

				     /**
				      * Average of the largest and
				      * smallest eigenvalue under
				      * consideration.
				      */
    double theta;

				     /**
				      * Half the interval length between
				      * the largest and smallest
				      * eigenvalue under consideration.
				      */
    double delta;

				     /**
				      * Stores whether the preconditioner
				      * has been set up.
				      */
    bool is_initialized;
};



/*@}*/
/* ---------------------------------- Inline functions ------------------- */

#ifndef DOXYGEN

template <class MATRIX>
inline void
PreconditionIdentity::initialize (
  const MATRIX&,
  const PreconditionIdentity::AdditionalData&)
{}


template<class VECTOR>
inline void
PreconditionIdentity::vmult (VECTOR &dst, const VECTOR &src) const
{
  dst = src;
}



template<class VECTOR>
inline void
PreconditionIdentity::Tvmult (VECTOR &dst, const VECTOR &src) const
{
  dst = src;
}

template<class VECTOR>
inline void
PreconditionIdentity::vmult_add (VECTOR &dst, const VECTOR &src) const
{
  dst.add(src);
}



template<class VECTOR>
inline void
PreconditionIdentity::Tvmult_add (VECTOR &dst, const VECTOR &src) const
{
  dst.add(src);
}

//---------------------------------------------------------------------------

inline
PreconditionRichardson::AdditionalData::AdditionalData (
  const double relaxation)
		:
		relaxation(relaxation)
{}


inline
PreconditionRichardson::PreconditionRichardson ()
		:
		relaxation(0)
{
  AdditionalData add_data;
  relaxation=add_data.relaxation;
}



inline void
PreconditionRichardson::initialize (
  const PreconditionRichardson::AdditionalData &parameters)
{
  relaxation = parameters.relaxation;
}



template <class MATRIX>
inline void
PreconditionRichardson::initialize (
  const MATRIX&,
  const PreconditionRichardson::AdditionalData &parameters)
{
  relaxation = parameters.relaxation;
}



template<class VECTOR>
inline void
PreconditionRichardson::vmult (VECTOR &dst, const VECTOR &src) const
{
  dst.equ(relaxation,src);
}



template<class VECTOR>
inline void
PreconditionRichardson::Tvmult (VECTOR &dst, const VECTOR &src) const
{
  dst.equ(relaxation,src);
}

template<class VECTOR>
inline void
PreconditionRichardson::vmult_add (VECTOR &dst, const VECTOR &src) const
{
  dst.add(relaxation,src);
}



template<class VECTOR>
inline void
PreconditionRichardson::Tvmult_add (VECTOR &dst, const VECTOR &src) const
{
  dst.add(relaxation,src);
}

//---------------------------------------------------------------------------

template <class MATRIX>
inline void
PreconditionRelaxation<MATRIX>::initialize (const MATRIX &rA,
					    const AdditionalData &parameters)
{
  A = &rA;
  relaxation = parameters.relaxation;
}


template <class MATRIX>
inline void
PreconditionRelaxation<MATRIX>::clear ()
{
  A = 0;
}


//---------------------------------------------------------------------------

template <class MATRIX>
template<class VECTOR>
inline void
PreconditionJacobi<MATRIX>::vmult (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  this->A->precondition_Jacobi (dst, src, this->relaxation);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionJacobi<MATRIX>::Tvmult (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  this->A->precondition_Jacobi (dst, src, this->relaxation);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionJacobi<MATRIX>::step (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  this->A->Jacobi_step (dst, src, this->relaxation);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionJacobi<MATRIX>::Tstep (VECTOR &dst, const VECTOR &src) const
{
  step (dst, src);
}



//---------------------------------------------------------------------------

template <class MATRIX>
template<class VECTOR>
inline void
PreconditionSOR<MATRIX>::vmult (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  this->A->precondition_SOR (dst, src, this->relaxation);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionSOR<MATRIX>::Tvmult (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  this->A->precondition_TSOR (dst, src, this->relaxation);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionSOR<MATRIX>::step (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  this->A->SOR_step (dst, src, this->relaxation);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionSOR<MATRIX>::Tstep (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  this->A->TSOR_step (dst, src, this->relaxation);
}



//---------------------------------------------------------------------------

template <class MATRIX>
inline void
PreconditionSSOR<MATRIX>::initialize (const MATRIX &rA,
				      const typename BaseClass::AdditionalData &parameters)
{
  this->PreconditionRelaxation<MATRIX>::initialize (rA, parameters);

				   // in case we have a SparseMatrix class,
				   // we can extract information about the
				   // diagonal.
  const SparseMatrix<typename MATRIX::value_type> * mat =
    dynamic_cast<const SparseMatrix<typename MATRIX::value_type> *>(&*this->A);

				   // calculate the positions first after
				   // the diagonal.
  if (mat != 0)
    {
      const std::size_t  * rowstart_ptr =
	mat->get_sparsity_pattern().get_rowstart_indices();
      const unsigned int * const colnums =
	mat->get_sparsity_pattern().get_column_numbers();
      const unsigned int n = this->A->n();
      pos_right_of_diagonal.resize(n);
      for (unsigned int row=0; row<n; ++row, ++rowstart_ptr)
	{
				       // find the first element in this line
				       // which is on the right of the diagonal.
				       // we need to precondition with the
				       // elements on the left only.
				       // note: the first entry in each
				       // line denotes the diagonal element,
				       // which we need not check.
	  pos_right_of_diagonal[row] =
	    Utilities::lower_bound (&colnums[*rowstart_ptr+1],
			      &colnums[*(rowstart_ptr+1)],
			      row)
	    - colnums;
	}
    }
}


template <class MATRIX>
template<class VECTOR>
inline void
PreconditionSSOR<MATRIX>::vmult (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  this->A->precondition_SSOR (dst, src, this->relaxation, pos_right_of_diagonal);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionSSOR<MATRIX>::Tvmult (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  this->A->precondition_SSOR (dst, src, this->relaxation, pos_right_of_diagonal);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionSSOR<MATRIX>::step (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  this->A->SSOR_step (dst, src, this->relaxation);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionSSOR<MATRIX>::Tstep (VECTOR &dst, const VECTOR &src) const
{
  step (dst, src);
}



//---------------------------------------------------------------------------

template <class MATRIX>
inline void
PreconditionPSOR<MATRIX>::initialize (
  const MATRIX &rA,
  const std::vector<unsigned int> &p,
  const std::vector<unsigned int> &ip,
  const typename PreconditionRelaxation<MATRIX>::AdditionalData &parameters)
{
  permutation = &p;
  inverse_permutation = &ip;
  PreconditionRelaxation<MATRIX>::initialize(rA, parameters);
}


template <class MATRIX>
template<class VECTOR>
inline void
PreconditionPSOR<MATRIX>::vmult (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  dst = src;
  this->A->PSOR (dst, *permutation, *inverse_permutation, this->relaxation);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionPSOR<MATRIX>::Tvmult (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  dst = src;
  this->A->TPSOR (dst, *permutation, *inverse_permutation, this->relaxation);
}


//---------------------------------------------------------------------------


template<class MATRIX, class VECTOR>
PreconditionUseMatrix<MATRIX,VECTOR>::PreconditionUseMatrix(const MATRIX       &M,
							    const function_ptr  method)
		:
		matrix(M), precondition(method)
{}



template<class MATRIX, class VECTOR>
void
PreconditionUseMatrix<MATRIX,VECTOR>::vmult (VECTOR &dst,
					     const VECTOR &src) const
{
  (matrix.*precondition)(dst, src);
}

//---------------------------------------------------------------------------

template<class MATRIX>
inline
PreconditionRelaxation<MATRIX>::AdditionalData::
AdditionalData (const double relaxation)
		:
		relaxation (relaxation)
{}



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
::initialize (SOLVER &s,
	      const MATRIX &m,
	      const PRECONDITION &p)
{
  solver = &s;
  matrix = &m;
  precondition = &p;
}


template<class SOLVER, class MATRIX, class PRECONDITION>
template<class VECTOR>
void
PreconditionLACSolver<SOLVER,MATRIX,PRECONDITION>::vmult (VECTOR &dst,
							  const VECTOR &src) const
{
  Assert (solver !=0 && matrix != 0 && precondition != 0,
	  ExcNotInitialized());

  solver->solve(*matrix, dst, src, *precondition);
}

//////////////////////////////////////////////////////////////////////


template<class MATRIX, class PRECOND, class VECTOR>
inline
PreconditionedMatrix<MATRIX, PRECOND, VECTOR>
::PreconditionedMatrix (const MATRIX & A,
			const PRECOND &P,
			VectorMemory<VECTOR> & mem):
		A(A), P(P), mem(mem)
{}


template<class MATRIX, class PRECOND, class VECTOR>
inline void
PreconditionedMatrix<MATRIX, PRECOND, VECTOR>
::vmult (VECTOR &dst,
	 const VECTOR &src) const
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
::Tvmult (VECTOR &dst,
	 const VECTOR &src) const
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
::residual (VECTOR &dst,
	    const VECTOR &src,
	    const VECTOR &rhs) const
{
  VECTOR* h = mem.alloc();
  h->reinit(src);
  A.vmult(*h, src);
  P.vmult(dst, *h);
  mem.free(h);
  dst.sadd(-1.,1.,rhs);
  return dst.l2_norm ();
}



//---------------------------------------------------------------------------

template <class MATRIX, class VECTOR>
inline
PreconditionChebyshev<MATRIX,VECTOR>::AdditionalData::
AdditionalData (const unsigned int degree,
		const double       smoothing_range,
		const bool         nonzero_starting,
		const unsigned int eig_cg_n_iterations,
		const double       eig_cg_residual)
                :
                degree  (degree),
		smoothing_range (smoothing_range),
		nonzero_starting (nonzero_starting),
		eig_cg_n_iterations (eig_cg_n_iterations),
		eig_cg_residual (eig_cg_residual)
{}



template <class MATRIX, class VECTOR>
inline
PreconditionChebyshev<MATRIX,VECTOR>::PreconditionChebyshev ()
  :
    is_initialized (false)
{}



template <class MATRIX, class VECTOR>
inline
void
PreconditionChebyshev<MATRIX,VECTOR>::initialize (const MATRIX &matrix,
						  const AdditionalData &additional_data)
{
  Assert (matrix.m() == matrix.n(), ExcMessage("Matrix not quadratic."));
  matrix_ptr = &matrix;
  data = additional_data;
  if (data.matrix_diagonal_inverse.size() != matrix.m())
    {
      data.matrix_diagonal_inverse.reinit(matrix.m());
      for (unsigned int i=0; i<matrix.m(); ++i)
	data.matrix_diagonal_inverse(i) = 1./matrix.el(i,i);
    }
  update1.reinit (data.matrix_diagonal_inverse, true);
  update2.reinit (data.matrix_diagonal_inverse, true);


				 // calculate largest eigenvalue using a
				 // hand-tuned CG iteration on the matrix
				 // weighted by its diagonal. we start
				 // with a vector that consists of ones
				 // only, weighted by the length.
				 //
				 // TODO: can we obtain this with the
				 // regular CG implementation? we would need
				 // to read the logfile in that case, which
				 // does not seem feasible.
  double max_eigenvalue, min_eigenvalue;
  {
    double eigen_beta_alpha = 0;

    std::vector<double> diagonal;
    std::vector<double> offdiagonal;

    VECTOR rhs, g;
    rhs.reinit(data.matrix_diagonal_inverse, true);
    rhs = 1./std::sqrt(static_cast<double>(matrix.m()));
    g.reinit(data.matrix_diagonal_inverse, true);

    unsigned int it=0;
    double res,gh,alpha,beta;

    g.equ(-1.,rhs);
    res = g.l2_norm();
    update2.equ (-1., g);
    gh = res*res;

    while (true)
      {
	it++;
	matrix.vmult (update1, update2);
	update1.scale (data.matrix_diagonal_inverse);
	alpha = update2 * update1;
	alpha = gh/alpha;
	g.add (alpha, update1);
	res = g.l2_norm();

	if (it > data.eig_cg_n_iterations || res < data.eig_cg_residual)
	  break;

	beta = gh;
	gh = res*res;
	beta = gh/beta;
	update2.sadd (beta, -1., g);

	diagonal.push_back (1./alpha + eigen_beta_alpha);
	eigen_beta_alpha = beta/alpha;
	offdiagonal.push_back(std::sqrt(beta)/alpha);
      }

    TridiagonalMatrix<double> T(diagonal.size(), true);
    for (unsigned int i=0;i<diagonal.size();++i)
      {
	T(i,i) = diagonal[i];
	if (i< diagonal.size()-1)
	  T(i,i+1) = offdiagonal[i];
      }
    T.compute_eigenvalues();
    min_eigenvalue = T.eigenvalue(0);
    max_eigenvalue = T.eigenvalue(T.n()-1);
  }

				 // include a safety factor since the CG
				 // method will in general not be converged
  const double beta = 1.2 * max_eigenvalue;
  const double alpha = (data.smoothing_range > 0 ?
			max_eigenvalue / data.smoothing_range :
			max_eigenvalue / min_eigenvalue);
  delta = (beta-alpha)*0.5;
  theta = (beta+alpha)*0.5;
  is_initialized = true;
}



template <class MATRIX, class VECTOR>
inline
void
PreconditionChebyshev<MATRIX,VECTOR>::vmult (VECTOR &dst,
					     const VECTOR &src) const
{
  Assert (is_initialized, ExcMessage("Preconditioner not initialized"));
  double rhok  = delta / theta,  sigma = theta / delta;
  if (data.nonzero_starting && !dst.all_zero())
    {
      matrix_ptr->vmult (update1, dst);
      update1 -= src;
      update1 /= theta;
      update1.scale (data.matrix_diagonal_inverse);
      dst -= update1;
    }
  else
    {
      dst.equ (1./theta, src);
      dst.scale (data.matrix_diagonal_inverse);
      update1.equ(-1.,dst);
    }

  for (unsigned int k=0; k<data.degree; ++k)
    {
      matrix_ptr->vmult (update2, dst);
      update2 -= src;
      update2.scale (data.matrix_diagonal_inverse);
      const double rhokp = 1./(2.*sigma-rhok);
      const double factor1 = rhokp * rhok, factor2 = 2.*rhokp/delta;
      rhok = rhokp;
      update1.sadd (factor1, factor2, update2);
      dst -= update1;
    }
}



template <class MATRIX, class VECTOR>
inline
void
PreconditionChebyshev<MATRIX,VECTOR>::Tvmult (VECTOR &dst,
					      const VECTOR &src) const
{
  Assert (is_initialized, ExcMessage("Preconditioner not initialized"));
  double rhok  = delta / theta,  sigma = theta / delta;
  if (data.nonzero_starting && !dst.all_zero())
    {
      matrix_ptr->Tvmult (update1, dst);
      update1 -= src;
      update1 /= theta;
      update1.scale (data.matrix_diagonal_inverse);
      dst -= update1;
    }
  else
    {
      dst.equ (1./theta, src);
      dst.scale (data.matrix_diagonal_inverse);
      update1.equ(-1.,dst);
    }

  for (unsigned int k=0; k<data.degree-1; ++k)
    {
      matrix_ptr->Tvmult (update2, dst);
      update2 -= src;
      update2.scale (data.matrix_diagonal_inverse);
      const double rhokp = 1./(2.*sigma-rhok);
      const double factor1 = rhokp * rhok, factor2 = 2.*rhokp/delta;
      rhok = rhokp;
      update1.sadd (factor1, factor2, update2);
      dst -= update1;
    }
}



template <class MATRIX, class VECTOR>
inline
void PreconditionChebyshev<MATRIX,VECTOR>::clear ()
{
  is_initialized = false;
  matrix_ptr = 0;
  data.matrix_diagonal_inverse.reinit(0);
  update1.reinit(0);
  update2.reinit(0);
}



#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

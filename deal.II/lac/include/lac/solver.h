/*----------------------------   solver.h     ---------------------------*/
/*      $Id$                 */
#ifndef __solver_H
#define __solver_H
/*----------------------------   solver.h     ---------------------------*/



// forward declaration
class SolverControl;
template <class Vector> class VectorMemory;



/**
 * Base class for iterative solvers.
 *
 * HAS TO BE UPDATED!
 *
 * This class defines possible
 * return states of linear solvers and provides interfaces to a memory
 * pool and the control object.
 *
 *
 * \subsection{Requirements for template classes}
 *
 * The class is templated to allow for different matrix and vector
 * classes, since iterative solvers do not rely on any special structure
 * of matrices or the format of storage. However, there are some common
 * requirements a matrix or vector type must fulfil to qualify as an
 * applicable type for the solvers in this hierarchy. These requirements
 * are listed following. The classes do not declare any concrete
 * class, they are rather intended to form a `signature' which a concrete
 * class has to conform to.
 *
 * \begin{verbatim}
 * class Matrix
 * {
 *   public:
 *                        // Application to a Vector
 *     void vmult (Vector& dst, const Vector& src) const;
 *
 *                        // Application of a preconditioner to
 *                        // a Vector, i.e. $dst=\tilde A^(-1) src$,
 *                        // where $\tilde A^(-1)$ is an approximation
 *                        // to the inverse if the matrix stored in
 *                        // this object.
 *     void precondition (Vector& dst, const Vector& src) const;
 * 
 *                        // Application of transpose to a Vector.
 *                        // Only used by special iterative methods.
 *     void T_vmult (Vector& dst, const Vector& src) const;
 *
 *                        // Application of a transposed preconditioner
 *                        // to a Vector. Only used by special
 *                        // iterative methods
 *    
 *     void T_precondition (Vector& dst, const Vector& src) const;
 * };
 *
 *
 * class Vector
 * {
 *   public:
 *                        // scalar product
 *     double operator * (const Vector& v) const;
 *
 *                        // addition of vectors
 *                        // $y = y + x$.
 *     void add (const Vector& x);
 *                        // $y = y + ax$.
 *     void add (double a, const Vector& x);
 *
 *                        // scaled addition of vectors
 *                        // $y = ay + x$.
 *     void sadd (double a,
 *                const Vector& x);
 *                        // $y = ay + bx$.
 *     void sadd (double a,
 *                double b, const Vector& x);
 *                        // $y = ay + bx + cz$.
 *     void sadd (double a,
 *                double b, const Vector& x,
 *                double c, const Vector& z);
 * 
 *                        // $y = ax$.
 *     void equ (double a, const Vector& x);
 *                        // $y = ax + bz$.
 *     void equ (double a, const Vector& x,
 *               double b, const Vector& z);
 * };
 * \end{verbatim}
 *
 * \subsection{AdditionalData}
 * Several solvers need additional data, like the damping parameter #omega# of the
 * #SolverRichardson# class or the maximum number of tmp vectors of the #SolverGMRES#.
 * To have a standardized constructor for each solver class the #struct AdditionalData#
 * has been introduced to each solver class. Some solver needs no additional data, like
 * the SolverCG or SolverBicgstab. For these solvers the struct #AdditionalData# is empty
 * and the calling
 * of the constructor has not change as a default #AdditionalData# is set by default. 
 *
 * Now the generating of a solver looks like
 * \begin{verbatim}
 *                               // GMRES with 50 tmp vectors
 * SolverGMRES solver_gmres (solver_control, vector_memory,
 *                           SolverGMRES::AdditionalData(50));
 *
 *                               // Richardson with omega=0.8
 * SolverRichardson solver_richardson (solver_control, vector_memory,
 *                                     SolverGMRES::AdditionalData(0.8));
 *
 *                               // CG with default AdditionalData
 * SolverCG solver_cg (solver_control, vector_memory);
 * \end{verbatim}
 */
template <class Matrix, class Vector>
class Solver
{
  public:
				     /**
				      * Declare possible return values of a
				      * solver object.
				      */
    enum ReturnState {
	  success=0, exceeded, breakdown
    };
    
    
    
				     /**
				      * Constructor. Assign a control
				      * object which stores the required
				      * precision and an object to provide
				      * memory.
				      *
				      * Of both objects, a reference is
				      * stored, so it is the user's
				      * responsibility to guarantee that the
				      * lifetime of the two arguments is at
				      * least as long as that of the solver
				      * object.
				      */
    Solver (SolverControl &, VectorMemory<Vector> &);

				     /**
				      * Solver procedure.  This is in
				      * fact only a template for the
				      * solve-function to be
				      * implemented in derived classes
				      */
    template<class Preconditioner>
    ReturnState solve (const Matrix &A,
		       Vector &x,
		       const Vector &b,
		       const Preconditioner& precondition) const;
    
				     /**
				      * Access to control object.
				      */
    SolverControl& control() const;
    
  protected:
				     /**
				      * Calculation of convergence
				      * criterion. To allow further
				      * flexibility in solvers, the
				      * convergence criterion can be
				      * implemented by the user. Each
				      * method has it's standard
				      * criterion (usually some kind of
				      * residual) already implemented.
				      */
    virtual double criterion() = 0;
    
				     /**
				      * Additional entry point for
				      * examination.  Allows access to
				      * internal variables of a solver
				      * at a given point dependent of
				      * the actual solver method. Does
				      * nothing by default and should
				      * only be implemented to test
				      * iterative methods.
				      */
    virtual void expert() {};

				     /**
				      * Control structure.
				      */
    SolverControl &cntrl;
    
				     /**
				      * Memory for auxilliary vectors.
				      */
    VectorMemory<Vector> &memory;
};





/**
 * Base class for non-symmetric linear solvers.  A second entry point
 * allows the simultaneous computation of a dual problem.
 */
template<class Matrix, class Vector>
class SolverDual : public Solver<Matrix, Vector>
{
  public:
				     /**
				      * Constructor.
				      */
    SolverDual(SolverControl&, VectorMemory<Vector>&);
    
				     /**
				      * Solve the original problem
				      * $Ax=b$.
				      */
    virtual ReturnState solve (const Matrix &A,
			       Vector       &x,
			       const Vector &b) = 0;

				     /**
				      * Solve the two problems
				      * $Ax=b1$ and $A^Tz=b2$ simultanously.
				      */
    virtual ReturnState solve (const Matrix &A,
			       Vector       &x,
			       const Vector &b1,
			       Vector       &z,
			       const Vector &b2) = 0;
};






/*-------------------------------- Inline functions ------------------------*/

template <class Matrix, class Vector>
inline
SolverControl & Solver<Matrix,Vector>::control() const {
  return cntrl;
};



template<class Matrix, class Vector>
inline
Solver<Matrix, Vector>::Solver(SolverControl &cn, VectorMemory<Vector> &mem)
		: cntrl(cn),
		  memory(mem)
{};





/*----------------------------   solver.h     ---------------------------*/
/* end of #ifndef __solver_H */
#endif
/*----------------------------   solver.h     ---------------------------*/

//----------------------------  solver.h  ---------------------------
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
//----------------------------  solver.h  ---------------------------
#ifndef __deal2__solver_h
#define __deal2__solver_h

template<typename number> class Vector;
template<class VECTOR>  class VectorMemory;
class SolverControl;


/**
 * Base class for iterative solvers.
 *
//TODO:[?] * HAS TO BE UPDATED!
 *
 * This class defines possible
 * return states of linear solvers and provides interfaces to a memory
 * pool and the control object.
 *
 *
 * @sect3{Requirements for template classes}
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
 * @begin{verbatim}
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
 * @end{verbatim}
 *
 *
 * @sect3{AdditionalData}
 *
 * Several solvers need additional data, like the damping parameter @p{omega} of the
 * @p{SolverRichardson} class or the maximum number of temporary vectors of the @p{SolverGMRES}.
 * To have a standardized constructor for each solver class the @p{struct AdditionalData}
 * has been introduced to each solver class. Some solvers need no additional data, like
 * @p{SolverCG} or @p{SolverBicgstab}. For these solvers the struct @p{AdditionalData} is empty
 * and calling the constructor may be done without giving the additional structure as
 * an argument as a default @p{AdditionalData} is set by default. 
 *
 * Now the generating of a solver looks like
 * @begin{verbatim}
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
 * @end{verbatim}
 *
 * Using a unified constructor parameter list for all solvers was introduced when the
 * @p{SolverSelector} class was written; the unified interface enabled us to use this
 * class unchanged even if the number of types of parameters to a certain solver
 * changes and it is still possible in a simple way to give these additional data to
 * the @p{SolverSelector} object for each solver which it may use.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, Ralf Hartmann, 1997-2001
 */
template <class Vector = ::Vector<double> >
class Solver
{
  public:
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
				      * Access to control object.
				      */
    SolverControl& control() const;
  
  protected:

				     /**
				      * Control structure.
				      */
    SolverControl &cntrl;
    
				     /**
				      * Memory for auxilliary vectors.
				      */
    VectorMemory<Vector> &memory;
};


/*-------------------------------- Inline functions ------------------------*/

template <class Vector>
inline
SolverControl &
Solver<Vector>::control() const
{
  return cntrl;
};


template<class Vector>
inline
Solver<Vector>::Solver(SolverControl &cn, VectorMemory<Vector> &mem)
		: cntrl(cn),
		  memory(mem)
{};


#endif

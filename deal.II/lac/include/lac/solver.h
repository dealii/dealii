//----------------------------  solver.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  solver.h  ---------------------------
#ifndef __deal2__solver_h
#define __deal2__solver_h

#include <base/config.h>
#include <base/subscriptor.h>

template <typename number> class Vector;
template <class VECTOR>    class VectorMemory;
class SolverControl;

/*!@addtogroup Solvers */
/*@{*/

//! Base class for iterative solvers.
/**
 * This class defines possible return states of linear solvers and
 * provides interfaces to a memory pool and the control object.
 *
 * @sect3{Requirements for template classes}
 *
 * Since iterative solvers do not rely on any special structure of
 * matrices or the format of storage, but only require that matrices
 * and vector define certain operations such as matrix-vector
 * products, or scalar products between vectors, this class as well as
 * the derived classes implementing concrete linear solvers are
 * templated on the types of matrices and vectors. However, there are
 * some common requirements a matrix or vector type must fulfill to
 * qualify as an applicable type for the solvers in this
 * hierarchy. These requirements are listed following. The listed
 * classes are not any concrete class, they are rather intended to
 * form a `signature' which a concrete class has to conform to. Note
 * that the matrix and vector classes within this library of course
 * conform to this interface.
 *
 * @begin{verbatim}
 * class Matrix
 * {
 *   public:
 *                        // Application of matrix to vector src.
 *                        // write result into dst
 *     void vmult (VECTOR &dst, const VECTOR &src) const;
 * 
 *                        // Application of transpose to a Vector.
 *                        // Only used by certain iterative methods.
 *     void Tvmult (VECTOR &dst, const VECTOR &src) const;
 * };
 *
 *
 * class VECTOR
 * {
 *   public:
 *                        // resize to have the same structure
 *                        // as the one provided and/or
 *                        // clear vector. note
 *                        // that the second argument must have
 *                        // a default value equal to false
 *     void reinit (const VECTOR&,
 *                  bool  leave_elements_uninitialized = false);
 *
 *                        // scalar product
 *     double operator * (const VECTOR &v) const;
 *
 *                        // addition of vectors
 *                        // $y = y + x$.
 *     void add (const VECTOR &x);
 *
 *                        // $y = y + ax$.
 *     void add (const double  a,
 *               const VECTOR &x);
 *
 *                        // $y = ay + bx$.
 *     void sadd (const double  a,
 *                const double  b,
 *                const VECTOR &x);
 * 
 *                        // $y = ax$.
 *     void equ (const double  a,
 *               const VECTOR &x);
 *
 *                        // scale the elements of the vector
 *                        // by a fixed value
 *     void scale (const double a);
 *
 *                        // return the l2 norm of the vector
 *     double l2_norm () const;
 * };
 * @end{verbatim}
 *
 * In addition, for some solvers there has to be a global function
 * @p{swap(VECTOR &a, VECTOR &b)} that exchanges the values of the two vectors.
 *
 * The preconditioners used must have the same interface as matrices,
 * i.e. in particular they have to provide a member function @p{vmult}
 * which denotes the application of the preconditioner.
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
template <class VECTOR = Vector<double> >
class Solver : public Subscriptor
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
    Solver (SolverControl &, VectorMemory<VECTOR> &);

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
    VectorMemory<VECTOR> &memory;
};

/*@}*/
/*-------------------------------- Inline functions ------------------------*/

template <class VECTOR>
inline
SolverControl &
Solver<VECTOR>::control() const
{
  return cntrl;
}


template<class VECTOR>
inline
Solver<VECTOR>::Solver(SolverControl &cn, VectorMemory<VECTOR> &mem)
		: cntrl(cn),
		  memory(mem)
{}


#endif

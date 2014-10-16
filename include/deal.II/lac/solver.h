// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__solver_h
#define __deal2__solver_h

#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/lac/vector_memory.h>

DEAL_II_NAMESPACE_OPEN

template <typename number> class Vector;
class SolverControl;

/**
 * A base class for iterative linear solvers. This class
 * provides interfaces to a memory pool and the objects that
 * determine whether a solver has converged.
 *
 *
 * <h3>Requirements common to derived solver classes</h3>
 *
 * Since iterative solvers do not rely on any special structure of
 * matrices or the format of storage but only require that matrices
 * and vectors define certain operations such as matrix-vector
 * products, or scalar products between vectors, this class as well as
 * the derived classes and their member functions implementing concrete
 * linear solvers are templated on the types of matrices and vectors.
 * However, there are some common requirements a matrix or vector type
 * must fulfill to qualify as an acceptable type for the solvers in this
 * hierarchy. These requirements are listed below.
 *
 * The classes we show below are not any concrete class. Rather, they are intended to
 * form a `signature' which a concrete class has to conform to. Note
 * that the matrix and vector classes within this library of course
 * conform to this interface; therefore, SparseMatrix and Vector are
 * good examples for these classes as they provide the necessary
 * signatures of member functions.
 *
 * @code
 * class Matrix
 * {
 *   public:
 *                        // Application of matrix to vector src.
 *                        // Write result into dst
 *     void vmult (VECTOR &dst,
 *                 const VECTOR &src) const;
 *
 *                        // Application of transpose to a vector.
 *                        // Only used by some iterative methods.
 *     void Tvmult (VECTOR &dst,
 *                  const VECTOR &src) const;
 * };
 *
 *
 * class Vector
 * {
 *   public:
 *                        // Resize the current object to have
 *                        // the same size and layout as the model_vector
 *                        // argument provided. The second argument
 *                        // indicates whether to clear the current
 *                        // object after resizing.
 *                        // The second argument must have
 *                        // a default value equal to false
 *     void reinit (const Vector &model_vector,
 *                  const bool  leave_elements_uninitialized = false);
 *
 *                        // Scalar product between the current object
 *                        // and the argument
 *     double operator * (const Vector &v) const;
 *
 *                        // Addition of vectors
 *     void add (const Vector &x);
 *
 *                        // Scaled addition of vectors
 *     void add (const double  a,
 *               const Vector &x);
 *
 *                        // Scaled addition of vectors
 *     void sadd (const double  a,
 *                const double  b,
 *                const Vector &x);
 *
 *                        // Scaled assignment of a vector
 *     void equ (const double  a,
 *               const Vector &x);
 *
 *                        // Multiply the elements of the current
 *                        // object by a fixed value
 *     Vector & operator *= (const double a);
 *
 *                        // Return the l2 norm of the vector
 *     double l2_norm () const;
 * };
 * @endcode
 *
 * In addition, for some solvers there has to be a global function
 * <tt>swap(VECTOR &a, VECTOR &b)</tt> that exchanges the values of the two vectors.
 *
 * The preconditioners used must have the same interface as matrices,
 * i.e. in particular they have to provide a member function @p vmult
 * which denotes the application of the preconditioner.
 *
 *
 * <h3>AdditionalData</h3>
 *
 * Several solvers need additional data, like the damping parameter @p omega
 * of the @p SolverRichardson class or the maximum number of temporary
 * vectors of @p SolverGMRES.  To have a standardized way of constructing solvers,
 * each solver class has a <tt>struct AdditionalData</tt> as a member, and constructors
 * of all solver classes take such an argument. Some solvers need no additional data,
 * or may not at the current time. For these solvers the struct @p AdditionalData is
 * empty and calling the constructor may be done without giving the additional
 * structure as an argument as a default @p AdditionalData is set by default.
 *
 * With this, creating a solver looks like
 * @code
 *                               // GMRES with restart every 50 iterations
 * SolverGMRES solver_gmres (solver_control, vector_memory,
 *                           SolverGMRES::AdditionalData(50));
 *
 *                               // Richardson with omega=0.8
 * SolverRichardson solver_richardson (solver_control, vector_memory,
 *                                     SolverGMRES::AdditionalData(0.8));
 *
 *                               // CG with default AdditionalData
 * SolverCG solver_cg (solver_control, vector_memory);
 * @endcode
 *
 * Using a unified constructor parameter list for all solvers supports
 * the @p SolverSelector class; the unified interface
 * enables us to use this class unchanged even if the number of types of
 * parameters to a certain solver changes and it is still possible in a simple
 * way to give these additional data to the @p SolverSelector object for each
 * solver which it may use.
 *
 * @ingroup Solvers
 * @author Wolfgang Bangerth, Guido Kanschat, Ralf Hartmann, 1997-2001, 2005
 */
template <class VECTOR = Vector<double> >
class Solver : public Subscriptor
{
public:
  /**
   * Constructor. Takes a control
   * object which evaluates the
   * conditions for convergence,
   * and an object that allows solvers to allocate
   * memory for temporary objects.
   *
   * Of both objects, a reference is
   * stored, so it is the user's
   * responsibility to guarantee that the
   * lifetime of the two arguments is at
   * least as long as that of the solver
   * object.
   */
  Solver (SolverControl        &solver_control,
          VectorMemory<VECTOR> &vector_memory);

  /**
   * Constructor. Takes a control
   * object which evaluates the
   * conditions for convergence. In
   * contrast to the other
   * constructor, this constructor
   * designates an internal object of
   * type GrowingVectorMemory to
   * allocate memory.
   *
   * A reference to the control
   * object is stored, so it is the
   * user's responsibility to
   * guarantee that the lifetime of
   * the argument is at least
   * as long as that of the solver
   * object.
   */
  Solver (SolverControl        &solver_control);

  /**
   * Return a reference to the object that controls
   * convergence.
   */
  SolverControl &control() const;

protected:
  /**
   * A static vector memory object
   * to be used whenever no such
   * object has been given to the
   * constructor.
   */
  mutable GrowingVectorMemory<VECTOR> static_vector_memory;

  /**
   * Reference to the object that determines convergence.
   */
  SolverControl &cntrl;

  /**
   * A reference to an object that provides memory for auxiliary vectors.
   */
  VectorMemory<VECTOR> &memory;
};


/*-------------------------------- Inline functions ------------------------*/

template<class VECTOR>
inline
Solver<VECTOR>::Solver (SolverControl        &solver_control,
                        VectorMemory<VECTOR> &vector_memory)
  :
  cntrl(solver_control),
  memory(vector_memory)
{}



template<class VECTOR>
inline
Solver<VECTOR>::Solver (SolverControl        &solver_control)
  :
  cntrl(solver_control),
  memory(static_vector_memory)
{}



template <class VECTOR>
inline
SolverControl &
Solver<VECTOR>::control() const
{
  return cntrl;
}




DEAL_II_NAMESPACE_CLOSE

#endif

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
#include <deal.II/lac/solver_control.h>

#include <boost/signals2.hpp>

DEAL_II_NAMESPACE_OPEN

template <typename number> class Vector;

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

protected:
  /**
   * A static vector memory object
   * to be used whenever no such
   * object has been given to the
   * constructor.
   */
  mutable GrowingVectorMemory<VECTOR> static_vector_memory;

  /**
   * A reference to an object that provides memory for auxiliary vectors.
   */
  VectorMemory<VECTOR> &memory;

private:
  /**
   * A class whose operator() combines two states indicating whether
   * we should continue iterating or stop, and returns a state that
   * dominates. The rules are:
   * - If one of the two states indicates failure, return failure.
   * - Otherwise, if one of the two states indicates to continue
   *   iterating, then continue iterating.
   * - Otherwise, return success.
   */
  struct StateCombiner
  {
    typedef SolverControl::State result_type;

    SolverControl::State operator() (const SolverControl::State state1,
                                     const SolverControl::State state2) const;

    template <typename Iterator>
    SolverControl::State operator() (const Iterator begin,
                                     const Iterator end) const;
  };

protected:
  /**
   * A signal that iterative solvers can execute at the end of every
   * iteration (or in an otherwise periodic fashion) to find out whether
   * we should continue iterating or not. The signal may call one or
   * more slots that each will make this determination by themselves,
   * and the result over all slots (function calls) will be determined
   * by the StateCombiner object.
   *
   * The arguments passed to the signal are (i) the number of the current
   * iteration; (ii) the value that is used to determine convergence (oftentimes
   * the residual, but in other cases other quantities may be used as long
   * as they converge to zero as the iterate approaches the solution of
   * the linear system); and (iii) a vector that corresponds to the current
   * best guess for the solution at the point where the signal is called.
   * Note that some solvers do not update the approximate solution in every
   * iteration but only after convergence or failure has been determined
   * (GMRES is an example); in such cases, the vector passed as the last
   * argument to the signal is simply the best approximate at the time
   * the signal is called, but not the vector that will be returned if
   * the signal's return value indicates that the iteration should be
   * terminated.
   */
  boost::signals2::signal<SolverControl::State (const unsigned int iteration,
                                                const double        check_value,
                                                const VECTOR       &current_iterate),
                                                      StateCombiner> iteration_status;
};


/*-------------------------------- Inline functions ------------------------*/


template <class VECTOR>
inline
SolverControl::State
Solver<VECTOR>::StateCombiner::operator ()(const SolverControl::State state1,
                                           const SolverControl::State state2) const
{
  if ((state1 == SolverControl::failure)
      ||
      (state2 == SolverControl::failure))
    return SolverControl::failure;
  else if ((state1 == SolverControl::iterate)
           ||
           (state2 == SolverControl::iterate))
    return SolverControl::iterate;
  else
    return SolverControl::success;
}


template <class VECTOR>
template <typename Iterator>
inline
SolverControl::State
Solver<VECTOR>::StateCombiner::operator ()(const Iterator begin,
                                           const Iterator end) const
{
  Assert (begin != end, ExcMessage ("You can't combine iterator states if no state is given."));

  // combine the first with all of the following states
  SolverControl::State state = *begin;
  Iterator p = begin;
  ++p;
  for (; p != end; ++p)
    state = this->operator()(state, *p);

  return state;
}


template<class VECTOR>
inline
Solver<VECTOR>::Solver (SolverControl        &solver_control,
                        VectorMemory<VECTOR> &vector_memory)
  :
  memory(vector_memory)
{
  // connect the solver control object to the signal. SolverControl::check
  // only takes two arguments, the iteration and the check_value, and so
  // we simply ignore the third argument that is passed in whenever the
  // signal is executed
  iteration_status.connect (std_cxx11::bind(&SolverControl::check,
                                            std_cxx11::ref(solver_control),
                                            std_cxx11::_1,
                                            std_cxx11::_2));
}



template<class VECTOR>
inline
Solver<VECTOR>::Solver (SolverControl        &solver_control)
  :
  // use the static memory object this class owns
  memory(static_vector_memory)
{
  // connect the solver control object to the signal. SolverControl::check
  // only takes two arguments, the iteration and the check_value, and so
  // we simply ignore the third argument that is passed in whenever the
  // signal is executed
  iteration_status.connect (std_cxx11::bind(&SolverControl::check,
                                            std_cxx11::ref(solver_control),
                                            std_cxx11::_1,
                                            std_cxx11::_2));
}




DEAL_II_NAMESPACE_CLOSE

#endif

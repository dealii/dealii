// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2015 by the deal.II authors
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

#ifndef dealii__solver_h
#define dealii__solver_h

#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/solver_control.h>

// Ignore deprecation warnings for auto_ptr.
DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/signals2.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

DEAL_II_NAMESPACE_OPEN

template <typename number> class Vector;

/**
 * A base class for iterative linear solvers. This class provides interfaces
 * to a memory pool and the objects that determine whether a solver has
 * converged.
 *
 *
 * <h3>Requirements common to derived solver classes</h3>
 *
 * Since iterative solvers do not rely on any special structure of matrices or
 * the format of storage but only require that matrices and vectors define
 * certain operations such as matrix-vector products, or scalar products
 * between vectors, this class as well as the derived classes and their member
 * functions implementing concrete linear solvers are templated on the types
 * of matrices and vectors. However, there are some common requirements a
 * matrix or vector type must fulfill to qualify as an acceptable type for the
 * solvers in this hierarchy. These requirements are listed below.
 *
 * The classes we show below are not any concrete class. Rather, they are
 * intended to form a `signature' which a concrete class has to conform to.
 * Note that the matrix and vector classes within this library of course
 * conform to this interface; therefore, SparseMatrix and Vector are good
 * examples for these classes as they provide the necessary signatures of
 * member functions.
 *
 * @code
 * class Matrix
 * {
 *   public:
 *                        // Application of matrix to vector src.
 *                        // Write result into dst
 *     void vmult (VectorType       &dst,
 *                 const VectorType &src) const;
 *
 *                        // Application of transpose to a vector.
 *                        // Only used by some iterative methods.
 *     void Tvmult (VectorType       &dst,
 *                  const VectorType &src) const;
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
 *                        // Inner product between the current object
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
 *                        // Combined scaled addition of vector x into
 *                        // the current object and subsequent inner
 *                        // product of the current object with v
 *     double add_and_dot (const double  a,
 *                         const Vector &x,
 *                         const Vector &v);
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
 * <tt>swap(VectorType &a, VectorType &b)</tt> that exchanges the values of
 * the two vectors.
 *
 * Finally, the solvers also expect an instantiation of
 * GrowingVectorMemory@<VectorType@>. These instantiations are provided by the
 * deal.II library for the built-in vector types, but must be explicitly added
 * for user-provided vector classes. Otherwise, the linker will complain that
 * it cannot find the constructors and destructors of GrowingVectorMemory that
 * happen in the @p Solver class.
 *
 * @code
 * // Definition and implementation of vector class
 * class UserVector { ... };
 *
 * // Create explicit instantiation for the vector class. If your project
 * // consists of multiple files, including header files, this instantiation
 * // must be put in a <code>.cc</code> file in order to instantiate only
 * // once.
 * #include <deal.II/lac/vector_memory.templates.h>
 *
 * template class VectorMemory<UserVector>;
 * template class GrowingVectorMemory<UserVector>;
 * @endcode
 *
 * The preconditioners used must have the same interface as matrices, i.e. in
 * particular they have to provide a member function @p vmult which denotes
 * the application of the preconditioner.
 *
 *
 * <h3>AdditionalData</h3>
 *
 * Several solvers need additional data, like the damping parameter @p omega
 * of the @p SolverRichardson class or the maximum number of temporary vectors
 * of @p SolverGMRES.  To have a standardized way of constructing solvers,
 * each solver class has a <tt>struct AdditionalData</tt> as a member, and
 * constructors of all solver classes take such an argument. Some solvers need
 * no additional data, or may not at the current time. For these solvers the
 * struct @p AdditionalData is empty and calling the constructor may be done
 * without giving the additional structure as an argument as a default @p
 * AdditionalData is set by default.
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
 * Using a unified constructor parameter list for all solvers supports the @p
 * SolverSelector class; the unified interface enables us to use this class
 * unchanged even if the number of types of parameters to a certain solver
 * changes and it is still possible in a simple way to give these additional
 * data to the @p SolverSelector object for each solver which it may use.
 *
 *
 * <h3>Observing the progress of linear solver iterations</h3>
 *
 * The Solver class, being the base class for all of the iterative solvers
 * such as SolverCG, SolverGMRES, etc, provides the facilities by which actual
 * solver implementations determine whether the iteration is converged, not
 * yet converged, or has failed. Typically, this is done using an object of
 * type SolverControl that is passed to the solver classes's constructors and
 * from them down to the constructor of this base class. Every one of the
 * tutorial programs that solves a linear problem (starting with step-3) uses
 * this method and it is described in detail there. However, the underlying
 * mechanism is more general and allows for many other uses to observe how the
 * linear solver iterations progress.
 *
 * The basic approach is that the iterative solvers invoke a <i>signal</i> at
 * the end of each iteration to determine whether the solution is converged. A
 * signal is a class that has, conceptually, a list of pointers to functions
 * and every time the signal is invoked, each of these functions are called.
 * In the language of signals, the functions called are called <i>slots</i>
 * and one can attach any number of slots to a signal. (The implementation of
 * signals and slots we use here is the one from the BOOST.signals2 library.)
 * A number of details may clarify what is happening underneath: - In reality,
 * the signal object does not store pointers to functions, but function
 * objects as slots. Each slot must conform to a particular signature: here,
 * it is an object that can be called with three arguments (the number of the
 * current linear iteration, the current residual, and the current iterate;
 * more specifics are discussed in the documentation of the connect()
 * function). A pointer to a function with this argument list satisfies the
 * requirements, but you can also pass a member function whose
 * <code>this</code> argument has been bound using the
 * <code>std_cxx11::bind</code> mechanism (see the example below). - Each of
 * the slots will return a value that indicates whether the iteration should
 * continue, should stop because it has succeeded, or stop because it has
 * failed. The return type of slots is therefore of type SolverControl::State.
 * The returned values from all of the slots will then have to be combined
 * before they are returned to the iterative solver that invoked the signal.
 * The way this works is that if at least one slot returned
 * SolverControl::failure, then the combined value is SolverControl::failure;
 * otherwise, if at least one slot returned SolverControl::iterate, then this
 * is going to be the return value of the signal; finally, only if all slots
 * return SolverControl::success will the signal's return value be
 * SolverControl::success. - It may of course be that a particular slot has
 * been connected to the signal only to observe how the solution or a specific
 * part of it converges, but has no particular opinion on whether the
 * iteration should continue or not. In such cases, the slot should just
 * return SolverControl::success, which is the weakest of all return values
 * according to the rules laid out above.
 *
 * Given all this, it should now be obvious how the SolverControl object fits
 * into this scheme: when a SolverControl object is passed to the constructor
 * of the current class, we simply connect the SolverControl::check() function
 * of that object as a slot to the signal we maintain here. In other words,
 * since a Solver object is always constructed using a SolverControl object,
 * there is always at least one slot associated with the signal, namely the
 * one that determines convergence.
 *
 * On the other hand, using the connect() member function, it is possible to
 * connect any number of other slots to the signal to observe whatever it is
 * you want to observe. The connect() function also returns an object that
 * describes the connection from the signal to the slot, and the corresponding
 * BOOST functions then allow you to disconnect the slot if you want.
 *
 * An example may illuminate these issues. In the step-3 tutorial program, let
 * us add a member function as follows to the main class:
 * @code
 *  SolverControl::State
 *  Step3::write_intermediate_solution (const unsigned int    iteration,
 *                                      const double          , //check_value
 *                                      const Vector<double> &current_iterate) const
 *    {
 *      DataOut<2> data_out;
 *      data_out.attach_dof_handler (dof_handler);
 *      data_out.add_data_vector (current_iterate, "solution");
 *      data_out.build_patches ();
 *
 *      std::ofstream output ((std::string("solution-")
 *                             + Utilities::int_to_string(iteration,4) + ".vtu").c_str());
 *      data_out.write_vtu (output);
 *
 *      return SolverControl::success;
 *    }
 * @endcode
 * The function satisfies the signature necessary to be a slot for the signal
 * discussed above, with the exception that it is a member function and
 * consequently requires a <code>this</code> pointer. What the function does
 * is to take the vector given as last argument and write it into a file in
 * VTU format with a file name derived from the number of the iteration.
 *
 * This function can then be hooked into the CG solver by modifying the
 * <code>Step3::solve()</code> function as follows:
 * @code
 * void Step3::solve ()
 * {
 *   SolverControl           solver_control (1000, 1e-12);
 *   SolverCG<>              solver (solver_control);
 *
 *   solver.connect (std_cxx11::bind (&Step3::write_intermediate_solution,
 *                                    this,
 *                                    std_cxx11::_1,
 *                                    std_cxx11::_2,
 *                                    std_cxx11::_3));
 *   solver.solve (system_matrix, solution, system_rhs,
 *                 PreconditionIdentity());
 * }
 * @endcode
 * The use of <code>std_cxx11::bind</code> here ensures that we convert the
 * member function with its three arguments plus the <code>this</code>
 * pointer, to a function that only takes three arguments, by fixing the
 * implicit <code>this</code> argument of the function to the
 * <code>this</code> pointer in the current function.
 *
 * It is well understood that the CG method is a smoothing iteration (in the
 * same way as the more commonly used Jacobi or SSOR iterations are
 * smoothers). The code above therefore allows to observe how the solution
 * becomes smoother and smoother in every iteration. This is best observed by
 * initializing the solution vector with randomly distributed numbers in
 * $[-1,1]$, using code such as
 * @code
 *   for (unsigned int i=0; i<solution.size(); ++i)
 *     solution(i) = 2.*rand()/RAND_MAX-1;
 * @endcode
 * Using this, the slot will then generate files that when visualized look
 * like this over the course of iterations zero to five: <table> <tr> <td>
 * @image html "cg-monitor-smoothing-0.png"
 * </td> <td>
 * @image html "cg-monitor-smoothing-1.png"
 * </td> <td>
 * @image html "cg-monitor-smoothing-2.png"
 * </td> </tr> <tr> <td>
 * @image html "cg-monitor-smoothing-3.png"
 * </td> <td>
 * @image html "cg-monitor-smoothing-4.png"
 * </td> <td>
 * @image html "cg-monitor-smoothing-5.png"
 * </td> </tr> </table>
 *
 * @ingroup Solvers
 * @author Wolfgang Bangerth, Guido Kanschat, Ralf Hartmann, 1997-2001, 2005,
 * 2014
 */
template <class VectorType = Vector<double> >
class Solver : public Subscriptor
{
public:
  /**
   * A typedef for the underlying vector type
   */
  typedef VectorType vector_type;

  /**
   * Constructor. Takes a control object which evaluates the conditions for
   * convergence, and an object that allows solvers to allocate memory for
   * temporary objects.
   *
   * Of both objects, a reference is stored, so it is the user's
   * responsibility to guarantee that the lifetime of the two arguments is at
   * least as long as that of the solver object.
   */
  Solver (SolverControl            &solver_control,
          VectorMemory<VectorType> &vector_memory);

  /**
   * Constructor. Takes a control object which evaluates the conditions for
   * convergence. In contrast to the other constructor, this constructor
   * designates an internal object of type GrowingVectorMemory to allocate
   * memory.
   *
   * A reference to the control object is stored, so it is the user's
   * responsibility to guarantee that the lifetime of the argument is at least
   * as long as that of the solver object.
   */
  Solver (SolverControl &solver_control);

  /**
   * Connect a function object that will be called periodically within
   * iterative solvers. This function is used to attach monitors to iterative
   * solvers, either to determine when convergence has happened, or simply to
   * observe the progress of an iteration. See the documentation of this class
   * for more information.
   *
   * @param slot A function object specified here will, with each call,
   * receive the number of the current iteration, the value that is used to
   * check for convergence (typically the residual of the current iterate with
   * respect to the linear system to be solved) and the currently best
   * available guess for the current iterate. Note that some solvers do not
   * update the approximate solution in every iteration but only after
   * convergence or failure has been determined (GMRES is an example); in such
   * cases, the vector passed as the last argument to the signal is simply the
   * best approximate at the time the signal is called, but not the vector
   * that will be returned if the signal's return value indicates that the
   * iteration should be terminated. The function object must return a
   * SolverControl::State value that indicates whether the iteration should
   * continue, has failed, or has succeeded. The results of all connected
   * functions will then be combined to determine what should happen with the
   * iteration.
   *
   * @return A connection object that represents the connection from the
   * signal to the function object. It can be used to disconnect the function
   * object again from the signal. See the documentation of the BOOST Signals2
   * library for more information on connection management.
   */
  boost::signals2::connection
  connect (const std_cxx11::function<SolverControl::State (const unsigned int iteration,
                                                           const double       check_value,
                                                           const VectorType   &current_iterate)> &slot);



protected:
  /**
   * A static vector memory object to be used whenever no such object has been
   * given to the constructor.
   */
  mutable GrowingVectorMemory<VectorType> static_vector_memory;

  /**
   * A reference to an object that provides memory for auxiliary vectors.
   */
  VectorMemory<VectorType> &memory;

private:
  /**
   * A class whose operator() combines two states indicating whether we should
   * continue iterating or stop, and returns a state that dominates. The rules
   * are:
   * - If one of the two states indicates failure, return failure.
   * - Otherwise, if one of the two states indicates to continue iterating, then
   * continue iterating.
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
   * A signal that iterative solvers can execute at the end of every iteration
   * (or in an otherwise periodic fashion) to find out whether we should
   * continue iterating or not. The signal may call one or more slots that
   * each will make this determination by themselves, and the result over all
   * slots (function calls) will be determined by the StateCombiner object.
   *
   * The arguments passed to the signal are (i) the number of the current
   * iteration; (ii) the value that is used to determine convergence
   * (oftentimes the residual, but in other cases other quantities may be used
   * as long as they converge to zero as the iterate approaches the solution
   * of the linear system); and (iii) a vector that corresponds to the current
   * best guess for the solution at the point where the signal is called. Note
   * that some solvers do not update the approximate solution in every
   * iteration but only after convergence or failure has been determined
   * (GMRES is an example); in such cases, the vector passed as the last
   * argument to the signal is simply the best approximate at the time the
   * signal is called, but not the vector that will be returned if the
   * signal's return value indicates that the iteration should be terminated.
   */
  boost::signals2::signal<SolverControl::State (const unsigned int iteration,
                                                const double       check_value,
                                                const VectorType   &current_iterate),
                                                      StateCombiner> iteration_status;
};


/*-------------------------------- Inline functions ------------------------*/


template <class VectorType>
inline
SolverControl::State
Solver<VectorType>::StateCombiner::operator ()(const SolverControl::State state1,
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


template <class VectorType>
template <typename Iterator>
inline
SolverControl::State
Solver<VectorType>::StateCombiner::operator ()(const Iterator begin,
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


template<class VectorType>
inline
Solver<VectorType>::Solver (SolverControl        &solver_control,
                            VectorMemory<VectorType> &vector_memory)
  :
  memory(vector_memory)
{
  // connect the solver control object to the signal. SolverControl::check
  // only takes two arguments, the iteration and the check_value, and so
  // we simply ignore the third argument that is passed in whenever the
  // signal is executed
  connect (std_cxx11::bind(&SolverControl::check,
                           std_cxx11::ref(solver_control),
                           std_cxx11::_1,
                           std_cxx11::_2));
}



template<class VectorType>
inline
Solver<VectorType>::Solver (SolverControl &solver_control)
  :
  // use the static memory object this class owns
  memory(static_vector_memory)
{
  // connect the solver control object to the signal. SolverControl::check
  // only takes two arguments, the iteration and the check_value, and so
  // we simply ignore the third argument that is passed in whenever the
  // signal is executed
  connect (std_cxx11::bind(&SolverControl::check,
                           std_cxx11::ref(solver_control),
                           std_cxx11::_1,
                           std_cxx11::_2));
}



template<class VectorType>
inline
boost::signals2::connection
Solver<VectorType>::
connect (const std_cxx11::function<SolverControl::State (const unsigned int iteration,
                                                         const double       check_value,
                                                         const VectorType   &current_iterate)> &slot)
{
  return iteration_status.connect (slot);
}



DEAL_II_NAMESPACE_CLOSE

#endif

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_scope_exit_h
#define dealii_scope_exit_h

#include <deal.II/base/config.h>

#include <functional>


DEAL_II_NAMESPACE_OPEN

/**
 * A class that stores a function object (typically a lambda function)
 * that is executed in the destructor of the current object. Such
 * a class is useful to execute an action *whenever* the function in
 * which the object is declared exits, whether that is by falling
 * off the end of the function, an explicit `return` statement, or
 * because an exception is thrown. In all of these cases, the
 * destructors of all local objects are run, and the destructor
 * of the ScopeExit object then executes the stored action.
 *
 * The typical use case is that we want to perform a clean-up
 * action. For example, let's say that we want to print to the
 * console that we're leaving a long-running and potentially
 * complicated function. This could be achieved as follows:
 * @code
 * void long_running_and_complex_function ()
 * {
 *   ScopeExit scope_exit_action ([]() { std::cout << "Leaving the long
 * running function!\n"; });
 *
 *   for (...)
 *   {
 *     // some 800 lines of complex code
 *   }
 * @endcode
 * The point here is that in order to have the message printed, we do not
 * need to go through the 800 lines of complex code and find all places
 * where either there is an explicit `return` statement, or where a
 * function that is called may throw an exception that we do not explicitly
 * catch: *Any* way the current function is exited will lead to the
 * message being printed.
 *
 * A perhaps more instructive example is the following (from some past
 * iteration of the implementation of SUNDIALS::KINSOL::solve()):
 * @code
 *    status = KINSol(kinsol_mem, solution, data.strategy, u_scale, f_scale);
 *
 *    if (pending_exception)
 *      {
 *        try
 *          {
 *            std::rethrow_exception(pending_exception);
 *          }
 *        catch (const RecoverableUserCallbackError &exc)
 *          {
 *            pending_exception = nullptr;
 *            if (status == 0)
 *              ; // just eat the exception
 *            else
 *              throw;
 *          }
 *        catch (...)
 *          {
 *            pending_exception = nullptr;
 *            throw;
 *          }
 *      }
 *    AssertKINSOL(status);
 *
 *    internal::copy(initial_guess_and_solution, solution);
 *    // Free the vectors which are no longer used.
 * #ifdef DEAL_II_WITH_MPI
 *    if (is_serial_vector<VectorType>::value == false)
 *      {
 *        N_VDestroy_Parallel(solution);
 *        N_VDestroy_Parallel(u_scale);
 *        N_VDestroy_Parallel(f_scale);
 *      }
 *    else
 * #endif
 *      {
 *        N_VDestroy_Serial(solution);
 *        N_VDestroy_Serial(u_scale);
 *        N_VDestroy_Serial(f_scale);
 *      }
 *    long nniters;
 *    status = KINGetNumNonlinSolvIters(kinsol_mem, &nniters);
 *    AssertKINSOL(status);
 *    if (J != nullptr)
 *      SUNMatDestroy(J);
 *    if (LS != nullptr)
 *      SUNLinSolFree(LS);
 *    KINFree(&kinsol_mem);
 *    return static_cast<unsigned int>(nniters);
 * @endcode
 *
 * This code has a bug: At the bottom, it contains clean-up code
 * that destroys a variety of objects created throughout the life of
 * the function we are looking at. Because SUNDIALS is a library
 * written in C, objects such as `solution`, `u_scale`, `J`, `LS`,
 * etc., are really just regular pointers that need to be explicitly
 * freed. But this does not happen in the places higher up where we
 * `throw` an exception and don't catch it again immediately. One
 * might be tempted to just copy the clean-up code to right before the
 * `std::rethrow_exception` statement, but there is a path through the
 * `catch` clauses marked with "just eat the exception" that ends up
 * returning control flow back to the regular path, and if we hit this
 * path, we would call the clean-up functions *twice* -- also not what
 * we want. The upshot of this kind of example is that there are cases
 * where it is quite difficult to say where exactly a function will
 * exit and where the clean-up code needs to be placed; it is far
 * easier to write it only once and let the compiler determine where
 * it needs to be called, via a class such as the current one.
 *
 * The function is conceptually very similar to the
 * [std::experimental::scope_exit](https://en.cppreference.com/w/cpp/experimental/scope_exit)
 * class template that may find its way into a future C++ standard.
 */
class ScopeExit
{
public:
  /**
   * Constructor. Takes a function object that is to be executed at the
   * place where the current object goes out of scope as argument.
   */
  explicit ScopeExit(const std::function<void()> &exit_function);

  /**
   * Copy constructor. These kinds of objects cannot be copied, so the
   * constructor is deleted.
   */
  ScopeExit(const ScopeExit &) = delete;

  /**
   * Destructor. Execute the stored action.
   */
  ~ScopeExit();

  /**
   * Copy operator. These kinds of objects cannot be copied, so the
   * operator is deleted.
   */
  ScopeExit &
  operator=(const ScopeExit &) = delete;

private:
  /**
   * A copy of the function to be executed.
   */
  const std::function<void()> exit_function;
};



inline ScopeExit::ScopeExit(const std::function<void()> &exit_function)
  : exit_function(exit_function)
{}



inline ScopeExit::~ScopeExit()
{
  // Actually trigger the stored function
  exit_function();
}


DEAL_II_NAMESPACE_CLOSE

#endif

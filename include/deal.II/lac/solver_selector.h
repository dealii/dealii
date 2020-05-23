// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_solver_selector_h
#define dealii_solver_selector_h


#include <deal.II/base/config.h>

#include <deal.II/base/smartpointer.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/vector.h>

DEAL_II_NAMESPACE_OPEN


/*!@addtogroup Solvers */
/*@{*/

/**
 * Selects a solver by changing a parameter.
 *
 * By calling the @p solve function of this @p SolverSelector, it selects the
 * @p solve function of that @p Solver that was specified in the constructor
 * of this class.
 *
 * <h3>Usage</h3> The simplest use of this class is the following:
 * @code
 * // generate a @p SolverSelector that calls the @p SolverCG
 * SolverControl control;
 * SolverSelector<Vector<double> > solver_selector ("cg", control);
 *
 * // generate e.g. a @p PreconditionRelaxation
 * PreconditionRelaxation<SparseMatrix<double>, Vector<double> >
 *   preconditioning(
 *     A, &SparseMatrix<double>::template precondition_SSOR<double>,0.8);
 *
 * // call the @p solve function with this preconditioning as last argument
 * solver_selector.solve(A,x,b,preconditioning);
 * @endcode
 * But the full usefulness of the @p SolverSelector class is not clear until
 * the presentation of the following example that assumes the user using the
 * @p ParameterHandler class and having declared a "solver" entry, e.g. with
 * @code
 * Parameter_Handler prm;
 * prm.declare_entry ("solver", "none",
 *                    Patterns::Selection(
 *                      SolverSelector<>::get_solver_names()));
 * ...
 * @endcode
 * Assuming that in the users parameter file there exists the line
 * @code
 * set solver = cg
 * @endcode
 * then the constructor call in the above example can be written as
 * @code
 * SolverSelector<SparseMatrix<double>, Vector<double> >
 *   solver_selector(prm.get("solver"), control);
 * @endcode
 *
 *
 * If at some time there exists a new solver "xyz" then the user does not need
 * to change their program. Only in the implementation of the @p SolverSelector
 * the calling of this solver has to be added and each user with program lines
 * quoted above only needs to 'set solver = xyz' in their parameter file to get
 * access to that new solver.
 *
 * @author Ralf Hartmann, 1999
 */
template <typename VectorType = Vector<double>>
class SolverSelector : public Subscriptor
{
public:
  /**
   * An alias for the underlying vector type
   */
  using vector_type = VectorType;

  /**
   * Constructor, filling in default values
   */
  SolverSelector() = default;

  /**
   * Constructor, selecting the solver @p name
   * and the SolverControl object @p control already.
   */
  SolverSelector(const std::string &name, SolverControl &control);

  /**
   * Destructor
   */
  virtual ~SolverSelector() override = default;

  /**
   * Solver procedure. Calls the @p solve function of the @p solver whose @p
   * SolverName was specified in the constructor.
   *
   */
  template <class Matrix, class Preconditioner>
  void
  solve(const Matrix &        A,
        VectorType &          x,
        const VectorType &    b,
        const Preconditioner &precond) const;

  /**
   * Select a new solver. Note that all solver names used in this class are
   * all lower case.
   */
  void
  select(const std::string &name);

  /**
   * Set a new SolverControl. This needs to be set before solving.
   */
  void
  set_control(SolverControl &ctrl);

  /**
   * Set the additional data. For more information see the @p Solver class.
   */
  void
  set_data(const typename SolverRichardson<VectorType>::AdditionalData &data);

  /**
   * Set the additional data. For more information see the @p Solver class.
   */
  void
  set_data(const typename SolverCG<VectorType>::AdditionalData &data);

  /**
   * Set the additional data. For more information see the @p Solver class.
   */
  void
  set_data(const typename SolverMinRes<VectorType>::AdditionalData &data);

  /**
   * Set the additional data. For more information see the @p Solver class.
   */
  void
  set_data(const typename SolverBicgstab<VectorType>::AdditionalData &data);

  /**
   * Set the additional data. For more information see the @p Solver class.
   */
  void
  set_data(const typename SolverGMRES<VectorType>::AdditionalData &data);

  /**
   * Set the additional data. For more information see the @p Solver class.
   */
  void
  set_data(const typename SolverFGMRES<VectorType>::AdditionalData &data);

  /**
   * Get the names of all implemented solvers. The list of possible
   * options includes:
   * <ul>
   * <li>  "richardson" </li>
   * <li>  "cg" </li>
   * <li>  "bicgstab" </li>
   * <li>  "gmres" </li>
   * <li>  "fgmres" </li>
   * <li>  "minres" </li>
   * </ul>
   */
  static std::string
  get_solver_names();

  /**
   * Exception.
   */
  DeclException1(ExcSolverDoesNotExist,
                 std::string,
                 << "Solver " << arg1 << " does not exist. Use one of "
                 << std::endl
                 << get_solver_names());



protected:
  /**
   * Stores the @p SolverControl that is needed in the constructor of each @p
   * Solver class. This can be changed with @p set_control().
   */
  SmartPointer<SolverControl, SolverSelector<VectorType>> control;

  /**
   * Stores the name of the solver.
   */
  std::string solver_name;

private:
  /**
   * Stores the additional data.
   */
  typename SolverRichardson<VectorType>::AdditionalData richardson_data;

  /**
   * Stores the additional data.
   */
  typename SolverCG<VectorType>::AdditionalData cg_data;

  /**
   * Stores the additional data.
   */
  typename SolverMinRes<VectorType>::AdditionalData minres_data;

  /**
   * Stores the additional data.
   */
  typename SolverBicgstab<VectorType>::AdditionalData bicgstab_data;

  /**
   * Stores the additional data.
   */
  typename SolverGMRES<VectorType>::AdditionalData gmres_data;

  /**
   * Stores the additional data.
   */
  typename SolverFGMRES<VectorType>::AdditionalData fgmres_data;
};

/*@}*/
/* --------------------- Inline and template functions ------------------- */


template <typename VectorType>
SolverSelector<VectorType>::SolverSelector(const std::string &name,
                                           SolverControl &    solver_control)
  : solver_name(name)
  , control(&solver_control)
{}



template <typename VectorType>
void
SolverSelector<VectorType>::select(const std::string &name)
{
  solver_name = name;
}



template <typename VectorType>
template <class Matrix, class Preconditioner>
void
SolverSelector<VectorType>::solve(const Matrix &        A,
                                  VectorType &          x,
                                  const VectorType &    b,
                                  const Preconditioner &precond) const
{
  if (solver_name == "richardson")
    {
      SolverRichardson<VectorType> solver(*control, richardson_data);
      solver.solve(A, x, b, precond);
    }
  else if (solver_name == "cg")
    {
      SolverCG<VectorType> solver(*control, cg_data);
      solver.solve(A, x, b, precond);
    }
  else if (solver_name == "minres")
    {
      SolverMinRes<VectorType> solver(*control, minres_data);
      solver.solve(A, x, b, precond);
    }
  else if (solver_name == "bicgstab")
    {
      SolverBicgstab<VectorType> solver(*control, bicgstab_data);
      solver.solve(A, x, b, precond);
    }
  else if (solver_name == "gmres")
    {
      SolverGMRES<VectorType> solver(*control, gmres_data);
      solver.solve(A, x, b, precond);
    }
  else if (solver_name == "fgmres")
    {
      SolverFGMRES<VectorType> solver(*control, fgmres_data);
      solver.solve(A, x, b, precond);
    }
  else
    Assert(false, ExcSolverDoesNotExist(solver_name));
}



template <typename VectorType>
void
SolverSelector<VectorType>::set_control(SolverControl &ctrl)
{
  control = &ctrl;
}



template <typename VectorType>
std::string
SolverSelector<VectorType>::get_solver_names()
{
  return "richardson|cg|bicgstab|gmres|fgmres|minres";
}



template <typename VectorType>
void
SolverSelector<VectorType>::set_data(
  const typename SolverGMRES<VectorType>::AdditionalData &data)
{
  gmres_data = data;
}



template <typename VectorType>
void
SolverSelector<VectorType>::set_data(
  const typename SolverFGMRES<VectorType>::AdditionalData &data)
{
  fgmres_data = data;
}



template <typename VectorType>
void
SolverSelector<VectorType>::set_data(
  const typename SolverRichardson<VectorType>::AdditionalData &data)
{
  richardson_data = data;
}



template <typename VectorType>
void
SolverSelector<VectorType>::set_data(
  const typename SolverCG<VectorType>::AdditionalData &data)
{
  cg_data = data;
}



template <typename VectorType>
void
SolverSelector<VectorType>::set_data(
  const typename SolverMinRes<VectorType>::AdditionalData &data)
{
  minres_data = data;
}



template <typename VectorType>
void
SolverSelector<VectorType>::set_data(
  const typename SolverBicgstab<VectorType>::AdditionalData &data)
{
  bicgstab_data = data;
}

DEAL_II_NAMESPACE_CLOSE

#endif

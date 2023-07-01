// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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

#ifndef DOXYGEN

#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_PETSC

#    include <deal.II/lac/petsc_block_sparse_matrix.h>
#    include <deal.II/lac/petsc_snes.templates.h>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  void
  NonlinearSolverData::add_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Running parameters");
    prm.add_parameter(
      "options prefix",
      options_prefix,
      "The string indicating the options prefix for command line customization.");
    prm.add_parameter("solver type",
                      snes_type,
                      "The string indicating the PETSc SNES type.");
    prm.add_parameter("linesearch type",
                      snes_linesearch_type,
                      "The string indicating the PETSc linesearch type.");
    prm.add_parameter("absolute error tolerance",
                      absolute_tolerance,
                      "Absolute error tolerance.");
    prm.add_parameter("relative error tolerance",
                      relative_tolerance,
                      "Relative error tolerance.");
    prm.add_parameter("step tolerance", step_tolerance, "Step tolerance.");
    prm.add_parameter("maximum iterations",
                      maximum_non_linear_iterations,
                      "Maximum number of iterations allowed.");
    prm.add_parameter("maximum function evaluations",
                      max_n_function_evaluations,
                      "Maximum number of function evaluations allowed.");
    prm.leave_subsection();
  }

} // namespace PETScWrappers

template class PETScWrappers::NonlinearSolver<>;
template class PETScWrappers::NonlinearSolver<PETScWrappers::MPI::Vector>;
template class PETScWrappers::NonlinearSolver<PETScWrappers::MPI::BlockVector>;
template class PETScWrappers::NonlinearSolver<PETScWrappers::MPI::Vector,
                                              PETScWrappers::MPI::SparseMatrix>;
template class PETScWrappers::NonlinearSolver<
  PETScWrappers::MPI::BlockVector,
  PETScWrappers::MPI::BlockSparseMatrix>;


DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_WITH_PETSC
#endif

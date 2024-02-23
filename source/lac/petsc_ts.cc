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

#ifndef DOXYGEN

#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_PETSC

#    include <deal.II/lac/petsc_block_sparse_matrix.h>
#    include <deal.II/lac/petsc_ts.templates.h>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  void
  TimeStepperData::add_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Running parameters");
    prm.add_parameter(
      "options prefix",
      options_prefix,
      "The string indicating the options prefix for command line customization.");
    prm.add_parameter("solver type",
                      ts_type,
                      "The string indicating the PETSc TS type.");
    prm.add_parameter("initial time",
                      initial_time,
                      "The value for the initial time.");
    prm.add_parameter("final time",
                      final_time,
                      "The value for the final time.");
    prm.add_parameter("initial step size",
                      initial_step_size,
                      "The value for the initial time step.");
    prm.add_parameter("maximum number of steps",
                      max_steps,
                      "Maximum number of time steps allowed.");
    prm.add_parameter(
      "match final time",
      match_step,
      "Whether or not to exactly stop at final time or step over it.");
    prm.leave_subsection();

    prm.enter_subsection("Error control");
    prm.add_parameter("adaptor type",
                      ts_adapt_type,
                      "The string for the TSAdapt type.");
    prm.add_parameter("minimum step size",
                      minimum_step_size,
                      "Minimum time step size allowed.");
    prm.add_parameter("maximum step size",
                      maximum_step_size,
                      "Maximum time step size allowed.");
    prm.add_parameter("absolute error tolerance",
                      absolute_tolerance,
                      "Absolute error tolerance.");
    prm.add_parameter("relative error tolerance",
                      relative_tolerance,
                      "Absolute error tolerance.");
    prm.add_parameter("ignore algebraic lte",
                      ignore_algebraic_lte,
                      "Indicate whether or not to suppress algebraic variables "
                      "in the local truncation error test.");
    prm.leave_subsection();
  }

} // namespace PETScWrappers

template class PETScWrappers::TimeStepper<>;
template class PETScWrappers::TimeStepper<PETScWrappers::MPI::Vector>;
template class PETScWrappers::TimeStepper<PETScWrappers::MPI::BlockVector>;
template class PETScWrappers::TimeStepper<PETScWrappers::MPI::Vector,
                                          PETScWrappers::MPI::SparseMatrix>;
template class PETScWrappers::TimeStepper<
  PETScWrappers::MPI::BlockVector,
  PETScWrappers::MPI::BlockSparseMatrix>;


DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_WITH_PETSC
#endif

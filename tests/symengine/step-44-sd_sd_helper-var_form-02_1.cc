// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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

// This is a modified version of step-44, which tests the implementation of
// cell-level symbolic-differentiation (via a variational formulation).
// The optimization/evaluation of symbolic expressions may be split across
// multiple optimizer instances.

// Optimizer mode: Lambda
// Symbolic expression optimizations: Disabled
// Load balancing: Disabled
// Batching: Disabled

#include "../tests.h"
#include "sd_common_tests/step-44-sd_sd_helper-var_form.h"

int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  using namespace dealii;
  using namespace Step44;
  try
    {
      const enum SD::OptimizerType     opt_method = SD::OptimizerType::lambda;
      const enum SD::OptimizationFlags opt_flags =
        SD::OptimizationFlags::optimize_default; // Increase test speed.
      const enum SD::LoadBalancing load_balancing = SD::LoadBalancing::none;
      const unsigned int           n_batches      = 1;

      const unsigned int dim = 2; // Increase test speed.
      Solid<dim, opt_method, opt_flags, load_balancing, n_batches> solid(
        SOURCE_DIR "/prm/parameters-step-44.prm");
      solid.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}

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

// #include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/petsc_ts.h>
#include <deal.II/lac/petsc_vector.h>

#include "../tests.h"

/**
 * Solves a dummy index-1 DAE to test the solver.
 */
using VectorType  = PETScWrappers::MPI::Vector;
using TimeStepper = PETScWrappers::TimeStepper<VectorType>;
using real_type   = TimeStepper::real_type;


class Dummy : public TimeStepper
{
public:
  Dummy()
  {
    // Customize solver: use BDF and basic adaptive time stepping
    PETScWrappers::TimeStepperData data;
    data.ts_type       = "bdf";
    data.final_time    = 1.0;
    data.ts_adapt_type = "basic";
    reinit(data);

    // Here we solve the two variables system:
    //   x_dot = x
    //   z = x
    // with y = (x,z)
    implicit_function = [&](const real_type   t,
                            const VectorType &y,
                            const VectorType &y_dot,
                            VectorType       &res) -> void {
      res(0) = y_dot(0) - y(0);
      res(1) = y(1) - y(0);
      res.compress(VectorOperation::insert);
    };

    // Return the index set representing the
    // algebraic components
    algebraic_components = [&]() -> IndexSet {
      IndexSet algebraic_is;

      algebraic_is.set_size(2);
      algebraic_is.add_index(1);
      return algebraic_is;
    };
  }

  unsigned int
  solve()
  {
    VectorType y(MPI_COMM_SELF, 2, 2);
    y(0) = 1.0;
    y(1) = y(0);
    y.compress(VectorOperation::insert);

    return TimeStepper::solve(y);
  }
};

int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  Dummy dae;
  dae.solve();
}

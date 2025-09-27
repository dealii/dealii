// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Make sure that we can initialize the SUNDIALS solvers
// for serial vector types without initializing MPI

#include <deal.II/lac/vector.h>

#include <deal.II/sundials/arkode.h>
#include <deal.II/sundials/ida.h>
#include <deal.II/sundials/kinsol.h>

#include "../tests.h"


int
main()
{
  initlog();
#if DEAL_II_SUNDIALS_VERSION_LT(4, 0, 0)
  SUNDIALS::IDA<Vector<double>> ida_solver;
  (void)ida_solver;
  deallog << "IDA OK" << std::endl;

  SUNDIALS::KINSOL<Vector<double>> kinsol_solver;
  (void)kinsol_solver;
  deallog << "KINSOL OK" << std::endl;
#endif

  SUNDIALS::ARKode<Vector<double>> arkode_solver;
  (void)arkode_solver;
  deallog << "ARKODE OK" << std::endl;

  return 0;
}

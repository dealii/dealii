//-----------------------------------------------------------
//
//    Copyright (C) 2017 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE at
//    the top level of the deal.II distribution.
//
//-----------------------------------------------------------

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

  SUNDIALS::IDA<Vector<double>> ida_solver;
  (void)ida_solver;
  deallog << "IDA OK" << std::endl;

  SUNDIALS::KINSOL<Vector<double>> kinsol_solver;
  (void)kinsol_solver;
  deallog << "KINSOL OK" << std::endl;

  SUNDIALS::ARKode<Vector<double>> arkode_solver;
  (void)arkode_solver;
  deallog << "ARKODE OK" << std::endl;

  return 0;
}

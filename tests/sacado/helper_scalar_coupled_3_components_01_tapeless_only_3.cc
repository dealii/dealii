// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Evaluation of a coupled system (tensor + vector + scalar components)
// using a helper class.
// This test is based off of helper_scalar_coupled_3_components_01.h, and checks
// that everything still works in tapeless only mode (i.e. when the
// start_recording_operations and stop_recording_operations calls are
// removed).
//
// AD number type: Sacado Rad

#include "../tests.h"

#include "../ad_common_tests/helper_scalar_coupled_3_components_01_tapeless_only.h"

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    mpi_log;

  deallog.push("Double");
  {
    test_tensor_vector_scalar_coupled<2, double, AD::NumberTypes::sacado_rad>();
    test_tensor_vector_scalar_coupled<3, double, AD::NumberTypes::sacado_rad>();
    deallog << "OK" << std::endl;
  }
  deallog.pop();

  deallog.push("Float");
  {
    test_tensor_vector_scalar_coupled<2, float, AD::NumberTypes::sacado_rad>();
    test_tensor_vector_scalar_coupled<3, float, AD::NumberTypes::sacado_rad>();
    deallog << "OK" << std::endl;
  }
  deallog.pop();

  deallog << "OK" << std::endl;
}

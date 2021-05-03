// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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


// Evaluation of a coupled system
//   (symmetric tensor + tensor + vector + scalar components)
// using a helper class
//
// AD number type: Sacado Rad-DFad

#include "../tests.h"

#include "../ad_common_tests/helper_scalar_coupled_4_components_01.h"

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    mpi_log;

  deallog.push("Double");
  {
    test_symmetric_tensor_tensor_vector_scalar_coupled<
      2,
      double,
      AD::NumberTypes::sacado_rad_dfad>();
    test_symmetric_tensor_tensor_vector_scalar_coupled<
      3,
      double,
      AD::NumberTypes::sacado_rad_dfad>();
    deallog << "OK" << std::endl;
  }
  deallog.pop();

  deallog.push("Float");
  {
    // test_symmetric_tensor_tensor_vector_scalar_coupled<
    //   2,
    //   float,
    //   AD::NumberTypes::sacado_rad_dfad>();
    // test_symmetric_tensor_tensor_vector_scalar_coupled<
    //   3,
    //   float,
    //   AD::NumberTypes::sacado_rad_dfad>();
    deallog << "OK" << std::endl;
  }
  deallog.pop();

  deallog << "OK" << std::endl;
}

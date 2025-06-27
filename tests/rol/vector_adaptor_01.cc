// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test to check Rol::VectorAdaptor::set() and Rol::VectorAdaptor::plus().

#include <deal.II/lac/generic_linear_algebra.h>

#include <deal.II/optimization/rol/vector_adaptor.h>

#include "../tests.h"


template <typename VectorType>
void
test(const VectorType &given_vector)
{
  ROL::Ptr<VectorType> given_vector_rcp =
    ROL::makePtr<VectorType>(given_vector);

  // --- Testing the constructor
  Rol::VectorAdaptor<VectorType> given_vector_rol(given_vector_rcp);
  AssertThrow(given_vector == *given_vector_rol.getVector(),
              ExcInternalError());


  ROL::Ptr<VectorType>           w_rcp = ROL::makePtr<VectorType>();
  Rol::VectorAdaptor<VectorType> w_rol(w_rcp);

  // --- Testing VectorAdaptor::set()
  {
    w_rol.set(given_vector_rol);
    AssertThrow(given_vector == *w_rol.getVector(), ExcInternalError());
  }

  // --- Testing VectorAdaptor::plus()
  {
    VectorType u;
    u = given_vector;
    u *= 2.;
    w_rol.plus(given_vector_rol);
    AssertThrow(u == *w_rol.getVector(), ExcInternalError());
  }

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  deallog.depth_console(10);

  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, dealii::numbers::invalid_unsigned_int);

  try
    {
      {
        LinearAlgebraTrilinos::MPI::Vector trilinos_vector;
        trilinos_vector.reinit(complete_index_set(100), MPI_COMM_WORLD);

        // set the first vector
        for (unsigned int i = 0; i < trilinos_vector.size(); ++i)
          trilinos_vector(i) = i;

        test(trilinos_vector);
      }
    }
  catch (const std::exception &exc)
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
    };
}

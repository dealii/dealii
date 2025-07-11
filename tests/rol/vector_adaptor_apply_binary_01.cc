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

// Check the Rol::VectorAdaptor's applyBinary function applied to a fully
// distributed (MPI) vector.

#include <deal.II/lac/generic_linear_algebra.h>

#include <deal.II/optimization/rol/vector_adaptor.h>

#include "../tests.h"


// Taken from deal.II's test: parallel_vector_07
template <typename VectorType>
void
prepare_vector(VectorType &v)
{
  const unsigned int myid =
                       dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD),
                     numproc =
                       dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  const unsigned int set = 10;
  AssertIndexRange(numproc, set - 2);
  const unsigned int local_size  = set - myid;
  unsigned int       global_size = 0;
  unsigned int       my_start    = 0;
  for (unsigned int i = 0; i < numproc; ++i)
    {
      global_size += set - i;
      if (i < myid)
        my_start += set - i;
    }
  // each processor owns some indices and all
  // are ghosting elements from three
  // processors (the second). some entries
  // are right around the border between two
  // processors
  IndexSet local_owned(global_size);
  local_owned.add_range(my_start, my_start + local_size);

  // --- Prepare vector.
  v.reinit(local_owned, MPI_COMM_WORLD);
}


template <typename VectorType>
void
test()
{
  VectorType a, b;
  prepare_vector(a);
  prepare_vector(b);

  Testing::srand(1);

  for (auto iterator = a.begin(); iterator != a.end(); ++iterator)
    *iterator = static_cast<double>(Testing::rand()) / RAND_MAX;

  for (auto iterator = b.begin(); iterator != b.end(); ++iterator)
    *iterator = static_cast<double>(Testing::rand()) / RAND_MAX;

  a.compress(VectorOperation::insert);
  b.compress(VectorOperation::insert);

  ROL::Ptr<VectorType> a_ptr = ROL::makePtr<VectorType>(a);
  ROL::Ptr<VectorType> b_ptr = ROL::makePtr<VectorType>(b);

  ROL::Elementwise::Multiply<double> mult;
  ROL::Elementwise::Plus<double>     plus;

  // --- Testing the constructor
  Rol::VectorAdaptor<VectorType> a_rol(a_ptr);
  Rol::VectorAdaptor<VectorType> b_rol(b_ptr);

  a_rol.print(std::cout);
  b_rol.print(std::cout);

  a_rol.applyBinary(mult, b_rol);
  a_rol.print(std::cout);

  a_rol.applyBinary(plus, b_rol);
  a_rol.print(std::cout);
}



int
main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  deallog.push(dealii::Utilities::int_to_string(myid));


  if (myid == 0)
    {
      deallog.depth_console(10); // initlog();
      deallog << std::setprecision(4);
    }

  try
    {
      test<LinearAlgebraTrilinos::MPI::Vector>();
      test<LinearAlgebra::distributed::Vector<double>>();
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

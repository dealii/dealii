// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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

// Check the Rol::VectorAdaptor with MPI fully distributed vectors
// using ROL::Vector's checkVector method.

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

  const unsigned int set = 200;
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
  VectorType a, b, c;
  prepare_vector(a);
  prepare_vector(b);
  prepare_vector(c);

  for (auto iterator = a.begin(); iterator != a.end(); iterator++)
    *iterator = static_cast<double>(Testing::rand()) / RAND_MAX;

  for (auto iterator = b.begin(); iterator != b.end(); iterator++)
    *iterator = static_cast<double>(Testing::rand()) / RAND_MAX;

  for (auto iterator = c.begin(); iterator != c.end(); iterator++)
    *iterator = static_cast<double>(Testing::rand()) / RAND_MAX;

  a.compress(VectorOperation::insert);
  b.compress(VectorOperation::insert);
  c.compress(VectorOperation::insert);

  Teuchos::RCP<VectorType> a_rcp(new VectorType(a));
  Teuchos::RCP<VectorType> b_rcp(new VectorType(b));
  Teuchos::RCP<VectorType> c_rcp(new VectorType(c));

  // --- Testing the constructor
  Rol::VectorAdaptor<VectorType> a_rol(a_rcp);
  Rol::VectorAdaptor<VectorType> b_rol(b_rcp);
  Rol::VectorAdaptor<VectorType> c_rol(c_rcp);

  Teuchos::RCP<std::ostream> out_stream;
  Teuchos::oblackholestream  bhs; // outputs nothing

  if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    out_stream = Teuchos::rcp(&std::cout, false);
  else
    out_stream = Teuchos::rcp(&bhs, false);

  a_rol.checkVector(b_rol, c_rol, true, *out_stream);
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
    };
}

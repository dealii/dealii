// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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



// test TrilinosVector::locally_owned_elements


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/trilinos_vector.h>

#include <fstream>
#include <sstream>



template <int dim>
void test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  const unsigned int n_processes = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  // create non-contiguous index set
  {
    Assert(n_processes == 2, ExcNotImplemented());
    IndexSet index (10);
    for (unsigned int i=0; i<10; i+=2)
      index.add_range(i+myid, i+myid+1);
    index.compress();

    TrilinosWrappers::MPI::Vector vec(index, MPI_COMM_WORLD);

    IndexSet index2 = vec.locally_owned_elements();
    Assert (index == index2, ExcInternalError());
  }

  // create contiguous index set
  {
    IndexSet index (10);
    index.add_range(5*myid, 5*(myid+1));
    index.compress();

    TrilinosWrappers::MPI::Vector vec(index, MPI_COMM_WORLD);

    IndexSet index2 = vec.locally_owned_elements();
    Assert (index == index2, ExcInternalError());
  }

  if (myid == 0)
    deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test<2>();
    }
  else
    {
      deallog.depth_console(0);
      test<2>();
    }

}

//----------------------------  trilinos_vector_reinit.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  trilinos_vector_reinit.cc  ---------------------------


// check correct behaviour of reinit of Trilinos vectors

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/trilinos_vector.h>
#include <fstream>
#include <iostream>
#include <vector>


void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0) deallog << "numproc=" << numproc << std::endl;


  TrilinosWrappers::MPI::Vector test1, test2;

  Assert (test1.vector_partitioner().SameAs(test2.vector_partitioner()),
	  ExcInternalError());

				   // first processor owns 2 indices, second
				   // processor owns none
  IndexSet local_owned(2);
  if (myid == 0)
    local_owned.add_range (0,2);

  test1.reinit (local_owned, MPI_COMM_WORLD);

				// reinit Trilinos vector from other vector
  test2.reinit (test1, true);

  Assert (test1.vector_partitioner().SameAs(test2.vector_partitioner()),
	  ExcInternalError());

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("trilinos_vector_reinit").c_str());
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();

}

//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// PETScWrappers: document bug when using GrowingVectorMemory
// at the end of a run.

#include "../tests.h"

#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_parallel_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <fstream>
//#include <mpi.h>

template<class v>
void testit()
{
  GrowingVectorMemory<v> m;
  m.free(m.alloc());
}



void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numprocs = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0)
    deallog << "Running on " << numprocs << " CPU(s)." << std::endl;

  testit<PETScWrappers::MPI::Vector>();
  testit<PETScWrappers::MPI::BlockVector>();
  testit<PETScWrappers::Vector>();
  testit<PETScWrappers::BlockVector>();

  if (myid==0)
    deallog << "done" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      std::ofstream logfile(output_file_for_mpi(JobIdentifier::base_name(__FILE__)).c_str());
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();
}

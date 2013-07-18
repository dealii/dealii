//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// Test correct copying  of ghosted vectors in PETScWrappers::mpi::vectors

#include "../tests.h"
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/base/index_set.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>


void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);
  
  if (myid==0)
    deallog << "numproc=" << numproc << std::endl;

				   // each processor owns 2 indices and all
				   // are ghosting Element 1 (the second)

  IndexSet local_active(numproc*2);
  local_active.add_range(myid*2,myid*2+2);
  IndexSet local_relevant(numproc*2);
  local_relevant.add_range(1,2);

  PETScWrappers::MPI::Vector vb(local_active, MPI_COMM_WORLD);
  PETScWrappers::MPI::Vector v(local_active, local_relevant, MPI_COMM_WORLD);
  PETScWrappers::MPI::Vector v2(local_active, local_relevant, MPI_COMM_WORLD);

  vb = 1.5;
  v2 = vb;
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "ghost: " << v2(1) << std::endl;
  Assert(v2(1) == 1.5, ExcInternalError());
  Assert(v2(myid*2) == 1.5, ExcInternalError());
  Assert(v2(myid*2+1) == 1.5, ExcInternalError());
  
				   // set local values
  vb(myid*2)=myid*2.0;
  vb(myid*2+1)=myid*2.0+1.0;

  vb.compress(VectorOperation::insert);
  vb*=2.0;
  v=vb;
  
  
				   // check local values
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      deallog << myid*2 << ":" << v(myid*2) << std::endl;
      deallog << myid*2+1 << ":" << v(myid*2+1) << std::endl;
    }
  
  Assert(v(myid*2) == myid*4.0, ExcInternalError());
  Assert(v(myid*2+1) == myid*4.0+2.0, ExcInternalError());
  

				   // check ghost values
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "ghost: " << v(1) << std::endl;
  Assert(v(1) == 2.0, ExcInternalError());

				   //assignment from ghosted to ghosted
  v2 = v;
  Assert(v2(1) == 2.0, ExcInternalError());
  Assert(v2(myid*2) == myid*4.0, ExcInternalError());
  Assert(v2(myid*2+1) == myid*4.0+2.0, ExcInternalError());
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "ghost: " << v2(1) << std::endl;
  
				   // done
  if (myid==0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("ghost_02").c_str());
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();
  
  
}

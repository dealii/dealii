//---------------------------------------------------------------------------
//    $Id: vector_assign_01.cc 17443 2008-10-31 19:25:47Z bangerth $
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


// Test correct handling of ghost elements in PETScWrappers::mpi::vectors

#include "../tests.h"
#include <lac/petsc_parallel_vector.h>
#include <base/index_set.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>


void test ()
{
  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::System::get_n_mpi_processes (MPI_COMM_WORLD);
  
  if (myid==0)
    deallog << "numproc=" << numproc << std::endl;

				   // each processor owns 2 indices and all
				   // are ghosting Element 1 (the second)

  IndexSet local_active(numproc*2);
  local_active.add_range(myid*2,myid*2+2);
  IndexSet local_relevant(numproc*2);
  local_relevant.add_range(1,2);

  PETScWrappers::MPI::Vector v_tmp(MPI_COMM_WORLD, local_active.size(), local_active.n_elements());
  PETScWrappers::MPI::Vector v(MPI_COMM_WORLD, local_active, local_relevant);

				   // set local values
  v(myid*2)=myid*2.0;
  v(myid*2+1)=myid*2.0+1.0;

  v.compress();
  v*=2.0;
  v.update_ghost_values();
  
  
				   // check local values
  if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      deallog << myid*2 << ":" << v(myid*2) << std::endl;
      deallog << myid*2+1 << ":" << v(myid*2+1) << std::endl;
    }
  
  Assert(v(myid*2) == myid*4.0, ExcInternalError());
  Assert(v(myid*2+1) == myid*4.0+2.0, ExcInternalError());
  

				   // check ghost values
  if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "ghost: " << v(1) << std::endl;
  Assert(v(1) == 2.0, ExcInternalError());

				   // done
  if (myid==0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  PetscInitialize(&argc,&argv,0,0);
  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);

  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("ghost_01").c_str());
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();
  
  PetscFinalize();
}

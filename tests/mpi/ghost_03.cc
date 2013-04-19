//---------------------------------------------------------------------------
//    $Id: ghost_01.cc 28403 2013-02-15 05:23:40Z heister $
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


// Test that ghosted vectors are read-only

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

  PETScWrappers::MPI::Vector vb(MPI_COMM_WORLD, local_active);
  PETScWrappers::MPI::Vector v(MPI_COMM_WORLD, local_active, local_relevant);

				   // set local values
  vb(myid*2)=myid*2.0;
  vb(myid*2+1)=myid*2.0+1.0;

  vb.compress(VectorOperation::insert);
  vb*=2.0;
  v=vb;

  Assert(!vb.has_ghost_elements(), ExcInternalError());
  Assert(v.has_ghost_elements(), ExcInternalError());

  try
    {
      v(0)=1.0;
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }
  try
    {
      v(0)*=2.0;
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }
  try
    {
      v=0.0;
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }
  try
    {
      v+=v;
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }

				   // done
  if (myid==0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  deal_II_exceptions::disable_abort_on_exception();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("ghost_03").c_str());
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();
}

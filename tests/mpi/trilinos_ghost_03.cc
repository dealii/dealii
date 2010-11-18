//----------------------------  trilinos_vector_equality_4.cc  ---------------------------
//    $Id: vector_equality_4.cc 17992 2008-12-18 03:03:20Z bangerth $
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  trilinos_vector_equality_4.cc  ---------------------------


// check correct behaviour of Trilinos ghosted vectors and l2_norm().  The
// test now checks that l2_norm() is disabled.

#include "../tests.h" 
#include <base/utilities.h>
#include <base/index_set.h>
#include <lac/trilinos_vector.h>    
#include <fstream>
#include <iostream>
#include <vector>


void test ()
{
  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::System::get_n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0) deallog << "numproc=" << numproc << std::endl;


				   // each processor owns 2 indices and all
                                   // are ghosting element 1 (the second)
  IndexSet local_active(numproc*2);
  local_active.add_range(myid*2,myid*2+2);
  IndexSet local_relevant(numproc*2);
  local_relevant = local_active;
  local_relevant.add_range(1,2);

  TrilinosWrappers::MPI::Vector v(local_active, MPI_COMM_WORLD);
  TrilinosWrappers::MPI::Vector v_tmp(local_relevant, MPI_COMM_WORLD);
  

                                     // set local values
  v(myid*2)=1.0;
  v(myid*2+1)=1.0;

  v.compress();

  v_tmp.reinit(v,false,true);

  deal_II_exceptions::disable_abort_on_exception();

  bool exc = false;
  try
    {
	v_tmp.l2_norm();
    }
  catch (TrilinosWrappers::VectorBase::ExcTrilinosError e)
    {
      exc = true;
    }
  
  Assert (exc == true, ExcInternalError());
  if (Utilities::System::get_this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv);

  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("trilinos_ghost_03").c_str());
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();

}

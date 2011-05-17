//----------------------------  trilinos_vector_equality_4.cc  ---------------------------
//    $Id$
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


// check correct behaviour of Trilinos ghosted vectors document bug: assigning
// local values into a ghosted vector does not transfer it to other machines.
// update: we now disallow changing vectors which have ghost elements. Check
// that we get that Assert().

#include "../tests.h" 
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/trilinos_vector.h>    
#include <fstream>
#include <iostream>
#include <vector>


void test ()
{
  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::System::get_n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0) deallog << "numproc=" << numproc << std::endl;

  unsigned int ghostel=(numproc>1)?2:1;

				   // each processor owns 2 indices and all
                                   // are ghosting one element
  IndexSet local_active(numproc*2);
  local_active.add_range(myid*2,myid*2+2);
  IndexSet local_relevant(numproc*2);
  local_relevant = local_active;
  local_relevant.add_range(ghostel,ghostel+1);

  TrilinosWrappers::MPI::Vector v(local_relevant, MPI_COMM_WORLD);
  

                                     // set local values
  deal_II_exceptions::disable_abort_on_exception();
  v(myid*2)=myid*2.0;
  v(myid*2+1)=myid*2.0+1.0;

  v.compress();
  v*=2.0;

  if (myid == 0)
    {
      deallog << myid*2 << ":" << v(myid*2) << std::endl;
      deallog << myid*2+1 << ":" << v(myid*2+1) << std::endl;
    }
  
  Assert(v(myid*2) == myid*4.0, ExcInternalError());
  Assert(v(myid*2+1) == myid*4.0+2.0, ExcInternalError());

                                  // check ghost values
  if (myid == 0)
    deallog << "ghost: " << v(ghostel) << std::endl;
  
  if (myid == 0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv);

  unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("trilinos_ghost_04").c_str());
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();

}

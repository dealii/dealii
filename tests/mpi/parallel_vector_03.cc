//--------------------------  parallel_vector_03.cc  -----------------------
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
//--------------------------  parallel_vector_03.cc  -----------------------

// check set of ghosts for parallel vector

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/parallel_vector.h>
#include <fstream>
#include <iostream>
#include <vector>


void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0) deallog << "numproc=" << numproc << std::endl;


                                   // each processor owns 2 indices and all
                                   // are ghosting element 1 (the second)
  IndexSet local_owned(numproc*2);
  local_owned.add_range(myid*2,myid*2+2);
  IndexSet local_relevant(numproc*2);
  local_relevant = local_owned;
  local_relevant.add_range(1,2);

  parallel::distributed::Vector<double> v(local_owned, local_relevant, MPI_COMM_WORLD);

                                     // set local values and check them
  v(myid*2)=myid*2.0;
  v(myid*2+1)=myid*2.0+1.0;

  v.compress(VectorOperation::insert);
  v*=2.0;

  Assert(v(myid*2) == myid*4.0, ExcInternalError());
  Assert(v(myid*2+1) == myid*4.0+2.0, ExcInternalError());

                                // set ghost dof on all processors, compress
                                // (insert mode)
  v(1) = 7;
  v.compress(VectorOperation::insert);

  if (myid == 0)
    {
      deallog << myid*2 << ":" << v(myid*2) << std::endl;
      deallog << myid*2+1 << ":" << v(myid*2+1) << std::endl;
    }

                                // import ghosts onto all procs
  v.update_ghost_values();
  Assert (v(1) == 7.0, ExcInternalError());

                                // check l2 norm
  const double l2_norm = v.l2_norm();
  if (myid == 0)
    deallog << "L2 norm: " << l2_norm << std::endl;

  if (myid == 0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("parallel_vector_03").c_str());
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();

}

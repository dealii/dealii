//---------------------------------------------------------------------------
//    $Id: trilinos_distribute_02.cc 29593 2013-05-25 12:29:14Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2009, 2010, 2012, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


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
      std::ofstream logfile(output_file_for_mpi("trilinos_locally_owned_elements_01").c_str());
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

//----------------------------  trilinos_vector_reinit.cc  ---------------------------
//    $Id: trilinos_vector_reinit.cc 23710 2011-05-17 04:50:10Z bangerth $
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


// document a bug in Trilinos regarding compress()
// this fails in 10.6.4, but works in 10.4.2!

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/constraint_matrix.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "../../deal.II/include/deal.II/base/index_set.h"
#include "../../deal.II/include/deal.II/lac/trilinos_vector.h"


void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0) deallog << "numproc=" << numproc << std::endl;

  IndexSet locally_owned(21);
  if (myid==0)
    locally_owned.add_range(0,9);
  else
    locally_owned.add_range(9,21);
// 
  TrilinosWrappers::MPI::Vector test1(locally_owned);
  if (myid==0)
  {
    unsigned int idx[]={0,1,2,3,4,5,6,7,8};
    double val[]={0,1,2,3,4,5,6,7,8};
    test1.add(9,idx,val);
  }
  else
  {
   {
      unsigned int idx[]={1,9,3,10,5,11,12,13,14};
    double val[]={1,9,3,10,5,11,12,13,14};
    test1.add(9,idx,val);
    }
  }
  
  test1.compress(Add);

  //TrilinosWrappers::MPI::Vector test(test1.vector_partitioner()); // works
  //TrilinosWrappers::MPI::Vector test(locally_owned); // works
  TrilinosWrappers::MPI::Vector test(test1); // fails

  test = 0;

  if (myid==0)
    test(locally_owned.nth_index_in_set(5))=7;  

  if (myid==0)
    deallog << "before compress: " << test(locally_owned.nth_index_in_set(5)) << endl;
      
   test.compress(Insert);
  
  if (myid==0)
    deallog << "after compress: " << test(locally_owned.nth_index_in_set(5)) << endl;
      
  // Trilinos produces a 0 instead of a 7 here. Why? 
  if (myid==0)
  {
    Assert(test(locally_owned.nth_index_in_set(5)) == 7, ExcInternalError());
  }

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("trilinos_compress_bug").c_str());
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();

}

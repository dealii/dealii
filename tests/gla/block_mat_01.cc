//---------------------------------------------------------------------------
//    $Id: ghost_01.cc 28440 2013-02-17 14:35:08Z heister $
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


// test LA::MPI::BlockSparseMatrix

#include "../tests.h"
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/base/index_set.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

#include "gla.h"

template <class LA> 
void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);
  
  if (myid==0)
    deallog << "numproc=" << numproc << std::endl;

  IndexSet block1(10);
  if (numproc==1)
    block1.add_range(0,10);
  
  if (myid==0)
    block1.add_range(0,7);
  if (myid==1)
    block1.add_range(7,10);

  IndexSet block2(5);
  if (numproc==1)
    block2.add_range(0,5);
  
  if (myid==0)
    block2.add_range(0,2);
  if (myid==1)
    block2.add_range(2,5);


  std::vector<IndexSet> partitioning;
  partitioning.push_back(block1);
  partitioning.push_back(block2);
  
  //LA::MPI::CompressedBlockSparsityPattern sp(partitioning);
  BlockCompressedSimpleSparsityPattern sp(partitioning);
  for (unsigned int i=0;i<15;++i)
    {
      sp.add(i,i);
      sp.add(i,1);
    }
  sp.compress();

  typename LA::MPI::BlockSparseMatrix matrix;  
  matrix.reinit (partitioning, sp, MPI_COMM_WORLD);

  matrix.add(1,1,1.3);

  matrix.compress(VectorOperation::add);

  if (myid==0)
    {
      deallog << "(0,0) = " << matrix(0,0) << std::endl;
      deallog << "(1,1) = " << matrix(1,1) << std::endl;
    }
 
				   // done
  if (myid==0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll log(__FILE__);

  {	
    deallog.push("PETSc");
    test<LA_PETSc>();
    deallog.pop();	
    deallog.push("Trilinos");
    test<LA_Trilinos>();
    deallog.pop();	
  }
  
  // compile, don't run
  //if (myid==9999)
    //  test<LA_Dummy>();
  

}

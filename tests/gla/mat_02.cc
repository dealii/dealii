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


// create LA::MPI::SparseMatrix where one CPU has no DoFs

// this is a bug in PETScWrappers right now!

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
  
  IndexSet local_active(10);
  if (myid==0)
    local_active.add_range(0,10);
  IndexSet local_relevant= local_active;
  if (myid==1)
    local_relevant.add_range(5,10);

  CompressedSimpleSparsityPattern csp (local_relevant);
  
  for (unsigned int i=0;i<10;++i)
    if (local_relevant.is_element(i))
      csp.add(i,i);

  if (myid==0)
    csp.add(0,1);
 
  typename LA::MPI::SparseMatrix mat;
  mat.reinit (local_active, local_active, csp, MPI_COMM_WORLD);

  Assert(mat.n()==10, ExcInternalError());
  Assert(mat.m()==10, ExcInternalError());

  mat.set(0,0,0.1);
  mat.set(1,1,0.1);  
  
  mat.compress(VectorOperation::insert);

  mat.add(0,1,1.0);
  
  mat.compress(VectorOperation::add);
  
				   // check local values
  if (myid==0)
    {
      deallog << "1,1 : " << mat(1,1) << std::endl;
      deallog << "0,1 : " << mat(0,1) << std::endl;
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

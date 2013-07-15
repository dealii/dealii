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


// assignment of ghost vectors

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

				   // each processor owns 2 indices and all
				   // are ghosting Element 1 (the second)

  IndexSet local_active(numproc*2);
  local_active.add_range(myid*2,myid*2+2);
  IndexSet local_relevant(numproc*2);
  local_relevant.add_range(1,2);

  typename LA::MPI::Vector vb(local_active);
  typename LA::MPI::Vector v(local_active, local_relevant);
  typename LA::MPI::Vector v2(local_active, local_relevant);

  vb = 1.0;
  v2 = vb;

				   // set local values
  vb(myid*2)=myid*2.0;
  vb(myid*2+1)=myid*2.0+1.0;

  vb.compress(VectorOperation::insert);
  vb*=2.0;
  v=vb;

  Assert(vb.size() == numproc*2, ExcInternalError());
  Assert(v.size() == numproc*2, ExcInternalError());

  Assert(!vb.has_ghost_elements(), ExcInternalError());
  Assert(v.has_ghost_elements(), ExcInternalError());
  
				   // check local values
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      deallog << myid*2 << ":" << v(myid*2) << std::endl;
      deallog << myid*2+1 << ":" << v(myid*2+1) << std::endl;
    }
  
				   //assignment from ghosted to ghosted
  v2 = v;
  Assert(v2(1) == 2.0, ExcInternalError());
  Assert(v2(myid*2) == myid*4.0, ExcInternalError());
  Assert(v2(myid*2+1) == myid*4.0+2.0, ExcInternalError());


  
  Assert(v(myid*2) == myid*4.0, ExcInternalError());
  Assert(v(myid*2+1) == myid*4.0+2.0, ExcInternalError());
  

				   // check ghost values
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "ghost: " << v(1) << std::endl;
  Assert(v(1) == 2.0, ExcInternalError());

				   // done
  if (myid==0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("vec_02").c_str());
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      {	
	deallog.push("PETSc");
	test<LA_PETSc>();
	deallog.pop();	
	deallog.push("Trilinos");
	test<LA_Trilinos>();
	deallog.pop();	
      }
      
    }
  else
      {	
	deallog.push("PETSc");
	test<LA_PETSc>();
	deallog.pop();	
	deallog.push("Trilinos");
	test<LA_Trilinos>();
	deallog.pop();	
      }

  if (myid==9999)
    test<LA_Dummy>();
  

}

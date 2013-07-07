//---------------------------------------------------------------------------
//    $Id: constraint_matrix_condense_01.cc 29076 2013-03-27 15:08:37Z heister $
//    Version: $Name$
//
//    Copyright (C) 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// check ConstraintMatrix::set_zero(Vector) for parallel vectors.
// this documents a bug introduced in r29678 for block vectors
// that was fixed in 29940.

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/petsc_parallel_block_vector.h>

#include <fstream>


void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);
 
  std::vector<IndexSet> local_active(2);

  // block 0:
  local_active[0].set_size(numproc);
  local_active[0].add_range(myid,myid+1);

  // block 1:
  local_active[1].set_size(2*numproc);
  local_active[1].add_range(myid*2,myid*2+2);

  PETScWrappers::MPI::BlockVector v(local_active, MPI_COMM_WORLD);

  for (unsigned int i=0;i<v.size();++i)
    v(i)=1.0+i;
  v.compress(VectorOperation::insert);

  IndexSet local_active_together(3*numproc);
  local_active_together.add_range(myid,myid+1);
  local_active_together.add_range(numproc+myid*2,numproc+myid*2+2);
  
  ConstraintMatrix cm(local_active_together);
  cm.add_line(numproc + myid*2);
  cm.close();

  deallog << "vector before:" << std::endl;  
  v.print(deallog.get_file_stream());

  deallog << "CM:" << std::endl;  
  cm.print(deallog.get_file_stream());
  
  cm.set_zero(v);
  
  deallog << "vector after:" << std::endl;  
  v.print(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll log(__FILE__);

  test();
  return 0;
}

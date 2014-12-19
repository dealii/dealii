// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// check SparsityTools::distribute_sparsity_pattern

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/base/index_set.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>

#include <fstream>


void test_mpi()
{
  Assert( Utilities::System::job_supports_mpi(), ExcInternalError());


  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  const unsigned int numprocs = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0)
    deallog << "Running on " << numprocs << " CPU(s)." << std::endl;

  unsigned int num_local=10;
  unsigned int n=numprocs*num_local;
  std::vector<types::global_dof_index> rows_per_cpu;
  for (unsigned int i=0; i<numprocs; ++i)
    rows_per_cpu.push_back(num_local);

  IndexSet locally_rel(n);
  locally_rel.add_range(myid*num_local, (myid+1)*num_local);
  if (myid>0)
    locally_rel.add_range((myid-1)*num_local, (myid+0)*num_local);
  if (myid<numprocs-1)
    locally_rel.add_range((myid+1)*num_local, (myid+2)*num_local);

  CompressedSimpleSparsityPattern csp(n,n, locally_rel);

  for (unsigned int i=0; i<n; ++i)
    csp.add(i, myid);

  SparsityTools::distribute_sparsity_pattern<>(csp,
                                               rows_per_cpu,
                                               MPI_COMM_WORLD,
                                               locally_rel);
  /*  {
      std::ofstream f((std::string("after")+Utilities::int_to_string(myid)).c_str());
      csp.print(f);
      }*/

  // checking...
  for (unsigned int r=0; r<num_local; ++r)
    {
      unsigned int indx=r+myid*num_local;
      unsigned int len=csp.row_length(indx);

//std::cout << "myid=" << myid << " idx=" << indx << " len=" << len <<std::endl;

      if (myid>0 && myid<numprocs-1)
        Assert(len==3, ExcInternalError());
      if (myid==0 || myid==numprocs-1)
        Assert(len==2, ExcInternalError());

      Assert(csp.exists(indx, myid), ExcInternalError());
      if (myid>0)
        Assert(csp.exists(indx, myid-1), ExcInternalError());
      if (myid<numprocs-1)
        Assert(csp.exists(indx, myid+1), ExcInternalError());
    }

  if (myid==0)
    deallog << "done" << std::endl;

}


int main(int argc, char *argv[])
{
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi (argc, argv);
#else
  (void)argc;
  (void)argv;
  compile_time_error;

#endif

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("mpi");
      test_mpi();
      deallog.pop();
    }
  else
    test_mpi();
}

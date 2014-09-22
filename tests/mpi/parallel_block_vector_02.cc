// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2013 by the deal.II authors
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


// check ghost handling on parallel block vectors

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/parallel_block_vector.h>
#include <deal.II/lac/parallel_vector.h>
#include <fstream>
#include <iostream>
#include <vector>


void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0) deallog << "numproc=" << numproc << std::endl;


  // each processor from processor 1 to 8 owns 2 indices (the other processors
  // do not own any dof), and all processors are ghosting element 1
  IndexSet local_owned(std::min(16U,numproc*2));
  if (myid < 8)
    local_owned.add_range(myid*2,myid*2+2);
  IndexSet local_relevant(numproc*2);
  local_relevant = local_owned;
  local_relevant.add_range(1,2);

  parallel::distributed::Vector<double> v(local_owned, local_relevant,
                                          MPI_COMM_WORLD);

  // set local values
  if (myid < 8)
    {
      v(myid*2)=myid*2.0;
      v(myid*2+1)=myid*2.0+1.0;
    }
  v.compress(VectorOperation::insert);

  parallel::distributed::BlockVector<double> w(3);
  for (unsigned int i=0; i<3; ++i)
    {
      w.block(i) = v;
      w.block(i) *= (i+1);
    }
  w.collect_sizes();

  // see if we can access the ghosted entry 1 in each block with the correct
  // value: the initialization should also update the ghost values, so all
  // processors should have i+1
  for (unsigned int i=0; i<3; ++i)
    AssertDimension(i+1,(unsigned int)w.block(i)(1));

  // import ghost values, all processors should still have i+1
  w.update_ghost_values();
  for (unsigned int i=0; i<3; ++i)
    AssertDimension(i+1,(unsigned int)w.block(i)(1));

  // zero out ghosts, now all processors except processor 1 should have 0.
  w.zero_out_ghosts();
  for (unsigned int i=0; i<3; ++i)
    if (myid == 0)
      {
        AssertDimension(i+1,(unsigned int)w.block(i)(1));
      }
    else
      {
        AssertDimension(0,(unsigned int)w.block(i)(1));
      }

  // create a vector copy that gets the entries from w. First, it should have
  // updated the ghosts because it is created from an empty state.
  parallel::distributed::BlockVector<double> x(w);
  Assert(x.has_ghost_elements() == true, ExcInternalError());
  for (unsigned int i=0; i<3; ++i)
    AssertDimension(i+1,(unsigned int)x.block(i)(1));

  // now we zero the vector, which should disable ghost elements
  x = 0;
  Assert(x.has_ghost_elements() == false, ExcInternalError());

  // we copy from w (i.e., the same vector but one that does not have ghosts
  // enabled) -> should not have ghosts enabled
  x = w;
  Assert(x.has_ghost_elements() == false, ExcInternalError());
  for (unsigned int i=0; i<3; ++i)
    if (myid == 0)
      {
        AssertDimension(i+1,(unsigned int)x.block(i)(1));
      }
    else
      {
        AssertDimension(0,(unsigned int)x.block(i)(1));
      }

  x.update_ghost_values();
  Assert(x.has_ghost_elements() == true, ExcInternalError());


  // add something to entry 1 on all processors
  w(1) += myid+1;
  w.compress(VectorOperation::add);
  if (myid == 0)
    AssertDimension((unsigned int)w(1), 1+(numproc*(numproc+1))/2);

  // add again and check if everything is still correct
  w(1+v.size()) += myid+1;
  w.compress(VectorOperation::add);
  if (myid == 0)
    AssertDimension((unsigned int)w(1), 1+(numproc*(numproc+1))/2);
  if (myid == 0)
    AssertDimension((unsigned int)w(v.size()+1), 2+(numproc*(numproc+1))/2);

  w.update_ghost_values();
  AssertDimension((unsigned int)w(1), 1+(numproc*(numproc+1))/2);
  AssertDimension((unsigned int)w(v.size()+1), 2+(numproc*(numproc+1))/2);

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
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();

}

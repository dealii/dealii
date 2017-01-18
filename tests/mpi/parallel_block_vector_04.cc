// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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


// Check global reduction operations (norms, operator==, operator!=) with the
// parallel block vector class. This is the same as parallel_vector_01 except
// that this test uses the version of the add() method that takes a std::vector
// of indices and a std::vector of values to build the block vector instead of
// vector assignment.

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <fstream>
#include <iostream>
#include <vector>


void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0) deallog << "numproc=" << numproc << std::endl;


  // each processor from processor 1 to 8
  // owns 2 indices (the other processors do
  // not own any dof), and all processors are
  // ghosting element 1 (the second)
  IndexSet local_owned(std::min(16U,numproc*2));
  if (myid < 8)
    local_owned.add_range(myid*2,myid*2+2);
  IndexSet local_relevant(numproc*2);
  local_relevant = local_owned;
  local_relevant.add_range(1,2);

  LinearAlgebra::distributed::Vector<double> v(local_owned, local_relevant,
                                               MPI_COMM_WORLD);
  std::vector<types::global_dof_index> indices;
  std::vector<double> values;

  // set local values
  if (myid < 8)
    {
      // unlike the first test, we scale values here since these are added to
      // w directly
      indices.push_back(myid*2);
      values.push_back(myid*4.0);
      indices.push_back(myid*2 + 1);
      values.push_back(myid*4.0 + 2.0);
    }

  LinearAlgebra::distributed::BlockVector<double> w(3);
  for (unsigned int i=0; i<3; ++i)
    w.block(i) = v;
  w.collect_sizes();

  // set up v for comparison purposes below
  v.add(indices, values);
  v.compress(VectorOperation::add);
  // in analogy to the first test, v is scaled differently than w

  if (myid < 8)
    {
      for (unsigned int i = 0; i<3; ++i)
        {
          // update the values in w
          w.add(indices, values);
          // update indices to point to the next block
          indices[0] += v.size();
          indices[1] += v.size();
        }
    }
  w.compress(VectorOperation::add);

  // check l2 norm
  {
    const double l2_norm = w.l2_norm();
    if (myid == 0)
      deallog << "l2 norm: " << l2_norm << std::endl;
    AssertThrow (std::abs(v.l2_norm()*std::sqrt(3.)-w.l2_norm()) < 1e-13,
                 ExcInternalError());
  }

  // check l1 norm
  {
    const double l1_norm = w.l1_norm();
    if (myid == 0)
      deallog << "l1 norm: " << l1_norm << std::endl;
    AssertThrow (std::abs(v.l1_norm()*3.-w.l1_norm()) < 1e-14,
                 ExcInternalError());
  }

  // check linfty norm
  {
    const double linfty_norm = w.linfty_norm();
    if (myid == 0)
      deallog << "linfty norm: " << linfty_norm << std::endl;
    AssertThrow (v.linfty_norm()==w.linfty_norm(),
                 ExcInternalError());
  }

  // check lp norm
  {
    const double lp_norm = w.lp_norm(2.2);
    if (myid == 0)
      deallog << "l2.2 norm: " << lp_norm << std::endl;

    AssertThrow (std::fabs (w.l2_norm() - w.lp_norm(2.0)) < 1e-14,
                 ExcInternalError());
  }

  // check mean value (should be equal to l1
  // norm divided by vector size here since we
  // have no negative entries)
  {
    const double mean = w.mean_value();
    if (myid == 0)
      deallog << "Mean value: " << mean << std::endl;

    AssertThrow (std::fabs (mean * w.size() - w.l1_norm()) < 1e-15,
                 ExcInternalError());
  }
  // check inner product
  {
    const double norm_sqr = w.l2_norm() * w.l2_norm();
    AssertThrow (std::fabs(w * w - norm_sqr) < 1e-12,
                 ExcInternalError());
    LinearAlgebra::distributed::BlockVector<double> w2;
    w2 = w;
    AssertThrow (std::fabs(w2 * w - norm_sqr) < 1e-12,
                 ExcInternalError());

    if (myid<8)
      w2.block(0).local_element(0) = -1;
    const double inner_prod = w * w2;
    if (myid == 0)
      deallog << "Inner product: " << inner_prod << std::endl;
  }

  // check all_zero
  {
    bool allzero = w.all_zero();
    if (myid == 0)
      deallog << " v==0 ? " << allzero << std::endl;
    LinearAlgebra::distributed::BlockVector<double> w2;
    w2.reinit (w);
    allzero = w2.all_zero();
    if (myid == 0)
      deallog << " v2==0 ? " << allzero << std::endl;

    // now change one element to nonzero
    if (myid == 0)
      w2.block(1).local_element(1) = 1;
    allzero = w2.all_zero();
    if (myid == 0)
      deallog << " v2==0 ? " << allzero << std::endl;
  }

  if (myid == 0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, testing_max_num_threads());

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();

}

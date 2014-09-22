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


// check global reduction operation (norms, operator==, operator!=) on
// parallel vector

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

  parallel::distributed::Vector<double> v(local_owned, local_owned, MPI_COMM_WORLD);

  // set local values
  if (myid < 8)
    {
      v(myid*2)=myid*2.0;
      v(myid*2+1)=myid*2.0+1.0;
    }
  v.compress(VectorOperation::insert);
  v*=2.0;
  if (myid < 8)
    {
      Assert(v(myid*2) == myid*4.0, ExcInternalError());
      Assert(v(myid*2+1) == myid*4.0+2.0, ExcInternalError());
    }

  // check l2 norm
  {
    const double l2_norm = v.l2_norm();
    if (myid == 0)
      deallog << "l2 norm: " << l2_norm << std::endl;
  }

  // check l1 norm
  {
    const double l1_norm = v.l1_norm();
    if (myid == 0)
      deallog << "l1 norm: " << l1_norm << std::endl;
  }

  // check linfty norm
  {
    const double linfty_norm = v.linfty_norm();
    if (myid == 0)
      deallog << "linfty norm: " << linfty_norm << std::endl;
  }

  // check lp norm
  {
    const double lp_norm = v.lp_norm(2.2);
    if (myid == 0)
      deallog << "l2.2 norm: " << lp_norm << std::endl;

    Assert (std::fabs (v.l2_norm() - v.lp_norm(2.0)) < 1e-14,
            ExcInternalError());
  }

  // check mean value (should be equal to l1
  // norm divided by vector size here since we
  // have no negative entries)
  {
    const double mean = v.mean_value();
    if (myid == 0)
      deallog << "Mean value: " << mean << std::endl;

    Assert (std::fabs (mean * v.size() - v.l1_norm()) < 1e-15,
            ExcInternalError());
  }
  // check inner product
  {
    const double norm_sqr = v.norm_sqr();
    Assert (std::fabs(v * v - norm_sqr) < 1e-15,
            ExcInternalError());
    parallel::distributed::Vector<double> v2;
    v2 = v;
    Assert (std::fabs(v2 * v - norm_sqr) < 1e-15,
            ExcInternalError());

    if (myid<8)
      v2.local_element(0) = -1;
    const double inner_prod = v * v2;
    if (myid == 0)
      deallog << "Inner product: " << inner_prod << std::endl;
  }

  // check operator ==
  {
    parallel::distributed::Vector<double> v2 (v);
    bool equal = (v2 == v);
    if (myid == 0)
      deallog << " v==v2 ? " << equal << std::endl;

    bool not_equal = (v2 != v);
    if (myid == 0)
      deallog << " v!=v2 ? " << not_equal << std::endl;

    // change v2 on one proc only
    if (myid == 0)
      v2.local_element(1) = 2.2212;

    equal = (v2 == v);
    if (myid == 0)
      deallog << " v==v2 ? " << equal << std::endl;
    not_equal = (v2 != v);
    if (myid == 0)
      deallog << " v!=v2 ? " << not_equal << std::endl;

    // reset
    v2 = v;
    equal = (v2 == v);
    if (myid == 0)
      deallog << " v==v2 ? " << equal << std::endl;
    not_equal = (v2 != v);
    if (myid == 0)
      deallog << " v!=v2 ? " << not_equal << std::endl;

    // change some value on all procs
    if (myid < 8)
      v2.local_element(0) = -1;
    equal = (v2 == v);
    if (myid == 0)
      deallog << " v==v2 ? " << equal << std::endl;
    not_equal = (v2 != v);
    if (myid == 0)
      deallog << " v!=v2 ? " << not_equal << std::endl;
  }

  // check all_zero
  {
    bool allzero = v.all_zero();
    if (myid == 0)
      deallog << " v==0 ? " << allzero << std::endl;
    parallel::distributed::Vector<double> v2;
    v2.reinit (v);
    allzero = v2.all_zero();
    if (myid == 0)
      deallog << " v2==0 ? " << allzero << std::endl;

    // now change one element to nonzero
    if (myid == 0)
      v2.local_element(1) = 1;
    allzero = v2.all_zero();
    if (myid == 0)
      deallog << " v2==0 ? " << allzero << std::endl;
  }


  // check all_non_negative
  {
    bool allnonneg = v.is_non_negative();
    if (myid == 0)
      deallog << " v>=0 ? " << allnonneg << std::endl;
    parallel::distributed::Vector<double> v2, v3;

    // vector where all processors have
    // non-negative entries
    v2 = v;
    if (myid < 8)
      v2.local_element(0) = -1;
    allnonneg = v2.is_non_negative();
    if (myid == 0)
      deallog << " v2>=0 ? " << allnonneg << std::endl;

    // zero vector
    v3.reinit (v2);
    allnonneg = v3.is_non_negative();
    if (myid == 0)
      deallog << " v3>=0 ? " << allnonneg << std::endl;

    // only one processor has non-negative entry
    v3 = v;
    if (myid == 1 || numproc==1)
      v3.local_element(0) = -1;
    allnonneg = v3.is_non_negative();
    if (myid == 0)
      deallog << " v3>=0 ? " << allnonneg << std::endl;
  }

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

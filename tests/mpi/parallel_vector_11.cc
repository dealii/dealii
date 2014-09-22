// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2013 by the deal.II authors
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


// check that add, sadd, equ, scale work correctly on a vector where some
// processor do not own any degrees of freedom

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

  // global size: 20, local_size: 3 as long as
  // less than 20
  const unsigned int local_size = 3;
  const unsigned int global_size = std::min(20U, local_size * numproc);
  const int my_start = std::min (local_size * myid, global_size);
  const int my_end   = std::min (local_size * (myid+1), global_size);
  const int actual_local_size = my_end-my_start;

  IndexSet local_owned (global_size);
  if (my_end > my_start)
    local_owned.add_range(static_cast<unsigned int>(my_start),
                          static_cast<unsigned int>(my_end));
  IndexSet local_relevant(global_size);
  local_relevant = local_owned;
  local_relevant.add_index (2);

  parallel::distributed::Vector<double> v(local_owned, local_relevant,
                                          MPI_COMM_WORLD);
  AssertDimension (static_cast<unsigned int>(actual_local_size), v.local_size());
  parallel::distributed::Vector<double> w (v), x(v), y(v);

  // set local elements
  for (int i=0; i<actual_local_size; ++i)
    {
      v.local_element(i) = i + my_start;
      w.local_element(i) = 1000 + 2 * (my_start + i);
      x.local_element(i) = 10000;
    }

  y = v;
  for (int i=0; i<actual_local_size; ++i)
    Assert (y.local_element(i) == i+my_start, ExcInternalError());

  if (myid==0) deallog << "Check add (scalar): ";
  y.add (42);
  for (int i=0; i<actual_local_size; ++i)
    Assert (y.local_element(i) == i+my_start+42, ExcInternalError());
  if (myid==0) deallog << "OK" << std::endl;

  if (myid==0) deallog << "Check add (vector): ";
  y.add (w);
  for (int i=0; i<actual_local_size; ++i)
    Assert (y.local_element(i) == 3*(i+my_start)+1042, ExcInternalError());
  if (myid==0) deallog << "OK" << std::endl;

  if (myid==0) deallog << "Check add (factor, vector): ";
  y.add (-1., w);
  for (int i=0; i<actual_local_size; ++i)
    Assert (y.local_element(i) == i+my_start+42, ExcInternalError());
  if (myid==0) deallog << "OK" << std::endl;

  if (myid==0) deallog << "Check add (factor, vector, factor, vector): ";
  y.add (2., w, -0.5, x);
  for (int i=0; i<actual_local_size; ++i)
    Assert (y.local_element(i) == 5*(i+my_start)+2042-5000,ExcInternalError());
  if (myid==0) deallog << "OK" << std::endl;

  if (myid==0) deallog << "Check sadd (factor, factor, vector): ";
  y = v;
  y.sadd (-3.,2.,v);
  for (int i=0; i<actual_local_size; ++i)
    Assert (y.local_element(i)==(-i-my_start), ExcInternalError());
  if (myid==0) deallog << "OK" << std::endl;

  if (myid==0) deallog << "Check sadd (factor, factor, vector, factor, vector): ";
  y.sadd (2.,3.,v, 2., w);
  for (int i=0; i<actual_local_size; ++i)
    Assert (y.local_element(i) == 5*(i+my_start)+2000, ExcInternalError());
  if (myid==0) deallog << "OK" << std::endl;

  if (myid==0) deallog << "Check sadd (factor, factor, vector, factor, vector, factor, vector): ";
  y.sadd (-1.,1.,v, 2., w, 2., x);
  for (int i=0; i<actual_local_size; ++i)
    Assert (y.local_element(i) == 20000, ExcInternalError());
  if (myid==0) deallog << "OK" << std::endl;

  if (myid==0) deallog << "Check add (factor, vector_1, factor, vector_1): ";
  y = 0;
  y.add (1.,v, 3., v);
  for (int i=0; i<actual_local_size; ++i)
    Assert (y.local_element(i) == 4*(i+my_start), ExcInternalError());
  if (myid==0) deallog << "OK" << std::endl;

  if (myid==0) deallog << "Check operator * (scalar): ";
  x *= 2.;
  for (int i=0; i<actual_local_size; ++i)
    Assert (x.local_element(i) == 20000., ExcInternalError());
  if (myid==0) deallog << "OK" << std::endl;

  if (myid==0) deallog << "Check operator / (scalar): ";
  x /= 2.;
  for (int i=0; i<actual_local_size; ++i)
    Assert (x.local_element(i) == 10000., ExcInternalError());
  if (myid==0) deallog << "OK" << std::endl;

  if (myid==0) deallog << "Check scale (vector): ";
  y.scale (x);
  for (int i=0; i<actual_local_size; ++i)
    Assert (y.local_element(i) == 40000.*(i+my_start), ExcInternalError());
  if (myid==0) deallog << "OK" << std::endl;

  if (myid==0) deallog << "Check equ (factor, vector): ";
  y. equ (10., x);
  for (int i=0; i<actual_local_size; ++i)
    Assert (y.local_element(i) == 100000., ExcInternalError());
  if (myid==0) deallog << "OK" << std::endl;

  if (myid==0) deallog << "Check equ (factor, vector, factor, vector): ";
  y. equ (10., v, -2., w);
  for (int i=0; i<actual_local_size; ++i)
    Assert (y.local_element(i) == 6.*(i+my_start)-2000, ExcInternalError());
  if (myid==0) deallog << "OK" << std::endl;

  if (myid==0) deallog << "Check equ (factor, vector, factor, vector, factor, vector): ";
  y. equ (10., v, -2., w, 3., x);
  for (int i=0; i<actual_local_size; ++i)
    Assert (y.local_element(i) == 6.*(i+my_start)+28000, ExcInternalError());
  if (myid==0) deallog << "OK" << std::endl;

  if (myid==0) deallog << "Check equ<float> (factor, vector): ";
  parallel::distributed::Vector<float> z;
  z = v;
  y.equ (1., z);
  for (int i=0; i<actual_local_size; ++i)
    Assert (y.local_element(i) == i+my_start, ExcInternalError());
  if (myid==0) deallog << "OK" << std::endl;
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

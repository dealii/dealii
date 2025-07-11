// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check that add, sadd, equ, scale work correctly on a vector where some
// processor do not own any degrees of freedom

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/read_write_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "numproc=" << numproc << std::endl;

  // global size: 20, local_size: 3 as long as
  // less than 20
  const unsigned int local_size  = 3;
  const unsigned int global_size = std::min(20U, local_size * numproc);
  const int          my_start    = std::min(local_size * myid, global_size);
  const int          my_end = std::min(local_size * (myid + 1), global_size);
  const int          actual_local_size = my_end - my_start;

  IndexSet local_owned(global_size);
  if (my_end > my_start)
    local_owned.add_range(static_cast<unsigned int>(my_start),
                          static_cast<unsigned int>(my_end));
  IndexSet local_relevant(global_size);
  local_relevant = local_owned;
  local_relevant.add_index(2);

  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> v(
    local_owned, local_relevant, MPI_COMM_WORLD);
  AssertDimension(static_cast<unsigned int>(actual_local_size),
                  v.locally_owned_size());
  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> w(v), x(v),
    y(v);

  // set local elements
  LinearAlgebra::ReadWriteVector<double> v_rw(local_owned);
  LinearAlgebra::ReadWriteVector<double> w_rw(local_owned);
  LinearAlgebra::ReadWriteVector<double> x_rw(local_owned);
  for (int i = 0; i < actual_local_size; ++i)
    {
      v_rw.local_element(i) = i + my_start;
      w_rw.local_element(i) = 1000 + 2 * (my_start + i);
      x_rw.local_element(i) = 10000;
    }
  v.import_elements(v_rw, VectorOperation::insert);
  w.import_elements(w_rw, VectorOperation::insert);
  x.import_elements(x_rw, VectorOperation::insert);

  y = v;
  LinearAlgebra::ReadWriteVector<double> y_rw(local_owned);
  y_rw.import_elements(y, VectorOperation::insert);
  for (int i = 0; i < actual_local_size; ++i)
    AssertThrow(y_rw.local_element(i) == i + my_start, ExcInternalError());

  if (myid == 0)
    deallog << "Check add (scalar): ";
  y.add(42);
  y_rw.import_elements(y, VectorOperation::insert);
  for (int i = 0; i < actual_local_size; ++i)
    AssertThrow(y_rw.local_element(i) == i + my_start + 42, ExcInternalError());
  if (myid == 0)
    deallog << "OK" << std::endl;

  if (myid == 0)
    deallog << "Check add (vector): ";
  y.add(1., w);
  y_rw.import_elements(y, VectorOperation::insert);
  for (int i = 0; i < actual_local_size; ++i)
    AssertThrow(y_rw.local_element(i) == 3 * (i + my_start) + 1042,
                ExcInternalError());
  if (myid == 0)
    deallog << "OK" << std::endl;

  if (myid == 0)
    deallog << "Check add (factor, vector): ";
  y.add(-1., w);
  y_rw.import_elements(y, VectorOperation::insert);
  for (int i = 0; i < actual_local_size; ++i)
    AssertThrow(y_rw.local_element(i) == i + my_start + 42, ExcInternalError());
  if (myid == 0)
    deallog << "OK" << std::endl;

  if (myid == 0)
    deallog << "Check add (factor, vector, factor, vector): ";
  y.add(2., w, -0.5, x);
  y_rw.import_elements(y, VectorOperation::insert);
  for (int i = 0; i < actual_local_size; ++i)
    AssertThrow(y_rw.local_element(i) == 5 * (i + my_start) + 2042 - 5000,
                ExcInternalError());
  if (myid == 0)
    deallog << "OK" << std::endl;

  if (myid == 0)
    deallog << "Check sadd (factor, factor, vector): ";
  y = v;
  y.sadd(-3., 2., v);
  y_rw.import_elements(y, VectorOperation::insert);
  for (int i = 0; i < actual_local_size; ++i)
    AssertThrow(y_rw.local_element(i) == (-i - my_start), ExcInternalError());
  if (myid == 0)
    deallog << "OK" << std::endl;

  if (myid == 0)
    deallog << "Check sadd (factor, factor, vector, factor, vector): ";
  y.sadd(2., 3., v);
  y.add(2., w);
  y_rw.import_elements(y, VectorOperation::insert);
  for (int i = 0; i < actual_local_size; ++i)
    {
      AssertThrow(y_rw.local_element(i) == 5 * (i + my_start) + 2000,
                  ExcInternalError());
    }
  if (myid == 0)
    deallog << "OK" << std::endl;

  if (myid == 0)
    deallog
      << "Check sadd (factor, factor, vector, factor, vector, factor, vector): ";
  y.sadd(-1., 1., v);
  y.add(2., w);
  y.add(2., x);
  y_rw.import_elements(y, VectorOperation::insert);
  for (int i = 0; i < actual_local_size; ++i)
    AssertThrow(y_rw.local_element(i) == 20000, ExcInternalError());
  if (myid == 0)
    deallog << "OK" << std::endl;

  if (myid == 0)
    deallog << "Check add (factor, vector_1, factor, vector_1): ";
  y = 0;
  y.add(1., v, 3., v);
  y_rw.import_elements(y, VectorOperation::insert);
  for (int i = 0; i < actual_local_size; ++i)
    AssertThrow(y_rw.local_element(i) == 4 * (i + my_start),
                ExcInternalError());
  if (myid == 0)
    deallog << "OK" << std::endl;

  if (myid == 0)
    deallog << "Check operator * (scalar): ";
  x *= 2.;
  x_rw.import_elements(x, VectorOperation::insert);
  for (int i = 0; i < actual_local_size; ++i)
    AssertThrow(x_rw.local_element(i) == 20000., ExcInternalError());
  if (myid == 0)
    deallog << "OK" << std::endl;

  if (myid == 0)
    deallog << "Check operator / (scalar): ";
  x /= 2.;
  x_rw.import_elements(x, VectorOperation::insert);
  for (int i = 0; i < actual_local_size; ++i)
    AssertThrow(x_rw.local_element(i) == 10000., ExcInternalError());
  if (myid == 0)
    deallog << "OK" << std::endl;

  if (myid == 0)
    deallog << "Check scale (vector): ";
  y.scale(x);
  y_rw.import_elements(y, VectorOperation::insert);
  for (int i = 0; i < actual_local_size; ++i)
    AssertThrow(y_rw.local_element(i) == 40000. * (i + my_start),
                ExcInternalError());
  if (myid == 0)
    deallog << "OK" << std::endl;

  if (myid == 0)
    deallog << "Check equ (factor, vector): ";
  y.equ(10., x);
  y_rw.import_elements(y, VectorOperation::insert);
  for (int i = 0; i < actual_local_size; ++i)
    AssertThrow(y_rw.local_element(i) == 100000., ExcInternalError());
  if (myid == 0)
    deallog << "OK" << std::endl;

  if (myid == 0)
    deallog << "Check equ (factor, vector, factor, vector): ";
  y.equ(10., v);
  y.add(-2., w);
  y_rw.import_elements(y, VectorOperation::insert);
  for (int i = 0; i < actual_local_size; ++i)
    AssertThrow(y_rw.local_element(i) == 6. * (i + my_start) - 2000,
                ExcInternalError());
  if (myid == 0)
    deallog << "OK" << std::endl;

  if (myid == 0)
    deallog << "Check equ (factor, vector, factor, vector, factor, vector): ";
  y.equ(10., v);
  y.add(-2., w);
  y.add(3., x);
  y_rw.import_elements(y, VectorOperation::insert);
  for (int i = 0; i < actual_local_size; ++i)
    AssertThrow(y_rw.local_element(i) == 6. * (i + my_start) + 28000,
                ExcInternalError());
  if (myid == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();
      deallog << std::setprecision(4);

      test();
    }
  else
    test();
}

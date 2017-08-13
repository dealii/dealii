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

// the same as timer_04.cc but this time with MPI.

#include "../tests.h"
#include <deal.II/base/timer.h>

// burn computer time

double s = 0.;
void burn (unsigned int n)
{
  for (unsigned int i=0 ; i<n ; ++i)
    {
      for (unsigned int j=1 ; j<100000 ; ++j)
        {
          s += 1./j * i;
        }
    }
}



void
test_timer(Timer &t)
{
  burn(50);

  const double old_wall_time = t.wall_time();
  AssertThrow(old_wall_time > 0., ExcInternalError());
  const double old_cpu_time = t.wall_time();
  AssertThrow(old_cpu_time > 0., ExcInternalError());
  AssertThrow(t.get_total_data().max == 0., ExcInternalError());

  burn(50);
  AssertThrow(t.stop() > 0., ExcInternalError());

  AssertThrow(t.wall_time() > old_wall_time, ExcInternalError());
  AssertThrow(t.cpu_time() > old_cpu_time, ExcInternalError());
  AssertThrow(t.last_wall_time() > 0., ExcInternalError());
  AssertThrow(t.last_cpu_time() > 0., ExcInternalError());
  AssertThrow(t.get_data().min > 0., ExcInternalError());
  AssertThrow(t.get_total_data().min > 0., ExcInternalError());

  t.reset();
  AssertThrow(t.wall_time() == 0., ExcInternalError());
  AssertThrow(t.cpu_time()== 0., ExcInternalError());
  AssertThrow(t.get_total_data().max == 0, ExcInternalError());

  deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();

  Timer t1(MPI_COMM_WORLD, false);
  Timer t2(MPI_COMM_WORLD, true);

  test_timer(t1);
  test_timer(t2);
}


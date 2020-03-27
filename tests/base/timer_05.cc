// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// the same as timer_04.cc but this time with MPI.

#include <deal.II/base/timer.h>

#include "../tests.h"

// burn computer time

double s = 0.;
void
burn(unsigned int n)
{
  for (unsigned int i = 0; i < n; ++i)
    {
      for (unsigned int j = 1; j < 100000; ++j)
        {
          s += 1. / j * i;
        }
    }
}

// check that the MinMaxAvg is set to a default, unpopulated state.
void
assert_min_max_avg_invalid(const Utilities::MPI::MinMaxAvg &data)
{
  // we want to explicitly check that some things are signaling NaNs, so
  // disable FP exceptions if they are on
#if defined(DEBUG) && defined(DEAL_II_HAVE_FP_EXCEPTIONS)
  fedisableexcept(FE_INVALID);
#endif

  AssertThrow(std::isnan(data.min), ExcInternalError());
  AssertThrow(std::isnan(data.max), ExcInternalError());
  AssertThrow(std::isnan(data.avg), ExcInternalError());
  AssertThrow(data.min_index == numbers::invalid_unsigned_int,
              ExcInternalError());
  AssertThrow(data.max_index == numbers::invalid_unsigned_int,
              ExcInternalError());
}

// check that the MinMaxAvg values are reasonable.
void
assert_min_max_avg_valid(const Utilities::MPI::MinMaxAvg &data)
{
  AssertThrow(data.min > 0., ExcInternalError());
  AssertThrow(data.max >= data.min, ExcInternalError());
  AssertThrow(data.avg >= data.min, ExcInternalError());
  AssertThrow(data.min_index != numbers::invalid_unsigned_int,
              ExcInternalError());
  AssertThrow(data.max_index != numbers::invalid_unsigned_int,
              ExcInternalError());
}


void
test_timer(Timer &t)
{
  burn(50);

  const double old_wall_time = t.wall_time();
  AssertThrow(old_wall_time > 0., ExcInternalError());
  const double old_cpu_time = t.cpu_time();
  AssertThrow(old_cpu_time > 0., ExcInternalError());
  assert_min_max_avg_invalid(t.get_last_lap_wall_time_data());
  assert_min_max_avg_invalid(t.get_accumulated_wall_time_data());

  burn(50);
  AssertThrow(t.stop() > 0., ExcInternalError());
  assert_min_max_avg_valid(t.get_last_lap_wall_time_data());
  assert_min_max_avg_valid(t.get_accumulated_wall_time_data());

  AssertThrow(t.wall_time() > old_wall_time, ExcInternalError());
  AssertThrow(t.cpu_time() > old_cpu_time, ExcInternalError());
  AssertThrow(t.last_wall_time() > 0., ExcInternalError());
  AssertThrow(t.last_cpu_time() > 0., ExcInternalError());

  t.reset();
  AssertThrow(t.wall_time() == 0., ExcInternalError());
  AssertThrow(t.cpu_time() == 0., ExcInternalError());
  assert_min_max_avg_invalid(t.get_last_lap_wall_time_data());
  assert_min_max_avg_invalid(t.get_accumulated_wall_time_data());

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();

  Timer t1(MPI_COMM_WORLD, false);
  Timer t2(MPI_COMM_WORLD, true);

  test_timer(t1);
  test_timer(t2);
}

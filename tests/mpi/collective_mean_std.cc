// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check Utilities::MPI::mean_and_standard_deviation()

#include <deal.II/base/mpi.templates.h>

#include "../tests.h"

/* MWE in Python
>>> import numpy as np

>>> a = np.array([-3, 55, 0, 1, 11, -12])
>>> np.mean(a)
8.6666666666666661
>>> np.std(a, ddof=1)
23.871880249923059

>>> b = np.array([-1+2.j, 3+7.j, -5.j,6])
>>> np.mean(b)
(2+1j)
>>> np.std(b,ddof=1)
5.8878405775518976

>>> c = np.array([1,1,1])
>>> np.mean(c)
1.0
>>> np.std(c,ddof=1)
0.0

>>> d = np.array([1,2])
>>> np.mean(d)
1.5
>>> np.std(d,ddof=1)
0.70710678118654757
*/

void
test()
{
  unsigned int     myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  std::vector<int> values;
  std::vector<std::complex<double>> val_c;
  std::vector<int>                  same;
  std::vector<int>                  empty;
  if (myid == 0)
    {
      values.resize(4);
      values[0] = -3;
      values[1] = 55;
      values[2] = 0;
      values[3] = 1;

      val_c.resize(1);
      val_c[0].real(-1);
      val_c[0].imag(2.);

      same.resize(2);
      same[0] = 1;
      same[1] = 1;

      empty.resize(2);
      empty[0] = 1;
      empty[1] = 2;
    }
  else
    {
      values.resize(2);
      values[0] = 11;
      values[1] = -12;

      val_c.resize(3);
      val_c[0].real(3);
      val_c[0].imag(7);
      val_c[1].real(0);
      val_c[1].imag(-5);
      val_c[2].real(6);
      val_c[2].imag(0);

      same.resize(1);
      same[0] = 1;
    }

  const auto pair = Utilities::MPI::mean_and_standard_deviation(values.begin(),
                                                                values.end(),
                                                                MPI_COMM_WORLD);

  const auto pair_c =
    Utilities::MPI::mean_and_standard_deviation<decltype(val_c.begin()),
                                                std::complex<double>>(
      val_c.begin(), val_c.end(), MPI_COMM_WORLD);

  const auto pair_same =
    Utilities::MPI::mean_and_standard_deviation(same.begin(),
                                                same.end(),
                                                MPI_COMM_WORLD);

  const auto pair_empty =
    Utilities::MPI::mean_and_standard_deviation(empty.begin(),
                                                empty.end(),
                                                MPI_COMM_WORLD);

  if (myid == 0)
    deallog << pair.first << ' ' << pair.second << std::endl
            << pair_c.first << ' ' << pair_c.second << std::endl
            << pair_same.first << ' ' << pair_same.second << std::endl
            << pair_empty.first << ' ' << pair_empty.second << std::endl;
}


int
main(int argc, char *argv[])
{
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
#else
  (void)argc;
  (void)argv;
  compile_time_error;

#endif

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      initlog();

      deallog.push("mpi");
      test();
      deallog.pop();
    }
  else
    test();
}

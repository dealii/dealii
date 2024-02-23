// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check that LinearAlgebra::distributed::Vector::add_and_dot works correctly
// for complex-valued vectors

#include <deal.II/base/numbers.h>

#include <deal.II/lac/la_parallel_vector.h>

#include "../tests.h"



template <typename number>
void
check()
{
  for (unsigned int test = 0; test < 5; ++test)
    {
      const unsigned int size         = 17 + test * 1101;
      const IndexSet     complete_set = complete_index_set(size);
      LinearAlgebra::distributed::Vector<std::complex<number>> v1(
        complete_set, MPI_COMM_SELF),
        v2(complete_set, MPI_COMM_SELF), v3(complete_set, MPI_COMM_SELF),
        check(complete_set, MPI_COMM_SELF);
      for (unsigned int i = 0; i < size; ++i)
        {
          v1(i) = std::complex<number>(0.1 + 0.005 * i, 1.234 + 12 * i);
          v2(i) = std::complex<number>(-5.2 + 0.18 * i, 42.4242 + 42 * i);
          v3(i) =
            std::complex<number>(numbers::PI + numbers::E / (1. + i), 13.);
        }
      check                             = v1;
      const std::complex<number> factor = std::complex<number>(0.01432);

      // do things by hand once
      v1.add(factor, v2);
      const std::complex<number> prod = v1 * v3;

      // then do it a second time with the add_and_dot function
      const std::complex<number> prod_check = check.add_and_dot(factor, v2, v3);
      if (test == 0 && std::is_same_v<number, double>)
        {
          deallog << "Vector add reference:   ";
          v1.print(deallog.get_file_stream());
          deallog << std::endl;
          deallog << "Vector check reference: ";
          check.print(deallog.get_file_stream());
          deallog << std::endl;
        }

      deallog << "Add and dot is ";
      // check tolerance with respect to the expected size of result which is
      // ~ size^2 including the roundoff error ~ sqrt(size) we expect
      const number tolerance = 4. * std::numeric_limits<number>::epsilon() *
                               std::sqrt(static_cast<number>(size)) * size *
                               size;
      if (std::abs(prod - prod_check) < tolerance)
        deallog << "correct" << std::endl;
      else
        deallog << "wrong by " << std::abs(prod - prod_check)
                << " with tolerance " << tolerance << "; should be "
                << prod / static_cast<number>(size) << ", is "
                << prod_check / static_cast<number>(size) << std::endl;
    }
}


int
main()
{
  initlog();
  deallog << std::setprecision(8);
  check<float>();
  check<double>();
  deallog << "OK" << std::endl;
}

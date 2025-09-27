// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check that the l2 norm is exactly the same for many runs on random vector
// size with random vector entries. this is to ensure that the result is
// reproducible also when parallel evaluation is done

#include <deal.II/lac/vector.h>

#include "../tests.h"



template <typename number>
void
check_norms()
{
  for (unsigned int test = 0; test < 20; ++test)
    {
      const unsigned int size = Testing::rand() % 100000;
      Vector<number>     vec(size);
      for (unsigned int i = 0; i < size; ++i)
        vec(i) = random_value<number>();
      const typename Vector<number>::real_type norm = vec.l2_norm();
      for (unsigned int i = 0; i < 30; ++i)
        AssertThrow(vec.l2_norm() == norm, ExcInternalError());

      Vector<number> vec2(vec);
      for (unsigned int i = 0; i < 10; ++i)
        AssertThrow(vec2.l2_norm() == norm, ExcInternalError());
    }
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(2);
  deallog.attach(logfile);

  check_norms<float>();
  check_norms<double>();
  check_norms<long double>();
  check_norms<std::complex<double>>();
  deallog << "OK" << std::endl;
}

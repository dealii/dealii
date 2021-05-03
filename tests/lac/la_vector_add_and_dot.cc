// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2019 by the deal.II authors
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


// check that LinearAlgebra::Vector::add_and_dot works correctly

#include <deal.II/lac/la_vector.h>

#include "../tests.h"



template <typename number>
void
check()
{
  for (unsigned int test = 0; test < 5; ++test)
    {
      const unsigned int            size = 17 + test * 1101;
      LinearAlgebra::Vector<number> v1(size), v2(size), v3(size), check(size);
      // Check that the assignment works
      v1 = 0.;
      for (unsigned int i = 0; i < size; ++i)
        {
          v1[i] = 0.1 + 0.005 * i;
          v2[i] = -5.2 + 0.18 * i;
          v3[i] = 3.14159 + 2.7183 / (1. + i);
        }
      check               = v1;
      const number factor = 0.01432;

      v1.add(factor, v2);
      const number prod       = v1 * v3;
      const number prod_check = check.add_and_dot(factor, v2, v3);
      if (test == 0 && std::is_same<number, double>::value)
        {
          deallog << "Vector add reference:   ";
          for (unsigned int i = 0; i < size; ++i)
            deallog << v1[i] << " ";
          deallog << std::endl;
          deallog << "Vector check reference: ";
          for (unsigned int i = 0; i < size; ++i)
            deallog << check[i] << " ";
          deallog << std::endl;

          const number constant = 1.;
          v1.add(constant);
          deallog << "Vector add constant:    ";
          for (unsigned int i = 0; i < size; ++i)
            deallog << v1[i] << " ";
          deallog << std::endl;
        }

      deallog << "Add and dot is ";
      if (std::abs(prod - prod_check) <
          4. *
            std::abs(
              std::numeric_limits<
                typename numbers::NumberTraits<number>::real_type>::epsilon()) *
            std::sqrt(static_cast<double>(size)) * size)
        deallog << "correct" << std::endl;
      else
        deallog << "wrong; should be " << prod / static_cast<number>(size)
                << ", is " << prod_check / static_cast<number>(size)
                << std::endl;
    }
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(10);
  deallog.attach(logfile);

  check<float>();
  check<double>();
#ifdef DEAL_II_WITH_COMPLEX_VALUES
  check<std::complex<double>>();
#endif
  deallog << "OK" << std::endl;
}

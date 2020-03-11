// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2019 by the deal.II authors
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


// check that the inner product is exactly the same down to roundoff for
// various vectors.
//
// It was previously possible to create an unaligned view (a VectorView
// object) into a Vector: an original goal of this test was to verify that
// unaligned dot products were identical to aligned ones, but since Vector now
// uses AlignedVector to store its data it is no longer possible to test this.

#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector.templates.h>

#include "../tests.h"



template <typename number>
void
check_norms()
{
  for (unsigned int test = 0; test < 5; ++test)
    {
      const unsigned int size = Testing::rand() % 20000;
      Vector<number> larger1(size + 8), larger2(size + 8), in1(size), in2(size);
      for (unsigned int i = 0; i < size; ++i)
        {
          in1(i) = random_value<number>();
          in2(i) = random_value<number>();
        }

      const number inner_product = in1 * in2;
      deallog << "Vector difference: ";
      for (unsigned int shift1 = 0; shift1 < 8; ++shift1)
        for (unsigned int shift2 = 0; shift2 < 8; ++shift2)
          {
            Vector<number> v1(larger1.begin() + shift1,
                              larger1.begin() + shift1 + size);
            Vector<number> v2(larger2.begin() + shift2,
                              larger2.begin() + shift2 + size);
            for (unsigned int i = 0; i < size; ++i)
              {
                v1(i) = in1(i);
                v2(i) = in2(i);
              }

            const number result = v1 * v2;
            deallog << static_cast<double>(std::abs(result - inner_product))
                    << " ";
          }
      deallog << std::endl;
    }
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  check_norms<float>();
  check_norms<double>();
  check_norms<long double>();
#ifdef DEAL_II_WITH_COMPLEX_VALUES
  check_norms<std::complex<double>>();
#endif
}

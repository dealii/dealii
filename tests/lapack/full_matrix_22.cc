// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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


// test LAPACKFullMatrix::reciprocal_condition_number() for triangular matrices

/* MWE for size=3 in Octave:
R = [1,2,3; 0, 5, 6; 0, 0, 9]
rcond(R)
ans =  0.055556
*/

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

#include <iostream>

#include "../tests.h"

#include "create_matrix.h"


template <typename NumberType>
void
test(const unsigned int size)
{
  LAPACKFullMatrix<NumberType> M(size);
  M.set_property(LAPACKSupport::upper_triangular);

  M                    = 0.;
  unsigned int counter = 1;
  for (unsigned int i = 0; i < size; ++i)
    for (unsigned int j = 0; j < size; ++j)
      {
        if (j >= i)
          M(i, j) = counter;

        counter++;
      }

  // M.print_formatted(deallog.get_file_stream(), 3, false, 8);
  const double rcond = M.reciprocal_condition_number();

  deallog << rcond << std::endl;
}


int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);

  const std::vector<unsigned int> sizes = {{1, 3, 11}};
  for (const auto &s : sizes)
    {
      deallog << "size=" << s << std::endl;
      // test<float>(s);
      test<double>(s);
    }
}

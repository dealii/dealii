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


// test LAPACKFullMatrix:: norms for non-symmetric matrices

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
  // Full matrix:
  FullMatrix<NumberType> F(size);
  create_random(F);

  // Lapack:
  LAPACKFullMatrix<NumberType> M(size);
  M = F;

  deallog << "l1:        " << (F.l1_norm() - M.l1_norm()) << std::endl
          << "linfty:    " << (F.linfty_norm() - M.linfty_norm()) << std::endl
          << "frobenius: " << (F.frobenius_norm() - M.frobenius_norm())
          << std::endl;
}


int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);

  const std::vector<unsigned int> sizes = {{1, 3, 11, 17, 32, 64, 200, 391}};
  for (const auto &s : sizes)
    {
      deallog << "size=" << s << std::endl;
      // test<float>(s);
      test<double>(s);
    }
}

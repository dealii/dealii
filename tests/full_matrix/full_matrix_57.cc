// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check FullMatrix::outer_product for correct functionality


#include "../tests.h"

#include "full_matrix_common.h"



const double ints[9] = {0, -1., 1., -2., 2., -3., 3., -4., 4.};

template <typename number>
void
check()
{
  deallog << std::fixed;
  deallog << std::setprecision(1);

  Vector<double>     First4(4), Second4(4), First9(9);
  FullMatrix<double> F4(4, 4), F9(9, 9);


  deallog << "Original Vector V" << std::endl;
  deallog << "=================" << std::endl;

  for (unsigned int i = 0; i < First4.size(); ++i)
    {
      First4(i) = ints[i];
      deallog << First4(i) << "  ";
    }
  deallog << std::endl;
  deallog << std::endl;

  deallog << "Original Vector W" << std::endl;
  deallog << "=================" << std::endl;

  for (unsigned int i = 0; i < Second4.size(); ++i)
    {
      Second4(i) = ints[i + 1];
      deallog << Second4(i) << "  ";
    }
  deallog << std::endl;
  deallog << std::endl;

  F4.outer_product(First4, Second4);

  deallog << "Outer_Product of V and W" << std::endl;
  deallog << "========================" << std::endl;
  display_matrix(F4);
  deallog << std::endl;
  deallog << std::endl;

  deallog << std::endl;

  F4.outer_product(First4, First4);

  deallog << "Outer_Product of V and V" << std::endl;
  deallog << "========================" << std::endl;
  display_matrix(F4);
  deallog << std::endl;
  deallog << std::endl;


  F4.outer_product(Second4, Second4);

  deallog << "Outer_Product of W and W" << std::endl;
  deallog << "========================" << std::endl;
  display_matrix(F4);
  deallog << std::endl;
  deallog << std::endl;
  deallog << std::endl;

  deallog << "Vector V" << std::endl;
  deallog << "========" << std::endl;

  for (unsigned int i = 0; i < First9.size(); ++i)
    {
      First9(i) = ints[i];
      deallog << First9(i) << "  ";
    }
  deallog << std::endl;
  deallog << std::endl;

  F9.outer_product(First9, First9);

  deallog << "Outer_Product of V and V" << std::endl;
  deallog << "========================" << std::endl;
  display_matrix(F9);
  deallog << std::endl;
  deallog << std::endl;

  deallog << std::endl;

  deallog << "Vector V" << std::endl;
  deallog << "========" << std::endl;

  for (unsigned int i = 0; i < First9.size(); ++i)
    {
      First9(i) = ints[8 - i];
      deallog << First9(i) << "  ";
    }
  deallog << std::endl;
  deallog << std::endl;

  F9.outer_product(First9, First9);

  deallog << "Outer_Product of V and V" << std::endl;
  deallog << "========================" << std::endl;
  display_matrix(F9);
  deallog << std::endl;
  deallog << std::endl;

  deallog << std::endl;
}

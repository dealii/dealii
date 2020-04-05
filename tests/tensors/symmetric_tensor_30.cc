// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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

// Verify that SymmetricTensor::operator() and SymmetricTensor::operator[] do
// what we think they should do.

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/table_indices.h>

#include "../tests.h"


template <int dim>
void
check_2()
{
  Testing::rand(true, 42);

  SymmetricTensor<2, dim> change_with_brackets;
  SymmetricTensor<2, dim> change_with_parentheses;
  for (unsigned int k = 0; k < dim; ++k)
    {
      for (unsigned int l = 0; l < dim; ++l)
        {
          const double    entry = double(Testing::rand());
          TableIndices<2> indices(k, l);

          // check assignment
          change_with_brackets[indices]    = entry;
          change_with_parentheses(indices) = entry;
          AssertThrow(change_with_brackets[k][l] == entry,
                      ExcMessage("Entries should match"));
          AssertThrow(change_with_parentheses[k][l] == entry,
                      ExcMessage("Entries should match"));

          // and access
          const double brackets_entry    = change_with_brackets[indices];
          const double parentheses_entry = change_with_parentheses(indices);
          AssertThrow(brackets_entry == entry,
                      ExcMessage("Entries should match"));
          AssertThrow(parentheses_entry == entry,
                      ExcMessage("Entries should match"));
        }
    }
}

template <int dim>
void
check_4()
{
  Testing::rand(true, 42);

  SymmetricTensor<4, dim> change_with_brackets;
  SymmetricTensor<4, dim> change_with_parentheses;
  for (unsigned int i = 0; i < dim; ++i)
    {
      for (unsigned int j = 0; j < dim; ++j)
        {
          for (unsigned int k = 0; k < dim; ++k)
            {
              for (unsigned int l = 0; l < dim; ++l)
                {
                  const double    entry = double(Testing::rand());
                  TableIndices<4> indices(i, j, k, l);

                  // check assignment
                  change_with_brackets[indices]    = entry;
                  change_with_parentheses(indices) = entry;
                  AssertThrow(change_with_brackets[i][j][k][l] == entry,
                              ExcMessage("Entries should match"));
                  AssertThrow(change_with_parentheses[i][j][k][l] == entry,
                              ExcMessage("Entries should match"));

                  // and access
                  const double brackets_entry = change_with_brackets[indices];
                  const double parentheses_entry =
                    change_with_parentheses(indices);
                  AssertThrow(brackets_entry == entry,
                              ExcMessage("Entries should match"));
                  AssertThrow(parentheses_entry == entry,
                              ExcMessage("Entries should match"));
                }
            }
        }
    }
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  deallog << "check rank 2 tensors" << std::endl;
  check_2<1>();
  check_2<2>();
  check_2<3>();

  deallog << "check rank 4 tensors" << std::endl;
  check_4<1>();
  check_4<2>();
  check_4<3>();

  deallog << "OK" << std::endl;
}

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check assignment between block vectors and regular vectors

#include <deal.II/lac/block_vector.h>

#include <algorithm>
#include <numeric>
#include <utility>
#include <vector>

#include "../tests.h"

namespace dealii::Testing
{
  template <typename Vector1, typename Vector2>
  bool
  compare_equal(const Vector1 &v1, const Vector2 &v2)
  {
    if (v1.size() != v2.size())
      return false;
    for (unsigned int i = 0; i < v1.size(); ++i)
      if (v1(i) != v2(i))
        return false;
    return true;
  }
} // namespace dealii::Testing


void
test()
{
  using namespace dealii::Testing;

  std::vector<types::global_dof_index> ivector(4);
  ivector[0] = 2;
  ivector[1] = 4;
  ivector[2] = 3;
  ivector[3] = 5;

  BlockVector<double> v1(ivector);
  Vector<double>      v2(v1.size());

  for (unsigned int i = 0; i < v1.size(); ++i)
    v1(i) = 1 + i * i;

  v2 = v1;
  AssertThrow(compare_equal(v1, v2), ExcInternalError());

  BlockVector<double> v3(ivector);
  v3 = v2;
  AssertThrow(compare_equal(v3, v2), ExcInternalError());
  AssertThrow(compare_equal(v3, v1), ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();
  deallog << std::setprecision(3) << std::fixed;

  try
    {
      test();
    }
  catch (const std::exception &e)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << e.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      // abort
      return 2;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      // abort
      return 3;
    };


  return 0;
}

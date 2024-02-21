// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check VectorTools::subtract_mean_value() for deal.II serial vectors

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <vector>

#include "../tests.h"

template <class VectorType>
void
test(VectorType &v)
{
  std::vector<bool> filter(v.size(), false);
  // set some elements of the vector
  for (unsigned int i = 0; i < v.size(); i += 1 + i)
    {
      filter[i] = true;
      v(i)      = i;
    }

  // then check the norm
  VectorTools::subtract_mean_value(v, filter);
  AssertThrow(std::fabs(v.mean_value()) < 1e-10 * v.l2_norm(),
              ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  try
    {
      {
        Vector<double> v(10);
        test(v);
      }

      {
        Vector<float> v(10);
        test(v);
      }

      {
        BlockVector<double> v(std::vector<types::global_dof_index>(1, 10));
        test(v);
      }

      {
        BlockVector<float> v(std::vector<types::global_dof_index>(1, 10));
        test(v);
      }
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
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
      return 1;
    };
}

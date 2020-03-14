// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2018 by the deal.II authors
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


// check assignment between block vectors and regular vectors

#include <deal.II/lac/block_vector.h>

#include <algorithm>
#include <numeric>
#include <utility>
#include <vector>

#include "../tests.h"

template <typename Vector1, typename Vector2>
bool
operator==(const Vector1 &v1, const Vector2 &v2)
{
  if (v1.size() != v2.size())
    return false;
  for (unsigned int i = 0; i < v1.size(); ++i)
    if (v1(i) != v2(i))
      return false;
  return true;
}


void
test()
{
  std::vector<types::global_dof_index> ivector(4);
  ivector[0] = 2;
  ivector[1] = 4;
  ivector[2] = 3;
  ivector[3] = 5;

  BlockVector<std::complex<double>> v1(ivector);
  Vector<std::complex<double>>      v2(v1.size());

  for (unsigned int i = 0; i < v1.size(); ++i)
    v1(i) = 1 + i * i;

  v2 = v1;
  AssertThrow(v1 == v2, ExcInternalError());

  BlockVector<std::complex<double>> v3(ivector);
  v3 = v2;
  AssertThrow(v3 == v2, ExcInternalError());
  AssertThrow(v3 == v1, ExcInternalError());

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
  catch (std::exception &e)
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

// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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


// check IdentityMatrix::vmult and friends


#include <deal.II/lac/identity_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"


template <typename number>
void
check_vmult()
{
  IdentityMatrix M(4);
  Vector<number> u(4);
  Vector<number> v(4);

  for (unsigned int i = 0; i < 4; ++i)
    u(i) = i + 1;

  M.vmult(v, u);
  Assert(v == u, ExcInternalError());
  for (unsigned int i = 0; i < v.size(); ++i)
    deallog << ' ' << v(i);
  deallog << std::endl;

  M.vmult_add(v, u);
  v /= 2;
  Assert(v == u, ExcInternalError());
  for (unsigned int i = 0; i < v.size(); ++i)
    deallog << ' ' << v(i);
  deallog << std::endl;

  M.Tvmult(v, u);
  Assert(v == u, ExcInternalError());
  for (unsigned int i = 0; i < v.size(); ++i)
    deallog << ' ' << v(i);
  deallog << std::endl;

  M.Tvmult_add(v, u);
  v /= 2;
  Assert(v == u, ExcInternalError());
  for (unsigned int i = 0; i < v.size(); ++i)
    deallog << ' ' << v(i);
  deallog << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(0) << std::fixed;

  check_vmult<double>();
  check_vmult<float>();
}

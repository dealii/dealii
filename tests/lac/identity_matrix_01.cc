// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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

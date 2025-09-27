// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check running over const iterators starting at the second line. like the
// full_matrix_* tests, but use complex-valued matrices and vectors; this time
// we actually store complex values in them


#include "../tests.h"

#include "full_matrix_common.h"



template <typename number>
void
check()
{
  FullMatrix<std::complex<number>> m;
  make_complex_matrix(m);


  for (typename FullMatrix<std::complex<number>>::const_iterator p = m.begin(1);
       p != m.end(1);
       ++p)
    deallog << p->row() << ' ' << p->column() << ' ' << p->value() << std::endl;
}

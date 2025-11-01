// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test AffineConstraints<double>::shift for a AffineConstraints<double> object
// initialized with an IndexSet object

#include <deal.II/base/index_set.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"


void
test()
{
  const int size = 20;

  dealii::IndexSet index_set(size);
  dealii::IndexSet set1(size);
  deallog << "Create index set with entries: ";
  index_set.add_index(2);
  index_set.add_index(5);
  index_set.add_index(8);
  index_set.print(deallog);

  deallog
    << "Create AffineConstraints<double> with constraints u(2)=.5*u(5), u(5)=.7*u(8)"
    << std::endl;
  dealii::AffineConstraints<double> constraints1(index_set, index_set);
  constraints1.add_line(5);
  constraints1.add_entry(5, 8, .7);
  constraints1.add_line(2);
  constraints1.add_entry(2, 5, .5);
  dealii::AffineConstraints<double> constraints2(constraints1);
  constraints1.print(deallog.get_file_stream());

  constraints1.shift(size / 2);
  deallog << "Shifted constraints" << std::endl;
  constraints1.print(deallog.get_file_stream());

  constraints1.merge(constraints2,
                     AffineConstraints<double>::no_conflicts_allowed,
                     true);
  deallog << "Shifted and merged constraints" << std::endl;
  constraints1.print(deallog.get_file_stream());

  constraints1.close();
  deallog << "Close" << std::endl;
  constraints1.print(deallog.get_file_stream());

  Vector<double> vec(size);
  for (unsigned int i = 0; i < size; ++i)
    vec(i) = i;
  constraints1.distribute(vec);
  vec.print(deallog.get_file_stream(), 3, true, false);
}

int
main()
{
  initlog();

  test();

  return 0;
}

// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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



// Merge two ConstraintMatrix objects initialized with IndexSets where entries
// have been added in a way that local_lines extends beyond the end of the
// IndexSet

#include <deal.II/lac/constraint_matrix.h>

#include "../tests.h"



void
merge_check()
{
  deallog << "Checking ConstraintMatrix::merge with localized lines"
          << std::endl;

  // set local lines to a very large range that surely triggers an error if
  // the implementation is wrong
  IndexSet local_lines(100000000);
  local_lines.add_range(99999800, 100000000);
  local_lines.compress();

  // works correctly
  ConstraintMatrix c1(local_lines), c2(local_lines);
  for (types::global_dof_index i = 99999800; i < local_lines.size(); ++i)
    if (i % 2 == 1)
      c1.add_line(i);

  c2.add_line(99999800);
  c2.add_line(99999801);
  c2.add_line(99999802);

  // now merge the two and print the results
  c2.merge(c1, ConstraintMatrix::right_object_wins);
  c2.print(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  merge_check();
}

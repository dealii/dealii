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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// This is the same as constraints_merge_08 but local_lines differs for
// the objects to be merged.

#include "../tests.h"
#include <deal.II/lac/constraint_matrix.h>



void
merge_check ()
{
  deallog << "Checking ConstraintMatrix::merge with localized lines" << std::endl;

  // set local lines to a very large range that
  // surely triggers an error if the
  // implementation is wrong
  IndexSet local_lines1 (100000000);
  local_lines1.add_range (99999890, 99999900);
  local_lines1.add_range (99999990, 99999992);
  local_lines1.compress();
  local_lines1.print(std::cout);

  IndexSet local_lines2 (100000000);
  local_lines2.add_range (99999893, 99999900);
  local_lines2.add_range (99999990, 100000000);
  local_lines2.compress();
  local_lines2.print(std::cout);

  // the test is the same as
  // constraints_merge_02, but we add very large
  // indices here
  const unsigned int index_0  = 99999890;
  const unsigned int index_1  = 99999891;
  const unsigned int index_3  = 99999893;
  const unsigned int index_4  = 99999894;
  const unsigned int index_10 = 99999990;
  const unsigned int index_11 = 99999991;
  const unsigned int index_12 = 99999992;
  const unsigned int index_13 = 99999993;

  // check twice, once with closed
  // objects, once with open ones
  for (unsigned int run=0; run<2; ++run)
    {
      deallog << "Checking with " << (run == 0 ? "open" : "closed")
              << " objects" << std::endl;

      // check that the `merge' function
      // works correctly
      ConstraintMatrix c1 (local_lines1), c2 (local_lines2);

      // enter simple line
      c1.add_line          (index_0);
      c1.add_entry         (index_0, index_11, 1.);
      c1.set_inhomogeneity (index_0, 42);

      // add more complex line
      c1.add_line          (index_1);
      c1.add_entry         (index_1, index_3,  0.5);
      c1.add_entry         (index_1, index_4,  0.5);
      c1.set_inhomogeneity (index_1, 100);

      // fill second constraints
      // object with one trivial line
      // and one which further
      // constrains one of the
      // entries in the first object
      c2.add_line          (index_10);
      c2.add_entry         (index_10, index_11, 1.);
      c2.set_inhomogeneity (index_10, 142);

      c2.add_line          (index_3);
      c2.add_entry         (index_3, index_12, 0.25);
      c2.add_entry         (index_3, index_13, 0.75);
      c2.set_inhomogeneity (index_3, 242);

      // in one of the two runs,
      // close the objects
      if (run == 1)
        {
          c1.close ();
          c2.close ();
        };

      // now merge the two and print the
      // results
      c1.merge (c2, ConstraintMatrix::no_conflicts_allowed, true);
      c1.print(deallog.get_file_stream());
    };
}


int
main ()
{
  initlog();

  merge_check ();
}


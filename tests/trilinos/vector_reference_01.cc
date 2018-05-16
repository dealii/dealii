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



// TrilinosWrappers::internal::VectorReference had its non-const copy
// operator return a const reference. That was non-intuitive but,
// because of the particular semantics of copying vector reference
// objects, made no difference. Either way, we now check these semantics.

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_vector.h>
#include <iostream>


void
test ()
{
  TrilinosWrappers::MPI::Vector v;
  v.reinit(complete_index_set(3), MPI_COMM_WORLD);
  v(0) = 0;
  v(1) = 1;
  v(2) = 2;

  TrilinosWrappers::internal::VectorReference a (v(0));
  TrilinosWrappers::internal::VectorReference b (v(1));
  TrilinosWrappers::internal::VectorReference c (v(2));

  // Copy the VectorReference objects. Note that operator= for these
  // objects does not copy the *content* of the reference object, but
  // rather assigns the value of the vector element referenced on the
  // right to the vector element referenced on the left.
  // So the applied operations result in the following actions:
  // 1. We first assign c to point to where a points (c points to entry 0)
  // 2. We next assign c to point to where b points (c points to entry 1)
  // 3. Lastly, we assign b to point to where a points (b points to entry 0)
  (c = a) = b;
  b = a;

  deallog << static_cast<TrilinosScalar>(a) << std::endl; // should point to v(0)
  deallog << static_cast<TrilinosScalar>(b) << std::endl; // should point to v(0)
  deallog << static_cast<TrilinosScalar>(c) << std::endl; // should point to v(1)
}



int
main (int argc,char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, testing_max_num_threads());
  initlog();;

  test ();
}

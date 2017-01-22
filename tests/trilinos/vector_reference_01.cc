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
// operator return a const reference. that was non-intuitive but,
// because of the particular semantics of copying vector reference
// objects, made no difference. either way, have a check for these
// semantics

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_vector.h>
#include <fstream>
#include <iostream>


void test ()
{
  TrilinosWrappers::Vector v(3);
  v(0) = 0;
  v(1) = 1;
  v(2) = 2;

  TrilinosWrappers::internal::VectorReference a (v(0));
  TrilinosWrappers::internal::VectorReference b (v(1));
  TrilinosWrappers::internal::VectorReference c (v(2));

  // copy VectorReference objects. note that operator= for these
  // objects does not copy the *content* of the reference object, but
  // assigns the value of the vector element referenced on the right
  // to the vector element referenced on the left
  (c = a) = b;
  b = a;

  deallog << static_cast<TrilinosScalar>(a) << std::endl;
  deallog << static_cast<TrilinosScalar>(b) << std::endl;
  deallog << static_cast<TrilinosScalar>(c) << std::endl;
}



int main (int argc,char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, testing_max_num_threads());
  MPILogInitAll log;
  
  test ();
}

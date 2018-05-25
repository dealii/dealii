// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2017 by the deal.II authors
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



// check TrilinosWrappers::MPI::Vector::operator() in set/add-mode
// alternatingly. this test doesn't really make sense any more -- at
// least one a single processor. on multiple processors, one has to
// call compress() and it used to be that we had code that called
// compress() either with VectorOperation::insert or
// VectorOperation::add automatically whenever we switched from one
// kind of operation to the other. this is no longer the case -- the
// user has to select which operation herself. on the other hand, as
// long as we're working on a single processor, it doesn't matter
// anyway.
//
// so the test doesn't make much sense any more for the original
// purpose, but it's still worth keeping it around just to test that
// set and write operations on a single processor go right through to
// the stored data, rather than ending up in separate caches that only
// get flushed when the right kind of compress() function is called.

#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test(TrilinosWrappers::MPI::Vector &v)
{
  // set only certain elements of the
  // vector. have a bit pattern of where we
  // actually wrote elements to
  std::vector<bool> pattern(v.size(), false);
  bool              flag = false;
  for (unsigned int i = 0; i < v.size(); i += 1 + i)
    {
      if (flag == true)
        v(i) += i;
      else
        v(i) = i;
      flag = !flag;

      pattern[i] = true;
    }

  // check that they are ok, and this time
  // all of them
  for (unsigned int i = 0; i < v.size(); ++i)
    AssertThrow((((pattern[i] == true) && (v(i) == i)) ||
                 ((pattern[i] == false) && (v(i) == 0))),
                ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());


  try
    {
      {
        TrilinosWrappers::MPI::Vector v;
        v.reinit(complete_index_set(100), MPI_COMM_WORLD);
        test(v);
      }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}

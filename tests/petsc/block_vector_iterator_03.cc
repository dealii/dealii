// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// this test is an adaptation of lac/block_vector_iterator for PETSc block
// vectors

#include <deal.II/lac/petsc_block_vector.h>

#include <algorithm>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>

#include "../tests.h"

template <typename number>
bool
operator==(const PETScWrappers::MPI::BlockVector &v1,
           const PETScWrappers::MPI::BlockVector &v2)
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

  // Check 1: initialization via
  // iterators
  if (true)
    {
      PETScWrappers::MPI::BlockVector v1(ivector, MPI_COMM_WORLD, ivector);
      PETScWrappers::MPI::BlockVector v2(ivector, MPI_COMM_WORLD, ivector);

      // initialize first vector with
      // simple loop
      for (unsigned int i = 0; i < v1.size(); ++i)
        v1(i) = i;
      // initialize other vector
      // through iterators
      PETScWrappers::MPI::BlockVector::iterator p2 = v2.begin();
      for (unsigned int i = 0; i < v1.size(); ++i, ++p2)
        *p2 = i;
      Assert(p2 == v2.end(), ExcInternalError());

      // check that the two vectors are equal
      deallog << "Check 1: " << (v1 == v2 ? "true" : "false") << std::endl;

      // print vectors
      v1.print(deallog.get_file_stream(), 10, true, false);

      // Extract the PETSc VECNEST and use print from PETScWrappers::VectorBase
      PETScWrappers::VectorBase v1b(v1.petsc_vector());
      v1b.print(deallog.get_file_stream(), 10, true, false);
    };

  // Check 2: loop forward and back
  // and check that things are the
  // same
  if (true)
    {
      PETScWrappers::MPI::BlockVector v1(ivector, MPI_COMM_WORLD, ivector);
      // initialize first vector with
      // simple loop
      for (unsigned int i = 0; i < v1.size(); ++i)
        v1(i) = i;

      PETScWrappers::MPI::BlockVector::iterator p1 = v1.begin();
      for (unsigned int i = 0; i < v1.size(); ++i, ++p1)
        AssertThrow(*p1 == i, ExcInternalError());

      Assert(p1 == v1.end(), ExcInternalError());

      // move back into allowable
      // region
      --p1;

      // check backwards
      for (unsigned int i = 0; i < v1.size(); ++i, --p1)
        AssertThrow(*p1 == v1.size() - i - 1, ExcInternalError());

      // if we came thus far,
      // everything is alright
      deallog << "Check 2: true" << std::endl;
    };


  // Check 3: same, but this time
  // with const iterators
  if (true)
    {
      PETScWrappers::MPI::BlockVector v1(ivector, MPI_COMM_WORLD, ivector);
      // initialize first vector with
      // simple loop
      for (unsigned int i = 0; i < v1.size(); ++i)
        v1(i) = i;

      PETScWrappers::MPI::BlockVector::const_iterator p1 = v1.begin();
      for (unsigned int i = 0; i < v1.size(); ++i, ++p1)
        AssertThrow(*p1 == i, ExcInternalError());

      Assert(p1 == v1.end(), ExcInternalError());

      // move back into allowable
      // region
      --p1;

      // check backwards
      for (unsigned int i = 0; i < v1.size(); ++i, --p1)
        {
          const double val = *p1;
          const double ref = v1.size() - i - 1;
          AssertThrow(val == ref, ExcInternalError());
        };

      // if we came thus far,
      // everything is alright
      deallog << "Check 3: true" << std::endl;
    };

  // Checks 4-13: use some standard
  // algorithms
  if (true)
    {
      PETScWrappers::MPI::BlockVector v1(ivector, MPI_COMM_WORLD, ivector);
      // initialize first vector with
      // simple loop
      for (unsigned int i = 0; i < v1.size(); ++i)
        v1(i) = i;

      // check std::distance
      // algorithm
      deallog << "Check 4: "
              << (std::distance(v1.begin(), v1.end()) ==
                      static_cast<signed int>(v1.size()) ?
                    "true" :
                    "false")
              << std::endl;

      // check std::copy
      PETScWrappers::MPI::BlockVector v2(ivector, MPI_COMM_WORLD, ivector);
      std::copy(v1.begin(), v1.end(), v2.begin());
      deallog << "Check 5: " << (v1 == v2 ? "true" : "false") << std::endl;

      // check std::transform
      std::transform(v1.begin(),
                     v1.end(),
                     v2.begin(),
                     std::bind(std::multiplies<double>(),
                               std::placeholders::_1,
                               2.0));
      v2.compress(VectorOperation::insert);
      v2 *= 1. / 2.;
      deallog << "Check 6: " << (v1 == v2 ? "true" : "false") << std::endl;


      // check operators +/-, +=/-=
      deallog << "Check 7: "
              << (std::distance(v1.begin(), v1.begin() + 3) == 3 ? "true" :
                                                                   "false")
              << std::endl;
      deallog << "Check 8: "
              << (std::distance(v1.end() - 6, v1.end()) == 6 ? "true" : "false")
              << std::endl;
      deallog << "Check 9: "
              << (std::distance(v1.begin(), v1.end()) == (signed)v1.size() ?
                    "true" :
                    "false")
              << std::endl;
      deallog << "Check 10: "
              << (std::distance(v1.begin(), (v1.begin() += 7)) == 7 ? "true" :
                                                                      "false")
              << std::endl;
      deallog << "Check 11: "
              << (std::distance((v1.end() -= 4), v1.end()) == 4 ? "true" :
                                                                  "false")
              << std::endl;

      // check advance
      PETScWrappers::MPI::BlockVector::iterator p2 = v1.begin();
      std::advance(p2, v1.size());
      deallog << "Check 12: " << (p2 == v1.end() ? "true" : "false")
              << std::endl;

      PETScWrappers::MPI::BlockVector::const_iterator p3 = v1.begin();
      std::advance(p3, v1.size());
      deallog << "Check 13: " << (p3 == v1.end() ? "true" : "false")
              << std::endl;
    };

  // Check 14: operator[]
  if (true)
    {
      PETScWrappers::MPI::BlockVector v1(ivector, MPI_COMM_WORLD, ivector);
      for (unsigned int i = 0; i < v1.size(); ++i)
        v1(i) = i;

      for (unsigned int i = 0; i < v1.size(); ++i)
        {
          const PETScWrappers::MPI::BlockVector::iterator p = (v1.begin() + i);
          for (unsigned int j = 0; j < v1.size(); ++j)
            AssertThrow(p[(signed)j - (signed)i] == j, ExcInternalError());
        };

      // if we came thus far,
      // everything is alright
      deallog << "Check 14: true" << std::endl;
    };
}



int
main(int argc, char **argv)
{
  std::ofstream logfile("output");
  logfile.setf(std::ios::fixed);
  logfile.precision(3);
  deallog.attach(logfile);

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        test();
      }
    }
  catch (const std::exception &e)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << e.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      // abort
      return 2;
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
      // abort
      return 3;
    };


  return 0;
}

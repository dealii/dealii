// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/utilities.h>

#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"

// Check LinearAlgebra::TpetraWrappers::Vector assignment and import

template <typename Number>
void
test()
{
  IndexSet     parallel_partitioner_1(10);
  IndexSet     parallel_partitioner_2(10);
  unsigned int rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (rank == 0)
    {
      parallel_partitioner_1.add_range(0, 5);
      parallel_partitioner_2.add_range(0, 3);
    }
  else
    {
      parallel_partitioner_1.add_range(5, 10);
      parallel_partitioner_2.add_range(3, 10);
    }
  parallel_partitioner_1.compress();
  parallel_partitioner_2.compress();
  LinearAlgebra::TpetraWrappers::Vector<Number> a;
  LinearAlgebra::TpetraWrappers::Vector<Number> b(parallel_partitioner_1,
                                                  MPI_COMM_WORLD);
  LinearAlgebra::TpetraWrappers::Vector<Number> c(b);

  AssertThrow(a.size() == 0, ExcMessage("Vector has the wrong size."));
  AssertThrow(b.size() == 10, ExcMessage("Vector has the wrong size."));
  AssertThrow(c.size() == 10, ExcMessage("Vector has the wrong size."));

  a.reinit(parallel_partitioner_2, MPI_COMM_WORLD);
  AssertThrow(a.size() == 10, ExcMessage("Vector has the wrong size."));

  AssertThrow(parallel_partitioner_1 == b.locally_owned_elements(),
              ExcMessage("IndexSet has been modified."));
  AssertThrow(parallel_partitioner_2 == a.locally_owned_elements(),
              ExcMessage("IndexSet has been modified."));

  IndexSet read_write_index_set(10);
  if (rank == 0)
    read_write_index_set.add_range(0, 6);
  else
    read_write_index_set.add_range(4, 10);
  read_write_index_set.compress();

  LinearAlgebra::ReadWriteVector<Number> read_write_1(read_write_index_set);
  LinearAlgebra::ReadWriteVector<Number> read_write_2(read_write_index_set);
  LinearAlgebra::ReadWriteVector<Number> read_write_3(read_write_index_set);
  if (rank == 0)
    {
      for (unsigned int i = 0; i < 6; ++i)
        {
          read_write_1[i] = i;
          read_write_2[i] = 5. + i;
        }
    }
  else
    {
      for (unsigned int i = 4; i < 10; ++i)
        {
          read_write_1[i] = i;
          read_write_2[i] = 5. + i;
        }
    }

  a.import_elements(read_write_2, VectorOperation::insert);
  AssertThrow(a.size() == 10, ExcMessage("Vector has the wrong size."));

  read_write_3.import_elements(a, VectorOperation::insert);
  if (rank == 0)
    {
      for (unsigned int i = 0; i < 6; ++i)
        {
          AssertThrow(read_write_2[i] == read_write_3[i],
                      ExcMessage("Vector a has been modified."));
        }
    }
  else
    {
      for (unsigned int i = 4; i < 10; ++i)
        AssertThrow(read_write_2[i] == read_write_3[i],
                    ExcMessage("Vector a has been modified."));
    }

  b.import_elements(read_write_1, VectorOperation::insert);
  AssertThrow(b.size() == 10, ExcMessage("Vector has the wrong size."));

  read_write_3.import_elements(b, VectorOperation::insert);
  if (rank == 0)
    {
      for (unsigned int i = 0; i < 6; ++i)
        AssertThrow(read_write_1[i] == read_write_3[i],
                    ExcMessage("Vector b has been modified."));
    }
  else
    {
      for (unsigned int i = 4; i < 10; ++i)
        AssertThrow(read_write_1[i] == read_write_3[i],
                    ExcMessage("Vector b has been modified."));
    }

  c.import_elements(read_write_2, VectorOperation::insert);
  AssertThrow(c.size() == 10, ExcMessage("Vector has the wrong size."));

  read_write_3.import_elements(c, VectorOperation::insert);
  if (rank == 0)
    {
      for (unsigned int i = 0; i < 6; ++i)
        AssertThrow(read_write_2[i] == read_write_3[i],
                    ExcMessage("Vector c has been modified."));
    }
  else
    {
      for (unsigned int i = 4; i < 10; ++i)
        AssertThrow(read_write_2[i] == read_write_3[i],
                    ExcMessage("Vector c has been modified."));
    }

  a *= 2;
  read_write_3.import_elements(a, VectorOperation::insert);
  if (rank == 0)
    {
      for (unsigned int i = 0; i < 6; ++i)
        AssertThrow(Number(2.) * read_write_2[i] == read_write_3[i],
                    ExcMessage("Problem in operator *=."));
    }
  else
    {
      for (unsigned int i = 4; i < 10; ++i)
        AssertThrow(Number(2.) * read_write_2[i] == read_write_3[i],
                    ExcMessage("Problem in operator *=."));
    }

  c /= 2.;
  read_write_3.import_elements(c, VectorOperation::insert);
  if (rank == 0)
    {
      for (unsigned int i = 0; i < 6; ++i)
        AssertThrow(Number(0.5) * read_write_2[i] == read_write_3[i],
                    ExcMessage("Problem in operator /=."));
    }
  else
    {
      for (unsigned int i = 4; i < 10; ++i)
        AssertThrow(Number(0.5) * read_write_2[i] == read_write_3[i],
                    ExcMessage("Problem in operator /=."));
    }

  b += a;
  read_write_3.import_elements(b, VectorOperation::insert);
  if (rank == 0)
    {
      for (unsigned int i = 0; i < 6; ++i)
        AssertThrow(Number(2.) * read_write_2[i] + read_write_1[i] ==
                      read_write_3[i],
                    ExcMessage("Problem in operator +=."));
    }
  else
    {
      for (unsigned int i = 4; i < 10; ++i)
        AssertThrow(Number(2.) * read_write_2[i] + read_write_1[i] ==
                      read_write_3[i],
                    ExcMessage("Problem in operator +=."));
    }

  b -= c;
  read_write_3.import_elements(b, VectorOperation::insert);
  if (rank == 0)
    {
      for (unsigned int i = 0; i < 6; ++i)
        AssertThrow(Number(1.5) * read_write_2[i] + read_write_1[i] ==
                      read_write_3[i],
                    ExcMessage("Problem in operator -=."));
    }
  else
    {
      for (unsigned int i = 4; i < 10; ++i)
        AssertThrow(Number(1.5) * read_write_2[i] + read_write_1[i] ==
                      read_write_3[i],
                    ExcMessage("Problem in operator -=."));
    }

  b.import_elements(read_write_1, VectorOperation::insert);
  c.import_elements(read_write_1, VectorOperation::insert);
  const Number val = b * c;
  AssertThrow(val == Number(285.), ExcMessage("Problem in operator *."));
}


int
main(int argc, char **argv)
{
  initlog();
  deallog.depth_console(0);

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  test<double>();

  deallog << "OK" << std::endl;

  return 0;
}

// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Check that a LinearAlgebra::distributed::Vector<double> in
// MemorySpace::Default can be assigned to a Vector<float> and vice versa.

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/read_write_vector.h>

#include "../tests.h"


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  initlog();

  const unsigned int size = 4;

  // fill a double vector on the device with values that are exact in float
  LinearAlgebra::ReadWriteVector<double> rw_double(size);
  for (unsigned int i = 0; i < size; ++i)
    rw_double(i) = 0.5 + 0.25 * i;

  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> v_double(
    size);
  v_double.import_elements(rw_double, VectorOperation::insert);

  // double -> float
  LinearAlgebra::distributed::Vector<float, MemorySpace::Default> v_float;
  v_float = v_double;

  LinearAlgebra::ReadWriteVector<float> rw_float(size);
  rw_float.import_elements(v_float, VectorOperation::insert);
  for (unsigned int i = 0; i < size; ++i)
    deallog << "float " << i << ": " << rw_float(i) << std::endl;

  // float -> double
  v_float *= 2.f;
  v_double = v_float;

  rw_double.import_elements(v_double, VectorOperation::insert);
  for (unsigned int i = 0; i < size; ++i)
    deallog << "double " << i << ": " << rw_double(i) << std::endl;
}

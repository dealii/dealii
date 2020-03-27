
// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2018 by the deal.II authors
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


// test the identity_operator LinearOperator and its payload initialization


#include <deal.II/lac/linear_operator_tools.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

#include "../testmatrix.h"


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  initlog();

  // Need to get the linear operator to do some nonsense work to verify
  // that it can indeed be used as expected
  const unsigned int size = 32;
  unsigned int       dim  = (size - 1) * (size - 1);

  // Make matrix
  FDMatrix testproblem(size, size);

  {
    SparsityPattern sp(dim, dim);
    testproblem.five_point_structure(sp);
    sp.compress();

    SparseMatrix<double> A;
    A.reinit(sp);
    testproblem.five_point(A);

    Vector<double> f(dim);
    Vector<double> u(dim);

    const auto lo_A = linear_operator<Vector<double>>(A);
    const auto lo_id =
      identity_operator<Vector<double>>(lo_A.reinit_range_vector);
    const auto lo_A_plus_id = lo_A + lo_id;

    u = lo_A_plus_id * f;
  }

  {
    DynamicSparsityPattern csp(dim, dim);
    testproblem.five_point_structure(csp);

    TrilinosWrappers::SparseMatrix A;
    A.reinit(csp);
    testproblem.five_point(A);

    TrilinosWrappers::MPI::Vector f;
    f.reinit(complete_index_set(dim));
    TrilinosWrappers::MPI::Vector u;
    u.reinit(complete_index_set(dim));

    A.compress(VectorOperation::insert);
    f.compress(VectorOperation::insert);
    u.compress(VectorOperation::insert);

    const auto lo_A    = linear_operator<TrilinosWrappers::MPI::Vector>(A);
    const auto lo_id_1 = identity_operator<TrilinosWrappers::MPI::Vector>(
      lo_A.reinit_range_vector);
    const auto lo_id_2 = identity_operator<TrilinosWrappers::MPI::Vector>(lo_A);
    const auto lo_A_plus_id_1 = lo_A + lo_id_1; // Not a good idea. See below.
    const auto lo_A_plus_id_2 = lo_A + lo_id_2;

    // u = lo_A_plus_id_1 * f; // This is not allowed -- different payloads
    u = lo_A_plus_id_2 * f;
  }

  deallog << "OK" << std::endl;

  return 0;
}

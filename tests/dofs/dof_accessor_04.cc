// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2023 by the deal.II authors
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

// Test DoFCellAccessor::distribute_local_to_global(
// local_matrix, local_vector, global_matrix, global_vector)

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

int
main()
{
  initlog();

  // create the triangulation
  Triangulation<2> tria;
  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global(1);

  // distribute dofs
  DoFHandler<2> dof_handler(tria);

  FE_Q<2> fe(1);
  dof_handler.distribute_dofs(fe);

  // create the global matrix and global vector
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  SparsityPattern sp;
  sp.copy_from(dsp);

  SparseMatrix<double> global_matrix(sp);
  Vector<double>       global_vector(dof_handler.n_dofs());

  // create local matrix and local vector
  FullMatrix<double> local_matrix(fe.dofs_per_cell, fe.dofs_per_cell);
  Vector<double>     local_vector(fe.dofs_per_cell);
  for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
    {
      local_vector(i) = i + 1;
      for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
        local_matrix(i, j) = i * fe.dofs_per_cell + j + 1;
    }

  // call function distribute_local_to_global()
  auto cell = dof_handler.begin_active();
  cell->distribute_local_to_global(local_matrix,
                                   local_vector,
                                   global_matrix,
                                   global_vector);

  // check if each entry of the local matrix/vector has
  // found its position in the global matrix/vector
  bool passed = true;
  for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
    {
      if (global_vector(i) != local_vector(i))
        passed = false;

      for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
        if (global_matrix(i, j) != local_matrix(i, j))
          {
            passed = false;
            break;
          }

      if (passed == false)
        break;
    }

  if (passed)
    deallog << "PASSED" << std::endl;
  else
    deallog << "NOT PASSED" << std::endl;

  return 0;
}

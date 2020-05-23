// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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


// Similar to dof_tools_04 but for a parallel::distributed::Triangulation
// instead of a (serial) Triangulation:
// check
//   DoFTools::extract_hanging_node_constraints
// using a slightly different refinement and less different FiniteElements

#include <deal.II/base/logstream.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_iterator.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <string>

#include "../tests.h"



template <int dim>
void
check_this(const DoFHandler<dim> &dof_handler)
{
  const types::global_dof_index n_dofs = dof_handler.n_dofs();

  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  const IndexSet is_hanging_node_constrained =
    DoFTools::extract_hanging_node_dofs(dof_handler);

  AffineConstraints<double> constraints(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  for (const auto &dof : locally_relevant_dofs)
    if (is_hanging_node_constrained.is_element(dof))
      AssertThrow(constraints.is_constrained(dof), ExcInternalError());

  AssertThrow(is_hanging_node_constrained.n_elements() ==
                constraints.n_constraints(),
              ExcInternalError());

  deallog << "OK" << std::endl;
}



template <int dim>
void
check(const FiniteElement<dim> &fe, const std::string &name)
{
  deallog << "Checking " << name << " in " << dim << "d:" << std::endl;

  // create tria and dofhandler
  // objects. set different boundary
  // and sub-domain ids
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global(1);
  for (unsigned int ref = 0; ref < 2; ++ref)
    {
      for (auto &cell : tria.active_cell_iterators())
        if (cell->is_locally_owned() && cell->center()(0) < .5 &&
            cell->center()(1) < .5)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  check_this(dof_handler);
}

#define CHECK(EL, deg, dim) \
  {                         \
    FE_##EL<dim> EL(deg);   \
    check(EL, #EL #deg);    \
  }

#define CHECK_ALL(EL, deg) \
  {                        \
    CHECK(EL, deg, 2);     \
    CHECK(EL, deg, 3);     \
  }

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();
  CHECK_ALL(Q, 1);

  CHECK_ALL(Q, 2);

  return 0;
}

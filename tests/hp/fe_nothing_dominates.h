// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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



// Check for continuity requirements of two neighboring finite elements
// and project a given function subject to these constraints.


#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/vector_tools_project.h>

#include "../tests.h"

#include "../test_grids.h"


template <int dim>
void
project(const hp::FECollection<dim> &fe_collection,
        const hp::QCollection<dim> & q_collection,
        const Function<dim> &        function)
{
#ifdef DEBUG
  Assert(fe_collection.size() == 2, ExcInternalError());
  Assert(q_collection.size() == 2, ExcInternalError());
  for (unsigned int f = 0; f < fe_collection.size(); ++f)
    Assert(fe_collection[f].n_components() == function.n_components,
           ExcInternalError());
#endif

  // setup
  // +---+---+
  // | 0 | 1 |
  // +---+---+
  Triangulation<dim> tria;
  TestGrids::hyper_line(tria, 2);

  DoFHandler<dim> dofh(tria);
  (++(dofh.begin_active()))->set_active_fe_index(1);
  dofh.distribute_dofs(fe_collection);

  // make constraints
  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dofh, constraints);
  constraints.close();
  deallog << "constraints:" << std::endl;
  constraints.print(deallog.get_file_stream());

  // project function with constraints
  Vector<double> solution(dofh.n_dofs());
  VectorTools::project(dofh, constraints, q_collection, function, solution);

  // verify output
  deallog << "dof values:" << std::endl;
  Vector<double> cell_values;
  for (const auto &cell : dofh.active_cell_iterators())
    {
      cell_values.reinit(cell->get_fe().n_dofs_per_cell());
      cell->get_dof_values(solution, cell_values);

      deallog << " cell " << cell->active_cell_index() << ":";
      for (const auto &value : cell_values)
        deallog << " " << value;
      deallog << std::endl;
    }

  deallog << "OK" << std::endl;
}

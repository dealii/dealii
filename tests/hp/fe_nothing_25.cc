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



// Check for continuity requirements for neighboring
// (FE_QxFE_Nothing) and (FE_NothingxFE_Q) elements.
// The twist: only one FE_Nothing element dominates.


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/numerics/vector_tools_project.h>

#include "../tests.h"

#include "../test_grids.h"



template <int dim>
void
test()
{
  // setup
  // +-----+-----+
  // | FEQ | FEN |
  // | FEN | FEQ |
  // +-----+-----+
  Triangulation<dim> tria;
  TestGrids::hyper_line(tria, 2);

  hp::FECollection<dim> fes;
  fes.push_back(FESystem<dim>(FE_Nothing<dim>(1, true), FE_Q<dim>(1)));
  fes.push_back(FESystem<dim>(FE_Q<dim>(1), FE_Nothing<dim>(1, false)));

  DoFHandler<dim> dofh(tria);
  dofh.begin_active()->set_active_fe_index(1);
  dofh.distribute_dofs(fes);

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dofh, constraints);
  constraints.close();
  deallog << "constraints:" << std::endl;
  constraints.print(deallog.get_file_stream());

  // init constant solution with all constraints
  Vector<double> solution(dofh.n_dofs());
  VectorTools::project(dofh,
                       constraints,
                       hp::QCollection<dim>(QGauss<dim>(2), QGauss<dim>(2)),
                       Functions::ConstantFunction<dim>(1., 2),
                       solution);

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


int
main()
{
  initlog();

  deallog.push("1d");
  test<1>();
  deallog.pop();
  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();
}

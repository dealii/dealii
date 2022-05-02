// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2020 by the deal.II authors
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

// test for bug #82
// (http://code.google.com/p/dealii/issues/detail?id=82) which
// demonstrates that under some circumstances we create cycles in
// constraints

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/matrix_out.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <sstream>
#include <vector>

#include "../tests.h"


class Deal2PeriodicBug
{
public:
  Deal2PeriodicBug();
  void
  run();

private:
  void
  makeGrid();
  void
  make_periodicity_constraints();
  void
  setup_system();

  Triangulation<2>          triangulation;
  FE_Q<2>                   fe;
  DoFHandler<2>             dof_handler;
  AffineConstraints<double> constraints;
};

Deal2PeriodicBug::Deal2PeriodicBug()
  : fe(2)
  , dof_handler(triangulation)
{}


void
Deal2PeriodicBug::run()
{
  makeGrid();
  setup_system();
}

void
Deal2PeriodicBug::make_periodicity_constraints()
{
  std::vector<bool> mask(1);
  mask[0] = true;
  ComponentMask cmask(mask);
  // Here we use the DoFTools function to place periodic constraints in the
  // constraints matrix We have set the boundary index for the left face of the
  // boundary to 0 and the right face index to 2 (Hence the 2nd and 3rd
  // arguments) The direction indicator is to specify that the DOFs are to be
  // matched in the y-direction (0 in the 4th argument since they are allowed to
  // differ in the x-direction, it compares each of the coordinates not
  // specified by the direction integer which here is 0) Since we wanted to only
  // specify the director and electric components for periodicity we use a
  // component mask as well
  DoFTools::make_periodicity_constraints(
    dof_handler, 2, 0, 0, constraints, cmask);
}

void
Deal2PeriodicBug::setup_system()
{
  dof_handler.distribute_dofs(fe);

  constraints.clear();
  deallog << "Making Constraint Matrix..." << std::endl;
  make_periodicity_constraints();
  deallog << "Constraint Matrix Complete" << std::endl;

  constraints.print(deallog.get_file_stream());

  constraints.close();
}

void
Deal2PeriodicBug::makeGrid()
{
  deallog << "Constructing the grid..." << std::endl;
  const Point<2> p1(0, 0), p2(1, 1);
  GridGenerator::hyper_rectangle(triangulation, p1, p2);
  triangulation.begin_active()->face(2)->set_boundary_id(1);
  triangulation.begin_active()->face(3)->set_boundary_id(1);
  triangulation.begin_active()->face(0)->set_boundary_id(0);
  triangulation.begin_active()->face(1)->set_boundary_id(2);
  triangulation.refine_global(1);

  Triangulation<2>::active_cell_iterator cell = triangulation.begin_active();
  (++(++cell))->set_refine_flag();
  triangulation.prepare_coarsening_and_refinement();
  triangulation.execute_coarsening_and_refinement();

  deallog << "Number of active cells: " << triangulation.n_active_cells()
          << std::endl;
  GridOut grid_out;
  grid_out.write_eps(triangulation, deallog.get_file_stream());
  deallog << "Grid construction complete..." << std::endl;
}



int
main()
{
  initlog();

  Deal2PeriodicBug Deal2Bug;
  Deal2Bug.run();
  return 0;
}

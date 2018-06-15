// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2017 by the deal.II authors
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



// no normal flux constraints on a hyper cube for all faces this caused
// ExcMessage (\"Cycle in constraints detected!\")" in 3d with a higher order
// mapping.  to make things even weirder, mappings of order <4 work.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools.templates.h>

#include "../tests.h"



template <int dim>
void
test_hyper_cube()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_rectangle(tr, Point<dim>(), Point<dim>(1, 1, 1), true);

  FESystem<dim> fe(FE_Q<dim>(2), dim);

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  deallog << "FE=" << fe.get_name() << std::endl;

  std::set<types::boundary_id> boundary_ids;
  boundary_ids.insert(1);
  boundary_ids.insert(3);

  ConstraintMatrix    cm;
  const MappingQ<dim> mapping(4);
  VectorTools::compute_no_normal_flux_constraints(
    dof, 0, boundary_ids, cm, mapping);
  cm.close();

  cm.print(deallog.get_file_stream());
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(2);
  deallog << std::fixed;
  deallog.attach(logfile);

  test_hyper_cube<3>();
}

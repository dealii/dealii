// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



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

  AffineConstraints<double> cm;
  const MappingQ<dim>       mapping(4);
  VectorTools::compute_no_normal_flux_constraints(
    dof, 0, boundary_ids, cm, mapping);
  cm.close();

  cm.print(deallog.get_file_stream());
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);
  deallog << std::fixed;

  test_hyper_cube<3>();
}

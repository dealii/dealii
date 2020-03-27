// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2018 by the deal.II authors
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

// bug report from mailing list from 11/15/2013 (simplified). no_normal_flux
// throws an ExcInternalError when handing it an FE with less than dim
// components. This is now fixed (throws ExcMessage).

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -1.0, 1.0);
  triangulation.begin_active()->face(1)->set_all_boundary_ids(1);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<double>    constraints;
  std::set<types::boundary_id> no_normal_flux_boundaries;
  no_normal_flux_boundaries.insert(1);
  deal_II_exceptions::disable_abort_on_exception();
  try
    {
      VectorTools::compute_no_normal_flux_constraints(dof_handler,
                                                      0,
                                                      no_normal_flux_boundaries,
                                                      constraints);
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }

  constraints.close();
  constraints.print(deallog.get_file_stream());
}


int
main()
{
  initlog();

  test<3>();
}

// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Similar to no_flux_14 for a higher order element including a high-order
// mapping. The purpose here is to test that the result of
// Manifold::normal_vector function inside of
// VectorTools::compute_no_normal_flux_constraints is accurate for higher
// order mappings.


#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim>
void
check()
{
  SphericalManifold<dim> spherical;
  Triangulation<dim>     tria;
  GridGenerator::hyper_shell(tria, Point<dim>(), 0.5, 1., 96, true);
  tria.set_all_manifold_ids(0);
  tria.set_manifold(0, spherical);

  AffineConstraints<double> cm;
  MappingQ<dim>             mapping(4);

  FESystem<dim>   fe(FE_Q<dim>(2), dim);
  DoFHandler<dim> dofh(tria);

  dofh.distribute_dofs(fe);

  const std::set<types::boundary_id> no_normal_flux_boundaries = {0, 1};
  VectorTools::compute_no_normal_flux_constraints(
    dofh, 0, no_normal_flux_boundaries, cm, mapping);

  cm.print(deallog.get_file_stream());
}



int
main()
{
  initlog();
  deallog.get_file_stream().precision(8);

  check<3>();
}

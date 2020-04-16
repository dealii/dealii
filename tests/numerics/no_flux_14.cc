// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2018 by the deal.II authors
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



// Copy of no_flux_09 for a higher order element. The purpose here is to test
// the Manifold::normal_vector function inside of
// VectorTools::compute_no_normal_flux_constraints for points other than the
// vertices.


#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/dof_handler.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
void
check()
{
  Triangulation<dim> tr;
  GridGenerator::quarter_hyper_shell(tr, Point<dim>(), 0.5, 1.0, 3, true);
  tr.reset_manifold(0);

  AffineConstraints<double> cm;
  MappingQ<dim>             mapping(1);

  FESystem<dim>   fe(FE_Q<dim>(3), dim);
  DoFHandler<dim> dofh(tr);

  dofh.distribute_dofs(fe);

  std::set<types::boundary_id> no_normal_flux_boundaries;
  no_normal_flux_boundaries.insert(1);
  //  no_normal_flux_boundaries.insert (2); // not required for the crash for
  //  now, please test with it later!
  no_normal_flux_boundaries.insert(3);
  no_normal_flux_boundaries.insert(4);
  VectorTools::compute_no_normal_flux_constraints(
    dofh, 0, no_normal_flux_boundaries, cm, mapping);

  cm.print(deallog.get_file_stream());
}



int
main()
{
  initlog();
  deallog.get_file_stream().precision(4);
  deallog.get_file_stream().setf(std::ios::fixed);

  check<3>();
}

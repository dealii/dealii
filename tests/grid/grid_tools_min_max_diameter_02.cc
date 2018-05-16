// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// test GridTools::minimal_cell_diameter and GridTools::maximal_cell_diameter with a
// MappingFEField

#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/fe/mapping_cartesian.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q_eulerian.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_bernstein.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/numerics/vector_tools.h>

#include <vector>
#include <string>
#include <sstream>

#define PRECISION 2

template <int dim, int spacedim>
void
test(const unsigned int degree)
{
  deallog << "dim = " << dim << ", spacedim = " << spacedim << std::endl;
  deallog << "degree = " << degree << std::endl;

  Triangulation<dim,spacedim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global (1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FE_Q<dim,spacedim> fe(degree);
  FESystem<dim,spacedim> fe_sys(fe, spacedim);

  DoFHandler<dim,spacedim> dof_sys(tria);
  dof_sys.distribute_dofs(fe_sys);

  Vector<double> euler;
  euler.reinit(dof_sys.n_dofs());
  const ComponentMask mask(spacedim, true);

  VectorTools::get_position_vector(dof_sys, euler, mask);

  MappingFEField<dim,spacedim> map_fe(dof_sys, euler, mask);

  // Make the grid bigger by a factor two.
  euler *= 2.0;

  deallog << "Min diameter        : " << GridTools::minimal_cell_diameter(tria) << std::endl
          << "Max diameter        : " << GridTools::maximal_cell_diameter(tria) << std::endl
          << "Min mapped diameter : " << GridTools::minimal_cell_diameter(tria,map_fe) << std::endl
          << "Max mapped diameter : " << GridTools::maximal_cell_diameter(tria,map_fe) << std::endl;
}

int
main()
{
  initlog();

  for (unsigned int d=1; d<4; ++d)
    {
      test<2,2>(d);
      test<2,3>(d);
      test<3,3>(d);
    }
}

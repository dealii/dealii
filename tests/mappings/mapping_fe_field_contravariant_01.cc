// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2019 by the deal.II authors
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


// Like mapping_fe_field_01, but add an update flag that led to an
// assertion that checked the size of the wrong field and consequently
// accidentally aborted an otherwise perfectly valid program

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_bernstein.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_cartesian.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <sstream>
#include <string>
#include <vector>

#include "../tests.h"

#define PRECISION 2

template <int dim, int spacedim>
void
test(const unsigned int degree)
{
  deallog << "dim = " << dim << ", spacedim = " << spacedim << std::endl;
  deallog << "degree = " << degree << std::endl;

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global(2);

  FE_Bernstein<dim, spacedim> fe(degree);
  // FE_Q<dim> fe(degree);
  FESystem<dim, spacedim> fe_sys(fe, spacedim);

  // DoFHandler<dim> dof(tria);
  DoFHandler<dim, spacedim> dof_sys(tria);
  // dof.distribute_dofs(fe);
  dof_sys.distribute_dofs(fe_sys);

  Vector<double> euler;
  euler.reinit(dof_sys.n_dofs());
  const ComponentMask mask(spacedim, true);

  VectorTools::get_position_vector(dof_sys, euler, mask);
  MappingFEField<dim, spacedim> map_fe(dof_sys, euler, mask);

  QIterated<dim> quadrature_formula(QTrapez<1>(), fe.degree);

  FEValues<dim, spacedim> fe_values(map_fe,
                                    fe_sys,
                                    quadrature_formula,
                                    update_values | update_JxW_values |
                                      update_volume_elements);

  typename DoFHandler<dim, spacedim>::active_cell_iterator cell =
                                                             dof_sys
                                                               .begin_active(),
                                                           endc = dof_sys.end();

  for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      deallog << "Cell " << cell << ": OK" << std::endl;
    }
}

int
main()
{
  initlog();

  for (unsigned int d = 1; d < 4; ++d)
    {
      test<2, 2>(d);
      test<2, 3>(d);
      test<3, 3>(d);
    }
}

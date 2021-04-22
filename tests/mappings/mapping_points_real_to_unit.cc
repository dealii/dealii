// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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



// on a test case similar to mapping_real_to_unit_q4_curved, check the
// implementation of the many-point interface
// Mapping::transform_points_real_to_unit_cell for both a MappingQGeneric and
// MappingFEField

#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim, int spacedim>
void
test_real_to_unit_cell(const Mapping<dim, spacedim> &mapping)
{
  // define a boundary that fits the the vertices of the hyper cube we're
  // going to create below
  SphericalManifold<dim, spacedim> boundary;

  Triangulation<dim, spacedim> triangulation;
  GridGenerator::hyper_cube(triangulation, -1, 1);

  // set the boundary indicator for one face of the single cell
  triangulation.set_manifold(1, boundary);
  triangulation.begin_active()->face(0)->set_boundary_id(1);

  const unsigned int      n_points = 5;
  std::vector<Point<dim>> unit_points(Utilities::pow(n_points, dim));

  switch (dim)
    {
      case 1:
        for (unsigned int x = 0; x < n_points; ++x)
          unit_points[x][0] = static_cast<double>(x) / n_points;
        break;

      case 2:
        for (unsigned int x = 0; x < n_points; ++x)
          for (unsigned int y = 0; y < n_points; ++y)
            {
              unit_points[y * n_points + x][0] =
                static_cast<double>(x) / n_points;
              unit_points[y * n_points + x][1] =
                static_cast<double>(y) / n_points;
            }
        break;

      case 3:
        for (unsigned int x = 0; x < n_points; ++x)
          for (unsigned int y = 0; y < n_points; ++y)
            for (unsigned int z = 0; z < n_points; ++z)
              {
                unit_points[z * n_points * n_points + y * n_points + x][0] =
                  static_cast<double>(x) / n_points;
                unit_points[z * n_points * n_points + y * n_points + x][1] =
                  static_cast<double>(y) / n_points;
                unit_points[z * n_points * n_points + y * n_points + x][2] =
                  static_cast<double>(z) / n_points;
              }
        break;
    }

  typename Triangulation<dim, spacedim>::active_cell_iterator cell =
    triangulation.begin_active();
  std::vector<Point<spacedim>> real_points(unit_points.size());
  for (unsigned int i = 0; i < unit_points.size(); ++i)
    real_points[i] = mapping.transform_unit_to_real_cell(cell, unit_points[i]);
  std::vector<Point<dim>> new_points(unit_points.size());
  mapping.transform_points_real_to_unit_cell(cell, real_points, new_points);
  for (unsigned int i = 0; i < unit_points.size(); ++i)
    {
      // for each of the points, verify that applying the forward map and
      // then pull back get the same point again
      AssertThrow(unit_points[i].distance(new_points[i]) < 1e-10,
                  ExcInternalError());
    }
  deallog << "OK" << std::endl;
}



template <int dim, int spacedim>
void
test()
{
  deallog << "dim=" << dim << ", spacedim=" << spacedim << std::endl;
  deallog << "MappingQ(1): ";
  test_real_to_unit_cell(MappingQGeneric<dim, spacedim>(1));
  deallog << "MappingQ(4): ";
  test_real_to_unit_cell(MappingQGeneric<dim, spacedim>(4));

  deallog << "MappingFEField(FESystem(FE_Q(4))): ";
  Triangulation<dim, spacedim> triangulation;
  GridGenerator::hyper_cube(triangulation, -1, 1);

  FESystem<dim, spacedim>   fe(FE_Q<dim, spacedim>(4), spacedim);
  DoFHandler<dim, spacedim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  const ComponentMask mask(spacedim, true);
  Vector<double>      location_vector(dof_handler.n_dofs());
  VectorTools::get_position_vector(dof_handler, location_vector, mask);
  MappingFEField<dim, spacedim> mapping(dof_handler, location_vector, mask);
  test_real_to_unit_cell(mapping);
}



int
main()
{
  initlog();

  test<1, 1>();
  test<2, 2>();
  test<3, 3>();

  test<1, 2>();
}

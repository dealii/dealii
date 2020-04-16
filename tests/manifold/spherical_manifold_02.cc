// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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

// Check SphericalManifold for get_intermediate_point and get_tangent_vector
// issues.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_manifold.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <memory>

#include "../tests.h"


struct MappingEnum
{
  enum type
  {
    MappingManifold,
    MappingQ
  };
};

void
test(MappingEnum::type mapping_name, unsigned int refinements = 1)
{
  deallog.depth_console(0);

  const unsigned int degree = 2; // Degree of shape functions

  Triangulation<2, 3> triangulation;

  FE_Q<2, 3>       fe(degree);
  DoFHandler<2, 3> dof_handler(triangulation);
  QGaussLobatto<2> cell_quadrature(degree + 1);



  const double radius = 1.0;
  Point<3>     center(0.0, 0.0, 0.0);
  GridGenerator::hyper_sphere(triangulation, center, radius);

  static const SphericalManifold<2, 3> sphere;
  triangulation.set_manifold(0, sphere);
  // static const RotatedSphericalManifold rotated_sphere;
  // triangulation.set_manifold (1, rotated_sphere);

  for (Triangulation<2, 3>::active_cell_iterator cell =
         triangulation.begin_active();
       cell != triangulation.end();
       ++cell)
    {
      cell->set_all_manifold_ids(0);
      // deallog << "Setting SphericalManifold\n";
    }

  triangulation.refine_global(refinements);
  dof_handler.distribute_dofs(fe);

  if (false) // reenable for visualization
    {
      GridOut       grid_out;
      std::ofstream grid_file("grid.vtk");
      grid_out.write_vtk(triangulation, grid_file);
      // deallog << "Grid has been saved into grid.vtk" << std::endl;
    }

  // deallog << "Surface mesh has " << triangulation.n_active_cells()
  //           << " cells."
  //           << std::endl;
  // deallog << "Surface mesh has " << dof_handler.n_dofs()
  //           << " degrees of freedom."
  //           << std::endl;

  std::shared_ptr<Mapping<2, 3>> mapping;
  switch (mapping_name)
    {
      case MappingEnum::MappingManifold:
        // deallog << " MappingManifold" << std::endl;
        mapping.reset(new MappingManifold<2, 3>());
        break;
      case MappingEnum::MappingQ:
        // deallog << " MappingQ" << std::endl;
        mapping.reset(new MappingQ<2, 3>(fe.degree));
        break;
    }

  FEValues<2, 3> fe_values(*mapping, fe, cell_quadrature, update_JxW_values);
  const unsigned int n_q_points = cell_quadrature.size();

  double surface_area = 0;
  for (DoFHandler<2, 3>::active_cell_iterator cell = dof_handler.begin_active(),
                                              endc = dof_handler.end();
       cell != endc;
       ++cell)
    {
      double patch_surface = 0;
      fe_values.reinit(cell);


      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        {
          patch_surface += fe_values.JxW(q_point);
          // deallog << "--> " << qp[q_point] << std::endl;
        }
      // deallog  << " Patch area       = "
      //            << patch_surface << std::endl;
      surface_area += patch_surface;
    }

  deallog << " Ref      = " << std::setw(5) << refinements;
  // deallog << " Surface area     = "
  //         << surface_area << std::endl;
  deallog << "  RelErr  = "
          << (surface_area - 4 * numbers::PI * radius * radius) /
               (4 * numbers::PI * radius * radius)
          << std::endl;

  return;
}

int
main()
{
  initlog();

  std::string bar(35, '-');

  deallog << bar << std::endl;
  for (unsigned int i = 1; i < 8; ++i)
    test(MappingEnum::MappingManifold, i);
  deallog << bar << std::endl;
  for (unsigned int i = 1; i < 8; ++i)
    test(MappingEnum::MappingQ, i);
  deallog << bar << std::endl;
}

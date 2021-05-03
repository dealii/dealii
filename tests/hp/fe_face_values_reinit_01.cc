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



// Test different versions of hp::FEFaceValues::reinit().


#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

#include "../tests.h"



template <int dim, int spacedim = dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);

  hp::FECollection<dim>      fe_collection(FE_Q<dim>(1));
  hp::QCollection<dim - 1>   quad_collection(QGauss<dim - 1>(2));
  hp::MappingCollection<dim> mapping_collection(MappingQ<dim, spacedim>(1));

  DoFHandler<dim> dof_handler(triangulation);

  dof_handler.distribute_dofs(fe_collection);



  hp::FEFaceValues<dim, spacedim> hp_fe_face_values(mapping_collection,
                                                    fe_collection,
                                                    quad_collection,
                                                    update_quadrature_points);

  for (const auto &cell : triangulation.active_cell_iterators())
    {
      for (const auto face : cell->face_indices())
        {
          hp_fe_face_values.reinit(cell, face);

          for (const auto &p : hp_fe_face_values.get_present_fe_values()
                                 .get_quadrature_points())
            deallog << p << " ";
          deallog << std::endl;
        }
      deallog << std::endl;
    }

  for (const auto &cell : triangulation.active_cell_iterators())
    {
      for (const auto face : cell->face_iterators())
        {
          hp_fe_face_values.reinit(cell, face);

          for (const auto &p : hp_fe_face_values.get_present_fe_values()
                                 .get_quadrature_points())
            deallog << p << " ";
          deallog << std::endl;
        }
      deallog << std::endl;
    }

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      for (const auto face : cell->face_indices())
        {
          hp_fe_face_values.reinit(cell, face);

          for (const auto &p : hp_fe_face_values.get_present_fe_values()
                                 .get_quadrature_points())
            deallog << p << " ";
          deallog << std::endl;
        }
      deallog << std::endl;
    }

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      for (const auto face : cell->face_iterators())
        {
          hp_fe_face_values.reinit(cell, face);

          for (const auto &p : hp_fe_face_values.get_present_fe_values()
                                 .get_quadrature_points())
            deallog << p << " ";
          deallog << std::endl;
        }
      deallog << std::endl;
    }
}



int
main()
{
  initlog();
  deallog.get_file_stream().precision(2);

  test<2>();
  test<3>();
}

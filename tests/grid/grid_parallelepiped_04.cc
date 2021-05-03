// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
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

// check boundary indicators of colored subdivided_parallelepiped

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <map>

#include "../tests.h"

template <int dim>
void
check_parallelepiped(bool colorize, bool log, const unsigned int (&subd)[dim])
{
  deallog << "* checking dim=" << dim << " subd=";
  for (unsigned int i = 0; i < dim; ++i)
    deallog << subd[i] << " ";
  deallog << std::endl;

  // Data structure defining dim coordinates that make up a
  // parallelepiped.
  Point<dim>(corners)[dim];

  // build corners for this particular dim:
  switch (dim)
    {
      case 1:
        corners[0] = Point<dim>(0.5);
        break;

      case 2:
        corners[0] = Point<dim>(0.0, 0.5);
        corners[1] = Point<dim>(0.5, 0.0);
        break;

      case 3:
        corners[0] = Point<dim>(0.0, 0.3, 0.5);
        corners[1] = Point<dim>(0.4, 0.0, 0.5);
        corners[2] = Point<dim>(0.4, 0.3, 0.0);
        break;

      default:
        Assert(false, ExcInternalError());
    }

  Triangulation<dim> triangulation;

  GridGenerator::subdivided_parallelepiped(triangulation,
                                           subd,
                                           corners,
                                           colorize);

  {
    std::map<unsigned int, unsigned int>              boundary_count;
    typename Triangulation<dim>::active_cell_iterator cell = triangulation
                                                               .begin_active(),
                                                      endc =
                                                        triangulation.end();
    for (; cell != endc; ++cell)
      {
        for (const unsigned int face : GeometryInfo<dim>::face_indices())
          {
            if (cell->face(face)->at_boundary())
              {
                boundary_count[cell->face(face)->boundary_id()]++;
                deallog << " center: " << cell->face(face)->center()
                        << " id: " << (int)cell->face(face)->boundary_id()
                        << std::endl;
              }
          }
      }

    deallog << " boundary indicators: ";
    for (std::map<unsigned int, unsigned int>::iterator it =
           boundary_count.begin();
         it != boundary_count.end();
         ++it)
      {
        deallog << it->first << "(" << it->second << " times) ";
      }
    deallog << std::endl;
  }

  GridOut grid_out;

  if (log)
    {
      FE_Q<dim>       fe(2);
      DoFHandler<dim> dh(triangulation);
      dh.distribute_dofs(fe);
      DataOut<dim> d_o;
      d_o.attach_dof_handler(dh);
      Vector<double>            vec(dh.n_dofs());
      AffineConstraints<double> constraints;
      for (unsigned int c = 0; c < 6; ++c)
        VectorTools::interpolate_boundary_values(
          dh, c, Functions::ConstantFunction<dim>(c), constraints);
      constraints.close();
      constraints.distribute(vec);

      d_o.add_data_vector(vec, "v");
      d_o.build_patches(2);
      char       fname[] = "0.vtk";
      static int counter = 0;
      fname[0] += counter;
      ++counter;

      std::ofstream ss(fname);

      // d_o.write_vtk(ss);
      d_o.write_vtk(deallog.get_file_stream());
    }
}

int
main()
{
  initlog(true);

  // check_parallelepiped<1> (false, true);
  // check_parallelepiped<2> (false, true);
  for (unsigned int subd = 1; subd <= 3; ++subd)
    {
      unsigned int subdivisions[3] = {subd, subd, subd};
      check_parallelepiped<3>(true, subd == 2, subdivisions);
    }
  for (unsigned int subd = 1; subd <= 3; ++subd)
    {
      unsigned int subdivisions[3] = {1, 2, subd};
      check_parallelepiped<3>(true, false, subdivisions);
    }
}

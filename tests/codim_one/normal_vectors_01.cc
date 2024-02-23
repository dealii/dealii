// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


//  Asking for face normal vectors in the codim-1 case led to aborts.


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
test()
{
  Triangulation<dim, dim> volume_mesh;
  GridGenerator::hyper_cube(volume_mesh);

  Triangulation<dim - 1, dim> tria;
  GridGenerator::extract_boundary_mesh(volume_mesh, tria);

  FE_Q<dim - 1, dim>       fe(1);
  DoFHandler<dim - 1, dim> dh(tria);
  MappingQ<dim - 1, dim>   mapping(1);

  dh.distribute_dofs(fe);

  FEFaceValues<dim - 1, dim> fe_face_values(mapping,
                                            fe,
                                            QTrapezoid<dim - 2>(),
                                            update_normal_vectors);

  for (typename DoFHandler<dim - 1, dim>::active_cell_iterator cell =
         dh.begin_active();
       cell != dh.end();
       ++cell)
    {
      deallog << "Face centered at " << cell->center() << std::endl;

      for (unsigned int f = 0; f < GeometryInfo<dim - 1>::faces_per_cell; ++f)
        {
          deallog << "  Edge centered at " << cell->face(f)->center()
                  << std::endl;

          fe_face_values.reinit(cell, f);
          for (unsigned int q = 0; q < fe_face_values.n_quadrature_points; ++q)
            deallog << "    normal_vector=" << fe_face_values.normal_vector(q)
                    << std::endl;
        }
    }
}



int
main()
{
  initlog();

  //  test<2> ();
  test<3>();

  return 0;
}

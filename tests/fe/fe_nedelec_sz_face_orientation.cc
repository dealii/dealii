// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test FE_NedelecSZ<3> for meshes with faces with non-standard orientation.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_nedelec_sz.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

void
create_reference_triangulation(Triangulation<3> &tria)
{
  std::vector<unsigned int> repetitions(3, 1);

  repetitions[0] = 2;
  GridGenerator::subdivided_hyper_rectangle(tria,
                                            repetitions,
                                            Point<3>(-1.0, 0.0, 0.0),
                                            Point<3>(1.0, 1.0, 1.0));
}

void
create_triangulation(Triangulation<3> &tria,
                     const bool        face_orientation,
                     const bool        face_flip,
                     const bool        face_rotation)
{
  std::vector<CellData<3>> cells(2);

  cells[0].vertices[0] = 0;
  cells[0].vertices[1] = 1;
  cells[0].vertices[2] = 6;
  cells[0].vertices[3] = 7;
  cells[0].vertices[4] = 3;
  cells[0].vertices[5] = 4;
  cells[0].vertices[6] = 9;
  cells[0].vertices[7] = 10;
  cells[0].material_id = 0;

  if (face_orientation)
    {
      if (face_flip)
        {
          if (face_rotation)
            {
              cells[1].vertices[0] = 7;
              cells[1].vertices[1] = 8;
              cells[1].vertices[2] = 10;
              cells[1].vertices[3] = 11;
              cells[1].vertices[4] = 1;
              cells[1].vertices[5] = 2;
              cells[1].vertices[6] = 4;
              cells[1].vertices[7] = 5;
            }

          else
            {
              cells[1].vertices[0] = 10;
              cells[1].vertices[1] = 11;
              cells[1].vertices[2] = 4;
              cells[1].vertices[3] = 5;
              cells[1].vertices[4] = 7;
              cells[1].vertices[5] = 8;
              cells[1].vertices[6] = 1;
              cells[1].vertices[7] = 2;
            }
        }

      else if (face_rotation)
        {
          cells[1].vertices[0] = 4;
          cells[1].vertices[1] = 5;
          cells[1].vertices[2] = 1;
          cells[1].vertices[3] = 2;
          cells[1].vertices[4] = 10;
          cells[1].vertices[5] = 11;
          cells[1].vertices[6] = 7;
          cells[1].vertices[7] = 8;
        }

      else
        {
          cells[1].vertices[0] = 2;
          cells[1].vertices[1] = 1;
          cells[1].vertices[2] = 5;
          cells[1].vertices[3] = 4;
          cells[1].vertices[4] = 8;
          cells[1].vertices[5] = 7;
          cells[1].vertices[6] = 11;
          cells[1].vertices[7] = 10;
        }
    }

  else if (face_flip)
    {
      if (face_rotation)
        {
          cells[1].vertices[0] = 8;
          cells[1].vertices[1] = 7;
          cells[1].vertices[2] = 2;
          cells[1].vertices[3] = 1;
          cells[1].vertices[4] = 11;
          cells[1].vertices[5] = 10;
          cells[1].vertices[6] = 5;
          cells[1].vertices[7] = 4;
        }

      else
        {
          cells[1].vertices[0] = 11;
          cells[1].vertices[1] = 10;
          cells[1].vertices[2] = 8;
          cells[1].vertices[3] = 7;
          cells[1].vertices[4] = 5;
          cells[1].vertices[5] = 4;
          cells[1].vertices[6] = 2;
          cells[1].vertices[7] = 1;
        }
    }

  else if (face_rotation)
    {
      cells[1].vertices[0] = 5;
      cells[1].vertices[1] = 4;
      cells[1].vertices[2] = 11;
      cells[1].vertices[3] = 10;
      cells[1].vertices[4] = 2;
      cells[1].vertices[5] = 1;
      cells[1].vertices[6] = 8;
      cells[1].vertices[7] = 7;
    }

  else
    {
      cells[1].vertices[0] = 1;
      cells[1].vertices[1] = 2;
      cells[1].vertices[2] = 7;
      cells[1].vertices[3] = 8;
      cells[1].vertices[4] = 4;
      cells[1].vertices[5] = 5;
      cells[1].vertices[6] = 10;
      cells[1].vertices[7] = 11;
    }

  cells[1].material_id = 0;

  std::vector<Point<3>> vertices(12);

  vertices[0]  = Point<3>(-1.0, 0.0, 0.0);
  vertices[1]  = Point<3>();
  vertices[2]  = Point<3>(1.0, 0.0, 0.0);
  vertices[3]  = Point<3>(-1.0, 0.0, 1.0);
  vertices[4]  = Point<3>(0.0, 0.0, 1.0);
  vertices[5]  = Point<3>(1.0, 0.0, 1.0);
  vertices[6]  = Point<3>(-1.0, 1.0, 0.0);
  vertices[7]  = Point<3>(0.0, 1.0, 0.0);
  vertices[8]  = Point<3>(1.0, 1.0, 0.0);
  vertices[9]  = Point<3>(-1.0, 1.0, 1.0);
  vertices[10] = Point<3>(0.0, 1.0, 1.0);
  vertices[11] = Point<3>(1.0, 1.0, 1.0);
  tria.create_triangulation(vertices, cells, SubCellData());
}

void
evaluate(const FiniteElement<3> &fe, const DoFHandler<3> &dof_handler)
{
  const FEValuesExtractors::Vector component(0);
  const Quadrature<3>              quadrature(Point<3>(0.5, 0.5, 0.5));
  FEValues<3>                      fe_values(fe,
                        quadrature,
                        update_quadrature_points | update_values |
                          update_gradients);

  for (DoFHandler<3>::active_cell_iterator cell = dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    {
      fe_values.reinit(cell);

      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
        {
          deallog << "DoF#" << i << ", value=["
                  << fe_values[component].value(i, 0) << "], curl=["
                  << fe_values[component].curl(i, 0) << ']' << std::endl;
        }
    }
}

void
run(const bool face_orientation, const bool face_flip, const bool face_rotation)
{
  //  Triangulation<3> tria_ref;
  //  create_reference_triangulation (tria_ref);

  FE_NedelecSZ<3> fe(1);
  //  DoFHandler<3> dof_handler_ref (tria_ref);
  //  dof_handler_ref.distribute_dofs (fe);

  Triangulation<3> tria;
  create_triangulation(tria, face_orientation, face_flip, face_rotation);

  DoFHandler<3> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  evaluate(fe, dof_handler);
}

int
main()
{
  initlog();

  for (int face_orientation = 0; face_orientation <= 1; ++face_orientation)
    for (int face_flip = 0; face_flip <= 1; ++face_flip)
      for (int face_rotation = 0; face_rotation <= 1; ++face_rotation)
        {
          deallog << face_orientation << face_flip << face_rotation
                  << std::endl;
          run(face_orientation, face_flip, face_rotation);
        }
}

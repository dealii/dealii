// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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

// Test FE_Nedelec<3> for meshes with faces with non-standard orientation.

#include "../tests.h"
#include<deal.II/base/quadrature_lib.h>
#include<deal.II/dofs/dof_handler.h>
#include<deal.II/fe/fe_nedelec.h>
#include<deal.II/fe/fe_values.h>
#include<deal.II/grid/grid_generator.h>
#include<deal.II/grid/tria.h>
#include<deal.II/lac/constraint_matrix.h>
#include<deal.II/lac/vector.h>
#include<deal.II/numerics/fe_field_function.h>
#include<deal.II/numerics/vector_tools.h>

void create_reference_triangulation (Triangulation<3> &tria)
{
  std::vector<unsigned int> repetitions (3, 1);

  repetitions[0] = 2;
  GridGenerator::subdivided_hyper_rectangle (tria, repetitions, Point<3> (-1.0, 0.0, 0.0), Point<3> (1.0, 1.0, 1.0));
}

void create_triangulation (Triangulation<3> &tria, const bool face_orientation, const bool face_flip, const bool face_rotation)
{
  std::vector<CellData<3> > cells (2);

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

  std::vector<Point<3> > vertices (12);

  vertices[0] = Point<3> (-1.0, 0.0, 0.0);
  vertices[1] = Point<3> ();
  vertices[2] = Point<3> (1.0, 0.0, 0.0);
  vertices[3] = Point<3> (-1.0, 0.0, 1.0);
  vertices[4] = Point<3> (0.0, 0.0, 1.0);
  vertices[5] = Point<3> (1.0, 0.0, 1.0);
  vertices[6] = Point<3> (-1.0, 1.0, 0.0);
  vertices[7] = Point<3> (0.0, 1.0, 0.0);
  vertices[8] = Point<3> (1.0, 1.0, 0.0);
  vertices[9] = Point<3> (-1.0, 1.0, 1.0);
  vertices[10] = Point<3> (0.0, 1.0, 1.0);
  vertices[11] = Point<3> (1.0, 1.0, 1.0);
  tria.create_triangulation (vertices, cells, SubCellData ());
}

void evaluate (const FE_Nedelec<3> &fe, const DoFHandler<3> &dof_handler_ref, const Vector<double> &u_ref, const DoFHandler<3> &dof_handler, const Vector<double> &u)
{
  const FEValuesExtractors::Vector component (0);
  const QGauss<3> quadrature (2);
  const unsigned int n_q_points = quadrature.size ();
  Functions::FEFieldFunction<3> fe_field_function (dof_handler, u);
  FEValues<3> fe_values (fe, quadrature, update_quadrature_points | update_values);
  std::vector<Vector<double> > values (n_q_points, Vector<double> (3));
  std::vector<Tensor<1, 3> > values_ref (n_q_points);

  for (DoFHandler<3>::active_cell_iterator cell = dof_handler_ref.begin_active (); cell != dof_handler_ref.end (); ++cell)
    {
      fe_values.reinit (cell);
      fe_values[component].get_function_values (u_ref, values_ref);
      fe_field_function.vector_value_list (fe_values.get_quadrature_points (), values);

      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        {
          for (unsigned int d = 0; d < 3; ++d)
            deallog << values_ref[q_point][d] - values[q_point] (d) << "  ";

          deallog << std::endl;
        }
    }
}

void set_reference_solution (Vector<double> &vector)
{
  for (unsigned int i = 0; i < vector.size (); ++i)
    vector (i) = 1.0;
}

void set_solution (Vector<double> &vector, const DoFHandler<3> &dof_handler, const DoFHandler<3> &dof_handler_ref, const Vector<double> &u_ref)
{
  ConstraintMatrix constraints;

  constraints.close ();

  Functions::FEFieldFunction<3> fe_field_function (dof_handler_ref, u_ref);

  VectorTools::project (dof_handler, constraints, QGauss<3> (2), fe_field_function, vector);
}

void run (const bool face_orientation, const bool face_flip, const bool face_rotation)
{
  Triangulation<3> tria_ref;

  create_reference_triangulation (tria_ref);

  FE_Nedelec<3> fe (0);
  DoFHandler<3> dof_handler_ref (tria_ref);

  dof_handler_ref.distribute_dofs (fe);

  Vector<double> u_ref (dof_handler_ref.n_dofs ());

  set_reference_solution (u_ref);

  Triangulation<3> tria;

  create_triangulation (tria, face_orientation, face_flip, face_rotation);

  DoFHandler<3> dof_handler (tria);

  dof_handler.distribute_dofs (fe);

  Vector<double> u (dof_handler.n_dofs ());

  set_solution (u, dof_handler, dof_handler_ref, u_ref);
  evaluate (fe, dof_handler_ref, u_ref, dof_handler, u);
}

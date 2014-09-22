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

void create_reference_triangulation (Triangulation<3> &triangulation)
{
  GridGenerator::hyper_cube (triangulation, -1.0, 1.0);
  triangulation.refine_global(1);
}

void create_triangulation (Triangulation<3> &triangulation)
{
  static const Point<3> vertices_parallelograms[] =
  {
    Point<3> (-1., -1., -1.),  //0
    Point<3> (0., -1., -1.),
    Point<3> (1., -1., -1.),

    Point<3> (-1., -1., 0.),   //3
    Point<3> (0., -1., 0.),
    Point<3> (1., -1., 0.),

    Point<3> (-1., -1., 1.),   //6
    Point<3> (0., -1., 1.),
    Point<3> (1., -1., 1.),

    Point<3> (-1., 0., -1.),   //9
    Point<3> (0., 0., -1.),
    Point<3> (1., 0., -1.),

    Point<3> (-1., 0., 0.),    //12
    Point<3> (0., 0., 0.),
    Point<3> (1., 0., 0.),

    Point<3> (-1., 0., 1.),    //15
    Point<3> (0., 0., 1.),
    Point<3> (1., 0., 1.),

    Point<3> (-1., 1., -1.),   //18
    Point<3> (0., 1., -1.),
    Point<3> (1., 1., -1.),

    Point<3> (-1., 1., 0.),    //21
    Point<3> (0., 1., 0.),
    Point<3> (1., 1., 0.),

    Point<3> (-1., 1., 1.),    //24
    Point<3> (0., 1., 1.),
    Point<3> (1., 1., 1.)
  };
  const unsigned n_vertices = sizeof(vertices_parallelograms) / sizeof(vertices_parallelograms[0]);

  const unsigned n_cells = 8;

  const std::vector<Point<3> > vertices (&vertices_parallelograms[0],
                                         &vertices_parallelograms[n_vertices]);

  // create grid with all possible combintations of face_flip, face_orientation and face_rotation flags
  static const int cell_vertices[][GeometryInfo<3>::vertices_per_cell] =
  {
    {0, 1, 9, 10, 3, 4, 12, 13},        // cell 1 standard
    {10, 11, 13, 14, 1, 2, 4, 5},       // cell 2 rotated by 270 deg
    {9, 10, 18, 19, 12, 13, 21, 22},    // cell 3 standard
    {13, 14, 10, 11, 22, 23, 19, 20}, // cell 4 rotated by 90 deg
    {3, 4, 12, 13, 6, 7, 15, 16},       // cell 5 standard
    {4, 5, 13, 14, 7, 8, 16, 17},       // cell 6 standard
    {24, 25, 15, 16, 21, 22, 12, 13}, // cell 7 rotated by 180 deg
    {13, 14, 22, 23, 16, 17, 25, 26}    // cell 8 standard
  };

  std::vector<CellData<3> > cells (n_cells, CellData<3>());
  for (unsigned i = 0; i<n_cells; ++i)
    {
      for (unsigned int j=0; j<GeometryInfo<3>::vertices_per_cell; ++j)
        cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    }

  triangulation.create_triangulation (vertices, cells, SubCellData());
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
      std::vector<types::global_dof_index> dof_indices (fe.dofs_per_cell);
      cell->get_dof_indices (dof_indices);
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

int main ()
{
  initlog();
  deallog.threshold_double (1.e-10);

  Triangulation<3> tria_ref;

  create_reference_triangulation (tria_ref);

  FE_Nedelec<3> fe (0);
  DoFHandler<3> dof_handler_ref (tria_ref);

  dof_handler_ref.distribute_dofs (fe);

  Vector<double> u_ref (dof_handler_ref.n_dofs ());

  set_reference_solution (u_ref);

  Triangulation<3> tria;

  create_triangulation (tria);

  DoFHandler<3> dof_handler (tria);

  dof_handler.distribute_dofs (fe);

  Vector<double> u (dof_handler.n_dofs ());

  set_solution (u, dof_handler, dof_handler_ref, u_ref);
  evaluate (fe, dof_handler_ref, u_ref, dof_handler, u);
}

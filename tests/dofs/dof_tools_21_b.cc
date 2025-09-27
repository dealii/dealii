// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>

#include <bitset>
#include <iostream>
#include <utility>

#include "../tests.h"

//
// Test
//   DoFTools::
//   make_periodicity_constraints(const FaceIterator        &,
//                                const FaceIterator        &,
//                                AffineConstraints<double> &,
//                                const ComponentMask       &,
//                                types::geometric_orientation)
// for correct behavior on non standard oriented meshes.
//


/*
 * Generate a grid consisting of two disjoint cells, colorize the two
 * outermost faces. They will be matched via collect_periodic_faces
 *
 * The integer orientation determines the orientation of the second cell
 * to get something else than the boring default orientation.
 */

/* The 2D case */
void
generate_grid(Triangulation<2> &triangulation, int orientation)
{
  Point<2> vertices_1[] = {
    Point<2>(-1., -3.),
    Point<2>(+1., -3.),
    Point<2>(-1., -1.),
    Point<2>(+1., -1.),
    Point<2>(-1., +1.),
    Point<2>(+1., +1.),
    Point<2>(-1., +3.),
    Point<2>(+1., +3.),
  };
  std::vector<Point<2>> vertices(&vertices_1[0], &vertices_1[8]);

  std::vector<CellData<2>> cells(2, CellData<2>());

  /* cell 0 */
  int cell_vertices_0[GeometryInfo<2>::vertices_per_cell] = {0, 1, 2, 3};

  /* cell 1 */
  int cell_vertices_1[2][GeometryInfo<2>::vertices_per_cell] = {
    {4, 5, 6, 7},
    {7, 6, 5, 4},
  };

  for (const unsigned int j : GeometryInfo<2>::vertex_indices())
    {
      cells[0].vertices[j] = cell_vertices_0[j];
      cells[1].vertices[j] = cell_vertices_1[orientation][j];
    }
  cells[0].material_id = 0;
  cells[1].material_id = 0;

  triangulation.create_triangulation(vertices, cells, SubCellData());

  Triangulation<2>::cell_iterator cell_1 = triangulation.begin();
  Triangulation<2>::cell_iterator cell_2 = cell_1++;
  Triangulation<2>::face_iterator face_1;
  Triangulation<2>::face_iterator face_2;

  // Look for the two outermost faces:
  for (const unsigned int j : GeometryInfo<2>::face_indices())
    {
      if (cell_1->face(j)->center()[1] > 2.9)
        face_1 = cell_1->face(j);
      if (cell_2->face(j)->center()[1] < -2.9)
        face_2 = cell_2->face(j);
    }
  face_1->set_boundary_id(42);
  face_2->set_boundary_id(43);

  triangulation.refine_global(0);
}

/* The 2D in 3D case */
void
generate_grid(Triangulation<2, 3> &triangulation, int orientation)
{
  Point<3> vertices_1[] = {
    Point<3>(-1., -3., 0.),
    Point<3>(+1., -3., 0.),
    Point<3>(-1., -1., 0.),
    Point<3>(+1., -1., 0.),
    Point<3>(-1., +1., 0.),
    Point<3>(+1., +1., 0.),
    Point<3>(-1., +3., 0.),
    Point<3>(+1., +3., 0.),
  };
  std::vector<Point<3>> vertices(&vertices_1[0], &vertices_1[8]);

  std::vector<CellData<2>> cells(2, CellData<2>());

  /* cell 0 */
  int cell_vertices_0[GeometryInfo<2>::vertices_per_cell] = {0, 1, 2, 3};

  /* cell 1 */
  int cell_vertices_1[2][GeometryInfo<2>::vertices_per_cell] = {
    {4, 5, 6, 7},
    {7, 6, 5, 4},
  };

  for (const unsigned int j : GeometryInfo<2>::vertex_indices())
    {
      cells[0].vertices[j] = cell_vertices_0[j];
      cells[1].vertices[j] = cell_vertices_1[orientation][j];
    }
  cells[0].material_id = 0;
  cells[1].material_id = 0;

  triangulation.create_triangulation(vertices, cells, SubCellData());

  Triangulation<2, 3>::cell_iterator cell_1 = triangulation.begin();
  Triangulation<2, 3>::cell_iterator cell_2 = cell_1++;
  Triangulation<2, 3>::face_iterator face_1;
  Triangulation<2, 3>::face_iterator face_2;

  // Look for the two outermost faces:
  for (const unsigned int j : GeometryInfo<2>::face_indices())
    {
      if (cell_1->face(j)->center()[1] > 2.9)
        face_1 = cell_1->face(j);
      if (cell_2->face(j)->center()[1] < -2.9)
        face_2 = cell_2->face(j);
    }
  face_1->set_boundary_id(42);
  face_2->set_boundary_id(43);

  triangulation.refine_global(0);
}

/* The 3D case */
void
generate_grid(Triangulation<3> &triangulation, int orientation)
{
  Point<3>              vertices_1[] = {Point<3>(-1., -1., -3.),
                                        Point<3>(+1., -1., -3.),
                                        Point<3>(-1., +1., -3.),
                                        Point<3>(+1., +1., -3.),
                                        Point<3>(-1., -1., -1.),
                                        Point<3>(+1., -1., -1.),
                                        Point<3>(-1., +1., -1.),
                                        Point<3>(+1., +1., -1.),
                                        Point<3>(-1., -1., +1.),
                                        Point<3>(+1., -1., +1.),
                                        Point<3>(-1., +1., +1.),
                                        Point<3>(+1., +1., +1.),
                                        Point<3>(-1., -1., +3.),
                                        Point<3>(+1., -1., +3.),
                                        Point<3>(-1., +1., +3.),
                                        Point<3>(+1., +1., +3.)};
  std::vector<Point<3>> vertices(&vertices_1[0], &vertices_1[16]);

  std::vector<CellData<3>> cells(2, CellData<3>());

  /* cell 0 */
  int cell_vertices_0[GeometryInfo<3>::vertices_per_cell] = {
    0, 1, 2, 3, 4, 5, 6, 7};

  /* cell 1 */
  int cell_vertices_1[8][GeometryInfo<3>::vertices_per_cell] = {
    {8, 9, 10, 11, 12, 13, 14, 15},
    {9, 11, 8, 10, 13, 15, 12, 14},
    {11, 10, 9, 8, 15, 14, 13, 12},
    {10, 8, 11, 9, 14, 12, 15, 13},
    {13, 12, 15, 14, 9, 8, 11, 10},
    {12, 14, 13, 15, 8, 10, 9, 11},
    {14, 15, 12, 13, 10, 11, 8, 9},
    {15, 13, 14, 12, 11, 9, 10, 8},
  };

  for (const unsigned int j : GeometryInfo<3>::vertex_indices())
    {
      cells[0].vertices[j] = cell_vertices_0[j];
      cells[1].vertices[j] = cell_vertices_1[orientation][j];
    }
  cells[0].material_id = 0;
  cells[1].material_id = 0;


  triangulation.create_triangulation(vertices, cells, SubCellData());

  Triangulation<3>::cell_iterator cell_1 = triangulation.begin();
  Triangulation<3>::cell_iterator cell_2 = cell_1++;
  Triangulation<3>::face_iterator face_1;
  Triangulation<3>::face_iterator face_2;

  // Look for the two outermost faces:
  for (const unsigned int j : GeometryInfo<3>::face_indices())
    {
      if (cell_1->face(j)->center()[2] > 2.9)
        face_1 = cell_1->face(j);
      if (cell_2->face(j)->center()[2] < -2.9)
        face_2 = cell_2->face(j);
    }
  face_1->set_boundary_id(42);
  face_2->set_boundary_id(43);

  triangulation.refine_global(0);
}


/*
 * Print out all face DoFs and support points as well as the actual
 * matching via make_periodicity_constraints()
 */
template <int dim, int spacedim>
void
print_matching(DoFHandler<dim, spacedim> &dof_handler,
               bool                       constrain_only_velocity = false)
{
  const FiniteElement<dim, spacedim> &fe = dof_handler.get_fe();
  MappingQ<dim, spacedim>             mapping(1);

  AffineConstraints<double>    constraint_matrix;
  AffineConstraints<double>    constraint_matrix_reverse;
  std::vector<Point<spacedim>> support_points(dof_handler.n_dofs());
  DoFTools::map_dofs_to_support_points<dim, spacedim>(mapping,
                                                      dof_handler,
                                                      support_points);

  const FEValuesExtractors::Vector v(0);
  const FEValuesExtractors::Scalar v_1(0);
  const FEValuesExtractors::Scalar v_2(1);
  const FEValuesExtractors::Scalar v_3(2);
  const FEValuesExtractors::Scalar pressure(3);

  ComponentMask velocity_mask;
  if (constrain_only_velocity)
    velocity_mask = fe.component_mask(v);

  // Look for the two outermost faces:
  typename DoFHandler<dim, spacedim>::face_iterator face_1;
  typename DoFHandler<dim, spacedim>::face_iterator face_2;
  for (typename DoFHandler<dim, spacedim>::cell_iterator cell =
         dof_handler.begin(0);
       cell != dof_handler.end(0);
       ++cell)
    {
      for (const unsigned int j : GeometryInfo<dim>::face_indices())
        {
          if (cell->face(j)->center()[dim == 2 ? 1 : 2] > 2.9)
            face_1 = cell->face(j);
          if (cell->face(j)->center()[dim == 2 ? 1 : 2] < -2.9)
            face_2 = cell->face(j);
        }
    }

  // Determine the orientation of the two faces:

  std::vector<types::global_dof_index> dofs_1(fe.dofs_per_face);
  std::vector<types::global_dof_index> dofs_2(fe.dofs_per_face);
  const unsigned int face_1_index = face_1->nth_active_fe_index(0);
  const unsigned int face_2_index = face_2->nth_active_fe_index(0);
  face_1->get_dof_indices(dofs_1, face_1_index);
  face_2->get_dof_indices(dofs_2, face_2_index);

  // Print out all DoF support points on the two faces:
  deallog << "DoFs of face_1:";
  for (unsigned int c = 0; c < fe.n_components(); ++c)
    {
      deallog << std::endl << " component " << c << ':';
      for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
        {
          if (fe.face_system_to_component_index(i).first == c)
            deallog << " (" << dofs_1[i] << " - " << support_points[dofs_1[i]]
                    << ')';
        }
    }
  deallog << std::endl;
  deallog << "DoFs of face_2:";
  for (unsigned int c = 0; c < fe.n_components(); ++c)
    {
      deallog << std::endl << " component " << c << ':';
      for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
        {
          if (fe.face_system_to_component_index(i).first == c)
            deallog << " (" << dofs_2[i] << " - " << support_points[dofs_2[i]]
                    << ')';
        }
    }
  deallog << std::endl;


  const auto orientation = GridTools::orthogonal_equality(
    face_1, face_2, dim == 2 ? 1 : 2, Tensor<1, spacedim>());
  AssertThrow(orientation, ExcInternalError());
  const auto o = *orientation;
  // Preserve old output by printing a bitset and also using the old
  // (orientation, flip, rotation) format
  std::bitset<3> bo(o);
  deallog << "Orientation: " << bo[0] << bo[2] << bo[1] << std::endl;

  DoFTools::make_periodicity_constraints(
    face_1, face_2, constraint_matrix, velocity_mask, o);
  deallog << "Matching:" << std::endl;
  constraint_matrix.print(deallog.get_file_stream());
  constraint_matrix.close();

  const auto reverse_orientation = GridTools::orthogonal_equality(
    face_2, face_1, dim == 2 ? 1 : 2, Tensor<1, spacedim>());
  AssertThrow(reverse_orientation, ExcInternalError());
  const auto ro = *reverse_orientation;
  DoFTools::make_periodicity_constraints(
    face_2, face_1, constraint_matrix_reverse, velocity_mask, ro);
  deallog << "Reverse Matching:" << std::endl;
  constraint_matrix_reverse.print(deallog.get_file_stream());
  constraint_matrix_reverse.close();
}



int
main()
{
  initlog();
  deallog << std::setprecision(4);
  deallog.get_file_stream() << std::setprecision(4);

  deallog << "Test for 2D, Q1:" << std::endl << std::endl;

  for (int i = 0; i < 2; ++i)
    {
      // Generate a triangulation and match:
      Triangulation<2> triangulation;
      FE_Q<2>          fe(1);
      DoFHandler<2>    dof_handler(triangulation);

      deallog << "Triangulation:" << i << std::endl;

      generate_grid(triangulation, i);
      dof_handler.distribute_dofs(fe);
      print_matching(dof_handler);
    }

  deallog << "Test for 2D in 3D, Q1:" << std::endl << std::endl;

  for (int i = 0; i < 2; ++i)
    {
      // Generate a triangulation and match:
      Triangulation<2, 3> triangulation;
      FE_Q<2, 3>          fe(1);
      DoFHandler<2, 3>    dof_handler(triangulation);

      deallog << "Triangulation:" << i << std::endl;

      generate_grid(triangulation, i);
      dof_handler.distribute_dofs(fe);
      print_matching(dof_handler);
    }

  deallog << "Test for 3D, Q1:" << std::endl << std::endl;

  for (int i = 0; i < 8; ++i)
    {
      // Generate a triangulation and match:
      Triangulation<3> triangulation;
      FE_Q<3>          fe(1);
      DoFHandler<3>    dof_handler(triangulation);

      deallog << "Triangulation:" << i << std::endl;

      generate_grid(triangulation, i);
      dof_handler.distribute_dofs(fe);
      print_matching(dof_handler);
    }


  deallog << "Test for 3D, Q1, correct subface iteration:" << std::endl
          << std::endl;

  for (int i = 0; i < 8; ++i)
    {
      // Generate a triangulation and match:
      Triangulation<3> triangulation;
      FE_Q<3>          fe(1);
      DoFHandler<3>    dof_handler(triangulation);

      deallog << "Triangulation:" << i << std::endl;

      generate_grid(triangulation, i);
      triangulation.refine_global(1);
      dof_handler.distribute_dofs(fe);
      print_matching(dof_handler);
    }


  deallog << "Test for 3D, Taylor-Hood with Component-Mask on v:" << std::endl
          << std::endl;

  for (int i = 0; i < 8; ++i)
    {
      // Generate a triangulation and match:
      Triangulation<3> triangulation;
      FE_Q<3>          u(2);
      FE_Q<3>          p(1);
      FESystem<3>      taylor_hood(u, 3, p, 1);

      DoFHandler<3> dof_handler(triangulation);

      deallog << "Triangulation:" << i << std::endl;

      generate_grid(triangulation, i);
      dof_handler.distribute_dofs(taylor_hood);
      print_matching(dof_handler, true);
    }

  return 0;
}

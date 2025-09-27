// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
// like _21_b_x but use a Q3 element. For this, the face_flip has a
// tricky component in that we need not only flip the two vertics of
// the line, but also revert the order of the dofs located on the line
// itself
//
// at the time of writing this test, this does not work correctly. the
// mesh looks like this:
//
//   17-25-24-16
//   |         |
//   22 .   . 20
//   |         |
//   23 .   . 21
//   |         |
//   19-27-26-18
//
//
//   2--10-11--3
//   |         |
//   5  .  .   7
//   |         |
//   4  .  .   6
//   |         |
//   0---8--9--1
//
// we try to match line 0->1 against line 16->17 (not the line
// orientation) against each other, with face_flip=true. this should
// produce constraints
//    16->1
//    17->0
//    24->9
//    25->8
// but at the time of writing the test, it gives
//    24->8
//    25->9
// instead


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>

#include <iostream>
#include <utility>

#include "../tests.h"

/*
 * Generate a grid consisting of two disjoint cells, colorize the two
 * outermost faces. They will be matched via collect_periodic_faces
 *
 * The integer orientation determines the orientation of the second cell
 * to get something else than the boring default orientation.
 */

/* The 2D case */
void
generate_grid(Triangulation<2> &triangulation)
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
  int cell_vertices_1[GeometryInfo<2>::vertices_per_cell] = {7, 6, 5, 4};

  for (const unsigned int j : GeometryInfo<2>::vertex_indices())
    {
      cells[0].vertices[j] = cell_vertices_0[j];
      cells[1].vertices[j] = cell_vertices_1[j];
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
}

/* The 2D in 3D case */
void
generate_grid(Triangulation<2, 3> &triangulation)
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
  int cell_vertices_1[GeometryInfo<2>::vertices_per_cell] = {7, 6, 5, 4};

  for (const unsigned int j : GeometryInfo<2>::vertex_indices())
    {
      cells[0].vertices[j] = cell_vertices_0[j];
      cells[1].vertices[j] = cell_vertices_1[j];
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
}


/*
 * Print out all face DoFs and support points as well as the actual
 * matching via make_periodicity_constraints()
 */
template <int dim, int spacedim>
void
print_matching(DoFHandler<dim, spacedim> &dof_handler)
{
  const FiniteElement<dim, spacedim> &fe = dof_handler.get_fe();
  MappingQ<dim, spacedim>             mapping(1);

  AffineConstraints<double>    constraint_matrix;
  std::vector<Point<spacedim>> support_points(dof_handler.n_dofs());
  DoFTools::map_dofs_to_support_points<dim, spacedim>(mapping,
                                                      dof_handler,
                                                      support_points);

  // Look for the two outermost faces:
  typename DoFHandler<dim, spacedim>::face_iterator face_1 =
    (std::next(dof_handler.begin(0)))->face(2);
  typename DoFHandler<dim, spacedim>::face_iterator face_2 =
    dof_handler.begin(0)->face(2);

  // Determine the orientation of the two faces:

  std::vector<types::global_dof_index> dofs_1(fe.dofs_per_face);
  std::vector<types::global_dof_index> dofs_2(fe.dofs_per_face);
  face_1->get_dof_indices(dofs_1);
  face_2->get_dof_indices(dofs_2);

  // Print out all DoF support points on the two faces:
  deallog << "DoFs of face_1:" << std::endl;
  for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
    deallog << dofs_1[i] << " is located at " << support_points[dofs_1[i]]
            << std::endl;
  deallog << "DoFs of face_2:" << std::endl;
  for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
    deallog << dofs_2[i] << " is located at " << support_points[dofs_2[i]]
            << std::endl;

  DoFTools::make_periodicity_constraints(face_1,
                                         face_2,
                                         constraint_matrix,
                                         ComponentMask(),
                                         numbers::reverse_line_orientation);
  constraint_matrix.print(deallog.get_file_stream());
  constraint_matrix.close();
  deallog << "Matching:" << std::endl;
}



int
main()
{
  initlog();
  deallog << std::setprecision(4);
  deallog.get_file_stream() << std::setprecision(4);

  {
    // Generate a triangulation and match:
    Triangulation<2> triangulation;
    FE_Q<2>          fe(3);
    DoFHandler<2>    dof_handler(triangulation);

    generate_grid(triangulation);
    dof_handler.distribute_dofs(fe);
    print_matching(dof_handler);
  }

  {
    // Generate a triangulation and match:
    Triangulation<2, 3> triangulation;
    FE_Q<2, 3>          fe(3);
    DoFHandler<2, 3>    dof_handler(triangulation);

    generate_grid(triangulation);
    dof_handler.distribute_dofs(fe);
    print_matching(dof_handler);
  }


  return 0;
}

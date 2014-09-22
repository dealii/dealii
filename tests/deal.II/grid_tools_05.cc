// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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



//
// check collect_periodic_faces(b_id) for correct return values
//


#include "../tests.h"
#include <deal.II/base/logstream.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <fstream>
std::ofstream logfile("output");

using namespace dealii;


/*
 * Generate a grid consisting of two disjoint cells, colorize the two
 * outermost faces. They will be matched via collect_periodic_faces.
 * Only default orientation is available for this function.
 *
 */

/* The 2D case */
void generate_grid(Triangulation<2> &triangulation)
{
  Point<2> vertices_1[]
  =
  {
    Point<2> (-1.,-3.),
    Point<2> (+1.,-3.),
    Point<2> (-1.,-1.),
    Point<2> (+1.,-1.),
    Point<2> (-1.,+1.),
    Point<2> (+1.,+1.),
    Point<2> (-1.,+3.),
    Point<2> (+1.,+3.),
  };
  std::vector<Point<2> > vertices (&vertices_1[0], &vertices_1[8]);

  std::vector<CellData<2> > cells (2, CellData<2>());

  /* cell 0 */
  int cell_vertices_0[GeometryInfo<2>::vertices_per_cell] = {0, 1, 2, 3};

  /* cell 1 */
  int cell_vertices_1[GeometryInfo<2>::vertices_per_cell] = {4, 5, 6, 7};

  for (unsigned int j=0; j<GeometryInfo<2>::vertices_per_cell; ++j)
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
  for (unsigned int j=0; j<GeometryInfo<2>::faces_per_cell; ++j)
    {
      if (cell_1->face(j)->center()(1) > 2.9)
        face_1 = cell_1->face(j);
      if (cell_2->face(j)->center()(1) < -2.9)
        face_2 = cell_2->face(j);
    }
  face_1->set_boundary_indicator(42);
  face_2->set_boundary_indicator(42);

  triangulation.refine_global(1);
}


/* The 3D case */
void generate_grid(Triangulation<3> &triangulation)
{
  Point<3> vertices_1[]
  =
  {
    Point<3> (-1.,-1.,-3.),
    Point<3> (+1.,-1.,-3.),
    Point<3> (-1.,+1.,-3.),
    Point<3> (+1.,+1.,-3.),
    Point<3> (-1.,-1.,-1.),
    Point<3> (+1.,-1.,-1.),
    Point<3> (-1.,+1.,-1.),
    Point<3> (+1.,+1.,-1.),
    Point<3> (-1.,-1.,+1.),
    Point<3> (+1.,-1.,+1.),
    Point<3> (-1.,+1.,+1.),
    Point<3> (+1.,+1.,+1.),
    Point<3> (-1.,-1.,+3.),
    Point<3> (+1.,-1.,+3.),
    Point<3> (-1.,+1.,+3.),
    Point<3> (+1.,+1.,+3.)
  };
  std::vector<Point<3> > vertices (&vertices_1[0], &vertices_1[16]);

  std::vector<CellData<3> > cells (2, CellData<3>());

  /* cell 0 */
  int cell_vertices_0[GeometryInfo<3>::vertices_per_cell] = {0, 1,  2,  3,  4,  5,  6,  7};

  /* cell 1 */
  int cell_vertices_1[GeometryInfo<3>::vertices_per_cell] = {8, 9, 10, 11, 12, 13, 14, 15};

  for (unsigned int j=0; j<GeometryInfo<3>::vertices_per_cell; ++j)
    {
      cells[0].vertices[j] = cell_vertices_0[j];
      cells[1].vertices[j] = cell_vertices_1[j];
    }
  cells[0].material_id = 0;
  cells[1].material_id = 0;


  triangulation.create_triangulation(vertices, cells, SubCellData());

  Triangulation<3>::cell_iterator cell_1 = triangulation.begin();
  Triangulation<3>::cell_iterator cell_2 = cell_1++;
  Triangulation<3>::face_iterator face_1;
  Triangulation<3>::face_iterator face_2;

  // Look for the two outermost faces:
  for (unsigned int j=0; j<GeometryInfo<3>::faces_per_cell; ++j)
    {
      if (cell_1->face(j)->center()(2) > 2.9)
        face_1 = cell_1->face(j);
      if (cell_2->face(j)->center()(2) < -2.9)
        face_2 = cell_2->face(j);
    }
  face_1->set_boundary_indicator(42);
  face_2->set_boundary_indicator(42);

  triangulation.refine_global(1);
}



/*
 * Print out the face vertices as well as the orientation of a match:
 */
template<typename FaceIterator>
void print_match(const FaceIterator &face_1,
                 const FaceIterator &face_2,
                 const std::bitset<3> &orientation)
{
  static const int dim = FaceIterator::AccessorType::dimension;

  deallog << "face 1";
  for (unsigned int j=0; j<GeometryInfo<dim>::vertices_per_face; ++j)
    deallog << " :: " << face_1->vertex(j);
  deallog << std::endl;

  deallog << "face 2";
  for (unsigned int j=0; j<GeometryInfo<dim>::vertices_per_face; ++j)
    deallog << " :: " << face_2->vertex(j);
  deallog << std::endl;

  deallog << "orientation: " << orientation[0]
          << "  flip: " << orientation[1]
          << "  rotation: " << orientation[2]
          << std::endl
          << std::endl;
}

int main()
{
  deallog << std::setprecision(4);
  logfile << std::setprecision(4);
  deallog.attach(logfile, false);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog << "Test for 2D: Hypercube" << std::endl << std::endl;

  // Generate a triangulation and match:
  Triangulation<2> triangulation2;

  generate_grid(triangulation2);

  typedef Triangulation<2>::cell_iterator CellIterator2;
  typedef std::vector<GridTools::PeriodicFacePair<CellIterator2> > FaceVector2;
  FaceVector2 test2;
  GridTools::collect_periodic_faces
    (triangulation2, 42, 1, test2, dealii::Tensor<1,2>());

  for (FaceVector2::iterator it = test2.begin(); it != test2.end(); ++it)
    print_match(it->cell[0]->face(it->face_idx[0]),
                it->cell[1]->face(it->face_idx[1]),
                it->orientation);

  deallog << "Test for 3D: Hypercube" << std::endl << std::endl;


  // Generate a triangulation and match:
  Triangulation<3> triangulation3;

  generate_grid(triangulation3);

  typedef Triangulation<3>::cell_iterator CellIterator3;
  typedef std::vector<GridTools::PeriodicFacePair<CellIterator3> > FaceVector3;
  FaceVector3 test3;
  GridTools::collect_periodic_faces
    (triangulation3, 42, 2, test3, dealii::Tensor<1,3>());

  for (FaceVector3::iterator it = test3.begin(); it != test3.end(); ++it)
    print_match(it->cell[0]->face(it->face_idx[0]),
                it->cell[1]->face(it->face_idx[1]),
                it->orientation);

  return 0;
}

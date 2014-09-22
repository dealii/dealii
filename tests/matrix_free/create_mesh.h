//-----------------------  create_mesh.h  -----------------------------
//    Version: $Name$
//
//-----------------------  create_mesh.h  -----------------------------


// this creates a mesh that contains cells of all different kinds detected in
// the MatrixFree class: Use a mesh that consists of a square cell, then a
// triangular cell, a parallelepiped cell and two trapezoidal cells, where one
// is very close to being Cartesian (by 1e-8).

#include <deal.II/base/point.h>
#include <deal.II/grid/tria.h>

#include <fstream>
#include <iostream>


void create_mesh (Triangulation<2> &tria,
                  const double scale_grid = 1.)
{
  const unsigned int dim = 2;
  std::vector<Point<dim> > points (12);

  // build the mesh layer by layer from points

  // 1. cube cell
  points[0] = Point<dim> (0, 0);
  points[1] = Point<dim> (0, 1);
  points[2] = Point<dim> (1 ,0);
  points[3] = Point<dim> (1 ,1);

  // 2. rectangular cell
  points[4] = Point<dim> (3., 0);
  points[5] = Point<dim> (3., 1);

  // 3. parallelogram cell
  points[6] = Point<dim> (5., 1.);
  points[7] = Point<dim> (5., 2.);

  // almost square cell (but trapezoidal by
  // 1e-8)
  points[8] = Point<dim> (6., 1.);
  points[9] = Point<dim> (6., 2.+1e-8);

  // apparently trapezoidal cell
  points[10] = Point<dim> (7., 1.4);
  points[11] = Point<dim> (7.5, numbers::PI);

  if (scale_grid != 1.)
    for (unsigned int i=0; i<points.size(); ++i)
      points[i] *= scale_grid;


  // connect the points to cells
  std::vector<CellData<dim> > cells(5);
  for (unsigned int i=0; i<5; ++i)
    {
      cells[i].vertices[0] = 0+2*i;
      cells[i].vertices[1] = 2+2*i;
      cells[i].vertices[2] = 1+2*i;
      cells[i].vertices[3] = 3+2*i;
      cells[i].material_id = 0;
    }

  tria.create_triangulation (points, cells, SubCellData());
}



void create_mesh (Triangulation<3> &tria,
                  const double scale_grid = 1.)
{
  const unsigned int dim = 3;
  std::vector<Point<dim> > points (24);

  // build the mesh layer by layer from points

  // 1. cube cell
  points[0] = Point<dim> (0,0,0);
  points[1] = Point<dim> (0,1.,0);
  points[2] = Point<dim> (0,0,1);
  points[3] = Point<dim> (0,1.,1);
  points[4] = Point<dim> (1.,0,0);
  points[5] = Point<dim> (1.,1.,0);
  points[6] = Point<dim> (1.,0,1);
  points[7] = Point<dim> (1.,1.,1);

  // 2. rectangular cell
  points[8] = Point<dim> (3., 0, 0);
  points[9] = Point<dim> (3., 1, 0);
  points[10] = Point<dim> (3., 0,1);
  points[11] = Point<dim> (3., 1,1);

  // 3. parallelogram cell
  points[12] = Point<dim> (5., 1., 1.);
  points[13] = Point<dim> (5., 2., 1.);
  points[14] = Point<dim> (5., 1., 2.);
  points[15] = Point<dim> (5., 2., 2.);

  // almost square cell (but trapezoidal by
  // 1e-8 in y-direction)
  points[16] = Point<dim> (6., 1., 1.);
  points[17] = Point<dim> (6., 2.+1e-8, 1.);
  points[18] = Point<dim> (6., 1., 2.);
  points[19] = Point<dim> (6., 2., 2.);

  // apparently trapezoidal cell
  points[20] = Point<dim> (7., 1.4, 1.2231);
  points[21] = Point<dim> (7.5, numbers::PI, 1.334);
  points[22] = Point<dim> (7., 1.5, 7.1);
  points[23] = Point<dim> (7.5, 3.8, 2.99);

  if (scale_grid != 1.)
    for (unsigned int i=0; i<points.size(); ++i)
      points[i] *= scale_grid;

  // connect the points to cells
  std::vector<CellData<dim> > cells(5);
  for (unsigned int i=0; i<5; ++i)
    {
      cells[i].vertices[0] = 0+4*i;
      cells[i].vertices[1] = 4+4*i;
      cells[i].vertices[2] = 1+4*i;
      cells[i].vertices[3] = 5+4*i;
      cells[i].vertices[4] = 2+4*i;
      cells[i].vertices[5] = 6+4*i;
      cells[i].vertices[6] = 3+4*i;
      cells[i].vertices[7] = 7+4*i;
      cells[i].material_id = 0;
    }
  tria.create_triangulation (points, cells, SubCellData());
}

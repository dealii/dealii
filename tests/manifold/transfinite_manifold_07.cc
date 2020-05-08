// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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


// This test verifies that we can use transfinite interpolation on 2D cells
// that have large aspect ratios: if the initial guess of a chart point is not
// accurate enough then we will fail to calculate the weighted average in
// chart coordinates. This test is based on a minimum working example provided
// by Juan Carlos Araujo Cabarcas which did not work with the changes made to
// transfinite interpolation in late 2017 but worked with an earlier version.

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <sstream>

#include "../tests.h"

//---------------------------- GRID DEFINITION ---------------------------------

struct Geom_parameters
{
  unsigned int n_balls;

  std::vector<Point<2>> ball_centers;
  std::vector<double>   radius;
};

void concentric_disks(Triangulation<2> &         tria,
                      const double               s,
                      const std::vector<double> &x,
                      Geom_parameters &          gp)
{
  double r = x[0], d = 0.5 * x[0],
         q = 1.0 / sqrt(2.0); // q: corner points factor

  const unsigned int nlay = x.size() - 1, vlayer = 8, lcells = 8,
                     n_vert = (2 + nlay) * vlayer + 1,
                     n_cell = 3 + (1 + nlay) * lcells + 1; // 19

  int fc = 4; // fv = 1,  f: fixed, displacement of index on the vertices and
  // cells ... inner objects

  // geometric information-----------------------
  gp.n_balls = x.size(); // the resonator and all other centered at the origin
  gp.ball_centers.resize(gp.n_balls);
  gp.radius.resize(gp.n_balls);
  gp.radius[0]       = x[0];
  gp.ball_centers[0] = Point<2>(s, 0.0);

  for (unsigned int k = 1; k < x.size(); k++)
    {
      gp.radius[k]       = x[k];
      gp.ball_centers[k] = Point<2>(0.0, 0.0);
    }

  std::vector<Point<2>> vertices(n_vert);

  vertices[0] = Point<2>(0.0 + s, 0.0);

  // inner 4 squares
  vertices[1] = Point<2>(d + s, 0.0);
  vertices[2] = Point<2>(.5 * d + s, .5 * d);
  vertices[3] = Point<2>(0.0 + s, d);
  vertices[4] = Point<2>(-.5 * d + s, .5 * d);
  vertices[5] = Point<2>(-d + s, 0.0);
  vertices[6] = Point<2>(-.5 * d + s, -.5 * d);
  vertices[7] = Point<2>(0.0 + s, -d);
  vertices[8] = Point<2>(.5 * d + s, -.5 * d);

  // touching circle
  for (unsigned int k = 0; k < x.size(); k++)
    {
      double z                     = (k == 0) ? s : 0.0;
      r                            = x[k];
      vertices[8 + k * vlayer + 1] = Point<2>(r + z, 0.0);
      vertices[8 + k * vlayer + 2] = Point<2>(r * q + z, r * q);
      vertices[8 + k * vlayer + 3] = Point<2>(0.0 + z, r);
      vertices[8 + k * vlayer + 4] = Point<2>(-r * q + z, r * q);
      vertices[8 + k * vlayer + 5] = Point<2>(-r + z, 0.0);
      vertices[8 + k * vlayer + 6] = Point<2>(-r * q + z, -r * q);
      vertices[8 + k * vlayer + 7] = Point<2>(0.0 + z, -r);
      vertices[8 + k * vlayer + 8] = Point<2>(r * q + z, -r * q);
    }

  std::vector<std::vector<int>> cell_v(n_cell, std::vector<int>(4));

  cell_v[0][0] = 0; // cell, i = 0
  cell_v[0][1] = 8;
  cell_v[0][2] = 2;
  cell_v[0][3] = 1;

  cell_v[1][0] = 0; // cell, i = 1
  cell_v[1][1] = 2;
  cell_v[1][2] = 4;
  cell_v[1][3] = 3;

  cell_v[2][0] = 0; // cell, i = 2
  cell_v[2][1] = 4;
  cell_v[2][2] = 6;
  cell_v[2][3] = 5;

  cell_v[3][0] = 0; // cell, i = 3
  cell_v[3][1] = 6;
  cell_v[3][2] = 8;
  cell_v[3][3] = 7;

  cell_v[4][0] = 8; // cell, i = 4
  cell_v[4][1] = 16;
  cell_v[4][2] = 1;
  cell_v[4][3] = 9;

  cell_v[5][0] = 2; // cell, i = 5
  cell_v[5][1] = 1;
  cell_v[5][2] = 10;
  cell_v[5][3] = 9;

  cell_v[6][0] = 2; // cell, i = 6
  cell_v[6][1] = 10;
  cell_v[6][2] = 3;
  cell_v[6][3] = 11;

  cell_v[7][0] = 4; // cell, i = 7
  cell_v[7][1] = 3;
  cell_v[7][2] = 12;
  cell_v[7][3] = 11;

  cell_v[8][0] = 4; // cell, i = 8
  cell_v[8][1] = 12;
  cell_v[8][2] = 5;
  cell_v[8][3] = 13;

  cell_v[9][0] = 6; // cell, i = 9
  cell_v[9][1] = 5;
  cell_v[9][2] = 14;
  cell_v[9][3] = 13;

  cell_v[10][0] = 6; // cell, i = 10
  cell_v[10][1] = 14;
  cell_v[10][2] = 7;
  cell_v[10][3] = 15;

  cell_v[11][0] = 8; // cell, i = 11
  cell_v[11][1] = 7;
  cell_v[11][2] = 16;
  cell_v[11][3] = 15;

  // layer cells
  for (unsigned int k = 1; k < x.size(); k++)
    {
      const unsigned int m = k + 1;

      // cell 12
      cell_v[fc + k * lcells + 0][0] = 7 + 8 * k + 1; // 16;  // cell, i = 12
      cell_v[fc + k * lcells + 0][1] = 7 + 8 * m + 1; // 24;
      cell_v[fc + k * lcells + 0][2] = 0 + 8 * k + 1; // 9;
      cell_v[fc + k * lcells + 0][3] = 0 + 8 * m + 1; // 17;

      // cell 13
      cell_v[fc + k * lcells + 1][0] = 1 + 8 * k + 1; // 10;  // cell, i = 13
      cell_v[fc + k * lcells + 1][1] = 0 + 8 * k + 1; // 9;
      cell_v[fc + k * lcells + 1][2] = 1 + 8 * m + 1; // 18;
      cell_v[fc + k * lcells + 1][3] = 0 + 8 * m + 1; // 17;

      // cell 14
      cell_v[fc + k * lcells + 2][0] = 1 + 8 * k + 1; // 10;  // cell, i = 14
      cell_v[fc + k * lcells + 2][1] = 1 + 8 * m + 1; // 18;
      cell_v[fc + k * lcells + 2][2] = 2 + 8 * k + 1; // 11;
      cell_v[fc + k * lcells + 2][3] = 2 + 8 * m + 1; // 19;

      // cell 15
      cell_v[fc + k * lcells + 3][0] = 3 + 8 * k + 1; // 12;  // cell, i = 15
      cell_v[fc + k * lcells + 3][1] = 2 + 8 * k + 1; // 11;
      cell_v[fc + k * lcells + 3][2] = 3 + 8 * m + 1; // 20;
      cell_v[fc + k * lcells + 3][3] = 2 + 8 * m + 1; // 19;

      // cell 16
      cell_v[fc + k * lcells + 4][0] = 3 + 8 * k + 1; // 12;  // cell, i = 16
      cell_v[fc + k * lcells + 4][1] = 3 + 8 * m + 1; // 20;
      cell_v[fc + k * lcells + 4][2] = 4 + 8 * k + 1; // 13;
      cell_v[fc + k * lcells + 4][3] = 4 + 8 * m + 1; // 21;

      // cell 17
      cell_v[fc + k * lcells + 5][0] = 5 + 8 * k + 1; // 14;  // cell, i = 17
      cell_v[fc + k * lcells + 5][1] = 4 + 8 * k + 1; // 13;
      cell_v[fc + k * lcells + 5][2] = 5 + 8 * m + 1; // 22;
      cell_v[fc + k * lcells + 5][3] = 4 + 8 * m + 1; // 21;

      // cell 18
      cell_v[fc + k * lcells + 6][0] = 5 + 8 * k + 1; // 14;  // cell, i = 18
      cell_v[fc + k * lcells + 6][1] = 5 + 8 * m + 1; // 22;
      cell_v[fc + k * lcells + 6][2] = 6 + 8 * k + 1; // 15;
      cell_v[fc + k * lcells + 6][3] = 6 + 8 * m + 1; // 23;

      // cell 19
      cell_v[fc + k * lcells + 7][0] = 7 + 8 * k + 1; // 16;  // cell, i =
      cell_v[fc + k * lcells + 7][1] = 6 + 8 * k + 1; // 15;
      cell_v[fc + k * lcells + 7][2] = 7 + 8 * m + 1; // 24;
      cell_v[fc + k * lcells + 7][3] = 6 + 8 * m + 1; // 23;
    }

  std::vector<CellData<2>> cells(n_cell, CellData<2>());

  unsigned int cell_idx = 0;
  for (unsigned int i = 0; i < n_cell; ++i)
    {
      for (unsigned int j = 0; j < 4; ++j)
        {
          cells[i].vertices[j] = cell_v[i][j];
        }
      cell_idx++;
    }

  tria.create_triangulation(vertices,
                            cells,
                            SubCellData()); // no boundary information

  double       eps   = 1e-5 * x[0];
  unsigned int label = 100;

  for (Triangulation<2>::active_cell_iterator cell = tria.begin_active();
       cell != tria.end();
       ++cell)
    {
      cell->set_all_manifold_ids(
        1); // not faces ... for Transfinite interpolation
      for (unsigned int k = 0; k < gp.n_balls; ++k)
        {
          for (const unsigned int f : GeometryInfo<2>::face_indices())
            {
              const Point<2> p0 = cell->face(f)->vertex(0),
                             p1 = cell->face(f)->vertex(1);
              const double d0   = p0.distance(gp.ball_centers[k]),
                           d1   = p1.distance(gp.ball_centers[k]);

              if ((std::abs(d0 - gp.radius[k]) < eps) &&
                  (std::abs(d1 - gp.radius[k]) < eps))
                {
                  cell->face(f)->set_all_manifold_ids(label + k);
                }
            }
        } // face
    }     // cell

  // end-colorizing
  // --------------------------------------------------------------------------------
}

void concentric_disks(Triangulation<2> &  tria,
                      std::vector<double> x,
                      Geom_parameters &   gp)
{
  concentric_disks(tria, 0.0, x, gp);
}

//--------------------- MAIN CLASS ---------------------------------------------

template <int dim>
class Mygrid
{
public:
  Mygrid(unsigned int r);
  void
  run();

private:
  void
  make_grid();

  Geom_parameters                       gp;
  std::vector<PolarManifold<dim>>       balls;
  TransfiniteInterpolationManifold<dim> inner_manifold;
  Triangulation<dim>                    triangulation;
  unsigned int                          refinement;
};

template <int dim>
Mygrid<dim>::Mygrid(unsigned int r)
  : refinement(r)
{}

template <int dim>
void
Mygrid<dim>::make_grid()
{
  const double        s = 0.1;
  std::vector<double> x{1.0, 1.5, 2.0, 2.5, 3.0};
  concentric_disks(triangulation, s, x, gp);
  for (unsigned int i = 0; i < gp.n_balls; i++)
    {
      balls.emplace_back(gp.ball_centers[i]);
    }

  // assigning manifolds, with layers 100, 101, 102, 103
  for (unsigned int i = 0; i < gp.n_balls; i++)
    {
      triangulation.set_manifold(100 + i, balls[i]);
    }

  const unsigned int not_curved = 1;
  inner_manifold.initialize(triangulation);
  triangulation.set_manifold(not_curved, inner_manifold);

  GridOut grid_out;
  grid_out.write_vtk(triangulation, deallog.get_file_stream());
  triangulation.refine_global(refinement);
  grid_out.write_vtk(triangulation, deallog.get_file_stream());
}

template <int dim>
void
Mygrid<dim>::run()
{
  make_grid();
}

int
main(int argc, char **argv)
{
  initlog();

  Mygrid<2> problem_2d(1);
  problem_2d.run();
}

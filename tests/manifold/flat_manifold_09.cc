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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// Test that Manifold::get_normal_vector also works for very small cells (9+
// refinement steps) instead of failing with
//
// "The Newton iteration to find the reference point did not converge in 10
// iterations"

// clang-format off
/*
Example output of a failing iteration:

  1 delta_xi: 0.707107 xi: 1 1 F: 0.449978 0.159888 0 p: 0.45 0.159912 0 delta: 0.5
  2 delta_xi: 1.26319e-12 xi: 1 1 F: 0.45 0.159912 0 p: 0.45 0.159912 0 delta: 8.45028e-13
  3 delta_xi: 1.26319e-12 xi: 1 1 F: 0.45 0.159912 0 p: 0.45 0.159912 0 delta: 8.45028e-13
  4 delta_xi: 1.26319e-12 xi: 1 1 F: 0.45 0.159912 0 p: 0.45 0.159912 0 delta: 8.45028e-13
  5 delta_xi: 1.26319e-12 xi: 1 1 F: 0.45 0.159912 0 p: 0.45 0.159912 0 delta: 8.45028e-13
  6 delta_xi: 1.26319e-12 xi: 1 1 F: 0.45 0.159912 0 p: 0.45 0.159912 0 delta: 8.45028e-13
  7 delta_xi: 1.26319e-12 xi: 1 1 F: 0.45 0.159912 0 p: 0.45 0.159912 0 delta: 8.45028e-13
  8 delta_xi: 1.26319e-12 xi: 1 1 F: 0.45 0.159912 0 p: 0.45 0.159912 0 delta: 8.45028e-13
  9 delta_xi: 1.26319e-12 xi: 1 1 F: 0.45 0.159912 0 p: 0.45 0.159912 0 delta: 8.45028e-13
  10 delta_xi: 1.26319e-12 xi: 1 1 F: 0.45 0.159912 0 p: 0.45 0.159912 0 delta: 8.45028e-13
*/
// clang-format on


#include "../tests.h"


// all include files you need here
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

template <int dim>
void
make_grid(Triangulation<dim> &triangulation)
{
  std::vector<std::vector<double>> spacing(3, std::vector<double>());

  for (int i = 0; i < 5; ++i)
    spacing[0].push_back(0.09);
  spacing[0].push_back(0.1);
  for (int i = 0; i < 10; ++i)
    spacing[0].push_back(0.095);
  for (int i = 0; i < 5; ++i)
    spacing[0].push_back(0.2);

  spacing[1].push_back(0.075);
  spacing[1].push_back(0.075);
  spacing[1].push_back(0.1);
  spacing[1].push_back(0.08);
  spacing[1].push_back(0.08);

  const unsigned int z_slices = 4;
  for (unsigned int z = 0; z < z_slices; ++z)
    spacing[2].push_back(0.41 / z_slices);

  Point<3>        p;
  TableIndices<3> sizes;
  sizes[0] = spacing[0].size();
  sizes[1] = spacing[1].size();
  sizes[2] = spacing[2].size();
  Table<3, types::material_id> material_id;
  material_id.TableBase<3, types::material_id>::reinit(sizes);

  for (unsigned int z = 0; z < z_slices; ++z)
    material_id(5, 2, z) = -1;

  GridGenerator::subdivided_hyper_rectangle(
    triangulation, spacing, p, material_id, false);
}

template <int dim>
void
check(Triangulation<dim> &tria)
{
  FE_Q<dim> fe(2);

  Quadrature<dim - 1> quadrature_formula(fe.get_unit_face_support_points());

  FEFaceValues<dim> fe_face_values(
    fe, quadrature_formula, update_normal_vectors | update_quadrature_points);

  Tensor<1, dim> n;
  for (const auto &cell : tria.active_cell_iterators())
    {
      for (const unsigned int face : GeometryInfo<dim>::face_indices())
        if (cell->face(face)->at_boundary())
          {
            fe_face_values.reinit(cell, face);
            for (unsigned int q_point = 0; q_point < quadrature_formula.size();
                 ++q_point)
              n += (cell->face(face)->get_manifold().normal_vector(
                cell->face(face), fe_face_values.quadrature_point(q_point)));
          }
    }
  std::cout << n << std::endl;
}

template <int dim>
void
test()
{
  Triangulation<dim> tria;
  make_grid(tria);
  tria.refine_global(1);


  Point<dim> p;
  p(0) = 0.45;
  p(1) = 0.16;
  for (unsigned int step = 0; step < 10; ++step)
    {
      std::cout << "step " << step << std::endl;
      for (auto &cell : tria.active_cell_iterators())
        {
          for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
            {
              if (p.distance(cell->vertex(v)) < 0.03 / (1 << step))
                {
                  cell->set_refine_flag();
                  break;
                }
            }
        }

      tria.execute_coarsening_and_refinement();

      check(tria);
    }


  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();


  test<3>();
  return 0;
}

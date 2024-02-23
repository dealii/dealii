// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Small test to analyse the equivalence from both sides of a face of
// the normal component for the Raviart-Thomas element
//
// Test adapted from a small program by Oliver Kayser-Herold.
//
// the test presently fails because of the issue with computing the
// normals using FEFaceValue, where FEFaceValue by accident uses the
// *face* mapping, not the *cell* mapping to compute the Piola
// transform (leading to a missing power of h in the determinant)


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q1_eulerian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


void
evaluate_normal(DoFHandler<2> &dof_handler, Vector<double> &solution)
{
  // This quadrature rule determines
  // the points, where the continuity
  // will be tested.
  QGauss<1>       quad(6);
  FEFaceValues<2> fe_v_face(
    dof_handler.get_fe(),
    quad,
    UpdateFlags(update_values | update_quadrature_points | update_gradients |
                update_normal_vectors | update_JxW_values));

  FEFaceValues<2> fe_v_face_n(
    dof_handler.get_fe(),
    quad,
    UpdateFlags(update_values | update_quadrature_points | update_gradients |
                update_normal_vectors | update_JxW_values));

  const unsigned int n_q_face     = quad.size();
  const unsigned int n_components = dof_handler.get_fe().n_components();

  // Cell iterators
  DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(),
                                      endc = dof_handler.end();

  for (; cell != endc; ++cell)
    {
      for (const unsigned int f : GeometryInfo<2>::face_indices())
        {
          if (!cell->face(f)->at_boundary())
            {
              deallog << "cell=" << cell << ", face=" << f << std::endl;

              fe_v_face.reinit(cell, f);

              const unsigned int neighbor = cell->neighbor_of_neighbor(f);
              fe_v_face_n.reinit(cell->neighbor(f), neighbor);

              // Get values from
              // solution vector (For
              // Trap.Rule)
              std::vector<Vector<double>> this_value(
                n_q_face, Vector<double>(n_components));
              fe_v_face.get_function_values(solution, this_value);

              // Same for neighbor cell
              std::vector<Vector<double>> this_value_n(
                n_q_face, Vector<double>(n_components));
              fe_v_face_n.get_function_values(solution, this_value_n);

              for (unsigned int q_point = 0; q_point < n_q_face; ++q_point)
                {
                  Tensor<1, 2> vn = fe_v_face.normal_vector(q_point);
                  double       nx = vn[0];
                  double       ny = vn[1];

                  double u = this_value[q_point](0);
                  double v = this_value[q_point](1);

                  double u_n = this_value_n[q_point](0);
                  double v_n = this_value_n[q_point](1);

                  deallog << "quadrature point " << q_point
                          << " (u-u_n)*nx+(v-v_n)*ny = "
                          << (u - u_n) * nx + (v - v_n) * ny
                          << ", u*nx+v*ny = " << u * nx + v * ny
                          << ", u_n*nx+v_n*ny = " << u_n * nx + v_n * ny
                          << std::endl;
                }
            }
        }
    }
}


int
main()
{
  initlog();

  Triangulation<2>          tria_test;
  Point<2>                  p1(0, 0), p2(1, 1);
  std::vector<unsigned int> sub_div;

  sub_div.push_back(1);
  sub_div.push_back(4);

  GridGenerator::subdivided_hyper_rectangle(tria_test, sub_div, p1, p2);
  tria_test.refine_global(2);

  FE_RaviartThomas<2> fe(1);
  DoFHandler<2>       dof_handler(tria_test);
  dof_handler.distribute_dofs(fe);

  // Fill solution vector with random
  // values between 0 and 1.
  Vector<double> solution(dof_handler.n_dofs());
  for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
    solution(i) = random_value<double>();

  // Now check if the function is
  // continuous in normal direction
  // on uniform mesh:
  deallog << "Uniform mesh test" << std::endl;
  evaluate_normal(dof_handler, solution);


  // Then test same on distorted mesh
  GridTools::distort_random(0.05, tria_test);
  deallog << "Distorted mesh test" << std::endl;
  evaluate_normal(dof_handler, solution);
}

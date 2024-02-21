// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/*
 * Small test to analyse the equivalence of the normal component
 * on the element edges for the Raviart-Thomas elements.
 */



#include "../tests.h"

#define PRECISION 8


#include <deal.II/base/function.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature.h>
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
#include <deal.II/grid/reference_cell.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

/*
 * Check if the normal component is continuous over element edges.
 */

void
evaluate_normal_component(const DoFHandler<2> &dof_handler,
                          Vector<double>      &solution)
{
  // This quadrature rule determines the points, where the
  // continuity will be tested.
  QGauss<1>     quad(6);
  Quadrature<2> qproject =
    QProjector<2>::project_to_all_faces(ReferenceCells::Quadrilateral, quad);

  FEFaceValues<2> fe_v_face(
    dof_handler.get_fe(),
    quad,
    UpdateFlags(update_values | update_quadrature_points | update_gradients |
                update_normal_vectors | update_JxW_values));

  FEValues<2> fe_v(dof_handler.get_fe(),
                   qproject,
                   UpdateFlags(update_values | update_quadrature_points |
                               update_gradients | update_JxW_values));

  FEValues<2> fe_v_n(dof_handler.get_fe(),
                     qproject,
                     UpdateFlags(update_values | update_quadrature_points |
                                 update_gradients | update_JxW_values));

  const unsigned int n_q_face      = quad.size();
  const unsigned int n_q_proj      = qproject.size();
  const unsigned int n_components  = dof_handler.get_fe().n_components();
  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

  deallog << "n_q_face=" << n_q_face << ", n_q_proj=" << n_q_proj << std::endl;

  // Cell iterators
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell->get_dof_indices(local_dof_indices);
      fe_v.reinit(cell);

      for (const unsigned int f : cell->face_indices())
        {
          if (!cell->face(f)->at_boundary())
            {
              deallog << "Testing face with center at "
                      << cell->face(f)->center() << std::endl;

              const QProjector<2>::DataSetDescriptor offset =
                (QProjector<2>::DataSetDescriptor::face(
                  ReferenceCells::Quadrilateral,
                  f,
                  cell->face_orientation(f),
                  cell->face_flip(f),
                  cell->face_rotation(f),
                  quad.size()));
              fe_v_face.reinit(cell, f);

              DoFHandler<2>::active_cell_iterator cell_n = cell->neighbor(f);
              const unsigned int neighbor = cell->neighbor_of_neighbor(f);
              fe_v_n.reinit(cell_n);

              const QProjector<2>::DataSetDescriptor offset_n =
                (QProjector<2>::DataSetDescriptor::face(
                  cell_n->reference_cell(),
                  neighbor,
                  cell_n->face_orientation(neighbor),
                  cell_n->face_flip(neighbor),
                  cell_n->face_rotation(neighbor),
                  quad.size()));

              // Get values from solution vector (For Trap.Rule)
              std::vector<Vector<double>> this_value(
                n_q_proj, Vector<double>(n_components));
              fe_v.get_function_values(solution, this_value);

              // Same for neighbor cell
              std::vector<Vector<double>> this_value_n(
                n_q_proj, Vector<double>(n_components));
              fe_v_n.get_function_values(solution, this_value_n);

              for (const auto q_point : fe_v_face.quadrature_point_indices())
                {
                  const Tensor<1, 2> n  = fe_v_face.normal_vector(q_point);
                  const double       nx = n[0];
                  const double       ny = n[1];

                  const Tensor<1, 2> u({this_value[q_point + offset](0),
                                        this_value[q_point + offset](1)});

                  const Tensor<1, 2> u_n({this_value_n[q_point + offset_n](0),
                                          this_value_n[q_point + offset_n](1)});

                  const double un1 = u * n;
                  const double un2 = u_n * n;

                  deallog << "  QP=" << q_point << ", error=" << (u - u_n) * n
                          << ", u.n=" << un1 << ", u_neighbor.n=" << un2
                          << std::endl;

                  Assert(std::fabs((u - u_n) * n) < 1e-12, ExcInternalError());
                }
            }
        }
    }
}



int
main()
{
  initlog();
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;

  Triangulation<2>                tria_test;
  const Point<2>                  p1(0, 0), p2(1, 1);
  const std::vector<unsigned int> sub_div = {1, 4};

  GridGenerator::subdivided_hyper_rectangle(tria_test, sub_div, p1, p2);
  tria_test.refine_global(2);
  GridTools::distort_random(0.05, tria_test);

  // Create a DoFHandler
  FE_RaviartThomas<2> fe(1);
  DoFHandler<2>       dof_handler(tria_test);
  dof_handler.distribute_dofs(fe);

  // Alloc some DoFs
  Vector<double> solution(dof_handler.n_dofs());

  // Fill solution vector with random values between 0 and 1.
  for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
    solution(i) = random_value<double>();

  // Now check if the function is continuous in normal
  // direction.
  evaluate_normal_component(dof_handler, solution);
}

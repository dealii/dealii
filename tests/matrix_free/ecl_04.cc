// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/tools.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


// Test FEFaceEvaluation::gather_evaluate() for values, ECL, shared-memory MPI,
// and different geometries (i.e., different orientations).

template <int dim>
class ExactSolution : public Function<dim>
{
public:
  ExactSolution()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int /*component*/ = 0) const
  {
    double result = p[0];

    for (unsigned i = 1; i < dim; ++i)
      result *= p[i];

    return result;
  }
};


template <int dim,
          int fe_degree,
          int n_points                 = fe_degree + 1,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
void
test(const unsigned int geometry, const MPI_Comm comm = MPI_COMM_SELF)
{
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  if (geometry == 0)
    {
      GridGenerator::hyper_cube(tria);

      if (dim == 2)
        tria.refine_global(4);
      else
        tria.refine_global(3);
    }
  else if (geometry == 1)
    {
      GridGenerator::hyper_shell(tria, Point<dim>(), 0.5, 1.0);

      if (dim == 2)
        tria.refine_global(4);
      else
        tria.refine_global(3);
    }
  else if (geometry == 2)
    {
      GridGenerator::hyper_ball(tria, Point<dim>(), 1.0);

      if (dim == 2)
        tria.refine_global(4);
      else
        tria.refine_global(3);
    }

  FE_DGQ<dim>     fe(fe_degree);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  MappingQ<dim> mapping(1);
  QGauss<1>     quad(n_points);

  AffineConstraints<Number> constraint;

  using MF = MatrixFree<dim, Number, VectorizedArrayType>;

  typename MF::AdditionalData additional_data;
  additional_data.mapping_update_flags = update_values | update_gradients;
  additional_data.mapping_update_flags_inner_faces =
    update_values | update_gradients;
  additional_data.mapping_update_flags_boundary_faces =
    update_values | update_gradients;
  additional_data.mapping_update_flags_faces_by_cells =
    update_values | update_gradients;
  additional_data.hold_all_faces_to_owned_cells = true;
  additional_data.communicator_sm               = comm;
  additional_data.tasks_parallel_scheme =
    MF::AdditionalData::TasksParallelScheme::none;

  MatrixFreeTools::categorize_by_boundary_ids(tria, additional_data);

  MF matrix_free;
  matrix_free.reinit(mapping, dof_handler, constraint, quad, additional_data);

  VectorType src, dst;

  matrix_free.initialize_dof_vector(src);
  matrix_free.initialize_dof_vector(dst);

  VectorTools::interpolate(dof_handler, ExactSolution<dim>(), src);

  dst = 0.0;

  /**
   * Element-centric loop
   */
  matrix_free.template loop_cell_centric<VectorType, VectorType>(
    [&](const auto &, auto &dst, const auto &src, const auto range) {
      FEFaceEvaluation<dim, fe_degree, n_points, 1, Number, VectorizedArrayType>
        phi_m(matrix_free, true);
      FEFaceEvaluation<dim, fe_degree, n_points, 1, Number, VectorizedArrayType>
        phi_p(matrix_free, false);

      for (unsigned int cell = range.first; cell < range.second; ++cell)
        {
          for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
               face++)
            {
              auto bids =
                matrix_free.get_faces_by_cells_boundary_id(cell, face);

              if (bids[0] != numbers::internal_face_boundary_id)
                continue;

              phi_m.reinit(cell, face);
              phi_p.reinit(cell, face);

              phi_m.gather_evaluate(src, EvaluationFlags::values);
              phi_p.gather_evaluate(src, EvaluationFlags::values);

              for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
                {
                  const auto u_minus = phi_m.get_value(q);
                  const auto u_plus  = phi_p.get_value(q);

                  for (unsigned int v = 0; v < VectorizedArray<double>::size();
                       ++v)
                    {
                      Assert(std::abs(u_minus[v] - u_plus[v]) < 1e-10,
                             ExcMessage("Entries do not match!"));
                    }
                }

              phi_m.gather_evaluate(src,
                                    EvaluationFlags::values |
                                      EvaluationFlags::gradients);
              phi_p.gather_evaluate(src,
                                    EvaluationFlags::values |
                                      EvaluationFlags::gradients);

              for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
                {
                  const auto u_minus = phi_m.get_value(q);
                  const auto u_plus  = phi_p.get_value(q);

                  const auto grad_u_minus = phi_m.get_gradient(q);
                  const auto grad_u_plus  = phi_p.get_gradient(q);

                  for (unsigned int v = 0; v < VectorizedArray<double>::size();
                       ++v)
                    {
                      Assert(std::abs(u_minus[v] - u_plus[v]) < 1e-10,
                             ExcMessage("Entries do not match!"));

                      if (false)
                        for (int d = 0; d < dim; ++d)
                          Assert(std::abs(grad_u_minus[d][v] -
                                          grad_u_plus[d][v]) < 1e-6,
                                 ExcMessage("Entries do not match!"));
                    }
                }
            }
        }
    },
    dst,
    src,
    false,
    MF::DataAccessOnFaces::values);

  deallog << "OK" << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  mpi_initlog();

  const auto rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  MPI_Comm subcommunicator;
  MPI_Comm_split_type(MPI_COMM_WORLD,
                      MPI_COMM_TYPE_SHARED,
                      rank,
                      MPI_INFO_NULL,
                      &subcommunicator);

  for (unsigned int i = 0; i < 3; ++i)
    test<2, 3, 4, double>(i, subcommunicator);

  for (unsigned int i = 0; i < 3; ++i)
    test<3, 3, 4, double>(i, subcommunicator);

  MPI_Comm_free(&subcommunicator);
}

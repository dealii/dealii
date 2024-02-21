// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/tools.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

// Test MatrixFreeTools::ElementActivationAndDeactivationMatrixFree by comparing
// the results with results obtained on fitted mesh.


template <int dim>
class Fu : public Function<dim>
{
public:
  double
  value(const Point<dim> &point, const unsigned int = 0) const
  {
    return std::sin(point[0] * 2 * numbers::PI) *
           std::sin(point[1] * 2 * numbers::PI);
  }
};

template <int dim>
void
test(const unsigned int n_refinements)
{
  using Number              = double;
  using VectorizedArrayType = VectorizedArray<Number>;
  using VectorType          = LinearAlgebra::distributed::Vector<Number>;
  using FECellIntegrator =
    FEEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType>;
  using FEFaceIntegrator =
    FEFaceEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType>;

  for (unsigned int i = 0; i < 2; ++i)
    {
      parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
      if (i == 0)
        GridGenerator::subdivided_hyper_rectangle(
          tria, {2, 2}, {0.0, 0.0}, {1.0, 1.0}, true);
      else
        GridGenerator::subdivided_hyper_rectangle(
          tria, {2, 1}, {0.0, 0.0}, {1.0, 0.5}, true);

      tria.refine_global(n_refinements);

      DoFHandler<dim> dof_handler(tria);

      for (const auto &cell :
           filter_iterators(dof_handler.active_cell_iterators(),
                            IteratorFilters::LocallyOwnedCell()))
        {
          if (cell->center()[1] < 0.5)
            cell->set_active_fe_index(0);
          else
            cell->set_active_fe_index(1);
        }

      dof_handler.distribute_dofs(
        hp::FECollection<dim>(FE_Q<dim>(1), FE_Nothing<dim>(1)));

      typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
        data;
      data.tasks_parallel_scheme =
        MatrixFree<dim, double>::AdditionalData::none;
      data.mapping_update_flags                = update_gradients;
      data.mapping_update_flags_boundary_faces = update_gradients;
      data.mapping_update_flags_inner_faces    = update_gradients;

      MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
      matrix_free.reinit(MappingQ1<dim>(),
                         dof_handler,
                         AffineConstraints<Number>(),
                         QGauss<dim>(2),
                         data);

      VectorType dst, src;
      matrix_free.initialize_dof_vector(dst);
      matrix_free.initialize_dof_vector(src);

      VectorTools::interpolate(dof_handler, Fu<dim>(), src);

      const auto cell_operation = [&](const auto &matrix_free,
                                      auto       &dst,
                                      const auto &src,
                                      const auto  range) {
        FECellIntegrator phi(matrix_free);

        for (unsigned int cell = range.first; cell < range.second; ++cell)
          {
            phi.reinit(cell);
            phi.gather_evaluate(src, EvaluationFlags::gradients);

            for (const auto q : phi.quadrature_point_indices())
              phi.submit_gradient(phi.get_gradient(q), q);

            phi.integrate_scatter(EvaluationFlags::gradients, dst);
          }
      };

      const auto face_operation =
        [&](const auto &, auto &, const auto &, const auto) {
          // nothing to do here but in the DG case
        };

      const auto boundary_operation = [&](const auto &matrix_free,
                                          auto &,
                                          const auto &,
                                          const auto range,
                                          const bool is_interior_face = true) {
        FEFaceIntegrator phi(matrix_free, is_interior_face);

        for (unsigned int face = range.first; face < range.second; ++face)
          {
            phi.reinit(face);

            const double scaling =
              phi.at_boundary() ? (phi.boundary_id() + 1.0) : (2.0 * dim);

            phi.gather_evaluate(src, EvaluationFlags::gradients);

            for (const auto q : phi.quadrature_point_indices())
              {
                auto gradient = phi.get_gradient(q);
                auto normal   = phi.normal_vector(q);

                if (is_interior_face == false) // fix sign!
                  normal *= -1.0;

                phi.submit_value(gradient * normal * scaling, q);
              }

            phi.integrate_scatter(EvaluationFlags::values, dst);
          }
      };

      if (i == 0)
        {
          MatrixFreeTools::ElementActivationAndDeactivationMatrixFree<
            dim,
            Number,
            VectorizedArrayType>
            element_birth_and_death;

          element_birth_and_death.reinit(matrix_free);

          element_birth_and_death.template cell_loop<VectorType, VectorType>(
            cell_operation, dst, src, true);

          deallog << dst.l2_norm() << std::endl;

          element_birth_and_death.template loop<VectorType, VectorType>(
            cell_operation, face_operation, boundary_operation, dst, src, true);

          deallog << dst.l2_norm() << std::endl;
        }
      else
        {
          matrix_free.template cell_loop<VectorType, VectorType>(cell_operation,
                                                                 dst,
                                                                 src,
                                                                 true);

          deallog << dst.l2_norm() << std::endl;

          matrix_free.template loop<VectorType, VectorType>(
            cell_operation, face_operation, boundary_operation, dst, src, true);

          deallog << dst.l2_norm() << std::endl;
        }
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  mpi_initlog();
  test<2>(2);
}

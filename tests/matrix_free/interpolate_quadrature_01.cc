// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
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

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



// Tests internal::FEFaceNormalEvaluationImpl::interpolate_quadrature()
// by comparing the results of FEFaceEvaluation::gather_evaluate().

template <int dim>
class ExactSolution : public Function<dim>
{
public:
  ExactSolution()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int /*component*/ = 0) const
  {
    return p[0] + p[1];
  }
};

template <int dim,
          int fe_degree,
          int n_points                 = fe_degree + 1,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
void
test(const unsigned int n_refinements = 1)
{
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0, 1, true);

  tria.refine_global(n_refinements);

  FE_DGQ<dim>     fe(fe_degree);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  MappingQ<dim> mapping(1);
  QGauss<1>     quad(n_points);

  AffineConstraints<Number> constraint;

  using MF = MatrixFree<dim, Number, VectorizedArrayType>;

  typename MF::AdditionalData additional_data;
  additional_data.mapping_update_flags                = update_values;
  additional_data.mapping_update_flags_inner_faces    = update_values;
  additional_data.mapping_update_flags_boundary_faces = update_values;
  additional_data.mapping_update_flags_faces_by_cells = update_values;
  additional_data.hold_all_faces_to_owned_cells       = true;

  MF matrix_free;
  matrix_free.reinit(mapping, dof_handler, constraint, quad, additional_data);

  VectorType src, dst;

  matrix_free.initialize_dof_vector(src);
  matrix_free.initialize_dof_vector(dst);

  VectorTools::interpolate(dof_handler, ExactSolution<dim>(), src);

  FEEvaluation<dim, fe_degree, n_points, 1, Number, VectorizedArrayType> phi(
    matrix_free);
  FEFaceEvaluation<dim, fe_degree, n_points, 1, Number, VectorizedArrayType>
    phi_m(matrix_free, true);

  matrix_free.template loop_cell_centric<VectorType, VectorType>(
    [&](const auto &, auto &, const auto &src, const auto range) {
      for (unsigned int cell = range.first; cell < range.second; ++cell)
        {
          phi.reinit(cell);

          phi.gather_evaluate(src, EvaluationFlags::values);

          for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
               face++)
            {
              phi_m.reinit(cell, face);
              phi_m.gather_evaluate(src, EvaluationFlags::values);

              AlignedVector<VectorizedArrayType> temp(
                phi_m.static_dofs_per_cell);

              internal::FEFaceNormalEvaluationImpl<dim,
                                                   n_points - 1,
                                                   VectorizedArrayType>::
                template interpolate_quadrature<true, false>(
                  1,
                  EvaluationFlags::values,
                  matrix_free.get_shape_info(),
                  phi.begin_values(),
                  temp.data(),
                  face);

              for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
                {
                  const auto u_cell = temp[q];
                  const auto u_face = phi_m.get_value(q);

                  for (unsigned int v = 0; v < VectorizedArray<double>::size();
                       ++v)
                    {
                      Assert(std::abs(u_cell[v] - u_face[v]) < 1e-10,
                             ExcMessage("Entries do not match!"));
                    }
                }
            }
        }
    },
    dst,
    src);
}

int
main()
{
  initlog();
  test<2, 1, 2, double, VectorizedArray<double>>();
  test<3, 1, 2, double, VectorizedArray<double>>();

  deallog << "OK!" << std::endl;
}

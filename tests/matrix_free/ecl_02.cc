// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


// Test FEFaceEvaluation::read_dof_values() and
// FEFaceEvaluation::gather_evaluate() for ECL for two cells.
//
// @note Since this program assumes that both cells are within the same
//   macro cell, this test is only run if vectorization is enabled.

template <int dim,
          int fe_degree,
          int n_points                 = fe_degree + 1,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
void
test(const unsigned int n_refinements = 1)
{
  if (VectorizedArrayType::size() == 1)
    return;

  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_rectangle(
    tria, {2, 1}, Point<dim>(0.0, 0.0), Point<dim>(2.0, 1.0), true);

  std::vector<dealii::GridTools::PeriodicFacePair<
    typename dealii::Triangulation<dim>::cell_iterator>>
    periodic_faces;

  if (dim >= 1)
    dealii::GridTools::collect_periodic_faces(tria, 0, 1, 0, periodic_faces);

  if (dim >= 2)
    dealii::GridTools::collect_periodic_faces(tria, 2, 3, 1, periodic_faces);

  tria.add_periodicity(periodic_faces);

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

  FEEvaluation<dim, fe_degree, n_points, 1, Number, VectorizedArrayType> phi(
    matrix_free);
  FEFaceEvaluation<dim, fe_degree, n_points, 1, Number, VectorizedArrayType>
    phi_m(matrix_free, true);
  FEFaceEvaluation<dim, fe_degree, n_points, 1, Number, VectorizedArrayType>
    phi_p(matrix_free, false);


  for (unsigned int i = 0; i < src.size() / 2; i++)
    src[i] = 1;

  for (unsigned int i = src.size() / 2; i < src.size(); i++)
    src[i] = 2;

  dst = 0.0;

  /**
   * Element-centric loop
   */
  matrix_free.template loop_cell_centric<VectorType, VectorType>(
    [&](const auto &, auto &dst, const auto &src, const auto range) {
      for (unsigned int cell = range.first; cell < range.second; ++cell)
        {
          for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
               face++)
            {
              phi_m.reinit(cell, face);
              phi_p.reinit(cell, face);

              phi_m.read_dof_values(src);

              for (unsigned int i = 0; i < phi_m.static_dofs_per_component; i++)
                deallog << static_cast<int>(phi_m.begin_dof_values()[i][0])
                        << " ";
              deallog << std::endl;
              phi_m.gather_evaluate(src, EvaluationFlags::values);
              for (unsigned int i = 0; i < phi_m.static_n_q_points; i++)
                deallog << static_cast<int>(phi_m.begin_values()[i][0]) << " ";
              deallog << std::endl;

              phi_p.read_dof_values(src);
              for (unsigned int i = 0; i < phi_p.static_dofs_per_component; i++)
                deallog << static_cast<int>(phi_p.begin_dof_values()[i][0])
                        << " ";
              deallog << std::endl;

              phi_p.gather_evaluate(src, EvaluationFlags::values);
              for (unsigned int i = 0; i < phi_p.static_n_q_points; i++)
                deallog << static_cast<int>(phi_p.begin_values()[i][0]) << " ";
              deallog << std::endl << std::endl;
            }
        }
    },
    dst,
    src,
    false,
    MF::DataAccessOnFaces::gradients);
}

int
main()
{
  initlog();
  test<2,
       1,
       2,
       double,
       VectorizedArray<double,
                       VectorizedArray<double>::size() < 2 ?
                         VectorizedArray<double>::size() :
                         2>>();
}

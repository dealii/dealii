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


// tests matrix-free read_cell_data for locally refined mesh

template <int dim, typename Number>
class SinusFunction : public Function<dim, Number>
{
public:
  Number
  value(const Point<dim> &p, const unsigned int component = 0) const
  {
    (void)component;
    return std::cos(2 * numbers::PI * p[0]);
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
  GridGenerator::hyper_cube(tria, 0.0, 1.0, true);

  std::vector<dealii::GridTools::PeriodicFacePair<
    typename dealii::Triangulation<dim>::cell_iterator>>
    periodic_faces;

  if (dim >= 1)
    dealii::GridTools::collect_periodic_faces(tria, 0, 1, 0, periodic_faces);

  if (dim >= 2)
    dealii::GridTools::collect_periodic_faces(tria, 2, 3, 1, periodic_faces);

  if (dim >= 3)
    dealii::GridTools::collect_periodic_faces(tria, 4, 5, 2, periodic_faces);

  tria.add_periodicity(periodic_faces);

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

  FEEvaluation<dim, fe_degree, n_points, 1, Number, VectorizedArrayType> phi(
    matrix_free);
  FEFaceEvaluation<dim, fe_degree, n_points, 1, Number, VectorizedArrayType>
    phi_m(matrix_free, true);
  FEFaceEvaluation<dim, fe_degree, n_points, 1, Number, VectorizedArrayType>
    phi_p(matrix_free, false);


  SinusFunction<dim, Number> function;
  VectorTools::interpolate(dof_handler, function, src);

  dst = 0.0;

  /**
   * Face-centric loop
   */
  matrix_free.template loop<VectorType, VectorType>(
    [&](const auto &, auto &dst, const auto &src, const auto range) {
      for (unsigned int cell = range.first; cell < range.second; ++cell)
        {
          phi.reinit(cell);
          phi.read_dof_values(src);
          phi.evaluate(false, true, false);
          for (unsigned int q = 0; q < phi.n_q_points; ++q)
            phi.submit_gradient(phi.get_gradient(q), q);
          phi.integrate(false, true);
          phi.set_dof_values(dst);
        }
    },
    [&](const auto &, auto &dst, const auto &src, const auto range) {
      for (unsigned int face = range.first; face < range.second; face++)
        {
          phi_m.reinit(face);
          phi_p.reinit(face);

          phi_m.read_dof_values(src);
          phi_m.evaluate(true, true);
          phi_p.read_dof_values(src);
          phi_p.evaluate(true, true);
          VectorizedArrayType sigmaF =
            (std::abs((phi_m.get_normal_vector(0) *
                       phi_m.inverse_jacobian(0))[dim - 1]) +
             std::abs((phi_m.get_normal_vector(0) *
                       phi_p.inverse_jacobian(0))[dim - 1])) *
            (Number)(std::max(fe_degree, 1) * (fe_degree + 1.0));

          for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
            {
              VectorizedArrayType average_value =
                (phi_m.get_value(q) - phi_p.get_value(q)) * 0.5;
              VectorizedArrayType average_valgrad =
                phi_m.get_normal_derivative(q) + phi_p.get_normal_derivative(q);
              average_valgrad =
                average_value * 2. * sigmaF - average_valgrad * 0.5;
              phi_m.submit_normal_derivative(-average_value, q);
              phi_p.submit_normal_derivative(-average_value, q);
              phi_m.submit_value(average_valgrad, q);
              phi_p.submit_value(-average_valgrad, q);
            }
          phi_m.integrate(true, true);
          phi_m.distribute_local_to_global(dst);
          phi_p.integrate(true, true);
          phi_p.distribute_local_to_global(dst);
        }
    },
    [&](const auto &, auto &, const auto &, const auto) {
      Assert(false, StandardExceptions::ExcMessage("Periodic BC!"));
    },
    dst,
    src,
    false,
    MF::DataAccessOnFaces::gradients,
    MF::DataAccessOnFaces::gradients);

  deallog << dst.l2_norm() << std::endl;
  dst.print(deallog.get_file_stream());

  dst = 0.0;

  /**
   * Element-centric loop
   */
  matrix_free.template loop_cell_centric<VectorType, VectorType>(
    [&](const auto &, auto &dst, const auto &src, const auto range) {
      for (unsigned int cell = range.first; cell < range.second; ++cell)
        {
          phi.reinit(cell);
          phi.read_dof_values(src);
          phi.evaluate(false, true, false);
          for (unsigned int q = 0; q < phi.n_q_points; ++q)
            phi.submit_gradient(phi.get_gradient(q), q);
          phi.integrate(false, true);

          for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
               face++)
            {
              auto bids =
                matrix_free.get_faces_by_cells_boundary_id(cell, face);

              // make sure that no boundaries are used
              AssertThrow(bids[0] == numbers::internal_face_boundary_id,
                          StandardExceptions::ExcMessage("Periodic BC!"));

              // make sure that all bids are the same
              for (unsigned int v = 1;
                   v < matrix_free.n_active_entries_per_cell_batch(cell);
                   v++)
                AssertThrow(bids[0] == bids[v],
                            ExcDimensionMismatch(bids[0], bids[v]));


              phi_m.reinit(cell, face);
              phi_p.reinit(cell, face);

              phi_m.read_dof_values(src);
              phi_m.evaluate(true, true);
              phi_p.read_dof_values(src);
              phi_p.evaluate(true, true);

              VectorizedArrayType sigmaF =
                (std::abs((phi_m.get_normal_vector(0) *
                           phi_m.inverse_jacobian(0))[dim - 1]) +
                 std::abs((phi_m.get_normal_vector(0) *
                           phi_p.inverse_jacobian(0))[dim - 1])) *
                (Number)(std::max(fe_degree, 1) * (fe_degree + 1.0));

              for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
                {
                  VectorizedArrayType average_value =
                    (phi_m.get_value(q) - phi_p.get_value(q)) * 0.5;
                  VectorizedArrayType average_valgrad =
                    phi_m.get_normal_derivative(q) +
                    phi_p.get_normal_derivative(q);
                  average_valgrad =
                    average_value * 2. * sigmaF - average_valgrad * 0.5;
                  phi_m.submit_normal_derivative(-average_value, q);
                  phi_m.submit_value(average_valgrad, q);
                }
              phi_m.integrate(true, true);

              for (unsigned int q = 0; q < phi.dofs_per_cell; ++q)
                phi.begin_dof_values()[q] += phi_m.begin_dof_values()[q];
            }

          phi.set_dof_values(dst);
        }
    },
    dst,
    src,
    false,
    MF::DataAccessOnFaces::gradients);

  deallog << dst.l2_norm() << std::endl;
  dst.print(deallog.get_file_stream());
}

int
main()
{
  initlog();
  test<2, 2, 3, double, VectorizedArray<double>>();
}
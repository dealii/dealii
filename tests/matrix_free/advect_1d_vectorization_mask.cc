// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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



// similar to  advect_1d but testing the mask for distribute_local_to_global

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim,
          int fe_degree,
          int n_q_points_1d   = fe_degree + 1,
          typename number     = double,
          typename VectorType = Vector<number>,
          int n_components    = 1>
class MatrixFreeAdvectionBasic
{
public:
  MatrixFreeAdvectionBasic(const MatrixFree<dim, number> &data,
                           const bool         zero_within_loop       = true,
                           const unsigned int start_vector_component = 0)
    : data(data)
    , zero_within_loop(zero_within_loop)
    , start_vector_component(start_vector_component)
  {
    for (unsigned int d = 0; d < dim; ++d)
      advection[d] = 0.4 + 0.12 * d;
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    if (!zero_within_loop)
      dst = 0;
    data.loop(&MatrixFreeAdvectionBasic::local_apply,
              &MatrixFreeAdvectionBasic::local_apply_face,
              &MatrixFreeAdvectionBasic::local_apply_boundary_face,
              this,
              dst,
              src,
              zero_within_loop,
              MatrixFree<dim, number>::DataAccessOnFaces::values,
              MatrixFree<dim, number>::DataAccessOnFaces::values);

    FEEvaluation<dim, fe_degree, fe_degree + 1, n_components, number> phi(data);

    const unsigned int dofs_per_cell = phi.dofs_per_cell;

    AlignedVector<VectorizedArray<number>> coefficients(phi.dofs_per_cell);
    MatrixFreeOperators::
      CellwiseInverseMassMatrix<dim, fe_degree, n_components, number>
        inverse(phi);

    for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values(dst);

        inverse.fill_inverse_JxW_values(coefficients);
        inverse.apply(coefficients,
                      n_components,
                      phi.begin_dof_values(),
                      phi.begin_dof_values());

        phi.set_dof_values(dst);
      }
  }

private:
  void
  local_apply(const MatrixFree<dim, number> &              data,
              VectorType &                                 dst,
              const VectorType &                           src,
              const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, n_q_points_1d, n_components, number> phi(
      data, 0, 0, start_vector_component);

    const unsigned int n_vect = VectorizedArray<number>::size();

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values(src);
        phi.evaluate(true, false);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_gradient(advection * phi.get_value(q), q);
        phi.integrate(false, true);
        for (unsigned int v = 0; v < n_vect; ++v)
          {
            std::bitset<n_vect> mask;
            mask[v] = true;
            phi.distribute_local_to_global(dst, 0, mask);
          }
      }
  }

  void
  local_apply_face(
    const MatrixFree<dim, number> &              data,
    VectorType &                                 dst,
    const VectorType &                           src,
    const std::pair<unsigned int, unsigned int> &face_range) const
  {
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, n_components, number> phi_m(
      data, true, 0, 0, start_vector_component);
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, n_components, number> phi_p(
      data, false, 0, 0, start_vector_component);
    typedef typename FEFaceEvaluation<dim,
                                      fe_degree,
                                      n_q_points_1d,
                                      n_components,
                                      number>::value_type value_type;

    const unsigned int n_vect = VectorizedArray<number>::size();

    for (unsigned int face = face_range.first; face < face_range.second; face++)
      {
        phi_m.reinit(face);
        phi_m.read_dof_values(src);
        phi_m.evaluate(true, false);
        phi_p.reinit(face);
        phi_p.read_dof_values(src);
        phi_p.evaluate(true, false);

        for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
          {
            value_type u_minus = phi_m.get_value(q),
                       u_plus  = phi_p.get_value(q);
            const VectorizedArray<number> normal_times_advection =
              advection * phi_m.get_normal_vector(q);
            const value_type flux_times_normal =
              make_vectorized_array<number>(0.5) *
              ((u_minus + u_plus) * normal_times_advection +
               std::abs(normal_times_advection) * (u_minus - u_plus));
            phi_m.submit_value(-flux_times_normal, q);
            phi_p.submit_value(flux_times_normal, q);
          }

        phi_m.integrate(true, false);
        for (unsigned int v = 0; v < n_vect; ++v)
          {
            std::bitset<n_vect> mask;
            mask[v] = true;
            phi_m.distribute_local_to_global(dst, 0, mask);
          }
        phi_p.integrate(true, false);
        for (unsigned int v = 0; v < n_vect; ++v)
          {
            std::bitset<n_vect> mask;
            mask[v] = true;
            phi_p.distribute_local_to_global(dst, 0, mask);
          }
      }
  }

  void
  local_apply_boundary_face(
    const MatrixFree<dim, number> &              data,
    VectorType &                                 dst,
    const VectorType &                           src,
    const std::pair<unsigned int, unsigned int> &face_range) const
  {
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, n_components, number>
                                                          fe_eval(data, true, 0, 0, start_vector_component);
    typedef typename FEFaceEvaluation<dim,
                                      fe_degree,
                                      n_q_points_1d,
                                      n_components,
                                      number>::value_type value_type;

    const unsigned int n_vect = VectorizedArray<number>::size();

    for (unsigned int face = face_range.first; face < face_range.second; face++)
      {
        fe_eval.reinit(face);
        fe_eval.read_dof_values(src);
        fe_eval.evaluate(true, false);

        for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
          {
            value_type                    u_minus = fe_eval.get_value(q);
            value_type                    u_plus  = -u_minus;
            const VectorizedArray<number> normal_times_advection =
              advection * fe_eval.get_normal_vector(q);
            const value_type flux_times_normal =
              make_vectorized_array<number>(0.5) *
              ((u_minus + u_plus) * normal_times_advection +
               std::abs(normal_times_advection) * (u_minus - u_plus));
            fe_eval.submit_value(-flux_times_normal, q);
          }

        fe_eval.integrate(true, false);
        for (unsigned int v = 0; v < n_vect; ++v)
          {
            std::bitset<n_vect> mask =
              std::bitset<VectorizedArray<number>::size()>();
            mask[v] = true;
            fe_eval.distribute_local_to_global(dst, 0, mask);
          }
      }
  }

  const MatrixFree<dim, number> &         data;
  const bool                              zero_within_loop;
  const unsigned int                      start_vector_component;
  Tensor<1, dim, VectorizedArray<number>> advection;
};



template <int dim>
class AnalyticFunction : public Function<dim>
{
public:
  static_assert(dim == 1, "Only 1D implemented");
  AnalyticFunction()
    : Function<dim>(1)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int) const override
  {
    return std::sin(3 * numbers::PI * p[0] / 0.8);
  }
};



template <int dim>
class AnalyticDerivative : public Function<dim>
{
public:
  static_assert(dim == 1, "Only 1D implemented");
  AnalyticDerivative()
    : Function<dim>(1)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int) const override
  {
    Tensor<1, dim> advection;
    for (unsigned int d = 0; d < dim; ++d)
      advection[d] = 0.4 + 0.12 * d;

    return -std::cos(3 * numbers::PI * p[0] / 0.8) * advection[0] * 3 *
           numbers::PI / 0.8;
  }
};



template <int dim, int fe_degree>
void
test(const unsigned int n_refine)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0, 0.8);
  tria.refine_global(n_refine);

  FE_DGQ<dim>     fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  if (n_refine == 3)
    {
      deallog << "Testing " << dof.get_fe().get_name();
      deallog << std::endl;
    }

  LinearAlgebra::distributed::Vector<double> in, out;

  const QGauss<1>                                  quad(fe_degree + 1);
  typename MatrixFree<dim, double>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim, double>::AdditionalData::none;
  data.mapping_update_flags_inner_faces =
    (update_gradients | update_JxW_values);
  data.mapping_update_flags_boundary_faces =
    (update_gradients | update_JxW_values);

  MatrixFree<dim, double> mf_data;
  mf_data.reinit(dof, constraints, quad, data);

  mf_data.initialize_dof_vector(in);
  mf_data.initialize_dof_vector(out);

  VectorTools::interpolate(dof, AnalyticFunction<dim>(), in);

  MatrixFreeAdvectionBasic<dim,
                           fe_degree,
                           fe_degree + 1,
                           double,
                           LinearAlgebra::distributed::Vector<double>>
    mf2(mf_data);
  mf2.vmult(out, in);

  VectorTools::interpolate(dof, AnalyticDerivative<dim>(), in);
  out -= in;

  double diff_norm = out.linfty_norm();
  deallog << "Norm of difference:          " << diff_norm << " " << std::endl;
}



int
main()
{
  initlog();

  for (unsigned int r = 3; r < 9; ++r)
    test<1, 2>(r);

  for (unsigned int r = 3; r < 9; ++r)
    test<1, 4>(r);

  for (unsigned int r = 3; r < 9; ++r)
    test<1, 5>(r);
}

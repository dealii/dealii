// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that FEInterfaceViews::get_function_*_from_local_dof_values
// works with different number types.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <fstream>
#include <type_traits>

#include "../tests.h"

template <typename NumberType, int dim, typename ExtractorType>
void
test_feiv_view(const Vector<double>          &solution,
               const FEInterfaceValues<dim>  &fe_iv,
               const unsigned int            &n_q_points,
               const ExtractorType           &extractor,
               const std::vector<NumberType> &local_dof_values);

// Scalar view
template <typename NumberType, int dim>
void
test_feiv_view(const Vector<double>             &solution,
               const FEInterfaceValues<dim>     &fe_iv,
               const unsigned int               &n_q_points,
               const FEValuesExtractors::Scalar &extractor,
               const std::vector<NumberType>    &local_dof_values)
{
  using View =
    std::remove_reference_t<std::remove_const_t<decltype(fe_iv[extractor])>>;
  const View &fe_iv_view = fe_iv[extractor];

  // Typedefs
  using value_type =
    typename ProductType<typename View::value_type, NumberType>::type;
  using gradient_type =
    typename ProductType<typename View::gradient_type, NumberType>::type;
  using hessian_type =
    typename ProductType<typename View::hessian_type, NumberType>::type;
  using laplacian_type =
    typename ProductType<typename View::value_type, NumberType>::type;
  using third_derivative_type =
    typename ProductType<typename View::third_derivative_type,
                         NumberType>::type;

  // Values
  std::vector<typename View::template solution_value_type<NumberType>>
    qp_jumps_local(n_q_points), qp_averages_local(n_q_points),
    qp_interface_values_local(n_q_points);
  std::vector<value_type> qp_jumps_global(n_q_points),
    qp_averages_global(n_q_points), qp_interface_values_global(n_q_points);

  fe_iv_view.get_jump_in_function_values_from_local_dof_values(local_dof_values,
                                                               qp_jumps_local);
  fe_iv_view.get_jump_in_function_values(solution, qp_jumps_global);

  fe_iv_view.get_average_of_function_values_from_local_dof_values(
    local_dof_values, qp_averages_local);
  fe_iv_view.get_average_of_function_values(solution, qp_averages_global);

  fe_iv_view.get_function_values_from_local_dof_values(
    true, local_dof_values, qp_interface_values_local);

  fe_iv_view.get_function_values(true, solution, qp_interface_values_global);

  // Gradients
  std::vector<typename View::template solution_gradient_type<NumberType>>
    qp_jump_grads_local(n_q_points), qp_average_grads_local(n_q_points);
  std::vector<gradient_type> qp_jump_grads_global(n_q_points),
    qp_average_grads_global(n_q_points);

  fe_iv_view.get_jump_in_function_gradients_from_local_dof_values(
    local_dof_values, qp_jump_grads_local);
  fe_iv_view.get_jump_in_function_gradients(solution, qp_jump_grads_global);

  fe_iv_view.get_average_of_function_gradients_from_local_dof_values(
    local_dof_values, qp_average_grads_local);
  fe_iv_view.get_average_of_function_gradients(solution,
                                               qp_average_grads_global);

  // Hessians
  std::vector<typename View::template solution_hessian_type<NumberType>>
    qp_jump_hess_local(n_q_points), qp_average_hess_local(n_q_points);
  std::vector<hessian_type> qp_jump_hess_global(n_q_points),
    qp_average_hess_global(n_q_points);

  fe_iv_view.get_jump_in_function_hessians_from_local_dof_values(
    local_dof_values, qp_jump_hess_local);
  fe_iv_view.get_jump_in_function_hessians(solution, qp_jump_hess_global);
  fe_iv_view.get_average_of_function_hessians_from_local_dof_values(
    local_dof_values, qp_average_hess_local);
  fe_iv_view.get_average_of_function_hessians(solution, qp_average_hess_global);

  // Third derivatives
  std::vector<
    typename View::template solution_third_derivative_type<NumberType>>
                                     qp_jump_third_deriv_local(n_q_points);
  std::vector<third_derivative_type> qp_jump_third_deriv_global(n_q_points);

  fe_iv_view.get_jump_in_function_third_derivatives_from_local_dof_values(
    local_dof_values, qp_jump_third_deriv_local);
  fe_iv_view.get_jump_in_function_third_derivatives(solution,
                                                    qp_jump_third_deriv_global);


  // Output
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      if (value_type(qp_interface_values_local[q]) !=
          value_type(qp_interface_values_global[q]))
        deallog << "Interface value NOT OK: Value @ " << q << std::endl;

      if (value_type(qp_jumps_local[q]) != value_type(qp_jumps_global[q]))
        deallog << "Jump NOT OK: Value @ " << q << std::endl;

      if (value_type(qp_averages_local[q]) != value_type(qp_averages_global[q]))
        deallog << "Average NOT OK: Value @ " << q << std::endl;

      if (gradient_type(qp_jump_grads_local[q]) !=
          gradient_type(qp_jump_grads_global[q]))
        deallog << "Jump NOT OK: Grad @ " << q << std::endl;

      if (gradient_type(qp_average_grads_local[q]) !=
          gradient_type(qp_average_grads_global[q]))
        deallog << "Average NOT OK: Grad @ " << q << std::endl;

      if (hessian_type(qp_jump_hess_local[q]) !=
          hessian_type(qp_jump_hess_global[q]))
        deallog << "Jump NOT OK: Hess @ " << q << std::endl;

      if (hessian_type(qp_average_hess_local[q]) !=
          hessian_type(qp_average_hess_global[q]))
        deallog << "Average NOT OK: Hess @ " << q << std::endl;

      if (third_derivative_type(qp_jump_third_deriv_local[q]) !=
          third_derivative_type(qp_jump_third_deriv_global[q]))
        deallog << "Jump NOT OK: 3rd der @ " << q << std::endl;
    }
}

// Vector view
template <typename NumberType, int dim>
void
test_feiv_view(const Vector<double>             &solution,
               const FEInterfaceValues<dim>     &fe_iv,
               const unsigned int               &n_q_points,
               const FEValuesExtractors::Vector &extractor,
               const std::vector<NumberType>    &local_dof_values)
{
  using View =
    std::remove_reference_t<std::remove_const_t<decltype(fe_iv[extractor])>>;
  const View &fe_iv_view = fe_iv[extractor];

  // Typedefs
  using value_type =
    typename ProductType<typename View::value_type, NumberType>::type;
  using gradient_type =
    typename ProductType<typename View::gradient_type, NumberType>::type;
  using hessian_type =
    typename ProductType<typename View::hessian_type, NumberType>::type;
  using third_derivative_type =
    typename ProductType<typename View::third_derivative_type,
                         NumberType>::type;

  // Values
  std::vector<typename View::template solution_value_type<NumberType>>
    qp_jumps_local(n_q_points), qp_averages_local(n_q_points),
    qp_interface_values_local(n_q_points);
  std::vector<value_type> qp_jumps_global(n_q_points),
    qp_averages_global(n_q_points), qp_interface_values_global(n_q_points);

  fe_iv_view.get_jump_in_function_values_from_local_dof_values(local_dof_values,
                                                               qp_jumps_local);
  fe_iv_view.get_jump_in_function_values(solution, qp_jumps_global);

  fe_iv_view.get_average_of_function_values_from_local_dof_values(
    local_dof_values, qp_averages_local);
  fe_iv_view.get_average_of_function_values(solution, qp_averages_global);

  fe_iv_view.get_function_values_from_local_dof_values(
    true, local_dof_values, qp_interface_values_local);

  fe_iv_view.get_function_values(true, solution, qp_interface_values_global);

  // Gradients
  std::vector<typename View::template solution_gradient_type<NumberType>>
    qp_jump_grads_local(n_q_points), qp_average_grads_local(n_q_points);
  std::vector<gradient_type> qp_jump_grads_global(n_q_points),
    qp_average_grads_global(n_q_points);

  fe_iv_view.get_jump_in_function_gradients_from_local_dof_values(
    local_dof_values, qp_jump_grads_local);
  fe_iv_view.get_jump_in_function_gradients(solution, qp_jump_grads_global);

  fe_iv_view.get_average_of_function_gradients_from_local_dof_values(
    local_dof_values, qp_average_grads_local);
  fe_iv_view.get_average_of_function_gradients(solution,
                                               qp_average_grads_global);

  // Hessians
  std::vector<typename View::template solution_hessian_type<NumberType>>
    qp_jump_hess_local(n_q_points), qp_average_hess_local(n_q_points);
  std::vector<hessian_type> qp_jump_hess_global(n_q_points),
    qp_average_hess_global(n_q_points);

  fe_iv_view.get_jump_in_function_hessians_from_local_dof_values(
    local_dof_values, qp_jump_hess_local);
  fe_iv_view.get_jump_in_function_hessians(solution, qp_jump_hess_global);
  fe_iv_view.get_average_of_function_hessians_from_local_dof_values(
    local_dof_values, qp_average_hess_local);
  fe_iv_view.get_average_of_function_hessians(solution, qp_average_hess_global);

  // Third derivatives
  std::vector<
    typename View::template solution_third_derivative_type<NumberType>>
                                     qp_jump_third_deriv_local(n_q_points);
  std::vector<third_derivative_type> qp_jump_third_deriv_global(n_q_points);

  fe_iv_view.get_jump_in_function_third_derivatives_from_local_dof_values(
    local_dof_values, qp_jump_third_deriv_local);
  fe_iv_view.get_jump_in_function_third_derivatives(solution,
                                                    qp_jump_third_deriv_global);


  // Output
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      if (value_type(qp_interface_values_local[q]) !=
          value_type(qp_interface_values_global[q]))
        deallog << "Interface value NOT OK: Value @ " << q << std::endl;

      if (value_type(qp_jumps_local[q]) != value_type(qp_jumps_global[q]))
        deallog << "Jump NOT OK: Value @ " << q << std::endl;

      if (value_type(qp_averages_local[q]) != value_type(qp_averages_global[q]))
        deallog << "Average NOT OK: Value @ " << q << std::endl;

      if (value_type(qp_averages_local[q]) != value_type(qp_averages_global[q]))
        deallog << "Average NOT OK: Value @ " << q << std::endl;

      if (gradient_type(qp_jump_grads_local[q]) !=
          gradient_type(qp_jump_grads_global[q]))
        deallog << "Jump NOT OK: Grad @ " << q << std::endl;

      if (gradient_type(qp_average_grads_local[q]) !=
          gradient_type(qp_average_grads_global[q]))
        deallog << "Average NOT OK: Grad @ " << q << std::endl;

      if (hessian_type(qp_jump_hess_local[q]) !=
          hessian_type(qp_jump_hess_global[q]))
        deallog << "Jump NOT OK: Hess @ " << q << std::endl;

      if (hessian_type(qp_average_hess_local[q]) !=
          hessian_type(qp_average_hess_global[q]))
        deallog << "Average NOT OK: Hess @ " << q << std::endl;

      if (third_derivative_type(qp_jump_third_deriv_local[q]) !=
          third_derivative_type(qp_jump_third_deriv_global[q]))
        deallog << "Jump NOT OK: 3rd der @ " << q << std::endl;
    }
}

template <typename NumberType, int dim, typename FEType, typename ExtractorType>
void
test_feiv_extractor(const FEType &fe, const ExtractorType &extractor)
{
  Triangulation<dim> tria;
  DoFHandler<dim>    dof_handler(tria);
  Vector<double>     solution;

  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global(1);
  dof_handler.distribute_dofs(fe);
  DoFRenumbering::random(dof_handler);
  solution.reinit(dof_handler.n_dofs());

  // Populate with non-trivial values
  for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
    {
      solution(i) = i + 1;
    }

  const QGauss<dim - 1>  face_quadrature(2);
  FEInterfaceValues<dim> fe_iv(fe,
                               face_quadrature,
                               update_values | update_gradients |
                                 update_quadrature_points |
                                 update_3rd_derivatives);

  typename DoFHandler<dim>::active_cell_iterator cell =
    dof_handler.begin_active();
  for (const auto f : cell->face_indices())
    if (!cell->face(f)->at_boundary())
      {
        {
          fe_iv.reinit(cell,
                       f,
                       numbers::invalid_unsigned_int,
                       cell->neighbor(f),
                       cell->neighbor_of_neighbor(f),
                       numbers::invalid_unsigned_int);
          const std::vector<types::global_dof_index> &interface_dof_indices =
            fe_iv.get_interface_dof_indices();

          // Convert the DoF values so that they are potentially of
          // a different number type
          std::vector<NumberType> interface_dof_values_other(
            fe_iv.n_current_interface_dofs());
          for (unsigned int i = 0; i < fe_iv.n_current_interface_dofs(); ++i)
            interface_dof_values_other[i] = solution(interface_dof_indices[i]);

          test_feiv_view(solution,
                         fe_iv,
                         face_quadrature.size(),
                         extractor,
                         interface_dof_values_other);
        }
      }

  deallog << "OK" << std::endl;
}


template <int dim = 2>
void
test_jump_function()
{
  FE_DGP<dim>                      fe(0);
  const FEValuesExtractors::Scalar extractor(0);

  Triangulation<dim> tria;
  DoFHandler<dim>    dof_handler(tria);
  Vector<double>     solution;

  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global(1);
  dof_handler.distribute_dofs(fe);

  solution.reinit(dof_handler.n_dofs());

  std::vector<types::global_dof_index> local_dof_indices(1);
  // Populate with non-trivial values
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell->get_dof_indices(local_dof_indices);
      solution(local_dof_indices[0]) =
        static_cast<double>(cell->active_cell_index());
    }

  const QGauss<dim - 1>  face_quadrature(2);
  FEInterfaceValues<dim> fe_iv(fe,
                               face_quadrature,
                               update_values | update_gradients |
                                 update_quadrature_points |
                                 update_3rd_derivatives);

  std::vector<double> qp_jumps_global(face_quadrature.size());

  typename DoFHandler<dim>::active_cell_iterator cell =
    dof_handler.begin_active();
  for (; cell != dof_handler.end(); ++cell)
    for (const auto f : cell->face_indices())
      if (!cell->face(f)->at_boundary())
        {
          {
            fe_iv.reinit(cell,
                         f,
                         numbers::invalid_unsigned_int,
                         cell->neighbor(f),
                         cell->neighbor_of_neighbor(f),
                         numbers::invalid_unsigned_int);

            fe_iv[extractor].get_jump_in_function_values(solution,
                                                         qp_jumps_global);

            double exact = cell->active_cell_index() * 1.0 -
                           cell->neighbor(f)->active_cell_index() * 1.0;

            for (unsigned int q = 0; q < face_quadrature.size(); ++q)
              Assert(std::fabs(qp_jumps_global[q] - exact) < 1e-15,
                     ExcNotImplemented());
          }
        }

  deallog << "OK" << std::endl;
}

template <typename NumberType, int dim = 2>
void
test()
{
  const unsigned int degree = 3; // Need third derivatives

  deallog.push("FEInterfaceViews Scalar");
  {
    FE_Q<dim>                        fe(degree);
    const FEValuesExtractors::Scalar extractor(0);
    test_feiv_extractor<NumberType, dim>(fe, extractor);
  }
  deallog.pop();

  deallog.push("FEInterfaceViews Vector");
  {
    FESystem<dim>                    fe(FE_Q<dim>(degree), dim);
    const FEValuesExtractors::Vector extractor(0);
    test_feiv_extractor<NumberType, dim>(fe, extractor);
  }
  deallog.pop();
}

int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  deallog.push("Float");
  {
    test<float>();
  }
  deallog.pop();

  deallog.push("Double");
  {
    test<double>();
  }
  deallog.pop();

  deallog.push("Test jump functions");
  {
    test_jump_function();
  }
  deallog.pop();

  deallog << "OK" << std::endl;
}

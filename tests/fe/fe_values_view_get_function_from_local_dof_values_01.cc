// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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



// Check that FEValuesViews::get_function_*_from_local_dof_values
// works with different number types.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

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
test_view(const Vector<double> &         solution,
          const FEValues<dim> &          fe_values,
          const unsigned int &           n_q_points,
          const ExtractorType &          extractor,
          const std::vector<NumberType> &local_dof_values);

// Scalar view
template <typename NumberType, int dim>
void
test_view(const Vector<double> &            solution,
          const FEValues<dim> &             fe_values,
          const unsigned int &              n_q_points,
          const FEValuesExtractors::Scalar &extractor,
          const std::vector<NumberType> &   local_dof_values)
{
  using View = typename std::remove_reference<
    typename std::remove_const<decltype(fe_values[extractor])>::type>::type;
  const View &fe_values_view = fe_values[extractor];

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
                          qp_values_local(n_q_points);
  std::vector<value_type> qp_values_global(n_q_points);
  fe_values_view.get_function_values_from_local_dof_values(local_dof_values,
                                                           qp_values_local);
  fe_values_view.get_function_values(solution, qp_values_global);

  // Gradients
  std::vector<typename View::template solution_gradient_type<NumberType>>
                             qp_grads_local(n_q_points);
  std::vector<gradient_type> qp_grads_global(n_q_points);
  fe_values_view.get_function_gradients_from_local_dof_values(local_dof_values,
                                                              qp_grads_local);
  fe_values_view.get_function_gradients(solution, qp_grads_global);

  // Hessians
  std::vector<typename View::template solution_hessian_type<NumberType>>
                            qp_hess_local(n_q_points);
  std::vector<hessian_type> qp_hess_global(n_q_points);
  fe_values_view.get_function_hessians_from_local_dof_values(local_dof_values,
                                                             qp_hess_local);
  fe_values_view.get_function_hessians(solution, qp_hess_global);

  // Laplacians
  std::vector<typename View::template solution_laplacian_type<NumberType>>
                              qp_laplace_local(n_q_points);
  std::vector<laplacian_type> qp_laplace_global(n_q_points);
  fe_values_view.get_function_laplacians_from_local_dof_values(
    local_dof_values, qp_laplace_local);
  fe_values_view.get_function_laplacians(solution, qp_laplace_global);

  // Third derivatives
  std::vector<
    typename View::template solution_third_derivative_type<NumberType>>
                                     qp_third_deriv_local(n_q_points);
  std::vector<third_derivative_type> qp_third_deriv_global(n_q_points);
  fe_values_view.get_function_third_derivatives_from_local_dof_values(
    local_dof_values, qp_third_deriv_local);
  fe_values_view.get_function_third_derivatives(solution,
                                                qp_third_deriv_global);

  // Output
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      if (value_type(qp_values_local[q]) != value_type(qp_values_global[q]))
        deallog << "NOT OK: Value @ " << q << std::endl;

      if (gradient_type(qp_grads_local[q]) != gradient_type(qp_grads_global[q]))
        deallog << "NOT OK: Grad @ " << q << std::endl;

      if (hessian_type(qp_hess_local[q]) != hessian_type(qp_hess_global[q]))
        deallog << "NOT OK: Hess @ " << q << std::endl;

      if (laplacian_type(qp_laplace_local[q]) !=
          laplacian_type(qp_laplace_global[q]))
        deallog << "NOT OK: Laplace @ " << q << std::endl;

      if (third_derivative_type(qp_third_deriv_local[q]) !=
          third_derivative_type(qp_third_deriv_global[q]))
        deallog << "NOT OK: 3rd der @ " << q << std::endl;
    }
}

// Vector view
template <typename NumberType, int dim>
void
test_view(const Vector<double> &            solution,
          const FEValues<dim> &             fe_values,
          const unsigned int &              n_q_points,
          const FEValuesExtractors::Vector &extractor,
          const std::vector<NumberType> &   local_dof_values)
{
  using View = typename std::remove_reference<
    typename std::remove_const<decltype(fe_values[extractor])>::type>::type;
  const View &fe_values_view = fe_values[extractor];

  // Typedefs
  using value_type =
    typename ProductType<typename View::value_type, NumberType>::type;
  using gradient_type =
    typename ProductType<typename View::gradient_type, NumberType>::type;
  using symmetric_gradient_type =
    typename ProductType<typename View::symmetric_gradient_type,
                         NumberType>::type;
  using divergence_type =
    typename ProductType<typename View::divergence_type, NumberType>::type;
  using curl_type =
    typename ProductType<typename View::curl_type, NumberType>::type;
  using hessian_type =
    typename ProductType<typename View::hessian_type, NumberType>::type;
  using laplacian_type =
    typename ProductType<typename View::value_type, NumberType>::type;
  using third_derivative_type =
    typename ProductType<typename View::third_derivative_type,
                         NumberType>::type;

  // Values
  std::vector<typename View::template solution_value_type<NumberType>>
                          qp_values_local(n_q_points);
  std::vector<value_type> qp_values_global(n_q_points);
  fe_values_view.get_function_values_from_local_dof_values(local_dof_values,
                                                           qp_values_local);
  fe_values_view.get_function_values(solution, qp_values_global);

  // Gradients
  std::vector<typename View::template solution_gradient_type<NumberType>>
                             qp_grads_local(n_q_points);
  std::vector<gradient_type> qp_grads_global(n_q_points);
  fe_values_view.get_function_gradients_from_local_dof_values(local_dof_values,
                                                              qp_grads_local);
  fe_values_view.get_function_gradients(solution, qp_grads_global);

  // Symmetric gradients
  std::vector<
    typename View::template solution_symmetric_gradient_type<NumberType>>
                                       qp_symm_grads_local(n_q_points);
  std::vector<symmetric_gradient_type> qp_symm_grads_global(n_q_points);
  fe_values_view.get_function_symmetric_gradients_from_local_dof_values(
    local_dof_values, qp_symm_grads_local);
  fe_values_view.get_function_symmetric_gradients(solution,
                                                  qp_symm_grads_global);

  // Divergences
  std::vector<typename View::template solution_divergence_type<NumberType>>
                               qp_divs_local(n_q_points);
  std::vector<divergence_type> qp_divs_global(n_q_points);
  fe_values_view.get_function_divergences_from_local_dof_values(
    local_dof_values, qp_divs_local);
  fe_values_view.get_function_divergences(solution, qp_divs_global);

  // Curls
  std::vector<typename View::template solution_curl_type<NumberType>>
                         qp_curls_local(n_q_points);
  std::vector<curl_type> qp_curls_global(n_q_points);
  fe_values_view.get_function_curls_from_local_dof_values(local_dof_values,
                                                          qp_curls_local);
  fe_values_view.get_function_curls(solution, qp_curls_global);

  // Hessians
  std::vector<typename View::template solution_hessian_type<NumberType>>
                            qp_hess_local(n_q_points);
  std::vector<hessian_type> qp_hess_global(n_q_points);
  fe_values_view.get_function_hessians_from_local_dof_values(local_dof_values,
                                                             qp_hess_local);
  fe_values_view.get_function_hessians(solution, qp_hess_global);

  // Laplacians
  std::vector<typename View::template solution_laplacian_type<NumberType>>
                              qp_laplace_local(n_q_points);
  std::vector<laplacian_type> qp_laplace_global(n_q_points);
  fe_values_view.get_function_laplacians_from_local_dof_values(
    local_dof_values, qp_laplace_local);
  fe_values_view.get_function_laplacians(solution, qp_laplace_global);

  // Third derivatives
  std::vector<
    typename View::template solution_third_derivative_type<NumberType>>
                                     qp_third_deriv_local(n_q_points);
  std::vector<third_derivative_type> qp_third_deriv_global(n_q_points);
  fe_values_view.get_function_third_derivatives_from_local_dof_values(
    local_dof_values, qp_third_deriv_local);
  fe_values_view.get_function_third_derivatives(solution,
                                                qp_third_deriv_global);

  // Output
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      if (value_type(qp_values_local[q]) != value_type(qp_values_global[q]))
        deallog << "NOT OK: Value @ " << q << std::endl;

      if (gradient_type(qp_grads_local[q]) != gradient_type(qp_grads_global[q]))
        deallog << "NOT OK: Grad @ " << q << std::endl;

      if (gradient_type(qp_symm_grads_local[q]) !=
          gradient_type(qp_symm_grads_global[q]))
        deallog << "NOT OK: Symm grad @ " << q << std::endl;

      if (divergence_type(qp_divs_local[q]) !=
          divergence_type(qp_divs_global[q]))
        deallog << "NOT OK: Div @ " << q << std::endl;

      // Note: FE_Q's are curl free: Should always be zero'd
      // So we are just checking that we don't hit an internal assert
      // when doing the above calls, rather than testing the values
      if (curl_type(qp_curls_local[q]) != curl_type(qp_curls_global[q]))
        deallog << "NOT OK: Curl @ " << q << std::endl;

      if (hessian_type(qp_hess_local[q]) != hessian_type(qp_hess_global[q]))
        deallog << "NOT OK: Hess @ " << q << std::endl;

      if (laplacian_type(qp_laplace_local[q]) !=
          laplacian_type(qp_laplace_global[q]))
        deallog << "NOT OK: Laplace @ " << q << std::endl;

      if (third_derivative_type(qp_third_deriv_local[q]) !=
          third_derivative_type(qp_third_deriv_global[q]))
        deallog << "NOT OK: 3rd der @ " << q << std::endl;
    }
}

// SymmetricTensor view
template <typename NumberType, int dim>
void
test_view(const Vector<double> &                        solution,
          const FEValues<dim> &                         fe_values,
          const unsigned int &                          n_q_points,
          const FEValuesExtractors::SymmetricTensor<2> &extractor,
          const std::vector<NumberType> &               local_dof_values)
{
  using View = typename std::remove_reference<
    typename std::remove_const<decltype(fe_values[extractor])>::type>::type;
  const View &fe_values_view = fe_values[extractor];

  // Typedefs
  using value_type =
    typename ProductType<typename View::value_type, NumberType>::type;
  using divergence_type =
    typename ProductType<typename View::divergence_type, NumberType>::type;

  // Values
  std::vector<typename View::template solution_value_type<NumberType>>
                          qp_values_local(n_q_points);
  std::vector<value_type> qp_values_global(n_q_points);
  fe_values_view.get_function_values_from_local_dof_values(local_dof_values,
                                                           qp_values_local);
  fe_values_view.get_function_values(solution, qp_values_global);

  // Divergences
  std::vector<typename View::template solution_divergence_type<NumberType>>
                               qp_divs_local(n_q_points);
  std::vector<divergence_type> qp_divs_global(n_q_points);
  fe_values_view.get_function_divergences_from_local_dof_values(
    local_dof_values, qp_divs_local);
  fe_values_view.get_function_divergences(solution, qp_divs_global);

  // Output
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      if (value_type(qp_values_local[q]) != value_type(qp_values_global[q]))
        deallog << "NOT OK: Value @ " << q << std::endl;

      if (divergence_type(qp_divs_local[q]) !=
          divergence_type(qp_divs_global[q]))
        deallog << "NOT OK: Div @ " << q << std::endl;
    }
}

// Tensor view
template <typename NumberType, int dim>
void
test_view(const Vector<double> &               solution,
          const FEValues<dim> &                fe_values,
          const unsigned int &                 n_q_points,
          const FEValuesExtractors::Tensor<2> &extractor,
          const std::vector<NumberType> &      local_dof_values)
{
  using View = typename std::remove_reference<
    typename std::remove_const<decltype(fe_values[extractor])>::type>::type;
  const View &fe_values_view = fe_values[extractor];

  // Typedefs
  using value_type =
    typename ProductType<typename View::value_type, NumberType>::type;
  using gradient_type =
    typename ProductType<typename View::gradient_type, NumberType>::type;
  using divergence_type =
    typename ProductType<typename View::divergence_type, NumberType>::type;

  // Values
  std::vector<typename View::template solution_value_type<NumberType>>
                          qp_values_local(n_q_points);
  std::vector<value_type> qp_values_global(n_q_points);
  fe_values_view.get_function_values_from_local_dof_values(local_dof_values,
                                                           qp_values_local);
  fe_values_view.get_function_values(solution, qp_values_global);

  // Divergences
  std::vector<typename View::template solution_divergence_type<NumberType>>
                               qp_divs_local(n_q_points);
  std::vector<divergence_type> qp_divs_global(n_q_points);
  fe_values_view.get_function_divergences_from_local_dof_values(
    local_dof_values, qp_divs_local);
  fe_values_view.get_function_divergences(solution, qp_divs_global);

  // Gradients
  std::vector<typename View::template solution_gradient_type<NumberType>>
                             qp_grads_local(n_q_points);
  std::vector<gradient_type> qp_grads_global(n_q_points);
  fe_values_view.get_function_gradients_from_local_dof_values(local_dof_values,
                                                              qp_grads_local);
  fe_values_view.get_function_gradients(solution, qp_grads_global);

  // Output
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      if (value_type(qp_values_local[q]) != value_type(qp_values_global[q]))
        deallog << "NOT OK: Value @ " << q << std::endl;

      if (divergence_type(qp_divs_local[q]) !=
          divergence_type(qp_divs_global[q]))
        deallog << "NOT OK: Div @ " << q << std::endl;

      if (gradient_type(qp_grads_local[q]) != gradient_type(qp_grads_global[q]))
        deallog << "NOT OK: Grad @ " << q << std::endl;
    }
}

template <typename NumberType, int dim, typename FEType, typename ExtractorType>
void
test_extractor(const FEType &fe, const ExtractorType &extractor)
{
  QGauss<dim> quadrature_formula(2);

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

  FEValues<dim>                        fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients | update_hessians |
                            update_3rd_derivatives);
  std::vector<types::global_dof_index> local_dof_indices(fe.dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell =
    dof_handler.begin_active();
  {
    fe_values.reinit(cell);
    cell->get_dof_indices(local_dof_indices);

    std::vector<double> local_dof_values(fe.dofs_per_cell);
    cell->get_dof_values(solution,
                         local_dof_values.begin(),
                         local_dof_values.end());

    // Convert the DoF values so that they are potentially of
    // a different number type
    std::vector<NumberType> local_dof_values_other(fe.dofs_per_cell);
    for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
      local_dof_values_other[i] = local_dof_values[i];

    test_view(solution,
              fe_values,
              quadrature_formula.size(),
              extractor,
              local_dof_values_other);
  }

  deallog << "OK" << std::endl;
}

template <typename NumberType, int dim = 2>
void
test()
{
  const unsigned int degree = 3; // Need third derivatives

  deallog.push("Scalar");
  {
    FE_Q<dim>                  fe(degree);
    FEValuesExtractors::Scalar extractor(0);
    test_extractor<NumberType, dim>(fe, extractor);
  }
  deallog.pop();

  deallog.push("Vector");
  {
    FESystem<dim>              fe(FE_Q<dim>(degree), dim);
    FEValuesExtractors::Vector extractor(0);
    test_extractor<NumberType, dim>(fe, extractor);
  }
  deallog.pop();

  deallog.push("SymmetricTensor");
  {
    FESystem<dim>                          fe(FE_Q<dim>(degree),
                     SymmetricTensor<2, dim>::n_independent_components);
    FEValuesExtractors::SymmetricTensor<2> extractor(0);
    test_extractor<NumberType, dim>(fe, extractor);
  }
  deallog.pop();

  deallog.push("Tensor");
  {
    FESystem<dim>                 fe(FE_Q<dim>(degree),
                     Tensor<2, dim>::n_independent_components);
    FEValuesExtractors::Tensor<2> extractor(0);
    test_extractor<NumberType, dim>(fe, extractor);
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

  deallog << "OK" << std::endl;
}

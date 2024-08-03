// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check that FEValuesViews::get_function_*_from_local_dof_values
// works with different number types.
// As opposed to fe_values_view_get_function_from_local_dof_values_01.cc, this
// test does not check the output types, but rather that all necessary
// instantiations are in place and that no assertions are triggered.
// This is because FEValuesViews::get_function_* cannot work with Sacado types.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/differentiation/ad/sacado_product_types.h>

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
test_view(const Vector<double>          &solution,
          const FEValues<dim>           &fe_values,
          const unsigned int            &n_q_points,
          const ExtractorType           &extractor,
          const std::vector<NumberType> &local_dof_values);

// Scalar view
template <typename NumberType, int dim>
void
test_view(const Vector<double>             &solution,
          const FEValues<dim>              &fe_values,
          const unsigned int               &n_q_points,
          const FEValuesExtractors::Scalar &extractor,
          const std::vector<NumberType>    &local_dof_values)
{
  using View = std::remove_reference_t<
    std::remove_const_t<decltype(fe_values[extractor])>>;
  const View &fe_values_view = fe_values[extractor];

  // Typedefs
  // Values
  std::vector<typename View::template solution_value_type<NumberType>>
    qp_values_local(n_q_points);
  fe_values_view.get_function_values_from_local_dof_values(local_dof_values,
                                                           qp_values_local);

  // Gradients
  std::vector<typename View::template solution_gradient_type<NumberType>>
    qp_grads_local(n_q_points);
  fe_values_view.get_function_gradients_from_local_dof_values(local_dof_values,
                                                              qp_grads_local);

  // Hessians
  std::vector<typename View::template solution_hessian_type<NumberType>>
    qp_hess_local(n_q_points);
  fe_values_view.get_function_hessians_from_local_dof_values(local_dof_values,
                                                             qp_hess_local);

  // Laplacians
  std::vector<typename View::template solution_laplacian_type<NumberType>>
    qp_laplace_local(n_q_points);
  fe_values_view.get_function_laplacians_from_local_dof_values(
    local_dof_values, qp_laplace_local);

  // Third derivatives
  std::vector<
    typename View::template solution_third_derivative_type<NumberType>>
    qp_third_deriv_local(n_q_points);
  fe_values_view.get_function_third_derivatives_from_local_dof_values(
    local_dof_values, qp_third_deriv_local);
}

// Vector view
template <typename NumberType, int dim>
void
test_view(const Vector<double>             &solution,
          const FEValues<dim>              &fe_values,
          const unsigned int               &n_q_points,
          const FEValuesExtractors::Vector &extractor,
          const std::vector<NumberType>    &local_dof_values)
{
  using View = std::remove_reference_t<
    std::remove_const_t<decltype(fe_values[extractor])>>;
  const View &fe_values_view = fe_values[extractor];

  // Typedefs
  // Values
  std::vector<typename View::template solution_value_type<NumberType>>
    qp_values_local(n_q_points);
  fe_values_view.get_function_values_from_local_dof_values(local_dof_values,
                                                           qp_values_local);

  // Gradients
  std::vector<typename View::template solution_gradient_type<NumberType>>
    qp_grads_local(n_q_points);
  fe_values_view.get_function_gradients_from_local_dof_values(local_dof_values,
                                                              qp_grads_local);

  // Symmetric gradients
  std::vector<
    typename View::template solution_symmetric_gradient_type<NumberType>>
    qp_symm_grads_local(n_q_points);
  fe_values_view.get_function_symmetric_gradients_from_local_dof_values(
    local_dof_values, qp_symm_grads_local);

  // Divergences
  std::vector<typename View::template solution_divergence_type<NumberType>>
    qp_divs_local(n_q_points);
  fe_values_view.get_function_divergences_from_local_dof_values(
    local_dof_values, qp_divs_local);

  // Curls
  std::vector<typename View::template solution_curl_type<NumberType>>
    qp_curls_local(n_q_points);
  fe_values_view.get_function_curls_from_local_dof_values(local_dof_values,
                                                          qp_curls_local);

  // Hessians
  std::vector<typename View::template solution_hessian_type<NumberType>>
    qp_hess_local(n_q_points);
  fe_values_view.get_function_hessians_from_local_dof_values(local_dof_values,
                                                             qp_hess_local);

  // Laplacians
  std::vector<typename View::template solution_laplacian_type<NumberType>>
    qp_laplace_local(n_q_points);
  fe_values_view.get_function_laplacians_from_local_dof_values(
    local_dof_values, qp_laplace_local);

  // Third derivatives
  std::vector<
    typename View::template solution_third_derivative_type<NumberType>>
    qp_third_deriv_local(n_q_points);
  fe_values_view.get_function_third_derivatives_from_local_dof_values(
    local_dof_values, qp_third_deriv_local);
}

// SymmetricTensor view
template <typename NumberType, int dim>
void
test_view(const Vector<double>                         &solution,
          const FEValues<dim>                          &fe_values,
          const unsigned int                           &n_q_points,
          const FEValuesExtractors::SymmetricTensor<2> &extractor,
          const std::vector<NumberType>                &local_dof_values)
{
  using View = std::remove_reference_t<
    std::remove_const_t<decltype(fe_values[extractor])>>;
  const View &fe_values_view = fe_values[extractor];

  // Typedefs
  using value_type =
    typename ProductType<typename View::value_type, NumberType>::type;
  using divergence_type =
    typename ProductType<typename View::divergence_type, NumberType>::type;

  // Values
  std::vector<typename View::template solution_value_type<NumberType>>
    qp_values_local(n_q_points);
  fe_values_view.get_function_values_from_local_dof_values(local_dof_values,
                                                           qp_values_local);

  // Divergences
  std::vector<typename View::template solution_divergence_type<NumberType>>
    qp_divs_local(n_q_points);
  fe_values_view.get_function_divergences_from_local_dof_values(
    local_dof_values, qp_divs_local);
}

// Tensor view
template <typename NumberType, int dim>
void
test_view(const Vector<double>                &solution,
          const FEValues<dim>                 &fe_values,
          const unsigned int                  &n_q_points,
          const FEValuesExtractors::Tensor<2> &extractor,
          const std::vector<NumberType>       &local_dof_values)
{
  using View = std::remove_reference_t<
    std::remove_const_t<decltype(fe_values[extractor])>>;
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
  fe_values_view.get_function_values_from_local_dof_values(local_dof_values,
                                                           qp_values_local);

  // Divergences
  std::vector<typename View::template solution_divergence_type<NumberType>>
    qp_divs_local(n_q_points);
  fe_values_view.get_function_divergences_from_local_dof_values(
    local_dof_values, qp_divs_local);

  // Gradients
  std::vector<typename View::template solution_gradient_type<NumberType>>
    qp_grads_local(n_q_points);
  fe_values_view.get_function_gradients_from_local_dof_values(local_dof_values,
                                                              qp_grads_local);
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
    FE_Q<dim>                        fe(degree);
    const FEValuesExtractors::Scalar extractor(0);
    test_extractor<NumberType, dim>(fe, extractor);
  }
  deallog.pop();

  deallog.push("Vector");
  {
    FESystem<dim>                    fe(FE_Q<dim>(degree), dim);
    const FEValuesExtractors::Vector extractor(0);
    test_extractor<NumberType, dim>(fe, extractor);
  }
  deallog.pop();

  deallog.push("SymmetricTensor");
  {
    FESystem<dim>                                fe(FE_Q<dim>(degree),
                     SymmetricTensor<2, dim>::n_independent_components);
    const FEValuesExtractors::SymmetricTensor<2> extractor(0);
    test_extractor<NumberType, dim>(fe, extractor);
  }
  deallog.pop();

  deallog.push("Tensor");
  {
    FESystem<dim>                       fe(FE_Q<dim>(degree),
                     Tensor<2, dim>::n_independent_components);
    const FEValuesExtractors::Tensor<2> extractor(0);
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

  deallog.push("Sacado::Fad::DFad<float>");
  {
    test<Sacado::Fad::DFad<float>>();
  }
  deallog.pop();

  deallog.push("Sacado::Fad::DFad<Sacado::Fad::DFad<float>>");
  {
    test<Sacado::Fad::DFad<Sacado::Fad::DFad<float>>>();
  }
  deallog.pop();

  deallog.push("Sacado::Fad::DFad<double>");
  {
    test<Sacado::Fad::DFad<double>>();
  }
  deallog.pop();

  deallog.push("Sacado::Fad::DFad<Sacado::Fad::DFad<double>>");
  {
    test<Sacado::Fad::DFad<Sacado::Fad::DFad<double>>>();
  }
  deallog.pop();

  deallog << "OK" << std::endl;
}

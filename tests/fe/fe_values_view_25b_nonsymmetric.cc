// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2018 by the deal.II authors
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



// similar to _25_nonsymmetric but also tests gradients as well as
// get_function_XYZ.
// Linear primitive elements, so for scalar 2D element those are:
//
//   *-----*
//   |     |
//   |     |
//   x-----*
//  N0 = (1-x1)(1-x2)
//  N1 =    x1 (1-x2)
//  N3 = (1-x1)   x2
//  N4 =    x1    x2
//  d/dx1 N0 = -(1-x2)
//  d/dx2 N0 = -(1-x1)
// and same for other support points
// To facilitate comparison, we use auxiliary scalar-valued FE
// and take its values and gradients. For primitive FEs we can
// deduce which coefficient of tensor-valued shape function is
// nonzero and compare values/gradient to the scalar-valued counterpart.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"



template <int dim>
void
test(const Triangulation<dim> &tr, const FiniteElement<dim> &fe)
{
  deallog << "FE=" << fe.get_name() << std::endl;

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  const QGauss<dim> quadrature(2);

  FE_Q<dim>                        fe_scalar(1);
  std::vector<std::vector<double>> scalar_values(
    fe_scalar.dofs_per_cell, std::vector<double>(quadrature.size()));
  std::vector<std::vector<Tensor<1, dim>>> scalar_gradients(
    fe_scalar.dofs_per_cell, std::vector<Tensor<1, dim>>(quadrature.size()));
  // fill-in scalar values and gradients
  {
    DoFHandler<dim> dof_scalar(tr);
    dof_scalar.distribute_dofs(fe_scalar);
    FEValues<dim> fe_values_scalar(fe_scalar,
                                   quadrature,
                                   update_values | update_gradients);
    fe_values_scalar.reinit(dof_scalar.begin_active());

    for (unsigned int i = 0; i < fe_scalar.dofs_per_cell; ++i)
      {
        const unsigned int i_node = fe_scalar.system_to_base_index(i).second;
        Assert(i_node < fe_scalar.dofs_per_cell, ExcInternalError());
        for (unsigned int q = 0; q < quadrature.size(); ++q)
          {
            scalar_values[i_node][q]    = fe_values_scalar.shape_value(i, q);
            scalar_gradients[i_node][q] = fe_values_scalar.shape_grad(i, q);
          }
      }
  }

  FEValues<dim> fe_values(fe,
                          quadrature,
                          update_values | update_gradients |
                            update_quadrature_points);
  fe_values.reinit(dof.begin_active());

  // let the FEValues object compute the
  // divergences at quadrature points
  FEValuesExtractors::Tensor<2> extractor(0);

  // also compare get_function_values/gradients/divergences
  // the manual evaluation
  Vector<double> fe_function(dof.n_dofs());
  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    fe_function(i) = (i + 1) * (3 + i);

  std::vector<Tensor<1, dim>> divergences(quadrature.size()),
    divergences_manual(quadrature.size());
  fe_values[extractor].get_function_divergences(fe_function, divergences);

  std::vector<Tensor<3, dim>> gradients(quadrature.size()),
    gradients_manual(quadrature.size());
  fe_values[extractor].get_function_gradients(fe_function, gradients);

  std::vector<Tensor<2, dim>> values(quadrature.size()),
    values_manual(quadrature.size());
  fe_values[extractor].get_function_values(fe_function, values);

  std::vector<types::global_dof_index> local_dof_indices(fe.dofs_per_cell);
  dof.begin_active()->get_dof_indices(local_dof_indices);

  for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
    {
      const double f_val = fe_function(local_dof_indices[i]);
      const auto   val   = fe_values[extractor].value(i, 0);
      // find out which component is non-zero
      TableIndices<2> nonzero_ind;
      for (unsigned int k = 0; k < Tensor<2, dim>::n_independent_components;
           ++k)
        {
          nonzero_ind =
            dealii::Tensor<2, dim>::unrolled_to_component_indices(k);
          if (std::abs(val[nonzero_ind]) > 1e-10)
            break;
        }

      // the support point (node) id
      const unsigned int i_node = fe.system_to_base_index(i).second;

      deallog << "i=" << i << " ii=" << nonzero_ind[0]
              << " jj=" << nonzero_ind[1] << " node=" << i_node << std::endl;

      for (unsigned int q = 0; q < quadrature.size(); ++q)
        {
          const auto val_q  = fe_values[extractor].value(i, q);
          const auto grad_q = fe_values[extractor].gradient(i, q);
          const auto div_q  = fe_values[extractor].divergence(i, q);
          deallog << "  q_point=" << fe_values.quadrature_point(q) << std::endl
                  << "    value= " << val_q << std::endl
                  << "    grad= " << grad_q << std::endl
                  << "    div= " << div_q << std::endl;

          values_manual[q] += val_q * f_val;
          gradients_manual[q] += grad_q * f_val;
          divergences_manual[q] += div_q * f_val;

          // check value and gradients:
          for (unsigned int k = 0; k < Tensor<2, dim>::n_independent_components;
               ++k)
            {
              const TableIndices<2> ind_k =
                dealii::Tensor<2, dim>::unrolled_to_component_indices(k);
              if (ind_k == nonzero_ind)
                {
                  AssertThrow(val_q[ind_k] == scalar_values[i_node][q],
                              ExcInternalError());
                  AssertThrow((grad_q[ind_k[0]][ind_k[1]] ==
                               scalar_gradients[i_node][q]),
                              ExcInternalError());
                }
              else
                {
                  AssertThrow(val_q[ind_k] == 0.,
                              ExcMessage(std::to_string(k) + " " +
                                         std::to_string(ind_k[0]) + " " +
                                         std::to_string(ind_k[1])));
                  AssertThrow((grad_q[ind_k[0]][ind_k[1]] == Tensor<1, dim>()),
                              ExcInternalError());
                }
            }
          // finally check consistency between gradient and divergence, namely
          // G : I = Div or
          // \sum_j G_{ijj} = Div_i
          for (unsigned int i = 0; i < dim; ++i)
            {
              double div_tmp = 0.;
              for (unsigned int j = 0; j < dim; ++j)
                div_tmp += grad_q[i][j][j];

              AssertThrow(div_tmp == div_q[i], ExcInternalError());
            }
        }
    }

  // compare function values/gradients/divergences:
  for (unsigned int q = 0; q < quadrature.size(); ++q)
    {
      AssertThrow(values_manual[q] == values[q], ExcInternalError());
      AssertThrow(gradients_manual[q] == gradients[q], ExcInternalError());
      AssertThrow(divergences_manual[q] == divergences[q], ExcInternalError());
    }
}



template <int dim>
void
test_hyper_cube()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);

  FESystem<dim> fe(FE_Q<dim>(1), Tensor<2, dim>::n_independent_components);
  test(tr, fe);
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test_hyper_cube<1>();
  test_hyper_cube<2>();
  test_hyper_cube<3>();
}

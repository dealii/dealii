// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <iostream>

using namespace dealii;


template <int dim>
class MixedElastoPlasticity
{
public:
  MixedElastoPlasticity(const unsigned int degree);
  void run();
private:
  void make_grid_and_dofs();
  void assemble_system();

  const unsigned int degree;
  const unsigned int n_stress_components; // components of stress
  const unsigned int n_gamma_components; // scalar plastic multiplier

  Triangulation<dim> triangulation;
  FESystem<dim> fe;
  DoFHandler<dim> dof_handler;

  BlockSparsityPattern sparsity_pattern;
  BlockSparseMatrix<double> system_matrix;
  BlockVector<double> solution;
};



template <int dim>
MixedElastoPlasticity<dim>::MixedElastoPlasticity(const unsigned int degree):
  degree(degree),
  n_stress_components((dim*(dim + 1)) / 2),
  n_gamma_components(1),
  fe( FE_Q<dim>(degree), n_stress_components,
      FE_Q<dim>(degree), n_gamma_components),
  dof_handler(triangulation)
{
}



template <int dim>
void MixedElastoPlasticity<dim>::make_grid_and_dofs()
{
  GridGenerator::hyper_cube(triangulation, 0, 1);
  triangulation.refine_global(0);

  dof_handler.distribute_dofs(fe);

  deallog << "Number of stress components: "
          << n_stress_components
          << "\tNumber of gamma components: "
          << n_gamma_components << std::endl;

  // stress -> 0 gamma -> 1
  std::vector<unsigned int> block_component(n_stress_components + n_gamma_components, 1);
  for (unsigned int ii = 0; ii < n_stress_components; ii++)
    block_component[ii] = 0;

  DoFRenumbering::component_wise(dof_handler);

  std::vector<types::global_dof_index> dofs_per_block(2);
  DoFTools::count_dofs_per_block(dof_handler, dofs_per_block, block_component);

  const unsigned int n_stress_dof = dofs_per_block[0];
  const unsigned int n_gamma_dof = dofs_per_block[1];

  deallog << "Number of active cells: "
          << triangulation.n_active_cells()
          << std::endl
          << "Total number of cells: "
          << triangulation.n_cells()
          << std::endl
          << "Number of degrees of freedom: "
          << dof_handler.n_dofs()
          << " = (" << n_stress_dof << " + "
          << n_gamma_dof << ")"
          << std::endl;

  {
    BlockCompressedSimpleSparsityPattern csp(2, 2);

    csp.block(0, 0).reinit(n_stress_dof, n_stress_dof);
    csp.block(1, 0).reinit(n_gamma_dof, n_stress_dof);
    csp.block(0, 1).reinit(n_stress_dof, n_gamma_dof);
    csp.block(1, 1).reinit(n_gamma_dof, n_gamma_dof);
    csp.collect_sizes();

    DoFTools::make_sparsity_pattern(dof_handler, csp);
    sparsity_pattern.copy_from(csp);
  }

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(2);
  solution.block(0).reinit(n_stress_dof);
  solution.block(1).reinit(n_gamma_dof);
  solution.collect_sizes();
}



template <int dim>
void MixedElastoPlasticity<dim>::assemble_system()
{
  QGauss<dim> quadrature_formula(1);

  FEValues<dim> fe_values(fe, quadrature_formula,
                          update_values | update_gradients
                          | update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  deallog << "dofs_per_cell: " << fe.dofs_per_cell << std::endl;

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  const FEValuesExtractors::SymmetricTensor<2> stress(0);
  const FEValuesExtractors::Scalar gamma(n_stress_components);

  deallog   << "fe.dofs_per_cell: " << fe.dofs_per_cell
            << "\tquadrature_formula.size(): " << quadrature_formula.size()
            << std::endl;

  // constant stress and constant gamma distribution accross cell
  const double stress_value = 125.0;
  const double gamma_value = 1.0;

  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
  {
    fe_values.reinit (cell);
    cell->get_dof_indices (local_dof_indices);
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        // 0 = stress (symmetric 2nd order tensor), 1 = Gamma (scalar) interpolation fields
        const unsigned int i_group = fe.system_to_base_index(i).first.first;
        // the index of the non-zero shape function
        const unsigned int i_index = fe.system_to_base_index(i).first.second;
        // the support point (node) id
        const unsigned int i_node = fe.system_to_base_index(i).second;

        deallog << "\t" << i << "\t" << i_group << "\t" << i_index << "\t" << i_node;


        if (i_group == 0) // if i corresponds to tensor stress
          solution(local_dof_indices[i]) = stress_value;
        else  // i corresponds to scalar gamma
          solution(local_dof_indices[i]) = gamma_value;

        deallog << "\t" << solution(local_dof_indices[i]) << std::endl;
      }


    std::vector<Tensor<1,dim> > local_divergences (quadrature_formula.size());
    std::vector<SymmetricTensor<2,dim> > local_values (quadrature_formula.size());
    std::vector<double > local_scalar_values (quadrature_formula.size());

    fe_values.reinit (dof_handler.begin_active());
    fe_values[stress].get_function_values (solution, local_values);
    fe_values[stress].get_function_divergences (solution, local_divergences);
    fe_values[gamma].get_function_values (solution, local_scalar_values);

    // expect the nodal stress value at the quadrature point
    // constant stress field, therefore expect zero divergence
    // constant gamma field, therefore nodal value should equal value at quadrature point
    for (unsigned int q=0; q<quadrature_formula.size(); ++q)
      {
        deallog << local_values[q]
                << std::endl
                << local_divergences[q]
                << std::endl
                << local_scalar_values[q]
                << std::endl;

        for (unsigned int m = 0; m < dim; m++)
          {
            for (unsigned int n = 0; n< dim; n++)
              Assert ( (local_values[q])[m][n] == stress_value, ExcInternalError());

            Assert ( (local_divergences[q])[m] == 0.0, ExcInternalError());
            Assert ( local_scalar_values[q] == gamma_value, ExcInternalError());

          }
      }
  }
}



template <int dim>
void MixedElastoPlasticity<dim>::run()
{
  make_grid_and_dofs();
  assemble_system();
}



int main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (3);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  MixedElastoPlasticity < 3 > elasto_plasticity(1);
  elasto_plasticity.run();

  return 0;
}


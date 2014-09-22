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



// a simple example code demonstrating the usage of tensor-valued nodal unknowns
// in deal.II

#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
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
  BlockVector<double> system_rhs;

};

template <int dim>
MixedElastoPlasticity<dim>::MixedElastoPlasticity(const unsigned int degree):
  degree(degree),
  n_stress_components(dim *dim),
  n_gamma_components(1),
  fe( FE_Q<dim>(degree), n_stress_components,
      FE_Q<dim>(degree), n_gamma_components),
  dof_handler(triangulation)
{
}

template <int dim>
void MixedElastoPlasticity<dim>::make_grid_and_dofs()
{
  GridGenerator::hyper_cube(triangulation, -1, 1);
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

  // total number of dof per block component
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

  // following step-22 use of simple
  // compressed block sparsity
  // pattern for efficiency
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

  system_rhs.reinit(2);
  system_rhs.block(0).reinit(n_stress_dof);
  system_rhs.block(1).reinit(n_gamma_dof);
  system_rhs.collect_sizes();

  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    system_rhs(i) = i;
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

  FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> local_rhs(dofs_per_cell);

  std::vector<unsigned int> local_dof_indices(dofs_per_cell);


  const FEValuesExtractors::Tensor<2> stress(0);
  const FEValuesExtractors::Scalar gamma(n_stress_components);

  deallog   << "fe.dofs_per_cell: " << fe.dofs_per_cell
            << "\tquadrature_formula.size(): " << quadrature_formula.size()
            << std::endl;

  std::vector<Tensor<1,dim> > local_divergences (quadrature_formula.size());
  std::vector<Tensor<2,dim> > local_values (quadrature_formula.size());


  fe_values.reinit (dof_handler.begin_active());
  fe_values[stress].get_function_values (system_rhs, local_values);
  fe_values[stress].get_function_divergences (system_rhs, local_divergences);

  for (unsigned int q=0; q<quadrature_formula.size(); ++q)
    deallog << local_values[q]
            << std::endl
            << local_divergences[q]
            << std::endl;
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
}


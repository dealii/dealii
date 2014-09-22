// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
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



// a complete testcase by Denis Davydov for FEValuesExtractors::Tensor

#include "../tests.h"
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <base/function.h>
#include <lac/block_vector.h>
#include <lac/full_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>
#include <numerics/vector_tools.h>
#include <numerics/matrix_tools.h>
#include <numerics/data_out.h>
#include <fe/fe_q.h>

#include <fstream>

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

  // following step-22 use of simple compressed block sparsity pattern for efficiency
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
}

template <int dim>
void MixedElastoPlasticity<dim>::assemble_system()
{
  QGauss<dim> quadrature_formula(1);

  FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  deallog << "dofs_per_cell: " << fe.dofs_per_cell << std::endl;
  //return;
  const unsigned int n_q_points = quadrature_formula.size();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);

  std::vector<unsigned int> local_dof_indices(dofs_per_cell);


  const FEValuesExtractors::Tensor<2> stress_extr(0);//rank2
  const FEValuesExtractors::Scalar gamma_extr(n_stress_components);

  deallog   << "fe.dofs_per_cell: " << fe.dofs_per_cell
            << "\tquadrature_formula.size(): " << quadrature_formula.size()
            << std::endl;

  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                 endc = dof_handler.end();


  std::vector<Tensor<1,dim> > div_values(n_q_points);
  std::vector<Tensor<2,dim> > stress_values(n_q_points);

  unsigned int cc = 0;
  for (; cell!=endc; ++cell) //loop over all cells
    {
      deallog<<++cc<<" ";
      cell_matrix = 0;
      cell_rhs = 0;

      fe_values.reinit (cell); //compute requested values for a given cell

      fe_values[stress_extr].get_function_values(solution,stress_values);
      fe_values[stress_extr].get_function_divergences(solution,div_values);
    }
}

template <int dim>
void MixedElastoPlasticity<dim>::run()
{
  make_grid_and_dofs();
  assemble_system();
  deallog << "here" << std::endl;

}

void check()
{
  static const unsigned int dim = 3;
  {
    for (unsigned int i = 0; i < dim; i++)
      {
        for (unsigned int j = 0; j < dim; j++)
          {
            TableIndices<2> indices(i,j);
            unsigned int unrolled = Tensor<2,dim>::component_to_unrolled_index(indices);
            deallog<<i<<" "<<j<<" -> "<<unrolled<<std::endl;
            indices = Tensor<2,dim>::unrolled_to_component_indices(unrolled);
            deallog<<unrolled<<" -> "<<indices[0]<<" "<<indices[1]<<std::endl;
          }
      }
  }

  MixedElastoPlasticity < 3 > elasto_plasticity(1);
  elasto_plasticity.run();

  deallog << "Analysis complete" << std::endl;

}


int main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (3);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  check ();
}

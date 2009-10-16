//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009 by Andrew McBride and the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


// a simple example code demonstrating the usage of tensor-valued nodal unknowns
// in deal.II

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
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <fe/fe_q.h>


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
  std::vector<unsigned int> dofs_per_block(2);
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


  const FEValuesExtractors::SymmetricTensor<2> stress(0);
  const FEValuesExtractors::Scalar gamma(n_stress_components);

  deallog   << "fe.dofs_per_cell: " << fe.dofs_per_cell
	    << "\tquadrature_formula.size(): " << quadrature_formula.size()
	    << std::endl;

  std::vector<Tensor<1,dim> > local_divergences (quadrature_formula.size());
  std::vector<SymmetricTensor<2,dim> > local_values (quadrature_formula.size());


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
  std::ofstream logfile ("fe_values_view_21/output");
  deallog << std::setprecision (3);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  MixedElastoPlasticity < 3 > elasto_plasticity(1);
  elasto_plasticity.run();
}


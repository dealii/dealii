// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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



// this function tests the correctness of the implementation of matrix free
// matrix-vector products by comparing with the result of deal.II sparse
// matrix. The mesh uses a hypershell mesh with hanging nodes and constraints
// between the vector components in the form of no-normal flux constraints on
// the Stokes equations.

#include "../tests.h"

std::ofstream logfile("output");

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>
#include <complex>



template <int dim, int degree_p, typename VectorType>
class MatrixFreeTest
{
public:
  typedef typename DoFHandler<dim>::active_cell_iterator CellIterator;
  typedef double Number;

  MatrixFreeTest(const MatrixFree<dim,Number> &data_in):
    data (data_in)
  {};

  void
  local_apply (const MatrixFree<dim,Number> &data,
               VectorType          &dst,
               const VectorType    &src,
               const std::pair<unsigned int,unsigned int> &cell_range) const
  {
    typedef VectorizedArray<Number> vector_t;
    FEEvaluation<dim,degree_p+1,degree_p+2,dim,Number> velocity (data, 0);
    FEEvaluation<dim,degree_p,  degree_p+2,1,  Number> pressure (data, 1);

    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        velocity.reinit (cell);
        velocity.read_dof_values (src.block(0));
        velocity.evaluate (false,true,false);
        pressure.reinit (cell);
        pressure.read_dof_values (src.block(1));
        pressure.evaluate (true,false,false);

        for (unsigned int q=0; q<velocity.n_q_points; ++q)
          {
            SymmetricTensor<2,dim,vector_t> sym_grad_u =
              velocity.get_symmetric_gradient (q);
            vector_t pres = pressure.get_value(q);
            vector_t div = -trace(sym_grad_u);
            pressure.submit_value   (div, q);

            // subtract p * I
            for (unsigned int d=0; d<dim; ++d)
              sym_grad_u[d][d] -= pres;

            velocity.submit_symmetric_gradient(sym_grad_u, q);
          }

        velocity.integrate (false,true);
        velocity.distribute_local_to_global (dst.block(0));
        pressure.integrate (true,false);
        pressure.distribute_local_to_global (dst.block(1));
      }
  }


  void vmult (VectorType &dst,
              const VectorType &src) const
  {
    dst = 0;
    data.cell_loop (&MatrixFreeTest<dim,degree_p,VectorType>::local_apply,
                    this, dst, src);
  };

private:
  const MatrixFree<dim,Number> &data;
};



template <int dim, int fe_degree>
void test ()
{
  Triangulation<dim>   triangulation;
  GridGenerator::hyper_shell (triangulation, Point<dim>(),
                              0.5, 1., 96, true);
  static HyperShellBoundary<dim> boundary;
  triangulation.set_boundary (0, boundary);
  triangulation.set_boundary (1, boundary);
  triangulation.begin_active()->set_refine_flag();
  triangulation.last()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();
  triangulation.refine_global (3-dim);
  triangulation.last()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();

  MappingQ<dim>        mapping (3);
  FE_Q<dim>            fe_u_scal (fe_degree+1);
  FESystem<dim>        fe_u (fe_u_scal,dim);
  FE_Q<dim>            fe_p (fe_degree);
  FESystem<dim>        fe (fe_u_scal, dim, fe_p, 1);
  DoFHandler<dim>      dof_handler_u (triangulation);
  DoFHandler<dim>      dof_handler_p (triangulation);
  DoFHandler<dim>      dof_handler (triangulation);

  MatrixFree<dim,double> mf_data;

  ConstraintMatrix     constraints, constraints_u, constraints_p;

  BlockSparsityPattern      sparsity_pattern;
  BlockSparseMatrix<double> system_matrix;

  BlockVector<double> solution;
  BlockVector<double> system_rhs;
  BlockVector<double> mf_solution;

  dof_handler.distribute_dofs (fe);
  dof_handler_u.distribute_dofs (fe_u);
  dof_handler_p.distribute_dofs (fe_p);
  std::vector<unsigned int> stokes_sub_blocks (dim+1,0);
  stokes_sub_blocks[dim] = 1;
  DoFRenumbering::component_wise (dof_handler, stokes_sub_blocks);

  std::set<unsigned char> no_normal_flux_boundaries;
  no_normal_flux_boundaries.insert (0);
  no_normal_flux_boundaries.insert (1);
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           constraints);
  VectorTools::compute_no_normal_flux_constraints (dof_handler, 0,
                                                   no_normal_flux_boundaries,
                                                   constraints, mapping);
  constraints.close ();
  DoFTools::make_hanging_node_constraints (dof_handler_u,
                                           constraints_u);
  VectorTools::compute_no_normal_flux_constraints (dof_handler_u, 0,
                                                   no_normal_flux_boundaries,
                                                   constraints_u, mapping);
  constraints_u.close ();
  DoFTools::make_hanging_node_constraints (dof_handler_p,
                                           constraints_p);
  constraints_p.close ();

  std::vector<types::global_dof_index> dofs_per_block (2);
  DoFTools::count_dofs_per_block (dof_handler, dofs_per_block,
                                  stokes_sub_blocks);

  //std::cout << "Number of active cells: "
  //          << triangulation.n_active_cells()
  //          << std::endl
  //          << "Number of degrees of freedom: "
  //          << dof_handler.n_dofs()
  //          << " (" << n_u << '+' << n_p << ')'
  //          << std::endl;

  {
    BlockCompressedSimpleSparsityPattern csp (2,2);

    for (unsigned int d=0; d<2; ++d)
      for (unsigned int e=0; e<2; ++e)
        csp.block(d,e).reinit (dofs_per_block[d], dofs_per_block[e]);

    csp.collect_sizes();

    DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
    sparsity_pattern.copy_from (csp);
  }

  system_matrix.reinit (sparsity_pattern);

  // this is from step-22
  {
    QGauss<dim>   quadrature_formula(fe_degree+2);

    FEValues<dim> fe_values (mapping, fe, quadrature_formula,
                             update_values    |
                             update_JxW_values |
                             update_gradients);

    const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
    const unsigned int   n_q_points      = quadrature_formula.size();

    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);

    std::vector<SymmetricTensor<2,dim> > phi_grads_u (dofs_per_cell);
    std::vector<double>                  div_phi_u   (dofs_per_cell);
    std::vector<double>                  phi_p       (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);
        local_matrix = 0;

        for (unsigned int q=0; q<n_q_points; ++q)
          {
            for (unsigned int k=0; k<dofs_per_cell; ++k)
              {
                phi_grads_u[k] = fe_values[velocities].symmetric_gradient (k, q);
                div_phi_u[k]   = fe_values[velocities].divergence (k, q);
                phi_p[k]       = fe_values[pressure].value (k, q);
              }

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              {
                for (unsigned int j=0; j<=i; ++j)
                  {
                    local_matrix(i,j) += (phi_grads_u[i] * phi_grads_u[j]
                                          - div_phi_u[i] * phi_p[j]
                                          - phi_p[i] * div_phi_u[j])
                                         * fe_values.JxW(q);
                  }
              }
          }
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=i+1; j<dofs_per_cell; ++j)
            local_matrix(i,j) = local_matrix(j,i);

        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global (local_matrix,
                                                local_dof_indices,
                                                system_matrix);
      }
  }


  solution.reinit (2);
  for (unsigned int d=0; d<2; ++d)
    solution.block(d).reinit (dofs_per_block[d]);
  solution.collect_sizes ();

  system_rhs.reinit (solution);
  mf_solution.reinit (solution);

  // fill system_rhs with random numbers
  for (unsigned int j=0; j<system_rhs.block(0).size(); ++j)
    if (constraints_u.is_constrained(j) == false)
      {
        const double val = -1 + 2.*(double)Testing::rand()/double(RAND_MAX);
        system_rhs.block(0)(j) = val;
      }
  for (unsigned int j=0; j<system_rhs.block(1).size(); ++j)
    if (constraints_p.is_constrained(j) == false)
      {
        const double val = -1 + 2.*(double)Testing::rand()/double(RAND_MAX);
        system_rhs.block(1)(j) = val;
      }

  // setup matrix-free structure
  {
    std::vector<const DoFHandler<dim>*> dofs;
    dofs.push_back(&dof_handler_u);
    dofs.push_back(&dof_handler_p);
    std::vector<const ConstraintMatrix *> constraints;
    constraints.push_back (&constraints_u);
    constraints.push_back (&constraints_p);
    QGauss<1> quad(fe_degree+2);
    // no parallelism
    mf_data.reinit (mapping, dofs, constraints, quad,
                    typename MatrixFree<dim>::AdditionalData
                    (MPI_COMM_WORLD,
                     MatrixFree<dim>::AdditionalData::none));
  }

  system_matrix.vmult (solution, system_rhs);

  MatrixFreeTest<dim,fe_degree,BlockVector<double> > mf (mf_data);
  mf.vmult (mf_solution, system_rhs);

  // Verification
  mf_solution -= solution;
  const double error = mf_solution.linfty_norm();
  const double relative = solution.linfty_norm();
  deallog << "Verification fe degree " << fe_degree  <<  ": "
          << error/relative << std::endl << std::endl;
}



int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  deallog << std::setprecision (3);

  {
    deallog << std::endl << "Test with doubles" << std::endl << std::endl;
    deallog.threshold_double(1.e-12);
    deallog.push("2d");
    test<2,1>();
    test<2,2>();
    test<2,3>();
    deallog.pop();
    deallog.push("3d");
    test<3,1>();
    deallog.pop();
  }
}


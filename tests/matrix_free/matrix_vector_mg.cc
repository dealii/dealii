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
// matrix for MG DoFHandler on a hyperball mesh with no hanging nodes but
// homogeneous Dirichlet conditions

#include "../tests.h"

std::ofstream logfile("output");

#include "matrix_vector_common.h"
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/base/function.h>



template <int dim, int fe_degree>
void test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_ball (tria);
  static const HyperBallBoundary<dim> boundary;
  tria.set_boundary (0, boundary);
  tria.refine_global(5-dim);

  FE_Q<dim> fe (fe_degree);

  // setup DoFs
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  dof.distribute_mg_dofs(fe);
  ConstraintMatrix constraints;
  VectorTools::interpolate_boundary_values (dof,
                                            0,
                                            ZeroFunction<dim>(),
                                            constraints);
  constraints.close ();

  //std::cout << "Number of cells: " << dof.get_tria().n_active_cells() << std::endl;
  //std::cout << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;

  // set up MatrixFree
  QGauss<1> quad (fe_degree+1);
  MatrixFree<dim> mf_data;
  mf_data.reinit (dof, constraints, quad);
  SparsityPattern sparsity;
  SparseMatrix<double> system_matrix;
  {
    CompressedSimpleSparsityPattern csp (dof.n_dofs(), dof.n_dofs());
    DoFTools::make_sparsity_pattern (static_cast<const DoFHandler<dim>&>(dof),
                                     csp, constraints, false);
    sparsity.copy_from (csp);
  }
  system_matrix.reinit (sparsity);

  // setup MG levels
  const unsigned int nlevels = tria.n_levels();
  typedef MatrixFree<dim> MatrixFreeTestType;
  MGLevelObject<MatrixFreeTestType>    mg_matrices;
  MGLevelObject<ConstraintMatrix>  mg_constraints;
  MGLevelObject<SparsityPattern>   mg_sparsities;
  MGLevelObject<SparseMatrix<double> > mg_ref_matrices;
  mg_matrices.resize (0, nlevels-1);
  mg_constraints.resize (0, nlevels-1);
  mg_sparsities.resize (0, nlevels-1);
  mg_ref_matrices.resize (0, nlevels-1);

  typename FunctionMap<dim>::type dirichlet_boundary;
  ZeroFunction<dim>               homogeneous_dirichlet_bc (1);
  dirichlet_boundary[0] = &homogeneous_dirichlet_bc;
  std::vector<std::set<types::global_dof_index> > boundary_indices(nlevels);
  MGTools::make_boundary_list (dof, dirichlet_boundary, boundary_indices);
  for (unsigned int level=0; level<nlevels; ++level)
    {
      std::set<types::global_dof_index>::iterator bc_it = boundary_indices[level].begin();
      for ( ; bc_it != boundary_indices[level].end(); ++bc_it)
        mg_constraints[level].add_line(*bc_it);
      mg_constraints[level].close();
      typename MatrixFree<dim>::AdditionalData data;
      data.level_mg_handler = level;
      mg_matrices[level].reinit(dof, mg_constraints[level], quad, data);

      CompressedSimpleSparsityPattern csp;
      csp.reinit (dof.n_dofs(level), dof.n_dofs(level));
      MGTools::make_sparsity_pattern (dof, csp, level);
      mg_sparsities[level].copy_from (csp);
      mg_ref_matrices[level].reinit (mg_sparsities[level]);
    }

  // assemble sparse matrix with (\nabla v,
  // \nabla u) + (v, 10 * u) on the actual
  // discretization and on all levels
  {
    QGauss<dim> quad (fe_degree+1);
    FEValues<dim> fe_values (fe, quad,
                             update_values    |  update_gradients |
                             update_JxW_values);
    const unsigned int n_quadrature_points = quad.size();
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof.begin_active(),
    endc = dof.end();
    for (; cell!=endc; ++cell)
      {
        cell_matrix = 0;
        fe_values.reinit (cell);

        for (unsigned int q_point=0; q_point<n_quadrature_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                cell_matrix(i,j) += ((fe_values.shape_grad(i,q_point) *
                                      fe_values.shape_grad(j,q_point) +
                                      10. * fe_values.shape_value(i,q_point) *
                                      fe_values.shape_value(j,q_point)) *
                                     fe_values.JxW(q_point));
            }
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global (cell_matrix,
                                                local_dof_indices,
                                                system_matrix);
      }

    // now to the MG assembly
    typename DoFHandler<dim>::cell_iterator
    cellm = dof.begin(),
    endcm = dof.end();
    for (; cellm!=endcm; ++cellm)
      {
        cell_matrix = 0;
        fe_values.reinit (cellm);

        for (unsigned int q_point=0; q_point<n_quadrature_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                cell_matrix(i,j) += ((fe_values.shape_grad(i,q_point) *
                                      fe_values.shape_grad(j,q_point) +
                                      10. * fe_values.shape_value(i,q_point) *
                                      fe_values.shape_value(j,q_point)) *
                                     fe_values.JxW(q_point));
            }
        cellm->get_mg_dof_indices (local_dof_indices);
        mg_constraints[cellm->level()]
        .distribute_local_to_global (cell_matrix,
                                     local_dof_indices,
                                     mg_ref_matrices[cellm->level()]);
      }
  }

  // fill a right hand side vector with random
  // numbers in unconstrained degrees of freedom
  Vector<double> src (dof.n_dofs());
  Vector<double> result_spmv(src), result_mf (src);

  for (unsigned int i=0; i<dof.n_dofs(); ++i)
    {
      if (constraints.is_constrained(i) == false)
        src(i) = (double)Testing::rand()/RAND_MAX;
    }

  // now perform matrix-vector product and check
  // its correctness
  system_matrix.vmult (result_spmv, src);
  MatrixFreeTest<dim,fe_degree,double> mf (mf_data);
  mf.vmult (result_mf, src);

  result_mf -= result_spmv;
  const double diff_norm = result_mf.linfty_norm();
  deallog << "Norm of difference active: " << diff_norm << std::endl;

  for (unsigned int level=0; level<nlevels; ++level)
    {
      Vector<double> src (dof.n_dofs(level));
      Vector<double> result_spmv(src), result_mf (src);

      for (unsigned int i=0; i<dof.n_dofs(level); ++i)
        {
          if (mg_constraints[level].is_constrained(i) == false)
            src(i) = (double)Testing::rand()/RAND_MAX;
        }

      // now perform matrix-vector product and check
      // its correctness
      mg_ref_matrices[level].vmult (result_spmv, src);
      MatrixFreeTest<dim,fe_degree,double> mf_lev (mg_matrices[level]);
      mf_lev.vmult (result_mf, src);

      result_mf -= result_spmv;
      const double diff_norm = result_mf.linfty_norm();
      deallog << "Norm of difference MG level " << level
              << ": " << diff_norm << std::endl;
    }
  deallog << std::endl;
}


// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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


// Test the correctness of an implementation for the Stokes operator
// based on MatrixFreeOperators::Base by comparing with a matrix version.

#include <deal.II/base/multithread_info.h>

#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/operators.h>

#include "../tests.h"


template <int dim, int degree_p, typename BlockVectorType>
class MatrixFreeTest : public MatrixFreeOperators::Base<dim, BlockVectorType>
{
public:
  typedef typename BlockVectorType::value_type                     Number;
  typedef typename MatrixFreeOperators::Base<dim, BlockVectorType> Base;

  void
  compute_diagonal()
  {
    AssertThrow(false, ExcNotImplemented());
  }

protected:
  void
  apply_add(BlockVectorType &dst, const BlockVectorType &src) const
  {
    Base::data->cell_loop(&MatrixFreeTest::local_apply_cell, this, dst, src);
  }

  void
  local_apply_cell(
    const dealii::MatrixFree<dim, Number> &      data,
    BlockVectorType &                            dst,
    const BlockVectorType &                      src,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    typedef VectorizedArray<Number>                            vector_t;
    FEEvaluation<dim, degree_p + 1, degree_p + 2, dim, Number> velocity(data,
                                                                        0);
    FEEvaluation<dim, degree_p, degree_p + 2, 1, Number> pressure(data, 1);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        velocity.reinit(cell);
        velocity.read_dof_values(src.block(0));
        velocity.evaluate(false, true, false);
        pressure.reinit(cell);
        pressure.read_dof_values(src.block(1));
        pressure.evaluate(true, false, false);

        for (unsigned int q = 0; q < velocity.n_q_points; ++q)
          {
            SymmetricTensor<2, dim, vector_t> sym_grad_u =
              velocity.get_symmetric_gradient(q);
            vector_t pres = pressure.get_value(q);
            vector_t div  = -trace(sym_grad_u);
            pressure.submit_value(div, q);

            // subtract p * I
            for (unsigned int d = 0; d < dim; ++d)
              sym_grad_u[d][d] -= pres;

            velocity.submit_symmetric_gradient(sym_grad_u, q);
          }

        velocity.integrate(false, true);
        velocity.distribute_local_to_global(dst.block(0));
        pressure.integrate(true, false);
        pressure.distribute_local_to_global(dst.block(1));
      }
  }
};



template <int dim, int degree>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FE_Q<dim>     fe_u_scal(degree + 1);
  FESystem<dim> fe_u(fe_u_scal, dim);
  FE_Q<dim>     fe_p(degree);
  FESystem<dim> fe(fe_u_scal, dim, fe_p, 1);

  DoFHandler<dim> dof(tria);
  DoFHandler<dim> dof_u(tria);
  DoFHandler<dim> dof_p(tria);

  dof.distribute_dofs(fe);
  dof_u.distribute_dofs(fe_u);
  dof_p.distribute_dofs(fe_p);

  std::shared_ptr<MatrixFree<dim, double>> mf_data;

  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints, constraints_u, constraints_p;

  BlockSparsityPattern      sparsity_pattern;
  BlockSparseMatrix<double> system_matrix;

  BlockVector<double>                             solution;
  BlockVector<double>                             system_rhs;
  LinearAlgebra::distributed::BlockVector<double> mf_system_rhs;
  LinearAlgebra::distributed::BlockVector<double> mf_solution;

  std::vector<unsigned int> stokes_sub_blocks(dim + 1, 0);
  stokes_sub_blocks[dim] = 1;
  DoFRenumbering::component_wise(dof, stokes_sub_blocks);

  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();
  DoFTools::make_hanging_node_constraints(dof_u, constraints_u);
  constraints_u.close();
  DoFTools::make_hanging_node_constraints(dof_p, constraints_p);
  constraints_p.close();

  const std::vector<types::global_dof_index> dofs_per_block =
    DoFTools::count_dofs_per_fe_block(dof, stokes_sub_blocks);
  {
    BlockDynamicSparsityPattern csp(2, 2);

    for (unsigned int d = 0; d < 2; ++d)
      {
        for (unsigned int e = 0; e < 2; ++e)
          csp.block(d, e).reinit(dofs_per_block[d], dofs_per_block[e]);
      }

    csp.collect_sizes();

    DoFTools::make_sparsity_pattern(dof, csp, constraints, false);
    sparsity_pattern.copy_from(csp);
  }

  system_matrix.reinit(sparsity_pattern);

  // this is from step-22
  {
    QGauss<dim> quadrature_formula(degree + 2);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_JxW_values |
                              update_gradients);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar pressure(dim);

    std::vector<SymmetricTensor<2, dim>> phi_grads_u(dofs_per_cell);
    std::vector<double>                  div_phi_u(dofs_per_cell);
    std::vector<double>                  phi_p(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
                                                   endc = dof.end();
    for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);
        local_matrix = 0;

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            for (unsigned int k = 0; k < dofs_per_cell; ++k)
              {
                phi_grads_u[k] = fe_values[velocities].symmetric_gradient(k, q);
                div_phi_u[k]   = fe_values[velocities].divergence(k, q);
                phi_p[k]       = fe_values[pressure].value(k, q);
              }

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  {
                    local_matrix(i, j) +=
                      (phi_grads_u[i] * phi_grads_u[j] -
                       div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j]) *
                      fe_values.JxW(q);
                  }
              }
          }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(local_matrix,
                                               local_dof_indices,
                                               system_matrix);
      }
  }

  solution.reinit(2);
  mf_solution.reinit(2);
  for (unsigned int d = 0; d < 2; ++d)
    {
      solution.block(d).reinit(dofs_per_block[d]);
      mf_solution.block(d).reinit(dofs_per_block[d]);
    }
  solution.collect_sizes();
  mf_solution.collect_sizes();
  system_rhs.reinit(solution);
  mf_system_rhs.reinit(mf_solution);

  // fill system_rhs with random numbers
  for (unsigned int j = 0; j < system_rhs.block(0).size(); ++j)
    if (constraints_u.is_constrained(j) == false)
      {
        const double val          = -1 + 2. * random_value<double>();
        system_rhs.block(0)(j)    = val;
        mf_system_rhs.block(0)(j) = val;
      }
  for (unsigned int j = 0; j < system_rhs.block(1).size(); ++j)
    if (constraints_p.is_constrained(j) == false)
      {
        const double val          = -1 + 2. * random_value<double>();
        system_rhs.block(1)(j)    = val;
        mf_system_rhs.block(1)(j) = val;
      }

  mf_data =
    std::shared_ptr<MatrixFree<dim, double>>(new MatrixFree<dim, double>());
  // setup matrix-free structure
  {
    std::vector<const DoFHandler<dim> *> dofs;
    dofs.push_back(&dof_u);
    dofs.push_back(&dof_p);
    std::vector<const AffineConstraints<double> *> constraints;
    constraints.push_back(&constraints_u);
    constraints.push_back(&constraints_p);
    QGauss<1> quad(degree + 2);
    // no parallelism
    mf_data->reinit(dofs,
                    constraints,
                    quad,
                    typename MatrixFree<dim>::AdditionalData(
                      MatrixFree<dim>::AdditionalData::none));
  }
  system_matrix.vmult(solution, system_rhs);

  MatrixFreeTest<dim, degree, LinearAlgebra::distributed::BlockVector<double>>
    mf;
  mf.initialize(mf_data);
  mf.vmult(mf_solution, mf_system_rhs);

  // Verification
  for (unsigned int i = 0; i < mf_solution.size(); ++i)
    mf_solution(i) -= solution(i);
  const double error    = mf_solution.linfty_norm();
  const double relative = solution.linfty_norm();
  deallog << "Verification fe degree " << degree << ": " << error / relative
          << std::endl
          << std::endl;
}



int
main(int /*argc*/, char ** /*argv*/)
{
  initlog();
  deallog.push("2D");
  test<2, 1>();
  test<2, 2>();
  deallog.pop();
  deallog.push("3D");
  test<3, 1>();
  test<3, 2>();
  deallog.pop();

  return 0;
}

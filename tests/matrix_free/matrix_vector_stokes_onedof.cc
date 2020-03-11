// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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



// correctness matrix free matrix-vector products by comparing with the result
// of a deal.II sparse matrix. Similar to matrix_vector_stokes_noflux but
// putting all degrees of freedom into a single DoFHandler, where the
// selection is done through FEEvaluation

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include <complex>
#include <iostream>

#include "../tests.h"



template <int dim, int degree_p, typename VectorType>
class MatrixFreeTest
{
public:
  typedef typename DoFHandler<dim>::active_cell_iterator CellIterator;
  typedef double                                         Number;

  MatrixFreeTest(const MatrixFree<dim, Number> &data_in)
    : data(data_in){};

  void
  local_apply(const MatrixFree<dim, Number> &              data,
              VectorType &                                 dst,
              const VectorType &                           src,
              const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    typedef VectorizedArray<Number>                            vector_t;
    FEEvaluation<dim, degree_p + 1, degree_p + 2, dim, Number> velocity(data,
                                                                        0,
                                                                        0,
                                                                        0);
    FEEvaluation<dim, degree_p, degree_p + 2, 1, Number>       pressure(data,
                                                                  0,
                                                                  0,
                                                                  dim);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        velocity.reinit(cell);
        velocity.read_dof_values(src);
        velocity.evaluate(false, true, false);
        pressure.reinit(cell);
        pressure.read_dof_values(src);
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
        velocity.distribute_local_to_global(dst);
        pressure.integrate(true, false);
        pressure.distribute_local_to_global(dst);
      }
  }


  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    dst = 0;
    data.cell_loop(&MatrixFreeTest<dim, degree_p, VectorType>::local_apply,
                   this,
                   dst,
                   src);
  };

private:
  const MatrixFree<dim, Number> &data;
};



template <int dim, int fe_degree>
void
test(const FESystem<dim> &fe)
{
  SphericalManifold<dim> manifold;
  Triangulation<dim>     triangulation;
  GridGenerator::hyper_shell(triangulation, Point<dim>(), 0.5, 1., 96, true);
  triangulation.set_all_manifold_ids(0);
  triangulation.set_manifold(0, manifold);
  triangulation.begin_active()->set_refine_flag();
  triangulation.last()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();
  triangulation.refine_global(3 - dim);
  triangulation.last()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();

  MappingQ<dim>   mapping(3);
  DoFHandler<dim> dof_handler(triangulation);

  MatrixFree<dim, double> mf_data;

  AffineConstraints<double> constraints, constraints_u, constraints_p;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> solution;
  Vector<double> system_rhs;
  Vector<double> mf_solution;

  dof_handler.distribute_dofs(fe);

  std::set<types::boundary_id> no_normal_flux_boundaries;
  no_normal_flux_boundaries.insert(0);
  no_normal_flux_boundaries.insert(1);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  if (fe.dofs_per_vertex > 0)
    VectorTools::compute_no_normal_flux_constraints(
      dof_handler, 0, no_normal_flux_boundaries, constraints, mapping);
  constraints.close();

  // std::cout << "Number of active cells: "
  //          << triangulation.n_active_cells()
  //          << std::endl
  //          << "Number of degrees of freedom: "
  //          << dof_handler.n_dofs()
  //          << " (" << n_u << '+' << n_p << ')'
  //          << std::endl;

  {
    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());

    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    sparsity_pattern.copy_from(dsp);
  }

  system_matrix.reinit(sparsity_pattern);

  // this is from step-22
  {
    QGauss<dim> quadrature_formula(fe_degree + 2);

    FEValues<dim> fe_values(mapping,
                            fe,
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

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                   endc = dof_handler.end();
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
                for (unsigned int j = 0; j <= i; ++j)
                  {
                    local_matrix(i, j) +=
                      (phi_grads_u[i] * phi_grads_u[j] -
                       div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j]) *
                      fe_values.JxW(q);
                  }
              }
          }
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = i + 1; j < dofs_per_cell; ++j)
            local_matrix(i, j) = local_matrix(j, i);

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(local_matrix,
                                               local_dof_indices,
                                               system_matrix);
      }
  }

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(solution);
  mf_solution.reinit(solution);

  // fill system_rhs with random numbers
  for (unsigned int j = 0; j < system_rhs.size(); ++j)
    if (constraints.is_constrained(j) == false)
      {
        const double val = -1 + 2. * (double)Testing::rand() / double(RAND_MAX);
        system_rhs(j)    = val;
      }

  // setup matrix-free structure
  {
    QGauss<1> quad(fe_degree + 2);
    // no parallelism
    mf_data.reinit(mapping,
                   dof_handler,
                   constraints,
                   quad,
                   typename MatrixFree<dim>::AdditionalData(
                     MatrixFree<dim>::AdditionalData::none));
  }

  system_matrix.vmult(solution, system_rhs);

  MatrixFreeTest<dim, fe_degree, Vector<double>> mf(mf_data);
  mf.vmult(mf_solution, system_rhs);

  // Verification
  mf_solution -= solution;
  const double error    = mf_solution.linfty_norm();
  const double relative = solution.linfty_norm();
  deallog << "Verification " << fe.get_name() << ": " << error / relative
          << std::endl
          << std::endl;
}



int
main()
{
  initlog();

  {
    deallog << std::endl << "Test with doubles" << std::endl << std::endl;
    deallog.push("2d");
    test<2, 1>(FESystem<2>(FE_Q<2>(2), 2, FE_Q<2>(1), 1));
    test<2, 2>(FESystem<2>(FE_Q<2>(3), 2, FE_Q<2>(2), 1));
    test<2, 3>(FESystem<2>(FE_Q<2>(4), 2, FE_Q<2>(3), 1));
    test<2, 1>(FESystem<2>(FE_DGQ<2>(2), 2, FE_DGQ<2>(1), 1));
    test<2, 0>(FESystem<2>(FE_DGQ<2>(1), 2, FE_DGQ<2>(0), 1));
    test<2, 3>(FESystem<2>(FE_DGQ<2>(4), 2, FE_DGQ<2>(3), 1));
    deallog.pop();
    deallog.push("3d");
    test<3, 1>(FESystem<3>(FE_Q<3>(2), 3, FE_Q<3>(1), 1));
    test<3, 1>(FESystem<3>(FE_DGQ<3>(2), 3, FE_DGQ<3>(1), 1));
    deallog.pop();
  }
}

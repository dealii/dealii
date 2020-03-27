/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2008 - 2018 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 */

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/packaged_operation.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/schur_complement.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <sstream>

#include "../tests.h"

namespace Step22
{
  template <int dim>
  class StokesProblem
  {
  public:
    StokesProblem(const unsigned int degree);
    void
    run();

  private:
    void
    setup_dofs();
    void
    assemble_system();
    void
    solve();
    void
                              refine_mesh();
    const unsigned int        degree;
    Triangulation<dim>        triangulation;
    FESystem<dim>             fe;
    DoFHandler<dim>           dof_handler;
    AffineConstraints<double> constraints;
    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;
    BlockVector<double>       solution;
    BlockVector<double>       system_rhs;
  };

  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    BoundaryValues()
      : Function<dim>(dim + 1)
    {}
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const;
    virtual void
    vector_value(const Point<dim> &p, Vector<double> &value) const;
  };

  template <int dim>
  double
  BoundaryValues<dim>::value(const Point<dim> & p,
                             const unsigned int component) const
  {
    Assert(component < this->n_components,
           ExcIndexRange(component, 0, this->n_components));
    if (component == 0)
      return (p[0] < 0 ? -1 : (p[0] > 0 ? 1 : 0));
    return 0;
  }

  template <int dim>
  void
  BoundaryValues<dim>::vector_value(const Point<dim> &p,
                                    Vector<double> &  values) const
  {
    for (unsigned int c = 0; c < this->n_components; ++c)
      values(c) = BoundaryValues<dim>::value(p, c);
  }

  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide()
      : Function<dim>(dim + 1)
    {}
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const;
    virtual void
    vector_value(const Point<dim> &p, Vector<double> &value) const;
  };

  template <int dim>
  double
  RightHandSide<dim>::value(const Point<dim> & /*p*/,
                            const unsigned int /*component*/) const
  {
    return 0;
  }

  template <int dim>
  void
  RightHandSide<dim>::vector_value(const Point<dim> &p,
                                   Vector<double> &  values) const
  {
    for (unsigned int c = 0; c < this->n_components; ++c)
      values(c) = RightHandSide<dim>::value(p, c);
  }

  template <int dim>
  StokesProblem<dim>::StokesProblem(const unsigned int degree)
    : degree(degree)
    , triangulation(Triangulation<dim>::maximum_smoothing)
    , fe(FE_Q<dim>(degree + 1), dim, FE_Q<dim>(degree), 1)
    , dof_handler(triangulation)
  {}

  template <int dim>
  void
  StokesProblem<dim>::setup_dofs()
  {
    system_matrix.clear();
    dof_handler.distribute_dofs(fe);
    DoFRenumbering::Cuthill_McKee(dof_handler);
    std::vector<unsigned int> block_component(dim + 1, 0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, block_component);
    {
      constraints.clear();
      FEValuesExtractors::Vector velocities(0);
      DoFTools::make_hanging_node_constraints(dof_handler, constraints);
      VectorTools::interpolate_boundary_values(dof_handler,
                                               1,
                                               BoundaryValues<dim>(),
                                               constraints,
                                               fe.component_mask(velocities));
    }
    constraints.close();
    const std::vector<types::global_dof_index> dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
    const unsigned int n_u = dofs_per_block[0], n_p = dofs_per_block[1];
    deallog << "   Number of active cells: " << triangulation.n_active_cells()
            << std::endl
            << "   Number of degrees of freedom: " << dof_handler.n_dofs()
            << " (" << n_u << '+' << n_p << ')' << std::endl;
    {
      BlockDynamicSparsityPattern dsp(2, 2);
      dsp.block(0, 0).reinit(n_u, n_u);
      dsp.block(1, 0).reinit(n_p, n_u);
      dsp.block(0, 1).reinit(n_u, n_p);
      dsp.block(1, 1).reinit(n_p, n_p);
      dsp.collect_sizes();
      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
      sparsity_pattern.copy_from(dsp);
    }
    system_matrix.reinit(sparsity_pattern);
    solution.reinit(2);
    solution.block(0).reinit(n_u);
    solution.block(1).reinit(n_p);
    solution.collect_sizes();
    system_rhs.reinit(2);
    system_rhs.block(0).reinit(n_u);
    system_rhs.block(1).reinit(n_p);
    system_rhs.collect_sizes();
  }

  template <int dim>
  void
  StokesProblem<dim>::assemble_system()
  {
    system_matrix = 0;
    system_rhs    = 0;
    QGauss<dim>        quadrature_formula(degree + 2);
    FEValues<dim>      fe_values(fe,
                            quadrature_formula,
                            update_values | update_quadrature_points |
                              update_JxW_values | update_gradients);
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();
    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     local_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    const RightHandSide<dim>             right_hand_side;
    std::vector<Vector<double>> rhs_values(n_q_points, Vector<double>(dim + 1));
    const FEValuesExtractors::Vector     velocities(0);
    const FEValuesExtractors::Scalar     pressure(dim);
    std::vector<SymmetricTensor<2, dim>> symgrad_phi_u(dofs_per_cell);
    std::vector<double>                  div_phi_u(dofs_per_cell);
    std::vector<double>                  phi_p(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);
        local_matrix = 0;
        local_rhs    = 0;
        right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
                                          rhs_values);
        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            for (unsigned int k = 0; k < dofs_per_cell; ++k)
              {
                symgrad_phi_u[k] =
                  fe_values[velocities].symmetric_gradient(k, q);
                div_phi_u[k] = fe_values[velocities].divergence(k, q);
                phi_p[k]     = fe_values[pressure].value(k, q);
              }
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j <= i; ++j)
                  {
                    local_matrix(i, j) +=
                      (2 * (symgrad_phi_u[i] * symgrad_phi_u[j]) -
                       div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j] +
                       phi_p[i] * phi_p[j]) *
                      fe_values.JxW(q);
                  }
                const unsigned int component_i =
                  fe.system_to_component_index(i).first;
                local_rhs(i) += fe_values.shape_value(i, q) *
                                rhs_values[q](component_i) * fe_values.JxW(q);
              }
          }
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = i + 1; j < dofs_per_cell; ++j)
            local_matrix(i, j) = local_matrix(j, i);
        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(local_matrix,
                                               local_rhs,
                                               local_dof_indices,
                                               system_matrix,
                                               system_rhs);
      }
  }

  template <int dim>
  void
  StokesProblem<dim>::solve()
  {
    // Linear operators
    const auto A = linear_operator(system_matrix.block(0, 0));
    const auto B = linear_operator(system_matrix.block(0, 1));
    const auto C = linear_operator(system_matrix.block(1, 0));
    const auto M = linear_operator(
      system_matrix.block(1, 1)); // Mass matrix stored in this block
    const auto D0 = null_operator(M);

    // Inverse of A
    SparseILU<double> preconditioner_A;
    preconditioner_A.initialize(system_matrix.block(0, 0),
                                SparseILU<double>::AdditionalData());
    ReductionControl solver_control_A(system_matrix.block(0, 0).m(),
                                      1e-10,
                                      1e-6);
    SolverCG<>       solver_A(solver_control_A);
    const auto       A_inv = inverse_operator(A, solver_A, preconditioner_A);

    // Inverse of mass matrix stored in block "D"
    SparseILU<double> preconditioner_M;
    preconditioner_M.initialize(system_matrix.block(1, 1),
                                SparseILU<double>::AdditionalData());
    ReductionControl solver_control_M(system_matrix.block(1, 1).m(),
                                      1e-10,
                                      1e-6);
    SolverCG<>       solver_M(solver_control_M);
    const auto       M_inv = inverse_operator(M, solver_M, preconditioner_M);

    // Schur complement
    const auto S = schur_complement(A_inv, B, C, D0);

    // Inverse of Schur complement
    ReductionControl solver_control_S(system_matrix.block(1, 1).m(),
                                      1e-10,
                                      1e-6);
    SolverCG<>       solver_S(solver_control_S);
    const auto       S_inv = inverse_operator(S, solver_S, M_inv);

    Vector<double> &      x   = solution.block(0);
    Vector<double> &      y   = solution.block(1);
    const Vector<double> &f   = system_rhs.block(0);
    const Vector<double> &g   = system_rhs.block(1);
    auto                  rhs = condense_schur_rhs(A_inv, C, f, g);
    y                         = S_inv * rhs;
    x                         = postprocess_schur_solution(A_inv, B, y, f);

    constraints.distribute(solution);
    deallog << "  " << solver_control_S.last_step()
            << " outer CG Schur complement iterations for pressure"
            << std::endl;
  }

  template <int dim>
  void
  StokesProblem<dim>::refine_mesh()
  {
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
    FEValuesExtractors::Scalar pressure(dim);
    KellyErrorEstimator<dim>::estimate(
      dof_handler,
      QGauss<dim - 1>(degree + 1),
      std::map<types::boundary_id, const Function<dim> *>(),
      solution,
      estimated_error_per_cell,
      fe.component_mask(pressure));
    GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                    estimated_error_per_cell,
                                                    0.3,
                                                    0.0);
    triangulation.execute_coarsening_and_refinement();
  }

  template <int dim>
  void
  StokesProblem<dim>::run()
  {
    {
      std::vector<unsigned int> subdivisions(dim, 1);

      const Point<dim> bottom_left =
        (dim == 2 ? Point<dim>(-2, -1) : Point<dim>(-2, 0, -1));
      const Point<dim> top_right =
        (dim == 2 ? Point<dim>(2, 0) : Point<dim>(2, 1, 0));
      GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                subdivisions,
                                                bottom_left,
                                                top_right);
    }
    for (typename Triangulation<dim>::active_cell_iterator cell =
           triangulation.begin_active();
         cell != triangulation.end();
         ++cell)
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        if (cell->face(f)->center()[dim - 1] == 0)
          cell->face(f)->set_all_boundary_ids(1);
    triangulation.refine_global(4 - dim);
    for (unsigned int refinement_cycle = 0; refinement_cycle < 2;
         ++refinement_cycle)
      {
        deallog << "Refinement cycle " << refinement_cycle << std::endl;
        if (refinement_cycle > 0)
          refine_mesh();
        setup_dofs();
        deallog << "   Assembling..." << std::endl;
        assemble_system();
        deallog << "   Solving..." << std::endl;
        solve();
        deallog << std::endl;
      }
  }
} // namespace Step22
int
main()
{
  initlog();
  deallog.depth_file(1);

  using namespace Step22;
  StokesProblem<2> flow_problem(1);
  flow_problem.run();
}

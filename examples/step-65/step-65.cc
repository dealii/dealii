// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

// This tutorial program was contributed by Martin Kronbichler

#include <deal.II/base/timer.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_cache.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>



namespace step65
{
  using namespace dealii;

  template <int dim>
  class ExactSolution : public Function<dim>
  {
  public:
    virtual double value(const Point<dim> &p,
                         const unsigned int /*component*/ = 0) const override
    {
      if (p.norm_square() < 0.25)
        return p.norm_square();
      else
        return 0.1 * p.norm_square() + (0.25 - 0.025);
    }

    virtual Tensor<1, dim>
    gradient(const Point<dim> &p,
             const unsigned int /*component*/ = 0) const override
    {
      if (p.norm_square() < 0.25)
        return 2. * p;
      else
        return 0.2 * p;
    }
  };


  template <int dim>
  class Coefficient : public Function<dim>
  {
  public:
    virtual double value(const Point<dim> &p,
                         const unsigned int /*component*/ = 0) const override
    {
      if (p.norm_square() < 0.25)
        return 0.5;
      else
        return 5.;
    }
  };



  template <int dim>
  class PoissonProblem
  {
  public:
    PoissonProblem();
    void run();

  private:
    void create_grid();
    void setup_system(const Mapping<dim> &mapping);
    void assemble_system(const Mapping<dim> &mapping);
    void solve();
    void postprocess(const Mapping<dim> &mapping);
    void output_results(const Mapping<dim> &mapping);

    Triangulation<dim> triangulation;
    FE_Q<dim>          fe;
    DoFHandler<dim>    dof_handler;

    AffineConstraints<double> constraints;
    SparsityPattern           sparsity_pattern;
    SparseMatrix<double>      system_matrix;
    Vector<double>            solution;
    Vector<double>            system_rhs;

    TimerOutput timer;
  };



  template <int dim>
  PoissonProblem<dim>::PoissonProblem()
    : fe(3)
    , dof_handler(triangulation)
    , timer(std::cout, TimerOutput::never, TimerOutput::wall_times)
  {}



  template <int dim>
  void PoissonProblem<dim>::create_grid()
  {
    Triangulation<dim> tria_outer, tria_inner;
    GridGenerator::hyper_shell(
      tria_outer, Point<dim>(), 0.5, std::sqrt(dim), 2 * dim);
    GridGenerator::hyper_ball(tria_inner, Point<dim>(), 0.5);
    GridGenerator::merge_triangulations(tria_inner, tria_outer, triangulation);
    triangulation.reset_all_manifolds();
    triangulation.set_all_manifold_ids(0);
    for (const auto &cell : triangulation.cell_iterators())
      {
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
          {
            bool face_at_sphere_boundary = true;
            for (unsigned int v = 0;
                 v < GeometryInfo<dim - 1>::vertices_per_cell;
                 ++v)
              if (std::abs(cell->face(f)->vertex(v).norm_square() - 0.25) >
                  1e-12)
                face_at_sphere_boundary = false;
            if (face_at_sphere_boundary)
              cell->face(f)->set_all_manifold_ids(1);
          }
        if (cell->center().norm_square() < 0.25)
          cell->set_material_id(1);
        else
          cell->set_material_id(0);
      }
    triangulation.set_manifold(1, SphericalManifold<dim>());
    TransfiniteInterpolationManifold<dim> transfinite_manifold;
    transfinite_manifold.initialize(triangulation);
    triangulation.set_manifold(0, transfinite_manifold);

    triangulation.refine_global(9 - 2 * dim);
  }



  template <int dim>
  void PoissonProblem<dim>::setup_system(const Mapping<dim> &mapping)
  {
    dof_handler.distribute_dofs(fe);
    std::cout << "   Number of active cells:       "
              << triangulation.n_global_active_cells() << std::endl;
    std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;

    {
      TimerOutput::Scope scope(timer, "Compute constraints");
      constraints.clear();
      DoFTools::make_hanging_node_constraints(dof_handler, constraints);
      VectorTools::interpolate_boundary_values(
        mapping, dof_handler, 0, ExactSolution<dim>(), constraints);
      constraints.close();
    }

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);

    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }



  template <int dim>
  void PoissonProblem<dim>::assemble_system(const Mapping<dim> &mapping)
  {
    TimerOutput::Scope scope(timer, "Assemble linear system");
    QGauss<dim>        quadrature_formula(fe.degree + 1);
    FEValues<dim>      fe_values(mapping,
                            fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);
    FullMatrix<double> partial_matrix(dofs_per_cell, dim * n_q_points);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_rhs = 0.;
        fe_values.reinit(cell);
        for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
          {
            const double current_coefficient =
              Coefficient<dim>().value(fe_values.quadrature_point(q_index));
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int d = 0; d < dim; ++d)
                  partial_matrix(i, q_index * dim + d) =
                    std::sqrt(fe_values.JxW(q_index) * current_coefficient) *
                    fe_values.shape_grad(i, q_index)[d];
                cell_rhs(i) +=
                  (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                   (-dim) *                            // f(x_q)
                   fe_values.JxW(q_index));            // dx
              }
          }
        partial_matrix.mTmult(cell_matrix, partial_matrix);
        cell->get_dof_indices(local_dof_indices);

        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
  }



  template <int dim>
  void PoissonProblem<dim>::solve()
  {
    TimerOutput::Scope   scope(timer, "Solve linear system");
    SolverControl        solver_control(1000, 1e-12);
    SolverCG<>           solver(solver_control);
    PreconditionJacobi<> preconditioner;
    preconditioner.initialize(system_matrix);
    solver.solve(system_matrix, solution, system_rhs, preconditioner);
    constraints.distribute(solution);
    std::cout << "   Number of solver iterations:  "
              << solver_control.last_step() << std::endl;
  }



  template <int dim>
  void PoissonProblem<dim>::postprocess(const Mapping<dim> &mapping)
  {
    {
      TimerOutput::Scope scope(timer, "Write output");
      Timer              time;
      DataOut<dim>       data_out;
      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(solution, "solution");

      DataOutBase::VtkFlags flags;
      flags.write_higher_order_cells = true;
      data_out.set_flags(flags);

      Vector<double> material_ids(triangulation.n_active_cells());
      for (const auto &cell : triangulation.active_cell_iterators())
        material_ids[cell->active_cell_index()] = cell->material_id();
      data_out.add_data_vector(material_ids, "material_ids");

      data_out.build_patches(mapping,
                             fe.degree,
                             DataOut<dim>::curved_inner_cells);

      std::ofstream file(
        ("solution-" +
         std::to_string(triangulation.n_global_levels() - 10 + 2 * dim) +
         ".vtu")
          .c_str());
      data_out.write_vtu(file);
    }

    {
      TimerOutput::Scope scope(timer, "Compute error norms");
      Vector<double>     norm_per_cell_p(triangulation.n_active_cells());

      VectorTools::integrate_difference(mapping,
                                        dof_handler,
                                        solution,
                                        ExactSolution<dim>(),
                                        norm_per_cell_p,
                                        QGauss<dim>(fe.degree + 2),
                                        VectorTools::L2_norm);
      std::cout << "   L2 error vs exact solution:   "
                << norm_per_cell_p.l2_norm() << std::endl;

      VectorTools::integrate_difference(mapping,
                                        dof_handler,
                                        solution,
                                        ExactSolution<dim>(),
                                        norm_per_cell_p,
                                        QGauss<dim>(fe.degree + 2),
                                        VectorTools::H1_norm);
      std::cout << "   H1 error vs exact solution:   "
                << norm_per_cell_p.l2_norm() << std::endl;
    }

    {
      TimerOutput::Scope scope(timer, "Compute error estimator");
      Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
      KellyErrorEstimator<dim>::estimate(
        mapping,
        dof_handler,
        QGauss<dim - 1>(fe.degree + 1),
        std::map<types::boundary_id, const Function<dim> *>(),
        solution,
        estimated_error_per_cell);
    }
  }



  template <int dim>
  void PoissonProblem<dim>::output_results(const Mapping<dim> &mapping)
  {
    TimerOutput::Scope scope(timer, "Write output");
    Timer              time;
    DataOut<dim>       data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");

    DataOutBase::VtkFlags flags;
    flags.write_higher_order_cells = true;
    data_out.set_flags(flags);

    Vector<double> material_ids(triangulation.n_active_cells());
    for (const auto &cell : triangulation.active_cell_iterators())
      material_ids[cell->active_cell_index()] = cell->material_id();
    data_out.add_data_vector(material_ids, "material_ids");

    data_out.build_patches(mapping,
                           fe.degree,
                           DataOut<dim>::curved_inner_cells);

    std::ofstream file(
      ("solution-" +
       std::to_string(triangulation.n_global_levels() - 10 + 2 * dim) + ".vtu")
        .c_str());
    data_out.write_vtu(file);
  }



  template <int dim>
  void PoissonProblem<dim>::run()
  {
    create_grid();
    {
      std::cout << std::endl
                << "====== Running with the basic MappingQGeneric class ====== "
                << std::endl
                << std::endl;

      MappingQGeneric<dim> mapping(fe.degree + 1);
      setup_system(mapping);
      assemble_system(mapping);
      solve();
      postprocess(mapping);
      output_results(mapping);

      timer.print_summary();
    }

    timer.reset();

    {
      std::cout
        << "====== Running with the optimized MappingQCache class ====== "
        << std::endl
        << std::endl;

      MappingQCache<dim> mapping(fe.degree + 1);
      {
        TimerOutput::Scope scope(timer, "Initialize mapping cache");
        mapping.initialize(triangulation, MappingQGeneric<dim>(fe.degree + 1));
      }
      std::cout << "   Memory consumption cache:     "
                << 1e-6 * mapping.memory_consumption() << " MB" << std::endl;

      setup_system(mapping);
      assemble_system(mapping);
      solve();
      postprocess(mapping);
      output_results(mapping);

      timer.print_summary();
    }
  }

} // namespace step65


int main()
{
  step65::PoissonProblem<3> test_program;
  test_program.run();
  return 0;
}

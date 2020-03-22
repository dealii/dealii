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


// A lightly adapted version of the step-4 tutorial program. Compared
// to step-4, this test uses hp::DoFHandler and
// parallel::distributed::Triangulation. In the upper half of the domain
// FENothing is used to save effort, as these elements would only be used later.
// The test checks whether the hanging-node constraints are setup correctly.

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/dof_handler.h>

#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"

namespace Step4
{
  template <int dim>
  class Step4
  {
  public:
    Step4();
    void
    run();

  private:
    void
    make_grid();
    void
    setup_system();
    void
    assemble_system();
    void
    solve();
    void
    output_results() const;

    MPI_Comm communicator;

    parallel::distributed::Triangulation<dim> triangulation;
    hp::FECollection<dim>                     fe;
    hp::DoFHandler<dim>                       dof_handler;
    hp::QCollection<dim>                      q_collection;

    IndexSet locally_owned_dofs;

    IndexSet locally_relevant_dofs;

    TrilinosWrappers::SparseMatrix system_matrix;

    TrilinosWrappers::MPI::Vector solution, rhs;

    AffineConstraints<double> constraints;

    ConditionalOStream pcout;
  };



  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide()
      : Function<dim>()
    {}

    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;
  };



  template <int dim>
  double
  RightHandSide<dim>::value(const Point<dim> &p,
                            const unsigned int /*component*/) const
  {
    double return_value = 0.0;
    for (unsigned int i = 0; i < dim; ++i)
      return_value += 4.0 * std::pow(p(i), 4.0);

    return return_value;
  }



  template <int dim>
  Step4<dim>::Step4()
    : communicator(MPI_COMM_WORLD)
    , triangulation(communicator,
                    typename Triangulation<dim>::MeshSmoothing(
                      Triangulation<dim>::smoothing_on_refinement |
                      Triangulation<dim>::smoothing_on_coarsening))
    , dof_handler(triangulation)
    , pcout(std::cout, (Utilities::MPI::this_mpi_process(communicator) == 0))
  {
    fe.push_back(FE_Q<dim>(2));
    fe.push_back(FE_Nothing<dim>());

    q_collection.push_back(QGauss<dim>(fe[0].degree + 1));
  }



  template <int dim>
  void
  Step4<dim>::make_grid()
  {
    GridGenerator::hyper_cube(triangulation, -1, 1);

    if (dim == 3)
      triangulation.refine_global(3);
    else
      triangulation.refine_global(4);

    for (auto &cell : triangulation.active_cell_iterators())
      if (cell->is_locally_owned())
        if ((cell->center()[dim - 1] < 0.) && (cell->center()[dim - 1] > -0.2))
          cell->set_refine_flag();

    triangulation.execute_coarsening_and_refinement();
  }



  template <int dim>
  void
  Step4<dim>::setup_system()
  {
    for (auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          Point<dim> center = cell->center();

          if (center[dim - 1] >= 0)
            cell->set_active_fe_index(1);

          else
            cell->set_active_fe_index(0);
        }
    dof_handler.distribute_dofs(fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs.clear();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
    locally_relevant_dofs.compress();

    constraints.clear();

    constraints.reinit(locally_relevant_dofs);

    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    DoFTools::make_zero_boundary_constraints(dof_handler, constraints);

#ifdef DEBUG
    // We did not think about hp constraints on ghost cells yet.
    // Thus, we are content with verifying their consistency for now.
    IndexSet locally_active_dofs;
    DoFTools::extract_locally_active_dofs(dof_handler, locally_active_dofs);
    AssertThrow(constraints.is_consistent_in_parallel(
                  dof_handler.locally_owned_dofs_per_processor(),
                  locally_active_dofs,
                  communicator,
                  /*verbose=*/true),
                ExcMessage(
                  "AffineConstraints object contains inconsistencies!"));
#endif

    constraints.close();

    DynamicSparsityPattern dsp(locally_relevant_dofs);

    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);

    SparsityTools::distribute_sparsity_pattern(dsp,
                                               locally_owned_dofs,
                                               communicator,
                                               locally_relevant_dofs);

    system_matrix.reinit(locally_owned_dofs,
                         locally_owned_dofs,
                         dsp,
                         communicator);

    solution.reinit(locally_relevant_dofs, communicator);

    rhs.reinit(locally_owned_dofs, communicator);
  }



  template <int dim>
  void
  Step4<dim>::assemble_system()
  {
    const RightHandSide<dim> right_hand_side;

    system_matrix = 0;
    rhs           = 0;

    hp::FEValues<dim> hp_fe_v(dof_handler.get_fe_collection(),
                              q_collection,
                              update_JxW_values | update_gradients |
                                update_values | update_quadrature_points);

    std::vector<types::global_dof_index> local_dof_indices;

    FullMatrix<double> cell_matrix;
    Vector<double>     cell_rhs;

    for (auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned() && cell->active_fe_index() == 0)
          {
            const unsigned int active_fe_index = cell->active_fe_index();

            if (active_fe_index == 0)
              {
                hp_fe_v.reinit(cell);

                const FEValues<dim> &fe_values =
                  hp_fe_v.get_present_fe_values();

                cell_matrix.reinit(fe[active_fe_index].dofs_per_cell,
                                   fe[active_fe_index].dofs_per_cell);

                cell_rhs.reinit(fe[active_fe_index].dofs_per_cell);

                local_dof_indices.resize(fe[active_fe_index].dofs_per_cell);

                cell->get_dof_indices(local_dof_indices);

                const unsigned int dofs_per_cell =
                  fe[active_fe_index].dofs_per_cell;

                const unsigned int n_q_points =
                  q_collection[active_fe_index].size();

                for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    {
                      for (unsigned int j = 0; j < dofs_per_cell; ++j)
                        cell_matrix(i, j) += (fe_values.shape_grad(i, q_index) *
                                              fe_values.shape_grad(j, q_index) *
                                              fe_values.JxW(q_index));

                      const auto x_q = fe_values.quadrature_point(q_index);
                      cell_rhs(i) +=
                        (fe_values.shape_value(i, q_index) *
                         right_hand_side.value(x_q) * fe_values.JxW(q_index));
                    }

                constraints.distribute_local_to_global(
                  cell_matrix, cell_rhs, local_dof_indices, system_matrix, rhs);
              }
          }
      }
    system_matrix.compress(VectorOperation::add);

    rhs.compress(VectorOperation::add);
  }



  template <int dim>
  void
  Step4<dim>::solve()
  {
    SolverControl solver_control(system_matrix.m(), 1e-12 * rhs.l2_norm());

    TrilinosWrappers::SolverDirect::AdditionalData add_data(false,
                                                            "Amesos_Mumps");

    TrilinosWrappers::SolverDirect direct_solver(solver_control, add_data);

    TrilinosWrappers::MPI::Vector completely_distributed_solution(
      locally_owned_dofs, communicator);

    direct_solver.solve(system_matrix, completely_distributed_solution, rhs);

    constraints.distribute(completely_distributed_solution);

    solution = completely_distributed_solution;
  }



  template <int dim>
  void
  Step4<dim>::run()
  {
    pcout << "Solving problem in " << dim << " space dimensions." << std::endl;

    make_grid();

    setup_system();

    pcout << "   Number of active cells:       "
          << triangulation.n_global_active_cells() << std::endl;

    pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl;

    assemble_system();
    solve();
  }
} // namespace Step4



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  Step4::Step4<2> laplace_problem_2d;

  laplace_problem_2d.run();

  Step4::Step4<3> laplace_problem_3d;

  laplace_problem_3d.run();

  return 0;
}

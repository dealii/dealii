// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Solve iteratively the domain-decomposition problem.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_evaluate.h>

#include "../tests.h"


template <typename MeshType>
MPI_Comm
get_mpi_comm(const MeshType &mesh)
{
  const auto *tria_parallel = dynamic_cast<
    const parallel::TriangulationBase<MeshType::dimension,
                                      MeshType::space_dimension> *>(
    &(mesh.get_triangulation()));

  return tria_parallel != nullptr ? tria_parallel->get_mpi_communicator() :
                                    MPI_COMM_SELF;
}

template <int dim, int spacedim>
std::shared_ptr<const Utilities::MPI::Partitioner>
create_partitioner(const DoFHandler<dim, spacedim> &dof_handler)
{
  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  return std::make_shared<const Utilities::MPI::Partitioner>(
    dof_handler.locally_owned_dofs(),
    locally_relevant_dofs,
    get_mpi_comm(dof_handler));
}


template <int dim>
void
print(const Mapping<dim>                               &mapping,
      const DoFHandler<dim>                            &dof_handler,
      const LinearAlgebra::distributed::Vector<double> &result,
      const unsigned int                                counter)
{
  result.update_ghost_values();
  DataOutBase::VtkFlags flags;
  if (counter == 0)
    flags.write_higher_order_cells = true;

  DataOut<dim> data_out;
  data_out.set_flags(flags);
  data_out.attach_dof_handler(dof_handler);

  const auto &tria = dof_handler.get_triangulation();

  Vector<double> ranks(tria.n_active_cells());
  for (const auto &cell :
       tria.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
    ranks(cell->active_cell_index()) = cell->subdomain_id();
  data_out.add_data_vector(ranks, "rank");
  data_out.add_data_vector(result, "result");

  data_out.build_patches(
    mapping,
    std::min<unsigned int>(2, dof_handler.get_fe().tensor_degree() + 1));
  data_out.write_vtu_with_pvtu_record(
    "./", "example-2", counter, MPI_COMM_WORLD, 1, 1);

  result.zero_out_ghost_values();
}

template <int dim>
class PoissonProblem
{
public:
  PoissonProblem(const Triangulation<dim> &tria,
                 const unsigned int        id,
                 const Mapping<dim>       &mapping,
                 const FiniteElement<dim> &fe,
                 const Quadrature<dim>    &quad)
    : id(id)
    , mapping(mapping)
    , fe(fe)
    , quad(quad)
    , dof_handler(tria)
  {
    dof_handler.distribute_dofs(fe);
    partitioner = create_partitioner(dof_handler);
    solution.reinit(partitioner);

    print(mapping, dof_handler, solution, id);
  }

  void
  solve(const DoFHandler<dim>                            &dof_handler_other,
        const LinearAlgebra::distributed::Vector<double> &solution_other,
        const Mapping<dim>                               &mapping_other)
  {
    AffineConstraints<double> constraints(
      dof_handler.locally_owned_dofs(),
      DoFTools::extract_locally_relevant_dofs(dof_handler));
    {
      std::map<types::global_dof_index, Point<dim>> support_points;

      DoFTools::map_dofs_to_support_points(mapping,
                                           dof_handler,
                                           support_points);

      std::vector<types::global_dof_index> global_ids;

      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (cell->is_locally_owned() == false)
            continue;

          for (const auto &face : cell->face_iterators())
            {
              if (face->at_boundary() == false)
                continue;

              std::vector<types::global_dof_index> dofs_on_face(
                dof_handler.get_fe().n_dofs_per_face());
              face->get_dof_indices(dofs_on_face);

              global_ids.insert(global_ids.end(),
                                dofs_on_face.begin(),
                                dofs_on_face.end());
            }
        }

      std::sort(global_ids.begin(), global_ids.end());
      global_ids.erase(unique(global_ids.begin(), global_ids.end()),
                       global_ids.end());

      std::vector<Point<dim>> evaluation_points;
      evaluation_points.reserve(global_ids.size());

      for (const auto &global_id : global_ids)
        evaluation_points.push_back(support_points[global_id]);

      solution_other.update_ghost_values();
      Utilities::MPI::RemotePointEvaluation<dim> evaluation_cache;
      const auto                                 evaluation_point_results =
        VectorTools::point_values<1>(mapping_other,
                                     dof_handler_other,
                                     solution_other,
                                     evaluation_points,
                                     evaluation_cache);
      solution_other.zero_out_ghost_values();
      for (unsigned int i = 0; i < evaluation_points.size(); ++i)
        {
          if (global_ids[i] == numbers::invalid_size_type)
            continue;

          constraints.add_constraint(global_ids[i],
                                     {},
                                     evaluation_point_results[i]);
        }
    }
    constraints.close();

    {
      // initialize vectors and system matrix
      LinearAlgebra::distributed::Vector<double> b;
      b.reinit(partitioner);

      TrilinosWrappers::SparseMatrix A;

      TrilinosWrappers::SparsityPattern dsp(dof_handler.locally_owned_dofs(),
                                            get_mpi_comm(dof_handler));
      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
      dsp.compress();
      A.reinit(dsp);

      // assemble right-hand side and system matrix
      FEValues<dim> fe_values(mapping,
                              fe,
                              quad,
                              update_values | update_gradients |
                                update_JxW_values);

      FullMatrix<double>                   cell_matrix;
      Vector<double>                       cell_rhs;
      std::vector<types::global_dof_index> local_dof_indices;

      // loop over all cells
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (cell->is_locally_owned() == false)
            continue;

          fe_values.reinit(cell);

          const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
          cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
          cell_rhs.reinit(dofs_per_cell);

          // loop over cell dofs
          for (const auto q : fe_values.quadrature_point_indices())
            {
              for (const auto i : fe_values.dof_indices())
                for (const auto j : fe_values.dof_indices())
                  cell_matrix(i, j) +=
                    (fe_values.shape_grad(i, q) * // grad phi_i(x_q)
                     fe_values.shape_grad(j, q) * // grad phi_j(x_q)
                     fe_values.JxW(q));           // dx

              for (const unsigned int i : fe_values.dof_indices())
                cell_rhs(i) += (fe_values.shape_value(i, q) * // phi_i(x_q)
                                1. *                          // f(x_q)
                                fe_values.JxW(q));            // dx
            }

          local_dof_indices.resize(cell->get_fe().dofs_per_cell);
          cell->get_dof_indices(local_dof_indices);

          constraints.distribute_local_to_global(
            cell_matrix, cell_rhs, local_dof_indices, A, b);
        }

      b.compress(VectorOperation::values::add);
      A.compress(VectorOperation::values::add);

      // solve linear equation system
      ReductionControl reduction_control(1000, 1e-12, 1e-10);
      SolverCG<LinearAlgebra::distributed::Vector<double>> solver(
        reduction_control);
      solver.solve(A, solution, b, PreconditionIdentity());

      if (false && Utilities::MPI::this_mpi_process(
                     get_mpi_comm(dof_handler.get_triangulation())) == 0)
        deallog << "Solved in " << reduction_control.last_step()
                << " iterations." << std::endl;

      constraints.distribute(solution);
    }

    print(mapping, dof_handler, solution, id);
  }

public:
  const unsigned int                                 id;
  const Mapping<dim>                                &mapping;
  const FiniteElement<dim>                          &fe;
  const Quadrature<dim>                             &quad;
  DoFHandler<dim>                                    dof_handler;
  std::shared_ptr<const Utilities::MPI::Partitioner> partitioner;
  LinearAlgebra::distributed::Vector<double>         solution;
};

void
test()
{
  const unsigned int dim             = 2;
  const unsigned int n_refinements_1 = 4;
  const unsigned int n_refinements_2 = 4;

  const MappingQ<dim> mapping_1(1);
  const FE_Q<dim>     fe_1(2);
  const QGauss<dim>   quad_1(3);

#if false
  const MappingFE<dim>     mapping_2(Simplex::FE_P<dim>(1));
  const Simplex::FE_P<dim> fe_2(2);
  const Simplex::QGauss<dim>        quad_2(3);
#else
  const MappingQ<dim> mapping_2(1);
  const FE_Q<dim>     fe_2(2);
  const QGauss<dim>   quad_2(3);
#endif

  parallel::distributed::Triangulation<dim> tria_1(MPI_COMM_WORLD);
  GridGenerator::hyper_ball(tria_1);
  tria_1.refine_global(n_refinements_1);

  parallel::shared::Triangulation<dim> tria_2(MPI_COMM_WORLD);
#if false
  GridIn<dim> grid_in;
  grid_in.attach_triangulation(tria_2);
  std::ifstream input_file("bar.msh");
  grid_in.read_msh(input_file);
#else
  GridGenerator::hyper_rectangle(tria_2, {0.0, -0.5}, {2.0, 0.5});
  tria_2.refine_global(n_refinements_2);
#endif

  PoissonProblem<dim> pp_1(tria_1, 0, mapping_1, fe_1, quad_1);
  PoissonProblem<dim> pp_2(tria_2, 1, mapping_2, fe_2, quad_2);

  for (unsigned int i = 0; i < 10; ++i)
    {
      pp_1.solve(pp_2.dof_handler, pp_2.solution, mapping_2);
      pp_2.solve(pp_1.dof_handler, pp_1.solution, mapping_1);

      // break;
    }

  pp_1.solution.print(deallog.get_file_stream());
  pp_2.solution.print(deallog.get_file_stream());
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.depth_file(1);

  test();
}

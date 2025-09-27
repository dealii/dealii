// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Same as assemble_matrix_parallel_02, but using block matrices

#include <deal.II/base/function.h>
#include <deal.II/base/graph_coloring.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <complex>
#include <iostream>

#include "../tests.h"

namespace Assembly
{
  namespace Scratch
  {
    template <int dim>
    struct Data
    {
      Data(const FiniteElement<dim> &fe, const Quadrature<dim> &quadrature)
        : fe_values(fe,
                    quadrature,
                    update_values | update_gradients |
                      update_quadrature_points | update_JxW_values)
      {}

      Data(const Data &data)
        : fe_values(data.fe_values.get_mapping(),
                    data.fe_values.get_fe(),
                    data.fe_values.get_quadrature(),
                    data.fe_values.get_update_flags())
      {}

      FEValues<dim> fe_values;
    };
  } // namespace Scratch

  namespace Copy
  {
    struct Data
    {
      Data(const bool assemble_reference)
        : assemble_reference(assemble_reference)
      {}
      std::vector<types::global_dof_index> local_dof_indices;
      FullMatrix<double>                   local_matrix;
      Vector<double>                       local_rhs;
      const bool                           assemble_reference;
    };
  } // namespace Copy
} // namespace Assembly

template <int dim>
class LaplaceProblem
{
public:
  LaplaceProblem();
  ~LaplaceProblem();

  void
  run();

private:
  void
  setup_system();
  void
  test_equality();
  void
  assemble_reference();
  void
  assemble_test();
  void
  solve();
  void
  create_coarse_grid();
  void
  postprocess();

  void
  local_assemble(
    const FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
                                 &cell,
    Assembly::Scratch::Data<dim> &scratch,
    Assembly::Copy::Data         &data);
  void
  copy_local_to_global(const Assembly::Copy::Data &data);

  std::vector<types::global_dof_index>
  get_conflict_indices(
    const FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
      &cell) const;

  parallel::distributed::Triangulation<dim> triangulation;

  DoFHandler<dim> dof_handler;
  FESystem<dim>   fe;
  QGauss<dim>     quadrature;

  AffineConstraints<double> constraints;

  TrilinosWrappers::BlockSparseMatrix reference_matrix;
  TrilinosWrappers::BlockSparseMatrix test_matrix;

  TrilinosWrappers::MPI::BlockVector reference_rhs;
  TrilinosWrappers::MPI::BlockVector test_rhs;

  std::vector<std::vector<
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>>>
    graph;
};



template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  BoundaryValues()
    : Function<dim>(2)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component) const;
};


template <int dim>
double
BoundaryValues<dim>::value(const Point<dim> &p,
                           const unsigned int /*component*/) const
{
  double sum = 0;
  for (unsigned int d = 0; d < dim; ++d)
    sum += std::sin(numbers::PI * p[d]);
  return sum;
}


template <int dim>
class RightHandSide : public Function<dim>
{
public:
  RightHandSide()
    : Function<dim>()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component) const;
};


template <int dim>
double
RightHandSide<dim>::value(const Point<dim> &p,
                          const unsigned int /*component*/) const
{
  double product = 1;
  for (unsigned int d = 0; d < dim; ++d)
    product *= (p[d] + 1);
  return product;
}


template <int dim>
LaplaceProblem<dim>::LaplaceProblem()
  : triangulation(MPI_COMM_WORLD)
  , dof_handler(triangulation)
  , fe(FE_Q<dim>(1), 1, FE_Q<dim>(2), 1)
  , quadrature(3)
{}


template <int dim>
LaplaceProblem<dim>::~LaplaceProblem()
{
  dof_handler.clear();
}



template <int dim>
std::vector<types::global_dof_index>
LaplaceProblem<dim>::get_conflict_indices(
  const FilteredIterator<typename DoFHandler<dim>::active_cell_iterator> &cell)
  const
{
  std::vector<types::global_dof_index> local_dof_indices(
    cell->get_fe().dofs_per_cell);
  cell->get_dof_indices(local_dof_indices);

  constraints.resolve_indices(local_dof_indices);
  return local_dof_indices;
}

template <int dim>
void
LaplaceProblem<dim>::setup_system()
{
  reference_matrix.clear();
  test_matrix.clear();
  dof_handler.distribute_dofs(fe);
  std::vector<unsigned int> blocks(2, 0);
  blocks[1] = 1;
  DoFRenumbering::component_wise(dof_handler, blocks);

  constraints.clear();

  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  // add boundary conditions as inhomogeneous constraints here, do it after
  // having added the hanging node constraints in order to be consistent and
  // skip dofs that are already constrained (i.e., are hanging nodes on the
  // boundary in 3D). In contrast to step-27, we choose a sine function.
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           BoundaryValues<dim>(),
                                           constraints);
  constraints.close();

  using CellFilter =
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>;
  CellFilter begin(IteratorFilters::LocallyOwnedCell(),
                   dof_handler.begin_active());
  CellFilter end(IteratorFilters::LocallyOwnedCell(), dof_handler.end());
  graph = GraphColoring::make_graph_coloring(
    begin,
    end,
    static_cast<std::function<std::vector<types::global_dof_index>(
      const FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
        &)>>(std::bind(&LaplaceProblem<dim>::get_conflict_indices,
                       this,
                       std::placeholders::_1)));

  TrilinosWrappers::BlockSparsityPattern csp(2, 2);
  const IndexSet &locally_owned_total = dof_handler.locally_owned_dofs();
  const IndexSet  relevant_total =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  const std::vector<types::global_dof_index> dofs_per_block =
    DoFTools::count_dofs_per_fe_block(dof_handler, blocks);

  std::vector<IndexSet> locally_owned(2), relevant_set(2);
  locally_owned[0] = locally_owned_total.get_view(0, dofs_per_block[0]);
  locally_owned[1] =
    locally_owned_total.get_view(dofs_per_block[0], dof_handler.n_dofs());
  relevant_set[0] = relevant_total.get_view(0, dofs_per_block[0]);
  relevant_set[1] =
    relevant_total.get_view(dofs_per_block[0], dof_handler.n_dofs());

  {
    csp.reinit(locally_owned, locally_owned, MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, csp, constraints, false);
    csp.compress();
    reference_matrix.reinit(csp);
    reference_rhs.reinit(locally_owned, MPI_COMM_WORLD);
  }
  {
    csp.reinit(locally_owned, locally_owned, relevant_set, MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, csp, constraints, false);
    csp.compress();
    test_matrix.reinit(csp);
    test_rhs.reinit(locally_owned, relevant_set, MPI_COMM_WORLD, true);
  }
}



template <int dim>
void
LaplaceProblem<dim>::local_assemble(
  const FilteredIterator<typename DoFHandler<dim>::active_cell_iterator> &cell,
  Assembly::Scratch::Data<dim> &scratch,
  Assembly::Copy::Data         &data)
{
  const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

  data.local_matrix.reinit(dofs_per_cell, dofs_per_cell);
  data.local_matrix = 0;

  data.local_rhs.reinit(dofs_per_cell);
  data.local_rhs = 0;

  scratch.fe_values.reinit(cell);

  const FEValues<dim> &fe_values = scratch.fe_values;

  const RightHandSide<dim> rhs_function;

  // this does not make a lot of sense physically but it serves the purpose of
  // the test well
  for (unsigned int q_point = 0; q_point < fe_values.n_quadrature_points;
       ++q_point)
    {
      const double rhs_value =
        rhs_function.value(fe_values.quadrature_point(q_point), 0);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            data.local_matrix(i, j) +=
              (fe_values.shape_grad(i, q_point) *
               fe_values.shape_grad(j, q_point) * fe_values.JxW(q_point));

          data.local_rhs(i) += (fe_values.shape_value(i, q_point) * rhs_value *
                                fe_values.JxW(q_point));
        }
    }

  data.local_dof_indices.resize(dofs_per_cell);
  cell->get_dof_indices(data.local_dof_indices);
}



template <int dim>
void
LaplaceProblem<dim>::copy_local_to_global(const Assembly::Copy::Data &data)
{
  if (data.assemble_reference)
    constraints.distribute_local_to_global(data.local_matrix,
                                           data.local_rhs,
                                           data.local_dof_indices,
                                           reference_matrix,
                                           reference_rhs);
  else
    constraints.distribute_local_to_global(data.local_matrix,
                                           data.local_rhs,
                                           data.local_dof_indices,
                                           test_matrix,
                                           test_rhs);
}



template <int dim>
void
LaplaceProblem<dim>::assemble_reference()
{
  reference_matrix = 0;
  reference_rhs    = 0;

  Assembly::Copy::Data         copy_data(true);
  Assembly::Scratch::Data<dim> assembly_data(fe, quadrature);

  for (unsigned int color = 0; color < graph.size(); ++color)
    for (typename std::vector<FilteredIterator<
           typename DoFHandler<dim>::active_cell_iterator>>::const_iterator p =
           graph[color].begin();
         p != graph[color].end();
         ++p)
      {
        local_assemble(*p, assembly_data, copy_data);
        copy_local_to_global(copy_data);
      }
  reference_matrix.compress(VectorOperation::add);
  reference_rhs.compress(VectorOperation::add);
}



template <int dim>
void
LaplaceProblem<dim>::assemble_test()
{
  test_matrix = 0;
  test_rhs    = 0;

  WorkStream::run(graph,
                  std::bind(&LaplaceProblem<dim>::local_assemble,
                            this,
                            std::placeholders::_1,
                            std::placeholders::_2,
                            std::placeholders::_3),
                  std::bind(&LaplaceProblem<dim>::copy_local_to_global,
                            this,
                            std::placeholders::_1),
                  Assembly::Scratch::Data<dim>(fe, quadrature),
                  Assembly::Copy::Data(false),
                  2 * MultithreadInfo::n_threads(),
                  1);
  test_matrix.compress(VectorOperation::add);
  test_rhs.compress(VectorOperation::add);

  test_matrix.add(-1, reference_matrix);

  // there should not even be roundoff difference between matrices
  double frobenius_norm = 0;
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < 2; ++j)
      frobenius_norm += numbers::NumberTraits<double>::abs_square(
        test_matrix.block(i, j).frobenius_norm());
  deallog << "error in matrix: " << std::sqrt(frobenius_norm) << std::endl;
  test_rhs.add(-1., reference_rhs);
  deallog << "error in vector: " << test_rhs.l2_norm() << std::endl;
}



template <int dim>
void
LaplaceProblem<dim>::postprocess()
{
  Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
  for (unsigned int i = 0; i < estimated_error_per_cell.size(); ++i)
    estimated_error_per_cell(i) = i;

  GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                  estimated_error_per_cell,
                                                  0.3,
                                                  0.03);
  triangulation.execute_coarsening_and_refinement();
}



template <int dim>
void
LaplaceProblem<dim>::run()
{
  for (unsigned int cycle = 0; cycle < 3; ++cycle)
    {
      if (cycle == 0)
        {
          GridGenerator::hyper_shell(
            triangulation, Point<dim>(), 0.5, 1., (dim == 3) ? 96 : 12, false);
          if constexpr (running_in_debug_mode())
            {
              triangulation.refine_global(3);
            }
          else
            {
              triangulation.refine_global(5);
            }
        }

      setup_system();

      assemble_reference();
      assemble_test();

      if (cycle < 2)
        postprocess();
    }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
  mpi_initlog();
  deallog << std::setprecision(2);

  {
    deallog.push("2d");
    LaplaceProblem<2> laplace_problem;
    laplace_problem.run();
    deallog.pop();
  }
}

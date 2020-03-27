// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2018 by the deal.II authors
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



// similar as parallel_multigrid but using adaptive meshes with hanging nodes
// (doing local smoothing)

#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim,
          int fe_degree,
          int n_q_points_1d = fe_degree + 1,
          typename number   = double>
class LaplaceOperator : public Subscriptor
{
public:
  LaplaceOperator(){};


  void
  initialize(const Mapping<dim> &     mapping,
             const DoFHandler<dim> &  dof_handler,
             const MGConstrainedDoFs &mg_constrained_dofs,
             const std::map<types::boundary_id, const Function<dim> *>
               &                dirichlet_boundary,
             const unsigned int level = numbers::invalid_unsigned_int)
  {
    const QGauss<1>                                  quad(n_q_points_1d);
    typename MatrixFree<dim, number>::AdditionalData addit_data;
    addit_data.tasks_parallel_scheme =
      MatrixFree<dim, number>::AdditionalData::none;
    addit_data.tasks_block_size = 3;
    addit_data.mg_level         = level;
    AffineConstraints<double> constraints;
    if (level == numbers::invalid_unsigned_int)
      {
        IndexSet relevant_dofs;
        DoFTools::extract_locally_relevant_dofs(dof_handler, relevant_dofs);
        constraints.reinit(relevant_dofs);
        DoFTools::make_hanging_node_constraints(dof_handler, constraints);
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 dirichlet_boundary,
                                                 constraints);
      }
    else
      {
        IndexSet relevant_dofs;
        DoFTools::extract_locally_relevant_level_dofs(dof_handler,
                                                      level,
                                                      relevant_dofs);
        constraints.reinit(relevant_dofs);
        constraints.add_lines(mg_constrained_dofs.get_boundary_indices(level));

        std::vector<types::global_dof_index> interface_indices;
        mg_constrained_dofs.get_refinement_edge_indices(level)
          .fill_index_vector(interface_indices);
        edge_constrained_indices.clear();
        edge_constrained_indices.reserve(interface_indices.size());
        edge_constrained_values.resize(interface_indices.size());
        const IndexSet &locally_owned =
          dof_handler.locally_owned_mg_dofs(level);
        for (unsigned int i = 0; i < interface_indices.size(); ++i)
          if (locally_owned.is_element(interface_indices[i]))
            edge_constrained_indices.push_back(
              locally_owned.index_within_set(interface_indices[i]));
        have_interface_matrices =
          Utilities::MPI::max((unsigned int)edge_constrained_indices.size(),
                              MPI_COMM_WORLD) > 0;
      }
    constraints.close();

    data.reinit(mapping, dof_handler, constraints, quad, addit_data);

    if (level != numbers::invalid_unsigned_int)
      compute_inverse_diagonal();
  }

  void
  vmult(LinearAlgebra::distributed::Vector<number> &      dst,
        const LinearAlgebra::distributed::Vector<number> &src) const
  {
    dst = 0;
    vmult_add(dst, src);
  }

  void
  Tvmult(LinearAlgebra::distributed::Vector<number> &      dst,
         const LinearAlgebra::distributed::Vector<number> &src) const
  {
    dst = 0;
    vmult_add(dst, src);
  }

  void
  Tvmult_add(LinearAlgebra::distributed::Vector<number> &      dst,
             const LinearAlgebra::distributed::Vector<number> &src) const
  {
    vmult_add(dst, src);
  }

  void
  vmult_add(LinearAlgebra::distributed::Vector<number> &      dst,
            const LinearAlgebra::distributed::Vector<number> &src) const
  {
    Assert(src.partitioners_are_globally_compatible(
             *data.get_dof_info(0).vector_partitioner),
           ExcInternalError());
    Assert(dst.partitioners_are_globally_compatible(
             *data.get_dof_info(0).vector_partitioner),
           ExcInternalError());

    // set zero Dirichlet values on the input vector (and remember the src and
    // dst values because we need to reset them at the end)
    for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
      {
        edge_constrained_values[i] = std::pair<number, number>(
          src.local_element(edge_constrained_indices[i]),
          dst.local_element(edge_constrained_indices[i]));
        const_cast<LinearAlgebra::distributed::Vector<double> &>(src)
          .local_element(edge_constrained_indices[i]) = 0.;
      }

    data.cell_loop(&LaplaceOperator::local_apply, this, dst, src);

    const std::vector<unsigned int> &constrained_dofs =
      data.get_constrained_dofs();
    for (unsigned int i = 0; i < constrained_dofs.size(); ++i)
      dst.local_element(constrained_dofs[i]) +=
        src.local_element(constrained_dofs[i]);

    // reset edge constrained values, multiply by unit matrix and add into
    // destination
    for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
      {
        const_cast<LinearAlgebra::distributed::Vector<double> &>(src)
          .local_element(edge_constrained_indices[i]) =
          edge_constrained_values[i].first;
        dst.local_element(edge_constrained_indices[i]) =
          edge_constrained_values[i].second + edge_constrained_values[i].first;
      }
  }

  void
  vmult_interface_down(
    LinearAlgebra::distributed::Vector<double> &      dst,
    const LinearAlgebra::distributed::Vector<double> &src) const
  {
    Assert(src.partitioners_are_globally_compatible(
             *data.get_dof_info(0).vector_partitioner),
           ExcInternalError());
    Assert(dst.partitioners_are_globally_compatible(
             *data.get_dof_info(0).vector_partitioner),
           ExcInternalError());

    dst = 0;

    if (!have_interface_matrices)
      return;

    // set zero Dirichlet values on the input vector (and remember the src and
    // dst values because we need to reset them at the end)
    for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
      {
        const double src_val = src.local_element(edge_constrained_indices[i]);
        const_cast<LinearAlgebra::distributed::Vector<double> &>(src)
          .local_element(edge_constrained_indices[i]) = 0.;
        edge_constrained_values[i] = std::pair<number, number>(
          src_val, dst.local_element(edge_constrained_indices[i]));
      }

    data.cell_loop(&LaplaceOperator::local_apply, this, dst, src);

    unsigned int c = 0;
    for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
      {
        for (; c < edge_constrained_indices[i]; ++c)
          dst.local_element(c) = 0.;
        ++c;

        // reset the src values
        const_cast<LinearAlgebra::distributed::Vector<double> &>(src)
          .local_element(edge_constrained_indices[i]) =
          edge_constrained_values[i].first;
      }
    for (; c < dst.local_size(); ++c)
      dst.local_element(c) = 0.;
  }

  void
  vmult_interface_up(
    LinearAlgebra::distributed::Vector<double> &      dst,
    const LinearAlgebra::distributed::Vector<double> &src) const
  {
    Assert(src.partitioners_are_globally_compatible(
             *data.get_dof_info(0).vector_partitioner),
           ExcInternalError());
    Assert(dst.partitioners_are_globally_compatible(
             *data.get_dof_info(0).vector_partitioner),
           ExcInternalError());

    dst = 0;

    if (!have_interface_matrices)
      return;

    LinearAlgebra::distributed::Vector<double> src_cpy(src);
    unsigned int                               c = 0;
    for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
      {
        for (; c < edge_constrained_indices[i]; ++c)
          src_cpy.local_element(c) = 0.;
        ++c;
      }
    for (; c < src_cpy.local_size(); ++c)
      src_cpy.local_element(c) = 0.;

    data.cell_loop(&LaplaceOperator::local_apply, this, dst, src_cpy);
    for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
      {
        dst.local_element(edge_constrained_indices[i]) = 0.;
      }
  }

  types::global_dof_index
  m() const
  {
    return data.get_vector_partitioner()->size();
  }

  types::global_dof_index
  n() const
  {
    return data.get_vector_partitioner()->size();
  }

  number
  el(const unsigned int row, const unsigned int col) const
  {
    AssertThrow(false,
                ExcMessage("Matrix-free does not allow for entry access"));
    return number();
  }

  void
  initialize_dof_vector(
    LinearAlgebra::distributed::Vector<number> &vector) const
  {
    if (!vector.partitioners_are_compatible(
          *data.get_dof_info(0).vector_partitioner))
      data.initialize_dof_vector(vector);
    Assert(vector.partitioners_are_globally_compatible(
             *data.get_dof_info(0).vector_partitioner),
           ExcInternalError());
  }

  const LinearAlgebra::distributed::Vector<number> &
  get_matrix_diagonal_inverse() const
  {
    Assert(inverse_diagonal_entries.size() > 0, ExcNotInitialized());
    return inverse_diagonal_entries;
  }


private:
  void
  local_apply(const MatrixFree<dim, number> &                   data,
              LinearAlgebra::distributed::Vector<number> &      dst,
              const LinearAlgebra::distributed::Vector<number> &src,
              const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, n_q_points_1d, 1, number> phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values(src);
        phi.evaluate(false, true, false);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_gradient(phi.get_gradient(q), q);
        phi.integrate(false, true);
        phi.distribute_local_to_global(dst);
      }
  }

  void
  compute_inverse_diagonal()
  {
    data.initialize_dof_vector(inverse_diagonal_entries);
    unsigned int dummy;
    data.cell_loop(&LaplaceOperator::local_diagonal_cell,
                   this,
                   inverse_diagonal_entries,
                   dummy);

    const std::vector<unsigned int> &constrained_dofs =
      data.get_constrained_dofs();
    for (unsigned int i = 0; i < constrained_dofs.size(); ++i)
      inverse_diagonal_entries.local_element(constrained_dofs[i]) = 1.;
    for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
      {
        inverse_diagonal_entries.local_element(edge_constrained_indices[i]) =
          1.;
      }


    for (unsigned int i = 0; i < inverse_diagonal_entries.local_size(); ++i)
      if (std::abs(inverse_diagonal_entries.local_element(i)) > 1e-10)
        inverse_diagonal_entries.local_element(i) =
          1. / inverse_diagonal_entries.local_element(i);
      else
        inverse_diagonal_entries.local_element(i) = 1.;
  }

  void
  local_diagonal_cell(
    const MatrixFree<dim, number> &             data,
    LinearAlgebra::distributed::Vector<number> &dst,
    const unsigned int &,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, n_q_points_1d, 1, number> phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);

        VectorizedArray<number> local_diagonal_vector[phi.tensor_dofs_per_cell];
        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < phi.dofs_per_cell; ++j)
              phi.begin_dof_values()[j] = VectorizedArray<number>();
            phi.begin_dof_values()[i] = 1.;
            phi.evaluate(false, true, false);
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              phi.submit_gradient(phi.get_gradient(q), q);
            phi.integrate(false, true);
            local_diagonal_vector[i] = phi.begin_dof_values()[i];
          }
        for (unsigned int i = 0; i < phi.tensor_dofs_per_cell; ++i)
          phi.begin_dof_values()[i] = local_diagonal_vector[i];
        phi.distribute_local_to_global(dst);
      }
  }

  MatrixFree<dim, number>                        data;
  LinearAlgebra::distributed::Vector<number>     inverse_diagonal_entries;
  std::vector<unsigned int>                      edge_constrained_indices;
  mutable std::vector<std::pair<number, number>> edge_constrained_values;
  bool                                           have_interface_matrices;
};



template <typename LAPLACEOPERATOR>
class MGInterfaceMatrix : public Subscriptor
{
public:
  void
  initialize(const LAPLACEOPERATOR &laplace)
  {
    this->laplace = &laplace;
  }

  void
  vmult(LinearAlgebra::distributed::Vector<double> &      dst,
        const LinearAlgebra::distributed::Vector<double> &src) const
  {
    laplace->vmult_interface_down(dst, src);
  }

  void
  Tvmult(LinearAlgebra::distributed::Vector<double> &      dst,
         const LinearAlgebra::distributed::Vector<double> &src) const
  {
    laplace->vmult_interface_up(dst, src);
  }

private:
  SmartPointer<const LAPLACEOPERATOR> laplace;
};



template <typename LAPLACEOPERATOR>
class MGTransferMF
  : public MGTransferPrebuilt<LinearAlgebra::distributed::Vector<double>>
{
public:
  MGTransferMF(const MGLevelObject<LAPLACEOPERATOR> &laplace,
               const MGConstrainedDoFs &             mg_constrained_dofs)
    : MGTransferPrebuilt<LinearAlgebra::distributed::Vector<double>>(
        mg_constrained_dofs)
    , laplace_operator(laplace)
  {}

  /**
   * Overload copy_to_mg from MGTransferPrebuilt to get the vectors compatible
   * with MatrixFree and bypass the crude vector initialization in
   * MGTransferPrebuilt
   */
  template <int dim, class InVector, int spacedim>
  void
  copy_to_mg(const DoFHandler<dim, spacedim> &mg_dof_handler,
             MGLevelObject<LinearAlgebra::distributed::Vector<double>> &dst,
             const InVector &src) const
  {
    for (unsigned int level = dst.min_level(); level <= dst.max_level();
         ++level)
      laplace_operator[level].initialize_dof_vector(dst[level]);
    MGTransferPrebuilt<LinearAlgebra::distributed::Vector<double>>::copy_to_mg(
      mg_dof_handler, dst, src);
  }

private:
  const MGLevelObject<LAPLACEOPERATOR> &laplace_operator;
};



template <typename MatrixType, typename Number>
class MGCoarseIterative
  : public MGCoarseGridBase<LinearAlgebra::distributed::Vector<Number>>
{
public:
  MGCoarseIterative()
  {}

  void
  initialize(const MatrixType &matrix)
  {
    coarse_matrix = &matrix;
  }

  virtual void
  operator()(const unsigned int                                level,
             LinearAlgebra::distributed::Vector<double> &      dst,
             const LinearAlgebra::distributed::Vector<double> &src) const
  {
    ReductionControl solver_control(1e4, 1e-50, 1e-10);
    SolverCG<LinearAlgebra::distributed::Vector<double>> solver_coarse(
      solver_control);
    solver_coarse.solve(*coarse_matrix, dst, src, PreconditionIdentity());
  }

  const MatrixType *coarse_matrix;
};



template <int dim, int fe_degree, int n_q_points_1d, typename number>
void
do_test(const DoFHandler<dim> &dof)
{
  if (std::is_same<number, float>::value == true)
    {
      deallog.push("float");
    }
  else
    {}

  deallog << "Testing " << dof.get_fe().get_name();
  deallog << std::endl;
  deallog << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;

  AffineConstraints<double> hanging_node_constraints;
  IndexSet                  locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof, locally_relevant_dofs);
  hanging_node_constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof, hanging_node_constraints);
  hanging_node_constraints.close();

  MGConstrainedDoFs                                   mg_constrained_dofs;
  Functions::ZeroFunction<dim>                        zero_function;
  std::map<types::boundary_id, const Function<dim> *> dirichlet_boundary;
  dirichlet_boundary[0] = &zero_function;
  mg_constrained_dofs.initialize(dof, dirichlet_boundary);

  MappingQ<dim>                                          mapping(fe_degree + 1);
  LaplaceOperator<dim, fe_degree, n_q_points_1d, number> fine_matrix;
  fine_matrix.initialize(mapping,
                         dof,
                         mg_constrained_dofs,
                         dirichlet_boundary,
                         numbers::invalid_unsigned_int);

  LinearAlgebra::distributed::Vector<number> in, sol;
  fine_matrix.initialize_dof_vector(in);
  fine_matrix.initialize_dof_vector(sol);

  // set constant rhs vector
  for (unsigned int i = 0; i < in.local_size(); ++i)
    if (!hanging_node_constraints.is_constrained(
          in.get_partitioner()->local_to_global(i)))
      in.local_element(i) = 1.;

  // set up multigrid in analogy to step-37
  typedef LaplaceOperator<dim, fe_degree, n_q_points_1d, number>
    LevelMatrixType;

  MGLevelObject<LevelMatrixType> mg_matrices;
  mg_matrices.resize(0, dof.get_triangulation().n_global_levels() - 1);
  for (unsigned int level = 0;
       level < dof.get_triangulation().n_global_levels();
       ++level)
    {
      mg_matrices[level].initialize(
        mapping, dof, mg_constrained_dofs, dirichlet_boundary, level);
    }
  MGLevelObject<MGInterfaceMatrix<LevelMatrixType>> mg_interface_matrices;
  mg_interface_matrices.resize(0,
                               dof.get_triangulation().n_global_levels() - 1);
  for (unsigned int level = 0;
       level < dof.get_triangulation().n_global_levels();
       ++level)
    mg_interface_matrices[level].initialize(mg_matrices[level]);

  MGTransferMF<LevelMatrixType> mg_transfer(mg_matrices, mg_constrained_dofs);
  mg_transfer.build(dof);

  MGCoarseIterative<LevelMatrixType, number> mg_coarse;
  mg_coarse.initialize(mg_matrices[0]);

  typedef PreconditionChebyshev<LevelMatrixType,
                                LinearAlgebra::distributed::Vector<number>>
    SMOOTHER;
  MGSmootherPrecondition<LevelMatrixType,
                         SMOOTHER,
                         LinearAlgebra::distributed::Vector<number>>
    mg_smoother;

  MGLevelObject<typename SMOOTHER::AdditionalData> smoother_data;
  smoother_data.resize(0, dof.get_triangulation().n_global_levels() - 1);
  for (unsigned int level = 0;
       level < dof.get_triangulation().n_global_levels();
       ++level)
    {
      smoother_data[level].smoothing_range     = 15.;
      smoother_data[level].degree              = 5;
      smoother_data[level].eig_cg_n_iterations = 15;
      smoother_data[level].preconditioner.reset(
        new DiagonalMatrix<LinearAlgebra::distributed::Vector<number>>());
      smoother_data[level].preconditioner->get_vector() =
        mg_matrices[level].get_matrix_diagonal_inverse();
    }
  mg_smoother.initialize(mg_matrices, smoother_data);

  mg::Matrix<LinearAlgebra::distributed::Vector<double>> mg_matrix(mg_matrices);
  mg::Matrix<LinearAlgebra::distributed::Vector<double>> mg_interface(
    mg_interface_matrices);

  Multigrid<LinearAlgebra::distributed::Vector<double>> mg(
    mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);
  mg.set_edge_matrices(mg_interface, mg_interface);
  PreconditionMG<dim,
                 LinearAlgebra::distributed::Vector<double>,
                 MGTransferMF<LevelMatrixType>>
    preconditioner(dof, mg, mg_transfer);

  {
    // avoid output from inner (coarse-level) solver
    deallog.depth_file(3);

    ReductionControl control(30, 1e-20, 1e-7);
    SolverCG<LinearAlgebra::distributed::Vector<double>> solver(control);
    solver.solve(fine_matrix, sol, in, preconditioner);
  }

  if (std::is_same<number, float>::value == true)
    deallog.pop();
}



template <int dim, int fe_degree>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(6 - dim);
  const unsigned int n_runs = fe_degree == 1 ? 6 - dim : 5 - dim;
  for (unsigned int i = 0; i < n_runs; ++i)
    {
      for (typename Triangulation<dim>::active_cell_iterator cell =
             tria.begin_active();
           cell != tria.end();
           ++cell)
        if (cell->is_locally_owned() &&
            ((cell->center().norm() < 0.5 &&
              (cell->level() < 5 || cell->center().norm() > 0.45)) ||
             (dim == 2 && cell->center().norm() > 1.2)))
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
      FE_Q<dim>       fe(fe_degree);
      DoFHandler<dim> dof(tria);
      dof.distribute_dofs(fe);
      dof.distribute_mg_dofs();

      do_test<dim, fe_degree, fe_degree + 1, double>(dof);
    }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
  mpi_initlog();
  deallog << std::setprecision(4);

  {
    deallog.push("2d");
    test<2, 1>();
    deallog.pop();
    deallog.push("3d");
    test<3, 1>();
    deallog.pop();
  }
}

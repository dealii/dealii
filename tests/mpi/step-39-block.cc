// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2018 by the deal.II authors
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

/*
 * Test RelaxationBlockJacobi in parallel
 */

#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/integrators/laplace.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/relaxation_block.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/loop.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"

namespace Step39
{
  Functions::SlitSingularityFunction<2> exact_solution;



  template <int dim>
  class MatrixIntegrator : public MeshWorker::LocalIntegrator<dim>
  {
  public:
    void
    cell(MeshWorker::DoFInfo<dim> &                 dinfo,
         typename MeshWorker::IntegrationInfo<dim> &info) const;
    void
    boundary(MeshWorker::DoFInfo<dim> &                 dinfo,
             typename MeshWorker::IntegrationInfo<dim> &info) const;
    void
    face(MeshWorker::DoFInfo<dim> &                 dinfo1,
         MeshWorker::DoFInfo<dim> &                 dinfo2,
         typename MeshWorker::IntegrationInfo<dim> &info1,
         typename MeshWorker::IntegrationInfo<dim> &info2) const;
  };


  template <int dim>
  void
  MatrixIntegrator<dim>::cell(
    MeshWorker::DoFInfo<dim> &                 dinfo,
    typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    LocalIntegrators::Laplace::cell_matrix(dinfo.matrix(0, false).matrix,
                                           info.fe_values());
  }


  template <int dim>
  void
  MatrixIntegrator<dim>::boundary(
    MeshWorker::DoFInfo<dim> &                 dinfo,
    typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    const unsigned int deg = info.fe_values(0).get_fe().tensor_degree();
    LocalIntegrators::Laplace::nitsche_matrix(
      dinfo.matrix(0, false).matrix,
      info.fe_values(0),
      LocalIntegrators::Laplace::compute_penalty(dinfo, dinfo, deg, deg));
  }

  template <int dim>
  void
  MatrixIntegrator<dim>::face(
    MeshWorker::DoFInfo<dim> &                 dinfo1,
    MeshWorker::DoFInfo<dim> &                 dinfo2,
    typename MeshWorker::IntegrationInfo<dim> &info1,
    typename MeshWorker::IntegrationInfo<dim> &info2) const
  {
    const unsigned int deg = info1.fe_values(0).get_fe().tensor_degree();
    LocalIntegrators::Laplace::ip_matrix(
      dinfo1.matrix(0, false).matrix,
      dinfo1.matrix(0, true).matrix,
      dinfo2.matrix(0, true).matrix,
      dinfo2.matrix(0, false).matrix,
      info1.fe_values(0),
      info2.fe_values(0),
      LocalIntegrators::Laplace::compute_penalty(dinfo1, dinfo2, deg, deg));
  }

  template <int dim>
  class RHSIntegrator : public MeshWorker::LocalIntegrator<dim>
  {
  public:
    void
    cell(MeshWorker::DoFInfo<dim> &                 dinfo,
         typename MeshWorker::IntegrationInfo<dim> &info) const;
    void
    boundary(MeshWorker::DoFInfo<dim> &                 dinfo,
             typename MeshWorker::IntegrationInfo<dim> &info) const;
    void
    face(MeshWorker::DoFInfo<dim> &                 dinfo1,
         MeshWorker::DoFInfo<dim> &                 dinfo2,
         typename MeshWorker::IntegrationInfo<dim> &info1,
         typename MeshWorker::IntegrationInfo<dim> &info2) const;
  };


  template <int dim>
  void
  RHSIntegrator<dim>::cell(MeshWorker::DoFInfo<dim> &,
                           typename MeshWorker::IntegrationInfo<dim> &) const
  {}


  template <int dim>
  void
  RHSIntegrator<dim>::boundary(
    MeshWorker::DoFInfo<dim> &                 dinfo,
    typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    const FEValuesBase<dim> &fe           = info.fe_values();
    Vector<double> &         local_vector = dinfo.vector(0).block(0);

    std::vector<double> boundary_values(fe.n_quadrature_points);
    exact_solution.value_list(fe.get_quadrature_points(), boundary_values);

    const unsigned int deg = fe.get_fe().tensor_degree();
    const double       penalty =
      2. * deg * (deg + 1) * dinfo.face->measure() / dinfo.cell->measure();

    for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
        local_vector(i) +=
          (-fe.shape_value(i, k) * penalty * boundary_values[k] +
           (fe.normal_vector(k) * fe.shape_grad(i, k)) * boundary_values[k]) *
          fe.JxW(k);
  }


  template <int dim>
  void
  RHSIntegrator<dim>::face(MeshWorker::DoFInfo<dim> &,
                           MeshWorker::DoFInfo<dim> &,
                           typename MeshWorker::IntegrationInfo<dim> &,
                           typename MeshWorker::IntegrationInfo<dim> &) const
  {}


  template <int dim>
  class Estimator : public MeshWorker::LocalIntegrator<dim>
  {
  public:
    void
    cell(MeshWorker::DoFInfo<dim> &                 dinfo,
         typename MeshWorker::IntegrationInfo<dim> &info) const;
    void
    boundary(MeshWorker::DoFInfo<dim> &                 dinfo,
             typename MeshWorker::IntegrationInfo<dim> &info) const;
    void
    face(MeshWorker::DoFInfo<dim> &                 dinfo1,
         MeshWorker::DoFInfo<dim> &                 dinfo2,
         typename MeshWorker::IntegrationInfo<dim> &info1,
         typename MeshWorker::IntegrationInfo<dim> &info2) const;
  };


  template <int dim>
  void
  Estimator<dim>::cell(MeshWorker::DoFInfo<dim> &                 dinfo,
                       typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    const FEValuesBase<dim> &fe = info.fe_values();

    const std::vector<Tensor<2, dim>> &DDuh = info.hessians[0][0];
    for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
      {
        const double t = dinfo.cell->diameter() * trace(DDuh[k]);
        dinfo.value(0) += t * t * fe.JxW(k);
      }
    dinfo.value(0) = std::sqrt(dinfo.value(0));
  }

  template <int dim>
  void
  Estimator<dim>::boundary(
    MeshWorker::DoFInfo<dim> &                 dinfo,
    typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    const FEValuesBase<dim> &fe = info.fe_values();

    std::vector<double> boundary_values(fe.n_quadrature_points);
    exact_solution.value_list(fe.get_quadrature_points(), boundary_values);

    const std::vector<double> &uh = info.values[0][0];

    const unsigned int deg = fe.get_fe().tensor_degree();
    const double       penalty =
      2. * deg * (deg + 1) * dinfo.face->measure() / dinfo.cell->measure();

    for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
      dinfo.value(0) += penalty * (boundary_values[k] - uh[k]) *
                        (boundary_values[k] - uh[k]) * fe.JxW(k);
    dinfo.value(0) = std::sqrt(dinfo.value(0));
  }


  template <int dim>
  void
  Estimator<dim>::face(MeshWorker::DoFInfo<dim> &                 dinfo1,
                       MeshWorker::DoFInfo<dim> &                 dinfo2,
                       typename MeshWorker::IntegrationInfo<dim> &info1,
                       typename MeshWorker::IntegrationInfo<dim> &info2) const
  {
    const FEValuesBase<dim> &          fe   = info1.fe_values();
    const std::vector<double> &        uh1  = info1.values[0][0];
    const std::vector<double> &        uh2  = info2.values[0][0];
    const std::vector<Tensor<1, dim>> &Duh1 = info1.gradients[0][0];
    const std::vector<Tensor<1, dim>> &Duh2 = info2.gradients[0][0];

    const unsigned int deg = fe.get_fe().tensor_degree();
    const double       penalty1 =
      deg * (deg + 1) * dinfo1.face->measure() / dinfo1.cell->measure();
    const double penalty2 =
      deg * (deg + 1) * dinfo2.face->measure() / dinfo2.cell->measure();
    const double penalty = penalty1 + penalty2;
    const double h       = dinfo1.face->measure();

    for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
      {
        double diff1 = uh1[k] - uh2[k];
        double diff2 =
          fe.normal_vector(k) * Duh1[k] - fe.normal_vector(k) * Duh2[k];
        dinfo1.value(0) +=
          (penalty * diff1 * diff1 + h * diff2 * diff2) * fe.JxW(k);
      }
    dinfo1.value(0) = std::sqrt(dinfo1.value(0));
    dinfo2.value(0) = dinfo1.value(0);
    // do not fill values if cells are ghost cells because we don't communicate
    if (!dinfo1.cell->is_locally_owned())
      dinfo1.value(0) = 0.0;
    if (!dinfo2.cell->is_locally_owned())
      dinfo2.value(0) = 0.0;
  }



  template <int dim>
  class ErrorIntegrator : public MeshWorker::LocalIntegrator<dim>
  {
  public:
    void
    cell(MeshWorker::DoFInfo<dim> &                 dinfo,
         typename MeshWorker::IntegrationInfo<dim> &info) const;
    void
    boundary(MeshWorker::DoFInfo<dim> &                 dinfo,
             typename MeshWorker::IntegrationInfo<dim> &info) const;
    void
    face(MeshWorker::DoFInfo<dim> &                 dinfo1,
         MeshWorker::DoFInfo<dim> &                 dinfo2,
         typename MeshWorker::IntegrationInfo<dim> &info1,
         typename MeshWorker::IntegrationInfo<dim> &info2) const;
  };


  template <int dim>
  void
  ErrorIntegrator<dim>::cell(
    MeshWorker::DoFInfo<dim> &                 dinfo,
    typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    const FEValuesBase<dim> &   fe = info.fe_values();
    std::vector<Tensor<1, dim>> exact_gradients(fe.n_quadrature_points);
    std::vector<double>         exact_values(fe.n_quadrature_points);

    exact_solution.gradient_list(fe.get_quadrature_points(), exact_gradients);
    exact_solution.value_list(fe.get_quadrature_points(), exact_values);

    const std::vector<Tensor<1, dim>> &Duh = info.gradients[0][0];
    const std::vector<double> &        uh  = info.values[0][0];

    for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
      {
        double sum = 0;
        for (unsigned int d = 0; d < dim; ++d)
          {
            const double diff = exact_gradients[k][d] - Duh[k][d];
            sum += diff * diff;
          }
        const double diff = exact_values[k] - uh[k];
        dinfo.value(0) += sum * fe.JxW(k);
        dinfo.value(1) += diff * diff * fe.JxW(k);
      }
    dinfo.value(0) = std::sqrt(dinfo.value(0));
    dinfo.value(1) = std::sqrt(dinfo.value(1));
  }


  template <int dim>
  void
  ErrorIntegrator<dim>::boundary(
    MeshWorker::DoFInfo<dim> &                 dinfo,
    typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    const FEValuesBase<dim> &fe = info.fe_values();

    std::vector<double> exact_values(fe.n_quadrature_points);
    exact_solution.value_list(fe.get_quadrature_points(), exact_values);

    const std::vector<double> &uh = info.values[0][0];

    const unsigned int deg = fe.get_fe().tensor_degree();
    const double       penalty =
      2. * deg * (deg + 1) * dinfo.face->measure() / dinfo.cell->measure();

    for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
      {
        const double diff = exact_values[k] - uh[k];
        dinfo.value(0) += penalty * diff * diff * fe.JxW(k);
      }
    dinfo.value(0) = std::sqrt(dinfo.value(0));
  }


  template <int dim>
  void
  ErrorIntegrator<dim>::face(
    MeshWorker::DoFInfo<dim> &                 dinfo1,
    MeshWorker::DoFInfo<dim> &                 dinfo2,
    typename MeshWorker::IntegrationInfo<dim> &info1,
    typename MeshWorker::IntegrationInfo<dim> &info2) const
  {
    const FEValuesBase<dim> &  fe  = info1.fe_values();
    const std::vector<double> &uh1 = info1.values[0][0];
    const std::vector<double> &uh2 = info2.values[0][0];

    const unsigned int deg = fe.get_fe().tensor_degree();
    const double       penalty1 =
      deg * (deg + 1) * dinfo1.face->measure() / dinfo1.cell->measure();
    const double penalty2 =
      deg * (deg + 1) * dinfo2.face->measure() / dinfo2.cell->measure();
    const double penalty = penalty1 + penalty2;

    for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
      {
        double diff = uh1[k] - uh2[k];
        dinfo1.value(0) += (penalty * diff * diff) * fe.JxW(k);
      }
    dinfo1.value(0) = std::sqrt(dinfo1.value(0));
    dinfo2.value(0) = dinfo1.value(0);
  }



  template <int dim>
  class InteriorPenaltyProblem
  {
  public:
    typedef MeshWorker::IntegrationInfo<dim> CellInfo;

    InteriorPenaltyProblem(const FiniteElement<dim> &fe);

    void
    run(unsigned int n_steps);

  private:
    void
    setup_system();
    void
    assemble_matrix();
    void
    assemble_mg_matrix();
    void
    assemble_right_hand_side();
    void
    error();
    double
    estimate();
    void
    solve();
    void
    output_results(const unsigned int cycle) const;

    parallel::distributed::Triangulation<dim> triangulation;
    const MappingQGeneric<dim>                mapping;
    const FiniteElement<dim> &                fe;
    DoFHandler<dim>                           dof_handler;

    IndexSet locally_relevant_set;

    TrilinosWrappers::SparseMatrix matrix;
    TrilinosWrappers::MPI::Vector  solution;
    TrilinosWrappers::MPI::Vector  right_hand_side;
    BlockVector<double>            estimates;

    MGLevelObject<TrilinosWrappers::SparseMatrix> mg_matrix;

    MGLevelObject<TrilinosWrappers::SparseMatrix> mg_matrix_dg_down;
    MGLevelObject<TrilinosWrappers::SparseMatrix> mg_matrix_dg_up;
  };


  template <int dim>
  InteriorPenaltyProblem<dim>::InteriorPenaltyProblem(
    const FiniteElement<dim> &fe)
    : triangulation(MPI_COMM_WORLD,
                    Triangulation<dim>::limit_level_difference_at_vertices,
                    parallel::distributed::Triangulation<
                      dim>::construct_multigrid_hierarchy)
    , mapping(1)
    , fe(fe)
    , dof_handler(triangulation)
    , estimates(1)
  {
    GridGenerator::hyper_L(triangulation, -1, 1);
  }


  template <int dim>
  void
  InteriorPenaltyProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    dof_handler.distribute_mg_dofs();

    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_set);
    solution.reinit(dof_handler.locally_owned_dofs(), MPI_COMM_WORLD);
    right_hand_side.reinit(dof_handler.locally_owned_dofs(), MPI_COMM_WORLD);

    DynamicSparsityPattern c_sparsity(dof_handler.n_dofs(),
                                      dof_handler.n_dofs());
    DoFTools::make_flux_sparsity_pattern(dof_handler, c_sparsity);
    matrix.reinit(dof_handler.locally_owned_dofs(),
                  c_sparsity,
                  MPI_COMM_WORLD,
                  true);

    const unsigned int n_levels = triangulation.n_global_levels();
    mg_matrix.resize(0, n_levels - 1);
    mg_matrix.clear_elements();
    mg_matrix_dg_up.resize(0, n_levels - 1);
    mg_matrix_dg_up.clear_elements();
    mg_matrix_dg_down.resize(0, n_levels - 1);
    mg_matrix_dg_down.clear_elements();

    for (unsigned int level = mg_matrix.min_level();
         level <= mg_matrix.max_level();
         ++level)
      {
        DynamicSparsityPattern c_sparsity(dof_handler.n_dofs(level));
        MGTools::make_flux_sparsity_pattern(dof_handler, c_sparsity, level);
        mg_matrix[level].reinit(dof_handler.locally_owned_mg_dofs(level),
                                dof_handler.locally_owned_mg_dofs(level),
                                c_sparsity,
                                MPI_COMM_WORLD,
                                true);

        if (level > 0)
          {
            DynamicSparsityPattern ci_sparsity;
            ci_sparsity.reinit(dof_handler.n_dofs(level - 1),
                               dof_handler.n_dofs(level));
            MGTools::make_flux_sparsity_pattern_edge(dof_handler,
                                                     ci_sparsity,
                                                     level);

            mg_matrix_dg_up[level].reinit(
              dof_handler.locally_owned_mg_dofs(level - 1),
              dof_handler.locally_owned_mg_dofs(level),
              ci_sparsity,
              MPI_COMM_WORLD,
              true);
            mg_matrix_dg_down[level].reinit(
              dof_handler.locally_owned_mg_dofs(level - 1),
              dof_handler.locally_owned_mg_dofs(level),
              ci_sparsity,
              MPI_COMM_WORLD,
              true);
          }
      }
  }


  template <int dim>
  void
  InteriorPenaltyProblem<dim>::assemble_matrix()
  {
    MeshWorker::IntegrationInfoBox<dim> info_box;
    UpdateFlags update_flags = update_values | update_gradients;
    info_box.add_update_flags_all(update_flags);
    info_box.initialize(fe, mapping);

    MeshWorker::DoFInfo<dim> dof_info(dof_handler);

    MeshWorker::Assembler::MatrixSimple<TrilinosWrappers::SparseMatrix>
      assembler;
    assembler.initialize(matrix);

    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator> begin(
      IteratorFilters::LocallyOwnedCell(), dof_handler.begin_active());
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator> end(
      IteratorFilters::LocallyOwnedCell(), dof_handler.end());

    MatrixIntegrator<dim> integrator;
    MeshWorker::integration_loop<dim, dim>(
      begin, end, dof_info, info_box, integrator, assembler);

    matrix.compress(VectorOperation::add);
  }


  template <int dim>
  void
  InteriorPenaltyProblem<dim>::assemble_mg_matrix()
  {
    MeshWorker::IntegrationInfoBox<dim> info_box;
    UpdateFlags update_flags = update_values | update_gradients;
    info_box.add_update_flags_all(update_flags);
    info_box.initialize(fe, mapping);

    MeshWorker::DoFInfo<dim> dof_info(dof_handler);

    MeshWorker::Assembler::MGMatrixSimple<TrilinosWrappers::SparseMatrix>
      assembler;
    assembler.initialize(mg_matrix);
    assembler.initialize_fluxes(mg_matrix_dg_up, mg_matrix_dg_down);

    FilteredIterator<typename DoFHandler<dim>::level_cell_iterator> begin(
      IteratorFilters::LocallyOwnedLevelCell(), dof_handler.begin());
    FilteredIterator<typename DoFHandler<dim>::level_cell_iterator> end(
      IteratorFilters::LocallyOwnedLevelCell(), dof_handler.end());

    MatrixIntegrator<dim> integrator;
    MeshWorker::integration_loop<dim, dim>(
      begin, end, dof_info, info_box, integrator, assembler);

    for (unsigned int level = mg_matrix.min_level();
         level <= mg_matrix.max_level();
         ++level)
      {
        mg_matrix[level].compress(VectorOperation::add);
        if (level > mg_matrix.min_level())
          {
            mg_matrix_dg_up[level].compress(VectorOperation::add);
            mg_matrix_dg_down[level].compress(VectorOperation::add);
          }
      }
  }


  template <int dim>
  void
  InteriorPenaltyProblem<dim>::assemble_right_hand_side()
  {
    MeshWorker::IntegrationInfoBox<dim> info_box;
    UpdateFlags                         update_flags =
      update_quadrature_points | update_values | update_gradients;
    info_box.add_update_flags_all(update_flags);
    info_box.initialize(fe, mapping);

    MeshWorker::DoFInfo<dim> dof_info(dof_handler);

    MeshWorker::Assembler::ResidualSimple<TrilinosWrappers::MPI::Vector>
            assembler;
    AnyData data;
    data.add<TrilinosWrappers::MPI::Vector *>(&right_hand_side, "RHS");
    assembler.initialize(data);

    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator> begin(
      IteratorFilters::LocallyOwnedCell(), dof_handler.begin_active());
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator> end(
      IteratorFilters::LocallyOwnedCell(), dof_handler.end());


    RHSIntegrator<dim> integrator;
    MeshWorker::integration_loop<dim, dim>(
      begin, end, dof_info, info_box, integrator, assembler);

    right_hand_side.compress(VectorOperation::add);
    right_hand_side *= -1.;
  }


  template <int dim>
  void
  InteriorPenaltyProblem<dim>::solve()
  {
    SolverControl                           control(1000, 1.e-12);
    SolverCG<TrilinosWrappers::MPI::Vector> solver(control);

    MGTransferPrebuilt<TrilinosWrappers::MPI::Vector> mg_transfer;
    mg_transfer.build(dof_handler);

    SolverControl coarse_solver_control(1000, 1e-10, false, false);
    SolverCG<TrilinosWrappers::MPI::Vector> coarse_solver(
      coarse_solver_control);
    PreconditionIdentity            identity;
    TrilinosWrappers::SparseMatrix &coarse_matrix = mg_matrix[0];
    MGCoarseGridLACIteration<SolverCG<TrilinosWrappers::MPI::Vector>,
                             TrilinosWrappers::MPI::Vector>
      coarse_grid_solver(coarse_solver, coarse_matrix, identity);

    typedef RelaxationBlockJacobi<TrilinosWrappers::SparseMatrix,
                                  double,
                                  TrilinosWrappers::MPI::Vector>
      Smoother;

    MGLevelObject<typename Smoother::AdditionalData> smoother_data;
    smoother_data.resize(0, triangulation.n_levels() - 1);
    mg::SmootherRelaxation<Smoother, TrilinosWrappers::MPI::Vector> mg_smoother;

    MGLevelObject<TrilinosWrappers::MPI::Vector> temp_vectors(
      0, triangulation.n_levels() - 1);


    for (unsigned int l = smoother_data.min_level() + 1;
         l <= smoother_data.max_level();
         ++l)
      {
        DoFTools::make_cell_patches(smoother_data[l].block_list,
                                    this->dof_handler,
                                    l);
        if (smoother_data[l].block_list.n_rows() > 0)
          smoother_data[l].block_list.compress();
        smoother_data[l].relaxation = 0.7;
        smoother_data[l].inversion  = PreconditionBlockBase<double>::svd;
        TrilinosWrappers::MPI::Vector *ghost = &(temp_vectors[l]);
        IndexSet                       relevant_dofs;
        DoFTools::extract_locally_relevant_level_dofs(dof_handler,
                                                      l,
                                                      relevant_dofs);
        ghost->reinit(dof_handler.locally_owned_mg_dofs(l),
                      relevant_dofs,
                      MPI_COMM_WORLD);
        smoother_data[l].temp_ghost_vector = ghost;
      }

    mg_smoother.initialize(mg_matrix, smoother_data);
    mg_smoother.set_steps(2);
    mg_smoother.set_variable(false);

    mg::Matrix<TrilinosWrappers::MPI::Vector> mgmatrix(mg_matrix);
    mg::Matrix<TrilinosWrappers::MPI::Vector> mgdown(mg_matrix_dg_down);
    mg::Matrix<TrilinosWrappers::MPI::Vector> mgup(mg_matrix_dg_up);

    Multigrid<TrilinosWrappers::MPI::Vector> mg(
      mgmatrix, coarse_grid_solver, mg_transfer, mg_smoother, mg_smoother);
    mg.set_edge_flux_matrices(mgdown, mgup);

    PreconditionMG<dim,
                   TrilinosWrappers::MPI::Vector,
                   MGTransferPrebuilt<TrilinosWrappers::MPI::Vector>>
      preconditioner(dof_handler, mg, mg_transfer);
    solver.solve(matrix, solution, right_hand_side, preconditioner);
  }


  template <int dim>
  double
  InteriorPenaltyProblem<dim>::estimate()
  {
    TrilinosWrappers::MPI::Vector ghost;
    ghost.reinit(locally_relevant_set, MPI_COMM_WORLD);
    ghost = solution;

    std::vector<unsigned int> old_user_indices;
    triangulation.save_user_indices(old_user_indices);

    estimates.block(0).reinit(triangulation.n_active_cells());
    unsigned int i = 0;
    for (typename Triangulation<dim>::active_cell_iterator cell =
           triangulation.begin_active();
         cell != triangulation.end();
         ++cell, ++i)
      cell->set_user_index(i);

    MeshWorker::IntegrationInfoBox<dim> info_box;
    const unsigned int                  n_gauss_points =
      dof_handler.get_fe().tensor_degree() + 1;
    info_box.initialize_gauss_quadrature(n_gauss_points,
                                         n_gauss_points + 1,
                                         n_gauss_points);

    AnyData solution_data;
    solution_data.add<TrilinosWrappers::MPI::Vector *>(&ghost, "solution");

    info_box.cell_selector.add("solution", false, false, true);
    info_box.boundary_selector.add("solution", true, true, false);
    info_box.face_selector.add("solution", true, true, false);

    info_box.add_update_flags_boundary(update_quadrature_points);
    info_box.initialize(fe, mapping, solution_data, solution);

    MeshWorker::DoFInfo<dim> dof_info(dof_handler);

    MeshWorker::Assembler::CellsAndFaces<double> assembler;
    AnyData                                      out_data;
    out_data.add<BlockVector<double> *>(&estimates, "cells");
    assembler.initialize(out_data, false);

    Estimator<dim>          integrator;
    MeshWorker::LoopControl lctrl;
    // assemble all faces adjacent to ghost cells to get the full
    // information for all own cells without communication
    lctrl.faces_to_ghost = MeshWorker::LoopControl::both;

    MeshWorker::integration_loop<dim, dim>(dof_handler.begin_active(),
                                           dof_handler.end(),
                                           dof_info,
                                           info_box,
                                           integrator,
                                           assembler,
                                           lctrl);

    triangulation.load_user_indices(old_user_indices);
    // estimates is a BlockVector<double> (so serial) on each processor
    // with one entry per active cell. Note that only the locally owned
    // cells are !=0, so summing the contributions of l2_norm() over all
    // processors is the right way to do this.
    double local_norm = estimates.block(0).l2_norm();
    local_norm *= local_norm;
    return std::sqrt(Utilities::MPI::sum(local_norm, MPI_COMM_WORLD));
  }

  template <int dim>
  void
  InteriorPenaltyProblem<dim>::run(unsigned int n_steps)
  {
    deallog << "Element: " << fe.get_name() << std::endl;
    for (unsigned int s = 0; s < n_steps; ++s)
      {
        deallog << "Step " << s << std::endl;
        if (estimates.block(0).size() == 0)
          triangulation.refine_global(1);
        else
          {
            parallel::distributed::GridRefinement::
              refine_and_coarsen_fixed_fraction(triangulation,
                                                estimates.block(0),
                                                0.5,
                                                0.0);
            triangulation.execute_coarsening_and_refinement();
          }

        deallog << "Triangulation " << triangulation.n_global_active_cells()
                << " cells, " << triangulation.n_global_levels() << " levels"
                << std::endl;

        setup_system();
        deallog << "DoFHandler " << dof_handler.n_dofs() << " dofs, level dofs";
        for (unsigned int l = 0; l < triangulation.n_global_levels(); ++l)
          deallog << ' ' << dof_handler.n_dofs(l);
        deallog << std::endl;

        deallog << "Assemble matrix" << std::endl;
        assemble_matrix();
        deallog << "Assemble multilevel matrix" << std::endl;
        assemble_mg_matrix();
        deallog << "Assemble right hand side" << std::endl;
        assemble_right_hand_side();
        deallog << "Solve" << std::endl;
        solve();
        // error();
        deallog << "Estimate " << estimate() << std::endl;
        // output_results(s);
      }
  }
} // namespace Step39



int
main(int argc, char *argv[])
{
  using namespace Step39;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
  MPILogInitAll log;

  FE_DGQ<2>                 fe1(2);
  InteriorPenaltyProblem<2> test1(fe1);
  test1.run(6);
  return 0;
}

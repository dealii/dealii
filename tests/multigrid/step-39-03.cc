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


// Multigrid with continuous and discontinuous elements works, if we
// enforce continuity at refinement edges through interior penalty

// TH: the test was failing because the solver did not converge. It turns out
//     that switching from CG to GMRES makes everything work. I have no idea
//     if the problem is not SPD...

#include "../tests.h"
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/block_vector.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/multigrid/mg_dof_handler.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/loop.h>

#include <deal.II/integrators/laplace.h>

#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>

#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <iostream>
#include <fstream>

namespace Step39
{
  using namespace dealii;

  Functions::SlitSingularityFunction<2> exact_solution;




  template <int dim>
  class MatrixIntegrator : public MeshWorker::LocalIntegrator<dim>
  {
  public:
    void cell(MeshWorker::DoFInfo<dim> &dinfo,
              typename MeshWorker::IntegrationInfo<dim> &info) const;
    void boundary(MeshWorker::DoFInfo<dim> &dinfo,
                  typename MeshWorker::IntegrationInfo<dim> &info) const;
    void face(MeshWorker::DoFInfo<dim> &dinfo1,
              MeshWorker::DoFInfo<dim> &dinfo2,
              typename MeshWorker::IntegrationInfo<dim> &info1,
              typename MeshWorker::IntegrationInfo<dim> &info2) const;
  };


  template <int dim>
  void MatrixIntegrator<dim>::cell(
    MeshWorker::DoFInfo<dim> &dinfo,
    typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    LocalIntegrators::Laplace::cell_matrix(dinfo.matrix(0,false).matrix, info.fe_values());
  }


  template <int dim>
  void MatrixIntegrator<dim>::boundary(
    MeshWorker::DoFInfo<dim> &dinfo,
    typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    const unsigned int deg = info.fe_values(0).get_fe().tensor_degree();
    LocalIntegrators::Laplace::nitsche_matrix(
      dinfo.matrix(0,false).matrix, info.fe_values(0),
      LocalIntegrators::Laplace::compute_penalty(dinfo, dinfo, deg, deg));
  }

  template <int dim>
  void MatrixIntegrator<dim>::face(
    MeshWorker::DoFInfo<dim> &dinfo1,
    MeshWorker::DoFInfo<dim> &dinfo2,
    typename MeshWorker::IntegrationInfo<dim> &info1,
    typename MeshWorker::IntegrationInfo<dim> &info2) const
  {
    const unsigned int deg = info1.fe_values(0).get_fe().tensor_degree();
    LocalIntegrators::Laplace::ip_matrix(
      dinfo1.matrix(0,false).matrix, dinfo1.matrix(0,true).matrix,
      dinfo2.matrix(0,true).matrix, dinfo2.matrix(0,false).matrix,
      info1.fe_values(0), info2.fe_values(0),
      LocalIntegrators::Laplace::compute_penalty(dinfo1, dinfo2, deg, deg));
  }

  template <int dim>
  class RHSIntegrator : public MeshWorker::LocalIntegrator<dim>
  {
  public:
    void cell(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const;
    void boundary(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const;
    void face(MeshWorker::DoFInfo<dim> &dinfo1,
              MeshWorker::DoFInfo<dim> &dinfo2,
              typename MeshWorker::IntegrationInfo<dim> &info1,
              typename MeshWorker::IntegrationInfo<dim> &info2) const;
  };


  template <int dim>
  void RHSIntegrator<dim>::cell(MeshWorker::DoFInfo<dim> &, typename MeshWorker::IntegrationInfo<dim> &) const
  {}


  template <int dim>
  void RHSIntegrator<dim>::boundary(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    const FEValuesBase<dim> &fe = info.fe_values();
    Vector<double> &local_vector = dinfo.vector(0).block(0);

    std::vector<double> boundary_values(fe.n_quadrature_points);
    exact_solution.value_list(fe.get_quadrature_points(), boundary_values);

    const unsigned int deg = fe.get_fe().tensor_degree();
    const double penalty = 2. * deg * (deg+1) * dinfo.face->measure() / dinfo.cell->measure();

    for (unsigned k=0; k<fe.n_quadrature_points; ++k)
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        local_vector(i) += (- fe.shape_value(i,k) * penalty * boundary_values[k]
                            + (fe.normal_vector(k) * fe.shape_grad(i,k)) * boundary_values[k])
                           * fe.JxW(k);
  }


  template <int dim>
  void RHSIntegrator<dim>::face(MeshWorker::DoFInfo<dim> &,
                                MeshWorker::DoFInfo<dim> &,
                                typename MeshWorker::IntegrationInfo<dim> &,
                                typename MeshWorker::IntegrationInfo<dim> &) const
  {}


  template <int dim>
  class Estimator : public MeshWorker::LocalIntegrator<dim>
  {
  public:
    void cell(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const;
    void boundary(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const;
    void face(MeshWorker::DoFInfo<dim> &dinfo1,
              MeshWorker::DoFInfo<dim> &dinfo2,
              typename MeshWorker::IntegrationInfo<dim> &info1,
              typename MeshWorker::IntegrationInfo<dim> &info2) const;
  };


  template <int dim>
  void Estimator<dim>::cell(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    const FEValuesBase<dim> &fe = info.fe_values();

    const std::vector<Tensor<2,dim> > &DDuh = info.hessians[0][0];
    for (unsigned k=0; k<fe.n_quadrature_points; ++k)
      {
        const double t = dinfo.cell->diameter() * trace(DDuh[k]);
        dinfo.value(0) +=  t*t * fe.JxW(k);
      }
    dinfo.value(0) = std::sqrt(dinfo.value(0));
  }

  template <int dim>
  void Estimator<dim>::boundary(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    const FEValuesBase<dim> &fe = info.fe_values();

    std::vector<double> boundary_values(fe.n_quadrature_points);
    exact_solution.value_list(fe.get_quadrature_points(), boundary_values);

    const std::vector<double> &uh = info.values[0][0];

    const unsigned int deg = fe.get_fe().tensor_degree();
    const double penalty = 2. * deg * (deg+1) * dinfo.face->measure() / dinfo.cell->measure();

    for (unsigned k=0; k<fe.n_quadrature_points; ++k)
      dinfo.value(0) += penalty * (boundary_values[k] - uh[k]) * (boundary_values[k] - uh[k])
                        * fe.JxW(k);
    dinfo.value(0) = std::sqrt(dinfo.value(0));
  }


  template <int dim>
  void Estimator<dim>::face(MeshWorker::DoFInfo<dim> &dinfo1,
                            MeshWorker::DoFInfo<dim> &dinfo2,
                            typename MeshWorker::IntegrationInfo<dim> &info1,
                            typename MeshWorker::IntegrationInfo<dim> &info2) const
  {
    const FEValuesBase<dim> &fe = info1.fe_values();
    const std::vector<double> &uh1 = info1.values[0][0];
    const std::vector<double> &uh2 = info2.values[0][0];
    const std::vector<Tensor<1,dim> > &Duh1 = info1.gradients[0][0];
    const std::vector<Tensor<1,dim> > &Duh2 = info2.gradients[0][0];

    const unsigned int deg = fe.get_fe().tensor_degree();
    const double penalty1 = deg * (deg+1) * dinfo1.face->measure() / dinfo1.cell->measure();
    const double penalty2 = deg * (deg+1) * dinfo2.face->measure() / dinfo2.cell->measure();
    const double penalty = penalty1 + penalty2;
    const double h = dinfo1.face->measure();

    for (unsigned k=0; k<fe.n_quadrature_points; ++k)
      {
        double diff1 = uh1[k] - uh2[k];
        double diff2 = fe.normal_vector(k) * Duh1[k] - fe.normal_vector(k) * Duh2[k];
        dinfo1.value(0) += (penalty * diff1*diff1 + h * diff2*diff2)
                           * fe.JxW(k);
      }
    dinfo1.value(0) = std::sqrt(dinfo1.value(0));
    dinfo2.value(0) = dinfo1.value(0);
  }



  template <int dim>
  class ErrorIntegrator : public MeshWorker::LocalIntegrator<dim>
  {
  public:
    void cell(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const;
    void boundary(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const;
    void face(MeshWorker::DoFInfo<dim> &dinfo1,
              MeshWorker::DoFInfo<dim> &dinfo2,
              typename MeshWorker::IntegrationInfo<dim> &info1,
              typename MeshWorker::IntegrationInfo<dim> &info2) const;
  };


  template <int dim>
  void ErrorIntegrator<dim>::cell(
    MeshWorker::DoFInfo<dim> &dinfo,
    typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    const FEValuesBase<dim> &fe = info.fe_values();
    std::vector<Tensor<1,dim> > exact_gradients(fe.n_quadrature_points);
    std::vector<double> exact_values(fe.n_quadrature_points);

    exact_solution.gradient_list(fe.get_quadrature_points(), exact_gradients);
    exact_solution.value_list(fe.get_quadrature_points(), exact_values);

    const std::vector<Tensor<1,dim> > &Duh = info.gradients[0][0];
    const std::vector<double> &uh = info.values[0][0];

    for (unsigned k=0; k<fe.n_quadrature_points; ++k)
      {
        double sum = 0;
        for (unsigned int d=0; d<dim; ++d)
          {
            const double diff = exact_gradients[k][d] - Duh[k][d];
            sum += diff*diff;
          }
        const double diff = exact_values[k] - uh[k];
        dinfo.value(0) +=  sum * fe.JxW(k);
        dinfo.value(1) +=  diff*diff * fe.JxW(k);
      }
    dinfo.value(0) = std::sqrt(dinfo.value(0));
    dinfo.value(1) = std::sqrt(dinfo.value(1));
  }


  template <int dim>
  void ErrorIntegrator<dim>::boundary(
    MeshWorker::DoFInfo<dim> &dinfo,
    typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    const FEValuesBase<dim> &fe = info.fe_values();

    std::vector<double> exact_values(fe.n_quadrature_points);
    exact_solution.value_list(fe.get_quadrature_points(), exact_values);

    const std::vector<double> &uh = info.values[0][0];

    const unsigned int deg = fe.get_fe().tensor_degree();
    const double penalty = 2. * deg * (deg+1) * dinfo.face->measure() / dinfo.cell->measure();

    for (unsigned k=0; k<fe.n_quadrature_points; ++k)
      {
        const double diff = exact_values[k] - uh[k];
        dinfo.value(0) += penalty * diff * diff * fe.JxW(k);
      }
    dinfo.value(0) = std::sqrt(dinfo.value(0));
  }


  template <int dim>
  void ErrorIntegrator<dim>::face(
    MeshWorker::DoFInfo<dim> &dinfo1,
    MeshWorker::DoFInfo<dim> &dinfo2,
    typename MeshWorker::IntegrationInfo<dim> &info1,
    typename MeshWorker::IntegrationInfo<dim> &info2) const
  {
    const FEValuesBase<dim> &fe = info1.fe_values();
    const std::vector<double> &uh1 = info1.values[0][0];
    const std::vector<double> &uh2 = info2.values[0][0];

    const unsigned int deg = fe.get_fe().tensor_degree();
    const double penalty1 = deg * (deg+1) * dinfo1.face->measure() / dinfo1.cell->measure();
    const double penalty2 = deg * (deg+1) * dinfo2.face->measure() / dinfo2.cell->measure();
    const double penalty = penalty1 + penalty2;

    for (unsigned k=0; k<fe.n_quadrature_points; ++k)
      {
        double diff = uh1[k] - uh2[k];
        dinfo1.value(0) += (penalty * diff*diff)
                           * fe.JxW(k);
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

    void run(unsigned int n_steps);

  private:
    void setup_system ();
    void assemble_matrix ();
    void assemble_mg_matrix ();
    void assemble_right_hand_side ();
    void error ();
    double estimate ();
    void solve ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>        triangulation;
    const MappingQ1<dim>      mapping;
    const FiniteElement<dim> &fe;
    MGDoFHandler<dim>         mg_dof_handler;
    DoFHandler<dim>          &dof_handler;

    SparsityPattern      sparsity;
    SparseMatrix<double> matrix;
    Vector<double>       solution;
    Vector<double>       right_hand_side;
    BlockVector<double>  estimates;

    MGLevelObject<SparsityPattern> mg_sparsity;
    MGLevelObject<SparseMatrix<double> > mg_matrix;

    MGLevelObject<SparsityPattern> mg_sparsity_dg_interface;
    MGLevelObject<SparseMatrix<double> > mg_matrix_dg_down;
    MGLevelObject<SparseMatrix<double> > mg_matrix_dg_up;
    MGLevelObject<SparseMatrix<double> > mg_matrix_in_out;
  };


  template <int dim>
  InteriorPenaltyProblem<dim>::InteriorPenaltyProblem(const FiniteElement<dim> &fe)
    :
    mapping(),
    fe(fe),
    mg_dof_handler(triangulation),
    dof_handler(mg_dof_handler),
    estimates(1)
  {
    GridGenerator::hyper_cube_slit(triangulation, -1, 1);
  }


  template <int dim>
  void
  InteriorPenaltyProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    types::global_dof_index n_dofs = dof_handler.n_dofs();
    solution.reinit(n_dofs);
    right_hand_side.reinit(n_dofs);

    CompressedSparsityPattern c_sparsity(n_dofs);
    DoFTools::make_flux_sparsity_pattern(dof_handler, c_sparsity);
    sparsity.copy_from(c_sparsity);
    matrix.reinit(sparsity);

    const unsigned int n_levels = triangulation.n_levels();
    mg_matrix.resize(0, n_levels-1);
    mg_matrix.clear();
    mg_matrix_dg_up.resize(0, n_levels-1);
    mg_matrix_dg_up.clear();
    mg_matrix_dg_down.resize(0, n_levels-1);
    mg_matrix_dg_down.clear();
    mg_matrix_in_out.resize(0, n_levels-1);
    mg_matrix_in_out.clear();
    mg_sparsity.resize(0, n_levels-1);
    mg_sparsity_dg_interface.resize(0, n_levels-1);

    for (unsigned int level=mg_sparsity.min_level();
         level<=mg_sparsity.max_level(); ++level)
      {
        CompressedSparsityPattern c_sparsity(mg_dof_handler.n_dofs(level));
        MGTools::make_flux_sparsity_pattern(mg_dof_handler, c_sparsity, level);
        mg_sparsity[level].copy_from(c_sparsity);
        mg_matrix[level].reinit(mg_sparsity[level]);
        mg_matrix_in_out[level].reinit(mg_sparsity[level]);

        if (level>0)
          {
            CompressedSparsityPattern ci_sparsity;
            ci_sparsity.reinit(mg_dof_handler.n_dofs(level-1), mg_dof_handler.n_dofs(level));
            MGTools::make_flux_sparsity_pattern_edge(mg_dof_handler, ci_sparsity, level);
            mg_sparsity_dg_interface[level].copy_from(ci_sparsity);
            mg_matrix_dg_up[level].reinit(mg_sparsity_dg_interface[level]);
            mg_matrix_dg_down[level].reinit(mg_sparsity_dg_interface[level]);
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

    MeshWorker::Assembler::MatrixSimple<SparseMatrix<double> > assembler;
    assembler.initialize(matrix);

    MatrixIntegrator<dim> integrator;
    MeshWorker::integration_loop<dim, dim>(
      dof_handler.begin_active(), dof_handler.end(),
      dof_info, info_box,
      integrator, assembler);
  }


  template <int dim>
  void
  InteriorPenaltyProblem<dim>::assemble_mg_matrix()
  {
    MeshWorker::IntegrationInfoBox<dim> info_box;
    UpdateFlags update_flags = update_values | update_gradients;
    info_box.add_update_flags_all(update_flags);
    info_box.initialize(fe, mapping);

    MeshWorker::DoFInfo<dim> dof_info(mg_dof_handler);

    MeshWorker::Assembler::MGMatrixSimple<SparseMatrix<double> > assembler;
    assembler.initialize(mg_matrix);
    assembler.initialize_interfaces(mg_matrix_in_out, mg_matrix_in_out);
    assembler.initialize_fluxes(mg_matrix_dg_up, mg_matrix_dg_down);

    MatrixIntegrator<dim> integrator;
    MeshWorker::integration_loop<dim, dim> (
      mg_dof_handler.begin(), mg_dof_handler.end(),
      dof_info, info_box,
      integrator, assembler);

    for (unsigned int level=mg_matrix_in_out.min_level();
         level<=mg_matrix_in_out.min_level(); ++level)
      if (mg_matrix_in_out[level].frobenius_norm() != 0.)
        deallog << "Oops!" << std::endl;
  }


  template <int dim>
  void
  InteriorPenaltyProblem<dim>::assemble_right_hand_side()
  {
    MeshWorker::IntegrationInfoBox<dim> info_box;
    UpdateFlags update_flags = update_quadrature_points | update_values | update_gradients;
    info_box.add_update_flags_all(update_flags);
    info_box.initialize(fe, mapping);

    MeshWorker::DoFInfo<dim> dof_info(dof_handler);

    MeshWorker::Assembler::ResidualSimple<Vector<double> > assembler;
    NamedData<Vector<double>* > data;
    Vector<double> *rhs = &right_hand_side;
    data.add(rhs, "RHS");
    assembler.initialize(data);

    RHSIntegrator<dim> integrator;
    MeshWorker::integration_loop<dim, dim>(
      dof_handler.begin_active(), dof_handler.end(),
      dof_info, info_box,
      integrator, assembler);

    right_hand_side *= -1.;
  }


  template <int dim>
  void
  InteriorPenaltyProblem<dim>::solve()
  {
    SolverControl control(1000, 1.e-12);
    SolverGMRES<Vector<double> > solver(control);

    MGTransferPrebuilt<Vector<double> > mg_transfer;
    mg_transfer.build_matrices(mg_dof_handler);

    FullMatrix<double> coarse_matrix;
    coarse_matrix.copy_from (mg_matrix[0]);
    MGCoarseGridHouseholder<double, Vector<double> > mg_coarse;
    mg_coarse.initialize(coarse_matrix);

    GrowingVectorMemory<Vector<double> > mem;
    typedef PreconditionSOR<SparseMatrix<double> > RELAXATION;
    mg::SmootherRelaxation<RELAXATION, Vector<double> >
    mg_smoother;
    RELAXATION::AdditionalData smoother_data(1.);
    mg_smoother.initialize(mg_matrix, smoother_data);

    mg_smoother.set_steps(2);
    mg_smoother.set_symmetric(true);
    mg_smoother.set_variable(false);

    MGMatrix<SparseMatrix<double>, Vector<double> > mgmatrix(&mg_matrix);
    MGMatrix<SparseMatrix<double>, Vector<double> > mgdown(&mg_matrix_dg_down);
    MGMatrix<SparseMatrix<double>, Vector<double> > mgup(&mg_matrix_dg_up);
    MGMatrix<SparseMatrix<double>, Vector<double> > mgedge(&mg_matrix_in_out);

    Multigrid<Vector<double> > mg(mg_dof_handler, mgmatrix,
                                  mg_coarse, mg_transfer,
                                  mg_smoother, mg_smoother);
    mg.set_edge_flux_matrices(mgdown, mgup);
    mg.set_edge_matrices(mgedge, mgedge);

    PreconditionMG<dim, Vector<double>,
                   MGTransferPrebuilt<Vector<double> > >
                   preconditioner(mg_dof_handler, mg, mg_transfer);
    solver.solve(matrix, solution, right_hand_side, preconditioner);
  }


  template <int dim>
  double
  InteriorPenaltyProblem<dim>::estimate()
  {
    std::vector<unsigned int> old_user_indices;
    triangulation.save_user_indices(old_user_indices);

    estimates.block(0).reinit(triangulation.n_active_cells());
    unsigned int i=0;
    for (typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell,++i)
      cell->set_user_index(i);

    MeshWorker::IntegrationInfoBox<dim> info_box;
    const unsigned int n_gauss_points = dof_handler.get_fe().tensor_degree()+1;
    info_box.initialize_gauss_quadrature(n_gauss_points, n_gauss_points+1, n_gauss_points);

    NamedData<Vector<double>* > solution_data;
    solution_data.add(&solution, "solution");

    info_box.cell_selector.add("solution", false, false, true);
    info_box.boundary_selector.add("solution", true, true, false);
    info_box.face_selector.add("solution", true, true, false);

    info_box.add_update_flags_boundary(update_quadrature_points);
    info_box.initialize(fe, mapping, solution_data);

    MeshWorker::DoFInfo<dim> dof_info(dof_handler);

    MeshWorker::Assembler::CellsAndFaces<double> assembler;
    NamedData<BlockVector<double>* > out_data;
    BlockVector<double> *est = &estimates;
    out_data.add(est, "cells");
    assembler.initialize(out_data, false);

    Estimator<dim> integrator;
    MeshWorker::integration_loop<dim, dim> (
      dof_handler.begin_active(), dof_handler.end(),
      dof_info, info_box,
      integrator, assembler);

    triangulation.load_user_indices(old_user_indices);
    return estimates.block(0).l2_norm();
  }


  template <int dim>
  void
  InteriorPenaltyProblem<dim>::error()
  {
    BlockVector<double> errors(2);
    errors.block(0).reinit(triangulation.n_active_cells());
    errors.block(1).reinit(triangulation.n_active_cells());
    unsigned int i=0;
    for (typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell,++i)
      cell->set_user_index(i);

    MeshWorker::IntegrationInfoBox<dim> info_box;
    const unsigned int n_gauss_points = dof_handler.get_fe().tensor_degree()+1;
    info_box.initialize_gauss_quadrature(n_gauss_points, n_gauss_points+1, n_gauss_points);

    NamedData<Vector<double>* > solution_data;
    solution_data.add(&solution, "solution");

    info_box.cell_selector.add("solution", true, true, false);
    info_box.boundary_selector.add("solution", true, false, false);
    info_box.face_selector.add("solution", true, false, false);

    info_box.add_update_flags_cell(update_quadrature_points);
    info_box.add_update_flags_boundary(update_quadrature_points);
    info_box.initialize(fe, mapping, solution_data);

    MeshWorker::DoFInfo<dim> dof_info(dof_handler);

    MeshWorker::Assembler::CellsAndFaces<double> assembler;
    NamedData<BlockVector<double>* > out_data;
    BlockVector<double> *est = &errors;
    out_data.add(est, "cells");
    assembler.initialize(out_data, false);

    ErrorIntegrator<dim> integrator;
    MeshWorker::integration_loop<dim, dim> (
      dof_handler.begin_active(), dof_handler.end(),
      dof_info, info_box,
      integrator, assembler);

    deallog << "energy-error: " << errors.block(0).l2_norm() << std::endl;
    deallog << "L2-error:     " << errors.block(1).l2_norm() << std::endl;
  }


  template <int dim>
  void InteriorPenaltyProblem<dim>::output_results (const unsigned int cycle) const
  {
    char *fn = new char[100];
    sprintf(fn, "sol-%02d", cycle);

    std::string filename(fn);
    filename += ".gnuplot";
    deallog << "Writing solution to <" << filename << ">..."
            << std::endl << std::endl;
    std::ofstream gnuplot_output (filename.c_str());

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "u");
    data_out.add_data_vector (estimates.block(0), "est");

    data_out.build_patches ();

    data_out.write_gnuplot(gnuplot_output);
  }

  template <int dim>
  void
  InteriorPenaltyProblem<dim>::run(unsigned int n_steps)
  {
    deallog << "Element: " << fe.get_name() << std::endl;
    for (unsigned int s=0; s<n_steps; ++s)
      {
        deallog << "Step " << s << std::endl;
        if (estimates.block(0).size() == 0)
          triangulation.refine_global(1);
        else
          {
            GridRefinement::refine_and_coarsen_fixed_fraction (triangulation,
                                                               estimates.block(0),
                                                               0.5, 0.0);
            triangulation.execute_coarsening_and_refinement ();
          }

        deallog << "Triangulation "
                << triangulation.n_active_cells() << " cells, "
                << triangulation.n_levels() << " levels" << std::endl;

        setup_system();
        deallog << "DoFHandler " << dof_handler.n_dofs() << " dofs, level dofs";
        for (unsigned int l=0; l<triangulation.n_levels(); ++l)
          deallog << ' ' << mg_dof_handler.n_dofs(l);
        deallog << std::endl;

        deallog << "Assemble matrix" << std::endl;
        assemble_matrix();
        deallog << "Assemble multilevel matrix" << std::endl;
        assemble_mg_matrix();
        deallog << "Assemble right hand side" << std::endl;
        assemble_right_hand_side();
        deallog << "Solve" << std::endl;
        solve();
        error();
        deallog << "Estimate " << estimate() << std::endl;
        output_results(s);
      }
  }
}



int main()
{
  try
    {
      using namespace dealii;
      using namespace Step39;
      initlog(__FILE__);

      FESystem<2> fe1(FE_DGQ<2>(2), 1, FE_Q<2>(2), 1);
      InteriorPenaltyProblem<2> test1(fe1);
      test1.run(6);
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}

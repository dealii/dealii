/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the deal.II authors
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

 *
 * Author: Thomas C. Clevenger, Clemson University
 *         Timo Heister, Clemson University
 *         Guido Kanschat, Heidelberg University
 *         Martin Kronbichler, TU Munich
 */

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/meshworker/mesh_loop.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>


// uncomment the following #define if you have PETSc and Trilinos installed
// and you prefer using Trilinos in this example:
#define FORCE_USE_OF_TRILINOS

namespace LA
{
#if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \
  !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
  using namespace dealii::LinearAlgebraPETSc;
#  define USE_PETSC_LA
#elif defined(DEAL_II_WITH_TRILINOS)
  using namespace dealii::LinearAlgebraTrilinos;
#else
#  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required
#endif
} // namespace LA


using namespace dealii;

template <int dim>
class RightHandSide : public Function<dim>
{
public:
  virtual double value(const Point<dim> & /*p*/,
                       const unsigned int /*component*/ = 0) const override
  {
    return 1.0;
  }
};



template <int dim>
class Coefficient : public Function<dim>
{
public:
  virtual double value(const Point<dim> & p,
                       const unsigned int component = 0) const override;
};


template <int dim>
double Coefficient<dim>::value(const Point<dim> &p, const unsigned int) const
{
  for (int d = 0; d < dim; ++d)
    {
      if (p[d] < -0.5)
        return 100.0;
    }
  return 1.0;
}



void average(std::vector<double> &values)
{
  double sum = 0.0;
  for (unsigned int i = 0; i < values.size(); ++i)
    sum += values[i];
  sum /= values.size();

  for (unsigned int i = 0; i < values.size(); ++i)
    values[i] = sum;
}



struct Settings
{
  bool try_parse(const std::string &prm_filename);

  enum AssembleEnum
  {
    gmg,
    amg
  } assembler;
  std::string assembler_text;

  int          dimension;
  double       smoother_dampen;
  unsigned int smoother_steps;
  unsigned int n_steps;
  bool         output;
};

template <int dim>
class LaplaceProblem
{
  typedef LA::MPI::SparseMatrix       MatrixType;
  typedef LA::MPI::Vector             VectorType;
  typedef LA::MPI::PreconditionAMG    PreconditionAMG;
  typedef LA::MPI::PreconditionJacobi PreconditionJacobi;

public:
  LaplaceProblem(const Settings &settings);
  void run();

private:
  void setup_system();
  void setup_multigrid();
  void assemble_system();
  void assemble_multigrid();
  void solve();
  void estimate();
  void refine_grid();
  void output_results(const unsigned int cycle);

  Settings settings;

  MPI_Comm           mpi_communicator;
  ConditionalOStream pcout;

  parallel::distributed::Triangulation<dim> triangulation;
  const MappingQ1<dim>                      mapping;
  FE_Q<dim>                                 fe;

  DoFHandler<dim> dof_handler;

  IndexSet                  locally_owned_set;
  IndexSet                  locally_relevant_set;
  AffineConstraints<double> constraints;

  MatrixType     system_matrix;
  VectorType     solution;
  VectorType     right_hand_side;
  Vector<double> estimate_vector;

  MGLevelObject<MatrixType> mg_matrix;
  MGLevelObject<MatrixType> mg_interface_in;
  MGConstrainedDoFs         mg_constrained_dofs;

  TimerOutput computing_timer;
};


template <int dim>
LaplaceProblem<dim>::LaplaceProblem(const Settings &settings)
  : settings(settings)
  , mpi_communicator(MPI_COMM_WORLD)
  , pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
  , triangulation(mpi_communicator,
                  Triangulation<dim>::limit_level_difference_at_vertices,
                  (settings.assembler == Settings::amg) ?
                    parallel::distributed::Triangulation<dim>::default_setting :
                    parallel::distributed::Triangulation<
                      dim>::construct_multigrid_hierarchy)
  , mapping()
  , fe(2)
  , dof_handler(triangulation)
  , computing_timer(pcout, TimerOutput::summary, TimerOutput::wall_times)
{
  GridGenerator::hyper_L(triangulation, -1, 1, /*colorize*/ false);
  triangulation.refine_global(1);
}


bool Settings::try_parse(const std::string &prm_filename)
{
  ParameterHandler prm;
  prm.declare_entry("dim", "2", Patterns::Integer(), "The problem dimension.");
  prm.declare_entry("n_steps",
                    "20",
                    Patterns::Integer(0),
                    "Number of adaptive refinement steps.");
  prm.declare_entry("smoother dampen",
                    "1.0",
                    Patterns::Double(0.0),
                    "Dampen factor for the smoother.");
  prm.declare_entry("smoother steps",
                    "2",
                    Patterns::Integer(1),
                    "Number of smoother steps.");
  prm.declare_entry("assembler",
                    "GMG",
                    Patterns::Selection("GMG|AMG"),
                    "Switch between GMG and AMG.");
  prm.declare_entry("output",
                    "false",
                    Patterns::Bool(),
                    "Output graphical results.");

  try
    {
      prm.parse_input(prm_filename);
    }
  catch (...)
    {
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        prm.print_parameters(std::cout, ParameterHandler::Text);
      return false;
    }

  if (prm.get("assembler") == "GMG")
    this->assembler = gmg;
  else if (prm.get("assembler") == "AMG")
    this->assembler = amg;
  else
    AssertThrow(false, ExcNotImplemented());
  this->assembler_text = prm.get("assembler");

  this->dimension       = prm.get_integer("dim");
  this->n_steps         = prm.get_integer("n_steps");
  this->smoother_dampen = prm.get_double("smoother dampen");
  this->smoother_steps  = prm.get_integer("smoother steps");

  this->output = prm.get_bool("output");

  return true;
}


template <int dim>
void LaplaceProblem<dim>::setup_system()
{
  TimerOutput::Scope timing(computing_timer, "Setup");

  dof_handler.distribute_dofs(fe);

  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_set);
  locally_owned_set = dof_handler.locally_owned_dofs();

  solution.reinit(locally_owned_set, mpi_communicator);
  right_hand_side.reinit(locally_owned_set, mpi_communicator);
  constraints.reinit(locally_relevant_set);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  VectorTools::interpolate_boundary_values(
    mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraints);
  constraints.close();

#ifdef USE_PETSC_LA
  DynamicSparsityPattern dsp(locally_relevant_set);
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);

  SparsityTools::distribute_sparsity_pattern(dsp,
                                             locally_owned_set,
                                             mpi_communicator,
                                             locally_relevant_set);

  system_matrix.reinit(locally_owned_set,
                       locally_owned_set,
                       dsp,
                       mpi_communicator);
#else
  TrilinosWrappers::SparsityPattern dsp(locally_owned_set,
                                        locally_owned_set,
                                        locally_relevant_set,
                                        MPI_COMM_WORLD);
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);
  dsp.compress();
  system_matrix.reinit(dsp);
#endif
}


template <int dim>
void LaplaceProblem<dim>::setup_multigrid()
{
  TimerOutput::Scope timing(computing_timer, "Setup multigrid");

  dof_handler.distribute_mg_dofs();

  mg_constrained_dofs.clear();
  mg_constrained_dofs.initialize(dof_handler);
  std::set<types::boundary_id> bset;
  bset.insert(0);
  mg_constrained_dofs.make_zero_boundary_constraints(dof_handler, bset);

  const unsigned int n_levels = triangulation.n_global_levels();
  mg_matrix.resize(0, n_levels - 1);
  mg_matrix.clear_elements();
  mg_interface_in.resize(0, n_levels - 1);
  mg_interface_in.clear_elements();

  for (unsigned int level = 0; level < n_levels; ++level)
    {
      IndexSet dofset;
      DoFTools::extract_locally_relevant_level_dofs(dof_handler, level, dofset);

      {
#ifdef USE_PETSC_LA
        DynamicSparsityPattern dsp(dofset);
        MGTools::make_sparsity_pattern(dof_handler, dsp, level);
        dsp.compress();
        SparsityTools::distribute_sparsity_pattern(
          dsp,
          dof_handler.locally_owned_mg_dofs(level),
          mpi_communicator,
          dofset);

        mg_matrix[level].reinit(dof_handler.locally_owned_mg_dofs(level),
                                dof_handler.locally_owned_mg_dofs(level),
                                dsp,
                                mpi_communicator);
#else
        TrilinosWrappers::SparsityPattern dsp(
          dof_handler.locally_owned_mg_dofs(level),
          dof_handler.locally_owned_mg_dofs(level),
          dofset,
          mpi_communicator);
        MGTools::make_sparsity_pattern(dof_handler, dsp, level);

        dsp.compress();
        mg_matrix[level].reinit(dsp);
#endif
      }

      {
#ifdef USE_PETSC_LA
        DynamicSparsityPattern dsp(dofset);
        MGTools::make_interface_sparsity_pattern(dof_handler,
                                                 mg_constrained_dofs,
                                                 dsp,
                                                 level);
        dsp.compress();
        SparsityTools::distribute_sparsity_pattern(
          dsp,
          dof_handler.locally_owned_mg_dofs(level),
          mpi_communicator,
          dofset);

        mg_interface_in[level].reinit(dof_handler.locally_owned_mg_dofs(level),
                                      dof_handler.locally_owned_mg_dofs(level),
                                      dsp,
                                      mpi_communicator);
#else
        TrilinosWrappers::SparsityPattern dsp(
          dof_handler.locally_owned_mg_dofs(level),
          dof_handler.locally_owned_mg_dofs(level),
          dofset,
          mpi_communicator);

        MGTools::make_interface_sparsity_pattern(dof_handler,
                                                 mg_constrained_dofs,
                                                 dsp,
                                                 level);
        dsp.compress();
        mg_interface_in[level].reinit(dsp);
#endif
      }
    }
}


template <int dim>
void LaplaceProblem<dim>::assemble_system()
{
  TimerOutput::Scope timing(computing_timer, "Assemble");

  const QGauss<dim> quadrature_formula(fe.degree + 1);

  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  const Coefficient<dim> coefficient;
  std::vector<double>    coefficient_values(n_q_points);
  RightHandSide<dim>     rhs;
  std::vector<double>    rhs_values(n_q_points);

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        cell_matrix = 0;
        cell_rhs    = 0;

        fe_values.reinit(cell);

        coefficient.value_list(fe_values.get_quadrature_points(),
                               coefficient_values);
        average(coefficient_values);
        const double coefficient_value = coefficient_values[0];

        rhs.value_list(fe_values.get_quadrature_points(), rhs_values);

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                cell_matrix(i, j) +=
                  (coefficient_value * fe_values.shape_grad(i, q_point) *
                   fe_values.shape_grad(j, q_point)) *
                  fe_values.JxW(q_point);

              cell_rhs(i) +=
                (fe_values.shape_value(i, q_point) * rhs_values[q_point]) *
                fe_values.JxW(q_point);
            }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix,
                                               cell_rhs,
                                               local_dof_indices,
                                               system_matrix,
                                               right_hand_side);
      }

  system_matrix.compress(VectorOperation::add);
  right_hand_side.compress(VectorOperation::add);
}


template <int dim>
void LaplaceProblem<dim>::assemble_multigrid()
{
  TimerOutput::Scope timing(computing_timer, "Assemble multigrid");

  QGauss<dim> quadrature_formula(1 + fe.degree);

  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  const Coefficient<dim> coefficient;
  std::vector<double>    coefficient_values(n_q_points);

  std::vector<AffineConstraints<double>> boundary_constraints(
    triangulation.n_global_levels());
  for (unsigned int level = 0; level < triangulation.n_global_levels(); ++level)
    {
      IndexSet dofset;
      DoFTools::extract_locally_relevant_level_dofs(dof_handler, level, dofset);
      boundary_constraints[level].reinit(dofset);
      boundary_constraints[level].add_lines(
        mg_constrained_dofs.get_refinement_edge_indices(level));
      boundary_constraints[level].add_lines(
        mg_constrained_dofs.get_boundary_indices(level));

      boundary_constraints[level].close();
    }

  for (const auto &cell : dof_handler.cell_iterators())
    if (cell->level_subdomain_id() == triangulation.locally_owned_subdomain())
      {
        cell_matrix = 0;
        fe_values.reinit(cell);

        coefficient.value_list(fe_values.get_quadrature_points(),
                               coefficient_values);
        average(coefficient_values);
        const double coefficient_value = coefficient_values[0];

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              cell_matrix(i, j) +=
                (coefficient_value * fe_values.shape_grad(i, q_point) *
                 fe_values.shape_grad(j, q_point)) *
                fe_values.JxW(q_point);

        cell->get_mg_dof_indices(local_dof_indices);

        boundary_constraints[cell->level()].distribute_local_to_global(
          cell_matrix, local_dof_indices, mg_matrix[cell->level()]);

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            if (mg_constrained_dofs.is_interface_matrix_entry(
                  cell->level(), local_dof_indices[i], local_dof_indices[j]))
              mg_interface_in[cell->level()].add(local_dof_indices[i],
                                                 local_dof_indices[j],
                                                 cell_matrix(i, j));
      }

  for (unsigned int i = 0; i < triangulation.n_global_levels(); ++i)
    {
      mg_matrix[i].compress(VectorOperation::add);
      mg_interface_in[i].compress(VectorOperation::add);
    }
}


template <int dim>
void LaplaceProblem<dim>::solve()
{
  TimerOutput::Scope timing(computing_timer, "Solve");

  SolverControl solver_control(1000, 1.e-10 * right_hand_side.l2_norm());
  solver_control.enable_history_data();
  SolverCG<VectorType> solver(solver_control);

  solution = 0.;

  if (settings.assembler == Settings::amg)
    {
      computing_timer.enter_subsection("Solve: AMG preconditioner setup");

      PreconditionAMG                 prec;
      PreconditionAMG::AdditionalData Amg_data;

#ifdef USE_PETSC_LA
      Amg_data.symmetric_operator = true;
#else
      Amg_data.elliptic              = true;
      Amg_data.smoother_type         = "Jacobi";
      Amg_data.higher_order_elements = true;
      Amg_data.smoother_sweeps       = settings.smoother_steps;
      Amg_data.aggregation_threshold = 0.02;
#endif

      Amg_data.output_details = false;

      prec.initialize(system_matrix, Amg_data);
      computing_timer.leave_subsection("Solve: AMG preconditioner setup");

      {
        TimerOutput::Scope timing(computing_timer, "Solve: 1 AMG vcycle");
        prec.vmult(solution, right_hand_side);
      }
      solution = 0.;

      {
        TimerOutput::Scope timing(computing_timer, "Solve: CG");
        solver.solve(system_matrix, solution, right_hand_side, prec);
      }
      constraints.distribute(solution);
    }
  else
    {
      computing_timer.enter_subsection("Solve: GMG preconditioner setup");

      MGTransferPrebuilt<VectorType> mg_transfer(mg_constrained_dofs);
      mg_transfer.build(dof_handler);

      MatrixType &         coarse_matrix = mg_matrix[0];
      SolverControl        coarse_solver_control(1000, 1e-12, false, false);
      SolverCG<VectorType> coarse_solver(coarse_solver_control);
      PreconditionIdentity identity;

      MGCoarseGridIterativeSolver<VectorType,
                                  SolverCG<VectorType>,
                                  MatrixType,
                                  PreconditionIdentity>
        coarse_grid_solver(coarse_solver, coarse_matrix, identity);

      typedef LA::MPI::PreconditionJacobi                      Smoother;
      MGSmootherPrecondition<MatrixType, Smoother, VectorType> smoother;

#ifdef USE_PETSC_LA
      smoother.initialize(mg_matrix);
      Assert(
        settings.smoother_dampen == 1.0,
        ExcNotImplemented(
          "PETSc's PreconditionJacobi has no support for a damping parameter."));
#else
      smoother.initialize(mg_matrix, settings.smoother_dampen);
#endif

      smoother.set_steps(settings.smoother_steps);

      mg::Matrix<VectorType> mg_m(mg_matrix);
      mg::Matrix<VectorType> mg_in(mg_interface_in);
      mg::Matrix<VectorType> mg_out(mg_interface_in);

      Multigrid<VectorType> mg(
        mg_m, coarse_grid_solver, mg_transfer, smoother, smoother);
      mg.set_edge_matrices(mg_out, mg_in);


      PreconditionMG<dim, VectorType, MGTransferPrebuilt<VectorType>>
        preconditioner(dof_handler, mg, mg_transfer);

      computing_timer.leave_subsection("Solve: GMG preconditioner setup");

      {
        TimerOutput::Scope timing(computing_timer, "Solve: 1 GMG vcycle");
        preconditioner.vmult(solution, right_hand_side);
      }
      solution = 0.;

      {
        TimerOutput::Scope timing(computing_timer, "Solve: CG");
        solver.solve(system_matrix, solution, right_hand_side, preconditioner);
      }

      constraints.distribute(solution);
    }

  double rate = solver_control.final_reduction();
  {
    double r0 = right_hand_side.l2_norm();
    double rn = solver_control.last_value();
    rate      = 1.0 / solver_control.last_step() * log(r0 / rn) / log(10);
  }

  pcout << "   CG iterations: " << solver_control.last_step()
        << ", iters: " << 10.0 / rate << ", rate: " << rate << std::endl;
}



template <int dim>
struct ScratchData
{
  ScratchData(const Mapping<dim> &      mapping,
              const FiniteElement<dim> &fe,
              const unsigned int        quadrature_degree,
              const UpdateFlags         update_flags,
              const UpdateFlags         interface_update_flags)
    : fe_values(mapping, fe, QGauss<dim>(quadrature_degree), update_flags)
    , fe_interface_values(mapping,
                          fe,
                          QGauss<dim - 1>(quadrature_degree),
                          interface_update_flags)
  {}


  ScratchData(const ScratchData<dim> &scratch_data)
    : fe_values(scratch_data.fe_values.get_mapping(),
                scratch_data.fe_values.get_fe(),
                scratch_data.fe_values.get_quadrature(),
                scratch_data.fe_values.get_update_flags())
    , fe_interface_values(scratch_data.fe_values.get_mapping(),
                          scratch_data.fe_values.get_fe(),
                          scratch_data.fe_interface_values.get_quadrature(),
                          scratch_data.fe_interface_values.get_update_flags())
  {}

  FEValues<dim>          fe_values;
  FEInterfaceValues<dim> fe_interface_values;
};



struct CopyData
{
  CopyData()
    : cell_index(numbers::invalid_unsigned_int)
    , value(0.)
  {}

  CopyData(const CopyData &) = default;

  struct FaceData
  {
    unsigned int cell_indices[2];
    double       values[2];
  };

  unsigned int          cell_index;
  double                value;
  std::vector<FaceData> face_data;
};


template <int dim>
void LaplaceProblem<dim>::estimate()
{
  TimerOutput::Scope timing(computing_timer, "Estimate");

  VectorType temp_solution;
  temp_solution.reinit(locally_owned_set,
                       locally_relevant_set,
                       mpi_communicator);
  temp_solution = solution;

  Coefficient<dim> coefficient;

  estimate_vector.reinit(triangulation.n_active_cells());

  using Iterator = typename DoFHandler<dim>::active_cell_iterator;

  auto cell_worker = [&](const Iterator &  cell,
                         ScratchData<dim> &scratch_data,
                         CopyData &        copy_data) {
    // assemble cell residual $h^2 \| f + \epsilon \triangle u \|_K^2$

    FEValues<dim> &fe_values = scratch_data.fe_values;
    fe_values.reinit(cell);

    RightHandSide<dim> rhs;
    const double       rhs_value = rhs.value(cell->center());

    const double nu = coefficient.value(cell->center());

    std::vector<Tensor<2, dim>> hessians(fe_values.n_quadrature_points);
    fe_values.get_function_hessians(temp_solution, hessians);

    copy_data.cell_index = cell->active_cell_index();

    double value = 0.;
    for (unsigned k = 0; k < fe_values.n_quadrature_points; ++k)
      {
        const double res =
          cell->diameter() * (rhs_value + nu * trace(hessians[k]));
        value += res * res * fe_values.JxW(k);
      }
    copy_data.value = std::sqrt(value);
  };

  auto face_worker = [&](const Iterator &    cell,
                         const unsigned int &f,
                         const unsigned int &sf,
                         const Iterator &    ncell,
                         const unsigned int &nf,
                         const unsigned int &nsf,
                         ScratchData<dim> &  scratch_data,
                         CopyData &          copy_data) {
    // face term $\sum_F h_F \| [ \epsilon \nabla u \cdot n ] \|_F^2$

    FEInterfaceValues<dim> &fe_interface_values =
      scratch_data.fe_interface_values;
    fe_interface_values.reinit(cell, f, sf, ncell, nf, nsf);

    copy_data.face_data.emplace_back();
    CopyData::FaceData &copy_data_face = copy_data.face_data.back();

    copy_data_face.cell_indices[0] = cell->active_cell_index();
    copy_data_face.cell_indices[1] = ncell->active_cell_index();

    const double nu1 = coefficient.value(cell->center());
    const double nu2 = coefficient.value(ncell->center());
    const double h   = cell->face(f)->measure(); // TODO: FEIV.measure

    std::vector<Tensor<1, dim>> grad_u[2];

    for (unsigned int i = 0; i < 2; ++i)
      {
        grad_u[i].resize(fe_interface_values.n_quadrature_points);
        fe_interface_values.get_fe_face_values(i).get_function_gradients(
          temp_solution, grad_u[i]);
      }

    double value = 0.;

    for (unsigned int qpoint = 0;
         qpoint < fe_interface_values.n_quadrature_points;
         ++qpoint)
      {
        const double jump =
          nu1 * grad_u[0][qpoint] * fe_interface_values.normal(qpoint) -
          nu2 * grad_u[1][qpoint] * fe_interface_values.normal(qpoint);

        value += h * jump * jump * fe_interface_values.JxW(qpoint);
      }

    copy_data_face.values[0] = 0.5 * std::sqrt(value);
    copy_data_face.values[1] = copy_data_face.values[0];
  };

  auto copier = [&](const CopyData &copy_data) {
    if (copy_data.cell_index != numbers::invalid_unsigned_int)
      estimate_vector[copy_data.cell_index] += copy_data.value;

    for (auto &cdf : copy_data.face_data)
      for (unsigned int j = 0; j < 2; ++j)
        estimate_vector[cdf.cell_indices[j]] += cdf.values[j];
  };

  const unsigned int n_gauss_points = dof_handler.get_fe().degree + 1;
  ScratchData<dim>   scratch_data(mapping,
                                fe,
                                n_gauss_points,
                                update_hessians | update_quadrature_points |
                                  update_JxW_values,
                                update_values | update_gradients |
                                  update_JxW_values | update_normal_vectors);

  CopyData copy_data;
  MeshWorker::mesh_loop(dof_handler.begin_active(),
                        dof_handler.end(),
                        cell_worker,
                        copier,
                        scratch_data,
                        copy_data,
                        MeshWorker::assemble_own_cells |
                          MeshWorker::assemble_ghost_faces_both |
                          MeshWorker::assemble_own_interior_faces_once,
                        nullptr /*boundary_worker*/,
                        face_worker);
}



template <int dim>
void LaplaceProblem<dim>::refine_grid()
{
  TimerOutput::Scope timing(computing_timer, "Refine grid");

  const double refinement_fraction = 1. / (std::pow(2.0, dim) - 1.);
  parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
    triangulation, estimate_vector, refinement_fraction, 0.0);

  triangulation.execute_coarsening_and_refinement();
}



template <int dim>
void LaplaceProblem<dim>::output_results(const unsigned int cycle)
{
  TimerOutput::Scope timing(computing_timer, "Output results");

  DataOut<dim> data_out;

  VectorType temp_solution;
  temp_solution.reinit(locally_owned_set,
                       locally_relevant_set,
                       mpi_communicator);
  temp_solution = solution;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(temp_solution, "solution");
  Vector<float> subdomain(triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain();
  data_out.add_data_vector(subdomain, "subdomain");

  Vector<float> level(triangulation.n_active_cells());
  for (const auto &cell : triangulation.active_cell_iterators())
    level(cell->active_cell_index()) = cell->level();
  data_out.add_data_vector(level, "level");

  if (estimate_vector.size() > 0)
    data_out.add_data_vector(estimate_vector, "estimator");

  data_out.build_patches(0);

  const std::string filename =
    ("solution-" + Utilities::int_to_string(cycle, 5) + "." +
     Utilities::int_to_string(triangulation.locally_owned_subdomain(), 4) +
     ".vtu");
  std::ofstream output(filename.c_str());
  data_out.write_vtu(output);

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::vector<std::string> filenames;
      for (unsigned int i = 0;
           i < Utilities::MPI::n_mpi_processes(mpi_communicator);
           ++i)
        filenames.push_back(std::string("solution-") +
                            Utilities::int_to_string(cycle, 5) + "." +
                            Utilities::int_to_string(i, 4) + ".vtu");
      const std::string pvtu_master_filename =
        ("solution-" + Utilities::int_to_string(cycle, 5) + ".pvtu");
      std::ofstream pvtu_master(pvtu_master_filename.c_str());
      data_out.write_pvtu_record(pvtu_master, filenames);

      const std::string visit_master_filename =
        ("solution-" + Utilities::int_to_string(cycle, 5) + ".visit");
      std::ofstream visit_master(visit_master_filename.c_str());
      DataOutBase::write_visit_record(visit_master, filenames);
    }
}


template <int dim>
void LaplaceProblem<dim>::run()
{
  for (unsigned int cycle = 0; cycle < settings.n_steps; ++cycle)
    {
      pcout << "Cycle " << cycle << ':' << std::endl;
      if (cycle > 0)
        refine_grid();

      pcout << "   Number of active cells:       "
            << triangulation.n_global_active_cells();
      if (settings.assembler == Settings::gmg)
        pcout << " (" << triangulation.n_global_levels() << " global levels)"
              << std::endl
              << "   Workload imbalance:           "
              << MGTools::workload_imbalance(triangulation);
      pcout << std::endl;

      setup_system();
      if (settings.assembler == Settings::gmg)
        setup_multigrid();

      pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs();
      if (settings.assembler != Settings::amg)
        {
          pcout << " (by level: ";
          for (unsigned int level = 0; level < triangulation.n_global_levels();
               ++level)
            pcout << dof_handler.n_dofs(level)
                  << (level == triangulation.n_global_levels() - 1 ? ")" :
                                                                     ", ");
        }
      pcout << std::endl;

      assemble_system();
      if (settings.assembler == Settings::gmg)
        assemble_multigrid();

      solve();
      estimate();

      if (settings.output)
        output_results(cycle);

      computing_timer.print_summary();
      computing_timer.reset();
    }
}


int main(int argc, char *argv[])
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  using namespace dealii;

  Settings settings;
  if (!settings.try_parse((argc > 1) ? (argv[1]) : ""))
    return 0;

  try
    {
      if (settings.dimension == 2)
        {
          LaplaceProblem<2> test(settings);
          test.run();
        }
      else if (settings.dimension == 3)
        {
          LaplaceProblem<3> test(settings);
          test.run();
        }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 2);
      return 1;
    }

  return 0;
}

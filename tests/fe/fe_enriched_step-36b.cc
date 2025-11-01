// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test FE_Enriched in real-life application on eigenvalue problem similar
// to Step-36. That involves assembly (shape values and gradients) and
// error estimator (Kelly - > face gradients) and MPI run.
//
// same as fe_enriched_step-36, but using FESystem only. Note that the results
// of KellyErrorEstimator won't be the same, thus this tests does only a single
// refinement cycle where the refinement path is the same.

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_enriched.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/slepc_solver.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

const unsigned int dim = 3;


/**
 * Coulomb potential
 */
template <int dim>
class PotentialFunction : public dealii::Function<dim>
{
public:
  PotentialFunction()
    : dealii::Function<dim>(1)
  {}

  virtual double
  value(const dealii::Point<dim> &point,
        const unsigned int        component = 0) const;
};

template <int dim>
double
PotentialFunction<dim>::value(const dealii::Point<dim> &p,
                              const unsigned int) const
{
  return -1.0 / std::sqrt(p.square());
}

template <int dim>
class EnrichmentFunction : public Function<dim>
{
public:
  EnrichmentFunction(const Point<dim> &origin,
                     const double     &Z,
                     const double     &radius)
    : Function<dim>(1)
    , origin(origin)
    , Z(Z)
    , radius(radius)
  {}

  virtual double
  value(const Point<dim> &point, const unsigned int component = 0) const
  {
    Tensor<1, dim> dist = point - origin;
    const double   r    = dist.norm();
    return std::exp(-Z * r);
  }

  bool
  is_enriched(const Point<dim> &point) const
  {
    if (origin.distance(point) < radius)
      return true;
    else
      return false;
  }

  virtual Tensor<1, dim>
  gradient(const Point<dim> &p, const unsigned int component = 0) const
  {
    Tensor<1, dim> dist = p - origin;
    const double   r    = dist.norm();
    dist /= r;
    return -Z * std::exp(-Z * r) * dist;
  }

private:
  /**
   * origin
   */
  const Point<dim> origin;

  /**
   * charge
   */
  const double Z;

  /**
   * enrichment radius
   */
  const double radius;
};

namespace Step36
{
  /**
   * Main class
   */
  template <int dim>
  class EigenvalueProblem
  {
  public:
    EigenvalueProblem();
    virtual ~EigenvalueProblem();
    void
    run();

  private:
    bool
    cell_is_pou(const typename DoFHandler<dim>::cell_iterator &cell) const;

    std::pair<unsigned int, unsigned int>
    setup_system();
    void
    assemble_system();
    std::pair<unsigned int, double>
    solve();
    void
    constrain_pou_dofs();
    void
    estimate_error();
    void
    refine_grid();
    void
    output_results(const unsigned int cycle) const;

    Triangulation<dim>    triangulation;
    DoFHandler<dim>       dof_handler;
    hp::FECollection<dim> fe_collection;
    hp::QCollection<dim>  q_collection;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    std::vector<PETScWrappers::MPI::Vector> eigenfunctions;
    std::vector<PETScWrappers::MPI::Vector> eigenfunctions_locally_relevant;
    std::vector<PetscScalar>                eigenvalues;
    PETScWrappers::MPI::SparseMatrix        stiffness_matrix, mass_matrix;

    AffineConstraints<double> constraints;

    MPI_Comm           mpi_communicator;
    const unsigned int n_mpi_processes;
    const unsigned int this_mpi_process;

    ConditionalOStream pcout;

    const unsigned int number_of_eigenvalues;

    PotentialFunction<dim> potential;

    EnrichmentFunction<dim> enrichment;

    const FEValuesExtractors::Scalar fe_extractor;
    const unsigned int               fe_group;
    const unsigned int               fe_fe_index;
    const unsigned int               fe_material_id;
    const FEValuesExtractors::Scalar pou_extractor;
    const unsigned int               pou_group;
    const unsigned int               pou_fe_index;
    const unsigned int               pou_material_id;

    std::vector<Vector<float>> vec_estimated_error_per_cell;
    Vector<float>              estimated_error_per_cell;
  };

  template <int dim>
  EigenvalueProblem<dim>::EigenvalueProblem()
    : dof_handler(triangulation)
    , mpi_communicator(MPI_COMM_WORLD)
    , n_mpi_processes(dealii::Utilities::MPI::n_mpi_processes(mpi_communicator))
    , this_mpi_process(
        dealii::Utilities::MPI::this_mpi_process(mpi_communicator))
    , pcout(std::cout, (this_mpi_process == 0))
    , number_of_eigenvalues(1)
    , enrichment(Point<dim>(),
                 /*Z*/ 1.0,
                 /*radius*/ 2.5)
    , // radius is set such that 8 cells are marked as enriched
    fe_extractor(/*dofs start at...*/ 0)
    , fe_group(/*in FE*/ 0)
    , fe_fe_index(0)
    , fe_material_id(0)
    , pou_extractor(/*dofs start at (scalar fields!)*/ 1)
    , pou_group(/*FE*/ 1)
    , pou_fe_index(1)
    , pou_material_id(1)
  {
    dealii::GridGenerator::hyper_cube(triangulation, -10, 10);
    triangulation.refine_global(2); // 64 cells

    for (auto cell = triangulation.begin_active(); cell != triangulation.end();
         ++cell)
      if (std::sqrt(cell->center().square()) < 5.0)
        cell->set_refine_flag();

    triangulation.execute_coarsening_and_refinement(); // 120 cells

    q_collection.push_back(QGauss<dim>(4));
    q_collection.push_back(QGauss<dim>(10));

    // usual elements (active_fe_index ==0):
    fe_collection.push_back(
      FESystem<dim>(FE_Q<dim>(2), 1, FE_Nothing<dim>(), 1));

    // enriched elements (active_fe_index==1):
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(2), 1, FE_Q<dim>(1), 1));

    // assign FE index in the constructor so that FE/FE+POU is determined
    // on the coarsest mesh to avoid issues like
    // +---------+----+----+
    // |         | fe |    |
    // |    fe   +----+----+
    // |         | pou|    |
    // +---------+----+----+
    // see discussion in Step46.
    for (typename DoFHandler<dim>::cell_iterator cell =
           dof_handler.begin_active();
         cell != dof_handler.end();
         ++cell)
      if (enrichment.is_enriched(cell->center()))
        cell->set_material_id(pou_material_id);
      else
        cell->set_material_id(fe_material_id);
  }

  template <int dim>
  std::pair<unsigned int, unsigned int>
  EigenvalueProblem<dim>::setup_system()
  {
    for (typename DoFHandler<dim>::active_cell_iterator cell =
           dof_handler.begin_active();
         cell != dof_handler.end();
         ++cell)
      {
        if (cell->material_id() == fe_material_id)
          cell->set_active_fe_index(fe_fe_index);
        else if (cell->material_id() == pou_material_id)
          cell->set_active_fe_index(pou_fe_index);
        else
          DEAL_II_NOT_IMPLEMENTED();
      }

    GridTools::partition_triangulation(n_mpi_processes, triangulation);
    dof_handler.distribute_dofs(fe_collection);
    DoFRenumbering::subdomain_wise(dof_handler);
    std::vector<IndexSet> locally_owned_dofs_per_processor =
      DoFTools::locally_owned_dofs_per_subdomain(dof_handler);
    locally_owned_dofs = locally_owned_dofs_per_processor[this_mpi_process];
    locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);

    constraints.clear();
    constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(2),
                                             constraints);
    constrain_pou_dofs();
    constraints.close();

    // Initialise the stiffness and mass matrices
    DynamicSparsityPattern csp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler,
                                    csp,
                                    constraints,
                                    /* keep constrained dofs */ false);

    std::vector<types::global_dof_index> n_locally_owned_dofs(n_mpi_processes);
    for (unsigned int i = 0; i < n_mpi_processes; ++i)
      n_locally_owned_dofs[i] =
        locally_owned_dofs_per_processor[i].n_elements();

    SparsityTools::distribute_sparsity_pattern(csp,
                                               n_locally_owned_dofs,
                                               mpi_communicator,
                                               locally_relevant_dofs);

    stiffness_matrix.reinit(locally_owned_dofs,
                            locally_owned_dofs,
                            csp,
                            mpi_communicator);

    mass_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       csp,
                       mpi_communicator);

    // reinit vectors
    eigenfunctions.resize(number_of_eigenvalues);
    eigenfunctions_locally_relevant.resize(number_of_eigenvalues);
    vec_estimated_error_per_cell.resize(number_of_eigenvalues);
    for (unsigned int i = 0; i < eigenfunctions.size(); ++i)
      {
        eigenfunctions[i].reinit(locally_owned_dofs,
                                 mpi_communicator); // without ghost dofs
        eigenfunctions_locally_relevant[i].reinit(locally_owned_dofs,
                                                  locally_relevant_dofs,
                                                  mpi_communicator);

        vec_estimated_error_per_cell[i].reinit(triangulation.n_active_cells());
      }

    eigenvalues.resize(eigenfunctions.size());

    estimated_error_per_cell.reinit(triangulation.n_active_cells());

    unsigned int n_pou_cells = 0, n_fem_cells = 0;

    for (typename DoFHandler<dim>::cell_iterator cell =
           dof_handler.begin_active();
         cell != dof_handler.end();
         ++cell)
      if (cell_is_pou(cell))
        n_pou_cells++;
      else
        n_fem_cells++;

    return std::make_pair(n_fem_cells, n_pou_cells);
  }

  template <int dim>
  bool
  EigenvalueProblem<dim>::cell_is_pou(
    const typename DoFHandler<dim>::cell_iterator &cell) const
  {
    return cell->material_id() == pou_material_id;
  }

  template <int dim>
  void
  EigenvalueProblem<dim>::constrain_pou_dofs()
  {
    std::vector<types::global_dof_index> local_face_dof_indices(
      fe_collection[pou_fe_index].dofs_per_face);
    for (typename DoFHandler<dim>::active_cell_iterator cell =
           dof_handler.begin_active();
         cell != dof_handler.end();
         ++cell)
      if (cell_is_pou(cell))
        for (const unsigned int f : GeometryInfo<dim>::face_indices())
          if (!cell->at_boundary(f))
            {
              bool face_is_on_interface = false;
              if ((cell->neighbor(f)->has_children() ==
                   false /* => it is active */) &&
                  (!cell_is_pou(cell->neighbor(f))))
                face_is_on_interface = true;
              else if (cell->neighbor(f)->has_children() == true)
                {
                  // The neighbor does
                  // have
                  // children. See if
                  // any of the cells
                  // on the other
                  // side are vanilla FEM
                  for (unsigned int sf = 0; sf < cell->face(f)->n_children();
                       ++sf)
                    if (!cell_is_pou(cell->neighbor_child_on_subface(f, sf)))
                      {
                        face_is_on_interface = true;
                        break;
                      }
                }
              // add constraints
              if (face_is_on_interface)
                {
                  cell->face(f)->get_dof_indices(local_face_dof_indices,
                                                 pou_fe_index);
                  for (unsigned int i = 0; i < local_face_dof_indices.size();
                       ++i)
                    if (fe_collection[pou_fe_index]
                          .face_system_to_base_index(i)
                          .first.first == pou_group)
                      // if
                      // (fe_collection[1].face_system_to_component_index(i).first
                      // /*component*/ > 0)
                      if (constraints.is_constrained(
                            local_face_dof_indices[i]) == false)
                        constraints.constrain_dof_to_zero(
                          local_face_dof_indices[i]);
                }
            }
  }

  template <int dim>
  void
  EigenvalueProblem<dim>::assemble_system()
  {
    stiffness_matrix = 0;
    mass_matrix      = 0;

    dealii::FullMatrix<double>           cell_stiffness_matrix;
    dealii::FullMatrix<double>           cell_mass_matrix;
    std::vector<types::global_dof_index> local_dof_indices;
    std::vector<double>                  potential_values;
    std::vector<double>                  enrichment_values;
    std::vector<Tensor<1, dim>>          enrichment_gradients;

    hp::FEValues<dim> fe_values_hp(fe_collection,
                                   q_collection,
                                   update_values | update_gradients |
                                     update_quadrature_points |
                                     update_JxW_values);


    typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell != endc; ++cell)
      if (cell->subdomain_id() == this_mpi_process)
        {
          fe_values_hp.reinit(cell);
          const FEValues<dim> &fe_values = fe_values_hp.get_present_fe_values();

          const unsigned int &dofs_per_cell = cell->get_fe().dofs_per_cell;
          const unsigned int &n_q_points    = fe_values.n_quadrature_points;

          potential_values.resize(n_q_points);
          enrichment_values.resize(n_q_points);
          enrichment_gradients.resize(n_q_points);

          local_dof_indices.resize(dofs_per_cell);
          cell_stiffness_matrix.reinit(dofs_per_cell, dofs_per_cell);
          cell_mass_matrix.reinit(dofs_per_cell, dofs_per_cell);

          cell_stiffness_matrix = 0;
          cell_mass_matrix      = 0;

          potential.value_list(fe_values.get_quadrature_points(),
                               potential_values);

          if (cell->active_fe_index() == 0) // plain FE
            {
              for (const auto q_point : fe_values.quadrature_point_indices())
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  for (unsigned int j = i; j < dofs_per_cell; ++j)
                    {
                      cell_stiffness_matrix(i, j) +=
                        (fe_values[fe_extractor].gradient(i, q_point) *
                           fe_values[fe_extractor].gradient(j, q_point) * 0.5 +
                         potential_values[q_point] *
                           fe_values[fe_extractor].value(i, q_point) *
                           fe_values[fe_extractor].value(j, q_point)) *
                        fe_values.JxW(q_point);

                      cell_mass_matrix(i, j) +=
                        (fe_values[fe_extractor].value(i, q_point) *
                         fe_values[fe_extractor].value(j, q_point)) *
                        fe_values.JxW(q_point);
                    }
            }
          else // POUFEM
            {
              Assert(cell->active_fe_index() == 1, ExcInternalError());

              enrichment.value_list(fe_values.get_quadrature_points(),
                                    enrichment_values);
              enrichment.gradient_list(fe_values.get_quadrature_points(),
                                       enrichment_gradients);

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  const unsigned int i_group =
                    cell->get_fe().system_to_base_index(i).first.first;
                  for (unsigned int j = i; j < dofs_per_cell; ++j)
                    {
                      const unsigned int j_group =
                        cell->get_fe().system_to_base_index(j).first.first;

                      if (i_group == fe_group && j_group == fe_group) // fe - fe
                        {
                          for (const auto q_point :
                               fe_values.quadrature_point_indices())
                            {
                              cell_stiffness_matrix(i, j) +=
                                (fe_values[fe_extractor].gradient(i, q_point) *
                                   fe_values[fe_extractor].gradient(j,
                                                                    q_point) *
                                   0.5 +
                                 potential_values[q_point] *
                                   fe_values[fe_extractor].value(i, q_point) *
                                   fe_values[fe_extractor].value(j, q_point)) *
                                fe_values.JxW(q_point);

                              cell_mass_matrix(i, j) +=
                                (fe_values[fe_extractor].value(i, q_point) *
                                 fe_values[fe_extractor].value(j, q_point)) *
                                fe_values.JxW(q_point);
                            }
                        }
                      else if (i_group == fe_group &&
                               j_group == pou_group) // fe - pou
                        {
                          for (const auto q_point :
                               fe_values.quadrature_point_indices())
                            {
                              cell_stiffness_matrix(i, j) +=
                                (fe_values[fe_extractor].gradient(i, q_point) *
                                   (fe_values[pou_extractor].gradient(j,
                                                                      q_point) *
                                      enrichment_values[q_point] +
                                    fe_values[pou_extractor].value(j, q_point) *
                                      enrichment_gradients[q_point]) *
                                   0.5 +
                                 potential_values[q_point] *
                                   fe_values[fe_extractor].value(i, q_point) *
                                   fe_values[pou_extractor].value(j, q_point) *
                                   enrichment_values[q_point]) *
                                fe_values.JxW(q_point);

                              cell_mass_matrix(i, j) +=
                                (fe_values[fe_extractor].value(i, q_point) *
                                 fe_values[pou_extractor].value(j, q_point) *
                                 enrichment_values[q_point]) *
                                fe_values.JxW(q_point);
                            }
                        }
                      else if (i_group == pou_group &&
                               j_group == fe_group) // pou - fe
                        {
                          for (const auto q_point :
                               fe_values.quadrature_point_indices())
                            {
                              cell_stiffness_matrix(i, j) +=
                                ((fe_values[pou_extractor].gradient(i,
                                                                    q_point) *
                                    enrichment_values[q_point] +
                                  fe_values[pou_extractor].value(i, q_point) *
                                    enrichment_gradients[q_point]) *
                                   fe_values[fe_extractor].gradient(j,
                                                                    q_point) *
                                   0.5 +
                                 potential_values[q_point] *
                                   fe_values[pou_extractor].value(i, q_point) *
                                   enrichment_values[q_point] *
                                   fe_values[fe_extractor].value(j, q_point)) *
                                fe_values.JxW(q_point);

                              cell_mass_matrix(i, j) +=
                                (fe_values[pou_extractor].value(i, q_point) *
                                 enrichment_values[q_point] *
                                 fe_values[fe_extractor].value(j, q_point)) *
                                fe_values.JxW(q_point);
                            }
                        }
                      else // pou - pou
                        {
                          Assert(i_group == pou_group && j_group == pou_group,
                                 ExcInternalError());

                          for (const auto q_point :
                               fe_values.quadrature_point_indices())
                            {
                              cell_stiffness_matrix(i, j) +=
                                ((fe_values[pou_extractor].gradient(i,
                                                                    q_point) *
                                    enrichment_values[q_point] +
                                  fe_values[pou_extractor].value(i, q_point) *
                                    enrichment_gradients[q_point]) *
                                   (fe_values[pou_extractor].gradient(j,
                                                                      q_point) *
                                      enrichment_values[q_point] +
                                    fe_values[pou_extractor].value(j, q_point) *
                                      enrichment_gradients[q_point]) *
                                   0.5 +
                                 potential_values[q_point] *
                                   fe_values[pou_extractor].value(i, q_point) *
                                   enrichment_values[q_point] *
                                   fe_values[pou_extractor].value(j, q_point) *
                                   enrichment_values[q_point]) *
                                fe_values.JxW(q_point);

                              cell_mass_matrix(i, j) +=
                                (fe_values[pou_extractor].value(i, q_point) *
                                 enrichment_values[q_point] *
                                 fe_values[pou_extractor].value(j, q_point) *
                                 enrichment_values[q_point]) *
                                fe_values.JxW(q_point);
                            }
                        }
                    }
                }
            }

          // exploit symmetry
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = i; j < dofs_per_cell; ++j)
              {
                cell_stiffness_matrix(j, i) = cell_stiffness_matrix(i, j);
                cell_mass_matrix(j, i)      = cell_mass_matrix(i, j);
              }

          cell->get_dof_indices(local_dof_indices);

          constraints.distribute_local_to_global(cell_stiffness_matrix,
                                                 local_dof_indices,
                                                 stiffness_matrix);
          constraints.distribute_local_to_global(cell_mass_matrix,
                                                 local_dof_indices,
                                                 mass_matrix);
        }

    stiffness_matrix.compress(dealii::VectorOperation::add);
    mass_matrix.compress(dealii::VectorOperation::add);
  }

  template <int dim>
  std::pair<unsigned int, double>
  EigenvalueProblem<dim>::solve()
  {
    dealii::SolverControl solver_control(dof_handler.n_dofs(),
                                         1e-9,
                                         false,
                                         false);
    dealii::SLEPcWrappers::SolverKrylovSchur eigensolver(solver_control,
                                                         mpi_communicator);

    eigensolver.set_which_eigenpairs(EPS_SMALLEST_REAL);
    eigensolver.set_problem_type(EPS_GHEP);

    eigensolver.solve(stiffness_matrix,
                      mass_matrix,
                      eigenvalues,
                      eigenfunctions,
                      eigenfunctions.size());

    for (unsigned int i = 0; i < eigenfunctions.size(); ++i)
      {
        constraints.distribute(eigenfunctions[i]);
        eigenfunctions_locally_relevant[i] = eigenfunctions[i];
      }

    return std::make_pair(solver_control.last_step(),
                          solver_control.last_value());
  }

  template <int dim>
  EigenvalueProblem<dim>::~EigenvalueProblem()
  {
    dof_handler.clear();
  }

  template <int dim>
  void
  EigenvalueProblem<dim>::estimate_error()
  {
    {
      std::vector<const ReadVector<PetscScalar> *> sol(number_of_eigenvalues);
      std::vector<Vector<float> *>                 error(number_of_eigenvalues);

      for (unsigned int i = 0; i < number_of_eigenvalues; ++i)
        {
          sol[i]   = &eigenfunctions_locally_relevant[i];
          error[i] = &vec_estimated_error_per_cell[i];
        }

      hp::QCollection<dim - 1> face_quadrature_formula;
      face_quadrature_formula.push_back(dealii::QGauss<dim - 1>(3));
      face_quadrature_formula.push_back(dealii::QGauss<dim - 1>(3));

      ArrayView<const ReadVector<PetscScalar> *> sol_view =
        make_array_view(sol);
      ArrayView<Vector<float> *> error_view = make_array_view(error);
      KellyErrorEstimator<dim>::estimate(
        dof_handler,
        face_quadrature_formula,
        std::map<types::boundary_id, const Function<dim> *>(),
        sol_view,
        error_view);
    }

    // sum up for a global:
    for (unsigned int c = 0; c < estimated_error_per_cell.size(); ++c)
      {
        double er = 0.0;
        for (unsigned int i = 0; i < number_of_eigenvalues; ++i)
          er += vec_estimated_error_per_cell[i][c] *
                vec_estimated_error_per_cell[i][c];

        estimated_error_per_cell[c] = sqrt(er);
      }
  }

  template <int dim>
  void
  EigenvalueProblem<dim>::refine_grid()
  {
    const double threshold = 0.9 * estimated_error_per_cell.linfty_norm();
    GridRefinement::refine(triangulation, estimated_error_per_cell, threshold);

    triangulation.prepare_coarsening_and_refinement();
    triangulation.execute_coarsening_and_refinement();
  }

  template <int dim>
  class Postprocessor : public DataPostprocessorScalar<dim>
  {
  public:
    Postprocessor(const EnrichmentFunction<dim> &enrichment);

    virtual void
    compute_derived_quantities_vector(
      const std::vector<Vector<double>>              &solution_values,
      const std::vector<std::vector<Tensor<1, dim>>> &solution_gradients,
      const std::vector<std::vector<Tensor<2, dim>>> &solution_hessians,
      const std::vector<Point<dim>>                  &normals,
      const std::vector<Point<dim>>                  &evaluation_points,
      std::vector<Vector<double>> &computed_quantities) const;

  private:
    const EnrichmentFunction<dim> &enrichment;
  };

  template <int dim>
  Postprocessor<dim>::Postprocessor(const EnrichmentFunction<dim> &enrichment)
    : DataPostprocessorScalar<dim>("total_solution",
                                   update_values | update_quadrature_points)
    , enrichment(enrichment)
  {}

  template <int dim>
  void
  Postprocessor<dim>::compute_derived_quantities_vector(
    const std::vector<Vector<double>> &solution_values,
    const std::vector<std::vector<Tensor<1, dim>>> & /*solution_gradients*/,
    const std::vector<std::vector<Tensor<2, dim>>> & /*solution_hessians*/,
    const std::vector<Point<dim>> & /*normals*/,
    const std::vector<Point<dim>> &evaluation_points,
    std::vector<Vector<double>>   &computed_quantities) const
  {
    const unsigned int n_quadrature_points = solution_values.size();
    Assert(computed_quantities.size() == n_quadrature_points,
           ExcInternalError());
    for (unsigned int q = 0; q < n_quadrature_points; ++q)
      {
        Assert(solution_values[q].size() == 2,
               ExcDimensionMismatch(solution_values[q].size(),
                                    2)); // FESystem with 2 components

        computed_quantities[q](0) =
          (solution_values[q](0) +
           solution_values[q](1) *
             enrichment.value(
               evaluation_points[q])); // for FE_Nothing solution_values[q](1)
                                       // will be zero
      }
  }

  template <int dim>
  void
  EigenvalueProblem<dim>::output_results(const unsigned int cycle) const
  {
    dealii::Vector<float> fe_index(triangulation.n_active_cells());
    {
      typename dealii::DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
      for (unsigned int index = 0; cell != endc; ++cell, ++index)
        fe_index(index) = cell->active_fe_index();
    }

    Assert(cycle < 10, ExcNotImplemented());
    if (this_mpi_process == 0)
      {
        std::string filename = "solution-";
        filename += ('0' + cycle);
        filename += ".vtk";
        std::ofstream output(filename);

        Postprocessor<dim> postprocessor(
          enrichment); // has to live until the DataOut object is destroyed;
                       // objects are destroyed in reverse order of declaration
        DataOut<dim> data_out;
        data_out.attach_dof_handler(dof_handler);
        data_out.add_data_vector(eigenfunctions_locally_relevant[0],
                                 "solution");
        data_out.add_data_vector(eigenfunctions_locally_relevant[0],
                                 postprocessor);
        data_out.build_patches(6);
        data_out.write_vtk(output);
      }

    // second output without making mesh finer
    if (this_mpi_process == 0)
      {
        std::string filename = "mesh-";
        filename += ('0' + cycle);
        filename += ".vtk";
        std::ofstream output(filename);

        DataOut<dim> data_out;
        data_out.attach_dof_handler(dof_handler);
        data_out.add_data_vector(fe_index, "fe_index");
        data_out.add_data_vector(estimated_error_per_cell, "estimated_error");
        data_out.build_patches();
        data_out.write_vtk(output);
      }

    // scalar data for plotting
    // output scalar data (eigenvalues, energies, ndofs, etc)
    if (this_mpi_process == 0)
      {
        const std::string scalar_fname = "scalar-data.txt";

        std::ofstream output(scalar_fname,
                             std::ios::out |
                               (cycle == 0 ? std::ios::trunc : std::ios::app));

        std::streamsize max_precision = std::numeric_limits<double>::digits10;

        std::string sep("  ");
        output << cycle << sep << std::setprecision(max_precision)
               << triangulation.n_active_cells() << sep << dof_handler.n_dofs()
               << sep << std::scientific;

        for (unsigned int i = 0; i < eigenvalues.size(); ++i)
          output << eigenvalues[i] << sep;

        output << std::endl;
        output.close();
      } // end scope
  }


  template <int dim>
  void
  EigenvalueProblem<dim>::run()
  {
    for (unsigned int cycle = 0; cycle < 2; ++cycle)
      {
        pcout << "Cycle " << cycle << std::endl;
        const std::pair<unsigned int, unsigned int> n_cells = setup_system();

        pcout << "   Number of active cells:       "
              << triangulation.n_active_cells() << std::endl
              << "     FE / POUFE :                " << n_cells.first << " / "
              << n_cells.second << std::endl
              << "   Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;

        assemble_system();

        const std::pair<unsigned int, double> res = solve();
        AssertThrow(res.second < 5e-8, ExcInternalError());

        estimate_error();
        // output_results(cycle);
        refine_grid();

        pcout << std::endl;
        for (unsigned int i = 0; i < eigenvalues.size(); ++i)
          pcout << "      Eigenvalue " << i << " : " << eigenvalues[i]
                << std::endl;
      }
  }
} // namespace Step36

int
main(int argc, char **argv)
{
  try
    {
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc,
                                                                  argv,
                                                                  1);
      {
        Step36::EigenvalueProblem<dim> step36;
        dealii::PETScWrappers::set_option_value("-eps_target", "-1.0");
        dealii::PETScWrappers::set_option_value("-st_type", "sinvert");
        dealii::PETScWrappers::set_option_value("-st_ksp_type", "cg");
        dealii::PETScWrappers::set_option_value("-st_pc_type", "jacobi");
        step36.run();
      }
    }
  catch (const std::exception &exc)
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
      return 1;
    };
}

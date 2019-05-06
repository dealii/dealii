/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2018 - 2019 by the deal.II authors
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
 * Authors: Thomas Clevenger, Clemson University
 *          Timo Heister, University of Utah
 */

// @note: This is work in progress and will be an example for block smoothers
// in geometric multigrid.
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/std_cxx14/memory.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/path_search.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/relaxation_block.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_matrix.h>

#include <fstream>
#include <iostream>
#include <random>


namespace Step63
{
  using namespace dealii;


  template <int dim>
  struct ScratchData
  {
    ScratchData(const FiniteElement<dim> &fe,
                const unsigned int        quadrature_degree)
      : fe_values(fe,
                  QGauss<dim>(quadrature_degree),
                  update_values | update_gradients | update_hessians |
                    update_quadrature_points | update_JxW_values)
    {}

    ScratchData(const ScratchData<dim> &scratch_data)
      : fe_values(scratch_data.fe_values.get_fe(),
                  scratch_data.fe_values.get_quadrature(),
                  update_values | update_gradients | update_hessians |
                    update_quadrature_points | update_JxW_values)
    {}

    FEValues<dim> fe_values;
  };



  struct CopyData
  {
    CopyData() = default;

    unsigned int level;
    unsigned int dofs_per_cell;

    FullMatrix<double>                   cell_matrix;
    Vector<double>                       cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;
  };


  struct Settings
  {
    enum DoFRenumberingStrategy
    {
      none,
      downstream,
      upstream,
      random
    };

    void get_parameters(const std::string &prm_filename);

    unsigned int           fe_degree;
    std::string            smoother_type;
    DoFRenumberingStrategy dof_renumbering;
    bool                   with_streamline_diffusion;
    bool                   output;
  };

  void Settings::get_parameters(const std::string &prm_filename)
  {
    ParameterHandler prm;

    prm.declare_entry("Fe degree",
                      "1",
                      Patterns::Integer(0),
                      "Finite Element degree");
    prm.declare_entry("Smoother type",
                      "block SOR",
                      Patterns::Selection("sor|Jacobi|block SOR|block Jacobi"),
                      "Smoother Type: SOR|Jacobi|block SOR|block Jacobi");
    prm.declare_entry("DoF renumbering",
                      "downstream",
                      Patterns::Selection("none|downstream|upstream|random"),
                      "DoF renumbering: none|downstream|upstream|random");
    prm.declare_entry("With streamline diffusion",
                      "true",
                      Patterns::Bool(),
                      "With streamline diffusion: true|false");
    prm.declare_entry("Output",
                      "true",
                      Patterns::Bool(),
                      "Generate graphical output: true|false");

    if (prm_filename.empty())
      {
        prm.print_parameters(std::cout, ParameterHandler::Text);
        AssertThrow(
          false, ExcMessage("please pass a .prm file as the first argument!"));
      }

    prm.parse_input(prm_filename);

    fe_degree     = prm.get_integer("Fe degree");
    smoother_type = prm.get("Smoother type");

    std::string renumbering = prm.get("DoF renumbering");
    if (renumbering == "none")
      dof_renumbering = DoFRenumberingStrategy::none;
    else if (renumbering == "downstream")
      dof_renumbering = DoFRenumberingStrategy::downstream;
    else if (renumbering == "upstream")
      dof_renumbering = DoFRenumberingStrategy::upstream;
    else if (renumbering == "random")
      dof_renumbering = DoFRenumberingStrategy::random;

    with_streamline_diffusion = prm.get_bool("With streamline diffusion");
    output                    = prm.get_bool("Output");
  }


  // Functions for creating permutation of cells for output and Block
  // smoothers
  template <int dim>
  std::vector<unsigned int>
  create_downstream_cell_ordering(const DoFHandler<dim> &dof,
                                  const Tensor<1, dim>   direction,
                                  const unsigned int     level)
  {
    std::vector<typename DoFHandler<dim>::level_cell_iterator> ordered_cells;
    ordered_cells.reserve(dof.get_triangulation().n_cells(level));
    const DoFRenumbering::
      CompareDownstream<typename DoFHandler<dim>::level_cell_iterator, dim>
        comparator(direction);

    for (const auto &cell : dof.cell_iterators_on_level(level))
      ordered_cells.push_back(cell);

    std::sort(ordered_cells.begin(), ordered_cells.end(), comparator);

    std::vector<unsigned> ordered_indices;
    ordered_indices.reserve(dof.get_triangulation().n_cells(level));

    for (const auto &cell : ordered_cells)
      ordered_indices.push_back(cell->index());

    return ordered_indices;
  }

  template <int dim>
  std::vector<unsigned int>
  create_downstream_cell_ordering(const DoFHandler<dim> &dof,
                                  const Tensor<1, dim>   direction)
  {
    std::vector<typename DoFHandler<dim>::active_cell_iterator> ordered_cells;
    ordered_cells.reserve(dof.get_triangulation().n_active_cells());
    const DoFRenumbering::
      CompareDownstream<typename DoFHandler<dim>::active_cell_iterator, dim>
        comparator(direction);

    for (const auto &cell : dof.active_cell_iterators())
      ordered_cells.push_back(cell);

    std::sort(ordered_cells.begin(), ordered_cells.end(), comparator);

    std::vector<unsigned int> ordered_indices;
    ordered_indices.reserve(dof.get_triangulation().n_active_cells());

    for (const auto &cell : ordered_cells)
      ordered_indices.push_back(cell->index());

    return ordered_indices;
  }



  template <int dim>
  std::vector<unsigned int>
  create_random_cell_ordering(const DoFHandler<dim> &dof,
                              const unsigned int     level)
  {
    const unsigned int n_cells = dof.get_triangulation().n_cells(level);

    std::vector<unsigned int> ordered_cells;
    ordered_cells.reserve(n_cells);

    for (const auto &cell : dof.cell_iterators_on_level(level))
      ordered_cells.push_back(cell->index());

    // shuffle the elements
    std::mt19937 random_number_generator;
    for (unsigned int i = 1; i < n_cells; ++i)
      {
        // get a random number between 0 and i (inclusive)
        const unsigned int j =
          std::uniform_int_distribution<>(0, i)(random_number_generator);

        // if possible, swap the elements
        if (i != j)
          std::swap(ordered_cells[i], ordered_cells[j]);
      }

    return ordered_cells;
  }

  template <int dim>
  std::vector<unsigned int>
  create_random_cell_ordering(const DoFHandler<dim> &dof)
  {
    const unsigned int n_cells = dof.get_triangulation().n_active_cells();

    std::vector<unsigned int> ordered_cells;
    ordered_cells.reserve(n_cells);

    for (const auto &cell : dof.active_cell_iterators())
      ordered_cells.push_back(cell->index());

    // shuffle the elements
    std::mt19937 random_number_generator;
    for (unsigned int i = 1; i < n_cells; ++i)
      {
        // get a random number between 0 and i (inclusive)
        const unsigned int j =
          std::uniform_int_distribution<>(0, i)(random_number_generator);

        // if possible, swap the elements
        if (i != j)
          std::swap(ordered_cells[i], ordered_cells[j]);
      }

    return ordered_cells;
  }



  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide()
      : Function<dim>()
    {}

    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;

    virtual void value_list(const std::vector<Point<dim>> &points,
                            std::vector<double> &          values,
                            const unsigned int component = 0) const override;
  };

  template <int dim>
  double RightHandSide<dim>::value(const Point<dim> &,
                                   const unsigned int component) const
  {
    Assert(component == 0, ExcIndexRange(component, 0, 1));
    (void)component;

    return 0.0;
  }



  template <int dim>
  void RightHandSide<dim>::value_list(const std::vector<Point<dim>> &points,
                                      std::vector<double> &          values,
                                      const unsigned int component) const
  {
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));

    for (unsigned int i = 0; i < points.size(); ++i)
      values[i] = RightHandSide<dim>::value(points[i], component);
  }



  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    BoundaryValues()
      : Function<dim>()
    {}

    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;

    virtual void value_list(const std::vector<Point<dim>> &points,
                            std::vector<double> &          values,
                            const unsigned int component = 0) const override;
  };


  template <int dim>
  double BoundaryValues<dim>::value(const Point<dim> & p,
                                    const unsigned int component) const
  {
    Assert(component == 0, ExcIndexRange(component, 0, 1));
    (void)component;

    if (std::fabs(p[0] - 1) < 1e-8                     // x == 1
        || (std::fabs(p[1] + 1) < 1e-8 && p[0] >= 0.5) // y == -1, x > 0.5
    )
      {
        return 1.0;
      }
    else
      {
        return 0.0;
      }
  }



  template <int dim>
  void BoundaryValues<dim>::value_list(const std::vector<Point<dim>> &points,
                                       std::vector<double> &          values,
                                       const unsigned int component) const
  {
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));

    for (unsigned int i = 0; i < points.size(); ++i)
      values[i] = BoundaryValues<dim>::value(points[i], component);
  }

  template <int dim>
  double compute_stabilization_delta(const double         hk,
                                     const double         eps,
                                     const Tensor<1, dim> dir,
                                     const double         pk)
  {
    // Value defined in 'On discontinuity–capturing methods for
    // convection–diffusion equations'
    const double Peclet = dir.norm() * hk / (2.0 * eps * pk);
    const double coth =
      (1.0 + std::exp(-2.0 * Peclet)) / (1.0 - std::exp(-2.0 * Peclet));

    return hk / (2.0 * dir.norm() * pk) * (coth - 1.0 / Peclet);
  }


  template <int dim>
  class AdvectionProblem
  {
  public:
    AdvectionProblem(Settings settings);
    void run();

  private:
    void setup_system();

    template <class IteratorType>
    void assemble_cell(const IteratorType &cell,
                       ScratchData<dim> &  scratch_data,
                       CopyData &          copy_data);
    void assemble_system_and_multigrid();

    void setup_smoother();

    void solve();
    void refine_grid();
    void output_results(const unsigned int cycle) const;

    Triangulation<dim> triangulation;
    DoFHandler<dim>    dof_handler;

    const FE_Q<dim>     fe;
    const MappingQ<dim> mapping;

    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;

    MGLevelObject<SparsityPattern> mg_sparsity_patterns;
    MGLevelObject<SparsityPattern> mg_interface_sparsity_patterns;

    MGLevelObject<SparseMatrix<double>> mg_matrices;
    MGLevelObject<SparseMatrix<double>> mg_interface_in;
    MGLevelObject<SparseMatrix<double>> mg_interface_out;

    mg::Matrix<Vector<double>> mg_matrix;
    mg::Matrix<Vector<double>> mg_interface_matrix_in;
    mg::Matrix<Vector<double>> mg_interface_matrix_out;

    using SmootherType =
      RelaxationBlock<SparseMatrix<double>, double, Vector<double>>;
    using SmootherAdditionalDataType = SmootherType::AdditionalData;
    std::unique_ptr<MGSmoother<Vector<double>>> mg_smoother;

    MGLevelObject<SmootherAdditionalDataType> smoother_data;

    MGConstrainedDoFs mg_constrained_dofs;

    const double   epsilon;
    Tensor<1, dim> advection_direction;

    const Settings settings;
  };



  template <int dim>
  AdvectionProblem<dim>::AdvectionProblem(Settings settings)
    : triangulation(Triangulation<dim>::limit_level_difference_at_vertices)
    , dof_handler(triangulation)
    , fe(settings.fe_degree)
    , mapping(settings.fe_degree)
    , epsilon(0.005)
    , settings(settings)
  {
    advection_direction[0] = -std::sin(numbers::PI / 6.0);
    if (dim > 1)
      advection_direction[1] = std::cos(numbers::PI / 6.0);
    if (dim > 2)
      advection_direction[2] = std::sin(numbers::PI / 6.0);
  }



  template <int dim>
  void AdvectionProblem<dim>::setup_system()
  {
    const unsigned int n_levels = triangulation.n_levels();

    dof_handler.distribute_dofs(fe);

    // We could renumber the active dofs with DoFRenumbering::downstream()
    // here, but the smoothers only act on multigrid levels and as such, this
    // wouldn't matter. Instead, we will renumber the DoFs on each multigrid
    // level below.

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    VectorTools::interpolate_boundary_values(
      mapping, dof_handler, 0, BoundaryValues<dim>(), constraints);
    VectorTools::interpolate_boundary_values(
      mapping, dof_handler, 1, BoundaryValues<dim>(), constraints);
    constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ false);

    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);


    // Setup GMG DoFs
    dof_handler.distribute_mg_dofs();

    // Renumber DoFs on each level in downstream or upstream direction if
    // needed. This is only necessary for point smoothers (SOR and Jacobi) as
    // the block smoothers operate on cells (see create_smoother()):
    if (settings.smoother_type == "SOR" || settings.smoother_type == "Jacobi")
      {
        if (settings.dof_renumbering ==
              Settings::DoFRenumberingStrategy::downstream ||
            settings.dof_renumbering ==
              Settings::DoFRenumberingStrategy::upstream)
          {
            const Tensor<1, dim> direction =
              (settings.dof_renumbering ==
                   Settings::DoFRenumberingStrategy::upstream ?
                 -1.0 :
                 1.0) *
              advection_direction;

            for (unsigned int level = 0; level < n_levels; ++level)
              DoFRenumbering::downstream(dof_handler,
                                         level,
                                         direction,
                                         /*dof_wise_renumbering = */ true);
          }
        else if (settings.dof_renumbering ==
                 Settings::DoFRenumberingStrategy::random)
          {
            for (unsigned int level = 0; level < n_levels; ++level)
              DoFRenumbering::random(dof_handler, level);
          }
        else
          Assert(false, ExcNotImplemented());
      }

    mg_constrained_dofs.clear();
    mg_constrained_dofs.initialize(dof_handler);

    mg_constrained_dofs.make_zero_boundary_constraints(dof_handler, {0, 1});

    mg_matrices.resize(0, n_levels - 1);
    mg_matrices.clear_elements();
    mg_interface_in.resize(0, n_levels - 1);
    mg_interface_in.clear_elements();
    mg_interface_out.resize(0, n_levels - 1);
    mg_interface_out.clear_elements();
    mg_sparsity_patterns.resize(0, n_levels - 1);
    mg_interface_sparsity_patterns.resize(0, n_levels - 1);

    for (unsigned int level = 0; level < n_levels; ++level)
      {
        {
          DynamicSparsityPattern dsp(dof_handler.n_dofs(level),
                                     dof_handler.n_dofs(level));
          MGTools::make_sparsity_pattern(dof_handler, dsp, level);
          mg_sparsity_patterns[level].copy_from(dsp);
          mg_matrices[level].reinit(mg_sparsity_patterns[level]);
        }
        {
          DynamicSparsityPattern dsp(dof_handler.n_dofs(level),
                                     dof_handler.n_dofs(level));
          MGTools::make_interface_sparsity_pattern(dof_handler,
                                                   mg_constrained_dofs,
                                                   dsp,
                                                   level);
          mg_interface_sparsity_patterns[level].copy_from(dsp);

          // We need both interface in and out matrices since our problem is not
          // symmetric
          mg_interface_in[level].reinit(mg_interface_sparsity_patterns[level]);
          mg_interface_out[level].reinit(mg_interface_sparsity_patterns[level]);
        }
      }
  }


  template <int dim>
  template <class IteratorType>
  void AdvectionProblem<dim>::assemble_cell(const IteratorType &cell,
                                            ScratchData<dim> &  scratch_data,
                                            CopyData &          copy_data)
  {
    copy_data.level = cell->level();

    const unsigned int dofs_per_cell =
      scratch_data.fe_values.get_fe().dofs_per_cell;
    copy_data.dofs_per_cell = dofs_per_cell;
    copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);

    const unsigned int n_q_points =
      scratch_data.fe_values.get_quadrature().size();

    if (cell->is_level_cell() == false)
      copy_data.cell_rhs.reinit(dofs_per_cell);

    copy_data.local_dof_indices.resize(dofs_per_cell);
    cell->get_active_or_mg_dof_indices(copy_data.local_dof_indices);

    scratch_data.fe_values.reinit(cell);

    const RightHandSide<dim> right_hand_side;
    std::vector<double>      rhs_values(n_q_points);

    right_hand_side.value_list(scratch_data.fe_values.get_quadrature_points(),
                               rhs_values);

    const double delta = settings.with_streamline_diffusion ?
                           compute_stabilization_delta(cell->diameter(),
                                                       epsilon,
                                                       advection_direction,
                                                       settings.fe_degree) :
                           0.0;

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              copy_data.cell_matrix(i, j) +=
                (epsilon * scratch_data.fe_values.shape_grad(j, q_point) *
                 scratch_data.fe_values.shape_grad(i, q_point) *
                 scratch_data.fe_values.JxW(q_point)) +
                ((advection_direction *
                  scratch_data.fe_values.shape_grad(j, q_point)) *
                 scratch_data.fe_values.shape_value(i, q_point)) *
                  scratch_data.fe_values.JxW(q_point);

              if (settings.with_streamline_diffusion)
                copy_data.cell_matrix(i, j) +=
                  delta *
                    (advection_direction *
                     scratch_data.fe_values.shape_grad(j, q_point)) *
                    (advection_direction *
                     scratch_data.fe_values.shape_grad(i, q_point)) *
                    scratch_data.fe_values.JxW(q_point) -
                  delta * epsilon *
                    trace(scratch_data.fe_values.shape_hessian(j, q_point)) *
                    (advection_direction *
                     scratch_data.fe_values.shape_grad(i, q_point)) *
                    scratch_data.fe_values.JxW(q_point);
            }
          if (cell->is_level_cell() == false)
            {
              copy_data.cell_rhs(i) +=
                scratch_data.fe_values.shape_value(i, q_point) *
                rhs_values[q_point] * scratch_data.fe_values.JxW(q_point);
              if (settings.with_streamline_diffusion)
                copy_data.cell_rhs(i) +=
                  delta * rhs_values[q_point] * advection_direction *
                  scratch_data.fe_values.shape_grad(i, q_point) *
                  scratch_data.fe_values.JxW(q_point);
            }
        }
  }


  template <int dim>
  void AdvectionProblem<dim>::assemble_system_and_multigrid()
  {
    auto cell_worker_active =
      [&](const decltype(dof_handler.begin_active()) &cell,
          ScratchData<dim> &                          scratch_data,
          CopyData &                                  copy_data) {
        this->assemble_cell(cell, scratch_data, copy_data);
      };



    auto copier_active = [&](const CopyData &copy_data) {
      constraints.distribute_local_to_global(copy_data.cell_matrix,
                                             copy_data.cell_rhs,
                                             copy_data.local_dof_indices,
                                             system_matrix,
                                             system_rhs);
    };


    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker_active,
                          copier_active,
                          ScratchData<dim>(fe, fe.degree + 1),
                          CopyData(),
                          MeshWorker::assemble_own_cells);

    std::vector<AffineConstraints<double>> boundary_constraints(
      triangulation.n_global_levels());
    for (unsigned int level = 0; level < triangulation.n_global_levels();
         ++level)
      {
        IndexSet locally_owned_level_dof_indices;
        DoFTools::extract_locally_relevant_level_dofs(
          dof_handler, level, locally_owned_level_dof_indices);
        boundary_constraints[level].reinit(locally_owned_level_dof_indices);
        boundary_constraints[level].add_lines(
          mg_constrained_dofs.get_refinement_edge_indices(level));
        boundary_constraints[level].add_lines(
          mg_constrained_dofs.get_boundary_indices(level));
        boundary_constraints[level].close();
      }

    auto cell_worker_mg = [&](const decltype(dof_handler.begin_mg()) &cell,
                              ScratchData<dim> &scratch_data,
                              CopyData &        copy_data) {
      this->assemble_cell(cell, scratch_data, copy_data);
    };

    auto copier_mg = [&](const CopyData &copy_data) {
      boundary_constraints[copy_data.level].distribute_local_to_global(
        copy_data.cell_matrix,
        copy_data.local_dof_indices,
        mg_matrices[copy_data.level]);

      // If (i,j) is an interface_out dof pair, then (j,i) is an interface_in
      // dof pair. Note: for interface_in, we load the transpose of the
      // interface entries, i.e., the entry for dof pair (j,i) is stored in
      // interface_in(i,j).
      for (unsigned int i = 0; i < copy_data.dofs_per_cell; ++i)
        for (unsigned int j = 0; j < copy_data.dofs_per_cell; ++j)
          if (mg_constrained_dofs.is_interface_matrix_entry(
                copy_data.level,
                copy_data.local_dof_indices[i],
                copy_data.local_dof_indices[j]))
            {
              mg_interface_out[copy_data.level].add(
                copy_data.local_dof_indices[i],
                copy_data.local_dof_indices[j],
                copy_data.cell_matrix(i, j));
              mg_interface_in[copy_data.level].add(
                copy_data.local_dof_indices[i],
                copy_data.local_dof_indices[j],
                copy_data.cell_matrix(j, i));
            }
    };

    MeshWorker::mesh_loop(dof_handler.begin_mg(),
                          dof_handler.end_mg(),
                          cell_worker_mg,
                          copier_mg,
                          ScratchData<dim>(fe, fe.degree + 1),
                          CopyData(),
                          MeshWorker::assemble_own_cells);
  }


  template <int dim>
  void AdvectionProblem<dim>::setup_smoother()
  {
    if (settings.smoother_type == "SOR")
      {
        using Smoother = PreconditionSOR<SparseMatrix<double>>;

        auto smoother =
          std_cxx14::make_unique<MGSmootherPrecondition<SparseMatrix<double>,
                                                        Smoother,
                                                        Vector<double>>>();
        smoother->initialize(mg_matrices,
                             Smoother::AdditionalData(fe.degree == 1 ? 1.0 :
                                                                       0.62));
        smoother->set_steps(2);
        mg_smoother = std::move(smoother);
      }
    else if (settings.smoother_type == "Jacobi")
      {
        using Smoother = PreconditionJacobi<SparseMatrix<double>>;
        auto smoother =
          std_cxx14::make_unique<MGSmootherPrecondition<SparseMatrix<double>,
                                                        Smoother,
                                                        Vector<double>>>();
        smoother->initialize(mg_matrices,
                             Smoother::AdditionalData(fe.degree == 1 ? 0.6667 :
                                                                       0.47));
        smoother->set_steps(4);
        mg_smoother = std::move(smoother);
      }
    else if (settings.smoother_type == "block SOR" ||
             settings.smoother_type == "block Jacobi")
      {
        smoother_data.resize(0, triangulation.n_levels() - 1);

        for (unsigned int level = 0; level < triangulation.n_levels(); ++level)
          {
            DoFTools::make_cell_patches(smoother_data[level].block_list,
                                        dof_handler,
                                        level);

            smoother_data[level].relaxation =
              (settings.smoother_type == "block SOR" ? 1.0 : 0.25);
            smoother_data[level].inversion = PreconditionBlockBase<double>::svd;

            std::vector<unsigned int> ordered_indices;
            switch (settings.dof_renumbering)
              {
                case Settings::DoFRenumberingStrategy::downstream:
                  ordered_indices =
                    create_downstream_cell_ordering(dof_handler,
                                                    advection_direction,
                                                    level);
                  break;

                case Settings::DoFRenumberingStrategy::upstream:
                  ordered_indices =
                    create_downstream_cell_ordering(dof_handler,
                                                    -1.0 * advection_direction,
                                                    level);
                  break;

                case Settings::DoFRenumberingStrategy::random:
                  ordered_indices =
                    create_random_cell_ordering(dof_handler, level);
                  break;

                case Settings::DoFRenumberingStrategy::none: // Do nothing
                  break;

                default:
                  AssertThrow(false, ExcNotImplemented());
                  break;
              }

            smoother_data[level].order =
              std::vector<std::vector<unsigned int>>(1, ordered_indices);
          }

        if (settings.smoother_type == "block SOR")
          {
            auto smoother = std_cxx14::make_unique<MGSmootherPrecondition<
              SparseMatrix<double>,
              RelaxationBlockSOR<SparseMatrix<double>, double, Vector<double>>,
              Vector<double>>>();
            smoother->initialize(mg_matrices, smoother_data);
            smoother->set_steps(1);
            mg_smoother = std::move(smoother);
          }
        else if (settings.smoother_type == "block Jacobi")
          {
            auto smoother = std_cxx14::make_unique<
              MGSmootherPrecondition<SparseMatrix<double>,
                                     RelaxationBlockJacobi<SparseMatrix<double>,
                                                           double,
                                                           Vector<double>>,
                                     Vector<double>>>();
            smoother->initialize(mg_matrices, smoother_data);
            smoother->set_steps(2);
            mg_smoother = std::move(smoother);
          }
      }
    else
      AssertThrow(false, ExcNotImplemented());
  }


  template <int dim>
  void AdvectionProblem<dim>::solve()
  {
    const unsigned int max_iters       = 200;
    const double       solve_tolerance = 1e-8 * system_rhs.l2_norm();
    SolverControl      solver_control(max_iters, solve_tolerance, true, true);
    solver_control.enable_history_data();

    using Transfer = MGTransferPrebuilt<Vector<double>>;
    Transfer mg_transfer(mg_constrained_dofs);
    mg_transfer.build_matrices(dof_handler);

    FullMatrix<double> coarse_matrix;
    coarse_matrix.copy_from(mg_matrices[0]);
    MGCoarseGridHouseholder<double, Vector<double>> coarse_grid_solver;
    coarse_grid_solver.initialize(coarse_matrix);

    setup_smoother();

    mg_matrix.initialize(mg_matrices);
    mg_interface_matrix_in.initialize(mg_interface_in);
    mg_interface_matrix_out.initialize(mg_interface_out);

    Multigrid<Vector<double>> mg(
      mg_matrix, coarse_grid_solver, mg_transfer, *mg_smoother, *mg_smoother);
    mg.set_edge_matrices(mg_interface_matrix_out, mg_interface_matrix_in);

    PreconditionMG<dim, Vector<double>, Transfer> preconditioner(dof_handler,
                                                                 mg,
                                                                 mg_transfer);

    std::cout << "     Solving with GMRES to tol " << solve_tolerance << "..."
              << std::endl;
    SolverGMRES<> solver(solver_control);

    Timer time;
    time.start();
    solver.solve(system_matrix, solution, system_rhs, preconditioner);
    time.stop();

    std::cout << "          converged in " << solver_control.last_step()
              << " iterations"
              << " in " << time.last_wall_time() << " seconds " << std::endl;

    constraints.distribute(solution);

    mg_smoother.release();
  }



  template <int dim>
  void AdvectionProblem<dim>::output_results(const unsigned int cycle) const
  {
    // Here we generate an index for each cell to visualize the ordering used
    // by the smoothers. Note that we do this only for the active cells
    // instead of the levels, where the smoothers are actually used. For the
    // point smoothers we renumber DoFs instead of cells, so this is only an
    // approximation of what happens in reality. Finally, the random ordering
    // is not the random ordering we actually use (see create_smoother() for
    // that).
    const unsigned int n_active_cells = triangulation.n_active_cells();
    Vector<double>     cell_indices(n_active_cells);
    {
      // First generate a permutation vector for the cell indices:
      std::vector<unsigned int> ordered_indices;
      switch (settings.dof_renumbering)
        {
          case Settings::DoFRenumberingStrategy::downstream:
            ordered_indices =
              create_downstream_cell_ordering(dof_handler, advection_direction);
            break;

          case Settings::DoFRenumberingStrategy::upstream:
            ordered_indices =
              create_downstream_cell_ordering(dof_handler,
                                              -1.0 * advection_direction);
            break;

          case Settings::DoFRenumberingStrategy::random:
            ordered_indices = create_random_cell_ordering(dof_handler);
            break;

          case Settings::DoFRenumberingStrategy::none:
            ordered_indices.resize(n_active_cells);
            for (unsigned int i = 0; i < n_active_cells; ++i)
              ordered_indices[i] = i;
            break;

          default:
            AssertThrow(false, ExcNotImplemented());
            break;
        }

      // Then copy the permutation in ordered_indices into an output vector:
      for (unsigned int i = 0; i < n_active_cells; ++i)
        cell_indices(ordered_indices[i]) = static_cast<double>(i);
    }

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");
    data_out.add_data_vector(cell_indices, "cell_index");
    data_out.build_patches();

    std::string filename =
      "solution-" + Utilities::int_to_string(cycle) + ".vtu";
    std::ofstream output(filename.c_str());
    data_out.write_vtu(output);
  }


  template <int dim>
  void AdvectionProblem<dim>::run()
  {
    for (unsigned int cycle = 0; cycle < (settings.fe_degree == 1 ? 7 : 5);
         ++cycle)
      {
        std::cout << "  Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          {
            GridGenerator::hyper_cube_with_cylindrical_hole(
              triangulation, 0.3, 1.0, 0.5, 1, false);
            static const SphericalManifold<dim> manifold_description(
              Point<dim>(0, 0));
            triangulation.set_manifold(1, manifold_description);
          }

        triangulation.refine_global();

        setup_system();

        std::cout << "     Number of active cells:       "
                  << triangulation.n_active_cells() << " ("
                  << triangulation.n_levels() << " levels)" << std::endl;
        std::cout << "     Number of degrees of freedom: "
                  << dof_handler.n_dofs() << std::endl;

        assemble_system_and_multigrid();

        solve();

        if (settings.output)
          output_results(cycle);

        std::cout << std::endl;
      }
  }
} // namespace Step63



int main(int argc, char *argv[])
{
  try
    {
      Step63::Settings settings;
      settings.get_parameters((argc > 1) ? (argv[1]) : "");

      Step63::AdvectionProblem<2> advection_problem_2d(settings);
      advection_problem_2d.run();
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
    }

  return 0;
}

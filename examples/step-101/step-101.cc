/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2025 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * Author: Michał Wichrowski, Heidelberg University, 2025
 */



#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/function_signed_distance.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <deal.II/non_matching/fe_values.h>
#include <deal.II/non_matching/mesh_classifier.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <vector>



namespace Step85
{
  using namespace dealii;


  template <int dim>
  class LaplaceSolver
  {
  public:
    LaplaceSolver();

    void run();

  private:
    void make_grid();

    void setup_discrete_level_set();

    void distribute_dofs();

    void initialize_matrices();

    void assemble_system();

    void solve();

    void output_results() const;

    double compute_L2_error() const;


    const unsigned int fe_degree;

    FunctionParser<dim> rhs_function;
    FunctionParser<dim> boundary_condition;
    FunctionParser<dim> analytical_solution;

    Triangulation<dim> triangulation;

    const FE_Q<dim> fe_level_set;
    DoFHandler<dim> level_set_dof_handler;
    Vector<double>  level_set;

    hp::FECollection<dim> fe_collection;
    DoFHandler<dim>       dof_handler;
    Vector<double>        solution;

    NonMatching::MeshClassifier<dim> mesh_classifier;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> stiffness_matrix;
    Vector<double>       rhs;
  };



  template <int dim>
  LaplaceSolver<dim>::LaplaceSolver()
    : fe_degree(1)
    , rhs_function("2 * cos(x) * sin(y)")
    , boundary_condition("cos(x)*sin(y)")
    , analytical_solution("cos(x)*sin(y)")
    , fe_level_set(fe_degree)
    , level_set_dof_handler(triangulation)
    , dof_handler(triangulation)
    , mesh_classifier(level_set_dof_handler, level_set)
  {}



  template <int dim>
  void LaplaceSolver<dim>::make_grid()
  {
    std::cout << "Creating background mesh" << std::endl;

    GridGenerator::hyper_cube(triangulation, -1.21, 1.21);
    triangulation.refine_global(2);
  }



  template <int dim>
  void LaplaceSolver<dim>::setup_discrete_level_set()
  {
    std::cout << "Setting up discrete level set function" << std::endl;

    level_set_dof_handler.distribute_dofs(fe_level_set);
    level_set.reinit(level_set_dof_handler.n_dofs());

    const Functions::SignedDistance::Sphere<dim> signed_distance_sphere;
    VectorTools::interpolate(level_set_dof_handler,
                             signed_distance_sphere,
                             level_set);
  }



  enum ActiveFEIndex
  {
    lagrange = 0,
    nothing  = 1
  };

  template <int dim>
  void LaplaceSolver<dim>::distribute_dofs()
  {
    std::cout << "Distributing degrees of freedom" << std::endl;

    fe_collection.push_back(FE_Q<dim>(fe_degree));
    fe_collection.push_back(FE_Nothing<dim>());

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        const NonMatching::LocationToLevelSet cell_location =
          mesh_classifier.location_to_level_set(cell);

        if (cell_location != NonMatching::LocationToLevelSet::inside)
          cell->set_active_fe_index(ActiveFEIndex::nothing);
        else
          cell->set_active_fe_index(ActiveFEIndex::lagrange);
      }

    dof_handler.distribute_dofs(fe_collection);
  }



  template <int dim>
  void LaplaceSolver<dim>::initialize_matrices()
  {
    std::cout << "Initializing matrices" << std::endl;

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());

    const AffineConstraints<double> constraints;
    const bool                      keep_constrained_dofs = true;

    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    keep_constrained_dofs);
    sparsity_pattern.copy_from(dsp);

    stiffness_matrix.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs());
    rhs.reinit(dof_handler.n_dofs());
  }



  template <int dim>
  void LaplaceSolver<dim>::assemble_system()
  {
    std::cout << "Assembling" << std::endl;

    std::vector<types::global_dof_index> level_set_dof_indices(
      fe_level_set.dofs_per_cell);
    std::vector<double> dof_values_level_set(fe_level_set.dofs_per_cell);

    const unsigned int n_dofs_per_cell = fe_collection[0].dofs_per_cell;
    FullMatrix<double> local_stiffness(n_dofs_per_cell, n_dofs_per_cell);
    Vector<double>     local_rhs(n_dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(n_dofs_per_cell);

    std::vector<Point<dim>> face_quad_points(n_dofs_per_cell);
    std::vector<Point<dim>> shifted_points(n_dofs_per_cell);
    std::vector<double>     shifted_shape_value(n_dofs_per_cell);

    const double nitsche_parameter = 5 * (fe_degree + 1) * fe_degree;

    const QGauss<dim - 1> face_quadrature(fe_degree + 1);
    const QGauss<dim>     cell_quadrature(fe_degree + 1);

    const FiniteElement<dim> &fe_lagrange =
      fe_collection[ActiveFEIndex::lagrange];

    FEValues<dim> fe_values(fe_lagrange,
                            cell_quadrature,
                            update_values | update_gradients |
                              update_JxW_values | update_quadrature_points);

    FEFaceValues<dim> surface_fe_values(fe_lagrange,
                                        face_quadrature,
                                        update_values | update_gradients |
                                          update_normal_vectors |
                                          update_JxW_values |
                                          update_quadrature_points);

    // We collect shift to output them later for debugging purposes.
    // Since we are working on a sphere, the shift can also be computed
    // analytically, so we can compare the computed shifts with the exact ones.
    std::vector<std::pair<Point<dim>, Point<dim>>> shifts;
    std::vector<std::pair<Point<dim>, Point<dim>>> exact_shifts;

    for (const auto &cell :
         dof_handler.active_cell_iterators() |
           IteratorFilters::ActiveFEIndexEqualTo(ActiveFEIndex::lagrange))
      {
        local_stiffness = 0;
        local_rhs       = 0;

        const double cell_side_length = cell->minimum_vertex_distance();
        fe_values.reinit(cell);

        for (const unsigned int q : fe_values.quadrature_point_indices())
          {
            const Point<dim> &point = fe_values.quadrature_point(q);
            for (const unsigned int i : fe_values.dof_indices())
              {
                for (const unsigned int j : fe_values.dof_indices())
                  {
                    local_stiffness(i, j) += fe_values.shape_grad(i, q) *
                                             fe_values.shape_grad(j, q) *
                                             fe_values.JxW(q);
                  }
                local_rhs(i) += rhs_function.value(point) *
                                fe_values.shape_value(i, q) * fe_values.JxW(q);
              }
          }

        for (const auto &face : GeometryInfo<dim>::face_indices())
          {
            // not at boundary
            if (cell->neighbor(face)->active_fe_index() ==
                ActiveFEIndex::lagrange)
              continue;
            surface_fe_values.reinit(cell, face);

            auto neighbour = cell->neighbor(face);
            typename DoFHandler<dim>::cell_iterator neighbours_lvl_set_cell(
              &neighbour->get_triangulation(),
              neighbour->level(),
              neighbour->index(),
              &level_set_dof_handler);

            neighbours_lvl_set_cell->get_dof_indices(level_set_dof_indices);
            level_set.extract_subvector_to(level_set_dof_indices.begin(),
                                           level_set_dof_indices.end(),
                                           dof_values_level_set.begin());



            for (const unsigned int q :
                 surface_fe_values.quadrature_point_indices())
              {
                const Point<dim> &point = surface_fe_values.quadrature_point(q);

                const Point<dim> closest_boundary_point =
                  point * (1. / point.norm());

                // Computing shifts from levels set is PR:18680

                const Point<dim> unit_shifted_point =
                  surface_fe_values.get_mapping().transform_real_to_unit_cell(
                    cell, closest_boundary_point);

                const Point<dim> real_shifted_point =
                  surface_fe_values.get_mapping().transform_unit_to_real_cell(
                    cell, unit_shifted_point);



                shifts.push_back(std::make_pair(point, real_shifted_point));
                exact_shifts.push_back(
                  std::make_pair(point, point * (1. / point.norm())));

                // Get the shifted shape values
                for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                  shifted_shape_value[j] =
                    fe_lagrange.shape_value(j, unit_shifted_point);

                // Only use mesh normal according to:
                // https://doi.org/10.1016/j.cma.2022.114885
                const Tensor<1, dim> &mesh_normal =
                  surface_fe_values.normal_vector(q);

                for (const unsigned int i : surface_fe_values.dof_indices())
                  {
                    for (const unsigned int j : surface_fe_values.dof_indices())
                      {
                        local_stiffness(i, j) +=
                          (-mesh_normal * surface_fe_values.shape_grad(i, q) *
                             shifted_shape_value[j] +
                           -mesh_normal * surface_fe_values.shape_grad(j, q) *
                             surface_fe_values.shape_value(i, q) +
                           nitsche_parameter / cell_side_length *
                             shifted_shape_value[i] * shifted_shape_value[j]) *
                          surface_fe_values.JxW(q);
                      }
                    local_rhs(i) +=
                      boundary_condition.value(closest_boundary_point) *
                      (nitsche_parameter / cell_side_length *
                         shifted_shape_value[i] -
                       mesh_normal * surface_fe_values.shape_grad(i, q)) *
                      surface_fe_values.JxW(q);
                  }
              }
          }

        cell->get_dof_indices(local_dof_indices);

        stiffness_matrix.add(local_dof_indices, local_stiffness);
        rhs.add(local_dof_indices, local_rhs);
      }
    // Here we output the shifts, not implemented (PR:18741)
    //  export_line_segments("shifts", shifts);
    //  export_line_segments("shifts_exact", exact_shifts);
  }


  template <int dim>
  void LaplaceSolver<dim>::solve()
  {
    std::cout << "Solving system" << std::endl;

    SparseDirectUMFPACK solver_direct;
    solver_direct.initialize(stiffness_matrix);
    solver_direct.vmult(solution, rhs);
  }



  template <int dim>
  void LaplaceSolver<dim>::output_results() const
  {
    std::cout << "Writing vtu file" << std::endl;

    DataOut<dim> data_out;
    data_out.add_data_vector(dof_handler, solution, "solution");
    data_out.add_data_vector(level_set_dof_handler, level_set, "level_set");

    data_out.set_cell_selection(
      [this](const typename Triangulation<dim>::cell_iterator &cell) {
        return cell->is_active() &&
               mesh_classifier.location_to_level_set(cell) !=
                 NonMatching::LocationToLevelSet::outside;
      });

    data_out.build_patches();
    std::ofstream output("sbm.vtu");
    data_out.write_vtu(output);
  }



  template <int dim>
  class AnalyticalSolution : public Function<dim>
  {
  public:
    double value(const Point<dim>  &point,
                 const unsigned int component = 0) const override;
  };



  template <int dim>
  double AnalyticalSolution<dim>::value(const Point<dim>  &point,
                                        const unsigned int component) const
  {
    AssertIndexRange(component, this->n_components);
    (void)component;

    return 1. - 2. / dim * (point.norm_square() - 1.);
  }



  template <int dim>
  double LaplaceSolver<dim>::compute_L2_error() const
  {
    std::cout << "Computing L2 error" << std::endl;

    const QGauss<dim> quadrature(fe_degree + 2);

    FEValues<dim> fe_values(fe_collection[0],
                            quadrature,
                            update_values | update_JxW_values |
                              update_quadrature_points);

    // AnalyticalSolution<dim> analytical_solution;
    double error_L2_squared = 0;

    for (const auto &cell :
         dof_handler.active_cell_iterators() |
           IteratorFilters::ActiveFEIndexEqualTo(ActiveFEIndex::lagrange))
      {
        fe_values.reinit(cell);



        std::vector<double> solution_values(fe_values.n_quadrature_points);
        fe_values.get_function_values(solution, solution_values);

        for (const unsigned int q : fe_values.quadrature_point_indices())
          {
            const Point<dim> &point = fe_values.quadrature_point(q);
            const double      error_at_point =
              solution_values.at(q) - analytical_solution.value(point);
            error_L2_squared +=
              Utilities::fixed_power<2>(error_at_point) * fe_values.JxW(q);
          }
      }

    return std::sqrt(error_L2_squared);
  }



  template <int dim>
  void LaplaceSolver<dim>::run()
  {
    ConvergenceTable   convergence_table;
    const unsigned int n_refinements = 5;

    make_grid();
    for (unsigned int cycle = 0; cycle <= n_refinements; cycle++)
      {
        std::cout << "Refinement cycle " << cycle << std::endl;
        triangulation.refine_global(1);
        setup_discrete_level_set();
        std::cout << "Classifying cells" << std::endl;
        mesh_classifier.reclassify();
        distribute_dofs();
        initialize_matrices();
        assemble_system();
        solve();
        // if (cycle == 1)
        output_results();
        const double error_L2 = compute_L2_error();
        const double cell_side_length =
          triangulation.begin_active()->minimum_vertex_distance();

        convergence_table.add_value("Cycle", cycle);
        convergence_table.add_value("Mesh size", cell_side_length);
        convergence_table.add_value("L2-Error", error_L2);

        convergence_table.evaluate_convergence_rates(
          "L2-Error", ConvergenceTable::reduction_rate_log2);
        convergence_table.set_scientific("L2-Error", true);

        std::cout << std::endl;
        convergence_table.write_text(std::cout);
        std::cout << std::endl;
      }
  }

} // namespace Step85



int main()
{
  const int                  dim = 2;
  Step85::LaplaceSolver<dim> laplace_solver;
  laplace_solver.run();
}
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

// @sect3{Include files}

// The first include files have all been treated in previous examples.
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

// These headers are needed for the shifted boundary method implementation
// There will be more once other PR get accepted.
#include <deal.II/non_matching/closest_surface_point.h>
#include <deal.II/non_matching/mesh_classifier.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <vector>


// @sect3{The LaplaceSolver class template}

namespace Step101
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

    // We use FunctionParser to define the right-hand side, boundary condition,
    // and analytical solution. If combined with  parameter parsing, this
    // allows for easy modification of the test problem without recompiling.
    FunctionParser<dim> rhs_function;
    FunctionParser<dim> boundary_condition;
    FunctionParser<dim> neumann_condition;
    FunctionParser<dim> analytical_solution;

    // The triangulation and level set function setup (same as in step-85)
    // followed by DoFHandler, finite element collection for the solution
    // and standard linear algebra objects
    Triangulation<dim> triangulation;

    const FE_Q<dim> fe_level_set;
    DoFHandler<dim> level_set_dof_handler;
    Vector<double>  level_set;

    hp::FECollection<dim> fe_collection;
    DoFHandler<dim>       dof_handler;
    Vector<double>        solution;


    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> stiffness_matrix;
    Vector<double>       rhs;

    // The mesh classifier helps us determine which cells are inside,
    // outside, or intersected by the level set
    NonMatching::MeshClassifier<dim> mesh_classifier;
  };



  // @sect4{Constructor}
  // In the constructor, we set up the finite element degree and define
  // the problem functions. Due to the naive form of the Neumann condition,
  // we loose one order of convergence, so we use a quadratic
  // finite element space to achieve a second-order convergence rate.
  template <int dim>
  LaplaceSolver<dim>::LaplaceSolver()
    : fe_degree(2)
    , rhs_function("2 * cos(x) * sin(y)")
    , boundary_condition("cos(x)*sin(y)")
    , neumann_condition("y * cos(x)* cos(y) - x * sin(x) * sin(y)")
    , analytical_solution("cos(x)*sin(y)")
    , fe_level_set(fe_degree)
    , level_set_dof_handler(triangulation)
    , dof_handler(triangulation)
    , mesh_classifier(level_set_dof_handler, level_set)
  {}



  // @sect4{Setting up the background mesh}
  // We create a Cartesian background mesh that covers the domain.
  // Since our domain is the unit disk, we need the mesh to extend
  // slightly beyond [-1,1]^dim.
  template <int dim>
  void LaplaceSolver<dim>::make_grid()
  {
    std::cout << "Creating background mesh" << std::endl;

    GridGenerator::hyper_cube(triangulation, -1.21, 1.21);
    triangulation.refine_global(2);
  }



  // @sect4{Setting up the discrete level set function}
  // The level set function @f$\phi(x)@f$ describes the geometry: @f$\phi < 0@f$
  // inside @f$\Omega@f$, @f$\phi = 0@f$ on the boundary @f$\Gamma@f$, and
  // @f$\phi > 0@f$ outside @f$\Omega@f$. We use a signed distance function to
  // the unit sphere for this purpose.
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


  // @sect4{Distributing degrees of freedom}
  // We use hp finite elements: FE_Q for cells inside the domain
  // and FE_Nothing for cells outside the domain.
  enum ActiveFEIndex
  {
    lagrange = 0,
    nothing  = 1
  };

  // The key difference from step-85 is that we only use finite elements
  // on cells that are completely inside the domain. Cells that are
  // intersected or outside get FE_Nothing elements.
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



  // @sect4{Initialize matrices and vectors}
  // Since we don't have flux terms between cells (unlike step-85 with its
  // ghost penalty), we can use the standard sparsity pattern.
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



  // @sect4{System assembly}
  // This is where the Shifted Boundary Method is implemented. The key idea
  // is to shift Nitsche penalty terms from the curved boundary to mesh faces.
  template <int dim>
  void LaplaceSolver<dim>::assemble_system()
  {
    std::cout << "Assembling" << std::endl;

    // Standard local assembly variables and  Nitsche penalty parameter (similar
    // to step-85)
    const unsigned int n_dofs_per_cell = fe_collection[0].dofs_per_cell;
    FullMatrix<double> local_stiffness(n_dofs_per_cell, n_dofs_per_cell);
    Vector<double>     local_rhs(n_dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(n_dofs_per_cell);
    const double nitsche_parameter = 5 * (fe_degree + 1) * fe_degree;


    // In the Shifted Boundary Method trial quantities are evaluated at
    // a shifted point on the true boundary while test functions remain
    // evaluated at the mesh quadrature point. As a result, when assembling
    // the Nitsche/penalty and flux terms we need shape-function values and
    // their gradients at the shifted location. The arrays below store these
    // precomputed values for all local dofs so they can be reused when
    // accumulating the face contributions (consistency, symmetry and penalty
    // terms). Gradients also need to be transformed appropriately for the
    // mapping (here handled approximately by scaling with a characteristic
    // cell length).
    std::vector<Point<dim>>     face_quad_points(n_dofs_per_cell);
    std::vector<Point<dim>>     shifted_points(n_dofs_per_cell);
    std::vector<double>         shifted_shape_value(n_dofs_per_cell);
    std::vector<Tensor<1, dim>> shifted_shape_grad(n_dofs_per_cell);

    // For debugging: we collect the shifts to output them later.
    // This allows us to verify that the shifting is working correctly.
    std::vector<std::pair<Point<dim>, Point<dim>>> shifts;
    std::vector<std::pair<Point<dim>, Point<dim>>> exact_shifts;


    // Quadrature rules and finite element objects
    const QGauss<dim - 1> face_quadrature(fe_degree + 1);
    const QGauss<dim>     cell_quadrature(fe_degree + 1);

    const FiniteElement<dim> &fe_lagrange =
      fe_collection[ActiveFEIndex::lagrange];

    // Contrary to step-85, we use regular FEValues for the volume integrals
    // since we only integrate over cells that are fully inside the domain.
    // No non-matching data is needed for this part.
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

    MappingCartesian<dim>                         cartesian_mapping;
    NonMatching::ClosestSurfacePoint<dim, double> closest_surface_point(
      level_set, level_set_dof_handler, cartesian_mapping);


    for (const auto &cell :
         dof_handler.active_cell_iterators() |
           IteratorFilters::ActiveFEIndexEqualTo(ActiveFEIndex::lagrange))
      {
        local_stiffness = 0;
        local_rhs       = 0;

        const double cell_side_length = cell->minimum_vertex_distance();
        fe_values.reinit(cell);

        // Standard volume integration for the Laplace operator
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

        // Now we loop over faces. If a face borders a cell with FE_Nothing
        // (i.e., a cell outside the domain), we apply shifted boundary
        // conditions.
        for (const auto &face : GeometryInfo<dim>::face_indices())
          {
            // Skip faces that border cells with active finite elements
            if (cell->neighbor(face)->active_fe_index() ==
                ActiveFEIndex::lagrange)
              continue;

            surface_fe_values.reinit(cell, face);

            // We need the neighboring cell as the "searching" cell for the
            // closest-surface-point routine. In the SBM we handle faces that
            // separate an active (inside-domain) cell from a FE_Nothing
            // neighbor (outside-domain). To find the shifted (closest) points
            // on the true boundary we call the closest_surface_point helper
            // with three pieces of information:
            //   - the cell to search from (here the neighbor, which typically
            //     lies outside the domain and contains the level-set zero set),
            //   - the reference cell that provides the face quadrature points
            //     (the active cell), and
            //   - the quadrature points on that face.
            //
            // The routine returns two arrays:
            //   - real_points: physical coordinates of the closest points on
            //   the
            //                  true boundary corresponding to each face
            //                  quadrature point, and
            //   - unit_points: the corresponding coordinates in the reference
            //                  cell of the searching cell (so we can evaluate
            //                  shape values/gradients there).
            //
            // We then use these returned points to evaluate trial functions at
            // the shifted locations and to transfer Dirichlet/Neumann data
            // from the true boundary back to the surrogate mesh face.
            auto neighbour = cell->neighbor(face);
            auto [real_points, unit_points] =
              closest_surface_point.compute_closest_surface_points(
                neighbour, // searching cell
                cell,      // reference cell for face quadrature points
                surface_fe_values.get_quadrature_points());

            // Loop over quadrature points on the face
            for (const unsigned int q :
                 surface_fe_values.quadrature_point_indices())
              {
                const Point<dim> &point = surface_fe_values.quadrature_point(q);

                // @sect5{Computing the shift}
                // For a point on the mesh boundary, we find the closest point
                // on the true boundary Γ. For a sphere, this is simply the
                // point projected onto the unit circle/sphere.
                const Point<dim> closest_boundary_point =
                  point * (1. / point.norm());


                // Transform the shifted point to the reference cell coordinates
                const Point<dim> &unit_shifted_point = unit_points[q];
                const Point<dim> &real_shifted_point = real_points[q];


                // Store shifts for debugging (will be exported using
                // functionality in PR #18741)
                shifts.push_back(std::make_pair(point, real_shifted_point));
                exact_shifts.push_back(
                  std::make_pair(point, point * (1. / point.norm())));

                // @sect5{Evaluate shifted shape functions}
                // The key ingredient of SBM: we evaluate shape functions at
                // the shifted point rather than the quadrature point
                // Since the mesh is cartesian, we have to scale the gradients
                // with the mesh size.
                // We will need both the mesh normal and  the normal to the true
                // boundary (for Neumann BC). The Dirichlet BC follow the
                // analysis in: https://doi.org/10.1016/j.cma.2022.114885,
                // while the Neumann BC are treated as in: TODO
                for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                  {
                    shifted_shape_value[j] =
                      fe_lagrange.shape_value(j, unit_shifted_point);
                    shifted_shape_grad[j] =
                      fe_lagrange.shape_grad(j, unit_shifted_point);

                    shifted_shape_grad[j] /= cell_side_length;
                  }

                const Tensor<1, dim> &mesh_normal =
                  surface_fe_values.normal_vector(q);

                // fixme: this only work for unit spehre.
                // We need functionality to compute normals from the level set.
                const Tensor<1, dim> &surface_normal = closest_boundary_point;

                // @sect5{Shifted Nitsche method assembly}
                // Assemble the shifted Nitsche penalty terms. The difference
                // from standard Nitsche is that test functions are evaluated
                // at quadrature points while trial functions are evaluated
                // at shifted points.
                if (cell->center()(0) < 0.0)
                  // Dirichlet condition on left half of boundary.
                  // We assemble the Nitsche-type integrals over the surrogate
                  // boundary (tilde-Gamma), i.e., the mesh faces that lie
                  // between an active cell and a FE_Nothing neighbor. All face
                  // integrals are evaluated with the chosen face quadrature and
                  // multiplied by the usual JxW. Note the key SBM feature: test
                  // functions and their gradients are evaluated at the
                  // quadrature point on the mesh face, whereas trial functions
                  // are evaluated at the shifted (closest-point) location on
                  // the true boundary. This is why we precompute
                  // shifted_shape_value and shifted_shape_grad above and scale
                  // the gradients with the characteristic cell length. The
                  // penalty term uses nitsche_parameter/cell_side_length and
                  // enforces the Dirichlet data transferred from the true
                  // boundary (see the introduction doc for the mathematical
                  // formulation).
                  for (const unsigned int i : surface_fe_values.dof_indices())
                    {
                      for (const unsigned int j :
                           surface_fe_values.dof_indices())
                        {
                          local_stiffness(i, j) +=
                            (-mesh_normal * surface_fe_values.shape_grad(i, q) *
                               shifted_shape_value[j] +
                             -mesh_normal * surface_fe_values.shape_grad(j, q) *
                               surface_fe_values.shape_value(i, q) +
                             nitsche_parameter / cell_side_length *
                               shifted_shape_value[i] *
                               shifted_shape_value[j]) *
                            surface_fe_values.JxW(q);
                        }
                      local_rhs(i) +=
                        boundary_condition.value(closest_boundary_point) *
                        (nitsche_parameter / cell_side_length *
                           shifted_shape_value[i] -
                         mesh_normal * surface_fe_values.shape_grad(i, q)) *
                        surface_fe_values.JxW(q);
                    }
                else
                  {
                    // Neumann condition on right half.
                    const double n_ntilde =
                      scalar_product(surface_normal, mesh_normal);

                    for (const unsigned int i : surface_fe_values.dof_indices())
                      {
                        for (const unsigned int j :
                             surface_fe_values.dof_indices())
                          {
                            local_stiffness(i, j) +=
                              (-mesh_normal *
                                 surface_fe_values.shape_grad(j, q) *
                                 surface_fe_values.shape_value(i, q)

                               + n_ntilde * surface_normal *
                                   shifted_shape_grad[j] *
                                   surface_fe_values.shape_value(i, q)

                                 ) *
                              surface_fe_values.JxW(q);
                          }
                        local_rhs(i) +=
                          n_ntilde *
                          neumann_condition.value(closest_boundary_point) *
                          surface_fe_values.shape_value(i, q) *
                          surface_fe_values.JxW(q);
                      }
                  }
              }
          }

        // Add local contributions to global system
        cell->get_dof_indices(local_dof_indices);
        stiffness_matrix.add(local_dof_indices, local_stiffness);
        rhs.add(local_dof_indices, local_rhs);
      }

    // Future work: Export the shifts for visualization and debugging
    // This functionality will be implemented in PR #18741
    // export_line_segments("shifts", shifts);
    // export_line_segments("shifts_exact", exact_shifts);
  }


  // @sect4{Solving the linear system}
  // We use a direct solver for simplicity. For larger problems,
  // iterative solvers would be more appropriate. The system matrix is
  // non-symmetric due to the Nitsche terms. For iterative solvers, GMRES is
  // always a good choice for non-symmetric systems, and in practice BiCGStab
  // also seems to work well for this problem.
  template <int dim>
  void LaplaceSolver<dim>::solve()
  {
    std::cout << "Solving system" << std::endl;

    SparseDirectUMFPACK solver_direct;
    solver_direct.initialize(stiffness_matrix);
    solver_direct.vmult(solution, rhs);
  }



  // @sect4{Output results}
  // We output both the solution and the level set function.
  // We exclude cells that are outside the domain from the output.
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


  // @sect4{Computing the L2 error}
  // We compute the L2 error only over cells that are inside the domain.
  // This gives us a measure of how well the SBM approximates the solution.
  // Similarly as in step-85 we loop only over cells with active finite
  // elements.
  template <int dim>
  double LaplaceSolver<dim>::compute_L2_error() const
  {
    std::cout << "Computing L2 error" << std::endl;

    const QGauss<dim> quadrature(fe_degree + 2);

    FEValues<dim> fe_values(fe_collection[0],
                            quadrature,
                            update_values | update_JxW_values |
                              update_quadrature_points);

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



  // @sect4{The main run function}
  // We perform a convergence study to verify that the SBM achieves
  // the expected convergence rates.
  template <int dim>
  void LaplaceSolver<dim>::run()
  {
    ConvergenceTable   convergence_table;
    const unsigned int n_refinements = 4;

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

} // namespace Step101



// @sect3{The main function}
int main()
{
  const int                   dim = 2;
  Step101::LaplaceSolver<dim> laplace_solver;
  laplace_solver.run();
}

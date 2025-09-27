// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// A combination of step-27 from 8.4 with corrected k-vectors, that is 2\pi*k
// instead of \pi*k and a new step-27 from 8.5 which use FESeries namespace. By
// default, the new version is used, but the blessed output file is obtained
// using the modified 8.4 version.


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_series.h>

#include <deal.II/grid/cell_data.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/hp/refinement.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/smoothness_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <complex>
#include <fstream>
#include <iostream>

#include "../tests.h"

namespace Step27
{
  using namespace dealii;



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
    assemble_system();
    void
    solve();
    void
    create_coarse_grid();
    void
    postprocess(const unsigned int cycle);

    Triangulation<dim> triangulation;

    DoFHandler<dim>          dof_handler;
    hp::FECollection<dim>    fe_collection;
    hp::QCollection<dim>     quadrature_collection;
    hp::QCollection<dim - 1> face_quadrature_collection;

    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;

    const unsigned int                      max_degree;
    hp::QCollection<dim>                    fourier_q_collection;
    std::unique_ptr<FESeries::Fourier<dim>> fourier;
  };



  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide()
      : Function<dim>()
    {}

    virtual double
    value(const Point<dim> &p, const unsigned int component) const override;
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
    : dof_handler(triangulation)
    , max_degree(dim <= 2 ? 7 : 5)
  {
    for (unsigned int degree = 2; degree <= max_degree; ++degree)
      {
        fe_collection.push_back(FE_Q<dim>(degree));
        quadrature_collection.push_back(QGauss<dim>(degree + 1));
        face_quadrature_collection.push_back(QGauss<dim - 1>(degree + 1));
      }

    QGauss<1>      base_quadrature(2);
    QIterated<dim> quadrature(base_quadrature, max_degree);
    for (unsigned int i = 0; i < fe_collection.size(); ++i)
      fourier_q_collection.push_back(quadrature);

    const std::vector<unsigned int> n_coefficients_per_direction(
      fe_collection.size(), max_degree);
    fourier =
      std::make_unique<FESeries::Fourier<dim>>(n_coefficients_per_direction,
                                               fe_collection,
                                               fourier_q_collection);
  }



  template <int dim>
  LaplaceProblem<dim>::~LaplaceProblem()
  {
    dof_handler.clear();
  }


  template <int dim>
  void
  LaplaceProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe_collection);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             constraints);
    constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
  }



  template <int dim>
  void
  LaplaceProblem<dim>::assemble_system()
  {
    hp::FEValues<dim> hp_fe_values(fe_collection,
                                   quadrature_collection,
                                   update_values | update_gradients |
                                     update_quadrature_points |
                                     update_JxW_values);

    const RightHandSide<dim> rhs_function;

    FullMatrix<double> cell_matrix;
    Vector<double>     cell_rhs;

    std::vector<types::global_dof_index> local_dof_indices;

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell != endc; ++cell)
      {
        const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

        cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
        cell_matrix = 0;

        cell_rhs.reinit(dofs_per_cell);
        cell_rhs = 0;

        hp_fe_values.reinit(cell);

        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

        std::vector<double> rhs_values(fe_values.n_quadrature_points);
        rhs_function.value_list(fe_values.get_quadrature_points(), rhs_values);

        for (unsigned int q_point = 0; q_point < fe_values.n_quadrature_points;
             ++q_point)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                cell_matrix(i, j) +=
                  (fe_values.shape_grad(i, q_point) *
                   fe_values.shape_grad(j, q_point) * fe_values.JxW(q_point));

              cell_rhs(i) += (fe_values.shape_value(i, q_point) *
                              rhs_values[q_point] * fe_values.JxW(q_point));
            }

        local_dof_indices.resize(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);

        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
  }



  template <int dim>
  void
  LaplaceProblem<dim>::solve()
  {
    SolverControl solver_control(system_rhs.size(),
                                 1e-8 * system_rhs.l2_norm());
    SolverCG<>    cg(solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    constraints.distribute(solution);
  }



  template <int dim>
  void
  LaplaceProblem<dim>::postprocess(const unsigned int cycle)
  {
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
    KellyErrorEstimator<dim>::estimate(
      dof_handler,
      face_quadrature_collection,
      std::map<types::boundary_id, const Function<dim> *>(),
      solution,
      estimated_error_per_cell);

    GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                    estimated_error_per_cell,
                                                    0.3,
                                                    0.03);

    Vector<float> smoothness_indicators;
    SmoothnessEstimator::Fourier::coefficient_decay(
      *fourier,
      dof_handler,
      solution,
      smoothness_indicators,
      /*regression_strategy=*/VectorTools::Linfty_norm,
      /*smallest_abs_coefficient=*/1e-10,
      /*only_flagged_cells=*/true);

    hp::Refinement::p_adaptivity_from_relative_threshold(dof_handler,
                                                         smoothness_indicators,
                                                         0.5,
                                                         0);
    hp::Refinement::choose_p_over_h(dof_handler);

    // Output to VTK
    if (false)
      {
        Vector<float> fe_degrees(triangulation.n_active_cells());
        {
          typename DoFHandler<dim>::active_cell_iterator
            cell = dof_handler.begin_active(),
            endc = dof_handler.end();
          for (; cell != endc; ++cell)
            fe_degrees(cell->active_cell_index()) =
              fe_collection[cell->active_fe_index()].degree;
        }

        DataOut<dim> data_out;

        data_out.attach_dof_handler(dof_handler);
        data_out.add_data_vector(solution, "solution");
        data_out.add_data_vector(estimated_error_per_cell, "error");
        data_out.add_data_vector(smoothness_indicators, "smoothness");
        data_out.add_data_vector(fe_degrees, "fe_degree");
        data_out.build_patches();

        const std::string filename =
          "solution-" + Utilities::int_to_string(cycle, 2) + ".vtk";
        std::ofstream output(filename);
        data_out.write_vtk(output);
      }

    triangulation.execute_coarsening_and_refinement();
  }



  template <>
  void
  LaplaceProblem<2>::create_coarse_grid()
  {
    const unsigned int dim = 2;

    static const Point<2> vertices_1[] = {
      Point<2>(-1., -1.),      Point<2>(-1. / 2, -1.),
      Point<2>(0., -1.),       Point<2>(+1. / 2, -1.),
      Point<2>(+1, -1.),

      Point<2>(-1., -1. / 2.), Point<2>(-1. / 2, -1. / 2.),
      Point<2>(0., -1. / 2.),  Point<2>(+1. / 2, -1. / 2.),
      Point<2>(+1, -1. / 2.),

      Point<2>(-1., 0.),       Point<2>(-1. / 2, 0.),
      Point<2>(+1. / 2, 0.),   Point<2>(+1, 0.),

      Point<2>(-1., 1. / 2.),  Point<2>(-1. / 2, 1. / 2.),
      Point<2>(0., 1. / 2.),   Point<2>(+1. / 2, 1. / 2.),
      Point<2>(+1, 1. / 2.),

      Point<2>(-1., 1.),       Point<2>(-1. / 2, 1.),
      Point<2>(0., 1.),        Point<2>(+1. / 2, 1.),
      Point<2>(+1, 1.)};
    const unsigned int n_vertices = sizeof(vertices_1) / sizeof(vertices_1[0]);
    const std::vector<Point<dim>> vertices(&vertices_1[0],
                                           &vertices_1[n_vertices]);
    static const int cell_vertices[][GeometryInfo<dim>::vertices_per_cell] = {
      {0, 1, 5, 6},
      {1, 2, 6, 7},
      {2, 3, 7, 8},
      {3, 4, 8, 9},
      {5, 6, 10, 11},
      {8, 9, 12, 13},
      {10, 11, 14, 15},
      {12, 13, 17, 18},
      {14, 15, 19, 20},
      {15, 16, 20, 21},
      {16, 17, 21, 22},
      {17, 18, 22, 23}};
    const unsigned int n_cells =
      sizeof(cell_vertices) / sizeof(cell_vertices[0]);

    std::vector<CellData<dim>> cells(n_cells, CellData<dim>());
    for (unsigned int i = 0; i < n_cells; ++i)
      {
        for (const unsigned int j : GeometryInfo<dim>::vertex_indices())
          cells[i].vertices[j] = cell_vertices[i][j];
        cells[i].material_id = 0;
      }

    triangulation.create_triangulation(vertices, cells, SubCellData());
    triangulation.refine_global(3);
  }



  template <int dim>
  void
  LaplaceProblem<dim>::run()
  {
    for (unsigned int cycle = 0; cycle < 6; ++cycle)
      {
        std::cout << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          create_coarse_grid();

        setup_system();

        std::cout << "   Number of active cells:       "
                  << triangulation.n_active_cells() << std::endl
                  << "   Number of degrees of freedom: " << dof_handler.n_dofs()
                  << std::endl
                  << "   Number of constraints       : "
                  << constraints.n_constraints() << std::endl;

        assemble_system();
        solve();
        postprocess(cycle);
      }
  }
} // namespace Step27



int
main()
{
  try
    {
      using namespace dealii;
      using namespace Step27;

      LaplaceProblem<2> laplace_problem;
      laplace_problem.run();
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
    }

  return 0;
}

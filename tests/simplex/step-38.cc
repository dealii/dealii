/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2021 - 2025 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

// Test step-38 on a simplex mesh.
// Note that currently the quadratic mapping for simplices is not implemented.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h> // new header
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

// simplex
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_wedge_p.h>

#include <deal.II/grid/grid_generator.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

#define USE_SIMPLEX

namespace Step38
{
  using namespace dealii;

  template <int spacedim>
  class LaplaceBeltramiProblem
  {
  public:
    LaplaceBeltramiProblem(const unsigned degree = 2);
    void
    run();

  private:
    static constexpr unsigned int dim = spacedim - 1;

    void
    make_grid_and_dofs();
    void
    assemble_system();
    void
    solve();
    void
    output_results() const;
    void
    compute_error() const;

    Triangulation<dim, spacedim> triangulation;
    DoFHandler<dim, spacedim>    dof_handler;

#ifdef USE_SIMPLEX
    FE_SimplexP<dim, spacedim> fe;
    MappingFE<dim, spacedim>   mapping;
#else
    FE_Q<dim, spacedim>     fe;
    MappingQ<dim, spacedim> mapping;
#endif

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;
  };

  template <int dim>
  class Solution : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    virtual Tensor<1, dim>
    gradient(const Point<dim>  &p,
             const unsigned int component = 0) const override;
  };

  template <>
  double
  Solution<2>::value(const Point<2> &p, const unsigned int) const
  {
    return (-2. * p[0] * p[1]);
  }

  template <>
  Tensor<1, 2>
  Solution<2>::gradient(const Point<2> &p, const unsigned int) const
  {
    Tensor<1, 2> return_value;
    return_value[0] = -2. * p[1] * (1 - 2. * p[0] * p[0]);
    return_value[1] = -2. * p[0] * (1 - 2. * p[1] * p[1]);

    return return_value;
  }

  template <>
  double
  Solution<3>::value(const Point<3> &p, const unsigned int) const
  {
    return (std::sin(numbers::PI * p[0]) * std::cos(numbers::PI * p[1]) *
            exp(p[2]));
  }

  template <>
  Tensor<1, 3>
  Solution<3>::gradient(const Point<3> &p, const unsigned int) const
  {
    using numbers::PI;

    Tensor<1, 3> return_value;

    return_value[0] = PI * cos(PI * p[0]) * cos(PI * p[1]) * exp(p[2]);
    return_value[1] = -PI * sin(PI * p[0]) * sin(PI * p[1]) * exp(p[2]);
    return_value[2] = sin(PI * p[0]) * cos(PI * p[1]) * exp(p[2]);

    return return_value;
  }

  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;
  };

  template <>
  double
  RightHandSide<2>::value(const Point<2> &p,
                          const unsigned int /*component*/) const
  {
    return (-8. * p[0] * p[1]);
  }

  template <>
  double
  RightHandSide<3>::value(const Point<3> &p,
                          const unsigned int /*component*/) const
  {
    using numbers::PI;

    Tensor<2, 3> hessian;

    hessian[0][0] = -PI * PI * sin(PI * p[0]) * cos(PI * p[1]) * exp(p[2]);
    hessian[1][1] = -PI * PI * sin(PI * p[0]) * cos(PI * p[1]) * exp(p[2]);
    hessian[2][2] = sin(PI * p[0]) * cos(PI * p[1]) * exp(p[2]);

    hessian[0][1] = -PI * PI * cos(PI * p[0]) * sin(PI * p[1]) * exp(p[2]);
    hessian[1][0] = -PI * PI * cos(PI * p[0]) * sin(PI * p[1]) * exp(p[2]);

    hessian[0][2] = PI * cos(PI * p[0]) * cos(PI * p[1]) * exp(p[2]);
    hessian[2][0] = PI * cos(PI * p[0]) * cos(PI * p[1]) * exp(p[2]);

    hessian[1][2] = -PI * sin(PI * p[0]) * sin(PI * p[1]) * exp(p[2]);
    hessian[2][1] = -PI * sin(PI * p[0]) * sin(PI * p[1]) * exp(p[2]);

    Tensor<1, 3> gradient;
    gradient[0] = PI * cos(PI * p[0]) * cos(PI * p[1]) * exp(p[2]);
    gradient[1] = -PI * sin(PI * p[0]) * sin(PI * p[1]) * exp(p[2]);
    gradient[2] = sin(PI * p[0]) * cos(PI * p[1]) * exp(p[2]);

    Point<3> normal = p;
    normal /= p.norm();

    return (-trace(hessian) + 2 * (gradient * normal) +
            (hessian * normal) * normal);
  }

  template <int spacedim>
  LaplaceBeltramiProblem<spacedim>::LaplaceBeltramiProblem(
    const unsigned degree)
    : dof_handler(triangulation)
    , fe(degree)
#ifdef USE_SIMPLEX
    // Quadratic mapping for simplex is still not implemented.
    , mapping(FE_SimplexP<dim, spacedim>(degree))
#else
    , mapping(degree)
#endif
  {}

  template <int spacedim>
  void
  LaplaceBeltramiProblem<spacedim>::make_grid_and_dofs()
  {
    {
      Triangulation<spacedim> volume_mesh;
      GridGenerator::half_hyper_ball(volume_mesh);

      const std::set<types::boundary_id> boundary_ids = {0};

#ifdef USE_SIMPLEX
      Triangulation<dim, spacedim> tria_temp;
      GridGenerator::extract_boundary_mesh(volume_mesh,
                                           tria_temp,
                                           boundary_ids);
      tria_temp.set_all_manifold_ids(0);
      tria_temp.set_manifold(0, SphericalManifold<dim, spacedim>());
      GridGenerator::convert_hypercube_to_simplex_mesh(tria_temp,
                                                       triangulation);
      for (const auto i : tria_temp.get_manifold_ids())
        if (i != numbers::flat_manifold_id)
          triangulation.set_manifold(i, tria_temp.get_manifold(i));

#else
      GridGenerator::extract_boundary_mesh(volume_mesh,
                                           triangulation,
                                           boundary_ids);
      triangulation.set_all_manifold_ids(0);
      triangulation.set_manifold(0, SphericalManifold<dim, spacedim>());
#endif
      triangulation.refine_global(4);
    }

    std::cout << "Surface mesh has " << triangulation.n_active_cells()
              << " cells." << std::endl;

    dof_handler.distribute_dofs(fe);

    std::cout << "Surface mesh has " << dof_handler.n_dofs()
              << " degrees of freedom." << std::endl;

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }

  template <int spacedim>
  void
  LaplaceBeltramiProblem<spacedim>::assemble_system()
  {
    system_matrix = 0;
    system_rhs    = 0;
#ifdef USE_SIMPLEX
    const QGaussSimplex<dim> quadrature_formula(fe.degree + 1);
#else
    const QGauss<dim> quadrature_formula(2 * fe.degree);
#endif
    FEValues<dim, spacedim> fe_values(mapping,
                                      fe,
                                      quadrature_formula,
                                      update_values | update_gradients |
                                        update_quadrature_points |
                                        update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<double>                  rhs_values(n_q_points);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    RightHandSide<spacedim> rhs;

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_matrix = 0;
        cell_rhs    = 0;

        fe_values.reinit(cell);

        rhs.value_list(fe_values.get_quadrature_points(), rhs_values);

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
              cell_matrix(i, j) += fe_values.shape_grad(i, q_point) *
                                   fe_values.shape_grad(j, q_point) *
                                   fe_values.JxW(q_point);

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            cell_rhs(i) += fe_values.shape_value(i, q_point) *
                           rhs_values[q_point] * fe_values.JxW(q_point);

        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              system_matrix.add(local_dof_indices[i],
                                local_dof_indices[j],
                                cell_matrix(i, j));

            system_rhs(local_dof_indices[i]) += cell_rhs(i);
          }
      }

    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(
      mapping, dof_handler, 0, Solution<spacedim>(), boundary_values);

    MatrixTools::apply_boundary_values(
      boundary_values, system_matrix, solution, system_rhs, false);
  }

  template <int spacedim>
  void
  LaplaceBeltramiProblem<spacedim>::solve()
  {
    SolverControl solver_control(solution.size(), 1e-7 * system_rhs.l2_norm());
    SolverCG<Vector<double>> cg(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve(system_matrix, solution, system_rhs, preconditioner);
  }

  template <int spacedim>
  void
  LaplaceBeltramiProblem<spacedim>::output_results() const
  {
    DataOut<dim, spacedim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution,
                             "solution",
                             DataOut<dim, spacedim>::type_dof_data);
    data_out.build_patches(mapping, mapping.get_degree());

    const std::string filename =
      "solution-" + std::to_string(spacedim) + "d.vtk";
    std::ofstream output(filename);
    data_out.write_vtk(output);
  }

  template <int spacedim>
  void
  LaplaceBeltramiProblem<spacedim>::compute_error() const
  {
    Vector<float> difference_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(mapping,
                                      dof_handler,
                                      solution,
                                      Solution<spacedim>(),
                                      difference_per_cell,
#ifdef USE_SIMPLEX
                                      QGaussSimplex<dim>(fe.degree + 1),
#else
                                      QGauss<dim>(2 * fe.degree +
                                                  1), // This also works on
                                                      // triangular meshes
#endif
                                      VectorTools::H1_norm);

    double h1_error = VectorTools::compute_global_error(triangulation,
                                                        difference_per_cell,
                                                        VectorTools::H1_norm);
    std::cout << "H1 error = " << h1_error << std::endl;
  }

  template <int spacedim>
  void
  LaplaceBeltramiProblem<spacedim>::run()
  {
    make_grid_and_dofs();
    assemble_system();
    solve();
    output_results();
    compute_error();
  }
} // namespace Step38

int
main()
{
  try
    {
      using namespace Step38;

      LaplaceBeltramiProblem<3> laplace_beltrami;
      laplace_beltrami.run();
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

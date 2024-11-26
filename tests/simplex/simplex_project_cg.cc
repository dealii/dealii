/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2023 - 2024 by the deal.II authors
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

// Verify convergence rates for continuous simplex elements.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_creator.h>
#include <deal.II/numerics/vector_tools_integrate_difference.h>
#include <deal.II/numerics/vector_tools_project.h>
#include <deal.II/numerics/vector_tools_rhs.h>

#include "../tests.h"

#define SIMPLEX

template <int dim>
class Solution : public Function<dim>
{
public:
  double
  value(const Point<dim> &p, const unsigned int = 0) const override
  {
    double u = 1.0;
    for (int d = 0; d < dim; ++d)
      u *= std::sin(numbers::PI * p[d]);

    return u;
  }
};

template <int dim>
class ForcingH1 : public Function<dim>
{
public:
  double
  value(const Point<dim> &p, const unsigned int = 0) const override
  {
    double u = 1.0;
    for (int d = 0; d < dim; ++d)
      u *= std::sin(numbers::PI * p[d]);

    return (1 + 0.0 * dim * numbers::PI * numbers::PI) * u;
  }
};

template <int dim>
void
test(const unsigned int degree, const bool use_amr = false)
{
#ifdef SIMPLEX
  FE_SimplexP<dim>   fe(degree);
  QGaussSimplex<dim> quadrature(degree + 1);
#else
  FE_Q<dim>   fe(degree);
  QGauss<dim> quadrature(degree + 1);
#endif
  deallog << "FE = " << fe.get_name() << std::endl;
  deallog << std::setprecision(4);
  if (use_amr)
    deallog << "using AMR" << std::endl;

  double previous_error = 1.0;

  for (unsigned int r = 0; r < (use_amr ? 3 : 4); ++r)
    {
      Triangulation<dim> tria_hex, tria;
      GridGenerator::hyper_cube(tria_hex);
      tria_hex.refine_global(r);
#ifdef SIMPLEX
      GridGenerator::convert_hypercube_to_simplex_mesh(tria_hex, tria);
#else
      tria.copy_triangulation(tria_hex);
#endif
      if (use_amr)
        {
          for (auto &cell : tria.active_cell_iterators())
            if (cell->index() % 3 == 0)
              cell->set_refine_flag();
          tria.execute_coarsening_and_refinement();
        }

      ReferenceCell   reference_cell = tria.begin_active()->reference_cell();
      DoFHandler<dim> dof_handler(tria);
      dof_handler.distribute_dofs(fe);

      Vector<double>            cell_errors(tria.n_active_cells());
      Vector<double>            solution(dof_handler.n_dofs());
      Solution<dim>             function;
      AffineConstraints<double> constraints;
      const auto               &mapping =
        reference_cell.template get_default_linear_mapping<dim>();
      DoFTools::make_hanging_node_constraints(dof_handler, constraints);
      constraints.close();

      SparsityPattern sparsity_pattern;
      {
        DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
        DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, true);
        sparsity_pattern.copy_from(dsp);
      }

      SparseMatrix<double> h1_matrix(sparsity_pattern);
      SparseMatrix<double> laplace_matrix(sparsity_pattern);

      MatrixCreator::create_mass_matrix(
        mapping, dof_handler, quadrature, h1_matrix, {}, constraints);
      MatrixCreator::create_laplace_matrix(
        mapping, dof_handler, quadrature, laplace_matrix, {}, constraints);

      h1_matrix.add(0.0, laplace_matrix);

      Vector<double> rhs(solution.size());
      VectorTools::create_right_hand_side(
        mapping, dof_handler, quadrature, ForcingH1<dim>(), rhs, constraints);

      SolverControl solver_control(1000, 1e-12 * rhs.l2_norm(), false, false);
      SolverCG<Vector<double>> cg(solver_control);

      cg.solve(h1_matrix, solution, rhs, PreconditionIdentity());
      constraints.distribute(solution);

      VectorTools::integrate_difference(mapping,
                                        dof_handler,
                                        solution,
                                        function,
                                        cell_errors,
                                        quadrature,
                                        VectorTools::L2_norm);
      const double L2_error =
        VectorTools::compute_global_error(tria,
                                          cell_errors,
                                          VectorTools::L2_norm);

      deallog << "L2 error = " << L2_error << std::endl;
      if (L2_error != 0.0)
        deallog << "ratio = " << previous_error / L2_error << std::endl;
      previous_error = L2_error;

#if 0
      if (dim == 2)
        {
          DataOut<dim> data_out;
          data_out.attach_dof_handler(dof_handler);
          data_out.add_data_vector(solution, "u");
          data_out.build_patches(2);

          std::ofstream output("out-" + std::to_string(degree) + "-" +
                               std::to_string(r) + ".vtu");
          data_out.write_vtu(output);
        }
#endif
    }
}

int
main()
{
  initlog();

  test<2>(1);
  test<2>(2);
  test<2>(3);

  test<3>(1);
  test<3>(2);
  test<3>(3);

  test<2>(1, true);
  test<2>(2, true);
  test<2>(3, true);
}

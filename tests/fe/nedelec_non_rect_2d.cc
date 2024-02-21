// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------
//
// By Ross Kynch
//
// Test to confirm that FE_Nedelec works on deformed elements in 2D when using
// Dirichlet boundary conditions (n x E = f, type). This is handled by the
// function project_boundary_values_curl_conforming_l2().
//
// This test solves the real valued curl-curl equation in 2D:
//
// curl(curl(E)) + E = Js
//
// where the solution is:
// E = (0, x^2)
//
// so, Js = (0,-2) + E.
//
// i.e the solution should be "exact" at p=2.
//
// The domain is a distorted quad after 1 refinement. The distortion is
// performed by GridTools::distort_random.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <sstream>

#include "../tests.h"



namespace polytest
{
  template <int dim>
  class SimplePolynomial : public Function<dim>
  {
  public:
    SimplePolynomial()
      : Function<dim>(2)
    {}
    void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>>   &values) const;

    void
    rhs_value_list(const std::vector<Point<dim>> &points,
                   std::vector<Vector<double>>   &values) const;
  };
  template <int dim>
  void
  SimplePolynomial<dim>::vector_value_list(
    const std::vector<Point<dim>> &points,
    std::vector<Vector<double>>   &values) const
  {
    Assert(dim == 2, ExcNotImplemented());
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));
    for (unsigned int i = 0; i < points.size(); ++i)
      {
        const Point<dim> &p = points[i];
        // non-zero curl-curl:
        values[i][0] = 0.0;
        values[i][1] = p[0] * p[0];
      }
  }

  template <int dim>
  void
  SimplePolynomial<dim>::rhs_value_list(
    const std::vector<Point<dim>> &points,
    std::vector<Vector<double>>   &values) const
  {
    Assert(dim == 2, ExcNotImplemented());
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));

    for (unsigned int i = 0; i < points.size(); ++i)
      {
        const Point<dim> &p = points[i];
        // non-zero curl-curl:
        values[i][0] = 0.0;
        values[i][1] = -2.0 + p[0] * p[0];
      }
  }

  template <int dim>
  class polytest
  {
  public:
    polytest(unsigned int degree);
    ~polytest();

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
    output_error();

    unsigned int p_order;
    unsigned int quad_order;

    Triangulation<dim>        tria;
    DoFHandler<dim>           dof_handler;
    FE_Nedelec<dim>           fe;
    AffineConstraints<double> constraints;
    SparsityPattern           sparsity_pattern;
    SparseMatrix<double>      system_matrix;
    Vector<double>            solution;
    Vector<double>            system_rhs;
  };

  template <int dim>
  polytest<dim>::polytest(unsigned int degree)
    : p_order(degree)
    , quad_order(2 * degree + 3)
    , dof_handler(tria)
    , fe(degree)
  {}
  template <int dim>
  polytest<dim>::~polytest()
  {
    dof_handler.clear();
  }

  template <int dim>
  void
  polytest<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());


    constraints.clear();
    SimplePolynomial<dim> boundary_function;

    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler,
      0,
      boundary_function,
      0,
      constraints,
      StaticMappingQ1<dim>::mapping);
    constraints.close();
    DynamicSparsityPattern c_sparsity(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    c_sparsity,
                                    constraints,
                                    false);

    sparsity_pattern.copy_from(c_sparsity);
    system_matrix.reinit(sparsity_pattern);
  }

  template <int dim>
  void
  polytest<dim>::assemble_system()
  {
    const QGauss<dim> test_quad(quad_order);
    FEValues<dim>     fe_values_test(fe,
                                 test_quad,
                                 update_values | update_gradients |
                                   update_quadrature_points |
                                   update_JxW_values);

    const QGauss<dim>  quadrature_formula(quad_order);
    const unsigned int n_q_points = quadrature_formula.size();

    const QGauss<dim - 1> face_quadrature_formula(quad_order);
    const unsigned int    n_face_q_points = face_quadrature_formula.size();

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    FEFaceValues<dim> fe_face_values(fe,
                                     face_quadrature_formula,
                                     update_values | update_quadrature_points |
                                       update_normal_vectors |
                                       update_JxW_values);

    const FEValuesExtractors::Vector vec(0);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;


    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    SimplePolynomial<dim>       right_hand_side;
    std::vector<Vector<double>> rhs_value_list(
      n_q_points, Vector<double>(fe.n_components()));

    typename DoFHandler<dim>::active_cell_iterator cell, endc;
    endc = dof_handler.end();
    cell = dof_handler.begin_active();
    for (; cell != endc; ++cell)
      {
        fe_values_test.reinit(cell);
        fe_values.reinit(cell);
        cell_matrix = 0;
        cell_rhs    = 0;

        right_hand_side.rhs_value_list(fe_values.get_quadrature_points(),
                                       rhs_value_list);
        for (const auto q : fe_values.quadrature_point_indices())
          {
            Tensor<1, dim> rhs_value;
            for (unsigned int d = 0; d < dim; ++d)
              {
                rhs_value[d] = rhs_value_list[q](d);
              }
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  {
                    cell_matrix(i, j) +=
                      (fe_values[vec].curl(i, q) * fe_values[vec].curl(j, q) +
                       fe_values[vec].value(i, q) *
                         fe_values[vec].value(j, q)) *
                      fe_values.JxW(q);
                  }
                cell_rhs(j) +=
                  rhs_value * fe_values[vec].value(j, q) * fe_values.JxW(q);
              }
          }
        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
  }

  template <int dim>
  void
  polytest<dim>::solve()
  {
    SparseDirectUMFPACK direct;
    direct.initialize(system_matrix);

    direct.vmult(solution, system_rhs);
    constraints.distribute(solution);
  }

  template <int dim>
  void
  polytest<dim>::output_error()
  {
    SimplePolynomial<dim> exact_solution;
    Vector<double>        diff_per_cell(tria.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      exact_solution,
                                      diff_per_cell,
                                      QGauss<dim>(quad_order),
                                      VectorTools::L2_norm);
    const double L2_error = diff_per_cell.l2_norm();

    deallog << "p=" << p_order << " L2_error: " << L2_error << std::endl;
  }

  template <int dim>
  void
  polytest<dim>::run()
  {
    GridGenerator::hyper_cube(tria, -1.2, 1.3);
    // REFINE
    tria.refine_global(1);

    // DISTORT ALL:
    GridTools::distort_random(0.2, tria, false);

    setup_system();
    assemble_system();
    solve();
    output_error();
  }
} // namespace polytest
int
main()
{
  const unsigned int dim(2);

  initlog();

  for (unsigned int p = 0; p < 3; ++p)
    {
      polytest::polytest<dim> poly(p);
      poly.run();
    }
}

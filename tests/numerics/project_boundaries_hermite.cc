// ---------------------------------------------------------------------
//
// Copyright (C) 2023 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/base/config.h>

#define PRECISION 8


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_hermite.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_cartesian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <sstream>
#include <vector>



/**
 * Test of Hermite finite elements with the Laplace equation on a regular grid
 * in 1,2,3D. <code>FE_Hermite<dim>(poly_degree)<\code> should be able to
 * perfectly
 * represent any polynomial function up to degree @p poly_degree, including
 * on the boundaries. If all basis functions are correctly scaled according to
 * element size, then solving the Laplace equation with a polynomial solution
 * in the Hermite FE space will produce negligible pointwise errors for
 * non-homogeneous Dirichlet boundary conditions.
 */
using namespace dealii;

// Define a function to project onto the domain [-1,1]^dim
template <int dim>
class Solution : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p, unsigned int c = 0) const override
  {
    (void)c;
    if (dim == 0)
      Assert(false, ExcNotImplemented());
    double temp = p(0) * (1.0 + p(0) * (0.5 - p(0)));
    if (dim > 1)
      temp *= 1.0 - p(1) * p(1);
    if (dim == 3)
      temp *= p(2);
    return temp;
  }



  std::string
  get_function_string()
  {
    switch (dim)
      {
        case 1:
          return "X + 0.5 X^2 - X^3";
        case 2:
          return "(X + 0.5 X^2 - X^3)(1 - Y^2)";
        case 3:
          return "(X + 0.5 X^2 - X^3)(1 - Y^2)Z";
        default:
          Assert(false, ExcNotImplemented());
          return "";
      }
  }
};



// This function should be the analytical Laplacian of the solution
// function above.
template <int dim>
class RHSFunction : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p, unsigned int c = 0) const override
  {
    (void)c;
    if (dim == 0)
      Assert(false, ExcNotImplemented());
    double temp = -1.0 + 6 * p(0);
    if (dim > 1)
      {
        temp *= -p(1) * p(1);
        temp -= 1.0 - 8.0 * p(0) - p(0) * p(0) + 2.0 * p(0) * p(0) * p(0);
      }
    if (dim == 3)
      temp *= p(2);
    return temp;
  }
};



template <int dim>
void
test_fe_on_domain(const unsigned int regularity)
{
  deallog << std::endl;
  char fname[50];
  sprintf(fname, "Cell-%dd-Hermite-%d", dim, regularity);
  deallog.push(fname);

  deallog.depth_file(2);

  deallog << "Test polynomial:" << std::endl;

  Triangulation<dim> tr;
  DoFHandler<dim>    dof(tr);

  double left = -1.0, right = 1.0;
  GridGenerator::subdivided_hyper_cube(tr, 4, left, right);

  FE_Hermite<dim> herm(2 * regularity + 1);
  dof.distribute_dofs(herm);

  MappingCartesian<dim> mapping;

  QGauss<dim> quadr(2 * regularity + 2);

  Vector<double> sol(dof.n_dofs());
  Vector<double> rhs(dof.n_dofs());

  Solution<dim>    sol_object;
  RHSFunction<dim> rhs_object;

  deallog << sol_object.get_function_string() << std::endl;
  deallog << std::endl;

  AffineConstraints<double> constraints;
  constraints.close();

  DynamicSparsityPattern dsp(dof.n_dofs());
  DoFTools::make_sparsity_pattern(dof, dsp);
  SparsityPattern sp;
  sp.copy_from(dsp);

  SparseMatrix<double> stiffness_matrix;
  stiffness_matrix.reinit(sp);
  MatrixCreator::create_laplace_matrix(mapping, dof, quadr, stiffness_matrix);

  FEValues<dim> fe_herm(mapping,
                        herm,
                        quadr,
                        update_values | update_gradients |
                          update_quadrature_points | update_JxW_values);

  std::vector<types::global_dof_index> local_to_global(herm.n_dofs_per_cell());

  for (const auto &cell : dof.active_cell_iterators())
    {
      fe_herm.reinit(cell);
      cell->get_dof_indices(local_to_global);
      for (const unsigned int i : fe_herm.dof_indices())
        {
          double rhs_temp = 0;
          for (const unsigned int q : fe_herm.quadrature_point_indices())
            rhs_temp += fe_herm.shape_value(i, q) *
                        rhs_object.value(fe_herm.quadrature_point(q)) *
                        fe_herm.JxW(q);

          rhs(local_to_global[i]) += rhs_temp;
        }
    }

  std::map<types::global_dof_index, double> bound_vals;

  std::map<types::boundary_id, const Function<dim, double> *> bound_map;
  bound_map.emplace(std::make_pair(0U, &sol_object));
  if (dim == 1)
    bound_map.emplace(std::make_pair(1U, &sol_object));

  // The following is the main function being tested here
  VectorTools::project_hermite_boundary_values(mapping,
                                               dof,
                                               bound_map,
                                               QGauss<dim - 1>(2 * regularity +
                                                               2),
                                               0,
                                               bound_vals);

  MatrixTools::apply_boundary_values(bound_vals, stiffness_matrix, sol, rhs);

  IterationNumberControl solver_control_values(8000, 1e-11);
  SolverCG<>             solver(solver_control_values);

  solver.solve(stiffness_matrix, sol, rhs, PreconditionIdentity());

  double err_sq = 0;

  for (auto &cell : dof.active_cell_iterators())
    {
      fe_herm.reinit(cell);
      cell->get_dof_indices(local_to_global);
      for (const unsigned int q : fe_herm.quadrature_point_indices())
        {
          double sol_at_point = 0;
          for (const unsigned int i : fe_herm.dof_indices())
            sol_at_point += fe_herm.shape_value(i, q) * sol(local_to_global[i]);

          sol_at_point -= sol_object.value(fe_herm.quadrature_point(q));
          err_sq += sol_at_point * sol_at_point * fe_herm.JxW(q);
        }
    }

  err_sq = std::sqrt(err_sq);

  deallog.depth_file(2);

  deallog << "Interpolation error:" << std::endl;
  deallog << err_sq << "\n\n" << std::endl;
  deallog.pop();
}



int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);

  test_fe_on_domain<1>(0);
  test_fe_on_domain<1>(1);
  test_fe_on_domain<1>(2);
  test_fe_on_domain<1>(3);

  test_fe_on_domain<2>(0);
  test_fe_on_domain<2>(1);
  test_fe_on_domain<2>(2);
  test_fe_on_domain<2>(3);

  test_fe_on_domain<3>(0);
  test_fe_on_domain<3>(1);

  return 0;
}

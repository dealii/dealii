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
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_cartesian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <cmath>
#include <fstream>
#include <iterator>
#include <sstream>
#include <vector>



/*
 * Test case for Hermite with a time-stepping method for $u_{tt}=u_{xx}$
 * to check that the stability conditions have transferred correctly. Uses a
 * manufactured solution of a plane wave entering the domain at $x=0$
 * and leaving at $x=3$. The solution uses $y = x-t-0.5$ as a local
 * coordinate and has the following Gaussian bell shape to avoid
 * discontinuities in derivatives:
 *
 * @f{align*}{
 *  u(x,t) = exp(-a y^2),
 * @f}
 *
 * where $a$ is some predefined constant. The right-most tail of the wave
 * is already in the domain at $t=0$. The solution is chosen so that initial
 * conditions are not uniformly zero and the wave is currently entering
 * the domain through the boundary. This allows both the use of
 * boundary projection and initial projection for both value and
 * derivatives to be tested.
 */
using namespace dealii;

template <int dim>
class Solution : public Function<dim>
{
public:
  Solution(const Solution &sol) = default;
  Solution(double initial_time = 0, double alpha = 40, bool boundary = false)
    : alpha(alpha)
    , current_time(initial_time)
    , zero_boundary(boundary){};



  double inline update_time(double time_step)
  {
    return this->current_time += time_step;
  }



  double inline get_time() const
  {
    return this->current_time;
  }



  bool inline is_boundary_zero() const
  {
    return this->zero_boundary;
  }



  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    (void)component;
    if (zero_boundary)
      {
        double return_val = 1.;
        for (unsigned int i = 0; i < dim; ++i)
          return_val *= std::sin(numbers::PI * p(i));
        return_val *= std::cos(numbers::PI * std::sqrt(dim) * current_time);
        return return_val;
      }
    else
      {
        double y = p(0) + 0.5 - current_time;
        return std::exp(-alpha * y * y);
      }
  }



  double alpha;

private:
  double current_time;
  bool   zero_boundary;
};



template <int dim>
class SolutionDerivative : public Function<dim>
{
public:
  SolutionDerivative() = delete;
  SolutionDerivative(const Solution<dim> &parent)
    : parent_solution(parent){};



  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    (void)component;
    if (parent_solution.is_boundary_zero())
      {
        double return_val = -numbers::PI * std::sqrt(dim);
        for (unsigned int i = 0; i < dim; ++i)
          return_val *= std::sin(numbers::PI * p(i));
        return_val *= std::sin(numbers::PI * std::sqrt(dim) *
                               this->parent_solution.get_time());
        return return_val;
      }
    else
      {
        double y = p(0) + 0.5 - this->parent_solution.get_time();
        return 2 * (parent_solution.alpha) * y *
               std::exp(-(parent_solution.alpha) * y * y);
      }
  }



private:
  Solution<dim> parent_solution;
};



template <int dim>
double
get_cfl_number(const unsigned int regularity)
{
  double c_eff = 1.0;
  return c_eff / (std::sqrt(12.0 * dim) * (regularity + 1));
}



template <int dim>
double
l2_error(const Mapping<dim>    &mapping,
         const DoFHandler<dim> &dof,
         const Quadrature<dim> &quadrature,
         const Function<dim>   &solution,
         const Vector<double>  &fe_solution)
{
  FEValues<dim>                        herm_vals(mapping,
                          dof.get_fe(),
                          quadrature,
                          update_values | update_quadrature_points |
                            update_JxW_values);
  std::vector<types::global_dof_index> local_to_global(
    dof.get_fe().n_dofs_per_cell());
  double err_sq = 0;

  for (const auto &cell : dof.active_cell_iterators())
    {
      herm_vals.reinit(cell);
      cell->get_dof_indices(local_to_global);

      for (const unsigned int q : herm_vals.quadrature_point_indices())
        {
          double diff = solution.value(herm_vals.quadrature_point(q));

          for (const unsigned int i : herm_vals.dof_indices())
            diff -=
              fe_solution(local_to_global[i]) * herm_vals.shape_value(i, q);

          err_sq += herm_vals.JxW(q) * diff * diff;
        }
    }

  return std::sqrt(err_sq);
}



template <int dim>
void
test_wave_solver(const double initial_time, const unsigned int regularity)
{
  Triangulation<dim> tr;
  DoFHandler<dim>    dof(tr);

  double x_left = 0.0, x_right = 3.0;
  int    divisions  = 8;
  double dx         = (x_right - x_left) / divisions;
  double dt         = dx * get_cfl_number<dim>(regularity);
  double final_time = 4.0;

  GridGenerator::subdivided_hyper_cube(tr, divisions, x_left, x_right);

  MappingCartesian<dim> mapping_h;
  FE_Hermite<dim>       fe_h(2 * regularity + 1);
  dof.distribute_dofs(fe_h);

  AffineConstraints<double> constraints;
  constraints.close();

  QGauss<dim>     quadrature(2 * regularity + 2);
  QGauss<dim - 1> face_quadrature(2 * regularity + 2);

  Solution<dim>           wave(initial_time);
  SolutionDerivative<dim> wave_tdev(wave);

  Vector<double> initial_conds_value(dof.n_dofs());
  Vector<double> iniital_conds_tdev(dof.n_dofs());

  Vector<double> sol_prev(dof.n_dofs());
  Vector<double> sol_curr(dof.n_dofs());
  Vector<double> sol_next(dof.n_dofs());

  Vector<double> temp(dof.n_dofs());

  SparsityPattern sp;
  {
    DynamicSparsityPattern dsp(dof.n_dofs());
    DoFTools::make_sparsity_pattern(dof, dsp);
    sp.copy_from(dsp);
  }

  SparseMatrix<double> mass, mass_solve, tstep_noninvert;
  mass.reinit(sp), mass_solve.reinit(sp), tstep_noninvert.reinit(sp);

  MatrixCreator::create_mass_matrix(mapping_h, dof, quadrature, mass);
  MatrixCreator::create_laplace_matrix(mapping_h,
                                       dof,
                                       quadrature,
                                       tstep_noninvert);
  tstep_noninvert *= -0.5 * dt * dt;
  tstep_noninvert.add(1.0, mass);

  VectorTools::project_hermite(mapping_h,
                               dof,
                               constraints,
                               quadrature,
                               wave,
                               sol_curr,
                               false,
                               face_quadrature,
                               true);
  VectorTools::project_hermite(mapping_h,
                               dof,
                               constraints,
                               quadrature,
                               wave_tdev,
                               sol_prev,
                               false,
                               face_quadrature,
                               true);

  std::map<types::global_dof_index, double>           boundary_values;
  std::map<types::boundary_id, const Function<dim> *> boundary_functions;

  boundary_functions.emplace(std::make_pair(0U, &wave));
  if (dim == 1)
    boundary_functions.emplace(std::make_pair(1U, &wave));

  std::vector<double> l2_errors;
  double              max_error;
  l2_errors.emplace_back(
    l2_error<dim>(mapping_h, dof, quadrature, wave, sol_curr));
  max_error = l2_errors[0];

  // Make an initial time-step using
  // u(t+dt) = u(t) + dt u_t(t) +
  //              0.5 dt^2 M^(-1)A u(t) + O(dt^3, dt^2 dx^(p+1))
  tstep_noninvert.vmult(temp, sol_curr);
  mass.vmult(sol_next, sol_prev);
  sol_next *= dt;
  temp += sol_next;
  sol_next = 0;

  wave.update_time(dt);
  VectorTools::project_hermite_boundary_values(
    mapping_h, dof, boundary_functions, face_quadrature, 0, boundary_values);

  mass_solve.copy_from(mass);
  MatrixTools::apply_boundary_values(
    boundary_values, mass_solve, sol_next, temp, true);

  SparseDirectUMFPACK mass_inv;
  mass_inv.initialize(mass_solve);
  mass_inv.vmult(sol_next, temp);

  sol_prev = sol_curr;
  sol_curr = sol_next;
  sol_next = 0;

  boundary_values.clear();

  l2_errors.emplace_back(
    l2_error<dim>(mapping_h, dof, quadrature, wave, sol_curr));
  max_error = std::max(max_error, *l2_errors.cend());

  while (wave.get_time() < final_time)
    {
      mass.vmult(sol_next, sol_prev);
      tstep_noninvert.vmult(temp, sol_curr);
      temp *= 2.0;
      temp -= sol_next;
      sol_next = 0;

      wave.update_time(dt);
      VectorTools::project_hermite_boundary_values(mapping_h,
                                                   dof,
                                                   boundary_functions,
                                                   face_quadrature,
                                                   0,
                                                   boundary_values);

      mass_solve.copy_from(mass);
      MatrixTools::apply_boundary_values(
        boundary_values, mass_solve, sol_next, temp, true);
      mass_inv.vmult(sol_next, temp);

      sol_prev = sol_curr;
      sol_curr = sol_next;
      sol_next = 0;
      boundary_values.clear();

      l2_errors.emplace_back(
        l2_error<dim>(mapping_h, dof, quadrature, wave, sol_curr));
      max_error = std::max(max_error, *(--l2_errors.cend()));
    }

  deallog << std::endl;
  char fname[50];
  sprintf(fname, "Cell-%dd-Hermite-%d", dim, regularity);
  deallog.push(fname);

  deallog << "initial_time: " << initial_time << std::endl;
  deallog << "Final time: " << final_time << std::endl;
  deallog << std::endl;

  deallog << "Error at final time: " << l2_errors.back() << std::endl;
  deallog << "Largest error: " << max_error << std::endl;
  deallog << std::endl;

  double t = initial_time - dt;
  deallog << "Error time series:";
  for (const auto &it : l2_errors)
    deallog << std::endl << (t += dt) << "\t" << it;
  deallog << "\n\n" << std::endl;

  deallog.pop();
}



int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);

  test_wave_solver<1>(-0.3, 1);

  test_wave_solver<1>(0.0, 1);
  test_wave_solver<1>(0.0, 2);
  test_wave_solver<1>(0.0, 3);

  test_wave_solver<2>(0.0, 1);
  test_wave_solver<2>(0.0, 2);
  test_wave_solver<2>(0.0, 3);

  return 0;
}

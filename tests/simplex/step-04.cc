/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2021 - 2024 by the deal.II authors
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
 * Author: Wolfgang Bangerth, University of Heidelberg, 1999
 */


// Step-04 on a simplex mesh. Following incompatible modifications had to be
// made:
//  - Change the FE_Q to FE_SimplexP
//  - Put the MappingFE as a class member and use as an argument instead of
//  default mapping
//  - Change QGauss to QGaussSimplex
//  - Convert triangulation to a triangulation based on simplices


#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "../tests.h"



template <int dim>
class Step4
{
public:
  Step4();
  void
  run();

private:
  void
  make_grid();
  void
  setup_system();
  void
  assemble_system();
  void
  solve();
  void
  output_results() const;

  Triangulation<dim> triangulation;
  FE_SimplexP<dim>   fe;
  DoFHandler<dim>    dof_handler;

  MappingFE<dim> mapping;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> solution;
  Vector<double> system_rhs;
};


template <int dim>
class RightHandSide : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override;
};


template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override;
};


template <int dim>
double
RightHandSide<dim>::value(const Point<dim> &p,
                          const unsigned int /*component*/) const
{
  double return_value = 0.0;
  for (unsigned int i = 0; i < dim; ++i)
    return_value += 4.0 * std::pow(p[i], 4.0);

  return return_value;
}


template <int dim>
double
BoundaryValues<dim>::value(const Point<dim> &p,
                           const unsigned int /*component*/) const
{
  return p.square();
}



template <int dim>
Step4<dim>::Step4()
  : fe(1)
  , dof_handler(triangulation)
  , mapping(fe)
{}


template <int dim>
void
Step4<dim>::make_grid()
{
  Triangulation<dim, dim> temp_tria;
  GridGenerator::subdivided_hyper_cube(temp_tria, 8, -1., 1., false);
  GridGenerator::convert_hypercube_to_simplex_mesh<dim, dim>(temp_tria,
                                                             triangulation);

  deallog << "   Number of active cells: " << triangulation.n_active_cells()
          << std::endl
          << "   Total number of cells: " << triangulation.n_cells()
          << std::endl;
}


template <int dim>
void
Step4<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  deallog << "   Number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl;

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}


template <int dim>
void
Step4<dim>::assemble_system()
{
  QGaussSimplex<dim> quadrature_formula(fe.degree + 1);

  RightHandSide<dim> right_hand_side;

  const UpdateFlags flag = update_JxW_values | update_values |
                           update_gradients | update_quadrature_points;

  FEValues<dim> fe_values(mapping, fe, quadrature_formula, flag);


  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell_matrix = 0;
      cell_rhs    = 0;

      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        for (const unsigned int i : fe_values.dof_indices())
          {
            for (const unsigned int j : fe_values.dof_indices())
              cell_matrix(i, j) +=
                (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                 fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                 fe_values.JxW(q_index));           // dx

            const auto x_q = fe_values.quadrature_point(q_index);
            cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                            right_hand_side.value(x_q) *        // f(x_q)
                            fe_values.JxW(q_index));            // dx
          }

      cell->get_dof_indices(local_dof_indices);
      for (const unsigned int i : fe_values.dof_indices())
        {
          for (const unsigned int j : fe_values.dof_indices())
            system_matrix.add(local_dof_indices[i],
                              local_dof_indices[j],
                              cell_matrix(i, j));

          system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    }

  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(
    mapping, dof_handler, 0, BoundaryValues<dim>(), boundary_values);
  MatrixTools::apply_boundary_values(boundary_values,
                                     system_matrix,
                                     solution,
                                     system_rhs);
}



template <int dim>
void
Step4<dim>::solve()
{
  SolverControl            solver_control(1000, 1e-12);
  SolverCG<Vector<double>> solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());

  deallog << "   " << solver_control.last_step()
          << " CG iterations needed to obtain convergence." << std::endl;
}



template <int dim>
void
Step4<dim>::output_results() const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");

  data_out.build_patches(mapping, 0);

  std::ofstream output(dim == 2 ? "solution-2d.vtk" : "solution-3d.vtk");
  // data_out.write_vtk(output);
}


template <int dim>
void
Step4<dim>::run()
{
  deallog << "Solving problem in " << dim << " space dimensions." << std::endl;

  make_grid();
  setup_system();
  assemble_system();
  solve();
  output_results();
}



int
main()
{
  initlog();
  deallog.depth_console(2);
  {
    Step4<2> laplace_problem_2d;
    laplace_problem_2d.run();
  }

  {
    Step4<3> laplace_problem_3d;
    laplace_problem_3d.run();
  }

  return 0;
}

/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2021 by the deal.II authors
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
 * Authors: Wolfgang Bangerth, 1999,
 *          Guido Kanschat, 2011,
 *          Peter Munch, 2021
 */


// @sect3{Include files}

// Include files, as used in step-3:
#include <deal.II/base/function.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <fstream>
#include <iostream>

// Include files, as added in step-3simplex:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/grid/grid_in.h>

// Include files that we need in this tutorial to be able to deal with
// collections of finite elements, quadrature rules, mapping objects, and
// FEValues.
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

using namespace dealii;

// @sect3{The <code>Step3</code> class}
//
// This is the main class of the tutorial. Since it is very similar to the
// version from step-3 and step-3simplex, we will only point out and explain
// the relevant differences that allow to perform simulations on mixed meshes.

class Step3
{
public:
  Step3();

  void run();

private:
  void make_grid();
  void setup_system();
  void assemble_system();
  void solve();
  void output_results() const;

  Triangulation<2> triangulation;

  // As already explained, we are not working with mapping objects, finite
  // elements, and quadrature rules directly but with collections of them.
  const hp::MappingCollection<2> mapping;
  const hp::FECollection<2>      fe;
  const hp::QCollection<2>       quadrature_formula;

  DoFHandler<2> dof_handler;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> solution;
  Vector<double> system_rhs;
};


// @sect4{Step3::Step3}
//
// In the constructor of the Step3 class, we fill the collections. Here, we
// position the objects related to triangles in the first place (index 0) and
// the ones related to quadrilaterals in the second place (index 1).
Step3::Step3()
  : mapping(MappingFE<2>(FE_SimplexP<2>(1)), MappingFE<2>(FE_Q<2>(1)))
  , fe(FE_SimplexP<2>(2), FE_Q<2>(2))
  , quadrature_formula(QGaussSimplex<2>(3), QGauss<2>(3))
  , dof_handler(triangulation)
{}


// @sect4{Step3::make_grid}
//
// Read the external mesh file "box_2D_mixed.msh" as in step-3simplex.
void Step3::make_grid()
{
  GridIn<2>(triangulation).read("box_2D_mixed.msh");

  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;
}


// @sect4{Step3::setup_system}
//
// In contrast to step-3 and step-3simplex, we need here a preprocessing step
// that assigns an active_fe_index to each cell consistently according to the
// indices in the collections and the cell type.
void Step3::setup_system()
{
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->reference_cell() == ReferenceCells::Triangle)
        cell->set_active_fe_index(0);
      else if (cell->reference_cell() == ReferenceCells::Quadrilateral)
        cell->set_active_fe_index(1);
      else
        Assert(false, ExcNotImplemented());
    }

  dof_handler.distribute_dofs(fe);
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}


// @sect4{Step3::assemble_system}
//
// The following function looks similar to the version in step-3 and
// step-3simplex with the following two differences:
//  - We do not work with FEValues directly but with the collection class
//    hp::FEValues. It gives us - after it has been initialized with the current
//    cell - a reference to the right FEValues object (constructed
//    with the correct mapping object, finite element, and quadrature rule),
//    which can be used as usual to compute the cell integrals.
//  - The cell-local stiffness matrix and the right-hand-side vector have
//    different sizes depending on the cell type (6 DoFs vs. 9 DoFs) so that
//    they might need to be resized for each cell.
//
// Apart from these two changes, the code has not changed. Not even, the
// cell integrals have been changed depending on whether one operates on
// hypercube, simplex, or mixed meshes.
void Step3::assemble_system()
{
  hp::FEValues<2> hp_fe_values(mapping,
                               fe,
                               quadrature_formula,
                               update_values | update_gradients |
                                 update_JxW_values);

  FullMatrix<double>                   cell_matrix;
  Vector<double>                       cell_rhs;
  std::vector<types::global_dof_index> local_dof_indices;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      hp_fe_values.reinit(cell);

      const auto &fe_values = hp_fe_values.get_present_fe_values();

      const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();
      cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
      cell_rhs.reinit(dofs_per_cell);
      local_dof_indices.resize(dofs_per_cell);

      cell_matrix = 0;
      cell_rhs    = 0;

      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
          for (const unsigned int i : fe_values.dof_indices())
            for (const unsigned int j : fe_values.dof_indices())
              cell_matrix(i, j) +=
                (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                 fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                 fe_values.JxW(q_index));           // dx

          for (const unsigned int i : fe_values.dof_indices())
            cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                            1. *                                // f(x_q)
                            fe_values.JxW(q_index));            // dx
        }
      cell->get_dof_indices(local_dof_indices);

      for (const unsigned int i : fe_values.dof_indices())
        for (const unsigned int j : fe_values.dof_indices())
          system_matrix.add(local_dof_indices[i],
                            local_dof_indices[j],
                            cell_matrix(i, j));

      for (const unsigned int i : fe_values.dof_indices())
        system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }


  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(
    mapping, dof_handler, 0, Functions::ZeroFunction<2>(), boundary_values);
  MatrixTools::apply_boundary_values(boundary_values,
                                     system_matrix,
                                     solution,
                                     system_rhs);
}


// @sect4{Step3::solve}
//
// Nothing has changed here.
void Step3::solve()
{
  SolverControl            solver_control(1000, 1e-12);
  SolverCG<Vector<double>> solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
}


// @sect4{Step3::output_results}
//
// Nothing has changed here.
void Step3::output_results() const
{
  DataOut<2> data_out;

  DataOutBase::VtkFlags flags;
  flags.write_higher_order_cells = true;
  data_out.set_flags(flags);

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches(mapping, 2);
  std::ofstream output("solution.vtk");
  data_out.write_vtk(output);
}


// @sect4{Step3::run}
//
// Nothing has changed here.
void Step3::run()
{
  make_grid();
  setup_system();
  assemble_system();
  solve();
  output_results();
}


// @sect3{The <code>main</code> function}
//
// Nothing has changed here.
int main()
{
  deallog.depth_console(2);

  Step3 laplace_problem;
  laplace_problem.run();

  return 0;
}

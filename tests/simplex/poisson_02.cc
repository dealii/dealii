// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Solve Poisson problem on a tet mesh with DG.

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

// #define HEX


template <int dim>
struct ScratchData
{
  ScratchData(const Mapping<dim>        &mapping,
              const FiniteElement<dim>  &fe,
              const Quadrature<dim>     &quad,
              const Quadrature<dim - 1> &quad_face,
              const UpdateFlags          update_flags = update_values |
                                               update_gradients |
                                               update_quadrature_points |
                                               update_JxW_values,
              const UpdateFlags interface_update_flags =
                update_values | update_gradients | update_quadrature_points |
                update_JxW_values | update_normal_vectors)
    : fe_values(mapping, fe, quad, update_flags)
    , fe_interface_values(mapping, fe, quad_face, interface_update_flags)
  {}


  ScratchData(const ScratchData<dim> &scratch_data)
    : fe_values(scratch_data.fe_values.get_mapping(),
                scratch_data.fe_values.get_fe(),
                scratch_data.fe_values.get_quadrature(),
                scratch_data.fe_values.get_update_flags())
    , fe_interface_values(scratch_data.fe_values.get_mapping(),
                          scratch_data.fe_values.get_fe(),
                          scratch_data.fe_interface_values.get_quadrature(),
                          scratch_data.fe_interface_values.get_update_flags())
  {}

  FEValues<dim>          fe_values;
  FEInterfaceValues<dim> fe_interface_values;
};



struct CopyDataFace
{
  FullMatrix<double>                   cell_matrix;
  std::vector<types::global_dof_index> joint_dof_indices;
};



struct CopyData
{
  FullMatrix<double>                   cell_matrix;
  Vector<double>                       cell_rhs;
  std::vector<types::global_dof_index> local_dof_indices;
  std::vector<CopyDataFace>            face_data;

  template <class Iterator>
  void
  reinit(const Iterator &cell, unsigned int dofs_per_cell)
  {
    cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
    cell_rhs.reinit(dofs_per_cell);

    local_dof_indices.resize(dofs_per_cell);
    cell->get_dof_indices(local_dof_indices);
  }
};

template <int dim>
class RightHandSideFunction : public Function<dim>
{
public:
  RightHandSideFunction()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int /*component*/ = 0) const
  {
    if (dim == 2)
      return -2. * M_PI * M_PI * std::sin(M_PI * p[0]) * std::sin(M_PI * p[1]);
    else /* if(dim == 3)*/
      return -3. * M_PI * M_PI * std::sin(M_PI * p[0]) * std::sin(M_PI * p[1]) *
             std::sin(M_PI * p[2]);
  }
};

template <int dim>
class DGHeat
{
public:
  DGHeat(const bool           hex,
         FiniteElement<dim>  *fe,
         Mapping<dim>        *mapping,
         Quadrature<dim>     *quad,
         Quadrature<dim - 1> *face_quad,
         unsigned int         initial_refinement,
         unsigned int         number_refinement)
    : hex(hex)
    , fe(fe)
    , mapping(mapping)
    , quad(quad)
    , face_quad(face_quad)
    , dof_handler(triangulation)
    , initial_refinement_level(initial_refinement)
    , number_refinement(number_refinement)
  {}

  static std::unique_ptr<DGHeat<dim>>
  HEX(unsigned int degree,
      unsigned int initial_refinement,
      unsigned int number_refinement)
  {
    return std::make_unique<DGHeat<dim>>(true,
                                         new FE_DGQ<dim>(degree),
                                         new MappingQ<dim>(1),
                                         new QGauss<dim>(degree + 1),
                                         new QGauss<dim - 1>(degree + 1),
                                         initial_refinement,
                                         number_refinement);
  }

  static std::unique_ptr<DGHeat<dim>>
  TET(unsigned int degree,
      unsigned int initial_refinement,
      unsigned int number_refinement)
  {
    return std::make_unique<DGHeat<dim>>(false,
                                         new FE_SimplexDGP<dim>(degree),
                                         new MappingFE<dim>(
                                           FE_SimplexP<dim>(1)),
                                         new QGaussSimplex<dim>(degree + 1),
                                         new QGaussSimplex<dim - 1>(degree + 1),
                                         initial_refinement,
                                         number_refinement);
  }



  void
  run();



private:
  void
  make_grid(int refinements = -1);
  void
  setup_system();
  void
  assemble_system();
  void
  solve();
  void
  output_results(unsigned int it) const;
  void
  calculateL2Error();

  Triangulation<dim> triangulation;

  bool hex;

  std::unique_ptr<FiniteElement<dim>>  fe;
  const std::unique_ptr<Mapping<dim>>  mapping;
  std::unique_ptr<Quadrature<dim>>     quad;
  std::unique_ptr<Quadrature<dim - 1>> face_quad;

  DoFHandler<dim> dof_handler;


  RightHandSideFunction<dim> right_hand_side;


  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> solution;
  Vector<double> system_rhs;
  Point<dim>     center;

  ConvergenceTable error_table;

  unsigned int initial_refinement_level;
  unsigned int number_refinement;
};


template <int dim>
void
DGHeat<dim>::make_grid(int refinements)
{
  triangulation.clear();

  const unsigned int ref =
    refinements == -1 ? initial_refinement_level : refinements;

  if (hex)
    GridGenerator::subdivided_hyper_cube(triangulation,
                                         Utilities::pow(2, ref),
                                         -1.0,
                                         +1.0);
  else
    GridGenerator::subdivided_hyper_cube_with_simplices(triangulation,
                                                        Utilities::pow(2, ref),
                                                        -1.0,
                                                        +1.0);

  // deallog << "   Number of active cells: " <<
  // triangulation.n_active_cells()
  //          << std::endl
  //          << "   Total number of cells: " << triangulation.n_cells()
  //          << std::endl;
}


template <int dim>
void
DGHeat<dim>::setup_system()
{
  dof_handler.distribute_dofs(*fe);

  // deallog << "   Number of degrees of freedom: " << dof_handler.n_dofs()
  //          << std::endl;

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}


template <int dim>
void
DGHeat<dim>::assemble_system()
{
  using Iterator = typename DoFHandler<dim>::active_cell_iterator;

  auto cell_worker = [&](const Iterator   &cell,
                         ScratchData<dim> &scratch_data,
                         CopyData         &copy_data) {
    const unsigned int n_dofs = scratch_data.fe_values.get_fe().dofs_per_cell;
    copy_data.reinit(cell, n_dofs);
    scratch_data.fe_values.reinit(cell);

    const auto &q_points = scratch_data.fe_values.get_quadrature_points();

    const FEValues<dim>       &fe_v = scratch_data.fe_values;
    const std::vector<double> &JxW  = fe_v.get_JxW_values();

    std::vector<double> f(q_points.size());
    right_hand_side.value_list(q_points, f);

    for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point)
      {
        for (unsigned int i = 0; i < n_dofs; ++i)
          {
            for (unsigned int j = 0; j < n_dofs; ++j)
              {
                copy_data.cell_matrix(i, j) +=
                  fe_v.shape_grad(i, point)   // \nabla \phi_i
                  * fe_v.shape_grad(j, point) // \nabla \phi_j
                  * JxW[point];               // dx
              }

            // Right Hand Side
            copy_data.cell_rhs(i) +=
              (fe_v.shape_value(i, point) * f[point] * JxW[point]);
          }
      }
  };

  auto boundary_worker = [&](const Iterator     &cell,
                             const unsigned int &face_no,
                             ScratchData<dim>   &scratch_data,
                             CopyData           &copy_data) {
    scratch_data.fe_interface_values.reinit(cell, face_no);

    const FEFaceValuesBase<dim> &fe_face =
      scratch_data.fe_interface_values.get_fe_face_values(0);

    const auto        &q_points     = fe_face.get_quadrature_points();
    const unsigned int n_facet_dofs = fe_face.get_fe().n_dofs_per_cell();
    const std::vector<double> &JxW  = fe_face.get_JxW_values();

    const std::vector<Tensor<1, dim>> &normals = fe_face.get_normal_vectors();

    double h;
    if (dim == 2)
      {
        if (hex)
          h = std::sqrt(4. * cell->measure() / M_PI);
        else
          h = std::sqrt(4. * (4.0 / triangulation.n_cells()) / M_PI);
      }
    else if (dim == 3)
      {
        if (hex)
          h = pow(6 * cell->measure() / M_PI, 1. / 3.);
        else
          h = pow(6 * (8.0 / triangulation.n_cells()) / M_PI, 1. / 3.);
      }



    const double beta = 10.;

    for (unsigned int point = 0; point < q_points.size(); ++point)
      for (unsigned int i = 0; i < n_facet_dofs; ++i)
        for (unsigned int j = 0; j < n_facet_dofs; ++j)
          {
            copy_data.cell_matrix(i, j) +=
              -normals[point] * fe_face.shape_grad(i, point) // n*\nabla \phi_i
              * fe_face.shape_value(j, point)                // \phi_j
              * JxW[point];                                  // dx

            copy_data.cell_matrix(i, j) +=
              -fe_face.shape_value(i, point)                  // \phi_i
              * fe_face.shape_grad(j, point) * normals[point] // n*\nabla \phi_j
              * JxW[point];                                   // dx

            copy_data.cell_matrix(i, j) +=
              beta * 1. / h * fe_face.shape_value(i, point) // \phi_i
              * fe_face.shape_value(j, point) * JxW[point]; // dx
          }
  };

  auto face_worker = [&](const Iterator     &cell,
                         const unsigned int &f,
                         const unsigned int &sf,
                         const Iterator     &ncell,
                         const unsigned int &nf,
                         const unsigned int &nsf,
                         ScratchData<dim>   &scratch_data,
                         CopyData           &copy_data) {
    FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values;

    fe_iv.reinit(cell, f, sf, ncell, nf, nsf);

    const auto &q_points = fe_iv.get_quadrature_points();

    copy_data.face_data.emplace_back();
    CopyDataFace &copy_data_face = copy_data.face_data.back();

    const unsigned int n_dofs        = fe_iv.n_current_interface_dofs();
    copy_data_face.joint_dof_indices = fe_iv.get_interface_dof_indices();

    copy_data_face.cell_matrix.reinit(n_dofs, n_dofs);

    const std::vector<double>         &JxW     = fe_iv.get_JxW_values();
    const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();


    double h;
    if (dim == 2)
      {
        if (hex)
          h = std::sqrt(4. * cell->measure() / M_PI);
        else
          h = std::sqrt(4. * (4.0 / triangulation.n_cells()) / M_PI);
      }
    else if (dim == 3)
      {
        if (hex)
          h = pow(6 * cell->measure() / M_PI, 1. / 3.);
        else
          h = pow(6 * (8.0 / triangulation.n_cells()) / M_PI, 1. / 3.);
      }

    const double beta = 10.;

    for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
      {
        for (unsigned int i = 0; i < n_dofs; ++i)
          {
            for (unsigned int j = 0; j < n_dofs; ++j)
              {
                copy_data_face.cell_matrix(i, j) +=
                  -normals[qpoint] *
                  fe_iv.average_of_shape_gradients(i, qpoint) *
                  fe_iv.jump_in_shape_values(j, qpoint) * JxW[qpoint];

                copy_data_face.cell_matrix(i, j) +=
                  -fe_iv.jump_in_shape_values(i, qpoint) // \phi_i
                  * fe_iv.average_of_shape_gradients(j, qpoint) *
                  normals[qpoint] // n*\nabla \phi_j
                  * JxW[qpoint];  // dx

                copy_data_face.cell_matrix(i, j) +=
                  beta * 1. / h * fe_iv.jump_in_shape_values(i, qpoint) *
                  fe_iv.jump_in_shape_values(j, qpoint) * JxW[qpoint];
              }
          }
      }
  };

  AffineConstraints<double> constraints;

  auto copier = [&](const CopyData &c) {
    constraints.distribute_local_to_global(c.cell_matrix,
                                           c.cell_rhs,
                                           c.local_dof_indices,
                                           system_matrix,
                                           system_rhs);

    for (auto &cdf : c.face_data)
      {
        constraints.distribute_local_to_global(cdf.cell_matrix,
                                               cdf.joint_dof_indices,
                                               system_matrix);
      }
  };


  ScratchData<dim> scratch_data(*mapping, *fe, *quad, *face_quad);
  CopyData         copy_data;

  MeshWorker::mesh_loop(dof_handler.begin_active(),
                        dof_handler.end(),
                        cell_worker,
                        copier,
                        scratch_data,
                        copy_data,
                        MeshWorker::assemble_own_cells |
                          MeshWorker::assemble_boundary_faces |
                          MeshWorker::assemble_own_interior_faces_once,
                        boundary_worker,
                        face_worker);
}

template <int dim>
void
DGHeat<dim>::solve()
{
  SolverControl solver_control(10000, 1e-8);
  SolverCG<>    solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());

  // We have made one addition, though: since we suppress output from the
  // linear solvers, we have to print the number of iterations by hand.
  // deallog << "   " << solver_control.last_step()
  //          << " CG iterations needed to obtain convergence." << std::endl;

  // error_table.add_value("iterations", solver_control.last_step());
}

template <int dim>
void
DGHeat<dim>::output_results(unsigned int it) const
{
  return;

  std::string type = hex ? "hex" : "tet";

  std::string dimension(dim == 2 ? "solution-2d-" + type + "-case-" :
                                   "solution-3d-" + type + "-case-");

  std::string fname = dimension + Utilities::int_to_string(it) + ".vtk";

  deallog << "  Writing solution to <" << fname << '>' << std::endl;

  std::ofstream output(fname);

  if (false)
    {
      DataOut<dim> data_out;

      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(solution, "solution");

      data_out.build_patches(*mapping);
      data_out.write_vtk(output);
    }
}

// Find the l2 norm of the error between the finite element sol'n and the exact
// sol'n
template <int dim>
void
DGHeat<dim>::calculateL2Error()
{
  FEValues<dim> fe_values(*mapping,
                          *fe,
                          *quad,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell =
    fe->dofs_per_cell; // This gives you dofs per cell
  std::vector<types::global_dof_index> local_dof_indices(
    dofs_per_cell); //  Local connectivity

  const unsigned int n_q_points = quad->size();

  double l2error = 0.;

  // loop over elements
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);

      cell->get_dof_indices(local_dof_indices);

      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          const double u_exact =
            dim == 2 ? -std::sin(M_PI * fe_values.quadrature_point(q)[0]) *
                         std::sin(M_PI * fe_values.quadrature_point(q)[1]) :
                       -std::sin(M_PI * fe_values.quadrature_point(q)[0]) *
                         std::sin(M_PI * fe_values.quadrature_point(q)[1]) *
                         std::sin(M_PI * fe_values.quadrature_point(q)[2]);

          double u_sim = 0;

          // Find the values of x and u_h (the finite element solution) at the
          // quadrature points
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              u_sim +=
                fe_values.shape_value(i, q) * solution[local_dof_indices[i]];
            }
          l2error += (u_sim - u_exact) * (u_sim - u_exact) * fe_values.JxW(q);
          //       deallog << " x = " << x << " y = " << y <<  " r = " << r <<
          //       "   u_exact = " << u_exact << "   u_sim=" << u_sim <<
          //       std::endl;
        }
    }


  // deallog << "L2Error is : " << std::sqrt(l2error) << std::endl;
  error_table.add_value("error", std::sqrt(l2error));
  error_table.add_value("cells", triangulation.n_global_active_cells());
  error_table.add_value("dofs", dof_handler.n_dofs());
}



template <int dim>
void
DGHeat<dim>::run()
{
  for (unsigned int it = 0; it < number_refinement; ++it)
    {
      make_grid(initial_refinement_level + it);
      setup_system();
      assemble_system();
      solve();
      output_results(it);
      calculateL2Error();
    }

  // error_table.omit_column_from_convergence_rate_evaluation("iterations");
  error_table.omit_column_from_convergence_rate_evaluation("cells");
  error_table.evaluate_all_convergence_rates(
    ConvergenceTable::reduction_rate_log2);

  error_table.set_scientific("error", true);

  error_table.write_text(deallog.get_file_stream());
  deallog << std::endl;
}

int
main()
{
  initlog();

  deallog.depth_file(1);


  {
    auto problem = DGHeat<2>::TET(1 /*=degree*/, 2, 3);
    problem->run();
  }
  {
    auto problem = DGHeat<2>::TET(2 /*=degree*/, 2, 3);
    problem->run();
  }
  {
    auto problem = DGHeat<3>::TET(1 /*=degree*/, 2, 2);
    problem->run();
  }
  {
    auto problem = DGHeat<3>::TET(2 /*=degree*/, 2, 2);
    problem->run();
  }

  {
    auto problem = DGHeat<2>::HEX(1 /*=degree*/, 2, 3);
    problem->run();
  }
  {
    auto problem = DGHeat<2>::HEX(2 /*=degree*/, 2, 3);
    problem->run();
  }
  {
    auto problem = DGHeat<3>::HEX(1 /*=degree*/, 2, 2);
    problem->run();
  }
  {
    auto problem = DGHeat<3>::HEX(2 /*=degree*/, 2, 2);
    problem->run();
  }

  return 0;
}

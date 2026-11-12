// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2024 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// solves a 2D Poisson equation for linear FE_DGP elements (SIP
// discretization) with MueLu preconditioner

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/trilinos_tpetra_precondition.h>
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



// Symmetric interior penalty discretization of the Laplacian, assembled
// with mesh_loop and FEInterfaceValues.
template <int dim>
struct ScratchData
{
  ScratchData(const Mapping<dim>        &mapping,
              const FiniteElement<dim>  &fe,
              const Quadrature<dim>     &quadrature,
              const Quadrature<dim - 1> &face_quadrature)
    : mapping(&mapping)
    , fe(&fe)
    , quadrature(quadrature)
    , face_quadrature(face_quadrature)
    , fe_values(mapping, fe, quadrature, update_gradients | update_JxW_values)
    , fe_face_values(mapping,
                     fe,
                     face_quadrature,
                     update_values | update_gradients | update_JxW_values |
                       update_normal_vectors)
    , fe_interface_values(mapping,
                          fe,
                          face_quadrature,
                          update_values | update_gradients | update_JxW_values |
                            update_normal_vectors)
  {}

  ScratchData(const ScratchData &other)
    : ScratchData(*other.mapping,
                  *other.fe,
                  other.quadrature,
                  other.face_quadrature)
  {}

  const Mapping<dim>       *mapping;
  const FiniteElement<dim> *fe;
  Quadrature<dim>           quadrature;
  Quadrature<dim - 1>       face_quadrature;
  FEValues<dim>             fe_values;
  FEFaceValues<dim>         fe_face_values;
  FEInterfaceValues<dim>    fe_interface_values;
};



struct CopyData
{
  struct FaceData
  {
    FullMatrix<double>                   matrix;
    std::vector<types::global_dof_index> dof_indices;
  };

  FullMatrix<double>                   cell_matrix;
  std::vector<types::global_dof_index> cell_dof_indices;
  std::vector<FaceData>                face_data;
};



template <int dim>
double
penalty_in_direction(const typename DoFHandler<dim>::active_cell_iterator &cell,
                     const unsigned int face_no,
                     const unsigned int degree)
{
  const unsigned int deg_sq = (degree == 0) ? 1 : degree * (degree + 1);
  return deg_sq / cell->extent_in_direction(
                    GeometryInfo<dim>::unit_normal_direction[face_no]);
}



template <int dim, typename MatrixType>
void
assemble_laplace_matrix(const DoFHandler<dim> &dof_handler,
                        const Mapping<dim>    &mapping,
                        MatrixType            &matrix)
{
  const FiniteElement<dim> &fe  = dof_handler.get_fe();
  const unsigned int        deg = fe.degree;
  using Iterator = typename DoFHandler<dim>::active_cell_iterator;

  const auto cell_worker =
    [&](const Iterator &cell, ScratchData<dim> &scratch, CopyData &copy) {
      scratch.fe_values.reinit(cell);
      const unsigned int n_dofs = fe.n_dofs_per_cell();
      copy.cell_matrix.reinit(n_dofs, n_dofs);
      for (unsigned int q = 0; q < scratch.fe_values.n_quadrature_points; ++q)
        for (unsigned int i = 0; i < n_dofs; ++i)
          for (unsigned int j = 0; j < n_dofs; ++j)
            copy.cell_matrix(i, j) += scratch.fe_values.shape_grad(i, q) *
                                      scratch.fe_values.shape_grad(j, q) *
                                      scratch.fe_values.JxW(q);
      copy.cell_dof_indices.resize(n_dofs);
      cell->get_dof_indices(copy.cell_dof_indices);
    };

  // Nitsche boundary terms: 2 gamma u v - d_n u v - u d_n v
  const auto boundary_worker = [&](const Iterator    &cell,
                                   const unsigned int face_no,
                                   ScratchData<dim>  &scratch,
                                   CopyData          &copy) {
    scratch.fe_face_values.reinit(cell, face_no);
    const auto        &fev     = scratch.fe_face_values;
    const unsigned int n_dofs  = fe.n_dofs_per_cell();
    const double       penalty = penalty_in_direction<dim>(cell, face_no, deg);
    auto              &face_data = copy.face_data.emplace_back();
    face_data.matrix.reinit(n_dofs, n_dofs);
    face_data.dof_indices.resize(n_dofs);
    cell->get_dof_indices(face_data.dof_indices);
    for (unsigned int q = 0; q < fev.n_quadrature_points; ++q)
      for (unsigned int i = 0; i < n_dofs; ++i)
        for (unsigned int j = 0; j < n_dofs; ++j)
          face_data.matrix(i, j) +=
            (2. * penalty * fev.shape_value(i, q) * fev.shape_value(j, q) -
             (fev.normal_vector(q) * fev.shape_grad(i, q)) *
               fev.shape_value(j, q) -
             (fev.normal_vector(q) * fev.shape_grad(j, q)) *
               fev.shape_value(i, q)) *
            fev.JxW(q);
  };

  // Interior penalty: gamma [u][v] - {d_n u}[v] - [u]{d_n v}
  const auto face_worker = [&](const Iterator    &cell,
                               const unsigned int f,
                               const unsigned int sf,
                               const Iterator    &ncell,
                               const unsigned int nf,
                               const unsigned int nsf,
                               ScratchData<dim>  &scratch,
                               CopyData          &copy) {
    auto &fiv = scratch.fe_interface_values;
    fiv.reinit(cell, f, sf, ncell, nf, nsf);

    const double penalty = 0.5 * (penalty_in_direction<dim>(cell, f, deg) +
                                  penalty_in_direction<dim>(ncell, nf, deg));

    const unsigned int n_interface_dofs = fiv.n_current_interface_dofs();
    auto              &face_data        = copy.face_data.emplace_back();
    face_data.matrix.reinit(n_interface_dofs, n_interface_dofs);
    face_data.dof_indices = fiv.get_interface_dof_indices();

    const std::vector<Tensor<1, dim>> &normals = fiv.get_normal_vectors();
    for (unsigned int q = 0; q < fiv.n_quadrature_points; ++q)
      for (unsigned int i = 0; i < n_interface_dofs; ++i)
        for (unsigned int j = 0; j < n_interface_dofs; ++j)
          face_data.matrix(i, j) +=
            (penalty * fiv.jump_in_shape_values(i, q) *
               fiv.jump_in_shape_values(j, q) -
             (fiv.average_of_shape_gradients(j, q) * normals[q]) *
               fiv.jump_in_shape_values(i, q) -
             (fiv.average_of_shape_gradients(i, q) * normals[q]) *
               fiv.jump_in_shape_values(j, q)) *
            fiv.JxW(q);
  };

  const auto copier = [&](const CopyData &copy) {
    for (unsigned int i = 0; i < copy.cell_dof_indices.size(); ++i)
      for (unsigned int j = 0; j < copy.cell_dof_indices.size(); ++j)
        matrix.add(copy.cell_dof_indices[i],
                   copy.cell_dof_indices[j],
                   copy.cell_matrix(i, j));

    for (const auto &face : copy.face_data)
      for (unsigned int i = 0; i < face.dof_indices.size(); ++i)
        for (unsigned int j = 0; j < face.dof_indices.size(); ++j)
          matrix.add(face.dof_indices[i],
                     face.dof_indices[j],
                     face.matrix(i, j));
  };

  ScratchData<dim> scratch(mapping,
                           fe,
                           QGauss<dim>(deg + 1),
                           QGauss<dim - 1>(deg + 1));
  CopyData         copy;
  MeshWorker::mesh_loop(dof_handler.begin_active(),
                        dof_handler.end(),
                        cell_worker,
                        copier,
                        scratch,
                        copy,
                        MeshWorker::assemble_own_cells |
                          MeshWorker::assemble_boundary_faces |
                          MeshWorker::assemble_own_interior_faces_once,
                        boundary_worker,
                        face_worker);
}


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
  solve();

  Triangulation<dim> triangulation;
  FE_DGP<dim>        fe;
  DoFHandler<dim>    dof_handler;

  LinearAlgebra::TpetraWrappers::SparseMatrix<double> system_matrix;

  LinearAlgebra::TpetraWrappers::Vector<double> solution;
  LinearAlgebra::TpetraWrappers::Vector<double> system_rhs;
};



template <int dim>
Step4<dim>::Step4()
  : fe(1)
  , dof_handler(triangulation)
{}


template <int dim>
void
Step4<dim>::make_grid()
{
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(4);
}



template <int dim>
void
Step4<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  DynamicSparsityPattern c_sparsity(dof_handler.n_dofs());
  DoFTools::make_flux_sparsity_pattern(dof_handler, c_sparsity);
  system_matrix.reinit(c_sparsity);

  solution.reinit(dof_handler.locally_owned_dofs());
  system_rhs.reinit(dof_handler.locally_owned_dofs());

  MappingQ<dim> mapping(1);
  assemble_laplace_matrix(dof_handler, mapping, system_matrix);

  system_matrix.compress(VectorOperation::add);

  for (unsigned int i = 0; i < system_rhs.size(); ++i)
    system_rhs(i) = 0.01 * i - 0.000001 * i * i;
}



template <int dim>
void
Step4<dim>::solve()
{
  deallog.push(Utilities::int_to_string(dof_handler.n_dofs(), 5));
  LinearAlgebra::TpetraWrappers::PreconditionAMGMueLu<double> preconditioner;
  LinearAlgebra::TpetraWrappers::PreconditionAMGMueLu<double>::AdditionalData
    data;
  data.constant_modes =
    DoFTools::extract_constant_modes(dof_handler, ComponentMask(1, true));
  data.smoother_sweeps = 2;
  {
    solution = 0;
    SolverControl solver_control(1000, 1e-10);
    SolverCG<LinearAlgebra::TpetraWrappers::Vector<double>> solver(
      solver_control);
    preconditioner.initialize(system_matrix, data);

    check_solver_within_range(
      solver.solve(system_matrix, solution, system_rhs, preconditioner),
      solver_control.last_step(),
      18,
      32);
  }
  deallog.pop();
}



template <int dim>
void
Step4<dim>::run()
{
  for (unsigned int cycle = 0; cycle < 2; ++cycle)
    {
      if (cycle == 0)
        make_grid();
      else
        triangulation.refine_global(1);

      setup_system();
      solve();
    }
}


int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  try
    {
      Step4<2> test;
      test.run();
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}

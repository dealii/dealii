// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

#include <deal.II/base/function_lib.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <deal.II/non_matching/coupling.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


// Test that a coupling matrix can be constructed for each pair of dimension
// and immersed dimension, and check that quadratic functions are correctly
// projected.

template <int dim, int spacedim>
void
test()
{
  deallog << "dim: " << dim << ", spacedim: " << spacedim << std::endl;

  const auto &comm = MPI_COMM_WORLD;

  parallel::shared::Triangulation<dim, spacedim>           tria(comm);
  parallel::distributed::Triangulation<spacedim, spacedim> space_tria(comm);

  GridGenerator::hyper_cube(tria, -.4, .3);
  GridGenerator::hyper_cube(space_tria, -1, 1);

  tria.refine_global(3);
  space_tria.refine_global(3);

  FE_Q<dim, spacedim>      fe(2);
  FE_Q<spacedim, spacedim> space_fe(2);

  deallog << "FE      : " << fe.get_name() << std::endl
          << "Space FE: " << space_fe.get_name() << std::endl;

  DoFHandler<dim, spacedim>      dh(tria);
  DoFHandler<spacedim, spacedim> space_dh(space_tria);

  dh.distribute_dofs(fe);
  space_dh.distribute_dofs(space_fe);

  auto space_locally_owned_dofs = space_dh.locally_owned_dofs();
  auto locally_owned_dofs       = dh.locally_owned_dofs();


  deallog << "Dofs      : " << dh.n_dofs() << std::endl
          << "Space dofs: " << space_dh.n_dofs() << std::endl;

  deallog << "Local dofs      : " << locally_owned_dofs.n_elements()
          << std::endl;
  deallog << "Local space dofs: " << space_locally_owned_dofs.n_elements()
          << std::endl;
  QGauss<dim> quad(3); // Quadrature for coupling

  GridTools::Cache<spacedim, spacedim> cache(space_tria);

  TrilinosWrappers::SparsityPattern sparsity(space_locally_owned_dofs,
                                             locally_owned_dofs,
                                             comm);
  NonMatching::create_coupling_sparsity_pattern(
    cache, space_dh, dh, quad, sparsity);
  sparsity.compress();

  TrilinosWrappers::SparseMatrix coupling(sparsity);
  NonMatching::create_coupling_mass_matrix(cache, space_dh, dh, quad, coupling);
  coupling.compress(VectorOperation::add);

  TrilinosWrappers::SparsityPattern mass_sparsity(locally_owned_dofs, comm);
  DoFTools::make_sparsity_pattern(dh, mass_sparsity);
  mass_sparsity.compress();

  TrilinosWrappers::SparseMatrix mass_matrix(mass_sparsity);

  {
    QGauss<dim>             quad(4);
    FEValues<dim, spacedim> fev(fe, quad, update_values | update_JxW_values);
    std::vector<types::global_dof_index> dofs(fe.dofs_per_cell);
    FullMatrix<double>        cell_matrix(fe.dofs_per_cell, fe.dofs_per_cell);
    AffineConstraints<double> constraints;

    for (auto &cell : dh.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          cell_matrix = 0;
          fev.reinit(cell);
          cell->get_dof_indices(dofs);
          for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
            for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
              for (unsigned int q = 0; q < quad.size(); ++q)
                cell_matrix(i, j) +=
                  fev.shape_value(i, q) * fev.shape_value(j, q) * fev.JxW(q);
          constraints.distribute_local_to_global(cell_matrix,
                                                 dofs,
                                                 mass_matrix);
        }
    mass_matrix.compress(VectorOperation::add);
  }

  // now take the square function in space, project them onto the immersed
  // space, get back ones, and check for the error.
  TrilinosWrappers::MPI::Vector space_square(space_locally_owned_dofs, comm);
  TrilinosWrappers::MPI::Vector squares(locally_owned_dofs, comm);
  TrilinosWrappers::MPI::Vector Mprojected_squares(locally_owned_dofs, comm);
  TrilinosWrappers::MPI::Vector projected_squares(locally_owned_dofs, comm);

  VectorTools::interpolate(space_dh,
                           Functions::SquareFunction<spacedim>(),
                           space_square);
  VectorTools::interpolate(dh, Functions::SquareFunction<spacedim>(), squares);

  coupling.Tvmult(Mprojected_squares, space_square);

  SolverControl                     cn(100, 1e-12, false, false);
  TrilinosWrappers::SolverCG        solver(cn);
  TrilinosWrappers::PreconditionILU prec;
  prec.initialize(mass_matrix);

  solver.solve(mass_matrix, projected_squares, Mprojected_squares, prec);

  deallog << "Squares norm    : " << projected_squares.l2_norm() << std::endl;

  projected_squares -= squares;

  deallog << "Error on squares: " << projected_squares.l2_norm() << std::endl;
}



int
main(int argc, char **argv)
{
  auto          init = Utilities::MPI::MPI_InitFinalize(argc, argv, 1);
  MPILogInitAll log(true);
  test<1, 2>();
  test<2, 2>();
  test<2, 3>();
  test<3, 3>();
}

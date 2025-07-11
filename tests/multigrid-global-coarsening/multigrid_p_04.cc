// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/**
 * Similar to test multigrid_01 but the DoFs are renumbered in the DoFHandler
 * used by the outer solver.
 */

#include <deal.II/dofs/dof_renumbering.h>

#include "multigrid_util.h"

template <int dim, typename Number = double>
void
test(const unsigned int n_refinements,
     const unsigned int fe_degree_fine,
     const bool         do_simplex_mesh,
     const unsigned int mesh_type)
{
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  Triangulation<dim> tria;
  if (do_simplex_mesh)
    GridGenerator::subdivided_hyper_cube_with_simplices(tria, 2);
  else
    GridGenerator::subdivided_hyper_cube(tria, 2);

  if (mesh_type == 0)
    {
      tria.refine_global(n_refinements);
    }
  else if (mesh_type == 1)
    {
      for (unsigned int i = 1; i < n_refinements; ++i)
        {
          for (auto cell : tria.active_cell_iterators())
            if (cell->is_locally_owned())
              {
                bool flag = true;
                for (int d = 0; d < dim; ++d)
                  if (cell->center()[d] > 0.5)
                    flag = false;
                if (flag)
                  cell->set_refine_flag();
              }
          tria.execute_coarsening_and_refinement();
        }
    }
  else
    AssertThrow(false, ExcNotImplemented());

  ////////////////////////////////////////////////////////////////////////

  std::unique_ptr<FiniteElement<dim>> fe;
  std::unique_ptr<Quadrature<dim>>    quad;
  std::unique_ptr<Mapping<dim>>       mapping;

  if (do_simplex_mesh)
    {
      fe      = std::make_unique<FE_SimplexP<dim>>(fe_degree_fine);
      quad    = std::make_unique<QGaussSimplex<dim>>(fe_degree_fine + 1);
      mapping = std::make_unique<MappingFE<dim>>(FE_SimplexP<dim>(1));
    }
  else
    {
      fe      = std::make_unique<FE_Q<dim>>(fe_degree_fine);
      quad    = std::make_unique<QGauss<dim>>(fe_degree_fine + 1);
      mapping = std::make_unique<MappingFE<dim>>(FE_Q<dim>(1));
    }

  DoFHandler<dim> fine_dof_handler(tria);
  fine_dof_handler.distribute_dofs(*fe);
  DoFRenumbering::Cuthill_McKee(fine_dof_handler);

  AffineConstraints<Number> fine_constraints;
  const IndexSet            locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(fine_dof_handler);
  fine_constraints.reinit(fine_dof_handler.locally_owned_dofs(),
                          locally_relevant_dofs);
  VectorTools::interpolate_boundary_values(*mapping,
                                           fine_dof_handler,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           fine_constraints);
  DoFTools::make_hanging_node_constraints(fine_dof_handler, fine_constraints);
  fine_constraints.close();

  // set up operator
  Operator<dim, 1, Number> fine_operator;
  fine_operator.reinit(*mapping, fine_dof_handler, *quad, fine_constraints);

  ////////////////////////////////////////////////////////////////////////

  const auto level_degrees =
    MGTransferGlobalCoarseningTools::create_polynomial_coarsening_sequence(
      fe_degree_fine,
      MGTransferGlobalCoarseningTools::PolynomialCoarseningSequenceType::
        bisect);

  const unsigned int min_level = 0;
  const unsigned int max_level = level_degrees.size() - 1;

  MGLevelObject<DoFHandler<dim>> dof_handlers(min_level, max_level, tria);
  MGLevelObject<AffineConstraints<Number>> constraints(min_level, max_level);
  MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> transfers(min_level,
                                                               max_level);
  MGLevelObject<Operator<dim, 1, Number>> operators(min_level, max_level);

  std::unique_ptr<Mapping<dim>> mapping_;

  // set up levels
  for (auto l = min_level; l <= max_level; ++l)
    {
      auto &dof_handler = dof_handlers[l];
      auto &constraint  = constraints[l];
      auto &op          = operators[l];

      std::unique_ptr<FiniteElement<dim>> fe;
      std::unique_ptr<Quadrature<dim>>    quad;
      std::unique_ptr<Mapping<dim>>       mapping;

      if (do_simplex_mesh)
        {
          fe      = std::make_unique<FE_SimplexP<dim>>(level_degrees[l]);
          quad    = std::make_unique<QGaussSimplex<dim>>(level_degrees[l] + 1);
          mapping = std::make_unique<MappingFE<dim>>(FE_SimplexP<dim>(1));
        }
      else
        {
          fe      = std::make_unique<FE_Q<dim>>(level_degrees[l]);
          quad    = std::make_unique<QGauss<dim>>(level_degrees[l] + 1);
          mapping = std::make_unique<MappingFE<dim>>(FE_Q<dim>(1));
        }

      if (l == max_level)
        mapping_ = mapping->clone();

      // set up dofhandler
      dof_handler.distribute_dofs(*fe);

      // set up constraints
      const IndexSet locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(dof_handler);
      constraint.reinit(dof_handler.locally_owned_dofs(),
                        locally_relevant_dofs);
      VectorTools::interpolate_boundary_values(
        *mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraint);
      DoFTools::make_hanging_node_constraints(dof_handler, constraint);
      constraint.close();

      // set up operator
      op.reinit(*mapping, dof_handler, *quad, constraint);
    }

  // set up transfer operator
  for (unsigned int l = min_level; l < max_level; ++l)
    transfers[l + 1].reinit(dof_handlers[l + 1],
                            dof_handlers[l],
                            constraints[l + 1],
                            constraints[l]);

  MGTransferGlobalCoarsening<dim, VectorType> transfer;
  transfer.initialize_two_level_transfers(transfers);
  transfer.build(fine_dof_handler, [&](const auto l, auto &vec) {
    operators[l].initialize_dof_vector(vec);
  });

  GMGParameters mg_data; // TODO

  VectorType dst, src;
  fine_operator.initialize_dof_vector(dst);
  fine_operator.initialize_dof_vector(src);
  fine_operator.rhs(src);

  ReductionControl solver_control(
    mg_data.maxiter, mg_data.abstol, mg_data.reltol, false, false);

  mg_solve(solver_control,
           dst,
           src,
           mg_data,
           fine_dof_handler,
           fine_operator,
           operators,
           transfer);

  fine_constraints.distribute(dst);

  deallog << dim << ' ' << fe_degree_fine << ' ' << n_refinements << ' '
          << (do_simplex_mesh ? "tri " : "quad") << ' '
          << solver_control.last_step() << std::endl;

  static unsigned int counter = 0;

  MGLevelObject<VectorType> results(min_level, max_level);

  transfer.interpolate_to_mg(dof_handlers[max_level], results, dst);

  for (unsigned int l = min_level; l <= max_level; ++l)
    {
      deallog << "Norm interpolated solution on level " << l << ": "
              << results[l].l2_norm() << std::endl;
      DataOut<dim> data_out;

      data_out.attach_dof_handler(dof_handlers[l]);
      data_out.add_data_vector(
        results[l],
        "solution",
        DataOut_DoFData<dim, dim>::DataVectorType::type_dof_data);
      data_out.build_patches(*mapping_, 2);

      std::ofstream output("test." + std::to_string(dim) + "." +
                           std::to_string(counter) + "." + std::to_string(l) +
                           ".vtk");
      data_out.write_vtk(output);
    }

  counter++;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.precision(8);

  for (unsigned int mesh_type = 0; mesh_type < 2; ++mesh_type)
    {
      for (unsigned int n_refinements = 2; n_refinements <= 4; ++n_refinements)
        for (unsigned int degree = 2; degree <= 4; ++degree)
          test<2>(n_refinements, degree, false /*quadrilateral*/, mesh_type);

      for (unsigned int n_refinements = 2; n_refinements <= 4; ++n_refinements)
        for (unsigned int degree = 2; degree <= 2; ++degree)
          test<2>(n_refinements, degree, true /*triangle*/, mesh_type);
    }
}

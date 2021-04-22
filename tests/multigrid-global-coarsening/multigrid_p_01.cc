// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


/**
 * Test p-multigrid for a uniformly and locally refined meshes both for
 * simplex and hypercube mesh.
 */

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
      for (unsigned int i = 1; i < n_refinements; i++)
        {
          for (auto cell : tria.active_cell_iterators())
            if (cell->is_locally_owned())
              {
                bool flag = true;
                for (int d = 0; d < dim; d++)
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
  MGLevelObject<Operator<dim, Number>> operators(min_level, max_level);

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
      IndexSet locally_relevant_dofs;
      DoFTools::extract_locally_relevant_dofs(dof_handler,
                                              locally_relevant_dofs);
      constraint.reinit(locally_relevant_dofs);
      VectorTools::interpolate_boundary_values(
        *mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraint);
      DoFTools::make_hanging_node_constraints(dof_handler, constraint);
      constraint.close();

      // set up operator
      op.reinit(*mapping, dof_handler, *quad, constraint);
    }

  // set up transfer operator
  for (unsigned int l = min_level; l < max_level; ++l)
    transfers[l + 1].reinit_polynomial_transfer(dof_handlers[l + 1],
                                                dof_handlers[l],
                                                constraints[l + 1],
                                                constraints[l]);

  MGTransferGlobalCoarsening<dim, VectorType> transfer(
    transfers,
    [&](const auto l, auto &vec) { operators[l].initialize_dof_vector(vec); });

  GMGParameters mg_data; // TODO

  VectorType dst, src;
  operators[max_level].initialize_dof_vector(dst);
  operators[max_level].initialize_dof_vector(src);

  operators[max_level].rhs(src);

  ReductionControl solver_control(
    mg_data.maxiter, mg_data.abstol, mg_data.reltol, false, false);

  mg_solve(solver_control,
           dst,
           src,
           mg_data,
           dof_handlers[max_level],
           operators[max_level],
           operators,
           transfer);

  constraints[max_level].distribute(dst);

  deallog << dim << " " << fe_degree_fine << " " << n_refinements << " "
          << (do_simplex_mesh ? "tri " : "quad") << " "
          << solver_control.last_step() << std::endl;

  static unsigned int counter = 0;

  MGLevelObject<VectorType> results(min_level, max_level);

  transfer.interpolate_to_mg(dof_handlers[max_level], results, dst);

  for (unsigned int l = min_level; l <= max_level; ++l)
    {
      DataOut<dim> data_out;

      data_out.attach_dof_handler(dof_handlers[l]);
      data_out.add_data_vector(
        results[l],
        "solution",
        DataOut_DoFData<DoFHandler<dim>, dim>::DataVectorType::type_dof_data);
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

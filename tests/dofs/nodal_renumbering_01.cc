// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/function_lib.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec_sz.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools_interpolate.h>

#include <iostream>
#include <map>

#include "../tests.h"

#include "../simplex/simplex_grids.h"

// Test DoFRenumbering::support_point_wise(). The vector position should contain
// (x0, y0, z0, x1, y1, z1, ...) and solution1 should contain (u0, v0, w0, u1,
// v1, w1, ...). Verify that they match when we interpolate.


class Test : public Function<2>
{
public:
  Test(const unsigned int n_components)
    : Function<2>(n_components)
  {}

  double
  value(const Point<2> &p, const unsigned int component = 0) const override
  {
    switch (component)
      {
        case 0:
          return std::sin(p[0]) * std::sin(2.0 * p[1]) + 1.0;
        case 1:
          return std::cos(p[0]) * std::cos(3.0 * p[1]) + 2.0;
        default:
          DEAL_II_NOT_IMPLEMENTED();
      }
    return 0.0;
  }
};

void
test(DoFHandler<2> &dof_handler, const hp::MappingCollection<2> &mappings)
{
  DoFRenumbering::support_point_wise(dof_handler);
  const MPI_Comm comm = dof_handler.get_mpi_communicator();

  const IndexSet &local_dofs = dof_handler.locally_owned_dofs();
  deallog << "new case with locally owned dofs = ";
  local_dofs.print(deallog);
  deallog << std::endl;
  const IndexSet relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);
  LinearAlgebra::distributed::Vector<double> position(local_dofs,
                                                      relevant_dofs,
                                                      comm);
  LinearAlgebra::distributed::Vector<double> solution1(local_dofs,
                                                       relevant_dofs,
                                                       comm);
  LinearAlgebra::distributed::Vector<double> solution2(local_dofs,
                                                       relevant_dofs,
                                                       comm);

  const unsigned int n_components =
    dof_handler.get_fe_collection().n_components();
  // It doesn't make sense to treat position as a vector of coordinates unless
  // we have enough components
  if (n_components == 2)
    {
      VectorTools::interpolate(mappings,
                               dof_handler,
                               Functions::IdentityFunction<2>(),
                               position);

      Test test(n_components);
      VectorTools::interpolate(mappings, dof_handler, test, solution1);

      // Verify that we get the same results when we interpolate either manually
      // or by reading off position data and evaluating the interpolated
      // function
      for (unsigned int node_n = 0;
           node_n < position.locally_owned_size() / n_components;
           ++node_n)
        {
          const auto i0 = node_n * n_components;
          const auto i1 = node_n * n_components + 1;
          Point<2>   p(position.local_element(i0), position.local_element(i1));
          solution2.local_element(i0) = test.value(p, 0);
          solution2.local_element(i1) = test.value(p, 1);
        }
      LinearAlgebra::distributed::Vector<double> difference = solution1;
      difference -= solution2;
      deallog << "difference norm = " << difference.l2_norm() << std::endl;
    }

#if 0
  DataOut<2> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(position, "X");
  data_out.add_data_vector(solution1, "U1");
  data_out.add_data_vector(solution2, "U2");
  data_out.build_patches(4);
  data_out.write_vtu_with_pvtu_record(
    "./", "solution", 0, comm, 2, 8);
#endif
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll  all;
  const MPI_Comm comm = MPI_COMM_WORLD;

  // Test with p::s::T, mixed FE, multiple components
  {
    parallel::shared::Triangulation<2> tria(
      comm,
      dealii::Triangulation<2>::none,
      true,
      parallel::shared::Triangulation<2>::partition_zorder);
    GridGenerator::cube_and_pyramid(tria);

    hp::FECollection<2>      fe(FESystem<2>(FE_Q<2>(1), 2),
                           FESystem<2>(FE_SimplexP<2>(1), 2));
    hp::MappingCollection<2> mappings(MappingQ<2>(1),
                                      MappingFE<2>(FE_SimplexP<2>(1)));
    DoFHandler<2>            dof_handler(tria);
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            if (cell->reference_cell() == ReferenceCells::Quadrilateral)
              cell->set_active_fe_index(0);
            else
              cell->set_active_fe_index(1);
          }
      }
    dof_handler.distribute_dofs(fe);

    test(dof_handler, mappings);
  }

  // Try discontinuous elements
  {
    parallel::shared::Triangulation<2> tria(
      comm,
      dealii::Triangulation<2>::none,
      true,
      parallel::shared::Triangulation<2>::partition_zorder);
    GridGenerator::hyper_ball(tria);

    hp::FECollection<2>      fe(FESystem<2>(FE_DGQ<2>(4), 2));
    hp::MappingCollection<2> mappings(MappingQ<2>(1));
    DoFHandler<2>            dof_handler(tria);
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        cell->set_active_fe_index(0);
    dof_handler.distribute_dofs(fe);

    test(dof_handler, mappings);
  }

  // Test with p::d::T, hp-FE, multiple components
  {
    parallel::distributed::Triangulation<2> tria(comm);
    GridGenerator::hyper_ball(tria);
    tria.refine_global(3);
    for (const auto &cell : tria.active_cell_iterators())
      if (cell->is_locally_owned() && std::abs(cell->center()[0]) < 0.1)
        cell->set_refine_flag();
    tria.execute_coarsening_and_refinement();

    hp::FECollection<2>      fe(FESystem<2>(FE_Q<2>(2), 2),
                           FESystem<2>(FE_Q<2>(4), 2),
                           FESystem<2>(FE_Nothing<2>(), 2));
    hp::MappingCollection<2> mappings(MappingQ<2>(1),
                                      MappingQ<2>(1),
                                      MappingQ<2>(1));
    DoFHandler<2>            dof_handler(tria);
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        cell->set_active_fe_index(cell->index() % 3);
    dof_handler.distribute_dofs(fe);

    test(dof_handler, mappings);
  }

  // Test with a scalar FE
  {
    parallel::shared::Triangulation<2> tria(
      comm,
      dealii::Triangulation<2>::none,
      true,
      parallel::shared::Triangulation<2>::partition_zorder);
    GridGenerator::cube_and_pyramid(tria);

    hp::FECollection<2>      fe(FE_Q<2>(1), FE_SimplexP<2>(1));
    hp::MappingCollection<2> mappings(MappingQ<2>(1),
                                      MappingFE<2>(FE_SimplexP<2>(1)));
    DoFHandler<2>            dof_handler(tria);
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            if (cell->reference_cell() == ReferenceCells::Quadrilateral)
              cell->set_active_fe_index(0);
            else
              cell->set_active_fe_index(1);
          }
      }
    dof_handler.distribute_dofs(fe);

    test(dof_handler, mappings);
  }

  // Test with another scalar FE
  {
    parallel::shared::Triangulation<2> tria(
      comm,
      dealii::Triangulation<2>::none,
      true,
      parallel::shared::Triangulation<2>::partition_zorder);
    GridGenerator::hyper_cube(tria);
    tria.refine_global(2);

    FE_Q<2>                  fe(5);
    hp::MappingCollection<2> mappings(MappingQ<2>(1));
    DoFHandler<2>            dof_handler(tria);
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        cell->set_active_fe_index(0);
    dof_handler.distribute_dofs(fe);

    test(dof_handler, mappings);
  }

  // Test with another vector FE + grid refinement
  {
    parallel::distributed::Triangulation<2> tria(comm);
    GridGenerator::hyper_cube(tria);
    tria.refine_global(3);

    FESystem<2>   fe(FE_Q<2>(5), 2);
    DoFHandler<2> dof_handler(tria);
    dof_handler.distribute_dofs(fe);
    hp::MappingCollection<2> mappings(MappingQ<2>(1));
    for (const auto &cell : dof_handler.active_cell_iterators())
      cell->set_active_fe_index(0);

    test(dof_handler, mappings);
  }

  // Make sure that we fail correctly when we do weird things
  deal_II_exceptions::disable_abort_on_exception();
  {
    try
      {
        parallel::distributed::Triangulation<2> tria(comm);
        GridGenerator::hyper_cube(tria);
        tria.refine_global(3);

        FESystem<2>   fe(FE_NedelecSZ<2>(0), 1, FE_Q<2>(1), 2);
        DoFHandler<2> dof_handler(tria);
        dof_handler.distribute_dofs(fe);
        hp::MappingCollection<2> mappings(MappingQ<2>(1));
        for (const auto &cell : dof_handler.active_cell_iterators())
          cell->set_active_fe_index(0);

        test(dof_handler, mappings);
      }
    catch (const ExceptionBase &exc)
      {
        deallog << "FE_NedelecSZ nodal renumbering failed successfully."
                << std::endl;
      }
  }

  {
    try
      {
        parallel::distributed::Triangulation<2> tria(comm);
        GridGenerator::hyper_cube(tria);
        tria.refine_global(3);

        FE_NedelecSZ<2> fe(0);
        DoFHandler<2>   dof_handler(tria);
        dof_handler.distribute_dofs(fe);
        hp::MappingCollection<2> mappings(MappingQ<2>(1));
        for (const auto &cell : dof_handler.active_cell_iterators())
          cell->set_active_fe_index(0);

        test(dof_handler, mappings);
      }
    catch (const ExceptionBase &exc)
      {
        deallog << "FE_NedelecSZ nodal renumbering failed successfully."
                << std::endl;
      }
  }
}

// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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



// Test SolutionSerialization for a parallel simulation.

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/solution_serialization.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim>
class SolutionFunction : public Function<dim>
{
public:
  SolutionFunction()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const
  {
    (void)component;

    return p[0];
  }
};



template <int dim>
void
test(const bool use_ss)
{
  const MPI_Comm comm = MPI_COMM_WORLD;

  using VectorType = LinearAlgebra::distributed::Vector<double>;

  MappingQ1<dim> mapping;
  FE_Q<dim>      fe(2);
  QGauss<dim>    quad(3);

  {
    const unsigned int no_refinement = 4;

    parallel::distributed::Triangulation<dim> tria(comm);

    // create coarse grid
    GridGenerator::subdivided_hyper_cube(tria, 2);

    // refine mesh
    for (unsigned int i = 0; i < no_refinement; ++i)
      {
        for (const auto &cell : tria.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              const auto center = cell->center();

              bool flag = true;

              for (unsigned int d = 0; d < dim; ++d)
                flag &= center[d] <= 0.5;

              if (flag)
                cell->set_refine_flag();
            }
        tria.execute_coarsening_and_refinement();
      }

    // create vector to be serialized
    DoFHandler<dim> dof_handler(tria);
    dof_handler.distribute_dofs(fe);

    VectorType vector(std::make_shared<const Utilities::MPI::Partitioner>(
      dof_handler.locally_owned_dofs(),
      DoFTools::extract_locally_active_dofs(dof_handler),
      comm));

    VectorTools::interpolate(mapping,
                             dof_handler,
                             SolutionFunction<dim>(),
                             vector);

    if (use_ss)
      {
        // save vector
        SolutionSerialization<dim, VectorType> ss(dof_handler);
        ss.add_vector(vector);
        ss.save("vector");

        // save mesh
        tria.save("mesh");
      }
    else
      {
        // save vector
        vector.update_ghost_values();

        parallel::distributed::SolutionTransfer<dim, VectorType>
          solution_transfer(dof_handler);
        solution_transfer.prepare_for_coarsening_and_refinement(vector);

        // save mesh
        tria.save("mesh");
      }
  }

  {
    parallel::distributed::Triangulation<dim> tria(comm);

    // create coarse grid
    GridGenerator::subdivided_hyper_cube(tria, 2);

    // load mesh
    tria.load("mesh");

    if (use_ss)
      tria.repartition();

    // create dof-handler and empty vector
    DoFHandler<dim> dof_handler(tria);
    dof_handler.distribute_dofs(fe);

    VectorType vector(std::make_shared<const Utilities::MPI::Partitioner>(
      dof_handler.locally_owned_dofs(),
      DoFTools::extract_locally_active_dofs(dof_handler),
      comm));

    // load vector
    if (use_ss)
      {
        SolutionSerialization<dim, VectorType> ss(dof_handler);
        ss.add_vector(vector);
        ss.load("vector");
      }
    else
      {
        parallel::distributed::SolutionTransfer<dim, VectorType>
          solution_transfer(dof_handler);
        solution_transfer.deserialize(vector);
      }

    // check correctness
    vector.update_ghost_values();

    Vector<float> norm_per_cell(tria.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                      vector,
                                      SolutionFunction<dim>(),
                                      norm_per_cell,
                                      quad,
                                      VectorTools::L2_norm);
    const double error_L2_norm =
      VectorTools::compute_global_error(tria,
                                        norm_per_cell,
                                        VectorTools::L2_norm);

    if (error_L2_norm < 1e-10)
      deallog << "OK!" << std::endl;
  }
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();

  test<2>(true);  // use SolutionSerialization
  test<2>(false); // use SolutionTransfer
}

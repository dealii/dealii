// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// PETSc up to at least 3.5 has a bug in which creating a ghosted vector
// properly zeros out the locally owned elements but leaves the ghost elements
// uninitialized.
//
// this test, originally by Michal Wichrowski, first does all sorts of complex
// things to prime the memory, and then creates a PETSc MPI vector to check
// whether the ghost elements are zero
//
// the workaround for the bug in PETSc is conditional on PETSC_VERSION <= 3.5,
// which at the time of writing the test was the current version. if the test
// fails with newer PETSc versions, then this means that the PETSc folks
// didn't apply the fix in time
//
// reference: https://code.google.com/p/dealii/issues/detail?id=233


#include "../tests.h"

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/point.h>
#include <deal.II/base/synchronous_iterator.h>
#include <deal.II/base/function.h>

#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/solver_control.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <iostream>
#include <fstream>


void test ()
{
  const unsigned int this_mpi_process=Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int n_mpi_processes=Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  MPI_Comm mpi_communicator= MPI_COMM_WORLD;
  parallel::distributed::Triangulation<2> fluid_triangulation(mpi_communicator,
                                                              typename Triangulation<2>::MeshSmoothing
                                                              (Triangulation<2>::smoothing_on_refinement |
                                                                  Triangulation<2>::smoothing_on_coarsening));

  {
    Triangulation<2> tmp_tria_fluid;
    Triangulation<2> tmp_tria_buffer;
    Triangulation<2> tmp_rectangle;
    Triangulation<2> tmp_tria_solid;
    //-------Generate and colorize triangulation:-------
    //Triangulation parameters
    const unsigned int h_cells=11;
    const unsigned int outflow_cells=25;
    const unsigned int inflow_cells=5;
    const unsigned int solid_cells_x=20;
    const unsigned int solid_cells_y=1;
    const unsigned int v_space_cells=5;;
    const double cell_size_x=0.1;
    const double cell_size_y=0.1;
    const double lenght = solid_cells_x*cell_size_x,
                 height= static_cast<double>(h_cells*cell_size_y),
                 h_distance=static_cast<double>(cell_size_x*inflow_cells),
                 v_distance=static_cast<double>(cell_size_y*v_space_cells),
                 solid_top=static_cast<double>(v_distance+solid_cells_y*cell_size_y),
                 total_length=(inflow_cells+solid_cells_x+outflow_cells)*cell_size_x;


    std::vector< std::vector< double > > step_sizes(2);
    //----------------INFLOW ---------------
    for (unsigned int i =0; i<inflow_cells; ++i)
      step_sizes[0].push_back(cell_size_x);
    for (unsigned int i =0; i<h_cells; ++i)
      step_sizes[1].push_back(cell_size_y);

    GridGenerator::subdivided_hyper_rectangle(tmp_tria_fluid,
                                              step_sizes,
                                              Point<2>(0,0),
                                              Point<2>(h_distance,height),
                                              false);

    //------------------------LOWER-------------
    step_sizes[0].clear();
    step_sizes[1].clear();
    for (unsigned int i =0; i<solid_cells_x; ++i)
      step_sizes[0].push_back(cell_size_x);
    for (unsigned int i =0; i<v_space_cells; ++i)
      step_sizes[1].push_back(cell_size_y);

    GridGenerator::subdivided_hyper_rectangle(tmp_rectangle,
                                              step_sizes,
                                              Point<2>(h_distance,0),
                                              Point<2>(lenght+h_distance,v_distance),
                                              false);

    tmp_tria_buffer.copy_triangulation(tmp_tria_fluid);
    GridGenerator::merge_triangulations(tmp_tria_buffer,
                                        tmp_rectangle,
                                        tmp_tria_fluid);

    tmp_tria_buffer.clear();
    tmp_rectangle.clear();

    //---------------------------UPPER--------


    step_sizes[1].clear();
    for (unsigned int i =0; i<(h_cells-v_space_cells-solid_cells_y); ++i)
      step_sizes[1].push_back(cell_size_y);

    GridGenerator::subdivided_hyper_rectangle(tmp_rectangle,
                                              step_sizes,
                                              Point<2>(h_distance,solid_top),
                                              Point<2>(h_distance+lenght,height),
                                              false);

    tmp_tria_buffer.copy_triangulation(tmp_tria_fluid);
    GridGenerator::merge_triangulations(tmp_tria_buffer,
                                        tmp_rectangle,
                                        tmp_tria_fluid);
    tmp_tria_buffer.clear();
    tmp_rectangle.clear();
    //
    //----------------------------outflow----------
    step_sizes[0].clear();
    step_sizes[1].clear();
    for (unsigned int i =0; i<outflow_cells; ++i)
      step_sizes[0].push_back(cell_size_x);
    for (unsigned int i =0; i<h_cells; ++i)
      step_sizes[1].push_back(cell_size_y);
    GridGenerator::subdivided_hyper_rectangle(tmp_rectangle,
                                              step_sizes,
                                              Point<2>(h_distance+lenght,0),
                                              Point<2>(total_length,height),
                                              false);

    tmp_tria_buffer.copy_triangulation(tmp_tria_fluid);
    GridGenerator::merge_triangulations(tmp_tria_buffer,
                                        tmp_rectangle,
                                        tmp_tria_fluid);

    //-----------COPY to distributed
    fluid_triangulation.copy_triangulation(tmp_tria_fluid);
  }


  fluid_triangulation.refine_global(1);
  //     ---------END GENERATING TRIA ------------

  // call data_out.build_patches() once but then destroy the object.
  {
    FE_Q<2> fe(1);
    DoFHandler<2> dof_handler(fluid_triangulation);

    DataOut<2> data_out;
    dof_handler.distribute_dofs(fe);
    data_out.attach_dof_handler (dof_handler);

    Vector<float> subdomain (fluid_triangulation.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = fluid_triangulation.locally_owned_subdomain();
    data_out.add_data_vector (subdomain, "subdomain"  );

    data_out.build_patches ();
  }

  // do it again. the bug is that we get wrong data because of the call above
  {
    FE_Q<2> fe(2);
    DoFHandler<2> handler(fluid_triangulation);
    handler.distribute_dofs(fe);

    IndexSet locally_owned_dofs = handler.locally_owned_dofs();
    IndexSet locally_relevant_dofs = handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(handler,
                                            locally_relevant_dofs);


    PETScWrappers::MPI::Vector vector;
    vector.reinit(locally_owned_dofs,
                  locally_relevant_dofs,
                  mpi_communicator);

    // ensure that all elements we can access locally are initially at
    // zero. check this first for the locally owned elements...
    for (unsigned int i=0; i<handler.n_dofs(); ++i)
      if (locally_owned_dofs.is_element(i))
        AssertThrow (get_real_assert_zero_imag(vector(i)) == 0,
                     ExcInternalError());
    // ...end then also for the ghost elements
    for (unsigned int i=0; i<handler.n_dofs(); ++i)
      if (locally_relevant_dofs.is_element(i)
          &&
          !locally_owned_dofs.is_element(i))
        AssertThrow (get_real_assert_zero_imag(vector(i)) == 0,
                     ExcInternalError());
  }
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.threshold_double(1.e-10);

      test();

      // if we got here, everything was fine:
      deallog << "OK" << std::endl;
    }
  else
    test();
}

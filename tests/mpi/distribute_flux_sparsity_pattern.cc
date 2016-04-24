// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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

/*
 * Verify that distributed flux sparsity patterns are built correctly. There
 * used to be a bug where coupling between an active, locally owned cell and
 * an active ghost cell may be erroneously skipped when building a sparsity
 * pattern with flux terms.
 */


#include "../tests.h"

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <fstream>
#include <iostream>

namespace LinearAdvectionTest
{
  using namespace dealii;

  template<int dim>
  class AdvectionProblem
  {
  public:
    AdvectionProblem();
    void run();

  private:
    void setup_system();
    void assemble_system();
    void calculate_flux_terms
    (const TriaActiveIterator<DoFCellAccessor<DoFHandler<dim>, false> > &current_cell,
     FEFaceValues<dim> &current_face_values,
     const TriaIterator<DoFCellAccessor<DoFHandler<dim>, false> > &neighbor_cell,
     FEFaceValues<dim> &neighbor_face_values,
     FullMatrix<double> &current_to_current_flux,
     FullMatrix<double> &current_to_neighbor_flux,
     FullMatrix<double> &neighbor_to_current_flux,
     FullMatrix<double> &neighbor_to_neighbor_flux);

    parallel::distributed::Triangulation<dim> triangulation;
    FE_DGQ<dim>        fe;
    DoFHandler<dim>    dof_handler;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    TrilinosWrappers::SparseMatrix system_matrix;
  };



  template<int dim>
  AdvectionProblem<dim>::AdvectionProblem()
    : triangulation(MPI_COMM_WORLD),
      fe(1),
      dof_handler(triangulation)
  {
    std::vector<unsigned int> repetitions(2);
    repetitions[0] = 2;
    repetitions[1] = 1;

    const Point<2> p0(0.0, 0.0);
    const Point<2> p1(2.0, 1.0);
    GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions,
                                              p0, p1);
  }

  template<int dim>
  void AdvectionProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(fe);
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    DynamicSparsityPattern dynamic_sparsity_pattern(locally_relevant_dofs);
    Table<2, DoFTools::Coupling> cell_integral_mask(fe.n_components(), fe.n_components());
    Table<2, DoFTools::Coupling> flux_integral_mask(fe.n_components(), fe.n_components());
    for (unsigned int i = 0; i < fe.n_components(); ++i)
      {
        for (unsigned int j = 0; j < fe.n_components(); ++j)
          {
            cell_integral_mask(i, j) = DoFTools::Coupling::always;
            flux_integral_mask(i, j) = DoFTools::Coupling::nonzero;
          }
      }

    DoFTools::make_flux_sparsity_pattern(dof_handler, dynamic_sparsity_pattern,
                                         cell_integral_mask, flux_integral_mask);
    SparsityTools::distribute_sparsity_pattern
    (dynamic_sparsity_pattern, dof_handler.n_locally_owned_dofs_per_processor(),
     MPI_COMM_WORLD, locally_relevant_dofs);

    system_matrix.reinit(locally_owned_dofs,
                         locally_owned_dofs,
                         dynamic_sparsity_pattern,
                         MPI_COMM_WORLD);
  }



  template<int dim>
  void AdvectionProblem<dim>::calculate_flux_terms
  (const TriaActiveIterator<DoFCellAccessor<DoFHandler<dim>, false> > &current_cell,
   FEFaceValues<dim> &current_face_values,
   const TriaIterator<DoFCellAccessor<DoFHandler<dim>, false> > &neighbor_cell,
   FEFaceValues<dim> &neighbor_face_values,
   FullMatrix<double> &current_to_current_flux,
   FullMatrix<double> &current_to_neighbor_flux,
   FullMatrix<double> &neighbor_to_current_flux,
   FullMatrix<double> &neighbor_to_neighbor_flux)
  {
    for (unsigned int q_point_n = 0;
         q_point_n < neighbor_face_values.n_quadrature_points;
         ++q_point_n)
      {
        /*
         * Note that if we are integrating the neighbor's test functions, then
         * we are really performing the integration on the neighbor's face.
         *
         * Do integration on both cells (rather than use something like
         * upwinding) since this is just a test.
         */
        for (unsigned int i = 0; i < neighbor_face_values.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < neighbor_face_values.dofs_per_cell; ++j)
              {
                current_to_current_flux(i, j) +=
                  current_face_values.shape_value(i, q_point_n)
                  * current_face_values.shape_value(j, q_point_n)
                  * current_face_values.JxW(q_point_n);

                current_to_neighbor_flux(i, j) -=
                  neighbor_face_values.shape_value(i, q_point_n)
                  * current_face_values.shape_value(j, q_point_n)
                  * current_face_values.JxW(q_point_n);

                neighbor_to_current_flux(i, j) +=
                  current_face_values.shape_value(i, q_point_n)
                  * neighbor_face_values.shape_value(j, q_point_n)
                  * current_face_values.JxW(q_point_n);

                neighbor_to_neighbor_flux(i, j) -=
                  neighbor_face_values.shape_value(i, q_point_n)
                  * neighbor_face_values.shape_value(j, q_point_n)
                  * neighbor_face_values.JxW(q_point_n);
              }
          }
      }
    std::vector<types::global_dof_index> current_dofs
    (current_face_values.dofs_per_cell);
    current_cell->get_dof_indices(current_dofs);

    std::vector<types::global_dof_index> neighbor_dofs
    (neighbor_face_values.dofs_per_cell);
    neighbor_cell->get_dof_indices(neighbor_dofs);

    system_matrix.add(current_dofs, current_to_current_flux);
    system_matrix.add(neighbor_dofs, current_dofs, current_to_neighbor_flux);
    system_matrix.add(neighbor_dofs, neighbor_to_neighbor_flux);
    system_matrix.add(current_dofs, neighbor_dofs, neighbor_to_current_flux);
  }



  template<int dim>
  void AdvectionProblem<dim>::assemble_system()
  {
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    FullMatrix<double> current_to_current_flux(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> current_to_neighbor_flux(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> neighbor_to_current_flux(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> neighbor_to_neighbor_flux(dofs_per_cell, dofs_per_cell);

    const QGauss<dim - 1> face_quadrature(3);

    const UpdateFlags update_flags =
      update_values | update_quadrature_points | update_JxW_values;

    FEFaceValues<dim> current_face_values(fe, face_quadrature, update_flags
                                          | update_normal_vectors);

    FEFaceValues<dim> neighbor_face_values(fe, face_quadrature, update_flags);
    FESubfaceValues<dim> neighbor_subface_values(fe, face_quadrature, update_flags);

    typename DoFHandler<dim>::active_cell_iterator
    current_cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; current_cell != endc; ++current_cell)
      {
        if (current_cell->is_locally_owned())
          {
            for (unsigned int face_n = 0; face_n < GeometryInfo<dim>::faces_per_cell;
                 ++face_n)
              {
                const int neighbor_index = current_cell->neighbor_index(face_n);
                if (neighbor_index != -1) // interior face
                  {
                    typename DoFHandler<dim>::active_cell_iterator neighbor_cell =
                      neighbor_cell = current_cell->neighbor(face_n);

                    bool do_face_integration = false;
                    bool neighbor_is_level_lower = false;

                    /*
                     * Always integrate if the current cell is more refined
                     * than the neighbor.
                     */
                    if (current_cell->level() > neighbor_cell->level())
                      {
                        do_face_integration = true;
                        neighbor_is_level_lower = true;
                      }
                    // If the neighbor is not active, then it is at a higher
                    // refinement level (so we do not need to integrate now)
                    if (neighbor_cell->active())
                      {
                        if (neighbor_cell->is_locally_owned())
                          {
                            if (neighbor_cell < current_cell)
                              {
                                do_face_integration = true;
                              }
                          }
                        else
                          {
                            Assert(neighbor_cell->is_ghost(),
                                   ExcMessage("All neighbors should be locally "
                                              "owned or ghost cells."));
                            if (current_cell->level() == neighbor_cell->level()
                                && current_cell->subdomain_id() < neighbor_cell->subdomain_id())
                              {
                                do_face_integration = true;
                              }
                          }
                      }

                    if (do_face_integration)
                      {
                        const unsigned int neighbor_face_n =
                          current_cell->neighbor_face_no(face_n);
                        AssertThrow(!neighbor_is_level_lower, ExcInternalError());
                        neighbor_face_values.reinit(neighbor_cell, neighbor_face_n);

                        current_face_values.reinit(current_cell, face_n);
                        current_to_current_flux = 0.0;
                        current_to_neighbor_flux = 0.0;
                        neighbor_to_current_flux = 0.0;
                        neighbor_to_neighbor_flux = 0.0;

                        calculate_flux_terms
                        (current_cell, current_face_values, neighbor_cell, neighbor_face_values,
                         current_to_current_flux, current_to_neighbor_flux,
                         neighbor_to_current_flux, neighbor_to_neighbor_flux);
                      }
                  }
                else // boundary face: skip for the purposes of this test
                  {
                  }
              }
          }
      }

    system_matrix.compress(VectorOperation::add);
  }



  template<int dim>
  void AdvectionProblem<dim>::run()
  {
    setup_system();
    assemble_system();
  }
}

int main(int argc, char **argv)
{
  using namespace dealii;
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  Utilities::MPI::MPI_InitFinalize mpi_initialization
  (argc, argv, 1);
  LinearAdvectionTest::AdvectionProblem<2> advection_problem;
  advection_problem.run();
  deallog << "OK" << std::endl;
}

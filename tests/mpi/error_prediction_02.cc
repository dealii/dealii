// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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



// Test to check if ErrorPredictor works in parallel with hp::DoFHandler.
// This tests is based on hp/error_prediction.cc


#include <deal.II/distributed/cell_data_transfer.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/hp/refinement.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"

#include "../test_grids.h"


template <int dim>
void
test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  // ------ setup ------
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  TestGrids::hyper_line(tria, 4);

  for (auto cell = tria.begin(0); cell != tria.end(0); ++cell)
    if (cell->id().to_string() == "0_0:" || cell->id().to_string() == "1_0:")
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  hp::FECollection<dim> fes;
  for (unsigned int d = 1; d <= 3; ++d)
    fes.push_back(FE_Q<dim>(d));

  DoFHandler<dim> dh(tria);
  for (const auto &cell : dh.active_cell_iterators())
    {
      // set active FE index
      if (cell->is_locally_owned())
        cell->set_active_fe_index(1);
    }
  for (auto cell = dh.begin(0); cell != dh.end(0); ++cell)
    {
      // set refinement/coarsening flags
      if (cell->id().to_string() == "0_0:")
        {
          // h-coarsening and p-refinement
          for (unsigned int i = 0; i < cell->n_children(); ++i)
            if (cell->child(i)->is_locally_owned())
              {
                cell->child(i)->set_coarsen_flag();
                cell->child(i)->set_future_fe_index(2);
              }
        }
      else if (cell->id().to_string() == "1_0:")
        {
          // h-coarsening and p-coarsening
          for (unsigned int i = 0; i < cell->n_children(); ++i)
            if (cell->child(i)->is_locally_owned())
              {
                cell->child(i)->set_coarsen_flag();
                cell->child(i)->set_future_fe_index(0);
              }
        }
      else if (cell->id().to_string() == "2_0:")
        {
          // h-refinement and p-refinement
          if (cell->is_locally_owned())
            {
              cell->set_refine_flag();
              cell->set_future_fe_index(2);
            }
        }
      else if (cell->id().to_string() == "3_0:")
        {
          // h-refinement and p-coarsening
          if (cell->is_locally_owned())
            {
              cell->set_refine_flag();
              cell->set_future_fe_index(0);
            }
        }
    }
  dh.distribute_dofs(fes);

  // ----- prepare error indicators -----
  Vector<float> error_indicators(tria.n_active_cells());
  for (unsigned int i = 0; i < error_indicators.size(); ++i)
    error_indicators(i) = 10.;

  // ----- connect error predictor -----
  Vector<float> predicted_errors;
  tria.signals.post_p4est_refinement.connect([&]() {
    const parallel::distributed::TemporarilyMatchRefineFlags<dim>
      refine_modifier(tria);
    predicted_errors.reinit(tria.n_active_cells());
    hp::Refinement::predict_error(dh,
                                  error_indicators,
                                  predicted_errors,
                                  /*gamma_p=*/0.5,
                                  /*gamma_h=*/1.,
                                  /*gamma_n=*/1.);
  });

  // ----- verify ------
  deallog << "pre_adaptation" << std::endl;
  for (const auto &cell : dh.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        deallog << " cell:" << cell->id().to_string()
                << " fe_deg:" << cell->get_fe().degree
                << " error:" << error_indicators[cell->active_cell_index()];

        if (cell->coarsen_flag_set())
          deallog << " coarsening";
        else if (cell->refine_flag_set())
          deallog << " refining";

        if (cell->future_fe_index_set())
          deallog << " future_fe_deg:" << fes[cell->future_fe_index()].degree;

        deallog << std::endl;
      }

  // ----- execute adaptation -----
  parallel::distributed::CellDataTransfer<dim, dim, Vector<float>>
  data_transfer(tria,
                /*transfer_variable_size_data=*/false,
                &AdaptationStrategies::Refinement::l2_norm<dim, dim, float>,
                &AdaptationStrategies::Coarsening::l2_norm<dim, dim, float>);

  data_transfer.prepare_for_coarsening_and_refinement(predicted_errors);
  tria.execute_coarsening_and_refinement();

  predicted_errors.reinit(tria.n_active_cells());
  data_transfer.unpack(predicted_errors);

  // ------ verify ------
  deallog << "post_adaptation" << std::endl;
  for (const auto &cell : dh.active_cell_iterators())
    if (cell->is_locally_owned())
      deallog << " cell:" << cell->id().to_string()
              << " predicted:" << predicted_errors(cell->active_cell_index())
              << std::endl;

  // make sure no processor is hanging
  MPI_Barrier(MPI_COMM_WORLD);

  deallog << "OK" << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();
}

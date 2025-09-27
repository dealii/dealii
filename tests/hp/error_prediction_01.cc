// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// validate combination of error prediction and cell data transfer algorithms
// for hp-adaptive methods
// test either h- or p-adaptation


#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/tria.h>

#include <deal.II/hp/refinement.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/cell_data_transfer.h>

#include "../tests.h"

#include "../test_grids.h"


template <int dim>
void
test()
{
  // ----- setup -----
  Triangulation<dim> tria;
  TestGrids::hyper_line(tria, 4);

  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  hp::FECollection<dim> fes;
  for (unsigned int d = 1; d <= 3; ++d)
    fes.push_back(FE_Q<dim>(d));

  DoFHandler<dim> dh(tria);
  for (const auto &cell : dh.active_cell_iterators())
    {
      // set active FE index
      cell->set_active_fe_index(1);
    }
  for (auto cell = dh.begin(0); cell != dh.end(0); ++cell)
    {
      if (cell->id().to_string() == "0_0:")
        {
          // h-coarsening
          for (unsigned int i = 0; i < cell->n_children(); ++i)
            cell->child(i)->set_coarsen_flag();
        }
      else if (cell->id().to_string() == "1_0:")
        {
          // h-refinement
          cell->set_refine_flag();
        }
      else if (cell->id().to_string() == "2_0:")
        {
          // p-refinement
          cell->set_future_fe_index(2);
        }
      else if (cell->id().to_string() == "3_0:")
        {
          // p-coarsening
          cell->set_future_fe_index(0);
        }
    }
  dh.distribute_dofs(fes);

  // ----- predict -----
  Vector<float> error_indicators(tria.n_active_cells());
  for (unsigned int i = 0; i < tria.n_active_cells(); ++i)
    error_indicators[i] = 10.;

  Vector<float> predicted_error_indicators(tria.n_active_cells());
  hp::Refinement::predict_error(dh,
                                error_indicators,
                                predicted_error_indicators,
                                /*gamma_p=*/0.5,
                                /*gamma_h=*/1.,
                                /*gamma_n=*/1.);

  // ----- verify ------
  deallog << "pre_adaptation" << std::endl
          << " ncells:" << tria.n_active_cells() << std::endl;
  for (const auto &cell : dh.active_cell_iterators())
    {
      deallog << " cell:" << cell->id().to_string()
              << " fe_deg:" << cell->get_fe().degree
              << " error:" << error_indicators[cell->active_cell_index()]
              << " predicted:"
              << predicted_error_indicators[cell->active_cell_index()];

      if (cell->refine_flag_set())
        deallog << " refining";
      else if (cell->coarsen_flag_set())
        deallog << " coarsening";

      if (cell->future_fe_index_set())
        deallog << " future_fe_deg:" << fes[cell->future_fe_index()].degree;

      deallog << std::endl;
    }

  // ----- execute adaptation -----
  CellDataTransfer<dim, dim, Vector<float>> cell_data_transfer(
    tria,
    &AdaptationStrategies::Refinement::l2_norm<dim, dim, float>,
    &AdaptationStrategies::Coarsening::l2_norm<dim, dim, float>);
  cell_data_transfer.prepare_for_coarsening_and_refinement();

  tria.execute_coarsening_and_refinement();

  Vector<float> transferred_indicators(tria.n_active_cells());
  cell_data_transfer.unpack(predicted_error_indicators, transferred_indicators);

  // ----- verify -----
  deallog << "post_adaptation" << std::endl
          << " ncells:" << tria.n_active_cells() << std::endl;
  for (const auto &cell : dh.active_cell_iterators())
    deallog << " cell:" << cell->id().to_string()
            << " fe_deg:" << cell->get_fe().degree << " transferred:"
            << transferred_indicators[cell->active_cell_index()] << std::endl;

  // ----- verify norms -----
  const double predicted_error_pre  = predicted_error_indicators.l2_norm(),
               predicted_error_post = transferred_indicators.l2_norm();

  deallog << "predicted_error_norms" << std::endl
          << " pre_adaptation:" << predicted_error_pre << std::endl
          << " post_adaptation:" << predicted_error_post << std::endl;

  Assert(predicted_error_pre == predicted_error_post,
         ExcMessage("Transfer failed - Results not similar."));

  deallog << "OK" << std::endl;
}



int
main(int argc, char *argv[])
{
  initlog();

  deallog.push("1d");
  test<1>();
  deallog.pop();
  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();
}

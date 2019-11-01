// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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



// validate error prediction algorithm for hp adaptive methods


#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/refinement.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"



template <int dim>
void
test()
{
  const unsigned int n_cells = 4;

  // ----- setup -----
  Triangulation<dim>        tria;
  std::vector<unsigned int> rep(dim, 1);
  rep[0] = n_cells;
  Point<dim> p1, p2;
  for (unsigned int d = 0; d < dim; ++d)
    {
      p1[d] = 0;
      p2[d] = (d == 0) ? n_cells : 1;
    }
  GridGenerator::subdivided_hyper_rectangle(tria, rep, p1, p2);

  hp::FECollection<dim> fes;
  for (unsigned int d = 1; d <= 3; ++d)
    fes.push_back(FE_Q<dim>(d));

  hp::DoFHandler<dim> dh(tria);
  dh.set_fe(fes);
  for (const auto &cell : dh.active_cell_iterators())
    {
      if (cell->id().to_string() == "0_0:")
        cell->set_refine_flag();
      else if (cell->id().to_string() == "1_0:")
        cell->set_coarsen_flag();
      else if (cell->id().to_string() == "2_0:")
        cell->set_future_fe_index(2);
    }

  // ----- predict -----
  Vector<float> error_indicators, predicted_errors;
  error_indicators.reinit(tria.n_active_cells());
  predicted_errors.reinit(tria.n_active_cells());
  for (unsigned int i = 0; i < tria.n_active_cells(); ++i)
    error_indicators[i] = 10;

  hp::Refinement::predict_error(dh,
                                error_indicators,
                                predicted_errors,
                                /*gamma_p=*/0.5,
                                /*gamma_h=*/1.,
                                /*gamma_n=*/1.);

  // ----- verify ------
  deallog << "ncells:" << tria.n_active_cells() << std::endl;
  for (const auto &cell : dh.active_cell_iterators())
    {
      deallog << " cell:" << cell->id().to_string()
              << " fe_deg:" << cell->get_fe().degree
              << " error:" << error_indicators[cell->active_cell_index()]
              << " predict:" << predicted_errors[cell->active_cell_index()];

      if (cell->refine_flag_set())
        deallog << " refining";
      else if (cell->coarsen_flag_set())
        deallog << " coarsening";
      else if (cell->future_fe_index_set())
        deallog << " future_fe_deg:"
                << dh.get_fe_collection()[cell->future_fe_index()].degree;

      deallog << std::endl;
    }

  // ----- check feature -----
  hp::Refinement::p_adaptivity_from_reference(
    dh,
    error_indicators,
    predicted_errors,
    /*compare_refine=*/std::less<float>(),
    /*compare_coarsen=*/std::less<float>());

  deallog << "OK" << std::endl;
}



int
main(int argc, char *argv[])
{
  initlog();
  deallog << std::setprecision(3);

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

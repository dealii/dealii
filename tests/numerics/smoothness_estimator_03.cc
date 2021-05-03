// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2021 by the deal.II authors
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



// This test is used to make sure that FESeries::Fourier/Legendre
// work with non-primitive elements


#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_series.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/smoothness_estimator.h>

#include "../tests.h"

template <int dim>
void
test()
{
  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FE_Nedelec<dim>(0));
  fe_collection.push_back(FE_Nedelec<dim>(1));

  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube(tria, 2);

  DoFHandler<dim> dh(tria);
  dh.distribute_dofs(fe_collection);

  Vector<double> solution(dh.n_dofs());
  Vector<float>  smoothness(tria.n_active_cells());

  for (unsigned int i = 0; i < dim; ++i)
    {
      smoothness = 0.0;
      FESeries::Fourier<dim> fourier =
        SmoothnessEstimator::Fourier::default_fe_series(fe_collection, i);
      SmoothnessEstimator::Fourier::coefficient_decay(fourier,
                                                      dh,
                                                      solution,
                                                      smoothness);
    }


  for (unsigned int i = 0; i < dim; ++i)
    {
      smoothness = 0.0;
      FESeries::Legendre<dim> legendre =
        SmoothnessEstimator::Legendre::default_fe_series(fe_collection, i);
      SmoothnessEstimator::Legendre::coefficient_decay(legendre,
                                                       dh,
                                                       solution,
                                                       smoothness);
    }

  deallog << "Ok" << std::endl;
}

int
main()
{
  initlog();

  test<2>();
  test<3>();

  return 0;
}

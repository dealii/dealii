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



// check serialization for FESeries::Legendre

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_series.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>

#include "serialization.h"

template <int dim>
void
test()
{
  // setup
  std::vector<unsigned int> n_modes;
  hp::FECollection<dim>     hp_fe;
  hp::QCollection<dim>      hp_q;

  const unsigned int min_degree = 1, max_degree = 2;
  const QGauss<dim>  quadrature(max_degree + 1);
  const QSorted<dim> quadrature_sorted(quadrature);
  for (unsigned int p = min_degree; p <= max_degree; ++p)
    {
      n_modes.push_back(max_degree + 1);
      hp_fe.push_back(FE_Q<dim>(p));
      hp_q.push_back(quadrature_sorted);
    }

  FESeries::Legendre<dim> legendre_save(n_modes, hp_fe, hp_q);
  FESeries::Legendre<dim> legendre_load(n_modes, hp_fe, hp_q);

  // create transformation matrices
  legendre_save.precalculate_all_transformation_matrices();

  // save series expansion object
  std::ostringstream            oss;
  boost::archive::text_oarchive oa(oss, boost::archive::no_header);
  legendre_save.save_transformation_matrices(oa, 0);
  deallog << oss.str() << std::endl;

  // load series expansion object
  std::istringstream            iss(oss.str());
  boost::archive::text_iarchive ia(iss, boost::archive::no_header);
  legendre_load.load_transformation_matrices(ia, 0);
  AssertThrow(compare(legendre_save, legendre_load), ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main()
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

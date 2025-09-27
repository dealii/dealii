// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check serialization for FESeries::Fourier

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

  const unsigned int   min_degree = 1, max_degree = 2;
  const QGauss<1>      base_quadrature(4);
  const QIterated<dim> quadrature(base_quadrature, max_degree);
  const QSorted<dim>   quadrature_sorted(quadrature);
  for (unsigned int p = min_degree; p <= max_degree; ++p)
    {
      n_modes.push_back(max_degree + 1);
      hp_fe.push_back(FE_Q<dim>(p));
      hp_q.push_back(quadrature_sorted);
    }

  FESeries::Fourier<dim> fourier_save(n_modes, hp_fe, hp_q);
  FESeries::Fourier<dim> fourier_load(n_modes, hp_fe, hp_q);

  // create transformation matrices
  fourier_save.precalculate_all_transformation_matrices();

  // save series expansion object
  std::ostringstream            oss;
  boost::archive::text_oarchive oa(oss, boost::archive::no_header);
  fourier_save.save_transformation_matrices(oa, 0);
  deallog << oss.str() << std::endl;

  // load expansion object
  std::istringstream            iss(oss.str());
  boost::archive::text_iarchive ia(iss, boost::archive::no_header);
  fourier_load.load_transformation_matrices(ia, 0);
  AssertThrow(compare(fourier_save, fourier_load), ExcInternalError());

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

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



// precalculate hp::FEValues objects


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include "../tests.h"



template <int dim>
void
test()
{
  const unsigned int max_degree = 3;

  hp::FECollection<dim> fes;
  for (unsigned int degree = 1; degree <= max_degree; ++degree)
    fes.push_back(FE_Q<dim>(degree));

  // provide quadrature and mapping collections corresponding to the size of the
  // finite element collection
  {
    hp::MappingCollection<dim> mappings;
    hp::QCollection<dim>       quadratures;

    for (unsigned int i = 0; i < fes.size(); ++i)
      {
        mappings.push_back(MappingQ<dim>(fes[i].degree));
        quadratures.push_back(QGauss<dim>(fes[i].degree + 1));
      }

    hp::FEValues<dim> hp_fe_values(mappings, fes, quadratures, update_values);
    hp_fe_values.precalculate_fe_values();
  }

  // provide just one quadrature rule and the default mapping object
  {
    hp::QCollection<dim> quadratures(QGauss<dim>(max_degree + 1));

    hp::FEValues<dim> hp_fe_values(fes, quadratures, update_values);
    hp_fe_values.precalculate_fe_values();
  }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  deallog.push("2d");
  test<2>();
  deallog.pop();
}

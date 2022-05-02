// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

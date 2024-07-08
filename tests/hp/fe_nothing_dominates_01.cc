// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check for continuity requirements for neighboring FE_Nothing and FE_Q
// elements.
//
// Note: output corresponds to current results, which are not correct.


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>

#include "../tests.h"

#include "fe_nothing_dominates.h"


template <int dim>
void
test(const bool fe_nothing_dominates)
{
  Functions::ConstantFunction<dim> function(1.);

  {
    deallog << "(FE_Nothing, FE_Q)" << std::endl;
    hp::FECollection<dim> fe_collection(FE_Nothing<dim>(/*n_components=*/1,
                                                        fe_nothing_dominates),
                                        FE_Q<dim>(1));
    hp::QCollection<dim>  q_collection(
      Quadrature<dim>(std::vector<Point<dim>>{Point<dim>()},
                      std::vector<double>{0.}),
      QGauss<dim>(2));
    project(fe_collection, q_collection, function);
  }

  {
    deallog << "(FE_Q, FE_Nothing)" << std::endl;
    hp::FECollection<dim> fe_collection(FE_Q<dim>(1),
                                        FE_Nothing<dim>(/*n_components=*/1,
                                                        fe_nothing_dominates));
    hp::QCollection<dim>  q_collection(
      QGauss<dim>(2),
      Quadrature<dim>(std::vector<Point<dim>>{Point<dim>()},
                      std::vector<double>{0.}));
    project(fe_collection, q_collection, function);
  }
}


int
main()
{
  initlog();

  deallog << std::boolalpha;
  for (const bool fe_nothing_dominates : {false, true})
    {
      deallog << "FE_Nothing dominates: " << fe_nothing_dominates << std::endl;
      deallog.push("1d");
      test<1>(fe_nothing_dominates);
      deallog.pop();
      deallog.push("2d");
      test<2>(fe_nothing_dominates);
      deallog.pop();
      deallog.push("3d");
      test<3>(fe_nothing_dominates);
      deallog.pop();
    }
}

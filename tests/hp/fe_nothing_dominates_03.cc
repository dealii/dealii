// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check for continuity requirements for neighboring
// (FE_QxFE_Nothing) and (FE_NothingxFE_Q) elements.
// The twist: only one FE_Nothing element dominates.
//
// Note: output corresponds to current results, which are not correct.


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>

#include "../tests.h"

#include "fe_nothing_dominates.h"


template <int dim>
void
test()
{
  Functions::ConstantFunction<dim> function(1., 2);
  hp::QCollection<dim>             q_collection(QGauss<dim>(2), QGauss<dim>(2));

  {
    deallog << "(FE_Nothing(true)xFE_Q, FE_QxFE_Nothing(false))" << std::endl;
    hp::FECollection<dim> fe_collection(
      FESystem<dim>(FE_Nothing<dim>(1, true), FE_Q<dim>(1)),
      FESystem<dim>(FE_Q<dim>(1), FE_Nothing<dim>(1, false)));
    project(fe_collection, q_collection, function);
  }
  {
    deallog << "(FE_QxFE_Nothing(false), FE_Nothing(true)xFE_Q)" << std::endl;
    hp::FECollection<dim> fe_collection(
      FESystem<dim>(FE_Q<dim>(1), FE_Nothing<dim>(1, false)),
      FESystem<dim>(FE_Nothing<dim>(1, true), FE_Q<dim>(1)));
    project(fe_collection, q_collection, function);
  }
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

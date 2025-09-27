// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// call VectorTools::integrate_difference with fe's distributed in the
// same random way as in hp/random


#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
void
test()
{
  deallog << "dim=" << dim << std::endl;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.refine_global(4 - dim);

  hp::FECollection<dim> fe_collection;
  hp::QCollection<dim>  q_collection;
  for (unsigned int i = 1; i <= 4; ++i)
    {
      fe_collection.push_back(FE_Q<dim>(QIterated<1>(QTrapezoid<1>(), i)));
      q_collection.push_back(QGauss<dim>(i + 2));
    }


  DoFHandler<dim> dof_handler(tria);

  for (typename DoFHandler<dim>::active_cell_iterator cell =
         dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    cell->set_active_fe_index(Testing::rand() % fe_collection.size());

  dof_handler.distribute_dofs(fe_collection);

  Vector<double> vec(dof_handler.n_dofs());
  for (unsigned int i = 0; i < vec.size(); ++i)
    vec(i) = i;

  Vector<float> diff(tria.n_active_cells());

  VectorTools::NormType norms[] = {VectorTools::mean,
                                   VectorTools::L1_norm,
                                   VectorTools::L2_norm,
                                   VectorTools::Linfty_norm,
                                   VectorTools::H1_seminorm,
                                   VectorTools::W1p_seminorm};
  for (unsigned int i = 0; i < sizeof(norms) / sizeof(norms[0]); ++i)
    {
      VectorTools::integrate_difference(dof_handler,
                                        vec,
                                        Functions::SquareFunction<dim>(),
                                        diff,
                                        q_collection,
                                        norms[i]);
      deallog << "i=" << i << ", diff=" << diff.l2_norm() << std::endl;
    }
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  test<1>();
  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}

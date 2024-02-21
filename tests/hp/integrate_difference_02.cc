// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2020 by the deal.II authors
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
// same random way as in hp/random and on a function that is d-linear


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
      fe_collection.push_back(FE_Q<dim>(i));
      q_collection.push_back(QGauss<dim>(i + 2));
    }


  DoFHandler<dim> dof_handler(tria);

  for (typename DoFHandler<dim>::active_cell_iterator cell =
         dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    cell->set_active_fe_index(Testing::rand() % fe_collection.size());

  dof_handler.distribute_dofs(fe_collection);

  // interpolate a linear function
  Vector<double> vec(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler,
                           Functions::Monomial<dim>(
                             dim == 1 ? Point<dim>(1.) :
                                        (dim == 2 ? Point<dim>(1., 0.) :
                                                    Point<dim>(1., 0., 0.))),
                           vec);

  Vector<float> diff(tria.n_active_cells());

  // L1 norm. the function is u(x)=x, so its
  // L1 norm should be equal to 1/2
  {
    VectorTools::integrate_difference(dof_handler,
                                      vec,
                                      Functions::ZeroFunction<dim>(),
                                      diff,
                                      q_collection,
                                      VectorTools::L1_norm);
    deallog << "L1, diff=" << diff.l1_norm() << std::endl;
  }

  // H1 seminorm. the function is u(x)=x, so
  // its H1 seminorm should be equal to 1
  {
    VectorTools::integrate_difference(dof_handler,
                                      vec,
                                      Functions::ZeroFunction<dim>(),
                                      diff,
                                      q_collection,
                                      VectorTools::H1_seminorm);
    deallog << "H1 seminorm, diff=" << diff.l2_norm() << std::endl;
  }

  // W1infty seminorm. the function is
  // u(x)=x, so the norm must be equal to 1
  // on every cell
  {
    VectorTools::integrate_difference(dof_handler,
                                      vec,
                                      Functions::ZeroFunction<dim>(),
                                      diff,
                                      q_collection,
                                      VectorTools::W1infty_seminorm);
    deallog << "W1infty semi, diff=" << diff.linfty_norm() << std::endl;
    // also ensure that we indeed get the
    // same value on every cell
    diff.add(-1);
    AssertThrow(diff.l2_norm() == 0, ExcInternalError());
  }

  // W1infty norm. the Linfty norm is one, so
  // the W1infty norm must be two. but not on
  // every cell
  {
    VectorTools::integrate_difference(dof_handler,
                                      vec,
                                      Functions::ZeroFunction<dim>(),
                                      diff,
                                      q_collection,
                                      VectorTools::W1infty_norm);
    deallog << "W1infty, diff=" << diff.linfty_norm() << std::endl;
    diff.add(-2);
    AssertThrow(diff.l1_norm() > 0.5, ExcInternalError());
  }
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);
  deallog.get_file_stream().precision(2);

  test<1>();
  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}

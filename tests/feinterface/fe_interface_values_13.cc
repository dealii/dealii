// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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


// Set up a two cells with different elements, such that neither element
// dominates. Create an FEInterfaceValues object and call reinit on it.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/tria.h>

#include "../tests.h"

#include "../test_grids.h"


template <int dim>
void
test()
{
  Triangulation<dim> tria;
  TestGrids::hyper_line<dim>(tria, 2);

  const FESystem<dim>   fe_left(FE_Q<dim>(1), FE_Q<dim>(2));
  const QGauss<dim - 1> q_left(3);

  const FESystem<dim>   fe_right(FE_Q<dim>(3), FE_Q<dim>(1));
  const QGauss<dim - 1> q_right(4);

  const hp::FECollection<dim>    fe_collection(fe_left, fe_right);
  const hp::QCollection<dim - 1> q_collection(q_left, q_right);


  const unsigned int face_index = 1;


  DoFHandler<dim> dofh(tria);

  auto cell     = dofh.begin();
  auto neighbor = cell->neighbor(face_index);

  cell->set_active_fe_index(0);
  neighbor->set_active_fe_index(1);

  dofh.distribute_dofs(fe_collection);


  const UpdateFlags update_flags = update_values | update_gradients |
                                   update_quadrature_points | update_JxW_values;
  FEInterfaceValues<dim> fiv(fe_collection, q_collection, update_flags);

  try
    {
      fiv.reinit(cell,
                 face_index,
                 numbers::invalid_unsigned_int,
                 neighbor,
                 cell->neighbor_of_neighbor(face_index),
                 numbers::invalid_unsigned_int);
    }
  catch (ExceptionBase & /*exc*/)
    {
      deallog << "reinit() failed correctly." << std::endl;
    }

  deallog << "OK" << std::endl;
}



int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();

  test<2>();
}

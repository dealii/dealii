// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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


// evaluate jump_hessian and jump_third_derivative of FEInterfaceValues

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

template <int dim>
void
make_2_cells(Triangulation<dim> &tria);



template <>
void make_2_cells<2>(Triangulation<2> &tria)
{
  const unsigned int        dim         = 2;
  std::vector<unsigned int> repetitions = {2, 1};
  Point<dim>                p1;
  Point<dim>                p2(2.0, 1.0);

  GridGenerator::subdivided_hyper_rectangle(tria, repetitions, p1, p2);
}



template <>
void make_2_cells<3>(Triangulation<3> &tria)
{
  const unsigned int        dim         = 3;
  std::vector<unsigned int> repetitions = {2, 1, 1};
  Point<dim>                p1;
  Point<dim>                p2(2.0, 1.0, 1.0);

  GridGenerator::subdivided_hyper_rectangle(tria, repetitions, p1, p2);
}



template <int dim>
void
print_norm_of_average_over_quadrature_points(const FEInterfaceValues<dim> &fiv)
{
  const unsigned int n_dofs = fiv.n_current_interface_dofs();
  Vector<double>     cell_vector(n_dofs);

  cell_vector = 0.0;
  for (unsigned int qpoint = 0; qpoint < fiv.n_quadrature_points; ++qpoint)
    for (unsigned int i = 0; i < n_dofs; ++i)
      cell_vector(i) +=
        fiv.jump_hessian(i, qpoint).norm() * fiv.get_JxW_values()[qpoint];
  deallog << "jump_hessian.norm(): " << cell_vector << std::endl;

  cell_vector = 0.0;
  for (unsigned int qpoint = 0; qpoint < fiv.n_quadrature_points; ++qpoint)
    for (unsigned int i = 0; i < n_dofs; ++i)
      cell_vector(i) += fiv.jump_3rd_derivative(i, qpoint).norm() *
                        fiv.get_JxW_values()[qpoint];
  deallog << "jump_3rd_derivative.norm(): " << cell_vector << std::endl;
}



template <int dim>
void
test(const FiniteElement<dim> &fe)
{
  Triangulation<dim> tria;
  make_2_cells(tria);

  DoFHandler<dim> dofh(tria);
  deallog << fe.get_name() << std::endl;
  dofh.distribute_dofs(fe);

  MappingQ<dim>     mapping(1);
  const UpdateFlags update_flags =
    update_hessians | update_3rd_derivatives | update_JxW_values;

  FEInterfaceValues<dim> fiv(mapping,
                             fe,
                             QGauss<dim - 1>(fe.degree + 1),
                             update_flags);

  auto cell = dofh.begin();

  // Print jump in Hessian and third derivative for a face at the boundary.
  {
    unsigned int face = 0;
    Assert(cell->at_boundary(face), ExcInternalError());

    fiv.reinit(cell, face);

    print_norm_of_average_over_quadrature_points(fiv);
  }

  // Print jump in Hessian and third derivative for a face between two cells.
  {
    const unsigned int face = 1;
    Assert(!cell->at_boundary(face), ExcInternalError());

    fiv.reinit(cell,
               face,
               numbers::invalid_unsigned_int,
               cell->neighbor(face),
               cell->neighbor_of_neighbor(face),
               numbers::invalid_unsigned_int);

    print_norm_of_average_over_quadrature_points(fiv);
  }
}



int
main()
{
  initlog();
  // Test the lowest order cg and dg elements which have a non-zero third
  // derivative.
  test<2>(FE_Q<2>(2));
  test<2>(FE_DGQ<2>(2));
  test<3>(FE_Q<3>(1));
  test<3>(FE_DGQ<3>(1));
}

// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2021 by the deal.II authors
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

// Test that we can use FE_Bernstein::get_face_interpolation_matrix to create
// an interpolation matrix between FE_Bernstein and an element that has support
// points.

#include <deal.II/fe/fe_bernstein.h>
#include <deal.II/fe/fe_q.h>

#include <vector>

#include "../tests.h"


// Set up an FE_Q element and an FE_Bernstein element with the incoming order.
// Call get_face_interpolation_matrix to create the interpolation matrix
// between these and write it to deallog.
template <int dim>
void
create_and_print_face_interpolation_matrix(const unsigned int order)
{
  deallog << "dim = " << dim << ", order = " << order << std::endl;

  const FE_Q<dim>         fe_q(order);
  const FE_Bernstein<dim> fe_bernstein(order);

  FullMatrix<double> interpolation_matrix(fe_q.n_dofs_per_face(),
                                          fe_bernstein.n_dofs_per_face());
  const unsigned int face_index = 0;


  fe_bernstein.get_face_interpolation_matrix(fe_q,
                                             interpolation_matrix,
                                             face_index);

  interpolation_matrix.print(deallog);
  deallog << std::endl;
}



int
main()
{
  initlog();

  const int dim = 2;

  const std::vector<unsigned int> orders = {1, 2};
  for (const unsigned int order : orders)
    create_and_print_face_interpolation_matrix<dim>(order);
}

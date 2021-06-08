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
// an interpolation matrix when the incoming element is either FE_Bernstein or
// FE_Q.

#include <deal.II/fe/fe_bernstein.h>
#include <deal.II/fe/fe_q.h>

#include <vector>

#include "../tests.h"


// Call get_face_interpolation_matrix to create the interpolation matrix
// between the two incoming elements and write it to deallog.
template <int dim>
void
create_and_print_face_interpolation_matrix(
  const FiniteElement<dim> &fe_destination,
  const FiniteElement<dim> &fe_source)
{
  deallog << fe_source.get_name() << " to " << fe_destination.get_name()
          << std::endl;
  FullMatrix<double> interpolation_matrix(fe_destination.n_dofs_per_face(),
                                          fe_source.n_dofs_per_face());
  const unsigned int face_index = 0;

  fe_source.get_face_interpolation_matrix(fe_destination,
                                          interpolation_matrix,
                                          face_index);

  interpolation_matrix.print(deallog);
  deallog << std::endl;
}



template <int dim>
void
run_test()
{
  const std::vector<unsigned int> orders = {1, 2};
  for (const unsigned int order : orders)
    {
      const FE_Q<dim>         fe_q(order);
      const FE_Bernstein<dim> fe_bernstein(order);
      create_and_print_face_interpolation_matrix<dim>(fe_q, fe_bernstein);
      create_and_print_face_interpolation_matrix<dim>(fe_bernstein,
                                                      fe_bernstein);
    }
}



int
main()
{
  initlog();

  run_test<1>();
  run_test<2>();
}

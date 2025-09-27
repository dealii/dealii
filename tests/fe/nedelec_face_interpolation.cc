// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test by Alexander Grayver: The FE_Nedelec did not correctly compute face
// interpolation matrices from lower to higher order elements. In particular,
// the matrix sometimes did not have full column rank, meaning that a nonzero
// function on the lower-p side would be interpolated to zero on the higher-p
// side.


#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>

#include <iostream>

#include "../tests.h"


template <int dim>
void
test(unsigned p1, unsigned p2)
{
  // create a mesh like this (viewed
  // from top, if in 3d):
  // *----*----*
  // | p1 | p2 |
  // *----*----*
  Triangulation<dim>        triangulation;
  std::vector<unsigned int> subdivisions(dim, 1);
  subdivisions[0] = 2;
  GridGenerator::subdivided_hyper_rectangle(triangulation,
                                            subdivisions,
                                            Point<dim>(),
                                            (dim == 3 ? Point<dim>(2, 1, 1) :
                                                        Point<dim>(2, 1)));
  (std::next(triangulation.begin_active()))->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();

  hp::FECollection<dim> fe;
  fe.push_back(FE_Nedelec<dim>(p1));
  fe.push_back(FE_Nedelec<dim>(p2));

  deallog << "Testing " << fe[0].get_name() << " vs. " << fe[1].get_name()
          << std::endl;

  FullMatrix<double> face_int_matrix(fe[1].dofs_per_face, fe[0].dofs_per_face);

  fe[0].get_face_interpolation_matrix(fe[1], face_int_matrix);

  bool is_zero_column = true;
  for (unsigned int i = 0; i < face_int_matrix.n(); ++i)
    {
      is_zero_column = true;
      for (unsigned int j = 0; j < face_int_matrix.m(); ++j)
        {
          if (fabs(face_int_matrix(j, i)) > 1e-10)
            {
              is_zero_column = false;
              break;
            }
        }

      if (is_zero_column)
        {
          deallog << "Column " << i << " has no non-zero elements."
                  << std::endl;
        }
    }

  if (!is_zero_column)
    deallog << "OK" << std::endl;
  else
    deallog << "Failed" << std::endl;
}

int
main()
{
  initlog();

  test<3>(0, 1);
  test<3>(1, 2);
  // test<3>(2, 3);
  // test<3>(3, 4);
  return 0;
}

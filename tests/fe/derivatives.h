// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_tests_derivatives_h
#define dealii_tests_derivatives_h

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "../tests.h"

char fname[50];



//////////////////////////////////////////////////////////////////////////////////////////////////////
// Plot derivatives of the shape functions at the corner nodes of [0,1]^d
//
// Output values on each line take the form
//
//  x   (y)   (z)   xderiv_index   (yderiv_index)   (zderiv_index)   value[0]+1.
//  value[1]+1.   ...
//////////////////////////////////////////////////////////////////////////////////////////////////////

template <int dim>
inline void
plot_function_derivatives(Mapping<dim>       &mapping,
                          FiniteElement<dim> &finel,
                          const char         *name)
{
  Triangulation<dim> tr;
  DoFHandler<dim>    dof(tr);
  GridGenerator::hyper_cube(tr, 0., 1.);
  dof.distribute_dofs(finel);

  typename DoFHandler<dim>::cell_iterator c = dof.begin();

  const unsigned int sz     = finel.n_dofs_per_cell();
  const unsigned int degree = finel.tensor_degree();
  const unsigned int max_test_deriv =
    std::min(((dim == 1) ? 3 : 1), static_cast<int>(degree - 1) / 2);
  const unsigned int div = 1;

  QTrapezoid<dim> q;
  FEValues<dim>   fe(mapping,
                   finel,
                   q,
                   UpdateFlags(update_values | update_gradients |
                               update_hessians | update_3rd_derivatives));

  sprintf(fname, "Cell-%dd-%s", dim, name);
  deallog.push(fname);

  fe.reinit(c);

  unsigned int k = 0;

  // iterate over vertices
  for (unsigned int mz = 0; mz <= ((dim > 2) ? div : 0); ++mz)
    {
      for (unsigned int my = 0; my <= ((dim > 1) ? div : 0); ++my)
        {
          for (unsigned int mx = 0; mx <= div; ++mx)
            {
              // iterate over derivatives
              for (unsigned int dz = 0; dz <= ((dim > 2) ? max_test_deriv : 0);
                   ++dz)
                {
                  for (unsigned int dy = 0;
                       dy <= ((dim > 1) ? max_test_deriv : 0);
                       ++dy)
                    {
                      for (unsigned int dx = 0; dx <= max_test_deriv; ++dx)
                        {
                          deallog << q.point(k) << " " << dx;
                          if (dim > 1)
                            deallog << " " << dy;
                          if (dim > 2)
                            deallog << " " << dz;

                          for (unsigned int i = 0; i < sz; ++i)
                            {
                              unsigned int deriv_level = dx + dy + dz;
                              unsigned int entry_1, entry_2, entry_3;

                              switch (deriv_level)
                                {
                                  case 0:
                                    deallog << " " << fe.shape_value(i, k) + 1.;
                                    break;

                                  case 1:
                                    entry_1 = dy + 2 * dz;
                                    deallog
                                      << " "
                                      << fe.shape_grad(i, k)[entry_1] + 1.;
                                    break;

                                  case 2:
                                    entry_1 = (dx == 0) ? 1 : 0;
                                    entry_2 = 2 * dz - dy * dz + dy;
                                    deallog
                                      << " "
                                      << fe.shape_hessian(i,
                                                          k)[entry_1][entry_2] +
                                           1.;
                                    break;

                                  case 3:
                                    if (dz == 1)
                                      entry_1 = 0, entry_2 = 1,
                                      entry_3 = 2; // only occurs for dx=dy=dz=1
                                    else if (dx * dy == 0)
                                      entry_1 = entry_2 = entry_3 = dy / 3;
                                    else
                                      entry_1 = 0, entry_2 = dy - 1,
                                      entry_3 = 1;
                                    deallog
                                      << " "
                                      << fe.shape_3rd_derivative(
                                           i, k)[entry_1][entry_2][entry_3] +
                                           1.;
                                    break;

                                  default:
                                    ExcMessage(
                                      "Requested order of derivative is too high for current implementation.");
                                    break;
                                }
                            }
                          deallog << std::endl;
                        }
                    }
                }
              deallog << std::endl;
              ++k;
            }
          deallog << std::endl;
        }
      deallog << std::endl;
    }
  deallog.pop();
}

#endif

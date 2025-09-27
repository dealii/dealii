// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/logstream.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_bdm.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_nedelec_sz.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_bubbles.h>
#include <deal.II/fe/fe_q_dg0.h>
#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../tests.h"

template <typename number>
void
print_formatted(const FullMatrix<number> &A,
                const unsigned int        precision,
                const unsigned int        width)
{
  for (unsigned int i = 0; i < A.m(); ++i)
    {
      for (unsigned int j = 0; j < A.n(); ++j)
        {
          if (A(i, j) != 0)
            deallog << std::setw(width) << std::setprecision(precision)
                    << A(i, j);
          else
            deallog << std::setw(width) << std::setprecision(precision) << "~";
          deallog << ' ';
        };
      deallog << std::endl;
    };
}


template <int dim>
inline void
check_prolongation(FiniteElement<dim> &fe, const char *name)
{
  deallog << name << '<' << dim << '>' << " constraint " << std::endl;
  print_formatted(fe.constraints(), 8, 10);

  for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_cell; ++i)
    {
      deallog << name << '<' << dim << '>' << " prolongation " << i
              << std::endl;
      if (fe.isotropic_prolongation_is_implemented())
        print_formatted(fe.get_prolongation_matrix(i), 8, 8);
    }
}


#define CHECK_ALL(EL, deg, dim)       \
  {                                   \
    FE_##EL<dim> EL(deg);             \
    check_prolongation(EL, #EL #deg); \
  }
#define CHECK_SYS1(sub1, N1, dim)     \
  {                                   \
    FESystem<dim> q(sub1, N1);        \
    check_prolongation(q, #sub1 #N1); \
  }
#define CHECK_SYS2(sub1, N1, sub2, N2, dim)     \
  {                                             \
    FESystem<dim> q(sub1, N1, sub2, N2);        \
    check_prolongation(q, #sub1 #N1 #sub2 #N2); \
  }
#define CHECK_SYS3(sub1, N1, sub2, N2, sub3, N3, dim)     \
  {                                                       \
    FESystem<dim> q(sub1, N1, sub2, N2, sub3, N3);        \
    check_prolongation(q, #sub1 #N1 #sub2 #N2 #sub3 #N3); \
  }

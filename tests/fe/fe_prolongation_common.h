// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_dg0.h>
#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_bdm.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

template <typename number>
void print_formatted (const FullMatrix<number> &A,
                      const unsigned int        precision,
                      const unsigned int        width)
{
  for (unsigned int i=0; i<A.m(); ++i)
    {
      for (unsigned int j=0; j<A.n(); ++j)
        {
          if (A(i,j) != 0)
            deallog << std::setw(width) << std::setprecision(precision)
                    << A(i,j);
          else
            deallog << std::setw(width) << std::setprecision(precision)
                    << "~";
          deallog << ' ';
        };
      deallog << std::endl;
    };
}


template <int dim>
inline void
check_prolongation (FiniteElement<dim> &fe, const char *name)
{
  deallog << name << '<' << dim << '>' << " constraint " << std::endl;
  print_formatted (fe.constraints(), 8, 10);

  for (unsigned int i=0; i<GeometryInfo<dim>::max_children_per_cell; ++i)
    {
      deallog << name << '<' << dim << '>' << " prolongation " << i << std::endl;
      if (fe.isotropic_prolongation_is_implemented())
        print_formatted (fe.get_prolongation_matrix(i), 8, 8);
    }
}


#define CHECK_ALL(EL,deg,dim) { FE_ ## EL<dim> EL(deg); check_prolongation(EL, #EL #deg); }
#define CHECK_SYS1(sub1,N1,dim) { FESystem<dim> q(sub1, N1); check_prolongation(q, #sub1 #N1); }
#define CHECK_SYS2(sub1,N1,sub2,N2,dim) { FESystem<dim> q(sub1, N1, sub2, N2); \
    check_prolongation(q, #sub1 #N1 #sub2 #N2); }
#define CHECK_SYS3(sub1,N1,sub2,N2,sub3,N3,dim) { FESystem<dim> q(sub1, N1, sub2, N2, sub3, N3); \
    check_prolongation(q, #sub1 #N1 #sub2 #N2 #sub3 #N3); }



// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2015 by the deal.II authors
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


// Test the cell matrices generated in FETools and the local renumbering vector.

#include "../tests.h"
#include <iostream>
#include <iomanip>
#include <fstream>

#include <deal.II/base/logstream.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_tools.h>


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



template<int dim>
void test_embedding (const FiniteElement<dim> &fe)
{
  const std::string refine_case_names[8]=
  {
    "no_refinement",
    "cut_x",
    "cut_y",
    "cut_xy",
    "cut_z",
    "cut_xz",
    "cut_yz",
    "cut_xyz"
  };

  const unsigned int n = fe.dofs_per_cell;

  std::vector<std::vector<FullMatrix<double> > > P;
  P.resize(RefinementCase<dim>::isotropic_refinement);
  for (unsigned int ref_case=RefinementCase<dim>::cut_x;
       ref_case<RefinementCase<dim>::isotropic_refinement+1;
       ++ref_case)
    for (unsigned int c=0;
         c<GeometryInfo<dim>::n_children(RefinementCase<dim>(ref_case));
         ++c)
      P[ref_case-1].push_back(FullMatrix<double>(n));

  FETools::compute_embedding_matrices(fe, P);

  for (unsigned int ref_case=RefinementCase<dim>::cut_x;
       ref_case<RefinementCase<dim>::isotropic_refinement+1;
       ++ref_case)
    for (unsigned int c=0;
         c<GeometryInfo<dim>::n_children(RefinementCase<dim>(ref_case));
         ++c)
      {
        deallog << fe.get_name() << " embedding, RefinementCase<dim>:: " << refine_case_names[ref_case] << ", child " << c << std::endl;
        print_formatted(P[ref_case-1][c], 8, 6);
      }
}


template<int dim>
void test_projection (const FiniteElement<dim> &fe1,
                      const FiniteElement<dim> &fe2,
                      std::ostream &out)
{
  out << fe1.get_name() << " -> "
      << fe2.get_name() << std::endl;

  const unsigned int n1 = fe1.dofs_per_cell;
  const unsigned int n2 = fe2.dofs_per_cell;

  FullMatrix<double> P(n2,n1);

  FETools::get_projection_matrix(fe1, fe2, P);
  for (unsigned int i=0; i<n2; ++i)
    for (unsigned int j=0; j<n1; ++j)
      if (std::fabs(P(i,j)) < 1e-12)
        P(i,j) = 0;
  P.print_formatted(out, 8, false, 5);
}


template<int dim>
void test_projection (std::ostream &out)
{
  FE_DGQ<dim> q0(0);
  FE_DGQ<dim> q1(1);
  FE_DGQ<dim> q2(2);
  FE_DGQArbitraryNodes<dim> q3(QIterated<1>(QTrapez<1>(),3));
  FE_DGQArbitraryNodes<dim> q4(QIterated<1>(QTrapez<1>(),4));

  FE_DGP<dim> p0(0);
  FE_DGP<dim> p1(1);
  FE_DGP<dim> p2(2);
  FE_DGP<dim> p3(3);
  FE_DGP<dim> p4(4);

  test_embedding(q0);
  test_embedding(q1);
  test_embedding(q2);
  test_embedding(q3);
  test_embedding(p1);
  test_embedding(p2);
  test_embedding(p3);

  test_projection(p1,p0, out);
  test_projection(p0,p1, out);
  test_projection(p2,p1, out);
  test_projection(p1,p2, out);
  test_projection(p2,p0, out);
  test_projection(p0,p2, out);
  test_projection(p3,p2, out);
  test_projection(p2,p3, out);
  test_projection(p4,p3, out);
  test_projection(p3,p4, out);

  test_projection(q1,q0, out);
  test_projection(q2,q0, out);
  test_projection(q3,q0, out);
  test_projection(q4,q0, out);
  test_projection(q2,q1, out);
  test_projection(q1,q2, out);
  test_projection(q3,q2, out);
  test_projection(q4,q3, out);
}





int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog << std::setprecision(8);
  deallog.threshold_double(1.e-8);

  test_projection<2>(logfile);
  test_projection<3>(logfile);
}

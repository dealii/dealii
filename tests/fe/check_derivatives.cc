// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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


// at a number of quadrature points, evaluate the gradients of shape functions
// and compare it to a finite difference approximation computed using the
// shape values

#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgp_nonparametric.h>
#include <deal.II/fe/fe_dgp_monomial.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <vector>
#include <fstream>
#include <string>


const double delta_x = 1e-8;


template <int dim>
void test (const FiniteElement<dim> &fe,
           const Quadrature<dim>    &quadrature)
{
  deallog << fe.get_name() << ' ' << fe.dofs_per_cell << ' ';

  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
    for (unsigned int q=0; q<quadrature.size(); ++q)
      for (unsigned int c=0; c<fe.n_components(); ++c)
        {
          const Point<dim> point = quadrature.point(q);

          const Tensor<1,dim> gradient
            = fe.shape_grad_component (i, point, c);

          Tensor<1,dim> fd_grad;
          for (unsigned int d=0; d<dim; ++d)
            {
              Point<dim> point_plus_dx = point;
              point_plus_dx[d] += delta_x;
              fd_grad[d] = (fe.shape_value_component (i, point_plus_dx, c) -
                            fe.shape_value_component (i, point, c)) /
                           delta_x;
            }

          Assert ((gradient-fd_grad).norm () <= 2e-5,
                  ExcInternalError());
        }
  deallog << "OK" << std::endl;
}



template <template <int,int> class FE>
void check (const unsigned int min_degree,
            const unsigned int max_degree)
{
  for (unsigned int degree=min_degree; degree<=max_degree; ++degree)
    {
      FE<1,1> fe1(degree);
      test<1> (fe1, QGauss<1> (degree+1));
      FE<2,2> fe2(degree);
      test<2> (fe2, QGauss<2> (degree+1));
      FE<3,3> fe3(degree);
      test<3> (fe3, QGauss<3> (degree+1));
    }
}


template <template <int> class FE>
void check1 (const unsigned int min_degree,
             const unsigned int max_degree)
{
  for (unsigned int degree=min_degree; degree<=max_degree; ++degree)
    {
      FE<1> fe1(degree);
      test<1> (fe1, QGauss<1> (degree+1));
      FE<2> fe2(degree);
      test<2> (fe2, QGauss<2> (degree+1));
      FE<3> fe3(degree);
      test<3> (fe3, QGauss<3> (degree+1));
    }
}


// Nedelec exists only in 2d/3d
template <>
void check1<FE_Nedelec> (const unsigned int min_degree,
                         const unsigned int max_degree)
{
  for (unsigned int degree=min_degree; degree<=max_degree; ++degree)
    {
      test<2> (FE_Nedelec<2> (degree), QGauss<2> (degree+1));
      test<3> (FE_Nedelec<3> (degree), QGauss<3> (degree+1));
    }
}


// Raviart-Thomas doesn't exists 1d. so does the nodal variant of it. the
// former is also not implemented in 3d
template <>
void check1<FE_RaviartThomas> (const unsigned int min_degree,
                               const unsigned int max_degree)
{
  for (unsigned int degree=min_degree; degree<=max_degree; ++degree)
    {
      test<2> (FE_RaviartThomas<2> (degree), QGauss<2> (degree+1));
    }
}

template <>
void check1<FE_RaviartThomasNodal> (const unsigned int min_degree,
                                    const unsigned int max_degree)
{
  for (unsigned int degree=min_degree; degree<=max_degree; ++degree)
    {
      test<2> (FE_RaviartThomasNodal<2> (degree), QGauss<2> (degree+1));
      test<3> (FE_RaviartThomasNodal<3> (degree), QGauss<3> (degree+1));
    }
}



int
main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision(2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);

  check<FE_Q> (1,4);
  check1<FE_Q_Hierarchical> (1,4);
  check<FE_DGQ> (0,4);
  check<FE_DGP> (0,4);
  check<FE_DGPNonparametric> (0,4);
  check1<FE_DGPMonomial> (0,3);

  check1<FE_Nedelec> (0,1);
  check1<FE_RaviartThomas> (0,4);
  check1<FE_RaviartThomasNodal> (0,2);

  return 0;
}




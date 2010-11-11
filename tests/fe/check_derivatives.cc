//----------------------------  check_derivatives.cc  ---------------------------
//    check_derivatives.cc,v 1.3 2003/06/09 16:00:38 wolf Exp
//    Version: 
//
//    Copyright (C) 2005, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  check_derivatives.cc  ---------------------------

// at a number of quadrature points, evaluate the gradients of shape functions
// and compare it to a finite difference approximation computed using the
// shape values

#include "../tests.h"
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <fe/fe_q.h>
#include <fe/fe_q_hierarchical.h>
#include <fe/fe_dgq.h>
#include <fe/fe_dgp.h>
#include <fe/fe_dgp_nonparametric.h>
#include <fe/fe_dgp_monomial.h>
#include <fe/fe_nedelec.h>
#include <fe/fe_raviart_thomas.h>
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
  std::ofstream logfile ("check_derivatives/output");
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




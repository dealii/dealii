// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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


#include <deal.II/base/polynomials_adini.h>

#define ENTER_COEFFICIENTS(koefs,z,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11)       \
  koefs(0, z)= a0; koefs(1, z)=a1; koefs(2, z)=a2; koefs(3, z)=a3; koefs( 4, z)=a4 ; koefs( 5, z)=a5 ; \
  koefs(6, z)= a6; koefs(7, z)=a7; koefs(8, z)=a8; koefs(9, z)=a9; koefs(10, z)=a10; koefs(11, z)=a11;


DEAL_II_NAMESPACE_OPEN


PolynomialsAdini::PolynomialsAdini() :
  coef(12,12),
  dx(12,12),
  dy(12,12),
  dxx(12,12),
  dyy(12,12),
  dxy(12,12)
{
  //                       1  x  y  xx yy xy 3x 3y xyy xxy 3xy x3y
  //                       0  1  2  3  4  5  6  7  8  9 10 11
  ENTER_COEFFICIENTS(coef,  0, 1, 0, 0,-3,-3,-1, 2, 2, 3, 3,-2,-2);
  ENTER_COEFFICIENTS(coef,  1, 0, 1, 0,-2, 0,-1, 1, 0, 0, 2,-1, 0);
  ENTER_COEFFICIENTS(coef,  2, 0, 0, 1, 0,-2,-1, 0, 1, 2, 0, 0,-1);
  ENTER_COEFFICIENTS(coef,  3, 0, 0, 0, 3, 0, 1,-2, 0,-3,-3, 2, 2);
  ENTER_COEFFICIENTS(coef,  4, 0, 0, 0,-1, 0, 0, 1, 0, 0, 1,-1, 0);
  ENTER_COEFFICIENTS(coef,  5, 0, 0, 0, 0, 0, 1, 0, 0,-2, 0, 0, 1);
  ENTER_COEFFICIENTS(coef,  6, 0, 0, 0, 0, 3, 1, 0,-2,-3,-3, 2, 2);
  ENTER_COEFFICIENTS(coef,  7, 0, 0, 0, 0, 0, 1, 0, 0, 0,-2, 1, 0);
  ENTER_COEFFICIENTS(coef,  8, 0, 0, 0, 0,-1, 0, 0, 1, 1, 0, 0,-1);
  ENTER_COEFFICIENTS(coef,  9, 0, 0, 0, 0, 0,-1, 0, 0, 3, 3,-2,-2);
  ENTER_COEFFICIENTS(coef, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0);
  ENTER_COEFFICIENTS(coef, 11, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 1);

  ENTER_COEFFICIENTS( dx,  0, 0,-6,-1, 6, 3, 6, 0,-2, 0,-6, 0, 0);
  ENTER_COEFFICIENTS( dx,  1, 1,-4,-1, 3, 0, 4, 0, 0, 0,-3, 0, 0);
  ENTER_COEFFICIENTS( dx,  2, 0, 0,-1, 0, 2, 0, 0,-1, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dx,  3, 0, 6, 1,-6,-3,-6, 0, 2, 0, 6, 0, 0);
  ENTER_COEFFICIENTS( dx,  4, 0,-2, 0, 3, 0, 2, 0, 0, 0,-3, 0, 0);
  ENTER_COEFFICIENTS( dx,  5, 0, 0, 1, 0,-2, 0, 0, 1, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dx,  6, 0, 0, 1, 0,-3,-6, 0, 2, 0, 6, 0, 0);
  ENTER_COEFFICIENTS( dx,  7, 0, 0, 1, 0, 0,-4, 0, 0, 0, 3, 0, 0);
  ENTER_COEFFICIENTS( dx,  8, 0, 0, 0, 0, 1, 0, 0,-1, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dx,  9, 0, 0,-1, 0, 3, 6, 0,-2, 0,-6, 0, 0);
  ENTER_COEFFICIENTS( dx, 10, 0, 0, 0, 0, 0,-2, 0, 0, 0, 3, 0, 0);
  ENTER_COEFFICIENTS( dx, 11, 0, 0, 0, 0,-1, 0, 0, 1, 0, 0, 0, 0);

  ENTER_COEFFICIENTS( dy,  0, 0,-1,-6, 3, 6, 6,-2, 0,-6, 0, 0, 0);
  ENTER_COEFFICIENTS( dy,  1, 0,-1, 0, 2, 0, 0,-1, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dy,  2, 1,-1,-4, 0, 3, 4, 0, 0,-3, 0, 0, 0);
  ENTER_COEFFICIENTS( dy,  3, 0, 1, 0,-3, 0,-6, 2, 0, 6, 0, 0, 0);
  ENTER_COEFFICIENTS( dy,  4, 0, 0, 0, 1, 0, 0,-1, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dy,  5, 0, 1, 0, 0, 0,-4, 0, 0, 3, 0, 0, 0);
  ENTER_COEFFICIENTS( dy,  6, 0, 1, 6,-3,-6,-6, 2, 0, 6, 0, 0, 0);
  ENTER_COEFFICIENTS( dy,  7, 0, 1, 0,-2, 0, 0, 1, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dy,  8, 0, 0,-2, 0, 3, 2, 0, 0,-3, 0, 0, 0);
  ENTER_COEFFICIENTS( dy,  9, 0,-1, 0, 3, 0, 6,-2, 0,-6, 0, 0, 0);
  ENTER_COEFFICIENTS( dy, 10, 0, 0, 0,-1, 0, 0, 1, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dy, 11, 0, 0, 0, 0, 0,-2, 0, 0, 3, 0, 0, 0);

  //                       0  1  2  3  4  5  6  7  8  9 10 11
  ENTER_COEFFICIENTS( dxx,  0,-6, 12, 6, 0, 0,-12, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxx,  1,-4,  6, 4, 0, 0, -6, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxx,  2, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxx,  3, 6,-12,-6, 0, 0, 12, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxx,  4,-2,  6, 2, 0, 0, -6, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxx,  5, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxx,  6, 0,  0,-6, 0, 0, 12, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxx,  7, 0,  0,-4, 0, 0,  6, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxx,  8, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxx,  9, 0,  0, 6, 0, 0,-12, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxx, 10, 0,  0,-2, 0, 0,  6, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxx, 11, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0);

  ENTER_COEFFICIENTS( dyy,  0,-6, 6, 12, 0, 0,-12, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dyy,  1, 0, 0,  0, 0, 0,  0, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dyy,  2,-4, 4,  6, 0, 0, -6, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dyy,  3, 0,-6,  0, 0, 0, 12, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dyy,  4, 0, 0,  0, 0, 0,  0, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dyy,  5, 0,-4,  0, 0, 0,  6, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dyy,  6, 6,-6,-12, 0, 0, 12, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dyy,  7, 0, 0,  0, 0, 0,  0, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dyy,  8,-2, 2,  6, 0, 0, -6, 0,-0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dyy,  9, 0, 6,  0, 0, 0,-12, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dyy, 10, 0, 0,  0, 0, 0,  0, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dyy, 11, 0,-2,  0, 0, 0,  6, 0, 0, 0, 0, 0, 0);

  ENTER_COEFFICIENTS( dxy,  0,-1, 6, 6,-6,-6, 0, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxy,  1,-1, 4, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxy,  2,-1, 0, 4, 0,-3, 0, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxy,  3, 1,-6,-6, 6, 6, 0, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxy,  4, 0, 2, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxy,  5, 1, 0,-4, 0, 3, 0, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxy,  6, 1,-6,-6, 6, 6, 0, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxy,  7, 1,-4, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxy,  8, 0, 0, 2, 0,-3, 0, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxy,  9,-1, 6, 6,-6,-6, 0, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxy, 10, 0,-2, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0);
  ENTER_COEFFICIENTS( dxy, 11, 0, 0,-2, 0, 3, 0, 0, 0, 0, 0, 0, 0);
}

void
PolynomialsAdini::compute (const Point<2> &unit_point,
                           std::vector<double> &values,
                           std::vector<Tensor<1,2> > &grads,
                           std::vector< Tensor<2,2> > &grad_grads) const
{
  if (values.empty() == false) //nur dann was tun
    {
      for (unsigned int i=0; i<values.size(); ++i)
        {
          values[i] = compute_value(i, unit_point);
        }
    }

  if (grads.empty() == false) //nur dann was tun
    {
      for (unsigned int i=0; i<grads.size(); ++i)
        {
          grads[i] = compute_grad(i, unit_point);
        }
    }

  if (grad_grads.empty() == false) //nur dann was tun
    {
      for (unsigned int i=0; i<grad_grads.size(); ++i)
        {
          grad_grads[i] = compute_grad_grad(i, unit_point);
        }
    }
  return;
}

double PolynomialsAdini::compute_value (const unsigned int i, const Point<2> &p) const
{
  const double x = p(0);
  const double y = p(1);
  return coef( 0,i)           +
         coef( 1,i)*x         +
         coef( 2,i)*y         +
         coef( 3,i)*x*x       +
         coef( 4,i)*y*y       +
         coef( 5,i)*x*y       +
         coef( 6,i)*x*x*x     +
         coef( 7,i)*y*y*y     +
         coef( 8,i)*x*y*y     +
         coef( 9,i)*x*x*y     +
         coef(10,i)*x*x*x*y  +
         coef(11,i)*x*y*y*y;
}

Tensor<1,2> PolynomialsAdini::compute_grad (const unsigned int i, const Point<2> &p) const
{
  const double x = p(0);
  const double y = p(1);
  Tensor<1,2> tensor;
  tensor[0]= dx( 0,i)          +
             dx( 1,i)*x        +
             dx( 2,i)*y        +
             dx( 3,i)*x*x      +
             dx( 4,i)*y*y      +
             dx( 5,i)*x*y      +
             dx( 6,i)*x*x*x    +
             dx( 7,i)*y*y*y    +
             dx( 8,i)*x*y*y    +
             dx( 9,i)*x*x*y    +
             dx(10,i)*x*x*x*y +
             dx(11,i)*x*y*y*y;

  tensor[1]= dy( 0,i)          +
             dy( 1,i)*x        +
             dy( 2,i)*y        +
             dy( 3,i)*x*x      +
             dy( 4,i)*y*y      +
             dy( 5,i)*x*y      +
             dy( 6,i)*x*x*x    +
             dy( 7,i)*y*y*y    +
             dy( 8,i)*x*y*y    +
             dy( 9,i)*x*x*y    +
             dy(10,i)*x*x*x*y +
             dy(11,i)*x*y*y*y;
  return tensor;
}

Tensor<2,2> PolynomialsAdini::compute_grad_grad (const unsigned int i, const Point<2> &p) const
{
  const double x = p(0);
  const double y = p(1);
  Tensor<2,2> tensor;
  tensor[0][0]= dxx( 0,i)          +
                dxx( 1,i)*x        +
                dxx( 2,i)*y        +
                dxx( 3,i)*x*x      +
                dxx( 4,i)*y*y      +
                dxx( 5,i)*x*y      +
                dxx( 6,i)*x*x*x    +
                dxx( 7,i)*y*y*y    +
                dxx( 8,i)*x*y*y    +
                dxx( 9,i)*x*x*y    +
                dxx(10,i)*x*x*x*y +
                dxx(11,i)*x*y*y*y;
  tensor[0][1]= dxy( 0,i)          +
                dxy( 1,i)*x        +
                dxy( 2,i)*y        +
                dxy( 3,i)*x*x      +
                dxy( 4,i)*y*y      +
                dxy( 5,i)*x*y      +
                dxy( 6,i)*x*x*x    +
                dxy( 7,i)*y*y*y    +
                dxy( 8,i)*x*y*y    +
                dxy( 9,i)*x*x*y    +
                dxy(10,i)*x*x*x*y +
                dxy(11,i)*x*y*y*y;
  tensor[1][0]= tensor[0][1];
  tensor[1][1]= dyy( 0,i)          +
                dyy( 1,i)*x        +
                dyy( 2,i)*y        +
                dyy( 3,i)*x*x      +
                dyy( 4,i)*y*y      +
                dyy( 5,i)*x*y      +
                dyy( 6,i)*x*x*x    +
                dyy( 7,i)*y*y*y    +
                dyy( 8,i)*x*y*y    +
                dyy( 9,i)*x*x*y    +
                dyy(10,i)*x*x*x*y +
                dyy(11,i)*x*y*y*y;
  return tensor;
}


DEAL_II_NAMESPACE_CLOSE

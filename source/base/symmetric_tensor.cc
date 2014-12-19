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

#include <deal.II/base/symmetric_tensor.h>

DEAL_II_NAMESPACE_OPEN


template <>
SymmetricTensor<4,3,double>
invert<3,double> (const SymmetricTensor<4,3,double> &t)
{
  SymmetricTensor<4,3,double> tmp = t;

  // this function follows the exact same
  // scheme as the 2d case, except that
  // hardcoding the inverse of a 6x6 matrix
  // is pretty wasteful. instead, we use the
  // Gauss-Jordan algorithm implemented for
  // FullMatrix; the following code is copied
  // from there because using the FullMatrix
  // class would introduce circular
  // references between libbase and liblac
  const unsigned int N = 6;

  // first get an estimate of the
  // size of the elements of this
  // matrix, for later checks whether
  // the pivot element is large
  // enough, or whether we have to
  // fear that the matrix is not
  // regular
  double diagonal_sum = 0;
  for (unsigned int i=0; i<N; ++i)
    diagonal_sum += std::fabs(tmp.data[i][i]);
  const double typical_diagonal_element = diagonal_sum/N;

  unsigned int p[N];
  for (unsigned int i=0; i<N; ++i)
    p[i] = i;

  for (unsigned int j=0; j<N; ++j)
    {
      // pivot search: search that
      // part of the line on and
      // right of the diagonal for
      // the largest element
      double       max = std::fabs(tmp.data[j][j]);
      unsigned int r   = j;
      for (unsigned int i=j+1; i<N; ++i)
        if (std::fabs(tmp.data[i][j]) > max)
          {
            max = std::fabs(tmp.data[i][j]);
            r = i;
          }
      // check whether the pivot is
      // too small
      Assert(max > 1.e-16*typical_diagonal_element,
             ExcMessage("This tensor seems to be noninvertible"));

      // row interchange
      if (r>j)
        {
          for (unsigned int k=0; k<N; ++k)
            std::swap (tmp.data[j][k], tmp.data[r][k]);

          std::swap (p[j], p[r]);
        }

      // transformation
      const double hr = 1./tmp.data[j][j];
      tmp.data[j][j] = hr;
      for (unsigned int k=0; k<N; ++k)
        {
          if (k==j) continue;
          for (unsigned int i=0; i<N; ++i)
            {
              if (i==j) continue;
              tmp.data[i][k] -= tmp.data[i][j]*tmp.data[j][k]*hr;
            }
        }
      for (unsigned int i=0; i<N; ++i)
        {
          tmp.data[i][j] *= hr;
          tmp.data[j][i] *= -hr;
        }
      tmp.data[j][j] = hr;
    }
  // column interchange
  double hv[N];
  for (unsigned int i=0; i<N; ++i)
    {
      for (unsigned int k=0; k<N; ++k)
        hv[p[k]] = tmp.data[i][k];
      for (unsigned int k=0; k<N; ++k)
        tmp.data[i][k] = hv[k];
    }

  // scale rows and columns. the mult matrix
  // here is diag[1, 1, 1, 1/2, 1/2, 1/2]
  for (unsigned int i=3; i<6; ++i)
    for (unsigned int j=0; j<3; ++j)
      tmp.data[i][j] /= 2;

  for (unsigned int i=0; i<3; ++i)
    for (unsigned int j=3; j<6; ++j)
      tmp.data[i][j] /= 2;

  for (unsigned int i=3; i<6; ++i)
    for (unsigned int j=3; j<6; ++j)
      tmp.data[i][j] /= 4;

  return tmp;
}

DEAL_II_NAMESPACE_CLOSE

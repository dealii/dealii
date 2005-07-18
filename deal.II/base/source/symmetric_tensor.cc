//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/symmetric_tensor.h>
#include <lac/full_matrix.h>

template <>
SymmetricTensor<4,3>
invert (const SymmetricTensor<4,3> &t)
{
  SymmetricTensor<4,3> tmp;

                                   // this function follows the exact same
                                   // scheme as the 2d case, except that
                                   // hardcoding the inverse of a 6x6 matrix
                                   // is pretty wasteful. instead, we use the
                                   // Gauss-Jordan algorithm implemented for
                                   // FullMatrix
  FullMatrix<double> m(6,6);
  for (unsigned int i=0; i<6; ++i)
    for (unsigned int j=0; j<6; ++j)
      m(i,j) = t.data[i][j];
  m.gauss_jordan ();

                                   // copy back and scale rows and
                                   // columns. the mult matrix here is diag[1,
                                   // 1, 1, 1/2, 1/2, 1/2]
  for (unsigned int i=0; i<3; ++i)
    for (unsigned int j=0; j<3; ++j)
      tmp.data[i][j] = m(i,j);

  for (unsigned int i=3; i<6; ++i)
    for (unsigned int j=0; j<3; ++j)
      tmp.data[i][j] = m(i,j) / 2;

  for (unsigned int i=0; i<3; ++i)
    for (unsigned int j=3; j<6; ++j)
      tmp.data[i][j] = m(i,j) / 2;

  for (unsigned int i=3; i<6; ++i)
    for (unsigned int j=3; j<6; ++j)
      tmp.data[i][j] = m(i,j) / 4;
  
  return tmp;
}

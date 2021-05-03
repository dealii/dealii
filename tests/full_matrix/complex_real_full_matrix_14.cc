// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// check FullMatrix::matrix_scalar_product. like the full_matrix_* tests, but
// use complex-valued matrices and vectors, even though we only store real
// values in them


#include "../tests.h"

#include "full_matrix_common.h"



template <typename number>
void
check()
{
  FullMatrix<std::complex<number>> m;
  make_matrix(m);
  Vector<std::complex<number>> v, w;
  make_range_vector(v);
  make_domain_vector(w);
  for (unsigned int i = 0; i < w.size(); ++i)
    w(i) = w(i) + std::complex<number>(1.);

  deallog << m.matrix_scalar_product(v, w) << std::endl;
}

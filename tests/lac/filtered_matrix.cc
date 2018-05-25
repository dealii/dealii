// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2017 by the deal.II authors
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


// Test properties of FilteredMatrix and iterators

#include <deal.II/lac/filtered_matrix.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/matrix_lib.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"


template <typename VectorType>
void
test(const FilteredMatrix<VectorType> &M)
{
  deallog << "Iterator";

  unsigned int max = 0;
  for (typename FilteredMatrix<VectorType>::const_iterator i = M.begin();
       i != M.end();
       ++i)
    {
      Assert(i->row() == i->column(), ExcInternalError());
      deallog << ' ' << i->row() << ':' << i->value();
      max = i->row();
    }
  VectorType v(max + 1);
  VectorType w(max + 1);

  for (unsigned int i = 0; i < v.size(); ++i)
    v(i) = 31 + i;

  deallog << std::endl << "vmult ";

  w = 0.;
  M.vmult(w, v);
  for (unsigned int i = 0; i < v.size(); ++i)
    deallog << ' ' << w(i);

  deallog << std::endl << "Tvmult";

  w = 0.;
  M.Tvmult(w, v);
  for (unsigned int i = 0; i < v.size(); ++i)
    deallog << ' ' << w(i);

  deallog << std::endl << "vmult_add";

  M.vmult_add(w, v);
  for (unsigned int i = 0; i < v.size(); ++i)
    deallog << ' ' << w(i);

  deallog << std::endl << "boundary";

  M.apply_constraints(w, true);
  for (unsigned int i = 0; i < v.size(); ++i)
    deallog << ' ' << w(i);

  deallog << std::endl << "boundary";

  M.apply_constraints(w, false);
  for (unsigned int i = 0; i < v.size(); ++i)
    deallog << ' ' << w(i);

  deallog << std::endl;
}

int
main()
{
  initlog();

  PreconditionIdentity           identity;
  FilteredMatrix<Vector<double>> f;
  f.initialize(identity, false);
  test(f);

  for (unsigned int i = 0; i < 5; ++i)
    f.add_constraint(i * i, i / 2.);
  test(f);

  f.initialize(identity, true);
  test(f);
}


// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2017 by the deal.II authors
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



// check DynamicSparsityPattern::compute_mmult_pattern(). Test if multiplication
// with each combination of SparsityPatternTypes is possible.

#include "../tests.h"
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>


void
test ()
{
  // create different (dynamic) sparsity pattern and add entries at
  // destinct locations. Check if mmultiplication creates entries at the right
  // rows and columns.
  const unsigned int M = 100;
  const unsigned int N = 50;
  const unsigned int O = 75;

  const unsigned int m_entries = 4;
  const unsigned int n_entries = 5;
  const unsigned int o_entries = 6;

  bool test_failed = false;

  DynamicSparsityPattern dyn_left (M,N),
                         dyn_right (N,O);

  DynamicSparsityPattern out;

  SparsityPattern left,
                  right;

  // add manually entries to pattern
  for (unsigned int i=0; i<M; i+=m_entries)
    for (unsigned int j=0; j<N; j+=n_entries)
      dyn_left.add (i,j);
  // add manually entries to pattern
  for (unsigned int i=0; i<N; i+=n_entries)
    for (unsigned int j=0; j<O; j+=o_entries)
      dyn_right.add (i,j);
  // copy DynamicSparsityPattern to SparsityPattern
  left.copy_from(dyn_left);
  right.copy_from(dyn_right);
  // test DynamicSparsityPattern and DynamicSparsityPattern
  out.compute_mmult_pattern(dyn_left,dyn_right);
  for (unsigned int i=0; i<M; i+=m_entries)
    for (unsigned int j=0; j<O; j+=o_entries)
      if (!out.exists(i,j))
        test_failed = true;
  // test SparsityPattern and DynamicSparsityPattern
  out.compute_mmult_pattern(left,dyn_right);
  for (unsigned int i=0; i<M; i+=m_entries)
    for (unsigned int j=0; j<O; j+=o_entries)
      if (!out.exists(i,j))
        test_failed = true;
  // test DynamicSparsityPattern and SparsityPattern
  out.compute_mmult_pattern(dyn_left,right);
  for (unsigned int i=0; i<M; i+=m_entries)
    for (unsigned int j=0; j<O; j+=o_entries)
      if (!out.exists(i,j))
        test_failed = true;
  // test SparsityPattern and SparsityPattern
  out.compute_mmult_pattern(left,right);
  for (unsigned int i=0; i<M; i+=m_entries)
    for (unsigned int j=0; j<O; j+=o_entries)
      if (!out.exists(i,j))
        test_failed = true;

  Assert(!test_failed, ExcInternalError());
  deallog << "OK" << std::endl;
}



int
main ()
{
  initlog();

  test ();
  return 0;
}

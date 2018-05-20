// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2017 by the deal.II authors
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

#ifndef dealii_matrix_lib_templates_h
#define dealii_matrix_lib_templates_h

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/matrix_lib.h>
#include <deal.II/lac/vector.h>

DEAL_II_NAMESPACE_OPEN

template <typename number>
void
MeanValueFilter::filter(Vector<number>& v) const
{
  number mean = v.mean_value();

  for(size_type i = 0; i < v.size(); ++i)
    v(i) -= mean;
}

template <typename number>
void
MeanValueFilter::vmult(Vector<number>& dst, const Vector<number>& src) const
{
  Assert(dst.size() == src.size(),
         ExcDimensionMismatch(dst.size(), src.size()));

  number mean = src.mean_value();

  for(size_type i = 0; i < dst.size(); ++i)
    dst(i) = src(i) - mean;
}

template <typename number>
void
MeanValueFilter::vmult_add(Vector<number>& dst, const Vector<number>& src) const
{
  Assert(dst.size() == src.size(),
         ExcDimensionMismatch(dst.size(), src.size()));

  number mean = src.mean_value();

  for(size_type i = 0; i < dst.size(); ++i)
    dst(i) += src(i) - mean;
}

template <typename number>
void
MeanValueFilter::filter(BlockVector<number>& v) const
{
  Assert(component != numbers::invalid_unsigned_int, ExcNotInitialized());

  for(unsigned int i = 0; i < v.n_blocks(); ++i)
    if(i == component)
      vmult(v.block(i), v.block(i));
}

template <typename number>
void
MeanValueFilter::vmult(BlockVector<number>&       dst,
                       const BlockVector<number>& src) const
{
  Assert(component != numbers::invalid_unsigned_int, ExcNotInitialized());

  Assert(dst.n_blocks() == src.n_blocks(),
         ExcDimensionMismatch(dst.n_blocks(), src.n_blocks()));

  for(unsigned int i = 0; i < dst.n_blocks(); ++i)
    if(i == component)
      vmult(dst.block(i), src.block(i));
    else
      dst.block(i) = src.block(i);
}

template <typename number>
void
MeanValueFilter::vmult_add(BlockVector<number>&       dst,
                           const BlockVector<number>& src) const
{
  Assert(component != numbers::invalid_unsigned_int, ExcNotInitialized());

  Assert(dst.n_blocks() == src.n_blocks(),
         ExcDimensionMismatch(dst.n_blocks(), src.n_blocks()));

  for(unsigned int i = 0; i < dst.n_blocks(); ++i)
    if(i == component)
      vmult_add(dst.block(i), src.block(i));
    else
      dst.block(i) += src.block(i);
}

DEAL_II_NAMESPACE_CLOSE

#endif

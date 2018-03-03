// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

#include <deal.II/base/std_cxx14/memory.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_banded_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/lapack_support.h>
#include <deal.II/lac/lapack_templates.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/utilities.h>
#include <deal.II/lac/vector.h>

DEAL_II_NAMESPACE_OPEN

// ----------------------------------------------------------------------------
// LAPACKBandedMatrixData
// ----------------------------------------------------------------------------

LAPACKBandedMatrixData::LAPACKBandedMatrixData() :
  LAPACKBandedMatrixData(0, 0, 0, 0)
{}



LAPACKBandedMatrixData::LAPACKBandedMatrixData(
  const size_type n_rows,
  const size_type n_cols,
  const size_type n_subdiagonals,
  const size_type n_superdiagonals) :
  n_rows(n_rows),
  n_cols(n_cols),
  n_subdiagonals(n_subdiagonals),
  n_superdiagonals(n_superdiagonals),
  n_allocated_superdiagonals(n_superdiagonals + n_subdiagonals),
  n_allocated_rows(n_allocated_superdiagonals + n_subdiagonals + 1)
{}

// ----------------------------------------------------------------------------
// construction, copying, swapping, and reinitialization
// ----------------------------------------------------------------------------

template <typename Number>
LAPACKBandedMatrix<Number>::LAPACKBandedMatrix(
  const size_type n_rows,
  const size_type n_cols,
  const size_type n_subdiagonals,
  const size_type n_superdiagonals) :
  LAPACKBandedMatrixData(n_rows, n_cols, n_subdiagonals, n_superdiagonals)
{
  values.resize(n_allocated_rows * n_cols);
}



template <typename Number>
void
LAPACKBandedMatrix<Number>::reinit(const size_type n_rows,
                                   const size_type n_cols,
                                   const size_type n_subdiagonals,
                                   const size_type n_superdiagonals)
{
  LAPACKBandedMatrix<Number> other(
    n_rows, n_cols, n_subdiagonals, n_superdiagonals);
  swap(other);
}



template <typename Number>
void
LAPACKBandedMatrix<Number>::swap(LAPACKBandedMatrix<Number> &other)
{
  std::swap(*this, other);
}



#include "lapack_banded_matrix.inst"


DEAL_II_NAMESPACE_CLOSE

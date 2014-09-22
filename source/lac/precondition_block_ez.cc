// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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

#include <deal.II/lac/precondition_block.templates.h>
#include <deal.II/lac/sparse_matrix_ez.h>

DEAL_II_NAMESPACE_OPEN


// explicit instantiations for "float" PreconditionBlock
template class PreconditionBlock<SparseMatrixEZ<float>, float>;

// the instantiation for class PreconditionBlock<SparseMatrixEZ<float>, double> is skipped
// because it does not make sense to have inverse block matrices with
// higher precision than the matrix itself


// explicit instantiations for "double" PreconditionBlock
template class PreconditionBlock<SparseMatrixEZ<double>, float>;

template class PreconditionBlock<SparseMatrixEZ<double>, double>;


/*--------------------- PreconditionBlockJacobi -----------------------*/


// explicit instantiations for "float" PreconditionBlock
template class PreconditionBlockJacobi<SparseMatrixEZ<float>, float>;

template void PreconditionBlockJacobi<SparseMatrixEZ<float>, float>::vmult<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrixEZ<float>, float>::vmult<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<SparseMatrixEZ<float>, float>::Tvmult<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrixEZ<float>, float>::Tvmult<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<SparseMatrixEZ<float>, float>::vmult_add<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrixEZ<float>, float>::vmult_add<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<SparseMatrixEZ<float>, float>::Tvmult_add<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrixEZ<float>, float>::Tvmult_add<double>
(Vector<double> &, const Vector<double> &) const;

template class PreconditionBlockJacobi<SparseMatrixEZ<double>, float>;

template void PreconditionBlockJacobi<SparseMatrixEZ<double>, float>::vmult<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrixEZ<double>, float>::vmult<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<SparseMatrixEZ<double>, float>::Tvmult<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrixEZ<double>, float>::Tvmult<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<SparseMatrixEZ<double>, float>::vmult_add<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrixEZ<double>, float>::vmult_add<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<SparseMatrixEZ<double>, float>::Tvmult_add<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrixEZ<double>, float>::Tvmult_add<double>
(Vector<double> &, const Vector<double> &) const;

template class PreconditionBlockJacobi<SparseMatrixEZ<double>, double>;

template void PreconditionBlockJacobi<SparseMatrixEZ<double>, double>::vmult<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrixEZ<double>, double>::vmult<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<SparseMatrixEZ<double>, double>::Tvmult<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrixEZ<double>, double>::Tvmult<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<SparseMatrixEZ<double>, double>::vmult_add<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrixEZ<double>, double>::vmult_add<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<SparseMatrixEZ<double>, double>::Tvmult_add<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrixEZ<double>, double>::Tvmult_add<double>
(Vector<double> &, const Vector<double> &) const;

/*--------------------- PreconditionBlockGaussSeidel -----------------------*/


// explicit instantiations for "float" PreconditionBlock
template class PreconditionBlockSOR<SparseMatrixEZ<float>, float>;

template void PreconditionBlockSOR<SparseMatrixEZ<float>, float>::vmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<SparseMatrixEZ<float>, float>::vmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSOR<SparseMatrixEZ<float>, float>::Tvmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<SparseMatrixEZ<float>, float>::Tvmult<double> (
  Vector<double> &, const Vector<double> &) const;


// the instantiation for class PreconditionBlockSOR<SparseMatrixEZ<float>, double> is skipped
// because it does not make sense to have inverse block matrices with
// higher precision than the matrix itself


// explicit instantiations for "double" PreconditionBlockSOR
template class PreconditionBlockSOR<SparseMatrixEZ<double>, float>;


template void PreconditionBlockSOR<SparseMatrixEZ<double>, float>::vmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<SparseMatrixEZ<double>, float>::vmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSOR<SparseMatrixEZ<double>, float>::Tvmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<SparseMatrixEZ<double>, float>::Tvmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSOR<SparseMatrixEZ<double>, float>::vmult_add<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<SparseMatrixEZ<double>, float>::vmult_add<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSOR<SparseMatrixEZ<double>, float>::Tvmult_add<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<SparseMatrixEZ<double>, float>::Tvmult_add<double> (
  Vector<double> &, const Vector<double> &) const;

template class PreconditionBlockSOR<SparseMatrixEZ<double>, double>;

template void PreconditionBlockSOR<SparseMatrixEZ<double>, double>::vmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<SparseMatrixEZ<double>, double>::vmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSOR<SparseMatrixEZ<double>, double>::Tvmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<SparseMatrixEZ<double>, double>::Tvmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSOR<SparseMatrixEZ<double>, double>::vmult_add<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<SparseMatrixEZ<double>, double>::vmult_add<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSOR<SparseMatrixEZ<double>, double>::Tvmult_add<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<SparseMatrixEZ<double>, double>::Tvmult_add<double> (
  Vector<double> &, const Vector<double> &) const;


/*--------------------- PreconditionBlockSSOR -----------------------*/


// explicit instantiations for "float" PreconditionBlock
template class PreconditionBlockSSOR<SparseMatrixEZ<float>, float>;

template void PreconditionBlockSSOR<SparseMatrixEZ<float>, float>::vmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSSOR<SparseMatrixEZ<float>, float>::vmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSSOR<SparseMatrixEZ<float>, float>::Tvmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSSOR<SparseMatrixEZ<float>, float>::Tvmult<double> (
  Vector<double> &, const Vector<double> &) const;


// the instantiation for class PreconditionBlockSSOR<SparseMatrixEZ<float>, double> is skipped
// because it does not make sense to have inverse block matrices with
// higher precision than the matrix itself


// explicit instantiations for "double" PreconditionBlockSSOR
template class PreconditionBlockSSOR<SparseMatrixEZ<double>, float>;


template void PreconditionBlockSSOR<SparseMatrixEZ<double>, float>::vmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSSOR<SparseMatrixEZ<double>, float>::vmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSSOR<SparseMatrixEZ<double>, float>::Tvmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSSOR<SparseMatrixEZ<double>, float>::Tvmult<double> (
  Vector<double> &, const Vector<double> &) const;

template class PreconditionBlockSSOR<SparseMatrixEZ<double>, double>;

template void PreconditionBlockSSOR<SparseMatrixEZ<double>, double>::vmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSSOR<SparseMatrixEZ<double>, double>::vmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSSOR<SparseMatrixEZ<double>, double>::Tvmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSSOR<SparseMatrixEZ<double>, double>::Tvmult<double> (
  Vector<double> &, const Vector<double> &) const;

DEAL_II_NAMESPACE_CLOSE

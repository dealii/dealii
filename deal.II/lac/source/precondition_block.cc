//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


#include <lac/precondition_block.templates.h>
#include <lac/sparse_matrix.h>


// explicit instantiations for "float" PreconditionBlock
template class PreconditionBlock<SparseMatrix<float>, float>;

// the instantiation for class PreconditionBlock<SparseMatrix<float>, double> is skipped
// because it does not make sense to have inverse block matrices with
// higher precision than the matrix itself


// explicit instantiations for "double" PreconditionBlock
template class PreconditionBlock<SparseMatrix<double>, float>;

template class PreconditionBlock<SparseMatrix<double>, double>;


/*--------------------- PreconditionBlockJacobi -----------------------*/


// explicit instantiations for "float" PreconditionBlock
template class PreconditionBlockJacobi<SparseMatrix<float>, float>;

template void PreconditionBlockJacobi<SparseMatrix<float>, float>::vmult<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrix<float>, float>::vmult<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<SparseMatrix<float>, float>::Tvmult<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrix<float>, float>::Tvmult<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<SparseMatrix<float>, float>::vmult_add<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrix<float>, float>::vmult_add<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<SparseMatrix<float>, float>::Tvmult_add<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrix<float>, float>::Tvmult_add<double>
(Vector<double> &, const Vector<double> &) const;

template class PreconditionBlockJacobi<SparseMatrix<double>, float>;

template void PreconditionBlockJacobi<SparseMatrix<double>, float>::vmult<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrix<double>, float>::vmult<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<SparseMatrix<double>, float>::Tvmult<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrix<double>, float>::Tvmult<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<SparseMatrix<double>, float>::vmult_add<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrix<double>, float>::vmult_add<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<SparseMatrix<double>, float>::Tvmult_add<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrix<double>, float>::Tvmult_add<double>
(Vector<double> &, const Vector<double> &) const;

template class PreconditionBlockJacobi<SparseMatrix<double>, double>;

template void PreconditionBlockJacobi<SparseMatrix<double>, double>::vmult<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrix<double>, double>::vmult<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<SparseMatrix<double>, double>::Tvmult<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrix<double>, double>::Tvmult<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<SparseMatrix<double>, double>::vmult_add<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrix<double>, double>::vmult_add<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<SparseMatrix<double>, double>::Tvmult_add<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<SparseMatrix<double>, double>::Tvmult_add<double>
(Vector<double> &, const Vector<double> &) const;

/*--------------------- PreconditionBlockGaussSeidel -----------------------*/


// explicit instantiations for "float" PreconditionBlock
template class PreconditionBlockSOR<SparseMatrix<float>, float>;

template void PreconditionBlockSOR<SparseMatrix<float>, float>::vmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<SparseMatrix<float>, float>::vmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSOR<SparseMatrix<float>, float>::Tvmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<SparseMatrix<float>, float>::Tvmult<double> (
  Vector<double> &, const Vector<double> &) const;


// the instantiation for class PreconditionBlockSOR<SparseMatrix<float>, double> is skipped
// because it does not make sense to have inverse block matrices with
// higher precision than the matrix itself


// explicit instantiations for "double" PreconditionBlockSOR
template class PreconditionBlockSOR<SparseMatrix<double>, float>;


template void PreconditionBlockSOR<SparseMatrix<double>, float>::vmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<SparseMatrix<double>, float>::vmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSOR<SparseMatrix<double>, float>::Tvmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<SparseMatrix<double>, float>::Tvmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSOR<SparseMatrix<double>, float>::vmult_add<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<SparseMatrix<double>, float>::vmult_add<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSOR<SparseMatrix<double>, float>::Tvmult_add<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<SparseMatrix<double>, float>::Tvmult_add<double> (
  Vector<double> &, const Vector<double> &) const;

template class PreconditionBlockSOR<SparseMatrix<double>, double>;

template void PreconditionBlockSOR<SparseMatrix<double>, double>::vmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<SparseMatrix<double>, double>::vmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSOR<SparseMatrix<double>, double>::Tvmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<SparseMatrix<double>, double>::Tvmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSOR<SparseMatrix<double>, double>::vmult_add<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<SparseMatrix<double>, double>::vmult_add<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSOR<SparseMatrix<double>, double>::Tvmult_add<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<SparseMatrix<double>, double>::Tvmult_add<double> (
  Vector<double> &, const Vector<double> &) const;


/*--------------------- PreconditionBlockSSOR -----------------------*/


// explicit instantiations for "float" PreconditionBlock
template class PreconditionBlockSSOR<SparseMatrix<float>, float>;

template void PreconditionBlockSSOR<SparseMatrix<float>, float>::vmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSSOR<SparseMatrix<float>, float>::vmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSSOR<SparseMatrix<float>, float>::Tvmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSSOR<SparseMatrix<float>, float>::Tvmult<double> (
  Vector<double> &, const Vector<double> &) const;


// the instantiation for class PreconditionBlockSSOR<SparseMatrix<float>, double> is skipped
// because it does not make sense to have inverse block matrices with
// higher precision than the matrix itself


// explicit instantiations for "double" PreconditionBlockSSOR
template class PreconditionBlockSSOR<SparseMatrix<double>, float>;


template void PreconditionBlockSSOR<SparseMatrix<double>, float>::vmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSSOR<SparseMatrix<double>, float>::vmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSSOR<SparseMatrix<double>, float>::Tvmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSSOR<SparseMatrix<double>, float>::Tvmult<double> (
  Vector<double> &, const Vector<double> &) const;

template class PreconditionBlockSSOR<SparseMatrix<double>, double>;

template void PreconditionBlockSSOR<SparseMatrix<double>, double>::vmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSSOR<SparseMatrix<double>, double>::vmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSSOR<SparseMatrix<double>, double>::Tvmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSSOR<SparseMatrix<double>, double>::Tvmult<double> (
  Vector<double> &, const Vector<double> &) const;

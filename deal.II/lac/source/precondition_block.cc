//----------------------------  precondition_block.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  precondition_block.cc  ---------------------------


#include <lac/precondition_block.templates.h>


// explicit instantiations for "float" PreconditionBlock
template class PreconditionBlock<float, float>;

// the instantiation for class PreconditionBlock<float, double> is skipped
// because it does not make sense to have inverse block matrices with
// higher precision than the matrix itself


// explicit instantiations for "double" PreconditionBlock
template class PreconditionBlock<double, float>;

template class PreconditionBlock<double, double>;


/*--------------------- PreconditionBlockJacobi -----------------------*/


// explicit instantiations for "float" PreconditionBlock
template class PreconditionBlockJacobi<float, float>;

template void PreconditionBlockJacobi<float, float>::vmult<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<float, float>::vmult<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<float, float>::Tvmult<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<float, float>::Tvmult<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<float, float>::vmult_add<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<float, float>::vmult_add<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<float, float>::Tvmult_add<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<float, float>::Tvmult_add<double>
(Vector<double> &, const Vector<double> &) const;

template class PreconditionBlockJacobi<double, float>;

template void PreconditionBlockJacobi<double, float>::vmult<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<double, float>::vmult<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<double, float>::Tvmult<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<double, float>::Tvmult<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<double, float>::vmult_add<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<double, float>::vmult_add<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<double, float>::Tvmult_add<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<double, float>::Tvmult_add<double>
(Vector<double> &, const Vector<double> &) const;

template class PreconditionBlockJacobi<double, double>;

template void PreconditionBlockJacobi<double, double>::vmult<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<double, double>::vmult<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<double, double>::Tvmult<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<double, double>::Tvmult<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<double, double>::vmult_add<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<double, double>::vmult_add<double>
(Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockJacobi<double, double>::Tvmult_add<float>
(Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockJacobi<double, double>::Tvmult_add<double>
(Vector<double> &, const Vector<double> &) const;

/*--------------------- PreconditionBlockGaussSeidel -----------------------*/


// explicit instantiations for "float" PreconditionBlock
template class PreconditionBlockSOR<float, float>;

template void PreconditionBlockSOR<float, float>::vmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<float, float>::vmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSOR<float, float>::Tvmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<float, float>::Tvmult<double> (
  Vector<double> &, const Vector<double> &) const;


// the instantiation for class PreconditionBlockSOR<float, double> is skipped
// because it does not make sense to have inverse block matrices with
// higher precision than the matrix itself


// explicit instantiations for "double" PreconditionBlockSOR
template class PreconditionBlockSOR<double, float>;


template void PreconditionBlockSOR<double, float>::vmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<double, float>::vmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSOR<double, float>::Tvmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<double, float>::Tvmult<double> (
  Vector<double> &, const Vector<double> &) const;

template class PreconditionBlockSOR<double, double>;

template void PreconditionBlockSOR<double, double>::vmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<double, double>::vmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSOR<double, double>::Tvmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<double, double>::Tvmult<double> (
  Vector<double> &, const Vector<double> &) const;


/*--------------------- PreconditionBlockSSOR -----------------------*/


// explicit instantiations for "float" PreconditionBlock
template class PreconditionBlockSSOR<float, float>;

template void PreconditionBlockSSOR<float, float>::vmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSSOR<float, float>::vmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSSOR<float, float>::Tvmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSSOR<float, float>::Tvmult<double> (
  Vector<double> &, const Vector<double> &) const;


// the instantiation for class PreconditionBlockSSOR<float, double> is skipped
// because it does not make sense to have inverse block matrices with
// higher precision than the matrix itself


// explicit instantiations for "double" PreconditionBlockSSOR
template class PreconditionBlockSSOR<double, float>;


template void PreconditionBlockSSOR<double, float>::vmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSSOR<double, float>::vmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSSOR<double, float>::Tvmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSSOR<double, float>::Tvmult<double> (
  Vector<double> &, const Vector<double> &) const;

template class PreconditionBlockSSOR<double, double>;

template void PreconditionBlockSSOR<double, double>::vmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSSOR<double, double>::vmult<double> (
  Vector<double> &, const Vector<double> &) const;
template void PreconditionBlockSSOR<double, double>::Tvmult<float> (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSSOR<double, double>::Tvmult<double> (
  Vector<double> &, const Vector<double> &) const;

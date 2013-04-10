//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <deal.II/matrix_free/matrix_free.templates.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>

#include <iostream>

DEAL_II_NAMESPACE_OPEN

#include "matrix_free.inst"

template struct internal::MatrixFreeFunctions::ShapeInfo<double>;
template struct internal::MatrixFreeFunctions::ShapeInfo<float>;

DEAL_II_NAMESPACE_CLOSE

// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

#include <deal.II/base/config.h>

// This header undefines all macros that are not namespaced (i.e., macros that
// do not start with deal or DEAL).

#ifdef Assert
#  undef Assert
#endif // #ifdef Assert

#ifdef AssertARKode
#  undef AssertARKode
#endif // #ifdef AssertARKode

#ifdef AssertCuda
#  undef AssertCuda
#endif // #ifdef AssertCuda

#ifdef AssertCudaKernel
#  undef AssertCudaKernel
#endif // #ifdef AssertCudaKernel

#ifdef AssertCusolver
#  undef AssertCusolver
#endif // #ifdef AssertCusolver

#ifdef AssertCusparse
#  undef AssertCusparse
#endif // #ifdef AssertCusparse

#ifdef AssertDimension
#  undef AssertDimension
#endif // #ifdef AssertDimension

#ifdef AssertIDA
#  undef AssertIDA
#endif // #ifdef AssertIDA

#ifdef AssertIndexRange
#  undef AssertIndexRange
#endif // #ifdef AssertIndexRange

#ifdef AssertIsFinite
#  undef AssertIsFinite
#endif // #ifdef AssertIsFinite

#ifdef AssertKINSOL
#  undef AssertKINSOL
#endif // #ifdef AssertKINSOL

#ifdef AssertNothrow
#  undef AssertNothrow
#endif // #ifdef AssertNothrow

#ifdef AssertNothrowCuda
#  undef AssertNothrowCuda
#endif // #ifdef AssertNothrowCuda

#ifdef AssertNothrowCusparse
#  undef AssertNothrowCusparse
#endif // #ifdef AssertNothrowCusparse

#ifdef AssertThrow
#  undef AssertThrow
#endif // #ifdef AssertThrow

#ifdef AssertThrowMPI
#  undef AssertThrowMPI
#endif // #ifdef AssertThrowMPI

#ifdef AssertThrowExodusII
#  undef AssertThrowExodusII
#endif // #ifdef AssertThrowExodusII

#ifdef AssertVectorVectorDimension
#  undef AssertVectorVectorDimension
#endif // #ifdef AssertVectorVectorDimension

#ifdef DeclException0
#  undef DeclException0
#endif // #ifdef DeclException0

#ifdef DeclException1
#  undef DeclException1
#endif // #ifdef DeclException1

#ifdef DeclException2
#  undef DeclException2
#endif // #ifdef DeclException2

#ifdef DeclException3
#  undef DeclException3
#endif // #ifdef DeclException3

#ifdef DeclException4
#  undef DeclException4
#endif // #ifdef DeclException4

#ifdef DeclException5
#  undef DeclException5
#endif // #ifdef DeclException5

#ifdef DeclExceptionMsg
#  undef DeclExceptionMsg
#endif // #ifdef DeclExceptionMsg

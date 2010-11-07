//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/logstream.h>
#include <base/tensor.h>
#include <base/tensor_base.h>
#include <lac/full_matrix.templates.h>

DEAL_II_NAMESPACE_OPEN

#include "full_matrix.inst"

// do a few functions that currently don't fit the scheme because they have
// two template arguments that need to be different (the case of same
// arguments is covered by the default copy constructor and copy operator that
// is declared separately)

#define TEMPL_OP_EQ(S1,S2)			      \
  template FullMatrix<S1>& FullMatrix<S1>::operator = \
  (const FullMatrix<S2>&)

TEMPL_OP_EQ(double,float);
TEMPL_OP_EQ(float,double);

TEMPL_OP_EQ(long double,double);
TEMPL_OP_EQ(double,long double);

TEMPL_OP_EQ(long double,float);
TEMPL_OP_EQ(float,long double);


TEMPL_OP_EQ(std::complex<double>,std::complex<float>);
TEMPL_OP_EQ(std::complex<float>,std::complex<double>);

TEMPL_OP_EQ(std::complex<long double>,std::complex<double>);
TEMPL_OP_EQ(std::complex<double>,std::complex<long double>);

TEMPL_OP_EQ(std::complex<long double>,std::complex<float>);
TEMPL_OP_EQ(std::complex<float>,std::complex<long double>);

#undef TEMPL_OP_EQ

DEAL_II_NAMESPACE_CLOSE

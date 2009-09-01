//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/vector.templates.h>

DEAL_II_NAMESPACE_OPEN

#include "vector.inst"

// do a few functions that currently don't fit the scheme because they have
// two template arguments that need to be different (the case of same
// arguments is covered by the default copy constructor and copy operator that
// is declared separately)

#define TEMPL_COPY_CONSTRUCTOR(S1,S2)			\
  template Vector<S1>::Vector (const Vector<S2> &)

#ifndef DEAL_II_EXPLICIT_CONSTRUCTOR_BUG
TEMPL_COPY_CONSTRUCTOR(double,float);
TEMPL_COPY_CONSTRUCTOR(float,double);

TEMPL_COPY_CONSTRUCTOR(long double,double);
TEMPL_COPY_CONSTRUCTOR(double,long double);

TEMPL_COPY_CONSTRUCTOR(long double,float);
TEMPL_COPY_CONSTRUCTOR(float,long double);


TEMPL_COPY_CONSTRUCTOR(std::complex<double>,std::complex<float>);
TEMPL_COPY_CONSTRUCTOR(std::complex<float>,std::complex<double>);

TEMPL_COPY_CONSTRUCTOR(std::complex<long double>,std::complex<double>);
TEMPL_COPY_CONSTRUCTOR(std::complex<double>,std::complex<long double>);

TEMPL_COPY_CONSTRUCTOR(std::complex<long double>,std::complex<float>);
TEMPL_COPY_CONSTRUCTOR(std::complex<float>,std::complex<long double>);

#endif

#undef TEMPL_COPY_CONSTRUCTOR


#define TEMPL_OP_EQ(S1,S2) \
  template \
  Vector<S1>& \
  Vector<S1>::DEAL_II_MEMBER_OP_TEMPLATE_INST \
  operator=<S2>(const Vector<S2>&);		      \
  template void Vector<S1>::scale (const Vector<S2>&);	\
  template void Vector<S1>::equ (const S1, const Vector<S2>&)

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

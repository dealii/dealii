//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <base/qprojector.h>
#include <base/tensor_product_polynomials.h>
#include <base/polynomials_p.h>
#include <fe/fe_poly.h>
#include <fe/fe_values.h>

#include <fe/fe_poly.templates.h>

DEAL_II_NAMESPACE_OPEN

template class FE_Poly<TensorProductPolynomials<deal_II_dimension>, deal_II_dimension>;
template class FE_Poly<PolynomialSpace<deal_II_dimension>, deal_II_dimension>;
template class FE_Poly<PolynomialsP<deal_II_dimension>, deal_II_dimension>;

#if deal_II_dimension != 3

template class FE_Poly<TensorProductPolynomials<deal_II_dimension>, deal_II_dimension, deal_II_dimension+1>;
template class FE_Poly<PolynomialSpace<deal_II_dimension>, deal_II_dimension, deal_II_dimension+1>;
//template class FE_Poly<PolynomialsP<deal_II_dimension>, deal_II_dimension, deal_II_dimension+1>;

#endif
DEAL_II_NAMESPACE_CLOSE

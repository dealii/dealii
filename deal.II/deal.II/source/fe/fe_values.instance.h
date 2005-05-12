//-----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------

// Instantiations of functions in FEValuesBase

// Definition of vector type macros is in fe_values.cc


template
void FEValuesBase<deal_II_dimension>::get_function_values<IN>
(const IN&, std::vector<double>&) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values<IN>
(const IN&, std::vector<float>&) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values<IN>
(const IN&, const VectorSlice<const std::vector<unsigned int> >&, std::vector<double>&) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values<IN>
(const IN&, const VectorSlice<const std::vector<unsigned int> >&, std::vector<float>&) const;

template
void FEValuesBase<deal_II_dimension>::get_function_values<IN>
(const IN&, std::vector<Vector<double> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values<IN>
(const IN&, std::vector<Vector<float> > &) const;

template
void FEValuesBase<deal_II_dimension>::get_function_values<IN>
(const IN&, const VectorSlice<const std::vector<unsigned int> >&,
 std::vector<Vector<double> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values<IN>
(const IN&, const VectorSlice<const std::vector<unsigned int> >&,
 std::vector<Vector<float> > &) const;

template
void FEValuesBase<deal_II_dimension>::get_function_values<IN>
(const IN&, const VectorSlice<const std::vector<unsigned int> >&,
 std::vector<std::vector<double> > &, bool) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values<IN>
(const IN&, const VectorSlice<const std::vector<unsigned int> >&,
 std::vector<std::vector<float> > &, bool) const;

template
void FEValuesBase<deal_II_dimension>::get_function_grads<IN>
(const IN&, std::vector<Tensor<1,deal_II_dimension> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_grads<IN>
(const IN&, const VectorSlice<const std::vector<unsigned int> >&,
 std::vector<Tensor<1,deal_II_dimension> > &) const;

template
void FEValuesBase<deal_II_dimension>::get_function_grads<IN>
(const IN&, std::vector<std::vector<Tensor<1,deal_II_dimension> > > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_grads<IN>
(const IN&, const VectorSlice<const std::vector<unsigned int> >&,
 std::vector<std::vector<Tensor<1,deal_II_dimension> > > &, bool) const;

template
void FEValuesBase<deal_II_dimension>::get_function_2nd_derivatives<IN>
(const IN&, std::vector<Tensor<2,deal_II_dimension> > &) const;

template
void FEValuesBase<deal_II_dimension>::get_function_2nd_derivatives<IN>
(const IN&, std::vector<std::vector<Tensor<2,deal_II_dimension> > > &, bool) const;

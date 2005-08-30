//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

// Instantiations of functions in vectors.cc

// vectors.cc defines the vector type VEC

template
void VectorTools::interpolate<deal_II_dimension>
(const Mapping<deal_II_dimension>&,
 const DoFHandler<deal_II_dimension>&,
 const Function<deal_II_dimension>&,
 VEC&);
template
void VectorTools::interpolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension>&,
 const Function<deal_II_dimension>&,
 VEC&);

// Should these be instantiated for every combination of two types?
template
void VectorTools::interpolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension>&,
 const DoFHandler<deal_II_dimension>&,
 const FullMatrix<double>&,
 const VEC&,
 VEC&);


template
void VectorTools::integrate_difference<deal_II_dimension>
(const DoFHandler<deal_II_dimension>&,
 const VEC&,
 const Function<deal_II_dimension>&,
 Vector<double>&,
 const Quadrature<deal_II_dimension>&,
 const NormType&,
 const Function<deal_II_dimension>*,
 const double);
template
void VectorTools::integrate_difference<deal_II_dimension>
(const DoFHandler<deal_II_dimension>&,
 const VEC&,
 const Function<deal_II_dimension>&,
 Vector<float>&,
 const Quadrature<deal_II_dimension>&,
 const NormType&,
 const Function<deal_II_dimension>*,
 const double);
template
void VectorTools::integrate_difference<deal_II_dimension>
(const Mapping<deal_II_dimension>&,
 const DoFHandler<deal_II_dimension>&,
 const VEC&,
 const Function<deal_II_dimension>&,
 Vector<double>&,
 const Quadrature<deal_II_dimension>&,
 const NormType&,
 const Function<deal_II_dimension>*,
 const double);
template
void VectorTools::integrate_difference<deal_II_dimension>
(const Mapping<deal_II_dimension>&,
 const DoFHandler<deal_II_dimension>&,
 const VEC&,
 const Function<deal_II_dimension>&,
 Vector<float>&,
 const Quadrature<deal_II_dimension>&,
 const NormType&,
 const Function<deal_II_dimension>*,
 const double);

template
void VectorTools::point_difference<deal_II_dimension> (
  const DoFHandler<deal_II_dimension>&,
  const VEC&,
  const Function<deal_II_dimension>&,
  Vector<double>&,
  const Point<deal_II_dimension>&);

template
double VectorTools::compute_mean_value<deal_II_dimension>
(const Mapping<deal_II_dimension>&,
 const DoFHandler<deal_II_dimension>&,
 const Quadrature<deal_II_dimension>&,
 const VEC&,
 const unsigned int);
template
double VectorTools::compute_mean_value<deal_II_dimension>
(const DoFHandler<deal_II_dimension>&,
 const Quadrature<deal_II_dimension>&,
 const VEC&,
 const unsigned int);

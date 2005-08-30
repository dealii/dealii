//---------------------------------------------------------------------------
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
//---------------------------------------------------------------------------

#include<numerics/vectors.templates.h>

// explicit instantiations

#define VEC Vector<double>
#include "vectors.instance.h"
#undef VEC

#define VEC Vector<float>
#include "vectors.instance.h"
#undef VEC

#define VEC BlockVector<double>
#include "vectors.instance.h"
#undef VEC

#define VEC BlockVector<float>
#include "vectors.instance.h"
#undef VEC

template
void VectorTools::project<deal_II_dimension>
(const DoFHandler<deal_II_dimension>   &,
 const ConstraintMatrix                &,
 const Quadrature<deal_II_dimension>   &,
 const Function<deal_II_dimension>     &,
 Vector<double>                        &,
 const bool,
 const Quadrature<deal_II_dimension-1> &,
 const bool);

template
void VectorTools::create_right_hand_side<deal_II_dimension>
(const Mapping<deal_II_dimension>    &,
 const DoFHandler<deal_II_dimension> &,
 const Quadrature<deal_II_dimension> &,
 const Function<deal_II_dimension>   &,
 Vector<double>                      &);
template
void VectorTools::create_right_hand_side<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &,
 const Quadrature<deal_II_dimension> &,
 const Function<deal_II_dimension>   &,
 Vector<double>                      &);

#if deal_II_dimension != 1
template
void
VectorTools::create_boundary_right_hand_side<deal_II_dimension>
(const Mapping<deal_II_dimension>    &,
 const DoFHandler<deal_II_dimension> &,
 const Quadrature<deal_II_dimension-1> &,
 const Function<deal_II_dimension>   &,
 Vector<double>                      &,
 const std::set<unsigned char> &);
#endif

template
void
VectorTools::create_boundary_right_hand_side<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &,
 const Quadrature<deal_II_dimension-1> &,
 const Function<deal_II_dimension>   &,
 Vector<double>                      &,
 const std::set<unsigned char> &);
template
void VectorTools::interpolate_boundary_values<deal_II_dimension> (
  const DoFHandler<deal_II_dimension> &,
  const unsigned char,
  const Function<deal_II_dimension>   &,
  std::map<unsigned int,double>       &,
  const std::vector<bool>    &);

#if deal_II_dimension != 1
template
void VectorTools::project_boundary_values<deal_II_dimension>
(const Mapping<deal_II_dimension>     &,
 const DoFHandler<deal_II_dimension>  &,
 const FunctionMap<deal_II_dimension>::type &,
 const Quadrature<deal_II_dimension-1>&,
 std::map<unsigned int,double>        &);
#endif

template
void VectorTools::project_boundary_values<deal_II_dimension>
(const DoFHandler<deal_II_dimension>  &,
 const FunctionMap<deal_II_dimension>::type &,
 const Quadrature<deal_II_dimension-1>&,
 std::map<unsigned int,double>        &);


// the following two functions are not derived from a template in 1d
// and thus need no explicit instantiation
#if deal_II_dimension > 1
template
void VectorTools::interpolate_boundary_values<deal_II_dimension>
(const Mapping<deal_II_dimension>    &,
 const DoFHandler<deal_II_dimension> &,
 const unsigned char,
 const Function<deal_II_dimension>   &,
 std::map<unsigned int,double>       &,
 const std::vector<bool>    &);

template
void VectorTools::project<deal_II_dimension>
(const Mapping<deal_II_dimension>      &,
 const DoFHandler<deal_II_dimension>   &,
 const ConstraintMatrix                &,
 const Quadrature<deal_II_dimension>   &,
 const Function<deal_II_dimension>     &,
 Vector<double>                        &,
 const bool,
 const Quadrature<deal_II_dimension-1> &,
 const bool);
#endif


template
void
VectorTools::interpolate_boundary_values<deal_II_dimension>
(const DoFHandler<deal_II_dimension>         &,
 const FunctionMap<deal_II_dimension>::type &,
 std::map<unsigned int,double> &,
 const std::vector<bool>       &);


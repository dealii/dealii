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

#include <dofs/dof_handler.h>
#include <dofs/hp_dof_handler.h>

namespace 
{
  // Functions for the DoFHandler
  template <int dim>
  unsigned int
  max_dofs_per_cell (const DoFHandler<dim> &dh)
  {
    return dh.get_fe().dofs_per_cell;
  }


  template <int dim>
  unsigned int
  max_dofs_per_face (const DoFHandler<dim> &dh) 
  {
    return dh.get_fe().dofs_per_face;
  }


  template <int dim>
  unsigned int
  get_n_components (const DoFHandler<dim> &dh) 
  {
    return dh.get_fe().n_components();
  }


  // Functions for the hp::DoFHandler
  template <int dim>
  unsigned int
  max_dofs_per_cell (const hp::DoFHandler<dim> &dh) 
  {
    return dh.get_fe().max_dofs_per_cell ();
  }


  template <int dim>
  unsigned int
  max_dofs_per_face (const hp::DoFHandler<dim> &dh) 
  {
    return dh.get_fe().max_dofs_per_face ();
  }



  template <int dim>
  unsigned int
  get_n_components (const hp::DoFHandler<dim> &dh) 
  {
    return dh.get_fe().n_components();
  }
}


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



// Due to introducing the DoFHandler as a template parameter,
// the following instantiations are required in 1d
#if deal_II_dimension == 1
template
void VectorTools::interpolate_boundary_values<deal_II_dimension> 
(const Mapping<1>         &,
 const DoFHandler<1>      &,
 const unsigned char,
 const Function<1>        &,
 std::map<unsigned int,double> &,
 const std::vector<bool>       &);
#endif

// the following two functions are not derived from a template in 1d
// and thus need no explicit instantiation
#if deal_II_dimension > 1
template
void VectorTools::interpolate_boundary_values<deal_II_dimension>
(const Mapping<deal_II_dimension>    &,
 const DoFHandler<deal_II_dimension> &,
 const FunctionMap<deal_II_dimension>::type &,
 std::map<unsigned int,double>       &,
 const std::vector<bool>    &);

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


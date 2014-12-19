// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/numerics/data_postprocessor.h>

DEAL_II_NAMESPACE_OPEN



// -------------------------- DataPostprocessor ---------------------------

template <int dim>
DataPostprocessor<dim>::~DataPostprocessor()
{}



template <int dim>
void
DataPostprocessor<dim>::
compute_derived_quantities_scalar (const std::vector<double>         &/*uh*/,
                                   const std::vector<Tensor<1,dim> > &/*duh*/,
                                   const std::vector<Tensor<2,dim> > &/*dduh*/,
                                   const std::vector<Point<dim> >    &/*normals*/,
                                   std::vector<Vector<double> >      &computed_quantities) const
{
  computed_quantities.clear();
  AssertThrow(false,ExcPureFunctionCalled());
}


template <int dim>
void
DataPostprocessor<dim>::
compute_derived_quantities_scalar (const std::vector<double>         &uh,
                                   const std::vector<Tensor<1,dim> > &duh,
                                   const std::vector<Tensor<2,dim> > &dduh,
                                   const std::vector<Point<dim> >    &normals,
                                   const std::vector<Point<dim> >    &/*evaluation_points*/,
                                   std::vector<Vector<double> >      &computed_quantities) const
{
  compute_derived_quantities_scalar(uh, duh, dduh, normals, computed_quantities);
}



template <int dim>
void
DataPostprocessor<dim>::
compute_derived_quantities_vector (const std::vector<Vector<double> > &/*uh*/,
                                   const std::vector<std::vector<Tensor<1,dim> > > &/*duh*/,
                                   const std::vector<std::vector<Tensor<2,dim> > > &/*dduh*/,
                                   const std::vector<Point<dim> >                  &/*normals*/,
                                   std::vector<Vector<double> >                    &computed_quantities) const
{
  computed_quantities.clear();
  AssertThrow(false,ExcPureFunctionCalled());
}



template <int dim>
void
DataPostprocessor<dim>::
compute_derived_quantities_vector (const std::vector<Vector<double> > &uh,
                                   const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                   const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                   const std::vector<Point<dim> >                  &normals,
                                   const std::vector<Point<dim> >                  &/*evaluation_points*/,
                                   std::vector<Vector<double> >                    &computed_quantities) const
{
  compute_derived_quantities_vector(uh, duh, dduh, normals, computed_quantities);
}



template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
DataPostprocessor<dim>::get_data_component_interpretation () const
{
  // default implementation assumes that all
  // components are independent scalars
  return
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    (get_names().size(),
     DataComponentInterpretation::component_is_scalar);
}


// -------------------------- DataPostprocessorScalar ---------------------------

template <int dim>
DataPostprocessorScalar<dim>::
DataPostprocessorScalar (const std::string &name,
                         const UpdateFlags  update_flags)
  :
  name (name),
  update_flags (update_flags)
{}



template <int dim>
std::vector<std::string>
DataPostprocessorScalar<dim>::
get_names () const
{
  return std::vector<std::string> (1, name);
}



template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
DataPostprocessorScalar<dim>::
get_data_component_interpretation () const
{
  return
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    (1, DataComponentInterpretation::component_is_scalar);
}


template <int dim>
UpdateFlags
DataPostprocessorScalar<dim>::
get_needed_update_flags () const
{
  return update_flags;
}



// -------------------------- DataPostprocessorVector ---------------------------

template <int dim>
DataPostprocessorVector<dim>::
DataPostprocessorVector (const std::string &name,
                         const UpdateFlags  update_flags)
  :
  name (name),
  update_flags (update_flags)
{}



template <int dim>
std::vector<std::string>
DataPostprocessorVector<dim>::
get_names () const
{
  return std::vector<std::string> (dim, name);
}



template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
DataPostprocessorVector<dim>::
get_data_component_interpretation () const
{
  return
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    (dim, DataComponentInterpretation::component_is_part_of_vector);
}


template <int dim>
UpdateFlags
DataPostprocessorVector<dim>::
get_needed_update_flags () const
{
  return update_flags;
}


// explicit instantiation
#include "data_postprocessor.inst"


DEAL_II_NAMESPACE_CLOSE

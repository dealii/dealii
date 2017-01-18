// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2015 by the deal.II authors
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
evaluate_scalar_field (const DataPostprocessorInputs::Scalar<dim> &inputs,
                       std::vector<Vector<double> >               &computed_quantities) const
{
  // for backward compatibility, call the old function.
  // this also requires converting the accidental use
  // of Point<dim> for normal vectors
  std::vector<Point<dim> > normals (inputs.normals.begin(),
                                    inputs.normals.end());
  compute_derived_quantities_scalar(inputs.solution_values,
                                    inputs.solution_gradients,
                                    inputs.solution_hessians,
                                    normals,
                                    inputs.evaluation_points,
                                    computed_quantities);
}



template <int dim>
void
DataPostprocessor<dim>::
compute_derived_quantities_scalar (const std::vector<double>         &/*solution_values*/,
                                   const std::vector<Tensor<1,dim> > &/*solution_gradients*/,
                                   const std::vector<Tensor<2,dim> > &/*solution_hessians*/,
                                   const std::vector<Point<dim> >    &/*normals*/,
                                   const std::vector<Point<dim> >    &/*evaluation_points*/,
                                   std::vector<Vector<double> >      &computed_quantities) const
{
  computed_quantities.clear();
  AssertThrow(false,ExcPureFunctionCalled());
}



template <int dim>
void
DataPostprocessor<dim>::
evaluate_vector_field (const DataPostprocessorInputs::Vector<dim> &inputs,
                       std::vector<Vector<double> >               &computed_quantities) const
{
  // for backward compatibility, call the old function.
  // this also requires converting the accidental use
  // of Point<dim> for normal vectors
  std::vector<Point<dim> > normals (inputs.normals.begin(),
                                    inputs.normals.end());
  compute_derived_quantities_vector(inputs.solution_values,
                                    inputs.solution_gradients,
                                    inputs.solution_hessians,
                                    normals,
                                    inputs.evaluation_points,
                                    computed_quantities);
}



template <int dim>
void
DataPostprocessor<dim>::
compute_derived_quantities_vector (const std::vector<Vector<double> > &/*solution_values*/,
                                   const std::vector<std::vector<Tensor<1,dim> > > &/*solution_gradients*/,
                                   const std::vector<std::vector<Tensor<2,dim> > > &/*solution_hessians*/,
                                   const std::vector<Point<dim> >                  &/*normals*/,
                                   const std::vector<Point<dim> >                  &/*evaluation_points*/,
                                   std::vector<Vector<double> >                    &computed_quantities) const
{
  computed_quantities.clear();
  AssertThrow(false,ExcPureFunctionCalled());
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

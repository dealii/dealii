//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2007, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#include <numerics/data_postprocessor.h>

DEAL_II_NAMESPACE_OPEN


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
    (n_output_variables(),
     DataComponentInterpretation::component_is_scalar);
}



// explicit instantiation
#include "data_postprocessor.inst"


DEAL_II_NAMESPACE_CLOSE

//----------------------------  fe_values.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_values.cc  ---------------------------


#include <fe/fe.h>
#include <fe/mapping_q1.h>
#include <fe/fe_values.h>
#include <base/memory_consumption.h>
#include <base/quadrature.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_boundary.h>
#include <dofs/dof_accessor.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/petsc_vector.h>

#include <iomanip>


template <int dim>
void
FEValuesData<dim>::initialize (const unsigned int        n_quadrature_points,
			       const FiniteElement<dim> &fe,
			       const UpdateFlags         flags)
{
  this->update_flags = flags;

				   // initialize the table mapping
				   // from shape function number to
				   // the rows in the tables denoting
				   // its first non-zero
				   // component. with this also count
				   // the total number of non-zero
				   // components accumulated over all
				   // shape functions
  this->shape_function_to_row_table.resize (fe.dofs_per_cell);
  unsigned int row = 0;
  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
    {
      this->shape_function_to_row_table[i] = row;
      row += fe.n_nonzero_components (i);
    };
	 
  const unsigned int n_nonzero_shape_components = row;
  Assert (n_nonzero_shape_components >= fe.dofs_per_cell,
	  ExcInternalError());

				   // with the number of rows now
				   // known, initialize those fields
				   // that we will need to their
				   // correct size
  if (flags & update_values)
    this->shape_values.reinit(n_nonzero_shape_components,
			      n_quadrature_points);

  if (flags & update_gradients)
    this->shape_gradients.reinit(n_nonzero_shape_components,
				 n_quadrature_points);

  if (flags & update_second_derivatives)
    this->shape_2nd_derivatives.reinit(n_nonzero_shape_components,
				       n_quadrature_points);
  
  if (flags & update_q_points)
    this->quadrature_points.resize(n_quadrature_points);

  if (flags & update_JxW_values)
    this->JxW_values.resize(n_quadrature_points);

  if (flags & update_boundary_forms)
    this->boundary_forms.resize(n_quadrature_points);

  if (flags & update_normal_vectors)
    this->normal_vectors.resize(n_quadrature_points);
}


/*------------------------------- FEValuesBase ---------------------------*/


template <int dim>
FEValuesBase<dim>::FEValuesBase (const unsigned int n_q_points,
				 const unsigned int dofs_per_cell,
				 const UpdateFlags flags,
				 const Mapping<dim>       &mapping,
				 const FiniteElement<dim> &fe)
		:
		n_quadrature_points (n_q_points),
		dofs_per_cell (dofs_per_cell),
		mapping(&mapping),
		fe(&fe),
		mapping_data(0),
		fe_data(0)
{
  this->update_flags = flags;
}



template <int dim>
FEValuesBase<dim>::~FEValuesBase ()
{
				   // delete those fields that were
				   // created by the mapping and
				   // finite element objects,
				   // respectively, but of which we
				   // have assumed ownership
  if (fe_data != 0)
    {
      typename Mapping<dim>::InternalDataBase *tmp1=fe_data;
      fe_data=0;
      delete tmp1;
    }

  if (mapping_data != 0)
    {
      typename Mapping<dim>::InternalDataBase *tmp1=mapping_data;
      mapping_data=0;
      delete tmp1;
    }
}



template <int dim>
template <class InputVector, typename number>
void FEValuesBase<dim>::get_function_values (const InputVector            &fe_function,
					     std::vector<number> &values) const
{
  Assert (this->update_flags & update_values, ExcAccessToUninitializedField());
  Assert (fe->n_components() == 1,
	  ExcWrongNoOfComponents());
  Assert (values.size() == n_quadrature_points,
	  ExcWrongVectorSize(values.size(), n_quadrature_points));
  Assert (fe_function.size() == present_cell->get_dof_handler().n_dofs(),
	  ExcWrongVectorSize(fe_function.size(), present_cell->get_dof_handler().n_dofs()));

				   // get function values of dofs
				   // on this cell
  Vector<typename InputVector::value_type> dof_values (dofs_per_cell);
  if (present_cell->active())
    present_cell->get_dof_values (fe_function, dof_values);
  else
    present_cell->get_interpolated_dof_values(fe_function, dof_values);

				   // initialize with zero
  std::fill_n (values.begin(), n_quadrature_points, 0);

				   // add up contributions of trial
				   // functions. note that here we
				   // deal with scalar finite
				   // elements, so no need to check
				   // for non-primitivity of shape
				   // functions
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
      values[point] += (dof_values(shape_func) *
			this->shape_value(shape_func, point));
}



template <int dim>
template <class InputVector, typename number>
void FEValuesBase<dim>::get_function_values (const InputVector                     &fe_function,
					     std::vector<Vector<number> > &values) const
{
  Assert (n_quadrature_points == values.size(),
	  ExcWrongVectorSize(values.size(), n_quadrature_points));

  const unsigned int n_components = fe->n_components();
  for (unsigned i=0;i<values.size();++i)
    Assert (values[i].size() == n_components, ExcWrongNoOfComponents());

  Assert (this->update_flags & update_values, ExcAccessToUninitializedField());
  Assert (fe_function.size() == present_cell->get_dof_handler().n_dofs(),
	  ExcWrongVectorSize(fe_function.size(), present_cell->get_dof_handler().n_dofs()));
    
				   // get function values of dofs
				   // on this cell
  Vector<typename InputVector::value_type> dof_values (dofs_per_cell);
  if (present_cell->active())
    present_cell->get_dof_values (fe_function, dof_values);
  else
    present_cell->get_interpolated_dof_values(fe_function, dof_values);
  
				   // initialize with zero
  for (unsigned i=0;i<values.size();++i)
    std::fill_n (values[i].begin(), values[i].size(), 0);

				   // add up contributions of trial
				   // functions. now check whether the
				   // shape function is primitive or
				   // not. if it is, then set its only
				   // non-zero component, otherwise
				   // loop over components
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
      if (fe->is_primitive(shape_func))
	values[point](fe->system_to_component_index(shape_func).first)
	  += (dof_values(shape_func) * shape_value(shape_func, point));
      else
	for (unsigned int c=0; c<n_components; ++c)
	  values[point](c) += (dof_values(shape_func) *
			       shape_value_component(shape_func, point, c));
}



template <int dim>
const std::vector<Point<dim> > &
FEValuesBase<dim>::get_quadrature_points () const
{
  Assert (this->update_flags & update_q_points, ExcAccessToUninitializedField());
  return this->quadrature_points;
}



template <int dim>
const std::vector<double> &
FEValuesBase<dim>::get_JxW_values () const
{
  Assert (this->update_flags & update_JxW_values, ExcAccessToUninitializedField());
  return this->JxW_values;
}



template <int dim>
template <class InputVector>
void
FEValuesBase<dim>::
get_function_grads (const InputVector           &fe_function,
		    std::vector<Tensor<1,dim> > &gradients) const
{
  Assert (this->update_flags & update_gradients, ExcAccessToUninitializedField());

  Assert (fe->n_components() == 1,
	  ExcWrongNoOfComponents());
  Assert (gradients.size() == n_quadrature_points,
	  ExcWrongVectorSize(gradients.size(), n_quadrature_points));
  Assert (fe_function.size() == present_cell->get_dof_handler().n_dofs(),
	  ExcWrongVectorSize(fe_function.size(), present_cell->get_dof_handler().n_dofs()));

				   // get function values of dofs
				   // on this cell
  Vector<typename InputVector::value_type> dof_values (dofs_per_cell);
  if (present_cell->active())
    present_cell->get_dof_values (fe_function, dof_values);
  else
    present_cell->get_interpolated_dof_values(fe_function, dof_values);

				   // initialize with zero
  std::fill_n (gradients.begin(), n_quadrature_points, Tensor<1,dim>());

				   // add up contributions of trial
				   // functions. note that here we
				   // deal with scalar finite
				   // elements, so no need to check
				   // for non-primitivity of shape
				   // functions
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
      {
	Tensor<1,dim> tmp = this->shape_grad(shape_func,point);
	tmp *= dof_values(shape_func);
	gradients[point] += tmp;
      };
}



template <int dim>
template <class InputVector>
void
FEValuesBase<dim>::
get_function_grads (const InputVector                         &fe_function,
		    std::vector<std::vector<Tensor<1,dim> > > &gradients) const
{
  Assert (n_quadrature_points == gradients.size(),
	  ExcWrongNoOfComponents());

  const unsigned int n_components = fe->n_components();
  for (unsigned i=0; i<gradients.size(); ++i)
    Assert (gradients[i].size() == n_components, ExcWrongNoOfComponents());

  Assert (this->update_flags & update_gradients, ExcAccessToUninitializedField());
  Assert (fe_function.size() == present_cell->get_dof_handler().n_dofs(),
	  ExcWrongVectorSize(fe_function.size(), present_cell->get_dof_handler().n_dofs()));

				   // get function values of dofs
				   // on this cell
  Vector<typename InputVector::value_type> dof_values (dofs_per_cell);
  if (present_cell->active())
    present_cell->get_dof_values (fe_function, dof_values);
  else
    present_cell->get_interpolated_dof_values(fe_function, dof_values);

				   // initialize with zero
  for (unsigned i=0;i<gradients.size();++i)
    std::fill_n (gradients[i].begin(), gradients[i].size(), Tensor<1,dim>());

				   // add up contributions of trial
				   // functions. now check whether the
				   // shape function is primitive or
				   // not. if it is, then set its only
				   // non-zero component, otherwise
				   // loop over components
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
      if (fe->is_primitive (shape_func))
	{
	  Tensor<1,dim> tmp = this->shape_grad(shape_func,point);
	  tmp *= dof_values(shape_func);
	  gradients[point][fe->system_to_component_index(shape_func).first]
	    += tmp;
	}
      else
	for (unsigned int c=0; c<n_components; ++c)
	  {
	    Tensor<1,dim> tmp = this->shape_grad_component(shape_func,point,c);
	    tmp *= dof_values(shape_func);
	    gradients[point][c] += tmp;
	  };
}



template <int dim>
template <class InputVector>
void
FEValuesBase<dim>::
get_function_2nd_derivatives (const InputVector           &fe_function,
			      std::vector<Tensor<2,dim> > &second_derivatives) const
{
  Assert (fe->n_components() == 1,
	  ExcWrongNoOfComponents());
  Assert (second_derivatives.size() == n_quadrature_points,
	  ExcWrongVectorSize(second_derivatives.size(), n_quadrature_points));
  Assert (this->update_flags & update_second_derivatives, ExcAccessToUninitializedField());
  Assert (fe_function.size() == present_cell->get_dof_handler().n_dofs(),
	  ExcWrongVectorSize(fe_function.size(), present_cell->get_dof_handler().n_dofs()));

				   // get function values of dofs
				   // on this cell
  Vector<typename InputVector::value_type> dof_values (dofs_per_cell);
  if (present_cell->active())
    present_cell->get_dof_values (fe_function, dof_values);
  else
    present_cell->get_interpolated_dof_values(fe_function, dof_values);

				   // initialize with zero
  std::fill_n (second_derivatives.begin(), n_quadrature_points, Tensor<2,dim>());

				   // add up contributions of trial
				   // functions. note that here we
				   // deal with scalar finite
				   // elements, so no need to check
				   // for non-primitivity of shape
				   // functions
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
      {
	Tensor<2,dim> tmp = this->shape_2nd_derivative(shape_func,point);
	tmp *= dof_values(shape_func);
	second_derivatives[point] += tmp;
      };
}



template <int dim>
template <class InputVector>
void
FEValuesBase<dim>::
get_function_2nd_derivatives (const InputVector                         &fe_function,
			      std::vector<std::vector<Tensor<2,dim> > > &second_derivs) const
{
  Assert (n_quadrature_points == second_derivs.size(),
	  ExcWrongNoOfComponents());

  const unsigned int n_components = fe->n_components();
  for (unsigned i=0;i<second_derivs.size();++i)
    Assert (second_derivs[i].size() == n_components, ExcWrongNoOfComponents());

  Assert (this->update_flags & update_second_derivatives, ExcAccessToUninitializedField());
  Assert (fe_function.size() == present_cell->get_dof_handler().n_dofs(),
	  ExcWrongVectorSize(fe_function.size(), present_cell->get_dof_handler().n_dofs()));

				   // get function values of dofs
				   // on this cell
  Vector<typename InputVector::value_type> dof_values (dofs_per_cell);
  if (present_cell->active())
    present_cell->get_dof_values (fe_function, dof_values);
  else
    present_cell->get_interpolated_dof_values(fe_function, dof_values);

				   // initialize with zero
  for (unsigned i=0;i<second_derivs.size();++i)
    std::fill_n (second_derivs[i].begin(), second_derivs[i].size(), Tensor<2,dim>());

				   // add up contributions of trial
				   // functions
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
      if (fe->is_primitive(shape_func))
	{
	  Tensor<2,dim> tmp(shape_2nd_derivative(shape_func,point));
	  tmp *= dof_values(shape_func);
	  second_derivs[point][fe->system_to_component_index(shape_func).first]
	    += tmp;
	}
      else
	for (unsigned int c=0; c<n_components; ++c)
	  {
	    Tensor<2,dim> tmp = this->shape_2nd_derivative_component(shape_func,point,c);
	    tmp *= dof_values(shape_func);
	    second_derivs[point][c] += tmp;
	  };
}



template <int dim>
const Point<dim> &
FEValuesBase<dim>::quadrature_point (const unsigned int i) const
{
  Assert (this->update_flags & update_q_points, ExcAccessToUninitializedField());
  Assert (i<this->quadrature_points.size(), ExcIndexRange(i, 0, this->quadrature_points.size()));
  
  return this->quadrature_points[i];
}




template <int dim>
double FEValuesBase<dim>::JxW (const unsigned int i) const
{
  Assert (this->update_flags & update_JxW_values, ExcAccessToUninitializedField());
  Assert (i<this->JxW_values.size(), ExcIndexRange(i, 0, this->JxW_values.size()));
  
  return this->JxW_values[i];
}



template <int dim>
unsigned int
FEValuesBase<dim>::memory_consumption () const
{
  return (MemoryConsumption::memory_consumption (this->shape_values) +
	  MemoryConsumption::memory_consumption (this->shape_gradients) +
	  MemoryConsumption::memory_consumption (this->shape_2nd_derivatives) +
	  MemoryConsumption::memory_consumption (this->JxW_values) +
	  MemoryConsumption::memory_consumption (this->quadrature_points) +
	  MemoryConsumption::memory_consumption (this->normal_vectors) +
	  MemoryConsumption::memory_consumption (this->boundary_forms) +
	  sizeof(this->update_flags) +
	  MemoryConsumption::memory_consumption (present_cell) +
	  MemoryConsumption::memory_consumption (n_quadrature_points) +
	  MemoryConsumption::memory_consumption (dofs_per_cell) +
	  MemoryConsumption::memory_consumption (mapping) +
	  MemoryConsumption::memory_consumption (fe) +
	  MemoryConsumption::memory_consumption (mapping_data) +
	  MemoryConsumption::memory_consumption (*mapping_data) +
	  MemoryConsumption::memory_consumption (fe_data) +
	  MemoryConsumption::memory_consumption (*fe_data) +
	  MemoryConsumption::memory_consumption (this->shape_function_to_row_table));
}



template <int dim>
UpdateFlags
FEValuesBase<dim>::compute_update_flags (const UpdateFlags update_flags) const
{
  
				   // first find out which objects
				   // need to be recomputed on each
				   // cell we visit. this we have to
				   // ask the finite element and mapping.
				   // elements are first since they
				   // might require update in mapping
  UpdateFlags flags = update_flags
		      | fe->update_once (update_flags)
		      | fe->update_each (update_flags);
  flags |= mapping->update_once (flags)
	   | mapping->update_each (flags);

  return flags;
}



template <int dim>
const Mapping<dim> &
FEValuesBase<dim>::get_default_mapping ()
{
  static const MappingQ1<dim> mapping_q1;
  return mapping_q1;
}


/*------------------------------- FEValues -------------------------------*/



template <int dim>
FEValues<dim>::FEValues (const Mapping<dim>       &mapping,
			 const FiniteElement<dim> &fe,
			 const Quadrature<dim>    &q,
			 const UpdateFlags         update_flags)
		:
		FEValuesBase<dim> (q.n_quadrature_points,
				   fe.dofs_per_cell,
				   update_default,
				   mapping,
				   fe),
  quadrature (q)
{
  initialize (update_flags);
}



template <int dim>
FEValues<dim>::FEValues (const FiniteElement<dim> &fe,
			 const Quadrature<dim>    &q,
			 const UpdateFlags         update_flags)
		:
		FEValuesBase<dim> (q.n_quadrature_points,
				   fe.dofs_per_cell,
				   update_default,
				   this->get_default_mapping(),
				   fe),
  quadrature (q)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  initialize (update_flags);
}



template <int dim>
void
FEValues<dim>::initialize (const UpdateFlags update_flags)
{
				   // you can't compute normal vectors
				   // on cells, only on faces
  Assert ((update_flags & update_normal_vectors) == false,
	  typename FEValuesBase<dim>::ExcInvalidUpdateFlag());

  const UpdateFlags flags = this->compute_update_flags (update_flags);
  
				   // then get objects into which the
				   // FE and the Mapping can store
				   // intermediate data used across
				   // calls to reinit
  this->mapping_data = this->mapping->get_data(flags, quadrature);
  this->fe_data      = this->fe->get_data(flags, *this->mapping, quadrature);

				   // set up objects within this class
  FEValuesData<dim>::initialize (this->n_quadrature_points, *this->fe, flags);
}



template <int dim>
void FEValues<dim>::reinit (const typename DoFHandler<dim>::cell_iterator &cell)
{
				   // assert that the finite elements
				   // passed to the constructor and
				   // used by the DoFHandler used by
				   // this cell, are the same
  Assert (static_cast<const FiniteElementData<dim>&>(*this->fe) ==
	  static_cast<const FiniteElementData<dim>&>(cell->get_dof_handler().get_fe()),
	  typename FEValuesBase<dim>::ExcFEDontMatch());

  this->present_cell = cell;

  this->get_mapping().fill_fe_values(cell,
				     quadrature,
				     *this->mapping_data,
				     this->quadrature_points,
				     this->JxW_values);
  
  this->get_fe().fill_fe_values(this->get_mapping(),
				cell,
				quadrature,
				*this->mapping_data,
				*this->fe_data,
				*this);

  this->fe_data->clear_first_cell ();
  this->mapping_data->clear_first_cell ();
}



template <int dim>
unsigned int
FEValues<dim>::memory_consumption () const
{
  return (FEValuesBase<dim>::memory_consumption () +
	  MemoryConsumption::memory_consumption (quadrature));
}


/*------------------------------- FEFaceValuesBase --------------------------*/


template <int dim>
FEFaceValuesBase<dim>::FEFaceValuesBase (const unsigned int n_q_points,
					 const unsigned int dofs_per_cell,
					 const UpdateFlags,
					 const Mapping<dim> &mapping,      
					 const FiniteElement<dim> &fe,
					 const Quadrature<dim-1>& quadrature)
		:
		FEValuesBase<dim> (n_q_points,
				   dofs_per_cell,
				   update_default,
				   mapping,
				   fe),
                quadrature(quadrature)
{}



template <int dim>
const std::vector<Point<dim> > &
FEFaceValuesBase<dim>::get_normal_vectors () const
{
  Assert (this->update_flags & update_normal_vectors,
	  typename FEValuesBase<dim>::ExcAccessToUninitializedField());
  return this->normal_vectors;
}



template <int dim>
const std::vector<Tensor<1,dim> > &
FEFaceValuesBase<dim>::get_boundary_forms () const
{
  Assert (this->update_flags & update_boundary_forms,
	  typename FEValuesBase<dim>::ExcAccessToUninitializedField());
  return this->boundary_forms;
}



template <int dim>
unsigned int
FEFaceValuesBase<dim>::memory_consumption () const
{
  return (FEValuesBase<dim>::memory_consumption () +
	  MemoryConsumption::memory_consumption (quadrature));
}


/*------------------------------- FEFaceValues -------------------------------*/


template <int dim>
FEFaceValues<dim>::FEFaceValues (const Mapping<dim>       &mapping,
				 const FiniteElement<dim> &fe,
				 const Quadrature<dim-1>  &quadrature,
				 const UpdateFlags         update_flags)
		:
		FEFaceValuesBase<dim> (quadrature.n_quadrature_points,
				       fe.dofs_per_cell,
				       update_flags,
				       mapping,
				       fe, quadrature)
{
  initialize (update_flags);
}



template <int dim>
FEFaceValues<dim>::FEFaceValues (const FiniteElement<dim> &fe,
				 const Quadrature<dim-1>  &quadrature,
				 const UpdateFlags         update_flags)
		:
		FEFaceValuesBase<dim> (quadrature.n_quadrature_points,
				       fe.dofs_per_cell,
				       update_flags,
				       this->get_default_mapping(),
				       fe, quadrature)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  initialize (update_flags);
}



template <int dim>
void
FEFaceValues<dim>::initialize (const UpdateFlags update_flags)
{
  const UpdateFlags flags = this->compute_update_flags (update_flags);
  
				   // then get objects into which the
				   // FE and the Mapping can store
				   // intermediate data used across
				   // calls to reinit
  this->mapping_data = this->mapping->get_face_data(flags, this->quadrature);
  this->fe_data      = this->fe->get_face_data(flags, *this->mapping, this->quadrature);

				   // set up objects within this class
  FEValuesData<dim>::initialize(this->n_quadrature_points, *this->fe, flags);
}



template <int dim>
void FEFaceValues<dim>::reinit (const typename DoFHandler<dim>::cell_iterator &cell,
				const unsigned int              face_no)
{
				   // assert that the finite elements
				   // passed to the constructor and
				   // used by the DoFHandler used by
				   // this cell, are the same
  Assert (static_cast<const FiniteElementData<dim>&>(*this->fe) ==
	  static_cast<const FiniteElementData<dim>&>(cell->get_dof_handler().get_fe()),
	  typename FEValuesBase<dim>::ExcFEDontMatch());

  this->present_cell = cell;
  
  this->get_mapping().fill_fe_face_values(cell, face_no,
					  this->quadrature,
					  *this->mapping_data,
					  this->quadrature_points,
					  this->JxW_values,
					  this->boundary_forms,
					  this->normal_vectors);
  
  this->get_fe().fill_fe_face_values(this->get_mapping(),
				     cell, face_no,
				     this->quadrature,
				     *this->mapping_data,
				     *this->fe_data,
				     *this);

  this->fe_data->clear_first_cell ();
  this->mapping_data->clear_first_cell ();
}


/*------------------------------- FESubFaceValues -------------------------------*/


template <int dim>
FESubfaceValues<dim>::FESubfaceValues (const Mapping<dim>       &mapping,
				       const FiniteElement<dim> &fe,
				       const Quadrature<dim-1>  &quadrature,
				       const UpdateFlags         update_flags)
		:
		FEFaceValuesBase<dim> (quadrature.n_quadrature_points,
				       fe.dofs_per_cell,
				       update_flags,
				       mapping,
				       fe, quadrature)
{
  initialize (update_flags);
}



template <int dim>
FESubfaceValues<dim>::FESubfaceValues (const FiniteElement<dim> &fe,
				       const Quadrature<dim-1>  &quadrature,
				       const UpdateFlags         update_flags)
		:
		FEFaceValuesBase<dim> (quadrature.n_quadrature_points,
				       fe.dofs_per_cell,
				       update_flags,
				       this->get_default_mapping(),
				       fe, quadrature)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  initialize (update_flags);
}



template <int dim>
void
FESubfaceValues<dim>::initialize (const UpdateFlags update_flags)
{
  const UpdateFlags flags = this->compute_update_flags (update_flags);
  
				   // then get objects into which the
				   // FE and the Mapping can store
				   // intermediate data used across
				   // calls to reinit
  this->mapping_data = this->mapping->get_subface_data(flags, this->quadrature);
  this->fe_data      = this->fe->get_subface_data(flags,
						  *this->mapping,
						  this->quadrature);

				   // set up objects within this class
  FEValuesData<dim>::initialize(this->n_quadrature_points, *this->fe, flags);
}



template <int dim>
void FESubfaceValues<dim>::reinit (const typename DoFHandler<dim>::cell_iterator &cell,
				   const unsigned int         face_no,
				   const unsigned int         subface_no)
{
				   // assert that the finite elements
				   // passed to the constructor and
				   // used by the DoFHandler used by
				   // this cell, are the same
  Assert (static_cast<const FiniteElementData<dim>&>(*this->fe) ==
	  static_cast<const FiniteElementData<dim>&>(cell->get_dof_handler().get_fe()),
	  typename FEValuesBase<dim>::ExcFEDontMatch());
  Assert (face_no < GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (face_no, 0, GeometryInfo<dim>::faces_per_cell));
  Assert (subface_no < GeometryInfo<dim>::subfaces_per_face,
	  ExcIndexRange (subface_no, 0, GeometryInfo<dim>::subfaces_per_face));
  Assert (cell->has_children() == false,
          ExcMessage ("You can't use subface data for cells that are "
                      "already refined. Iterate over their children "
                      "instead in these cases."));
  
  this->present_cell  = cell;

  this->get_mapping().fill_fe_subface_values(cell, face_no, subface_no,
					     this->quadrature,
					     *this->mapping_data,
					     this->quadrature_points,
					     this->JxW_values,
					     this->boundary_forms,
					     this->normal_vectors);
  
  this->get_fe().fill_fe_subface_values(this->get_mapping(),
					cell, face_no, subface_no,
					this->quadrature,
					*this->mapping_data,
					*this->fe_data,
					*this);

  this->fe_data->clear_first_cell ();
  this->mapping_data->clear_first_cell ();
}


/*------------------------------- Explicit Instantiations -------------*/

template class FEValuesData<deal_II_dimension>;
template class FEValuesBase<deal_II_dimension>;
template class FEValues<deal_II_dimension>;

#if deal_II_dimension >= 2
template class FEFaceValuesBase<deal_II_dimension>;
template class FEFaceValues<deal_II_dimension>;
template class FESubfaceValues<deal_II_dimension>;
#endif


//-----------------------------------------------------------------------------

template
void FEValuesBase<deal_II_dimension>::get_function_values<Vector<double> >
(const Vector<double>&,
 std::vector<double>&) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values<Vector<float> >
(const Vector<float>&,
 std::vector<double>&) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values<Vector<float> >
(const Vector<float>&,
 std::vector<float>&) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values<BlockVector<double> >
(const BlockVector<double>&,
 std::vector<double>&) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values<BlockVector<float> >
(const BlockVector<float>&,
 std::vector<double>&) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values<BlockVector<float> >
(const BlockVector<float>&,
 std::vector<float>&) const;

#ifdef DEAL_II_USE_PETSC
template
void FEValuesBase<deal_II_dimension>::get_function_values<PETScWrappers::Vector>
(const PETScWrappers::Vector &,
 std::vector<double>&) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values<PETScWrappers::Vector>
(const PETScWrappers::Vector &,
 std::vector<float>&) const;
#endif

//-----------------------------------------------------------------------------

template
void FEValuesBase<deal_II_dimension>::get_function_values<Vector<double> >
(const Vector<double> &,
 std::vector<Vector<double> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values<Vector<float> >
(const Vector<float> &,
 std::vector<Vector<double> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values<Vector<float> >
(const Vector<float> &,
 std::vector<Vector<float> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values<BlockVector<double> >
(const BlockVector<double> &,
 std::vector<Vector<double> >     &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values<BlockVector<float> >
(const BlockVector<float> &,
 std::vector<Vector<double> >     &) const;

#ifdef DEAL_II_USE_PETSC
template
void FEValuesBase<deal_II_dimension>::get_function_values<PETScWrappers::Vector>
(const PETScWrappers::Vector &,
 std::vector<Vector<double> >     &) const;
#endif

//-----------------------------------------------------------------------------

template
void FEValuesBase<deal_II_dimension>::get_function_grads<Vector<double> >
(const Vector<double> &,
 std::vector<Tensor<1,deal_II_dimension> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_grads<Vector<float> >
(const Vector<float> &,
 std::vector<Tensor<1,deal_II_dimension> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_grads<BlockVector<double> >
(const BlockVector<double> &,
 std::vector<Tensor<1,deal_II_dimension> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_grads<BlockVector<float> >
(const BlockVector<float> &,
 std::vector<Tensor<1,deal_II_dimension> > &) const;

#ifdef DEAL_II_USE_PETSC
template
void FEValuesBase<deal_II_dimension>::get_function_grads<PETScWrappers::Vector>
(const PETScWrappers::Vector &,
 std::vector<Tensor<1,deal_II_dimension> > &) const;
#endif

//-----------------------------------------------------------------------------

template
void FEValuesBase<deal_II_dimension>::get_function_grads<Vector<double> >
(const Vector<double> &,
 std::vector<std::vector<Tensor<1,deal_II_dimension> > > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_grads<Vector<float> >
(const Vector<float> &,
 std::vector<std::vector<Tensor<1,deal_II_dimension> > > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_grads<BlockVector<double> >
(const BlockVector<double> &,
 std::vector<std::vector<Tensor<1,deal_II_dimension> > > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_grads<BlockVector<float> >
(const BlockVector<float> &,
 std::vector<std::vector<Tensor<1,deal_II_dimension> > > &) const;

#ifdef DEAL_II_USE_PETSC
template
void FEValuesBase<deal_II_dimension>::get_function_grads<PETScWrappers::Vector>
(const PETScWrappers::Vector &,
 std::vector<std::vector<Tensor<1,deal_II_dimension> > > &) const;
#endif

//-----------------------------------------------------------------------------

template
void FEValuesBase<deal_II_dimension>::get_function_2nd_derivatives<Vector<double> >
(const Vector<double> &,
 std::vector<Tensor<2,deal_II_dimension> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_2nd_derivatives<Vector<float> >
(const Vector<float> &,
 std::vector<Tensor<2,deal_II_dimension> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_2nd_derivatives<BlockVector<double> >
(const BlockVector<double> &,
 std::vector<Tensor<2,deal_II_dimension> > &) const;

#ifdef DEAL_II_USE_PETSC
template
void FEValuesBase<deal_II_dimension>::get_function_2nd_derivatives<PETScWrappers::Vector>
(const PETScWrappers::Vector &,
 std::vector<Tensor<2,deal_II_dimension> > &) const;
#endif

//-----------------------------------------------------------------------------

template
void FEValuesBase<deal_II_dimension>::get_function_2nd_derivatives<Vector<double> >
(const Vector<double> &,
 std::vector<std::vector<Tensor<2,deal_II_dimension> > > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_2nd_derivatives<Vector<float> >
(const Vector<float> &,
 std::vector<std::vector<Tensor<2,deal_II_dimension> > > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_2nd_derivatives<BlockVector<double> >
(const BlockVector<double> &,
 std::vector<std::vector<Tensor<2,deal_II_dimension> > > &) const;

#ifdef DEAL_II_USE_PETSC
template
void FEValuesBase<deal_II_dimension>::get_function_2nd_derivatives<PETScWrappers::Vector>
(const PETScWrappers::Vector &,
 std::vector<std::vector<Tensor<2,deal_II_dimension> > > &) const;
#endif

//----------------------------  fe_values.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
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

#include <iomanip>

// if necessary try to work around a bug in the IBM xlC compiler
#ifdef XLC_WORK_AROUND_STD_BUG
using namespace std;
#endif


// TODO:[RH,GK] replace this by non-global object
static const MappingQ1<deal_II_dimension> mapping_q1;

template <int dim>
void
FEValuesData<dim>::initialize (const unsigned int n_quadrature_points,
			       const unsigned int n_shapes,
			       const UpdateFlags  flags)
{
  update_flags = flags;
  
  if (flags & update_values)
    shape_values.reinit(n_shapes, n_quadrature_points);

  if (flags & update_gradients)
    {
      shape_gradients.resize(n_shapes);
      for (unsigned int i=0;i<n_shapes;++i)
	shape_gradients[i].resize(n_quadrature_points);
    }

  if (flags & update_second_derivatives)
    {      
      shape_2nd_derivatives.resize(n_shapes);
      for (unsigned int i=0;i<n_shapes;++i)
	shape_2nd_derivatives[i].resize(n_quadrature_points);
    }
  
  if (flags & update_q_points)
    quadrature_points.resize(n_quadrature_points);

  if (flags & update_JxW_values)
    JxW_values.resize(n_quadrature_points);

  if (flags & update_boundary_forms)
    boundary_forms.resize(n_quadrature_points);

  if (flags & update_normal_vectors)
    normal_vectors.resize(n_quadrature_points);

//Todo:[?] support points missing
}


/*------------------------------- FEValuesBase ---------------------------*/


template <int dim>
FEValuesBase<dim>::FEValuesBase (const unsigned int n_q_points,
				 const unsigned int dofs_per_cell,
				 const unsigned int,
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
  update_flags = flags;
}



template <int dim>
FEValuesBase<dim>::~FEValuesBase ()
{
  if (fe_data)
    {
      typename Mapping<dim>::InternalDataBase *tmp1=fe_data;
      fe_data=0;
      delete tmp1;
    }

  if (mapping_data)
    {
      typename Mapping<dim>::InternalDataBase *tmp1=mapping_data;
      mapping_data=0;
      delete tmp1;
    }
}



template <int dim>
template <class InputVector, typename number>
void FEValuesBase<dim>::get_function_values (const InputVector            &fe_function,
					     typename std::vector<number> &values) const
{
  Assert (update_flags & update_values, ExcAccessToUninitializedField());
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
  fill_n (values.begin(), n_quadrature_points, 0);

				   // add up contributions of trial
				   // functions
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
      values[point] += (dof_values(shape_func) *
			shape_values(shape_func, point));
};



template <int dim>
template <class InputVector, typename number>
void FEValuesBase<dim>::get_function_values (const InputVector                     &fe_function,
					     typename std::vector<Vector<number> > &values) const
{
  Assert (n_quadrature_points == values.size(),
	  ExcWrongVectorSize(values.size(), n_quadrature_points));
  for (unsigned i=0;i<values.size();++i)
    Assert (values[i].size() == fe->n_components(),
	    ExcWrongNoOfComponents());
  Assert (update_flags & update_values, ExcAccessToUninitializedField());
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
				   // functions
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
      values[point](fe->system_to_component_index(shape_func).first)
	+= (dof_values(shape_func) * shape_values(shape_func, point));
};



template <int dim>
const typename FEValuesData<dim>::ShapeVector &
FEValuesBase<dim>::get_shape_values () const
{
  Assert (update_flags & update_values, ExcAccessToUninitializedField());
  return shape_values;
};



template <int dim>
const typename FEValuesData<dim>::GradientVector &
FEValuesBase<dim>::get_shape_grads () const
{
  Assert (update_flags & update_gradients, ExcAccessToUninitializedField());
  return shape_gradients;
};




template <int dim>
const typename std::vector<typename std::vector<Tensor<2,dim> > > &
FEValuesBase<dim>::get_shape_2nd_derivatives () const
{
  Assert (update_flags & update_second_derivatives, ExcAccessToUninitializedField());
  return shape_2nd_derivatives;
};



template <int dim>
const typename std::vector<Point<dim> > &
FEValuesBase<dim>::get_quadrature_points () const
{
  Assert (update_flags & update_q_points, ExcAccessToUninitializedField());
  return quadrature_points;
};



template <int dim>
const std::vector<double> &
FEValuesBase<dim>::get_JxW_values () const
{
  Assert (update_flags & update_JxW_values, ExcAccessToUninitializedField());
  return JxW_values;
}



template <int dim>
template <class InputVector>
void
FEValuesBase<dim>::
get_function_grads (const InputVector                    &fe_function,
		    typename std::vector<Tensor<1,dim> > &gradients) const
{
  Assert (update_flags & update_gradients, ExcAccessToUninitializedField());

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
  fill_n (gradients.begin(), n_quadrature_points, Tensor<1,dim>());

				   // add up contributions of trial
				   // functions
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
      {
	Tensor<1,dim> tmp(shape_gradients[shape_func][point]);
	tmp *= dof_values(shape_func);
	gradients[point] += tmp;
      };
};



template <int dim>
template <class InputVector>
void
FEValuesBase<dim>::
get_function_grads (const InputVector                                           &fe_function,
		    typename std::vector<typename std::vector<Tensor<1,dim> > > &gradients) const
{
  Assert (n_quadrature_points == gradients.size(),
	  ExcWrongNoOfComponents());
  for (unsigned i=0;i<gradients.size();++i)
    Assert (gradients[i].size() == fe->n_components(),
	    ExcWrongVectorSize(gradients[i].size(), fe->n_components()));
  Assert (update_flags & update_gradients, ExcAccessToUninitializedField());
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
    fill_n (gradients[i].begin(), gradients[i].size(), Tensor<1,dim>());

				   // add up contributions of trial
				   // functions
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
      {
	Tensor<1,dim> tmp(shape_gradients[shape_func][point]);
	tmp *= dof_values(shape_func);
	gradients[point][fe->system_to_component_index(shape_func).first]
	  += tmp;
      };
};



template <int dim>
template <class InputVector>
void
FEValuesBase<dim>::
get_function_2nd_derivatives (const InputVector                    &fe_function,
			      typename std::vector<Tensor<2,dim> > &second_derivatives) const
{
  Assert (fe->n_components() == 1,
	  ExcWrongNoOfComponents());
  Assert (second_derivatives.size() == n_quadrature_points,
	  ExcWrongVectorSize(second_derivatives.size(), n_quadrature_points));
  Assert (update_flags & update_second_derivatives, ExcAccessToUninitializedField());
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
				   // functions
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
      {
	Tensor<2,dim> tmp(shape_2nd_derivatives[shape_func][point]);
	tmp *= dof_values(shape_func);
	second_derivatives[point] += tmp;
      };
};



template <int dim>
template <class InputVector>
void
FEValuesBase<dim>::
get_function_2nd_derivatives (const InputVector                                           &fe_function,
			      typename std::vector<typename std::vector<Tensor<2,dim> > > &second_derivs) const
{
  Assert (n_quadrature_points == second_derivs.size(),
	  ExcWrongNoOfComponents());
  for (unsigned i=0;i<second_derivs.size();++i)
    Assert (second_derivs[i].size() == fe->n_components(),
	    ExcWrongVectorSize(second_derivs[i].size(), fe->n_components()));
  Assert (update_flags & update_second_derivatives, ExcAccessToUninitializedField());
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
    fill_n (second_derivs[i].begin(), second_derivs[i].size(), Tensor<2,dim>());

				   // add up contributions of trial
				   // functions
  for (unsigned int point=0; point<n_quadrature_points; ++point)
    for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
      {
	Tensor<2,dim> tmp(shape_2nd_derivatives[shape_func][point]);
	tmp *= dof_values(shape_func);
	second_derivs[point][fe->system_to_component_index(shape_func).first]
	  += tmp;
      };
};



template <int dim>
const Point<dim> &
FEValuesBase<dim>::quadrature_point (const unsigned int i) const
{
  Assert (update_flags & update_q_points, ExcAccessToUninitializedField());
  Assert (i<quadrature_points.size(), ExcIndexRange(i, 0, quadrature_points.size()));
  
  return quadrature_points[i];
};




template <int dim>
double FEValuesBase<dim>::JxW (const unsigned int i) const
{
  Assert (update_flags & update_JxW_values, ExcAccessToUninitializedField());
  Assert (i<JxW_values.size(), ExcIndexRange(i, 0, JxW_values.size()));
  
  return JxW_values[i];
};



template <int dim>
unsigned int
FEValuesBase<dim>::memory_consumption () const
{
//TODO:[WB] check whether this is up to date  
  return (MemoryConsumption::memory_consumption (shape_values) +
	  MemoryConsumption::memory_consumption (shape_gradients) +
	  MemoryConsumption::memory_consumption (shape_2nd_derivatives) +
	  MemoryConsumption::memory_consumption (JxW_values) +
	  MemoryConsumption::memory_consumption (quadrature_points) +
	  sizeof(update_flags) +
	  MemoryConsumption::memory_consumption (present_cell) +
	  MemoryConsumption::memory_consumption (fe));
};



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
};



/*------------------------------- FEValues -------------------------------*/



template <int dim>
FEValues<dim>::FEValues (const Mapping<dim>       &mapping,
			 const FiniteElement<dim> &fe,
			 const Quadrature<dim>    &q,
			 const UpdateFlags         update_flags)
		:
		FEValuesBase<dim> (q.n_quadrature_points,
				   fe.dofs_per_cell,
				   1,
				   update_default,
				   mapping,
				   fe),
		quadrature (q)
{
  initialize (update_flags);
};



template <int dim>
FEValues<dim>::FEValues (const FiniteElement<dim> &fe,
			 const Quadrature<dim>    &q,
			 const UpdateFlags         update_flags)
		:
		FEValuesBase<dim> (q.n_quadrature_points,
				   fe.dofs_per_cell,
				   1,
				   update_default,
				   mapping_q1,
				   fe),
                quadrature (q)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  initialize (update_flags);
};



template <int dim>
void
FEValues<dim>::initialize (const UpdateFlags update_flags)
{
				   // you can't compute normal vectors
				   // on cells, only on faces
  Assert ((update_flags & update_normal_vectors) == false,
	  typename FEValuesBase<dim>::ExcInvalidUpdateFlag());

  const UpdateFlags flags = compute_update_flags (update_flags);
  
				   // then get objects into which the
				   // FE and the Mapping can store
				   // intermediate data used across
				   // calls to reinit
  mapping_data = mapping->get_data(flags, quadrature);
  fe_data      = fe->get_data(flags, *mapping, quadrature);

				   // set up objects within this class
  FEValuesData<dim>::initialize(n_quadrature_points,
				dofs_per_cell,
				flags);
};



template <int dim>
void FEValues<dim>::reinit (const typename DoFHandler<dim>::cell_iterator &cell)
{
				   // assert that the finite elements
				   // passed to the constructor and
				   // used by the DoFHandler used by
				   // this cell, are the same
  Assert (static_cast<const FiniteElementData<dim>&>(*fe) ==
	  static_cast<const FiniteElementData<dim>&>(cell->get_dof_handler().get_fe()),
	  typename FEValuesBase<dim>::ExcFEDontMatch());

  present_cell = cell;

  get_mapping().fill_fe_values(cell,
			       quadrature,
			       *mapping_data,
			       quadrature_points,
			       JxW_values);
  
  get_fe().fill_fe_values(get_mapping(),
			  cell,
			  quadrature,
			  *mapping_data,
			  *fe_data,
			  *this);
}



template <int dim>
unsigned int
FEValues<dim>::memory_consumption () const
{
//TODO:[WB] check whether this is up to date  
  return FEValuesBase<dim>::memory_consumption ();
};


/*------------------------------- FEFaceValuesBase --------------------------*/


template <int dim>
FEFaceValuesBase<dim>::FEFaceValuesBase (const unsigned int n_q_points,
					 const unsigned int dofs_per_cell,
					 const unsigned int n_faces_or_subfaces,
					 const UpdateFlags,
					 const Mapping<dim> &mapping,      
					 const FiniteElement<dim> &fe,
					 const Quadrature<dim-1>& quadrature)
		:
		FEValuesBase<dim> (n_q_points,
				   dofs_per_cell,
				   n_faces_or_subfaces,
				   update_default,
				   mapping,
				   fe),
                quadrature(quadrature)
{};



template <int dim>
const std::vector<Point<dim> > &
FEFaceValuesBase<dim>::get_normal_vectors () const
{
  Assert (update_flags & update_normal_vectors,
	  typename FEValuesBase<dim>::ExcAccessToUninitializedField());
  return normal_vectors;
};



template <int dim>
const std::vector<Tensor<1,dim> > &
FEFaceValuesBase<dim>::get_boundary_forms () const
{
  Assert (update_flags & update_boundary_forms,
	  FEValuesBase<dim>::ExcAccessToUninitializedField());
  return boundary_forms;
};


/*------------------------------- FEFaceValues -------------------------------*/


template <int dim>
FEFaceValues<dim>::FEFaceValues (const Mapping<dim>       &mapping,
				 const FiniteElement<dim> &fe,
				 const Quadrature<dim-1>  &quadrature,
				 const UpdateFlags         update_flags)
		:
		FEFaceValuesBase<dim> (quadrature.n_quadrature_points,
				       fe.dofs_per_cell,
				       GeometryInfo<dim>::faces_per_cell,
				       update_flags,
				       mapping,
				       fe, quadrature)
{
  initialize (update_flags);
};



template <int dim>
FEFaceValues<dim>::FEFaceValues (const FiniteElement<dim> &fe,
				 const Quadrature<dim-1>  &quadrature,
				 const UpdateFlags         update_flags)
		:
		FEFaceValuesBase<dim> (quadrature.n_quadrature_points,
				       fe.dofs_per_cell,
				       GeometryInfo<dim>::faces_per_cell,
				       update_flags,
				       mapping_q1,
				       fe, quadrature)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  initialize (update_flags);
};



template <int dim>
void
FEFaceValues<dim>::initialize (const UpdateFlags update_flags)
{
  const UpdateFlags flags = compute_update_flags (update_flags);
  
				   // then get objects into which the
				   // FE and the Mapping can store
				   // intermediate data used across
				   // calls to reinit
  mapping_data = mapping->get_face_data(flags, quadrature);
  fe_data      = fe->get_face_data(flags, *mapping, quadrature);

				   // set up objects within this class
  FEValuesData<dim>::initialize(n_quadrature_points,
				dofs_per_cell,
				flags);
};



template <int dim>
void FEFaceValues<dim>::reinit (const typename DoFHandler<dim>::cell_iterator &cell,
				const unsigned int              face_no)
{
				   // assert that the finite elements
				   // passed to the constructor and
				   // used by the DoFHandler used by
				   // this cell, are the same
  Assert (static_cast<const FiniteElementData<dim>&>(*fe) ==
	  static_cast<const FiniteElementData<dim>&>(cell->get_dof_handler().get_fe()),
	  typename FEValuesBase<dim>::ExcFEDontMatch());

  present_cell = cell;
  present_face = cell->face(face_no);

  get_mapping().fill_fe_face_values(cell, face_no,
				    quadrature,
				    *mapping_data,
				    quadrature_points,
				    JxW_values,
				    boundary_forms,
				    normal_vectors);
  
  get_fe().fill_fe_face_values(get_mapping(),
			       cell, face_no,
			       quadrature,
			       *mapping_data,
			       *fe_data,
			       *this);
};


/*------------------------------- FESubFaceValues -------------------------------*/


template <int dim>
FESubfaceValues<dim>::FESubfaceValues (const Mapping<dim>       &mapping,
				       const FiniteElement<dim> &fe,
				       const Quadrature<dim-1>  &quadrature,
				       const UpdateFlags         update_flags)
		:
		FEFaceValuesBase<dim> (quadrature.n_quadrature_points,
				       fe.dofs_per_cell,
				       GeometryInfo<dim>::faces_per_cell *
				       GeometryInfo<dim>::subfaces_per_face,
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
				       GeometryInfo<dim>::faces_per_cell *
				       GeometryInfo<dim>::subfaces_per_face,
				       update_flags,
				       mapping_q1,
				       fe, quadrature)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  initialize (update_flags);
};



template <int dim>
void
FESubfaceValues<dim>::initialize (const UpdateFlags update_flags)
{
  const UpdateFlags flags = compute_update_flags (update_flags);
  
				   // then get objects into which the
				   // FE and the Mapping can store
				   // intermediate data used across
				   // calls to reinit
  mapping_data = mapping->get_subface_data(flags, quadrature);
  fe_data      = fe->get_subface_data(flags, *mapping, quadrature);

				   // set up objects within this class
  FEValuesData<dim>::initialize(n_quadrature_points,
				dofs_per_cell,
				flags);
};



template <int dim>
void FESubfaceValues<dim>::reinit (const typename DoFHandler<dim>::cell_iterator &cell,
				   const unsigned int         face_no,
				   const unsigned int         subface_no)
{
				   // assert that the finite elements
				   // passed to the constructor and
				   // used by the DoFHandler used by
				   // this cell, are the same
  Assert (static_cast<const FiniteElementData<dim>&>(*fe) ==
	  static_cast<const FiniteElementData<dim>&>(cell->get_dof_handler().get_fe()),
	  typename FEValuesBase<dim>::ExcFEDontMatch());
  Assert (face_no < GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (face_no, 0, GeometryInfo<dim>::faces_per_cell));
  Assert (subface_no < GeometryInfo<dim>::subfaces_per_face,
	  ExcIndexRange (subface_no, 0, GeometryInfo<dim>::subfaces_per_face));

//TODO:[GK] Reinsert this assertion? It tests a necessary, not a sufficient condition
  //  Assert (cell->face(face_no)->at_boundary() == false,
  //	  ExcReinitCalledWithBoundaryFace());
  
  present_cell  = cell;
  present_face  = cell->face(face_no);

  get_mapping().fill_fe_subface_values(cell, face_no, subface_no,
				       quadrature,
				       *mapping_data,
				       quadrature_points,
				       JxW_values,
				       boundary_forms,
				       normal_vectors);
  
  get_fe().fill_fe_subface_values(get_mapping(),
				  cell, face_no, subface_no,
				  quadrature,
				  *mapping_data,
				  *fe_data,
				  *this);
};


/*------------------------------- Explicit Instantiations -------------*/

template class FEValuesBase<deal_II_dimension>;
template class FEValues<deal_II_dimension>;

#if deal_II_dimension >= 2
template class FEFaceValuesBase<deal_II_dimension>;
template class FEFaceValues<deal_II_dimension>;
template class FESubfaceValues<deal_II_dimension>;
#endif


//-----------------------------------------------------------------------------

template
void FEValuesBase<deal_II_dimension>::get_function_values (const Vector<double> &,
							   std::vector<double>       &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values (const Vector<float> &,
							   std::vector<float>      &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values (const BlockVector<double> &,
							   std::vector<double>      &) const;

//-----------------------------------------------------------------------------

template
void FEValuesBase<deal_II_dimension>::get_function_values (const Vector<double> &,
							   std::vector<Vector<double> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values (const Vector<float> &,
							   std::vector<Vector<float> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_values (const BlockVector<double> &,
							   std::vector<Vector<double> >     &) const;

//-----------------------------------------------------------------------------

template
void FEValuesBase<deal_II_dimension>::get_function_grads (const Vector<double> &,
							  std::vector<Tensor<1,deal_II_dimension> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_grads (const Vector<float> &,
							  std::vector<Tensor<1,deal_II_dimension> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_grads (const BlockVector<double> &,
							  std::vector<Tensor<1,deal_II_dimension> > &) const;

//-----------------------------------------------------------------------------

template
void FEValuesBase<deal_II_dimension>::get_function_grads (const Vector<double> &,
							  std::vector<std::vector<Tensor<1,deal_II_dimension> > > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_grads (const Vector<float> &,
							  std::vector<std::vector<Tensor<1,deal_II_dimension> > > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_grads (const BlockVector<double> &,
							  std::vector<std::vector<Tensor<1,deal_II_dimension> > > &) const;

//-----------------------------------------------------------------------------

template
void FEValuesBase<deal_II_dimension>::get_function_2nd_derivatives (const Vector<double> &,
								    std::vector<Tensor<2,deal_II_dimension> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_2nd_derivatives (const Vector<float> &,
								    std::vector<Tensor<2,deal_II_dimension> > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_2nd_derivatives (const BlockVector<double> &,
								    std::vector<Tensor<2,deal_II_dimension> > &) const;
//-----------------------------------------------------------------------------

template
void FEValuesBase<deal_II_dimension>::get_function_2nd_derivatives (const Vector<double> &,
								    std::vector<std::vector<Tensor<2,deal_II_dimension> > > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_2nd_derivatives (const Vector<float> &,
								    std::vector<std::vector<Tensor<2,deal_II_dimension> > > &) const;
template
void FEValuesBase<deal_II_dimension>::get_function_2nd_derivatives (const BlockVector<double> &,
								    std::vector<std::vector<Tensor<2,deal_II_dimension> > > &) const;

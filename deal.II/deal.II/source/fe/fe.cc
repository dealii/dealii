//----------------------------  fe.cc  ---------------------------
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
//----------------------------  fe.cc  ---------------------------


#include <fe/fe.h>
#include <base/memory_consumption.h>
#include <fe/mapping.h>
#include <fe/mapping_q1.h>
#include <fe/fe_values.h>
#include <base/quadrature.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_boundary.h>

// if necessary try to work around a bug in the IBM xlC compiler
#ifdef XLC_WORK_AROUND_STD_BUG
using namespace std;
#endif

/*------------------------------- FiniteElementBase ----------------------*/


template <int dim>
void
FiniteElementBase<dim>::
InternalDataBase::initialize_2nd (const FiniteElement<dim> *element,
				  const Mapping<dim>       &mapping,
				  const Quadrature<dim>    &quadrature)
{
				   // if we shall compute second
				   // derivatives, then we do so by
				   // finite differencing the
				   // gradients. that we do by
				   // evaluating the gradients of
				   // shape values at points shifted
				   // star-like a little in each
				   // coordinate direction around each
				   // quadrature point.
				   //
				   // therefore generate 2*dim (the
				   // number of evaluation points)
				   // FEValues objects with slightly
				   // shifted positions
  std::vector<Point<dim> > diff_points (quadrature.n_quadrature_points);
  std::vector<double>      diff_weights (quadrature.n_quadrature_points, 0);
  
  differences.resize(2*dim);
  for (unsigned int d=0; d<dim; ++d)
    {
      Point<dim> shift;
//TODO: unify the places where the finite differencing step length is used      
      shift (d) = 1.e-6;

				       // generate points and FEValues
				       // objects shifted in
				       // plus-direction. note that
				       // they only need to compute
				       // gradients, not more
      for (unsigned int i=0; i<diff_points.size(); ++i)
	diff_points[i] = quadrature.point(i) + shift;
      const Quadrature<dim> plus_quad (diff_points, diff_weights);
      differences[d] = new FEValues<dim> (mapping, *element,
					  plus_quad, update_gradients);

				       // now same in minus-direction
      for (unsigned int i=0; i<diff_points.size(); ++i)
	diff_points[i] = quadrature.point(i) - shift;
      const Quadrature<dim> minus_quad (diff_points, diff_weights);
      differences[d+dim] = new FEValues<dim> (mapping, *element,
					      minus_quad, update_gradients); 
    }
}




template <int dim>
FiniteElementBase<dim>::InternalDataBase::~InternalDataBase ()
{
  for (unsigned int i=0; i<differences.size (); ++i)
    if (differences[i] != 0)
      {
					 // delete pointer and set it
					 // to zero to avoid
					 // inadvertant use
	delete differences[i];
	differences[i] = 0;
      };
}




template <int dim>
FiniteElementBase<dim>::FiniteElementBase (const FiniteElementData<dim> &fe_data,
					   const std::vector<bool> &restriction_is_additive_flags)
		:
		FiniteElementData<dim> (fe_data),
                system_to_component_table(dofs_per_cell),
                face_system_to_component_table(dofs_per_face),
                component_to_system_table(components, std::vector<unsigned>(dofs_per_cell)),
		face_component_to_system_table(components, std::vector<unsigned>(dofs_per_face)),
		component_to_base_table(components),
		restriction_is_additive_flags(restriction_is_additive_flags)
{
  Assert(restriction_is_additive_flags.size()==fe_data.components,
	 ExcDimensionMismatch(restriction_is_additive_flags.size(),fe_data.components));

  for (unsigned int i=0; i<GeometryInfo<dim>::children_per_cell; ++i) 
    {
      restriction[i].reinit (dofs_per_cell, dofs_per_cell);
      prolongation[i].reinit (dofs_per_cell, dofs_per_cell);
    };

				   // first set sizes of some
				   // matrices. they will be filled by
				   // derived classes, or re-set to
				   // zero size of not defined
  switch (dim)
    {
      case 1:
	    Assert ((interface_constraints.m() == 0) &&
		    (interface_constraints.n() == 0),
		    ExcInternalError());
	    break;
	    
      case 2:
	    interface_constraints.reinit (dofs_per_vertex+2*dofs_per_line,
					  dofs_per_face);
	    break;

      case 3:
	    interface_constraints.reinit (5*dofs_per_vertex +
					  12*dofs_per_line  +
					  4*dofs_per_quad,
					  dofs_per_face);
	    break;

      default:
	    Assert (false, ExcNotImplemented());
    };
	    
  				   // this is the default way, if there is only
				   // one component; if there are several, then
				   // the constructor of the derived class needs
				   // to fill these arrays
  for (unsigned int j=0 ; j<dofs_per_cell ; ++j)
    {
      system_to_component_table[j] = std::pair<unsigned,unsigned>(0,j);
      component_to_system_table[0][j] = j;
    }
  for (unsigned int j=0 ; j<dofs_per_face ; ++j)
    {
      face_system_to_component_table[j] = std::pair<unsigned,unsigned>(0,j);
      face_component_to_system_table[0][j] = j;
    }
};



template <int dim>
const FullMatrix<double> &
FiniteElementBase<dim>::restrict (const unsigned int child) const
{
  Assert (child<GeometryInfo<dim>::children_per_cell,
	  ExcIndexRange(child, 0, GeometryInfo<dim>::children_per_cell));
  Assert (restriction[child].n() != 0, ExcProjectionVoid());
  return restriction[child];
};



template <int dim>
const FullMatrix<double> &
FiniteElementBase<dim>::prolongate (const unsigned int child) const
{
  Assert (child<GeometryInfo<dim>::children_per_cell,
	  ExcIndexRange(child, 0, GeometryInfo<dim>::children_per_cell));
  Assert (prolongation[child].n() != 0, ExcEmbeddingVoid());
  return prolongation[child];
};



template <int dim>
const FullMatrix<double> &
FiniteElementBase<dim>::constraints () const
{
  Assert ((dofs_per_face  == 0) || (interface_constraints.m() != 0),
	  ExcConstraintsVoid());
  
  if (dim==1)
    Assert ((interface_constraints.m()==0) && (interface_constraints.n()==0),
	    ExcWrongInterfaceMatrixSize(interface_constraints.m(),
					interface_constraints.n()));
  
  return interface_constraints;
};



template <int dim>
bool FiniteElementBase<dim>::operator == (const FiniteElementBase<dim> &f) const
{
  return ((static_cast<const FiniteElementData<dim>&>(*this) ==
	   static_cast<const FiniteElementData<dim>&>(f)) &&
	  (interface_constraints == f.interface_constraints));
};



template <int dim>
const std::vector<Point<dim> > &
FiniteElementBase<dim>::get_unit_support_points () const
{
				   // a finite element may define
				   // support points, but only if
				   // there are as many as there are
				   // degrees of freedom
  Assert ((unit_support_points.size() == 0) ||
	  (unit_support_points.size() == dofs_per_cell),
	  ExcInternalError());
  return unit_support_points;
};



template <int dim>
bool
FiniteElementBase<dim>::has_support_points () const
{
  return (unit_support_points.size() != 0);
};



template <int dim>
const std::vector<Point<dim-1> > &
FiniteElementBase<dim>::get_unit_face_support_points () const
{
				   // a finite element may define
				   // support points, but only if
				   // there are as many as there are
				   // degrees of freedom on a face
  Assert ((unit_face_support_points.size() == 0) ||
	  (unit_face_support_points.size() == dofs_per_face),
	  ExcInternalError());
  return unit_face_support_points;
};



template <int dim>
bool
FiniteElementBase<dim>::has_face_support_points () const
{
  return (unit_face_support_points.size() != 0);
};



template <int dim>
unsigned int
FiniteElementBase<dim>::memory_consumption () const
{
  return (sizeof(FiniteElementData<dim>) +
	  MemoryConsumption::
	  memory_consumption<FullMatrix<double>, sizeof(restriction)/sizeof(restriction[0])>
	  (restriction)+
	  MemoryConsumption::memory_consumption
	  <FullMatrix<double>, sizeof(prolongation)/sizeof(prolongation[0])>
	  (prolongation) +
	  MemoryConsumption::memory_consumption (interface_constraints) +
	  MemoryConsumption::memory_consumption (system_to_component_table) +
	  MemoryConsumption::memory_consumption (face_system_to_component_table) +
	  MemoryConsumption::memory_consumption (component_to_system_table) +
	  MemoryConsumption::memory_consumption (face_component_to_system_table) +
	  MemoryConsumption::memory_consumption (component_to_base_table) +
	  MemoryConsumption::memory_consumption (restriction_is_additive_flags));
};



template <int dim>
void
FiniteElementBase<dim>::
compute_2nd (const Mapping<dim>                   &mapping,
	     const DoFHandler<dim>::cell_iterator &cell,
	     const unsigned int                    offset,
	     Mapping<dim>::InternalDataBase       &mapping_internal,
	     InternalDataBase                     &fe_internal,
	     FEValuesData<dim>                    &data) const
{
  Assert ((fe_internal.update_each | fe_internal.update_once)
	  & update_second_derivatives,
	  ExcInternalError());
  Assert (data.shape_2nd_derivatives.size() == dofs_per_cell,
	  ExcInternalError());
				   // Number of quadrature points
  const unsigned int n_q_points = data.shape_2nd_derivatives[0].size();

				   // first reinit the fe_values
				   // objects used for the finite
				   // differencing stuff
  for (unsigned int d=0; d<dim; ++d)
    {
      fe_internal.differences[d]->reinit(cell);
      fe_internal.differences[d+dim]->reinit(cell);
    }

				   // collection of difference
				   // quotients of gradients in each
				   // direction (first index) and at
				   // all q-points (second index)
  std::vector<std::vector<Tensor<1,dim> > > diff_quot (dim, std::vector<Tensor<1,dim> >(n_q_points));
  std::vector<Tensor<1,dim> >               diff_quot2 (n_q_points);

				   // for all shape functions at all
				   // quadrature points and difference
				   // quotients in all directions:
  for (unsigned int shape=0; shape<dofs_per_cell; ++shape)
    {
      for (unsigned int d1=0; d1<dim; ++d1)
	for (unsigned int q=0; q<n_q_points; ++q)
	  {
					     // get gradient at points
					     // shifted slightly to
					     // the right and to the
					     // left in the present
					     // coordinate direction
	    const Tensor<1,dim>& right
	      = fe_internal.differences[d1]->shape_grad(shape, q);
	    const Tensor<1,dim>& left
	      = fe_internal.differences[d1+dim]->shape_grad(shape, q);

//TODO: unify the places where the finite differencing step length is used      
					     // compute the second
					     // derivative from a
					     // symmetric difference
					     // approximation
	    for (unsigned int d=0; d<dim; ++d)
	      diff_quot[d][q][d1] = 1./(2*1.e-6) * (right[d]-left[d]);
	  }

				       // up to now we still have
				       // difference quotients on the
				       // unit cell, so transform it
				       // to something on the real
				       // cell
      for (unsigned int d=0; d<dim; ++d)
	{
	  mapping.transform_covariant (diff_quot2, diff_quot[d],
				       mapping_internal, offset);

	  for (unsigned int q=0; q<n_q_points; ++q)
	    for (unsigned int d1=0; d1<dim; ++d1)
	      data.shape_2nd_derivatives[shape][q][d][d1] = diff_quot2[q][d1];
	}
    }
}


/*------------------------------- FiniteElement ----------------------*/

template <int dim>
FiniteElement<dim>::FiniteElement (const FiniteElementData<dim> &fe_data,
				   const std::vector<bool> &restriction_is_additive_flags) :
		FiniteElementBase<dim> (fe_data,
					restriction_is_additive_flags)
{}



template <int dim>
FiniteElement<dim>::~FiniteElement ()
{}



template <int dim>
Mapping<dim>::InternalDataBase*
FiniteElement<dim>::get_face_data (const UpdateFlags       flags,
				   const Mapping<dim>      &mapping,
				   const Quadrature<dim-1> &quadrature) const
{
  QProjector<dim> q(quadrature, false);
  return get_data (flags, mapping, q);
}



template <int dim>
Mapping<dim>::InternalDataBase*
FiniteElement<dim>::get_subface_data (const UpdateFlags        flags,
				      const Mapping<dim>      &mapping,
				      const Quadrature<dim-1> &quadrature) const
{
  QProjector<dim> q(quadrature, true);
  return get_data (flags, mapping, q);
}



template <int dim>
unsigned int
FiniteElement<dim>::n_base_elements() const
{
				   // default implementation
  return 1;
}




template <int dim>
unsigned int
FiniteElement<dim>::memory_consumption () const
{
  return FiniteElementBase<dim>::memory_consumption ();
}




template <int dim>
const FiniteElement<dim>&
FiniteElement<dim>::base_element(unsigned index) const
{
  Assert (index==0, ExcIndexRange(index,0,1));
  return *this;
}

/*------------------------------- Explicit Instantiations -------------*/

template class FiniteElementBase<deal_II_dimension>;
template class FiniteElement<deal_II_dimension>;



//----------------------------  fe.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
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
#include <fe/fe_values.h>
#include <base/quadrature.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_boundary.h>

#include <algorithm>
#include <functional>
#include <numeric>



/*------------------------------- FiniteElementBase ----------------------*/


template <int dim>
const double FiniteElementBase<dim>::fd_step_length = 1.0e-6;


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
  
  differences.resize(2*dim);
  for (unsigned int d=0; d<dim; ++d)
    {
      Point<dim> shift;
      shift (d) = fd_step_length;

				       // generate points and FEValues
				       // objects shifted in
				       // plus-direction. note that
				       // they only need to compute
				       // gradients, not more
      for (unsigned int i=0; i<diff_points.size(); ++i)
	diff_points[i] = quadrature.point(i) + shift;
      const Quadrature<dim> plus_quad (diff_points);
      differences[d] = new FEValues<dim> (mapping, *element,
					  plus_quad, update_gradients);

				       // now same in minus-direction
      for (unsigned int i=0; i<diff_points.size(); ++i)
	diff_points[i] = quadrature.point(i) - shift;
      const Quadrature<dim> minus_quad (diff_points);
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
FiniteElementBase<dim>::
FiniteElementBase (const FiniteElementData<dim> &fe_data,
                   const std::vector<bool> &restriction_is_additive_flags,
                   const std::vector<std::vector<bool> > &nonzero_components)
		:
		FiniteElementData<dim> (fe_data),
                system_to_component_table (this->dofs_per_cell),
                face_system_to_component_table(this->dofs_per_face),
                system_to_base_table(this->dofs_per_cell),
                face_system_to_base_table(this->dofs_per_face),		
                component_to_base_table (this->components,
                                         std::make_pair(0U, 0U)),
                restriction_is_additive_flags(restriction_is_additive_flags),
		nonzero_components (nonzero_components),
		n_nonzero_components_table (compute_n_nonzero_components(nonzero_components)),
		cached_primitivity (std::find_if (n_nonzero_components_table.begin(),
                                                  n_nonzero_components_table.end(),
                                                  std::bind2nd(std::not_equal_to<unsigned int>(),
                                                               1U))
                                    ==
                                    n_nonzero_components_table.end())
{
  Assert (restriction_is_additive_flags.size() == this->dofs_per_cell,
	  ExcDimensionMismatch(restriction_is_additive_flags.size(),
			       this->dofs_per_cell));
  Assert (nonzero_components.size() == this->dofs_per_cell,
	  ExcInternalError());
  for (unsigned int i=0; i<nonzero_components.size(); ++i)
    {
      Assert (nonzero_components[i].size() == this->n_components(),
	      ExcInternalError());
      Assert (std::count (nonzero_components[i].begin(),
			  nonzero_components[i].end(),
			  true)
	      >= 1,
	      ExcInternalError());
      Assert (n_nonzero_components_table[i] >= 1,
	      ExcInternalError());
      Assert (n_nonzero_components_table[i] <= this->n_components(),
	      ExcInternalError());      
    };
  
                                   // initialize some tables in the
				   // default way, i.e. if there is
				   // only one (vector-)component; if
				   // there are several, then the
				   // constructor of the derived class
				   // needs to overwrite these arrays
  for (unsigned int j=0 ; j<this->dofs_per_cell ; ++j)
    {
      system_to_component_table[j] = std::pair<unsigned,unsigned>(0,j);
      system_to_base_table[j] = std::make_pair(std::make_pair(0U,0U),j);      
    }
  for (unsigned int j=0 ; j<this->dofs_per_face ; ++j)
    {
      face_system_to_component_table[j] = std::pair<unsigned,unsigned>(0,j);
      face_system_to_base_table[j] = std::make_pair(std::make_pair(0U,0U),j);      
    }
}


template <int dim>
double
FiniteElementBase<dim>::shape_value (const unsigned int,
				     const Point<dim> &) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return 0.;
}



template <int dim>
double
FiniteElementBase<dim>::shape_value_component (const unsigned int,
					       const Point<dim> &,
					       const unsigned int) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return 0.;
}



template <int dim>
Tensor<1,dim>
FiniteElementBase<dim>::shape_grad (const unsigned int,
				    const Point<dim> &) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return Tensor<1,dim> ();
}



template <int dim>
Tensor<1,dim>
FiniteElementBase<dim>::shape_grad_component (const unsigned int,
					      const Point<dim> &,
					      const unsigned int) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return Tensor<1,dim> ();
}



template <int dim>
Tensor<2,dim>
FiniteElementBase<dim>::shape_grad_grad (const unsigned int,
					 const Point<dim> &) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return Tensor<2,dim> ();
}



template <int dim>
Tensor<2,dim>
FiniteElementBase<dim>::shape_grad_grad_component (const unsigned int,
						   const Point<dim> &,
						   const unsigned int) const
{
  AssertThrow(false, ExcUnitShapeValuesDoNotExist());
  return Tensor<2,dim> ();
}


template <int dim>
const FullMatrix<double> &
FiniteElementBase<dim>::restrict (const unsigned int child) const
{
  Assert (child<GeometryInfo<dim>::children_per_cell,
	  ExcIndexRange(child, 0, GeometryInfo<dim>::children_per_cell));
  Assert (restriction[child].n() != 0, ExcProjectionVoid());
  return restriction[child];
}



template <int dim>
const FullMatrix<double> &
FiniteElementBase<dim>::prolongate (const unsigned int child) const
{
  Assert (child<GeometryInfo<dim>::children_per_cell,
	  ExcIndexRange(child, 0, GeometryInfo<dim>::children_per_cell));
  Assert (prolongation[child].n() != 0, ExcEmbeddingVoid());
  return prolongation[child];
}


template <int dim>
bool
FiniteElementBase<dim>::prolongation_is_implemented () const
{
  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
    {
      Assert ((prolongation[c].m() == this->dofs_per_cell) ||
              (prolongation[c].m() == 0),
              ExcInternalError());
      Assert ((prolongation[c].n() == this->dofs_per_cell) ||
              (prolongation[c].n() == 0),
              ExcInternalError());
      if ((prolongation[c].m() == 0) ||
          (prolongation[c].n() == 0))
        return false;
    }
  return true;
}



template <int dim>
bool
FiniteElementBase<dim>::constraints_are_implemented () const
{
  return (this->dofs_per_face  == 0) || (interface_constraints.m() != 0);
}



template <int dim>
const FullMatrix<double> &
FiniteElementBase<dim>::constraints () const
{
  Assert ((this->dofs_per_face  == 0) || (interface_constraints.m() != 0),
          ExcConstraintsVoid());
  
  if (dim==1)
    Assert ((interface_constraints.m()==0) && (interface_constraints.n()==0),
	    ExcWrongInterfaceMatrixSize(interface_constraints.m(),
					interface_constraints.n()));
  
  return interface_constraints;
}



template <int dim>
TableIndices<2>
FiniteElementBase<dim>::interface_constraints_size () const 
{
  switch (dim)
    {
      case 1:
            return TableIndices<2> (0U, 0U);
      case 2:
            return TableIndices<2> (this->dofs_per_vertex +
                                    2*this->dofs_per_line,
                                    this->dofs_per_face);
      case 3:
            return TableIndices<2> (5*this->dofs_per_vertex +
                                    12*this->dofs_per_line  +
                                    4*this->dofs_per_quad,
                                    this->dofs_per_face);
      default:
            Assert (false, ExcNotImplemented());
    };
  return TableIndices<2> (static_cast<unsigned int>(-1),
                          static_cast<unsigned int>(-1));
}



template <int dim>
void
FiniteElementBase<dim>::
get_interpolation_matrix (const FiniteElementBase<dim> &,
			  FullMatrix<double>           &) const
{
				   // by default, no interpolation
				   // implemented. so throw exception,
				   // as documentation says
  AssertThrow (false,
               typename FiniteElementBase<dim>::
               ExcInterpolationNotImplemented());
}

                                   


template <int dim>
bool
FiniteElementBase<dim>::operator == (const FiniteElementBase<dim> &f) const
{
  return ((static_cast<const FiniteElementData<dim>&>(*this) ==
	   static_cast<const FiniteElementData<dim>&>(f)) &&
	  (interface_constraints == f.interface_constraints));
}



template <int dim>
const std::vector<Point<dim> > &
FiniteElementBase<dim>::get_unit_support_points () const
{
				   // a finite element may define
				   // support points, but only if
				   // there are as many as there are
				   // degrees of freedom
  Assert ((unit_support_points.size() == 0) ||
	  (unit_support_points.size() == this->dofs_per_cell),
	  ExcInternalError());
  return unit_support_points;
}



template <int dim>
bool
FiniteElementBase<dim>::has_support_points () const
{
  return (unit_support_points.size() != 0);
}



template <int dim>
Point<dim>
FiniteElementBase<dim>::unit_support_point (const unsigned index) const
{
  Assert (index < this->dofs_per_cell,
          ExcIndexRange (index, 0, this->dofs_per_cell));
  Assert (unit_support_points.size() == this->dofs_per_cell,
          ExcFEHasNoSupportPoints ());
  return unit_support_points[index];
}



template <int dim>
const std::vector<Point<dim-1> > &
FiniteElementBase<dim>::get_unit_face_support_points () const
{
				   // a finite element may define
				   // support points, but only if
				   // there are as many as there are
				   // degrees of freedom on a face
  Assert ((unit_face_support_points.size() == 0) ||
	  (unit_face_support_points.size() == this->dofs_per_face),
	  ExcInternalError());
  return unit_face_support_points;
}



template <int dim>
bool
FiniteElementBase<dim>::has_face_support_points () const
{
  return (unit_face_support_points.size() != 0);
}



template <int dim>
Point<dim-1>
FiniteElementBase<dim>::unit_face_support_point (const unsigned index) const
{
  Assert (index < this->dofs_per_face,
          ExcIndexRange (index, 0, this->dofs_per_face));
  Assert (unit_face_support_points.size() == this->dofs_per_face,
          ExcFEHasNoSupportPoints ());
  return unit_face_support_points[index];
}



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
	  MemoryConsumption::memory_consumption (system_to_base_table) +
	  MemoryConsumption::memory_consumption (face_system_to_base_table) +	  
	  MemoryConsumption::memory_consumption (component_to_base_table) +
	  MemoryConsumption::memory_consumption (restriction_is_additive_flags) +
	  MemoryConsumption::memory_consumption (nonzero_components) +
	  MemoryConsumption::memory_consumption (n_nonzero_components_table));
}



template <int dim>
void
FiniteElementBase<dim>::
compute_2nd (const Mapping<dim>                   &mapping,
	     const typename DoFHandler<dim>::cell_iterator &cell,
	     const unsigned int,
	     typename Mapping<dim>::InternalDataBase &mapping_internal,
	     InternalDataBase                     &fe_internal,
	     FEValuesData<dim>                    &data) const
{
  Assert ((fe_internal.update_each | fe_internal.update_once)
	  & update_second_derivatives,
	  ExcInternalError());

  const unsigned int total_nonzero_components
    = std::accumulate (n_nonzero_components_table.begin(),
                       n_nonzero_components_table.end(),
                       0U);
  
  Assert (data.shape_2nd_derivatives.n_rows() == total_nonzero_components,
	  ExcInternalError());
				   // Number of quadrature points
  const unsigned int n_q_points = data.shape_2nd_derivatives.n_cols();

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
  Table<2,Tensor<1,dim> >    diff_quot (dim, n_q_points);
  std::vector<Tensor<1,dim> > diff_quot2 (n_q_points);

				   // for all nonzero components of
				   // all shape functions at all
				   // quadrature points and difference
				   // quotients in all directions:
  unsigned int total_index = 0;
  for (unsigned int shape_index=0; shape_index<this->dofs_per_cell; ++shape_index)
    for (unsigned int n=0; n<n_nonzero_components(shape_index); ++n, ++total_index)
      {
        for (unsigned int d1=0; d1<dim; ++d1)
          for (unsigned int q=0; q<n_q_points; ++q)
            {
                                               // get gradient at points
                                               // shifted slightly to
                                               // the right and to the
                                               // left in the present
                                               // coordinate direction
                                               //
                                               // note that things
                                               // might be more
                                               // difficult if the
                                               // shape function has
                                               // more than one
                                               // non-zero component,
                                               // so find out about
                                               // the actual component
                                               // if necessary
              Tensor<1,dim> right, left;
              if (is_primitive(shape_index))
                {
                  right = fe_internal.differences[d1]->shape_grad(shape_index, q);
                  left  = fe_internal.differences[d1+dim]->shape_grad(shape_index, q);
                }
              else
                {
                                                   // get the
                                                   // component index
                                                   // of the n-th
                                                   // nonzero
                                                   // compoment
                  unsigned int component=0;
                  for (unsigned int nonzero_comp=0; component<this->n_components();
                       ++component)
                    if (nonzero_components[shape_index][component] == true)
                      {
                        ++nonzero_comp;
                                                         // check
                                                         // whether we
                                                         // have found
                                                         // the
                                                         // component
                                                         // we are
                                                         // looking
                                                         // for. note
                                                         // that
                                                         // nonzero_comp
                                                         // is 1-based
                                                         // by the way
                                                         // we compute
                                                         // it
                        if (nonzero_comp == n+1)
                          break;
                      }
                  Assert (component < this->n_components(),
                          ExcInternalError());

                  right = fe_internal.differences[d1]
                          ->shape_grad_component(shape_index, q, component);
                  left  = fe_internal.differences[d1+dim]
                          ->shape_grad_component(shape_index, q, component);
                };
              
                                               // compute the second
                                               // derivative from a
                                               // symmetric difference
                                               // approximation
              for (unsigned int d=0; d<dim; ++d)
                diff_quot[d][q][d1] = 1./(2*fd_step_length) * (right[d]-left[d]);
            }
        
                                         // up to now we still have
                                         // difference quotients on the
                                         // unit cell, so transform it
                                         // to something on the real
                                         // cell
        for (unsigned int d=0; d<dim; ++d)
          {
            Assert (diff_quot2.size() <=
                    diff_quot[d].size(),
                    ExcInternalError());
            mapping.transform_covariant (&*diff_quot2.begin(), &*diff_quot2.end(),
                                         diff_quot[d].begin(),
                                         mapping_internal);
            
            for (unsigned int q=0; q<n_q_points; ++q)
              for (unsigned int d1=0; d1<dim; ++d1)
                data.shape_2nd_derivatives[total_index][q][d][d1] = diff_quot2[q][d1];
          }
      }
}



template <int dim>
std::vector<unsigned int>
FiniteElementBase<dim>::
compute_n_nonzero_components (const std::vector<std::vector<bool> > &nonzero_components)
{
  std::vector<unsigned int> retval (nonzero_components.size());
  for (unsigned int i=0; i<nonzero_components.size(); ++i)
    retval[i] = std::count (nonzero_components[i].begin(),
			    nonzero_components[i].end(),
			    true);
  return retval;
}



/*------------------------------- FiniteElement ----------------------*/

template <int dim>
FiniteElement<dim>::FiniteElement (const FiniteElementData<dim> &fe_data,
				   const std::vector<bool> &restriction_is_additive_flags,
				   const std::vector<std::vector<bool> > &nonzero_components) :
		FiniteElementBase<dim> (fe_data,
					restriction_is_additive_flags,
					nonzero_components)
{}



template <int dim>
FiniteElement<dim>::~FiniteElement ()
{}



template <int dim>
typename Mapping<dim>::InternalDataBase *
FiniteElement<dim>::get_face_data (const UpdateFlags       flags,
				   const Mapping<dim>      &mapping,
				   const Quadrature<dim-1> &quadrature) const
{
  return get_data (flags, mapping,
		   QProjector<dim>::project_to_all_faces(quadrature));
}



template <int dim>
typename Mapping<dim>::InternalDataBase *
FiniteElement<dim>::get_subface_data (const UpdateFlags        flags,
				      const Mapping<dim>      &mapping,
				      const Quadrature<dim-1> &quadrature) const
{
  return get_data (flags, mapping,
		   QProjector<dim>::project_to_all_subfaces(quadrature));
}



template <int dim>
unsigned int
FiniteElement<dim>::memory_consumption () const
{
  return FiniteElementBase<dim>::memory_consumption ();
}




template <int dim>
const FiniteElement<dim>&
FiniteElement<dim>::base_element(const unsigned index) const
{
  Assert (index==0, ExcIndexRange(index,0,1));
  return *this;
}

/*------------------------------- Explicit Instantiations -------------*/

template class FiniteElementBase<deal_II_dimension>;
template class FiniteElement<deal_II_dimension>;

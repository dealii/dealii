//----------------------------  fe_system.cc  ---------------------------
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
//----------------------------  fe_system.cc  ---------------------------


#include <base/memory_consumption.h>
#include <base/quadrature.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/mapping.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>


// if necessary try to work around a bug in the IBM xlC compiler
#ifdef XLC_WORK_AROUND_STD_BUG
using namespace std;
#endif



/* ----------------------- FESystem::InternalData ------------------- */


template <int dim>
FESystem<dim>::InternalData::InternalData(const unsigned int n_base_elements):
		base_fe_datas(n_base_elements),
		base_fe_values_datas(n_base_elements)
{}



template <int dim>
FESystem<dim>::InternalData::~InternalData()
{
				   // delete pointers and set them to
				   // zero to avoid inadvertant use
  for (unsigned int i=0; i<base_fe_datas.size(); ++i)
    if (base_fe_datas[i])
      {
	delete base_fe_datas[i];
	base_fe_datas[i] = 0;
      };
  
  for (unsigned int i=0; i<base_fe_values_datas.size(); ++i)
    if (base_fe_values_datas[i])
      {
	delete base_fe_values_datas[i];
	base_fe_values_datas[i] = 0;
      };
};



template <int dim>
inline
typename FiniteElementBase<dim>::InternalDataBase &
FESystem<dim>::
InternalData::get_fe_data (const unsigned int base_no) const
{
  Assert(base_no<base_fe_datas.size(),
	 ExcIndexRange(base_no,0,base_fe_datas.size()));
  return *base_fe_datas[base_no];
};



template <int dim>
inline void
FESystem<dim>::
InternalData::set_fe_data (const unsigned int base_no,
			   typename FiniteElementBase<dim>::InternalDataBase *ptr)
{
  Assert(base_no<base_fe_datas.size(),
	 ExcIndexRange(base_no,0,base_fe_datas.size()));
  base_fe_datas[base_no]=ptr;
};



template <int dim>
inline FEValuesData<dim> &
FESystem<dim>::
InternalData::get_fe_values_data (const unsigned int base_no) const
{
  Assert(base_no<base_fe_values_datas.size(),
	 ExcIndexRange(base_no,0,base_fe_values_datas.size()));
  Assert(base_fe_values_datas[base_no]!=0, ExcInternalError());
  return *base_fe_values_datas[base_no];
};



template <int dim>
inline void
FESystem<dim>::
InternalData::set_fe_values_data (const unsigned int base_no,
				  FEValuesData<dim> *ptr)
{
  Assert(base_no<base_fe_values_datas.size(),
	 ExcIndexRange(base_no,0,base_fe_values_datas.size()));
  base_fe_values_datas[base_no]=ptr;
};



template <int dim>
inline void
FESystem<dim>::
InternalData::delete_fe_values_data (const unsigned int base_no)
{
  Assert(base_no<base_fe_values_datas.size(),
	 ExcIndexRange(base_no,0,base_fe_values_datas.size()));
  Assert(base_fe_values_datas[base_no]!=0, ExcInternalError());
  delete base_fe_values_datas[base_no];
  base_fe_values_datas[base_no]=0;
}



/* ---------------------------------- FESystem ------------------- */


template <int dim>
const unsigned int FESystem<dim>::invalid_face_number;


template <int dim>
FESystem<dim>::FESystem (const FiniteElement<dim> &fe, const unsigned int n_elements) :
		FiniteElement<dim> (multiply_dof_numbers(fe, n_elements),
				    compute_restriction_is_additive_flags (fe, n_elements)),
                base_elements(1)
{
  base_elements[0] = ElementPair(fe.clone(), n_elements);
  base_elements[0].first->subscribe ();
  initialize ();
};



template <int dim>
FESystem<dim>::FESystem (const FiniteElement<dim> &fe1, const unsigned int n1,
			 const FiniteElement<dim> &fe2, const unsigned int n2) :
		FiniteElement<dim> (multiply_dof_numbers(fe1, n1, fe2, n2),
				    compute_restriction_is_additive_flags (fe1, n1,
									   fe2, n2)),
                base_elements(2)
{
  base_elements[0] = ElementPair(fe1.clone(), n1);
  base_elements[0].first->subscribe ();
  base_elements[1] = ElementPair(fe2.clone(), n2);
  base_elements[1].first->subscribe ();
  initialize ();
};



template <int dim>
FESystem<dim>::FESystem (const FiniteElement<dim> &fe1, const unsigned int n1,
			 const FiniteElement<dim> &fe2, const unsigned int n2,
			 const FiniteElement<dim> &fe3, const unsigned int n3) :
		FiniteElement<dim> (multiply_dof_numbers(fe1, n1,
							 fe2, n2,
							 fe3, n3),
				    compute_restriction_is_additive_flags (fe1, n1,
									   fe2, n2,
									   fe3, n3)),
                base_elements(3)
{
  base_elements[0] = ElementPair(fe1.clone(), n1);  
  base_elements[0].first->subscribe ();
  base_elements[1] = ElementPair(fe2.clone(), n2);
  base_elements[1].first->subscribe ();
  base_elements[2] = ElementPair(fe3.clone(), n3);
  base_elements[2].first->subscribe ();
  initialize ();
};



template <int dim>
FESystem<dim>::~FESystem ()
{
				   // delete base elements created in
				   // the constructor
  for (unsigned i=0; i<base_elements.size(); ++i)
    {
      base_elements[i].first->unsubscribe();
      delete base_elements[i].first;
      base_elements[i].first = 0;
    }
};



template <int dim>
FiniteElement<dim>*
FESystem<dim>::clone() const
{
  switch (n_base_elements())
    {
      case 1:
	    return new FESystem(base_element(0),
				element_multiplicity(0));
      case 2:
	    return new FESystem(base_element(0),
				element_multiplicity(0),
				base_element(1),
				element_multiplicity(1));
      case 3:
	    return new FESystem(base_element(0),
				element_multiplicity(0),
				base_element(1),
				element_multiplicity(1),
				base_element(2),
				element_multiplicity(2));
      default:
	    Assert(false, ExcNotImplemented());
    }
  return 0;
}



template <int dim>
double
FESystem<dim>::shape_value (const unsigned int i,
			    const Point<dim> &p) const
{
  Assert((i<dofs_per_cell), ExcIndexRange(i, 0, dofs_per_cell));

  std::pair<unsigned,unsigned> comp = system_to_component_index(i);
  
  return base_element(component_to_base_table[comp.first])
    .shape_value(comp.second, p);
}



template <int dim>
Tensor<1,dim>
FESystem<dim>::shape_grad (const unsigned int i,
			   const Point<dim> &p) const
{
  Assert((i<dofs_per_cell), ExcIndexRange(i, 0, dofs_per_cell));

  std::pair<unsigned,unsigned> comp = system_to_component_index(i);
  
  return base_element(component_to_base_table[comp.first])
    .shape_grad(comp.second, p);
}



template <int dim>
Tensor<2,dim>
FESystem<dim>::shape_grad_grad (const unsigned int i,
				const Point<dim> &p) const
{
  Assert((i<dofs_per_cell), ExcIndexRange(i, 0, dofs_per_cell));


  std::pair<unsigned,unsigned> comp = system_to_component_index(i);
  
  return base_element(component_to_base_table[comp.first])
    .shape_grad_grad(comp.second, p);
}



template <int dim>
UpdateFlags
FESystem<dim>::update_once (const UpdateFlags flags) const
{
  UpdateFlags out = update_default;
				   // generate maximal set of flags
				   // that are necessary
  for (unsigned int base_no=0; base_no<n_base_elements(); ++base_no)
    out |= base_element(base_no).update_once(flags);
  return out;
}



template <int dim>
UpdateFlags
FESystem<dim>::update_each (const UpdateFlags flags) const
{
  UpdateFlags out = update_default;
				   // generate maximal set of flags
				   // that are necessary
  for (unsigned int base_no=0; base_no<n_base_elements(); ++base_no)
    out |= base_element(base_no).update_each(flags);

				   // second derivatives are handled
				   // by the top-level finite element,
				   // rather than by the base elements
				   // since it is generated by finite
				   // differencing. if second
				   // derivatives are requested, we
				   // therefore have to set the
				   // respective flag since the base
				   // elements don't have them
  if (flags & update_second_derivatives)
    out |= update_second_derivatives | update_covariant_transformation;
    
  return out;
}



template <int dim>
typename Mapping<dim>::InternalDataBase *
FESystem<dim>::get_data (UpdateFlags      flags,
			 const Mapping<dim>    &mapping,
			 const Quadrature<dim> &quadrature) const
{
  InternalData* data = new InternalData(n_base_elements());

  data->update_once = update_once (flags);
  data->update_each = update_each (flags);
  flags = data->update_once | data->update_each;
  
  UpdateFlags sub_flags = flags;

				   // if second derivatives through
				   // finite differencing is required,
				   // then initialize some objects for
				   // that
  data->compute_second_derivatives = flags & update_second_derivatives;
  if (data->compute_second_derivatives)
    {
				       // delete
				       // update_second_derivatives
				       // from flags list
      sub_flags = UpdateFlags (sub_flags ^ update_second_derivatives);
      data->initialize_2nd (this, mapping, quadrature);
    }
  

				   // get data objects from each of
				   // the base elements and store them
  for (unsigned int base_no=0; base_no<n_base_elements(); ++base_no)
    {
      typename Mapping<dim>::InternalDataBase *base_fe_data_base =
	base_element(base_no).get_data(sub_flags, mapping, quadrature);

      typename FiniteElementBase<dim>::InternalDataBase *base_fe_data =
	dynamic_cast<typename FiniteElementBase<dim>::InternalDataBase *>
	(base_fe_data_base);
      
      data->set_fe_data(base_no, base_fe_data);

				       // make sure that *we* compute
				       // second derivatives, base
				       // elements should not do it
      Assert (!(base_fe_data->update_each & update_second_derivatives),
	      ExcInternalError());
      Assert (!(base_fe_data->update_once & update_second_derivatives),
	      ExcInternalError());
      
				       // The FEValuesData @p{data}
				       // given to the
				       // @p{fill_fe_values} function
				       // includes the FEValuesDatas
				       // of the FESystem. Here the
				       // FEValuesDatas @p{*base_data}
				       // needs to be created that
				       // later will be given to the
				       // @p{fill_fe_values} functions
				       // of the base
				       // elements. @p{base_data->initialize}
				       // cannot be called earlier as
				       // in the @p{fill_fe_values}
				       // function called for the
				       // first cell. This is because
				       // the initialize function
				       // needs the update flags as
				       // argument.
				       //
				       // The pointers @p{base_data}
				       // are stored into the
				       // FESystem::InternalData
				       // @p{data}, similar to the
				       // storing of the
				       // @p{base_fe_data}s.
      FEValuesData<dim> *base_data = new FEValuesData<dim>();
      data->set_fe_values_data(base_no, base_data);
    }
  data->update_flags=data->update_once | data->update_each;
  return data;
}



template <int dim>
void
FESystem<dim>::fill_fe_values (const Mapping<dim>                   &mapping,
			       const typename DoFHandler<dim>::cell_iterator &cell,
			       const Quadrature<dim>                &quadrature,
			       typename Mapping<dim>::InternalDataBase       &mapping_data,
			       typename Mapping<dim>::InternalDataBase       &fe_data,
			       FEValuesData<dim>                    &data) const
{
  compute_fill(mapping, cell, invalid_face_number, invalid_face_number,
	       quadrature, mapping_data, fe_data, data);
};



template <int dim>
void
FESystem<dim>::fill_fe_face_values (const Mapping<dim>                   &mapping,
				    const typename DoFHandler<dim>::cell_iterator &cell,
				    const unsigned int                    face_no,
				    const Quadrature<dim-1>              &quadrature,
				    typename Mapping<dim>::InternalDataBase       &mapping_data,
				    typename Mapping<dim>::InternalDataBase       &fe_data,
				    FEValuesData<dim>                    &data) const
{
  compute_fill(mapping, cell, face_no, invalid_face_number,
	       quadrature, mapping_data, fe_data, data);
};




template <int dim>
void
FESystem<dim>::fill_fe_subface_values (const Mapping<dim>                   &mapping,
				       const typename DoFHandler<dim>::cell_iterator &cell,
				       const unsigned int                    face_no,
				       const unsigned int                    sub_no,
				       const Quadrature<dim-1>              &quadrature,
				       typename Mapping<dim>::InternalDataBase       &mapping_data,
				       typename Mapping<dim>::InternalDataBase       &fe_data,
				       FEValuesData<dim>                    &data) const
{
  compute_fill(mapping, cell, face_no, sub_no,
	       quadrature, mapping_data, fe_data, data);
}



template <int dim>
template <int dim_1>
void
FESystem<dim>::compute_fill (const Mapping<dim>                   &mapping,
			     const typename DoFHandler<dim>::cell_iterator &cell,
			     const unsigned int                    face_no,
			     const unsigned int                    sub_no,
			     const Quadrature<dim_1>              &quadrature,
			     typename Mapping<dim>::InternalDataBase       &mapping_data,
			     typename Mapping<dim>::InternalDataBase       &fedata,
			     FEValuesData<dim>                    &data) const
{       
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  InternalData & fe_data = dynamic_cast<InternalData&> (fedata);
  
				   // Either dim_1==dim (fill_fe_values)
				   // or dim_1==dim-1 (fill_fe_(sub)face_values)
  Assert(dim_1==dim || dim_1==dim-1, ExcInternalError());
  const UpdateFlags flags(dim_1==dim ?
			  fe_data.current_update_flags() :
			  fe_data.update_flags);


  if (flags & (update_values | update_gradients))
    {
      if (fe_data.first_cell)
	{
					   // Initialize the FEValuesDatas
					   // for the base elements.
					   // Originally this is the task
					   // of FEValues::FEValues() but
					   // the latter initializes
					   // the FEValuesDatas only of the
					   // FESystem but not the
					   // FEValuesDatas needed by the
					   // base elements.
	  for (unsigned int base_no=0; base_no<n_base_elements(); ++base_no)
	    {
					       // Pointer needed to get
					       // the update flags of the
					       // base element
	      typename Mapping<dim>::InternalDataBase &
		base_fe_data = fe_data.get_fe_data(base_no);
	      
					       // compute update flags ...
	      const UpdateFlags base_update_flags(mapping_data.update_flags
						  | base_fe_data.update_flags);
	      
					       // Initialize the FEValuesDatas
					       // for the base elements.
	      FEValuesData<dim> &base_data=fe_data.get_fe_values_data(base_no);
	      const FiniteElement<dim> &base_fe=base_element(base_no);
	      base_data.initialize(quadrature.n_quadrature_points,
				   base_fe.dofs_per_cell,
				   base_update_flags);
	    }
	}
      
				       // fill_fe_face_values needs
				       // argument Quadrature<dim-1>
				       // for both cases
				       // dim_1==dim-1 and
				       // dim_1=dim. Hence the
				       // following workaround
      const Quadrature<dim>   *cell_quadrature = 0;
      const Quadrature<dim-1> *face_quadrature = 0;
      
				       // static cast to the
				       // common base class of
				       // quadrature being either
				       // Quadrature<dim> or
				       // Quadrature<dim-1>:
      const Subscriptor* quadrature_base_pointer = &quadrature;

      if (face_no==invalid_face_number)
	{
	  Assert(dim_1==dim, 
		 typename FiniteElementData<dim>::
		 ExcSpaceDimensionMismatch(dim_1,dim));
	  cell_quadrature=dynamic_cast<const Quadrature<dim> *>(quadrature_base_pointer);
	  Assert (cell_quadrature != 0, ExcInternalError());
	}
      else
	{
	  Assert(dim_1==dim-1, 
		 typename FiniteElementData<dim>::
		 ExcSpaceDimensionMismatch(dim_1,dim-1));
	  face_quadrature=dynamic_cast<const Quadrature<dim-1> *>(quadrature_base_pointer);
	  Assert (face_quadrature != 0, ExcInternalError());
	}

      
      for (unsigned int base_no=0, comp=0; base_no<n_base_elements(); ++base_no)
	{
	  const FiniteElement<dim> &base_fe=base_element(base_no);
	  typename FiniteElementBase<dim>::InternalDataBase &
	    base_fe_data = fe_data.get_fe_data(base_no);

					   // Make sure that in the
					   // case of fill_fe_values
					   // the data is only copied
					   // from base_data to data
					   // if base_data is
					   // changed. therefore use
					   // fe_fe_data.current_update_flags()
	  
					   // for the case of
					   // fill_fe_(sub)face_values
					   // the data needs to be
					   // copied from base_data to
					   // data on each face,
					   // therefore use
					   // base_fe_data.update_flags.
	  
					   // Store these flags into
					   // base_flags before
					   // calling
					   // base_fe.fill_fe_([sub]face_)values
					   // as the latter changes
					   // the return value of
					   // base_fe_data.current_update_flags()
	  const UpdateFlags base_flags(dim_1==dim ?
				       base_fe_data.current_update_flags() :
				       base_fe_data.update_flags);
	  
	  FEValuesData<dim> &base_data=fe_data.get_fe_values_data(base_no);

	  if (face_no==invalid_face_number)
	    base_fe.fill_fe_values(mapping, cell,
				   *cell_quadrature, mapping_data, base_fe_data, base_data);
	  else if (sub_no==invalid_face_number)
	    base_fe.fill_fe_face_values(mapping, cell, face_no,
					*face_quadrature, mapping_data, base_fe_data, base_data);
	  else
	    base_fe.fill_fe_subface_values(mapping, cell, face_no, sub_no,
					   *face_quadrature, mapping_data, base_fe_data, base_data);

	  for (unsigned int m=0; m<element_multiplicity(base_no); ++m, ++comp)
	    for (unsigned int point=0; point<quadrature.n_quadrature_points; ++point)
	      for (unsigned int k=0; k<base_fe.dofs_per_cell; ++k)
		{
		  const unsigned int system_index=component_to_system_index(comp,k);
		  if (base_flags & update_values)
		    data.shape_values(system_index, point)=
		      base_data.shape_values(k,point);
		  
		  if (base_flags & update_gradients)
		    data.shape_gradients[system_index][point]=
		      base_data.shape_gradients[k][point];
		
		  if (base_flags & update_second_derivatives)
		    data.shape_2nd_derivatives[system_index][point]=
		      base_data.shape_2nd_derivatives[k][point];
		}
	}
  

      if (fe_data.first_cell)
	{
	  fe_data.first_cell = false;
	  for (unsigned int base_no=0; base_no<n_base_elements(); ++base_no)
	    Assert(fe_data.get_fe_data(base_no).first_cell==false, ExcInternalError());
	  
					   // delete FEValuesDatas that
					   // are not needed any more
	  for (unsigned int base_no=0; base_no<n_base_elements(); ++base_no)
	    {
					       // Pointer needed to get
					       // the update flags of the
					       // base element
	      typename Mapping<dim>::InternalDataBase &base_fe_data=
		fe_data.get_fe_data(base_no);
	      
					       // compute update flags ...
	      UpdateFlags base_flags_each(
		dim_1==dim ?
		mapping_data.update_each | base_fe_data.update_each :
		mapping_data.update_flags | base_fe_data.update_flags);
	      
	      if (base_flags_each==update_default)
		fe_data.delete_fe_values_data(base_no);
	    }
	}
    }
  
  if (fe_data.compute_second_derivatives)
    {
      unsigned int offset = 0;
      if (face_no != invalid_face_number)
	offset = (sub_no == invalid_face_number)
		 ? face_no * quadrature.n_quadrature_points
		 :(face_no * GeometryInfo<dim>::subfaces_per_face
		   + sub_no) * quadrature.n_quadrature_points;  
      compute_2nd (mapping, cell, offset, mapping_data, fe_data, data);
    }
}



template <int dim>
void
FESystem<dim>::build_cell_table()
{
  unsigned total_index = 0;
  for (unsigned int base=0 ; base < n_base_elements() ; ++base)
    for (unsigned int m = 0; m < element_multiplicity(base); ++m)
      component_to_base_table[total_index++] = base;

				   // Initialize index table
				   // Multi-component base elements have
				   // to be thought of.
  
				   // 1. Vertices
  total_index = 0;
  for (unsigned int vertex_number= 0 ;
       vertex_number < GeometryInfo<dim>::vertices_per_cell ;
       ++vertex_number)
    {
      unsigned comp_start = 0;
      for(unsigned int base = 0; base < n_base_elements() ;
	  ++base)
	{
	  for (unsigned int m = 0; m < element_multiplicity(base); ++m)
	    {
	      for (unsigned int local_index = 0 ;
		   local_index < base_element(base).dofs_per_vertex ;
		   ++local_index)
		{
		  system_to_component_table[total_index++]
		    = std::make_pair (comp_start+m,
				      vertex_number*base_element(base).dofs_per_vertex
				      +local_index);
		}
	    }
	  comp_start += element_multiplicity(base);
	}
    }
  
				   // 2. Lines
  for (unsigned int line_number= 0 ; ((line_number != GeometryInfo<dim>::lines_per_cell) &&
				      (GeometryInfo<dim>::lines_per_cell > 0));
       ++line_number)
    {
      unsigned comp_start = 0;
      for(unsigned int base = 0; base < n_base_elements() ;
	  ++base)
	{
	  for (unsigned int m = 0; m < element_multiplicity(base); ++m)
	    {
	      for (unsigned int local_index = 0 ;
		   local_index < base_element(base).dofs_per_line ;
		   ++local_index)
		{
		  system_to_component_table[total_index++]
		    = std::pair<unsigned,unsigned>
		    (comp_start+m,
		     line_number*base_element(base).dofs_per_line
		     +local_index+base_element(base).first_line_index);
		}
	    }
	  comp_start += element_multiplicity(base);
	}
    }
  
				   // 3. Quads
  for (unsigned int quad_number= 0 ;
       ((quad_number != GeometryInfo<dim>::quads_per_cell) &&
	(GeometryInfo<dim>::quads_per_cell > 0));
       ++quad_number)
    {
      unsigned int comp_start = 0;
      for(unsigned int base = 0; base < n_base_elements() ;
	  ++base)
	{
	  for (unsigned int m = 0; m < element_multiplicity(base); ++m)
	    {
	      for (unsigned int local_index = 0 ;
		   local_index < base_element(base).dofs_per_quad ;
		   ++local_index)
		{
		  system_to_component_table[total_index++]
		    = std::make_pair (comp_start+m,
				      quad_number*base_element(base).dofs_per_quad
				      +local_index+base_element(base).first_quad_index);
		}
	    }
	  comp_start += element_multiplicity(base);
	}
    }
  
				   // 4. Hex
  for (unsigned int hex_number= 0 ;
       ((hex_number != GeometryInfo<dim>::hexes_per_cell) &&
	(GeometryInfo<dim>::hexes_per_cell > 0));
       ++hex_number)
    {
      unsigned int comp_start = 0;
      for(unsigned int base = 0; base < n_base_elements() ;
	  ++base)
	{
	  for (unsigned int m = 0; m < element_multiplicity(base); ++m)
	    {
	      for (unsigned int local_index = 0 ;
		   local_index < base_element(base).dofs_per_hex ;
		   ++local_index)
		{
		  system_to_component_table[total_index++]
		    = std::make_pair (comp_start+m,
				      hex_number*base_element(base).dofs_per_hex
				      +local_index+base_element(base).first_hex_index);
		}
	    }
	  comp_start += element_multiplicity(base);
	  
	}
    }
				   // Initialize mapping from component
				   // to base element
				   // Initialize mapping from components to
				   // linear index. Fortunately, this is
				   // the inverse of what we just did.
  for (unsigned int comp=0 ; comp<n_components() ; ++comp)
    component_to_system_table[comp]
      .resize(base_element(component_to_base_table[comp]).dofs_per_cell);

  for (unsigned int sys=0 ; sys < dofs_per_cell ; ++sys)
    component_to_system_table[system_to_component_table[sys].first]
      [system_to_component_table[sys].second] = sys;
}



template <int dim>
void
FESystem<dim>::build_face_table()
{
  unsigned total_index = 0;
  for (unsigned int base=0 ; base < n_base_elements() ; ++base)
    for (unsigned int m = 0; m < element_multiplicity(base); ++m)
      component_to_base_table[total_index++] = base;

				   // Initialize index table
				   // Multi-component base elements have
				   // to be thought of.
  
				   // 1. Vertices
  total_index = 0;
  for (unsigned int vertex_number= 0 ; vertex_number < GeometryInfo<dim>::vertices_per_face ;
       ++vertex_number)
    {
      unsigned int comp_start = 0;
      for(unsigned int base = 0; base < n_base_elements() ;
	  ++base)
	{
	  for (unsigned int m = 0; m < element_multiplicity(base); ++m)
	    {
	      for (unsigned int local_index = 0 ;
		   local_index < base_element(base).dofs_per_vertex ;
		   ++local_index)
		{
		  face_system_to_component_table[total_index++]
		    = std::pair<unsigned,unsigned>
		    (comp_start+m,
		     vertex_number*base_element(base).dofs_per_vertex+local_index);
		}
	    }
	  comp_start += element_multiplicity(base);
	}
    }
  Assert (total_index <= face_system_to_component_table.size(),
	  ExcInternalError());
  
				   // 2. Lines
  for (unsigned line_number= 0 ; ((line_number != GeometryInfo<dim>::lines_per_face) &&
				  (GeometryInfo<dim>::lines_per_cell > 0));
       ++line_number)
    {
      unsigned comp_start = 0;
      for(unsigned base = 0; base < n_base_elements() ;
	  ++base)
	{
	  for (unsigned m = 0; m < element_multiplicity(base); ++m)
	    {
	      for (unsigned local_index = 0 ;
		   local_index < base_element(base).dofs_per_line ;
		   ++local_index)
		{
		  face_system_to_component_table[total_index++]
		    = std::pair<unsigned,unsigned>
		    (comp_start+m,
		     line_number*base_element(base).dofs_per_line
		     +local_index+base_element(base).first_face_line_index);
		}
	    }
	  comp_start += element_multiplicity(base);
	}
    }
  Assert (total_index <= face_system_to_component_table.size(),
	  ExcInternalError());
  
				   // 3. Quads
  for (unsigned quad_number= 0 ; ((quad_number != GeometryInfo<dim>::quads_per_face) &&
				  (GeometryInfo<dim>::quads_per_cell > 0));
       ++quad_number)
    {
      unsigned comp_start = 0;
      for(unsigned base = 0; base < n_base_elements() ;
	  ++base)
	{
	  for (unsigned m = 0; m < element_multiplicity(base); ++m)
	    {
	      for (unsigned local_index = 0 ;
		   local_index < base_element(base).dofs_per_quad ;
		   ++local_index)
		{
		  face_system_to_component_table[total_index++]
		    = std::pair<unsigned,unsigned>
		    (comp_start+m,
		     quad_number*base_element(base).dofs_per_quad
		     +local_index+base_element(base).first_face_quad_index);
		}
	    }
	  comp_start += element_multiplicity(base);
	}
    }
  Assert (total_index <= face_system_to_component_table.size(),
	  ExcInternalError());
  
				   // Initialize mapping from component
				   // to base element
				   // Initialize mapping from components to
				   // linear index. Fortunately, this is
				   // the inverse of what we just did.
  for (unsigned comp=0 ; comp<n_components() ; ++comp)
    face_component_to_system_table[comp]
      .resize(base_element(component_to_base_table[comp]).dofs_per_cell);

  for (unsigned sys=0 ; sys < dofs_per_face ; ++sys)
    face_component_to_system_table[face_system_to_component_table[sys].first]
      [face_system_to_component_table[sys].second] = sys;
}



template <int dim>
void FESystem<dim>::build_interface_constraints () 
{
				   // the layout of the constraints matrix is
				   // described in the FiniteElement class. you
				   // may want to look there first before trying
				   // to understand the following, especially
				   // the mapping of the @p{n} index.
				   //
				   // in order to map it to the fe-system class,
				   // we have to know which base element a
				   // degree of freedom within a vertex, line,
				   // etc belongs to. this can be accomplished
				   // by the system_to_component_index
				   // function in conjunction with the
				   // numbers first_{line,quad,...}_index
  for (unsigned int n=0; n<interface_constraints.n(); ++n)
    for (unsigned int m=0; m<interface_constraints.m(); ++m)
      {
					 // for the pair (n,m) find out which
					 // component they belong to and
					 // the number therein
					 //
					 // first value in pair is component,
					 // second is index
	const std::pair<unsigned int, unsigned int> n_index
	  = face_system_to_component_index (n);

	std::pair<unsigned int, unsigned int> m_index;
	switch (dim)
	  {
	    case 1:
	    {
					       // we should never get here!
					       // (in 1d, the constraints matrix
					       // should be of size zero)
	      Assert (false, ExcInternalError());
	      break;
	    };

	    case 2:
	    {
					       // the indices m=0..d_v-1 are
					       // from the center vertex.
					       // their order is the same
					       // as for the first vertex
					       // of the whole cell, so we
					       // can use the
					       // system_to_component_index
					       // function (using the
					       // face_s_t_c_i function would
					       // yield the same)
	      if (m < dofs_per_vertex)
		m_index = system_to_component_index (m);
	      else
						 // then come the two sets of
						 // line indices
		{
		  const unsigned int index_in_line
		    = (m-dofs_per_vertex) % dofs_per_line;
		  const unsigned int sub_line
		    = (m-dofs_per_vertex) / dofs_per_line;
		  Assert (sub_line < 2, ExcInternalError());
		  
						   // get the component by
						   // asking s_t_c_index and
						   // tweaking the index a bit
		  m_index.first = system_to_component_index
				  (GeometryInfo<2>::vertices_per_cell * dofs_per_vertex
				   + index_in_line).first;
						   // first find out the how-many'th
						   // line index of that component
						   // this was
		  m_index.second = (system_to_component_index
				    (GeometryInfo<2>::vertices_per_cell * dofs_per_vertex
				     + index_in_line).second
				    - base_element (component_to_base_table[m_index.first]).first_line_index)
								    // then add the number of dofs
								    // per vertex to get the index
								    // on the first line
				   + base_element(component_to_base_table[m_index.first]).dofs_per_vertex
								    // if on the second line: add
								    // some more
				   + base_element(component_to_base_table[m_index.first]).dofs_per_line * sub_line;
		};
	      break;
	    };

	    case 3:
	    {
					       // the indices
					       // m=0..5*d_v-1 are
					       // from the center and
					       // the four subline
					       // vertices.  their
					       // order is the same as
					       // for the first vertex
					       // of the whole cell,
					       // so we can use the
					       // simple arithmetic
	      if (m < 5*dofs_per_vertex)
		{
		  m_index.first  = system_to_component_index(m % dofs_per_vertex).first;
		  m_index.second = m / dofs_per_vertex;
		}
	      else
						 // then come the 12 sets of
						 // line indices
		if (m < 5*dofs_per_vertex + 12*dofs_per_line)
		  {   
		    const unsigned int index_in_line
		      = (m-5*dofs_per_vertex) % dofs_per_line;
		    const unsigned int sub_line
		      = (m-5*dofs_per_vertex) / dofs_per_line;
		    Assert (sub_line < 12, ExcInternalError());
		  
						     // get the component by
						     // asking s_t_c_index and
						     // tweaking the index a bit
		    m_index.first = system_to_component_index
				    (GeometryInfo<3>::vertices_per_cell * dofs_per_vertex
				     + index_in_line).first;
		    
						     // first find out the how-many'th
						     // line index of that component
						     // this was
		    m_index.second = (system_to_component_index
				      (GeometryInfo<3>::vertices_per_cell * dofs_per_vertex
				       + index_in_line).second
				      - base_element (component_to_base_table[m_index.first]).first_line_index)
								      // then add the number of dofs
								      // for the five vertices to get
								      // the index on the first line
				     + 5*base_element(component_to_base_table[m_index.first]).dofs_per_vertex
								      // and correct for the
								      // how-many'th line
				     + base_element(component_to_base_table[m_index.first]).dofs_per_line * sub_line;
		  }
		else
						   // on one of the four sub-quads
		  {   
		    const unsigned int index_in_quad
		      = (m-5*dofs_per_vertex-12*dofs_per_line) % dofs_per_line;
		    const unsigned int sub_quad
		      = (m-5*dofs_per_vertex-12*dofs_per_line) / dofs_per_line;
		    Assert (sub_quad < 4, ExcInternalError());
		  
						     // get the component by
						     // asking s_t_c_index and
						     // tweaking the index a bit
		    m_index.first = system_to_component_index
				    (GeometryInfo<3>::vertices_per_cell * dofs_per_vertex
				     + GeometryInfo<3>::lines_per_cell * dofs_per_line
				     + index_in_quad).first;
		    
						     // first find out the how-many'th
						     // quad index of that component
						     // this was
		    m_index.second = (system_to_component_index
				      (GeometryInfo<3>::vertices_per_cell * dofs_per_vertex
				       + GeometryInfo<3>::lines_per_cell * dofs_per_line
				       + index_in_quad).second
				      - base_element (component_to_base_table[m_index.first]).first_quad_index)
								      // then add the number of dofs
								      // for the five vertices and 12 lines
								      // to get the index on the first quad
				     + 5*base_element(component_to_base_table[m_index.first]).dofs_per_vertex
				     + 12*base_element(component_to_base_table[m_index.first]).dofs_per_line
								      // and correct for the
								      // how-many'th line
				     + base_element(component_to_base_table[m_index.first]).dofs_per_quad * sub_quad;
		  };
	      
	      break;
	    };
		  
	    default:
		  Assert (false, ExcNotImplemented());
	  };

					 // now that we gathered all information:
					 // use it to build the matrix. note
					 // that if n and m belong to different
					 // components, there definitely will be
					 // no coupling
	if (n_index.first == m_index.first)
	  interface_constraints(m,n)
	    = (base_element(component_to_base_table[n_index.first])
	       .constraints()(m_index.second,
			      n_index.second));
      };
};



template <int dim>
void FESystem<dim>::initialize ()
{
  build_cell_table();
  build_face_table();
  
				   // Check if some of the matrices of
				   // the base elements are void.
  bool do_restriction = true;
  bool do_prolongation = true;

  for (unsigned int i=0; i<n_base_elements(); ++i)
    {
      if (base_element(i).restriction[0].n() == 0)
	do_restriction = false;
      if (base_element(i).prolongation[0].n() == 0)
	do_prolongation = false;
    }
  
				   // if we encountered void matrices,
				   // disable them for the composite
				   // element as well
  if (!do_restriction)
    for (unsigned int i=0;i<GeometryInfo<dim>::children_per_cell;++i)
      restriction[i].reinit(0,0);
  if (!do_prolongation)
    for (unsigned int i=0;i<GeometryInfo<dim>::children_per_cell;++i)
      prolongation[i].reinit(0,0);
	
				   // distribute the matrices of the base
				   // finite elements to the matrices of
				   // this object
  for (unsigned int component=0; component<n_components(); ++component)
				     // transform restriction and
				     // prolongation matrices
    for (unsigned int i=0; i<base_element(component_to_base_table[component]).dofs_per_cell; ++i)
      for (unsigned int j=0; j<base_element(component_to_base_table[component]).dofs_per_cell; ++j)
					 // only fill block diagonals, no
					 // intermixing of subelements
	for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
	  {
	    if (do_restriction)
	      restriction[child] (component_to_system_index (component,i),
				  component_to_system_index (component, j))
		= base_element(component_to_base_table[component]).restrict(child)(i,j);
	    if (do_prolongation)
	      prolongation[child] (component_to_system_index (component,i),
				   component_to_system_index (component, j))
		= base_element(component_to_base_table[component]).prolongate(child)(i,j);
	  };


				   // now set up the interface constraints.
				   // this is kind'o hairy, so don't try
				   // to do it dimension independent
  build_interface_constraints ();

				   // finally fill in support points
				   // on cell and face
  initialize_unit_support_points ();
  initialize_unit_face_support_points ();
};



template <int dim>
FiniteElementData<dim>
FESystem<dim>::multiply_dof_numbers (const FiniteElementData<dim> &fe_data,
				     const unsigned int            N)
{
  std::vector<unsigned int> dpo;
  dpo.push_back(fe_data.dofs_per_vertex * N);
  dpo.push_back(fe_data.dofs_per_line * N);
  if (dim>1) dpo.push_back(fe_data.dofs_per_quad * N);
  if (dim>2) dpo.push_back(fe_data.dofs_per_hex * N);
  
  return FiniteElementData<dim> (dpo, fe_data.n_components() * N);
};




template <int dim>
void
FESystem<dim>::
initialize_unit_support_points ()
{
				       // if one of the base elements
				       // has no support points, then
				       // it makes no sense to define
				       // support points for the
				       // composed element, so return
				       // an empty array to
				       // demonstrate that
				       // fact
  for (unsigned int base_el=0 ; base_el<n_base_elements(); ++base_el)
    if (!base_element(base_el).has_support_points())
      {
	unit_support_points.resize(0);
	return;
      };

				   // generate unit support points
				   // from unit support points of sub
				   // elements
  unit_support_points.resize(dofs_per_cell);
  
  unsigned int comp = 0;
  for (unsigned int base_el=0; base_el<n_base_elements(); ++base_el)
    {
				       // we know that there are
				       // support points on the cell,
				       // so collect them
      const unsigned int
	base_element_dofs_per_cell = base_element(base_el).dofs_per_cell;
 
      const typename std::vector<Point<dim> >
	& base_unit_support_points = base_element(base_el).get_unit_support_points ();
      
				       // otherwise distribute the
				       // support points of this base
				       // element to all degrees of
				       // freedom contributed by this
				       // base element
      Assert(base_unit_support_points.size()==base_element_dofs_per_cell,
	     ExcInternalError());
      for (unsigned int n=0; n<element_multiplicity(base_el); ++n, ++comp)
	for (unsigned int i=0; i<base_element_dofs_per_cell; ++i)
	  unit_support_points[component_to_system_index(comp,i)]
	    = base_unit_support_points[i];
    }
}


#if deal_II_dimension == 1

template <>
void
FESystem<1>::
initialize_unit_face_support_points ()
{
				   // no faces no work
};

#endif


template <int dim>
void
FESystem<dim>::
initialize_unit_face_support_points ()
{
				       // if one of the base elements
				       // has no support points, then
				       // it makes no sense to define
				       // support points for the
				       // composed element, so return
				       // an empty array to
				       // demonstrate that fact (note
				       // that we ask whether the base
				       // element has no support
				       // points at all, not only none
				       // on the face!)
  for (unsigned int base_el=0 ; base_el<n_base_elements(); ++base_el)
    if (!base_element(base_el).has_support_points())
      {
	unit_face_support_points.resize(0);
	return;
      };
  

				   // generate unit face support points
				   // from unit support points of sub
				   // elements
  unit_face_support_points.resize(dofs_per_face);
  
  unsigned int comp = 0;
  for (unsigned int base_el=0 ; base_el<n_base_elements(); ++base_el)
    {
				       // in some cases, finite
				       // elements have support points
				       // (we have made sure that they
				       // have above) but don't have
				       // any on the face (e.g. DG
				       // elements). in that case,
				       // don't even bother with this
				       // base element and directly go
				       // to the next one:
      if (!base_element(base_el).has_face_support_points())
	{
	  comp += element_multiplicity(base_el);
	  continue;
	};

				       // otherwise, we know that
				       // there are support points on
				       // the face, so collect them
      const unsigned int
	base_element_dofs_per_face = base_element(base_el).dofs_per_face;
 
      const typename std::vector<Point<dim-1> > &
	base_unit_support_points = base_element(base_el).get_unit_face_support_points ();
      
				       // distribute the support
				       // points of this base element
				       // to all degrees of freedom
				       // contributed by this base
				       // element
      Assert(base_unit_support_points.size()==base_element_dofs_per_face,
	     ExcNotImplemented());
      for (unsigned int n=0; n<element_multiplicity(base_el); ++n, ++comp)
	for (unsigned int i=0; i<base_element_dofs_per_face; ++i)
	  unit_face_support_points[face_component_to_system_index(comp,i)]
	    = base_unit_support_points[i];
    }
}



template <int dim>
FiniteElementData<dim>
FESystem<dim>::multiply_dof_numbers (const FiniteElementData<dim> &fe1,
				     const unsigned int            N1,
				     const FiniteElementData<dim> &fe2,
				     const unsigned int            N2)
{
  std::vector<unsigned int> dpo;
  dpo.push_back(fe1.dofs_per_vertex * N1 + fe2.dofs_per_vertex * N2);
  dpo.push_back(fe1.dofs_per_line * N1 + fe2.dofs_per_line * N2);
  if (dim>1) dpo.push_back(fe1.dofs_per_quad * N1 + fe2.dofs_per_quad * N2);
  if (dim>2) dpo.push_back(fe1.dofs_per_hex * N1 + fe2.dofs_per_hex * N2);
  
  return FiniteElementData<dim> (dpo,
				 fe1.n_components() * N1 +
				 fe2.n_components() * N2);
};



template <int dim>
FiniteElementData<dim>
FESystem<dim>::multiply_dof_numbers (const FiniteElementData<dim> &fe1,
				     const unsigned int            N1,
				     const FiniteElementData<dim> &fe2,
				     const unsigned int            N2,
				     const FiniteElementData<dim> &fe3,
				     const unsigned int            N3)
{
  std::vector<unsigned int> dpo;
  dpo.push_back(fe1.dofs_per_vertex * N1 +
		fe2.dofs_per_vertex * N2 +
		fe3.dofs_per_vertex * N3);
  dpo.push_back(fe1.dofs_per_line * N1 +
		fe2.dofs_per_line * N2 +
		fe3.dofs_per_line * N3);
  if (dim>1) dpo.push_back(fe1.dofs_per_quad * N1 +
			   fe2.dofs_per_quad * N2 +
			   fe3.dofs_per_quad * N3);
  if (dim>2) dpo.push_back(fe1.dofs_per_hex * N1 +
			   fe2.dofs_per_hex * N2 +
			   fe3.dofs_per_hex * N3);
  
  return FiniteElementData<dim> (dpo,
				 fe1.n_components() * N1 +
				 fe2.n_components() * N2 +
				 fe3.n_components() * N3);
};



template <int dim>
std::vector<bool>
FESystem<dim>::compute_restriction_is_additive_flags (const FiniteElement<dim> &fe,
						      const unsigned int n_elements) 
{
  std::vector<bool> tmp;
  for (unsigned int i=0; i<n_elements; ++i)
    for (unsigned int component=0; component<fe.n_components(); ++component)
      tmp.push_back (fe.restriction_is_additive (component));
  return tmp;
};



template <int dim>
std::vector<bool>
FESystem<dim>::compute_restriction_is_additive_flags (const FiniteElement<dim> &fe1,
						      const unsigned int        N1,
						      const FiniteElement<dim> &fe2,
						      const unsigned int        N2) 
{
  std::vector<bool> tmp;
  for (unsigned int i=0; i<N1; ++i)
    for (unsigned int component=0; component<fe1.n_components(); ++component)
      tmp.push_back (fe1.restriction_is_additive (component));
  for (unsigned int i=0; i<N2; ++i)
    for (unsigned int component=0; component<fe2.n_components(); ++component)
      tmp.push_back (fe2.restriction_is_additive (component));
  return tmp;
};



template <int dim>
std::vector<bool>
FESystem<dim>::compute_restriction_is_additive_flags (const FiniteElement<dim> &fe1,
						      const unsigned int        N1,
						      const FiniteElement<dim> &fe2,
						      const unsigned int        N2,
						      const FiniteElement<dim> &fe3,
						      const unsigned int        N3) 
{
  std::vector<bool> tmp;
  for (unsigned int i=0; i<N1; ++i)
    for (unsigned int component=0; component<fe1.n_components(); ++component)
      tmp.push_back (fe1.restriction_is_additive (component));
  for (unsigned int i=0; i<N2; ++i)
    for (unsigned int component=0; component<fe2.n_components(); ++component)
      tmp.push_back (fe2.restriction_is_additive (component));
  for (unsigned int i=0; i<N3; ++i)
    for (unsigned int component=0; component<fe3.n_components(); ++component)
      tmp.push_back (fe3.restriction_is_additive (component));
  return tmp;
};



template <int dim>
unsigned int
FESystem<dim>::memory_consumption () const 
{
                                 // neglect size of data stored in
                                 // @p{base_elements} due to some
                                 // problems with teh
                                 // compiler. should be neglectable
                                 // after all, considering the size
                                 // of the data of the subelements
  unsigned int mem = (FiniteElement<dim>::memory_consumption () +
		      sizeof (base_elements));
  for (unsigned int i=0; i<base_elements.size(); ++i)
    mem += MemoryConsumption::memory_consumption (*base_elements[i].first);
  return mem;
};




// explicit instantiations
template class FESystem<deal_II_dimension>;



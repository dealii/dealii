//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/quadrature_lib.h>
#include <base/memory_consumption.h>
#include <base/thread_management.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/petsc_vector.h>
#include <lac/petsc_block_vector.h>
#include <numerics/data_out.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <hp/dof_handler.h>
#include <grid/tria_iterator.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <hp/fe_values.h>
#include <fe/mapping_q1.h>

#include <sstream>

DEAL_II_NAMESPACE_OPEN


template <class DH, int patch_dim, int patch_space_dim>
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
DataEntryBase::DataEntryBase (const std::vector<std::string> &names_in,
			      const std::vector<DataComponentInterpretation::DataComponentInterpretation> &data_component_interpretation)
		:
		names(names_in),
		data_component_interpretation (data_component_interpretation),
		postprocessor(0),
		n_output_variables(names.size())
{
  Assert (names.size() == data_component_interpretation.size(),
	  ExcDimensionMismatch(data_component_interpretation.size(),
			       names.size()));

				   // check that the names use only allowed
				   // characters
				   // check names for invalid characters
  for (unsigned int i=0; i<names.size(); ++i)
    Assert (names[i].find_first_not_of("abcdefghijklmnopqrstuvwxyz"
				       "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
				       "0123456789_<>()") == std::string::npos,
	    ExcInvalidCharacter (names[i],
				 names[i].find_first_not_of("abcdefghijklmnopqrstuvwxyz"
							    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
							    "0123456789_<>()")));  
}



template <class DH, int patch_dim, int patch_space_dim>
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
DataEntryBase::DataEntryBase (const DataPostprocessor<DH::space_dimension> *data_postprocessor)
		:
		names(data_postprocessor->get_names()),
		data_component_interpretation (data_postprocessor->get_data_component_interpretation()),
		postprocessor(data_postprocessor),
		n_output_variables(data_postprocessor->n_output_variables())
{
				   // if there is a post processor, then we
				   // should have gotten the names from the
				   // postprocessor. check that the number of
				   // elements in the names vector is
				   // correct. otherwise there is nothing for
				   // us to check
  Assert(data_postprocessor->n_output_variables()==names.size(),
	 ExcDimensionMismatch(data_postprocessor->n_output_variables(),
			      names.size()));
  Assert (names.size() == data_component_interpretation.size(),
	  ExcDimensionMismatch(data_component_interpretation.size(),
			       names.size()));

				   // check that the names use only allowed
				   // characters
				   // check names for invalid characters
  for (unsigned int i=0; i<names.size(); ++i)
    Assert (names[i].find_first_not_of("abcdefghijklmnopqrstuvwxyz"
				       "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
				       "0123456789_<>()") == std::string::npos,
	    ExcInvalidCharacter (names[i],
				 names[i].find_first_not_of("abcdefghijklmnopqrstuvwxyz"
							    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
							    "0123456789_<>()")));  
}



template <class DH, int patch_dim, int patch_space_dim>
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
DataEntryBase::~DataEntryBase ()
{}



template <class DH, int patch_dim, int patch_space_dim>
template <typename VectorType>
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
DataEntry<VectorType>::
DataEntry (const VectorType                       *data,
	   const std::vector<std::string>         &names,
	   const std::vector<DataComponentInterpretation::DataComponentInterpretation> &data_component_interpretation)
		:
                DataEntryBase (names, data_component_interpretation),
		vector (data)
{}



template <class DH, int patch_dim, int patch_space_dim>
template <typename VectorType>
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
DataEntry<VectorType>::
DataEntry (const VectorType                       *data,
	   const DataPostprocessor<DH::space_dimension> *data_postprocessor)
		:
                DataEntryBase (data_postprocessor),
		vector (data)
{}



template <class DH, int patch_dim, int patch_space_dim>
template <typename VectorType>
double
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
DataEntry<VectorType>::
get_cell_data_value (const unsigned int cell_number) const
{
  return (*vector)(cell_number);
}



template <class DH, int patch_dim, int patch_space_dim>
template <typename VectorType>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
DataEntry<VectorType>::
get_function_values (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                     std::vector<Vector<double> >    &patch_values_system) const
{
  fe_patch_values.get_function_values (*vector, patch_values_system);
}



template <class DH,
	  int patch_dim, int patch_space_dim>
template <typename VectorType>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
DataEntry<VectorType>::
get_function_values (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                     std::vector<double>             &patch_values) const
{
  fe_patch_values.get_function_values (*vector, patch_values);
}



template <class DH, int patch_dim, int patch_space_dim>
template <typename VectorType>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
DataEntry<VectorType>::
get_function_gradients (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
			std::vector<std::vector<Tensor<1,DH::space_dimension> > >   &patch_gradients_system) const
{
  fe_patch_values.get_function_grads (*vector, patch_gradients_system);
}



template <class DH,
	  int patch_dim, int patch_space_dim>
template <typename VectorType>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
DataEntry<VectorType>::
get_function_gradients (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
			std::vector<Tensor<1,DH::space_dimension> >       &patch_gradients) const
{
  fe_patch_values.get_function_grads (*vector, patch_gradients);
}



template <class DH, int patch_dim, int patch_space_dim>
template <typename VectorType>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
DataEntry<VectorType>::
get_function_hessians (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
		       std::vector<std::vector<Tensor<2,DH::space_dimension> > >   &patch_hessians_system) const
{
  fe_patch_values.get_function_2nd_derivatives (*vector, patch_hessians_system);
}



template <class DH,
	  int patch_dim, int patch_space_dim>
template <typename VectorType>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
DataEntry<VectorType>::
get_function_hessians (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
		       std::vector<Tensor<2,DH::space_dimension> >       &patch_hessians) const
{
  fe_patch_values.get_function_2nd_derivatives (*vector, patch_hessians);
}



template <class DH,
	  int patch_dim, int patch_space_dim>
template <typename VectorType>
unsigned int
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
DataEntry<VectorType>::memory_consumption () const
{
  return (sizeof (vector) +
	  MemoryConsumption::memory_consumption (this->names));
}



template <class DH, int patch_dim, int patch_space_dim>
template <typename VectorType>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
DataEntry<VectorType>::clear ()
{
  vector = 0;
}




template <class DH,
	  int patch_dim, int patch_space_dim>
DataOut_DoFData<DH,patch_dim,patch_space_dim>::DataOut_DoFData ()
                :
		dofs(0,typeid(*this).name())
{}



template <class DH, int patch_dim, int patch_space_dim>
DataOut_DoFData<DH,patch_dim,patch_space_dim>::~DataOut_DoFData ()
{
  clear ();
}



template <class DH, int patch_dim, int patch_space_dim>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
attach_dof_handler (const DH &d)
{
  Assert (dof_data.size() == 0, ExcOldDataStillPresent());
  Assert (cell_data.size() == 0, ExcOldDataStillPresent());
  
  dofs = &d;
}




template <class DH,
	  int patch_dim, int patch_space_dim>
template <class VECTOR>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
add_data_vector (const VECTOR                             &vec,
		 const std::string                        &name,
		 const DataVectorType                      type,
		 const std::vector<DataComponentInterpretation::DataComponentInterpretation> &data_component_interpretation)
{
  Assert (dofs != 0, ExcNoDoFHandlerSelected ());
  const unsigned int n_components = dofs->get_fe().n_components ();

  std::vector<std::string> names;
				   // if only one component or vector
				   // is cell vector: we only need one
				   // name
  if ((n_components == 1) ||
      (vec.size() == dofs->get_tria().n_active_cells()))
    {
      names.resize (1, name);
    }
  else
				     // otherwise append _i to the
				     // given name
    {
      names.resize (n_components);
      for (unsigned int i=0; i<n_components; ++i)
	{
	  std::ostringstream namebuf;
	  namebuf << '_' << i;
  	  names[i] = name + namebuf.str();
  	}
    }
  
  add_data_vector (vec, names, type, data_component_interpretation);
}



template <class DH,
	  int patch_dim, int patch_space_dim>
template <class VECTOR>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
add_data_vector (const VECTOR                             &vec,
		 const std::vector<std::string>           &names,
		 const DataVectorType                      type,
		 const std::vector<DataComponentInterpretation::DataComponentInterpretation> &data_component_interpretation_)
{
  Assert (dofs != 0, ExcNoDoFHandlerSelected ());

  const std::vector<DataComponentInterpretation::DataComponentInterpretation> &
    data_component_interpretation
    = (data_component_interpretation_.size() != 0
       ?
       data_component_interpretation_
       :
       std::vector<DataComponentInterpretation::DataComponentInterpretation>
       (names.size(), DataComponentInterpretation::component_is_scalar));
  
				   // either cell data and one name,
				   // or dof data and n_components names
  DataVectorType actual_type = type;
  if (type == type_automatic)
    {
      if (vec.size() == dofs->get_tria().n_active_cells())
	actual_type = type_cell_data;
      else
	actual_type = type_dof_data;
    }
  
  switch (actual_type)
    {
      case type_cell_data:
	    Assert (vec.size() == dofs->get_tria().n_active_cells(),
		    ExcInvalidVectorSize (vec.size(),
					  dofs->n_dofs(),
					  dofs->get_tria().n_active_cells()));
	    Assert (names.size() == 1,
		    ExcInvalidNumberOfNames (names.size(), 1));
	    break;
      
      case type_dof_data:
	    Assert (vec.size() == dofs->n_dofs(),
		    ExcInvalidVectorSize (vec.size(),
					  dofs->n_dofs(),
					  dofs->get_tria().n_active_cells()));
	    Assert (names.size() == dofs->get_fe().n_components(),
		    ExcInvalidNumberOfNames (names.size(), dofs->get_fe().n_components()));
	    break;

      case type_automatic:
					     // this case should have
					     // been handled above...
	    Assert (false, ExcInternalError());
    }

  DataEntryBase * new_entry = new DataEntry<VECTOR>(&vec, names,
						    data_component_interpretation);
  if (actual_type == type_dof_data)
    dof_data.push_back (std_cxx0x::shared_ptr<DataEntryBase>(new_entry));
  else
    cell_data.push_back (std_cxx0x::shared_ptr<DataEntryBase>(new_entry));
}



template <class DH,
	  int patch_dim, int patch_space_dim>
template <class VECTOR>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
add_data_vector (const VECTOR                           &vec,
		 const DataPostprocessor<DH::space_dimension> &data_postprocessor)
{
				   // this is a specialized version of the
				   // other function where we have a
				   // postprocessor. if we do, we know that we
				   // have type_dof_data, which makes things a
				   // bit simpler, we also don't need to deal
				   // with some of the other stuff and use a
				   // different constructor of DataEntry
  
  Assert (dofs != 0, ExcNoDoFHandlerSelected ());

  Assert (vec.size() == dofs->n_dofs(),
	  ExcInvalidVectorSize (vec.size(),
				dofs->n_dofs(),
				dofs->get_tria().n_active_cells()));

  DataEntryBase * new_entry = new DataEntry<VECTOR>(&vec, &data_postprocessor);
  dof_data.push_back (std_cxx0x::shared_ptr<DataEntryBase>(new_entry));
}



template <class DH,
	  int patch_dim, int patch_space_dim>
void DataOut_DoFData<DH,patch_dim,patch_space_dim>::clear_data_vectors ()
{
  dof_data.erase (dof_data.begin(), dof_data.end());
  cell_data.erase (cell_data.begin(), cell_data.end());

				   // delete patches
  std::vector<Patch> dummy;
  patches.swap (dummy);
}



template <class DH,
	  int patch_dim, int patch_space_dim>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
clear_input_data_references ()
{
  for (unsigned int i=0; i<dof_data.size(); ++i)
    dof_data[i]->clear ();
  
  for (unsigned int i=0; i<cell_data.size(); ++i)
    cell_data[i]->clear ();

  if (dofs != 0)
    dofs = 0;
}



template <class DH,
          int patch_dim, int patch_space_dim>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::clear ()
{
  dof_data.erase (dof_data.begin(), dof_data.end());
  cell_data.erase (cell_data.begin(), cell_data.end());

  if (dofs != 0)
    dofs = 0;

				   // delete patches
  std::vector<Patch> dummy;
  patches.swap (dummy);
}



template <class DH,
	  int patch_dim, int patch_space_dim>
std::vector<std::string>
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
get_dataset_names () const 
{
  std::vector<std::string> names;
				   // collect the names of dof
				   // and cell data
  typedef
    typename std::vector<std_cxx0x::shared_ptr<DataEntryBase> >::const_iterator
    data_iterator;
  
  for (data_iterator  d=dof_data.begin();
       d!=dof_data.end(); ++d)
    for (unsigned int i=0; i<(*d)->names.size(); ++i)
      names.push_back ((*d)->names[i]);
  for (data_iterator d=cell_data.begin(); d!=cell_data.end(); ++d)
    {
      Assert ((*d)->names.size() == 1, ExcInternalError());
      names.push_back ((*d)->names[0]);
    }

  return names;
}



template <class DH,
	  int patch_dim, int patch_space_dim>
std::vector<std_cxx0x::tuple<unsigned int, unsigned int, std::string> >
DataOut_DoFData<DH,patch_dim,patch_space_dim>::get_vector_data_ranges () const
{
  std::vector<std_cxx0x::tuple<unsigned int, unsigned int, std::string> >
    ranges;
  
				   // collect the ranges of dof
				   // and cell data
  typedef
    typename std::vector<std_cxx0x::shared_ptr<DataEntryBase> >::const_iterator
    data_iterator;

  unsigned int output_component = 0;
  for (data_iterator  d=dof_data.begin();
       d!=dof_data.end(); ++d)
    for (unsigned int i=0; i<(*d)->n_output_variables;
	 ++i, ++output_component)
				       // see what kind of data we have
				       // here. note that for the purpose of
				       // the current function all we care
				       // about is vector data
      if ((*d)->data_component_interpretation[i] ==
	  DataComponentInterpretation::component_is_part_of_vector)
	{
					   // ensure that there is a
					   // continuous number of next
					   // space_dim components that all
					   // deal with vectors
	  Assert (i+patch_space_dim <=
		  (*d)->n_output_variables,
		  ExcInvalidVectorDeclaration (i,
					       (*d)->names[i]));
	  for (unsigned int dd=1; dd<patch_space_dim; ++dd)
	    Assert ((*d)->data_component_interpretation[i+dd]
		    ==
		    DataComponentInterpretation::component_is_part_of_vector,
		    ExcInvalidVectorDeclaration (i,
						 (*d)->names[i]));

					   // all seems alright, so figure out
					   // whether there is a common name
					   // to these components. if not,
					   // leave the name empty and let the
					   // output format writer decide what
					   // to do here
	  std::string name = (*d)->names[i];
	  for (unsigned int dd=1; dd<patch_space_dim; ++dd)
	    if (name != (*d)->names[i+dd])
	      {
		name = "";
		break;
	      }

					   // finally add a corresponding
					   // range
	  std_cxx0x::tuple<unsigned int, unsigned int, std::string>
	    range (output_component,
		   output_component+patch_space_dim-1,
		   name);
	  
	  ranges.push_back (range);

					   // increase the 'component' counter
					   // by the appropriate amount, same
					   // for 'i', since we have already
					   // dealt with all these components
	  output_component += patch_space_dim-1;
	  i += patch_space_dim-1;
	}

				   // note that we do not have to traverse the
				   // list of cell data here because cell data
				   // is one value per (logical) cell and
				   // therefore cannot be a vector

				   // as a final check, the 'component'
				   // counter should be at the total number of
				   // components added up now
#ifdef DEBUG
  unsigned int n_output_components = 0;
  for (data_iterator  d=dof_data.begin();
       d!=dof_data.end(); ++d)
    n_output_components += (*d)->n_output_variables;
  Assert (output_component == n_output_components,
	  ExcInternalError());
#endif
  
  return ranges;
}



template <class DH,
	  int patch_dim, int patch_space_dim>
const std::vector< dealii::DataOutBase::Patch<patch_dim, patch_space_dim> > &
DataOut_DoFData<DH,patch_dim,patch_space_dim>::get_patches () const
{
  return patches;
}



template <class DH,
	  int patch_dim, int patch_space_dim>
unsigned int
DataOut_DoFData<DH,patch_dim,patch_space_dim>::memory_consumption () const
{
  return (DataOutInterface<patch_dim,patch_space_dim>::memory_consumption () +
	  MemoryConsumption::memory_consumption (dofs) +
	  MemoryConsumption::memory_consumption (patches));
}



/* ---------------------------------------------------------------------- */



template <int dim, class DH>
void DataOut<dim,DH>::build_some_patches (internal::DataOut::ParallelData<DH::dimension, DH::space_dimension> &data)
{
				   // Check consistency of redundant
				   // template parameter
  Assert (dim==DH::dimension, ExcDimensionMismatch(dim, DH::dimension));
  
  QTrapez<1>     q_trapez;
  QIterated<DH::dimension> patch_points (q_trapez, data.n_subdivisions);

// We use the mapping to transform the vertex coordinates and the shape
// functions (necessary for example for Raviart-Thomas elements). On the
// boundary, general mappings do not reduce to a MappingQ1, therefore the mapped
// (quadrature) points are stored in the patch, whereas for cells in the
// interior of the domain these points are obtained by a dim-linear mapping and
// can be recovered from the vertices later on, thus they need not to be stored.

				   // create collection objects from
				   // single quadratures,
				   // mappings and finite elements. if we have
				   // an hp DoFHandler,
				   // dof_handler.get_fe() returns a
				   // collection of which we do a
				   // shallow copy instead
  const hp::QCollection<DH::dimension>       q_collection (patch_points);
  const hp::FECollection<DH::dimension,DH::space_dimension>      fe_collection(this->dofs->get_fe());
  const hp::MappingCollection<DH::dimension,DH::space_dimension> mapping_collection(*(data.mapping));

  UpdateFlags update_flags=update_values;
  if (curved_cell_region != no_curved_cells)
    update_flags |= update_quadrature_points;
  
  for (unsigned int i=0; i<this->dof_data.size(); ++i)
    if (this->dof_data[i]->postprocessor)
      update_flags |= this->dof_data[i]->postprocessor->get_needed_update_flags();
				   // perhaps update_normal_vectors is present,
				   // which would only be useful on faces, but
				   // we may not use it here.
  Assert (!(update_flags & update_normal_vectors),
	  ExcMessage("The update of normal vectors may not be requested for evaluation of data on cells via DataPostprocessor."));
  
  hp::FEValues<DH::dimension,DH::space_dimension> x_fe_patch_values (mapping_collection,
						 fe_collection,
						 q_collection,
						 update_flags);

  const unsigned int n_q_points = patch_points.n_quadrature_points;
  
  typename std::vector< dealii::DataOutBase::Patch<DH::dimension,DH::space_dimension> >::iterator
    patch = this->patches.begin();
  cell_iterator cell=first_cell();

				   // keep track of the index of the
				   // current cell so we can
				   // efficiently evaluate cell-based
				   // data (as opposed to DoF-based
				   // data). we do so only if
				   // this->cell_data.size() != 0
  unsigned int cell_index = (this->cell_data.size() != 0
			     ?
			     std::distance (this->dofs->begin_active(),
					    active_cell_iterator (cell))
			     :
			     numbers::invalid_unsigned_int);
  
  
				   // get first cell in this thread
  for (unsigned int i=0; (i<data.this_thread)&&(cell != this->dofs->end()); ++i)
    {
      ++patch;

      const cell_iterator old_cell = cell;
      
      cell = next_cell(cell);

      if (this->cell_data.size() != 0)
	cell_index += std::distance (active_cell_iterator(old_cell),
				     active_cell_iterator(cell));
    }

  				   // now loop over all cells and
				   // actually create the patches
  for (; cell != this->dofs->end();)
    {
      Assert (patch != this->patches.end(), ExcInternalError());

				       // use ucd_to_deal map as patch
				       // vertices are in the old,
				       // unnatural ordering
      for (unsigned int vertex=0; vertex<GeometryInfo<DH::dimension>::vertices_per_cell; ++vertex)
	patch->vertices[vertex] = data.mapping->transform_unit_to_real_cell
				  (cell, GeometryInfo<DH::dimension>::unit_cell_vertex (vertex));
      
      if (data.n_datasets > 0)
	{
          x_fe_patch_values.reinit (cell);
          const FEValues<DH::dimension,DH::space_dimension> &fe_patch_values
            = x_fe_patch_values.get_present_fe_values ();
	  
					   // depending on the requested output
					   // of curved cells, if necessary
					   // append the quadrature points to
					   // the last rows of the patch->data
					   // member. THis is the case if we
					   // want to produce curved cells at
					   // the boundary and this cell
					   // actually is at the boundary, or
					   // else if we want to produce curved
					   // cells everywhere
	  if (curved_cell_region==curved_inner_cells ||
	      (curved_cell_region==curved_boundary && cell->at_boundary()))
	    {
	      Assert(patch->space_dim==dim, ExcInternalError());
	      const std::vector<Point<DH::space_dimension> > & q_points=fe_patch_values.get_quadrature_points();
					       // resize the patch->data member
					       // in order to have enough memory
					       // for the quadrature points as
					       // well
	      patch->data.reinit(patch->data.size(0)+dim,patch->data.size(1));
					       // set the flag indicating that
					       // for this cell the points are
					       // explicitly given
	      patch->points_are_available=true;
					       // copy points to patch->data
	      for (unsigned int i=0; i<dim; ++i)
		for (unsigned int q=0; q<n_q_points; ++q)
		  patch->data(patch->data.size(0)-dim+i,q)=q_points[q][i];
	    }

					   // counter for data records
	  unsigned int offset=0;
	  
					   // first fill dof_data
	  for (unsigned int dataset=0; dataset<this->dof_data.size(); ++dataset)
	    {
	      const DataPostprocessor<DH::space_dimension> *postprocessor=this->dof_data[dataset]->postprocessor;
	      
	      if (postprocessor != 0)
		{
						   // we have to postprocess the
						   // data, so determine, which
						   // fields have to be updated
		  const UpdateFlags update_flags=postprocessor->get_needed_update_flags();
		  if (data.n_components == 1)
		    {
						       // at each point there is
						       // only one component of
						       // value, gradient etc.
		      if (update_flags & update_values)
			this->dof_data[dataset]->get_function_values (fe_patch_values,
								      data.patch_values);
		      if (update_flags & update_gradients)
			this->dof_data[dataset]->get_function_gradients (fe_patch_values,
									 data.patch_gradients);
		      if (update_flags & update_hessians)
			this->dof_data[dataset]->get_function_hessians (fe_patch_values,
									data.patch_hessians);
		      postprocessor->
			compute_derived_quantities_scalar(data.patch_values,
							  data.patch_gradients,
							  data.patch_hessians,
							  data.dummy_normals,
							  data.postprocessed_values[dataset]);
		    }
		  else
		    {
						       // at each point there is
						       // a vector valued
						       // function and its
						       // derivative...
		      if (update_flags & update_values)
			this->dof_data[dataset]->get_function_values (fe_patch_values,
								      data.patch_values_system);
		      if (update_flags & update_gradients)
			this->dof_data[dataset]->get_function_gradients (fe_patch_values,
									 data.patch_gradients_system);
		      if (update_flags & update_hessians)
			this->dof_data[dataset]->get_function_hessians (fe_patch_values,
									data.patch_hessians_system);
		      postprocessor->
			compute_derived_quantities_vector(data.patch_values_system,
							  data.patch_gradients_system,
							  data.patch_hessians_system,
							  data.dummy_normals,
							  data.postprocessed_values[dataset]);
		    }
		  
		  for (unsigned int q=0; q<n_q_points; ++q)
		    for (unsigned int component=0; component<this->dof_data[dataset]->n_output_variables;++component)
		      patch->data(offset+component,q)= data.postprocessed_values[dataset][q](component);
		}
	      else
						 // now we use the given data
						 // vector without
						 // modifications. again, we
						 // treat single component
						 // functions separately for
						 // efficiency reasons.
		if (data.n_components == 1)
		  {
		    this->dof_data[dataset]->get_function_values (fe_patch_values,
								  data.patch_values);
		    for (unsigned int q=0; q<n_q_points; ++q)
		      patch->data(offset,q) = data.patch_values[q];
		  }
		else
		  {
		    this->dof_data[dataset]->get_function_values (fe_patch_values,
								  data.patch_values_system);
		    for (unsigned int component=0; component<data.n_components;
			 ++component)
		      for (unsigned int q=0; q<n_q_points; ++q)
			patch->data(offset+component,q) =
			  data.patch_values_system[q](component);
		  }
					       // increment the counter for the
					       // actual data record
	      offset+=this->dof_data[dataset]->n_output_variables;
	    }

					   // then do the cell data. only
					   // compute the number of a cell if
					   // needed; also make sure that we
					   // only access cell data if the
					   // first_cell/next_cell functions
					   // only return active cells
          if (this->cell_data.size() != 0)
            {
              Assert (!cell->has_children(), ExcNotImplemented());
              
              for (unsigned int dataset=0; dataset<this->cell_data.size(); ++dataset)
                {
                  const double value
                    = this->cell_data[dataset]->get_cell_data_value (cell_index);
                  for (unsigned int q=0; q<n_q_points; ++q)
                    patch->data(offset+dataset,q) =
                      value;
                }
            }
	}


      for (unsigned int f=0; f<GeometryInfo<DH::dimension>::faces_per_cell; ++f)
        {
                                           // let's look up whether
                                           // the neighbor behind that
                                           // face is noted in the
                                           // table of cells which we
                                           // treat. this can only
                                           // happen if the neighbor
                                           // exists, and is on the
                                           // same level as this cell,
                                           // but it may also happen
                                           // that the neighbor is not
                                           // a member of the range of
                                           // cells over which we
                                           // loop, in which case the
                                           // respective entry in the
                                           // cell_to_patch_index_map
                                           // will have the value
                                           // no_neighbor. (note that
                                           // since we allocated only
                                           // as much space in this
                                           // array as the maximum
                                           // index of the cells we
                                           // loop over, not every
                                           // neighbor may have its
                                           // space in it, so we have
                                           // to assume that it is
                                           // extended by values
                                           // no_neighbor)
          if (cell->at_boundary(f)
              ||
              (cell->neighbor(f)->level() != cell->level()))
            continue;

          const cell_iterator neighbor = cell->neighbor(f);
          Assert (static_cast<unsigned int>(neighbor->level()) <
                  data.cell_to_patch_index_map->size(),
                  ExcInternalError());
          if ((static_cast<unsigned int>(neighbor->index()) >=
               (*data.cell_to_patch_index_map)[neighbor->level()].size())
              ||
              ((*data.cell_to_patch_index_map)[neighbor->level()][neighbor->index()]
               ==
               dealii::DataOutBase::Patch<DH::dimension>::no_neighbor))
            continue;

                                           // now, there is a
                                           // neighbor, so get its
                                           // patch number and set it
                                           // for the neighbor index
          patch->neighbors[f] = this->patches[(*data.cell_to_patch_index_map)
					      [neighbor->level()][neighbor->index()]].patch_index;
        }
      
      				       // next cell (patch) in this
      				       // thread
      for (unsigned int i=0;
	   (i<data.n_threads)&&(cell != this->dofs->end()); ++i)
	{
	  ++patch;

	  const cell_iterator old_cell = cell;
	  
          cell = next_cell(cell);

	  if (this->cell_data.size() != 0)
	    cell_index += std::distance (active_cell_iterator(old_cell),
					 active_cell_iterator(cell));
	}
    }
}



template <int dim, class DH>
void DataOut<dim,DH>::build_patches (const unsigned int n_subdivisions,
				     const unsigned int n_threads_) 
{
  build_patches (StaticMappingQ1<DH::dimension,DH::space_dimension>::mapping,
		 n_subdivisions, n_threads_, no_curved_cells);
}


template <int dim, class DH>
void DataOut<dim,DH>::build_patches (const Mapping<DH::dimension,DH::space_dimension> &mapping,
				     const unsigned int nnnn_subdivisions,
				     const unsigned int n_threads_,
				     const CurvedCellRegion curved_region) 
{
  unsigned int n_subdivisions = (nnnn_subdivisions != 0)
				? nnnn_subdivisions
				: this->default_subdivisions;
				   // store the region in which cells shall be
				   // curved. If only one subdivision is
				   // requested then there is no need to do this
				   // at all
  curved_cell_region=curved_region;
  if (n_subdivisions<2)
    curved_cell_region=no_curved_cells;
  
  Assert (n_subdivisions >= 1,
	  ExcInvalidNumberOfSubdivisions(n_subdivisions));

  typedef DataOut_DoFData<DH, DH::dimension, DH::space_dimension> BaseClass;
  Assert (this->dofs != 0, typename BaseClass::ExcNoDoFHandlerSelected());

  Assert (!DEAL_II_USE_MT || (n_threads_ >= 1),
	  ExcMessage ("Must run with at least one thread!"));
  
  const unsigned int n_threads = (DEAL_II_USE_MT ? n_threads_ : 1);

 				   // before we start the loop:
 				   // create a quadrature rule that
 				   // actually has the points on this
 				   // patch
  QTrapez<1>     q_trapez;
  QIterated<DH::dimension> patch_points (q_trapez, n_subdivisions);
  
  const unsigned int n_q_points     = patch_points.n_quadrature_points;
  const unsigned int n_components   = this->dofs->get_fe().n_components();
  unsigned int n_datasets=this->cell_data.size();
  for (unsigned int i=0; i<this->dof_data.size(); ++i)
    n_datasets+= this->dof_data[i]->n_output_variables;
  
				   // clear the patches array
  if (true)
    {
      std::vector< dealii::DataOutBase::Patch<DH::dimension, DH::space_dimension> > dummy;
      this->patches.swap (dummy);
    };
  
				   // first count the cells we want to
				   // create patches of. also fill the
				   // object that maps the cell
				   // indices to the patch numbers, as
				   // this will be needed for
				   // generation of neighborship
				   // information
  std::vector<std::vector<unsigned int> > cell_to_patch_index_map;
  cell_to_patch_index_map.resize (this->dofs->get_tria().n_levels());
  for (unsigned int l=0; l<this->dofs->get_tria().n_levels(); ++l) 
    {
      unsigned int max_index = 0;
      for (cell_iterator cell=first_cell(); cell != this->dofs->end();
           cell = next_cell(cell))
        if (static_cast<unsigned int>(cell->level()) == l)
          max_index = std::max (max_index,
                                static_cast<unsigned int>(cell->index()));
      
      cell_to_patch_index_map[l].resize (max_index+1,
                                         dealii::DataOutBase::Patch<DH::dimension,DH::space_dimension>::no_neighbor);
    };
                                  
  unsigned int n_patches = 0;
  for (cell_iterator cell=first_cell(); cell != this->dofs->end();
       cell = next_cell(cell))
    {
      Assert (static_cast<unsigned int>(cell->level()) <
              cell_to_patch_index_map.size(),
              ExcInternalError());
      Assert (static_cast<unsigned int>(cell->index()) <
              cell_to_patch_index_map[cell->level()].size(),
              ExcInternalError());
      
      cell_to_patch_index_map[cell->level()][cell->index()] = n_patches;
      ++n_patches;
    };

				   // create the patches with default
				   // values. allocate as many patches
				   // as are needed, as this reduces
				   // expensive copying when push_back
				   // or similar operations are used
				   // which would regularly overflow
				   // the allocated amount of memory
                                   //
                                   // then number the patches
                                   // consecutively
  dealii::DataOutBase::Patch<DH::dimension,DH::space_dimension>  default_patch;
  default_patch.n_subdivisions = n_subdivisions;
  default_patch.data.reinit (n_datasets, n_q_points);
  this->patches.insert (this->patches.end(), n_patches, default_patch);

  for (unsigned int i=0; i<this->patches.size(); ++i)
    this->patches[i].patch_index = i;
  

				   // init data for the threads    
  std::vector<internal::DataOut::ParallelData<DH::dimension, DH::space_dimension> > thread_data(n_threads);
  for (unsigned int i=0;i<n_threads;++i)
    {
      thread_data[i].n_threads      = n_threads;
      thread_data[i].this_thread    = i;
      thread_data[i].n_components   = n_components;
      thread_data[i].n_datasets     = n_datasets;
      thread_data[i].n_subdivisions = n_subdivisions;
      thread_data[i].patch_values.resize (n_q_points);
      thread_data[i].patch_values_system.resize (n_q_points);
      thread_data[i].patch_gradients.resize (n_q_points);
      thread_data[i].patch_gradients_system.resize (n_q_points);
      thread_data[i].patch_hessians.resize (n_q_points);
      thread_data[i].patch_hessians_system.resize (n_q_points);
      thread_data[i].dummy_normals.clear();
      thread_data[i].postprocessed_values.resize(this->dof_data.size());
      thread_data[i].mapping        = &mapping;

      thread_data[i].cell_to_patch_index_map = &cell_to_patch_index_map;
      
      for (unsigned int k=0; k<n_q_points; ++k)
	{
	  thread_data[i].patch_values_system[k].reinit(n_components);
	  thread_data[i].patch_gradients_system[k].resize(n_components);
	  thread_data[i].patch_hessians_system[k].resize(n_components);
	}

      for (unsigned int dataset=0; dataset<this->dof_data.size(); ++dataset)
	if (this->dof_data[dataset]->postprocessor)
	  thread_data[i].postprocessed_values[dataset].resize(n_q_points,Vector<double>(this->dof_data[dataset]->n_output_variables));
    }

  Threads::ThreadGroup<> threads;  
  for (unsigned int l=0; l<n_threads; ++l)
    threads += Threads::spawn (*this, &DataOut<dim,DH>::build_some_patches)(thread_data[l]);
  threads.join_all();
}



template <int dim, class DH>
typename DataOut<dim,DH>::cell_iterator
DataOut<dim,DH>::first_cell () 
{
  return this->dofs->begin_active ();
}



template <int dim, class DH>
typename DataOut<dim,DH>::cell_iterator
DataOut<dim,DH>::next_cell (const typename DataOut<dim,DH>::cell_iterator &cell) 
{
				   // convert the iterator to an
				   // active_iterator and advance
				   // this to the next active cell
  typename DH::active_cell_iterator active_cell = cell;
  ++active_cell;
  return active_cell;
}


// explicit instantiations 

#define INSTANTIATE(DH,D0,D1,D2,D3,V)		\
template void \
DataOut_DoFData<DH<D0,D1>,D2,D3>::	    \
add_data_vector<V > (const V             &, \
                     const std::string   &, \
                     const DataVectorType,  \
		     const std::vector<DataComponentInterpretation::DataComponentInterpretation> &); \
\
template void \
 DataOut_DoFData<DH<D0,D1>,D2,D3>::		       \
add_data_vector<V > (const V                        &, \
                     const std::vector<std::string> &, \
                     const DataVectorType,  \
		     const std::vector<DataComponentInterpretation::DataComponentInterpretation> &); \
\
template void \
 DataOut_DoFData<DH<D0,D1>,D2,D3>::		 \
add_data_vector<V > (const V                  &, \
		     const DataPostprocessor<DH<D0,D1>::space_dimension> &)

#ifndef DEAL_II_USE_PETSC
# define INSTANTIATE_VECTORS(DH,D0,D1,D2,D3)	\
    INSTANTIATE(DH,D0,D1,D2,D3,Vector<double>); \
    INSTANTIATE(DH,D0,D1,D2,D3,Vector<float>); \
    INSTANTIATE(DH,D0,D1,D2,D3,BlockVector<double>) ; \
    INSTANTIATE(DH,D0,D1,D2,D3,BlockVector<float>)
#else
# define INSTANTIATE_VECTORS(DH,D0,D1,D2,D3)	 \
    INSTANTIATE(DH,D0,D1,D2,D3,Vector<double>); \
    INSTANTIATE(DH,D0,D1,D2,D3,Vector<float>); \
    INSTANTIATE(DH,D0,D1,D2,D3,BlockVector<double>) ; \
    INSTANTIATE(DH,D0,D1,D2,D3,BlockVector<float>); \
    INSTANTIATE(DH,D0,D1,D2,D3,PETScWrappers::Vector);		\
    INSTANTIATE(DH,D0,D1,D2,D3,PETScWrappers::BlockVector)
#endif

// now do actual instantiations, first for DoFHandler...
template class DataOut_DoFData<DoFHandler<deal_II_dimension>,deal_II_dimension>;
template class DataOut_DoFData<DoFHandler<deal_II_dimension>,deal_II_dimension+1>;
INSTANTIATE_VECTORS(DoFHandler,deal_II_dimension,deal_II_dimension,deal_II_dimension,deal_II_dimension);
INSTANTIATE_VECTORS(DoFHandler,deal_II_dimension,deal_II_dimension,deal_II_dimension+1,deal_II_dimension+1);

// for 3d, also generate an extra class
#if deal_II_dimension >= 2
template class DataOut_DoFData<DoFHandler<deal_II_dimension>,deal_II_dimension-1,deal_II_dimension>;
INSTANTIATE_VECTORS(DoFHandler,deal_II_dimension,deal_II_dimension,deal_II_dimension-1,deal_II_dimension);
#endif

template class DataOut<deal_II_dimension, DoFHandler<deal_II_dimension> >;


// ...and now for hp::DoFHandler
template class DataOut_DoFData<hp::DoFHandler<deal_II_dimension>,deal_II_dimension>;
template class DataOut_DoFData<hp::DoFHandler<deal_II_dimension>,deal_II_dimension+1>;
INSTANTIATE_VECTORS(hp::DoFHandler,deal_II_dimension,deal_II_dimension,deal_II_dimension,deal_II_dimension);
INSTANTIATE_VECTORS(hp::DoFHandler,deal_II_dimension,deal_II_dimension,deal_II_dimension+1,deal_II_dimension+1);

#if deal_II_dimension >= 2
template class DataOut_DoFData<hp::DoFHandler<deal_II_dimension>,deal_II_dimension-1,deal_II_dimension>;
INSTANTIATE_VECTORS(hp::DoFHandler,deal_II_dimension,deal_II_dimension,deal_II_dimension-1,deal_II_dimension);
#endif

template class DataOut<deal_II_dimension, hp::DoFHandler<deal_II_dimension> >;

#if deal_II_dimension == 2 || deal_II_dimension ==1
// now do actual instantiations, first for DoFHandler...
template class DataOut_DoFData<DoFHandler<deal_II_dimension,deal_II_dimension+1>,deal_II_dimension, deal_II_dimension+1>;
template class DataOut_DoFData<DoFHandler<deal_II_dimension,deal_II_dimension+1>,deal_II_dimension+1, deal_II_dimension+1>;

INSTANTIATE_VECTORS(DoFHandler,deal_II_dimension,deal_II_dimension+1,deal_II_dimension,deal_II_dimension+1);
INSTANTIATE_VECTORS(DoFHandler,deal_II_dimension,deal_II_dimension+1,deal_II_dimension+1,deal_II_dimension+1);

template class DataOut<deal_II_dimension, DoFHandler<deal_II_dimension,deal_II_dimension+1> >;
#endif


#undef INSTANTIATE
#undef INSTANTIATE_VECTORS

DEAL_II_NAMESPACE_CLOSE

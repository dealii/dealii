//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/quadrature_lib.h>
#include <base/work_stream.h>
#include <base/memory_consumption.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/petsc_vector.h>
#include <lac/petsc_block_vector.h>
#include <lac/trilinos_vector.h>
#include <lac/trilinos_block_vector.h>
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


namespace internal
{
  namespace DataOut
  {
    template <int dim, int spacedim>
    template <class FE>
    ParallelData<dim,spacedim>::
    ParallelData (const Quadrature<dim> &quadrature,
		  const unsigned int n_components,
		  const unsigned int n_datasets,
		  const unsigned int n_subdivisions,
		  const std::vector<unsigned int> &n_postprocessor_outputs,
		  const Mapping<dim,spacedim> &mapping,
		  const std::vector<std::vector<unsigned int> > &cell_to_patch_index_map,
		  const FE &finite_elements,
		  const UpdateFlags update_flags)
		    :
		    ParallelDataBase<dim,spacedim> (n_components,
						    n_datasets,
						    n_subdivisions,
						    quadrature.n_quadrature_points,
						    n_postprocessor_outputs,
						    finite_elements),
		    q_collection (quadrature),
		    mapping_collection (mapping),
		    x_fe_values (this->mapping_collection,
				 this->fe_collection,
				 q_collection,
				 update_flags),
		    cell_to_patch_index_map (&cell_to_patch_index_map)
    {}




				     /**
				      * In a WorkStream context, use
				      * this function to append the
				      * patch computed by the parallel
				      * stage to the array of patches.
				      */
    template <int dim, int spacedim>
    void
    append_patch_to_list (const DataOutBase::Patch<dim,spacedim> &patch,
			  std::vector<DataOutBase::Patch<dim,spacedim> > &patches)
    {
      patches.push_back (patch);
      patches.back().patch_index = patches.size()-1;
    }
  }
}


template <class DH, int patch_dim, int patch_space_dim>
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
DataEntryBase::DataEntryBase (const std::vector<std::string> &names_in,
			      const std::vector<DataComponentInterpretation::DataComponentInterpretation> &data_component_interpretation)
		:
		names(names_in),
		data_component_interpretation (data_component_interpretation),
		postprocessor(0, typeid(*this).name()),
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
		postprocessor(data_postprocessor, typeid(*this).name()),
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
    dof_data.push_back (std_cxx1x::shared_ptr<DataEntryBase>(new_entry));
  else
    cell_data.push_back (std_cxx1x::shared_ptr<DataEntryBase>(new_entry));
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
  dof_data.push_back (std_cxx1x::shared_ptr<DataEntryBase>(new_entry));
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
    typename std::vector<std_cxx1x::shared_ptr<DataEntryBase> >::const_iterator
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
std::vector<std_cxx1x::tuple<unsigned int, unsigned int, std::string> >
DataOut_DoFData<DH,patch_dim,patch_space_dim>::get_vector_data_ranges () const
{
  std::vector<std_cxx1x::tuple<unsigned int, unsigned int, std::string> >
    ranges;

				   // collect the ranges of dof
				   // and cell data
  typedef
    typename std::vector<std_cxx1x::shared_ptr<DataEntryBase> >::const_iterator
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
	  std_cxx1x::tuple<unsigned int, unsigned int, std::string>
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
void
DataOut<dim,DH>::
build_one_patch (const std::pair<cell_iterator, unsigned int> *cell_and_index,
		 internal::DataOut::ParallelData<DH::dimension, DH::space_dimension> &data,
		 DataOutBase::Patch<DH::dimension, DH::space_dimension> &patch,
		 const CurvedCellRegion curved_cell_region)
{
				   // use ucd_to_deal map as patch vertices
				   // are in the old, unnatural ordering. if
				   // the mapping does not preserve locations
				   // (e.g. MappingQEulerian), we need to
				   // compute the offset of the vertex for the
				   // graphical output. Otherwise, we can just
				   // use the vertex info.
  for (unsigned int vertex=0; vertex<GeometryInfo<DH::dimension>::vertices_per_cell; ++vertex)
    if (data.mapping_collection[0].preserves_vertex_locations())
      patch.vertices[vertex] = cell_and_index->first->vertex(vertex);
    else
      patch.vertices[vertex] = data.mapping_collection[0].transform_unit_to_real_cell
			       (cell_and_index->first,
				GeometryInfo<DH::dimension>::unit_cell_vertex (vertex));

  if (data.n_datasets > 0)
    {
      data.x_fe_values.reinit (cell_and_index->first);
      const FEValues<DH::dimension,DH::space_dimension> &fe_patch_values
	= data.x_fe_values.get_present_fe_values ();

      const unsigned int n_q_points = fe_patch_values.n_quadrature_points;

				       // depending on the requested output
				       // of curved cells, if necessary
				       // append the quadrature points to
				       // the last rows of the patch.data
				       // member. THis is the case if we
				       // want to produce curved cells at
				       // the boundary and this cell
				       // actually is at the boundary, or
				       // else if we want to produce curved
				       // cells everywhere
      if (curved_cell_region==curved_inner_cells ||
	  (curved_cell_region==curved_boundary && cell_and_index->first->at_boundary()))
	{
	  Assert(patch.space_dim==dim, ExcInternalError());
	  const std::vector<Point<DH::space_dimension> > & q_points=fe_patch_values.get_quadrature_points();
					   // resize the patch.data member
					   // in order to have enough memory
					   // for the quadrature points as
					   // well
	  patch.data.reinit (data.n_datasets+DH::space_dimension, n_q_points);
					   // set the flag indicating that
					   // for this cell the points are
					   // explicitly given
	  patch.points_are_available=true;
					   // copy points to patch.data
	  for (unsigned int i=0; i<DH::space_dimension; ++i)
	    for (unsigned int q=0; q<n_q_points; ++q)
	      patch.data(patch.data.size(0)-dim+i,q)=q_points[q][i];
	}
      else
	{
	  patch.data.reinit(data.n_datasets, n_q_points);
	  patch.points_are_available = false;
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
		  std::vector<Point<DH::space_dimension> > dummy_normals;
		  postprocessor->
		    compute_derived_quantities_scalar(data.patch_values,
						      data.patch_gradients,
						      data.patch_hessians,
						      dummy_normals,
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
		  std::vector<Point<DH::space_dimension> > dummy_normals;
		  postprocessor->
		    compute_derived_quantities_vector(data.patch_values_system,
						      data.patch_gradients_system,
						      data.patch_hessians_system,
						      dummy_normals,
						      data.postprocessed_values[dataset]);
		}

	      for (unsigned int q=0; q<n_q_points; ++q)
		for (unsigned int component=0;
		     component<this->dof_data[dataset]->n_output_variables;
		     ++component)
		  patch.data(offset+component,q)
		    = data.postprocessed_values[dataset][q](component);
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
		  patch.data(offset,q) = data.patch_values[q];
	      }
	    else
	      {
		this->dof_data[dataset]->get_function_values (fe_patch_values,
							      data.patch_values_system);
		for (unsigned int component=0; component<data.n_components;
		     ++component)
		  for (unsigned int q=0; q<n_q_points; ++q)
		    patch.data(offset+component,q) =
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
	  Assert (!cell_and_index->first->has_children(), ExcNotImplemented());

	  for (unsigned int dataset=0; dataset<this->cell_data.size(); ++dataset)
	    {
	      const double value
		= this->cell_data[dataset]->get_cell_data_value (cell_and_index->second);
	      for (unsigned int q=0; q<n_q_points; ++q)
		patch.data(offset+dataset,q) = value;
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
      if (cell_and_index->first->at_boundary(f)
	  ||
	  (cell_and_index->first->neighbor(f)->level() != cell_and_index->first->level()))
	{
	  patch.neighbors[f] = numbers::invalid_unsigned_int;
	  continue;
	}

      const cell_iterator neighbor = cell_and_index->first->neighbor(f);
      Assert (static_cast<unsigned int>(neighbor->level()) <
	      data.cell_to_patch_index_map->size(),
	      ExcInternalError());
      if ((static_cast<unsigned int>(neighbor->index()) >=
	   (*data.cell_to_patch_index_map)[neighbor->level()].size())
	  ||
	  ((*data.cell_to_patch_index_map)[neighbor->level()][neighbor->index()]
	   ==
	   dealii::DataOutBase::Patch<DH::dimension>::no_neighbor))
	{
	  patch.neighbors[f] = numbers::invalid_unsigned_int;
	  continue;
	}

				       // now, there is a
				       // neighbor, so get its
				       // patch number and set it
				       // for the neighbor index
      patch.neighbors[f]
	= (*data.cell_to_patch_index_map)[neighbor->level()][neighbor->index()];
    }
}



template <int dim, class DH>
void DataOut<dim,DH>::build_patches (const unsigned int n_subdivisions)
{
  build_patches (StaticMappingQ1<DH::dimension,DH::space_dimension>::mapping,
		 n_subdivisions, no_curved_cells);
}


template <int dim, class DH>
void DataOut<dim,DH>::build_patches (const Mapping<DH::dimension,DH::space_dimension> &mapping,
				     const unsigned int nnnn_subdivisions,
				     const CurvedCellRegion curved_region)
{
				   // Check consistency of redundant
				   // template parameter
  Assert (dim==DH::dimension, ExcDimensionMismatch(dim, DH::dimension));

  typedef DataOut_DoFData<DH, DH::dimension, DH::space_dimension> BaseClass;
  Assert (this->dofs != 0, typename BaseClass::ExcNoDoFHandlerSelected());

  const unsigned int n_subdivisions = (nnnn_subdivisions != 0)
				      ? nnnn_subdivisions
				      : this->default_subdivisions;
  Assert (n_subdivisions >= 1,
	  ExcInvalidNumberOfSubdivisions(n_subdivisions));

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
    }

  std::vector<std::pair<cell_iterator, unsigned int> > all_cells;
  {
				     // set the index of the first
				     // cell. if first_cell/next_cell
				     // returns non-active cells, then
				     // the index is not usable
				     // anyway, but otherwise we
				     // should keep track where we are
    unsigned int index;
    if (first_cell()->has_children())
      index = 0;
    else
      index = std::distance (this->dofs->begin_active(),
			     active_cell_iterator(first_cell()));
    for (cell_iterator cell=first_cell(); cell != this->dofs->end();
	 cell = next_cell(cell))
      {
	Assert (static_cast<unsigned int>(cell->level()) <
		cell_to_patch_index_map.size(),
		ExcInternalError());
	Assert (static_cast<unsigned int>(cell->index()) <
		cell_to_patch_index_map[cell->level()].size(),
		ExcInternalError());

	cell_to_patch_index_map[cell->level()][cell->index()] = all_cells.size();

	all_cells.push_back (std::make_pair(cell, index));

					 // if both this and the next
					 // cell are active, then
					 // increment the index that
					 // keeps track on which
					 // active cell we are sitting
					 // correctly. if one of the
					 // cells is not active, then
					 // this index doesn't mean
					 // anything anyway, so just
					 // ignore it. same if we are
					 // at the end of the range
	if (!cell->has_children() &&
	    next_cell(cell) != this->dofs->end() &&
	    !next_cell(cell)->has_children())
	  index += std::distance (active_cell_iterator(cell),
				  active_cell_iterator(next_cell(cell)));
      }
  }

  this->patches.clear ();
  this->patches.reserve (all_cells.size());
  Assert (this->patches.size() == 0, ExcInternalError());

				   // now create a default patch and a
				   // default object for the
				   // WorkStream object to work with
  const QTrapez<1>     q_trapez;
  const QIterated<DH::dimension> patch_points (q_trapez, n_subdivisions);

  const unsigned int n_components   = this->dofs->get_fe().n_components();
  unsigned int n_datasets=this->cell_data.size();
  for (unsigned int i=0; i<this->dof_data.size(); ++i)
    n_datasets += this->dof_data[i]->n_output_variables;

  std::vector<unsigned int> n_postprocessor_outputs (this->dof_data.size());
  for (unsigned int dataset=0; dataset<this->dof_data.size(); ++dataset)
    if (this->dof_data[dataset]->postprocessor)
      n_postprocessor_outputs[dataset] = this->dof_data[dataset]->n_output_variables;
    else
      n_postprocessor_outputs[dataset] = 0;

  const CurvedCellRegion curved_cell_region
    = (n_subdivisions<2 ? no_curved_cells : curved_region);

  UpdateFlags update_flags = update_values;
  if (curved_cell_region != no_curved_cells)
    update_flags |= update_quadrature_points;

  for (unsigned int i=0; i<this->dof_data.size(); ++i)
    if (this->dof_data[i]->postprocessor)
      update_flags |= this->dof_data[i]->postprocessor->get_needed_update_flags();
				   // perhaps update_normal_vectors is present,
				   // which would only be useful on faces, but
				   // we may not use it here.
  Assert (!(update_flags & update_normal_vectors),
	  ExcMessage("The update of normal vectors may not be requested for evaluation of "
		     "data on cells via DataPostprocessor."));

  internal::DataOut::ParallelData<DH::dimension, DH::space_dimension>
    thread_data (patch_points,
		 n_components, n_datasets, n_subdivisions,
		 n_postprocessor_outputs,
		 mapping,
		 cell_to_patch_index_map,
		 this->dofs->get_fe(),
		 update_flags);

  DataOutBase::Patch<DH::dimension, DH::space_dimension> sample_patch;
  sample_patch.n_subdivisions = n_subdivisions;
  sample_patch.data.reinit (n_datasets, patch_points.n_quadrature_points);



				   // now build the patches in parallel
  WorkStream::run (&all_cells[0],
		   &all_cells[0]+all_cells.size(),
		   std_cxx1x::bind(&DataOut<dim,DH>::build_one_patch,
				   *this, _1, _2, _3,
				   curved_cell_region),
		   std_cxx1x::bind(&internal::DataOut::append_patch_to_list<dim,DH::space_dimension>,
				   _1, std_cxx1x::ref(this->patches)),
		   thread_data,
		   sample_patch);
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
#include "data_out.inst"

DEAL_II_NAMESPACE_CLOSE

//----------------------------  data_out_stack.cc  ---------------------------
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
//----------------------------  data_out_stack.cc  ---------------------------


#include <numerics/data_out_stack.h>
#include <base/quadrature_lib.h>
#include <base/memory_consumption.h>
#include <lac/vector.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <fe/mapping_q1.h>



// if necessary try to work around a bug in the IBM xlC compiler
#ifdef XLC_WORK_AROUND_STD_BUG
using namespace std;
#endif


//TODO:[RH,GK] Replace global by local object; better: have two functions, or by default arg
static const MappingQ1<deal_II_dimension> mapping;


template <int dim>
unsigned int
DataOutStack<dim>::DataVector::memory_consumption () const
{
  return (MemoryConsumption::memory_consumption (data) +
	  MemoryConsumption::memory_consumption (names));
};



template <int dim>
DataOutStack<dim>::~DataOutStack ()
{};


template <int dim>
void DataOutStack<dim>::new_parameter_value (const double p,
					     const double dp)
{
  parameter      = p;
  parameter_step = dp;

				   // check whether the user called @p{finish_...}
				   // at the end of the previous parameter step
				   //
				   // this is to prevent serious waste of
				   // memory
  for (typename std::vector<DataVector>::const_iterator i=dof_data.begin();
       i!=dof_data.end(); ++i)
    Assert (i->data.size() == 0,
	    ExcDataNotCleared ());
  for (typename std::vector<DataVector>::const_iterator i=cell_data.begin();
       i!=cell_data.end(); ++i)
    Assert (i->data.size() == 0,
	    ExcDataNotCleared ());
  
};


template <int dim>
void DataOutStack<dim>::attach_dof_handler (const DoFHandler<dim> &dof) 
{
  dof_handler = &dof;
};


template <int dim>
void DataOutStack<dim>::declare_data_vector (const std::string &name,
					     const VectorType   vector_type)
{
  std::vector<std::string> names;
  names.push_back (name);
  declare_data_vector (names, vector_type);
};


template <int dim>
void DataOutStack<dim>::declare_data_vector (const std::vector<std::string> &names,
					     const VectorType    vector_type)
{
				   // make sure this function is
				   // not called after some parameter
				   // values have already been
				   // processed
  Assert (patches.size() == 0, ExcDataAlreadyAdded());

				   // also make sure that no name is
				   // used twice
  for (std::vector<std::string>::const_iterator name=names.begin(); name!=names.end(); ++name)
    {
      for (typename std::vector<DataVector>::const_iterator data_set=dof_data.begin();
	   data_set!=dof_data.end(); ++data_set)
	for (unsigned int i=0; i<data_set->names.size(); ++i)
	  Assert (*name != data_set->names[i], ExcNameAlreadyUsed(*name));

      for (typename std::vector<DataVector>::const_iterator data_set=cell_data.begin();
	   data_set!=cell_data.end(); ++data_set)
	for (unsigned int i=0; i<data_set->names.size(); ++i)
	  Assert (*name != data_set->names[i], ExcNameAlreadyUsed(*name));
    };
  
  switch (vector_type)
    {
      case dof_vector:
	    dof_data.push_back (DataVector());
	    dof_data.back().names = names;
	    break;

      case cell_vector:
	    cell_data.push_back (DataVector());
	    cell_data.back().names = names;
	    break;
    };
};


template <int dim>
template <typename number>
void DataOutStack<dim>::add_data_vector (const Vector<number> &vec,
					 const std::string    &name)
{
  std::vector<std::string> names;
  names.push_back (name);
  add_data_vector (vec, names);
};


template <int dim>
template <typename number>
void DataOutStack<dim>::add_data_vector (const Vector<number> &vec,
					 const std::vector<std::string> &names)
{
  Assert (dof_handler != 0, ExcNoDoFHandlerSelected ());
				   // either cell data and one name,
				   // or dof data and n_components names
  Assert (((vec.size() == dof_handler->get_tria().n_active_cells()) &&
	   (names.size() == 1))
	  ||
	  ((vec.size() == dof_handler->n_dofs()) &&
	   (names.size() == dof_handler->get_fe().n_components())),
	  ExcInvalidNumberOfNames (names.size(), dof_handler->get_fe().n_components()));
  for (unsigned int i=0; i<names.size(); ++i)
    Assert (names[i].find_first_not_of("abcdefghijklmnopqrstuvwxyz"
				       "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
				       "0123456789_<>()") == std::string::npos,
	    ExcInvalidCharacter (names[i]));
  
  if (vec.size() == dof_handler->n_dofs())
    {
      typename std::vector<DataVector>::iterator data_vector=dof_data.begin();
      for (; data_vector!=dof_data.end(); ++data_vector)
	if (data_vector->names == names)
	  {
	    data_vector->data.reinit (vec.size());
	    std::copy (vec.begin(), vec.end(),
		       data_vector->data.begin());
	    break;
	  };
      Assert (data_vector != dof_data.end(),
	      ExcVectorNotDeclared (names[0]));
    }
  else
    {
      typename std::vector<DataVector>::iterator data_vector=cell_data.begin();
      for (; data_vector!=cell_data.end(); ++data_vector)
	if (data_vector->names == names)
	  {
	    data_vector->data.reinit (vec.size());
	    std::copy (vec.begin(), vec.end(),
		       data_vector->data.begin());
	    break;
	  };
      Assert (data_vector != cell_data.end(),
	      ExcVectorNotDeclared (names[0]));
    };
};


template <int dim>
void DataOutStack<dim>::build_patches (const unsigned int n_subdivisions) 
{
				   // this is mostly copied from the
				   // DataOut class
  Assert (n_subdivisions >= 1,
	  ExcInvalidNumberOfSubdivisions(n_subdivisions));  
  Assert (dof_handler != 0, ExcNoDoFHandlerSelected());
  
  const unsigned int n_components   = dof_handler->get_fe().n_components();
  const unsigned int n_datasets     = dof_data.size() * n_components +
				      cell_data.size();
  
				   // first count the cells we want to
				   // create patches of and make sure
				   // there is enough memory for that
  unsigned int n_patches = 0;
  for (DoFHandler<dim>::active_cell_iterator cell=dof_handler->begin_active();
       cell != dof_handler->end(); ++cell)
    ++n_patches;


				   // before we start the loop:
				   // create a quadrature rule that
				   // actually has the points on this
				   // patch, and an object that
				   // extracts the data on each
				   // cell to these points
  QTrapez<1>     q_trapez;
  QIterated<dim> patch_points (q_trapez, n_subdivisions);
  FEValues<dim>  fe_patch_values (mapping, dof_handler->get_fe(),
				  patch_points,
				  update_values);
  const unsigned int n_q_points = patch_points.n_quadrature_points;
  std::vector<double>          patch_values (n_q_points);
  std::vector<Vector<double> > patch_values_system (n_q_points,
						    Vector<double>(n_components));

				   // add the required number of patches
  DataOutBase::Patch<dim+1>  default_patch;
  default_patch.n_subdivisions = n_subdivisions;
  default_patch.data.reinit (n_datasets, n_q_points*(n_subdivisions+1));
  patches.insert (patches.end(), n_patches, default_patch);

				   // now loop over all cells and
				   // actually create the patches
  typename std::vector<DataOutBase::Patch<dim+1> >::iterator
    patch = patches.begin() + (patches.size()-n_patches);
  unsigned int cell_number = 0;
  for (typename DoFHandler<dim>::active_cell_iterator cell=dof_handler->begin_active();
       cell != dof_handler->end(); ++cell, ++patch, ++cell_number)
    {
      Assert (patch != patches.end(), ExcInternalError());

				       // first fill in the vertices of the patch
      switch (dim)
	{
	  case 1:
		patch->vertices[0] = Point<dim+1>(cell->vertex(0)(0),
						  parameter-parameter_step);
		patch->vertices[3] = Point<dim+1>(cell->vertex(0)(0),
						  parameter);
		patch->vertices[1] = Point<dim+1>(cell->vertex(1)(0),
						  parameter-parameter_step);
		patch->vertices[2] = Point<dim+1>(cell->vertex(1)(0),
						  parameter);
		break;
		
	  default:
		Assert (false, ExcNotImplemented());
	};


				       // now fill in the the data values.
				       // note that the required order is
				       // with highest coordinate running
				       // fastest, we need to enter each
				       // value (n_subdivisions+1) times
				       // in succession
      if (n_datasets > 0)
	{
	  fe_patch_values.reinit (cell);
	  
					   // first fill dof_data
	  for (unsigned int dataset=0; dataset<dof_data.size(); ++dataset)
	    {
	      if (n_components == 1)
		{
		  fe_patch_values.get_function_values (dof_data[dataset].data,
						       patch_values);
		  for (unsigned int q=0; q<n_q_points; ++q)
		    for (unsigned int i=0; i<n_subdivisions+1; ++i)
		      patch->data(dataset,q*(n_subdivisions+1)+i) = patch_values[q];
		}
	      else
						 // system of components
		{
		  fe_patch_values.get_function_values (dof_data[dataset].data,
						       patch_values_system);
		  for (unsigned int component=0; component<n_components; ++component)
		    for (unsigned int q=0; q<n_q_points; ++q)
		      for (unsigned int i=0; i<n_subdivisions+1; ++i)
			patch->data(dataset*n_components+component,
				    q*(n_subdivisions+1)+i)
			  = patch_values_system[q](component);
		};
	    };

					   // then do the cell data
	  for (unsigned int dataset=0; dataset<cell_data.size(); ++dataset)
	    {
	      const double value = cell_data[dataset].data(cell_number);
	      for (unsigned int q=0; q<n_q_points; ++q)
		for (unsigned int i=0; i<n_subdivisions+1; ++i)
		  patch->data(dataset+dof_data.size()*n_components,
			      q*(n_subdivisions+1)+i) = value;
	    };
	};
    };
};


template <int dim>
void DataOutStack<dim>::finish_parameter_value ()
{
				   // release lock on dof handler
  dof_handler = 0;
  for (typename std::vector<DataVector>::iterator i=dof_data.begin();
       i!=dof_data.end(); ++i)
    i->data.reinit (0);
  
  for (typename std::vector<DataVector>::iterator i=cell_data.begin();
       i!=cell_data.end(); ++i)
    i->data.reinit (0);
};



template <int dim>
unsigned int
DataOutStack<dim>::memory_consumption () const
{
  return (DataOutInterface<dim+1>::memory_consumption () +
	  MemoryConsumption::memory_consumption (parameter) +
	  MemoryConsumption::memory_consumption (parameter_step) +
	  MemoryConsumption::memory_consumption (dof_handler) +
	  MemoryConsumption::memory_consumption (patches) +
	  MemoryConsumption::memory_consumption (dof_data) +
	  MemoryConsumption::memory_consumption (cell_data));
};



template <int dim>
const typename std::vector<typename DataOutBase::Patch<dim+1> > &
DataOutStack<dim>::get_patches () const
{
  return patches;
};



template <int dim>
std::vector<std::string> DataOutStack<dim>::get_dataset_names () const
{
  std::vector<std::string> names;
  for (typename std::vector<DataVector>::const_iterator dataset=dof_data.begin();
       dataset!=dof_data.end(); ++dataset)
    names.insert (names.end(), dataset->names.begin(), dataset->names.end());
  for (typename std::vector<DataVector>::const_iterator dataset=cell_data.begin();
       dataset!=cell_data.end(); ++dataset)
    names.insert (names.end(), dataset->names.begin(), dataset->names.end());
  
  return names;
};


// explicit instantiations
template class DataOutStack<deal_II_dimension>;
template void DataOutStack<deal_II_dimension>::add_data_vector (const Vector<double> &vec,
								const std::string    &name);
template void DataOutStack<deal_II_dimension>::add_data_vector (const Vector<float>  &vec,
								const std::string    &name);


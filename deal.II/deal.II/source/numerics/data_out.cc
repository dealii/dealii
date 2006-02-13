//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006 by the deal.II authors
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
#include <dofs/hp_dof_handler.h>
#include <grid/tria_iterator.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <fe/hp_fe_values.h>
#include <fe/mapping_q1.h>
#include <fe/select_fe_values.h>

#include <sstream>



template <int dof_handler_dim, template <int> class DH,
	  int patch_dim, int patch_space_dim>
DataOut_DoFData<dof_handler_dim,DH,patch_dim,patch_space_dim>::
DataEntryBase::DataEntryBase (const std::vector<std::string> &names)
		:
		names(names)
{}



template <int dof_handler_dim, template <int> class DH,
	  int patch_dim, int patch_space_dim>
DataOut_DoFData<dof_handler_dim,DH,patch_dim,patch_space_dim>::
DataEntryBase::~DataEntryBase ()
{}



template <int dof_handler_dim,  template <int> class DH,
	  int patch_dim, int patch_space_dim>
template <typename VectorType>
DataOut_DoFData<dof_handler_dim,DH,patch_dim,patch_space_dim>::
DataEntry<VectorType>::
DataEntry (const VectorType               *data,
	   const std::vector<std::string> &names)
		:
                DataEntryBase (names),
		vector (data)
{}


template <int dof_handler_dim, template <int> class DH,
	  int patch_dim, int patch_space_dim>
template <typename VectorType>
double
DataOut_DoFData<dof_handler_dim,DH,patch_dim,patch_space_dim>::
DataEntry<VectorType>::
get_cell_data_value (const unsigned int cell_number) const
{
  return (*vector)(cell_number);
}



template <int dof_handler_dim, template <int> class DH,
	  int patch_dim, int patch_space_dim>
template <typename VectorType>
void
DataOut_DoFData<dof_handler_dim,DH,patch_dim,patch_space_dim>::
DataEntry<VectorType>::
get_function_values (const FEValuesBase<dof_handler_dim> &fe_patch_values,
                     std::vector<Vector<double> >    &patch_values_system) const
{
  fe_patch_values.get_function_values (*vector, patch_values_system);
}



template <int dof_handler_dim, template <int> class DH,
	  int patch_dim, int patch_space_dim>
template <typename VectorType>
void
DataOut_DoFData<dof_handler_dim,DH,patch_dim,patch_space_dim>::
DataEntry<VectorType>::
get_function_values (const FEValuesBase<dof_handler_dim> &fe_patch_values,
                     std::vector<double>             &patch_values) const
{
  fe_patch_values.get_function_values (*vector, patch_values);
}



template <int dof_handler_dim, template <int> class DH,
	  int patch_dim, int patch_space_dim>
template <typename VectorType>
unsigned int
DataOut_DoFData<dof_handler_dim,DH,patch_dim,patch_space_dim>::
DataEntry<VectorType>::memory_consumption () const
{
  return (sizeof (vector) +
	  MemoryConsumption::memory_consumption (this->names));
}



template <int dof_handler_dim,  template <int> class DH,
	  int patch_dim, int patch_space_dim>
template <typename VectorType>
void
DataOut_DoFData<dof_handler_dim,DH,patch_dim,patch_space_dim>::
DataEntry<VectorType>::clear ()
{
  vector = 0;
}




template <int dof_handler_dim, template <int> class DH,
	  int patch_dim, int patch_space_dim>
DataOut_DoFData<dof_handler_dim,DH,patch_dim,patch_space_dim>::DataOut_DoFData ()
                :
		dofs(0,typeid(*this).name())
{}



template <int dof_handler_dim, template <int> class DH,
	  int patch_dim, int patch_space_dim>
DataOut_DoFData<dof_handler_dim,DH,patch_dim,patch_space_dim>::~DataOut_DoFData ()
{
  clear ();
}



template <int dof_handler_dim, template <int> class DH,
	  int patch_dim, int patch_space_dim>
void
DataOut_DoFData<dof_handler_dim,DH,patch_dim,patch_space_dim>::
attach_dof_handler (const DH<dof_handler_dim> &d)
{
  Assert (dof_data.size() == 0, ExcOldDataStillPresent());
  Assert (cell_data.size() == 0, ExcOldDataStillPresent());
  
  dofs = &d;
}



template <int dof_handler_dim, template <int> class DH,
	  int patch_dim, int patch_space_dim>
template <class VECTOR>
void
DataOut_DoFData<dof_handler_dim,DH,patch_dim,patch_space_dim>::
add_data_vector (const VECTOR                   &vec,
		 const std::vector<std::string> &names,
		 const DataVectorType            type)
{
  Assert (dofs != 0, ExcNoDoFHandlerSelected ());

				   // check for allowed sizes of
				   // vectors
  Assert ((vec.size() == dofs->get_tria().n_active_cells()) ||
	  (vec.size() == dofs->n_dofs()),
	  ExcInvalidVectorSize(vec.size(),
			       dofs->n_dofs(),
			       dofs->get_tria().n_active_cells()));
  
				   // either cell data and one name,
				   // or dof data and n_components names
  DataVectorType actual_type = type;
  if (type == type_automatic)
    {
      if (vec.size() == dofs->get_tria().n_active_cells())
	actual_type = type_cell_data;
      else
	actual_type = type_dof_data;
    };
  
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
    };
  
  for (unsigned int i=0; i<names.size(); ++i)
    Assert (names[i].find_first_not_of("abcdefghijklmnopqrstuvwxyz"
				       "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
				       "0123456789_<>()") == std::string::npos,
	    ExcInvalidCharacter (names[i],
				 names[i].find_first_not_of("abcdefghijklmnopqrstuvwxyz"
							    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
							    "0123456789_<>()")));
  
  DataEntryBase * new_entry = new DataEntry<VECTOR>(&vec, names);
  if (actual_type == type_dof_data)
    dof_data.push_back (boost::shared_ptr<DataEntryBase>(new_entry));
  else
    cell_data.push_back (boost::shared_ptr<DataEntryBase>(new_entry));
}



template <int dof_handler_dim, template <int> class DH,
	  int patch_dim, int patch_space_dim>
template <class VECTOR>
void
DataOut_DoFData<dof_handler_dim,DH,patch_dim,patch_space_dim>::
add_data_vector (const VECTOR         &vec,
		 const std::string    &name,
		 const DataVectorType  type)
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
  
  add_data_vector (vec, names, type);
}



template <int dof_handler_dim, template <int> class DH,
	  int patch_dim, int patch_space_dim>
void DataOut_DoFData<dof_handler_dim,DH,patch_dim,patch_space_dim>::clear_data_vectors ()
{
  dof_data.erase (dof_data.begin(), dof_data.end());
  cell_data.erase (cell_data.begin(), cell_data.end());

				   // delete patches
  std::vector<Patch> dummy;
  patches.swap (dummy);
}



template <int dof_handler_dim, template <int> class DH,
	  int patch_dim, int patch_space_dim>
void
DataOut_DoFData<dof_handler_dim,DH,patch_dim,patch_space_dim>::
clear_input_data_references ()
{
  for (unsigned int i=0; i<dof_data.size(); ++i)
    dof_data[i]->clear ();
  
  for (unsigned int i=0; i<cell_data.size(); ++i)
    cell_data[i]->clear ();

  if (dofs != 0)
    dofs = 0;
}



template <int dof_handler_dim, template <int> class DH,
          int patch_dim, int patch_space_dim>
void
DataOut_DoFData<dof_handler_dim,DH,patch_dim,patch_space_dim>::clear ()
{
  dof_data.erase (dof_data.begin(), dof_data.end());
  cell_data.erase (cell_data.begin(), cell_data.end());

  if (dofs != 0)
    dofs = 0;

				   // delete patches
  std::vector<Patch> dummy;
  patches.swap (dummy);
}



template <int dof_handler_dim, template <int> class DH,
	  int patch_dim, int patch_space_dim>
std::vector<std::string>
DataOut_DoFData<dof_handler_dim,DH,patch_dim,patch_space_dim>::
get_dataset_names () const 
{
  std::vector<std::string> names;
				   // collect the names of dof
				   // and cell data
  typedef
    typename std::vector<boost::shared_ptr<DataEntryBase> >::const_iterator
    data_iterator;
  
  for (data_iterator  d=dof_data.begin();
       d!=dof_data.end(); ++d)
    for (unsigned int i=0; i<(*d)->names.size(); ++i)
      names.push_back ((*d)->names[i]);
  for (data_iterator d=cell_data.begin(); d!=cell_data.end(); ++d)
    {
      Assert ((*d)->names.size() == 1, ExcInternalError());
      names.push_back ((*d)->names[0]);
    };

  return names;
}



template <int dof_handler_dim, template <int> class DH,
	  int patch_dim, int patch_space_dim>
const std::vector< ::DataOutBase::Patch<patch_dim, patch_space_dim> > &
DataOut_DoFData<dof_handler_dim,DH,patch_dim,patch_space_dim>::get_patches () const
{
  return patches;
}



template <int dof_handler_dim, template <int> class DH,
	  int patch_dim, int patch_space_dim>
unsigned int
DataOut_DoFData<dof_handler_dim,DH,patch_dim,patch_space_dim>::memory_consumption () const
{
  return (DataOutInterface<patch_dim,patch_space_dim>::memory_consumption () +
	  MemoryConsumption::memory_consumption (dofs) +
	  MemoryConsumption::memory_consumption (patches));
}



/* ---------------------------------------------------------------------- */



template <int dim, template <int> class DH>
void DataOut<dim,DH>::build_some_patches (Data data)
{
  QTrapez<1>     q_trapez;
  QIterated<dim> patch_points (q_trapez, data.n_subdivisions);

//TODO[RH]: This is strange -- Data has a member 'mapping' that should
//be used here, but it isn't. Rather, up until version 1.94, we were
//actually initializing a local mapping object and used that...
  typename SelectFEValues<DH<dim> >::FEValues
    x_fe_patch_values (this->dofs->get_fe(),
                       patch_points, update_values);

  const unsigned int n_q_points = patch_points.n_quadrature_points;
  
  typename std::vector< ::DataOutBase::Patch<dim> >::iterator patch = this->patches.begin();
  cell_iterator cell=first_cell();
  
				   // get first cell in this thread
  for (unsigned int i=0; (i<data.this_thread)&&(cell != this->dofs->end()); ++i)
    {
      ++patch;
      cell = next_cell(cell);
    }

  				   // now loop over all cells and
				   // actually create the patches
  for (;cell != this->dofs->end();)
    {
      Assert (patch != this->patches.end(), ExcInternalError());

				       // use ucd_to_deal map as patch
				       // vertices are in the old,
				       // unnatural ordering
      for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
	patch->vertices[vertex] = data.mapping->transform_unit_to_real_cell
				  (cell, GeometryInfo<dim>::unit_cell_vertex (vertex));
      
      if (data.n_datasets > 0)
	{
          x_fe_patch_values.reinit (cell);
          const FEValues<dim> &fe_patch_values
            = x_fe_patch_values.get_present_fe_values ();
	  
					   // first fill dof_data
	  for (unsigned int dataset=0; dataset<this->dof_data.size(); ++dataset)
	    {
	      if (data.n_components == 1)
		{
                  this->dof_data[dataset]->get_function_values (fe_patch_values,
                                                                data.patch_values);

		  for (unsigned int q=0; q<n_q_points; ++q)
		    patch->data(dataset,q) = data.patch_values[q];
		}
	      else
						 // system of components
		{
                  this->dof_data[dataset]->get_function_values (fe_patch_values,
                                                                data.patch_values_system);

		  for (unsigned int component=0; component<data.n_components;
		       ++component)
		    for (unsigned int q=0; q<n_q_points; ++q)
		      patch->data(dataset*data.n_components+component,q) =
			data.patch_values_system[q](component);
		};
	    };

					   // then do the cell data. only
					   // compute the number of a cell if
					   // needed; also make sure that we
					   // only access cell data if the
					   // first_cell/next_cell functions
					   // only return active cells
          if (this->cell_data.size() != 0)
            {
              Assert (!cell->has_children(), ExcNotImplemented());
//TODO: This introduces a quadratic component into the algorithm. We may be doing better if we always kept track of the cell_number when 'cell' is increased, but we have to make sure that we only do this if cell_number is actually used, since we don't want to enforce the activeness constraints otherwise
              const unsigned int cell_number
                = std::distance (this->dofs->begin_active(),
                                 active_cell_iterator (cell));
              
              for (unsigned int dataset=0; dataset<this->cell_data.size(); ++dataset)
                {
                  const double value
                    = this->cell_data[dataset]->get_cell_data_value (cell_number);
                  for (unsigned int q=0; q<n_q_points; ++q)
                    patch->data(dataset+this->dof_data.size()*data.n_components,q) =
                      value;
                }
            }
	}


      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
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
               ::DataOutBase::Patch<dim>::no_neighbor))
            continue;

                                           // now, there is a
                                           // neighbor, so get its
                                           // patch number and set it
                                           // for the neighbor index
          patch->neighbors[f] = this->patches[(*data.cell_to_patch_index_map)
					     [neighbor->level()][neighbor->index()]].patch_index;
        };
      
      				       // next cell (patch) in this
      				       // thread
      for (unsigned int i=0;
	   (i<data.n_threads)&&(cell != this->dofs->end()); ++i)
	{
	  ++patch;
          cell = next_cell(cell);
	}
    }
}



template <int dim, template <int> class DH>
void DataOut<dim,DH>::build_patches (const unsigned int n_subdivisions,
				     const unsigned int n_threads_) 
{
  static const MappingQ1<dim> mapping;
  build_patches (mapping, n_subdivisions, n_threads_);
}


template <int dim, template <int> class DH>
void DataOut<dim,DH>::build_patches (const Mapping<dim> &mapping,
				     const unsigned int n_subdivisions,
				     const unsigned int n_threads_) 
{
  Assert (n_subdivisions >= 1,
	  ExcInvalidNumberOfSubdivisions(n_subdivisions));

  typedef DataOut_DoFData<dim,DH,dim> BaseClass;
  Assert (this->dofs != 0, typename BaseClass::ExcNoDoFHandlerSelected());

  const unsigned int n_threads = (DEAL_II_USE_MT ? n_threads_ : 1);

 				   // before we start the loop:
 				   // create a quadrature rule that
 				   // actually has the points on this
 				   // patch
  QTrapez<1>     q_trapez;
  QIterated<dim> patch_points (q_trapez, n_subdivisions);
  
  const unsigned int n_q_points     = patch_points.n_quadrature_points;
  const unsigned int n_components   = this->dofs->get_fe().n_components();
  const unsigned int n_datasets     = this->dof_data.size() * n_components +
				      this->cell_data.size();
  
				   // clear the patches array
  if (true)
    {
      std::vector< ::DataOutBase::Patch<dim> > dummy;
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
                                         ::DataOutBase::Patch<dim>::no_neighbor);
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
  ::DataOutBase::Patch<dim>  default_patch;
  default_patch.n_subdivisions = n_subdivisions;
  default_patch.data.reinit (n_datasets, n_q_points);
  this->patches.insert (this->patches.end(), n_patches, default_patch);

  for (unsigned int i=0; i<this->patches.size(); ++i)
    this->patches[i].patch_index = i;
  

				   // init data for the threads    
  std::vector<Data> thread_data(n_threads);
  for (unsigned int i=0;i<n_threads;++i)
    {
      thread_data[i].n_threads      = n_threads;
      thread_data[i].this_thread    = i;
      thread_data[i].n_components   = n_components;
      thread_data[i].n_datasets     = n_datasets;
      thread_data[i].n_subdivisions = n_subdivisions;
      thread_data[i].patch_values.resize (n_q_points);
      thread_data[i].patch_values_system.resize (n_q_points);
      thread_data[i].mapping        = &mapping;

      thread_data[i].cell_to_patch_index_map = &cell_to_patch_index_map;
      
      for (unsigned int k=0; k<n_q_points; ++k)
	thread_data[i].patch_values_system[k].reinit(n_components);
    }

  Threads::ThreadGroup<> threads;  
  for (unsigned int l=0; l<n_threads; ++l)
    threads += Threads::spawn (*this, &DataOut<dim,DH>::build_some_patches)(thread_data[l]);
  threads.join_all();
}



template <int dim, template <int> class DH>
typename DataOut<dim,DH>::cell_iterator
DataOut<dim,DH>::first_cell () 
{
  return this->dofs->begin_active ();
}



template <int dim, template <int> class DH>
typename DataOut<dim,DH>::cell_iterator
DataOut<dim,DH>::next_cell (const typename DataOut<dim,DH>::cell_iterator &cell) 
{
				   // convert the iterator to an
				   // active_iterator and advance
				   // this to the next active cell
  typename DH<dim>::active_cell_iterator active_cell = cell;
  ++active_cell;
  return active_cell;
}


// explicit instantiations 

#define INSTANTIATE(D1,DH,D2,D3,V) \
template void \
DataOut_DoFData<D1,DH,D2,D3>:: \
add_data_vector<V > (const V             &, \
                     const std::string   &, \
                     const DataVectorType); \
\
template void \
DataOut_DoFData<D1,DH,D2,D3>:: \
add_data_vector<V > (const V                        &, \
                     const std::vector<std::string> &, \
                     const DataVectorType)

#ifndef DEAL_II_USE_PETSC
# define INSTANTIATE_VECTORS(D1,DH,D2,D3) \
    INSTANTIATE(D1,DH,D2,D3,Vector<double>); \
    INSTANTIATE(D1,DH,D2,D3,Vector<float>); \
    INSTANTIATE(D1,DH,D2,D3,BlockVector<double>) ; \
    INSTANTIATE(D1,DH,D2,D3,BlockVector<float>)
#else
# define INSTANTIATE_VECTORS(D1,DH,D2,D3) \
    INSTANTIATE(D1,DH,D2,D3,Vector<double>); \
    INSTANTIATE(D1,DH,D2,D3,Vector<float>); \
    INSTANTIATE(D1,DH,D2,D3,BlockVector<double>) ; \
    INSTANTIATE(D1,DH,D2,D3,BlockVector<float>); \
    INSTANTIATE(D1,DH,D2,D3,PETScWrappers::Vector); \
    INSTANTIATE(D1,DH,D2,D3,PETScWrappers::BlockVector)
#endif

// now do actual instantiations, first for DoFHandler...
template class DataOut_DoFData<deal_II_dimension,DoFHandler,deal_II_dimension>;
template class DataOut_DoFData<deal_II_dimension,DoFHandler,deal_II_dimension+1>;
INSTANTIATE_VECTORS(deal_II_dimension,DoFHandler,deal_II_dimension,deal_II_dimension);
INSTANTIATE_VECTORS(deal_II_dimension,DoFHandler,deal_II_dimension+1,deal_II_dimension+1);

// for 3d, also generate an extra class
#if deal_II_dimension >= 2
template class DataOut_DoFData<deal_II_dimension,DoFHandler,
			       deal_II_dimension-1,deal_II_dimension>;
INSTANTIATE_VECTORS(deal_II_dimension,DoFHandler,deal_II_dimension-1,deal_II_dimension);
#endif

template class DataOut<deal_II_dimension,DoFHandler>;


// ...and now for hp::DoFHandler
template class DataOut_DoFData<deal_II_dimension,hp::DoFHandler,deal_II_dimension>;
template class DataOut_DoFData<deal_II_dimension,hp::DoFHandler,deal_II_dimension+1>;
INSTANTIATE_VECTORS(deal_II_dimension,hp::DoFHandler,deal_II_dimension,deal_II_dimension);
INSTANTIATE_VECTORS(deal_II_dimension,hp::DoFHandler,deal_II_dimension+1,deal_II_dimension+1);

#if deal_II_dimension >= 2
template class DataOut_DoFData<deal_II_dimension,hp::DoFHandler,
			       deal_II_dimension-1,deal_II_dimension>;
INSTANTIATE_VECTORS(deal_II_dimension,hp::DoFHandler,deal_II_dimension-1,deal_II_dimension);
#endif

template class DataOut<deal_II_dimension,hp::DoFHandler>;

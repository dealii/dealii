//----------------------------  data_out.cc  ---------------------------
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
//----------------------------  data_out.cc  ---------------------------


#include <base/quadrature_lib.h>
#include <base/memory_consumption.h>
#include <lac/vector.h>
#include <numerics/data_out.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <fe/fe.h>
#include <fe/fe_values.h>

#ifdef DEAL_II_USE_MT
#include <base/thread_management.h>
#endif

#include <strstream>


template <int dof_handler_dim, int patch_dim, int patch_space_dim>
DataOut_DoFData<dof_handler_dim,patch_dim,patch_space_dim>::DataEntry::
DataEntry (const Vector<double> *data,
	   const std::vector<std::string> &names) :
		data(data),
		names(names)
{};



template <int dof_handler_dim, int patch_dim, int patch_space_dim>
unsigned int
DataOut_DoFData<dof_handler_dim,patch_dim,patch_space_dim>::
DataEntry::memory_consumption () const
{
  return (sizeof (data) +
	  MemoryConsumption::memory_consumption (names) +
	  MemoryConsumption::memory_consumption (units));
};





template <int dof_handler_dim, int patch_dim, int patch_space_dim>
DataOut_DoFData<dof_handler_dim,patch_dim,patch_space_dim>::DataOut_DoFData () :
		dofs(0)
{};



template <int dof_handler_dim, int patch_dim, int patch_space_dim>
DataOut_DoFData<dof_handler_dim,patch_dim,patch_space_dim>::~DataOut_DoFData ()
{
  clear ();
};



template <int dof_handler_dim, int patch_dim, int patch_space_dim>
void
DataOut_DoFData<dof_handler_dim,patch_dim,patch_space_dim>::
attach_dof_handler (const DoFHandler<dof_handler_dim> &d)
{
  Assert (dof_data.size() == 0, ExcOldDataStillPresent());
  Assert (cell_data.size() == 0, ExcOldDataStillPresent());
  
  if (dofs != 0)
    dofs->unsubscribe ();
  
  dofs = &d;
  if (dofs != 0)
    dofs->subscribe ();
};



template <int dof_handler_dim, int patch_dim, int patch_space_dim>
void
DataOut_DoFData<dof_handler_dim,patch_dim,patch_space_dim>::add_data_vector (const Vector<double> &vec,
							     const std::vector<std::string> &names)
{
  Assert (dofs != 0, ExcNoDoFHandlerSelected ());

				   // check for allowed sizes of
				   // vectors
  Assert ((vec.size() == dofs->get_tria().n_active_cells()) ||
	  (vec.size() == dofs->n_dofs()),
	  ExcInvalidVectorSize(vec.size(), dofs->get_tria().n_active_cells(), dofs->n_dofs()));
  
				   // either cell data and one name,
				   // or dof data and n_components names
  if (vec.size() == dofs->get_tria().n_active_cells())
    Assert (names.size() == 1,
	    ExcInvalidNumberOfNames (names.size(), 1));
  if (vec.size() == dofs->n_dofs())
    Assert (names.size() == dofs->get_fe().n_components(),
	    ExcInvalidNumberOfNames (names.size(), dofs->get_fe().n_components()));
  for (unsigned int i=0; i<names.size(); ++i)
    Assert (names[i].find_first_not_of("abcdefghijklmnopqrstuvwxyz"
				       "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
				       "0123456789_<>()") == std::string::npos,
	    ExcInvalidCharacter (names[i]));
  
  DataEntry new_entry (&vec, names);
  if (vec.size() == dofs->n_dofs())
    dof_data.push_back (new_entry);
  else
    if (vec.size() == dofs->get_tria().n_active_cells())
      cell_data.push_back (new_entry);
    else
      Assert (false,
	      ExcInvalidVectorSize (vec.size(),
				    dofs->n_dofs(),
				    dofs->get_tria().n_active_cells()));
};


template <int dof_handler_dim, int patch_dim, int patch_space_dim>
void
DataOut_DoFData<dof_handler_dim,patch_dim,patch_space_dim>::add_data_vector (const Vector<double> &vec,
							     const std::string         &name)
{
  unsigned int n_components = dofs->get_fe().n_components ();

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
	  std::ostrstream namebuf;
	  namebuf << name << '_' << i << std::ends;
  	  names[i] = namebuf.str();
  	};
    };
  
  add_data_vector (vec, names);
}



template <int dof_handler_dim, int patch_dim, int patch_space_dim>
void DataOut_DoFData<dof_handler_dim,patch_dim,patch_space_dim>::clear_data_vectors ()
{
  dof_data.erase (dof_data.begin(), dof_data.end());
  cell_data.erase (cell_data.begin(), cell_data.end());

				   // delete patches
  std::vector<DataOutBase::Patch<patch_dim,patch_space_dim> > dummy;
  patches.swap (dummy);
}



template <int dof_handler_dim, int patch_dim, int patch_space_dim>
void
DataOut_DoFData<dof_handler_dim,patch_dim,patch_space_dim>::clear_input_data_references ()
{
  for (unsigned int i=0; i<dof_data.size(); ++i)
    dof_data[i].data = 0;
  
  for (unsigned int i=0; i<cell_data.size(); ++i)
    cell_data[i].data = 0;

  if (dofs != 0)
    {
      dofs->unsubscribe ();
      dofs = 0;
    };
};



template <int dof_handler_dim, int patch_dim, int patch_space_dim>
void
DataOut_DoFData<dof_handler_dim,patch_dim,patch_space_dim>::clear ()
{
  dof_data.erase (dof_data.begin(), dof_data.end());
  cell_data.erase (cell_data.begin(), cell_data.end());

  if (dofs != 0)
    {
      dofs->unsubscribe ();
      dofs = 0;
    };

				   // delete patches
  std::vector<DataOutBase::Patch<patch_dim,patch_space_dim> > dummy;
  patches.swap (dummy);
}


template <int dof_handler_dim, int patch_dim, int patch_space_dim>
std::vector<std::string>
DataOut_DoFData<dof_handler_dim,patch_dim,patch_space_dim>::get_dataset_names () const 
{
  std::vector<std::string> names;
				   // collect the names of dof
				   // and cell data
  for (typename std::vector<DataEntry>::const_iterator d=dof_data.begin();
       d!=dof_data.end(); ++d)
    for (unsigned int i=0; i<d->names.size(); ++i)
      names.push_back (d->names[i]);
  for (typename std::vector<DataEntry>::const_iterator d=cell_data.begin();
       d!=cell_data.end(); ++d)
    {
      Assert (d->names.size() == 1, ExcInternalError());
      names.push_back (d->names[0]);
    };

  return names;
};



template <int dof_handler_dim, int patch_dim, int patch_space_dim>
const std::vector<typename DataOutBase::Patch<patch_dim, patch_space_dim> > &
DataOut_DoFData<dof_handler_dim,patch_dim,patch_space_dim>::get_patches () const
{
  return patches;
};



template <int dof_handler_dim, int patch_dim, int patch_space_dim>
unsigned int
DataOut_DoFData<dof_handler_dim,patch_dim,patch_space_dim>::memory_consumption () const
{
  return (DataOutInterface<patch_dim,patch_space_dim>::memory_consumption () +
	  MemoryConsumption::memory_consumption (dofs) +
	  MemoryConsumption::memory_consumption (dof_data) +
	  MemoryConsumption::memory_consumption (cell_data) +
	  MemoryConsumption::memory_consumption (patches));
};



/* ---------------------------------------------------------------------- */



template <int dim>
void DataOut<dim>::build_some_patches (Data data)
{
  QTrapez<1>     q_trapez;
  QIterated<dim> patch_points (q_trapez, data.n_subdivisions);
  
  FEValues<dim> fe_patch_values(dofs->get_fe(),
				patch_points,
				update_values);

  const unsigned int n_q_points = patch_points.n_quadrature_points;
  
  unsigned int cell_number = 0;
  typename std::vector<DataOutBase::Patch<dim> >::iterator patch = patches.begin();
  DoFHandler<dim>::cell_iterator cell=first_cell();

				   // get first cell in this thread
  for (unsigned int i=0; (i<data.this_thread)&&(cell != dofs->end()); ++i)
    {
      ++patch;
      ++cell_number;
      cell=next_cell(cell);
    }

  				   // now loop over all cells and
				   // actually create the patches
  for (;cell != dofs->end();)
    {
      Assert (patch != patches.end(), ExcInternalError());
      
      for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
	patch->vertices[vertex] = cell->vertex(vertex);
      
      if (data.n_datasets > 0)
	{
	  fe_patch_values.reinit (cell);
	  
					   // first fill dof_data
	  for (unsigned int dataset=0; dataset<dof_data.size(); ++dataset)
	    {
	      if (data.n_components == 1)
		{
		  fe_patch_values.get_function_values (*dof_data[dataset].data,
						       data.patch_values);
		  for (unsigned int q=0; q<n_q_points; ++q)
		    patch->data(dataset,q) = data.patch_values[q];
		}
	      else
						 // system of components
		{
		  fe_patch_values.get_function_values (*dof_data[dataset].data,
						       data.patch_values_system);
		  for (unsigned int component=0; component<data.n_components;
		       ++component)
		    for (unsigned int q=0; q<n_q_points; ++q)
		      patch->data(dataset*data.n_components+component,q) =
			data.patch_values_system[q](component);
		};
	    };

					   // then do the cell data
	  for (unsigned int dataset=0; dataset<cell_data.size(); ++dataset)
	    {
	      const double value = (*cell_data[dataset].data)(cell_number);
	      for (unsigned int q=0; q<n_q_points; ++q)
		patch->data(dataset+dof_data.size()*data.n_components,q) =
		  value;
	    };
	};
      				       // next cell (patch) in this thread
      for (unsigned int i=0;
	   (i<data.n_threads)&&(cell != dofs->end()); ++i)
	{
	  ++patch;
	  ++cell_number;
	  cell=next_cell(cell);
	}
    };
}


template <int dim>
void DataOut<dim>::build_patches (const unsigned int n_subdivisions,
				  const unsigned int n_threads_) 
{
  Assert (n_subdivisions >= 1,
	  ExcInvalidNumberOfSubdivisions(n_subdivisions));

  typedef DataOut_DoFData<dim,dim> BaseClass;
  Assert (dofs != 0, typename BaseClass::ExcNoDoFHandlerSelected());

#ifdef DEAL_II_USE_MT
  const unsigned int n_threads = n_threads_;
#else
				   // access this variable to avoid
				   // compiler warning about unused
				   // var:
  (void)n_threads_;
  const unsigned int n_threads = 1;
#endif


				   // before we start the loop:
				   // create a quadrature rule that
				   // actually has the points on this
				   // patch
  QTrapez<1>     q_trapez;
  QIterated<dim> patch_points (q_trapez, n_subdivisions);

  const unsigned int n_q_points     = patch_points.n_quadrature_points;
  const unsigned int n_components   = dofs->get_fe().n_components();
  const unsigned int n_datasets     = dof_data.size() * n_components +
				      cell_data.size();
  
				   // clear the patches array
  if (true)
    {
      std::vector<DataOutBase::Patch<dim> > dummy;
      patches.swap (dummy);
    };
  
				   // first count the cells we want to
				   // create patches of and make sure
				   // there is enough memory for that
  unsigned int n_patches = 0;
  for (DoFHandler<dim>::cell_iterator cell=first_cell();
       cell != dofs->end();
       cell = next_cell(cell))
    ++n_patches;

  std::vector<Data> thread_data(n_threads);

				   // init data for the threads
  for (unsigned int i=0;i<n_threads;++i)
    {
      thread_data[i].n_threads      = n_threads;
      thread_data[i].this_thread    = i;
      thread_data[i].n_components   = n_components;
      thread_data[i].n_datasets     = n_datasets;
      thread_data[i].n_subdivisions = n_subdivisions;
      thread_data[i].patch_values.resize (n_q_points);
      thread_data[i].patch_values_system.resize (n_q_points);
      
      for (unsigned int k=0; k<n_q_points; ++k)
	thread_data[i].patch_values_system[k].reinit(n_components);
    }

				   // create the patches with
				   // default values
  DataOutBase::Patch<dim>  default_patch;
  default_patch.n_subdivisions = n_subdivisions;
  default_patch.data.reinit (n_datasets, n_q_points);
  patches.insert (patches.end(), n_patches, default_patch);

#ifdef DEAL_II_USE_MT

  Threads::ThreadManager thread_manager;  
  for (unsigned int l=0;l<n_threads;++l)
    Threads::spawn (thread_manager,
		    Threads::encapsulate (&DataOut<dim>::build_some_patches)
		    .collect_args (this, thread_data[l]));
  thread_manager.wait();
  
				   // just one thread
#else
  build_some_patches(thread_data[0]);
#endif
};


template <int dim>
typename DoFHandler<dim>::cell_iterator
DataOut<dim>::first_cell () 
{
  return dofs->begin_active ();
};


template <int dim>
typename DoFHandler<dim>::cell_iterator
DataOut<dim>::next_cell (const typename DoFHandler<dim>::cell_iterator &cell) 
{
				   // convert the iterator to an
				   // active_iterator and advance
				   // this to the next active cell
  typename DoFHandler<dim>::active_cell_iterator active_cell = cell;
  ++active_cell;
  return active_cell;
};


// explicit instantiations
template class DataOut_DoFData<deal_II_dimension,deal_II_dimension>;
template class DataOut_DoFData<deal_II_dimension,deal_II_dimension+1>;
template class DataOut<deal_II_dimension>;


// for 3d, also generate an extra class
#if deal_II_dimension >= 3
template class DataOut_DoFData<deal_II_dimension,deal_II_dimension-1,deal_II_dimension>;
#endif

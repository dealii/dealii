/* $Id$ */


#include <base/quadrature_lib.h>
#include <lac/vector.h>
#include <numerics/data_out.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <fe/fe.h>
#include <fe/fe_values.h>

#ifdef DEAL_II_USE_MT
#include <base/thread_manager.h>
#endif

template <int dim>
DataOut_DoFData<dim>::DataEntry::DataEntry (const Vector<double> *data,
					    const vector<string> &names) :
		data(data),
		names(names)
{};

template <int dim>
DataOut_DoFData<dim>::DataOut_DoFData () :
		dofs(0)
{};

template <int dim>
DataOut_DoFData<dim>::~DataOut_DoFData ()
{
  clear ();
};

template <int dim>
void DataOut_DoFData<dim>::attach_dof_handler (const DoFHandler<dim> &d)
{
  Assert (dof_data.size() == 0, ExcOldDataStillPresent());
  Assert (cell_data.size() == 0, ExcOldDataStillPresent());
  
  if (dofs != 0)
    dofs->unsubscribe ();
  
  dofs = &d;
  if (dofs != 0)
    dofs->subscribe ();
};

template <int dim>
void DataOut_DoFData<dim>::add_data_vector (const Vector<double> &vec,
					    const vector<string> &names)
{
  Assert (dofs != 0, ExcNoDoFHandlerSelected ());
				   // either cell data and one name,
				   // or dof data and n_components names
  Assert (((vec.size() == dofs->get_tria().n_active_cells()) &&
	   (names.size() == 1))
	  ||
	  ((vec.size() == dofs->n_dofs()) &&
	   (names.size() == dofs->get_fe().n_components())),
	  ExcInvalidNumberOfNames (names.size(), dofs->get_fe().n_components()));
  for (unsigned int i=0; i<names.size(); ++i)
    Assert (names[i].find_first_not_of("abcdefghijklmnopqrstuvwxyz"
				       "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
				       "0123456789_<>()") == string::npos,
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



template <int dim>
void DataOut_DoFData<dim>::add_data_vector (const Vector<double> &vec,
					    const string         &name)
{
  add_data_vector (vec, vector<string>(1,name));
};




template <int dim>
void DataOut_DoFData<dim>::clear ()
{
  dof_data.erase (dof_data.begin(), dof_data.end());
  cell_data.erase (cell_data.begin(), cell_data.end());

  if (dofs != 0)
    {
      dofs->unsubscribe ();
      dofs = 0;
    };

				   // delete patches
  vector<DataOutBase::Patch<dim> > dummy;
  patches.swap (dummy);
};



template <int dim>
vector<string> DataOut_DoFData<dim>::get_dataset_names () const 
{
  vector<string> names;
				   // collect the names of dof
				   // and cell data
  for (vector<DataEntry>::const_iterator d=dof_data.begin(); d!=dof_data.end(); ++d)
    for (unsigned int i=0; i<d->names.size(); ++i)
      names.push_back (d->names[i]);
  for (vector<DataEntry>::const_iterator d=cell_data.begin(); d!=cell_data.end(); ++d)
    {
      Assert (d->names.size() == 1, ExcInternalError());
      names.push_back (d->names[0]);
    };

  return names;
};



template <int dim>
const vector<typename DataOutBase::Patch<dim> > &
DataOut_DoFData<dim>::get_patches () const
{
  return patches;
};


template <int dim>
void * DataOut<dim>::build_some_patches (Data data)
{
  QTrapez<1>     q_trapez;
  QIterated<dim> patch_points (q_trapez, data.n_subdivisions);
  
  FEValues<dim> fe_patch_values(dofs->get_fe(),
				patch_points,
				update_default);

  const unsigned int n_q_points = patch_points.n_quadrature_points;
  
  unsigned int cell_number = 0;
  vector<DataOutBase::Patch<dim> >::iterator patch = patches.begin();
  DoFHandler<dim>::cell_iterator cell=first_cell();

				   // get first cell in this thread
  for (unsigned int i=0;i<data.this_thread;++i,++patch,++cell_number,
		 cell=next_cell(cell));
  
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
      for (unsigned int i=0;((i<data.n_threads)&&(cell != dofs->end()));
	   ++i,++patch,++cell_number,cell=next_cell(cell));
      
    };
  
  return 0;
}




template <int dim>
void DataOut<dim>::build_patches (const unsigned int n_subdivisions,
				  const unsigned int n_threads_) 
{
  unsigned int n_threads = n_threads_;

  Assert (dofs != 0, ExcNoDoFHandlerSelected());

				   // if not in multithread mode
				   // set nr of threads to one
#ifndef DEAL_II_USE_MT
  n_threads = 1;
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
      vector<DataOutBase::Patch<dim> > dummy;
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

  vector<Data> thread_data(n_threads);

				   // init data for the threads
  for (unsigned int i=0;i<n_threads;++i)
    {
      thread_data[i].n_threads      = n_threads;
      thread_data[i].this_thread    = i;
      thread_data[i].n_components   = n_components;
      thread_data[i].n_datasets     = n_datasets;
      thread_data[i].n_subdivisions = n_subdivisions;
      thread_data[i].patch_values.resize(n_q_points);
      thread_data[i].patch_values_system.resize(n_q_points);
      
      for (unsigned int k=0;k<n_components;++k)
	thread_data[i].patch_values_system[k].reinit(n_components);
    }

				   // create the patches with
				   // default values
  DataOutBase::Patch<dim>  default_patch;
  default_patch.n_subdivisions = n_subdivisions;
  default_patch.data.reinit (n_datasets, n_q_points);
  patches.insert (patches.end(), n_patches, default_patch);

#ifdef DEAL_II_USE_MT

  ThreadManager thread_manager;
  
  typedef ThreadManager::Mem_Fun_Data1
    <DataOut<dim>, Data > MemFunData;
  
				   // One struct of this type for every thread
  vector<MemFunData> mem_fun_data (n_threads,
				   MemFunData (this,
					       thread_data[0],
					       &DataOut<dim>::build_some_patches));
  
				   // get start cells for each thread
  for (unsigned int l=0;l<n_threads;++l)
    {
      mem_fun_data[l].arg=thread_data[l];
      thread_manager.spawn(&mem_fun_data[l],THR_SCOPE_SYSTEM | THR_DETACHED);
    };
				   // wait for all threads to return
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
template class DataOut_DoFData<deal_II_dimension>;
template class DataOut<deal_II_dimension>;




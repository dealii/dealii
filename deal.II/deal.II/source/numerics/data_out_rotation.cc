//----------------------------  data_out_rotation.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  data_out_rotation.cc  ---------------------------


#include <base/quadrature_lib.h>
#include <lac/vector.h>
#include <numerics/data_out_rotation.h>
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
#include <cmath>

// if necessary try to work around a bug in the IBM xlC compiler
#ifdef XLC_WORK_AROUND_STD_BUG
using namespace std;
#endif




template <int dim>
void DataOutRotation<dim>::build_some_patches (Data data)
{
  QTrapez<1>     q_trapez;
  QIterated<dim> patch_points (q_trapez, data.n_subdivisions);
  
  FEValues<dim> fe_patch_values(dofs->get_fe(),
				patch_points,
				update_values);

  const unsigned int n_patches_per_circle = data.n_patches_per_circle;

				   // another abbreviation denoting
				   // the number of q_points in each
				   // direction
  const unsigned int n_points = data.n_subdivisions+1;
  
				   // set up an array that holds the
				   // directions in the plane of
				   // rotation in which we will put
				   // points in the whole domain (not
				   // the rotationally reduced one in
				   // which the computation took
				   // place. for simplicity add the
				   // initial direction at the end
				   // again
  const double pi = 3.14159265358979323846;
  std::vector<Point<dim+1> > angle_directions (n_patches_per_circle+1);
  for (unsigned int i=0; i<=n_patches_per_circle; ++i)
    {
      angle_directions[i][0] = cos(2*pi*i/n_patches_per_circle);
      angle_directions[i][1] = sin(2*pi*i/n_patches_per_circle);
    };
  
  
  unsigned int cell_number = 0;
  typename std::vector<DataOutBase::Patch<dim+1> >::iterator patch = patches.begin();
  typename DoFHandler<dim>::cell_iterator cell=first_cell();

				   // get first cell in this thread
  for (unsigned int i=0; (i<data.this_thread)&&(cell != dofs->end()); ++i)
    {
      advance (patch, n_patches_per_circle);
      ++cell_number;
      cell=next_cell(cell);
    }

  				   // now loop over all cells and
				   // actually create the patches
  for (; cell != dofs->end(); )
    {
      for (unsigned int angle=0; angle<n_patches_per_circle; ++angle, ++patch)
	{
	  Assert (patch != patches.end(), ExcInternalError());
	  

					   // first compute the
					   // vertices of the
					   // patch. note that they
					   // will have to be computed
					   // from the vertices of the
					   // cell, which has one
					   // dimension less, however.
	  switch (dim)
	    {
	      case 1:
	      {
		const double r1 = cell->vertex(0)(0),
			     r2 = cell->vertex(1)(0);
		Assert (r1 >= 0, ExcRadialVariableHasNegativeValues(r1));	
		Assert (r2 >= 0, ExcRadialVariableHasNegativeValues(r2));
		
		patch->vertices[0] = r1*angle_directions[angle];
		patch->vertices[1] = r2*angle_directions[angle];
		patch->vertices[2] = r2*angle_directions[angle+1];
		patch->vertices[3] = r1*angle_directions[angle+1];

		break;
	      };
	       
	      case 2:
	      {
		for (unsigned int vertex=0;
		     vertex<GeometryInfo<dim>::vertices_per_cell;
		     ++vertex)
		  {
		    const Point<dim> v = cell->vertex(vertex);
		    
						     // make sure that the
						     // radial variable does
						     // attain negative
						     // values
		    Assert (v(0) >= 0, ExcRadialVariableHasNegativeValues(v(0)));
		    
						     // now set the vertices
						     // of the patch
		    patch->vertices[vertex] = v(0) * angle_directions[angle];
		    patch->vertices[vertex][2] = v(1);
		    
		    patch->vertices[vertex+GeometryInfo<dim>::vertices_per_cell]
		      = v(0) * angle_directions[angle+1];
		    patch->vertices[vertex+GeometryInfo<dim>::vertices_per_cell][2]
		      = v(1);
		  };
		
		break;
	      };

	      default:
		    Assert (false, ExcNotImplemented());
	    };
	  
	       
					   // then fill in data
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
		      switch (dim)
			{
			  case 1:
				for (unsigned int x=0; x<n_points; ++x)
				  for (unsigned int y=0; y<n_points; ++y)
				    patch->data(dataset,
						x*n_points + y)
					= data.patch_values[x];
				break;
				
			  case 2:
				for (unsigned int x=0; x<n_points; ++x)
				  for (unsigned int y=0; y<n_points; ++y)
				    for (unsigned int z=0; z<n_points; ++z)
				      patch->data(dataset,
						  x*n_points*n_points +
						  y*n_points +
						  z)
					= data.patch_values[x*n_points+z];
				break;
				
			  default:
				Assert (false, ExcNotImplemented());
			};
		    }
		  else
						     // system of components
		    {
		      fe_patch_values.get_function_values (*dof_data[dataset].data,
							   data.patch_values_system);
		      for (unsigned int component=0; component<data.n_components;
			   ++component)
			{
			  switch (dim)
			    {
			      case 1:
				    for (unsigned int x=0; x<n_points; ++x)
				      for (unsigned int y=0; y<n_points; ++y)
					patch->data(dataset*data.n_components+component,
						    x*n_points + y)
					    = data.patch_values_system[x](component);
				    break;

			      case 2:
				    for (unsigned int x=0; x<n_points; ++x)
				      for (unsigned int y=0; y<n_points; ++y)
					for (unsigned int z=0; z<n_points; ++z)
					  patch->data(dataset*data.n_components+component,
						      x*n_points*n_points +
						      y*n_points +
						      z)
					    = data.patch_values_system[x*n_points+z](component);
				    break;

			      default:
				    Assert (false, ExcNotImplemented());
			    };
			};
		    };
		};
		  
					       // then do the cell data
	      for (unsigned int dataset=0; dataset<cell_data.size(); ++dataset)
		{
		  const double value = (*cell_data[dataset].data)(cell_number);
		  switch (dim)
		    {
		      case 1:
			    for (unsigned int x=0; x<n_points; ++x)
			      for (unsigned int y=0; y<n_points; ++y)
				patch->data(dataset+dof_data.size()*data.n_components,
					    x*n_points +
					    y)
				  = value;
			    break;
			    
		      case 2:
			    for (unsigned int x=0; x<n_points; ++x)
			      for (unsigned int y=0; y<n_points; ++y)
				for (unsigned int z=0; z<n_points; ++z)
				  patch->data(dataset+dof_data.size()*data.n_components,
					      x*n_points*n_points +
					      y*n_points +
					      z)
				    = value;
			    break;

		      default:
			    Assert (false, ExcNotImplemented());
		    };
		};
	    };
	};
      
				       // next cell (patch) in this
				       // thread. note that we have
				       // already advanced the patches
				       // for the present cell,
				       // i.e. we only have to skip
				       // the cells belonging to other
				       // threads, not the ones
				       // belonging to this thread.
      const int skip_threads = static_cast<signed int>(data.n_threads)-1;
      for (int i=0; (i<skip_threads) && (cell != dofs->end()); ++i)
	advance (patch, n_patches_per_circle);

				       // however, cell and cell
				       // number have not yet been
				       // increased
      for (unsigned int i=0; (i<data.n_threads) && (cell != dofs->end()); ++i)
	{
	  ++cell_number;
	  cell=next_cell(cell);
	};
    };
};



#if deal_II_dimension == 3

template <>
void DataOutRotation<3>::build_some_patches (Data)
{
				   // would this function make any
				   // sense after all? who would want
				   // to output/compute in four space
				   // dimensions?
  Assert (false, ExcNotImplemented());
};

#endif



template <int dim>
void DataOutRotation<dim>::build_patches (const unsigned int n_patches_per_circle,
					  const unsigned int n_subdivisions,
					  const unsigned int n_threads_) 
{
  Assert (n_subdivisions >= 1,
	  ExcInvalidNumberOfSubdivisions(n_subdivisions));

  typedef DataOut_DoFData<dim,dim+1> BaseClass;
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
      std::vector<DataOutBase::Patch<dim+1> > dummy;
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
				   // then also take into account that
				   // we want more than one patch to
				   // come out of every cell, as they
				   // are repeated around the axis of
				   // rotation
  n_patches *= n_patches_per_circle;

  std::vector<Data> thread_data(n_threads);

				   // init data for the threads
  for (unsigned int i=0;i<n_threads;++i)
    {
      thread_data[i].n_threads            = n_threads;
      thread_data[i].this_thread          = i;
      thread_data[i].n_components         = n_components;
      thread_data[i].n_datasets           = n_datasets;
      thread_data[i].n_patches_per_circle = n_patches_per_circle;
      thread_data[i].n_subdivisions       = n_subdivisions;
      thread_data[i].patch_values.resize (n_q_points);
      thread_data[i].patch_values_system.resize (n_q_points);
      
      for (unsigned int k=0; k<n_q_points; ++k)
	thread_data[i].patch_values_system[k].reinit(n_components);
    }

				   // create the patches with default
				   // values. note that the evaluation
				   // points on the cell have to be
				   // repeated in angular direction
  DataOutBase::Patch<dim+1>  default_patch;
  default_patch.n_subdivisions = n_subdivisions;
  default_patch.data.reinit (n_datasets,
			     n_q_points * (n_subdivisions+1));
  patches.insert (patches.end(), n_patches, default_patch);

#ifdef DEAL_II_USE_MT

  Threads::ThreadManager thread_manager;  
  for (unsigned int l=0;l<n_threads;++l)
    Threads::spawn (thread_manager,
		    Threads::encapsulate (&DataOutRotation<dim>::build_some_patches)
		    .collect_args (this, thread_data[l]));
  thread_manager.wait();
  
				   // just one thread
#else
  build_some_patches(thread_data[0]);
#endif
};


template <int dim>
typename DoFHandler<dim>::cell_iterator
DataOutRotation<dim>::first_cell () 
{
  return dofs->begin_active ();
};


template <int dim>
typename DoFHandler<dim>::cell_iterator
DataOutRotation<dim>::next_cell (const typename DoFHandler<dim>::cell_iterator &cell) 
{
				   // convert the iterator to an
				   // active_iterator and advance
				   // this to the next active cell
  typename DoFHandler<dim>::active_cell_iterator active_cell = cell;
  ++active_cell;
  return active_cell;
};


// explicit instantiations
template class DataOutRotation<deal_II_dimension>;



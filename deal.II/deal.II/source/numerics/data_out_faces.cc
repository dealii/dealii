//----------------------------  data_out_faces.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  data_out_faces.cc  ---------------------------


#include <base/quadrature_lib.h>
#include <base/thread_management.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <numerics/data_out_faces.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <fe/mapping_q1.h>



template <int dim>
void DataOutFaces<dim>::build_some_patches (Data data)
{
  QTrapez<1>        q_trapez;
  QIterated<dim-1>  patch_points (q_trapez, data.n_subdivisions);
  
				   // since most output formats can't
				   // handle cells that are not
				   // transformed using a Q1 mapping,
				   // we don't support anything else
				   // as well
  static const MappingQ1<dim> mapping;
  FEFaceValues<dim> fe_patch_values (mapping, this->dofs->get_fe(),
				     patch_points, update_values);

  const unsigned int n_q_points = patch_points.n_quadrature_points;
  
  unsigned int face_number = 0;
  typename std::vector< ::DataOutBase::Patch<dim-1,dim> >::iterator patch = this->patches.begin();
  FaceDescriptor face=first_face();

				   // get first face in this thread
  for (unsigned int i=0; (i<data.this_thread)&&(face.first != this->dofs->end()); ++i)
    {
      ++patch;
      ++face_number;
      face=next_face(face);
    }

  				   // now loop over all cells and
				   // actually create the patches
  for (; face.first != this->dofs->end();)
    {
      Assert (patch != this->patches.end(), ExcInternalError());
      
      for (unsigned int vertex=0; vertex<GeometryInfo<dim-1>::vertices_per_cell; ++vertex)
	patch->vertices[vertex] = face.first->face(face.second)->vertex(vertex);
      
      if (data.n_datasets > 0)
	{
	  fe_patch_values.reinit (face.first, face.second);
	  
					   // first fill dof_data
	  for (unsigned int dataset=0; dataset<this->dof_data.size(); ++dataset)
	    {
	      if (data.n_components == 1)
		{
		  if (this->dof_data[dataset].has_block)
		    fe_patch_values.get_function_values (*this->dof_data[dataset].block_data,
							 data.patch_values);
		  else
		    fe_patch_values.get_function_values (*this->dof_data[dataset].single_data,
							 data.patch_values);

		  for (unsigned int q=0; q<n_q_points; ++q)
		    patch->data(dataset,q) = data.patch_values[q];
		}
	      else
						 // system of components
		{
		  if (this->dof_data[dataset].has_block)
		    fe_patch_values.get_function_values (*this->dof_data[dataset].block_data,
							 data.patch_values_system);
		  else
		    fe_patch_values.get_function_values (*this->dof_data[dataset].single_data,
							 data.patch_values_system);

		  for (unsigned int component=0; component<data.n_components;
		       ++component)
		    for (unsigned int q=0; q<n_q_points; ++q)
		      patch->data(dataset*data.n_components+component,q) =
			data.patch_values_system[q](component);
		};
	    };

					   // then do the cell data
	  for (unsigned int dataset=0; dataset<this->cell_data.size(); ++dataset)
	    {
	      if (this->cell_data[dataset].has_block)
		{
		  const double value = (*this->cell_data[dataset].block_data)(face_number);
		  for (unsigned int q=0; q<n_q_points; ++q)
		    patch->data(dataset+this->dof_data.size()*data.n_components,q) =
		      value;
		} else {
		  const double value = (*this->cell_data[dataset].single_data)(face_number);
		  for (unsigned int q=0; q<n_q_points; ++q)
		    patch->data(dataset+this->dof_data.size()*data.n_components,q) =
		      value;
		} 
	    };
	};
      				       // next cell (patch) in this thread
      for (unsigned int i=0;
	   (i<data.n_threads)&&(face.first != this->dofs->end()); ++i)
	{
	  ++patch;
	  ++face_number;
	  face=next_face(face);
	}
    };
}




template <int dim>
void DataOutFaces<dim>::build_patches (const unsigned int n_subdivisions,
				       const unsigned int n_threads_) 
{
  Assert (n_subdivisions >= 1,
	  ExcInvalidNumberOfSubdivisions(n_subdivisions));

  typedef DataOut_DoFData<dim,dim+1> BaseClass;
  Assert (this->dofs != 0, typename BaseClass::ExcNoDoFHandlerSelected());

  const unsigned int n_threads = (DEAL_II_USE_MT ? n_threads_ : 1);

				   // before we start the loop:
				   // create a quadrature rule that
				   // actually has the points on this
				   // patch
  QTrapez<1>       q_trapez;
  QIterated<dim-1> patch_points (q_trapez, n_subdivisions);

  const unsigned int n_q_points     = patch_points.n_quadrature_points;
  const unsigned int n_components   = this->dofs->get_fe().n_components();
  const unsigned int n_datasets     = this->dof_data.size() * n_components +
				      this->cell_data.size();
  
				   // clear the patches array
  if (true)
    {
      std::vector< ::DataOutBase::Patch<dim-1,dim> > dummy;
      this->patches.swap (dummy);
    };
  
				   // first count the cells we want to
				   // create patches of and make sure
				   // there is enough memory for that
  unsigned int n_patches = 0;
  for (FaceDescriptor face=first_face();
       face.first != this->dofs->end();
       face = next_face(face))
    ++n_patches;

  std::vector<Data> thread_data(n_threads);

				   // init data for the threads
  for (unsigned int i=0;i<n_threads;++i)
    {
      thread_data[i].n_threads            = n_threads;
      thread_data[i].this_thread          = i;
      thread_data[i].n_components         = n_components;
      thread_data[i].n_datasets           = n_datasets;
      thread_data[i].n_subdivisions       = n_subdivisions;
      thread_data[i].patch_values.resize (n_q_points);
      thread_data[i].patch_values_system.resize (n_q_points);
      
      for (unsigned int k=0; k<n_q_points; ++k)
	thread_data[i].patch_values_system[k].reinit(n_components);
    }

				   // create the patches with default
				   // values. note that the evaluation
				   // points on the face have to be
				   // repeated in angular direction
  ::DataOutBase::Patch<dim-1,dim>  default_patch;
  default_patch.n_subdivisions = n_subdivisions;
  default_patch.data.reinit (n_datasets, n_q_points);
  this->patches.insert (this->patches.end(), n_patches, default_patch);

  if (DEAL_II_USE_MT)
    {
      Threads::ThreadManager thread_manager;  
      for (unsigned int l=0;l<n_threads;++l)
        Threads::spawn (thread_manager,
                        Threads::encapsulate (&DataOutFaces<dim>::build_some_patches)
                        .collect_args (this, thread_data[l]));
      thread_manager.wait();
    }
  else
				   // just one thread
    build_some_patches(thread_data[0]);
}



template <int dim>
typename DataOutFaces<dim>::FaceDescriptor
DataOutFaces<dim>::first_face () 
{
				   // simply find first active cell
				   // with a face on the boundary
  typename DoFHandler<dim>::active_cell_iterator cell = this->dofs->begin_active();
  for (; cell != this->dofs->end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if (cell->face(f)->at_boundary())
	return FaceDescriptor(cell, f);

				   // ups, triangulation has no
				   // boundary? impossible!
  Assert (false, ExcInternalError());
  
  return FaceDescriptor();
}



template <int dim>
typename DataOutFaces<dim>::FaceDescriptor
DataOutFaces<dim>::next_face (const FaceDescriptor &old_face)
{
  FaceDescriptor face = old_face;
  
				   // first check whether the present
				   // cell has more faces on the
				   // boundary
  for (unsigned int f=face.second+1; f<GeometryInfo<dim>::faces_per_cell; ++f)
    if (face.first->face(f)->at_boundary())
				       // yup, that is so, so return it
      {
	face.second = f;
	return face;
      };

				   // otherwise find the next active
				   // cell that has a face on the
				   // boundary
  
				   // convert the iterator to an
				   // active_iterator and advance
				   // this to the next active cell
  typename DoFHandler<dim>::active_cell_iterator active_cell = face.first;

				   // increase face pointer by one
  ++active_cell;

				   // while there are active cells
  while (active_cell != this->dofs->end())
    {
				       // check all the faces of this
				       // active cell
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
	if (active_cell->face(f)->at_boundary())
	  {
	    face.first  = active_cell;
	    face.second = f;
	    return face;
	  };
				       // the present cell had no
				       // faces on the boundary, so
				       // check next cell
      ++active_cell;
    };

				   // we fell off the edge, so return
				   // with invalid pointer
  face.first  = this->dofs->end();
  face.second = 0;
  return face;
}


// explicit instantiations
// don't instantiate anything for the 1d and 2d cases
#if deal_II_dimension >=2
template class DataOutFaces<deal_II_dimension>;
#endif

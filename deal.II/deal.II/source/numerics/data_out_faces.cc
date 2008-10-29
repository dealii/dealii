//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


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
#include <hp/fe_values.h>
#include <fe/mapping_q1.h>

DEAL_II_NAMESPACE_OPEN



template <int dim, class DH>
void DataOutFaces<dim,DH>::build_some_patches (internal::DataOut::ParallelData<DH::dimension, DH::dimension> &data)
{
				   // Check consistency of redundant
				   // template parameter
  Assert (dim==DH::dimension, ExcDimensionMismatch(dim, DH::dimension));
  
  QTrapez<1>        q_trapez;
  QIterated<DH::dimension-1>  patch_points (q_trapez, data.n_subdivisions);
  
//TODO[?]: This is strange -- Data has a member 'mapping' that should
//be used here, but it isn't. Rather, up until version 1.94, we were
//actually initializing a local mapping object and used that... While
//we use the mapping to transform the vertex coordinates, we do not
//use the mapping to transform the shape functions (necessary for
//example for Raviart-Thomas elements). This could lead to trouble
//when someone tries to use MappingEulerian with such elements

				   // create collection objects from
				   // single quadratures,
				   // and finite elements. if we have
				   // an hp DoFHandler,
				   // dof_handler.get_fe() returns a
				   // collection of which we do a
				   // shallow copy instead
  const hp::QCollection<DH::dimension-1>     q_collection (patch_points);
  const hp::FECollection<DH::dimension>      fe_collection(this->dofs->get_fe());
  
  UpdateFlags update_flags=update_values;
  for (unsigned int i=0; i<this->dof_data.size(); ++i)
    if (this->dof_data[i]->postprocessor)
      update_flags |= this->dof_data[i]->postprocessor->get_needed_update_flags();
  
  hp::FEFaceValues<DH::dimension> x_fe_patch_values (fe_collection, q_collection,
						     update_flags);

  const unsigned int n_q_points = patch_points.n_quadrature_points;
  
  typename std::vector< dealii::DataOutBase::Patch<DH::dimension-1,DH::dimension> >::iterator patch = this->patches.begin();
  FaceDescriptor face=first_face();

				   // get first face in this thread
  for (unsigned int i=0; (i<data.this_thread)&&(face.first != this->dofs->end()); ++i)
    {
      ++patch;
      face=next_face(face);
    }

  				   // now loop over all cells and
				   // actually create the patches
  for (; face.first != this->dofs->end();)
    {
      Assert (patch != this->patches.end(), ExcInternalError());
      
      for (unsigned int vertex=0; vertex<GeometryInfo<DH::dimension-1>::vertices_per_cell; ++vertex)
	patch->vertices[vertex] = face.first->face(face.second)->vertex(vertex);
      
      if (data.n_datasets > 0)
	{
	  x_fe_patch_values.reinit (face.first, face.second);
          const FEFaceValues<DH::dimension> &fe_patch_values
            = x_fe_patch_values.get_present_fe_values ();
	  
					   // counter for data records
	  unsigned int offset=0;
	  
					   // first fill dof_data
	  for (unsigned int dataset=0; dataset<this->dof_data.size(); ++dataset)
	    {
	      const DataPostprocessor<dim> *postprocessor=this->dof_data[dataset]->postprocessor;
	      if (postprocessor != 0)
		{
						   // we have to postprocess the
						   // data, so determine, which
						   // fields have to be updated
		  const UpdateFlags update_flags=postprocessor->get_needed_update_flags();
		  
						   // get normals, if
						   // needed. this is a
						   // geometrical information
						   // and thus does not depend
						   // on the number of
						   // components of the data
						   // vector
		  if (update_flags & update_normal_vectors)
		    data.patch_normals=fe_patch_values.get_normal_vectors();

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
							  data.patch_normals,
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
							  data.patch_normals,
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

					   // then do the cell data
	  for (unsigned int dataset=0; dataset<this->cell_data.size(); ++dataset)
	    {
					       // we need to get at
					       // the number of the
					       // cell to which this
					       // face belongs in
					       // order to access the
					       // cell data. this is
					       // not readily
					       // available, so choose
					       // the following rather
					       // inefficient way:
	      Assert (face.first->active(), ExcCellNotActiveForCellData());
	      const unsigned int cell_number
		= std::distance (this->dofs->begin_active(),
				 typename DH::active_cell_iterator(face.first));

              const double value
                = this->cell_data[dataset]->get_cell_data_value (cell_number);
              for (unsigned int q=0; q<n_q_points; ++q)
                patch->data(dataset+offset,q) =
                  value;
	    }
	}
      				       // next cell (patch) in this thread
      for (unsigned int i=0;
	   (i<data.n_threads)&&(face.first != this->dofs->end()); ++i)
	{
	  ++patch;
	  face=next_face(face);
	}
    }
}




template <int dim, class DH>
void DataOutFaces<dim,DH>::build_patches (const unsigned int nnnn_subdivisions,
					  const unsigned int n_threads_) 
{
  unsigned int n_subdivisions = (nnnn_subdivisions != 0)
				? nnnn_subdivisions
				: this->default_subdivisions;
  
  Assert (n_subdivisions >= 1,
	  ExcInvalidNumberOfSubdivisions(n_subdivisions));

  typedef DataOut_DoFData<DH,DH::dimension+1> BaseClass;
  Assert (this->dofs != 0, typename BaseClass::ExcNoDoFHandlerSelected());

  Assert (!DEAL_II_USE_MT || (n_threads_ >= 1),
	  ExcMessage ("Must run with at least one thread!"));
  const unsigned int n_threads = (DEAL_II_USE_MT ? n_threads_ : 1);

				   // before we start the loop:
				   // create a quadrature rule that
				   // actually has the points on this
				   // patch
  QTrapez<1>       q_trapez;
  QIterated<DH::dimension-1> patch_points (q_trapez, n_subdivisions);

  const unsigned int n_q_points     = patch_points.n_quadrature_points;
  const unsigned int n_components   = this->dofs->get_fe().n_components();
  unsigned int n_datasets     = this->cell_data.size();
  for (unsigned int i=0; i<this->dof_data.size(); ++i)
    n_datasets+= this->dof_data[i]->n_output_variables;
  
				   // clear the patches array
  if (true)
    {
      std::vector< dealii::DataOutBase::Patch<DH::dimension-1,DH::dimension> > dummy;
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

  std::vector<internal::DataOut::ParallelData<DH::dimension, DH::dimension> > thread_data(n_threads);

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
      thread_data[i].patch_gradients.resize (n_q_points);
      thread_data[i].patch_gradients_system.resize (n_q_points);
      thread_data[i].patch_hessians.resize (n_q_points);
      thread_data[i].patch_hessians_system.resize (n_q_points);
      thread_data[i].patch_normals.resize (n_q_points);
      thread_data[i].postprocessed_values.resize (this->dof_data.size());
      
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

				   // create the patches with default
				   // values. note that the evaluation
				   // points on the face have to be
				   // repeated in angular direction
  dealii::DataOutBase::Patch<DH::dimension-1,DH::dimension>  default_patch;
  default_patch.n_subdivisions = n_subdivisions;
  default_patch.data.reinit (n_datasets, n_q_points);
  this->patches.insert (this->patches.end(), n_patches, default_patch);

  if (DEAL_II_USE_MT)
    {
      Threads::ThreadGroup<> threads;  
      for (unsigned int l=0;l<n_threads;++l)
        threads += Threads::spawn (*this, &DataOutFaces<dim,DH>::build_some_patches)(thread_data[l]);
      threads.join_all();
    }
  else
				     // just one thread
    build_some_patches(thread_data[0]);
}



template <int dim, class DH>
typename DataOutFaces<dim,DH>::FaceDescriptor
DataOutFaces<dim,DH>::first_face () 
{
				   // simply find first active cell
				   // with a face on the boundary
  typename DH::active_cell_iterator cell = this->dofs->begin_active();
  for (; cell != this->dofs->end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<DH::dimension>::faces_per_cell; ++f)
      if (cell->face(f)->at_boundary())
	return FaceDescriptor(cell, f);

				   // ups, triangulation has no
				   // boundary? impossible!
  Assert (false, ExcInternalError());
  
  return FaceDescriptor();
}



template <int dim, class DH>
typename DataOutFaces<dim,DH>::FaceDescriptor
DataOutFaces<dim,DH>::next_face (const FaceDescriptor &old_face)
{
  FaceDescriptor face = old_face;
  
				   // first check whether the present
				   // cell has more faces on the
				   // boundary
  for (unsigned int f=face.second+1; f<GeometryInfo<DH::dimension>::faces_per_cell; ++f)
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
  typename DH::active_cell_iterator active_cell = face.first;

				   // increase face pointer by one
  ++active_cell;

				   // while there are active cells
  while (active_cell != this->dofs->end())
    {
				       // check all the faces of this
				       // active cell
      for (unsigned int f=0; f<GeometryInfo<DH::dimension>::faces_per_cell; ++f)
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
template class DataOutFaces<deal_II_dimension, DoFHandler<deal_II_dimension> >;
template class DataOutFaces<deal_II_dimension, hp::DoFHandler<deal_II_dimension> >;
#endif

DEAL_II_NAMESPACE_CLOSE

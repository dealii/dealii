/* $Id$ */


#include <base/quadrature_lib.h>
#include <lac/vector.h>
#include <basic/data_out.h>
#include <grid/tria.h>
#include <grid/dof.h>
#include <grid/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <fe/fe.h>
#include <fe/fe_values.h>


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
	   (names.size() == dofs->get_fe().n_components)),
	  ExcInvalidNumberOfNames (names.size(), dofs->get_fe().n_components));
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
void DataOut<dim>::build_patches (const unsigned int n_subdivisions) 
{
  Assert (dofs != 0, ExcNoDoFHandlerSelected());
  
  const unsigned int n_components   = dofs->get_fe().n_components;
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
  for (DoFHandler<dim>::active_cell_iterator cell=first_cell();
       cell != dofs->end();
       cell = next_cell(cell))
    ++n_patches;


				   // before we start the loop:
				   // create a quadrature rule that
				   // actually has the points on this
				   // patch, and an object that
				   // extracts the data on each
				   // cell to these points
  QTrapez<1>     q_trapez;
  QIterated<dim> patch_points (q_trapez, n_subdivisions);
  FEValues<dim>  fe_patch_values (dofs->get_fe(),
				  patch_points,
				  update_default);
  const unsigned int n_q_points = patch_points.n_quadrature_points;
  vector<double>          patch_values (n_q_points);
  vector<Vector<double> > patch_values_system (n_q_points,
					       Vector<double>(n_components));
  
  DataOutBase::Patch<dim>  default_patch;
  default_patch.n_subdivisions = n_subdivisions;
  default_patch.data.reinit (n_datasets, n_q_points);
  patches.insert (patches.end(), n_patches, default_patch);

				   // now loop over all cells and
				   // actually create the patches
  vector<DataOutBase::Patch<dim> >::iterator patch = patches.begin();
  unsigned int                               cell_number = 0;
  for (DoFHandler<dim>::active_cell_iterator cell=dofs->begin_active();
       cell != dofs->end(); ++cell, ++patch, ++cell_number)
    {
      Assert (patch != patches.end(), ExcInternalError());
      
      for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
	patch->vertices[vertex] = cell->vertex(vertex);

      if (n_datasets > 0)
	{
	  fe_patch_values.reinit (cell);
	  
					   // first fill dof_data
	  for (unsigned int dataset=0; dataset<dof_data.size(); ++dataset)
	    {
	      if (n_components == 1)
		{
		  fe_patch_values.get_function_values (*dof_data[dataset].data,
						       patch_values);
		  for (unsigned int q=0; q<n_q_points; ++q)
		    patch->data(dataset,q) = patch_values[q];
		}
	      else
						 // system of components
		{
		  fe_patch_values.get_function_values (*dof_data[dataset].data,
						       patch_values_system);
		  for (unsigned int component=0; component<n_components; ++component)
		    for (unsigned int q=0; q<n_q_points; ++q)
		      patch->data(dataset*n_components+component,q) = patch_values_system[q](component);
		};
	    };

					   // then do the cell data
	  for (unsigned int dataset=0; dataset<cell_data.size(); ++dataset)
	    {
	      const double value = (*cell_data[dataset].data)(cell_number);
	      for (unsigned int q=0; q<n_q_points; ++q)
		patch->data(dataset+dof_data.size()*n_components,q) = value;
	    };
	};
    };
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




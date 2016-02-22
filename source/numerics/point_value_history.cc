// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/parallel_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <deal.II/numerics/point_value_history.h>

#include <algorithm>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace PointValueHistory
  {
/// Only a constructor needed for this class (a struct really)
    template <int dim>
    PointGeometryData<dim>
    ::PointGeometryData (const Point <dim> &new_requested_location,
                         const std::vector <Point <dim> > &new_locations,
                         const std::vector <types::global_dof_index> &new_sol_indices)
    {
      requested_location = new_requested_location;
      support_point_locations = new_locations;
      solution_indices = new_sol_indices;
    }
  }
}



template <int dim>
PointValueHistory<dim>
::PointValueHistory (const unsigned int n_independent_variables) :
  n_indep (n_independent_variables)
{
  closed = false;
  cleared = false;
  triangulation_changed = false;
  have_dof_handler = false;

  // make a vector for keys
  dataset_key = std::vector <double> (); // initialize the std::vector

  // make a vector of independent values
  independent_values
    = std::vector<std::vector <double> > (n_indep, std::vector <double> (0));
  indep_names = std::vector <std::string> ();
}



template <int dim>
PointValueHistory<dim>::PointValueHistory (const DoFHandler<dim> &dof_handler,
                                           const unsigned int n_independent_variables) :
  dof_handler (&dof_handler),
  n_indep (n_independent_variables)
{
  closed = false;
  cleared = false;
  triangulation_changed = false;
  have_dof_handler = true;

  // make a vector to store keys
  dataset_key = std::vector <double> (); // initialize the std::vector

  // make a vector for the independent values
  independent_values
    = std::vector<std::vector <double> > (n_indep, std::vector <double> (0));
  indep_names = std::vector <std::string> ();

  tria_listener = dof_handler.get_triangulation().signals.any_change.connect (std_cxx11::bind (&PointValueHistory<dim>::tria_change_listener,
                  std_cxx11::ref(*this)));
}



template <int dim>
PointValueHistory<dim>::PointValueHistory (const PointValueHistory &point_value_history)
{
  dataset_key = point_value_history.dataset_key;
  independent_values = point_value_history.independent_values;
  indep_names = point_value_history.indep_names;
  data_store = point_value_history.data_store;
  component_mask = point_value_history.component_mask;
  component_names_map = point_value_history.component_names_map;
  point_geometry_data = point_value_history.point_geometry_data;

  closed = point_value_history.closed;
  cleared = point_value_history.cleared;

  dof_handler = point_value_history.dof_handler;

  triangulation_changed = point_value_history.triangulation_changed;
  have_dof_handler = point_value_history.have_dof_handler;
  n_indep = point_value_history.n_indep;

  // What to do with tria_listener?
  // Presume subscribe new instance?
  if (have_dof_handler)
    {
      tria_listener = dof_handler->get_triangulation().signals.any_change.connect (std_cxx11::bind     (&PointValueHistory<dim>::tria_change_listener,
                      std_cxx11::ref(*this)));
    }
}



template <int dim>
PointValueHistory<dim> &
PointValueHistory<dim>::operator= (const PointValueHistory &point_value_history)
{
  dataset_key = point_value_history.dataset_key;
  independent_values = point_value_history.independent_values;
  indep_names = point_value_history.indep_names;
  data_store = point_value_history.data_store;
  component_mask = point_value_history.component_mask;
  component_names_map = point_value_history.component_names_map;
  point_geometry_data = point_value_history.point_geometry_data;

  closed = point_value_history.closed;
  cleared = point_value_history.cleared;

  dof_handler = point_value_history.dof_handler;

  triangulation_changed = point_value_history.triangulation_changed;
  have_dof_handler = point_value_history.have_dof_handler;
  n_indep = point_value_history.n_indep;

  // What to do with tria_listener?
  // Presume subscribe new instance?
  if (have_dof_handler)
    {
      tria_listener = dof_handler->get_triangulation().signals.any_change.connect (std_cxx11::bind     (&PointValueHistory<dim>::tria_change_listener,
                      std_cxx11::ref(*this)));
    }

  return * this;
}



template <int dim>
PointValueHistory<dim>
::~PointValueHistory ()
{
  if (have_dof_handler)
    {
      tria_listener.disconnect ();
    }
}



template <int dim>
void PointValueHistory<dim>
::add_point (const Point <dim> &location)
{
  // can't be closed to add additional points
  // or vectors
  AssertThrow (!closed, ExcInvalidState ());
  AssertThrow (!cleared, ExcInvalidState ());
  AssertThrow (have_dof_handler, ExcDoFHandlerRequired ());
  AssertThrow (!triangulation_changed, ExcDoFHandlerChanged ());

  // Implementation assumes that support
  // points locations are dofs locations
  AssertThrow (dof_handler->get_fe ().has_support_points (), ExcNotImplemented ());

  // FEValues object to extract quadrature
  // points from
  std::vector <Point <dim> >
  unit_support_points = dof_handler->get_fe ().get_unit_support_points ();

  // While in general quadrature points seems
  // to refer to Gauss quadrature points, in
  // this case the quadrature points are
  // forced to be the support points of the
  // FE.
  Quadrature<dim>
  support_point_quadrature (dof_handler->get_fe ().get_unit_support_points ());
  FEValues<dim> fe_values (dof_handler->get_fe (),
                           support_point_quadrature,
                           update_quadrature_points);
  unsigned int n_support_points
    = dof_handler->get_fe ().get_unit_support_points ().size ();
  unsigned int n_components
    = dof_handler->get_fe ().n_components ();

  // set up a loop over all the cells in the
  // DoFHandler
  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler->begin_active ();
  typename DoFHandler<dim>::active_cell_iterator
  endc = dof_handler->end ();

  // default values to be replaced as closer
  // points are found however they need to be
  // consistent in case they are actually
  // chosen
  typename DoFHandler<dim>::active_cell_iterator current_cell = cell;
  std::vector <unsigned int> current_fe_index (n_components, 0); // need one index per component
  fe_values.reinit (cell);
  std::vector <Point <dim> > current_points (n_components, Point <dim > ());
  for (unsigned int support_point = 0;
       support_point < n_support_points; support_point++)
    {
      // setup valid data in the empty
      // vectors
      unsigned int component
        = dof_handler->get_fe ().system_to_component_index (support_point).first;
      current_points [component] = fe_values.quadrature_point (support_point);
      current_fe_index [component] = support_point;
    }

  // check each cell to find a suitable
  // support points
  // GridTools::find_active_cell_around_point
  // is an alternative. That method is not
  // used here mostly because of the history
  // of the class. The algorithm used in
  // add_points below may be slightly more
  // efficient than find_active_cell_around_point
  // because it operates on a set of points.

  for (; cell != endc; cell++)
    {
      fe_values.reinit (cell);

      for (unsigned int support_point = 0;
           support_point < n_support_points; support_point++)
        {
          unsigned int component
            = dof_handler->get_fe ().system_to_component_index (support_point).first;
          Point<dim> test_point
            = fe_values.quadrature_point (support_point);

          if (location.distance (test_point) <
              location.distance (current_points [component]))
            {
              // save the data
              current_points [component] = test_point;
              current_cell = cell;
              current_fe_index [component] = support_point;
            }
        }
    }


  std::vector<types::global_dof_index>
  local_dof_indices (dof_handler->get_fe ().dofs_per_cell);
  std::vector <types::global_dof_index> new_solution_indices;
  current_cell->get_dof_indices (local_dof_indices);
  // there is an implicit assumption here
  // that all the closest support point to
  // the requested point for all finite
  // element components lie in the same cell.
  // this could possibly be violated if
  // components use different fe orders,
  // requested points are on the edge or
  // vertex of a cell and we are unlucky with
  // floating point rounding. Worst case
  // scenario however is that the point
  // selected isn't the closest possible, it
  // will still lie within one cell distance.
  // calling
  // GridTools::find_active_cell_around_point
  // to obtain a cell to search is an
  // option for these methods, but currently
  // the GridTools function does not cater for
  // a vector of points, and does not seem to
  // be intrinsicly faster than this method.
  for (unsigned int component = 0;
       component < dof_handler->get_fe ().n_components (); component++)
    {
      new_solution_indices
      .push_back (local_dof_indices[current_fe_index [component]]);
    }

  internal::PointValueHistory::PointGeometryData<dim>
  new_point_geometry_data (location, current_points, new_solution_indices);
  point_geometry_data.push_back (new_point_geometry_data);

  std::map <std::string, std::vector <std::vector <double> > >::iterator
  data_store_begin = data_store.begin ();
  for (; data_store_begin != data_store.end (); data_store_begin++)
    {
      // add an extra row to each vector
      // entry
      const ComponentMask &current_mask = (component_mask.find (data_store_begin->first))->second;
      unsigned int n_stored = current_mask.n_selected_components();
      for (unsigned int component = 0; component < n_stored; component++)
        {
          data_store_begin->second.push_back (std::vector<double> (0));
        }
    }
}



template <int dim>
void PointValueHistory<dim>
::add_points (const std::vector <Point <dim> > &locations)
{
  // This algorithm adds points in the same
  // order as they appear in the vector
  // locations and users may depend on this
  // so do not change order added!

  // can't be closed to add additional points or vectors
  AssertThrow (!closed, ExcInvalidState ());
  AssertThrow (!cleared, ExcInvalidState ());
  AssertThrow (have_dof_handler, ExcDoFHandlerRequired ());
  AssertThrow (!triangulation_changed, ExcDoFHandlerChanged ());


  // Implementation assumes that support
  // points locations are dofs locations
  AssertThrow (dof_handler->get_fe ().has_support_points (), ExcNotImplemented ());

  // FEValues object to extract quadrature
  // points from
  std::vector <Point <dim> > unit_support_points = dof_handler->get_fe ().get_unit_support_points ();

  // While in general quadrature points seems
  // to refer to Gauss quadrature points, in
  // this case the quadrature points are
  // forced to be the support points of the
  // FE.
  Quadrature<dim> support_point_quadrature (dof_handler->get_fe ().get_unit_support_points ());
  FEValues<dim> fe_values (dof_handler->get_fe (), support_point_quadrature, update_quadrature_points);
  unsigned int n_support_points = dof_handler->get_fe ().get_unit_support_points ().size ();
  unsigned int n_components = dof_handler->get_fe ().n_components ();

  // set up a loop over all the cells in the
  // DoFHandler
  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler->begin_active ();
  typename DoFHandler<dim>::active_cell_iterator endc = dof_handler->end ();

  // default values to be replaced as closer
  // points are found however they need to be
  // consistent in case they are actually
  // chosen vector <vector>s defined where
  // previously single vectors were used

  // need to store one value per point per component
  std::vector <typename DoFHandler<dim>::active_cell_iterator > current_cell (locations.size (), cell);

  fe_values.reinit (cell);
  std::vector <Point <dim> > temp_points (n_components, Point <dim > ());
  std::vector <unsigned int> temp_fe_index (n_components, 0);
  for (unsigned int support_point = 0; support_point < n_support_points; support_point++)
    {
      // setup valid data in the empty
      // vectors
      unsigned int component = dof_handler->get_fe ().system_to_component_index (support_point).first;
      temp_points [component] = fe_values.quadrature_point (support_point);
      temp_fe_index [component] = support_point;
    }
  std::vector <std::vector <Point <dim> > > current_points (locations.size (), temp_points); // give a valid start point
  std::vector <std::vector <unsigned int> > current_fe_index (locations.size (), temp_fe_index);

  // check each cell to find suitable support
  // points
  // GridTools::find_active_cell_around_point
  // is an alternative. That method is not
  // used here mostly because of the history
  // of the class. The algorithm used here
  // may be slightly more
  // efficient than find_active_cell_around_point
  // because it operates on a set of points.
  for (; cell != endc; cell++)
    {
      fe_values.reinit (cell);
      for (unsigned int support_point = 0; support_point < n_support_points; support_point++)
        {
          unsigned int component = dof_handler->get_fe ().system_to_component_index (support_point).first;
          Point<dim> test_point = fe_values.quadrature_point (support_point);

          for (unsigned int point = 0; point < locations.size (); point++)
            {
              if (locations[point].distance (test_point) < locations[point].distance (current_points[point][component]))
                {
                  // save the data
                  current_points[point][component] = test_point;
                  current_cell[point] = cell;
                  current_fe_index[point][component] = support_point;
                }
            }
        }
    }

  std::vector<types::global_dof_index> local_dof_indices (dof_handler->get_fe ().dofs_per_cell);
  for (unsigned int point = 0; point < locations.size (); point++)
    {
      current_cell[point]->get_dof_indices (local_dof_indices);
      std::vector<types::global_dof_index> new_solution_indices;

      for (unsigned int component = 0; component < dof_handler->get_fe ().n_components (); component++)
        {
          new_solution_indices.push_back (local_dof_indices[current_fe_index[point][component]]);
        }

      internal::PointValueHistory::PointGeometryData<dim> new_point_geometry_data (locations[point], current_points[point], new_solution_indices);

      point_geometry_data.push_back (new_point_geometry_data);

      std::map <std::string, std::vector <std::vector <double> > >::iterator
      data_store_begin = data_store.begin ();
      for (; data_store_begin != data_store.end (); data_store_begin++)
        {
          // add an extra row to each vector
          // entry
          const ComponentMask current_mask = (component_mask.find (data_store_begin->first))->second;
          unsigned int n_stored = current_mask.n_selected_components();
          for (unsigned int component = 0; component < n_stored; component++)
            {
              data_store_begin->second.push_back (std::vector<double> (0));
            }
        }
    }
}






template <int dim>
void PointValueHistory<dim>
::add_field_name (const std::string &vector_name,
                  const ComponentMask &mask)
{
  // can't be closed to add additional points
  // or vectors
  AssertThrow (!closed, ExcInvalidState ());
  AssertThrow (!cleared, ExcInvalidState ());
  AssertThrow (have_dof_handler, ExcDoFHandlerRequired ());
  AssertThrow (!triangulation_changed, ExcDoFHandlerChanged ());

  // insert a component mask that is always of the right size
  if (mask.represents_the_all_selected_mask() == false)
    component_mask.insert (std::make_pair (vector_name, mask));
  else
    component_mask.insert (std::make_pair (vector_name,
                                           ComponentMask(std::vector<bool>(dof_handler->get_fe().n_components(), true))));

  // insert an empty vector of strings
  // to ensure each field has an entry
  // in the map
  std::pair <std::string, std::vector <std::string> >
  empty_names (vector_name, std::vector <std::string> ());
  component_names_map.insert (empty_names);

  // make and add a new vector
  // point_geometry_data.size() long
  std::pair<std::string, std::vector <std::vector <double> > > pair_data;
  pair_data.first = vector_name;
  const unsigned int n_stored = (mask.represents_the_all_selected_mask() == false
                                 ?
                                 mask.n_selected_components()
                                 :
                                 dof_handler->get_fe().n_components());

  int n_datastreams = point_geometry_data.size () * n_stored; // each point has n_stored sub parts
  std::vector < std::vector <double> > vector_size (n_datastreams,
                                                    std::vector <double> (0));
  pair_data.second = vector_size;
  data_store.insert (pair_data);
}


template <int dim>
void PointValueHistory<dim>
::add_field_name(const std::string &vector_name, const unsigned int n_components)
{
  std::vector <bool> temp_mask (n_components, true);
  add_field_name (vector_name, temp_mask);
}


template <int dim>
void PointValueHistory<dim>
::add_component_names(const std::string &vector_name,
                      const std::vector <std::string> &component_names)
{
  typename std::map <std::string, std::vector <std::string> >::iterator names = component_names_map.find(vector_name);
  Assert (names != component_names_map.end(), ExcMessage("vector_name not in class"));

  typename std::map <std::string, ComponentMask>::iterator mask = component_mask.find(vector_name);
  Assert (mask != component_mask.end(), ExcMessage("vector_name not in class"));
  unsigned int n_stored = mask->second.n_selected_components();
  (void)n_stored;
  Assert (component_names.size() == n_stored, ExcDimensionMismatch (component_names.size(), n_stored));

  names->second = component_names;
}


template <int dim>
void PointValueHistory<dim>
::add_independent_names(const std::vector <std::string> &independent_names)
{
  Assert (independent_names.size() == n_indep, ExcDimensionMismatch (independent_names.size(), n_indep));

  indep_names = independent_names;
}


template <int dim>
void PointValueHistory<dim>
::close ()
{
  closed = true;
}



template <int dim>
void PointValueHistory<dim>
::clear ()
{
  cleared = true;
  dof_handler = 0;
  have_dof_handler = false;
}

// Need to test that the internal data has a full and complete dataset for
// each key. That is that the data has not got 'out of sync'. Testing that
// dataset_key is within 1 of independent_values is cheap and is done in all
// three methods. Evaluate_field will check that its vector_name is within 1
// of dataset_key. However this leaves the possibility that the user has
// neglected to call evaluate_field on one vector_name consistently. To catch
// this case start_new_dataset will call bool deap_check () which will test
// all vector_names and return a bool. This can be called from an Assert
// statement.



template <int dim>
template <typename VectorType>
void PointValueHistory<dim>
::evaluate_field (const std::string &vector_name, const VectorType &solution)
{
  // must be closed to add data to internal
  // members.
  Assert (closed, ExcInvalidState ());
  Assert (!cleared, ExcInvalidState ());
  AssertThrow (have_dof_handler, ExcDoFHandlerRequired ());
  AssertThrow (!triangulation_changed, ExcDoFHandlerChanged ());

  if (n_indep != 0) // hopefully this will get optimized, can't test independent_values[0] unless n_indep > 0
    {
      Assert (std::abs ((int) dataset_key.size () - (int) independent_values[0].size ()) < 2, ExcDataLostSync ());
    }
  // Look up the field name and get an
  // iterator for the map. Doing this
  // up front means that it only needs
  // to be done once and also allows us
  // to check vector_name is in the map.
  typename std::map <std::string, std::vector <std::vector <double> > >::iterator data_store_field = data_store.find(vector_name);
  Assert (data_store_field != data_store.end(), ExcMessage("vector_name not in class"));
  // Repeat for component_mask
  typename std::map <std::string, ComponentMask>::iterator mask = component_mask.find(vector_name);
  Assert (mask != component_mask.end(), ExcMessage("vector_name not in class"));

  unsigned int n_stored = mask->second.n_selected_components(dof_handler->get_fe ().n_components ());

  typename std::vector <internal::PointValueHistory::PointGeometryData <dim> >::iterator point = point_geometry_data.begin ();
  for (unsigned int data_store_index = 0; point != point_geometry_data.end (); point++, data_store_index++)
    {
      // Look up the components to add
      // in the component_mask, and
      // access the data associated with
      // those components

      for (unsigned int store_index = 0, comp = 0; comp < dof_handler->get_fe ().n_components (); comp++)
        {
          if (mask->second[comp])
            {
              unsigned int solution_index = point->solution_indices[comp];
              data_store_field->second[data_store_index * n_stored + store_index].push_back (solution (solution_index));
              store_index++;
            }
        }
    }
}





template <int dim>
template <typename VectorType>
void PointValueHistory<dim>
::evaluate_field (const std::vector <std::string> &vector_names,
                  const VectorType                &solution,
                  const DataPostprocessor< dim>   &data_postprocessor,
                  const Quadrature<dim>           &quadrature)
{
  // must be closed to add data to internal
  // members.
  Assert (closed, ExcInvalidState ());
  Assert (!cleared, ExcInvalidState ());
  AssertThrow (have_dof_handler, ExcDoFHandlerRequired ());
  if (n_indep != 0) // hopefully this will get optimized, can't test independent_values[0] unless n_indep > 0
    {
      Assert (std::abs ((int) dataset_key.size () - (int) independent_values[0].size ()) < 2, ExcDataLostSync ());
    }

  // Make an FEValues object
  const UpdateFlags update_flags = data_postprocessor.get_needed_update_flags() | update_quadrature_points;
  Assert (!(update_flags & update_normal_vectors),
          ExcMessage("The update of normal vectors may not be requested for evaluation of "
                     "data on cells via DataPostprocessor."));
  FEValues<dim> fe_values (dof_handler->get_fe (), quadrature, update_flags);
  unsigned int n_components = dof_handler->get_fe ().n_components ();
  unsigned int n_quadrature_points = quadrature.size();

  unsigned int n_output_variables = data_postprocessor.get_names().size();

  // Loop over points and find correct cell
  typename std::vector <internal::PointValueHistory::PointGeometryData <dim> >::iterator point = point_geometry_data.begin ();
  for (unsigned int data_store_index = 0; point != point_geometry_data.end (); point++, data_store_index++)
    {
      // we now have a point to query,
      // need to know what cell it is in
      Point <dim> requested_location = point->requested_location;
      typename DoFHandler<dim>::active_cell_iterator cell = GridTools::find_active_cell_around_point (StaticMappingQ1<dim>::mapping, *dof_handler, requested_location).first;


      fe_values.reinit (cell);
      std::vector< Vector< double > > computed_quantities (1, Vector <double> (n_output_variables)); // just one point needed

      // The case of a scalar FE
      if (n_components == 1)
        {
          // Extract data for the
          // PostProcessor object
          std::vector< typename VectorType::value_type > uh (n_quadrature_points, 0.0);
          std::vector< Tensor< 1, dim, typename VectorType::value_type > > duh (n_quadrature_points, Tensor <1, dim, typename VectorType::value_type> ());
          std::vector< Tensor< 2, dim, typename VectorType::value_type > > dduh (n_quadrature_points, Tensor <2, dim, typename VectorType::value_type> ());
          std::vector<Point<dim> > dummy_normals (1, Point<dim> ());
          std::vector<Point<dim> > evaluation_points;
          // at each point there is
          // only one component of
          // value, gradient etc.
          if (update_flags & update_values)
            fe_values.get_function_values (solution,
                                           uh);
          if (update_flags & update_gradients)
            fe_values.get_function_gradients (solution,
                                              duh);
          if (update_flags & update_hessians)
            fe_values.get_function_hessians (solution,
                                             dduh);

          // find the closest quadrature point
          evaluation_points = fe_values.get_quadrature_points();
          double distance = cell->diameter ();
          unsigned int selected_point = 0;
          for (unsigned int q_point = 0; q_point < n_quadrature_points; q_point++)
            {
              if (requested_location.distance (evaluation_points[q_point]) < distance)
                {
                  selected_point = q_point;
                  distance = requested_location.distance (evaluation_points[q_point]);
                }
            }

          // Call compute_derived_quantities_vector
          // or compute_derived_quantities_scalar
          // TODO this function should also operate with typename VectorType::value_type
          data_postprocessor.
          compute_derived_quantities_scalar(std::vector< double > (1, uh[selected_point]),
                                            std::vector< Tensor< 1, dim > > (1, Tensor< 1, dim >(duh[selected_point]) ),
                                            std::vector< Tensor< 2, dim > > (1, Tensor< 2, dim >(dduh[selected_point]) ),
                                            dummy_normals,
                                            std::vector<Point<dim> > (1, evaluation_points[selected_point]),
                                            computed_quantities);

        }
      else     // The case of a vector FE
        {
          // Extract data for the PostProcessor object
          std::vector< Vector< typename VectorType::value_type > > uh (n_quadrature_points, Vector <typename VectorType::value_type> (n_components));
          std::vector< std::vector< Tensor< 1, dim, typename VectorType::value_type > > > duh (n_quadrature_points, std::vector< Tensor< 1, dim, typename VectorType::value_type > > (n_components,  Tensor< 1, dim, typename VectorType::value_type >()));
          std::vector< std::vector< Tensor< 2, dim, typename VectorType::value_type > > > dduh (n_quadrature_points, std::vector< Tensor< 2, dim, typename VectorType::value_type > > (n_components,  Tensor< 2, dim, typename VectorType::value_type >()));
          std::vector<Point<dim> > dummy_normals  (1, Point<dim> ());
          std::vector<Point<dim> > evaluation_points;
          // at each point there is
          // a vector valued
          // function and its
          // derivative...
          if (update_flags & update_values)
            fe_values.get_function_values (solution,
                                           uh);
          if (update_flags & update_gradients)
            fe_values.get_function_gradients (solution,
                                              duh);
          if (update_flags & update_hessians)
            fe_values.get_function_hessians (solution,
                                             dduh);

          // find the closest quadrature point
          evaluation_points = fe_values.get_quadrature_points();
          double distance = cell->diameter ();
          unsigned int selected_point = 0;
          for (unsigned int q_point = 0; q_point < n_quadrature_points; q_point++)
            {
              if (requested_location.distance (evaluation_points[q_point]) < distance)
                {
                  selected_point = q_point;
                  distance = requested_location.distance (evaluation_points[q_point]);
                }
            }

          // FIXME: We need tmp vectors below because the data
          // postprocessors are not equipped to deal with anything but
          // doubles (scalars and tensors).
          const Vector< typename VectorType::value_type >                        &uh_s   = uh[selected_point];
          const std::vector< Tensor< 1, dim, typename VectorType::value_type > > &duh_s  = duh[selected_point];
          const std::vector< Tensor< 2, dim, typename VectorType::value_type > > &dduh_s = dduh[selected_point];
          std::vector< Tensor< 1, dim > > tmp_d (duh_s.size());
          for (unsigned int i = 0; i < duh_s.size(); i++)
            tmp_d[i] = duh_s[i];

          std::vector< Tensor< 2, dim > > tmp_dd (dduh_s.size());
          for (unsigned int i = 0; i < dduh_s.size(); i++)
            tmp_dd[i] = dduh_s[i];

          Vector< double > tmp(uh_s.size());
          for (unsigned int i = 0; i < uh_s.size(); i++)
            tmp[i] = uh_s[i];
          // Call compute_derived_quantities_vector
          // or compute_derived_quantities_scalar
          data_postprocessor.
          compute_derived_quantities_vector(std::vector< Vector< double > > (1, tmp),
                                            std::vector< std::vector< Tensor< 1, dim > > > (1, tmp_d),
                                            std::vector< std::vector< Tensor< 2, dim > > > (1, tmp_dd),
                                            dummy_normals,
                                            std::vector<Point<dim> > (1, evaluation_points[selected_point]),
                                            computed_quantities);
        }


      // we now have the data and need to save it
      // loop over data names
      typename std::vector<std::string>::const_iterator name = vector_names.begin();
      for (; name != vector_names.end(); name++)
        {
          typename std::map <std::string, std::vector <std::vector <double> > >::iterator data_store_field = data_store.find(*name);
          Assert (data_store_field != data_store.end(), ExcMessage("vector_name not in class"));
          // Repeat for component_mask
          typename std::map <std::string, ComponentMask>::iterator mask = component_mask.find(*name);
          Assert (mask != component_mask.end(), ExcMessage("vector_name not in class"));

          unsigned int n_stored = mask->second.n_selected_components(n_output_variables);

          // Push back computed quantities according
          // to the component_mask.
          for (unsigned int store_index = 0, comp = 0; comp < n_output_variables; comp++)
            {
              if (mask->second[comp])
                {
                  data_store_field->second[data_store_index * n_stored + store_index].push_back (computed_quantities[0](comp));
                  store_index++;
                }
            }
        }
    } // end of loop over points
}


template <int dim>
template <typename VectorType>
void PointValueHistory<dim>
::evaluate_field (const std::string            &vector_name,
                  const VectorType             &solution,
                  const DataPostprocessor<dim> &data_postprocessor,
                  const Quadrature<dim>        &quadrature)
{
  std::vector <std::string> vector_names;
  vector_names.push_back (vector_name);
  evaluate_field (vector_names, solution, data_postprocessor, quadrature);
}



template <int dim>
template <typename VectorType>
void PointValueHistory<dim>
::evaluate_field_at_requested_location (const std::string &vector_name,
                                        const VectorType  &solution)
{
  typedef typename VectorType::value_type number;
  // must be closed to add data to internal
  // members.
  Assert (closed, ExcInvalidState ());
  Assert (!cleared, ExcInvalidState ());
  AssertThrow (have_dof_handler, ExcDoFHandlerRequired ());

  if (n_indep != 0) // hopefully this will get optimized, can't test independent_values[0] unless n_indep > 0
    {
      Assert (std::abs ((int) dataset_key.size () - (int) independent_values[0].size ()) < 2, ExcDataLostSync ());
    }
  // Look up the field name and get an
  // iterator for the map. Doing this
  // up front means that it only needs
  // to be done once and also allows us
  // to check vector_name is in the map.
  typename std::map <std::string, std::vector <std::vector <double> > >::iterator data_store_field = data_store.find(vector_name);
  Assert (data_store_field != data_store.end(), ExcMessage("vector_name not in class"));
  // Repeat for component_mask
  typename std::map <std::string, ComponentMask>::iterator mask = component_mask.find(vector_name);
  Assert (mask != component_mask.end(), ExcMessage("vector_name not in class"));

  unsigned int n_stored = mask->second.n_selected_components(dof_handler->get_fe ().n_components ());

  typename std::vector <internal::PointValueHistory::PointGeometryData <dim> >::iterator point = point_geometry_data.begin ();
  Vector <number> value (dof_handler->get_fe().n_components());
  for (unsigned int data_store_index = 0; point != point_geometry_data.end (); point++, data_store_index++)
    {
      // Make a Vector <double> for the value
      // at the point. It will have as many
      // components as there are in the fe.
      VectorTools::point_value (*dof_handler, solution, point->requested_location, value);

      // Look up the component_mask and add
      // in components according to that mask
      for (unsigned int store_index = 0, comp = 0; comp < mask->second.size(); comp++)
        {
          if (mask->second[comp])
            {
              data_store_field->second[data_store_index * n_stored + store_index].push_back (value (comp));
              store_index++;
            }
        }
    }
}


template <int dim>
void PointValueHistory<dim>
::start_new_dataset (double key)
{
  // must be closed to add data to internal
  // members.
  Assert (closed, ExcInvalidState ());
  Assert (!cleared, ExcInvalidState ());
  Assert (deep_check (false), ExcDataLostSync ());

  dataset_key.push_back (key);
}



template <int dim>
void PointValueHistory<dim>
::push_back_independent (const std::vector <double> &indep_values)
{
  // must be closed to add data to internal
  // members.
  Assert (closed, ExcInvalidState ());
  Assert (!cleared, ExcInvalidState ());
  Assert (indep_values.size () == n_indep, ExcDimensionMismatch (indep_values.size (), n_indep));
  Assert (n_indep != 0, ExcNoIndependent ());
  Assert (std::abs ((int) dataset_key.size () - (int) independent_values[0].size ()) < 2, ExcDataLostSync ());

  for (unsigned int component = 0; component < n_indep; component++)
    independent_values[component].push_back (indep_values[component]);
}



template <int dim>
void PointValueHistory<dim>
::write_gnuplot (const std::string &base_name, const std::vector <Point <dim> > postprocessor_locations)
{
  AssertThrow (closed, ExcInvalidState ());
  AssertThrow (!cleared, ExcInvalidState ());
  AssertThrow (deep_check (true), ExcDataLostSync ());

  // write inputs to a file
  if (n_indep != 0)
    {
      std::string filename = base_name + "_indep.gpl";
      std::ofstream to_gnuplot (filename.c_str ());

      to_gnuplot << "# Data independent of mesh location\n";

      // write column headings
      to_gnuplot << "# <Key> ";

      if (indep_names.size() > 0)
        {
          for (unsigned int name = 0; name < indep_names.size(); name++)
            {
              to_gnuplot << "<" << indep_names [name] << "> ";
            }
          to_gnuplot << "\n";
        }
      else
        {
          for (unsigned int component = 0; component < n_indep; component++)
            {
              to_gnuplot << "<Indep_" << component << "> ";
            }
          to_gnuplot << "\n";
        }
      // write general data stored
      for (unsigned int key = 0; key < dataset_key.size (); key++)
        {
          to_gnuplot << dataset_key[key];

          for (unsigned int component = 0; component < n_indep; component++)
            {
              to_gnuplot << " " << independent_values[component][key];
            }
          to_gnuplot << "\n";
        }

      to_gnuplot.close ();
    }



  // write points to a file
  if (have_dof_handler)
    {
      AssertThrow (have_dof_handler, ExcDoFHandlerRequired ());
      AssertThrow (postprocessor_locations.size() == 0 || postprocessor_locations.size() == point_geometry_data.size(), ExcDimensionMismatch (postprocessor_locations.size(), point_geometry_data.size()));
      // We previously required the
      // number of dofs to remain the
      // same to provide some sort of
      // test on the relevance of the
      // support point indices stored.
      // We now relax that to allow
      // adaptive refinement strategies
      // to make use of the
      // evaluate_field_requested_locations
      // method. Note that the support point
      // information is not meaningful if
      // the number of dofs has changed.
      //AssertThrow (!triangulation_changed, ExcDoFHandlerChanged ());

      typename std::vector <internal::PointValueHistory::PointGeometryData <dim> >::iterator point = point_geometry_data.begin ();
      for (unsigned int data_store_index = 0; point != point_geometry_data.end (); point++, data_store_index++)
        {
          // for each point, open a file to
          // be written to
          std::string filename = base_name + "_" + Utilities::int_to_string (data_store_index, 2) + ".gpl"; // store by order pushed back
          // due to
          // Utilities::int_to_string(data_store_index,
          // 2) call, can handle up to 100
          // points
          std::ofstream to_gnuplot (filename.c_str ());

          // put helpful info about the
          // support point into the file as
          // comments
          to_gnuplot << "# Requested location: " << point->requested_location << "\n";
          to_gnuplot << "# DoF_index : Support location (for each component)\n";
          for (unsigned int component = 0; component < dof_handler->get_fe ().n_components (); component++)
            {
              to_gnuplot << "# " << point->solution_indices[component] << " : " << point->support_point_locations [component] << "\n";
            }
          if (triangulation_changed)
            to_gnuplot << "# (Original components and locations, may be invalidated by mesh change.)\n";

          if (postprocessor_locations.size() != 0)
            {
              to_gnuplot << "# Postprocessor location: " << postprocessor_locations[data_store_index];
              if (triangulation_changed)
                to_gnuplot << " (may be approximate)\n";
            }
          to_gnuplot << "#\n";


          // write column headings
          to_gnuplot << "# <Key> ";

          if (indep_names.size() > 0)
            {
              for (unsigned int name = 0; name < indep_names.size(); name++)
                {
                  to_gnuplot << "<" << indep_names [name] << "> ";
                }
            }
          else
            {
              for (unsigned int component = 0; component < n_indep; component++)
                {
                  to_gnuplot << "<Indep_" << component << "> ";
                }
            }

          for (std::map <std::string, std::vector <std::vector <double> > >::iterator
               data_store_begin = data_store.begin (); data_store_begin != data_store.end (); ++data_store_begin)
            {
              typename std::map <std::string, ComponentMask>::iterator mask = component_mask.find(data_store_begin->first);
              unsigned int n_stored = mask->second.n_selected_components();
              std::vector <std::string> names = (component_names_map.find (data_store_begin->first))->second;

              if (names.size() > 0)
                {
                  AssertThrow (names.size() == n_stored, ExcDimensionMismatch (names.size(), n_stored));
                  for (unsigned int component = 0; component < names.size(); component++)
                    {
                      to_gnuplot << "<" << names[component] << "> ";
                    }
                }
              else
                {
                  for (unsigned int component = 0; component < n_stored; component++)
                    {
                      to_gnuplot << "<" << data_store_begin->first << "_" << component << "> ";
                    }
                }
            }
          to_gnuplot << "\n";

          // write data stored for the point
          for (unsigned int key = 0; key < dataset_key.size (); key++)
            {
              to_gnuplot << dataset_key[key];

              for (unsigned int component = 0; component < n_indep; component++)
                {
                  to_gnuplot << " " << independent_values[component][key];
                }

              for (std::map <std::string, std::vector <std::vector <double> > >::iterator
                   data_store_begin = data_store.begin ();
                   data_store_begin != data_store.end (); ++data_store_begin)
                {
                  typename std::map <std::string, ComponentMask>::iterator mask = component_mask.find(data_store_begin->first);
                  unsigned int n_stored = mask->second.n_selected_components();

                  for (unsigned int component = 0; component < n_stored; component++)
                    {
                      to_gnuplot << " " << (data_store_begin->second)[data_store_index * n_stored + component][key];
                    }
                }
              to_gnuplot << "\n";
            }

          to_gnuplot.close ();
        }
    }
}



template <int dim>
Vector<double> PointValueHistory<dim>
::mark_support_locations ()
{
  // a method to put a one at each point on
  // the grid where a location is defined
  AssertThrow (!cleared, ExcInvalidState ());
  AssertThrow (have_dof_handler, ExcDoFHandlerRequired ());
  AssertThrow (!triangulation_changed, ExcDoFHandlerChanged ());

  Vector<double> dof_vector (dof_handler->n_dofs ());

  typename std::vector <internal::PointValueHistory::PointGeometryData <dim> >::iterator point = point_geometry_data.begin ();
  for (; point != point_geometry_data.end (); point++)
    {
      for (unsigned int component = 0; component < dof_handler->get_fe ().n_components (); component++)
        {
          dof_vector (point->solution_indices[component]) = 1;
        }
    }
  return dof_vector;
}


template <int dim>
void PointValueHistory<dim>
::get_support_locations (std::vector <std::vector<Point <dim> > > &locations)
{
  AssertThrow (!cleared, ExcInvalidState ());
  AssertThrow (have_dof_handler, ExcDoFHandlerRequired ());
  AssertThrow (!triangulation_changed, ExcDoFHandlerChanged ());

  std::vector <std::vector <Point <dim> > > actual_points;
  typename std::vector <internal::PointValueHistory::PointGeometryData <dim> >::iterator point = point_geometry_data.begin ();

  for (; point != point_geometry_data.end (); point++)
    {
      actual_points.push_back (point->support_point_locations);
    }
  locations = actual_points;
}


template <int dim>
void PointValueHistory<dim>
::get_points (std::vector <std::vector<Point <dim> > > &locations)
{
  get_support_locations (locations);
}


template <int dim>
void PointValueHistory<dim>
::get_postprocessor_locations (const Quadrature<dim> &quadrature, std::vector<Point <dim> > &locations)
{
  Assert (!cleared, ExcInvalidState ());
  AssertThrow (have_dof_handler, ExcDoFHandlerRequired ());

  locations = std::vector<Point <dim> > ();

  FEValues<dim> fe_values (dof_handler->get_fe (), quadrature, update_quadrature_points);
  unsigned int n_quadrature_points = quadrature.size();
  std::vector<Point<dim> > evaluation_points;

  // Loop over points and find correct cell
  typename std::vector <internal::PointValueHistory::PointGeometryData <dim> >::iterator point = point_geometry_data.begin ();
  for (unsigned int data_store_index = 0; point != point_geometry_data.end (); point++, data_store_index++)
    {
      // we now have a point to query,
      // need to know what cell it is in
      Point <dim> requested_location = point->requested_location;
      typename DoFHandler<dim>::active_cell_iterator cell = GridTools::find_active_cell_around_point (StaticMappingQ1<dim>::mapping, *dof_handler, requested_location).first;
      fe_values.reinit (cell);

      evaluation_points = fe_values.get_quadrature_points();
      double distance = cell->diameter ();
      unsigned int selected_point = 0;

      for (unsigned int q_point = 0; q_point < n_quadrature_points; q_point++)
        {
          if (requested_location.distance (evaluation_points[q_point]) < distance)
            {
              selected_point = q_point;
              distance = requested_location.distance (evaluation_points[q_point]);
            }
        }

      locations.push_back (evaluation_points[selected_point]);
    }
}


template <int dim>
void PointValueHistory<dim>
::status (std::ostream &out)
{
  out << "***PointValueHistory status output***\n\n";
  out << "Closed: " << closed << "\n";
  out << "Cleared: " << cleared << "\n";
  out << "Triangulation_changed: " << triangulation_changed << "\n";
  out << "Have_dof_handler: " << have_dof_handler << "\n";
  out << "Geometric Data" << "\n";

  typename std::vector <internal::PointValueHistory::PointGeometryData <dim> >::iterator point = point_geometry_data.begin ();
  if (point == point_geometry_data.end ())
    {
      out << "No points stored currently\n";
    }
  else
    {
      if (!cleared)
        {
          for (; point != point_geometry_data.end (); point++)
            {
              out << "# Requested location: " << point->requested_location << "\n";
              out << "# DoF_index : Support location (for each component)\n";
              for (unsigned int component = 0; component < dof_handler->get_fe ().n_components (); component++)
                {
                  out << point->solution_indices[component] << " : " << point->support_point_locations [component] << "\n";
                }
              out << "\n";
            }
        }
      else
        {
          out << "#Cannot access DoF_indices once cleared\n";
        }
    }
  out << "\n";

  if (independent_values.size () != 0)
    {
      out << "Independent value(s): " << independent_values.size () << " : " << independent_values[0].size () << "\n";
      if (indep_names.size() > 0)
        {
          out << "Names: ";
          for (unsigned int name = 0; name < indep_names.size(); name++)
            {
              out << "<" << indep_names [name] << "> ";
            }
          out << "\n";
        }
    }
  else
    {
      out << "No independent values stored\n";
    }

  std::map <std::string, std::vector <std::vector <double> > >::iterator
  data_store_begin = data_store.begin ();
  if (data_store_begin != data_store.end())
    {
      out << "Mnemonic: data set size (mask size, n true components) : n data sets\n";
    }
  for (; data_store_begin != data_store.end (); data_store_begin++)
    {
      // Find field mnemonic
      std::string vector_name = data_store_begin->first;
      typename std::map <std::string, ComponentMask>::iterator mask = component_mask.find(vector_name);
      Assert (mask != component_mask.end(), ExcMessage("vector_name not in class"));
      typename std::map <std::string, std::vector <std::string> >::iterator component_names = component_names_map.find(vector_name);
      Assert (component_names != component_names_map.end(), ExcMessage("vector_name not in class"));

      if (data_store_begin->second.size () != 0)
        {
          out << data_store_begin->first << ": " << data_store_begin->second.size () << " (";
          out << mask->second.size() << ", " << mask->second.n_selected_components() << ") : ";
          out << (data_store_begin->second)[0].size () << "\n";
        }
      else
        {
          out << data_store_begin->first << ": " << data_store_begin->second.size () << " (";
          out << mask->second.size() << ", " << mask->second.n_selected_components() << ") : ";
          out << "No points added" << "\n";
        }
      // add names, if available
      if (component_names->second.size() > 0)
        {
          for (unsigned int name = 0; name < component_names->second.size(); name++)
            {
              out << "<" << component_names->second[name] << "> ";
            }
          out << "\n";
        }
    }
  out << "\n";
  out << "***end of status output***\n\n";
}



template <int dim>
bool PointValueHistory<dim>
::deep_check (const bool strict)
{
  // test ways that it can fail, if control
  // reaches last statement return true
  if (strict)
    {
      if (n_indep != 0)
        {
          if (dataset_key.size () != independent_values[0].size ())
            {
              return false;
            }
        }
      std::map <std::string, std::vector <std::vector <double> > >::iterator
      data_store_begin = data_store.begin ();
      if (have_dof_handler)
        {
          for (; data_store_begin != data_store.end (); data_store_begin++)
            {
              Assert (data_store_begin->second.size() > 0,
                      ExcInternalError());
              if ((data_store_begin->second)[0].size () != dataset_key.size ())
                return false;
              // this loop only tests one
              // member for each name,
              // i.e. checks the user it will
              // not catch internal errors
              // which do not update all
              // fields for a name.
            }
        }
      return true;
    }
  if (n_indep != 0)
    {
      if (std::abs ((int) dataset_key.size () - (int) independent_values[0].size ()) >= 2)
        {
          return false;
        }
    }

  if (have_dof_handler)
    {
      std::map <std::string, std::vector <std::vector <double> > >::iterator
      data_store_begin = data_store.begin ();
      for (; data_store_begin != data_store.end (); data_store_begin++)
        {
          Assert (data_store_begin->second.size() > 0,
                  ExcInternalError());

          if (std::abs ((int) (data_store_begin->second)[0].size () - (int) dataset_key.size ()) >= 2)
            return false;
          // this loop only tests one member
          // for each name, i.e. checks the
          // user it will not catch internal
          // errors which do not update all
          // fields for a name.
        }
    }
  return true;
}



template <int dim>
void PointValueHistory<dim>
::tria_change_listener ()
{
  // this function is called by the
  // Triangulation whenever something
  // changes, by virtue of having
  // attached the function to the
  // signal handler in the
  // triangulation object

  // we record the fact that the mesh
  // has changed. we need to take
  // this into account next time we
  // evaluate the solution
  triangulation_changed = true;
}


// explicit instantiations
#include "point_value_history.inst"


DEAL_II_NAMESPACE_CLOSE

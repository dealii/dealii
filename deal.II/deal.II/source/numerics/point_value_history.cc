//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009 by Michael Rapson and the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <numerics/point_value_history.h>



DEAL_II_NAMESPACE_OPEN

#include "point_value_history.inst"

namespace internal
{
  namespace PointValueHistory
  {
/// Only a constructor needed for this class (a struct really)
    template <int dim>
    PointGeometryData<dim>
    ::PointGeometryData (const std::vector <Point <dim> > &new_locations,
			 const std::vector <int> &new_sol_indices)
    {
      locations = new_locations;
      solution_indices = new_sol_indices;
    }
  }
}



template <int dim>
PointValueHistory<dim>
::PointValueHistory (const unsigned int n_independent_variables) :
		dof_handler (0, "No DoFHandler"),
		n_dofs (0),
		n_indep (n_independent_variables)
{
  closed = false;
  cleared = false;
				   // make a vector for keys
  dataset_key = std::vector <double> (); // initialize the std::vector

				   // make a vector of independent values
  independent_values
    = std::vector<std::vector <double> > (n_indep, std::vector <double> (0));
}



template <int dim>
PointValueHistory<dim>
::PointValueHistory (const DoFHandler<dim> & dof_handler,
		     const unsigned int n_independent_variables) :
		dof_handler (&dof_handler, typeid (*this).name ()),
		n_dofs (dof_handler.n_dofs ()),
		n_indep (n_independent_variables)
{
  closed = false;
  cleared = false;
				   // make a vector to store keys
  dataset_key = std::vector <double> (); // initialize the std::vector

				   // make a vector for the independent values
  independent_values
    = std::vector<std::vector <double> > (n_indep, std::vector <double> (0));
}



template <int dim>
PointValueHistory<dim>
::PointValueHistory (const PointValueHistory & point_value_history)
{
  dataset_key = point_value_history.dataset_key;
  independent_values = point_value_history.independent_values;
  data_store = point_value_history.data_store;
  point_geometry_data = point_value_history.point_geometry_data;
  pair_data = point_value_history.pair_data;
  closed = point_value_history.closed;
  cleared = point_value_history.cleared;

  dof_handler = point_value_history.dof_handler;

  n_dofs = point_value_history.n_dofs;
  n_indep = point_value_history.n_indep;
}



template <int dim>
PointValueHistory<dim> &
PointValueHistory<dim>::operator= (const PointValueHistory & point_value_history)
{
  dataset_key = point_value_history.dataset_key;
  independent_values = point_value_history.independent_values;
  data_store = point_value_history.data_store;
  point_geometry_data = point_value_history.point_geometry_data;
  pair_data = point_value_history.pair_data;
  closed = point_value_history.closed;
  cleared = point_value_history.cleared;

  dof_handler = point_value_history.dof_handler;

  n_dofs = point_value_history.n_dofs;
  n_indep = point_value_history.n_indep;
  return * this;
}



template <int dim>
PointValueHistory<dim>
::~PointValueHistory ()
{
}



template <int dim>
void PointValueHistory<dim>
::add_point (const Point <dim> & location)
{
				   // can't be closed to add additional points
				   // or vectors
  AssertThrow (!closed, ExcInvalidState ());
  AssertThrow (!cleared, ExcInvalidState ());
  AssertThrow (n_dofs != 0, ExcDoFHandlerRequired ());
  AssertThrow (n_dofs == dof_handler->n_dofs (),
	       ExcDoFHandlerChanged (n_dofs, dof_handler->n_dofs ()));

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


  std::vector<unsigned int>
    local_dof_indices (dof_handler->get_fe ().dofs_per_cell);
  std::vector <int> new_solution_indices;
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
				   // to obtain a cell to search may be an
				   // option for these methods, but currently
				   // the GridTools method does not cater for
				   // a vector of points, and does not seem to
				   // be intrinsicly faster than this method.
  for (unsigned int component = 0;
       component < dof_handler->get_fe ().n_components (); component++)
    {
      new_solution_indices
	.push_back (local_dof_indices[current_fe_index [component]]);
    }

  internal::PointValueHistory::PointGeometryData<dim>
    new_point_geometry_data (current_points, new_solution_indices);
  point_geometry_data.push_back (new_point_geometry_data);

  std::map <std::string, std::vector <std::vector <double> > >::iterator
    data_store_begin = data_store.begin ();
  for (; data_store_begin != data_store.end (); data_store_begin++)
    {
				       // add an extra row to each vector
				       // entry
      for (unsigned int component = 0;
	   component < dof_handler->get_fe ().n_components (); component++)
        {
          data_store_begin->second.push_back (std::vector<double> (0));
        }
    }
}



template <int dim>
void PointValueHistory<dim>
::add_points (const std::vector <Point <dim> > & locations)
{
				   // This algorithm adds points in the same
				   // order as they appear in the vector
				   // locations and users may depend on this
				   // so do not change order added!

				   // can't be closed to add additional points or vectors
  AssertThrow (!closed, ExcInvalidState ());
  AssertThrow (!cleared, ExcInvalidState ());
  AssertThrow (n_dofs != 0, ExcDoFHandlerRequired ());
  AssertThrow (n_dofs == dof_handler->n_dofs (), ExcDoFHandlerChanged (n_dofs, dof_handler->n_dofs ()));
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

  std::vector<unsigned int> local_dof_indices (dof_handler->get_fe ().dofs_per_cell);
  for (unsigned int point = 0; point < locations.size (); point++)
    {
      current_cell[point]->get_dof_indices (local_dof_indices);
      std::vector <int> new_solution_indices;

      for (unsigned int component = 0; component < dof_handler->get_fe ().n_components (); component++)
        {
          new_solution_indices.push_back (local_dof_indices[current_fe_index[point][component]]);
        }

      internal::PointValueHistory::PointGeometryData<dim> new_point_geometry_data (current_points[point], new_solution_indices);

      point_geometry_data.push_back (new_point_geometry_data);

      std::map <std::string, std::vector <std::vector <double> > >::iterator
	data_store_begin = data_store.begin ();
      for (; data_store_begin != data_store.end (); data_store_begin++)
        {
					   // add an extra row to each vector
					   // entry
          for (unsigned int component = 0; component < dof_handler->get_fe ().n_components (); component++)
            {
              data_store_begin->second.push_back (std::vector<double> (0));
            }
        }
    }
}






template <int dim>
void PointValueHistory<dim>
::add_field_name (const std::string &vector_name)
{
				   // can't be closed to add additional points
				   // or vectors
  AssertThrow (!closed, ExcInvalidState ());
  AssertThrow (n_dofs != 0, ExcDoFHandlerRequired ());
  AssertThrow (!cleared, ExcInvalidState ());
  AssertThrow (n_dofs == dof_handler->n_dofs (), ExcDoFHandlerChanged (n_dofs, dof_handler->n_dofs ()));


				   // make and add a new vector
				   // point_geometry_data.size() long
  pair_data.first = vector_name;
  int n_datastreams = point_geometry_data.size () * (dof_handler->get_fe ().n_components ()); // each point has n_components sub parts
  std::vector < std::vector <double> > vector_size (n_datastreams,
                                                    std::vector <double> (0));
  pair_data.second = vector_size;
  data_store.insert (pair_data);
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
  dof_handler = SmartPointer<const DoFHandler<dim> > (0, "Cleared");
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
template <class VECTOR>
void PointValueHistory<dim>
::evaluate_field (const std::string &vector_name, const VECTOR & solution)
{
				   // must be closed to add data to internal
				   // members.
  Assert (closed, ExcInvalidState ());
  Assert (!cleared, ExcInvalidState ());
  Assert (n_dofs != 0, ExcDoFHandlerRequired ());
  Assert (n_dofs == dof_handler->n_dofs (), ExcDoFHandlerChanged (n_dofs, dof_handler->n_dofs ()));
  if (n_indep != 0) // hopefully this will get optimized, can't test independent_values[0] unless n_indep > 0
    {
      Assert (std::abs ((int) dataset_key.size () - (int) independent_values[0].size ()) < 2, ExcDataLostSync ());
    }
  typename std::vector <internal::PointValueHistory::PointGeometryData <dim> >::iterator node = point_geometry_data.begin ();
  for (unsigned int data_store_index = 0; node != point_geometry_data.end (); node++, data_store_index++)
    {
				       // step through each node, to access
				       // the Solution_index and the vector
				       // index
      for (unsigned int component = 0; component < dof_handler->get_fe ().n_components (); component++)
        {
          unsigned int solution_index = node->solution_indices[component];
          (data_store[vector_name])[data_store_index * dof_handler->get_fe ().n_components () + component].push_back (solution (solution_index));
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
::write_gnuplot (const std::string &base_name)
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
      std::map <std::string, std::vector <std::vector <double> > >::iterator
	data_store_begin = data_store.begin ();
      to_gnuplot << "# <Key> ";

      for (unsigned int component = 0; component < n_indep; component++)
        {
          to_gnuplot << "<Indep_" << component << "> ";
        }
      to_gnuplot << "\n";

				       // write general data stored
      for (unsigned int key = 0; key < dataset_key.size (); key++)
        {
          std::map <std::string, std::vector <std::vector <double> > >::iterator
	    data_store_begin = data_store.begin ();
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
  if (n_dofs != 0)
    {
      AssertThrow (n_dofs != 0, ExcDoFHandlerRequired ());
      AssertThrow (n_dofs == dof_handler->n_dofs (), ExcDoFHandlerChanged (n_dofs, dof_handler->n_dofs ()));

      unsigned int n_components = dof_handler->get_fe ().n_components ();

      typename std::vector <internal::PointValueHistory::PointGeometryData <dim> >::iterator node = point_geometry_data.begin ();
      for (unsigned int data_store_index = 0; node != point_geometry_data.end (); node++, data_store_index++)
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
          to_gnuplot << "# DoF_index : Location (for each component)\n";
          for (unsigned int component = 0; component < n_components; component++)
            {
              to_gnuplot << "# " << node->solution_indices[component] << " : " << node->locations [component] << "\n";
            }
          to_gnuplot << "\n";

					   // write column headings
          std::map <std::string, std::vector <std::vector <double> > >::iterator
	    data_store_begin = data_store.begin ();
          to_gnuplot << "# <Key> ";

          for (unsigned int component = 0; component < n_indep; component++)
            {
              to_gnuplot << "<Indep_" << component << "> ";
            }

          for (; data_store_begin != data_store.end (); data_store_begin++)
            {
              for (unsigned int component = 0; component < n_components; component++)
                {
                  to_gnuplot << "<" << data_store_begin->first << "_" << component << "> ";
                }
            }
          to_gnuplot << "\n";

					   // write data stored for the node
          for (unsigned int key = 0; key < dataset_key.size (); key++)
            {
              std::map <std::string, std::vector <std::vector <double> > >::iterator
		data_store_begin = data_store.begin ();
              to_gnuplot << dataset_key[key];

              for (unsigned int component = 0; component < n_indep; component++)
                {
                  to_gnuplot << " " << independent_values[component][key];
                }

              for (; data_store_begin != data_store.end (); data_store_begin++)
                {
                  for (unsigned int component = 0; component < n_components; component++)
                    {
                      to_gnuplot << " " << (data_store_begin->second)[data_store_index * n_components + component][key];
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
::mark_locations ()
{
				   // a method to put a one at each point on
				   // the grid where a location is defined
  AssertThrow (n_dofs != 0, ExcDoFHandlerRequired ());
  AssertThrow (!cleared, ExcInvalidState ());
  AssertThrow (n_dofs == dof_handler->n_dofs (), ExcDoFHandlerChanged (n_dofs, dof_handler->n_dofs ()));

  Vector<double> dof_vector (dof_handler->n_dofs ());

  typename std::vector <internal::PointValueHistory::PointGeometryData <dim> >::iterator node = point_geometry_data.begin ();
  for (; node != point_geometry_data.end (); node++)
    {
      for (unsigned int component = 0; component < dof_handler->get_fe ().n_components (); component++)
        {
          dof_vector (node->solution_indices[component]) = 1;
        }
    }
  return dof_vector;
}



template <int dim>
void PointValueHistory<dim>
::get_points (std::vector <std::vector<Point <dim> > > & locations)
{
  AssertThrow (n_dofs != 0, ExcDoFHandlerRequired ());
  AssertThrow (!cleared, ExcInvalidState ());
  AssertThrow (n_dofs == dof_handler->n_dofs (), ExcDoFHandlerChanged (n_dofs, dof_handler->n_dofs ()));

  std::vector <std::vector <Point <dim> > > actual_points;
  typename std::vector <internal::PointValueHistory::PointGeometryData <dim> >::iterator node = point_geometry_data.begin ();

  for (; node != point_geometry_data.end (); node++)
    {
      actual_points.push_back (node->locations);
    }
  locations = actual_points;
}



template <int dim>
void PointValueHistory<dim>
::status (std::ostream &out)
{
  out << "***PointValueHistory status output***\n\n";
  out << "Closed: " << closed << "\n";
  out << "Cleared: " << cleared << "\n";
  out << "Geometric Data" << "\n";

  typename std::vector <internal::PointValueHistory::PointGeometryData <dim> >::iterator node = point_geometry_data.begin ();
  if (node == point_geometry_data.end ())
    {
      out << "No nodes stored currently\n";
    }
  else
    {
      if (!cleared)
	{
	  out << "# DoF_index : Location (for each component)\n";
	  for (; node != point_geometry_data.end (); node++)
	    {
	      for (unsigned int component = 0; component < dof_handler->get_fe ().n_components (); component++)
		{
		  out << node->solution_indices[component] << " : " << node->locations [component] << "\n";
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

  std::cout << dataset_key.size () << " data keys are stored \n";


  if (independent_values.size () != 0)
    {
      out << "Independent value(s): " << independent_values.size () << " : " << independent_values[0].size () << "\n";
    }
  else
    {
      out << "No independent values stored\n";
    }

  std::map <std::string, std::vector <std::vector <double> > >::iterator
    data_store_begin = data_store.begin ();
  for (; data_store_begin != data_store.end (); data_store_begin++)
    {
      if (data_store_begin->second.size () != 0)
        {
          out << data_store_begin->first << ": " << data_store_begin->second.size () << " : " << (data_store_begin->second)[0].size () << "\n";
        }
      else
        {
          out << data_store_begin->first << ": " << data_store_begin->second.size () << " : " << "No points added" << "\n";
        }
    }
  out << "\n";
  out << "***end of status output***\n\n";
}



template <int dim>
bool PointValueHistory<dim>
::deep_check (bool strict)
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
      if (n_dofs != 0)
        {
          for (; data_store_begin != data_store.end (); data_store_begin++)
            {
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

  if (n_dofs != 0)
    {
      std::map <std::string, std::vector <std::vector <double> > >::iterator
	data_store_begin = data_store.begin ();
      for (; data_store_begin != data_store.end (); data_store_begin++)
        {
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


// explicit instantiations
template class PointValueHistory<deal_II_dimension>;


DEAL_II_NAMESPACE_CLOSE


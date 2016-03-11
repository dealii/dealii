// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2016 by the deal.II authors
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


#ifndef dealii__point_value_history_h
#define dealii__point_value_history_h

#include <deal.II/base/point.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/numerics/data_postprocessor.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace PointValueHistory
  {
    /**
     * A class that stores the data needed to reference the support points
     * closest to one requested point.
     */
    template <int dim>
    class PointGeometryData
    {
    public:
      PointGeometryData(const Point <dim> &new_requested_location, const std::vector <Point <dim> > &new_locations,
                        const std::vector <types::global_dof_index> &new_sol_indices);
      Point <dim> requested_location;
      std::vector <Point <dim> > support_point_locations;
      std::vector <types::global_dof_index> solution_indices;
    };
  }
}



/**
 * PointValueHistory tackles the overhead of plotting time (or any other
 * iterative process) graphs of solution values at specific points on the
 * mesh. The user specifies the points which the solution should be monitored
 * at ahead of time, as well as giving each solution vector that they want to
 * record a mnemonic name. Then, for each step the user calls one of the three
 * available "evaluate field" methods to store the data from each time step,
 * and the class extracts data for the requested points to store it. Finally,
 * once the computation is finished, the user can request output files to be
 * generated; these files are in Gnuplot format but are basically just regular
 * text and can easily be imported into other programs well, for example into
 * spreadsheets.
 *
 * The user can store extra variables which do not relate to mesh location
 * specifying n_independent_variables. The class then expects a std::vector of
 * size n_independent_variables to be added during each step using the method
 * @p push_back_independent. This may be used for example for recording
 * external input, logging solver performance data such as time taken to solve
 * the step and solver steps before convergence, saving norms calculated, or
 * simply saving the time, number of time step, or number of nonlinear
 * iteration along with the data evaluated from the mesh.
 *
 * The three "evaluate field" methods each have different strengths and
 * weaknesses making each suitable for different contexts:
 * <ol>
 * <li>Firstly, the @p evaluate_field version that does not take a @p
 * DataPostprocessor object selects the nearest support point (see
 * @ref GlossSupport "this entry in the glossary"
 * ) to a given point to extract data from. This makes the code that needs to
 * be run at each time step very short, since looping over the mesh to extract
 * the needed dof_index can be done just once at the start. However, this
 * method is not suitable for FiniteElement objects that do not assign dofs to
 * actual mesh locations (i.e. FEs without
 * @ref GlossSupport "support points"
 * ) or if adaptive mesh refinement is used. The reason for the latter
 * restriction is that the location of the closest support point to a given
 * point may change upon mesh refinement. The class will throw an exception if
 * any change to the triangulation is made (Although the nearest support point
 * could be re- computed upon mesh refinement, the location of the support
 * point will most likely change slightly, making the interpretation of the
 * data difficult, hence this is not implemented currently.)
 *
 * <li> Secondly, @p evaluate_field_at_requested_location calls @p
 * VectorTools::point_value to compute values at the specific point requested.
 * This method is valid for any FE that is supported by @p
 * VectorTools::point_value. Specifically, this method can be called by codes
 * using adaptive mesh refinement.
 *
 * <li>Finally, the class offers a function @p evaluate_field that takes a @p
 * DataPostprocessor object. This method allows the deal.II data postprocessor
 * to be used to compute new quantities from the solution on the fly. The
 * values are located at the nearest quadrature point to the requested point.
 * If the mesh is refined between calls, this point will change, so care must
 * be taken when using this method in code using adaptive refinement, but as
 * the output will be meaningful (in the sense that the quadrature point
 * selected is guaranteed to remain in the same vicinity, the class does not
 * prevent the use of this method in adaptive codes. The class provides
 * warnings in the output files if the mesh has changed. Note that one can
 * reduce the error this procedure introduces by providing a quadrature
 * formula that has more points, at the expense of performing more work since
 * then the closest quadrature points is nearer to the point at which the
 * evaluation is really supposed to happen. (As a sidenote: Why not do the
 * evaluation at the requested point right away? The reason for this is that
 * it would require setting up a new quadrature point object on each cell that
 * has only a single point corresponding to the reference coordinates of the
 * point you really want; then initializing a FEValues object with it; then
 * evaluating the solution at this point; then handing the result to the
 * DataPostprocessor object. This sequence of things is expensive -- which is
 * the reason why VectorTools::point_value() is expensive. Using the same
 * quadrature formula on each cell on which we want to evaluate the solution
 * and only having to initialize a FEValue object once is a much cheaper
 * alternative, albeit of course at the expense of getting only an approximate
 * result.)
 * </ol>
 *
 * When recording a new mnemonic name, the user must supply a component_mask
 * (see
 * @ref GlossComponentMask "this glossary entry"
 * ) to indicate the
 * @ref GlossComponent "(vector) components"
 * to be extracted from the given input. If the user simply wants to extract
 * all the components, the mask need not be explicitly supplied to the @p
 * add_field_name method and the default value of the parameter is sufficient.
 * If the @p evaluate_field with a @p DataPostprocessor object is used, the
 * component_mask is interpreted as the mask of the @p DataPostprocessor
 * return vector. The size of this mask can be different to that of the FE
 * space, but must be provided when the @p add_field_name method is called.
 * One variant of the @p add_field_name method allows an unsigned int input to
 * construct a suitable mask, if all values from the @p DataPostprocessor are
 * desired.
 *
 * The class automatically generates names for the data stored based on the
 * mnemonics supplied. The methods @p add_component_names and @p
 * add_independent_names allow the user to provide lists of names to use
 * instead if desired.
 *
 * Following is a little code snippet that shows a common usage of this class:
 *
 * @code
 * #include <deal.II/numerics/point_value_history.h>
 * //....
 *
 * //... code to setup Triangulation, perform any refinement necessary
 * // and setup DoFHandler, sizing solution Vectors etc
 *
 * // call the constructor
 * unsigned int n_inputs = 1; // just one independent value, which happens to be an input
 * PointValueHistory<dim> node_monitor(dof_handler, n_inputs);
 *
 * // setup fields and points required
 * node_monitor.add_field_name("Solution");
 * std::vector <Point <dim> > point_vector(2);
 * point_vector[0] = Point <dim>(0, 0);
 * point_vector[1] = Point <dim>(0.25, 0);
 * node_monitor.add_points(point_vector); // multiple points at once
 * node_monitor.add_point(Point<dim>(1, 0.2)); // add a single point
 * node_monitor.close(); // close the class once the setup is complete
 * node_monitor.status(std::cout); // print out status to check if desired
 *
 * // ... more code ...
 *
 * // ... in an iterative loop ...
 * // double time, vector <double> with size 1 input_value,
 * // and Vector <double> solution calculated in the loop
 * node_monitor.start_new_dataset(time);
 * node_monitor.push_back_independent(input_value);
 * node_monitor.evaluate_field("Solution", solution);
 *
 * // ... end of iterative loop ...
 *
 * node_monitor.write_gnuplot("node"); // write out data files
 *
 * @endcode
 */
template <int dim>
class PointValueHistory
{
public:
  /**
   * Provide a stripped down instance of the class which does not support
   * adding points or mesh data.  This may be used for example for recording
   * external input or logging solver performance data.
   */
  PointValueHistory (const unsigned int n_independent_variables = 0);

  /**
   * Constructor linking the class to a specific @p DoFHandler. This class
   * reads specific data from the @p DoFHandler and stores it internally for
   * quick access (in particular dof indices of closest neighbors to requested
   * points) the class is fairly intolerant to changes to the @p DoFHandler if
   * data at support points is required. Mesh refinement and @p DoFRenumbering
   * methods should be performed before the @p add_points method is called and
   * adaptive grid refinement is only supported by some methods.
   *
   * The user can store extra variables which do not relate to mesh location
   * by specifying the number required using n_independent_variables and
   * making calls to @p push_back_independent as needed.  This may be used for
   * example for recording external input or logging solver performance data.
   */
  PointValueHistory (const DoFHandler<dim> &dof_handler,
                     const unsigned int n_independent_variables = 0);

  /**
   * Copy constructor. This constructor can be safely called with a @p
   * PointValueHistory object that contains data, but this could be expensive
   * and should be avoided.
   */
  PointValueHistory (const PointValueHistory &point_value_history);

  /**
   * Assignment operator. This assignment operator can be safely called once
   * the class is closed and data added, but this is provided primarily to
   * allow a @p PointValueHistory object declared in a class to be
   * reinitialized later in the class. Using the assignment operator when the
   * object contains data could be expensive.
   */
  PointValueHistory &operator=(const PointValueHistory &point_value_history);

  /**
   * Deconstructor.
   */
  ~PointValueHistory ();

  /**
   * Add a single point to the class. The support points (one per component)
   * in the mesh that are closest to that point are found and their details
   * stored for use when @p evaluate_field is called. If more than one point
   * is required rather use the @p add_points method since this minimizes
   * iterations over the mesh.
   */
  void add_point(const Point <dim> &location);

  /**
   * Add multiple points to the class. The support points (one per component)
   * in the mesh that are closest to that point is found and their details
   * stored for use when @p evaluate_field is called. If more than one point
   * is required, rather call this method as it is more efficient than the
   * add_point method since it minimizes iterations over the mesh. The points
   * are added to the internal database in the order they appear in the list
   * and there is always a one to one correspondence between the requested
   * point and the added point, even if a point is requested multiple times.
   */
  void add_points (const std::vector <Point <dim> > &locations);



  /**
   * Put another mnemonic string (and hence @p VectorType) into the class.
   * This method adds storage space for variables equal to the number of true
   * values in component_mask. This also adds extra entries for points that
   * are already in the class, so @p add_field_name and @p add_points can be
   * called in any order.
   */
  void add_field_name(const std::string &vector_name,
                      const ComponentMask &component_mask = ComponentMask());

  /**
   * Put another mnemonic string (and hence @p VectorType) into the class.
   * This method adds storage space for n_components variables. This also adds
   * extra entries for points that are already in the class, so @p
   * add_field_name and @p add_points can be called in any order. This method
   * generates a std::vector 0, ..., n_components-1 and calls the previous
   * function.
   */
  void add_field_name(const std::string &vector_name,
                      const unsigned int n_components);

  /**
   * Provide optional names for each component of a field. These names will be
   * used instead of names generated from the field name, if supplied.
   */
  void add_component_names(const std::string &vector_name,
                           const std::vector <std::string> &component_names);

  /**
   * Provide optional names for the independent values. These names will be
   * used instead of "Indep_...", if supplied.
   */
  void add_independent_names(const std::vector <std::string> &independent_names);



  /**
   * Extract values at the stored points from the VectorType supplied and add
   * them to the new dataset in vector_name. The component mask supplied when
   * the field was added is used to select components to extract. If a @p
   * DoFHandler is used, one (and only one) evaluate_field method must be
   * called for each dataset (time step, iteration, etc) for each vector_name,
   * otherwise a @p ExcDataLostSync error can occur.
   */
  template <class VectorType>
  void evaluate_field(const std::string &name,
                      const VectorType  &solution);


  /**
   * Compute values using a @p DataPostprocessor object with the @p VectorType
   * supplied and add them to the new dataset in vector_name. The
   * component_mask supplied when the field was added is used to select
   * components to extract from the @p DataPostprocessor return vector. This
   * method takes a vector of field names to process and is preferred if many
   * fields use the same @p DataPostprocessor object as each cell is only
   * located once. The quadrature object supplied is used for all components
   * of a vector field. Although this method will not throw an exception if
   * the mesh has changed. (No internal data structures are invalidated as the
   * quadrature points are repicked each time the function is called.)
   * Nevertheless the user must be aware that if the mesh changes the point
   * selected will also vary slightly, making interpretation of the data more
   * difficult. If a @p DoFHandler is used, one (and only one) evaluate_field
   * method must be called for each dataset (time step, iteration, etc) for
   * each vector_name, otherwise a @p ExcDataLostSync error can occur.
   */
  template <class VectorType>
  void evaluate_field(const std::vector <std::string> &names,
                      const VectorType                &solution,
                      const DataPostprocessor<dim>    &data_postprocessor,
                      const Quadrature<dim>           &quadrature);

  /**
   * Construct a std::vector <std::string> containing only vector_name and
   * call the above function. The above function is more efficient if multiple
   * fields use the same @p DataPostprocessor object.
   */
  template <class VectorType>
  void evaluate_field(const std::string            &name,
                      const VectorType             &solution,
                      const DataPostprocessor<dim> &data_postprocessor,
                      const Quadrature<dim>        &quadrature);


  /**
   * Extract values at the points actually requested from the VectorType
   * supplied and add them to the new dataset in vector_name. Unlike the other
   * evaluate_field methods this method does not care if the dof_handler has
   * been modified because it uses calls to @p VectorTools::point_value to
   * extract there data. Therefore, if only this method is used, the class is
   * fully compatible with adaptive refinement. The component_mask supplied
   * when the field was added is used to select components to extract. If a @p
   * DoFHandler is used, one (and only one) evaluate_field method must be
   * called for each dataset (time step, iteration, etc) for each vector_name,
   * otherwise a @p ExcDataLostSync error can occur.
   */
  template <class VectorType>
  void evaluate_field_at_requested_location(const std::string &name,
                                            const VectorType  &solution);


  /**
   * Add the key for the current dataset to the dataset. Although calling this
   * method first is sensible, the order in which this method, @p
   * evaluate_field and @p push_back_independent is not important. It is
   * however important that all the data for a give dataset is added to each
   * dataset and that it is added before a new data set is started. This
   * prevents a @p ExcDataLostSync.
   */
  void start_new_dataset (const double key);

  /**
   * If independent values have been set up, this method stores these values.
   * This should only be called once per dataset, and if independent values
   * are used it must be called for every dataset. A @p ExcDataLostSync
   * exception can be thrown if this method is not called.
   */
  void push_back_independent (const std::vector <double> &independent_values);


  /**
   * Write out a series of .gpl files named base_name + "-00.gpl", base_name +
   * "-01.gpl" etc. The data file gives info about where the support points
   * selected and interpreting the data. If @p n_indep != 0 an additional file
   * base_name + "_indep.gpl" containing key and independent data. The file
   * name agrees with the order the points were added to the class. The names
   * of the data columns can be supplied using the functions @p
   * add_component_names and @p add_independent_names. The support point
   * information is only meaningful if the dof_handler has not been changed.
   * Therefore, if adaptive mesh refinement has been used the support point
   * data should not be used. The optional parameter postprocessor_locations
   * is used to add the postprocessor locations to the output files. If this
   * is desired, the data should be obtained from a call to
   * get_postprocessor_locations while the dof_handler is usable. The default
   * parameter is an empty vector of strings, and will suppress postprocessor
   * locations output.
   */
  void write_gnuplot (const std::string &base_name,
                      const std::vector <Point <dim> > postprocessor_locations = std::vector <Point <dim> > ());


  /**
   * Return a @p Vector with the indices of selected points flagged with a 1.
   * This method is mainly for testing and verifying that the class is working
   * correctly. By passing this vector to a DataOut object, the user can
   * verify that the positions returned by @p get_points agree with the
   * positions that @p DataOut interprets from the @p Vector returned. The
   * code snippet below demonstrates how this could be done:
   * @code
   * // Make a DataOut object and attach the dof_handler
   * DataOut<dim> data_out;
   * data_out.attach_dof_handler(dof_handler);
   *
   * // Call the mark_locations method to get the vector with indices flagged
   * Vector<double> support_point_locations = node_monitor.mark_locations();
   *
   * // Add the vector to the data_out object and write out a file in the usual way
   * data_out.add_data_vector(support_point_locations, "Monitor_Locations");
   * data_out.build_patches(2);
   * std::ofstream output("locations.gpl");
   * data_out.write_gnuplot(output);
   * @endcode
   */
  Vector<double> mark_support_locations();

  /**
   * Stores the actual location of each support point selected by the @p
   * add_point(s) method.  This can be used to compare with the point
   * requested, for example by using the @p Point<dim>::distance function. For
   * convenience, location is resized to the correct number of points by the
   * method.
   */
  void get_support_locations (std::vector <std::vector<Point <dim> > > &locations);

  /**
   * @deprecated
   *
   * This function only exists for backward compatibility as this is the
   * interface provided by previous versions of the library. The function
   * get_support_locations replaces it and reflects the fact that the points
   * returned are actually the support points.
   */
  void get_points (std::vector <std::vector<Point <dim> > > &locations);

  /**
   * Stores the actual location of the points used by the data_postprocessor.
   * This can be used to compare with the points requested, for example by
   * using the @p Point<dim>::distance function. Unlike the support_locations,
   * these locations are computed every time the evaluate_field method is
   * called with a postprocessor. This method uses the same algorithm so can
   * will find the same points. For convenience, location is resized to the
   * correct number of points by the method.
   */
  void get_postprocessor_locations (const Quadrature<dim> &quadrature,
                                    std::vector<Point <dim> > &locations);

  /**
   * Once datasets have been added to the class, requests to add additional
   * points will make the data interpretation unclear. The boolean @p closed
   * defines a state of the class and ensures this does not happen. Additional
   * points or vectors can only be added while the class is not closed, and
   * the class must be closed before datasets can be added or written to file.
   * @p PointValueHistory::get_points and @p PointValueHistory::status do not
   * require the class to be closed. If a method that requires a class to be
   * open or close is called while in the wrong state a @p ExcInvalidState
   * exception is thrown.
   */
  void close();


  /**
   * Delete the lock this object has to the @p DoFHandler used the last time
   * the class was created.  This method should not normally need to be
   * called, but can be useful to ensure that the @p DoFHandler is released
   * before it goes out of scope if the @p PointValueHistory class might live
   * longer than it. Once this method has been called, the majority of methods
   * will throw a @p ExcInvalidState exception, so if used this method should
   * be the last call to the class.
   */
  void clear();

  /**
   * Print useful debugging information about the class, include details about
   * which support points were selected for each point and sizes of the data
   * stored.
   */
  void status(std::ostream &out);


  /**
   * Check the internal data sizes to test for a loss of data sync. This is
   * often used in @p Assert statements with the @p ExcDataLostSync exception.
   * If @p strict is @p false this method returns @p true if all sizes are
   * within 1 of each other (needed to allow data to be added), with @p strict
   * = @p true they must be exactly equal.
   */

  bool deep_check (const bool strict);

  /**
   * Exception
   */
  DeclExceptionMsg(ExcNoIndependent,
                   "A call has been made to push_back_independent() when "
                   "no independent values were requested.");

  /**
   * Exception
   */
  DeclExceptionMsg(ExcDataLostSync,
                   "This error is thrown to indicate that the data sets appear to be out of "
                   "sync. The class requires that the number of dataset keys is the same as "
                   "the number of independent values sets and mesh linked value sets. The "
                   "number of each of these is allowed to differ by one to allow new values "
                   "to be added with out restricting the order the user choses to do so. "
                   "Special cases of no FHandler and no independent values should not "
                   "trigger this error.");


  /**
   * Exception
   */
  DeclExceptionMsg(ExcDoFHandlerRequired,
                   "A method which requires access to a @p DoFHandler to be meaningful has "
                   "been called when have_dof_handler is false (most likely due to default "
                   "constructor being called). Only independent variables may be logged with "
                   "no DoFHandler.");

  /**
   * Exception
   */
  DeclExceptionMsg(ExcDoFHandlerChanged,
                   "The triangulation has been refined or coarsened in some way. This "
                   "suggests that the internal DoF indices stored by the current "
                   "object are no longer meaningful.");

private:
  /**
   * Stores keys, values on the abscissa. This will often be time, but
   * possibly time step, iteration etc.
   */
  std::vector <double> dataset_key;

  /**
   * Values that do not depend on grid location.
   */
  std::vector <std::vector <double> > independent_values;

  /**
   * Saves a vector listing component names associated with a
   * independent_values. This will be an empty vector if the user does not
   * supplies names.
   */
  std::vector<std::string> indep_names;

  /**
   * Saves data for each mnemonic entry. data_store: mnemonic ->
   * [point_0_components point_1_components ... point_n-1_components][key]
   * This format facilitates scalar mnemonics in a vector space, because
   * scalar mnemonics will only have one component per point. Vector
   * components are strictly FE.n_components () long.
   */
  std::map <std::string, std::vector <std::vector <double> > > data_store;

  /**
   * Saves a component mask for each mnemonic.
   */
  std::map <std::string, ComponentMask> component_mask;


  /**
   * Saves a vector listing component names associated with a mnemonic. This
   * will be an empty vector if the user does not supplies names.
   */
  std::map <std::string, std::vector<std::string> > component_names_map;

  /**
   * Saves the location and other mesh info about support points.
   */
  std::vector <internal::PointValueHistory::PointGeometryData <dim> >
  point_geometry_data;


  /**
   * Used to enforce @p closed state for some methods.
   */
  bool closed;

  /**
   * Used to enforce @p !cleared state for some methods.
   */
  bool cleared;


  /**
   * A smart pointer to the dof_handler supplied to the constructor. This can
   * be released by calling @p clear().
   */
  SmartPointer<const DoFHandler<dim>,PointValueHistory<dim> > dof_handler;


  /**
   * Variable to check if the triangulation has changed. If it has changed,
   * certain data is out of date (especially the
   * PointGeometryData::solution_indices.
   */
  bool triangulation_changed;

  /**
   * A boolean to record whether the class was initialized with a DoFHandler
   * or not.
   */
  bool have_dof_handler;

  /**
   * Used to detect signals from the Triangulation.
   */
  boost::signals2::connection tria_listener;

  /**
   * Stores the number of independent variables requested.
   */
  unsigned int n_indep;


  /**
   * A function that will be triggered through signals whenever the
   * triangulation is modified.
   *
   * It is currently used to check if the triangulation has changed,
   * invalidating precomputed values.
   */
  void tria_change_listener ();
};


DEAL_II_NAMESPACE_CLOSE
#endif /* dealii__point_value_history_h */

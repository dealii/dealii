//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2009, 2012 by Michael Rapson and the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __dealii__point_value_history_h
#define __dealii__point_value_history_h

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/utilities.h>
#include <numerics/data_postprocessor.h>
#include <grid/grid_tools.h>

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
    template <int dim> class PointGeometryData;
  }
}



namespace internal
{
  namespace PointValueHistory
  {
                                     /**
                                      *  A class that stores the data needed
                                      *  to reference the support points
                                      *  closest to one requested point.
                                      */
    template <int dim>
    class PointGeometryData
    {
      public:
        PointGeometryData(const Point <dim> &new_requested_location, const std::vector <Point <dim> > &new_locations,
                          const std::vector <int> &new_sol_indices);
        Point <dim> requested_location;
        std::vector <Point <dim> > support_point_locations;
        std::vector <int> solution_indices;
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
 * and the class extracts data for the requested
 * points to store it. Finally, once the computation is finished, the user can
 * request Gnuplot output files to be generated.
 *
 * The user can store extra variables which do not relate to mesh location by
 * the parameter by specifying n_independent_variables. The class then expects
 * a std::vector of size n_independent_variables to be added during each step
 * using the method @p push_back_independent. This may be used for example for
 * recording external input, logging solver performance data such as time
 * taken to solve the step and solver steps before convergence, or saving
 * norms calculated.
 *
 * The three "evaluate field" methods each have different strengths and
 * weaknesses making each suitable for different contexts:
 * <ol>
 * <li>Firstly, the @p evaluate_field version that does not take a @p DataPostprocessor object
 * selects the nearest support points to a given point to
 * extract data from. This makes the code that needs to be run at each time step very short,
 * since looping over the mesh to extract the needed dof_index can be done
 * just once at the start. However, this method is not suitable for
 * FiniteElement objects that do not assign dofs to actual mesh locations
 * (i.e. FEs without @ref GlossSupport "support points") or if adaptive mesh refinement is used.
 * The class will throw an exception if the finite element does not have support points
 * or any change in the @p dof_handler is identified. (The test is simply
 * whether the number of dofs changes.)
 *
 * <li> Secondly, @p evaluate_field_at_requested_location calls @p
 * VectorTools::point_value to compute values at the specific point requested.
 * This method is valid for any FE that is supported by @p VectorTools::point_value.
 * Specifically, this method can be called by codes using adaptive mesh refinement.
 *
 * <li>Finally, the class offers a function @p evaluate_field that takes a
 * @p DataPostprocessor object. This method allows the deal.II
 * data postprocessor to be used to compute new quantities from the
 * solution on the fly. The values are located at the nearest quadrature
 * point to the requested point. If the mesh is refined between calls, this
 * point will change, so care must be taken when using this method in
 * code using adaptive refinement, but as the output will be meaningful (in
 * the sense that the quadrature point selected is guaranteed to remain in the
 * same vicinity, the class does not prevent the use of this method in adaptive
 * codes.
 * </ol>
 *
 * When recording a new mnemonic name, the user must chose whether it should be
 * a strictly scalar field (@p add_scalar_field_name) or whether it should
 * depend on the @p FE (@p add_field_name). Calling @p add_field_name with a
 * @p dof_handler using a vector @p FE will result in a vector field. The vector
 * field will be of dimention @p FE.n_components(). Therefore, in this class the
 * distinction between vector and scalar is as related to the @p FE and not
 * whether a quantity is a logical vector.
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
                                      * Provide a stripped down instance of
                                      * the class which does not support
                                      * adding points or mesh data.  This may
                                      * be used for example for recording
                                      * external input or logging solver
                                      * performance data.
                                      */
    PointValueHistory (const unsigned int n_independent_variables = 0);

                                     /**
                                      * Constructor linking the class to a
                                      * specific @p DoFHandler. Since this
                                      * class reads specific data from the @p
                                      * DoFHandler and stores it internally
                                      * for quick access (in particular dof
                                      * indices of closest neighbors to
                                      * requested points) the class is fairly
                                      * intolerant to changes to the @p
                                      * DoFHandler.  Mesh refinement and @p
                                      * DoFRenumbering methods should be
                                      * performed before the @p add_points
                                      * method is called and adaptive grid
                                      * refinement is not supported (unless it
                                      * is completed before points are added).
                                      * Changes to the @p DoFHandler are
                                      * tested for by looking for a change in
                                      * the number of dofs, which may not
                                      * always occur so the user must be aware
                                      * of these limitations.
                                      *
                                      * The user can store extra variables
                                      * which do not relate to mesh location
                                      * by specifying the number required
                                      * using n_independent_variables and
                                      * making calls to @p
                                      * push_back_independent as needed.  This
                                      * may be used for example for recording
                                      * external input or logging solver
                                      * performance data.
                                      */
    PointValueHistory (const DoFHandler<dim> & dof_handler,
                       const unsigned int n_independent_variables = 0);

                                     /**
                                      * Copy constructor explicitly provided,
                                      * although default copy operator would
                                      * be sufficient.  This constructor can
                                      * be safely called with a @p
                                      * PointValueHistory object that contains
                                      * data, but this could be expensive and
                                      * should be avoided.
                                      */
    PointValueHistory (const PointValueHistory & point_value_history);

                                     /**
                                      * Assignment operator explicitly
                                      * provided, although default assignment
                                      * operator would be sufficient.  This
                                      * assignment operator can be safely
                                      * called once the class is closed and
                                      * data added, but this is provided
                                      * primarily to allow a @p
                                      * PointValueHistory object declared in a
                                      * class to be reinitialized later in the
                                      * class. Using the assignment operator
                                      * when the object contains data could be
                                      * expensive.
                                      */
    PointValueHistory & operator=(const PointValueHistory & point_value_history);

                                     /**
                                        Deconstructor.
                                     */
    ~PointValueHistory ();

                                     /**
                                      * Add a single point to the class. The
                                      * support points (one per component) in
                                      * the mesh that are closest to that
                                      * point is found and its details stored
                                      * for use when @p evaluate_field is
                                      * called. If more than one point is
                                      * required rather used the @p add_points
                                      * method since this minimizes iterations
                                      * over the mesh.
                                      */
    void add_point(const Point <dim> & location);

                                     /**
                                      * Add multiple points to the class. The
                                      * support points (one per component) in
                                      * the mesh that are closest to that
                                      * point is found and its details stored
                                      * for use when @p evaluate_field is
                                      * called. If more than one point is
                                      * required, rather call this method as it
                                      * is more efficient than the add_point
                                      * method since it minimizes iterations
                                      * over the mesh. The points are added to
                                      * the internal database in the order
                                      * they appear in the list and there is
                                      * always a one to one correspondence
                                      * between the requested point and the
                                      * added point, even if a point is
                                      * requested multiple times.
                                      */
    void add_points (const std::vector <Point <dim> > & locations);



                                     /**
                                      * Put another mnemonic string (and hence
                                      * @p VECTOR) into the class. This method
                                      * adds storage space for FE.n_components()
                                      * variables and should be used with methods
                                      * that extract all FE components from a vector.
                                      * This also adds extra entries for points
                                      * that are already in the class, so
                                      * @p add_field_name and @p add_points can
                                      * be called in any order.
                                      */
    void add_field_name(const std::string &vector_name);

                                     /**
                                      * Put another mnemonic string (and hence
                                      * @p VECTOR) into the class. This method
                                      * only adds storage space for a single
                                      * variable whether the FE is scalar or not
                                      * so it can be used with
                                      * @p DataPostprocessor objects that compute
                                      * a scalar value(s) from a vector field.
                                      * This also adds extra entries for points
                                      * that are already in the class, so
                                      * @p add_field_name and @p add_points can
                                      * be called in any order.
                                      */
    void add_scalar_field_name(const std::string &vector_name);

                                     /**
                                      * Extract values at the stored points
                                      * from the VECTOR supplied and add them
                                      * to the new dataset in vector_name.
                                      * The scalar_component input is only
                                      * used if the FE has multiple components
                                      * and the field to be evaluated is scalar.
                                      * In this case, the scalar_component input
                                      * is used to select the component used. If
                                      * a @p DoFHandler is used, one (and only
                                      * one) evaluate_field method
                                      * must be called for each dataset (time
                                      * step, iteration, etc) for each
                                      * vector_name, otherwise a @p
                                      * ExcDataLostSync error can occur.
                                      */
    template <class VECTOR>
    void evaluate_field(const std::string &vector_name, const VECTOR & solution, const unsigned int scalar_component = 0);


                                     /**
                                      * Compute values using a DataPostprocessor object
                                      * with the VECTOR supplied and add them
                                      * to the new dataset in vector_name.
                                      * The DataPostprocessor component_interpretation
                                      * is used to determine select scalar or vector data.
                                      * If the component_is_part_of_vector flag is used, then
                                      * the indicated vector must have the same number of
                                      * components as the FE associated with the dof_handler.
                                      * The DataPostprocessor object is queried for field names.
                                      * In the case of vector fields, the name of the first
                                      * component is used as the overall name of the vector.
                                      * The quadrature object supplied is used for all
                                      * components of a vector field.
                                      * If a @p DoFHandler is used, one (and only
                                      * one) evaluate_field method
                                      * must be called for each dataset (time
                                      * step, iteration, etc) for each
                                      * vector_name, otherwise a @p
                                      * ExcDataLostSync error can occur.
                                      */
    template <class VECTOR>
    void evaluate_field(const VECTOR & solution, const DataPostprocessor<dim> & data_postprocessor, const Quadrature<dim> & quadrature);


                                     /**
                                      * Extract values at the points actually
                                      * requested from the VECTOR supplied and
                                      * add them to the new dataset in vector_name.
                                      * Unlike the other evaluate_field methods
                                      * this method does not care if the dof_handler
                                      * has been modified because it uses calls to
                                      * @p VectorTools::point_value to extract there
                                      * data. Therefore, if only this method is used,
                                      * the class can be used with adaptive refinement.
                                      * Specifically, the vector form of
                                      * @p VectorTools::point_value that assumes a
                                      * Q1 mapping is used.
                                      * The scalar_component input is only
                                      * used if the FE has multiple components
                                      * and the field to be evaluated is scalar.
                                      * In this case, the scalar_component input
                                      * is used to select the component used. If
                                      * a @p DoFHandler is used, one (and only
                                      * one) evaluate_field method
                                      * must be called for each dataset (time
                                      * step, iteration, etc) for each
                                      * vector_name, otherwise a @p
                                      * ExcDataLostSync error can occur.
                                      */
    template <class VECTOR>
    void evaluate_field_at_requested_location(const std::string &vector_name, const VECTOR & solution, const unsigned int scalar_component = 0);

                                     /**
                                      * Add the key for the current dataset to
                                      * the dataset. Although calling this
                                      * method first is sensible, the order in
                                      * which this method, @p evaluate_field
                                      * and @p push_back_independent is not
                                      * important. It is however important
                                      * that all the data for a give dataset
                                      * is added to each dataset and that it
                                      * is added before a new data set is
                                      * started. This prevents a @p
                                      * ExcDataLostSync.
                                      */
    void start_new_dataset(double key);

                                     /**
                                      * If independent values have been set
                                      * up, this method stores these
                                      * values. This should only be called
                                      * once per dataset, and if independent
                                      * values are used it must be called for
                                      * every dataset. A @p ExcDataLostSync
                                      * exception can be thrown if this method
                                      * is not called.
                                      */
    void push_back_independent(const std::vector <double> &independent_values);


                                     /**
                                      * Write out a series of .gpl files named
                                      * base_name + "-00.gpl", base_name +
                                      * "-01.gpl" etc. The data file gives
                                      * info about where the support points
                                      * selected and interpreting the data.
                                      * If @p n_indep != 0 an additional file
                                      * base_name + "_indep.gpl" containing
                                      * key and independent data. The file
                                      * name agrees with the order the points
                                      * were added to the class.
                                      * The support point information is only
                                      * meaningful if the dof_handler has not
                                      * been changed. Therefore, if adaptive
                                      * mesh refinement has been used the
                                      * support point data should not be used.
                                      * The optional parameter
                                      * postprocessor_locations is used to
                                      * add the postprocessor locations to the
                                      * output files. If this is desired, the
                                      * data should be obtained from a call to
                                      * get_postprocessor_locations while the
                                      * dof_handler is usable. The default
                                      * parameter is an empty vector of strings,
                                      * and will suppress postprocessor
                                      * locations output.
                                      */
    void write_gnuplot(const std::string &base_name, const std::vector <Point <dim> > postprocessor_locations = std::vector <Point <dim> > ());


                                     /**
                                      * Return a @p Vector with the indices of
                                      * selected points flagged with a 1. This
                                      * method is mainly for testing and
                                      * verifying that the class is working
                                      * correctly. By passing this vector to a
                                      * DataOut object, the user can verify
                                      * that the positions returned by @p
                                      * PointValueHistory< dim >::get_points
                                      * agree with the positions that @p
                                      * DataOut interprets from the @p Vector
                                      * returned. The code snippet below
                                      * demonstrates how this could be done:
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
                                      * Depreciated:
                                      *
                                      * This function only exists for backward
                                      * compatibility as this is the interface
                                      * provided by previous versions of the library.
                                      * The function mark_support_locations replaces
                                      * it and reflects the fact that the locations
                                      * marked are actually the support points.
                                      */
    Vector<double> mark_locations();


                                     /**
                                      * Stores the actual location of each
                                      * support point selected by the @p
                                      * add_point(s) method.  This can be used
                                      * to compare with the point requested,
                                      * for example by using the @p
                                      * Point<dim>::distance function. For
                                      * convenience, location is resized to
                                      * the correct number of points by the
                                      * method.
                                      */
    void get_support_locations (std::vector <std::vector<Point <dim> > > & locations);

                                     /**
                                      * Depreciated:
                                      *
                                      * This function only exists for backward
                                      * compatibility as this is the interface
                                      * provided by previous versions of the library.
                                      * The function get_support_locations replaces
                                      * it and reflects the fact that the points
                                      * returned are actually the support points.
                                      */
    void get_points (std::vector <std::vector<Point <dim> > > & locations);

                                     /**
                                      * Stores the actual location of the
                                      * points used by the data_postprocessor.
                                      * This can be used to compare with the
                                      * points requested, for example by using
                                      * the @p Point<dim>::distance function.
                                      * Unlike the support_locations, these
                                      * locations are computed every time the
                                      * evaluate_field method is called with a
                                      * postprocessor. This method uses the
                                      * same algorithm so can will find the
                                      * same points.
                                      * For convenience, location is resized
                                      * to the correct number of points by the
                                      * method.
                                      */
    void get_postprocessor_locations (std::vector<Point <dim> > & locations, const Quadrature<dim> & quadrature);

                                     /**
                                      * Once datasets have been added to the
                                      * class, requests to add additional
                                      * points will make the data
                                      * interpretation unclear. The boolean @p
                                      * closed defines a state of the class
                                      * and ensures this does not
                                      * happen. Additional points or vectors
                                      * can only be added while the class is
                                      * not closed, and the class must be
                                      * closed before datasets can be added or
                                      * written to file. @p
                                      * PointValueHistory::get_points and @p
                                      * PointValueHistory::status do not
                                      * require the class to be closed. If a
                                      * method that requires a class to be
                                      * open or close is called while in the
                                      * wrong state a @p ExcInvalidState
                                      * exception is thrown.
                                      */
    void close();


                                     /**
                                      * Delete the lock this object has to the
                                      * @p DoFHandler used the last time the
                                      * class was created.  This method should
                                      * not normally need to be called, but
                                      * can be useful to ensure that the @p
                                      * DoFHandler is released before it goes
                                      * out of scope if the @p
                                      * PointValueHistory class might live
                                      * longer than it. Once this method has
                                      * been called, the majority of methods
                                      * will throw a @p ExcInvalidState
                                      * exception, so if used this method
                                      * should be the last call to the class.
                                      */
    void clear();

                                     /**
                                      * Print useful debugging information
                                      * about the class, include details about
                                      * which support points were selected for
                                      * each point and sizes of the data
                                      * stored.
                                      */
    void status(std::ostream &out);


                                     /**
                                      * Check the internal data sizes to test
                                      * for a loss of data sync. This is often
                                      * used in @p Assert statements with the
                                      * @p ExcDataLostSync exception. If @p
                                      * strict is @p false this method returns
                                      * @p true if all sizes are within 1 of
                                      * each other (needed to allow data to be
                                      * added), with @p strict = @p true they
                                      * must be exactly equal.
                                      */

    bool deep_check (bool strict);

                                     /**
                                      * A call has been made to @p
                                      * push_back_independent when no
                                      * independent values were requested.
                                      */
    DeclException0(ExcNoIndependent);

                                     /**
                                      * This error is thrown to indicate that
                                      * the data sets appear to be out of
                                      * sync.  The class requires that the
                                      * number of dataset keys is the same as
                                      * the number of independent values sets
                                      * and mesh linked value sets. The number
                                      * of each of these is allowed to differ
                                      * by one to allow new values to be added
                                      * with out restricting the order the
                                      * user choses to do so.  Special cases
                                      * of no @p DoFHandler and no independent
                                      * values should not trigger this error.
                                      */
    DeclException0(ExcDataLostSync);


                                     /**
                                      * A method which requires access to a @p
                                      * DoFHandler to be meaningful has been
                                      * called when @p n_dofs is reported to
                                      * be zero (most likely due to default
                                      * constructor being called).  Only
                                      * independent variables may be logged
                                      * with no DoFHandler.
                                      */
    DeclException0(ExcDoFHandlerRequired);

                                     /**
                                      * @p n_dofs for the @p DoFHandler has
                                      * changed since the constructor was
                                      * called. Possibly the mesh has been
                                      * refined in some way. This suggests
                                      * that the internal dof indices stored
                                      * are no longer meaningful.
                                      */
    DeclException2(ExcDoFHandlerChanged,
                   int,
                   int,
                   << "Original n_dofs (" << arg1 << ") != current n_dofs (" << arg2 << ")");

  private:
                                     /**
                                      * Stores keys, values on the abscissa.
                                      * This will often be time, but possibly
                                      * time step, iteration etc.
                                      */
    std::vector <double> dataset_key;

                                     /**
                                      * Values that do not depend on grid
                                      * location.
                                      */
    std::vector <std::vector <double> > independent_values;

                                     /**
                                      * Saves data for each mnemonic entry.
                                      * data_store: mnemonic ->
                                      * [point_0_components point_1_components
                                      * ... point_n-1_components][key]
                                      * This format facilitates scalar mnemonics
                                      * in a vector space, because scalar mnemonics
                                      * will only have one component per point.
                                      * Vector components are strictly
                                      * FE.n_components () long.
                                      */
    std::map <std::string, std::vector <std::vector <double> > > data_store;

                                     /**
                                      * Saves a boolean representing whether
                                      * the mnemonic links to a scalar field
                                      * (true) or not (false).
                                      */
    std::map <std::string, bool> scalar_field;

                                     /**
                                      * Saves the location and other mesh info
                                      * about support points.
                                      */
    std::vector <internal::PointValueHistory::PointGeometryData <dim> >
    point_geometry_data;

                                     /**
                                      * A tempory pair used when adding to
                                      * @p data_store.
                                      */
    std::pair<std::string, std::vector <std::vector <double> > > pair_data; // could possibly be removed

                                     /**
                                      * Used to enforce @p closed state for some
                                      * methods.
                                      */
    bool closed;

                                     /**
                                      * Used to enforce @p !cleared state for
                                      * some methods.
                                      */
    bool cleared;


                                     /**
                                      * A smart pointer to the dof_handler
                                      * supplied to the constructor. This can be
                                      * released by calling @p clear().
                                      */
    SmartPointer<const DoFHandler<dim>,PointValueHistory<dim> > dof_handler;

                                     /**
                                      * Variable to check that number of dofs
                                      * in the dof_handler does not change.
                                      *
				      * A cheap way to check that the
                                      * triangulation has not been refined in
                                      * anyway since the class was set up.
                                      * Refinement will invalidate stored dof
                                      * indices linked to support points.
                                      */
    unsigned int n_dofs;
                                     /**
                                      * Stores the number of independent
                                      * variables requested.
                                      */
    unsigned int n_indep;
};


DEAL_II_NAMESPACE_CLOSE
#endif /* __dealii__point_value_history_h */

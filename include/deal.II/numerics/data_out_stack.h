// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2014 by the deal.II authors
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

#ifndef __deal2__data_out_stack_h
#define __deal2__data_out_stack_h


#include <deal.II/base/config.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/lac/vector.h>

#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim> class DoFHandler;

/**
 * This class is used to stack the output from several computations
 * into one output file by stacking the data sets in another
 * co-ordinate direction orthogonal to the space directions. The most
 * common use is to stack the results of several time steps into one
 * space-time output file, or for example to connect the results of
 * solutions of a parameter dependent equation for several parameter
 * value together into one. The interface is mostly modelled after the
 * DataOut class, see there for some more documentation.
 *
 * We will explain the concept for a time dependent problem, but
 * instead of the time any parameter can be substituted. In our
 * example, a solution of an equation is computed for each discrete
 * time level. This is then added to an object of the present class
 * and after all time levels are added, a space-time plot will be
 * written in any of the output formats supported by the base
 * class. Upon output, the (spatial) solution on each time level is
 * extended into the time direction by writing it twice, once for the
 * time level itself and once for a time equal to the time level minus
 * a given time step. These two copies are connected, to form a
 * space-time slab, with constant values in time.
 *
 * Due to the piecewise constant output in time, the written solution
 * will in general be discontinuous at discrete time levels, but the
 * output is still sufficient in most cases. More sophisticated
 * interpolations in time may be added in the future.
 *
 *
 * <h3>Example of Use</h3>
 *
 * The following little example shall illustrate the different steps
 * of use of this class. It is assumed that the finite element used is
 * composed of two components, @p u and @p v, that the solution vector
 * is named @p solution and that a vector @p error is computed which
 * contains an error indicator for each spatial cell.
 *
 * Note that unlike for the DataOut class it is necessary to first
 * declare data vectors and the names of the components before first
 * use. This is because on all time levels the same data should be
 * present to produce reasonable time-space output. The output is
 * generated with two subdivisions in each space and time direction,
 * which is suitable for quadratic finite elements in space, for
 * example.
 *
 * @code
 *   DataOutStack<dim> data_out_stack;
 *
 *                                  // first declare the vectors
 *                                  // to be used later
 *   std::vector<std::string> solution_names;
 *   solution_names.push_back ("u");
 *   solution_names.push_back ("v");
 *   data_out_stack.declare_data_vector (solution_names,
 *                                       DataOutStack<dim>::dof_vector);
 *   data_out_stack.declare_data_vector ("error",
 *                                       DataOutStack<dim>::cell_vector);
 *
 *                                  // now do computations
 *   for (double parameter=0; ...)
 *     {
 *       DoFHandler<dim,spacedim> dof_handler;
 *       ...                        // compute something
 *
 *                                  // now for output
 *       data_out_stack.new_parameter_value (parameter,
 *                                           delta_parameter);
 *       data_out_stack.attach_dof_handler (dof_handler);
 *       data_out_stack.add_data_vector (solution, solution_names);
 *       data_out_stack.add_data_vector (error, "error");
 *       data_out_stack.build_patches (2);
 *       data_out_stack.finish_parameter_value ();
 *     };
 * @endcode
 *
 * @ingroup output
 * @author Wolfgang Bangerth, 1999
 */
template <int dim, int spacedim=dim, class DH = DoFHandler<dim,spacedim> >
class DataOutStack : public DataOutInterface<dim+1>
{
public:
  /**
   * Data type declaring the two types
   * of vectors which are used in this
   * class.
   */
  enum VectorType { cell_vector, dof_vector };

  /**
   * Destructor. Only declared to make
   * it @p virtual.
   */
  virtual ~DataOutStack ();

  /**
   * Start the next set of data for a
   * specific parameter value. The argument
   * @p parameter_step denotes the interval
   * (in backward direction, counted from
   * @p parameter_value) with which the
   * output will be extended in parameter
   * direction, i.e. orthogonal to the
   * space directions.
   */
  void new_parameter_value (const double parameter_value,
                            const double parameter_step);

  /**
   * Attach the DoF handler for the
   * grid and data associated with the
   * parameter previously set by
   * @p new_parameter_value.
   *
   * This has to happen before adding
   * data vectors for the present
   * parameter value.
   */
  void attach_dof_handler (const DH &dof_handler);

  /**
   * Declare a data vector. The @p vector_type
   * argument determines whether the data
   * vector will be considered as DoF or
   * cell data.
   *
   * This version may be called if the
   * finite element presently used by the
   * DoFHandler (and previously attached
   * to this object) has only one component
   * and therefore only one name needs to
   * be given.
   */
  void declare_data_vector (const std::string &name,
                            const VectorType   vector_type);

  /**
   * Declare a data vector. The @p vector_type
   * argument determines whether the data
   * vector will be considered as DoF or
   * cell data.
   *
   * This version must be called if the
   * finite element presently used by the
   * DoFHandler (and previously attached
   * to this object) has more than one
   * component and therefore more than one
   * name needs to be given. However, you
   * can also call this function with a
   * <tt>std::vector@<std::string@></tt> containing only one
   * element if the finite element has
   * only one component.
   */
  void declare_data_vector (const std::vector<std::string> &name,
                            const VectorType                vector_type);


  /**
   * Add a data vector for the presently
   * set value of the parameter.
   *
   * This version may be called if the
   * finite element presently used by the
   * DoFHandler (and previously attached
   * to this object) has only one component
   * and therefore only one name needs to
   * be given.
   *
   * If @p vec is a vector with
   * multiple components this
   * function will generate
   * distinct names for all
   * components by appending an
   * underscore and the number of
   * each component to @p name
   *
   * The data vector must have been
   * registered using the @p declare_data_vector
   * function before actually using it the
   * first time.
   *
   * Note that a copy of this vector is
   * stored until @p finish_parameter_value
   * is called the next time, so if you are
   * short of memory you may want to call
   * this function only after all
   * computations involving large matrices
   * are already done.
   */
  template <typename number>
  void add_data_vector (const Vector<number> &vec,
                        const std::string    &name);

  /**
   * Add a data vector for the presently
   * set value of the parameter.
   *
   * This version must be called if the
   * finite element presently used by the
   * DoFHandler (and previously attached
   * to this object) has more than one
   * component and therefore more than one
   * name needs to be given. However, you
   * can also call this function with a
   * <tt>std::vector@<std::string@></tt> containing only one
   * element if the finite element has
   * only one component.
   *
   * The data vector must have been
   * registered using the @p declare_data_vector
   * function before actually using it the
   * first time.
   *
   * Note that a copy of this vector is
   * stored until @p finish_parameter_value
   * is called the next time, so if you are
   * short of memory you may want to call
   * this function only after all
   * computations involving large matrices
   * are already done.
   */
  template <typename number>
  void add_data_vector (const Vector<number>           &vec,
                        const std::vector<std::string> &names);

  /**
   * Actually build the patches for output
   * by the base classes from the data
   * stored in the data vectors and using
   * the previously attached DoFHandler
   * object.
   *
   * By @p n_subdivisions you can decide
   * into how many subdivisions (in each
   * space and parameter direction) each
   * patch is divided. This is useful
   * if higher order elements are used.
   * Note however, that the number of
   * subdivisions in parameter direction
   * is always the same as the one is space
   * direction for technical reasons.
   */
  void build_patches (const unsigned int n_subdivisions = 0);

  /**
   * Release all data that is no
   * more needed once @p build_patches
   * was called and all other transactions
   * for a given parameter value are done.
   *
   * Couterpart of @p new_parameter_value.
   */
  void finish_parameter_value ();

  /**
   * Clear all data presently stored
   * in this object.
   */
  void clear ();

  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   */
  std::size_t memory_consumption () const;

  /**
   * Exception
   */
  DeclException0 (ExcNoDoFHandlerSelected);
  /**
   * Exception
   */
  DeclException3 (ExcInvalidVectorSize,
                  int, int, int,
                  << "The vector has size " << arg1
                  << " but the DoFHandler objects says there are " << arg2
                  << " degrees of freedom and there are " << arg3
                  << " active cells.");
  /**
   * Exception
   */
  DeclException2 (ExcInvalidCharacter,
                  std::string, size_t,
                  << "Please use only the characters [a-zA-Z0-9_<>()] for" << std::endl
                  << "description strings since some graphics formats will only accept these."
                  << std::endl
                  << "The string you gave was <" << arg1
                  << ">, the invalid character is <" << arg1[arg2]
                  << ">." << std::endl);
  /**
   * Exception
   */
  DeclException2 (ExcInvalidNumberOfNames,
                  int, int,
                  << "You have to give one name per component in your "
                  << "data vector. The number you gave was " << arg1
                  << ", but the number of components is " << arg2);
  /**
   * Exception
   */
  DeclException1 (ExcVectorNotDeclared,
                  std::string,
                  << "The data vector for which the first component has the name "
                  << arg1 << " has not been declared before.");
  /**
   * Exception
   */
  DeclException0 (ExcDataNotCleared);
  /**
   * Exception
   */
  DeclException0 (ExcDataAlreadyAdded);
  /**
   * Exception
   */
  DeclException1 (ExcNameAlreadyUsed,
                  std::string,
                  << "You tried to declare a component of a data vector with "
                  << "the name <" << arg1 << ">, but that name is already used.");
  /**
   * Exception
   */
  DeclException1 (ExcInvalidNumberOfSubdivisions,
                  int,
                  << "The number of subdivisions per patch, " << arg1
                  << ", is not valid.");

private:
  /**
   * Present parameter value.
   */
  double                               parameter;

  /**
   * Present parameter step, i.e.
   * length of the parameter interval
   * to be written next.
   */
  double                               parameter_step;

  /**
   * DoF handler to be used for the data
   * corresponding to the present parameter
   * value.
   */
  SmartPointer<const DH,DataOutStack<dim,spacedim,DH> > dof_handler;

  /**
   * List of patches of all past and
   * present parameter value data sets.
   */
  std::vector< dealii::DataOutBase::Patch<dim+1,dim+1> >   patches;

  /**
   * Structure holding data vectors
   * (cell and dof data) for the present
   * parameter value.
   */
  struct DataVector
  {
    /**
     * Data vector.
     */
    Vector<double> data;

    /**
     * Names of the different components
     * within each such data set.
     */
    std::vector<std::string> names;

    /**
     * Determine an estimate for
     * the memory consumption (in
     * bytes) of this object.
     */
    std::size_t memory_consumption () const;
  };

  /**
   * List of DoF data vectors.
   */
  std::vector<DataVector> dof_data;

  /**
   * List of cell data vectors.
   */
  std::vector<DataVector> cell_data;

  /**
   * This is the function through
   * which derived classes propagate
   * preprocessed data in the form of
   * Patch structures (declared in
   * the base class DataOutBase) to
   * the actual output function.
   */
  virtual const std::vector< dealii::DataOutBase::Patch<dim+1,dim+1> > & get_patches () const;


  /**
   * Virtual function through
   * which the names of data sets are
   * obtained by the output functions
   * of the base class.
   */
  virtual std::vector<std::string> get_dataset_names () const;
};


DEAL_II_NAMESPACE_CLOSE

#endif

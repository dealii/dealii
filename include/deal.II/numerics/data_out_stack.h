// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_data_out_stack_h
#define dealii_data_out_stack_h


#include <deal.II/base/config.h>

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/observer_pointer.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out_dof_data.h>

#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class DoFHandler;
#endif

/**
 * This class is used to stack the output from several computations into one
 * output file by stacking the data sets in another coordinate direction
 * orthogonal to the space directions. The most common use is to stack the
 * results of several time steps into one space-time output file, or for
 * example to connect the results of solutions of a parameter dependent
 * equation for several parameter value together into one. The interface is
 * mostly modelled after the DataOut class, see there for some more
 * documentation. The class is used in step-78.
 *
 * We will explain the concept for a time dependent problem, but instead of
 * the time any parameter can be substituted. In our example, a solution of an
 * equation is computed for each discrete time level. This is then added to an
 * object of the present class and after all time levels are added, a
 * space-time plot will be written in any of the output formats supported by the
 * base class. Upon output, the (spatial) solution on each time level is
 * extended into the time direction by writing it twice, once for the time
 * level itself and once for a time equal to the time level minus a given time
 * step. These two copies are connected, to form a space-time slab, with
 * constant values in time.
 *
 * Due to the piecewise constant output in time, the written solution will in
 * general be discontinuous at discrete time levels, but the output is still
 * sufficient in most cases. More sophisticated interpolations in time may be
 * added in the future.
 *
 *
 * <h3>Example of use</h3>
 *
 * The following little example shall illustrate the different steps of use of
 * this class. It is assumed that the finite element used is composed of two
 * components, @p u and @p v, that the solution vector is named @p solution
 * and that a vector @p error is computed which contains an error indicator
 * for each spatial cell.
 *
 * Note that unlike for the DataOut class it is necessary to first declare
 * data vectors and the names of the components before first use. This is
 * because on all time levels the same data should be present to produce
 * reasonable time-space output. The output is generated with two subdivisions
 * in each space and time direction, which is suitable for quadratic finite
 * elements in space, for example.
 *
 * @code
 *   DataOutStack<dim> data_out_stack;
 *
 *                                  // first declare the vectors
 *                                  // to be used later
 *   std::vector<std::string> solution_names;
 *   solution_names.emplace_back ("u");
 *   solution_names.emplace_back ("v");
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
 */
template <int dim, int spacedim = dim>
class DataOutStack : public DataOutInterface<dim + 1, spacedim + 1>
{
  static_assert(dim < 3,
                "Because this class stacks data into the (dim+1)st "
                "dimension to create graphical output, it only works for "
                "dim<3.");
  static_assert(dim == spacedim,
                "This class is not implemented for dim != spacedim.");

public:
  /**
   * Dimension parameters for the patches.
   */
  static constexpr int patch_dim      = dim + 1;
  static constexpr int patch_spacedim = spacedim + 1;

  /**
   * Data type declaring the two types of vectors which are used in this
   * class.
   */
  enum VectorType
  {
    /**
     * The data describes one value for each cell.
     */
    cell_vector,
    /**
     * The data describes one value for each DoF.
     */
    dof_vector
  };

  /**
   * Destructor. Only declared to make it @p virtual.
   */
  virtual ~DataOutStack() override = default;

  /**
   * Start the next set of data for a specific parameter value. The argument
   * @p parameter_step denotes the interval (in backward direction, counted
   * from @p parameter_value) with which the output will be extended in
   * parameter direction, i.e. orthogonal to the space directions.
   */
  void
  new_parameter_value(const double parameter_value,
                      const double parameter_step);

  /**
   * Attach the DoF handler for the grid and data associated with the
   * parameter previously set by @p new_parameter_value.
   *
   * This has to happen before adding data vectors for the present parameter
   * value.
   */
  void
  attach_dof_handler(const DoFHandler<dim, spacedim> &dof_handler);

  /**
   * Declare a data vector. The @p vector_type argument determines whether the
   * data vector will be considered as DoF or cell data.
   *
   * This version may be called if the finite element presently used by the
   * DoFHandler (and previously attached to this object) has only one
   * component and therefore only one name needs to be given.
   */
  void
  declare_data_vector(const std::string &name, const VectorType vector_type);

  /**
   * Declare a data vector. The @p vector_type argument determines whether the
   * data vector will be considered as DoF or cell data.
   *
   * This version must be called if the finite element presently used by the
   * DoFHandler (and previously attached to this object) has more than one
   * component and therefore more than one name needs to be given. However,
   * you can also call this function with a
   * <tt>std::vector@<std::string@></tt> containing only one element if the
   * finite element has only one component.
   */
  void
  declare_data_vector(const std::vector<std::string> &name,
                      const VectorType                vector_type);


  /**
   * Add a data vector for the presently set value of the parameter.
   *
   * This version may be called if the finite element presently used by the
   * DoFHandler (and previously attached to this object) has only one
   * component and therefore only one name needs to be given.
   *
   * If @p vec is a vector with multiple components this function will
   * generate distinct names for all components by appending an underscore and
   * the number of each component to @p name
   *
   * The data vector must have been registered using the @p
   * declare_data_vector function before actually using it the first time.
   *
   * Note that a copy of this vector is stored until @p finish_parameter_value
   * is called the next time, so if you are short of memory you may want to
   * call this function only after all computations involving large matrices
   * are already done.
   */
  template <typename number>
  void
  add_data_vector(const Vector<number> &vec, const std::string &name);

  /**
   * Add a data vector for the presently set value of the parameter.
   *
   * This version must be called if the finite element presently used by the
   * DoFHandler (and previously attached to this object) has more than one
   * component and therefore more than one name needs to be given. However,
   * you can also call this function with a
   * <tt>std::vector@<std::string@></tt> containing only one element if the
   * finite element has only one component.
   *
   * The data vector must have been registered using the @p
   * declare_data_vector function before actually using it the first time.
   *
   * Note that a copy of this vector is stored until @p finish_parameter_value
   * is called the next time, so if you are short of memory you may want to
   * call this function only after all computations involving large matrices
   * are already done.
   */
  template <typename number>
  void
  add_data_vector(const Vector<number>           &vec,
                  const std::vector<std::string> &names);

  /**
   * This is the central function of this class since it builds the list of
   * patches to be written by the low-level functions of the base class. A
   * patch is, in essence, some intermediate representation of the data on
   * each cell of a triangulation and DoFHandler object that can then be used
   * to write files in some format that is readable by visualization programs.
   *
   * You can find an overview of the use of this function in the general
   * documentation of this class. An example is also provided in the
   * documentation of this class's base class DataOut_DoFData.
   *
   * @param n_subdivisions See DataOut::build_patches() for an extensive
   * description of this parameter. The number of subdivisions is always one
   * in the direction of the time-like parameter used by this class.
   */
  void
  build_patches(const unsigned int n_subdivisions = 0);

  /**
   * Release all data that is no more needed once @p build_patches was called
   * and all other transactions for a given parameter value are done.
   *
   * Counterpart of @p new_parameter_value.
   */
  void
  finish_parameter_value();

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  std::size_t
  memory_consumption() const;

  /**
   * Exception
   */
  DeclException1(
    ExcVectorNotDeclared,
    std::string,
    << "The data vector for which the first component has the name " << arg1
    << " has not been added before.");
  /**
   * Exception
   */
  DeclExceptionMsg(ExcDataNotCleared,
                   "You cannot start a new time/parameter step before calling "
                   "finish_parameter_value() on the previous step.");
  /**
   * Exception
   */
  DeclExceptionMsg(
    ExcDataAlreadyAdded,
    "You cannot declare additional vectors after already calling "
    "build_patches(). All data vectors need to be declared "
    "before you call this function the first time.");
  /**
   * Exception
   */
  DeclException1(ExcNameAlreadyUsed,
                 std::string,
                 << "You tried to declare a component of a data vector with "
                 << "the name <" << arg1
                 << ">, but that name is already used.");

private:
  /**
   * Present parameter value.
   */
  double parameter;

  /**
   * Present parameter step, i.e. length of the parameter interval to be
   * written next.
   */
  double parameter_step;

  /**
   * DoF handler to be used for the data corresponding to the present
   * parameter value.
   */
  ObserverPointer<const DoFHandler<dim, spacedim>, DataOutStack<dim, spacedim>>
    dof_handler;

  /**
   * List of patches of all past and present parameter value data sets.
   */
  std::vector<dealii::DataOutBase::Patch<patch_dim, patch_spacedim>> patches;

  /**
   * Structure holding data vectors (cell and dof data) for the present
   * parameter value.
   */
  struct DataVector
  {
    /**
     * Data vector.
     */
    Vector<double> data;

    /**
     * Names of the different components within each such data set.
     */
    std::vector<std::string> names;

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    std::size_t
    memory_consumption() const;
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
   * This is the function through which derived classes propagate preprocessed
   * data in the form of Patch structures (declared in the base class
   * DataOutBase) to the actual output function.
   */
  virtual const std::vector<
    dealii::DataOutBase::Patch<DataOutStack<dim, spacedim>::patch_dim,
                               DataOutStack<dim, spacedim>::patch_spacedim>> &
  get_patches() const override;


  /**
   * Virtual function through which the names of data sets are obtained by the
   * output functions of the base class.
   */
  virtual std::vector<std::string>
  get_dataset_names() const override;
};


DEAL_II_NAMESPACE_CLOSE

#endif

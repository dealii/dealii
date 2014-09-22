// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2014 by the deal.II authors
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

#ifndef __deal2__data_postprocessor_h
#define __deal2__data_postprocessor_h



#include <deal.II/base/subscriptor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/lac/vector.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/numerics/data_component_interpretation.h>

#include <vector>
#include <string>

DEAL_II_NAMESPACE_OPEN


/**
 * This class provides an interface to compute derived quantities from
 * a solution that can then be output in graphical formats for visualization,
 * using facilities such as the DataOut class.
 *
 * For the (graphical) output of a FE solution one frequently wants to include
 * derived quantities, which are calculated from the values of the solution
 * and possibly the first and second derivates of the solution. Examples are
 * the calculation of Mach numbers from velocity and density in supersonic flow
 * computations, or the computation of the magnitude of a complex-valued
 * solution as demonstrated in step-29. Other uses are shown in step-32 and
 * step-33. This class offers the interface to perform such
 * postprocessing. Given the values and derivatives of the solution at those
 * points where we want to generated output, the functions of this class can
 * be overloaded to compute new quantities.
 *
 * A data vector and an object of a derived class can be given to the
 * DataOut::add_data_vector() function, which will write the derived
 * quantities instead of the provided data to the output file. Note, that the
 * DataPostprocessor has to live until DataOut::build_patches has been
 * called. DataOutFaces and DataOutRotation can be used as well.
 *
 * In order not to perform needless calculations, DataPostprocessor
 * has to provide information which input data is needed for the
 * calculation of the derived quantities, i.e. whether it needs the
 * values, the first derivative and/or the second derivative of the
 * provided data. DataPostprocessor objects which are used in
 * combination with a DataOutFaces object can also ask for the normal
 * vectors at each point. The information which data is needed has to
 * be provided via the UpdateFlags returned by the virtual function
 * get_needed_update_flags(). It is your responsibility to use only
 * those values which were updated in the calculation of derived
 * quantities. The DataOut object will provide references to the
 * requested data in the call to compute_derived_quantities_scalar()
 * or compute_derived_quantities_vector() (DataOut decides which of
 * the two functions to call depending on whether the finite element
 * in use has only a single, or multiple vector components; note that
 * this is only determined by the number of components in the finite
 * element in use, and not by whether the data computed by a class
 * derived from the current one is scalar or vector valued).
 *
 * Furthermore, derived classes have to implement the get_names()
 * function, where the number of output variables returned
 * by the latter function has to match the size of the vector returned by the
 * former. Furthermore, this number has to match the number of computed
 * quantities, of course.
 *
 *
 * <h3>Use in simpler cases</h3>
 *
 * Deriving from the current class allows to implement very general postprocessors.
 * For example, in the step-32 program, we implement a postprocessor that
 * takes a solution that consists of velocity, pressure and temperature (dim+2
 * components) and computes a variety of output quantities, some of which
 * are vector valued and some of which are scalar. On the other hand,
 * in step-29 we implement a postprocessor that only computes the magnitude
 * of a complex number given by a two-component finite element. It seems silly
 * to have to implement four virtual functions for this
 * (compute_derived_quantities_scalar() or compute_derived_quantities_vector(),
 * get_names(), get_update_flags() and get_data_component_interpretation()).
 *
 * To this end there are two classes DataPostprocessorScalar and
 * DataPostprocessorVector that are meant to be used if the output quantity
 * is either a single scalar or a single vector (here used meaning to have
 * exactly dim components). When using these classes, one only has to write a
 * constructor that passes the name of the output variable and the update
 * flags to the constructor of the base class and overload the function
 * that actually computes the results.
 *
 * @ingroup output
 * @author Tobias Leicht, 2007
 */
template <int dim>
class DataPostprocessor: public Subscriptor
{
public:
  /**
   * Virtual desctructor for safety. Does not
   * do anything.
   */
  virtual ~DataPostprocessor ();

  /**
   * @deprecated
   *
   * This function only exists for backward
   * compatibility as this is the interface
   * provided by previous versions of the
   * library. The default implementation of
   * the other function of same name below
   * calls this function by simply dropping
   * the argument that denotes the
   * evaluation points. Since this function
   * might at one point go away, you should
   * overload the other function instead.
   */
  virtual
  void
  compute_derived_quantities_scalar (const std::vector<double>         &uh,
                                     const std::vector<Tensor<1,dim> > &duh,
                                     const std::vector<Tensor<2,dim> > &dduh,
                                     const std::vector<Point<dim> >    &normals,
                                     std::vector<Vector<double> >      &computed_quantities) const DEAL_II_DEPRECATED;

  /**
   * This is the main function which actually
   * performs the postprocessing. The last
   * argument is a reference to the
   * postprocessed data which has correct
   * size already and must be filled by this
   * function. @p uh is a reference to a
   * vector of data values at all points, @p
   * duh the same for gradients, @p dduh for
   * second derivatives and @p normals is a
   * reference to the normal vectors. Note,
   * that the last four references will only
   * contain valid data, if the respective
   * flags are returned by @p
   * get_needed_update_flags, otherwise those
   * vectors will be in an unspecified
   * state. @p normals will always be an
   * empty vector when working on cells, not
   * on faces.
   *
   * This function is called when
   * the original data vector
   * represents scalar data,
   * i.e. the finite element in use
   * has only a single vector
   * component.
   */
  virtual
  void
  compute_derived_quantities_scalar (const std::vector<double>         &uh,
                                     const std::vector<Tensor<1,dim> > &duh,
                                     const std::vector<Tensor<2,dim> > &dduh,
                                     const std::vector<Point<dim> >    &normals,
                                     const std::vector<Point<dim> >    &evaluation_points,
                                     std::vector<Vector<double> >      &computed_quantities) const;

  /**
   * @deprecated
   *
   * This function only exists for backward
   * compatibility as this is the interface
   * provided by previous versions of the
   * library. The default implementation of
   * the other function of same name below
   * calls this function by simply dropping
   * the argument that denotes the
   * evaluation points. Since this function
   * might at one point go away, you should
   * overload the other function instead.
   */
  virtual
  void
  compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                     const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                     const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                     const std::vector<Point<dim> >                  &normals,
                                     std::vector<Vector<double> >                    &computed_quantities) const DEAL_II_DEPRECATED;

  /**
   * Same as the
   * compute_derived_quantities_scalar()
   * function, but this function is called
   * when the original data vector
   * represents vector data, i.e. the
   * finite element in use has multiple
   * vector components.
   */
  virtual
  void
  compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                     const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                     const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                     const std::vector<Point<dim> >                  &normals,
                                     const std::vector<Point<dim> >                  &evaluation_points,
                                     std::vector<Vector<double> >                    &computed_quantities) const;

  /**
   * Return the vector of strings describing
   * the names of the computed quantities.
   */
  virtual std::vector<std::string> get_names () const = 0;

  /**
   * This functions returns
   * information about how the
   * individual components of
   * output files that consist of
   * more than one data set are to
   * be interpreted.
   *
   * For example, if one has a finite
   * element for the Stokes equations in
   * 2d, representing components (u,v,p),
   * one would like to indicate that the
   * first two, u and v, represent a
   * logical vector so that later on when
   * we generate graphical output we can
   * hand them off to a visualization
   * program that will automatically know
   * to render them as a vector field,
   * rather than as two separate and
   * independent scalar fields.
   *
   * The default implementation of this
   * function returns a vector of values
   * DataComponentInterpretation::component_is_scalar,
   * indicating that all output components
   * are independent scalar
   * fields. However, if a derived class
   * produces data that represents vectors,
   * it may return a vector that contains
   * values
   * DataComponentInterpretation::component_is_part_of_vector. In
   * the example above, one would return a
   * vector with components
   * (DataComponentInterpretation::component_is_part_of_vector,
   * DataComponentInterpretation::component_is_part_of_vector,
   * DataComponentInterpretation::component_is_scalar)
   * for (u,v,p).
   */
  virtual
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  get_data_component_interpretation () const;

  /**
   * Return, which data has to be provided to
   * compute the derived quantities. This has
   * to be a combination of @p update_values,
   * @p update_gradients and @p
   * update_hessians. If the
   * DataPostprocessor is to be used in
   * combination with DataOutFaces, you may
   * also ask for a update of normals via the
   * @p update_normal_vectors flag.
   */
  virtual UpdateFlags get_needed_update_flags () const = 0;
};



/**
 * This class provides a simpler interface to the functionality offered by
 * the DataPostprocessor class in case one wants to compute only a
 * single scalar quantity from the finite element field passed to the
 * DataOut class. For this particular case, it is clear what the returned
 * value of DataPostprocessor::get_data_component_interpretation() should
 * be and we pass the values returned by get_names() and get_needed_update_flags()
 * to the constructor so that derived classes do not have to implement these
 * functions by hand.
 *
 * All derived classes have to do is implement a constructor and overload
 * either DataPostprocessor::compute_derived_quantities_scalar() or
 * DataPostprocessor::compute_derived_quantities_vector().
 *
 * An example of how this class can be used can be found in step-29.
 *
 * @ingroup output
 * @author Wolfgang Bangerth, 2011
 */
template <int dim>
class DataPostprocessorScalar : public DataPostprocessor<dim>
{
public:
  /**
   * Constructor. Take the name of the single scalar variable
   * computed by classes derived from the current one, as well
   * as the update flags necessary to compute this quantity.
   *
   * @param name The name by which the scalar variable
   *   computed by this class should be made available in
   *   graphical output files.
   * @param update_flags This has
   * to be a combination of @p update_values,
   * @p update_gradients and @p
   * update_hessians. If the
   * DataPostprocessor is to be used in
   * combination with DataOutFaces, you may
   * also ask for a update of normals via the
   * @p update_normal_vectors flag.
   **/
  DataPostprocessorScalar (const std::string &name,
                           const UpdateFlags  update_flags);

  /**
   * Return the vector of strings describing
   * the names of the computed quantities.
   * Given the purpose of this class, this
   * is a vector with a single entry equal
   * to the name given to the constructor.
   */
  virtual std::vector<std::string> get_names () const;

  /**
   * This functions returns
   * information about how the
   * individual components of
   * output files that consist of
   * more than one data set are to
   * be interpreted. Since the current
   * class is meant to be used for a
   * single scalar result variable,
   * the returned value is obviously
   * DataComponentInterpretation::component_is_scalar.
   */
  virtual
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  get_data_component_interpretation () const;

  /**
   * Return, which data has to be provided to
   * compute the derived quantities.
   * The flags returned here are the ones
   * passed to the constructor of this
   * class.
   */
  virtual UpdateFlags get_needed_update_flags () const;

private:
  /**
   * Copies of the two arguments given to the constructor of this
   * class.
   */
  const std::string name;
  const UpdateFlags update_flags;
};



/**
 * This class provides a simpler interface to the functionality offered by
 * the DataPostprocessor class in case one wants to compute only a
 * single vector quantity (defined as having exactly dim components)
 * from the finite element field passed to the
 * DataOut class. For this particular case, it is clear what the returned
 * value of DataPostprocessor::get_data_component_interpretation() should
 * be and we pass the values returned by get_names() and get_needed_update_flags()
 * to the constructor so that derived classes do not have to implement these
 * functions by hand.
 *
 * All derived classes have to do is implement a constructor and overload
 * either DataPostprocessor::compute_derived_quantities_scalar() or
 * DataPostprocessor::compute_derived_quantities_vector().
 *
 * An example of how the closely related class DataPostprocessorScalar is
 * used can be found in step-29.
 *
 * @ingroup output
 * @author Wolfgang Bangerth, 2011
 */
template <int dim>
class DataPostprocessorVector : public DataPostprocessor<dim>
{
public:
  /**
   * Constructor. Take the name of the single vector variable
   * computed by classes derived from the current one, as well
   * as the update flags necessary to compute this quantity.
   *
   * @param name The name by which the vector variable
   *   computed by this class should be made available in
   *   graphical output files.
   * @param update_flags This has
   * to be a combination of @p update_values,
   * @p update_gradients and @p
   * update_hessians. If the
   * DataPostprocessor is to be used in
   * combination with DataOutFaces, you may
   * also ask for a update of normals via the
   * @p update_normal_vectors flag.
   **/
  DataPostprocessorVector (const std::string &name,
                           const UpdateFlags  update_flags);

  /**
   * Return the vector of strings describing
   * the names of the computed quantities.
   * Given the purpose of this class, this
   * is a vector with dim entries all equal
   * to the name given to the constructor.
   */
  virtual std::vector<std::string> get_names () const;

  /**
   * This functions returns
   * information about how the
   * individual components of
   * output files that consist of
   * more than one data set are to
   * be interpreted. Since the current
   * class is meant to be used for a
   * single vector result variable,
   * the returned value is obviously
   * DataComponentInterpretation::component_is_part
   * repeated dim times.
   */
  virtual
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  get_data_component_interpretation () const;

  /**
   * Return which data has to be provided to
   * compute the derived quantities.
   * The flags returned here are the ones
   * passed to the constructor of this
   * class.
   */
  virtual UpdateFlags get_needed_update_flags () const;

private:
  /**
   * Copies of the two arguments given to the constructor of this
   * class.
   */
  const std::string name;
  const UpdateFlags update_flags;
};


DEAL_II_NAMESPACE_CLOSE

#endif

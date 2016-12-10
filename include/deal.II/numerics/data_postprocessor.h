// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2016 by the deal.II authors
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

#ifndef dealii__data_postprocessor_h
#define dealii__data_postprocessor_h



#include <deal.II/base/subscriptor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/numerics/data_component_interpretation.h>

#include <vector>
#include <string>

DEAL_II_NAMESPACE_OPEN


/**
 * A namespace for data structures that are going to be passed from
 * DataOut to the member functions of DataPostprocessor.
 */
namespace DataPostprocessorInputs
{
  /**
   * A structure that is used to pass information to
   * DataPostprocessor::evaluate_scalar_field(). It contains
   * the values and (if requested) derivatives of a scalar solution
   * variable at the evaluation points on a cell or face. If appropriate,
   * it also contains the normal vectors to the geometry on which output
   * is generated, at these evaluation points.
   */
  template <int spacedim>
  struct Scalar
  {
    /**
     * An array of values of the (scalar) solution at each of the evaluation
     * points used to create graphical output from one cell, face, or other
     * object.
     */
    std::vector<double>               solution_values;

    /**
     * An array of gradients of the (scalar) solution at each of the evaluation
     * points used to create graphical output from one cell, face, or other
     * object.
     *
     * This array is only filled if a user-derived class overloads the
     * DataPostprocessor::get_needed_update_flags(), and the function
     * returns (possibly among other flags)
     * UpdateFlags::update_gradients.  Alternatively, a class derived
     * from DataPostprocessorScalar or DataPostprocessorVector may
     * pass this flag to the constructor of DataPostprocessorScalar or
     * DataPostprocessorVector.
     */
    std::vector<Tensor<1, spacedim> > solution_gradients;

    /**
     * An array of second derivatives of the (scalar) solution at each of the evaluation
     * points used to create graphical output from one cell, face, or other
     * object.
     *
     * This array is only filled if a user-derived class overloads the
     * DataPostprocessor::get_needed_update_flags(), and the function
     * returns (possibly among other flags)
     * UpdateFlags::update_hessians.  Alternatively, a class derived
     * from DataPostprocessorScalar or DataPostprocessorVector may
     * pass this flag to the constructor of DataPostprocessorScalar or
     * DataPostprocessorVector.
     */
    std::vector<Tensor<2, spacedim> > solution_hessians;

    /**
     * An array of vectors normal to the faces of cells, evaluated at the points
     * at which we are generating graphical output. This array is only used by
     * the DataOutFaces class, and is left empty by all other classes for
     * which the DataPostprocessor framework can be used. In the case of
     * DataOutFaces, the array contains the outward normal vectors to the
     * face, seen from the interior of the cell.
     *
     * This array is only filled if a user-derived class overloads the
     * DataPostprocessor::get_needed_update_flags(), and the function
     * returns (possibly among other flags)
     * UpdateFlags::update_normal_vectors.  Alternatively, a class
     * derived from DataPostprocessorScalar or DataPostprocessorVector
     * may pass this flag to the constructor of
     * DataPostprocessorScalar or DataPostprocessorVector.
     */
    std::vector<Tensor<1, spacedim> > normals;

    /**
     * An array of coordinates corresponding to the locations at which
     * we are generating graphical output on one cell.
     *
     * This array is only filled if a user-derived class overloads the
     * DataPostprocessor::get_needed_update_flags(), and the function
     * returns (possibly among other flags)
     * UpdateFlags::update_quadrature_points.  Alternatively, a class
     * derived from DataPostprocessorScalar or DataPostprocessorVector
     * may pass this flag to the constructor of
     * DataPostprocessorScalar or DataPostprocessorVector.
     */
    std::vector<Point<spacedim> >     evaluation_points;
  };



  /**
   * A structure that is used to pass information to
   * DataPostprocessor::evaluate_vector_field(). It contains
   * the values and (if requested) derivatives of a vector-valued solution
   * variable at the evaluation points on a cell or face. If appropriate,
   * it also contains the normal vectors to the geometry on which output
   * is generated, at these evaluation points.
   */
  template <int spacedim>
  struct Vector
  {
    /**
     * An array of values of a vector-valued solution at each of the evaluation
     * points used to create graphical output from one cell, face, or other
     * object.
     *
     * The outer vector runs over the evaluation points, whereas the inner
     * vector runs over the components of the finite element field for which
     * output will be generated.
     */
    std::vector<dealii::Vector<double> >            solution_values;

    /**
     * An array of gradients of a vector-valued solution at each of the evaluation
     * points used to create graphical output from one cell, face, or other
     * object.
     *
     * The outer vector runs over the evaluation points, whereas the inner
     * vector runs over the components of the finite element field for which
     * output will be generated.
     *
     * This array is only filled if a user-derived class overloads the
     * DataPostprocessor::get_needed_update_flags(), and the function
     * returns (possibly among other flags)
     * UpdateFlags::update_gradients.  Alternatively, a class derived
     * from DataPostprocessorScalar or DataPostprocessorVector may
     * pass this flag to the constructor of DataPostprocessorScalar or
     * DataPostprocessorVector.
     */
    std::vector<std::vector<Tensor<1, spacedim> > > solution_gradients;

    /**
     * An array of second derivatives of a vector-valued solution at each of the evaluation
     * points used to create graphical output from one cell, face, or other
     * object.
     *
     * The outer vector runs over the evaluation points, whereas the inner
     * vector runs over the components of the finite element field for which
     * output will be generated.
     *
     * This array is only filled if a user-derived class overloads the
     * DataPostprocessor::get_needed_update_flags(), and the function
     * returns (possibly among other flags)
     * UpdateFlags::update_hessians.  Alternatively, a class derived
     * from DataPostprocessorScalar or DataPostprocessorVector may
     * pass this flag to the constructor of DataPostprocessorScalar or
     * DataPostprocessorVector.
     */
    std::vector<std::vector<Tensor<2, spacedim> > > solution_hessians;

    /**
     * An array of vectors normal to the faces of cells, evaluated at the points
     * at which we are generating graphical output. This array is only used by
     * the DataOutFaces class, and is left empty by all other classes for
     * which the DataPostprocessor framework can be used. In the case of
     * DataOutFaces, the array contains the outward normal vectors to the
     * face, seen from the interior of the cell.
     *
     * This array is only filled if a user-derived class overloads the
     * DataPostprocessor::get_needed_update_flags(), and the function
     * returns (possibly among other flags)
     * UpdateFlags::update_normal_vectors.  Alternatively, a class
     * derived from DataPostprocessorScalar or DataPostprocessorVector
     * may pass this flag to the constructor of
     * DataPostprocessorScalar or DataPostprocessorVector.
     */
    std::vector<Tensor<1, spacedim> >               normals;

    /**
     * An array of coordinates corresponding to the locations at which
     * we are generating graphical output on one cell.
     *
     * This array is only filled if a user-derived class overloads the
     * DataPostprocessor::get_needed_update_flags(), and the function
     * returns (possibly among other flags)
     * UpdateFlags::update_quadrature_points.  Alternatively, a class
     * derived from DataPostprocessorScalar or DataPostprocessorVector
     * may pass this flag to the constructor of
     * DataPostprocessorScalar or DataPostprocessorVector.
     */
    std::vector<Point<spacedim> >                   evaluation_points;
  };
}


/**
 * This class provides an interface to compute derived quantities from a
 * solution that can then be output in graphical formats for visualization,
 * using facilities such as the DataOut class.
 *
 * For the (graphical) output of a FE solution one frequently wants to include
 * derived quantities, which are calculated from the values of the solution
 * and possibly the first and second derivatives of the solution. Examples are
 * the calculation of Mach numbers from velocity and density in supersonic
 * flow computations, or the computation of the magnitude of a complex-valued
 * solution as demonstrated in step-29. Other uses are shown in step-32 and
 * step-33. This class offers the interface to perform such postprocessing.
 * Given the values and derivatives of the solution at those points where we
 * want to generated output, the functions of this class can be overloaded to
 * compute new quantities.
 *
 * A data vector and an object of a class derived from the current one can be
 * given to the DataOut::add_data_vector() function (and similarly for
 * DataOutRotation and DataOutFaces). This will cause DataOut::build_patches()
 * to compute the derived quantities instead of using the data provided by the
 * data vector (typically the solution vector). Note that the
 * DataPostprocessor object (i.e., in reality the object of your derived
 * class) has to live until the DataOut object is destroyed as the latter
 * keeps a pointer to the former and will complain if the object pointed to is
 * destroyed while the latter still has a pointer to it. If both the data
 * postprocessor and DataOut objects are local variables of a function (as
 * they are, for example, in step-29), then you can avoid this error by
 * declaring the data postprocessor variable before the DataOut variable as
 * objects are destroyed in reverse order of declaration.
 *
 * In order not to perform needless calculations, DataPostprocessor has to
 * provide information which input data is needed for the calculation of the
 * derived quantities, i.e. whether it needs the values, the first derivative
 * and/or the second derivative of the provided data. DataPostprocessor
 * objects which are used in combination with a DataOutFaces object can also
 * ask for the normal vectors at each point. The information which data is
 * needed has to be provided via the UpdateFlags returned by the virtual
 * function get_needed_update_flags(). It is your responsibility to use only
 * those values which were updated in the calculation of derived quantities.
 * The DataOut object will provide references to the requested data in the
 * call to evaluate_scalar_field() or
 * evaluate_vector_field() (DataOut decides which of the two
 * functions to call depending on whether the finite element in use has only a
 * single, or multiple vector components; note that this is only determined by
 * the number of components in the finite element in use, and not by whether
 * the data computed by a class derived from the current one is scalar or
 * vector valued).
 *
 * Furthermore, derived classes have to implement the get_names() function,
 * where the number of output variables returned by the latter function has to
 * match the size of the vector returned by the former. Furthermore, this
 * number has to match the number of computed quantities, of course.
 *
 *
 * <h3>Use in simpler cases</h3>
 *
 * Deriving from the current class allows to implement very general
 * postprocessors. For example, in the step-32 program, we implement a
 * postprocessor that takes a solution that consists of velocity, pressure and
 * temperature (dim+2 components) and computes a variety of output quantities,
 * some of which are vector valued and some of which are scalar. On the other
 * hand, in step-29 we implement a postprocessor that only computes the
 * magnitude of a complex number given by a two-component finite element. It
 * seems silly to have to implement four virtual functions for this
 * (evaluate_scalar_field() or evaluate_vector_field(), get_names(), get_update_flags() and
 * get_data_component_interpretation()).
 *
 * To this end there are two classes DataPostprocessorScalar and
 * DataPostprocessorVector that are meant to be used if the output quantity is
 * either a single scalar or a single vector (here used meaning to have
 * exactly dim components). When using these classes, one only has to write a
 * constructor that passes the name of the output variable and the update
 * flags to the constructor of the base class and overload the function that
 * actually computes the results.
 *
 * @ingroup output
 * @author Tobias Leicht, 2007, Wolfgang Bangerth, 2016
 */
template <int dim>
class DataPostprocessor: public Subscriptor
{
public:

  /**
   * Destructor. This function doesn't actually do anything but is marked as
   * virtual to ensure that data postprocessors can be destroyed through
   * pointers to the base class.
   */
  virtual ~DataPostprocessor ();

  /**
   * This is the main function which actually performs the postprocessing. The
   * second argument is a reference to the postprocessed data which already has
   * correct size and must be filled by this function.
   *
   * The function takes the values, gradients, and higher derivatives of the
   * solution at all evaluation points, as well as other data such as the
   * cell, via the first argument. Not all of the member vectors of this
   * argument will be filled with data -- in fact, derivatives and other
   * quantities will only be contain valid data if the corresponding flags
   * are returned by by an overloaded version of the get_needed_update_flags()
   * function (implemented in a user's derived class).
   * Otherwise those vectors will be in an unspecified state.
   *
   * This function is called when the finite element field that is being
   * converted into graphical data by DataOut or similar classes represents scalar
   * data, i.e. the finite element in use has only a single vector component.
   */
  virtual
  void
  evaluate_scalar_field (const DataPostprocessorInputs::Scalar<dim> &input_data,
                         std::vector<Vector<double> >               &computed_quantities) const;

  /**
   * @deprecated This function is deprecated. It has been superseded by
   * the evaluate_scalar_field() function that receives a superset of the
   * information provided to the current function through the members
   * of the structure it receives as the first argument.
   *
   * If a user class derived from the current class (or from
   * DataPostprocessorScalar) does not overload the function above,
   * but instead overloads the current (legacy) form of the function,
   * then the default implementation of the function above will simply
   * call the current function. However, not all elements of the
   * DataPostprocessorInputs::Scalar argument the function above
   * receives have corresponding function arguments in the current
   * function, and consequently not all information that function has
   * available is passed on to the current one. In other words, there
   * are pieces of information you may need in an implementation of a
   * postprocess that are available if you overload the new form of this
   * function above, but that are not available if you overload the old
   * form of the function here.
   */
  virtual
  void
  compute_derived_quantities_scalar (const std::vector<double>         &solution_values,
                                     const std::vector<Tensor<1,dim> > &solution_gradients,
                                     const std::vector<Tensor<2,dim> > &solution_hessians,
                                     const std::vector<Point<dim> >    &normals,
                                     const std::vector<Point<dim> >    &evaluation_points,
                                     std::vector<Vector<double> >      &computed_quantities) const DEAL_II_DEPRECATED;

  /**
   * Same as the evaluate_scalar_field() function, but this
   * function is called when the original data vector represents vector data,
   * i.e. the finite element in use has multiple vector components.
   */
  virtual
  void
  evaluate_vector_field (const DataPostprocessorInputs::Vector<dim> &input_data,
                         std::vector<Vector<double> >               &computed_quantities) const;

  /**
   * @deprecated This function is deprecated. It has been superseded by
   * the evaluate_vector_field() function that receives a superset of the
   * information provided to the current function through the members
   * of the structure it receives as the first argument.
   *
   * If a user class derived from the current class (or from
   * DataPostprocessorVector) does not overload the function above,
   * but instead overloads the current (legacy) form of the function,
   * then the default implementation of the function above will simply
   * call the current function. However, not all elements of the
   * DataPostprocessorInputs::Vector argument the function above
   * receives have corresponding function arguments in the current
   * function, and consequently not all information that function has
   * available is passed on to the current one. In other words, there
   * are pieces of information you may need in an implementation of a
   * postprocess that are available if you overload the new form of this
   * function above, but that are not available if you overload the old
   * form of the function here.
   */
  virtual
  void
  compute_derived_quantities_vector (const std::vector<Vector<double> >              &solution_values,
                                     const std::vector<std::vector<Tensor<1,dim> > > &solution_gradients,
                                     const std::vector<std::vector<Tensor<2,dim> > > &solution_hessians,
                                     const std::vector<Point<dim> >                  &normals,
                                     const std::vector<Point<dim> >                  &evaluation_points,
                                     std::vector<Vector<double> >                    &computed_quantities) const DEAL_II_DEPRECATED;

  /**
   * Return the vector of strings describing the names of the computed
   * quantities.
   */
  virtual std::vector<std::string> get_names () const = 0;

  /**
   * This functions returns information about how the individual components of
   * output files that consist of more than one data set are to be
   * interpreted.
   *
   * For example, if one has a finite element for the Stokes equations in 2d,
   * representing components (u,v,p), one would like to indicate that the
   * first two, u and v, represent a logical vector so that later on when we
   * generate graphical output we can hand them off to a visualization program
   * that will automatically know to render them as a vector field, rather
   * than as two separate and independent scalar fields.
   *
   * The default implementation of this function returns a vector of values
   * DataComponentInterpretation::component_is_scalar, indicating that all
   * output components are independent scalar fields. However, if a derived
   * class produces data that represents vectors, it may return a vector that
   * contains values DataComponentInterpretation::component_is_part_of_vector.
   * In the example above, one would return a vector with components
   * (DataComponentInterpretation::component_is_part_of_vector,
   * DataComponentInterpretation::component_is_part_of_vector,
   * DataComponentInterpretation::component_is_scalar) for (u,v,p).
   */
  virtual
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  get_data_component_interpretation () const;

  /**
   * Return, which data has to be provided to compute the derived quantities.
   * This has to be a combination of @p update_values, @p update_gradients and
   * @p update_hessians. If the DataPostprocessor is to be used in combination
   * with DataOutFaces, you may also ask for a update of normals via the @p
   * update_normal_vectors flag.
   */
  virtual UpdateFlags get_needed_update_flags () const = 0;
};



/**
 * This class provides a simpler interface to the functionality offered by the
 * DataPostprocessor class in case one wants to compute only a single scalar
 * quantity from the finite element field passed to the DataOut class. For
 * this particular case, it is clear what the returned value of
 * DataPostprocessor::get_data_component_interpretation() should be and we
 * pass the values returned by get_names() and get_needed_update_flags() to
 * the constructor so that derived classes do not have to implement these
 * functions by hand.
 *
 * All derived classes have to do is implement a constructor and overload
 * either DataPostprocessor::evaluate_scalar_field() or
 * DataPostprocessor::evaluate_vector_field().
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
   * Constructor. Take the name of the single scalar variable computed by
   * classes derived from the current one, as well as the update flags
   * necessary to compute this quantity.
   *
   * @param name The name by which the scalar variable computed by this class
   * should be made available in graphical output files.
   * @param update_flags This has to be a combination of @p update_values, @p
   * update_gradients and @p update_hessians. If the DataPostprocessor is to
   * be used in combination with DataOutFaces, you may also ask for a update
   * of normals via the @p update_normal_vectors flag.
   */
  DataPostprocessorScalar (const std::string &name,
                           const UpdateFlags  update_flags);

  /**
   * Return the vector of strings describing the names of the computed
   * quantities. Given the purpose of this class, this is a vector with a
   * single entry equal to the name given to the constructor.
   */
  virtual std::vector<std::string> get_names () const;

  /**
   * This functions returns information about how the individual components of
   * output files that consist of more than one data set are to be
   * interpreted. Since the current class is meant to be used for a single
   * scalar result variable, the returned value is obviously
   * DataComponentInterpretation::component_is_scalar.
   */
  virtual
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  get_data_component_interpretation () const;

  /**
   * Return, which data has to be provided to compute the derived quantities.
   * The flags returned here are the ones passed to the constructor of this
   * class.
   */
  virtual UpdateFlags get_needed_update_flags () const;

private:
  /**
   * Copies of the two arguments given to the constructor of this class.
   */
  const std::string name;
  const UpdateFlags update_flags;
};



/**
 * This class provides a simpler interface to the functionality offered by the
 * DataPostprocessor class in case one wants to compute only a single vector
 * quantity (defined as having exactly dim components) from the finite element
 * field passed to the DataOut class. For this particular case, it is clear
 * what the returned value of
 * DataPostprocessor::get_data_component_interpretation() should be and we
 * pass the values returned by get_names() and get_needed_update_flags() to
 * the constructor so that derived classes do not have to implement these
 * functions by hand.
 *
 * All derived classes have to do is implement a constructor and overload
 * either DataPostprocessor::evaluate_scalar_field() or
 * DataPostprocessor::evaluate_vector_field().
 *
 * An example of how the closely related class DataPostprocessorScalar is used
 * can be found in step-29.
 *
 * @ingroup output
 * @author Wolfgang Bangerth, 2011
 */
template <int dim>
class DataPostprocessorVector : public DataPostprocessor<dim>
{
public:
  /**
   * Constructor. Take the name of the single vector variable computed by
   * classes derived from the current one, as well as the update flags
   * necessary to compute this quantity.
   *
   * @param name The name by which the vector variable computed by this class
   * should be made available in graphical output files.
   * @param update_flags This has to be a combination of @p update_values, @p
   * update_gradients and @p update_hessians. If the DataPostprocessor is to
   * be used in combination with DataOutFaces, you may also ask for a update
   * of normals via the @p update_normal_vectors flag.
   */
  DataPostprocessorVector (const std::string &name,
                           const UpdateFlags  update_flags);

  /**
   * Return the vector of strings describing the names of the computed
   * quantities. Given the purpose of this class, this is a vector with dim
   * entries all equal to the name given to the constructor.
   */
  virtual std::vector<std::string> get_names () const;

  /**
   * This functions returns information about how the individual components of
   * output files that consist of more than one data set are to be
   * interpreted. Since the current class is meant to be used for a single
   * vector result variable, the returned value is obviously
   * DataComponentInterpretation::component_is_part repeated dim times.
   */
  virtual
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  get_data_component_interpretation () const;

  /**
   * Return which data has to be provided to compute the derived quantities.
   * The flags returned here are the ones passed to the constructor of this
   * class.
   */
  virtual UpdateFlags get_needed_update_flags () const;

private:
  /**
   * Copies of the two arguments given to the constructor of this class.
   */
  const std::string name;
  const UpdateFlags update_flags;
};


DEAL_II_NAMESPACE_CLOSE

#endif

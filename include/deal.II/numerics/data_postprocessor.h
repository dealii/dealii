// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2018 by the deal.II authors
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

#ifndef dealii_data_postprocessor_h
#define dealii_data_postprocessor_h



#include <deal.II/base/subscriptor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/numerics/data_component_interpretation.h>

#include <boost/any.hpp>

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
   * A base class containing common elements for the Scalar and Vector classes
   * that are passed as arguments to
   * DataPostprocessor::evaluate_scalar_field() and
   * DataPostprocessor::evaluate_vector_field(). This common base class
   * provides access to the points at which the solution is being evaluated,
   * and a few other fields, as described in the following.
   *
   * <h4>Normal vector access</h4>
   *
   * If appropriate, i.e., if the object that is currently being processed
   * is a face of a cell and the DataPostprocessor object is called from
   * DataOutFaces or a similar class, then the current object also
   * stores the normal vectors to the geometry on which output
   * is generated, at these evaluation points.
   *
   * On the other hand, if the solution is being evaluated on a cell,
   * then the @p normal_vectors member variable does not contain anything
   * useful.
   *
   * <h4>Cell access</h4>
   *
   * DataPostprocessor is typically called from classes such as DataOut
   * or DataOutFaces that evaluate solution fields on a cell-by-cell
   * basis. As a consequence, classes derived from DataPostprocessor
   * (or DataPostprocessorScalar, DataPostprocessorVector, or
   * DataPostprocessorTensor) sometimes
   * need to use which cell is currently under investigation. Consequently,
   * DataOut and similar classes pass the cell they are currently working
   * on to DataPostprocessor via the classes in this namespace (and
   * specifically the current base class).
   *
   * However, the situation is not so simple. This is because the current
   * class (and those derived from it) only knows the space dimension in
   * which the output lives. But this can come from many sources. First,
   * the cell may be a cell in a DoFHandler or hp::DoFHandler object.
   * Second, if we are in 3d, this may be because we are working on a
   * DoFHandler<3>, or a DoFHandler<2,3> (i.e., either a 3d mesh, or a
   * 2d meshes of a 2d surface embedded in 3d space). Finally, if one
   * considers classes such as DataOutRotation or DataOutStack, then
   * @p spacedim being equal to 3 might mean that we are actually
   * working on a DoFHandler<2> or hp::DoFHandler<2>.
   *
   * In other words, just because we know the value of the @p spacedim
   * template argument of the current class does not mean that the
   * data type of the cell iterator that is currently being worked on
   * is obvious.
   *
   * To make the cell iterator accessible nevertheless, this class uses
   * an object of type boost::any to store the cell iterator. You can
   * think of this as being a void pointer that can point to anything.
   * To use what is being used therefore requires the user to know the
   * data type of the thing being pointed to.
   *
   * To make this work, the DataOut and related classes store in objects
   * of the current type a representation of the cell. To get it back out,
   * you would use the get_cell() function that requires you to say,
   * as a template parameter, the DoFHandler type to which the cell that
   * is currently being processed belongs. This is knowledge you typically
   * have in an application: for example, if your application runs in
   * @p dim space dimensions, uses a hp::DoFHandler, and you are currently
   * using the DataOut class, then the cells that are worked on have data
   * type <code>DataOut<dim>::cell_iterator</code>. Consequently, in a
   * postprocessor, you can call
   * <code>inputs.get_cell@<hp::DoFHandler@<dim@> @> </code>. For technical
   * reasons, however, C++ will typically require you to write this as
   * <code>inputs.template get_cell@<DoFHandler@<dim@> @> </code>
   * because the member function we call here requires that we explicitly
   * provide the template argument.
   *
   * Let us consider a complete example of a postprocessor that computes
   * the fluid norm of the stress $\|\sigma\| = \|\eta \nabla u\|$ from the
   * viscosity $\eta$ and the gradient of the fluid velocity, $\nabla u$,
   * assuming that the viscosity is something that depends on the cell's material
   * id. This can be done using a class we derive from DataPostprocessorScalar
   * where we overload the DataPostprocessor::evaluate_vector_field() function
   * that receives the values and gradients of the velocity (plus of
   * other solution variables such as the pressure, but let's ignore those
   * for the moment). Then we could use code such as this, assuming that we
   * use a hp::DoFHandler:
   * @code
   *   template <int dim>
   *   class ComputeStress : public DataPostprocessorScalar<dim>
   *   {
   *     public:
   *       ... // overload other necessary member variables
   *       virtual
   *       void
   *       evaluate_vector_field (const DataPostprocessorInputs::Vector<dim> &input_data,
   *                              std::vector<Vector<double> >               &computed_quantities) const
   *       {
   *         const typename hp::DoFHandler<dim>::cell_iterator
   *           current_cell = input_data.template get_cell<hp::DoFHandler<dim> >();
   *         const viscosity = look_up_viscosity (current_cell->material_id());
   *
   *         for (unsigned int q=0; q<input_data.solution_gradients.size(); ++q)
   *           computed_quantities[q][0] = (viscosity * input_data.solution_gradients[q]).norm();
   *       }
   *   };
   * @endcode
   *
   * @author Wolfgang Bangerth, 2016
   */
  template <int spacedim>
  struct CommonInputs
  {
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
     * derived from DataPostprocessorScalar, DataPostprocessorVector,
     * or DataPostprocessorTensor may pass this flag to the constructor of
     * these three classes.
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
     * derived from DataPostprocessorScalar, DataPostprocessorVector,
     * or DataPostprocessorTensor may pass this flag to the constructor of
     * these three classes.
     */
    std::vector<Point<spacedim> >     evaluation_points;

    /**
     * Set the cell that is currently being used in evaluating the data
     * for which the DataPostprocessor object is being called.
     *
     * This function is not usually called from user space, but is instead
     * called by DataOut and similar classes when creating the object that
     * is then passed to DataPostprocessor.
     */
    template <typename DoFHandlerType>
    void
    set_cell (const typename DoFHandlerType::cell_iterator &cell);

    /**
     * Query the cell on which we currently produce graphical output.
     * See the documentation of the current class for an example on how
     * to use this function.
     */
    template <typename DoFHandlerType>
    typename DoFHandlerType::cell_iterator
    get_cell () const;

  private:
    /**
     * The place where set_cell() stores the cell. Since the actual data
     * type of the cell iterator can be many different things, the
     * interface uses boost::any here. This makes assignment in set_cell()
     * simple, but requires knowing the data type of the stored object in
     * get_cell().
     */
    boost::any cell;
  };

  /**
   * A structure that is used to pass information to
   * DataPostprocessor::evaluate_scalar_field(). It contains
   * the values and (if requested) derivatives of a scalar solution
   * variable at the evaluation points on a cell or face.
   *
   * Through the fields in the CommonInputs base class, this class also
   * makes available access to the locations of evaluations points,
   * normal vectors (if appropriate), and which cell data is currently
   * being evaluated on (also if appropriate).
   *
   * @author Wolfgang Bangerth, 2016
   */
  template <int spacedim>
  struct Scalar : public CommonInputs<spacedim>
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
     * UpdateFlags::update_gradients.  Alternatively, a class
     * derived from DataPostprocessorScalar, DataPostprocessorVector,
     * or DataPostprocessorTensor may pass this flag to the constructor of
     * these three classes.
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
     * UpdateFlags::update_hessians.  Alternatively, a class
     * derived from DataPostprocessorScalar, DataPostprocessorVector,
     * or DataPostprocessorTensor may pass this flag to the constructor of
     * these three classes.
     */
    std::vector<Tensor<2, spacedim> > solution_hessians;
  };



  /**
   * A structure that is used to pass information to
   * DataPostprocessor::evaluate_vector_field(). It contains
   * the values and (if requested) derivatives of a vector-valued solution
   * variable at the evaluation points on a cell or face.
   *
   * Through the fields in the CommonInputs base class, this class also
   * makes available access to the locations of evaluations points,
   * normal vectors (if appropriate), and which cell data is currently
   * being evaluated on (also if appropriate).
   *
   * @author Wolfgang Bangerth, 2016
   */
  template <int spacedim>
  struct Vector : public CommonInputs<spacedim>
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
     * UpdateFlags::update_gradients.  Alternatively, a class
     * derived from DataPostprocessorScalar, DataPostprocessorVector,
     * or DataPostprocessorTensor may pass this flag to the constructor of
     * these three classes.
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
     * UpdateFlags::update_hessians.  Alternatively, a class
     * derived from DataPostprocessorScalar, DataPostprocessorVector,
     * or DataPostprocessorTensor may pass this flag to the constructor of
     * these three classes.
     */
    std::vector<std::vector<Tensor<2, spacedim> > > solution_hessians;
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
 * To this end there are three classes DataPostprocessorScalar,
 * DataPostprocessorVector, and DataPostprocessorTensor that are meant to be
 * used if the output quantity is either a single scalar, a single vector
 * (here used meaning to have exactly @p dim components), or a single
 * tensor (here used meaning to have exactly <code>dim*dim</code> components).
 * When using these classes, one only has to write a
 * constructor that passes the name of the output variable and the update
 * flags to the constructor of the base class and overload the function that
 * actually computes the results.
 *
 * The DataPostprocessorVector and DataPostprocessorTensor class documentations
 * also contains a extensive examples of how they can be used. The step-29
 * tutorial program contains an example of using the DataPostprocessorScalar
 * class.
 *
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
  virtual ~DataPostprocessor () = default;

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
   * are returned by an overridden version of the get_needed_update_flags()
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
   * Same as the evaluate_scalar_field() function, but this
   * function is called when the original data vector represents vector data,
   * i.e. the finite element in use has multiple vector components.
   */
  virtual
  void
  evaluate_vector_field (const DataPostprocessorInputs::Vector<dim> &input_data,
                         std::vector<Vector<double> >               &computed_quantities) const;

  /**
   * Return the vector of strings describing the names of the computed
   * quantities.
   */
  virtual std::vector<std::string> get_names () const = 0;

  /**
   * This function returns information about how the individual components of
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
 * An example of how the closely related DataPostprocessorVector
 * class can be used is found in the documentation of that class.
 * The same is true for the DataPostprocessorTensor class.
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
   * This function returns information about how the individual components of
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
 * quantity (defined as having exactly @p dim components) from the finite element
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
 * can be found in step-29. An example of how the DataPostprocessorTensor
 * class can be used is found in the documentation of that class.
 *
 *
 * <h3> An example </h3>
 *
 * A common example of what one wants to do with postprocessors is to visualize
 * not just the value of the solution, but the gradient. Let's, for simplicity,
 * assume that you have only a scalar solution. In fact, because it's readily
 * available, let us simply take the step-6 solver to produce such a scalar
 * solution. The gradient is a vector (with exactly @p dim components), so the
 * current class fits the bill to produce the gradient through postprocessing.
 * Then, the following code snippet implements everything you need to have
 * to visualize the gradient:
 * @code
 * template <int dim>
 * class GradientPostprocessor : public DataPostprocessorVector<dim>
 * {
 * public:
 *   GradientPostprocessor ()
 *     :
 *     // call the constructor of the base class. call the variable to
 *     // be output "grad_u" and make sure that DataOut provides us
 *     // with the gradients:
 *     DataPostprocessorVector<dim> ("grad_u",
 *                                   update_gradients)
 *   {}
 *
 *   virtual
 *   void
 *   evaluate_scalar_field (const DataPostprocessorInputs::Scalar<dim> &input_data,
 *                          std::vector<Vector<double> >               &computed_quantities) const
 *   {
 *     // ensure that there really are as many output slots
 *     // as there are points at which DataOut provides the
 *     // gradients:
 *     AssertDimension (input_data.solution_gradients.size(),
 *                      computed_quantities.size());
 *
 *     // then loop over all of these inputs:
 *     for (unsigned int p=0; p<input_data.solution_gradients.size(); ++p)
 *       {
 *         // ensure that each output slot has exactly 'dim'
 *         // components (as should be expected, given that we
 *         // want to create vector-valued outputs), and copy the
 *         // gradients of the solution at the evaluation points
 *         // into the output slots:
 *         AssertDimension (computed_quantities[p].size(), dim);
 *         for (unsigned int d=0; d<dim; ++d)
 *           computed_quantities[p][d]
 *             = input_data.solution_gradients[p][d];
 *       }
 *   }
 * };
 * @endcode
 * The only thing that is necessary is to add another output to the call
 * of DataOut::add_vector() in the @p run() function of the @p Step6 class
 * of that example program. The corresponding code snippet would then look
 * like this (where we also use VTU as the file format to output the data):
 * @code
 *   GradientPostprocessor<dim> gradient_postprocessor;
 *
 *   DataOut<dim> data_out;
 *   data_out.attach_dof_handler (dof_handler);
 *   data_out.add_data_vector (solution, "solution");
 *   data_out.add_data_vector (solution, gradient_postprocessor);
 *   data_out.build_patches ();
 *
 *   std::ofstream output ("solution.vtu");
 *   data_out.write_vtu (output);
 * @endcode
 *
 * This leads to the following output for the solution and the
 * gradients (you may want to compare with the solution shown in the results
 * section of step-6; the current data is generated on a coarser mesh for
 * simplicity):
 *
 * @image html data_postprocessor_vector_0.png
 * @image html data_postprocessor_vector_1.png
 *
 * In the second image, the background color corresponds to the magnitude of the
 * gradient vector and the vector glyphs to the gradient itself. It may be surprising
 * at first to see that from each vertex, multiple vectors originate, going in
 * different directions. But that is because the solution is only continuous:
 * in general, the gradient is discontinuous across edges, and so the multiple
 * vectors originating from each vertex simply represent the differing
 * gradients of the solution at each adjacent cell.
 *
 * The output above -- namely, the gradient $\nabla u$ of the solution --
 * corresponds to the temperature gradient if one interpreted step-6 as solving
 * a steady-state heat transfer problem.
 * It is very small in the central part of the domain because in step-6 we are
 * solving an equation that has a coefficient $a(\mathbf x)$ that is large in
 * the central part and small on the outside. This can be thought as a material
 * that conducts heat well, and consequently the temperature gradient is small.
 * On the other hand, the "heat flux" corresponds to the quantity
 * $a(\mathbf x) \nabla u(\mathbf x)$. For the
 * solution of that equation, the flux should be continuous across the interface.
 * This is easily verified by the following modification of the postprocessor:
 * @code
 * template <int dim>
 * class HeatFluxPostprocessor : public DataPostprocessorVector<dim>
 * {
 * public:
 *   HeatFluxPostprocessor ()
 *     :
 *     // like above, but now also make sure that DataOut provides
 *     // us with coordinates of the evaluation points:
 *     DataPostprocessorVector<dim> ("heat_flux",
 *                                   update_gradients | update_quadrature_points)
 *   {}
 *
 *   virtual
 *   void
 *   evaluate_scalar_field (const DataPostprocessorInputs::Scalar<dim> &input_data,
 *                          std::vector<Vector<double> >               &computed_quantities) const
 *   {
 *     AssertDimension (input_data.solution_gradients.size(),
 *                      computed_quantities.size());
 *
 *     for (unsigned int p=0; p<input_data.solution_gradients.size(); ++p)
 *       {
 *         AssertDimension (computed_quantities[p].size(), dim);
 *         for (unsigned int d=0; d<dim; ++d)
 *           // like above, but also multiply the gradients with
 *           // the coefficient evaluated at the current point:
 *           computed_quantities[p][d]
 *             = coefficient (input_data.evaluation_points[p]) *
 *               input_data.solution_gradients[p][d];
 *       }
 *   }
 * };
 * @endcode
 * With this postprocessor, we get the following picture of the heat flux:
 *
 * @image html data_postprocessor_vector_2.png
 *
 * As the background color shows, the gradient times the coefficient is now
 * a continuous function. There are (large) vectors around the interface
 * where the coefficient jumps (at half the distance between the center of the
 * disk to the perimeter) that seem to point in the wrong direction; this is
 * an artifact of the fact that the solution has a discontinuous gradient
 * at these points and that the numerical solution on the current grid
 * does not adequately resolve this interface. This, however, is not
 * important to the current discussion.
 *
 *
 * <h3> Extension to the gradients of vector-valued problems </h3>
 *
 * The example above uses a scalar solution and its gradient as an example.
 * On the other hand, one may want to do something similar for the gradient
 * of a vector-valued displacement field (such as the strain or
 * stress of a displacement field, like those computed in step-8, step-17,
 * step-18, or step-44). In that case, the solution is already vector
 * valued and the stress is a (symmetric) tensor.
 *
 * deal.II does not currently support outputting tensor-valued quantities, but
 * they can of course be output as a collection of scalar-valued components of
 * the tensor. This can be facilitated using the DataPostprocessorTensor
 * class. The documentation of that class contains an example.
 *
 *
 * @ingroup output
 * @author Wolfgang Bangerth, 2011, 2017
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
   * This function returns information about how the individual components of
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



/**
 * This class provides a simpler interface to the functionality offered by the
 * DataPostprocessor class in case one wants to compute only a single tensor
 * quantity (defined as having exactly <code>dim*dim</code> components) from
 * the finite element field passed to the DataOut class.
 *
 * For this case, we would like to output all of these components as parts
 * of a tensor-valued quantity. Unfortunately, the various backends that
 * write DataOut data in graphical file formats (see the DataOutBase
 * namespace for what formats can be written) do not support tensor data
 * at the current time. In fact, neither does the DataComponentInterpretation
 * namespace that provides semantic information how individual components
 * of graphical data should be interpreted. Nevertheless, like
 * DataPostprocessorScalar and DataPostprocessorVector, this class helps
 * with setting up what the get_names() and get_needed_update_flags()
 * functions required by the DataPostprocessor base class should return,
 * and so the current class implements these based on information that
 * the constructor of the current class receives from further derived
 * classes.
 *
 * (In order to visualize this collection of scalar fields that, together,
 * are then supposed to be interpreted as a tensor, one has to (i) use a
 * visualization program that can visualize tensors, and (ii) teach it
 * how to re-combine the scalar fields into tensors. In the case of
 * Visit -- see https://wci.llnl.gov/simulation/computer-codes/visit/ --
 * this is done by creating a new "Expression": in essence, one creates
 * a variable, say "grad_u", that is tensor-valued and whose value is
 * given by the expression <code>{{grad_u_xx,grad_u_xy},
 * {grad_u_yx, grad_u_yy}}</code>, where the referenced variables are
 * the names of scalar fields that, here, are produced by the example
 * below. Visit is then able to visualize this "new" variable as a
 * tensor.)
 *
 * All derived classes have to do is implement a constructor and overload
 * either DataPostprocessor::evaluate_scalar_field() or
 * DataPostprocessor::evaluate_vector_field().
 *
 * An example of how the closely related class DataPostprocessorScalar is used
 * can be found in step-29. An example of how the DataPostprocessorVector
 * class can be used is found in the documentation of that class.
 *
 *
 * <h3> An example </h3>
 *
 * A common example of what one wants to do with postprocessors is to visualize
 * not just the value of the solution, but the gradient. This class is meant for
 * tensor-valued outputs, so we will start with a vector-valued solution: the
 * displacement field of step-8. The gradient is a rank-2 tensor (with exactly
 * <code>dim*dim</code> components), so the
 * current class fits the bill to produce the gradient through postprocessing.
 * Then, the following code snippet implements everything you need to have
 * to visualize the gradient:
 * @code
 *   template <int dim>
 *   class GradientPostprocessor : public DataPostprocessorTensor<dim>
 *   {
 *   public:
 *     GradientPostprocessor ()
 *       :
 *       DataPostprocessorTensor<dim> ("grad_u",
 *                                     update_gradients)
 *     {}
 *
 *     virtual
 *     void
 *     evaluate_vector_field (const DataPostprocessorInputs::Vector<dim> &input_data,
 *                            std::vector<Vector<double> >               &computed_quantities) const
 *     {
 *       // ensure that there really are as many output slots
 *       // as there are points at which DataOut provides the
 *       // gradients:
 *       AssertDimension (input_data.solution_gradients.size(),
 *                        computed_quantities.size());
 *
 *       for (unsigned int p=0; p<input_data.solution_gradients.size(); ++p)
 *         {
 *           // ensure that each output slot has exactly 'dim*dim'
 *           // components (as should be expected, given that we
 *           // want to create tensor-valued outputs), and copy the
 *           // gradients of the solution at the evaluation points
 *           // into the output slots:
 *           AssertDimension (computed_quantities[p].size(),
 *                            (Tensor<2,dim>::n_independent_components));
 *           for (unsigned int d=0; d<dim; ++d)
 *             for (unsigned int e=0; e<dim; ++e)
 *               computed_quantities[p][Tensor<2,dim>::component_to_unrolled_index(TableIndices<2>(d,e))]
 *                 = input_data.solution_gradients[p][d][e];
 *         }
 *     }
 *   };
 * @endcode
 * The only tricky part in this piece of code is how to sort the <code>dim*dim</code>
 * elements of the strain tensor into the one vector of computed output
 * quantities -- in other words, how to <i>unroll</i> the elements of
 * the tensor into the vector. This is facilitated by the
 * Tensor::component_to_unrolled_index() function that takes a
 * pair of indices that specify a particular element of the
 * tensor and returns a vector index that is then used in the code
 * above to fill the @p computed_quantities array.
 *
 * The last thing that is necessary is to add another output to the call
 * of DataOut::add_vector() in the @p output_results() function of the @p Step8
 * class of that example program. The corresponding code snippet would then look
 * like this:
 * @code
 *     GradientPostprocessor<dim> grad_u;
 *
 *     DataOut<dim> data_out;
 *     data_out.attach_dof_handler (dof_handler);
 *
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *     data_component_interpretation
 *     (dim, DataComponentInterpretation::component_is_part_of_vector);
 *     data_out.add_data_vector (solution, std::vector<std::string>(dim,"displacement"),
 *                               DataOut<dim>::type_dof_data,
 *                               data_component_interpretation);
 *     data_out.add_data_vector (solution, grad_u);
 *     data_out.build_patches ();
 *     data_out.write_vtk (output);
 * @endcode
 *
 * This leads to the following output for the displacement field (i.e., the
 * solution) and the gradients (you may want to compare with the solution shown
 * in the results section of step-8; the current data is generated on a uniform
 * mesh for simplicity):
 *
 * @image html data_postprocessor_tensor_0.png
 * @image html data_postprocessor_tensor_1.png
 *
 * These pictures show an ellipse representing the gradient tensor at, on average,
 * every tenth mesh point. You may want to read through the documentation of the
 * Visit visualization program (see https://wci.llnl.gov/simulation/computer-codes/visit/)
 * for an interpretation of how exactly tensors are visualizated.
 *
 * In elasticity, one is often interested not in the gradient of the displacement,
 * but in the "strain", i.e., the symmetrized version of the gradient
 * $\varepsilon=\frac 12 (\nabla u + \nabla u^T)$. This is easily facilitated with
 * the following minor modification:
 * @code
 *   template <int dim>
 *   class StrainPostprocessor : public DataPostprocessorTensor<dim>
 *   {
 *   public:
 *     StrainPostprocessor ()
 *       :
 *       DataPostprocessorTensor<dim> ("strain",
 *                                     update_gradients)
 *     {}
 *
 *     virtual
 *     void
 *     evaluate_vector_field (const DataPostprocessorInputs::Vector<dim> &input_data,
 *                            std::vector<Vector<double> >               &computed_quantities) const
 *     {
 *       AssertDimension (input_data.solution_gradients.size(),
 *                        computed_quantities.size());
 *
 *       for (unsigned int p=0; p<input_data.solution_gradients.size(); ++p)
 *         {
 *           AssertDimension (computed_quantities[p].size(),
 *                            (Tensor<2,dim>::n_independent_components));
 *           for (unsigned int d=0; d<dim; ++d)
 *             for (unsigned int e=0; e<dim; ++e)
 *               computed_quantities[p][Tensor<2,dim>::component_to_unrolled_index(TableIndices<2>(d,e))]
 *                 = (input_data.solution_gradients[p][d][e]
 *                    +
 *                    input_data.solution_gradients[p][e][d]) / 2;
 *         }
 *     }
 *   };
 * @endcode
 *
 * Using this class in step-8 leads to the following visualization:
 *
 * @image html data_postprocessor_tensor_2.png
 *
 * Given how easy it is to output the strain, it would also not be very
 * complicated to write a postprocessor that computes the <i>stress</i>
 * in the solution field as the stress is easily computed from the
 * strain by multiplication with either the strain-stress tensor or,
 * in simple cases, the Lam&eacute; constants.
 *
 * @ingroup output
 * @author Wolfgang Bangerth, 2017
 */
template <int dim>
class DataPostprocessorTensor : public DataPostprocessor<dim>
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
  DataPostprocessorTensor (const std::string &name,
                           const UpdateFlags  update_flags);

  /**
   * Return the vector of strings describing the names of the computed
   * quantities. Given the purpose of this class, this is a vector with dim
   * entries all equal to the name given to the constructor.
   */
  virtual std::vector<std::string> get_names () const;

  /**
   * This function returns information about how the individual components of
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



#ifndef DOXYGEN
// -------------------- template functions ----------------------

namespace DataPostprocessorInputs
{

  template <int spacedim>
  template <typename DoFHandlerType>
  void
  CommonInputs<spacedim>::set_cell (const typename DoFHandlerType::cell_iterator &new_cell)
  {
    // see if we had previously already stored a cell that has the same
    // data type; if so, reuse the memory location and avoid calling 'new'
    // inside boost::any
    if (typename DoFHandlerType::cell_iterator * storage_location
        = boost::any_cast<typename DoFHandlerType::cell_iterator>(&cell))
      *storage_location = new_cell;
    else
      // if we had nothing stored before, or if we had stored a different
      // data type, just let boost::any replace things
      cell = new_cell;
  }



  template <int spacedim>
  template <typename DoFHandlerType>
  typename DoFHandlerType::cell_iterator
  CommonInputs<spacedim>::get_cell () const
  {
    Assert(cell.empty() == false,
           ExcMessage ("You are trying to access the cell associated with a "
                       "DataPostprocessorInputs::Scalar object for which no cell has "
                       "been set."));
    Assert(boost::any_cast<typename DoFHandlerType::cell_iterator>(&cell) != nullptr,
           ExcMessage ("You are trying to access the cell associated with a "
                       "DataPostprocessorInputs::Scalar with a DoFHandler type that "
                       "is different from the type with which it has been set. For "
                       "example, if the cell for which output is currently being "
                       "generated belongs to a hp::DoFHandler<2,3> object, then you can "
                       "only call the current function with a template argument "
                       "equal to hp::DoFHandler<2,3>, but not with any other class "
                       "type or dimension template argument."));
    return boost::any_cast<typename DoFHandlerType::cell_iterator>(cell);
  }
}

#endif

DEAL_II_NAMESPACE_CLOSE

#endif

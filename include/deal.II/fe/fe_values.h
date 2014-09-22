// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2014 by the deal.II authors
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

#ifndef __deal2__fe_values_h
#define __deal2__fe_values_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/derivative_form.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/vector_slice.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/table.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/multigrid/mg_dof_handler.h>

#include <algorithm>
#include <memory>

// dummy include in order to have the
// definition of PetscScalar available
// without including other PETSc stuff
#ifdef DEAL_II_WITH_PETSC
#  include <petsc.h>
#endif

DEAL_II_NAMESPACE_OPEN

template <int dim>   class Quadrature;
template <int dim, int spacedim=dim> class FEValuesBase;

template <typename Number> class Vector;
template <typename Number> class BlockVector;


namespace internal
{
  /**
   * A class whose specialization is used to define what type the curl of a
   * vector valued function corresponds to.
   */
  template <int dim>
  struct CurlType;

  /**
   * A class whose specialization is used to define what type the curl of a
   * vector valued function corresponds to.
   *
   * In 1d, the curl is a scalar.
   */
  template <>
  struct CurlType<1>
  {
    typedef Tensor<1,1>     type;
  };

  /**
   * A class whose specialization is used to define what type the curl of a
   * vector valued function corresponds to.
   *
   * In 2d, the curl is a scalar.
   */
  template <>
  struct CurlType<2>
  {
    typedef Tensor<1,1>     type;
  };

  /**
   * A class whose specialization is used to define what type the curl of a
   * vector valued function corresponds to.
   *
   * In 3d, the curl is a vector.
   */
  template <>
  struct CurlType<3>
  {
    typedef Tensor<1,3>     type;
  };
}




/**
 * A namespace for "views" on a FEValues, FEFaceValues, or FESubfaceValues
 * object. A view represents only a certain part of the whole: whereas the
 * FEValues object represents <i>all</i> values, gradients, or second
 * derivatives of all components of a vector-valued element, views restrict
 * the attention to only a single component or a subset of components. You
 * typically get objects of classes defined in this namespace by applying
 * FEValuesExtractors objects to a FEValues, FEFaceValues or FESubfaceValues
 * objects using the square bracket operator.
 *
 * There are classes that present views for single scalar components, vector
 * components consisting of <code>dim</code> elements, and symmetric second
 * order tensor components consisting of <code>(dim*dim + dim)/2</code>
 * elements
 *
 * See the description of the @ref vector_valued module for examples how to
 * use the features of this namespace.
 *
 * @ingroup feaccess vector_valued
 */
namespace FEValuesViews
{
  /**
   * A class representing a view to a single scalar component of a possibly
   * vector-valued finite element. Views are discussed in the @ref
   * vector_valued module.
   *
   * You get an object of this type if you apply a FEValuesExtractors::Scalar
   * to an FEValues, FEFaceValues or FESubfaceValues object.
   *
   * @ingroup feaccess vector_valued
   */
  template <int dim, int spacedim=dim>
  class Scalar
  {
  public:
    /**
     * A typedef for the data type of values of the view this class
     * represents. Since we deal with a single components, the value type is a
     * scalar double.
     */
    typedef double        value_type;

    /**
     * A typedef for the type of gradients of the view this class
     * represents. Here, for a scalar component of the finite element, the
     * gradient is a <code>Tensor@<1,dim@></code>.
     */
    typedef dealii::Tensor<1,spacedim> gradient_type;

    /**
     * A typedef for the type of second derivatives of the view this class
     * represents. Here, for a scalar component of the finite element, the
     * Hessian is a <code>Tensor@<2,dim@></code>.
     */
    typedef dealii::Tensor<2,spacedim> hessian_type;

    /**
     * A structure where for each shape function we pre-compute a bunch of
     * data that will make later accesses much cheaper.
     */
    struct ShapeFunctionData
    {
      /**
       * For each shape function, store whether the selected vector component
       * may be nonzero. For primitive shape functions we know for sure
       * whether a certain scalar component of a given shape function is
       * nonzero, whereas for non-primitive shape functions this may not be
       * entirely clear (e.g. for RT elements it depends on the shape of a
       * cell).
       */
      bool is_nonzero_shape_function_component;

      /**
       * For each shape function, store the row index within the shape_values,
       * shape_gradients, and shape_hessians tables (the column index is the
       * quadrature point index). If the shape function is primitive, then we
       * can get this information from the shape_function_to_row_table of the
       * FEValues object; otherwise, we have to work a bit harder to compute
       * this information.
       */
      unsigned int row_index;
    };

    /**
     * Default constructor. Creates an invalid object.
     */
    Scalar ();

    /**
     * Constructor for an object that represents a single scalar component of
     * a FEValuesBase object (or of one of the classes derived from
     * FEValuesBase).
     */
    Scalar (const FEValuesBase<dim,spacedim> &fe_values_base,
            const unsigned int       component);

    /**
     * Copy operator. This is not a lightweight object so we don't allow
     * copying and generate an exception if this function is called.
     */
    Scalar &operator= (const Scalar<dim,spacedim> &);

    /**
     * Return the value of the vector component selected by this view, for the
     * shape function and quadrature point selected by the arguments.
     *
     * @param shape_function Number of the shape function to be
     * evaluated. Note that this number runs from zero to dofs_per_cell, even
     * in the case of an FEFaceValues or FESubfaceValues object.
     *
     * @param q_point Number of the quadrature point at which function is to
     * be evaluated.
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    value_type
    value (const unsigned int shape_function,
           const unsigned int q_point) const;

    /**
     * Return the gradient (a tensor of rank 1) of the vector component
     * selected by this view, for the shape function and quadrature point
     * selected by the arguments.
     *
     * @note The meaning of the arguments is as documented for the value()
     * function.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    gradient_type
    gradient (const unsigned int shape_function,
              const unsigned int q_point) const;

    /**
     * Return the Hessian (the tensor of rank 2 of all second derivatives) of
     * the vector component selected by this view, for the shape function and
     * quadrature point selected by the arguments.
     *
     * @note The meaning of the arguments is as documented for the value()
     * function.
     *
     * @dealiiRequiresUpdateFlags{update_hessians}
     */
    hessian_type
    hessian (const unsigned int shape_function,
             const unsigned int q_point) const;

    /**
     * Return the values of the selected scalar component of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_values function but it only works on the
     * selected scalar component.
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    template <class InputVector>
    void get_function_values (const InputVector &fe_function,
                              std::vector<value_type> &values) const;

    /**
     * Return the gradients of the selected scalar component of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_gradients function but it only works on the
     * selected scalar component.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <class InputVector>
    void get_function_gradients (const InputVector &fe_function,
                                 std::vector<gradient_type> &gradients) const;

    /**
     * Return the Hessians of the selected scalar component of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_hessians function but it only works on the
     * selected scalar component.
     *
     * @dealiiRequiresUpdateFlags{update_hessians}
     */
    template <class InputVector>
    void get_function_hessians (const InputVector &fe_function,
                                std::vector<hessian_type> &hessians) const;

    /**
     * Return the Laplacians of the selected scalar component of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called. The
     * Laplacians are the trace of the Hessians.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_laplacians function but it only works on the
     * selected scalar component.
     *
     * @dealiiRequiresUpdateFlags{update_hessians}
     */
    template <class InputVector>
    void get_function_laplacians (const InputVector &fe_function,
                                  std::vector<value_type> &laplacians) const;

  private:
    /**
     * A reference to the FEValuesBase object we operate on.
     */
    const FEValuesBase<dim,spacedim> &fe_values;

    /**
     * The single scalar component this view represents of the FEValuesBase
     * object.
     */
    const unsigned int component;

    /**
     * Store the data about shape functions.
     */
    std::vector<ShapeFunctionData> shape_function_data;
  };



  /**
   * A class representing a view to a set of <code>spacedim</code> components
   * forming a vector part of a vector-valued finite element. Views are
   * discussed in the @ref vector_valued module.
   *
   * Note that in the current context, a vector is meant in the sense physics
   * uses it: it has <code>spacedim</code> components that behave in specific
   * ways under coordinate system transformations. Examples include velocity
   * or displacement fields. This is opposed to how mathematics uses the word
   * "vector" (and how we use this word in other contexts in the library, for
   * example in the Vector class), where it really stands for a collection of
   * numbers. An example of this latter use of the word could be the set of
   * concentrations of chemical species in a flame; however, these are really
   * just a collection of scalar variables, since they do not change if the
   * coordinate system is rotated, unlike the components of a velocity vector,
   * and consequently, this class should not be used for this context.
   *
   * This class allows to query the value, gradient and divergence of
   * (components of) shape functions and solutions representing vectors. The
   * gradient of a vector $d_{k}, 0\le k<\text{dim}$ is defined as $S_{ij} =
   * \frac{\partial d_{i}}{\partial x_j}, 0\le i,j<\text{dim}$.
   *
   * You get an object of this type if you apply a FEValuesExtractors::Vector
   * to an FEValues, FEFaceValues or FESubfaceValues object.
   *
   * @ingroup feaccess vector_valued
   */
  template <int dim, int spacedim=dim>
  class Vector
  {
  public:
    /**
     * A typedef for the data type of values of the view this class
     * represents. Since we deal with a set of <code>dim</code> components,
     * the value type is a Tensor<1,spacedim>.
     */
    typedef dealii::Tensor<1,spacedim>          value_type;

    /**
     * A typedef for the type of gradients of the view this class
     * represents. Here, for a set of <code>dim</code> components of the
     * finite element, the gradient is a <code>Tensor@<2,spacedim@></code>.
     *
     * See the general documentation of this class for how exactly the
     * gradient of a vector is defined.
     */
    typedef dealii::Tensor<2,spacedim>          gradient_type;

    /**
     * A typedef for the type of symmetrized gradients of the view this class
     * represents. Here, for a set of <code>dim</code> components of the
     * finite element, the symmetrized gradient is a
     * <code>SymmetricTensor@<2,spacedim@></code>.
     *
     * The symmetric gradient of a vector field $\mathbf v$ is defined as
     * $\varepsilon(\mathbf v)=\frac 12 (\nabla \mathbf v + \nabla \mathbf
     * v^T)$.
     */
    typedef dealii::SymmetricTensor<2,spacedim> symmetric_gradient_type;

    /**
     * A typedef for the type of the divergence of the view this class
     * represents. Here, for a set of <code>dim</code> components of the
     * finite element, the divergence of course is a scalar.
     */
    typedef double                 divergence_type;

    /**
     * A typedef for the type of the curl of the view this class
     * represents. Here, for a set of <code>spacedim=2</code> components of
     * the finite element, the curl is a <code>Tensor@<1, 1@></code>. For
     * <code>spacedim=3</code> it is a <code>Tensor@<1, dim@></code>.
     */
    typedef typename dealii::internal::CurlType<spacedim>::type   curl_type;

    /**
     * A typedef for the type of second derivatives of the view this class
     * represents. Here, for a set of <code>dim</code> components of the
     * finite element, the Hessian is a <code>Tensor@<3,dim@></code>.
     */
    typedef dealii::Tensor<3,spacedim>          hessian_type;

    /**
     * A structure where for each shape function we pre-compute a bunch of
     * data that will make later accesses much cheaper.
     */
    struct ShapeFunctionData
    {
      /**
       * For each pair (shape function,component within vector), store whether
       * the selected vector component may be nonzero. For primitive shape
       * functions we know for sure whether a certain scalar component of a
       * given shape function is nonzero, whereas for non-primitive shape
       * functions this may not be entirely clear (e.g. for RT elements it
       * depends on the shape of a cell).
       */
      bool is_nonzero_shape_function_component[spacedim];

      /**
       * For each pair (shape function, component within vector), store the
       * row index within the shape_values, shape_gradients, and
       * shape_hessians tables (the column index is the quadrature point
       * index). If the shape function is primitive, then we can get this
       * information from the shape_function_to_row_table of the FEValues
       * object; otherwise, we have to work a bit harder to compute this
       * information.
       */
      unsigned int row_index[spacedim];

      /**
       * For each shape function say the following: if only a single entry in
       * is_nonzero_shape_function_component for this shape function is
       * nonzero, then store the corresponding value of row_index and
       * single_nonzero_component_index represents the index between 0 and dim
       * for which it is attained. If multiple components are nonzero, then
       * store -1. If no components are nonzero then store -2.
       */
      int          single_nonzero_component;
      unsigned int single_nonzero_component_index;
    };

    /**
     * Default constructor. Creates an invalid object.
     */
    Vector ();

    /**
     * Constructor for an object that represents dim components of a
     * FEValuesBase object (or of one of the classes derived from
     * FEValuesBase), representing a vector-valued variable.
     *
     * The second argument denotes the index of the first component of the
     * selected vector.
     */
    Vector (const FEValuesBase<dim,spacedim> &fe_values_base,
            const unsigned int first_vector_component);

    /**
     * Copy operator. This is not a lightweight object so we don't allow
     * copying and generate an exception if this function is called.
     */
    Vector &operator= (const Vector<dim,spacedim> &);

    /**
     * Return the value of the vector components selected by this view, for
     * the shape function and quadrature point selected by the
     * arguments. Here, since the view represents a vector-valued part of the
     * FEValues object with <code>dim</code> components, the return type is a
     * tensor of rank 1 with <code>dim</code> components.
     *
     * @param shape_function Number of the shape function to be
     * evaluated. Note that this number runs from zero to dofs_per_cell, even
     * in the case of an FEFaceValues or FESubfaceValues object.
     *
     * @param q_point Number of the quadrature point at which function is to
     * be evaluated.
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    value_type
    value (const unsigned int shape_function,
           const unsigned int q_point) const;

    /**
     * Return the gradient (a tensor of rank 2) of the vector component
     * selected by this view, for the shape function and quadrature point
     * selected by the arguments.
     *
     * See the general documentation of this class for how exactly the
     * gradient of a vector is defined.
     *
     * @note The meaning of the arguments is as documented for the value()
     * function.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    gradient_type
    gradient (const unsigned int shape_function,
              const unsigned int q_point) const;

    /**
     * Return the symmetric gradient (a symmetric tensor of rank 2) of the
     * vector component selected by this view, for the shape function and
     * quadrature point selected by the arguments.
     *
     * The symmetric gradient is defined as $\frac 12 [(\nabla \phi_i(x_q)) +
     * (\nabla \phi_i(x_q))^T]$, where $\phi_i$ represents the
     * <code>dim</code> components selected from the FEValuesBase object, and
     * $x_q$ is the location of the $q$-th quadrature point.
     *
     * @note The meaning of the arguments is as documented for the value()
     * function.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    symmetric_gradient_type
    symmetric_gradient (const unsigned int shape_function,
                        const unsigned int q_point) const;

    /**
     * Return the scalar divergence of the vector components selected by this
     * view, for the shape function and quadrature point selected by the
     * arguments.
     *
     * @note The meaning of the arguments is as documented for the value()
     * function.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    divergence_type
    divergence (const unsigned int shape_function,
                const unsigned int q_point) const;

    /**
     * Return the vector curl of the vector components selected by this view,
     * for the shape function and quadrature point selected by the
     * arguments. For 1d this function does not make any sense. Thus it is not
     * implemented for <code>spacedim=1</code>.  In 2d the curl is defined as
     * @f{equation*} \operatorname{curl}(u):=\frac{du_2}{dx} -\frac{du_1}{dy},
     * @f} whereas in 3d it is given by @f{equation*}
     * \operatorname{curl}(u):=\left( \begin{array}{c}
     * \frac{du_3}{dy}-\frac{du_2}{dz}\\ \frac{du_1}{dz}-\frac{du_3}{dx}\\
     * \frac{du_2}{dx}-\frac{du_1}{dy} \end{array} \right).  @f}
     *
     * @note The meaning of the arguments is as documented for the value()
     * function.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    curl_type
    curl (const unsigned int shape_function,
          const unsigned int q_point) const;

    /**
     * Return the Hessian (the tensor of rank 2 of all second derivatives) of
     * the vector components selected by this view, for the shape function and
     * quadrature point selected by the arguments.
     *
     * @note The meaning of the arguments is as documented for the value()
     * function.
     *
     * @dealiiRequiresUpdateFlags{update_hessians}
     */
    hessian_type
    hessian (const unsigned int shape_function,
             const unsigned int q_point) const;

    /**
     * Return the values of the selected vector components of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_values function but it only works on the
     * selected vector components.
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    template <class InputVector>
    void get_function_values (const InputVector &fe_function,
                              std::vector<value_type> &values) const;

    /**
     * Return the gradients of the selected vector components of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_gradients function but it only works on the
     * selected vector components.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <class InputVector>
    void get_function_gradients (const InputVector &fe_function,
                                 std::vector<gradient_type> &gradients) const;

    /**
     * Return the symmetrized gradients of the selected vector components of
     * the finite element function characterized by <tt>fe_function</tt> at
     * the quadrature points of the cell, face or subface selected the last
     * time the <tt>reinit</tt> function of the FEValues object was called.
     *
     * The symmetric gradient of a vector field $\mathbf v$ is defined as
     * $\varepsilon(\mathbf v)=\frac 12 (\nabla \mathbf v + \nabla \mathbf
     * v^T)$.
     *
     * @note There is no equivalent function such as
     * FEValuesBase::get_function_symmetric_gradients in the FEValues classes
     * but the information can be obtained from
     * FEValuesBase::get_function_gradients, of course.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <class InputVector>
    void
    get_function_symmetric_gradients (const InputVector &fe_function,
                                      std::vector<symmetric_gradient_type> &symmetric_gradients) const;

    /**
     * Return the divergence of the selected vector components of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * There is no equivalent function such as
     * FEValuesBase::get_function_divergences in the FEValues classes but the
     * information can be obtained from FEValuesBase::get_function_gradients,
     * of course.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <class InputVector>
    void get_function_divergences (const InputVector &fe_function,
                                   std::vector<divergence_type> &divergences) const;

    /**
     * Return the curl of the selected vector components of the finite element
     * function characterized by <tt>fe_function</tt> at the quadrature points
     * of the cell, face or subface selected the last time the <tt>reinit</tt>
     * function of the FEValues object was called.
     *
     * There is no equivalent function such as
     * FEValuesBase::get_function_curls in the FEValues classes but the
     * information can be obtained from FEValuesBase::get_function_gradients,
     * of course.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <class InputVector>
    void get_function_curls (const InputVector &fe_function,
                             std::vector<curl_type> &curls) const;

    /**
     * Return the Hessians of the selected vector components of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_hessians function but it only works on the
     * selected vector components.
     *
     * @dealiiRequiresUpdateFlags{update_hessians}
     */
    template <class InputVector>
    void get_function_hessians (const InputVector &fe_function,
                                std::vector<hessian_type> &hessians) const;

    /**
     * Return the Laplacians of the selected vector components of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called. The
     * Laplacians are the trace of the Hessians.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_laplacians function but it only works on the
     * selected vector components.
     *
     * @dealiiRequiresUpdateFlags{update_hessians}
     */
    template <class InputVector>
    void get_function_laplacians (const InputVector &fe_function,
                                  std::vector<value_type> &laplacians) const;

  private:
    /**
     * A reference to the FEValuesBase object we operate on.
     */
    const FEValuesBase<dim,spacedim> &fe_values;

    /**
     * The first component of the vector this view represents of the
     * FEValuesBase object.
     */
    const unsigned int first_vector_component;

    /**
     * Store the data about shape functions.
     */
    std::vector<ShapeFunctionData> shape_function_data;
  };


  template <int rank, int dim, int spacedim = dim>
  class SymmetricTensor;

  /**
   * A class representing a view to a set of <code>(dim*dim + dim)/2</code>
   * components forming a symmetric second-order tensor from a vector-valued
   * finite element. Views are discussed in the @ref vector_valued module.
   *
   * This class allows to query the value and divergence of (components of)
   * shape functions and solutions representing symmetric tensors. The
   * divergence of a symmetric tensor $S_{ij}, 0\le i,j<\text{dim}$ is defined
   * as $d_i = \sum_j \frac{\partial S_{ij}}{\partial x_j}, 0\le
   * i<\text{dim}$, which due to the symmetry of the tensor is also $d_i =
   * \sum_j \frac{\partial S_{ji}}{\partial x_j}$.  In other words, it due to
   * the symmetry of $S$ it does not matter whether we apply the nabla
   * operator by row or by column to get the divergence.
   *
   * You get an object of this type if you apply a
   * FEValuesExtractors::SymmetricTensor to an FEValues, FEFaceValues or
   * FESubfaceValues object.
   *
   * @ingroup feaccess vector_valued
   *
   * @author Andrew McBride, 2009
   */
  template <int dim, int spacedim>
  class SymmetricTensor<2,dim,spacedim>
  {
  public:
    /**
     * A typedef for the data type of values of the view this class
     * represents. Since we deal with a set of <code>(dim*dim + dim)/2</code>
     * components (i.e. the unique components of a symmetric second-order
     * tensor), the value type is a SymmetricTensor<2,spacedim>.
     */
    typedef dealii::SymmetricTensor<2, spacedim> value_type;

    /**
     * A typedef for the type of the divergence of the view this class
     * represents. Here, for a set of of <code>(dim*dim + dim)/2</code> unique
     * components of the finite element representing a symmetric second-order
     * tensor, the divergence of course is a * <code>Tensor@<1,dim@></code>.
     *
     * See the general discussion of this class for a definition of the
     * divergence.
     */
    typedef dealii::Tensor<1, spacedim> divergence_type;

    /**
     * A structure where for each shape function we pre-compute a bunch of
     * data that will make later accesses much cheaper.
     */
    struct ShapeFunctionData
    {
      /**
       * For each pair (shape function,component within vector), store whether
       * the selected vector component may be nonzero. For primitive shape
       * functions we know for sure whether a certain scalar component of a
       * given shape function is nonzero, whereas for non-primitive shape
       * functions this may not be entirely clear (e.g. for RT elements it
       * depends on the shape of a cell).
       */
      bool is_nonzero_shape_function_component[value_type::n_independent_components];

      /**
       * For each pair (shape function, component within vector), store the
       * row index within the shape_values, shape_gradients, and
       * shape_hessians tables (the column index is the quadrature point
       * index). If the shape function is primitive, then we can get this
       * information from the shape_function_to_row_table of the FEValues
       * object; otherwise, we have to work a bit harder to compute this
       * information.
       */
      unsigned int row_index[value_type::n_independent_components];

      /**
       * For each shape function say the following: if only a single entry in
       * is_nonzero_shape_function_component for this shape function is
       * nonzero, then store the corresponding value of row_index and
       * single_nonzero_component_index represents the index between 0 and
       * (dim^2 + dim)/2 for which it is attained. If multiple components are
       * nonzero, then store -1. If no components are nonzero then store -2.
       */
      int single_nonzero_component;
      unsigned int single_nonzero_component_index;
    };

    /**
     * Default constructor. Creates an invalid object.
     */
    SymmetricTensor();

    /**
     * Constructor for an object that represents <code>(dim*dim +
     * dim)/2</code> components of a FEValuesBase object (or of one of the
     * classes derived from FEValuesBase), representing the unique components
     * comprising a symmetric second- order tensor valued variable.
     *
     * The second argument denotes the index of the first component of the
     * selected symmetric second order tensor.
     */
    SymmetricTensor(const FEValuesBase<dim, spacedim> &fe_values_base,
                    const unsigned int first_tensor_component);

    /**
     * Copy operator. This is not a lightweight object so we don't allow
     * copying and generate an exception if this function is called.
     */
    SymmetricTensor &operator=(const SymmetricTensor<2, dim, spacedim> &);

    /**
     * Return the value of the vector components selected by this view, for
     * the shape function and quadrature point selected by the
     * arguments. Here, since the view represents a vector-valued part of the
     * FEValues object with <code>(dim*dim + dim)/2</code> components (the
     * unique components of a symmetric second-order tensor), the return type
     * is a symmetric tensor of rank 2.
     *
     * @param shape_function Number of the shape function to be
     * evaluated. Note that this number runs from zero to dofs_per_cell, even
     * in the case of an FEFaceValues or FESubfaceValues object.
     *
     * @param q_point Number of the quadrature point at which function is to
     * be evaluated.
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    value_type
    value (const unsigned int shape_function,
           const unsigned int q_point) const;


    /**
     * Return the vector divergence of the vector components selected by this
     * view, for the shape function and quadrature point selected by the
     * arguments.
     *
     * See the general discussion of this class for a definition of the
     * divergence.
     *
     * @note The meaning of the arguments is as documented for the value()
     * function.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    divergence_type
    divergence (const unsigned int shape_function,
                const unsigned int q_point) const;

    /**
     * Return the values of the selected vector components of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_values function but it only works on the
     * selected vector components.
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    template <class InputVector>
    void get_function_values (const InputVector &fe_function,
                              std::vector<value_type> &values) const;

    /**
     * Return the divergence of the selected vector components of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * There is no equivalent function such as
     * FEValuesBase::get_function_divergences in the FEValues classes but the
     * information can be obtained from FEValuesBase::get_function_gradients,
     * of course.
     *
     * See the general discussion of this class for a definition of the
     * divergence.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <class InputVector>
    void get_function_divergences (const InputVector &fe_function,
                                   std::vector<divergence_type> &divergences) const;

  private:
    /**
     * A reference to the FEValuesBase object we operate on.
     */
    const FEValuesBase<dim, spacedim> &fe_values;

    /**
     * The first component of the vector this view represents of the
     * FEValuesBase object.
     */
    const unsigned int first_tensor_component;

    /**
     * Store the data about shape functions.
     */
    std::vector<ShapeFunctionData> shape_function_data;
  };


  template <int rank, int dim, int spacedim = dim>
  class Tensor;

  /**
   * A class representing a view to a set of <code>dim*dim</code> components
   * forming a second-order tensor from a vector-valued finite element. Views
   * are discussed in the @ref vector_valued module.
   *
   * This class allows to query the value and divergence of (components of)
   * shape functions and solutions representing tensors. The divergence of a
   * tensor $T_{ij}, 0\le i,j<\text{dim}$ is defined as $d_i = \sum_j
   * \frac{\partial T_{ji}}{\partial x_j}, 0\le i<\text{dim}$.
   *
   * You get an object of this type if you apply a FEValuesExtractors::Tensor
   * to an FEValues, FEFaceValues or FESubfaceValues object.
   *
   * @ingroup feaccess vector_valued
   *
   * @author Denis Davydov, 2013
   */
  template <int dim, int spacedim>
  class Tensor<2,dim,spacedim>
  {
  public:

    /**
     * Data type for what you get when you apply an extractor of this kind to
     * a vector-valued finite element.
     */
    typedef dealii::Tensor<2, spacedim> value_type;

    /**
     * Data type for taking the divergence of a tensor: a vector.
     */
    typedef dealii::Tensor<1, spacedim> divergence_type;

    /**
     * A structure where for each shape function we pre-compute a bunch of
     * data that will make later accesses much cheaper.
     */
    struct ShapeFunctionData
    {
      /**
       * For each pair (shape function,component within vector), store whether
       * the selected vector component may be nonzero. For primitive shape
       * functions we know for sure whether a certain scalar component of a
       * given shape function is nonzero, whereas for non-primitive shape
       * functions this may not be entirely clear (e.g. for RT elements it
       * depends on the shape of a cell).
       */
      bool is_nonzero_shape_function_component[value_type::n_independent_components];

      /**
       * For each pair (shape function, component within vector), store the
       * row index within the shape_values, shape_gradients, and
       * shape_hessians tables (the column index is the quadrature point
       * index). If the shape function is primitive, then we can get this
       * information from the shape_function_to_row_table of the FEValues
       * object; otherwise, we have to work a bit harder to compute this
       * information.
       */
      unsigned int row_index[value_type::n_independent_components];

      /**
       * For each shape function say the following: if only a single entry in
       * is_nonzero_shape_function_component for this shape function is
       * nonzero, then store the corresponding value of row_index and
       * single_nonzero_component_index represents the index between 0 and
       * (dim^2) for which it is attained. If multiple components are nonzero,
       * then store -1. If no components are nonzero then store -2.
       */
      int single_nonzero_component;
      unsigned int single_nonzero_component_index;
    };

    /**
     * Default constructor. Creates an invalid object.
     */
    Tensor();


    /**
     * Constructor for an object that represents <code>(dim*dim)</code>
     * components of a FEValuesBase object (or of one of the classes derived
     * from FEValuesBase), representing the unique components comprising a
     * second-order tensor valued variable.
     *
     * The second argument denotes the index of the first component of the
     * selected symmetric second order tensor.
     */
    Tensor(const FEValuesBase<dim, spacedim> &fe_values_base,
           const unsigned int first_tensor_component);


    /**
     * Copy operator. This is not a lightweight object so we don't allow
     * copying and generate an exception if this function is called.
     */
    Tensor &operator=(const Tensor<2, dim, spacedim> &);

    /**
     * Return the value of the vector components selected by this view, for
     * the shape function and quadrature point selected by the
     * arguments. Here, since the view represents a vector-valued part of the
     * FEValues object with <code>(dim*dim)</code> components (the unique
     * components of a second-order tensor), the return type is a tensor of
     * rank 2.
     *
     * @param shape_function Number of the shape function to be
     * evaluated. Note that this number runs from zero to dofs_per_cell, even
     * in the case of an FEFaceValues or FESubfaceValues object.
     *
     * @param q_point Number of the quadrature point at which function is to
     * be evaluated.
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    value_type
    value (const unsigned int shape_function,
           const unsigned int q_point) const;

    /**
     * Return the vector divergence of the vector components selected by this
     * view, for the shape function and quadrature point selected by the
     * arguments.
     *
     * See the general discussion of this class for a definition of the
     * divergence.
     *
     * @note The meaning of the arguments is as documented for the value()
     * function.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    divergence_type
    divergence (const unsigned int shape_function,
                const unsigned int q_point) const;

    /**
     * Return the values of the selected vector components of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_values function but it only works on the
     * selected vector components.
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    template <class InputVector>
    void get_function_values (const InputVector &fe_function,
                              std::vector<value_type> &values) const;


    /**
     * Return the divergence of the selected vector components of the finite
     * element function characterized by <tt>fe_function</tt> at the
     * quadrature points of the cell, face or subface selected the last time
     * the <tt>reinit</tt> function of the FEValues object was called.
     *
     * There is no equivalent function such as
     * FEValuesBase::get_function_divergences in the FEValues classes but the
     * information can be obtained from FEValuesBase::get_function_gradients,
     * of course.
     *
     * See the general discussion of this class for a definition of the
     * divergence.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <class InputVector>
    void get_function_divergences (const InputVector &fe_function,
                                   std::vector<divergence_type> &divergences) const;

  private:
    /**
     * A reference to the FEValuesBase object we operate on.
     */
    const FEValuesBase<dim, spacedim> &fe_values;

    /**
     * The first component of the vector this view represents of the
     * FEValuesBase object.
     */
    const unsigned int first_tensor_component;

    /**
     * Store the data about shape functions.
     */
    std::vector<ShapeFunctionData> shape_function_data;
  };

}


namespace internal
{
  namespace FEValuesViews
  {
    /**
     * A class objects of which store a collection of FEValuesViews::Scalar,
     * FEValuesViews::Vector, etc object. The FEValuesBase class uses it to
     * generate all possible Views classes upon construction time; we do this
     * at construction time since the Views classes cache some information and
     * are therefore relatively expensive to create.
     */
    template <int dim, int spacedim>
    struct Cache
    {
      /**
       * Caches for scalar and vector, and symmetric second-order tensor
       * valued views.
       */
      std::vector<dealii::FEValuesViews::Scalar<dim,spacedim> > scalars;
      std::vector<dealii::FEValuesViews::Vector<dim,spacedim> > vectors;
      std::vector<dealii::FEValuesViews::SymmetricTensor<2,dim,spacedim> >
      symmetric_second_order_tensors;
      std::vector<dealii::FEValuesViews::Tensor<2,dim,spacedim> >
      second_order_tensors;

      /**
       * Constructor.
       */
      Cache (const FEValuesBase<dim,spacedim> &fe_values);
    };
  }
}



//TODO: Add access to mapping values to FEValuesBase
//TODO: Several FEValuesBase of a system should share Mapping

/**
 * Contains all data vectors for FEValues.  This class has been extracted from
 * FEValuesBase to be handed over to the fill functions of Mapping and
 * FiniteElement.
 *
 * @note All data fields are public, but this is not critical, because access
 * to this object is private in FEValues.
 *
 * The purpose of this class is discussed on the page on @ref
 * UpdateFlagsEssay.
 *
 * @ingroup feaccess
 * @author Guido Kanschat
 * @date 2000
 */
template <int dim, int spacedim=dim>
class FEValuesData
{
public:
  /**
   * Initialize all vectors to correct size.
   */
  void initialize (const unsigned int        n_quadrature_points,
                   const FiniteElement<dim,spacedim> &fe,
                   const UpdateFlags         flags);

  /**
   * Storage type for shape values. Each row in the matrix denotes the values
   * of a single shape function at the different points, columns are for a
   * single point with the different shape functions.
   *
   * If a shape function has more than one non-zero component (in deal.II
   * diction: it is non-primitive), then we allocate one row per non-zero
   * component, and shift subsequent rows backward.  Lookup of the correct row
   * for a shape function is thus simple in case the entire finite element is
   * primitive (i.e. all shape functions are primitive), since then the shape
   * function number equals the row number. Otherwise, use the
   * #shape_function_to_row_table array to get at the first row that belongs
   * to this particular shape function, and navigate among all the rows for
   * this shape function using the FiniteElement::get_nonzero_components()
   * function which tells us which components are non-zero and thus have a row
   * in the array presently under discussion.
   */
  typedef Table<2,double> ShapeVector;

  /**
   * Storage type for gradients. The layout of data is the same as for the
   * #ShapeVector data type.
   */
  typedef std::vector<std::vector<Tensor<1,spacedim> > > GradientVector;

  /**
   * Likewise for second order derivatives.
   */
  typedef std::vector<std::vector<Tensor<2,spacedim> > > HessianVector;

  /**
   * Store the values of the shape functions at the quadrature points. See the
   * description of the data type for the layout of the data in this field.
   */
  ShapeVector shape_values;

  /**
   * Store the gradients of the shape functions at the quadrature points. See
   * the description of the data type for the layout of the data in this
   * field.
   */
  GradientVector shape_gradients;

  /**
   * Store the 2nd derivatives of the shape functions at the quadrature
   * points.  See the description of the data type for the layout of the data
   * in this field.
   */
  HessianVector shape_hessians;

  /**
   * Store an array of weights times the Jacobi determinant at the quadrature
   * points. This function is reset each time reinit() is called. The Jacobi
   * determinant is actually the reciprocal value of the Jacobi matrices
   * stored in this class, see the general documentation of this class for
   * more information.
   *
   * However, if this object refers to an FEFaceValues or FESubfaceValues
   * object, then the JxW_values correspond to the Jacobian of the
   * transformation of the face, not the cell, i.e. the dimensionality is that
   * of a surface measure, not of a volume measure. In this case, it is
   * computed from the boundary forms, rather than the Jacobian matrix.
   */
  std::vector<double>       JxW_values;

  /**
   * Array of the Jacobian matrices at the quadrature points.
   */
  std::vector< DerivativeForm<1,dim,spacedim> > jacobians;

  /**
   * Array of the derivatives of the Jacobian matrices at the quadrature
   * points.
   */
  std::vector<DerivativeForm<2,dim,spacedim> >  jacobian_grads;

  /**
   * Array of the inverse Jacobian matrices at the quadrature points.
   */
  std::vector<DerivativeForm<1,spacedim,dim> > inverse_jacobians;

  /**
   * Array of quadrature points. This array is set up upon calling reinit()
   * and contains the quadrature points on the real element, rather than on
   * the reference element.
   */
  std::vector<Point<spacedim> >  quadrature_points;

  /**
   * List of outward normal vectors at the quadrature points. This field is
   * filled in by the finite element class.
   */
  std::vector<Point<spacedim> >  normal_vectors;

  /**
   * List of boundary forms at the quadrature points. This field is filled in
   * by the finite element class.
   */
  std::vector<Tensor<1,spacedim> >  boundary_forms;

  /**
  * When asked for the value (or gradient, or Hessian) of shape function i's
  * c-th vector component, we need to look it up in the #shape_values,
  * #shape_gradients and #shape_hessians arrays.  The question is where in
  * this array does the data for shape function i, component c reside. This is
  * what this table answers.
  *
  * The format of the table is as
  * follows:
  * - It has dofs_per_cell times
  *   n_components entries.
  * - The entry that corresponds to
  *   shape function i, component c
  *   is <code>i * n_components + c</code>.
  * - The value stored at this
  *   position indicates the row
  *   in #shape_values and the
  *   other tables where the
  *   corresponding datum is stored
  *   for all the quadrature points.
  *
  * In the general, vector-valued context, the number of components is larger
  * than one, but for a given shape function, not all vector components may be
  * nonzero (e.g., if a shape function is primitive, then exactly one vector
  * component is non-zero, while the others are all zero). For such zero
  * components, #shape_values and friends do not have a row. Consequently, for
  * vector components for which shape function i is zero, the entry in the
  * current table is numbers::invalid_unsigned_int.
  *
  * On the other hand, the table is guaranteed to have at least one valid
  * index for each shape function. In particular, for a primitive finite
  * element, each shape function has exactly one nonzero component and so for
  * each i, there is exactly one valid index within the range
  * <code>[i*n_components, (i+1)*n_components)</code>.
   */
  std::vector<unsigned int> shape_function_to_row_table;

  /**
   * Original update flags handed to the constructor of FEValues.
   */
  UpdateFlags          update_flags;
};


/**
 * FEValues, FEFaceValues and FESubfaceValues objects are interfaces to finite
 * element and mapping classes on the one hand side, to cells and quadrature
 * rules on the other side. They allow to evaluate values or derivatives of
 * shape functions at the quadrature points of a quadrature formula when
 * projected by a mapping from the unit cell onto a cell in real space. The
 * reason for this abstraction is possible optimization: Depending on the type
 * of finite element and mapping, some values can be computed once on the unit
 * cell. Others must be computed on each cell, but maybe computation of
 * several values at the same time offers ways for optimization. Since this
 * interlay may be complex and depends on the actual finite element, it cannot
 * be left to the applications programmer.
 *
 * FEValues, FEFaceValues and FESubfaceValues provide only data handling:
 * computations are left to objects of type Mapping and FiniteElement. These
 * provide functions <tt>get_*_data</tt> and <tt>fill_*_values</tt> which are
 * called by the constructor and <tt>reinit</tt> functions of
 * <tt>FEValues*</tt>, respectively.
 *
 * <h3>General usage</h3>
 *
 * Usually, an object of <tt>FEValues*</tt> is used in integration loops over
 * all cells of a triangulation (or faces of cells). To take full advantage of
 * the optimization features, it should be constructed before the loop so that
 * information that does not depend on the location and shape of cells can be
 * computed once and for all (this includes, for example, the values of shape
 * functions at quadrature points for the most common elements: we can
 * evaluate them on the unit cell and they will be the same when mapped to the
 * real cell). Then, in the loop over all cells, it must be re-initialized for
 * each grid cell to compute that part of the information that changes
 * depending on the actual cell (for example, the gradient of shape functions
 * equals the gradient on the unit cell -- which can be computed once and for
 * all -- times the Jacobian matrix of the mapping between unit and real cell,
 * which needs to be recomputed for each cell).
 *
 * A typical piece of code, adding up local contributions to the Laplace
 * matrix looks like this:
 *
 * @code
 * FEValues values (mapping, finite_element, quadrature, flags);
 * for (cell = dof_handler.begin_active();
 *      cell != dof_handler.end();
 *      ++cell)
 *   {
 *     values.reinit(cell);
 *     for (unsigned int q=0; q<quadrature.size(); ++q)
 *       for (unsigned int i=0; i<finite_element.dofs_per_cell; ++i)
 *         for (unsigned int j=0; j<finite_element.dofs_per_cell; ++j)
 *         A(i,j) += fe_values.shape_value(i,q) *
 *                   fe_values.shape_value(j,q) *
 *                   fe_values.JxW(q);
 *     ...
 *   }
 * @endcode
 *
 * The individual functions used here are described below. Note that by
 * design, the order of quadrature points used inside the FEValues object is
 * the same as defined by the quadrature formula passed to the constructor of
 * the FEValues object above.
 *
 *  <h3>Member functions</h3>
 *
 *  The functions of this class fall into different cathegories:
 *  <ul>
 *  <li> shape_value(), shape_grad(), etc: return one of the values
 *    of this object at a time. These functions are inlined, so this
 *    is the suggested access to all finite element values. There
 *    should be no loss in performance with an optimizing compiler. If
 *    the finite element is vector valued, then these functions return
 *    the only non-zero component of the requested shape
 *    function. However, some finite elements have shape functions
 *    that have more than one non-zero component (we call them
 *    non-"primitive"), and in this case this set of functions will
 *    throw an exception since they cannot generate a useful
 *    result. Rather, use the next set of functions.
 *
 *  <li> shape_value_component(), shape_grad_component(), etc:
 *    This is the same set of functions as above, except that for vector
 *    valued finite elements they return only one vector component. This
 *    is useful for elements of which shape functions have more than one
 *    non-zero component, since then the above functions cannot be used,
 *    and you have to walk over all (or only the non-zero) components of
 *    the shape function using this set of functions.
 *
 *  <li> get_function_values(), get_function_gradients(), etc.: Compute a
 *    finite element function or its derivative in quadrature points.
 *
 *  <li> reinit: initialize the FEValues object for a certain cell.
 *    This function is not in the present class but only in the derived
 *    classes and has a variable call syntax.
 *    See the docs for the derived classes for more information.
 * </ul>
 *
 *
 * <h3>UpdateFlags</h3>
 *
 * The UpdateFlags object handed to the constructor is used to determine which
 * of the data fields to compute. This way, it is possible to avoid expensive
 * computations of useless derivatives.  In the beginning, these flags are
 * processed through the functions Mapping::update_once(),
 * Mapping::update_each(), FiniteElement::update_once()
 * FiniteElement::update_each(). All the results are bit-wise or'd and
 * determine the fields actually computed. This enables Mapping and
 * FiniteElement to schedule auxiliary data fields for updating. Still, it is
 * recommended to give <b>all</b> needed update flags to FEValues.
 *
 * The mechanisms by which this class works is also discussed on the page on
 * @ref UpdateFlagsEssay.
 *
 * @ingroup feaccess
 * @author Wolfgang Bangerth, 1998, 2003, Guido Kanschat, 2001
 */
template <int dim, int spacedim>
class FEValuesBase : protected FEValuesData<dim,spacedim>,
  public Subscriptor
{
public:
  /**
   * Dimension in which this object operates.
   */
  static const unsigned int dimension = dim;

  /**
   * Dimension of the space in which this object operates.
   */
  static const unsigned int space_dimension = spacedim;

  /**
   * Number of quadrature points.
   */
  const unsigned int n_quadrature_points;

  /**
   * Number of shape functions per cell. If we use this base class to evaluate
   * a finite element on faces of cells, this is still the number of degrees
   * of freedom per cell, not per face.
   */
  const unsigned int dofs_per_cell;


  /**
   * Constructor. Set up the array sizes with <tt>n_q_points</tt> quadrature
   * points, <tt>dofs_per_cell</tt> trial functions per cell and with the
   * given pattern to update the fields when the <tt>reinit</tt> function of
   * the derived classes is called. The fields themselves are not set up, this
   * must happen in the constructor of the derived class.
   */
  FEValuesBase (const unsigned int n_q_points,
                const unsigned int dofs_per_cell,
                const UpdateFlags update_flags,
                const Mapping<dim,spacedim> &mapping,
                const FiniteElement<dim,spacedim> &fe);


  /**
   * Destructor.
   */
  ~FEValuesBase ();
  /// @name ShapeAccess Access to shape function values
  //@{

  /**
   * Value of a shape function at a quadrature point on the cell, face or
   * subface selected the last time the <tt>reinit</tt> function of the
   * derived class was called.
   *
   * If the shape function is vector-valued, then this returns the only
   * non-zero component. If the shape function has more than one non-zero
   * component (i.e. it is not primitive), then throw an exception of type
   * ExcShapeFunctionNotPrimitive. In that case, use the
   * shape_value_component() function.
   *
   * @param function_no Number of the shape function to be evaluated. Note
   * that this number runs from zero to dofs_per_cell, even in the case of an
   * FEFaceValues or FESubfaceValues object.
   *
   * @param point_no Number of the quadrature point at which function is to be
   * evaluated
   *
   * @dealiiRequiresUpdateFlags{update_values}
   */
  const double &shape_value (const unsigned int function_no,
                             const unsigned int point_no) const;

  /**
   * Compute one vector component of the value of a shape function at a
   * quadrature point. If the finite element is scalar, then only component
   * zero is allowed and the return value equals that of the shape_value()
   * function. If the finite element is vector valued but all shape functions
   * are primitive (i.e. they are non-zero in only one component), then the
   * value returned by shape_value() equals that of this function for exactly
   * one component. This function is therefore only of greater interest if the
   * shape function is not primitive, but then it is necessary since the other
   * function cannot be used.
   *
   * @param function_no Number of the shape function to be evaluated.
   *
   * @param point_no Number of the quadrature point at which function is to be
   * evaluated.
   *
   * @param component vector component to be evaluated.
   *
   * @dealiiRequiresUpdateFlags{update_values}
   */
  double shape_value_component (const unsigned int function_no,
                                const unsigned int point_no,
                                const unsigned int component) const;

  /**
   * Compute the gradient of the <tt>function_no</tt>th shape function at the
   * <tt>quadrature_point</tt>th quadrature point with respect to real cell
   * coordinates.  If you want to get the derivative in one of the coordinate
   * directions, use the appropriate function of the Tensor class to extract
   * one component of the Tensor returned by this function. Since only a
   * reference to the gradient's value is returned, there should be no major
   * performance drawback.
   *
   * If the shape function is vector-valued, then this returns the only
   * non-zero component. If the shape function has more than one non-zero
   * component (i.e. it is not primitive), then it will throw an exception of
   * type ExcShapeFunctionNotPrimitive. In that case, use the
   * shape_grad_component() function.
   *
   * The same holds for the arguments of this function as for the
   * shape_value() function.
   *
   * @dealiiRequiresUpdateFlags{update_gradients}
   */
  const Tensor<1,spacedim> &
  shape_grad (const unsigned int function_no,
              const unsigned int quadrature_point) const;

  /**
   * Return one vector component of the gradient of a shape function at a
   * quadrature point. If the finite element is scalar, then only component
   * zero is allowed and the return value equals that of the shape_grad()
   * function. If the finite element is vector valued but all shape functions
   * are primitive (i.e. they are non-zero in only one component), then the
   * value returned by shape_grad() equals that of this function for exactly
   * one component. This function is therefore only of greater interest if the
   * shape function is not primitive, but then it is necessary since the other
   * function cannot be used.
   *
   * The same holds for the arguments of this function as for the
   * shape_value_component() function.
   *
   * @dealiiRequiresUpdateFlags{update_gradients}
   */
  Tensor<1,spacedim>
  shape_grad_component (const unsigned int function_no,
                        const unsigned int point_no,
                        const unsigned int component) const;

  /**
   * Second derivatives of the <tt>function_no</tt>th shape function at the
   * <tt>point_no</tt>th quadrature point with respect to real cell
   * coordinates. If you want to get the derivatives in one of the coordinate
   * directions, use the appropriate function of the Tensor class to extract
   * one component. Since only a reference to the derivative values is
   * returned, there should be no major performance drawback.
   *
   * If the shape function is vector-valued, then this returns the only
   * non-zero component. If the shape function has more than one non-zero
   * component (i.e. it is not primitive), then throw an exception of type
   * ExcShapeFunctionNotPrimitive. In that case, use the
   * shape_grad_grad_component() function.
   *
   * The same holds for the arguments of this function as for the
   * shape_value() function.
   *
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  const Tensor<2,spacedim> &
  shape_hessian (const unsigned int function_no,
                 const unsigned int point_no) const;

  /**
   * @deprecated Wrapper for shape_hessian()
   */
  const Tensor<2,spacedim> &
  shape_2nd_derivative (const unsigned int function_no,
                        const unsigned int point_no) const DEAL_II_DEPRECATED;


  /**
   * Return one vector component of the gradient of a shape function at a
   * quadrature point. If the finite element is scalar, then only component
   * zero is allowed and the return value equals that of the shape_hessian()
   * function. If the finite element is vector valued but all shape functions
   * are primitive (i.e. they are non-zero in only one component), then the
   * value returned by shape_hessian() equals that of this function for
   * exactly one component. This function is therefore only of greater
   * interest if the shape function is not primitive, but then it is necessary
   * since the other function cannot be used.
   *
   * The same holds for the arguments of this function as for the
   * shape_value_component() function.
   *
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  Tensor<2,spacedim>
  shape_hessian_component (const unsigned int function_no,
                           const unsigned int point_no,
                           const unsigned int component) const;

  /**
   * @deprecated Wrapper for shape_hessian_component()
   */
  Tensor<2,spacedim>
  shape_2nd_derivative_component (const unsigned int function_no,
                                  const unsigned int point_no,
                                  const unsigned int component) const DEAL_II_DEPRECATED;


  //@}
  /// @name Access to values of global finite element fields
  //@{

  /**
   * Returns the values of a finite element function restricted to the current
   * cell, face or subface selected the last time the <tt>reinit</tt> function
   * of the derived class was called, at the quadrature points.
   *
   * If the present cell is not active then values are interpolated to the
   * current cell and point values are computed from that.
   *
   * This function may only be used if the finite element in use is a scalar
   * one, i.e. has only one vector component.  To get values of
   * multi-component elements, there is another get_function_values() below,
   * returning a vector of vectors of results.
   *
   * @param[in] fe_function A vector of values that describes (globally) the
   * finite element function that this function should evaluate at the
   * quadrature points of the current cell.
   *
   * @param[out] values The values of the function specified by fe_function at
   * the quadrature points of the current cell.  The object is assume to
   * already have the correct size.
   *
   * @post <code>values[q]</code> will contain the value of the field
   * described by fe_function at the $q$th quadrature point.
   *
   * @note The actual data type of the input vector may be either a
   * Vector&lt;T&gt;, BlockVector&lt;T&gt;, or one of the sequential PETSc or
   * Trilinos vector wrapper classes. It represents a global vector of DoF
   * values associated with the DofHandler object with which this FEValues
   * object was last initialized. Alternatively, if the vector argument is of
   * type IndexSet, then the function is represented as one that is either
   * zero or one, depending on whether a DoF index is in the set or not.
   *
   * @dealiiRequiresUpdateFlags{update_values}
   */
  template <class InputVector, typename number>
  void get_function_values (const InputVector &fe_function,
                            std::vector<number> &values) const;

  /**
   * This function does the same as the other get_function_values(), but
   * applied to multi-component (vector-valued) elements. The meaning of the
   * arguments is as explained there.
   *
   * @post <code>values[q]</code> is a vector of values of the field described
   * by fe_function at the $q$th quadrature point. The size of the vector
   * accessed by <code>values[q]</code> equals the number of components of the
   * finite element, i.e. <code>values[q](c)</code> returns the value of the
   * $c$th vector component at the $q$th quadrature point.
   *
   * @dealiiRequiresUpdateFlags{update_values}
   */
  template <class InputVector, typename number>
  void get_function_values (const InputVector       &fe_function,
                            std::vector<Vector<number> > &values) const;

  /**
   * Generate function values from an arbitrary vector.
   *
   * This function offers the possibility to extract function values in
   * quadrature points from vectors not corresponding to a whole
   * discretization.
   *
   * The vector <tt>indices</tt> corresponds to the degrees of freedom on a
   * single cell. Its length may even be a multiple of the number of dofs per
   * cell. Then, the vectors in <tt>value</tt> should allow for the same
   * multiple of the components of the finite element.
   *
   * You may want to use this function, if you want to access just a single
   * block from a BlockVector, if you have a multi-level vector or if you
   * already have a local representation of your finite element data.
   *
   * @dealiiRequiresUpdateFlags{update_values}
   */
  template <class InputVector, typename number>
  void get_function_values (const InputVector &fe_function,
                            const VectorSlice<const std::vector<types::global_dof_index> > &indices,
                            std::vector<number> &values) const;

  /**
   * Generate vector function values from an arbitrary vector.
   *
   * This function offers the possibility to extract function values in
   * quadrature points from vectors not corresponding to a whole
   * discretization.
   *
   * The vector <tt>indices</tt> corresponds to the degrees of freedom on a
   * single cell. Its length may even be a multiple of the number of dofs per
   * cell. Then, the vectors in <tt>value</tt> should allow for the same
   * multiple of the components of the finite element.
   *
   * You may want to use this function, if you want to access just a single
   * block from a BlockVector, if you have a multi-level vector or if you
   * already have a local representation of your finite element data.
   *
   * Since this function allows for fairly general combinations of argument
   * sizes, be aware that the checks on the arguments may not detect errors.
   *
   * @dealiiRequiresUpdateFlags{update_values}
   */
  template <class InputVector, typename number>
  void get_function_values (const InputVector &fe_function,
                            const VectorSlice<const std::vector<types::global_dof_index> > &indices,
                            std::vector<Vector<number> > &values) const;


  /**
   * Generate vector function values from an arbitrary vector.
   *
   * This function offers the possibility to extract function values in
   * quadrature points from vectors not corresponding to a whole
   * discretization.
   *
   * The vector <tt>indices</tt> corresponds to the degrees of freedom on a
   * single cell. Its length may even be a multiple of the number of dofs per
   * cell. Then, the vectors in <tt>value</tt> should allow for the same
   * multiple of the components of the finite element.
   *
   * Depending on the value of the last argument, the outer vector of
   * <tt>values</tt> has either the length of the quadrature rule
   * (<tt>quadrature_points_fastest == false</tt>) or the length of components
   * to be filled <tt>quadrature_points_fastest == true</tt>. If <tt>p</tt> is
   * the current quadrature point number and <tt>i</tt> is the vector
   * component of the solution desired, the access to <tt>values</tt> is
   * <tt>values[p][i]</tt> if <tt>quadrature_points_fastest == false</tt>, and
   * <tt>values[i][p]</tt> otherwise.
   *
   * You may want to use this function, if you want to access just a single
   * block from a BlockVector, if you have a multi-level vector or if you
   * already have a local representation of your finite element data.
   *
   * Since this function allows for fairly general combinations of argument
   * sizes, be aware that the checks on the arguments may not detect errors.
   *
   * @dealiiRequiresUpdateFlags{update_values}
   */
  template <class InputVector>
  void get_function_values (const InputVector &fe_function,
                            const VectorSlice<const std::vector<types::global_dof_index> > &indices,
                            VectorSlice<std::vector<std::vector<double> > > values,
                            const bool quadrature_points_fastest) const;

  //@}
  /// @name Access to derivatives of global finite element fields
  //@{

  /**
   * Compute the gradients of a finite element at the quadrature points of a
   * cell. This function is the equivalent of the corresponding
   * get_function_values() function (see there for more information) but
   * evaluates the finite element field's gradient instead of its value.
   *
   * This function may only be used if the finite element in use is a scalar
   * one, i.e. has only one vector component. There is a corresponding
   * function of the same name for vector-valued finite elements.
   *
   * @param[in] fe_function A vector of values that describes (globally) the
   * finite element function that this function should evaluate at the
   * quadrature points of the current cell.
   *
   * @param[out] gradients The gradients of the function specified by
   * fe_function at the quadrature points of the current cell.  The gradients
   * are computed in real space (as opposed to on the unit cell).  The object
   * is assume to already have the correct size.
   *
   * @post <code>gradients[q]</code> will contain the gradient of the field
   * described by fe_function at the $q$th quadrature
   * point. <code>gradients[q][d]</code> represents the derivative in
   * coordinate direction $d$ at quadrature point $q$.
   *
   * @note The actual data type of the input vector may be either a
   * Vector&lt;T&gt;, BlockVector&lt;T&gt;, or one of the sequential PETSc or
   * Trilinos vector wrapper classes. It represents a global vector of DoF
   * values associated with the DofHandler object with which this FEValues
   * object was last initialized. Alternatively, if the vector argument is of
   * type IndexSet, then the function is represented as one that is either
   * zero or one, depending on whether a DoF index is in the set or not.
   *
   * @dealiiRequiresUpdateFlags{update_gradients}
   */
  template <class InputVector>
  void get_function_gradients (const InputVector      &fe_function,
                               std::vector<Tensor<1,spacedim> > &gradients) const;

  /**
   * This function does the same as the other get_function_gradients(), but
   * applied to multi-component (vector-valued) elements. The meaning of the
   * arguments is as explained there.
   *
   * @post <code>gradients[q]</code> is a vector of gradients of the field
   * described by fe_function at the $q$th quadrature point. The size of the
   * vector accessed by <code>gradients[q]</code> equals the number of
   * components of the finite element, i.e. <code>gradients[q][c]</code>
   * returns the gradient of the $c$th vector component at the $q$th
   * quadrature point. Consequently, <code>gradients[q][c][d]</code> is the
   * derivative in coordinate direction $d$ of the $c$th vector component of
   * the vector field at quadrature point $q$ of the current cell.
   *
   * @dealiiRequiresUpdateFlags{update_gradients}
   */
  template <class InputVector>
  void get_function_gradients (const InputVector               &fe_function,
                               std::vector<std::vector<Tensor<1,spacedim> > > &gradients) const;

  /**
   * Function gradient access with more flexibility. See get_function_values()
   * with corresponding arguments.
   *
   * @dealiiRequiresUpdateFlags{update_gradients}
   */
  template <class InputVector>
  void get_function_gradients (const InputVector &fe_function,
                               const VectorSlice<const std::vector<types::global_dof_index> > &indices,
                               std::vector<Tensor<1,spacedim> > &gradients) const;

  /**
   * Function gradient access with more flexibility. See get_function_values()
   * with corresponding arguments.
   *
   * @dealiiRequiresUpdateFlags{update_gradients}
   */
  template <class InputVector>
  void get_function_gradients (const InputVector &fe_function,
                               const VectorSlice<const std::vector<types::global_dof_index> > &indices,
                               VectorSlice<std::vector<std::vector<Tensor<1,spacedim> > > > gradients,
                               bool quadrature_points_fastest = false) const;

  /**
   * @deprecated Use get_function_gradients() instead.
   */
  template <class InputVector>
  void get_function_grads (const InputVector      &fe_function,
                           std::vector<Tensor<1,spacedim> > &gradients) const DEAL_II_DEPRECATED;

  /**
   * @deprecated Use get_function_gradients() instead.
   */
  template <class InputVector>
  void get_function_grads (const InputVector               &fe_function,
                           std::vector<std::vector<Tensor<1,spacedim> > > &gradients) const DEAL_II_DEPRECATED;

  /**
   * @deprecated Use get_function_gradients() instead.
   */
  template <class InputVector>
  void get_function_grads (const InputVector &fe_function,
                           const VectorSlice<const std::vector<types::global_dof_index> > &indices,
                           std::vector<Tensor<1,spacedim> > &gradients) const DEAL_II_DEPRECATED;

  /**
   * @deprecated Use get_function_gradients() instead.
   */
  template <class InputVector>
  void get_function_grads (const InputVector &fe_function,
                           const VectorSlice<const std::vector<types::global_dof_index> > &indices,
                           std::vector<std::vector<Tensor<1,spacedim> > > &gradients,
                           bool quadrature_points_fastest = false) const DEAL_II_DEPRECATED;

  //@}
  /// @name Access to second derivatives (Hessian matrices and Laplacians) of global finite element fields
  //@{

  /**
   * Compute the tensor of second derivatives of a finite element at the
   * quadrature points of a cell. This function is the equivalent of the
   * corresponding get_function_values() function (see there for more
   * information) but evaluates the finite element field's second derivatives
   * instead of its value.
   *
   * This function may only be used if the finite element in use is a scalar
   * one, i.e. has only one vector component. There is a corresponding
   * function of the same name for vector-valued finite elements.
   *
   * @param[in] fe_function A vector of values that describes (globally) the
   * finite element function that this function should evaluate at the
   * quadrature points of the current cell.
   *
   * @param[out] hessians The Hessians of the function specified by
   * fe_function at the quadrature points of the current cell.  The Hessians
   * are computed in real space (as opposed to on the unit cell).  The object
   * is assume to already have the correct size.
   *
   * @post <code>hessians[q]</code> will contain the Hessian of the field
   * described by fe_function at the $q$th quadrature
   * point. <code>hessians[q][i][j]</code> represents the $(i,j)$th component
   * of the matrix of second derivatives at quadrature point $q$.
   *
   * @note The actual data type of the input vector may be either a
   * Vector&lt;T&gt;, BlockVector&lt;T&gt;, or one of the sequential PETSc or
   * Trilinos vector wrapper classes. It represents a global vector of DoF
   * values associated with the DofHandler object with which this FEValues
   * object was last initialized. Alternatively, if the vector argument is of
   * type IndexSet, then the function is represented as one that is either
   * zero or one, depending on whether a DoF index is in the set or not.
   *
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  template <class InputVector>
  void
  get_function_hessians (const InputVector &fe_function,
                         std::vector<Tensor<2,spacedim> > &hessians) const;

  /**
   * This function does the same as the other get_function_hessians(), but
   * applied to multi-component (vector-valued) elements. The meaning of the
   * arguments is as explained there.
   *
   * @post <code>hessians[q]</code> is a vector of Hessians of the field
   * described by fe_function at the $q$th quadrature point. The size of the
   * vector accessed by <code>hessians[q]</code> equals the number of
   * components of the finite element, i.e. <code>hessians[q][c]</code>
   * returns the Hessian of the $c$th vector component at the $q$th quadrature
   * point. Consequently, <code>hessians[q][c][i][j]</code> is the $(i,j)$th
   * component of the matrix of second derivatives of the $c$th vector
   * component of the vector field at quadrature point $q$ of the current
   * cell.
   *
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  template <class InputVector>
  void
  get_function_hessians (const InputVector      &fe_function,
                         std::vector<std::vector<Tensor<2,spacedim> > > &hessians,
                         bool quadrature_points_fastest = false) const;

  /**
   * Access to the second derivatives of a function with more flexibility. See
   * get_function_values() with corresponding arguments.
   */
  template <class InputVector>
  void get_function_hessians (
    const InputVector &fe_function,
    const VectorSlice<const std::vector<types::global_dof_index> > &indices,
    std::vector<Tensor<2,spacedim> > &hessians) const;

  /**
   * Access to the second derivatives of a function with more flexibility. See
   * get_function_values() with corresponding arguments.
   *
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  template <class InputVector>
  void get_function_hessians (
    const InputVector &fe_function,
    const VectorSlice<const std::vector<types::global_dof_index> > &indices,
    VectorSlice<std::vector<std::vector<Tensor<2,spacedim> > > > hessians,
    bool quadrature_points_fastest = false) const;

  /**
   * @deprecated Wrapper for get_function_hessians()
   */
  template <class InputVector>
  void
  get_function_2nd_derivatives (const InputVector &,
                                std::vector<Tensor<2,spacedim> > &) const DEAL_II_DEPRECATED;

  /**
   * @deprecated Wrapper for get_function_hessians()
   */
  template <class InputVector>
  void
  get_function_2nd_derivatives (const InputVector &,
                                std::vector<std::vector<Tensor<2,spacedim> > > &,
                                bool = false) const DEAL_II_DEPRECATED;

  /**
   * Compute the (scalar) Laplacian (i.e. the trace of the tensor of second
   * derivatives) of a finite element at the quadrature points of a cell. This
   * function is the equivalent of the corresponding get_function_values()
   * function (see there for more information) but evaluates the finite
   * element field's second derivatives instead of its value.
   *
   * This function may only be used if the finite element in use is a scalar
   * one, i.e. has only one vector component. There is a corresponding
   * function of the same name for vector-valued finite elements.
   *
   * @param[in] fe_function A vector of values that describes (globally) the
   * finite element function that this function should evaluate at the
   * quadrature points of the current cell.
   *
   * @param[out] laplacians The Laplacians of the function specified by
   * fe_function at the quadrature points of the current cell.  The Laplacians
   * are computed in real space (as opposed to on the unit cell).  The object
   * is assume to already have the correct size.
   *
   * @post <code>laplacians[q]</code> will contain the Laplacian of the field
   * described by fe_function at the $q$th quadrature
   * point. <code>gradients[q][i][j]</code> represents the $(i,j)$th component
   * of the matrix of second derivatives at quadrature point $q$.
   *
   * @post For each component of the output vector, there holds
   * <code>laplacians[q]=trace(hessians[q])</code>, where <tt>hessians</tt>
   * would be the output of the get_function_hessians() function.
   *
   * @note The actual data type of the input vector may be either a
   * Vector&lt;T&gt;, BlockVector&lt;T&gt;, or one of the sequential PETSc or
   * Trilinos vector wrapper classes. It represents a global vector of DoF
   * values associated with the DofHandler object with which this FEValues
   * object was last initialized. Alternatively, if the vector argument is of
   * type IndexSet, then the function is represented as one that is either
   * zero or one, depending on whether a DoF index is in the set or not.
   *
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  template <class InputVector, typename number>
  void
  get_function_laplacians (const InputVector &fe_function,
                           std::vector<number> &laplacians) const;

  /**
   * This function does the same as the other get_function_laplacians(), but
   * applied to multi-component (vector-valued) elements. The meaning of the
   * arguments is as explained there.
   *
   * @post <code>laplacians[q]</code> is a vector of Laplacians of the field
   * described by fe_function at the $q$th quadrature point. The size of the
   * vector accessed by <code>laplacians[q]</code> equals the number of
   * components of the finite element, i.e. <code>laplacians[q][c]</code>
   * returns the Laplacian of the $c$th vector component at the $q$th
   * quadrature point.
   *
   * @post For each component of the output vector, there holds
   * <code>laplacians[q][c]=trace(hessians[q][c])</code>, where
   * <tt>hessians</tt> would be the output of the get_function_hessians()
   * function.
   *
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  template <class InputVector, typename number>
  void
  get_function_laplacians (const InputVector      &fe_function,
                           std::vector<Vector<number> > &laplacians) const;

  /**
   * Access to the second derivatives of a function with more flexibility. See
   * get_function_values() with corresponding arguments.
   *
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  template <class InputVector, typename number>
  void get_function_laplacians (
    const InputVector &fe_function,
    const VectorSlice<const std::vector<types::global_dof_index> > &indices,
    std::vector<number> &laplacians) const;

  /**
   * Access to the second derivatives of a function with more flexibility. See
   * get_function_values() with corresponding arguments.
   *
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  template <class InputVector, typename number>
  void get_function_laplacians (
    const InputVector &fe_function,
    const VectorSlice<const std::vector<types::global_dof_index> > &indices,
    std::vector<Vector<number> > &laplacians) const;

  /**
   * Access to the second derivatives of a function with more flexibility. See
   * get_function_values() with corresponding arguments.
   *
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  template <class InputVector, typename number>
  void get_function_laplacians (
    const InputVector &fe_function,
    const VectorSlice<const std::vector<types::global_dof_index> > &indices,
    std::vector<std::vector<number> > &laplacians,
    bool quadrature_points_fastest = false) const;
  //@}

  /// @name Geometry of the cell
  //@{

  /**
   * Position of the <tt>q</tt>th quadrature point in real space.
   *
   * @dealiiRequiresUpdateFlags{update_quadrature_points}
   */
  const Point<spacedim> &
  quadrature_point (const unsigned int q) const;

  /**
   * Return a reference to the vector of quadrature points in real space.
   *
   * @dealiiRequiresUpdateFlags{update_quadrature_points}
   */
  const std::vector<Point<spacedim> > &get_quadrature_points () const;

  /**
   * Mapped quadrature weight. If this object refers to a volume evaluation
   * (i.e. the derived class is of type FEValues), then this is the Jacobi
   * determinant times the weight of the *<tt>i</tt>th unit quadrature point.
   *
   * For surface evaluations (i.e. classes FEFaceValues or FESubfaceValues),
   * it is the mapped surface element times the weight of the quadrature
   * point.
   *
   * You can think of the quantity returned by this function as the volume or
   * surface element $dx, ds$ in the integral that we implement here by
   * quadrature.
   *
   * @dealiiRequiresUpdateFlags{update_JxW_values}
   */
  double JxW (const unsigned int quadrature_point) const;

  /**
   * Return a reference to the array holding the values returned by JxW().
   */
  const std::vector<double> &get_JxW_values () const;

  /**
   * Return the Jacobian of the transformation at the specified quadrature
   * point, i.e.  $J_{ij}=dx_i/d\hat x_j$
   *
   * @dealiiRequiresUpdateFlags{update_jacobians}
   */
  const DerivativeForm<1,dim,spacedim> &jacobian (const unsigned int quadrature_point) const;

  /**
   * Return a reference to the array holding the values returned by jacobian().
   *
   * @dealiiRequiresUpdateFlags{update_jacobians}
   */
  const std::vector<DerivativeForm<1,dim,spacedim> > &get_jacobians () const;

  /**
   * Return the second derivative of the transformation from unit to real
   * cell, i.e. the first derivative of the Jacobian, at the specified
   * quadrature point, i.e. $G_{ijk}=dJ_{jk}/d\hat x_i$.
   *
   * @dealiiRequiresUpdateFlags{update_jacobian_grads}
   */
  const DerivativeForm<2,dim,spacedim> &jacobian_grad (const unsigned int quadrature_point) const;

  /**
   * Return a reference to the array holding the values returned by jacobian_grads().
   *
   * @dealiiRequiresUpdateFlags{update_jacobian_grads}
   */
  const std::vector<DerivativeForm<2,dim,spacedim> > &get_jacobian_grads () const;

  /**
   * Return the inverse Jacobian of the transformation at the specified
   * quadrature point, i.e.  $J_{ij}=d\hat x_i/dx_j$
   *
   * @dealiiRequiresUpdateFlags{update_inverse_jacobians}
   */
  const DerivativeForm<1,spacedim,dim> &inverse_jacobian (const unsigned int quadrature_point) const;

  /**
   * Return a reference to the array holding the values returned by inverse_jacobian().
   *
   * @dealiiRequiresUpdateFlags{update_inverse_jacobians}
   */
  const std::vector<DerivativeForm<1,spacedim,dim> > &get_inverse_jacobians () const;

  /**
   * For a face, return the outward normal vector to the cell at the
   * <tt>i</tt>th quadrature point.
   *
   * For a cell of codimension one, return the normal vector, as it is
   * specified by the numbering of the vertices.
   *
   * The length of the vector is normalized to one.
   *
   * @dealiiRequiresUpdateFlags{update_normal_vectors}
   */
  const Point<spacedim> &normal_vector (const unsigned int i) const;

  /**
   * Return the normal vectors at the quadrature points. For a face, these are
   * the outward normal vectors to the cell. For a cell of codimension one,
   * the orientation is given by the numbering of vertices.
   *
   * @dealiiRequiresUpdateFlags{update_normal_vectors}
   */
  const std::vector<Point<spacedim> > &get_normal_vectors () const;

  /**
   * Transform a set of vectors, one for each quadrature point. The
   * <tt>mapping</tt> can be any of the ones defined in MappingType.
   */
  void transform (std::vector<Tensor<1,spacedim> > &transformed,
                  const std::vector<Tensor<1,dim> > &original,
                  MappingType mapping) const;

  /**
   * @deprecated Use normal_vector() instead.
   *
   * Return the outward normal vector to the cell at the <tt>i</tt>th
   * quadrature point. The length of the vector is normalized to one.
   */
  const Point<spacedim> &cell_normal_vector (const unsigned int i) const DEAL_II_DEPRECATED;

  /**
   * @deprecated Use get_normal_vectors() instead.
   *
   * Returns the vectors normal to the cell in each of the quadrature points.
   */
  const std::vector<Point<spacedim> > &get_cell_normal_vectors () const DEAL_II_DEPRECATED;

  //@}

  /// @name Extractors Methods to extract individual components
  //@{

  /**
   * Create a view of the current FEValues object that represents a particular
   * scalar component of the possibly vector-valued finite element. The
   * concept of views is explained in the documentation of the namespace
   * FEValuesViews and in particular in the @ref vector_valued module.
   */
  const FEValuesViews::Scalar<dim,spacedim> &
  operator[] (const FEValuesExtractors::Scalar &scalar) const;

  /**
   * Create a view of the current FEValues object that represents a set of
   * <code>dim</code> scalar components (i.e. a vector) of the vector-valued
   * finite element. The concept of views is explained in the documentation of
   * the namespace FEValuesViews and in particular in the @ref vector_valued
   * module.
   */
  const FEValuesViews::Vector<dim,spacedim> &
  operator[] (const FEValuesExtractors::Vector &vector) const;

  /**
   * Create a view of the current FEValues object that represents a set of
   * <code>(dim*dim + dim)/2</code> scalar components (i.e. a symmetric 2nd
   * order tensor) of the vector-valued finite element. The concept of views
   * is explained in the documentation of the namespace FEValuesViews and in
   * particular in the @ref vector_valued module.
   */
  const FEValuesViews::SymmetricTensor<2,dim,spacedim> &
  operator[] (const FEValuesExtractors::SymmetricTensor<2> &tensor) const;


  /**
   * Create a view of the current FEValues object that represents a set of
   * <code>(dim*dim)</code> scalar components (i.e. a 2nd order tensor) of the
   * vector-valued finite element. The concept of views is explained in the
   * documentation of the namespace FEValuesViews and in particular in the
   * @ref vector_valued module.
   */
  const FEValuesViews::Tensor<2,dim,spacedim> &
  operator[] (const FEValuesExtractors::Tensor<2> &tensor) const;

  //@}

  /// @name Access to the raw data
  //@{

  /**
   * Constant reference to the selected mapping object.
   */
  const Mapping<dim,spacedim> &get_mapping () const;

  /**
   * Constant reference to the selected finite element object.
   */
  const FiniteElement<dim,spacedim> &get_fe () const;

  /**
   * Return the update flags set for this object.
   */
  UpdateFlags get_update_flags () const;

  /**
   * Return a triangulation iterator to the current cell.
   */
  const typename Triangulation<dim,spacedim>::cell_iterator get_cell () const;

  /**
   * Return the relation of the current cell to the previous cell. This allows
   * re-use of some cell data (like local matrices for equations with constant
   * coefficients) if the result is <tt>CellSimilarity::translation</tt>.
   */
  CellSimilarity::Similarity get_cell_similarity () const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  std::size_t memory_consumption () const;
  //@}


  /**
   * This exception is thrown if FEValuesBase is asked to return the value of
   * a field which was not required by the UpdateFlags for this FEValuesBase.
   *
   * @ingroup Exceptions
   */
  DeclException1 (ExcAccessToUninitializedField,
                  char *,
                  << ("You are requesting information from an FEValues/FEFaceValues/FESubfaceValues "
                      "object for which this kind of information has not been computed. What "
                      "information these objects compute is determined by the update_* flags you "
                      "pass to the constructor. Here, the operation you are attempting requires "
                      "the <")
                  << arg1
                  << "> flag to be set, but it was apparently not specified upon construction.");
  /**
   * @todo Document this
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcCannotInitializeField);
  /**
   * @todo Document this
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcInvalidUpdateFlag);
  /**
   * @todo Document this
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcFEDontMatch);
  /**
   * @todo Document this
   *
   * @ingroup Exceptions
   */
  DeclException1 (ExcShapeFunctionNotPrimitive,
                  int,
                  << "The shape function with index " << arg1
                  << " is not primitive, i.e. it is vector-valued and "
                  << "has more than one non-zero vector component. This "
                  << "function cannot be called for these shape functions. "
                  << "Maybe you want to use the same function with the "
                  << "_component suffix?");
  /**
   * @todo Document this
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcFENotPrimitive);

protected:
  /**
   * Objects of the FEValues class need to store a pointer (i.e. an iterator)
   * to the present cell in order to be able to extract the values of the
   * degrees of freedom on this cell in the get_function_values() and assorted
   * functions. On the other hand, this class should also work for different
   * iterators, as long as they have the same interface to extract the DoF
   * values (i.e., for example, they need to have a @p
   * get_interpolated_dof_values function).
   *
   * This calls for a common base class of iterator classes, and making the
   * functions we need here @p virtual. On the other hand, this is the only
   * place in the library where we need this, and introducing a base class of
   * iterators and making a function virtual penalizes <em>all</em> users of
   * the iterators, which are basically intended as very fast accessor
   * functions. So we do not want to do this. Rather, what we do here is
   * making the functions we need virtual only for use with <em>this
   * class</em>. The idea is the following: have a common base class which
   * declares some pure virtual functions, and for each possible iterator
   * type, we have a derived class which stores the iterator to the cell and
   * implements these functions. Since the iterator classes have the same
   * interface, we can make the derived classes a template, templatized on the
   * iterator type.
   *
   * This way, the use of virtual functions is restricted to only this class,
   * and other users of iterators do not have to bear the negative effects.
   *
   * @author Wolfgang Bangerth, 2003
   */
  class CellIteratorBase;

  /**
   * Forward declaration of classes derived from CellIteratorBase. Their
   * definition and implementation is given in the .cc file.
   */
  template <typename CI> class CellIterator;
  class TriaCellIterator;

  /**
   * Store the cell selected last time the reinit() function was called.  This
   * is necessary for the <tt>get_function_*</tt> functions as well as the
   * functions of same name in the extractor classes.
   */
  std::auto_ptr<const CellIteratorBase> present_cell;

  /**
   * A signal connection we use to ensure we get informed whenever the
   * triangulation changes. We need to know about that because it invalidates
   * all cell iterators and, as part of that, the 'present_cell' iterator we
   * keep around between subsequent calls to reinit() in order to compute the
   * cell similarity.
   */
  boost::signals2::connection tria_listener;

  /**
   * A function that is connected to the triangulation in order to reset the
   * stored 'present_cell' iterator to an invalid one whenever the
   * triangulation is changed and the iterator consequently becomes invalid.
   */
  void invalidate_present_cell ();

  /**
   * This function is called by the various reinit() functions in derived
   * classes. Given the cell indicated by the argument, test whether we have
   * to throw away the previously stored present_cell argument because it
   * would require us to compare cells from different triangulations. In
   * checking all this, also make sure that we have tria_listener connected to
   * the triangulation to which we will set present_cell right after calling
   * this function.
   */
  void
  maybe_invalidate_previous_present_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell);

  /**
   * Storage for the mapping object.
   */
  const SmartPointer<const Mapping<dim,spacedim>,FEValuesBase<dim,spacedim> > mapping;

  /**
   * Store the finite element for later use.
   */
  const SmartPointer<const FiniteElement<dim,spacedim>,FEValuesBase<dim,spacedim> > fe;


  /**
   * Internal data of mapping.
   */
  SmartPointer<typename Mapping<dim,spacedim>::InternalDataBase,FEValuesBase<dim,spacedim> > mapping_data;

  /**
   * Internal data of finite element.
   */
  SmartPointer<typename Mapping<dim,spacedim>::InternalDataBase,FEValuesBase<dim,spacedim> > fe_data;

  /**
   * Initialize some update flags. Called from the @p initialize functions of
   * derived classes, which are in turn called from their constructors.
   *
   * Basically, this function finds out using the finite element and mapping
   * object already stored which flags need to be set to compute everything
   * the user wants, as expressed through the flags passed as argument.
   */
  UpdateFlags compute_update_flags (const UpdateFlags update_flags) const;

  /**
   * An enum variable that can store different states of the current cell in
   * comparison to the previously visited cell. If wanted, additional states
   * can be checked here and used in one of the methods used during reinit.
   */
  CellSimilarity::Similarity cell_similarity;

  /**
   * A function that checks whether the new cell is similar to the one
   * previously used. Then, a significant amount of the data can be reused,
   * e.g. the derivatives of the basis functions in real space, shape_grad.
   */
  void
  check_cell_similarity (const typename Triangulation<dim,spacedim>::cell_iterator &cell);

private:
  /**
   * Copy constructor. Since objects of this class are not copyable, we make
   * it private, and also do not implement it.
   */
  FEValuesBase (const FEValuesBase &);

  /**
   * Copy operator. Since objects of this class are not copyable, we make it
   * private, and also do not implement it.
   */
  FEValuesBase &operator= (const FEValuesBase &);

  /**
   * A cache for all possible FEValuesViews objects.
   */
  dealii::internal::FEValuesViews::Cache<dim,spacedim> fe_values_views_cache;

  /**
   * Make the view classes friends of this class, since they access internal
   * data.
   */
  template <int, int> friend class FEValuesViews::Scalar;
  template <int, int> friend class FEValuesViews::Vector;
  template <int, int, int> friend class FEValuesViews::SymmetricTensor;
  template <int, int, int> friend class FEValuesViews::Tensor;
};



/**
 * Finite element evaluated in quadrature points of a cell.
 *
 * This function implements the initialization routines for FEValuesBase, if
 * values in quadrature points of a cell are needed. For further documentation
 * see this class.
 *
 * @ingroup feaccess
 * @author Wolfgang Bangerth, 1998, Guido Kanschat, 2001
 */
template <int dim, int spacedim=dim>
class FEValues : public FEValuesBase<dim,spacedim>
{
public:
  /**
   * Dimension of the object over which we integrate. For the present class,
   * this is equal to <code>dim</code>.
   */
  static const unsigned int integral_dimension = dim;

  /**
   * Constructor. Gets cell independent data from mapping and finite element
   * objects, matching the quadrature rule and update flags.
   */
  FEValues (const Mapping<dim,spacedim>       &mapping,
            const FiniteElement<dim,spacedim> &fe,
            const Quadrature<dim>             &quadrature,
            const UpdateFlags                  update_flags);

  /**
   * Constructor. Uses MappingQ1 implicitly.
   */
  FEValues (const FiniteElement<dim,spacedim> &fe,
            const Quadrature<dim>             &quadrature,
            const UpdateFlags                  update_flags);

  /**
   * Reinitialize the gradients, Jacobi determinants, etc for the given cell
   * of type "iterator into a DoFHandler object", and the finite element
   * associated with this object. It is assumed that the finite element used
   * by the given cell is also the one used by this FEValues object.
   */
  template <class DH, bool level_dof_access>
  void reinit (const TriaIterator<DoFCellAccessor<DH,level_dof_access> > cell);

  /**
   * Reinitialize the gradients, Jacobi determinants, etc for the given cell
   * of type "iterator into a Triangulation object", and the given finite
   * element. Since iterators into triangulation alone only convey information
   * about the geometry of a cell, but not about degrees of freedom possibly
   * associated with this cell, you will not be able to call some functions of
   * this class if they need information about degrees of freedom. These
   * functions are, above all, the
   * <tt>get_function_value/gradients/hessians/laplacians</tt> functions. If
   * you want to call these functions, you have to call the @p reinit variants
   * that take iterators into DoFHandler or other DoF handler type objects.
   */
  void reinit (const typename Triangulation<dim,spacedim>::cell_iterator &cell);

  /**
   * Return a reference to the copy of the quadrature formula stored by this
   * object.
   */
  const Quadrature<dim> &get_quadrature () const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  std::size_t memory_consumption () const;

  /**
   * Return a reference to this very object.
   *
   * Though it seems that it is not very useful, this function is there to
   * provide capability to the hpFEValues class, in which case it provides the
   * FEValues object for the present cell (remember that for hp finite
   * elements, the actual FE object used may change from cell to cell, so we
   * also need different FEValues objects for different cells; once you
   * reinitialize the hpFEValues object for a specific cell, it retrieves the
   * FEValues object for the FE on that cell and returns it through a function
   * of the same name as this one; this function here therefore only provides
   * the same interface so that one can templatize on FEValues/hpFEValues).
   */
  const FEValues<dim,spacedim> &get_present_fe_values () const;

private:
  /**
   * Store a copy of the quadrature formula here.
   */
  const Quadrature<dim> quadrature;

  /**
   * Do work common to the two constructors.
   */
  void initialize (const UpdateFlags update_flags);

  /**
   * The reinit() functions do only that part of the work that requires
   * knowledge of the type of iterator. After setting present_cell(), they
   * pass on to this function, which does the real work, and which is
   * independent of the actual type of the cell iterator.
   */
  void do_reinit ();
};


/**
 * Extend the interface of FEValuesBase to values that only make sense when
 * evaluating something on the surface of a cell. All the data that is
 * available in the interior of cells is also available here.
 *
 * See FEValuesBase
 *
 * @ingroup feaccess
 *  @author Wolfgang Bangerth, 1998, Guido Kanschat, 2000, 2001
 */
template <int dim, int spacedim=dim>
class FEFaceValuesBase : public FEValuesBase<dim,spacedim>
{
public:
  /**
   * Dimension of the object over which we integrate. For the present class,
   * this is equal to <code>dim-1</code>.
   */
  static const unsigned int integral_dimension = dim-1;

  /**
   * Constructor. Call the constructor of the base class and set up the arrays
   * of this class with the right sizes.  Actually filling these arrays is a
   * duty of the derived class's constructors.
   *
   * @p n_faces_or_subfaces is the number of faces or subfaces that this
   * object is to store. The actual number depends on the derived class, for
   * FEFaceValues it is <tt>2*dim</tt>, while for the FESubfaceValues class it
   * is <tt>2*dim*(1<<(dim-1))</tt>, i.e. the number of faces times the number
   * of subfaces per face.
   */
  FEFaceValuesBase (const unsigned int                 n_q_points,
                    const unsigned int                 dofs_per_cell,
                    const UpdateFlags                  update_flags,
                    const Mapping<dim,spacedim>       &mapping,
                    const FiniteElement<dim,spacedim> &fe,
                    const Quadrature<dim-1>&           quadrature);

  /**
   * Boundary form of the transformation of the cell at the <tt>i</tt>th
   * quadrature point.  See @ref GlossBoundaryForm .
   *
   * @dealiiRequiresUpdateFlags{update_boundary_forms}
   */
  const Tensor<1,spacedim> &boundary_form (const unsigned int i) const;

  /**
   * Return the list of outward normal vectors times the Jacobian of the
   * surface mapping.
   *
   * @dealiiRequiresUpdateFlags{update_boundary_forms}
   */
  const std::vector<Tensor<1,spacedim> > &get_boundary_forms () const;

  /**
   * Return the index of the face selected the last time the reinit() function
   * was called.
   */
  unsigned int get_face_index() const;

  /**
   * Return a reference to the copy of the quadrature formula stored by this
   * object.
   */
  const Quadrature<dim-1> & get_quadrature () const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  std::size_t memory_consumption () const;

protected:

  /**
   * Index of the face selected the last time the reinit() function was
   * called.
   */
  unsigned int present_face_index;

  /**
   * Store a copy of the quadrature formula here.
   */
  const Quadrature<dim-1> quadrature;
};



/**
 * Finite element evaluated in quadrature points on a face.
 *
 * This class adds the functionality of FEFaceValuesBase to FEValues; see
 * there for more documentation.
 *
 * Since finite element functions and their derivatives may be discontinuous
 * at cell boundaries, there is no restriction of this function to a mesh
 * face. But, there are limits of these values approaching the face from
 * either of the neighboring cells.
 *
 * @ingroup feaccess
 * @author Wolfgang Bangerth, 1998, Guido Kanschat, 2000, 2001
 */
template <int dim, int spacedim=dim>
class FEFaceValues : public FEFaceValuesBase<dim,spacedim>
{
public:
  /**
   * Dimension in which this object operates.
   */

  static const unsigned int dimension = dim;

  static const unsigned int space_dimension = spacedim;

  /**
   * Dimension of the object over which we integrate. For the present class,
   * this is equal to <code>dim-1</code>.
   */
  static const unsigned int integral_dimension = dim-1;

  /**
   * Constructor. Gets cell independent data from mapping and finite element
   * objects, matching the quadrature rule and update flags.
   */
  FEFaceValues (const Mapping<dim,spacedim>       &mapping,
                const FiniteElement<dim,spacedim> &fe,
                const Quadrature<dim-1>           &quadrature,
                const UpdateFlags                  update_flags);

  /**
   * Constructor. Uses MappingQ1 implicitly.
   */
  FEFaceValues (const FiniteElement<dim,spacedim> &fe,
                const Quadrature<dim-1>           &quadrature,
                const UpdateFlags                  update_flags);

  /**
   * Reinitialize the gradients, Jacobi determinants, etc for the face with
   * number @p face_no of @p cell and the given finite element.
   */
  template <class DH, bool level_dof_access>
  void reinit (const TriaIterator<DoFCellAccessor<DH,level_dof_access> > cell,
               const unsigned int face_no);

  /**
   * Reinitialize the gradients, Jacobi determinants, etc for the given face
   * on given cell of type "iterator into a Triangulation object", and the
   * given finite element. Since iterators into triangulation alone only
   * convey information about the geometry of a cell, but not about degrees of
   * freedom possibly associated with this cell, you will not be able to call
   * some functions of this class if they need information about degrees of
   * freedom. These functions are, above all, the
   * <tt>get_function_value/gradients/hessians</tt> functions. If you want to
   * call these functions, you have to call the @p reinit variants that take
   * iterators into DoFHandler or other DoF handler type objects.
   */
  void reinit (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
               const unsigned int                                         face_no);

  /**
   * Return a reference to this very object.
   *
   * Though it seems that it is not very useful, this function is there to
   * provide capability to the hpFEValues class, in which case it provides the
   * FEValues object for the present cell (remember that for hp finite
   * elements, the actual FE object used may change from cell to cell, so we
   * also need different FEValues objects for different cells; once you
   * reinitialize the hpFEValues object for a specific cell, it retrieves the
   * FEValues object for the FE on that cell and returns it through a function
   * of the same name as this one; this function here therefore only provides
   * the same interface so that one can templatize on FEValues/hpFEValues).
   */
  const FEFaceValues<dim,spacedim> &get_present_fe_values () const;
private:

  /**
   * Do work common to the two constructors.
   */
  void initialize (const UpdateFlags update_flags);

  /**
   * The reinit() functions do only that part of the work that requires
   * knowledge of the type of iterator. After setting present_cell(), they
   * pass on to this function, which does the real work, and which is
   * independent of the actual type of the cell iterator.
   */
  void do_reinit (const unsigned int face_no);
};


/**
 * Finite element evaluated in quadrature points on a face.
 *
 * This class adds the functionality of FEFaceValuesBase to FEValues; see
 * there for more documentation.
 *
 * This class is used for faces lying on a refinement edge. In this case, the
 * neighboring cell is refined. To be able to compute differences between
 * interior and exterior function values, the refinement of the neighboring
 * cell must be simulated on this cell. This is achieved by applying a
 * quadrature rule that simulates the refinement. The resulting data fields
 * are split up to reflect the refinement structure of the neighbor: a subface
 * number corresponds to the number of the child of the neighboring face.
 *
 * @ingroup feaccess
 * @author Wolfgang Bangerth, 1998, Guido Kanschat, 2000, 2001
 */
template <int dim, int spacedim=dim>
class FESubfaceValues : public FEFaceValuesBase<dim,spacedim>
{
public:
  /**
   * Dimension in which this object operates.
   */
  static const unsigned int dimension = dim;

  /**
   * Dimension of the space in which this object operates.
   */
  static const unsigned int space_dimension = spacedim;

  /**
   * Dimension of the object over which we integrate. For the present class,
   * this is equal to <code>dim-1</code>.
   */
  static const unsigned int integral_dimension = dim-1;

  /**
   * Constructor. Gets cell independent data from mapping and finite element
   * objects, matching the quadrature rule and update flags.
   */
  FESubfaceValues (const Mapping<dim,spacedim>       &mapping,
                   const FiniteElement<dim,spacedim> &fe,
                   const Quadrature<dim-1>  &face_quadrature,
                   const UpdateFlags         update_flags);

  /**
   * Constructor. Uses MappingQ1 implicitly.
   */
  FESubfaceValues (const FiniteElement<dim,spacedim> &fe,
                   const Quadrature<dim-1>  &face_quadrature,
                   const UpdateFlags         update_flags);

  /**
   * Reinitialize the gradients, Jacobi determinants, etc for the given cell
   * of type "iterator into a DoFHandler object", and the finite element
   * associated with this object. It is assumed that the finite element used
   * by the given cell is also the one used by this FESubfaceValues object.
   */
  template <class DH, bool level_dof_access>
  void reinit (const TriaIterator<DoFCellAccessor<DH,level_dof_access> > cell,
               const unsigned int                    face_no,
               const unsigned int                    subface_no);

  /**
   * Reinitialize the gradients, Jacobi determinants, etc for the given
   * subface on given cell of type "iterator into a Triangulation object", and
   * the given finite element. Since iterators into triangulation alone only
   * convey information about the geometry of a cell, but not about degrees of
   * freedom possibly associated with this cell, you will not be able to call
   * some functions of this class if they need information about degrees of
   * freedom. These functions are, above all, the
   * <tt>get_function_value/gradients/hessians</tt> functions. If you want to
   * call these functions, you have to call the @p reinit variants that take
   * iterators into DoFHandler or other DoF handler type objects.
   */
  void reinit (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
               const unsigned int                    face_no,
               const unsigned int                    subface_no);

  /**
   * Return a reference to this very object.
   *
   * Though it seems that it is not very useful, this function is there to
   * provide capability to the hpFEValues class, in which case it provides the
   * FEValues object for the present cell (remember that for hp finite
   * elements, the actual FE object used may change from cell to cell, so we
   * also need different FEValues objects for different cells; once you
   * reinitialize the hpFEValues object for a specific cell, it retrieves the
   * FEValues object for the FE on that cell and returns it through a function
   * of the same name as this one; this function here therefore only provides
   * the same interface so that one can templatize on FEValues/hpFEValues).
   */
  const FESubfaceValues<dim,spacedim> &get_present_fe_values () const;

  /**
   * @todo Document this
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcReinitCalledWithBoundaryFace);

  /**
   * @todo Document this
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcFaceHasNoSubfaces);

private:

  /**
   * Do work common to the two constructors.
   */
  void initialize (const UpdateFlags update_flags);

  /**
   * The reinit() functions do only that part of the work that requires
   * knowledge of the type of iterator. After setting present_cell(), they
   * pass on to this function, which does the real work, and which is
   * independent of the actual type of the cell iterator.
   */
  void do_reinit (const unsigned int face_no,
                  const unsigned int subface_no);
};


#ifndef DOXYGEN


/*------------------------ Inline functions: namespace FEValuesViews --------*/

namespace FEValuesViews
{
  template <int dim, int spacedim>
  inline
  typename Scalar<dim,spacedim>::value_type
  Scalar<dim,spacedim>::value (const unsigned int shape_function,
                               const unsigned int q_point) const
  {
    typedef FEValuesBase<dim,spacedim> FVB;
    Assert (shape_function < fe_values.fe->dofs_per_cell,
            ExcIndexRange (shape_function, 0, fe_values.fe->dofs_per_cell));
    Assert (fe_values.update_flags & update_values,
            typename FVB::ExcAccessToUninitializedField("update_values"));

    // an adaptation of the FEValuesBase::shape_value_component function
    // except that here we know the component as fixed and we have
    // pre-computed and cached a bunch of information. See the comments there.
    if (shape_function_data[shape_function].is_nonzero_shape_function_component)
      return fe_values.shape_values(shape_function_data[shape_function]
                                    .row_index,
                                    q_point);
    else
      return 0;
  }




  template <int dim, int spacedim>
  inline
  typename Scalar<dim,spacedim>::gradient_type
  Scalar<dim,spacedim>::gradient (const unsigned int shape_function,
                                  const unsigned int q_point) const
  {
    typedef FEValuesBase<dim,spacedim> FVB;
    Assert (shape_function < fe_values.fe->dofs_per_cell,
            ExcIndexRange (shape_function, 0, fe_values.fe->dofs_per_cell));
    Assert (fe_values.update_flags & update_gradients,
            typename FVB::ExcAccessToUninitializedField("update_gradients"));

    // an adaptation of the
    // FEValuesBase::shape_grad_component
    // function except that here we know the
    // component as fixed and we have
    // pre-computed and cached a bunch of
    // information. See the comments there.
    if (shape_function_data[shape_function].is_nonzero_shape_function_component)
      return fe_values.shape_gradients[shape_function_data[shape_function]
                                       .row_index][q_point];
    else
      return gradient_type();
  }



  template <int dim, int spacedim>
  inline
  typename Scalar<dim,spacedim>::hessian_type
  Scalar<dim,spacedim>::hessian (const unsigned int shape_function,
                                 const unsigned int q_point) const
  {
    typedef FEValuesBase<dim,spacedim> FVB;
    Assert (shape_function < fe_values.fe->dofs_per_cell,
            ExcIndexRange (shape_function, 0, fe_values.fe->dofs_per_cell));
    Assert (fe_values.update_flags & update_hessians,
            typename FVB::ExcAccessToUninitializedField("update_hessians"));

    // an adaptation of the
    // FEValuesBase::shape_grad_component
    // function except that here we know the
    // component as fixed and we have
    // pre-computed and cached a bunch of
    // information. See the comments there.
    if (shape_function_data[shape_function].is_nonzero_shape_function_component)
      return fe_values.shape_hessians[shape_function_data[shape_function].row_index][q_point];
    else
      return hessian_type();
  }



  template <int dim, int spacedim>
  inline
  typename Vector<dim,spacedim>::value_type
  Vector<dim,spacedim>::value (const unsigned int shape_function,
                               const unsigned int q_point) const
  {
    typedef FEValuesBase<dim,spacedim> FVB;
    Assert (shape_function < fe_values.fe->dofs_per_cell,
            ExcIndexRange (shape_function, 0, fe_values.fe->dofs_per_cell));
    Assert (fe_values.update_flags & update_values,
            typename FVB::ExcAccessToUninitializedField("update_values"));

    // same as for the scalar case except
    // that we have one more index
    const int snc = shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return value_type();
    else if (snc != -1)
      {
        value_type return_value;
        return_value[shape_function_data[shape_function].single_nonzero_component_index]
          = fe_values.shape_values(snc,q_point);
        return return_value;
      }
    else
      {
        value_type return_value;
        for (unsigned int d=0; d<dim; ++d)
          if (shape_function_data[shape_function].is_nonzero_shape_function_component[d])
            return_value[d]
              = fe_values.shape_values(shape_function_data[shape_function].row_index[d],q_point);

        return return_value;
      }
  }



  template <int dim, int spacedim>
  inline
  typename Vector<dim,spacedim>::gradient_type
  Vector<dim,spacedim>::gradient (const unsigned int shape_function,
                                  const unsigned int q_point) const
  {
    typedef FEValuesBase<dim,spacedim> FVB;
    Assert (shape_function < fe_values.fe->dofs_per_cell,
            ExcIndexRange (shape_function, 0, fe_values.fe->dofs_per_cell));
    Assert (fe_values.update_flags & update_gradients,
            typename FVB::ExcAccessToUninitializedField("update_gradients"));

    // same as for the scalar case except
    // that we have one more index
    const int snc = shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return gradient_type();
    else if (snc != -1)
      {
        gradient_type return_value;
        return_value[shape_function_data[shape_function].single_nonzero_component_index]
          = fe_values.shape_gradients[snc][q_point];
        return return_value;
      }
    else
      {
        gradient_type return_value;
        for (unsigned int d=0; d<dim; ++d)
          if (shape_function_data[shape_function].is_nonzero_shape_function_component[d])
            return_value[d]
              = fe_values.shape_gradients[shape_function_data[shape_function].row_index[d]][q_point];

        return return_value;
      }
  }



  template <int dim, int spacedim>
  inline
  typename Vector<dim,spacedim>::divergence_type
  Vector<dim,spacedim>::divergence (const unsigned int shape_function,
                                    const unsigned int q_point) const
  {
    // this function works like in
    // the case above
    typedef FEValuesBase<dim,spacedim> FVB;
    Assert (shape_function < fe_values.fe->dofs_per_cell,
            ExcIndexRange (shape_function, 0, fe_values.fe->dofs_per_cell));
    Assert (fe_values.update_flags & update_gradients,
            typename FVB::ExcAccessToUninitializedField("update_gradients"));

    // same as for the scalar case except
    // that we have one more index
    const int snc = shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return divergence_type();
    else if (snc != -1)
      return
        fe_values.shape_gradients[snc][q_point][shape_function_data[shape_function].single_nonzero_component_index];
    else
      {
        divergence_type return_value = 0;
        for (unsigned int d=0; d<dim; ++d)
          if (shape_function_data[shape_function].is_nonzero_shape_function_component[d])
            return_value
            += fe_values.shape_gradients[shape_function_data[shape_function].row_index[d]][q_point][d];

        return return_value;
      }
  }



  template <int dim, int spacedim>
  inline
  typename Vector<dim,spacedim>::curl_type
  Vector<dim,spacedim>::curl (const unsigned int shape_function, const unsigned int q_point) const
  {
    // this function works like in the case above
    typedef FEValuesBase<dim,spacedim> FVB;

    Assert (shape_function < fe_values.fe->dofs_per_cell,
            ExcIndexRange (shape_function, 0, fe_values.fe->dofs_per_cell));
    Assert (fe_values.update_flags & update_gradients,
            typename FVB::ExcAccessToUninitializedField("update_gradients"));
    // same as for the scalar case except that we have one more index
    const int snc = shape_function_data[shape_function].single_nonzero_component;

    if (snc == -2)
      return curl_type ();

    else
      switch (dim)
        {
        case 1:
        {
          Assert (false, ExcMessage("Computing the curl in 1d is not a useful operation"));
          return curl_type ();
        }

        case 2:
        {
          if (snc != -1)
            {
              curl_type return_value;

              // the single
              // nonzero component
              // can only be zero
              // or one in 2d
              if (shape_function_data[shape_function].single_nonzero_component_index == 0)
                return_value[0] = -1.0 * fe_values.shape_gradients[snc][q_point][1];
              else
                return_value[0] = fe_values.shape_gradients[snc][q_point][0];

              return return_value;
            }

          else
            {
              curl_type return_value;

              return_value[0] = 0.0;

              if (shape_function_data[shape_function].is_nonzero_shape_function_component[0])
                return_value[0]
                -= fe_values.shape_gradients[shape_function_data[shape_function].row_index[0]][q_point][1];

              if (shape_function_data[shape_function].is_nonzero_shape_function_component[1])
                return_value[0]
                += fe_values.shape_gradients[shape_function_data[shape_function].row_index[1]][q_point][0];

              return return_value;
            }
        }

        case 3:
        {
          if (snc != -1)
            {
              curl_type return_value;

              switch (shape_function_data[shape_function].single_nonzero_component_index)
                {
                case 0:
                {
                  return_value[0] = 0;
                  return_value[1] = fe_values.shape_gradients[snc][q_point][2];
                  return_value[2] = -1.0 * fe_values.shape_gradients[snc][q_point][1];
                  return return_value;
                }

                case 1:
                {
                  return_value[0] = -1.0 * fe_values.shape_gradients[snc][q_point][2];
                  return_value[1] = 0;
                  return_value[2] = fe_values.shape_gradients[snc][q_point][0];
                  return return_value;
                }

                default:
                {
                  return_value[0] = fe_values.shape_gradients[snc][q_point][1];
                  return_value[1] = -1.0 * fe_values.shape_gradients[snc][q_point][0];
                  return_value[2] = 0;
                  return return_value;
                }
                }
            }

          else
            {
              curl_type return_value;

              for (unsigned int i = 0; i < dim; ++i)
                return_value[i] = 0.0;

              if (shape_function_data[shape_function].is_nonzero_shape_function_component[0])
                {
                  return_value[1]
                  += fe_values.shape_gradients[shape_function_data[shape_function].row_index[0]][q_point][2];
                  return_value[2]
                  -= fe_values.shape_gradients[shape_function_data[shape_function].row_index[0]][q_point][1];
                }

              if (shape_function_data[shape_function].is_nonzero_shape_function_component[1])
                {
                  return_value[0]
                  -= fe_values.shape_gradients[shape_function_data[shape_function].row_index[1]][q_point][2];
                  return_value[2]
                  += fe_values.shape_gradients[shape_function_data[shape_function].row_index[1]][q_point][0];
                }

              if (shape_function_data[shape_function].is_nonzero_shape_function_component[2])
                {
                  return_value[0]
                  += fe_values.shape_gradients[shape_function_data[shape_function].row_index[2]][q_point][1];
                  return_value[1]
                  -= fe_values.shape_gradients[shape_function_data[shape_function].row_index[2]][q_point][0];
                }

              return return_value;
            }
        }
        }
    // should not end up here
    Assert (false, ExcInternalError());
    return curl_type();
  }

  template <int dim, int spacedim>
  inline
  typename Vector<dim,spacedim>::hessian_type
  Vector<dim,spacedim>::hessian (const unsigned int shape_function,
                                 const unsigned int q_point) const
  {
    // this function works like in
    // the case above
    typedef FEValuesBase<dim,spacedim> FVB;
    Assert (shape_function < fe_values.fe->dofs_per_cell,
            ExcIndexRange (shape_function, 0, fe_values.fe->dofs_per_cell));
    Assert (fe_values.update_flags & update_hessians,
            typename FVB::ExcAccessToUninitializedField("update_hessians"));

    // same as for the scalar case except
    // that we have one more index
    const int snc = shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return hessian_type();
    else if (snc != -1)
      {
        hessian_type return_value;
        return_value[shape_function_data[shape_function].single_nonzero_component_index]
          = fe_values.shape_hessians[snc][q_point];
        return return_value;
      }
    else
      {
        hessian_type return_value;
        for (unsigned int d=0; d<dim; ++d)
          if (shape_function_data[shape_function].is_nonzero_shape_function_component[d])
            return_value[d]
              = fe_values.shape_hessians[shape_function_data[shape_function].row_index[d]][q_point];

        return return_value;
      }
  }


  namespace
  {
    /**
     * Return the symmetrized version of a
     * tensor whose n'th row equals the
     * second argument, with all other rows
     * equal to zero.
     */
    inline
    dealii::SymmetricTensor<2,1>
    symmetrize_single_row (const unsigned int n,
                           const dealii::Tensor<1,1> &t)
    {
      Assert (n < 1, ExcIndexRange (n, 0, 1));
      (void)n; // removes -Wunused-parameter warning in optimized mode

      const double array[1] = { t[0] };
      return dealii::SymmetricTensor<2,1>(array);
    }


    inline
    dealii::SymmetricTensor<2,2>
    symmetrize_single_row (const unsigned int n,
                           const dealii::Tensor<1,2> &t)
    {
      switch (n)
        {
        case 0:
        {
          const double array[3] = { t[0], 0, t[1]/2 };
          return dealii::SymmetricTensor<2,2>(array);
        }
        case 1:
        {
          const double array[3] = { 0, t[1], t[0]/2 };
          return dealii::SymmetricTensor<2,2>(array);
        }
        default:
        {
          Assert (false, ExcIndexRange (n, 0, 2));
          return dealii::SymmetricTensor<2,2>();
        }
        }
    }


    inline
    dealii::SymmetricTensor<2,3>
    symmetrize_single_row (const unsigned int n,
                           const dealii::Tensor<1,3> &t)
    {
      switch (n)
        {
        case 0:
        {
          const double array[6] = { t[0], 0, 0, t[1]/2, t[2]/2, 0 };
          return dealii::SymmetricTensor<2,3>(array);
        }
        case 1:
        {
          const double array[6] = { 0, t[1], 0, t[0]/2, 0, t[2]/2 };
          return dealii::SymmetricTensor<2,3>(array);
        }
        case 2:
        {
          const double array[6] = { 0, 0, t[2], 0, t[0]/2, t[1]/2 };
          return dealii::SymmetricTensor<2,3>(array);
        }
        default:
        {
          Assert (false, ExcIndexRange (n, 0, 3));
          return dealii::SymmetricTensor<2,3>();
        }
        }
    }
  }


  template <int dim, int spacedim>
  inline
  typename Vector<dim,spacedim>::symmetric_gradient_type
  Vector<dim,spacedim>::symmetric_gradient (const unsigned int shape_function,
                                            const unsigned int q_point) const
  {
    typedef FEValuesBase<dim,spacedim> FVB;
    Assert (shape_function < fe_values.fe->dofs_per_cell,
            ExcIndexRange (shape_function, 0, fe_values.fe->dofs_per_cell));
    Assert (fe_values.update_flags & update_gradients,
            typename FVB::ExcAccessToUninitializedField("update_gradients"));

    // same as for the scalar case except
    // that we have one more index
    const int snc = shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return symmetric_gradient_type();
    else if (snc != -1)
      return symmetrize_single_row (shape_function_data[shape_function].single_nonzero_component_index,
                                    fe_values.shape_gradients[snc][q_point]);
    else
      {
        gradient_type return_value;
        for (unsigned int d=0; d<dim; ++d)
          if (shape_function_data[shape_function].is_nonzero_shape_function_component[d])
            return_value[d]
              = fe_values.shape_gradients[shape_function_data[shape_function].row_index[d]][q_point];

        return symmetrize(return_value);
      }
  }



  template <int dim, int spacedim>
  inline
  typename SymmetricTensor<2, dim, spacedim>::value_type
  SymmetricTensor<2, dim, spacedim>::value (const unsigned int shape_function,
                                            const unsigned int q_point) const
  {
    typedef FEValuesBase<dim,spacedim> FVB;
    Assert (shape_function < fe_values.fe->dofs_per_cell,
            ExcIndexRange (shape_function, 0, fe_values.fe->dofs_per_cell));
    Assert (fe_values.update_flags & update_values,
            typename FVB::ExcAccessToUninitializedField("update_values"));

    // similar to the vector case where we
    // have more then one index and we need
    // to convert between unrolled and
    // component indexing for tensors
    const int snc
      = shape_function_data[shape_function].single_nonzero_component;

    if (snc == -2)
      {
        // shape function is zero for the
        // selected components
        return value_type();

      }
    else if (snc != -1)
      {
        value_type return_value;
        const unsigned int comp =
          shape_function_data[shape_function].single_nonzero_component_index;
        return_value[value_type::unrolled_to_component_indices(comp)]
          = fe_values.shape_values(snc,q_point);
        return return_value;
      }
    else
      {
        value_type return_value;
        for (unsigned int d = 0; d < value_type::n_independent_components; ++d)
          if (shape_function_data[shape_function].is_nonzero_shape_function_component[d])
            return_value[value_type::unrolled_to_component_indices(d)]
              = fe_values.shape_values(shape_function_data[shape_function].row_index[d],q_point);
        return return_value;
      }
  }


  template <int dim, int spacedim>
  inline
  typename SymmetricTensor<2, dim, spacedim>::divergence_type
  SymmetricTensor<2, dim, spacedim>::divergence(const unsigned int shape_function,
                                                const unsigned int q_point) const
  {
    typedef FEValuesBase<dim,spacedim> FVB;
    Assert (shape_function < fe_values.fe->dofs_per_cell,
            ExcIndexRange (shape_function, 0, fe_values.fe->dofs_per_cell));
    Assert (fe_values.update_flags & update_gradients,
            typename FVB::ExcAccessToUninitializedField("update_gradients"));

    const int snc = shape_function_data[shape_function].single_nonzero_component;

    if (snc == -2)
      {
        // shape function is zero for the
        // selected components
        return divergence_type();
      }
    else if (snc != -1)
      {
        // we have a single non-zero component
        // when the symmetric tensor is
        // represented in unrolled form.
        // this implies we potentially have
        // two non-zero components when
        // represented in component form!  we
        // will only have one non-zero entry
        // if the non-zero component lies on
        // the diagonal of the tensor.
        //
        // the divergence of a second-order tensor
        // is a first order tensor.
        //
        // assume the second-order tensor is
        // A with components A_{ij}.  then
        // A_{ij} = A_{ji} and there is only
        // one (if diagonal) or two non-zero
        // entries in the tensorial
        // representation.  define the
        // divergence as:
        // b_i := \dfrac{\partial phi_{ij}}{\partial x_j}.
        // (which is incidentally also
        // b_j := \dfrac{\partial phi_{ij}}{\partial x_i}).
        // In both cases, a sum is implied.
        //
        // Now, we know the nonzero component
        // in unrolled form: it is indicated
        // by 'snc'. we can figure out which
        // tensor components belong to this:
        const unsigned int comp =
          shape_function_data[shape_function].single_nonzero_component_index;
        const unsigned int ii = value_type::unrolled_to_component_indices(comp)[0];
        const unsigned int jj = value_type::unrolled_to_component_indices(comp)[1];

        // given the form of the divergence
        // above, if ii=jj there is only a
        // single nonzero component of the
        // full tensor and the gradient
        // equals
        // b_ii := \dfrac{\partial phi_{ii,ii}}{\partial x_ii}.
        // all other entries of 'b' are zero
        //
        // on the other hand, if ii!=jj, then
        // there are two nonzero entries in
        // the full tensor and
        // b_ii := \dfrac{\partial phi_{ii,jj}}{\partial x_ii}.
        // b_jj := \dfrac{\partial phi_{ii,jj}}{\partial x_jj}.
        // again, all other entries of 'b' are
        // zero
        const dealii::Tensor<1, spacedim> phi_grad = fe_values.shape_gradients[snc][q_point];

        divergence_type return_value;
        return_value[ii] = phi_grad[jj];

        if (ii != jj)
          return_value[jj] = phi_grad[ii];

        return return_value;

      }
    else
      {
        Assert (false, ExcNotImplemented());
        divergence_type return_value;
        return return_value;
      }
  }

  template <int dim, int spacedim>
  inline
  typename Tensor<2, dim, spacedim>::value_type
  Tensor<2, dim, spacedim>::value (const unsigned int shape_function,
                                   const unsigned int q_point) const
  {
    typedef FEValuesBase<dim,spacedim> FVB;
    Assert (shape_function < fe_values.fe->dofs_per_cell,
            ExcIndexRange (shape_function, 0, fe_values.fe->dofs_per_cell));
    Assert (fe_values.update_flags & update_values,
            typename FVB::ExcAccessToUninitializedField("update_values"));

    // similar to the vector case where we
    // have more then one index and we need
    // to convert between unrolled and
    // component indexing for tensors
    const int snc
      = shape_function_data[shape_function].single_nonzero_component;

    if (snc == -2)
      {
        // shape function is zero for the
        // selected components
        return value_type();

      }
    else if (snc != -1)
      {
        value_type return_value;
        const unsigned int comp =
          shape_function_data[shape_function].single_nonzero_component_index;
        const TableIndices<2> indices = dealii::Tensor<2,spacedim>::unrolled_to_component_indices(comp);
        return_value[indices] = fe_values.shape_values(snc,q_point);
        return return_value;
      }
    else
      {
        value_type return_value;
        for (unsigned int d = 0; d < dim*dim; ++d)
          if (shape_function_data[shape_function].is_nonzero_shape_function_component[d])
            {
              const TableIndices<2> indices = dealii::Tensor<2,spacedim>::unrolled_to_component_indices(d);
              return_value[indices]
                = fe_values.shape_values(shape_function_data[shape_function].row_index[d],q_point);
            }
        return return_value;
      }
  }


  template <int dim, int spacedim>
  inline
  typename Tensor<2, dim, spacedim>::divergence_type
  Tensor<2, dim, spacedim>::divergence(const unsigned int shape_function,
                                       const unsigned int q_point) const
  {
    typedef FEValuesBase<dim,spacedim> FVB;
    Assert (shape_function < fe_values.fe->dofs_per_cell,
            ExcIndexRange (shape_function, 0, fe_values.fe->dofs_per_cell));
    Assert (fe_values.update_flags & update_gradients,
            typename FVB::ExcAccessToUninitializedField("update_gradients"));

    const int snc = shape_function_data[shape_function].single_nonzero_component;

    if (snc == -2)
      {
        // shape function is zero for the
        // selected components
        return divergence_type();
      }
    else if (snc != -1)
      {
        // we have a single non-zero component
        // when the tensor is
        // represented in unrolled form.
        //
        // the divergence of a second-order tensor
        // is a first order tensor.
        //
        // assume the second-order tensor is
        // A with components A_{ij}.
        // divergence as:
        // b_j := \dfrac{\partial phi_{ij}}{\partial x_i}.
        //
        // Now, we know the nonzero component
        // in unrolled form: it is indicated
        // by 'snc'. we can figure out which
        // tensor components belong to this:
        const unsigned int comp =
          shape_function_data[shape_function].single_nonzero_component_index;
        const TableIndices<2> indices = dealii::Tensor<2,spacedim>::unrolled_to_component_indices(comp);
        const unsigned int ii = indices[0];
        const unsigned int jj = indices[1];

        const dealii::Tensor<1, spacedim> phi_grad = fe_values.shape_gradients[snc][q_point];

        divergence_type return_value;
        return_value[jj] = phi_grad[ii];

        return return_value;

      }
    else
      {
        Assert (false, ExcNotImplemented());
        divergence_type return_value;
        return return_value;
      }
  }
}



/*------------------------ Inline functions: FEValuesBase ------------------------*/



template <int dim, int spacedim>
inline
const FEValuesViews::Scalar<dim,spacedim> &
FEValuesBase<dim,spacedim>::
operator[] (const FEValuesExtractors::Scalar &scalar) const
{
  Assert (scalar.component < fe_values_views_cache.scalars.size(),
          ExcIndexRange (scalar.component,
                         0, fe_values_views_cache.scalars.size()));

  return fe_values_views_cache.scalars[scalar.component];
}



template <int dim, int spacedim>
inline
const FEValuesViews::Vector<dim,spacedim> &
FEValuesBase<dim,spacedim>::
operator[] (const FEValuesExtractors::Vector &vector) const
{
  Assert (vector.first_vector_component <
          fe_values_views_cache.vectors.size(),
          ExcIndexRange (vector.first_vector_component,
                         0, fe_values_views_cache.vectors.size()));

  return fe_values_views_cache.vectors[vector.first_vector_component];
}

template <int dim, int spacedim>
inline
const FEValuesViews::SymmetricTensor<2,dim,spacedim> &
FEValuesBase<dim,spacedim>::
operator[] (const FEValuesExtractors::SymmetricTensor<2> &tensor) const
{
  Assert (tensor.first_tensor_component <
          fe_values_views_cache.symmetric_second_order_tensors.size(),
          ExcIndexRange (tensor.first_tensor_component,
                         0, fe_values_views_cache.symmetric_second_order_tensors.size()));

  return fe_values_views_cache.symmetric_second_order_tensors[tensor.first_tensor_component];
}

template <int dim, int spacedim>
inline
const FEValuesViews::Tensor<2,dim,spacedim> &
FEValuesBase<dim,spacedim>::
operator[] (const FEValuesExtractors::Tensor<2> &tensor) const
{
  Assert (tensor.first_tensor_component <
          fe_values_views_cache.second_order_tensors.size(),
          ExcIndexRange (tensor.first_tensor_component,
                         0, fe_values_views_cache.second_order_tensors.size()));

  return fe_values_views_cache.second_order_tensors[tensor.first_tensor_component];
}




template <int dim, int spacedim>
inline
const double &
FEValuesBase<dim,spacedim>::shape_value (const unsigned int i,
                                         const unsigned int j) const
{
  Assert (i < fe->dofs_per_cell,
          ExcIndexRange (i, 0, fe->dofs_per_cell));
  Assert (this->update_flags & update_values,
          ExcAccessToUninitializedField("update_values"));
  Assert (fe->is_primitive (i),
          ExcShapeFunctionNotPrimitive(i));

  // if the entire FE is primitive,
  // then we can take a short-cut:
  if (fe->is_primitive())
    return this->shape_values(i,j);
  else
    {
      // otherwise, use the mapping
      // between shape function
      // numbers and rows. note that
      // by the assertions above, we
      // know that this particular
      // shape function is primitive,
      // so we can call
      // system_to_component_index
      const unsigned int
      row = this->shape_function_to_row_table[i * fe->n_components() + fe->system_to_component_index(i).first];
      return this->shape_values(row, j);
    }
}



template <int dim, int spacedim>
inline
double
FEValuesBase<dim,spacedim>::shape_value_component (const unsigned int i,
                                                   const unsigned int j,
                                                   const unsigned int component) const
{
  Assert (i < fe->dofs_per_cell,
          ExcIndexRange (i, 0, fe->dofs_per_cell));
  Assert (this->update_flags & update_values,
          ExcAccessToUninitializedField("update_values"));
  Assert (component < fe->n_components(),
          ExcIndexRange(component, 0, fe->n_components()));

  // check whether the shape function
  // is non-zero at all within
  // this component:
  if (fe->get_nonzero_components(i)[component] == false)
    return 0;

  // look up the right row in the
  // table and take the data from
  // there
  const unsigned int
  row = this->shape_function_to_row_table[i * fe->n_components() + component];
  return this->shape_values(row, j);
}



template <int dim, int spacedim>
inline
const Tensor<1,spacedim> &
FEValuesBase<dim,spacedim>::shape_grad (const unsigned int i,
                                        const unsigned int j) const
{
  Assert (i < fe->dofs_per_cell,
          ExcIndexRange (i, 0, fe->dofs_per_cell));
  Assert (this->update_flags & update_gradients,
          ExcAccessToUninitializedField("update_gradients"));
  Assert (fe->is_primitive (i),
          ExcShapeFunctionNotPrimitive(i));
  Assert (i<this->shape_gradients.size(),
          ExcIndexRange (i, 0, this->shape_gradients.size()));
  Assert (j<this->shape_gradients[0].size(),
          ExcIndexRange (j, 0, this->shape_gradients[0].size()));

  // if the entire FE is primitive,
  // then we can take a short-cut:
  if (fe->is_primitive())
    return this->shape_gradients[i][j];
  else
    {
      // otherwise, use the mapping
      // between shape function
      // numbers and rows. note that
      // by the assertions above, we
      // know that this particular
      // shape function is primitive,
      // so we can call
      // system_to_component_index
      const unsigned int
      row = this->shape_function_to_row_table[i * fe->n_components() + fe->system_to_component_index(i).first];
      return this->shape_gradients[row][j];
    }
}



template <int dim, int spacedim>
inline
Tensor<1,spacedim>
FEValuesBase<dim,spacedim>::shape_grad_component (const unsigned int i,
                                                  const unsigned int j,
                                                  const unsigned int component) const
{
  Assert (i < fe->dofs_per_cell,
          ExcIndexRange (i, 0, fe->dofs_per_cell));
  Assert (this->update_flags & update_gradients,
          ExcAccessToUninitializedField("update_gradients"));
  Assert (component < fe->n_components(),
          ExcIndexRange(component, 0, fe->n_components()));

  // check whether the shape function
  // is non-zero at all within
  // this component:
  if (fe->get_nonzero_components(i)[component] == false)
    return Tensor<1,spacedim>();

  // look up the right row in the
  // table and take the data from
  // there
  const unsigned int
  row = this->shape_function_to_row_table[i * fe->n_components() + component];
  return this->shape_gradients[row][j];
}



template <int dim, int spacedim>
inline
const Tensor<2,spacedim> &
FEValuesBase<dim,spacedim>::shape_hessian (const unsigned int i,
                                           const unsigned int j) const
{
  Assert (i < fe->dofs_per_cell,
          ExcIndexRange (i, 0, fe->dofs_per_cell));
  Assert (this->update_flags & update_hessians,
          ExcAccessToUninitializedField("update_hessians"));
  Assert (fe->is_primitive (i),
          ExcShapeFunctionNotPrimitive(i));
  Assert (i<this->shape_hessians.size(),
          ExcIndexRange (i, 0, this->shape_hessians.size()));
  Assert (j<this->shape_hessians[0].size(),
          ExcIndexRange (j, 0, this->shape_hessians[0].size()));

  // if the entire FE is primitive,
  // then we can take a short-cut:
  if (fe->is_primitive())
    return this->shape_hessians[i][j];
  else
    {
      // otherwise, use the mapping
      // between shape function
      // numbers and rows. note that
      // by the assertions above, we
      // know that this particular
      // shape function is primitive,
      // so we can call
      // system_to_component_index
      const unsigned int
      row = this->shape_function_to_row_table[i * fe->n_components() + fe->system_to_component_index(i).first];
      return this->shape_hessians[row][j];
    }
}



template <int dim, int spacedim>
inline
const Tensor<2,spacedim> &
FEValuesBase<dim,spacedim>::shape_2nd_derivative (const unsigned int i,
                                                  const unsigned int j) const
{
  return shape_hessian(i,j);
}



template <int dim, int spacedim>
inline
Tensor<2,spacedim>
FEValuesBase<dim,spacedim>::shape_hessian_component (const unsigned int i,
                                                     const unsigned int j,
                                                     const unsigned int component) const
{
  Assert (i < fe->dofs_per_cell,
          ExcIndexRange (i, 0, fe->dofs_per_cell));
  Assert (this->update_flags & update_hessians,
          ExcAccessToUninitializedField("update_hessians"));
  Assert (component < fe->n_components(),
          ExcIndexRange(component, 0, fe->n_components()));

  // check whether the shape function
  // is non-zero at all within
  // this component:
  if (fe->get_nonzero_components(i)[component] == false)
    return Tensor<2,spacedim>();

  // look up the right row in the
  // table and take the data from
  // there
  const unsigned int
  row = this->shape_function_to_row_table[i * fe->n_components() + component];
  return this->shape_hessians[row][j];
}



template <int dim, int spacedim>
inline
Tensor<2,spacedim>
FEValuesBase<dim,spacedim>::shape_2nd_derivative_component (const unsigned int i,
                                                            const unsigned int j,
                                                            const unsigned int component) const
{
  return shape_hessian_component(i,j,component);
}



template <int dim, int spacedim>
inline
const FiniteElement<dim,spacedim> &
FEValuesBase<dim,spacedim>::get_fe () const
{
  return *fe;
}


template <int dim, int spacedim>
inline
const Mapping<dim,spacedim> &
FEValuesBase<dim,spacedim>::get_mapping () const
{
  return *mapping;
}



template <int dim, int spacedim>
inline
UpdateFlags
FEValuesBase<dim,spacedim>::get_update_flags () const
{
  return this->update_flags;
}



template <int dim, int spacedim>
inline
const std::vector<Point<spacedim> > &
FEValuesBase<dim,spacedim>::get_quadrature_points () const
{
  Assert (this->update_flags & update_quadrature_points,
          ExcAccessToUninitializedField("update_quadrature_points"));
  return this->quadrature_points;
}



template <int dim, int spacedim>
inline
const std::vector<double> &
FEValuesBase<dim,spacedim>::get_JxW_values () const
{
  Assert (this->update_flags & update_JxW_values,
          ExcAccessToUninitializedField("update_JxW_values"));
  return this->JxW_values;
}



template <int dim, int spacedim>
inline
const std::vector<DerivativeForm<1,dim,spacedim> > &
FEValuesBase<dim,spacedim>::get_jacobians () const
{
  Assert (this->update_flags & update_jacobians,
          ExcAccessToUninitializedField("update_jacobians"));
  return this->jacobians;
}



template <int dim, int spacedim>
inline
const std::vector<DerivativeForm<2,dim,spacedim> > &
FEValuesBase<dim,spacedim>::get_jacobian_grads () const
{
  Assert (this->update_flags & update_jacobian_grads,
          ExcAccessToUninitializedField("update_jacobians_grads"));
  return this->jacobian_grads;
}



template <int dim, int spacedim>
inline
const std::vector<DerivativeForm<1,spacedim,dim> > &
FEValuesBase<dim,spacedim>::get_inverse_jacobians () const
{
  Assert (this->update_flags & update_inverse_jacobians,
          ExcAccessToUninitializedField("update_inverse_jacobians"));
  return this->inverse_jacobians;
}



template <int dim, int spacedim>
inline
const Point<spacedim> &
FEValuesBase<dim,spacedim>::quadrature_point (const unsigned int i) const
{
  Assert (this->update_flags & update_quadrature_points,
          ExcAccessToUninitializedField("update_quadrature_points"));
  Assert (i<this->quadrature_points.size(), ExcIndexRange(i, 0, this->quadrature_points.size()));

  return this->quadrature_points[i];
}




template <int dim, int spacedim>
inline
double
FEValuesBase<dim,spacedim>::JxW (const unsigned int i) const
{
  Assert (this->update_flags & update_JxW_values,
          ExcAccessToUninitializedField("update_JxW_values"));
  Assert (i<this->JxW_values.size(), ExcIndexRange(i, 0, this->JxW_values.size()));

  return this->JxW_values[i];
}



template <int dim, int spacedim>
inline
const DerivativeForm<1,dim,spacedim> &
FEValuesBase<dim,spacedim>::jacobian (const unsigned int i) const
{
  Assert (this->update_flags & update_jacobians,
          ExcAccessToUninitializedField("update_jacobians"));
  Assert (i<this->jacobians.size(), ExcIndexRange(i, 0, this->jacobians.size()));

  return this->jacobians[i];
}



template <int dim, int spacedim>
inline
const DerivativeForm<2,dim,spacedim> &
FEValuesBase<dim,spacedim>::jacobian_grad (const unsigned int i) const
{
  Assert (this->update_flags & update_jacobian_grads,
          ExcAccessToUninitializedField("update_jacobians_grads"));
  Assert (i<this->jacobian_grads.size(), ExcIndexRange(i, 0, this->jacobian_grads.size()));

  return this->jacobian_grads[i];
}



template <int dim, int spacedim>
inline
const DerivativeForm<1,spacedim,dim> &
FEValuesBase<dim,spacedim>::inverse_jacobian (const unsigned int i) const
{
  Assert (this->update_flags & update_inverse_jacobians,
          ExcAccessToUninitializedField("update_inverse_jacobians"));
  Assert (i<this->inverse_jacobians.size(), ExcIndexRange(i, 0, this->inverse_jacobians.size()));

  return this->inverse_jacobians[i];
}


template <int dim, int spacedim>
template <class InputVector>
inline
void
FEValuesBase<dim,spacedim>::get_function_grads (const InputVector           &fe_function,
                                                std::vector<Tensor<1,spacedim> > &gradients) const
{
  get_function_gradients(fe_function, gradients);
}



template <int dim, int spacedim>
template <class InputVector>
inline
void
FEValuesBase<dim,spacedim>::get_function_grads (
  const InputVector &fe_function,
  const VectorSlice<const std::vector<types::global_dof_index> > &indices,
  std::vector<Tensor<1,spacedim> > &values) const
{
  get_function_gradients(fe_function, indices, values);
}



template <int dim, int spacedim>
template <class InputVector>
inline
void
FEValuesBase<dim,spacedim>::
get_function_grads (const InputVector                         &fe_function,
                    std::vector<std::vector<Tensor<1,spacedim> > > &gradients) const
{
  get_function_gradients(fe_function, gradients);
}



template <int dim, int spacedim>
template <class InputVector>
inline
void
FEValuesBase<dim,spacedim>::get_function_grads (
  const InputVector &fe_function,
  const VectorSlice<const std::vector<types::global_dof_index> > &indices,
  std::vector<std::vector<Tensor<1,spacedim> > > &values,
  bool q_points_fastest) const
{
  get_function_gradients(fe_function, indices, values, q_points_fastest);
}



template <int dim, int spacedim>
template <class InputVector>
inline
void
FEValuesBase<dim,spacedim>::
get_function_2nd_derivatives (const InputVector           &fe_function,
                              std::vector<Tensor<2,spacedim> > &hessians) const
{
  get_function_hessians(fe_function, hessians);
}



template <int dim, int spacedim>
template <class InputVector>
inline
void
FEValuesBase<dim,spacedim>::
get_function_2nd_derivatives (const InputVector                         &fe_function,
                              std::vector<std::vector<Tensor<2,spacedim> > > &hessians,
                              bool quadrature_points_fastest) const
{
  get_function_hessians(fe_function, hessians, quadrature_points_fastest);
}



template <int dim, int spacedim>
inline
const Point<spacedim> &
FEValuesBase<dim,spacedim>::normal_vector (const unsigned int i) const
{
  typedef FEValuesBase<dim,spacedim> FVB;
  Assert (this->update_flags & update_normal_vectors,
          typename FVB::ExcAccessToUninitializedField("update_normal_vectors"));
  Assert (i<this->normal_vectors.size(),
          ExcIndexRange(i, 0, this->normal_vectors.size()));

  return this->normal_vectors[i];
}



template <int dim, int spacedim>
inline
const Point<spacedim> &
FEValuesBase<dim,spacedim>::cell_normal_vector (const unsigned int i) const
{
  return this->normal_vector(i);
}




/*------------------------ Inline functions: FEValues ----------------------------*/


template <int dim, int spacedim>
inline
const Quadrature<dim> &
FEValues<dim,spacedim>::get_quadrature () const
{
  return quadrature;
}



template <int dim, int spacedim>
inline
const FEValues<dim,spacedim> &
FEValues<dim,spacedim>::get_present_fe_values () const
{
  return *this;
}


/*------------------------ Inline functions: FEFaceValuesBase --------------------*/


template <int dim, int spacedim>
inline
unsigned int
FEFaceValuesBase<dim,spacedim>::get_face_index () const
{
  return present_face_index;
}


/*------------------------ Inline functions: FE*FaceValues --------------------*/

template <int dim, int spacedim>
inline
const Quadrature<dim-1> &
FEFaceValuesBase<dim,spacedim>::get_quadrature () const
{
  return quadrature;
}



template <int dim, int spacedim>
inline
const FEFaceValues<dim,spacedim> &
FEFaceValues<dim,spacedim>::get_present_fe_values () const
{
  return *this;
}



template <int dim, int spacedim>
inline
const FESubfaceValues<dim,spacedim> &
FESubfaceValues<dim,spacedim>::get_present_fe_values () const
{
  return *this;
}



template <int dim, int spacedim>
inline
const Tensor<1,spacedim> &
FEFaceValuesBase<dim,spacedim>::boundary_form (const unsigned int i) const
{
  typedef FEValuesBase<dim,spacedim> FVB;
  Assert (i<this->boundary_forms.size(),
          ExcIndexRange(i, 0, this->boundary_forms.size()));
  Assert (this->update_flags & update_boundary_forms,
          typename FVB::ExcAccessToUninitializedField("update_boundary_forms"));

  return this->boundary_forms[i];
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

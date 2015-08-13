// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2015 by the deal.II authors
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

#ifndef dealii__fe_poly_h
#define dealii__fe_poly_h


#include <deal.II/fe/fe.h>
#include <deal.II/base/quadrature.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup febase */
/*@{*/

/**
 * This class gives a unified framework for the implementation of
 * FiniteElement classes based on polynomial spaces like the
 * TensorProductPolynomials or PolynomialSpace classes.
 *
 * Every class conforming to the following interface can be used as template
 * parameter POLY.
 *
 * @code
 * static const unsigned int dimension;
 *
 * double compute_value (const unsigned int i,
 *                       const Point<dim> &p) const;
 *
 * Tensor<1,dim> compute_grad (const unsigned int i,
 *                             const Point<dim> &p) const;
 *
 * Tensor<2,dim> compute_grad_grad (const unsigned int i,
 *                                  const Point<dim> &p) const;
 * @endcode
 * Example classes are TensorProductPolynomials, PolynomialSpace or
 * PolynomialsP.
 *
 * This class is not a fully implemented FiniteElement class. Instead there
 * are several pure virtual functions declared in the FiniteElement and
 * FiniteElement classes which cannot implemented by this class but are left
 * for implementation in derived classes.
 *
 * Furthermore, this class assumes that shape functions of the FiniteElement
 * under consideration do <em>not</em> depend on the actual shape of the cells
 * in real space, i.e. update_once() includes <tt>update_values</tt>. For
 * FiniteElements whose shape functions depend on the cells in real space, the
 * update_once() and update_each() functions must be overloaded.
 *
 * @todo Since nearly all functions for spacedim != dim are specialized, this
 * class needs cleaning up.
 *
 * @author Ralf Hartmann 2004, Guido Kanschat, 2009
 */

template <class POLY, int dim=POLY::dimension, int spacedim=dim>
class FE_Poly : public FiniteElement<dim,spacedim>
{
public:
  /**
   * Constructor.
   */
  FE_Poly (const POLY &poly_space,
           const FiniteElementData<dim> &fe_data,
           const std::vector<bool> &restriction_is_additive_flags,
           const std::vector<ComponentMask> &nonzero_components);

  /**
   * Return the polynomial degree of this finite element, i.e. the value
   * passed to the constructor.
   */
  unsigned int get_degree () const;

  /**
   * Return the numbering of the underlying polynomial space compared to
   * lexicographic ordering of the basis functions. Returns
   * POLY::get_numbering().
   */
  std::vector<unsigned int> get_poly_space_numbering() const;

  /**
   * Return the inverse numbering of the underlying polynomial space. Returns
   * POLY::get_numbering_inverse().
   */
  std::vector<unsigned int> get_poly_space_numbering_inverse() const;

  /**
   * Return the value of the <tt>i</tt>th shape function at the point
   * <tt>p</tt>. See the FiniteElement base class for more information about
   * the semantics of this function.
   */
  virtual double shape_value (const unsigned int i,
                              const Point<dim> &p) const;

  /**
   * Return the value of the <tt>component</tt>th vector component of the
   * <tt>i</tt>th shape function at the point <tt>p</tt>. See the
   * FiniteElement base class for more information about the semantics of this
   * function.
   *
   * Since this element is scalar, the returned value is the same as if the
   * function without the <tt>_component</tt> suffix were called, provided
   * that the specified component is zero.
   */
  virtual double shape_value_component (const unsigned int i,
                                        const Point<dim> &p,
                                        const unsigned int component) const;

  /**
   * Return the gradient of the <tt>i</tt>th shape function at the point
   * <tt>p</tt>. See the FiniteElement base class for more information about
   * the semantics of this function.
   */
  virtual Tensor<1,dim> shape_grad (const unsigned int  i,
                                    const Point<dim>   &p) const;

  /**
   * Return the gradient of the <tt>component</tt>th vector component of the
   * <tt>i</tt>th shape function at the point <tt>p</tt>. See the
   * FiniteElement base class for more information about the semantics of this
   * function.
   *
   * Since this element is scalar, the returned value is the same as if the
   * function without the <tt>_component</tt> suffix were called, provided
   * that the specified component is zero.
   */
  virtual Tensor<1,dim> shape_grad_component (const unsigned int i,
                                              const Point<dim> &p,
                                              const unsigned int component) const;

  /**
   * Return the tensor of second derivatives of the <tt>i</tt>th shape
   * function at point <tt>p</tt> on the unit cell. See the FiniteElement base
   * class for more information about the semantics of this function.
   */
  virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
                                         const Point<dim> &p) const;

  /**
   * Return the second derivative of the <tt>component</tt>th vector component
   * of the <tt>i</tt>th shape function at the point <tt>p</tt>. See the
   * FiniteElement base class for more information about the semantics of this
   * function.
   *
   * Since this element is scalar, the returned value is the same as if the
   * function without the <tt>_component</tt> suffix were called, provided
   * that the specified component is zero.
   */
  virtual Tensor<2,dim> shape_grad_grad_component (const unsigned int i,
                                                   const Point<dim> &p,
                                                   const unsigned int component) const;

protected:
  /**
  NOTE: The following function has its definition inlined into the class declaration
    * because we otherwise run into a compiler error with MS Visual Studio.
    */
  virtual
  typename FiniteElement<dim,spacedim>::InternalDataBase *
  get_data(const UpdateFlags update_flags,
           const Mapping<dim,spacedim> &mapping,
           const Quadrature<dim> &quadrature) const
  {
    // generate a new data object and
    // initialize some fields
    InternalData *data = new InternalData;

    // check what needs to be
    // initialized only once and what
    // on every cell/face/subface we
    // visit
    data->update_once = update_once(update_flags);
    data->update_each = update_each(update_flags);
    data->update_flags = data->update_once | data->update_each;

    const UpdateFlags flags(data->update_flags);
    const unsigned int n_q_points = quadrature.size();

    // some scratch arrays
    std::vector<double> values(0);
    std::vector<Tensor<1, dim> > grads(0);
    std::vector<Tensor<2, dim> > grad_grads(0);

    // initialize fields only if really
    // necessary. otherwise, don't
    // allocate memory
    if (flags & update_values)
      {
        values.resize(this->dofs_per_cell);
        data->shape_values.resize(this->dofs_per_cell,
                                  std::vector<double>(n_q_points));
      }

    if (flags & update_gradients)
      {
        grads.resize(this->dofs_per_cell);
        data->shape_gradients.resize(this->dofs_per_cell,
                                     std::vector<Tensor<1, dim> >(n_q_points));
      }

    // if second derivatives through
    // finite differencing is required,
    // then initialize some objects for
    // that
    if (flags & update_hessians)
      data->initialize_2nd(this, mapping, quadrature);

    // next already fill those fields
    // of which we have information by
    // now. note that the shape
    // gradients are only those on the
    // unit cell, and need to be
    // transformed when visiting an
    // actual cell
    if (flags & (update_values | update_gradients))
      for (unsigned int i = 0; i<n_q_points; ++i)
        {
          poly_space.compute(quadrature.point(i),
                             values, grads, grad_grads);

          if (flags & update_values)
            for (unsigned int k = 0; k<this->dofs_per_cell; ++k)
              data->shape_values[k][i] = values[k];

          if (flags & update_gradients)
            for (unsigned int k = 0; k<this->dofs_per_cell; ++k)
              data->shape_gradients[k][i] = grads[k];
        }
    return data;
  }

  virtual
  void
  fill_fe_values (const Mapping<dim,spacedim>                               &mapping,
                  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                  const Quadrature<dim>                                     &quadrature,
                  const typename Mapping<dim,spacedim>::InternalDataBase    &mapping_internal,
                  const typename FiniteElement<dim,spacedim>::InternalDataBase    &fe_internal,
                  const internal::FEValues::MappingRelatedData<dim,spacedim> &mapping_data,
                  internal::FEValues::FiniteElementRelatedData<dim,spacedim> &output_data,
                  const CellSimilarity::Similarity                           cell_similarity) const;

  virtual
  void
  fill_fe_face_values (const Mapping<dim,spacedim>                               &mapping,
                       const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                       const unsigned int                                         face_no,
                       const Quadrature<dim-1>                                   &quadrature,
                       const typename Mapping<dim,spacedim>::InternalDataBase    &mapping_internal,
                       const typename FiniteElement<dim,spacedim>::InternalDataBase    &fe_internal,
                       const internal::FEValues::MappingRelatedData<dim,spacedim> &mapping_data,
                       internal::FEValues::FiniteElementRelatedData<dim,spacedim> &output_data) const;

  virtual
  void
  fill_fe_subface_values (const Mapping<dim,spacedim>                               &mapping,
                          const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                          const unsigned int                                         face_no,
                          const unsigned int                                         sub_no,
                          const Quadrature<dim-1>                                   &quadrature,
                          const typename Mapping<dim,spacedim>::InternalDataBase    &mapping_internal,
                          const typename FiniteElement<dim,spacedim>::InternalDataBase    &fe_internal,
                          const internal::FEValues::MappingRelatedData<dim,spacedim> &mapping_data,
                          internal::FEValues::FiniteElementRelatedData<dim,spacedim> &output_data) const;

  /**
   * Determine the values that need to be computed on the unit cell to be able
   * to compute all values required by <tt>flags</tt>.
   *
   * For the purpose of this function, refer to the documentation in
   * FiniteElement.
   *
   * This class assumes that shape functions of this FiniteElement do
   * <em>not</em> depend on the actual shape of the cells in real space.
   * Therefore, the effect in this element is as follows: if
   * <tt>update_values</tt> is set in <tt>flags</tt>, copy it to the result.
   * All other flags of the result are cleared, since everything else must be
   * computed for each cell.
   */
  virtual UpdateFlags update_once (const UpdateFlags flags) const;

  /**
   * Determine the values that need to be computed on every cell to be able to
   * compute all values required by <tt>flags</tt>.
   *
   * For the purpose of this function, refer to the documentation in
   * FiniteElement.
   *
   * This class assumes that shape functions of this FiniteElement do
   * <em>not</em> depend on the actual shape of the cells in real space.
   *
   * The effect in this element is as follows:
   * <ul>
   *
   * <li> if <tt>update_gradients</tt> is set, the result will contain
   * <tt>update_gradients</tt> and <tt>update_covariant_transformation</tt>.
   * The latter is required to transform the gradient on the unit cell to the
   * real cell. Remark, that the action required by
   * <tt>update_covariant_transformation</tt> is actually performed by the
   * Mapping object used in conjunction with this finite element.
   *
   * <li> if <tt>update_hessians</tt> is set, the result will contain
   * <tt>update_hessians</tt> and <tt>update_covariant_transformation</tt>.
   * The rationale is the same as above and no higher derivatives of the
   * transformation are required, since we use difference quotients for the
   * actual computation.
   *
   * </ul>
   */
  virtual UpdateFlags update_each (const UpdateFlags flags) const;

  /**
   * Fields of cell-independent data.
   *
   * For information about the general purpose of this class, see the
   * documentation of the base class.
   */
  class InternalData : public FiniteElement<dim,spacedim>::InternalDataBase
  {
  public:
    /**
     * Array with shape function values in quadrature points. There is one row
     * for each shape function, containing values for each quadrature point.
     *
     * In this array, we store the values of the shape function in the
     * quadrature points on the unit cell. Since these values do not change
     * under transformation to the real cell, we only need to copy them over
     * when visiting a concrete cell.
     */
    std::vector<std::vector<double> > shape_values;

    /**
     * Array with shape function gradients in quadrature points. There is one
     * row for each shape function, containing values for each quadrature
     * point.
     *
     * We store the gradients in the quadrature points on the unit cell. We
     * then only have to apply the transformation (which is a matrix-vector
     * multiplication) when visiting an actual cell.
     */
    std::vector<std::vector<Tensor<1,dim> > > shape_gradients;

    /**
     * Array with shape function hessians in quadrature points. There is one
     * row for each shape function, containing values for each quadrature
     * point.
     *
     * We store the hessians in the quadrature points on the unit cell. We
     * then only have to apply the transformation when visiting an actual cell.
     */
    std::vector<std::vector<Tensor<2,dim> > > shape_hessians;

    /**
     * Scratch array to store temporary values during hessian calculations in
     * actual cells.
     */
    mutable std::vector<Tensor<2,dim> > untransformed_shape_hessians;
  };

  /**
   * Correct the hessian in the reference cell by subtracting the term corresponding
   * to the Jacobian gradient for one degree of freedom. The result being given by:
   * @f[
   * \frac{\partial^2 \phi_i}{\partial\hat{x}_J\partial\hat{x}_K}
   * - \frac{\partial \phi_i}{\partial {x}_l}
   * \left( \frac{\partial^2{x}_l}{\partial\hat{x}_J\partial\hat{x}_K} \right).
   * @f]
   * After this correction, the shape hessians are simply a mapping_covariant_gradient
   * transformation.
   */
  void
  correct_untransformed_hessians (VectorSlice< std::vector<Tensor<2, dim> > >                       uncorrected_shape_hessians,
                                  const internal::FEValues::MappingRelatedData<dim,spacedim>       &mapping_data,
                                  const internal::FEValues::FiniteElementRelatedData<dim,spacedim> &fevalues_data,
                                  const unsigned int                                                n_q_points,
                                  const unsigned int                                                dof) const;

  /**
   * The polynomial space. Its type is given by the template parameter POLY.
   */
  POLY poly_space;
};

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif

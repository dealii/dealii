// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_fe_poly_h
#define dealii_fe_poly_h


#include <deal.II/base/config.h>

#include <deal.II/base/quadrature.h>
#include <deal.II/base/scalar_polynomials_base.h>

#include <deal.II/fe/fe.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup febase */
/*@{*/

/**
 * This class gives a unified framework for the implementation of
 * FiniteElement classes based on scalar polynomial spaces like the
 * TensorProductPolynomials or PolynomialSpace classes. This
 * class has a corresponding class for tensor-valued finite
 * elements in the FE_PolyTensor class.
 *
 * Every class that has the following public member variables and
 * functions can be used as template parameter @p PolynomialType.
 *
 * @code
 *  static const unsigned int dimension;
 *
 *  void evaluate (const Point<dim>            &unit_point,
 *                 std::vector<double>         &values,
 *                 std::vector<Tensor<1,dim> > &grads,
 *                 std::vector<Tensor<2,dim> > &grad_grads,
 *                 std::vector<Tensor<3,dim> > &third_derivatives,
 *                 std::vector<Tensor<4,dim> > &fourth_derivatives) const;
 *
 *  double compute_value (const unsigned int i,
 *                        const Point<dim> &p) const;
 *
 *  template <int order>
 *  Tensor<order,dim> compute_derivative (const unsigned int i,
 *                                        const Point<dim> &p) const;
 * @endcode
 * Example classes are TensorProductPolynomials, PolynomialSpace or
 * PolynomialsP.
 *
 * This class is not a fully implemented FiniteElement class. Instead there
 * are several pure virtual functions declared in the FiniteElement and
 * FiniteElement classes which cannot be implemented by this class but are
 * left for implementation in derived classes.
 *
 * @todo Since nearly all functions for spacedim != dim are specialized, this
 * class needs cleaning up.
 *
 * @author Ralf Hartmann 2004, Guido Kanschat, 2009
 */

template <int dim, int spacedim = dim>
class FE_Poly : public FiniteElement<dim, spacedim>
{
public:
  /**
   * Constructor.
   */
  FE_Poly(const ScalarPolynomialsBase<dim> &poly_space,
          const FiniteElementData<dim> &    fe_data,
          const std::vector<bool> &         restriction_is_additive_flags,
          const std::vector<ComponentMask> &nonzero_components);

  /**
   * Copy constructor.
   */
  FE_Poly(const FE_Poly &fe);

  /**
   * Return the polynomial degree of this finite element, i.e. the value
   * passed to the constructor.
   */
  unsigned int
  get_degree() const;

  // for documentation, see the FiniteElement base class
  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

  /**
   * Return the underlying polynomial space.
   */
  const ScalarPolynomialsBase<dim> &
  get_poly_space() const;

  /**
   * Return the numbering of the underlying polynomial space compared to
   * lexicographic ordering of the basis functions. Returns
   * PolynomialType::get_numbering().
   */
  std::vector<unsigned int>
  get_poly_space_numbering() const;

  /**
   * Return the inverse numbering of the underlying polynomial space. Returns
   * PolynomialType::get_numbering_inverse().
   */
  std::vector<unsigned int>
  get_poly_space_numbering_inverse() const;

  /**
   * Return the value of the <tt>i</tt>th shape function at the point
   * <tt>p</tt>. See the FiniteElement base class for more information about
   * the semantics of this function.
   */
  virtual double
  shape_value(const unsigned int i, const Point<dim> &p) const override;

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
  virtual double
  shape_value_component(const unsigned int i,
                        const Point<dim> & p,
                        const unsigned int component) const override;

  /**
   * Return the gradient of the <tt>i</tt>th shape function at the point
   * <tt>p</tt>. See the FiniteElement base class for more information about
   * the semantics of this function.
   */
  virtual Tensor<1, dim>
  shape_grad(const unsigned int i, const Point<dim> &p) const override;

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
  virtual Tensor<1, dim>
  shape_grad_component(const unsigned int i,
                       const Point<dim> & p,
                       const unsigned int component) const override;

  /**
   * Return the tensor of second derivatives of the <tt>i</tt>th shape
   * function at point <tt>p</tt> on the unit cell. See the FiniteElement base
   * class for more information about the semantics of this function.
   */
  virtual Tensor<2, dim>
  shape_grad_grad(const unsigned int i, const Point<dim> &p) const override;

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
  virtual Tensor<2, dim>
  shape_grad_grad_component(const unsigned int i,
                            const Point<dim> & p,
                            const unsigned int component) const override;

  /**
   * Return the tensor of third derivatives of the <tt>i</tt>th shape function
   * at point <tt>p</tt> on the unit cell. See the FiniteElement base class
   * for more information about the semantics of this function.
   */
  virtual Tensor<3, dim>
  shape_3rd_derivative(const unsigned int i,
                       const Point<dim> & p) const override;

  /**
   * Return the third derivative of the <tt>component</tt>th vector component
   * of the <tt>i</tt>th shape function at the point <tt>p</tt>. See the
   * FiniteElement base class for more information about the semantics of this
   * function.
   *
   * Since this element is scalar, the returned value is the same as if the
   * function without the <tt>_component</tt> suffix were called, provided
   * that the specified component is zero.
   */
  virtual Tensor<3, dim>
  shape_3rd_derivative_component(const unsigned int i,
                                 const Point<dim> & p,
                                 const unsigned int component) const override;

  /**
   * Return the tensor of fourth derivatives of the <tt>i</tt>th shape
   * function at point <tt>p</tt> on the unit cell. See the FiniteElement base
   * class for more information about the semantics of this function.
   */
  virtual Tensor<4, dim>
  shape_4th_derivative(const unsigned int i,
                       const Point<dim> & p) const override;

  /**
   * Return the fourth derivative of the <tt>component</tt>th vector component
   * of the <tt>i</tt>th shape function at the point <tt>p</tt>. See the
   * FiniteElement base class for more information about the semantics of this
   * function.
   *
   * Since this element is scalar, the returned value is the same as if the
   * function without the <tt>_component</tt> suffix were called, provided
   * that the specified component is zero.
   */
  virtual Tensor<4, dim>
  shape_4th_derivative_component(const unsigned int i,
                                 const Point<dim> & p,
                                 const unsigned int component) const override;

  /**
   * Return an estimate (in bytes) for the memory consumption of this object.
   */
  virtual std::size_t
  memory_consumption() const override;

protected:
  /*
   * NOTE: The following function has its definition inlined into the class
   * declaration because we otherwise run into a compiler error with MS Visual
   * Studio.
   */


  virtual std::unique_ptr<
    typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_data(
    const UpdateFlags update_flags,
    const Mapping<dim, spacedim> & /*mapping*/,
    const Quadrature<dim> &quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override
  {
    // generate a new data object and
    // initialize some fields
    std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
          data_ptr   = std::make_unique<InternalData>();
    auto &data       = dynamic_cast<InternalData &>(*data_ptr);
    data.update_each = requires_update_flags(update_flags);

    const unsigned int n_q_points = quadrature.size();

    // initialize some scratch arrays. we need them for the underlying
    // polynomial to put the values and derivatives of shape functions
    // to put there, depending on what the user requested
    std::vector<double> values(
      update_flags & update_values ? this->dofs_per_cell : 0);
    std::vector<Tensor<1, dim>> grads(
      update_flags & update_gradients ? this->dofs_per_cell : 0);
    std::vector<Tensor<2, dim>> grad_grads(
      update_flags & update_hessians ? this->dofs_per_cell : 0);
    std::vector<Tensor<3, dim>> third_derivatives(
      update_flags & update_3rd_derivatives ? this->dofs_per_cell : 0);
    std::vector<Tensor<4, dim>>
      fourth_derivatives; // won't be needed, so leave empty

    // now also initialize fields the fields of this class's own
    // temporary storage, depending on what we need for the given
    // update flags.
    //
    // there is one exception from the rule: if we are dealing with
    // cells (i.e., if this function is not called via
    // get_(sub)face_data()), then we can already store things in the
    // final location where FEValues::reinit() later wants to see
    // things. we then don't need the intermediate space. we determine
    // whether we are on a cell by asking whether the number of
    // elements in the output array equals the number of quadrature
    // points (yes, it's a cell) or not (because in that case the
    // number of quadrature points we use here equals the number of
    // quadrature points summed over *all* faces or subfaces, whereas
    // the number of output slots equals the number of quadrature
    // points on only *one* face)
    if ((update_flags & update_values) &&
        !((output_data.shape_values.n_rows() > 0) &&
          (output_data.shape_values.n_cols() == n_q_points)))
      data.shape_values.reinit(this->dofs_per_cell, n_q_points);

    if (update_flags & update_gradients)
      data.shape_gradients.reinit(this->dofs_per_cell, n_q_points);

    if (update_flags & update_hessians)
      data.shape_hessians.reinit(this->dofs_per_cell, n_q_points);

    if (update_flags & update_3rd_derivatives)
      data.shape_3rd_derivatives.reinit(this->dofs_per_cell, n_q_points);

    // next already fill those fields of which we have information by
    // now. note that the shape gradients are only those on the unit
    // cell, and need to be transformed when visiting an actual cell
    if (update_flags & (update_values | update_gradients | update_hessians |
                        update_3rd_derivatives))
      for (unsigned int i = 0; i < n_q_points; ++i)
        {
          poly_space->evaluate(quadrature.point(i),
                               values,
                               grads,
                               grad_grads,
                               third_derivatives,
                               fourth_derivatives);

          // the values of shape functions at quadrature points don't change.
          // consequently, write these values right into the output array if
          // we can, i.e., if the output array has the correct size. this is
          // the case on cells. on faces, we already precompute data on *all*
          // faces and subfaces, but we later on copy only a portion of it
          // into the output object; in that case, copy the data from all
          // faces into the scratch object
          if (update_flags & update_values)
            if (output_data.shape_values.n_rows() > 0)
              {
                if (output_data.shape_values.n_cols() == n_q_points)
                  for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
                    output_data.shape_values[k][i] = values[k];
                else
                  for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
                    data.shape_values[k][i] = values[k];
              }

          // for everything else, derivatives need to be transformed,
          // so we write them into our scratch space and only later
          // copy stuff into where FEValues wants it
          if (update_flags & update_gradients)
            for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
              data.shape_gradients[k][i] = grads[k];

          if (update_flags & update_hessians)
            for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
              data.shape_hessians[k][i] = grad_grads[k];

          if (update_flags & update_3rd_derivatives)
            for (unsigned int k = 0; k < this->dofs_per_cell; ++k)
              data.shape_3rd_derivatives[k][i] = third_derivatives[k];
        }
    return data_ptr;
  }

  virtual void
  fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim> &                                     quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const Quadrature<dim - 1> &                                 quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  virtual void
  fill_fe_subface_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          sub_no,
    const Quadrature<dim - 1> &                                 quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  /**
   * Fields of cell-independent data.
   *
   * For information about the general purpose of this class, see the
   * documentation of the base class.
   */
  class InternalData : public FiniteElement<dim, spacedim>::InternalDataBase
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
    Table<2, double> shape_values;

    /**
     * Array with shape function gradients in quadrature points. There is one
     * row for each shape function, containing values for each quadrature
     * point.
     *
     * We store the gradients in the quadrature points on the unit cell. We
     * then only have to apply the transformation (which is a matrix-vector
     * multiplication) when visiting an actual cell.
     */
    Table<2, Tensor<1, dim>> shape_gradients;

    /**
     * Array with shape function hessians in quadrature points. There is one
     * row for each shape function, containing values for each quadrature
     * point.
     *
     * We store the hessians in the quadrature points on the unit cell. We
     * then only have to apply the transformation when visiting an actual
     * cell.
     */
    Table<2, Tensor<2, dim>> shape_hessians;

    /**
     * Array with shape function third derivatives in quadrature points. There
     * is one row for each shape function, containing values for each
     * quadrature point.
     *
     * We store the third derivatives in the quadrature points on the unit
     * cell. We then only have to apply the transformation when visiting an
     * actual cell.
     */
    Table<2, Tensor<3, dim>> shape_3rd_derivatives;
  };

  /**
   * Correct the shape Hessians by subtracting the terms corresponding to the
   * Jacobian pushed forward gradient.
   *
   * Before the correction, the Hessians would be given by
   * @f[
   * D_{ijk} = \frac{d^2\phi_i}{d \hat x_J d \hat x_K} (J_{jJ})^{-1}
   * (J_{kK})^{-1},
   * @f]
   * where $J_{iI}=\frac{d x_i}{d \hat x_I}$. After the correction, the
   * correct Hessians would be given by
   * @f[
   * \frac{d^2 \phi_i}{d x_j d x_k} = D_{ijk} - H_{mjk} \frac{d \phi_i}{d x_m},
   * @f]
   * where $H_{ijk}$ is the Jacobian pushed-forward derivative.
   */
  void
  correct_hessians(
    internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
      &output_data,
    const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &                mapping_data,
    const unsigned int n_q_points) const;

  /**
   * Correct the shape third derivatives by subtracting the terms
   * corresponding to the Jacobian pushed forward gradient and second
   * derivative.
   *
   * Before the correction, the third derivatives would be given by
   * @f[
   * D_{ijkl} = \frac{d^3\phi_i}{d \hat x_J d \hat x_K d \hat x_L} (J_{jJ})^{-1}
   * (J_{kK})^{-1} (J_{lL})^{-1},
   * @f]
   * where $J_{iI}=\frac{d x_i}{d \hat x_I}$. After the correction, the
   * correct third derivative would be given by
   * @f[
   * \frac{d^3\phi_i}{d x_j d x_k d x_l} = D_{ijkl} - H_{mjl} \frac{d^2
   * \phi_i}{d x_k d x_m}
   * - H_{mkl} \frac{d^2 \phi_i}{d x_j d x_m} - H_{mjk} \frac{d^2 \phi_i}{d x_l
   * d x_m}
   * - K_{mjkl} \frac{d \phi_i}{d x_m},
   * @f]
   * where $H_{ijk}$ is the Jacobian pushed-forward derivative and $K_{ijkl}$
   * is the Jacobian pushed-forward second derivative.
   */
  void
  correct_third_derivatives(
    internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
      &output_data,
    const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &                mapping_data,
    const unsigned int n_q_points) const;


  /**
   * The polynomial space.
   */
  const std::unique_ptr<ScalarPolynomialsBase<dim>> poly_space;
};

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif

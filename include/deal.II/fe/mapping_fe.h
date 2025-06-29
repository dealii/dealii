// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mapping_fe_h
#define dealii_mapping_fe_h


#include <deal.II/base/config.h>

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria_iterator.h>

#include <array>
#include <cmath>

DEAL_II_NAMESPACE_OPEN


/**
 * @addtogroup mapping
 * @{
 */

/**
 * This class consistently uses a user-provided finite element on all cells of a
 * triangulation to implement a polynomial mapping.
 *
 * If one initializes this class with the same FiniteElement as the
 * discretization, one obtains an iso-parametric mapping.
 *
 * If one initializes this class with an FE_Q(degree) object, then this class is
 * equivalent to MappingQ(degree). Please note that no optimizations
 * exploiting tensor-product structures of finite elements have been added here.
 *
 * @note Currently, only implemented for elements with tensor_degree==1 and
 *   n_components==1.
 *
 * Also see
 * @ref simplex "Simplex support".
 */
template <int dim, int spacedim = dim>
class MappingFE : public Mapping<dim, spacedim>
{
public:
  /**
   * Constructor.
   */
  explicit MappingFE(const FiniteElement<dim, spacedim> &fe);

  /**
   * Copy constructor.
   */
  MappingFE(const MappingFE<dim, spacedim> &mapping);

  // for documentation, see the Mapping base class
  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const override;

  /**
   * Return the degree of the mapping, i.e., the degree of the finite element
   * which was passed to the constructor.
   */
  unsigned int
  get_degree() const;

  // for documentation, see the Mapping base class
  virtual BoundingBox<spacedim>
  get_bounding_box(const typename Triangulation<dim, spacedim>::cell_iterator
                     &cell) const override;


  virtual bool
  is_compatible_with(const ReferenceCell &reference_cell) const override;

  /**
   * Always returns @p true because the default implementation of functions in
   * this class preserves vertex locations.
   */
  virtual bool
  preserves_vertex_locations() const override;

  /**
   * @name Mapping points between reference and real cells
   * @{
   */

  // for documentation, see the Mapping base class
  virtual Point<spacedim>
  transform_unit_to_real_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Point<dim> &p) const override;

  // for documentation, see the Mapping base class
  virtual Point<dim>
  transform_real_to_unit_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Point<spacedim> &p) const override;

  /**
   * @}
   */

  /**
   * @name Functions to transform tensors from reference to real coordinates
   * @{
   */

  // for documentation, see the Mapping base class
  virtual void
  transform(const ArrayView<const Tensor<1, dim>>                   &input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<1, spacedim>> &output) const override;

  // for documentation, see the Mapping base class
  virtual void
  transform(const ArrayView<const DerivativeForm<1, dim, spacedim>> &input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<2, spacedim>> &output) const override;

  // for documentation, see the Mapping base class
  virtual void
  transform(const ArrayView<const Tensor<2, dim>>                   &input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<2, spacedim>> &output) const override;

  // for documentation, see the Mapping base class
  virtual void
  transform(const ArrayView<const DerivativeForm<2, dim, spacedim>> &input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<3, spacedim>> &output) const override;

  // for documentation, see the Mapping base class
  virtual void
  transform(const ArrayView<const Tensor<3, dim>>                   &input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<3, spacedim>> &output) const override;

  /**
   * @}
   */

  /**
   * @name Interface with FEValues
   * @{
   */

  /**
   * Storage for internal data of polynomial mappings. See
   * Mapping::InternalDataBase for an extensive description.
   *
   * For the current class, the InternalData class stores data that is
   * computed once when the object is created (in get_data()) as well as data
   * the class wants to store from between the call to fill_fe_values(),
   * fill_fe_face_values(), or fill_fe_subface_values() until possible later
   * calls from the finite element to functions such as transform(). The
   * latter class of member variables are marked as 'mutable'.
   */
  class InternalData : public Mapping<dim, spacedim>::InternalDataBase
  {
  public:
    /**
     * Constructor.
     */
    InternalData(const FiniteElement<dim, spacedim> &fe);

    // Documentation see Mapping::InternalDataBase.
    virtual void
    reinit(const UpdateFlags      update_flags,
           const Quadrature<dim> &quadrature) override;

    /**
     * Initialize the object's member variables related to cell and face data
     * based on the given arguments. In order to initialize cell data, this
     * function calls reinit().
     */
    void
    initialize_face(const UpdateFlags      update_flags,
                    const Quadrature<dim> &quadrature,
                    const unsigned int     n_original_q_points);

    /**
     * Compute the values and/or derivatives of the shape functions used for
     * the mapping.
     */
    void
    compute_shape_function_values(const std::vector<Point<dim>> &unit_points);


    /**
     * Shape function at quadrature point. Shape functions are in tensor
     * product order, so vertices must be reordered to obtain transformation.
     */
    const double &
    shape(const unsigned int qpoint, const unsigned int shape_nr) const;

    /**
     * Shape function at quadrature point. See above.
     */
    double &
    shape(const unsigned int qpoint, const unsigned int shape_nr);

    /**
     * Gradient of shape function in quadrature point. See above.
     */
    const Tensor<1, dim> &
    derivative(const unsigned int qpoint, const unsigned int shape_nr) const;

    /**
     * Gradient of shape function in quadrature point. See above.
     */
    Tensor<1, dim> &
    derivative(const unsigned int qpoint, const unsigned int shape_nr);

    /**
     * Second derivative of shape function in quadrature point. See above.
     */
    const Tensor<2, dim> &
    second_derivative(const unsigned int qpoint,
                      const unsigned int shape_nr) const;

    /**
     * Second derivative of shape function in quadrature point. See above.
     */
    Tensor<2, dim> &
    second_derivative(const unsigned int qpoint, const unsigned int shape_nr);

    /**
     * third derivative of shape function in quadrature point. See above.
     */
    const Tensor<3, dim> &
    third_derivative(const unsigned int qpoint,
                     const unsigned int shape_nr) const;

    /**
     * third derivative of shape function in quadrature point. See above.
     */
    Tensor<3, dim> &
    third_derivative(const unsigned int qpoint, const unsigned int shape_nr);

    /**
     * fourth derivative of shape function in quadrature point. See above.
     */
    const Tensor<4, dim> &
    fourth_derivative(const unsigned int qpoint,
                      const unsigned int shape_nr) const;

    /**
     * fourth derivative of shape function in quadrature point. See above.
     */
    Tensor<4, dim> &
    fourth_derivative(const unsigned int qpoint, const unsigned int shape_nr);

    /**
     * Return an estimate (in bytes) for the memory consumption of this object.
     */
    virtual std::size_t
    memory_consumption() const override;

    /**
     * Values of shape functions. Access by function @p shape.
     *
     * Computed once.
     */
    std::vector<double> shape_values;

    /**
     * Values of shape function derivatives. Access by function @p derivative.
     *
     * Computed once.
     */
    std::vector<Tensor<1, dim>> shape_derivatives;

    /**
     * Values of shape function second derivatives. Access by function @p
     * second_derivative.
     *
     * Computed once.
     */
    std::vector<Tensor<2, dim>> shape_second_derivatives;

    /**
     * Values of shape function third derivatives. Access by function @p
     * second_derivative.
     *
     * Computed once.
     */
    std::vector<Tensor<3, dim>> shape_third_derivatives;

    /**
     * Values of shape function fourth derivatives. Access by function @p
     * second_derivative.
     *
     * Computed once.
     */
    std::vector<Tensor<4, dim>> shape_fourth_derivatives;

    /**
     * Unit tangential vectors. Used for the computation of boundary forms and
     * normal vectors.
     *
     * Filled once.
     */
    std::array<std::vector<Tensor<1, dim>>,
               ReferenceCells::max_n_faces<dim>() * 2>
      unit_tangentials;

    /**
     * Underlying finite element.
     */
    const FiniteElement<dim, spacedim> &fe;

    /**
     * The polynomial degree of the mapping.
     */
    const unsigned int polynomial_degree;

    /**
     * Number of shape functions.
     */
    const unsigned int n_shape_functions;

    /**
     * Tensors of covariant transformation at each of the quadrature points.
     * The matrix stored is the Jacobian * G^{-1}, where G = Jacobian^{t} *
     * Jacobian, is the first fundamental form of the map; if dim=spacedim
     * then it reduces to the transpose of the inverse of the Jacobian matrix,
     * which itself is stored in the @p contravariant field of this structure.
     *
     * Computed on each cell.
     */
    mutable std::vector<DerivativeForm<1, dim, spacedim>> covariant;

    /**
     * Tensors of contravariant transformation at each of the quadrature
     * points. The contravariant matrix is the Jacobian of the transformation,
     * i.e. $J_{ij}=dx_i/d\hat x_j$.
     *
     * Computed on each cell.
     */
    mutable std::vector<DerivativeForm<1, dim, spacedim>> contravariant;

    /**
     * Auxiliary vectors for internal use.
     */
    mutable std::vector<std::vector<Tensor<1, spacedim>>> aux;

    /**
     * Stores the support points of the mapping shape functions on the @p
     * cell_of_current_support_points.
     */
    mutable std::vector<Point<spacedim>> mapping_support_points;

    /**
     * Stores the cell of which the @p mapping_support_points are stored.
     */
    mutable typename Triangulation<dim, spacedim>::cell_iterator
      cell_of_current_support_points;

    /**
     * The determinant of the Jacobian in each quadrature point. Filled if
     * #update_volume_elements.
     */
    mutable std::vector<double> volume_elements;

    /**
     * Projected quadrature weights.
     */
    mutable std::vector<double> quadrature_weights;
  };


  // documentation can be found in Mapping::requires_update_flags()
  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

  // documentation can be found in Mapping::get_data()
  virtual std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
  get_data(const UpdateFlags, const Quadrature<dim> &quadrature) const override;

  using Mapping<dim, spacedim>::get_face_data;

  // documentation can be found in Mapping::get_face_data()
  virtual std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
  get_face_data(const UpdateFlags               flags,
                const hp::QCollection<dim - 1> &quadrature) const override;

  // documentation can be found in Mapping::get_subface_data()
  virtual std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
  get_subface_data(const UpdateFlags          flags,
                   const Quadrature<dim - 1> &quadrature) const override;

  // documentation can be found in Mapping::fill_fe_values()
  virtual CellSimilarity::Similarity
  fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim>                                      &quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const override;

  using Mapping<dim, spacedim>::fill_fe_face_values;

  // documentation can be found in Mapping::fill_fe_face_values()
  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const hp::QCollection<dim - 1>                             &quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const override;

  // documentation can be found in Mapping::fill_fe_subface_values()
  virtual void
  fill_fe_subface_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          subface_no,
    const Quadrature<dim - 1>                                  &quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const override;

  /**
   * @}
   */

protected:
  const std::unique_ptr<FiniteElement<dim, spacedim>> fe;

  /**
   * The degree of the polynomials used as shape functions for the mapping of
   * cells.
   */
  const unsigned int polynomial_degree;

  /**
   * Return the locations of support points for the mapping.
   */
  virtual std::vector<Point<spacedim>>
  compute_mapping_support_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell) const;

private:
  Table<2, double> mapping_support_point_weights;
};



/** @} */

/*----------------------------------------------------------------------*/

#ifndef DOXYGEN

template <int dim, int spacedim>
inline const double &
MappingFE<dim, spacedim>::InternalData::shape(const unsigned int qpoint,
                                              const unsigned int shape_nr) const
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr, shape_values.size());
  return shape_values[qpoint * n_shape_functions + shape_nr];
}



template <int dim, int spacedim>
inline double &
MappingFE<dim, spacedim>::InternalData::shape(const unsigned int qpoint,
                                              const unsigned int shape_nr)
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr, shape_values.size());
  return shape_values[qpoint * n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline const Tensor<1, dim> &
MappingFE<dim, spacedim>::InternalData::derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr) const
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_derivatives.size());
  return shape_derivatives[qpoint * n_shape_functions + shape_nr];
}



template <int dim, int spacedim>
inline Tensor<1, dim> &
MappingFE<dim, spacedim>::InternalData::derivative(const unsigned int qpoint,
                                                   const unsigned int shape_nr)
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_derivatives.size());
  return shape_derivatives[qpoint * n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline const Tensor<2, dim> &
MappingFE<dim, spacedim>::InternalData::second_derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr) const
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_second_derivatives.size());
  return shape_second_derivatives[qpoint * n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline Tensor<2, dim> &
MappingFE<dim, spacedim>::InternalData::second_derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr)
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_second_derivatives.size());
  return shape_second_derivatives[qpoint * n_shape_functions + shape_nr];
}

template <int dim, int spacedim>
inline const Tensor<3, dim> &
MappingFE<dim, spacedim>::InternalData::third_derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr) const
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_third_derivatives.size());
  return shape_third_derivatives[qpoint * n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline Tensor<3, dim> &
MappingFE<dim, spacedim>::InternalData::third_derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr)
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_third_derivatives.size());
  return shape_third_derivatives[qpoint * n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline const Tensor<4, dim> &
MappingFE<dim, spacedim>::InternalData::fourth_derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr) const
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_fourth_derivatives.size());
  return shape_fourth_derivatives[qpoint * n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline Tensor<4, dim> &
MappingFE<dim, spacedim>::InternalData::fourth_derivative(
  const unsigned int qpoint,
  const unsigned int shape_nr)
{
  AssertIndexRange(qpoint * n_shape_functions + shape_nr,
                   shape_fourth_derivatives.size());
  return shape_fourth_derivatives[qpoint * n_shape_functions + shape_nr];
}



template <int dim, int spacedim>
inline bool
MappingFE<dim, spacedim>::preserves_vertex_locations() const
{
  return true;
}



#endif // DOXYGEN

/* -------------- declaration of explicit specializations ------------- */


DEAL_II_NAMESPACE_CLOSE

#endif

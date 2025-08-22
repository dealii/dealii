// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mapping_p1_h
#define dealii_mapping_p1_h


#include <deal.II/base/config.h>

#include <deal.II/base/qprojector.h>

#include <deal.II/fe/mapping.h>

#include <cmath>


DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup mapping
 * @{
 */

/**
 * @brief Implementation of the classic affine transformation mapping used for
 * simplices.
 *
 * For simplex cells, the transformation mapping the unit to real cell may be
 * written as an affine transformation: that is
 *
 * @f[
 * \begin{pmatrix}
 *   x \\
 *   y
 * \end{pmatrix}
 * =
 * \begin{pmatrix}
 *   x_1 - x_0 & x_2 - x_0 \\
 *   y_1 - y_0 & y_2 - y_0
 * \end{pmatrix}
 * \begin{pmatrix}
 *   \hat{x} \\
 *   \hat{y}
 * \end{pmatrix}
 * +
 * \begin{pmatrix}
 *   x_0 \\
 *   y_0
 * \end{pmatrix}
 * @f]
 *
 * in which $(x_0, y_0)$, $(x_1, y_1)$ and $(x_2, y_2)$ are the vertices of a
 * triangle (in a codimension zero space). Unlike MappingQ1, this mapping is
 * affine, not bilinear or trilinear. Hence, this mapping's transformation
 * matrices (i.e., the covariant and contravariant transformations) and shift
 * vector are constant on each cell. This property makes this mapping
 * significantly more efficient than the equivalent MappingFE set up with a
 * linear FE_SimplexP.
 */
template <int dim, int spacedim = dim>
class MappingP1 : public Mapping<dim, spacedim>
{
public:
  /**
   * Default constructor.
   */
  MappingP1() = default;

  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const override;

  /**
   * Return @p true because MappingP1 preserves vertex locations.
   */
  virtual bool
  preserves_vertex_locations() const override;

  /**
   * This mapping is compatible with simplex reference cells.
   */
  virtual bool
  is_compatible_with(const ReferenceCell &reference_cell) const override;

  /**
   * @name Mapping points between reference and real cells
   * @{
   */

  virtual Point<spacedim>
  transform_unit_to_real_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Point<dim> &p) const override;

  virtual Point<dim>
  transform_real_to_unit_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Point<spacedim> &p) const override;

  virtual void
  transform_points_real_to_unit_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const ArrayView<const Point<spacedim>>                     &real_points,
    const ArrayView<Point<dim>> &unit_points) const override;

  /**
   * @}
   */

  /**
   * @name Functions to transform tensors from reference to real coordinates
   * @{
   */

  virtual void
  transform(const ArrayView<const Tensor<1, dim>>                   &input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<1, spacedim>> &output) const override;

  virtual void
  transform(const ArrayView<const DerivativeForm<1, dim, spacedim>> &input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<2, spacedim>> &output) const override;

  virtual void
  transform(const ArrayView<const Tensor<2, dim>>                   &input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<2, spacedim>> &output) const override;

  virtual void
  transform(const ArrayView<const DerivativeForm<2, dim, spacedim>> &input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<3, spacedim>> &output) const override;

  virtual void
  transform(const ArrayView<const Tensor<3, dim>>                   &input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<3, spacedim>> &output) const override;

  /**
   * @}
   */

  /**
   * As opposed to the other fill_fe_values() and fill_fe_face_values()
   * functions that rely on pre-computed information of InternalDataBase, this
   * function chooses the flexible evaluation path on the cell and points
   * passed in to the current function.
   *
   * @param[in] cell The cell where to evaluate the mapping
   *
   * @param[in] unit_points The points in reference coordinates where the
   * transformation (Jacobians, positions) should be computed.
   *
   * @param[in] update_flags The kind of information that should be computed.
   *
   * @param[out] output_data A struct containing the evaluated quantities such
   * as the Jacobian resulting from application of the mapping on the given
   * cell with its underlying manifolds.
   */
  void
  fill_mapping_data_for_generic_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const ArrayView<const Point<dim>>                          &unit_points,
    const UpdateFlags                                           update_flags,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const;

  /**
   * @name Interface with FEValues
   * @{
   */

  /**
   * Storage for internal data of the mapping. See Mapping::InternalDataBase
   * for an extensive description.
   *
   * This includes data that is computed once when the object is created (in
   * get_data()) as well as data the class wants to store from between the
   * call to fill_fe_values(), fill_fe_face_values(), or
   * fill_fe_subface_values() until possible later calls from the finite
   * element to functions such as transform(). The latter class of member
   * variables are marked as 'mutable'.
   */
  class InternalData : public Mapping<dim, spacedim>::InternalDataBase
  {
  public:
    /**
     * Constructor for use with arbitrary points.
     */
    InternalData(const ArrayView<const Point<dim>> &quadrature_points);

    /**
     * Constructor.
     */
    InternalData(const Quadrature<dim> &quadrature);

    virtual void
    reinit(const UpdateFlags      update_flags,
           const Quadrature<dim> &quadrature) override;

    /**
     * Return an estimate (in bytes) for the memory consumption of this object.
     */
    virtual std::size_t
    memory_consumption() const override;

    /**
     * Affine component of the transformation.
     */
    mutable Tensor<1, spacedim> affine_component;

    /**
     * Linear component of the transformation (the contravariant).
     */
    mutable DerivativeForm<1, dim, spacedim> linear_component;

    /**
     * Covariant form of the linear transformation.
     */
    mutable DerivativeForm<1, dim, spacedim> covariant;

    /**
     * Determinant of linear_component.
     */
    mutable double determinant;

    /**
     * Quadrature. May be an amalgamation of rules created by, e.g.,
     * QProjector::project_to_all_faces().
     */
    Quadrature<dim> quadrature;
  };

private:
  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

  virtual std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
  get_data(const UpdateFlags, const Quadrature<dim> &quadrature) const override;

  using Mapping<dim, spacedim>::get_face_data;

  virtual std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
  get_face_data(const UpdateFlags               flags,
                const hp::QCollection<dim - 1> &quadrature) const override;

  virtual std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
  get_subface_data(const UpdateFlags          flags,
                   const Quadrature<dim - 1> &quadrature) const override;

  virtual CellSimilarity::Similarity
  fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim>                                      &quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const override;

  using Mapping<dim, spacedim>::fill_fe_face_values;

  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const hp::QCollection<dim - 1>                             &quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const override;

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

  /**
   * Update the affine and linear components fields of the incoming InternalData
   * object.
   */
  void
  update_transformation(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const InternalData                                         &data) const;

  /**
   * Transform quadrature points from the unit cell to the real cell.
   */
  void
  transform_quadrature_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const InternalData                                         &data,
    const typename QProjector<dim>::DataSetDescriptor          &offset,
    std::vector<Point<spacedim>> &quadrature_points) const;

  /**
   * Compute the normal vectors if the UpdateFlags of the incoming InternalData
   * object say that they should be updated.
   */
  void
  maybe_update_normal_vectors(
    const unsigned int                face_no,
    const InternalData               &data,
    std::vector<Tensor<1, spacedim>> &normal_vectors) const;

  /**
   * Since the Jacobian is constant for this mapping all derivatives of the
   * Jacobian are identically zero. Fill these quantities with zeros if the
   * corresponding update flags say that they should be updated.
   */
  void
  maybe_update_jacobian_derivatives(
    const InternalData              &data,
    const CellSimilarity::Similarity cell_similarity,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const;

  /**
   * Compute the Jacobians if the UpdateFlags of the incoming
   * InternalData object say that they should be updated.
   */
  void
  maybe_update_jacobians(
    const InternalData              &data,
    const CellSimilarity::Similarity cell_similarity,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const;

  /**
   * Compute the inverse Jacobians if the UpdateFlags of the incoming
   * InternalData object say that they should be updated.
   */
  void
  maybe_update_inverse_jacobians(
    const InternalData              &data,
    const CellSimilarity::Similarity cell_similarity,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const;
};

/** @} */

DEAL_II_NAMESPACE_CLOSE

#endif

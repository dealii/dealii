// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mapping_cartesian_h
#define dealii_mapping_cartesian_h


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
 * A class providing a mapping from the reference cell to cells that are
 * axiparallel, i.e., that have the shape of rectangles (in 2d) or
 * boxes (in 3d) with edges parallel to the coordinate directions. The
 * class therefore provides functionality that is equivalent to what,
 * for example, MappingQ would provide for such cells. However, knowledge
 * of the shape of cells allows this class to be substantially more
 * efficient.
 *
 * Specifically, the mapping is meant for cells for which the mapping from
 * the reference to the real cell is a scaling along the coordinate
 * directions: The transformation from reference coordinates $\hat {\mathbf
 * x}$ to real coordinates $\mathbf x$ on each cell is of the form
 * @f{align*}{
 *   {\mathbf x}(\hat {\mathbf x})
 *   =
 *   \begin{pmatrix}
 *     h_x & 0 \\
 *     0 & h_y
 *   \end{pmatrix}
 *   \hat{\mathbf x}
 *   + {\mathbf v}_0
 * @f}
 * in 2d, and
 * @f{align*}{
 *   {\mathbf x}(\hat {\mathbf x})
 *   =
 *   \begin{pmatrix}
 *     h_x & 0 & 0 \\
 *     0 & h_y & 0 \\
 *     0 & 0 & h_z
 *   \end{pmatrix}
 *   \hat{\mathbf x}
 *   + {\mathbf v}_0
 * @f}
 * in 3d, where ${\mathbf v}_0$ is the bottom left vertex and $h_x,h_y,h_z$
 * are the extents of the cell along the axes.
 *
 * The class is intended for efficiency, and it does not do a whole lot of
 * error checking. If you apply this mapping to a cell that does not conform
 * to the requirements above, you will get strange results.
 */
template <int dim, int spacedim = dim>
class MappingCartesian : public Mapping<dim, spacedim>
{
public:
  // for documentation, see the Mapping base class
  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const override;

  /**
   * Return @p true because MappingCartesian preserves vertex
   * locations.
   */
  virtual bool
  preserves_vertex_locations() const override;

  virtual bool
  is_compatible_with(const ReferenceCell &reference_cell) const override;

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

  // for documentation, see the Mapping base class
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
     * Default constructor.
     */
    InternalData() = default;

    /**
     * Constructor that initializes the object with a quadrature.
     */
    InternalData(const Quadrature<dim> &quadrature);

    // Documentation see Mapping::InternalDataBase.
    virtual void
    reinit(const UpdateFlags      update_flags,
           const Quadrature<dim> &quadrature) override;

    /**
     * Return an estimate (in bytes) for the memory consumption of this object.
     */
    virtual std::size_t
    memory_consumption() const override;

    /**
     * Extents of the last cell we have seen in the coordinate directions,
     * i.e., <i>h<sub>x</sub></i>, <i>h<sub>y</sub></i>, <i>h<sub>z</sub></i>.
     */
    mutable Tensor<1, dim> cell_extents;

    /**
     * Reciprocal of the extents of the last cell we have seen in the
     * coordinate directions, i.e., <i>h<sub>x</sub></i>,
     * <i>h<sub>y</sub></i>, <i>h<sub>z</sub></i>.
     */
    mutable Tensor<1, dim> inverse_cell_extents;

    /**
     * The volume element
     */
    mutable double volume_element;

    /**
     * Location of quadrature points of faces or subfaces in 3d with all
     * possible orientations. Can be accessed with the correct offset provided
     * via QProjector::DataSetDescriptor. Not needed/used for cells.
     */
    std::vector<Point<dim>> quadrature_points;
  };

private:
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

  // documentation can be found in Mapping::fill_fe_immersed_surface_values()
  virtual void
  fill_fe_immersed_surface_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const NonMatching::ImmersedSurfaceQuadrature<dim>          &quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const override;

  /**
   * @}
   */

  /**
   * Update the cell_extents field of the incoming InternalData object with the
   * size of the incoming cell.
   */
  void
  update_cell_extents(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const InternalData                                         &data) const;

  /**
   * Compute the quadrature points if the UpdateFlags of the incoming
   * InternalData object say that they should be updated.
   *
   * Called from fill_fe_values.
   */
  void
  maybe_update_cell_quadrature_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const InternalData                                         &data,
    const ArrayView<const Point<dim>> &unit_quadrature_points,
    std::vector<Point<dim>>           &quadrature_points) const;

  /**
   * Compute the quadrature points if the UpdateFlags of the incoming
   * InternalData object say that they should be updated.
   *
   * Called from fill_fe_face_values.
   */
  void
  maybe_update_face_quadrature_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const InternalData                                         &data,
    std::vector<Point<dim>> &quadrature_points) const;

  /**
   * Compute the quadrature points if the UpdateFlags of the incoming
   * InternalData object say that they should be updated.
   *
   * Called from fill_fe_subface_values.
   */
  void
  maybe_update_subface_quadrature_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          sub_no,
    const InternalData                                         &data,
    std::vector<Point<dim>> &quadrature_points) const;

  /**
   * Compute the normal vectors if the UpdateFlags of the incoming InternalData
   * object say that they should be updated.
   */
  void
  maybe_update_normal_vectors(
    const unsigned int           face_no,
    const InternalData          &data,
    std::vector<Tensor<1, dim>> &normal_vectors) const;

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
   * Compute the volume elements if the UpdateFlags of the incoming
   * InternalData object say that they should be updated.
   */
  void
  maybe_update_volume_elements(const InternalData &data) const;

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

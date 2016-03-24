// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2016 by the deal.II authors
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

#ifndef dealii__mapping_cartesian_h
#define dealii__mapping_cartesian_h


#include <deal.II/base/config.h>
#include <deal.II/base/table.h>
#include <deal.II/fe/mapping.h>

#include <cmath>


DEAL_II_NAMESPACE_OPEN

/*!@addtogroup mapping */
/*@{*/

/**
 * A class providing a mapping from the reference cell to cells that are
 * axiparallel.
 *
 * This class maps the unit cell to a grid cell with surfaces parallel to the
 * coordinate lines/planes. It is specifically developed for Cartesian meshes.
 * In other words, the mapping is meant for cells for which the mapping from
 * the reference to the real cell is a scaling along the coordinate
 * directions: The transformation from reference coordinates $\hat {\mathbf
 * x}$ to real coordinates $\mathbf x$ on each cell is of the form
 * @f{align*}{
 * {\mathbf x}(\hat {\mathbf x}) = \begin{pmatrix} h_x & 0 \\ 0 & h_y
 * \end{pmatrix} \hat{\mathbf x} + {\mathbf v}_0
 * @f}
 * in 2d, and
 * @f{align*}{
 * {\mathbf x}(\hat {\mathbf x}) = \begin{pmatrix} h_x & 0 & 0 \\ 0 & h_y & 0
 * \\ 0 & 0 & h_z \end{pmatrix} \hat{\mathbf x} + {\mathbf v}_0
 * @f}
 * in 3d, where ${\mathbf v}_0$ is the bottom left vertex and $h_x,h_y,h_z$
 * are the extents of the cell along the axes.
 *
 * The class is intended for efficiency, and it does not do a whole lot of
 * error checking. If you apply this mapping to a cell that does not conform
 * to the requirements above, you will get strange results.
 *
 * @author Guido Kanschat, 2001; Ralf Hartmann, 2005
 */
template <int dim, int spacedim=dim>
class MappingCartesian : public Mapping<dim,spacedim>
{
public:

  // for documentation, see the Mapping base class
  virtual
  Mapping<dim, spacedim> *clone () const;

  /**
   * Always returns @p true because MappingCartesian preserves vertex
   * locations.
   */
  bool preserves_vertex_locations () const;

  /**
   * @name Mapping points between reference and real cells
   * @{
   */

  // for documentation, see the Mapping base class
  virtual
  Point<spacedim>
  transform_unit_to_real_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                               const Point<dim>                                          &p) const;

  // for documentation, see the Mapping base class
  virtual
  Point<dim>
  transform_real_to_unit_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                               const Point<spacedim>                                     &p) const;

  /**
   * @}
   */

  /**
   * @name Functions to transform tensors from reference to real coordinates
   * @{
   */

  // for documentation, see the Mapping base class
  virtual
  void
  transform (const ArrayView<const Tensor<1,dim> >                  &input,
             const MappingType                                       type,
             const typename Mapping<dim,spacedim>::InternalDataBase &internal,
             const ArrayView<Tensor<1,spacedim> >                   &output) const;

  // for documentation, see the Mapping base class
  virtual
  void
  transform (const ArrayView<const DerivativeForm<1, dim, spacedim> > &input,
             const MappingType                                         type,
             const typename Mapping<dim,spacedim>::InternalDataBase   &internal,
             const ArrayView<Tensor<2,spacedim> >                     &output) const;

  // for documentation, see the Mapping base class
  virtual
  void
  transform (const ArrayView<const Tensor<2, dim> >                 &input,
             const MappingType                                       type,
             const typename Mapping<dim,spacedim>::InternalDataBase &internal,
             const ArrayView<Tensor<2,spacedim> >                   &output) const;

  // for documentation, see the Mapping base class
  virtual
  void
  transform (const ArrayView<const DerivativeForm<2, dim, spacedim> > &input,
             const MappingType                                         type,
             const typename Mapping<dim,spacedim>::InternalDataBase   &internal,
             const ArrayView<Tensor<3,spacedim> >                     &output) const;

  // for documentation, see the Mapping base class
  virtual
  void
  transform (const ArrayView<const Tensor<3, dim> >                 &input,
             const MappingType                                       type,
             const typename Mapping<dim,spacedim>::InternalDataBase &internal,
             const ArrayView<Tensor<3,spacedim> >                   &output) const;

  /**
   * @}
   */


private:

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
     * Constructor.
     */
    InternalData (const Quadrature<dim> &quadrature);

    /**
     * Return an estimate (in bytes) or the memory consumption of this object.
     */
    virtual std::size_t memory_consumption () const;

    /**
     * Extents of the last cell we have seen in the coordinate directions,
     * i.e., <i>h<sub>x</sub></i>, <i>h<sub>y</sub></i>, <i>h<sub>z</sub></i>.
     */
    mutable Tensor<1,dim> cell_extents;

    /**
     * The volume element
     */
    mutable double volume_element;

    /**
     * Vector of all quadrature points. Especially, all points on all faces.
     */
    std::vector<Point<dim> > quadrature_points;
  };

  // documentation can be found in Mapping::requires_update_flags()
  virtual
  UpdateFlags
  requires_update_flags (const UpdateFlags update_flags) const;

  // documentation can be found in Mapping::get_data()
  virtual
  typename Mapping<dim, spacedim>::InternalDataBase *
  get_data (const UpdateFlags,
            const Quadrature<dim> &quadrature) const;

  // documentation can be found in Mapping::get_face_data()
  virtual
  typename Mapping<dim, spacedim>::InternalDataBase *
  get_face_data (const UpdateFlags flags,
                 const Quadrature<dim-1>& quadrature) const;

  // documentation can be found in Mapping::get_subface_data()
  virtual
  typename Mapping<dim, spacedim>::InternalDataBase *
  get_subface_data (const UpdateFlags flags,
                    const Quadrature<dim-1>& quadrature) const;

  // documentation can be found in Mapping::fill_fe_values()
  virtual
  CellSimilarity::Similarity
  fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                  const CellSimilarity::Similarity                           cell_similarity,
                  const Quadrature<dim>                                     &quadrature,
                  const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                  internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const;

  // documentation can be found in Mapping::fill_fe_face_values()
  virtual void
  fill_fe_face_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                       const unsigned int                                         face_no,
                       const Quadrature<dim-1>                                   &quadrature,
                       const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                       internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const;

  // documentation can be found in Mapping::fill_fe_subface_values()
  virtual void
  fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                          const unsigned int                                         face_no,
                          const unsigned int                                         subface_no,
                          const Quadrature<dim-1>                                   &quadrature,
                          const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                          internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const;

  /**
   * @}
   */



  /**
   * Do the computation for the <tt>fill_*</tt> functions.
   */
  void compute_fill (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                     const unsigned int face_no,
                     const unsigned int sub_no,
                     const CellSimilarity::Similarity cell_similarity,
                     const InternalData &data,
                     std::vector<Point<dim> > &quadrature_points,
                     std::vector<Tensor<1,dim> > &normal_vectors) const;

  /**
   * Value to indicate that a given face or subface number is invalid.
   */
  static const unsigned int invalid_face_number = numbers::invalid_unsigned_int;
};

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif

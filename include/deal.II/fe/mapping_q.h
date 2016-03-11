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

#ifndef dealii__mapping_q_h
#define dealii__mapping_q_h


#include <deal.II/base/config.h>
#include <deal.II/fe/mapping_q_generic.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup mapping */
/*@{*/

/**
 * A class that implements a polynomial mapping $Q_p$ of degree $p$ on cells
 * at the boundary of the domain (or, if requested in the constructor, for all
 * cells) and linear mappings for interior cells.
 *
 * The class is in fact poorly named since (unless explicitly specified during
 * the construction of the object, see below), it does not actually use
 * mappings of degree $p$ <i>everywhere</i>, but only on cells at the
 * boundary. This is in contrast to the MappingQGeneric class which indeed
 * does use a polynomial mapping $Q_p$ of degree $p$ everywhere. The point of
 * the current class is that in many situations, curved domains are only
 * provided with information about how exactly edges at the boundary are
 * shaped, but we do not know anything about internal edges. Thus, in the
 * absence of other information, we can only assume that internal edges are
 * straight lines, and in that case internal cells may as well be treated is
 * bilinear quadrilaterals or trilinear hexahedra. (An example of how such
 * meshes look is shown in step-1 already, but it is also discussed in the
 * "Results" section of step-6.) Because bi-/trilinear mappings are
 * significantly cheaper to compute than higher order mappings, it is
 * advantageous in such situations to use the higher order mapping only on
 * cells at the boundary of the domain. This class implements exactly this
 * behavior.
 *
 * There are a number of special cases worth considering:
 * - If you want to use a higher order mapping for all cells, you can
 * achieve this by setting the second argument to the constructor to true.
 * This only makes sense if you can actually provide information about how
 * interior edges and faces of the mesh should be curved. This is typically
 * done by associating a Manifold with interior cells and edges. A simple
 * example of this is discussed in the "Results" section of step-6; a full
 * discussion of manifolds is provided in step-53.
 * - If you pass true as the second argument to this class, then it
 * is in fact completely equivalent to generating a MappingQGeneric object
 * right away.
 * - This class is also entirely equivalent to MappingQGeneric if the
 * polynomial degree provided is one. This is because in that case, no
 * distinction between the mapping used on cells in the interior and on the
 * boundary of the domain can be made.
 * - If you are working on meshes embedded in higher space dimensions,
 * i.e., if dim!=spacedim, then every cell is considered to be at the boundary
 * of the domain and consequently a higher order mapping is used for all
 * cells; again this class is then equivalent to using MappingQGeneric right
 * away.
 *
 * @author Ralf Hartmann, 2000, 2001, 2005; Guido Kanschat 2000, 2001,
 * Wolfgang Bangerth, 2015
 */
template <int dim, int spacedim=dim>
class MappingQ : public Mapping<dim,spacedim>
{
public:
  /**
   * Constructor.  @p polynomial_degree denotes the polynomial degree of the
   * polynomials that are used to map cells boundary.
   *
   * The second argument determines whether the higher order mapping should
   * also be used on interior cells. If its value is <code>false</code> (the
   * default), then a lower order mapping is used in the interior. This is
   * sufficient for most cases where higher order mappings are only used to
   * better approximate the boundary. In that case, cells bounded by straight
   * lines are acceptable in the interior. However, there are cases where one
   * would also like to use a higher order mapping in the interior. The
   * MappingQEulerian class is one such case.
   *
   * The value of @p use_mapping_q_on_all_cells is ignored if @p dim is not
   * equal to @p spacedim, i.e., if we are considering meshes on surfaces
   * embedded into higher dimensional spaces.
   */
  MappingQ (const unsigned int polynomial_degree,
            const bool use_mapping_q_on_all_cells = false);

  /**
   * Copy constructor.
   */
  MappingQ (const MappingQ<dim,spacedim> &mapping);

  /**
   * Return the degree of the mapping, i.e. the value which was passed to the
   * constructor.
   */
  unsigned int get_degree () const;

  /**
   * Always returns @p true because the default implementation of functions in
   * this class preserves vertex locations.
   */
  virtual
  bool preserves_vertex_locations () const;

  /**
   * Transforms the point @p p on the unit cell to the point @p p_real on the
   * real cell @p cell and returns @p p_real.
   */
  virtual
  Point<spacedim>
  transform_unit_to_real_cell (
    const typename Triangulation<dim,spacedim>::cell_iterator &cell,
    const Point<dim>                                 &p) const;

  /**
   * Transforms the point @p p on the real cell to the point @p p_unit on the
   * unit cell @p cell and returns @p p_unit.
   *
   * Uses Newton iteration and the @p transform_unit_to_real_cell function.
   *
   * In the codimension one case, this function returns the normal projection
   * of the real point @p p on the curve or surface identified by the @p cell.
   *
   * @note Polynomial mappings from the reference (unit) cell coordinates to
   * the coordinate system of a real cell are not always invertible if the
   * point for which the inverse mapping is to be computed lies outside the
   * cell's boundaries.  In such cases, the current function may fail to
   * compute a point on the reference cell whose image under the mapping
   * equals the given point @p p.  If this is the case then this function
   * throws an exception of type Mapping::ExcTransformationFailed .  Whether
   * the given point @p p lies outside the cell can therefore be determined by
   * checking whether the return reference coordinates lie inside of outside
   * the reference cell (e.g., using GeometryInfo::is_inside_unit_cell) or
   * whether the exception mentioned above has been thrown.
   */
  virtual
  Point<dim>
  transform_real_to_unit_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                               const Point<spacedim>                                     &p) const;

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
   * Return a pointer to a copy of the present object. The caller of this copy
   * then assumes ownership of it.
   */
  virtual
  Mapping<dim,spacedim> *clone () const;


  /**
   * @name Interface with FEValues
   * @{
   */

protected:

  /**
   * Storage for internal data of this mapping. See Mapping::InternalDataBase
   * for an extensive description.
   *
   * This includes data that is computed once when the object is created (in
   * get_data()) as well as data the class wants to store from between the
   * call to fill_fe_values(), fill_fe_face_values(), or
   * fill_fe_subface_values() until possible later calls from the finite
   * element to functions such as transform(). The latter class of member
   * variables are marked as 'mutable'.
   *
   * The current class uses essentially the same fields for storage as the
   * MappingQGeneric class. Consequently, it inherits from
   * MappingQGeneric::InternalData, rather than from
   * Mapping::InternalDataBase. The principal difference to
   * MappingQGeneric::InternalData is that MappingQ switches between $Q_1$ and
   * $Q_p$ mappings depending on the cell we are on, so the internal data
   * object needs to also store a pointer to an InternalData object that
   * pertains to a $Q_1$ mapping.
   */
  class InternalData : public Mapping<dim,spacedim>::InternalDataBase
  {
  public:
    /**
     * Constructor.
     */
    InternalData ();


    /**
     * Return an estimate (in bytes) or the memory consumption of this object.
     */
    virtual std::size_t memory_consumption () const;

    /**
     * Flag that is set by the <tt>fill_fe_[[sub]face]_values</tt> function.
     *
     * If this flag is @p true we are on an interior cell and the @p
     * mapping_q1_data is used.
     */
    mutable bool use_mapping_q1_on_current_cell;

    /**
     * A pointer to a structure to store the information for the pure $Q_1$
     * mapping that is, by default, used on all interior cells.
     */
    std_cxx11::unique_ptr<typename MappingQGeneric<dim,spacedim>::InternalData> mapping_q1_data;

    /**
     * A pointer to a structure to store the information for the full $Q_p$
     * mapping that is, by default, used on all boundary cells.
     */
    std_cxx11::unique_ptr<typename MappingQGeneric<dim,spacedim>::InternalData> mapping_qp_data;
  };

protected:

  // documentation can be found in Mapping::requires_update_flags()
  virtual
  UpdateFlags
  requires_update_flags (const UpdateFlags update_flags) const;

  // documentation can be found in Mapping::get_data()
  virtual
  InternalData *
  get_data (const UpdateFlags,
            const Quadrature<dim> &quadrature) const;

  // documentation can be found in Mapping::get_face_data()
  virtual
  InternalData *
  get_face_data (const UpdateFlags flags,
                 const Quadrature<dim-1>& quadrature) const;

  // documentation can be found in Mapping::get_subface_data()
  virtual
  InternalData *
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

protected:

  /**
   * The polynomial degree of the cells to be used on all cells at the
   * boundary of the domain, or everywhere if so specified.
   */
  const unsigned int polynomial_degree;

  /**
   * If this flag is set @p true then @p MappingQ is used on all cells, not
   * only on boundary cells.
   */
  const bool use_mapping_q_on_all_cells;

  /**
   * Pointer to a Q1 mapping. This mapping is used on interior cells unless
   * use_mapping_q_on_all_cells was set in the call to the constructor. The
   * mapping is also used on any cell in the transform_real_to_unit_cell() to
   * compute a cheap initial guess for the position of the point before we
   * employ the more expensive Newton iteration using the full mapping.
   *
   * @note MappingQEulerian resets this pointer to an object of type
   * MappingQ1Eulerian to ensure that the Q1 mapping also knows about the
   * proper shifts and transformations of the Eulerian displacements. This
   * also means that we really need to store our own Q1 mapping here, rather
   * than simply resorting to StaticMappingQ1::mapping.
   *
   * @note If the polynomial degree used for the current object is one, then
   * the qp_mapping and q1_mapping variables point to the same underlying
   * object.
   */
  std_cxx11::shared_ptr<const MappingQGeneric<dim,spacedim> > q1_mapping;

  /**
   * Pointer to a Q_p mapping. This mapping is used on boundary cells unless
   * use_mapping_q_on_all_cells was set in the call to the constructor (in
   * which case it is used for all cells).
   *
   * @note MappingQEulerian and MappingC1 reset this pointer to an object of
   * their own implementation to ensure that the Q_p mapping also knows about
   * the proper shifts and transformations of the Eulerian displacements
   * (Eulerian case) and proper choice of support points (C1 case).
   *
   * @note If the polynomial degree used for the current object is one, then
   * the qp_mapping and q1_mapping variables point to the same underlying
   * object.
   */
  std_cxx11::shared_ptr<const MappingQGeneric<dim,spacedim> > qp_mapping;
};

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif

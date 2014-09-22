// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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

#ifndef __deal2__mapping_cartesian_h
#define __deal2__mapping_cartesian_h


#include <deal.II/base/config.h>
#include <deal.II/base/table.h>
#include <cmath>
#include <deal.II/fe/mapping.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup mapping */
/*@{*/

/**
 * A class providing a mapping from the reference cell to cells that are axiparallel.
 *
 * This class maps the unit cell to a grid cell with surfaces parallel
 * to the coordinate lines/planes. It is specifically developed for
 * Cartesian meshes. In other words, the mapping is meant for cells for which the
 * mapping from the reference to the real cell is a scaling
 * along the coordinate directions: The transformation from reference coordinates
 * $\hat {\mathbf x}$ to real coordinates $\mathbf x$ on each cell is of the
 * form
 * @f{align*}
 *   {\mathbf x}(\hat {\mathbf x})
 *   =
 *   \begin{pmatrix} h_x & 0 \\ 0 & h_y \end{pmatrix}
 *   \hat{\mathbf x} + {\mathbf v}_0
 * @f}
 * in 2d, and
 * @f{align*}
 *   {\mathbf x}(\hat {\mathbf x})
 *   =
 *   \begin{pmatrix} h_x & 0 & 0 \\ 0 & h_y & 0 \\ 0 & 0 & h_z \end{pmatrix}
 *   \hat{\mathbf x} + {\mathbf v}_0
 * @f}
 * in 3d, where ${\mathbf v}_0$ is the bottom left vertex and $h_x,h_y,h_z$ are
 * the extents of the cell along the axes.
 *
 * The class is intended for efficiency, and it does not do a whole lot of
 * error checking. If you apply this mapping to a cell that does not conform to
 * the requirements above, you will get strange results.
 *
 * @author Guido Kanschat, 2001; Ralf Hartmann, 2005
 */
template <int dim, int spacedim=dim>
class MappingCartesian : public Mapping<dim,spacedim>
{
public:
  virtual
  typename Mapping<dim, spacedim>::InternalDataBase *
  get_data (const UpdateFlags,
            const Quadrature<dim> &quadrature) const;

  virtual
  typename Mapping<dim, spacedim>::InternalDataBase *
  get_face_data (const UpdateFlags flags,
                 const Quadrature<dim-1>& quadrature) const;

  virtual
  typename Mapping<dim, spacedim>::InternalDataBase *
  get_subface_data (const UpdateFlags flags,
                    const Quadrature<dim-1>& quadrature) const;

  virtual void
  fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                  const Quadrature<dim>                                     &quadrature,
                  typename Mapping<dim, spacedim>::InternalDataBase         &mapping_data,
                  std::vector<Point<spacedim> >                             &quadrature_points,
                  std::vector<double>                                       &JxW_values,
                  std::vector<DerivativeForm<1,dim,spacedim> >      &jacobians,
                  std::vector<DerivativeForm<2,dim,spacedim> >     &jacobian_grads,
                  std::vector<DerivativeForm<1,spacedim,dim> >   &inverse_jacobians,
                  std::vector<Point<spacedim> > &,
                  CellSimilarity::Similarity                           &cell_similarity) const ;


  virtual void
  fill_fe_face_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                       const unsigned int face_no,
                       const Quadrature<dim-1>& quadrature,
                       typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
                       std::vector<Point<dim> >        &quadrature_points,
                       std::vector<double>             &JxW_values,
                       std::vector<Tensor<1,dim> >        &boundary_form,
                       std::vector<Point<spacedim> >        &normal_vectors) const ;
  virtual void
  fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                          const unsigned int face_no,
                          const unsigned int sub_no,
                          const Quadrature<dim-1>& quadrature,
                          typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
                          std::vector<Point<dim> >        &quadrature_points,
                          std::vector<double>             &JxW_values,
                          std::vector<Tensor<1,dim> >        &boundary_form,
                          std::vector<Point<spacedim> >        &normal_vectors) const ;

  virtual void
  transform (const VectorSlice<const std::vector<Tensor<1,dim> > > input,
             VectorSlice<std::vector<Tensor<1,spacedim> > > output,
             const typename Mapping<dim,spacedim>::InternalDataBase &internal,
             const MappingType type) const;

  virtual void
  transform (const VectorSlice<const std::vector<DerivativeForm<1, dim,spacedim> > > input,
             VectorSlice<std::vector<Tensor<2,spacedim> > > output,
             const typename Mapping<dim,spacedim>::InternalDataBase &internal,
             const MappingType type) const;


  virtual
  void
  transform (const VectorSlice<const std::vector<Tensor<2, dim> > >     input,
             VectorSlice<std::vector<Tensor<2,spacedim> > >             output,
             const typename Mapping<dim,spacedim>::InternalDataBase &internal,
             const MappingType type) const;

  virtual Point<spacedim>
  transform_unit_to_real_cell (
    const typename Triangulation<dim,spacedim>::cell_iterator &cell,
    const Point<dim>                                 &p) const;

  /**
   * Transforms the point @p p on
   * the real cell to the point
   * @p p_unit on the unit cell
   * @p cell and returns @p p_unit.
   *
   * Uses Newton iteration and the
   * @p transform_unit_to_real_cell
   * function.
   */
  virtual Point<dim>
  transform_real_to_unit_cell (
    const typename Triangulation<dim,spacedim>::cell_iterator &cell,
    const Point<spacedim>                            &p) const;


  /**
   * Return a pointer to a copy of the
   * present object. The caller of this
   * copy then assumes ownership of it.
   */
  virtual
  Mapping<dim, spacedim> *clone () const;

  /**
   * Always returns @p true because
   * MappingCartesian preserves vertex
   * locations.
   */
  bool preserves_vertex_locations () const;

protected:
  /**
   * Storage for internal data of
   * the scaling.
   */
  class InternalData : public Mapping<dim, spacedim>::InternalDataBase
  {
  public:
    /**
     * Constructor.
     */
    InternalData (const Quadrature<dim> &quadrature);

    /**
     * Return an estimate (in
     * bytes) or the memory
     * consumption of this
     * object.
     */
    virtual std::size_t memory_consumption () const;

    /**
     * Length of the cell in
     * different coordinate
     * directions, <i>h<sub>x</sub></i>,
     * <i>h<sub>y</sub></i>, <i>h<sub>z</sub></i>.
     */
    Tensor<1,dim> length;

    /**
     * The volume element
     */
    double volume_element;

    /**
     * Vector of all quadrature
     * points. Especially, all
     * points on all faces.
     */
    std::vector<Point<dim> > quadrature_points;
  };

  /**
   * Do the computation for the
   * <tt>fill_*</tt> functions.
   */
  void compute_fill (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                     const unsigned int face_no,
                     const unsigned int sub_no,
                     const CellSimilarity::Similarity cell_similarity,
                     InternalData &data,
                     std::vector<Point<dim> > &quadrature_points,
                     std::vector<Point<dim> > &normal_vectors) const;

private:
  virtual UpdateFlags update_once (const UpdateFlags) const;
  virtual UpdateFlags update_each (const UpdateFlags) const;

  /**
   * Value to indicate that a given
   * face or subface number is
   * invalid.
   */
  static const unsigned int invalid_face_number = numbers::invalid_unsigned_int;
};

/*@}*/

/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN

template <int dim, int spacedim>
inline
bool
MappingCartesian<dim,spacedim>::preserves_vertex_locations () const
{
  return true;
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

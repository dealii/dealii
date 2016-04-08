// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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

#ifndef dealii__mapping_manifold_h
#define dealii__mapping_manifold_h


#include <deal.II/base/config.h>
#include <deal.II/base/derivative_form.h>
#include <deal.II/base/table.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/mapping.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

template <int,int> class MappingQ;


/*!@addtogroup mapping */
/*@{*/


/**
 * This class implements the functionality for Manifold conforming
 * mappings. This Mapping computes the transformation between the
 * reference and real cell by exploiting the geometrical information
 * coming from the underlying Manifold object.
 *
 * Quadrature points computed using this mapping lie on the exact
 * geometrical objects, and tangent and normal vectors computed using
 * this class are tangent and normal to the underlying geometry. This
 * is in constrast with the MappingQ class, which approximates the
 * geometry using a polynomial of some order, and then computes the
 * normals and tangents using the approximated surface.
 *
 * @author Luca Heltai, Wolfgang Bangerth, Alberto Sartori 2016
 */
template <int dim, int spacedim=dim>
class MappingManifold : public Mapping<dim,spacedim>
{
public:
  /**
   * Constructor.
   */
  MappingManifold ();

  /**
   * Copy constructor.
   */
  MappingManifold (const MappingManifold<dim,spacedim> &mapping);

  // for documentation, see the Mapping base class
  virtual
  Mapping<dim,spacedim> *clone () const;

  /**
   * Always returns @p true because this class assumes that the
   * vertices always lies on the underlying Manifold.
   */
  virtual
  bool preserves_vertex_locations () const;

  /**
   * @name Mapping points between reference and real cells
   * @{
   */

  // for documentation, see the Mapping base class
  virtual
  Point<spacedim>
  transform_unit_to_real_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                               const Point<dim>                                 &p) const;

  // for documentation, see the Mapping base class
  virtual
  Point<dim>
  transform_real_to_unit_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                               const Point<spacedim>                            &p) const;

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

  /**
   * @name Interface with FEValues
   * @{
   */

public:
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
  class InternalData : public Mapping<dim,spacedim>::InternalDataBase
  {
  public:
    /**
     * Constructor.
     */
    InternalData();

    /**
     * Initialize the object's member variables related to cell data based on
     * the given arguments.
     *
     * The function also calls compute_shape_function_values() to actually set
     * the member variables related to the values and derivatives of the
     * mapping shape functions.
     */
    void
    initialize (const UpdateFlags      update_flags,
                const Quadrature<dim> &quadrature,
                const unsigned int     n_original_q_points);

    /**
     * Initialize the object's member variables related to cell and face data
     * based on the given arguments. In order to initialize cell data, this
     * function calls initialize().
     */
    void
    initialize_face (const UpdateFlags      update_flags,
                     const Quadrature<dim> &quadrature,
                     const unsigned int     n_original_q_points);


    /**
     * Compute the weights associated to the Manifold object, that
     * need to be passed when computing the location of the quadrature
     * points.
     */
    void
    compute_manifold_quadrature_weights(const Quadrature<dim> &quadrature);

    /**
     * Store vertices internally.
     */
    void
    store_vertices(const typename Triangulation<dim,spacedim>::cell_iterator &cell) const;

    /**
     * Return an estimate (in bytes) or the memory consumption of this object.
     */
    virtual std::size_t memory_consumption () const;

    /**
     * The current cell vertices.
     *
     * Computed each.
     */
    mutable std::vector<Point<spacedim> > vertices;

    /**
     * The current cell.
     *
     * Computed each.
     */
    mutable typename Triangulation<dim,spacedim>::cell_iterator cell;

    /**
     * The actual quadrature on the reference cell.
     *
     * Computed once.
     */
    Quadrature<dim> quad;


    /**
     * Values of quadrature weights for manifold quadrature
     * formulas.
     *
     * The Manifold class has a function (Manifold::get_new_point())
     * that returns new points according to a weighted average of some
     * surrounding points on the Manifold. For each quadrature point,
     * we call this function with a Quadrature formula constructed
     * using the vertices of the current cell, and the values of the
     * basis functions of an FE_Q(1) finite element evaluated at the
     * quadrature point itslef. While the vertices of the cell change
     * for every cell, the weights can be computed once for each
     * quadrature point. We store this information in the following
     * variable, where the first index runs through the quadrature
     * points, and the second index runs through the vertex indices.
     *
     * Computed once.
     */
    std::vector<std::vector<double> > cell_manifold_quadrature_weights;


    /**
     * A vector of weights for use in Manifold::get_new_point(). For
     * each point (interior to a cell), we compute the weight each
     * vertex has for this point. If the point lies at a vertex, then
     * this vertex has weight one and all others have weight zero. If
     * the point lies interior to a cell, then the weight every vertex
     * has is just the $d$-linear shape functions associated with each
     * vertex evaluated at that point.
     *
     * This array has size GeometryInfo<dim>::vertices_per_cell, but it
     * can't be converted into a fixed size array because it is used
     * as input for Manifold::get_new_point() which wants to see a
     * std::vector<double> for the weights.
     */
    mutable std::vector<double> vertex_weights;

    /**
     * Unit tangential vectors. Used for the computation of boundary forms and
     * normal vectors.
     *
     * This array has (dim-1)*GeometryInfo::faces_per_cell entries. The first
     * GeometryInfo::faces_per_cell contain the vectors in the first
     * tangential direction for each face; the second set of
     * GeometryInfo::faces_per_cell entries contain the vectors in the second
     * tangential direction (only in 3d, since there we have 2 tangential
     * directions per face), etc.
     *
     * Filled once.
     */
    std_cxx11::array<std::vector<Tensor<1,dim> >, GeometryInfo<dim>::faces_per_cell *(dim-1)> unit_tangentials;

    /**
     * Tensors of covariant transformation at each of the quadrature points.
     * The matrix stored is the Jacobian * G^{-1}, where G = Jacobian^{t} *
     * Jacobian, is the first fundamental form of the map; if dim=spacedim
     * then it reduces to the transpose of the inverse of the Jacobian matrix,
     * which itself is stored in the @p contravariant field of this structure.
     *
     * Computed on each cell.
     */
    mutable std::vector<DerivativeForm<1,dim, spacedim > >  covariant;

    /**
     * Tensors of contravariant transformation at each of the quadrature
     * points. The contravariant matrix is the Jacobian of the transformation,
     * i.e. $J_{ij}=dx_i/d\hat x_j$.
     *
     * Computed on each cell.
     */
    mutable std::vector<DerivativeForm<1,dim,spacedim> > contravariant;

    /**
     * Auxiliary vectors for internal use.
     */
    mutable std::vector<std::vector<Tensor<1,spacedim> > > aux;

    /**
     * The determinant of the Jacobian in each quadrature point. Filled if
     * #update_volume_elements.
     */
    mutable std::vector<double> volume_elements;

    /**
     * A pointer to the Manifold in use.
     *
     * Updated each.
     */
    mutable SmartPointer<const Manifold<dim,spacedim> > manifold;
  };


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
  fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator     &cell,
                  const CellSimilarity::Similarity                               cell_similarity,
                  const Quadrature<dim>                                         &quadrature,
                  const typename Mapping<dim,spacedim>::InternalDataBase        &internal_data,
                  dealii::internal::FEValues::MappingRelatedData<dim, spacedim> &output_data) const;

  // documentation can be found in Mapping::fill_fe_face_values()
  virtual void
  fill_fe_face_values (const typename Triangulation<dim,spacedim>::cell_iterator     &cell,
                       const unsigned int                                             face_no,
                       const Quadrature<dim-1>                                       &quadrature,
                       const typename Mapping<dim,spacedim>::InternalDataBase        &internal_data,
                       dealii::internal::FEValues::MappingRelatedData<dim, spacedim> &output_data) const;

  // documentation can be found in Mapping::fill_fe_subface_values()
  virtual void
  fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator     &cell,
                          const unsigned int                                             face_no,
                          const unsigned int                                             subface_no,
                          const Quadrature<dim-1>                                       &quadrature,
                          const typename Mapping<dim,spacedim>::InternalDataBase        &internal_data,
                          dealii::internal::FEValues::MappingRelatedData<dim, spacedim> &output_data) const;

  /**
   * @}
   */
};



/*@}*/

/*----------------------------------------------------------------------*/

#ifndef DOXYGEN

template<int dim, int spacedim>
inline
void
MappingManifold<dim,spacedim>::InternalData::store_vertices (const typename Triangulation<dim,spacedim>::cell_iterator &cell) const
{
  vertices.resize(GeometryInfo<dim>::vertices_per_cell);
  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
    vertices[i] = cell->vertex(i);
  this->cell = cell;
}


template <int dim, int spacedim>
inline
void
MappingManifold<dim,spacedim>::InternalData::compute_manifold_quadrature_weights (const Quadrature<dim> &quad)
{
  cell_manifold_quadrature_weights.resize(quad.size(), std::vector<double>(GeometryInfo<dim>::vertices_per_cell));
  for (unsigned int q=0; q<quad.size(); ++q)
    {
      for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
        {
          cell_manifold_quadrature_weights[q][i] = GeometryInfo<dim>::d_linear_shape_function(quad.point(q), i);
        }
    }
}



template <int dim, int spacedim>
inline
bool
MappingManifold<dim,spacedim>::preserves_vertex_locations () const
{
  return true;
}

#endif // DOXYGEN

/* -------------- declaration of explicit specializations ------------- */


DEAL_II_NAMESPACE_CLOSE

#endif

// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2016 by the deal.II authors
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


#include <deal.II/base/derivative_form.h>
#include <deal.II/base/config.h>
#include <deal.II/base/table.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/fe_q.h>

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
 * Quadrature points computed using this mapping lye on the exact
 * geometrical objects, and tangent and normal vectors computed using
 * this class are normal and tangent to the underlying geometry. This
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
   * Always returns @p true because the default implementation of functions in
   * this class preserves vertex locations.
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
     * Return an estimate (in bytes) or the memory consumption of this object.
     */
    virtual std::size_t memory_consumption () const;

    /**
     * Store the current cell.
     *
     * Computed each.
     */
    typename Triangulation<dim,spacedim>::cell_iterator current_cell;

    /**
     * Values of quadrature weights for manifold quadrature formulas.
     *
     * Computed once.
     */
    std::vector<std::vector<double> > cell_manifold_quadratures_weights;

    // /**
    //  * Values of shape function derivatives. Access by function @p derivative.
    //  *
    //  * Computed once.
    //  */
    // std::vector<Tensor<1,dim> > shape_derivatives;

    // /**
    //  * Values of shape function second derivatives. Access by function @p
    //  * second_derivative.
    //  *
    //  * Computed once.
    //  */
    // std::vector<Tensor<2,dim> > shape_second_derivatives;

    // /**
    //  * Values of shape function third derivatives. Access by function @p
    //  * second_derivative.
    //  *
    //  * Computed once.
    //  */
    // std::vector<Tensor<3,dim> > shape_third_derivatives;

    // /**
    //  * Values of shape function fourth derivatives. Access by function @p
    //  * second_derivative.
    //  *
    //  * Computed once.
    //  */
    // std::vector<Tensor<4,dim> > shape_fourth_derivatives;

    // /**
    //  * Unit tangential vectors. Used for the computation of boundary forms and
    //  * normal vectors.
    //  *
    //  * This vector has (dim-1)GeometryInfo::faces_per_cell entries. The first
    //  * GeometryInfo::faces_per_cell contain the vectors in the first
    //  * tangential direction for each face; the second set of
    //  * GeometryInfo::faces_per_cell entries contain the vectors in the second
    //  * tangential direction (only in 3d, since there we have 2 tangential
    //  * directions per face), etc.
    //  *
    //  * Filled once.
    //  */
    // std::vector<std::vector<Tensor<1,dim> > > unit_tangentials;

    // /**
    //  * The polynomial degree of the mapping. Since the objects here are also
    //  * used (with minor adjustments) by MappingQ, we need to store this.
    //  */
    // unsigned int polynomial_degree;

    // /**
    //  * Number of shape functions. If this is a Q1 mapping, then it is simply
    //  * the number of vertices per cell. However, since also derived classes
    //  * use this class (e.g. the Mapping_Q() class), the number of shape
    //  * functions may also be different.
    //  *
    //  * In general, it is $(p+1)^\text{dim}$, where $p$ is the polynomial
    //  * degree of the mapping.
    //  */
    // const unsigned int n_shape_functions;

    // /**
    //  * Tensors of covariant transformation at each of the quadrature points.
    //  * The matrix stored is the Jacobian * G^{-1}, where G = Jacobian^{t} *
    //  * Jacobian, is the first fundamental form of the map; if dim=spacedim
    //  * then it reduces to the transpose of the inverse of the Jacobian matrix,
    //  * which itself is stored in the @p contravariant field of this structure.
    //  *
    //  * Computed on each cell.
    //  */
    // mutable std::vector<DerivativeForm<1,dim, spacedim > >  covariant;

    // /**
    //  * Tensors of contravariant transformation at each of the quadrature
    //  * points. The contravariant matrix is the Jacobian of the transformation,
    //  * i.e. $J_{ij}=dx_i/d\hat x_j$.
    //  *
    //  * Computed on each cell.
    //  */
    // mutable std::vector< DerivativeForm<1,dim,spacedim> > contravariant;

    // /**
    //  * Auxiliary vectors for internal use.
    //  */
    // mutable std::vector<std::vector<Tensor<1,spacedim> > > aux;

    // /**
    //  * Stores the support points of the mapping shape functions on the @p
    //  * cell_of_current_support_points.
    //  */
    // mutable std::vector<Point<spacedim> > mapping_support_points;

    // /**
    //  * Stores the cell of which the @p mapping_support_points are stored.
    //  */
    // mutable typename Triangulation<dim,spacedim>::cell_iterator cell_of_current_support_points;

    // /**
    //  * The determinant of the Jacobian in each quadrature point. Filled if
    //  * #update_volume_elements.
    //  */
    // mutable std::vector<double> volume_elements;
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

protected:

  /**
   * An FE_Q object, used to compute weights for Manifold quadratures.
   */
  const FE_Q<dim,spacedim> fe_q;

  /**
   * A table of weights by which we multiply the locations of the support
   * points on the perimeter of a quad to get the location of interior support
   * points.
   *
   * Sizes: support_point_weights_on_quad.size()= number of inner
   * unit_support_points support_point_weights_on_quad[i].size()= number of
   * outer unit_support_points, i.e.  unit_support_points on the boundary of
   * the quad
   *
   * For the definition of this vector see equation (8) of the `mapping'
   * report.
   */
  Table<2,double> support_point_weights_on_quad;

  /**
   * A table of weights by which we multiply the locations of the support
   * points on the perimeter of a hex to get the location of interior support
   * points.
   *
   * For the definition of this vector see equation (8) of the `mapping'
   * report.
   */
  Table<2,double> support_point_weights_on_hex;

  /**
   * Return the locations of support points for the mapping. For example, for
   * $Q_1$ mappings these are the vertices, and for higher order polynomial
   * mappings they are the vertices plus interior points on edges, faces, and
   * the cell interior that are placed in consultation with the Manifold
   * description of the domain and its boundary. However, other classes may
   * override this function differently. In particular, the MappingQ1Eulerian
   * class does exactly this by not computing the support points from the
   * geometry of the current cell but instead evaluating an externally given
   * displacement field in addition to the geometry of the cell.
   *
   * The default implementation of this function is appropriate for most
   * cases. It takes the locations of support points on the boundary of the
   * cell from the underlying manifold. Interior support points (ie. support
   * points in quads for 2d, in hexes for 3d) are then computed using the
   * solution of a Laplace equation with the position of the outer support
   * points as boundary values, in order to make the transformation as smooth
   * as possible.
   *
   * The function works its way from the vertices (which it takes from the
   * given cell) via the support points on the line (for which it calls the
   * add_line_support_points() function) and the support points on the quad
   * faces (in 3d, for which it calls the add_quad_support_points() function).
   * It then adds interior support points that are either computed by
   * interpolation from the surrounding points using weights computed by
   * solving a Laplace equation, or if dim<spacedim, it asks the underlying
   * manifold for the locations of interior points.
   */
  // virtual
  // std::vector<Point<spacedim> >
  // compute_mapping_support_points (const typename Triangulation<dim,spacedim>::cell_iterator &cell) const;

  /**
   * Transforms the point @p p on the real cell to the corresponding point on
   * the unit cell @p cell by a Newton iteration.
   */
  // Point<dim>
  // transform_real_to_unit_cell_internal (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  //                                       const Point<spacedim> &p,
  //                                       const Point<dim> &initial_p_unit) const;

  /**
   * Make MappingQ a friend since it needs to call the fill_fe_values()
   * functions on its MappingManifold(1) sub-object.
   */
  template <int, int> friend class MappingQ;

};



/*@}*/

/*----------------------------------------------------------------------*/

#ifndef DOXYGEN

// template<int dim, int spacedim>
// inline
// const Tensor<1,dim> &
// MappingManifold<dim,spacedim>::InternalData::derivative (const unsigned int qpoint,
//                                                          const unsigned int shape_nr) const
// {
//   Assert(qpoint*n_shape_functions + shape_nr < shape_derivatives.size(),
//          ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
//                        shape_derivatives.size()));
//   return shape_derivatives [qpoint*n_shape_functions + shape_nr];
// }



// template<int dim, int spacedim>
// inline
// Tensor<1,dim> &
// MappingManifold<dim,spacedim>::InternalData::derivative (const unsigned int qpoint,
//                                                          const unsigned int shape_nr)
// {
//   Assert(qpoint*n_shape_functions + shape_nr < shape_derivatives.size(),
//          ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
//                        shape_derivatives.size()));
//   return shape_derivatives [qpoint*n_shape_functions + shape_nr];
// }


// template <int dim, int spacedim>
// inline
// const Tensor<2,dim> &
// MappingManifold<dim,spacedim>::InternalData::second_derivative (const unsigned int qpoint,
//     const unsigned int shape_nr) const
// {
//   Assert(qpoint*n_shape_functions + shape_nr < shape_second_derivatives.size(),
//          ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
//                        shape_second_derivatives.size()));
//   return shape_second_derivatives [qpoint*n_shape_functions + shape_nr];
// }


// template <int dim, int spacedim>
// inline
// Tensor<2,dim> &
// MappingManifold<dim,spacedim>::InternalData::second_derivative (const unsigned int qpoint,
//     const unsigned int shape_nr)
// {
//   Assert(qpoint*n_shape_functions + shape_nr < shape_second_derivatives.size(),
//          ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
//                        shape_second_derivatives.size()));
//   return shape_second_derivatives [qpoint*n_shape_functions + shape_nr];
// }

// template <int dim, int spacedim>
// inline
// const Tensor<3,dim> &
// MappingManifold<dim,spacedim>::InternalData::third_derivative (const unsigned int qpoint,
//     const unsigned int shape_nr) const
// {
//   Assert(qpoint*n_shape_functions + shape_nr < shape_third_derivatives.size(),
//          ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
//                        shape_third_derivatives.size()));
//   return shape_third_derivatives [qpoint*n_shape_functions + shape_nr];
// }


// template <int dim, int spacedim>
// inline
// Tensor<3,dim> &
// MappingManifold<dim,spacedim>::InternalData::third_derivative (const unsigned int qpoint,
//     const unsigned int shape_nr)
// {
//   Assert(qpoint*n_shape_functions + shape_nr < shape_third_derivatives.size(),
//          ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
//                        shape_third_derivatives.size()));
//   return shape_third_derivatives [qpoint*n_shape_functions + shape_nr];
// }


// template <int dim, int spacedim>
// inline
// const Tensor<4,dim> &
// MappingManifold<dim,spacedim>::InternalData::fourth_derivative (const unsigned int qpoint,
//     const unsigned int shape_nr) const
// {
//   Assert(qpoint*n_shape_functions + shape_nr < shape_fourth_derivatives.size(),
//          ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
//                        shape_fourth_derivatives.size()));
//   return shape_fourth_derivatives [qpoint*n_shape_functions + shape_nr];
// }


// template <int dim, int spacedim>
// inline
// Tensor<4,dim> &
// MappingManifold<dim,spacedim>::InternalData::fourth_derivative (const unsigned int qpoint,
//     const unsigned int shape_nr)
// {
//   Assert(qpoint*n_shape_functions + shape_nr < shape_fourth_derivatives.size(),
//          ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
//                        shape_fourth_derivatives.size()));
//   return shape_fourth_derivatives [qpoint*n_shape_functions + shape_nr];
// }


template <int dim, int spacedim>
inline
void
MappingManifold<dim,spacedim>::InternalData::compute_manifold_quadrature_weights (const Quadrature<dim> &quad)
{
  static FE_Q<dim> fe_q(1);
  cell_manifold_quadratures_weights.resize(quad.size(), std::vector<double>(GeometryInfo<dim>::vertices_per_cell));
  for (unsigned int q=0; q<quad.size(); ++q)
    {
      for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
        {
          cell_manifold_quadratures_weights[q][i] = fe_q.shape_value(i, quad.point(q));
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

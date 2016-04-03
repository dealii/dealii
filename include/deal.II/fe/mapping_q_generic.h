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

#ifndef dealii__mapping_q_generic_h
#define dealii__mapping_q_generic_h


#include <deal.II/base/config.h>
#include <deal.II/base/derivative_form.h>
#include <deal.II/base/table.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/base/std_cxx11/array.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

template <int,int> class MappingQ;


/*!@addtogroup mapping */
/*@{*/


/**
 * This class implements the functionality for polynomial mappings $Q_p$ of
 * polynomial degree $p$ that will be used on all cells of the mesh. The
 * MappingQ1 and MappingQ classes specialize this behavior slightly.
 *
 * The class is poorly named. It should really have been called MappingQ
 * because it consistently uses $Q_p$ mappings on all cells of a
 * triangulation. However, the name MappingQ was already taken when we rewrote
 * the entire class hierarchy for mappings. One might argue that one should
 * always use MappingQGeneric over the existing class MappingQ (which, unless
 * explicitly specified during the construction of the object, only uses
 * mappings of degree $p$ <i>on cells at the boundary of the domain</i>). On
 * the other hand, there are good reasons to use MappingQ in many situations:
 * in many situations, curved domains are only provided with information about
 * how exactly edges at the boundary are shaped, but we do not know anything
 * about internal edges. Thus, in the absence of other information, we can
 * only assume that internal edges are straight lines, and in that case
 * internal cells may as well be treated is bilinear quadrilaterals or
 * trilinear hexahedra. (An example of how such meshes look is shown in step-1
 * already, but it is also discussed in the "Results" section of step-6.)
 * Because bi-/trilinear mappings are significantly cheaper to compute than
 * higher order mappings, it is advantageous in such situations to use the
 * higher order mapping only on cells at the boundary of the domain -- i.e.,
 * the behavior of MappingQ. Of course, MappingQGeneric also uses bilinear
 * mappings for interior cells as long as it has no knowledge about curvature
 * of interior edges, but it implements this the expensive way: as a general
 * $Q_p$ mapping where the mapping support points just <i>happen</i> to be
 * arranged along linear or bilinear edges or faces.
 *
 * There are a number of special cases worth considering:
 * - If you really want to use a higher order mapping for all cells,
 * you can do this using the current class, but this only makes sense if you
 * can actually provide information about how interior edges and faces of the
 * mesh should be curved. This is typically done by associating a Manifold
 * with interior cells and edges. A simple example of this is discussed in the
 * "Results" section of step-6; a full discussion of manifolds is provided in
 * step-53.
 * - If you are working on meshes that describe a (curved) manifold
 * embedded in higher space dimensions, i.e., if dim!=spacedim, then every
 * cell is at the boundary of the domain you will likely already have attached
 * a manifold object to all cells that can then also be used by the mapping
 * classes for higher order mappings.
 *
 *
 * @author Wolfgang Bangerth, 2015
 */
template <int dim, int spacedim=dim>
class MappingQGeneric : public Mapping<dim,spacedim>
{
public:
  /**
   * Constructor.  @p polynomial_degree denotes the polynomial degree of the
   * polynomials that are used to map cells from the reference to the real
   * cell.
   */
  MappingQGeneric (const unsigned int polynomial_degree);

  /**
   * Copy constructor.
   */
  MappingQGeneric (const MappingQGeneric<dim,spacedim> &mapping);

  // for documentation, see the Mapping base class
  virtual
  Mapping<dim,spacedim> *clone () const;

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
     * Constructor. The argument denotes the polynomial degree of the mapping
     * to which this object will correspond.
     */
    InternalData(const unsigned int polynomial_degree);

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
     * Compute the values and/or derivatives of the shape functions used for
     * the mapping.
     *
     * Which values, derivatives, or higher order derivatives are computed is
     * determined by which of the member arrays have nonzero sizes. They are
     * typically set to their appropriate sizes by the initialize() and
     * initialize_face() functions, which indeed call this function
     * internally. However, it is possible (and at times useful) to do the
     * resizing by hand and then call this function directly. An example is in
     * a Newton iteration where we update the location of a quadrature point
     * (e.g., in MappingQ::transform_real_to_uni_cell()) and need to re-
     * compute the mapping and its derivatives at this location, but have
     * already sized all internal arrays correctly.
     */
    void compute_shape_function_values (const std::vector<Point<dim> > &unit_points);


    /**
     * Shape function at quadrature point. Shape functions are in tensor
     * product order, so vertices must be reordered to obtain transformation.
     */
    const double &shape (const unsigned int qpoint,
                         const unsigned int shape_nr) const;

    /**
     * Shape function at quadrature point. See above.
     */
    double &shape (const unsigned int qpoint,
                   const unsigned int shape_nr);

    /**
     * Gradient of shape function in quadrature point. See above.
     */
    const Tensor<1,dim> &derivative (const unsigned int qpoint,
                                     const unsigned int shape_nr) const;

    /**
     * Gradient of shape function in quadrature point. See above.
     */
    Tensor<1,dim> &derivative (const unsigned int qpoint,
                               const unsigned int shape_nr);

    /**
     * Second derivative of shape function in quadrature point. See above.
     */
    const Tensor<2,dim> &second_derivative (const unsigned int qpoint,
                                            const unsigned int shape_nr) const;

    /**
     * Second derivative of shape function in quadrature point. See above.
     */
    Tensor<2,dim> &second_derivative (const unsigned int qpoint,
                                      const unsigned int shape_nr);

    /**
     * third derivative of shape function in quadrature point. See above.
     */
    const Tensor<3,dim> &third_derivative (const unsigned int qpoint,
                                           const unsigned int shape_nr) const;

    /**
     * third derivative of shape function in quadrature point. See above.
     */
    Tensor<3,dim> &third_derivative (const unsigned int qpoint,
                                     const unsigned int shape_nr);

    /**
     * fourth derivative of shape function in quadrature point. See above.
     */
    const Tensor<4,dim> &fourth_derivative (const unsigned int qpoint,
                                            const unsigned int shape_nr) const;

    /**
     * fourth derivative of shape function in quadrature point. See above.
     */
    Tensor<4,dim> &fourth_derivative (const unsigned int qpoint,
                                      const unsigned int shape_nr);

    /**
     * Return an estimate (in bytes) or the memory consumption of this object.
     */
    virtual std::size_t memory_consumption () const;

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
    std::vector<Tensor<1,dim> > shape_derivatives;

    /**
     * Values of shape function second derivatives. Access by function @p
     * second_derivative.
     *
     * Computed once.
     */
    std::vector<Tensor<2,dim> > shape_second_derivatives;

    /**
     * Values of shape function third derivatives. Access by function @p
     * second_derivative.
     *
     * Computed once.
     */
    std::vector<Tensor<3,dim> > shape_third_derivatives;

    /**
     * Values of shape function fourth derivatives. Access by function @p
     * second_derivative.
     *
     * Computed once.
     */
    std::vector<Tensor<4,dim> > shape_fourth_derivatives;

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
     * The polynomial degree of the mapping. Since the objects here are also
     * used (with minor adjustments) by MappingQ, we need to store this.
     */
    unsigned int polynomial_degree;

    /**
     * Number of shape functions. If this is a Q1 mapping, then it is simply
     * the number of vertices per cell. However, since also derived classes
     * use this class (e.g. the Mapping_Q() class), the number of shape
     * functions may also be different.
     *
     * In general, it is $(p+1)^\text{dim}$, where $p$ is the polynomial
     * degree of the mapping.
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
    mutable std::vector<DerivativeForm<1,dim, spacedim > >  covariant;

    /**
     * Tensors of contravariant transformation at each of the quadrature
     * points. The contravariant matrix is the Jacobian of the transformation,
     * i.e. $J_{ij}=dx_i/d\hat x_j$.
     *
     * Computed on each cell.
     */
    mutable std::vector< DerivativeForm<1,dim,spacedim> > contravariant;

    /**
     * Auxiliary vectors for internal use.
     */
    mutable std::vector<std::vector<Tensor<1,spacedim> > > aux;

    /**
     * Stores the support points of the mapping shape functions on the @p
     * cell_of_current_support_points.
     */
    mutable std::vector<Point<spacedim> > mapping_support_points;

    /**
     * Stores the cell of which the @p mapping_support_points are stored.
     */
    mutable typename Triangulation<dim,spacedim>::cell_iterator cell_of_current_support_points;

    /**
     * The determinant of the Jacobian in each quadrature point. Filled if
     * #update_volume_elements.
     */
    mutable std::vector<double> volume_elements;
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
   * The degree of the polynomials used as shape functions for the mapping of
   * cells.
   */
  const unsigned int polynomial_degree;

  /*
   * The default line support points. These are used when computing
   * the location in real space of the support points on lines and
   * quads, which are asked to the Manifold<dim,spacedim> class.
   *
   * The number of quadrature points depends on the degree of this
   * class, and it matches the number of degrees of freedom of an
   * FE_Q<1>(this->degree).
   */
  QGaussLobatto<1> line_support_points;

  /**
   * An FE_Q object which is only needed in 3D, since it knows how to reorder
   * shape functions/DoFs on non-standard faces. This is used to reorder
   * support points in the same way.
   */
  const std_cxx11::unique_ptr<FE_Q<dim> > fe_q;

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
  virtual
  std::vector<Point<spacedim> >
  compute_mapping_support_points (const typename Triangulation<dim,spacedim>::cell_iterator &cell) const;

  /**
   * Transforms the point @p p on the real cell to the corresponding point on
   * the unit cell @p cell by a Newton iteration.
   */
  Point<dim>
  transform_real_to_unit_cell_internal (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                                        const Point<spacedim> &p,
                                        const Point<dim> &initial_p_unit) const;

  /**
   * For <tt>dim=2,3</tt>. Append the support points of all shape functions
   * located on bounding lines of the given cell to the vector @p a. Points
   * located on the vertices of a line are not included.
   *
   * Needed by the @p compute_support_points() function. For <tt>dim=1</tt>
   * this function is empty. The function uses the underlying manifold object
   * of the line (or, if none is set, of the cell) for the location of the
   * requested points.
   *
   * This function is made virtual in order to allow derived classes to choose
   * shape function support points differently than the present class, which
   * chooses the points as interpolation points on the boundary.
   */
  virtual
  void
  add_line_support_points (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                           std::vector<Point<spacedim> > &a) const;

  /**
   * For <tt>dim=3</tt>. Append the support points of all shape functions
   * located on bounding faces (quads in 3d) of the given cell to the vector
   * @p a. Points located on the vertices or lines of a quad are not included.
   *
   * Needed by the @p compute_support_points() function. For <tt>dim=1</tt>
   * and <tt>dim=2</tt> this function is empty. The function uses the
   * underlying manifold object of the quad (or, if none is set, of the cell)
   * for the location of the requested points.
   *
   * This function is made virtual in order to allow derived classes to choose
   * shape function support points differently than the present class, which
   * chooses the points as interpolation points on the boundary.
   */
  virtual
  void
  add_quad_support_points(const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                          std::vector<Point<spacedim> > &a) const;

  /**
   * Make MappingQ a friend since it needs to call the fill_fe_values()
   * functions on its MappingQGeneric(1) sub-object.
   */
  template <int, int> friend class MappingQ;
};



/*@}*/

/*----------------------------------------------------------------------*/

#ifndef DOXYGEN

template<int dim, int spacedim>
inline
const double &
MappingQGeneric<dim,spacedim>::InternalData::shape (const unsigned int qpoint,
                                                    const unsigned int shape_nr) const
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_values.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_values.size()));
  return shape_values [qpoint*n_shape_functions + shape_nr];
}



template<int dim, int spacedim>
inline
double &
MappingQGeneric<dim,spacedim>::InternalData::shape (const unsigned int qpoint,
                                                    const unsigned int shape_nr)
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_values.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_values.size()));
  return shape_values [qpoint*n_shape_functions + shape_nr];
}


template<int dim, int spacedim>
inline
const Tensor<1,dim> &
MappingQGeneric<dim,spacedim>::InternalData::derivative (const unsigned int qpoint,
                                                         const unsigned int shape_nr) const
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_derivatives.size()));
  return shape_derivatives [qpoint*n_shape_functions + shape_nr];
}



template<int dim, int spacedim>
inline
Tensor<1,dim> &
MappingQGeneric<dim,spacedim>::InternalData::derivative (const unsigned int qpoint,
                                                         const unsigned int shape_nr)
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_derivatives.size()));
  return shape_derivatives [qpoint*n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline
const Tensor<2,dim> &
MappingQGeneric<dim,spacedim>::InternalData::second_derivative (const unsigned int qpoint,
    const unsigned int shape_nr) const
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_second_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_second_derivatives.size()));
  return shape_second_derivatives [qpoint*n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline
Tensor<2,dim> &
MappingQGeneric<dim,spacedim>::InternalData::second_derivative (const unsigned int qpoint,
    const unsigned int shape_nr)
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_second_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_second_derivatives.size()));
  return shape_second_derivatives [qpoint*n_shape_functions + shape_nr];
}

template <int dim, int spacedim>
inline
const Tensor<3,dim> &
MappingQGeneric<dim,spacedim>::InternalData::third_derivative (const unsigned int qpoint,
    const unsigned int shape_nr) const
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_third_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_third_derivatives.size()));
  return shape_third_derivatives [qpoint*n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline
Tensor<3,dim> &
MappingQGeneric<dim,spacedim>::InternalData::third_derivative (const unsigned int qpoint,
    const unsigned int shape_nr)
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_third_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_third_derivatives.size()));
  return shape_third_derivatives [qpoint*n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline
const Tensor<4,dim> &
MappingQGeneric<dim,spacedim>::InternalData::fourth_derivative (const unsigned int qpoint,
    const unsigned int shape_nr) const
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_fourth_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_fourth_derivatives.size()));
  return shape_fourth_derivatives [qpoint*n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline
Tensor<4,dim> &
MappingQGeneric<dim,spacedim>::InternalData::fourth_derivative (const unsigned int qpoint,
    const unsigned int shape_nr)
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_fourth_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_fourth_derivatives.size()));
  return shape_fourth_derivatives [qpoint*n_shape_functions + shape_nr];
}



template <int dim, int spacedim>
inline
bool
MappingQGeneric<dim,spacedim>::preserves_vertex_locations () const
{
  return true;
}

#endif // DOXYGEN

/* -------------- declaration of explicit specializations ------------- */


DEAL_II_NAMESPACE_CLOSE

#endif

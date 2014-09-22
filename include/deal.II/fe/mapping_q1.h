// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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

#ifndef __deal2__mapping_q1_h
#define __deal2__mapping_q1_h


#include <deal.II/base/derivative_form.h>
#include <deal.II/base/config.h>
#include <deal.II/base/table.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/mapping.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup mapping */
/*@{*/


/**
 * Mapping of general quadrilateral/hexahedra by d-linear shape
 * functions.
 *
 * This function maps the unit cell to a general grid cell with
 * straight lines in $d$ dimensions (remark that in 3D the surfaces
 * may be curved, even if the edges are not). This is the well-known
 * mapping for polyhedral domains.
 *
 * Shape function for this mapping are the same as for the finite
 * element FE_Q of order 1. Therefore, coupling these two yields
 * an isoparametric element.
 *
 * For more information about the <tt>spacedim</tt> template parameter
 * check the documentation of FiniteElement or the one of
 * Triangulation.
 *
 * @author Guido Kanschat, 2000, 2001; Ralf Hartmann, 2000, 2001, 2005
 */
template <int dim, int spacedim=dim>
class MappingQ1 : public Mapping<dim,spacedim>
{
public:
  /**
   * Default constructor.
   */
  MappingQ1 ();

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
   *
   * In the codimension one case,
   * this function returns the
   * normal projection of the real
   * point @p p on the curve or
   * surface identified by the @p
   * cell.
   *
   * @note Polynomial mappings from
   * the reference (unit) cell coordinates
   * to the coordinate system of a real
   * cell are not always invertible if
   * the point for which the inverse
   * mapping is to be computed lies
   * outside the cell's boundaries.
   * In such cases, the current function
   * may fail to compute a point on
   * the reference cell whose image
   * under the mapping equals the given
   * point @p p.  If this is the case
   * then this function throws an
   * exception of type
   * Mapping::ExcTransformationFailed .
   * Whether the given point @p p lies
   * outside the cell can therefore be
   * determined by checking whether the
   * return reference coordinates lie
   * inside of outside the reference
   * cell (e.g., using
   * GeometryInfo::is_inside_unit_cell)
   * or whether the exception mentioned
   * above has been thrown.
   */
  virtual Point<dim>
  transform_real_to_unit_cell (
    const typename Triangulation<dim,spacedim>::cell_iterator &cell,
    const Point<spacedim>                            &p) const;

  virtual void
  transform (const VectorSlice<const std::vector<Tensor<1,dim> > > input,
             VectorSlice<std::vector<Tensor<1,spacedim> > > output,
             const typename Mapping<dim,spacedim>::InternalDataBase &internal,
             const MappingType type) const;

  virtual void
  transform (const VectorSlice<const std::vector<DerivativeForm<1, dim,spacedim> > >    input,
             VectorSlice<std::vector<Tensor<2,spacedim> > > output,
             const typename Mapping<dim,spacedim>::InternalDataBase &internal,
             const MappingType type) const;

  virtual
  void
  transform (const VectorSlice<const std::vector<Tensor<2, dim> > >     input,
             VectorSlice<std::vector<Tensor<2,spacedim> > >             output,
             const typename Mapping<dim,spacedim>::InternalDataBase &internal,
             const MappingType type) const;


protected:
  /**
     This function and the next allow to generate the transform require by
     the virtual transform() in mapping, but unfortunately in C++ one cannot
     declare a virtual template function.
  */
  template < int rank >
  void
  transform_fields(const VectorSlice<const std::vector<Tensor<rank,dim>      > > input,
                   VectorSlice<      std::vector<Tensor<rank,spacedim> > > output,
                   const typename Mapping<dim,spacedim>::InternalDataBase &internal,
                   const MappingType type) const;
  /**
     see doc in transform_fields
  */
  template < int rank >
  void
  transform_gradients(const VectorSlice<const std::vector<Tensor<rank,dim>      > > input,
                      VectorSlice<      std::vector<Tensor<rank,spacedim> > > output,
                      const typename Mapping<dim,spacedim>::InternalDataBase &internal,
                      const MappingType type) const;
  /**
     see doc in transform_fields
  */
  template < int rank >
  void
  transform_differential_forms(
    const VectorSlice<const std::vector<DerivativeForm<rank, dim, spacedim> > >    input,
    VectorSlice<std::vector<DerivativeForm<rank, spacedim, spacedim> > > output,
    const typename Mapping<dim,spacedim>::InternalDataBase &internal,
    const MappingType type) const;

public:




  /**
   * Return a pointer to a copy of the
   * present object. The caller of this
   * copy then assumes ownership of it.
   */
  virtual
  Mapping<dim,spacedim> *clone () const;

  /**
   * Storage for internal data of
   * d-linear transformation.
   */
  class InternalData : public Mapping<dim,spacedim>::InternalDataBase
  {
  public:
    /**
     * Constructor. Pass the
     * number of shape functions.
     */
    InternalData(const unsigned int n_shape_functions);

    /**
     * Shape function at quadrature
     * point. Shape functions are
     * in tensor product order, so
     * vertices must be reordered
     * to obtain transformation.
     */
    double shape (const unsigned int qpoint,
                  const unsigned int shape_nr) const;

    /**
     * Shape function at quadrature
     * point. See above.
     */
    double &shape (const unsigned int qpoint,
                   const unsigned int shape_nr);

    /**
     * Gradient of shape function
     * in quadrature point. See
     * above.
     */
    Tensor<1,dim> derivative (const unsigned int qpoint,
                              const unsigned int shape_nr) const;

    /**
     * Gradient of shape function
     * in quadrature point. See
     * above.
     */
    Tensor<1,dim> &derivative (const unsigned int qpoint,
                               const unsigned int shape_nr);

    /**
     * Second derivative of shape
     * function in quadrature
     * point. See above.
     */
    Tensor<2,dim> second_derivative (const unsigned int qpoint,
                                     const unsigned int shape_nr) const;

    /**
     * Second derivative of shape
     * function in quadrature
     * point. See above.
     */
    Tensor<2,dim> &second_derivative (const unsigned int qpoint,
                                      const unsigned int shape_nr);

    /**
     * Return an estimate (in
     * bytes) or the memory
     * consumption of this
     * object.
     */
    virtual std::size_t memory_consumption () const;

    /**
     * Values of shape
     * functions. Access by
     * function @p shape.
     *
     * Computed once.
     */
    std::vector<double> shape_values;

    /**
     * Values of shape function
     * derivatives. Access by
     * function @p derivative.
     *
     * Computed once.
     */
    std::vector<Tensor<1,dim> > shape_derivatives;

    /**
     * Values of shape function
     * second derivatives. Access
     * by function
     * @p second_derivative.
     *
     * Computed once.
     */
    std::vector<Tensor<2,dim> > shape_second_derivatives;

    /**
     * Tensors of covariant
     * transformation at each of
     * the quadrature points. The
     * matrix stored is the
     * Jacobian * G^{-1},
     * where G = Jacobian^{t} * Jacobian,
     * is the first fundamental
     * form of the map;
     * if dim=spacedim then
     * it reduces to the transpose of the
     * inverse of the Jacobian
     * matrix, which itself is
     * stored in the
     * @p contravariant field of
     * this structure.
     *
     * Computed on each cell.
     */
    std::vector<DerivativeForm<1,dim, spacedim > >  covariant;

    /**
     * Tensors of contravariant
     * transformation at each of
     * the quadrature points. The
     * contravariant matrix is
     * the Jacobian of the
     * transformation,
     * i.e. $J_{ij}=dx_i/d\hat x_j$.
     *
     * Computed on each cell.
     */
    std::vector< DerivativeForm<1,dim,spacedim> > contravariant;

    /**
     * Unit tangential vectors. Used
     * for the computation of
     * boundary forms and normal
     * vectors.
     *
     * This vector has
     * (dim-1)GeometryInfo::faces_per_cell
     * entries. The first
     * GeometryInfo::faces_per_cell
     * contain the vectors in the first
     * tangential direction for each
     * face; the second set of
     * GeometryInfo::faces_per_cell
     * entries contain the vectors in the
     * second tangential direction (only
     * in 3d, since there we have 2
     * tangential directions per face),
     * etc.
     *
     * Filled once.
     */
    std::vector<std::vector<Tensor<1,dim> > > unit_tangentials;

    /**
     * Auxiliary vectors for internal use.
     */
    std::vector<std::vector<Tensor<1,spacedim> > > aux;

    /**
     * Stores the support points of
     * the mapping shape functions on
     * the @p cell_of_current_support_points.
     */
    std::vector<Point<spacedim> > mapping_support_points;

    /**
     * Stores the cell of which the
     * @p mapping_support_points are
     * stored.
     */
    typename Triangulation<dim,spacedim>::cell_iterator cell_of_current_support_points;

    /**
     * Default value of this flag
     * is @p true. If <tt>*this</tt>
     * is an object of a derived
     * class, this flag is set to
     * @p false.
     */
    bool is_mapping_q1_data;

    /**
     * Number of shape
     * functions. If this is a Q1
     * mapping, then it is simply
     * the number of vertices per
     * cell. However, since also
     * derived classes use this
     * class (e.g. the
     * Mapping_Q() class),
     * the number of shape
     * functions may also be
     * different.
     */
    unsigned int n_shape_functions;
  };

  /**
   * Declare a convenience typedef
   * for the class that describes
   * offsets into quadrature
   * formulas projected onto faces
   * and subfaces.
   */
  typedef
  typename QProjector<dim>::DataSetDescriptor
  DataSetDescriptor;

  /**
   * Implementation of the interface in
   * Mapping.
   */
  virtual void
  fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                  const Quadrature<dim>                                     &quadrature,
                  typename Mapping<dim,spacedim>::InternalDataBase          &mapping_data,
                  typename std::vector<Point<spacedim> >                    &quadrature_points,
                  std::vector<double>                                       &JxW_values,
                  std::vector<DerivativeForm<1,dim,spacedim> >        &jacobians,
                  std::vector<DerivativeForm<2,dim,spacedim> >       &jacobian_grads,
                  std::vector<DerivativeForm<1,spacedim,dim> >      &inverse_jacobians,
                  std::vector<Point<spacedim> >                             &cell_normal_vectors,
                  CellSimilarity::Similarity                           &cell_similarity) const;

  /**
   * Implementation of the interface in
   * Mapping.
   */
  virtual void
  fill_fe_face_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                       const unsigned int                               face_no,
                       const Quadrature<dim-1>                          &quadrature,
                       typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
                       typename std::vector<Point<spacedim> >                &quadrature_points,
                       std::vector<double>                              &JxW_values,
                       typename std::vector<Tensor<1,spacedim> >             &boundary_form,
                       typename std::vector<Point<spacedim> >           &normal_vectors) const ;

  /**
   * Implementation of the interface in
   * Mapping.
   */
  virtual void
  fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                          const unsigned int face_no,
                          const unsigned int sub_no,
                          const Quadrature<dim-1>& quadrature,
                          typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
                          typename std::vector<Point<spacedim> >        &quadrature_points,
                          std::vector<double>             &JxW_values,
                          typename std::vector<Tensor<1,spacedim> >        &boundary_form,
                          typename std::vector<Point<spacedim> >        &normal_vectors) const ;

  /**
   * Compute shape values and/or
   * derivatives.
   *
   * Calls either the
   * @p compute_shapes_virtual of
   * this class or that of the
   * derived class, depending on
   * whether
   * <tt>data.is_mapping_q1_data</tt>
   * equals @p true or @p false.
   */
  void compute_shapes (const std::vector<Point<dim> > &unit_points,
                       InternalData &data) const;

  /**
   * Do the computations for the
   * @p get_data functions. Here,
   * the data vectors of
   * @p InternalData are
   * reinitialized to proper size
   * and shape values are computed.
   */
  void compute_data (const UpdateFlags flags,
                     const Quadrature<dim> &quadrature,
                     const unsigned int n_orig_q_points,
                     InternalData &data) const;

  /**
   * Do the computations for the
   * @p get_face_data
   * functions. Here, the data
   * vectors of @p InternalData
   * are reinitialized to proper
   * size and shape values and
   * derivatives are
   * computed. Furthermore
   * @p unit_tangential vectors of
   * the face are computed.
   */
  void compute_face_data (const UpdateFlags flags,
                          const Quadrature<dim> &quadrature,
                          const unsigned int n_orig_q_points,
                          InternalData &data) const;

  /**
   * Do the computation for the
   * <tt>fill_*</tt> functions.
   */
  void compute_fill (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                     const unsigned int      npts,
                     const DataSetDescriptor data_set,
                     const CellSimilarity::Similarity cell_similarity,
                     InternalData           &data,
                     std::vector<Point<spacedim> > &quadrature_points) const;

  /**
   * Do the computation for the
   * <tt>fill_*</tt> functions.
   */
  void compute_fill_face (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                          const unsigned int      face_no,
                          const unsigned int      subface_no,
                          const unsigned int      npts,
                          const DataSetDescriptor data_set,
                          const std::vector<double>   &weights,
                          InternalData           &mapping_data,
                          std::vector<Point<spacedim> >    &quadrature_points,
                          std::vector<double>         &JxW_values,
                          std::vector<Tensor<1,spacedim> > &boundary_form,
                          std::vector<Point<spacedim> > &normal_vectors) const;

  /**
   * Compute shape values and/or
   * derivatives.
   */
  virtual void compute_shapes_virtual (const std::vector<Point<dim> > &unit_points,
                                       InternalData &data) const;

  /**
   * Transforms a point @p p on
   * the unit cell to the point
   * @p p_real on the real cell
   * @p cell and returns @p p_real.
   *
   * This function is called by
   * @p transform_unit_to_real_cell
   * and multiple times (through the
   * Newton iteration) by
   * @p transform_real_to_unit_cell_internal.
   *
   * Takes a reference to an
   * @p InternalData that must
   * already include the shape
   * values at point @p p and the
   * mapping support points of the
   * cell.
   *
   * This @p InternalData argument
   * avoids multiple computations
   * of the shape values at point
   * @p p and especially multiple
   * computations of the mapping
   * support points.
   */
  Point<spacedim>
  transform_unit_to_real_cell_internal (const InternalData &mdata) const;

  /**
   * Transforms the point @p p on
   * the real cell to the corresponding
   * point on the unit cell
   * @p cell by a Newton
   * iteration.
   *
   * Takes a reference to an
   * @p InternalData that is
   * assumed to be previously
   * created by the @p get_data
   * function with @p UpdateFlags
   * including
   * @p update_transformation_values
   * and
   * @p update_transformation_gradients
   * and a one point Quadrature
   * that includes the given
   * initial guess for the
   * transformation
   * @p initial_p_unit.  Hence this
   * function assumes that
   * @p mdata already includes the
   * transformation shape values
   * and gradients computed at
   * @p initial_p_unit.
   *
   * @p mdata will be changed by
   * this function.
   */
  Point<dim>
  transform_real_to_unit_cell_internal (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                                        const Point<spacedim> &p,
                                        const Point<dim> &initial_p_unit,
                                        InternalData &mdata) const;

  /**
   * Always returns @p true because
   * MappingQ1 preserves vertex locations.
   */
  virtual
  bool preserves_vertex_locations () const;

protected:
  /* Trick to templatize transform_real_to_unit_cell<dim, dim+1> */
  template<int dim_>
  Point<dim_>
  transform_real_to_unit_cell_internal_codim1
  (const typename Triangulation<dim_,dim_+1>::cell_iterator &cell,
   const Point<dim_+1> &p,
   const Point<dim_>         &initial_p_unit,
   InternalData        &mdata) const;

  /**
     Compute an initial guess to pass to the Newton method in
     transform_real_to_unit_cell.

     For the initial guess we proceed in the following way:
     <ul>
     <li> find the least square dim-dimensional plane
          approximating the cell vertices, i.e. we find and affine
          map A x_hat + b from the reference cell to the real space.
     <li> Solve the equation A x_hat + b = p for x_hat
     <li> This x_hat is the initial solution used for the Newton Method.
     </ul>

     @note if dim<spacedim we first project p onto the plane.
     @note if dim==1 (for any spacedim) the initial guess is the exact solution
     and no Newton iteration is needed.


     Some details about how we compute the least square plane.
     We look for a  spacedim x (dim + 1) matrix  X such that

     X * M = Y

     where M is a (dim+1) x n_vertices  matrix and Y a spacedim x n_vertices.

     And:
     The i-th column of M is unit_vertex[i] and the last row all 1's.
     The i-th column of Y is real_vertex[i].

     If we split X=[A|b], the least square approx is A x_hat+b

     Classically  X = Y * (M^t (M M^t)^{-1})

     Let K = M^t * (M M^t)^{-1} = [KA Kb]
     this can be precomputed, and that is exactely
     what we do.

     Finally A = Y*KA  and  b = Y*Kb.
  */
  Point<dim>
  transform_real_to_unit_cell_initial_guess (const std::vector<Point<spacedim> > &vertex,
                                             const Point<spacedim>                            &p) const;


private:
  /**
   * Implementation of the interface in
   * Mapping.
   *
   * Description of effects:
   * <ul>
   * <li> if @p update_quadrature_points
   * is required, the output will
   * contain
   * @p update_transformation_values. This
   * computes the values of the
   * transformation basis
   * polynomials at the unit cell
   * quadrature points.
   * <li> if any of
   * @p update_covariant_transformation,
   * @p update_contravariant_transformation,
   * @p update_JxW_values,
   * @p update_boundary_forms,
   * @p update_normal_vectors is
   * required, the output will
   * contain
   * @p update_transformation_gradients
   * to compute derivatives of the
   * transformation basis
   * polynomials.
   * </ul>
   */
  virtual UpdateFlags update_once (const UpdateFlags flags) const;

  /**
   * Implementation of the interface in
   * Mapping.
   *
   * Description of effects if
   * @p flags contains:
   * <ul>
   * <li> <code>update_quadrature_points</code> is
   * copied to the output to
   * compute the quadrature points
   * on the real cell.
   * <li> <code>update_JxW_values</code> is
   * copied and requires
   * @p update_boundary_forms on
   * faces. The latter, because the
   * surface element is just the
   * norm of the boundary form.
   * <li> <code>update_normal_vectors</code>
   * is copied and requires
   * @p update_boundary_forms. The
   * latter, because the normal
   * vector is the normalized
   * boundary form.
   * <li>
   * <code>update_covariant_transformation</code>
   * is copied and requires
   * @p update_contravariant_transformation,
   * since it is computed as the
   * inverse of the latter.
   * <li> <code>update_JxW_values</code> is
   * copied and requires
   * <code>update_contravariant_transformation</code>,
   * since it is computed as one
   * over determinant of the
   * latter.
   * <li> <code>update_boundary_forms</code>
   * is copied and requires
   * <code>update_contravariant_transformation</code>,
   * since the boundary form is
   * computed as the contravariant
   * image of the normal vector to
   * the unit cell.
   * </ul>
   */
  virtual UpdateFlags update_each (const UpdateFlags flags) const;

  virtual
  typename Mapping<dim,spacedim>::InternalDataBase *
  get_data (const UpdateFlags,
            const Quadrature<dim> &quadrature) const;

  virtual
  typename Mapping<dim,spacedim>::InternalDataBase *
  get_face_data (const UpdateFlags flags,
                 const Quadrature<dim-1>& quadrature) const;

  virtual
  typename Mapping<dim,spacedim>::InternalDataBase *
  get_subface_data (const UpdateFlags flags,
                    const Quadrature<dim-1>& quadrature) const;

  /**
   * Computes the support points of
   * the mapping. For @p MappingQ1
   * these are the
   * vertices. However, other
   * classes may override this
   * function. In particular, the
   * MappingQ1Eulerian class does
   * exactly this by not computing
   * the support points from the
   * geometry of the current cell
   * but instead evaluating an
   * externally given displacement
   * field in addition to the
   * geometry of the cell.
   */
  virtual void compute_mapping_support_points(
    const typename Triangulation<dim,spacedim>::cell_iterator &cell,
    std::vector<Point<spacedim> > &a) const;

  /**
   * Number of shape functions. Is
   * simply the number of vertices
   * per cell for the Q1 mapping.
   */
  static const unsigned int n_shape_functions = GeometryInfo<dim>::vertices_per_cell;
};


// explicit specializations

template<>
Point<2>
MappingQ1<2,3>::
transform_real_to_unit_cell_internal
(const Triangulation<2,3>::cell_iterator &cell,
 const Point<3> &p,
 const Point<2> &initial_p_unit,
 InternalData    &mdata) const;

template<>
Point<1>
MappingQ1<1,2>::
transform_real_to_unit_cell_internal
(const Triangulation<1,2>::cell_iterator &cell,
 const Point<2> &p,
 const Point<1> &initial_p_unit,
 InternalData    &mdata) const;

template<>
Point<1>
MappingQ1<1,3>::
transform_real_to_unit_cell_internal
(const Triangulation<1,3>::cell_iterator &cell,
 const Point<3> &p,
 const Point<1> &initial_p_unit,
 InternalData    &mdata) const;


/**
 * In order to avoid creation of static MappingQ1 objects at several
 * places in the library (in particular in backward compatibility
 * functions), we define a static MappingQ1 objects once and for all
 * places where it is needed.
 */
template <int dim, int spacedim=dim>
struct StaticMappingQ1
{
  static MappingQ1<dim, spacedim> mapping;
};


/*@}*/

/*----------------------------------------------------------------------*/

#ifndef DOXYGEN

template<int dim, int spacedim>
inline
double
MappingQ1<dim,spacedim>::InternalData::shape (const unsigned int qpoint,
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
MappingQ1<dim,spacedim>::InternalData::shape (const unsigned int qpoint,
                                              const unsigned int shape_nr)
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_values.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_values.size()));
  return shape_values [qpoint*n_shape_functions + shape_nr];
}


template<int dim, int spacedim>
inline
Tensor<1,dim>
MappingQ1<dim,spacedim>::InternalData::derivative (const unsigned int qpoint,
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
MappingQ1<dim,spacedim>::InternalData::derivative (const unsigned int qpoint,
                                                   const unsigned int shape_nr)
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_derivatives.size()));
  return shape_derivatives [qpoint*n_shape_functions + shape_nr];
}


template <int dim, int spacedim>
inline
Tensor<2,dim>
MappingQ1<dim,spacedim>::InternalData::second_derivative (const unsigned int qpoint,
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
MappingQ1<dim,spacedim>::InternalData::second_derivative (const unsigned int qpoint,
                                                          const unsigned int shape_nr)
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_second_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_second_derivatives.size()));
  return shape_second_derivatives [qpoint*n_shape_functions + shape_nr];
}



template <int dim, int spacedim>
inline
bool
MappingQ1<dim,spacedim>::preserves_vertex_locations () const
{
  return true;
}

#endif // DOXYGEN

/* -------------- declaration of explicit specializations ------------- */


DEAL_II_NAMESPACE_CLOSE

#endif

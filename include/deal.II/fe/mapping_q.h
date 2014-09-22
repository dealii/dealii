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

#ifndef __deal2__mapping_q_h
#define __deal2__mapping_q_h


#include <deal.II/base/config.h>
#include <deal.II/base/table.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/manifold.h>

DEAL_II_NAMESPACE_OPEN

template <int dim, typename POLY> class TensorProductPolynomials;


/*!@addtogroup mapping */
/*@{*/

/**
 * Mapping class that uses Qp-mappings on boundary cells. The mapping shape
 * functions make use of tensor product polynomials with unit cell support
 * points equal to the points of the Gauss-Lobatto quadrature formula. These
 * points give a well-conditioned interpolation also for very high orders and
 * are therefore preferred over equidistant support points.
 *
 * For more details about Qp-mappings, see the `mapping' report at
 * <tt>deal.II/doc/reports/mapping_q/index.html</tt> in the `Reports'
 * section of `Documentation'.
 *
 * For more information about the <tt>spacedim</tt> template parameter
 * check the documentation of FiniteElement or the one of
 * Triangulation.
 *
 * @note Since the boundary description is closely tied to the unit cell
 * support points, new boundary descriptions need to explicitly use the
 * Gauss-Lobatto points.
 *
 * @author Ralf Hartmann, 2000, 2001, 2005; Guido Kanschat 2000, 2001
 */
template <int dim, int spacedim=dim>
class MappingQ : public MappingQ1<dim,spacedim>
{
public:
  /**
   * Constructor.  @p p gives the degree of mapping polynomials on boundary
   * cells.
   *
   * The second argument determines whether the higher order mapping should
   * also be used on interior cells. If its value is <code>false</code> (the
   * default), the a lower-order mapping is used in the interior. This is
   * sufficient for most cases where higher order mappings are only used to
   * better approximate the boundary. In that case, cells bounded by straight
   * lines are acceptable in the interior. However, there are cases where one
   * would also like to use a higher order mapping in the interior. The
   * MappingQEulerian class is one such case.
   */
  MappingQ (const unsigned int p,
            const bool use_mapping_q_on_all_cells = false);

  /**
   * Copy constructor. Performs a deep copy, i.e. duplicates what #tensor_pols
   * points to instead of simply copying the #tensor_pols pointer as done by a
   * default copy constructor.
   */
  MappingQ (const MappingQ<dim,spacedim> &mapping);

  /**
   * Destructor.
   */
  virtual ~MappingQ ();

  /**
   * Transforms the point @p p on the unit cell to the point @p p_real on the
   * real cell @p cell and returns @p p_real.
   */
  virtual Point<spacedim>
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
  transform (const VectorSlice<const std::vector<DerivativeForm<1, dim, spacedim> > >    input,
             VectorSlice<std::vector<Tensor<2,spacedim> > > output,
             const typename Mapping<dim,spacedim>::InternalDataBase &internal,
             const MappingType type) const;

  virtual
  void
  transform (const VectorSlice<const std::vector<Tensor<2, dim> > >     input,
             VectorSlice<std::vector<Tensor<2,spacedim> > >             output,
             const typename Mapping<dim,spacedim>::InternalDataBase &internal,
             const MappingType type) const;

  /**
   * Return the degree of the mapping, i.e. the value which was passed to the
   * constructor.
   */
  unsigned int get_degree () const;

  /**
   * Return a pointer to a copy of the present object. The caller of this copy
   * then assumes ownership of it.
   */
  virtual
  Mapping<dim,spacedim> *clone () const;

  /**
   * Storage for internal data of Q_degree transformation.
   */
  class InternalData : public MappingQ1<dim,spacedim>::InternalData
  {
  public:
    /**
     * Constructor.
     */
    InternalData (const unsigned int n_shape_functions);


    /**
     * Return an estimate (in bytes) or the memory consumption of this object.
     */
    virtual std::size_t memory_consumption () const;

    /**
     * Unit normal vectors. Used for the alternative computation of the normal
     * vectors. See doc of the @p alternative_normals_computation flag.
     *
     * Filled (hardcoded) once in @p get_face_data.
     */
    std::vector<std::vector<Point<dim> > > unit_normals;

    /**
     * Flag that is set by the <tt>fill_fe_[[sub]face]_values</tt> function.
     *
     * If this flag is @p true we are on an interior cell and the @p
     * mapping_q1_data is used.
     */
    bool use_mapping_q1_on_current_cell;

    /**
     * On interior cells @p MappingQ1 is used.
     */
    typename MappingQ1<dim,spacedim>::InternalData mapping_q1_data;
  };

protected:
  /**
   * Implementation of the interface in Mapping.
   */
  virtual void
  fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                  const Quadrature<dim>                                     &quadrature,
                  typename Mapping<dim,spacedim>::InternalDataBase          &mapping_data,
                  typename std::vector<Point<spacedim> >                    &quadrature_points,
                  std::vector<double>                                       &JxW_values,
                  std::vector<DerivativeForm<1,dim,spacedim> >       &jacobians,
                  std::vector<DerivativeForm<2,dim,spacedim> >       &jacobian_grads,
                  std::vector<DerivativeForm<1,spacedim,dim> >      &inverse_jacobians,
                  std::vector<Point<spacedim> >                             &cell_normal_vectors,
                  CellSimilarity::Similarity                           &cell_similarity) const ;

  /**
   * Implementation of the interface in Mapping.
   */
  virtual void
  fill_fe_face_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                       const unsigned int face_no,
                       const Quadrature<dim-1>& quadrature,
                       typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
                       typename std::vector<Point<spacedim> >        &quadrature_points,
                       std::vector<double>             &JxW_values,
                       typename std::vector<Tensor<1,spacedim> >        &exterior_form,
                       typename std::vector<Point<spacedim> >        &normal_vectors) const ;

  /**
   * Implementation of the interface in Mapping.
   */
  virtual void
  fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                          const unsigned int face_no,
                          const unsigned int sub_no,
                          const Quadrature<dim-1>& quadrature,
                          typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
                          typename std::vector<Point<spacedim> >        &quadrature_points,
                          std::vector<double>             &JxW_values,
                          typename std::vector<Tensor<1,spacedim> >        &exterior_form,
                          typename std::vector<Point<spacedim> >        &normal_vectors) const ;

  /**
   * For <tt>dim=2,3</tt>. Append the support points of all shape functions
   * located on bounding lines to the vector @p a. Points located on the line
   * but not on vertices are not included.
   *
   * Needed by the @p compute_support_points_laplace function . For
   * <tt>dim=1</tt> this function is empty.
   *
   * This function is made virtual in order to allow derived classes to choose
   * shape function support points differently than the present class, which
   * chooses the points as interpolation points on the boundary.
   */
  virtual void
  add_line_support_points (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                           std::vector<Point<spacedim> > &a) const;

  /**
   * For <tt>dim=3</tt>. Append the support points of all shape functions
   * located on bounding faces (quads in 3d) to the vector @p a. Points
   * located on the quad but not on vertices are not included.
   *
   * Needed by the @p compute_support_points_laplace function. For
   * <tt>dim=1</tt> and <tt>dim=2</tt> this function is empty.
   *
   * This function is made virtual in order to allow derived classes to choose
   * shape function support points differently than the present class, which
   * chooses the points as interpolation points on the boundary.
   */
  virtual void
  add_quad_support_points(const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                          std::vector<Point<spacedim> > &a) const;


private:
  /**
   * Ask the manifold descriptor to return intermediate points on
   * lines or faces. The function needs to return one or multiple
   * points (depending on the number of elements in the output vector
   * @p points that lie inside a line, quad or hex). Whether it is a
   * line, quad or hex doesn't really matter to this function but it
   * can be inferred from the number of input points in the @p
   * surrounding_points vector.
   */
  void get_intermediate_points(const Manifold<dim, spacedim> &manifold,
                               const std::vector<Point<spacedim> > &surrounding_points,
                               std::vector<Point<spacedim> > &points) const;


  /**
   * Ask the manifold descriptor to return intermediate points on the
   * object pointed to by the TriaIterator @p iter. This function
   * tries to be backward compatible with respect to the differences
   * between Boundary<dim,spacedim> and Manifold<dim,spacedim>,
   * querying the first whenever the passed @p manifold can be
   * upgraded to a Boundary<dim,spacedim>.
   */
  template <class TriaIterator>
  void get_intermediate_points_on_object(const Manifold<dim, spacedim> &manifold,
                                         const TriaIterator &iter,
                                         std::vector<Point<spacedim> > &points) const;


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
   * Compute shape values and/or derivatives.
   */
  virtual void
  compute_shapes_virtual (const std::vector<Point<dim> > &unit_points,
                          typename MappingQ1<dim,spacedim>::InternalData &data) const;

  /**
   * This function is needed by the constructor of
   * <tt>MappingQ<dim,spacedim></tt> for <tt>dim=</tt> 2 and 3.
   *
   * For <tt>degree<4</tt> this function sets the @p laplace_on_quad_vector to
   * the hardcoded data. For <tt>degree>=4</tt> and MappingQ<2> this vector is
   * computed.
   *
   * For the definition of the @p laplace_on_quad_vector please refer to
   * equation (8) of the `mapping' report.
   */
  void
  set_laplace_on_quad_vector(Table<2,double> &loqvs) const;

  /**
   * This function is needed by the constructor of <tt>MappingQ<3></tt>.
   *
   * For <tt>degree==2</tt> this function sets the @p laplace_on_hex_vector to
   * the hardcoded data. For <tt>degree>2</tt> this vector is computed.
   *
   * For the definition of the @p laplace_on_hex_vector please refer to
   * equation (8) of the `mapping' report.
   */
  void set_laplace_on_hex_vector(Table<2,double> &lohvs) const;

  /**
   * Computes the <tt>laplace_on_quad(hex)_vector</tt>.
   *
   * Called by the <tt>set_laplace_on_quad(hex)_vector</tt> functions if the
   * data is not yet hardcoded.
   *
   * For the definition of the <tt>laplace_on_quad(hex)_vector</tt> please
   * refer to equation (8) of the `mapping' report.
   */
  void compute_laplace_vector(Table<2,double> &lvs) const;

  /**
   * Takes a <tt>laplace_on_hex(quad)_vector</tt> and applies it to the vector
   * @p a to compute the inner support points as a linear combination of the
   * exterior points.
   *
   * The vector @p a initially contains the locations of the @p n_outer
   * points, the @p n_inner computed inner points are appended.
   *
   * See equation (7) of the `mapping' report.
   */
  void apply_laplace_vector(const Table<2,double>   &lvs,
                            std::vector<Point<spacedim> > &a) const;

  /**
   * Computes the support points of the mapping.
   */
  virtual void compute_mapping_support_points(
    const typename Triangulation<dim,spacedim>::cell_iterator &cell,
    std::vector<Point<spacedim> > &a) const;

  /**
   * Computes all support points of the mapping shape functions. The inner
   * support points (ie. support points in quads for 2d, in hexes for 3d) are
   * computed using the solution of a Laplace equation with the position of
   * the outer support points as boundary values, in order to make the
   * transformation as smooth as possible.
   */
  void compute_support_points_laplace(
    const typename Triangulation<dim,spacedim>::cell_iterator &cell,
    std::vector<Point<spacedim> > &a) const;

  /**
   * Needed by the @p laplace_on_quad function (for <tt>dim==2</tt>). Filled
   * by the constructor.
   *
   * Sizes:
   * laplace_on_quad_vector.size()= number of inner unit_support_points
   * laplace_on_quad_vector[i].size()= number of outer unit_support_points,
   *   i.e.  unit_support_points on the boundary of the quad
   *
   * For the definition of this vector see equation (8) of the `mapping'
   * report.
   */
  Table<2,double> laplace_on_quad_vector;

  /**
   * Needed by the @p laplace_on_hex function (for <tt>dim==3</tt>). Filled by
   * the constructor.
   *
   * For the definition of this vector see equation (8) of the `mapping'
   * report.
   */
  Table<2,double> laplace_on_hex_vector;

  /**
   * Exception.
   */
  DeclException1 (ExcLaplaceVectorNotSet,
                  int,
                  << "laplace_vector not set for degree=" << arg1 << ".");

  /**
   * Degree @p p of the polynomials used as shape functions for the Qp mapping
   * of cells at the boundary.
   */
  const unsigned int degree;

  /**
   * Number of inner mapping shape functions.
   */
  const unsigned int n_inner;

  /**
   * Number of mapping shape functions on the boundary.
   */
  const unsigned int n_outer;

  /**
   * Pointer to the @p dim-dimensional tensor product polynomials used as
   * shape functions for the Qp mapping of cells at the boundary.
   */
  const TensorProductPolynomials<dim> *tensor_pols;

  /**
   * Number of the Qp tensor product shape functions.
   */
  const unsigned int n_shape_functions;

  /**
   * Mapping from lexicographic to to the Qp shape function numbering. Its
   * size is @p dofs_per_cell.
   */
  const std::vector<unsigned int> renumber;

  /**
   * If this flag is set @p true then @p MappingQ is used on all cells, not
   * only on boundary cells.
   */
  const bool use_mapping_q_on_all_cells;

  /**
   * An FE_Q object which is only needed in 3D, since it knows how to reorder
   * shape functions/DoFs on non-standard faces. This is used to reorder
   * support points in the same way. We could make this a pointer to prevent
   * construction in 1D and 2D, but since memory and time requirements are not
   * particularly high this seems unnecessary at the moment.
   */
  const FE_Q<dim> feq;


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
   * Declare other MappingQ classes friends.
   */
  template <int,int> friend class MappingQ;
};

/*@}*/

/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN

template<> MappingQ<1>::MappingQ (const unsigned int,
                                  const bool);
template<> MappingQ<1>::~MappingQ ();

template<>
void MappingQ<1>::compute_shapes_virtual (const std::vector<Point<1> > &unit_points,
                                          MappingQ1<1>::InternalData   &data) const;
template <>
void MappingQ<1>::set_laplace_on_quad_vector(Table<2,double> &) const;

template <>
void MappingQ<3>::set_laplace_on_hex_vector(Table<2,double> &lohvs) const;

template <>
void MappingQ<1>::compute_laplace_vector(Table<2,double> &) const;


template<>
void MappingQ<3>::add_quad_support_points(const Triangulation<3>::cell_iterator &cell,
                                          std::vector<Point<3> >                &a) const;

// ---- Templated functions ---- //
template <int dim, int spacedim>
template <class TriaIterator>
void
MappingQ<dim,spacedim>::get_intermediate_points_on_object(const Manifold<dim, spacedim> &manifold,
                                                          const TriaIterator &iter,
                                                          std::vector<Point<spacedim> > &points) const
{
  const unsigned int structdim = TriaIterator::AccessorType::structure_dimension;
  // Try backward compatibility option.
  const Boundary<dim,spacedim> *boundary = dynamic_cast<const Boundary<dim,spacedim> *>(&manifold);
  if (boundary) // This is actually a boundary. Call old methods.
    switch (structdim)
      {
      case 1:
      {
        const typename Triangulation<dim,spacedim>::line_iterator line = iter;
        boundary->get_intermediate_points_on_line(line, points);
        return;
      }
      case 2:
      {
        const typename Triangulation<dim,spacedim>::quad_iterator quad = iter;
        boundary->get_intermediate_points_on_quad(quad, points);
        return;
      }
      default:
        Assert(false, ExcInternalError());
        return;
      }
  else
    {
      std::vector<Point<spacedim> > sp(GeometryInfo<structdim>::vertices_per_cell);
      for (unsigned int i=0; i<sp.size(); ++i)
        sp[i] = iter->vertex(i);
      get_intermediate_points(manifold, sp, points);
    }
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2015 by the deal.II authors
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

#ifndef dealii__mapping_h
#define dealii__mapping_h


#include <deal.II/base/config.h>
#include <deal.II/base/derivative_form.h>
#include <deal.II/base/std_cxx11/array.h>
#include <deal.II/base/vector_slice.h>
#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe_update_flags.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

template <int dim> class Quadrature;
template <int dim, int spacedim> class FEValues;
template <int dim, int spacedim> class FEValuesBase;
template <int dim, int spacedim> class FEValues;
template <int dim, int spacedim> class FEFaceValues;
template <int dim, int spacedim> class FESubfaceValues;


/**
 * The transformation type used for the Mapping::transform() functions.
 *
 * Special finite elements may need special Mapping from the reference cell to
 * the actual mesh cell. In order to be most flexible, this enum provides an
 * extensible interface for arbitrary transformations. Nevertheless, these
 * must be implemented in the transform() functions of inheriting classes in
 * order to work.
 *
 * @ingroup mapping
 */
enum MappingType
{
  /**
   * No mapping, i.e., shape functions are not mapped from a reference cell
   * but instead are defined right on the real-space cell.
   */
  mapping_none = 0x0000,

  /**
   * Covariant mapping (see Mapping::transform() for details).
   */
  mapping_covariant = 0x0001,

  /**
   * Contravariant mapping (see Mapping::transform() for details).
   */
  mapping_contravariant = 0x0002,

  /**
   * Mapping of the gradient of a covariant vector field (see Mapping::transform() for details).
   */
  mapping_covariant_gradient = 0x0003,

  /**
   * Mapping of the gradient of a contravariant vector field (see Mapping::transform() for details).
   */
  mapping_contravariant_gradient = 0x0004,

  /**
   * The Piola transform usually used for Hdiv elements. Piola transform is
   * the standard transformation of vector valued elements in H<sup>div</sup>.
   * It amounts to a contravariant transformation scaled by the inverse of the
   * volume element.
   */
  mapping_piola = 0x0100,

  /**
   * Transformation for the gradient of a vector field corresponding to a
   * mapping_piola transformation (see Mapping::transform() for details).
   */
  mapping_piola_gradient = 0x0101,

  /**
   * The mapping used for Nedelec elements.
   *
   * Curl-conforming elements are mapped as covariant vectors. Nevertheless,
   * we introduce a separate mapping type, such that we can use the same flag
   * for the vector and its gradient (see Mapping::transform() for details).
   */
  mapping_nedelec = 0x0200,

  /**
   * The mapping used for Raviart-Thomas elements.
   */
  mapping_raviart_thomas = 0x0300,

  /**
   * The mapping used for BDM elements.
   */
  mapping_bdm = mapping_raviart_thomas
};


/**
 * Abstract base class for mapping classes.
 *
 * The interface for filling the tables of FEValues is provided. Everything
 * else has to happen in derived classes.
 *
 * The following paragraph applies to the implementation of FEValues. Usage of
 * the class is as follows: first, call the functions @p update_once and @p
 * update_each with the update flags you need. This includes the flags needed
 * by the FiniteElement. Then call <tt>get_*_data</tt> and with the or'd
 * results.  This will initialize and return some internal data structures. On
 * the first cell, call <tt>fill_fe_*_values</tt> with the result of @p
 * update_once. Finally, on each cell, use <tt>fill_fe_*_values</tt> with the
 * result of @p update_each to compute values for a special cell.
 *
 * <h3>Mathematics of the mapping</h3>
 *
 * The mapping is a transformation $\mathbf x = \Phi(\mathbf{\hat x})$ which
 * maps the reference cell [0,1]<sup>dim</sup> to the actual grid cell in
 * R<sup>spacedim</sup>. In order to describe the application of the mapping
 * to different objects, we introduce the notation for the Jacobian
 * $J(\mathbf{\hat x}) = \nabla\Phi(\mathbf{\hat x})$. For instance, if
 * dim=spacedim=2, we have
 * @f[
 * J(\mathbf{\hat x}) = \left(\begin{matrix}
 * \frac{\partial x}{\partial \hat x} & \frac{\partial x}{\partial \hat y}
 * \\
 * \frac{\partial y}{\partial \hat x} & \frac{\partial y}{\partial \hat y}
 * \end{matrix}\right)
 * @f]
 *
 * <h4>Mapping of functions</h4>
 *
 * Functions are simply mapped such that
 * @f[
 * u(\mathbf x) = u\bigl(\Phi(\mathbf{\hat x})\bigr)
 * = \hat u(\mathbf{\hat x}).
 * @f]
 * Since finite element shape functions are usually defined on the reference
 * cell, nothing needs to be done for them. For a function defined on the
 * computational domain, the quadrature points need to be mapped, which is
 * done in fill_fe_values() if @p update_quadrature_points is set in the
 * update flags. The mapped quadrature points are then accessed through
 * FEValuesBase::quadrature_point().
 *
 * @todo Add a function <tt>transform_quadrature_points</tt> for this.
 *
 * <h4>Mapping of integrals</h4>
 *
 * The volume form $d\hat x$ is mapped such that for a grid cell <i>Z</i>
 * @f[
 *  \int_Z u(\mathbf x)\,d\mathbf x = \int_{\hat Z} \hat
 * u(\mathbf{\hat x}) \left|\text{det}J(\mathbf{\hat x})\right|
 * \,d\mathbf{\hat x}.
 * @f]
 *
 * The transformed quadrature weights $\left|\text{det}J(\mathbf{\hat
 * x})\right|$ are accessed through FEValuesBase::JxW() and computed in
 * fill_fe_values(), if @p update_JxW_values is set in the update flags.
 *
 * @todo Add a function <tt>transform_quadrature_weights</tt> for this.
 *
 * @todo Add documentation on the codimension-one case
 *
 * <h4>Mapping of vector fields, differential forms and gradients of vector
 * fields</h4>
 *
 * The transformation of vector fields, differential forms
 * (gradients/jacobians) and gradients of vector fields between the reference
 * cell and the actual grid cell follows the general form
 *
 * @f[
 * \mathbf v(\mathbf x) = \mathbf A(\mathbf{\hat x})
 * \mathbf{\hat v}(\mathbf{\hat x}),
 * \qquad
 * \mathbf T(\mathbf x) = \mathbf A(\mathbf{\hat x})
 * \mathbf{\hat T}(\mathbf{\hat x}) \mathbf B(\mathbf{\hat x}),
 * @f]
 *
 * where <b>v</b> is a vector field or a differential form and and <b>T</b> a
 * tensor field of gradients. The differential forms <b>A</b> and <b>B</b> are
 * determined by the MappingType enumerator. These transformations are
 * performed through the functions transform(). See the documentation there
 * for possible choices.
 *
 * <h3>Technical notes</h3>
 *
 * A hint to implementors: no function except the two functions @p
 * update_once and @p update_each may add any flags.
 *
 * For more information about the <tt>spacedim</tt> template parameter check
 * the documentation of FiniteElement or the one of Triangulation.
 *
 * <h3>References</h3>
 *
 * A general publication on differential geometry and finite elements is the
 * survey
 * <ul>
 * <li>Douglas N. Arnold, Richard S. Falk, and Ragnar Winther. <i>Finite
 * element exterior calculus: from Hodge theory to numerical stability.</i>
 * Bull. Amer. Math. Soc. (N.S.), 47:281-354, 2010. <a
 * href="http://dx.doi.org/10.1090/S0273-0979-10-01278-4">DOI:
 * 10.1090/S0273-0979-10-01278-4</a>.
 * </ul>
 *
 * The description of the Piola transform has been taken from the <a
 * href="http://www.math.uh.edu/~rohop/spring_05/downloads/">lecture notes</a>
 * by Ronald H. W. Hoppe, University of Houston, Chapter 7.
 *
 * @ingroup mapping
 * @author Guido Kanschat, Ralf Hartmann 2000, 2001
 */
template <int dim, int spacedim=dim>
class Mapping : public Subscriptor
{
public:

  /**
   * Virtual destructor.
   */
  virtual ~Mapping ();

  /**
   * Return the mapped vertices of a cell. These values are not equal to the
   * vertex coordinates stored by the triangulation for MappingQEulerian and
   * MappingQ1Eulerian.
   */
  virtual
  std_cxx11::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
  get_vertices (
    const typename Triangulation<dim,spacedim>::cell_iterator &cell) const;

  /**
   * Transforms the point @p p on the unit cell to the point @p p_real on the
   * real cell @p cell and returns @p p_real.
   */
  virtual Point<spacedim>
  transform_unit_to_real_cell (
    const typename Triangulation<dim,spacedim>::cell_iterator &cell,
    const Point<dim>                                 &p) const = 0;

  /**
   * Transforms the point @p p on the real @p cell to the corresponding point
   * on the unit cell, and return its coordinates.
   *
   * In the codimension one case, this function returns the normal projection
   * of the real point @p p on the curve or surface identified by the @p cell.
   *
   * @note Polynomial mappings from the reference (unit) cell coordinates to
   * the coordinate system of a real cell are not always invertible if the
   * point for which the inverse mapping is to be computed lies outside the
   * cell's boundaries. In such cases, the current function may fail to
   * compute a point on the reference cell whose image under the mapping
   * equals the given point @p p.  If this is the case then this function
   * throws an exception of type Mapping::ExcTransformationFailed . Whether
   * the given point @p p lies outside the cell can therefore be determined by
   * checking whether the return reference coordinates lie inside of outside
   * the reference cell (e.g., using GeometryInfo::is_inside_unit_cell) or
   * whether the exception mentioned above has been thrown.
   */
  virtual Point<dim>
  transform_real_to_unit_cell (
    const typename Triangulation<dim,spacedim>::cell_iterator &cell,
    const Point<spacedim>                            &p) const = 0;

  /**
   * Transforms the point @p p on the real @p cell to the corresponding point
   * on the unit cell, and then projects it to a dim-1  point on the face with
   * the given face number @p face_no. Ideally the point @p p is near the face
   * @ face_no, but any point in the cell can technically be projected.
   *
   * This function does not make physical sense when dim=1,
   * so it throws an exception in this case.
   */
  virtual Point<dim-1>
  transform_real_to_unit_projected_to_face (
    const typename Triangulation<dim,spacedim>::cell_iterator &cell,
    const unsigned int &face_no,
    const Point<spacedim> &p) const;

  /**
   * Base class for internal data of mapping objects. The
   * internal mechanism is that upon construction of a FEValues object, it
   * asks the mapping and finite element classes that are to be used to
   * allocate memory for their own purpose in which they may store data that
   * only needs to be computed once. For example, most finite elements will
   * store the values of the shape functions at the quadrature points in this
   * object, since they do not change from cell to cell and only need to be
   * computed once. The same may be true for Mapping classes that want to
   * only evaluate the shape functions used for mapping once at the quadrature
   * points.
   *
   * Since different FEValues objects using different
   * quadrature rules might access the same mapping object at the same
   * time, it is necessary to create one such object per FEValues object.
   * FEValues does this by calling Mapping::get_data(), or in reality the
   * implementation of the corresponding function in derived classes.
   * Ownership of the object created by Mapping::get_data() is then transferred
   * to the FEValues object,
   * but a reference to this object is passed to the mapping object every
   * time it is asked to compute information on a concrete cell. This
   * happens when FEValues::reinit() (or the corresponding classes in
   * FEFaceValues and FESubfaceValues) call Mapping::fill_fe_values()
   * (and similarly via Mapping::fill_fe_face_values() and
   * Mapping::fill_fe_subface_values()).
   *
   * The purpose of this class is for mapping objects to store information
   * that can be computed once at the beginning, on the reference cell,
   * and to access it later when computing information on a concrete cell.
   * As such, the object handed to Mapping::fill_fe_values() is marked as
   * <code>const</code>, because the assumption is that at the time this
   * information is used, it will not need to modified again. However,
   * classes derived from Mapping can also use such objects for two other
   * purposes:
   *
   * - To provide scratch space for computations that are done in
   *   Mapping::fill_fe_values() and similar functions. Some of the
   *   derived classes would like to use scratch arrays and it would
   *   be a waste of time to allocate these arrays every time this
   *   function is called, just to de-allocate it again at the end
   *   of the function. Rather, one could allocate this memory once
   *   as a member variable of the current class, and simply use
   *   it in Mapping::fill_fe_values().
   * - After calling Mapping::fill_fe_values(), FEValues::reinit()
   *   calls FiniteElement::fill_fe_values() where the finite element
   *   computes values, gradients, etc of the shape functions using
   *   both information computed once at the beginning using a mechanism
   *   similar to the one described here (see FiniteElement::InternalDataBase)
   *   as well as the data already computed by Mapping::fill_fe_values().
   *   As part of its work, some implementations of
   *   FiniteElement::fill_fe_values() need to transform shape function
   *   data, and they do so by calling Mapping::transform(). The call
   *   to the latter function also receives a reference to the
   *   Mapping::InternalDataBase object. Since Mapping::transform()
   *   may be called many times on each cell, it is sometimes worth
   *   for derived classes to compute some information only once
   *   in Mapping::fill_fe_values() and reuse it in
   *   Mapping::transform(). This information can also be stored in
   *   the classes that derived mapping classes derive from
   *   InternalDataBase.
   *
   * In both of these cases, the InternalDataBase object being passed
   * around is "morally const", i.e., no external observer can tell
   * whether a scratch array or some intermediate data for
   * Mapping::transform() is being modified by Mapping::fill_fe_values()
   * or not. Consequently, the InternalDataBase objects are always
   * passed around as <code>const</code> objects. Derived classes
   * that would like to make use of the two additional uses outlined
   * above therefore need to mark the member variables they want to
   * use for these purposes as <code>mutable</code> to allow for their
   * modification despite the fact that the surrounding object is
   * marked as <code>const</code>.
   */
  class InternalDataBase: public Subscriptor
  {
  private:
    /**
     * Copy constructor forbidden.
     */
    InternalDataBase (const InternalDataBase &);

  public:
    /**
     * Constructor. Sets update_flags to @p update_default and @p first_cell
     * to @p true.
     */
    InternalDataBase ();

    /**
     * Virtual destructor for derived classes
     */
    virtual ~InternalDataBase ();

    /**
     * Values updated by the constructor or by reinit.
     */
    UpdateFlags          update_flags;

    /**
     * Values computed by constructor.
     */
    UpdateFlags          update_once;

    /**
     * Values updated on each cell by reinit.
     */
    UpdateFlags          update_each;

    /**
     * If <tt>first_cell==true</tt> this function returns @p update_flags,
     * i.e. <tt>update_once|update_each</tt>. If <tt>first_cell==false</tt> it
     * returns @p update_each.
     */
    UpdateFlags  current_update_flags() const;

    /**
     * Return whether we are presently initializing data for the first cell.
     * The value of the field this function is returning is set to @p true in
     * the constructor, and cleared by the @p FEValues class after the first
     * cell has been initialized.
     *
     * This function is used to determine whether we need to use the @p
     * update_once flags for computing data, or whether we can use the @p
     * update_each flags.
     */
    bool is_first_cell () const;

    /**
     * Set the @p first_cell flag to @p false. Used by the @p FEValues class
     * to indicate that we have already done the work on the first cell.
     */
    virtual void clear_first_cell ();

    /**
     * Return an estimate (in bytes) or the memory consumption of this object.
     */
    virtual std::size_t memory_consumption () const;

    /**
     * The positions of the mapped (generalized) support points.
     */
    std::vector<Point<spacedim> > support_point_values;

    /*
     * The Jacobian of the
     * transformation in the
     * (generalized) support
     * points.
     */
    std::vector<Tensor<2,spacedim> > support_point_gradients;

    /*
     * The inverse of the
     * Jacobian of the
     * transformation in the
     * (generalized) support
     * points.
     */
    std::vector<Tensor<2,spacedim> > support_point_inverse_gradients;

  private:
    /**
     * The value returned by @p is_first_cell.
     */
    bool first_cell;
  };


  /**
   * Transform a field of vectors or 1-differential forms according to the
   * selected MappingType.
   *
   * @note Normally, this function is called by a finite element, filling
   * FEValues objects. For this finite element, there should be an alias
   * MappingType like @p mapping_bdm, @p mapping_nedelec, etc. This alias
   * should be preferred to using the types below.
   *
   * The mapping types currently implemented by derived classes are:
   * <ul>
   * <li> @p mapping_contravariant: maps a vector field on the reference cell
   * is to the physical cell through the Jacobian:
   * @f[
   * \mathbf u(\mathbf x) = J(\mathbf{\hat x})\mathbf{\hat u}(\mathbf{\hat x}).
   * @f]
   * In physics, this is usually referred to as the contravariant
   * transformation. Mathematically, it is the push forward of a vector field.
   *
   * <li> @p mapping_covariant: maps a field of one-forms on the reference
   * cell to a field of one-forms on the physical cell. (theoretically this
   * would refer to a DerivativeForm<1, dim, 1> but it canonically identified
   * with a Tensor<1,dim>). Mathematically, it is the pull back of the
   * differential form
   * @f[
   * \mathbf u(\mathbf x) = J(J^{T} J)^{-1}(\mathbf{\hat x})\mathbf{\hat
   * u}(\mathbf{\hat x}).
   * @f]
   * In the case when dim=spacedim the previous formula reduces to
   * @f[
   * \mathbf u(\mathbf x) = J^{-T}(\mathbf{\hat x})\mathbf{\hat
   * u}(\mathbf{\hat x}).
   * @f]
   * Gradients of scalar differentiable functions are transformed this way.
   *
   * <li> @p mapping_piola: A field of <i>n-1</i>-forms on the reference cell
   * is also represented by a vector field, but again transforms differently,
   * namely by the Piola transform
   * @f[
   *  \mathbf u(\mathbf x) = \frac{1}{\text{det}J(\mathbf x)}
   * J(\mathbf x) \mathbf{\hat u}(\mathbf x).
   * @f]
   * </ul>
   *
   * @todo What is n in mapping_piola description?
   */
  virtual
  void
  transform (const VectorSlice<const std::vector<Tensor<1,dim> > > input,
             VectorSlice<std::vector<Tensor<1,spacedim> > >        output,
             const InternalDataBase &internal,
             const MappingType type) const = 0;



  /**
   * Transform a field of differential forms from the reference cell to the
   * physical cell.  It is useful to think of $\mathbf{T} = D \mathbf u$ and
   * $\mathbf{\hat T} = \hat D \mathbf{\hat u}$, with $\mathbf u$ a vector
   * field.  The mapping types currently implemented by derived classes are:
   * <ul>
   * <li> @p mapping_covariant: maps a field of forms on the reference cell to
   * a field of forms on the physical cell. Mathematically, it is the pull
   * back of the differential form
   * @f[
   * \mathbf T(\mathbf x) = \mathbf{\hat T}(\mathbf{\hat x})
   *                        J*(J^{T} J)^{-1}(\mathbf{\hat x}).
   * @f]
   * n the case when dim=spacedim the previous formula reduces to
   * @f[
   * \mathbf T(\mathbf x) = \mathbf{\hat u}(\mathbf{\hat x})
   *                        J^{-1}(\mathbf{\hat x}).
   * @f]
   * Jacobians of spacedim-vector valued differentiable functions are
   * transformed this way.
   * </ul>
   * @note It would have been more reasonable to make this transform a
   * template function with the rank in <code>DerivativeForm@<1, dim,
   * rank@></code>. Unfortunately C++ does not allow templatized virtual
   * functions. This is why we identify <code>DerivativeForm@<1, dim,
   * 1@></code> with a <code>Tensor@<1,dim@></code> when using
   * mapping_covariant() in the function transform above this one.
   */
  virtual
  void
  transform (const VectorSlice<const std::vector< DerivativeForm<1, dim, spacedim> > > input,
             VectorSlice<std::vector<Tensor<2,spacedim> > >             output,
             const InternalDataBase &internal,
             const MappingType type) const = 0;



  /**
   * Transform a tensor field from the reference cell to the physical cell.
   * This tensors are most of times the jacobians in the reference cell of
   * vector fields that have been pulled back from the physical cell.  The
   * mapping types currently implemented by derived classes are:
   * <ul>
   * <li> @p mapping_contravariant_gradient, it assumes $\mathbf u(\mathbf x)
   * = J \mathbf{\hat u}$ so that
   * @f[
   * \mathbf T(\mathbf x) =
   * J(\mathbf{\hat x}) \mathbf{\hat T}(\mathbf{\hat x})
   * J^{-1}(\mathbf{\hat x}).
   * @f]
   * <li> @p mapping_covariant_gradient, it assumes $\mathbf u(\mathbf x) =
   * J^{-T} \mathbf{\hat u}$ so that
   * @f[
   * \mathbf T(\mathbf x) =
   * J^{-T}(\mathbf{\hat x}) \mathbf{\hat T}(\mathbf{\hat x})
   * J^{-1}(\mathbf{\hat x}).
   * @f]
   * <li> @p mapping_piola_gradient, it assumes $\mathbf u(\mathbf x) =
   * \frac{1}{\text{det}J(\mathbf x)} J(\mathbf x) \mathbf{\hat u}(\mathbf x)$
   * so that
   * @f[
   * \mathbf T(\mathbf x) =
   * \frac{1}{\text{det}J(\mathbf x)}
   * J(\mathbf{\hat x}) \mathbf{\hat T}(\mathbf{\hat x})
   * J^{-1}(\mathbf{\hat x}).
   * @f]
   * </ul>
   * @todo The formulas for mapping_covariant_gradient(),
   * mapping_contravariant_gradient() and mapping_piola_gradient() are only
   * true as stated for linear mappings. If, for example, the mapping is
   * bilinear then there is a missing term associated with the derivative of
   * J.
   */
  virtual
  void
  transform (const VectorSlice<const std::vector<Tensor<2, dim> > >     input,
             VectorSlice<std::vector<Tensor<2,spacedim> > >             output,
             const InternalDataBase &internal,
             const MappingType type) const = 0;

  /**
   * The transformed (generalized) support point.
   */
  const Point<spacedim> &support_point_value(
    const unsigned int index,
    const typename Mapping<dim,spacedim>::InternalDataBase &internal) const;

  /**
   * The Jacobian matrix of the transformation in the (generalized) support
   * point.
   */
  const Tensor<2,spacedim> &support_point_gradient(
    const unsigned int index,
    const typename Mapping<dim,spacedim>::InternalDataBase &internal) const;

  /**
   * The inverse Jacobian matrix of the transformation in the (generalized)
   * support point.
   */
  const Tensor<2,spacedim> &support_point_inverse_gradient(
    const unsigned int index,
    const typename Mapping<dim,spacedim>::InternalDataBase &internal) const;

  /**
   * Return a pointer to a copy of the present object. The caller of this copy
   * then assumes ownership of it.
   *
   * Since one can't create objects of class Mapping, this function of course
   * has to be implemented by derived classes.
   *
   * This function is mainly used by the hp::MappingCollection class.
   */
  virtual
  Mapping<dim,spacedim> *clone () const = 0;

  /**
   * Returns whether the mapping preserves vertex locations, i.e. whether the
   * mapped location of the reference cell vertices (given by
   * GeometryInfo::unit_cell_vertex()) equals the result of
   * <code>cell-@>vertex()</code>.
   *
   * For example, implementations in derived classes return @p true for
   * MappingQ, MappingQ1, MappingCartesian, but @p false for MappingQEulerian,
   * MappingQ1Eulerian.
   */
  virtual
  bool preserves_vertex_locations () const = 0;

  /**
   * Exception
   */
  DeclException0 (ExcInvalidData);


  /**
   * Computing the mapping between a real space point and a point in reference
   * space failed, typically because the given point lies outside the cell
   * where the inverse mapping is not unique.
   *
   * @ingroup Exceptions
   */
  DeclException0(ExcTransformationFailed);

  /**
   * deal.II assumes the Jacobian determinant to be positive. When the cell
   * geometry is distorted under the image of the mapping, the mapping becomes
   * invalid and this exception is thrown.
   *
   * @ingroup Exceptions
   */
  DeclException3 (ExcDistortedMappedCell,
                  Point<spacedim>, double, int,
                  << "The image of the mapping applied to cell with center ["
                  << arg1 << "] is distorted. The cell geometry or the "
                  << "mapping are invalid, giving a non-positive volume "
                  << "fraction of " << arg2 << " in quadrature point "
                  << arg3 << ".");

private:

  /**
   * Indicate fields to be updated in the constructor of FEValues. Especially,
   * fields not asked for by FEValues, but computed for efficiency reasons
   * will be notified here.
   *
   * See
   * @ref UpdateFlagsEssay.
   */
  virtual UpdateFlags update_once (const UpdateFlags) const = 0;

  /**
   * The same as update_once(), but for the flags to be updated for each grid
   * cell.
   *
   * See
   * @ref UpdateFlagsEssay.
   */
  virtual UpdateFlags update_each (const UpdateFlags) const = 0;

  /**
   * Prepare internal data structures and fill in values independent of the
   * cell. See the documentation of Mapping::InternalDataBase for more
   * information on the purpose of this function.
   */
  virtual InternalDataBase *
  get_data (const UpdateFlags,
            const Quadrature<dim> &quadrature) const = 0;

  /**
   * Prepare internal data structure for transformation of faces and fill in
   * values independent of the cell. See the documentation of
   * Mapping::InternalDataBase for more
   * information on the purpose of this function.
   */
  virtual InternalDataBase *
  get_face_data (const UpdateFlags flags,
                 const Quadrature<dim-1>& quadrature) const = 0;

  /**
   * Prepare internal data structure for transformation of children of faces
   * and fill in values independent of the cell. See the documentation
   * of Mapping::InternalDataBase for more
   * information on the purpose of this function.
   */
  virtual InternalDataBase *
  get_subface_data (const UpdateFlags flags,
                    const Quadrature<dim-1>& quadrature) const = 0;

  /**
   * Compute information about the mapping from the reference cell
   * to the real cell indicated by the first argument to this function.
   * Derived classes will have to implement this function based on the
   * kind of mapping they represent. It is called by FEValues::reinit().
   *
   * Conceptually, this function's represents the application of the
   * mapping $\mathbf x=\mathbf F_K(\hat {\mathbf x})$ from reference
   * coordinates $\mathbf\in [0,1]^d$ to real space coordinates
   * $\mathbf x$ for a given cell $K$. Its purpose is to compute the following
   * kinds of data:
   *
   * - Data that results from the application of the mapping itself, e.g.,
   *   computing the location $\mathbf x_q = \mathbf F_K(\hat{\mathbf x}_q)$
   *   of quadrature points on the real cell, and that is directly useful
   *   to users of FEValues, for example during assembly.
   * - Data that is necessary for finite element implementations to compute
   *   their shape functions on the real cell. To this end, the
   *   FEValues::reinit() function calls FiniteElement::fill_fe_values()
   *   after the current function, and the output of this function serves
   *   as input to FiniteElement::fill_fe_values(). Examples of
   *   information that needs to be computed here for use by the
   *   finite element classes is the Jacobian of the mapping,
   *   $\hat\nabla \mathbf F_K(\hat{\mathbf x})$ or its inverse,
   *   for example to transform the gradients of shape functions on
   *   the reference cell to the gradients of shape functions on
   *   the real cell.
   *
   * The information computed by this function is used to fill the various
   * member variables of the output argument of this function. Which of
   * the member variables of that structure should be filled is determined
   * by the update flags stored in the Mapping::InternalDataBase object
   * passed to this function.
   *
   * @param[in] cell The cell of the triangulation for which this function
   *   is to compute a mapping from the reference cell to.
   * @param[in] cell_similarity Whether or not the cell given as first
   *   argument is simply a translation, rotation, etc of the cell for
   *   which this function was called the most recent time. This
   *   information is computed simply by matching the vertices (as stored
   *   by the Triangulation) between the previous and the current cell.
   *   The value passed here may be modified by implementations of
   *   this function and should then be returned (see the discussion of the
   *   return value of this function).
   * @param[in] quadrature A reference to the quadrature formula in use
   *   for the current evaluation. This quadrature object is the same
   *   as the one used when creating the @p internal_data object. The
   *   object is used both to map the location of quadrature points,
   *   as well as to compute the JxW values for each quadrature
   *   point (which involves the quadrature weights).
   * @param[in] internal_data A reference to an object previously
   *   created by get_data() and that may be used to store information
   *   the mapping can compute once on the reference cell. See the
   *   documentation of the Mapping::InternalDataBase class for an
   *   extensive description of the purpose of these objects.
   * @param[out] output_data A reference to an object whose member
   *   variables should be computed. Not all of the members of this
   *   argument need to be filled; which ones need to be filled is
   *   determined by the update flags stored inside the
   *   @p internal_data object.
   * @return An updated value of the @p cell_similarity argument to
   *   this function. The returned value will be used for the corresponding
   *   argument when FEValues::reinit() calls
   *   FiniteElement::fill_fe_values(). In most cases, derived classes will
   *   simply want to return the value passed for @p cell_similarity.
   *   However, implementations of this function may downgrade the
   *   level of cell similarity. This is, for example, the case for
   *   classes that take not only into account the locations of the
   *   vertices of a cell (as reported by the Triangulation), but also
   *   other information specific to the mapping. The purpose is that
   *   FEValues::reinit() can compute whether a cell is similar to the
   *   previous one only based on the cell's vertices, whereas the
   *   mapping may also consider displacement fields (e.g., in the
   *   MappingQ1Eulerian and MappingFEField classes). In such cases,
   *   the mapping may conclude that the previously computed
   *   cell similarity is too optimistic, and invalidate it for
   *   subsequent use in FiniteElement::fill_fe_values() by
   *   returning a less optimistic cell similarity value.
   */
  virtual
  CellSimilarity::Similarity
  fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                  const CellSimilarity::Similarity                           cell_similarity,
                  const Quadrature<dim>                                     &quadrature,
                  const InternalDataBase                                    &internal_data,
                  internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const = 0;

  /**
   * This function is the equivalent to Mapping::fill_fe_values(),
   * but for faces of cells. See there for an extensive discussion
   * of its purpose. It is called by FEFaceValues::reinit().
   *
   * @param[in] cell The cell of the triangulation for which this function
   *   is to compute a mapping from the reference cell to.
   * @param[in] face_no The number of the face of the given cell for which
   *   information is requested.
   * @param[in] quadrature A reference to the quadrature formula in use
   *   for the current evaluation. This quadrature object is the same
   *   as the one used when creating the @p internal_data object. The
   *   object is used both to map the location of quadrature points,
   *   as well as to compute the JxW values for each quadrature
   *   point (which involves the quadrature weights).
   * @param[in] internal_data A reference to an object previously
   *   created by get_data() and that may be used to store information
   *   the mapping can compute once on the reference cell. See the
   *   documentation of the Mapping::InternalDataBase class for an
   *   extensive description of the purpose of these objects.
   * @param[out] output_data A reference to an object whose member
   *   variables should be computed. Not all of the members of this
   *   argument need to be filled; which ones need to be filled is
   *   determined by the update flags stored inside the
   *   @p internal_data object.
   */
  virtual void
  fill_fe_face_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                       const unsigned int                                         face_no,
                       const Quadrature<dim-1>                                   &quadrature,
                       const InternalDataBase                                    &internal_data,
                       internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const = 0;

  /**
   * This function is the equivalent to Mapping::fill_fe_values(),
   * but for subfaces (i.e., children of faces) of cells.
   * See there for an extensive discussion
   * of its purpose. It is called by FESubfaceValues::reinit().
   *
   * @param[in] cell The cell of the triangulation for which this function
   *   is to compute a mapping from the reference cell to.
   * @param[in] face_no The number of the face of the given cell for which
   *   information is requested.
   * @param[in] subface_no The number of the child of a face of the
   *   given cell for which information is requested.
   * @param[in] quadrature A reference to the quadrature formula in use
   *   for the current evaluation. This quadrature object is the same
   *   as the one used when creating the @p internal_data object. The
   *   object is used both to map the location of quadrature points,
   *   as well as to compute the JxW values for each quadrature
   *   point (which involves the quadrature weights).
   * @param[in] internal_data A reference to an object previously
   *   created by get_data() and that may be used to store information
   *   the mapping can compute once on the reference cell. See the
   *   documentation of the Mapping::InternalDataBase class for an
   *   extensive description of the purpose of these objects.
   * @param[out] output_data A reference to an object whose member
   *   variables should be computed. Not all of the members of this
   *   argument need to be filled; which ones need to be filled is
   *   determined by the update flags stored inside the
   *   @p internal_data object.
   */
  virtual void
  fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                          const unsigned int                                         face_no,
                          const unsigned int                                         subface_no,
                          const Quadrature<dim-1>                                   &quadrature,
                          const InternalDataBase                                    &internal_data,
                          internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const = 0;

  /**
   * Give class @p FEValues access to the private <tt>get_...data</tt> and
   * <tt>fill_fe_...values</tt> functions.
   */
  friend class FEValuesBase<dim,spacedim>;
  friend class FEValues<dim,spacedim>;
  friend class FEFaceValues<dim,spacedim>;
  friend class FESubfaceValues<dim,spacedim>;
};


/* ------------------------- inline functions ------------------------- */

#ifndef DOXYGEN

template <int dim, int spacedim>
inline
UpdateFlags
Mapping<dim,spacedim>::InternalDataBase::current_update_flags () const
{
  if (first_cell)
    {
      Assert (update_flags==(update_once|update_each),
              ExcInternalError());
      return update_flags;
    }
  else
    return update_each;
}



template <int dim, int spacedim>
inline
bool
Mapping<dim,spacedim>::InternalDataBase::is_first_cell () const
{
  return first_cell;
}



template <int dim, int spacedim>
inline
void
Mapping<dim,spacedim>::InternalDataBase::clear_first_cell ()
{
  first_cell = false;
}



template <int dim, int spacedim>
inline
const Point<spacedim> &
Mapping<dim,spacedim>::support_point_value(
  const unsigned int index,
  const typename Mapping<dim,spacedim>::InternalDataBase &internal) const
{
  AssertIndexRange(index, internal.support_point_values.size());
  return internal.support_point_values[index];
}


template <int dim, int spacedim>
inline
const Tensor<2,spacedim> &
Mapping<dim,spacedim>::support_point_gradient(
  const unsigned int index,
  const typename Mapping<dim,spacedim>::InternalDataBase &internal) const
{
  AssertIndexRange(index, internal.support_point_gradients.size());
  return internal.support_point_gradients[index];
}


template <int dim, int spacedim>
inline
const Tensor<2,spacedim> &
Mapping<dim,spacedim>::support_point_inverse_gradient(
  const unsigned int index,
  const typename Mapping<dim,spacedim>::InternalDataBase &internal) const
{
  AssertIndexRange(index, internal.support_point_inverse_gradients.size());
  return internal.support_point_inverse_gradients[index];
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

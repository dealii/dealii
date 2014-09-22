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

#ifndef __deal2__mapping_h
#define __deal2__mapping_h


#include <deal.II/base/config.h>
#include <deal.II/base/derivative_form.h>
#include <deal.II/base/vector_slice.h>
#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe_update_flags.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

template <int dim> class Quadrature;
template <int dim, int spacedim> class FEValuesData;
template <int dim, int spacedim> class FEValuesBase;
template <int dim, int spacedim> class FEValues;
template <int dim, int spacedim> class FEFaceValues;
template <int dim, int spacedim> class FESubfaceValues;

/**
 * The transformation type used
 * for the Mapping::transform() functions.
 *
 * Special finite elements may
 * need special Mapping from the
 * reference cell to the actual
 * mesh cell. In order to be most
 * flexible, this enum provides
 * an extensible interface for
 * arbitrary
 * transformations. Nevertheless,
 * these must be implemented in
 * the transform() functions of
 * inheriting classes in order to
 * work.
 *
 * @ingroup mapping
 */
enum MappingType
{
/// No mapping
  mapping_none = 0x0000,
/// Covariant mapping (see Mapping::transform() for details)
  mapping_covariant = 0x0001,
/// Contravariant mapping (see Mapping::transform() for details)
  mapping_contravariant = 0x0002,
/// Mapping of the gradient of a covariant vector field (see Mapping::transform() for details)
  mapping_covariant_gradient = 0x0003,
/// Mapping of the gradient of a contravariant vector field (see Mapping::transform() for details)
  mapping_contravariant_gradient = 0x0004,
  /**
   * The Piola transform usually used for Hdiv elements.
   * Piola transform is the
   * standard transformation of
   * vector valued elements in
   * H<sup>div</sup>. It amounts
   * to a contravariant
   * transformation scaled by the
   * inverse of the volume
   * element.
   */
  mapping_piola = 0x0100,
  /**
     transformation for the gradient of a vector field corresponding to a
     mapping_piola transformation (see Mapping::transform() for details).
  */

  mapping_piola_gradient = 0x0101,
/// The mapping used for Nedelec elements
  /**
   * curl-conforming elements are
   * mapped as covariant
   * vectors. Nevertheless, we
   * introduce a separate mapping
   * type, such that we can use
   * the same flag for the vector
   * and its gradient (see Mapping::transform() for details).
   */
  mapping_nedelec = 0x0200,
/// The mapping used for Raviart-Thomas elements
  mapping_raviart_thomas = 0x0300,
/// The mapping used for BDM elements
  mapping_bdm = mapping_raviart_thomas
};


/**
 * Abstract base class for mapping classes.
 *
 * The interface for filling the tables of FEValues is provided.
 * Everything else has to happen in derived classes.
 *
 * The following paragraph applies to the implementation of
 * FEValues. Usage of the class is as follows: first, call the
 * functions @p update_once and @p update_each with the update
 * flags you need. This includes the flags needed by the
 * FiniteElement. Then call <tt>get_*_data</tt> and with the or'd
 * results.  This will initialize and return some internal data
 * structures.  On the first cell, call <tt>fill_fe_*_values</tt> with the
 * result of @p update_once. Finally, on each cell, use
 * <tt>fill_fe_*_values</tt> with the result of @p update_each to compute
 * values for a special cell.
 *
 * <h3>Mathematics of the mapping</h3>
 *
 * The mapping is a transformation $\mathbf x = \Phi(\mathbf{\hat x})$
 * which maps the reference cell [0,1]<sup>dim</sup> to the actual
 * grid cell in R<sup>spacedim</sup>.
 * In order to describe the application of the mapping to
 * different objects, we introduce the notation for the Jacobian
 * $J(\mathbf{\hat x}) = \nabla\Phi(\mathbf{\hat x})$. For instance,
 * if dim=spacedim=2, we have
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
 * Since finite element shape functions are usually defined on the
 * reference cell, nothing needs to be done for them. For a function
 * defined on the computational domain, the quadrature points need to
 * be mapped, which is done in fill_fe_values() if
 * @p update_quadrature_points is set in the update flags. The mapped
 * quadrature points are then accessed through FEValuesBase::quadrature_point().
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
 * x})\right|$ are accessed through FEValuesBase::JxW() and
 * computed in fill_fe_values(), if @p update_JxW_values is set in the
 * update flags.
 *
 * @todo Add a function <tt>transform_quadrature_weights</tt> for
 * this.
 *
 * @todo Add documentation on the codimension-one case
 *
 * <h4>Mapping of vector fields, differential forms and gradients of vector fields</h4>
 *
 * The transfomation of vector fields, differential forms (gradients/jacobians) and
 * gradients of vector fields between the reference cell and the actual grid cell
 * follows the general form
 *
 * @f[
 * \mathbf v(\mathbf x) = \mathbf A(\mathbf{\hat x})
 * \mathbf{\hat v}(\mathbf{\hat x}),
 * \qquad
 * \mathbf T(\mathbf x) = \mathbf A(\mathbf{\hat x})
 * \mathbf{\hat T}(\mathbf{\hat x}) \mathbf B(\mathbf{\hat x}),
 * @f]
 *
 * where <b>v</b> is a vector field or a differential form and
 * and <b>T</b> a tensor field of gradients.
 * The differential forms <b>A</b> and <b>B</b> are
 * determined by the MappingType enumerator.
 * These transformations are performed through the functions
 * transform(). See the documentation there for possible
 * choices.
 *
 * <h3>Technical notes</h3>
 *
 * A hint to implementators: no function except the two functions
 * @p update_once and @p update_each may add any flags.
 *
 * For more information about the <tt>spacedim</tt> template parameter
 * check the documentation of FiniteElement or the one of
 * Triangulation.
 *
 * <h3>References</h3>
 *
 * A general publication on differential geometry and finite elements
 * is the survey
 * <ul>
 * <li>Douglas N. Arnold, Richard S. Falk, and
 * Ragnar Winther. <i>Finite element exterior calculus: from
 * Hodge theory to numerical stability.</i>
 * Bull. Amer. Math. Soc. (N.S.), 47:281-354, 2010. <a
 * href="http://dx.doi.org/10.1090/S0273-0979-10-01278-4">DOI:
 * 10.1090/S0273-0979-10-01278-4</a>.
 * </ul>
 *
 * The description of the Piola transform has been taken from the <a
 * href="http://www.math.uh.edu/~rohop/spring_05/downloads/">lecture
 * notes</a> by Ronald H. W. Hoppe, University of Houston, Chapter 7.
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
   * Transforms the point @p p on
   * the unit cell to the point
   * @p p_real on the real cell
   * @p cell and returns @p p_real.
   */
  virtual Point<spacedim>
  transform_unit_to_real_cell (
    const typename Triangulation<dim,spacedim>::cell_iterator &cell,
    const Point<dim>                                 &p) const = 0;

  /**
   * Transforms the point @p p on
   * the real @p cell to the corresponding
   * point on the unit cell, and
   * return its coordinates.
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
    const Point<spacedim>                            &p) const = 0;

  /**
   * Base class for internal data
   * of finite element and mapping
   * objects. The internal
   * mechanism is that upon
   * construction of a @p FEValues
   * objects, it asks the mapping
   * and finite element classes
   * that are to be used to
   * allocate memory for their own
   * purpose in which they may
   * store data that only needs to
   * be computed once. For example,
   * most finite elements will
   * store the values of the shape
   * functions at the quadrature
   * points in this object, since
   * they do not change from cell
   * to cell and only need to be
   * computed once. Since different
   * @p FEValues objects using
   * different quadrature rules
   * might access the same finite
   * element object at the same
   * time, it is necessary to
   * create one such object per
   * @p FEValues object. Ownership
   * of this object is then
   * transferred to the
   * @p FEValues object, but a
   * pointer to this object is
   * passed to the finite element
   * object every time it shall
   * compute some data so that it
   * has access to the precomputed
   * values stored there.
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
     * Constructor. Sets
     * @p UpdateFlags to
     * @p update_default and
     * @p first_cell to @p true.
     */
    InternalDataBase ();

    /**
     * Virtual destructor for
     * derived classes
     */
    virtual ~InternalDataBase ();

    /**
     * Values updated by the constructor or
     * by reinit.
     */
    UpdateFlags          update_flags;

    /**
     * Values computed by
     * constructor.
     */
    UpdateFlags          update_once;

    /**
     * Values updated on each
     * cell by reinit.
     */
    UpdateFlags          update_each;

    /**
     * If <tt>first_cell==true</tt>
     * this function returns
     * @p update_flags,
     * i.e. <tt>update_once|update_each</tt>.
     * If <tt>first_cell==false</tt>
     * it returns
     * @p update_each.
     */
    UpdateFlags  current_update_flags() const;

    /**
     * Return whether we are
     * presently initializing
     * data for the first
     * cell. The value of the
     * field this function is
     * returning is set to
     * @p true in the
     * constructor, and cleared
     * by the @p FEValues class
     * after the first cell has
     * been initialized.
     *
     * This function is used to
     * determine whether we need
     * to use the @p update_once
     * flags for computing data,
     * or whether we can use the
     * @p update_each flags.
     */
    bool is_first_cell () const;

    /**
     * Set the @p first_cell
     * flag to @p false. Used by
     * the @p FEValues class to
     * indicate that we have
     * already done the work on
     * the first cell.
     */
    virtual void clear_first_cell ();

    /**
     * Return an estimate (in
     * bytes) or the memory
     * consumption of this
     * object.
     */
    virtual std::size_t memory_consumption () const;

    /**
     * The determinant of the
     * Jacobian in each
     * quadrature point. Filled
     * if #update_volume_elements.
     */
    std::vector<double> volume_elements;

    /**
     * The positions of the
     * mapped (generalized)
     * support points.
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
     * The value returned by
     * @p is_first_cell.
     */
    bool first_cell;
  };

  /**
   * Transform a field of vectors or 1-differential forms according to the selected
   * MappingType.
   *
   * @note Normally, this function is called by a finite element,
   * filling FEValues objects. For this finite element, there should be
   * an alias MappingType like @p mapping_bdm, @p mapping_nedelec, etc. This
   * alias should be preferred to using the types below.
   *
   * The mapping types currently implemented by derived classes are:
   * <ul>
   * <li> @p mapping_contravariant: maps a vector field on the reference cell
   * is to the physical cell through the Jacobian:
   * @f[
   * \mathbf u(\mathbf x) = J(\mathbf{\hat x})\mathbf{\hat u}(\mathbf{\hat x}).
   * @f]
   * In physics, this is usually referred to as the contravariant
   * transformation. Mathematically, it is the push forward of a
   * vector field.
   *
   * <li> @p mapping_covariant: maps a field of one-forms on the reference cell
   * to a field of one-forms on the physical cell.
   * (theoretically this would refer to a DerivativeForm<1, dim, 1> but it
   * canonically identified with a Tensor<1,dim>).
   * Mathematically, it is the pull back of the differential form
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
   * <li> @p mapping_piola: A field of <i>n-1</i>-forms on the reference cell is also
   * represented by a vector field, but again transforms differently,
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
     Transform a field of differential forms from the reference cell to the physical cell.

     It is useful to think of $\mathbf{T} = D \mathbf u$ and
     $\mathbf{\hat T} = \hat D \mathbf{\hat u}$, with $\mathbf u$ a vector field.

     The mapping types currently implemented by derived classes are:
     <ul>
     <li> @p mapping_covariant: maps a field of forms on the reference cell
     to a field of forms on the physical cell.
     Mathematically, it is the pull back of the differential form
     @f[
     \mathbf T(\mathbf x) = \mathbf{\hat T}(\mathbf{\hat x})
                            J*(J^{T} J)^{-1}(\mathbf{\hat x}).
     @f]
     n the case when dim=spacedim the previous formula reduces to
     @f[
     \mathbf T(\mathbf x) = \mathbf{\hat u}(\mathbf{\hat x})
                            J^{-1}(\mathbf{\hat x}).
     @f]
    Jacobians of spacedim-vector valued differentiable functions are transformed this way.
     </ul>

     @note It would have been more reasonable to make this transform a template function
     with the rank in <code>DerivativeForm@<1, dim, rank@></code>. Unfortunately C++ does not
     allow templatized virtual functions. This is why we identify
     <code>DerivativeForm@<1, dim, 1@></code> with a <code>Tensor@<1,dim@></code>
     when using  mapping_covariant() in the function transform above this one.
  */

  virtual
  void
  transform (const VectorSlice<const std::vector< DerivativeForm<1, dim, spacedim> > > input,
             VectorSlice<std::vector<Tensor<2,spacedim> > >             output,
             const InternalDataBase &internal,
             const MappingType type) const = 0;



  /**
     Transform a tensor field from the reference cell to the physical cell.
     This tensors are most of times the jacobians in the reference cell of
     vector fields that have been pulled back from the physical cell.

     The mapping types currently implemented by derived classes are:
     <ul>

     <li> @p mapping_contravariant_gradient, it
     assumes $\mathbf u(\mathbf x) = J \mathbf{\hat u}$ so that
     @f[
     \mathbf T(\mathbf x) =
     J(\mathbf{\hat x}) \mathbf{\hat T}(\mathbf{\hat x})
     J^{-1}(\mathbf{\hat x}).
     @f]

     <li> @p mapping_covariant_gradient, it
     assumes $\mathbf u(\mathbf x) = J^{-T} \mathbf{\hat u}$ so that
     @f[
     \mathbf T(\mathbf x) =
     J^{-T}(\mathbf{\hat x}) \mathbf{\hat T}(\mathbf{\hat x})
     J^{-1}(\mathbf{\hat x}).
     @f]

     <li> @p mapping_piola_gradient, it
     assumes $\mathbf u(\mathbf x) = \frac{1}{\text{det}J(\mathbf x)}
     J(\mathbf x) \mathbf{\hat u}(\mathbf x)$
     so that
     @f[
     \mathbf T(\mathbf x) =
     \frac{1}{\text{det}J(\mathbf x)}
     J(\mathbf{\hat x}) \mathbf{\hat T}(\mathbf{\hat x})
     J^{-1}(\mathbf{\hat x}).
     @f]
     </ul>

     @todo The formulas for mapping_covariant_gradient(),
     mapping_contravariant_gradient() and mapping_piola_gradient()
     are only true as stated for linear mappings.
     If, for example, the mapping is bilinear then there is a missing
     term associated with the derivative of J.
  */
  virtual
  void
  transform (const VectorSlice<const std::vector<Tensor<2, dim> > >     input,
             VectorSlice<std::vector<Tensor<2,spacedim> > >             output,
             const InternalDataBase &internal,
             const MappingType type) const = 0;

  /**
   * @deprecated Use transform() instead.
   */
  void
  transform_covariant (const VectorSlice<const std::vector<Tensor<1,dim> > > input,
                       const unsigned int                                    offset,
                       VectorSlice<std::vector<Tensor<1,spacedim> > >        output,
                       const InternalDataBase &internal) const DEAL_II_DEPRECATED;

  /**
   * @deprecated Use transform() instead.
   */
  void
  transform_covariant (const VectorSlice<const std::vector<DerivativeForm<1, dim ,spacedim> > > input,
                       const unsigned int                 offset,
                       VectorSlice<std::vector<Tensor<2,spacedim> > >      output,
                       const InternalDataBase &internal) const DEAL_II_DEPRECATED;

  /**
   * @deprecated Use transform() instead.
   */
  void
  transform_contravariant (const VectorSlice<const std::vector<Tensor<1,dim> > > input,
                           const unsigned int                 offset,
                           VectorSlice<std::vector<Tensor<1,spacedim> > >      output,
                           const typename Mapping<dim,spacedim>::InternalDataBase &internal) const DEAL_II_DEPRECATED;

  /**
   * @deprecated Use transform() instead.
   */

  void
  transform_contravariant (const VectorSlice<const std::vector<DerivativeForm<1, dim,spacedim> > > input,
                           const unsigned int                 offset,
                           const VectorSlice<std::vector<Tensor<2,spacedim> > > output,
                           const typename Mapping<dim,spacedim>::InternalDataBase &internal) const DEAL_II_DEPRECATED;

  /**
   * The transformed (generalized)
   * support point.
   */
  const Point<spacedim> &support_point_value(
    const unsigned int index,
    const typename Mapping<dim,spacedim>::InternalDataBase &internal) const;

  /**
   * The Jacobian
   * matrix of the transformation
   * in the (generalized) support
   * point.
   */
  const Tensor<2,spacedim> &support_point_gradient(
    const unsigned int index,
    const typename Mapping<dim,spacedim>::InternalDataBase &internal) const;

  /**
   * The inverse Jacobian
   * matrix of the transformation
   * in the (generalized) support
   * point.
   */
  const Tensor<2,spacedim> &support_point_inverse_gradient(
    const unsigned int index,
    const typename Mapping<dim,spacedim>::InternalDataBase &internal) const;

  /**
   * Return a pointer to a copy of the
   * present object. The caller of this
   * copy then assumes ownership of it.
   *
   * Since one can't create
   * objects of class Mapping, this
   * function of course has to be
   * implemented by derived classes.
   *
   * This function is mainly used by the
   * hp::MappingCollection class.
   */
  virtual
  Mapping<dim,spacedim> *clone () const = 0;

  /**
   * Returns whether the mapping preserves
   * vertex locations, i.e. whether the
   * mapped location of the reference cell
   * vertices (given by
   * GeometryInfo::unit_cell_vertex())
   * equals the result of
   * <code>cell-@>vertex()</code>.
   *
   * For example, implementations in
   * derived classes return @p true for
   * MappingQ, MappingQ1, MappingCartesian,
   * but @p false for MappingQEulerian,
   * MappingQ1Eulerian.
   */
  virtual
  bool preserves_vertex_locations () const = 0;

  /**
   * Exception
   */
  DeclException0 (ExcInvalidData);


  /**
   * Computing the mapping between a
   * real space point and a point
   * in reference space failed, typically because the given point
   * lies outside the cell where the inverse mapping is not
   * unique.
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
   * Indicate fields to be updated
   * in the constructor of
   * FEValues. Especially,
   * fields not asked for by
   * FEValues, but computed
   * for efficiency reasons will be
   * notified here.
   *
   * See @ref UpdateFlagsEssay.
   */
  virtual UpdateFlags update_once (const UpdateFlags) const = 0;

  /**
   * The same as update_once(),
   * but for the flags to be updated for
   * each grid cell.
   *
   * See @ref UpdateFlagsEssay.
   */
  virtual UpdateFlags update_each (const UpdateFlags) const = 0;

  /**
   * Prepare internal data
   * structures and fill in values
   * independent of the cell.
   */
  virtual InternalDataBase *
  get_data (const UpdateFlags,
            const Quadrature<dim> &quadrature) const = 0;

  /**
   * Prepare internal data
   * structure for transformation
   * of faces and fill in values
   * independent of the cell.
   */
  virtual InternalDataBase *
  get_face_data (const UpdateFlags flags,
                 const Quadrature<dim-1>& quadrature) const = 0;

  /**
   * Prepare internal data
   * structure for transformation
   * of children of faces and fill
   * in values independent of the
   * cell.
   */
  virtual InternalDataBase *
  get_subface_data (const UpdateFlags flags,
                    const Quadrature<dim-1>& quadrature) const = 0;


  /**
   * Fill the transformation fields
   * of @p FEValues.  Given a grid
   * cell and the quadrature points
   * on the unit cell, it computes
   * all values specified by
   * @p flags. The arrays to be
   * filled have to have the
   * correct size.
   *
   * Values are split into two
   * groups: first,
   * @p quadrature_points and
   * @p JxW_values are
   * filled with the quadrature
   * rule transformed to the
   * cell in physical space.
   *
   * The second group contains the
   * matrices needed to transform
   * vector-valued functions,
   * namely
   * @p jacobians,
   * the derivatives
   * @p jacobian_grads,
   * and the inverse operations in
   * @p inverse_jacobians.
   */
  /*     virtual void */
  /*     fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell, */
  /*                  const Quadrature<dim>                         &quadrature, */
  /*                  InternalDataBase                              &internal, */
  /*                  std::vector<Point<spacedim> >                 &quadrature_points, */
  /*                  std::vector<double>                           &JxW_values) const = 0; */

  /** The function above adjusted
   * with the variable
   * cell_normal_vectors for the
   * case of codimension 1
   */
  virtual void
  fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                  const Quadrature<dim>                                     &quadrature,
                  InternalDataBase                                          &internal,
                  std::vector<Point<spacedim> >                             &quadrature_points,
                  std::vector<double>                                       &JxW_values,
                  std::vector<DerivativeForm<1,dim,spacedim>  >       &jacobians,
                  std::vector<DerivativeForm<2,dim,spacedim>  >       &jacobian_grads,
                  std::vector<DerivativeForm<1,spacedim,dim>  >       &inverse_jacobians,
                  std::vector<Point<spacedim> >                             &cell_normal_vectors,
                  CellSimilarity::Similarity                           &cell_similarity
                 ) const=0;



  /**
   * Performs the same as @p
   * fill_fe_values on a face.
   * Additionally, @p boundary_form
   * (see @ref GlossBoundaryForm)
   * and @p normal_vectors can be
   * computed on surfaces. Since
   * the boundary form already
   * contains the determinant of
   * the Jacobian of the
   * transformation, it is
   * sometimes more economic to use
   * the boundary form instead of
   * the product of the unit normal
   * and the transformed quadrature
   * weight.
   */
  virtual void
  fill_fe_face_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                       const unsigned int                                        face_no,
                       const Quadrature<dim-1>                                   &quadrature,
                       InternalDataBase                                          &internal,
                       std::vector<Point<spacedim> >                             &quadrature_points,
                       std::vector<double>                                       &JxW_values,
                       std::vector<Tensor<1,spacedim> >                          &boundary_form,
                       std::vector<Point<spacedim> >                             &normal_vectors) const = 0;

  /**
   * See above.
   */
  virtual void
  fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                          const unsigned int                        face_no,
                          const unsigned int                        sub_no,
                          const Quadrature<dim-1>                  &quadrature,
                          InternalDataBase                         &internal,
                          std::vector<Point<spacedim> >        &quadrature_points,
                          std::vector<double>                      &JxW_values,
                          std::vector<Tensor<1,spacedim> >     &boundary_form,
                          std::vector<Point<spacedim> >        &normal_vectors) const = 0;

  /**
   * Give class @p FEValues access
   * to the private <tt>get_...data</tt>
   * and <tt>fill_fe_...values</tt>
   * functions.
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

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

#ifndef dealii__mapping_h
#define dealii__mapping_h


#include <deal.II/base/config.h>
#include <deal.II/base/derivative_form.h>
#include <deal.II/base/std_cxx11/array.h>
#include <deal.II/base/array_view.h>
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
   * Mapping of the gradient of a covariant vector field (see
   * Mapping::transform() for details).
   */
  mapping_covariant_gradient = 0x0003,

  /**
   * Mapping of the gradient of a contravariant vector field (see
   * Mapping::transform() for details).
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
  mapping_bdm = mapping_raviart_thomas,

  /**
   * The mappings for 2-forms and third order tensors.
   *
   * These are mappings typpically applied to hessians transformed to the
   * reference cell.
   *
   * Mapping of the hessian of a covariant vector field (see
   * Mapping::transform() for details).
   */
  mapping_covariant_hessian,

  /**
   * Mapping of the hessian of a contravariant vector field (see
   * Mapping::transform() for details).
   */
  mapping_contravariant_hessian,

  /**
   * Mapping of the hessian of a piola vector field (see Mapping::transform()
   * for details).
   */
  mapping_piola_hessian
};


/**
 * @short Abstract base class for mapping classes.
 *
 * This class declares the interface for the functionality to describe
 * mappings from the reference (unit) cell to a cell in real space, as well as
 * for filling the information necessary to use the FEValues, FEFaceValues,
 * and FESubfaceValues classes. Concrete implementations of these interfaces
 * are provided in derived classes.
 *
 * <h3>Mathematics of the mapping</h3>
 *
 * The mapping is a transformation $\mathbf x = \mathbf F_K(\hat{\mathbf  x})$
 * which maps points $\hat{\mathbf x}$ in the reference cell
 * $[0,1]^\text{dim}$ to points $\mathbf x$ in the actual grid cell
 * $K\subset{\mathbb R}^\text{spacedim}$. Many of the applications of such
 * mappings require the Jacobian of this mapping, $J(\hat{\mathbf x}) =
 * \hat\nabla {\mathbf F}_K(\hat{\mathbf  x})$. For instance, if
 * dim=spacedim=2, we have
 * @f[
 * J(\hat{\mathbf  x}) = \left(\begin{matrix}
 * \frac{\partial x}{\partial \hat x} & \frac{\partial x}{\partial \hat y}
 * \\
 * \frac{\partial y}{\partial \hat x} & \frac{\partial y}{\partial \hat y}
 * \end{matrix}\right)
 * @f]
 *
 * <h4>%Mapping of scalar functions</h4>
 *
 * The shape functions of scalar finite elements are typically defined on a
 * reference cell and are then simply mapped according to the rule
 * @f[
 * \varphi(\mathbf x) = \varphi\bigl(\mathbf F_K(\hat{\mathbf  x})\bigr)
 * = \hat \varphi(\hat{\mathbf  x}).
 * @f]
 *
 *
 * <h4>%Mapping of integrals</h4>
 *
 * Using simply a change of variables, integrals of scalar functions over a
 * cell $K$ can be expressed as an integral over the reference cell $\hat K$.
 * Specifically, The volume form $d\hat x$ is transformed so that
 * @f[
 *  \int_K u(\mathbf x)\,dx = \int_{\hat K} \hat
 * u(\hat{\mathbf  x}) \left|\text{det}J(\hat{\mathbf  x})\right|
 * \,d\hat x.
 * @f]
 *
 * In expressions where such integrals are approximated by quadrature, this
 * then leads to terms of the form
 * @f[
 *  \int_K u(\mathbf x)\,dx
 *  \approx
 *  \sum_{q}
 *  \hat u(\hat{\mathbf  x}_q)
 *  \underbrace{\left|\text{det}J(\hat{\mathbf  x}_q)\right| w_q}_{=: \text{JxW}_q}.
 * @f]
 * Here, the weights $\text{JxW}_q$ of each quadrature point (where <i>JxW</i>
 * mnemonically stands for <i>Jacobian times Quadrature Weights</i>) take the
 * role of the $dx$ in the original integral. Consequently, they appear in all
 * code that computes integrals approximated by quadrature, and are accessed
 * by FEValues::JxW().
 *
 * @todo Document what happens in the codimension-1 case.
 *
 *
 * <h4>%Mapping of vector fields, differential forms and gradients of vector
 * fields</h4>
 *
 * The transformation of vector fields or differential forms (gradients of
 * scalar functions) $\mathbf v$, and gradients of vector fields $\mathbf T$
 * follows the general form
 *
 * @f[
 * \mathbf v(\mathbf x) = \mathbf A(\hat{\mathbf  x})
 * \hat{\mathbf  v}(\hat{\mathbf  x}),
 * \qquad
 * \mathbf T(\mathbf x) = \mathbf A(\hat{\mathbf  x})
 * \hat{\mathbf  T}(\hat{\mathbf  x}) \mathbf B(\hat{\mathbf  x}).
 * @f]
 * The differential forms <b>A</b> and <b>B</b> are determined by the kind of
 * object being transformed. These transformations are performed through the
 * transform() functions, and the type of object being transformed is
 * specified by their MappingType argument. See the documentation there for
 * possible choices.
 *
 * <h4>Derivatives of the mapping</h4>
 *
 * Some applications require the derivatives of the mapping, of which the
 * first order derivative is the mapping Jacobian, $J_{iJ}(\hat{\mathbf
 * x})=\frac{\partial x_i}{\partial \hat x_J}$, described above. Higher order
 * derivatives of the mapping are similarly defined, for example the Jacobian
 * derivative, $\hat H_{iJK}(\hat{\mathbf  x}) = \frac{\partial^2
 * x_i}{\partial \hat x_J \partial \hat x_K}$, and the Jacobian second
 * derivative, $\hat K_{iJKL}(\hat{\mathbf  x}) = \frac{\partial^3
 * x_i}{\partial \hat x_J \partial \hat x_K \partial \hat x_L}$. It is also
 * useful to define the "pushed-forward" versions of the higher order
 * derivatives: the Jacobian pushed-forward derivative, $H_{ijk}(\hat{\mathbf
 * x}) = \frac{\partial^2 x_i}{\partial \hat x_J \partial \hat
 * x_K}(J_{jJ})^{-1}(J_{kK})^{-1}$, and the Jacobian pushed-forward second
 * derivative, $K_{ijkl}(\hat{\mathbf  x}) = \frac{\partial^3 x_i}{\partial
 * \hat x_J \partial \hat x_K \partial \hat
 * x_L}(J_{jJ})^{-1}(J_{kK})^{-1}(J_{lL})^{-1}$. These pushed-forward versions
 * can be used to compute the higher order derivatives of functions defined on
 * the reference cell with respect to the real cell coordinates. For instance,
 * the Jacobian derivative with respect to the real cell coordinates is given
 * by:
 *
 * @f[
 * \frac{\partial}{\partial x_j}\left[J_{iJ}(\hat{\mathbf  x})\right] =
 * H_{ikn}(\hat{\mathbf  x})J_{nJ}(\hat{\mathbf  x}),
 * @f]
 * and the derivative of the Jacobian inverse with respect to the real cell
 * coordinates is similarly given by:
 * @f[
 * \frac{\partial}{\partial x_j}\left[\left(J_{iJ}(\hat{\mathbf  x})\right)^{-1}\right]
 * = -H_{nik}(\hat{\mathbf  x})\left(J_{nJ}(\hat{\mathbf  x})\right)^{-1}.
 * @f]
 *
 * In a similar fashion, higher order derivatives, with respect to the real
 * cell coordinates, of functions defined on the reference cell can be defined
 * using the Jacobian pushed-forward higher-order derivatives. For example,
 * the derivative, with respect to the real cell coordinates, of the Jacobian
 * pushed-forward derivative is given by:
 *
 * @f[
 * \frac{\partial}{\partial x_l}\left[H_{ijk}(\hat{\mathbf  x})\right] = K_{ijkl}(\hat{\mathbf  x})
 * -H_{mjl}(\hat{\mathbf  x})H_{imk}(\hat{\mathbf  x})-H_{mkl}(\hat{\mathbf  x})H_{imj}(\hat{\mathbf  x}).
 * @f]
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
   * Return a pointer to a copy of the present object. The caller of this copy
   * then assumes ownership of it.
   *
   * The function is declared abstract virtual in this base class, and derived
   * classes will have to implement it.
   *
   * This function is mainly used by the hp::MappingCollection class.
   */
  virtual
  Mapping<dim,spacedim> *clone () const = 0;

  /**
   * Return the mapped vertices of a cell.
   *
   * Most of the time, these values will simply be the coordinates of the
   * vertices of a cell as returned by <code>cell-@>vertex(v)</code> for
   * vertex <code>v</code>, i.e., information stored by the triangulation.
   * However, there are also mappings that add displacements or choose
   * completely different locations, e.g., MappingQEulerian,
   * MappingQ1Eulerian, or MappingFEField.
   *
   * The default implementation of this function simply returns the
   * information stored by the triangulation, i.e.,
   * <code>cell-@>vertex(v)</code>.
   */
  virtual
  std_cxx11::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
  get_vertices (const typename Triangulation<dim,spacedim>::cell_iterator &cell) const;

  /**
   * Returns whether the mapping preserves vertex locations. In other words,
   * this function returns whether the mapped location of the reference cell
   * vertices (given by GeometryInfo::unit_cell_vertex()) equals the result of
   * <code>cell-@>vertex()</code> (i.e., information stored by the
   * triangulation).
   *
   * For example, implementations in derived classes return @p true for
   * MappingQ, MappingQGeneric, MappingCartesian, but @p false for
   * MappingQEulerian, MappingQ1Eulerian, and MappingFEField.
   */
  virtual
  bool preserves_vertex_locations () const = 0;

  /**
   * @name Mapping points between reference and real cells
   * @{
   */

  /**
   * Maps the point @p p on the unit cell to the corresponding point on the
   * real cell @p cell.
   *
   * @param cell Iterator to the cell that will be used to define the mapping.
   * @param p Location of a point on the reference cell.
   * @return The location of the reference point mapped to real space using
   * the mapping defined by the class derived from the current one that
   * implements the mapping, and the coordinates of the cell identified by the
   * first argument.
   */
  virtual
  Point<spacedim>
  transform_unit_to_real_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                               const Point<dim>                                          &p) const = 0;

  /**
   * Maps the point @p p on the real @p cell to the corresponding point on the
   * unit cell, and return its coordinates. This function provides the inverse
   * of the mapping provided by transform_unit_to_real_cell().
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
   * checking whether the returned reference coordinates lie inside or outside
   * the reference cell (e.g., using GeometryInfo::is_inside_unit_cell()) or
   * whether the exception mentioned above has been thrown.
   *
   * @param cell Iterator to the cell that will be used to define the mapping.
   * @param p Location of a point on the given cell.
   * @return The reference cell location of the point that when mapped to real
   * space equals the coordinates given by the second argument. This mapping
   * uses the mapping defined by the class derived from the current one that
   * implements the mapping, and the coordinates of the cell identified by the
   * first argument.
   */
  virtual
  Point<dim>
  transform_real_to_unit_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                               const Point<spacedim>                                     &p) const = 0;

  /**
   * Transforms the point @p p on the real @p cell to the corresponding point
   * on the unit cell, and then projects it to a dim-1  point on the face with
   * the given face number @p face_no. Ideally the point @p p is near the face
   * @p face_no, but any point in the cell can technically be projected.
   *
   * This function does not make physical sense when dim=1, so it throws an
   * exception in this case.
   */
  Point<dim-1>
  project_real_point_to_unit_point_on_face (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                                            const unsigned int &face_no,
                                            const Point<spacedim> &p) const;

  /**
   * @}
   */


  /**
   * @name Exceptions
   * @{
   */

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
  DeclExceptionMsg (ExcTransformationFailed,
                    "Computing the mapping between a real space point and a point in reference "
                    "space failed, typically because the given point lies outside the cell "
                    "where the inverse mapping is not unique.");

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

  /**
   * @}
   */

  /**
   * @name Interface with FEValues
   * @{
   */

public:
  /**
   * Base class for internal data of mapping objects. The internal mechanism
   * is that upon construction of a FEValues object, it asks the mapping and
   * finite element classes that are to be used to allocate memory for their
   * own purpose in which they may store data that only needs to be computed
   * once. For example, most finite elements will store the values of the
   * shape functions at the quadrature points in this object, since they do
   * not change from cell to cell and only need to be computed once. The same
   * may be true for Mapping classes that want to only evaluate the shape
   * functions used for mapping once at the quadrature points.
   *
   * Since different FEValues objects using different quadrature rules might
   * access the same mapping object at the same time, it is necessary to
   * create one such object per FEValues object. FEValues does this by calling
   * Mapping::get_data(), or in reality the implementation of the
   * corresponding function in derived classes. Ownership of the object
   * created by Mapping::get_data() is then transferred to the FEValues
   * object, but a reference to this object is passed to the mapping object
   * every time it is asked to compute information on a concrete cell. This
   * happens when FEValues::reinit() (or the corresponding classes in
   * FEFaceValues and FESubfaceValues) call Mapping::fill_fe_values() (and
   * similarly via Mapping::fill_fe_face_values() and
   * Mapping::fill_fe_subface_values()).
   *
   * The purpose of this class is for mapping objects to store information
   * that can be computed once at the beginning, on the reference cell, and to
   * access it later when computing information on a concrete cell. As such,
   * the object handed to Mapping::fill_fe_values() is marked as
   * <code>const</code>, because the assumption is that at the time this
   * information is used, it will not need to modified again. However, classes
   * derived from Mapping can also use such objects for two other purposes:
   *
   * - To provide scratch space for computations that are done in
   * Mapping::fill_fe_values() and similar functions. Some of the derived
   * classes would like to use scratch arrays and it would be a waste of time
   * to allocate these arrays every time this function is called, just to de-
   * allocate it again at the end of the function. Rather, one could allocate
   * this memory once as a member variable of the current class, and simply
   * use it in Mapping::fill_fe_values().
   * - After calling Mapping::fill_fe_values(), FEValues::reinit()
   * calls FiniteElement::fill_fe_values() where the finite element computes
   * values, gradients, etc of the shape functions using both information
   * computed once at the beginning using a mechanism similar to the one
   * described here (see FiniteElement::InternalDataBase) as well as the data
   * already computed by Mapping::fill_fe_values(). As part of its work, some
   * implementations of FiniteElement::fill_fe_values() need to transform
   * shape function data, and they do so by calling Mapping::transform(). The
   * call to the latter function also receives a reference to the
   * Mapping::InternalDataBase object. Since Mapping::transform() may be
   * called many times on each cell, it is sometimes worth for derived classes
   * to compute some information only once in Mapping::fill_fe_values() and
   * reuse it in Mapping::transform(). This information can also be stored in
   * the classes that derived mapping classes derive from InternalDataBase.
   *
   * In both of these cases, the InternalDataBase object being passed around
   * is "morally const", i.e., no external observer can tell whether a scratch
   * array or some intermediate data for Mapping::transform() is being
   * modified by Mapping::fill_fe_values() or not. Consequently, the
   * InternalDataBase objects are always passed around as <code>const</code>
   * objects. Derived classes that would like to make use of the two
   * additional uses outlined above therefore need to mark the member
   * variables they want to use for these purposes as <code>mutable</code> to
   * allow for their modification despite the fact that the surrounding object
   * is marked as <code>const</code>.
   */
  class InternalDataBase
  {
  private:
    /**
     * Copy construction is forbidden.
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
     * A set of update flags specifying the kind of information that an
     * implementation of the Mapping interface needs to compute on each cell
     * or face, i.e., in Mapping::fill_fe_values() and friends.
     *
     * This set of flags is stored here by implementations of
     * Mapping::get_data(), Mapping::get_face_data(), or
     * Mapping::get_subface_data(), and is that subset of the update flags
     * passed to those functions that require re-computation on every cell.
     * (The subset of the flags corresponding to information that can be
     * computed once and for all already at the time of the call to
     * Mapping::get_data() -- or an implementation of that interface -- need
     * not be stored here because it has already been taken care of.)
     */
    UpdateFlags          update_each;

    /**
     * Return an estimate (in bytes) or the memory consumption of this object.
     */
    virtual std::size_t memory_consumption () const;
  };


protected:
  /**
   * Given a set of update flags, compute which other quantities <i>also</i>
   * need to be computed in order to satisfy the request by the given flags.
   * Then return the combination of the original set of flags and those just
   * computed.
   *
   * As an example, if @p update_flags contains update_JxW_values (i.e., the
   * product of the determinant of the Jacobian and the weights provided by
   * the quadrature formula), a mapping may require the computation of the
   * full Jacobian matrix in order to compute its determinant. They would then
   * return not just update_JxW_values, but also update_jacobians. (This is
   * not how it is actually done internally in the derived classes that
   * compute the JxW values -- they set update_contravariant_transformation
   * instead, from which the determinant can also be computed -- but this does
   * not take away from the instructiveness of the example.)
   *
   * An extensive discussion of the interaction between this function and
   * FEValues can be found in the
   * @ref FE_vs_Mapping_vs_FEValues
   * documentation module.
   *
   * @see UpdateFlags
   */
  virtual
  UpdateFlags
  requires_update_flags (const UpdateFlags update_flags) const = 0;

  /**
   * Create and return a pointer to an object into which mappings can store
   * data that only needs to be computed once but that can then be used
   * whenever the mapping is applied to a concrete cell (e.g., in the various
   * transform() functions, as well as in the fill_fe_values(),
   * fill_fe_face_values() and fill_fe_subface_values() that form the
   * interface of mappings with the FEValues class).
   *
   * Derived classes will return pointers to objects of a type derived from
   * Mapping::InternalDataBase (see there for more information) and may pre-
   * compute some information already (in accordance with what will be asked
   * of the mapping in the future, as specified by the update flags) and for
   * the given quadrature object. Subsequent calls to transform() or
   * fill_fe_values() and friends will then receive back the object created
   * here (with the same set of update flags and for the same quadrature
   * object). Derived classes can therefore pre-compute some information in
   * their get_data() function and store it in the internal data object.
   *
   * The mapping classes do not keep track of the objects created by this
   * function. Ownership will therefore rest with the caller.
   *
   * An extensive discussion of the interaction between this function and
   * FEValues can be found in the
   * @ref FE_vs_Mapping_vs_FEValues
   * documentation module.
   *
   * @param update_flags A set of flags that define what is expected of the
   * mapping class in future calls to transform() or the fill_fe_values()
   * group of functions. This set of flags may contain flags that mappings do
   * not know how to deal with (e.g., for information that is in fact computed
   * by the finite element classes, such as UpdateFlags::update_values).
   * Derived classes will need to store these flags, or at least that subset
   * of flags that will require the mapping to perform any actions in
   * fill_fe_values(), in InternalDataBase::update_each.
   * @param quadrature The quadrature object for which mapping information
   * will have to be computed. This includes the locations and weights of
   * quadrature points.
   * @return A pointer to a newly created object of type InternalDataBase (or
   * a derived class). Ownership of this object passes to the calling
   * function.
   *
   * @note C++ allows that virtual functions in derived classes may return
   * pointers to objects not of type InternalDataBase but in fact pointers to
   * objects of classes <i>derived</i> from InternalDataBase. (This feature is
   * called "covariant return types".) This is useful in some contexts where
   * the calling is within the derived class and will immediately make use of
   * the returned object, knowing its real (derived) type.
   */
  virtual
  InternalDataBase *
  get_data (const UpdateFlags      update_flags,
            const Quadrature<dim> &quadrature) const = 0;

  /**
   * Like get_data(), but in preparation for later calls to transform() or
   * fill_fe_face_values() that will need information about mappings from the
   * reference face to a face of a concrete cell.
   *
   * @param update_flags A set of flags that define what is expected of the
   * mapping class in future calls to transform() or the fill_fe_values()
   * group of functions. This set of flags may contain flags that mappings do
   * not know how to deal with (e.g., for information that is in fact computed
   * by the finite element classes, such as UpdateFlags::update_values).
   * Derived classes will need to store these flags, or at least that subset
   * of flags that will require the mapping to perform any actions in
   * fill_fe_values(), in InternalDataBase::update_each.
   * @param quadrature The quadrature object for which mapping information
   * will have to be computed. This includes the locations and weights of
   * quadrature points.
   * @return A pointer to a newly created object of type InternalDataBase (or
   * a derived class). Ownership of this object passes to the calling
   * function.
   *
   * @note C++ allows that virtual functions in derived classes may return
   * pointers to objects not of type InternalDataBase but in fact pointers to
   * objects of classes <i>derived</i> from InternalDataBase. (This feature is
   * called "covariant return types".) This is useful in some contexts where
   * the calling is within the derived class and will immediately make use of
   * the returned object, knowing its real (derived) type.
   */
  virtual
  InternalDataBase *
  get_face_data (const UpdateFlags        update_flags,
                 const Quadrature<dim-1> &quadrature) const = 0;

  /**
   * Like get_data() and get_face_data(), but in preparation for later calls
   * to transform() or fill_fe_subface_values() that will need information
   * about mappings from the reference face to a child of a face (i.e.,
   * subface) of a concrete cell.
   *
   * @param update_flags A set of flags that define what is expected of the
   * mapping class in future calls to transform() or the fill_fe_values()
   * group of functions. This set of flags may contain flags that mappings do
   * not know how to deal with (e.g., for information that is in fact computed
   * by the finite element classes, such as UpdateFlags::update_values).
   * Derived classes will need to store these flags, or at least that subset
   * of flags that will require the mapping to perform any actions in
   * fill_fe_values(), in InternalDataBase::update_each.
   * @param quadrature The quadrature object for which mapping information
   * will have to be computed. This includes the locations and weights of
   * quadrature points.
   * @return A pointer to a newly created object of type InternalDataBase (or
   * a derived class). Ownership of this object passes to the calling
   * function.
   *
   * @note C++ allows that virtual functions in derived classes may return
   * pointers to objects not of type InternalDataBase but in fact pointers to
   * objects of classes <i>derived</i> from InternalDataBase. (This feature is
   * called "covariant return types".) This is useful in some contexts where
   * the calling is within the derived class and will immediately make use of
   * the returned object, knowing its real (derived) type.
   */
  virtual
  InternalDataBase *
  get_subface_data (const UpdateFlags        update_flags,
                    const Quadrature<dim-1> &quadrature) const = 0;

  /**
   * Compute information about the mapping from the reference cell to the real
   * cell indicated by the first argument to this function. Derived classes
   * will have to implement this function based on the kind of mapping they
   * represent. It is called by FEValues::reinit().
   *
   * Conceptually, this function's represents the application of the mapping
   * $\mathbf x=\mathbf F_K(\hat {\mathbf x})$ from reference coordinates
   * $\mathbf\in [0,1]^d$ to real space coordinates $\mathbf x$ for a given
   * cell $K$. Its purpose is to compute the following kinds of data:
   *
   * - Data that results from the application of the mapping itself, e.g.,
   * computing the location $\mathbf x_q = \mathbf F_K(\hat{\mathbf x}_q)$ of
   * quadrature points on the real cell, and that is directly useful to users
   * of FEValues, for example during assembly.
   * - Data that is necessary for finite element implementations to compute
   * their shape functions on the real cell. To this end, the
   * FEValues::reinit() function calls FiniteElement::fill_fe_values() after
   * the current function, and the output of this function serves as input to
   * FiniteElement::fill_fe_values(). Examples of information that needs to be
   * computed here for use by the finite element classes is the Jacobian of
   * the mapping, $\hat\nabla \mathbf F_K(\hat{\mathbf x})$ or its inverse,
   * for example to transform the gradients of shape functions on the
   * reference cell to the gradients of shape functions on the real cell.
   *
   * The information computed by this function is used to fill the various
   * member variables of the output argument of this function. Which of the
   * member variables of that structure should be filled is determined by the
   * update flags stored in the Mapping::InternalDataBase object passed to
   * this function.
   *
   * An extensive discussion of the interaction between this function and
   * FEValues can be found in the
   * @ref FE_vs_Mapping_vs_FEValues
   * documentation module.
   *
   * @param[in] cell The cell of the triangulation for which this function is
   * to compute a mapping from the reference cell to.
   * @param[in] cell_similarity Whether or not the cell given as first
   * argument is simply a translation, rotation, etc of the cell for which
   * this function was called the most recent time. This information is
   * computed simply by matching the vertices (as stored by the Triangulation)
   * between the previous and the current cell. The value passed here may be
   * modified by implementations of this function and should then be returned
   * (see the discussion of the return value of this function).
   * @param[in] quadrature A reference to the quadrature formula in use for
   * the current evaluation. This quadrature object is the same as the one
   * used when creating the @p internal_data object. The object is used both
   * to map the location of quadrature points, as well as to compute the JxW
   * values for each quadrature point (which involves the quadrature weights).
   * @param[in] internal_data A reference to an object previously created by
   * get_data() and that may be used to store information the mapping can
   * compute once on the reference cell. See the documentation of the
   * Mapping::InternalDataBase class for an extensive description of the
   * purpose of these objects.
   * @param[out] output_data A reference to an object whose member variables
   * should be computed. Not all of the members of this argument need to be
   * filled; which ones need to be filled is determined by the update flags
   * stored inside the @p internal_data object.
   * @return An updated value of the @p cell_similarity argument to this
   * function. The returned value will be used for the corresponding argument
   * when FEValues::reinit() calls FiniteElement::fill_fe_values(). In most
   * cases, derived classes will simply want to return the value passed for @p
   * cell_similarity. However, implementations of this function may downgrade
   * the level of cell similarity. This is, for example, the case for classes
   * that take not only into account the locations of the vertices of a cell
   * (as reported by the Triangulation), but also other information specific
   * to the mapping. The purpose is that FEValues::reinit() can compute
   * whether a cell is similar to the previous one only based on the cell's
   * vertices, whereas the mapping may also consider displacement fields
   * (e.g., in the MappingQ1Eulerian and MappingFEField classes). In such
   * cases, the mapping may conclude that the previously computed cell
   * similarity is too optimistic, and invalidate it for subsequent use in
   * FiniteElement::fill_fe_values() by returning a less optimistic cell
   * similarity value.
   *
   * @note FEValues ensures that this function is always called with the same
   * pair of @p internal_data and @p output_data objects. In other words, if
   * an implementation of this function knows that it has written a piece of
   * data into the output argument in a previous call, then there is no need
   * to copy it there again in a later call if the implementation knows that
   * this is the same value.
   */
  virtual
  CellSimilarity::Similarity
  fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator    &cell,
                  const CellSimilarity::Similarity                              cell_similarity,
                  const Quadrature<dim>                                        &quadrature,
                  const typename Mapping<dim,spacedim>::InternalDataBase       &internal_data,
                  dealii::internal::FEValues::MappingRelatedData<dim,spacedim> &output_data) const = 0;

  /**
   * This function is the equivalent to Mapping::fill_fe_values(), but for
   * faces of cells. See there for an extensive discussion of its purpose. It
   * is called by FEFaceValues::reinit().
   *
   * @param[in] cell The cell of the triangulation for which this function is
   * to compute a mapping from the reference cell to.
   * @param[in] face_no The number of the face of the given cell for which
   * information is requested.
   * @param[in] quadrature A reference to the quadrature formula in use for
   * the current evaluation. This quadrature object is the same as the one
   * used when creating the @p internal_data object. The object is used both
   * to map the location of quadrature points, as well as to compute the JxW
   * values for each quadrature point (which involves the quadrature weights).
   * @param[in] internal_data A reference to an object previously created by
   * get_data() and that may be used to store information the mapping can
   * compute once on the reference cell. See the documentation of the
   * Mapping::InternalDataBase class for an extensive description of the
   * purpose of these objects.
   * @param[out] output_data A reference to an object whose member variables
   * should be computed. Not all of the members of this argument need to be
   * filled; which ones need to be filled is determined by the update flags
   * stored inside the @p internal_data object.
   */
  virtual void
  fill_fe_face_values (const typename Triangulation<dim,spacedim>::cell_iterator    &cell,
                       const unsigned int                                            face_no,
                       const Quadrature<dim-1>                                      &quadrature,
                       const typename Mapping<dim,spacedim>::InternalDataBase       &internal_data,
                       dealii::internal::FEValues::MappingRelatedData<dim,spacedim> &output_data) const = 0;

  /**
   * This function is the equivalent to Mapping::fill_fe_values(), but for
   * subfaces (i.e., children of faces) of cells. See there for an extensive
   * discussion of its purpose. It is called by FESubfaceValues::reinit().
   *
   * @param[in] cell The cell of the triangulation for which this function is
   * to compute a mapping from the reference cell to.
   * @param[in] face_no The number of the face of the given cell for which
   * information is requested.
   * @param[in] subface_no The number of the child of a face of the given cell
   * for which information is requested.
   * @param[in] quadrature A reference to the quadrature formula in use for
   * the current evaluation. This quadrature object is the same as the one
   * used when creating the @p internal_data object. The object is used both
   * to map the location of quadrature points, as well as to compute the JxW
   * values for each quadrature point (which involves the quadrature weights).
   * @param[in] internal_data A reference to an object previously created by
   * get_data() and that may be used to store information the mapping can
   * compute once on the reference cell. See the documentation of the
   * Mapping::InternalDataBase class for an extensive description of the
   * purpose of these objects.
   * @param[out] output_data A reference to an object whose member variables
   * should be computed. Not all of the members of this argument need to be
   * filled; which ones need to be filled is determined by the update flags
   * stored inside the @p internal_data object.
   */
  virtual void
  fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator     &cell,
                          const unsigned int                                             face_no,
                          const unsigned int                                             subface_no,
                          const Quadrature<dim-1>                                       &quadrature,
                          const typename Mapping<dim,spacedim>::InternalDataBase        &internal_data,
                          dealii::internal::FEValues::MappingRelatedData<dim, spacedim> &output_data) const = 0;

  /**
   * @}
   */

public:

  /**
   * @name Functions to transform tensors from reference to real coordinates
   * @{
   */

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
   * to the physical cell through the Jacobian:
   * @f[
   * \mathbf u(\mathbf x) = J(\hat{\mathbf  x})\hat{\mathbf  u}(\hat{\mathbf  x}).
   * @f]
   * In physics, this is usually referred to as the contravariant
   * transformation. Mathematically, it is the push forward of a vector field.
   *
   * <li> @p mapping_covariant: maps a field of one-forms on the reference
   * cell to a field of one-forms on the physical cell. (Theoretically this
   * would refer to a DerivativeForm<1,dim,1> but we canonically identify this
   * type with a Tensor<1,dim>). Mathematically, it is the pull back of the
   * differential form
   * @f[
   * \mathbf u(\mathbf x) = J(\hat{\mathbf  x})(J(\hat{\mathbf  x})^{T} J(\hat{\mathbf  x}))^{-1}\hat{\mathbf
   * u}(\hat{\mathbf  x}).
   * @f]
   * Gradients of scalar differentiable functions are transformed this way.
   *
   * In the case when dim=spacedim the previous formula reduces to
   * @f[
   * \mathbf u(\mathbf x) = J(\hat{\mathbf  x})^{-T}\hat{\mathbf
   * u}(\hat{\mathbf  x})
   * @f]
   * because we assume that the mapping $\mathbf F_K$ is always invertible,
   * and consequently its Jacobian $J$ is an invertible matrix.
   *
   * <li> @p mapping_piola: A field of <i>dim-1</i>-forms on the reference
   * cell is also represented by a vector field, but again transforms
   * differently, namely by the Piola transform
   * @f[
   *  \mathbf u(\mathbf x) = \frac{1}{\text{det}\;J(\mathbf x)}
   * J(\mathbf x) \hat{\mathbf  u}(\mathbf x).
   * @f]
   * </ul>
   *
   * @param[in] input An array (or part of an array) of input objects that
   * should be mapped.
   * @param[in] type The kind of mapping to be applied.
   * @param[in] internal A pointer to an object of type
   * Mapping::InternalDataBase that contains information previously stored by
   * the mapping. The object pointed to was created by the get_data(),
   * get_face_data(), or get_subface_data() function, and will have been
   * updated as part of a call to fill_fe_values(), fill_fe_face_values(), or
   * fill_fe_subface_values() for the current cell, before calling the current
   * function. In other words, this object also represents with respect to
   * which cell the transformation should be applied to.
   * @param[out] output An array (or part of an array) into which the
   * transformed objects should be placed. (Note that the array view is @p
   * const, but the tensors it points to are not.)
   */
  virtual
  void
  transform (const ArrayView<const Tensor<1,dim> >                  &input,
             const MappingType                                       type,
             const typename Mapping<dim,spacedim>::InternalDataBase &internal,
             const ArrayView<Tensor<1,spacedim> >                   &output) const = 0;

  /**
   * Transform a field of differential forms from the reference cell to the
   * physical cell.  It is useful to think of $\mathbf{T} = \nabla \mathbf u$
   * and $\hat{\mathbf  T} = \hat \nabla \hat{\mathbf  u}$, with $\mathbf u$ a
   * vector field.  The mapping types currently implemented by derived classes
   * are:
   * <ul>
   * <li> @p mapping_covariant: maps a field of forms on the reference cell to
   * a field of forms on the physical cell. Mathematically, it is the pull
   * back of the differential form
   * @f[
   * \mathbf T(\mathbf x) = \hat{\mathbf  T}(\hat{\mathbf  x})
   *                        J(\hat{\mathbf  x})(J(\hat{\mathbf  x})^{T} J(\hat{\mathbf  x}))^{-1}.
   * @f]
   * Jacobians of spacedim-vector valued differentiable functions are
   * transformed this way.
   *
   * In the case when dim=spacedim the previous formula reduces to
   * @f[
   * \mathbf T(\mathbf x) = \hat{\mathbf  u}(\hat{\mathbf  x})
   *                        J(\hat{\mathbf  x})^{-1}.
   * @f]
   * </ul>
   *
   * @note It would have been more reasonable to make this transform a
   * template function with the rank in <code>DerivativeForm@<1, dim,
   * rank@></code>. Unfortunately C++ does not allow templatized virtual
   * functions. This is why we identify <code>DerivativeForm@<1, dim,
   * 1@></code> with a <code>Tensor@<1,dim@></code> when using
   * mapping_covariant() in the function transform() above this one.
   *
   * @param[in] input An array (or part of an array) of input objects that
   * should be mapped.
   * @param[in] type The kind of mapping to be applied.
   * @param[in] internal A pointer to an object of type
   * Mapping::InternalDataBase that contains information previously stored by
   * the mapping. The object pointed to was created by the get_data(),
   * get_face_data(), or get_subface_data() function, and will have been
   * updated as part of a call to fill_fe_values(), fill_fe_face_values(), or
   * fill_fe_subface_values() for the current cell, before calling the current
   * function. In other words, this object also represents with respect to
   * which cell the transformation should be applied to.
   * @param[out] output An array (or part of an array) into which the
   * transformed objects should be placed. (Note that the array view is @p
   * const, but the tensors it points to are not.)
   */
  virtual
  void
  transform (const ArrayView<const DerivativeForm<1, dim, spacedim> > &input,
             const MappingType                                         type,
             const typename Mapping<dim,spacedim>::InternalDataBase   &internal,
             const ArrayView<Tensor<2,spacedim> >                     &output) const = 0;

  /**
   * Transform a tensor field from the reference cell to the physical cell.
   * These tensors are usually the Jacobians in the reference cell of vector
   * fields that have been pulled back from the physical cell.  The mapping
   * types currently implemented by derived classes are:
   * <ul>
   * <li> @p mapping_contravariant_gradient: it assumes $\mathbf u(\mathbf x)
   * = J \hat{\mathbf  u}$ so that
   * @f[
   * \mathbf T(\mathbf x) =
   * J(\hat{\mathbf  x}) \hat{\mathbf  T}(\hat{\mathbf  x})
   * J(\hat{\mathbf  x})^{-1}.
   * @f]
   * <li> @p mapping_covariant_gradient: it assumes $\mathbf u(\mathbf x) =
   * J^{-T} \hat{\mathbf  u}$ so that
   * @f[
   * \mathbf T(\mathbf x) =
   * J(\hat{\mathbf  x})^{-T} \hat{\mathbf  T}(\hat{\mathbf  x})
   * J(\hat{\mathbf  x})^{-1}.
   * @f]
   * <li> @p mapping_piola_gradient: it assumes $\mathbf u(\mathbf x) =
   * \frac{1}{\text{det}\;J(\mathbf x)} J(\mathbf x) \hat{\mathbf  u}(\mathbf
   * x)$ so that
   * @f[
   * \mathbf T(\mathbf x) =
   * \frac{1}{\text{det}\;J(\mathbf x)}
   * J(\hat{\mathbf  x}) \hat{\mathbf  T}(\hat{\mathbf  x})
   * J(\hat{\mathbf  x})^{-1}.
   * @f]
   * </ul>
   *
   * @todo The formulas for mapping_covariant_gradient,
   * mapping_contravariant_gradient and mapping_piola_gradient are only true
   * as stated for linear mappings. If, for example, the mapping is bilinear
   * (or has a higher order polynomial degree) then there is a missing term
   * associated with the derivative of $J$.
   *
   * @param[in] input An array (or part of an array) of input objects that
   * should be mapped.
   * @param[in] type The kind of mapping to be applied.
   * @param[in] internal A pointer to an object of type
   * Mapping::InternalDataBase that contains information previously stored by
   * the mapping. The object pointed to was created by the get_data(),
   * get_face_data(), or get_subface_data() function, and will have been
   * updated as part of a call to fill_fe_values(), fill_fe_face_values(), or
   * fill_fe_subface_values() for the current cell, before calling the current
   * function. In other words, this object also represents with respect to
   * which cell the transformation should be applied to.
   * @param[out] output An array (or part of an array) into which the
   * transformed objects should be placed. (Note that the array view is @p
   * const, but the tensors it points to are not.)
   */
  virtual
  void
  transform (const ArrayView<const Tensor<2, dim> >                 &input,
             const MappingType                                       type,
             const typename Mapping<dim,spacedim>::InternalDataBase &internal,
             const ArrayView<Tensor<2,spacedim> >                   &output) const = 0;

  /**
   * Transform a tensor field from the reference cell to the physical cell.
   * This tensors are most of times the hessians in the reference cell of
   * vector fields that have been pulled back from the physical cell.
   *
   * The mapping types currently implemented by derived classes are:
   * <ul>
   * <li> @p mapping_covariant_gradient: maps a field of forms on the
   * reference cell to a field of forms on the physical cell. Mathematically,
   * it is the pull back of the differential form
   * @f[
   * \mathbf T_{ijk}(\mathbf x) = \hat{\mathbf  T}_{iJK}(\hat{\mathbf  x}) J_{jJ}^{\dagger} J_{kK}^{\dagger}@f],
   *
   * where @f[ J^{\dagger} = J(\hat{\mathbf  x})(J(\hat{\mathbf  x})^{T} J(\hat{\mathbf  x}))^{-1}.
   * @f]
   * </ul>
   *
   * Hessians of spacedim-vector valued differentiable functions are
   * transformed this way (After subtraction of the product of the derivative
   * with the Jacobian gradient).
   *
   * In the case when dim=spacedim the previous formula reduces to
   * @f[J^{\dagger} = J^{-1}@f]
   *
   * @param[in] input An array (or part of an array) of input objects that
   * should be mapped.
   * @param[in] type The kind of mapping to be applied.
   * @param[in] internal A pointer to an object of type
   * Mapping::InternalDataBase that contains information previously stored by
   * the mapping. The object pointed to was created by the get_data(),
   * get_face_data(), or get_subface_data() function, and will have been
   * updated as part of a call to fill_fe_values(), fill_fe_face_values(), or
   * fill_fe_subface_values() for the current cell, before calling the current
   * function. In other words, this object also represents with respect to
   * which cell the transformation should be applied to.
   * @param[out] output An array (or part of an array) into which the
   * transformed objects should be placed. (Note that the array view is @p
   * const, but the tensors it points to are not.)
   */
  virtual
  void
  transform (const ArrayView<const DerivativeForm<2, dim, spacedim> > &input,
             const MappingType                                         type,
             const typename Mapping<dim,spacedim>::InternalDataBase   &internal,
             const ArrayView<Tensor<3,spacedim> >                     &output) const = 0;

  /**
   * Transform a field of 3-differential forms from the reference cell to the
   * physical cell.  It is useful to think of $\mathbf{T}_{ijk} = D^2_{jk}
   * \mathbf u_i$ and $\mathbf{\hat T}_{IJK} = \hat D^2_{JK} \mathbf{\hat
   * u}_I$, with $\mathbf u_i$ a vector field.
   *
   * The mapping types currently implemented by derived classes are:
   * <ul>
   * <li> @p mapping_contravariant_hessian: it assumes $\mathbf u_i(\mathbf x)
   * = J_{iI} \hat{\mathbf  u}_I$ so that
   * @f[
   * \mathbf T_{ijk}(\mathbf x) =
   * J_{iI}(\hat{\mathbf  x}) \hat{\mathbf  T}_{IJK}(\hat{\mathbf  x})
   * J_{jJ}(\hat{\mathbf  x})^{-1} J_{kK}(\hat{\mathbf  x})^{-1}.
   * @f]
   * <li> @p mapping_covariant_hessian: it assumes $\mathbf u_i(\mathbf x) =
   * J_{iI}^{-T} \hat{\mathbf  u}_I$ so that
   * @f[
   * \mathbf T_{ijk}(\mathbf x) =
   * J_iI(\hat{\mathbf  x})^{-1} \hat{\mathbf  T}_{IJK}(\hat{\mathbf  x})
   * J_{jJ}(\hat{\mathbf  x})^{-1} J_{kK}(\hat{\mathbf  x})^{-1}.
   * @f]
   * <li> @p mapping_piola_hessian: it assumes $\mathbf u_i(\mathbf x) =
   * \frac{1}{\text{det}\;J(\mathbf x)} J_{iI}(\mathbf x) \hat{\mathbf
   * u}(\mathbf x)$ so that
   * @f[
   * \mathbf T_{ijk}(\mathbf x) =
   * \frac{1}{\text{det}\;J(\mathbf x)}
   * J_{iI}(\hat{\mathbf  x}) \hat{\mathbf  T}_{IJK}(\hat{\mathbf  x})
   * J_{jJ}(\hat{\mathbf  x})^{-1} J_{kK}(\hat{\mathbf  x})^{-1}.
   * @f]
   * </ul>
   *
   * @param[in] input An array (or part of an array) of input objects that
   * should be mapped.
   * @param[in] type The kind of mapping to be applied.
   * @param[in] internal A pointer to an object of type
   * Mapping::InternalDataBase that contains information previously stored by
   * the mapping. The object pointed to was created by the get_data(),
   * get_face_data(), or get_subface_data() function, and will have been
   * updated as part of a call to fill_fe_values(), fill_fe_face_values(), or
   * fill_fe_subface_values() for the current cell, before calling the current
   * function. In other words, this object also represents with respect to
   * which cell the transformation should be applied to.
   * @param[out] output An array (or part of an array) into which the
   * transformed objects should be placed.
   */
  virtual
  void
  transform (const ArrayView<const Tensor<3, dim> >                 &input,
             const MappingType                                       type,
             const typename Mapping<dim,spacedim>::InternalDataBase &internal,
             const ArrayView<Tensor<3,spacedim> >                   &output) const = 0;

  /**
   * @}
   */


  /**
   * Give class @p FEValues access to the private <tt>get_...data</tt> and
   * <tt>fill_fe_...values</tt> functions.
   */
  friend class FEValuesBase<dim,spacedim>;
  friend class FEValues<dim,spacedim>;
  friend class FEFaceValues<dim,spacedim>;
  friend class FESubfaceValues<dim,spacedim>;
};


DEAL_II_NAMESPACE_CLOSE

#endif

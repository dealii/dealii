// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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

#ifndef __deal2__fe_poly_h
#define __deal2__fe_poly_h


#include <deal.II/fe/fe.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup febase */
/*@{*/

/**
 * This class gives a unified framework for the implementation of
 * FiniteElement classes based on polynomial spaces like the
 * TensorProductPolynomials or PolynomialSpace classes.
 *
 * Every class conforming to the following interface can be used as
 * template parameter POLY.
 *
 * @code
 * static const unsigned int dimension;
 *
 * double compute_value (const unsigned int i,
 *                       const Point<dim> &p) const;
 *
 * Tensor<1,dim> compute_grad (const unsigned int i,
 *                             const Point<dim> &p) const;
 *
 * Tensor<2,dim> compute_grad_grad (const unsigned int i,
 *                                  const Point<dim> &p) const;
 * @endcode
 * Example classes are TensorProductPolynomials, PolynomialSpace or
 * PolynomialsP.
 *
 * This class is not a fully implemented FiniteElement class. Instead
 * there are several pure virtual functions declared in the
 * FiniteElement and FiniteElement classes which cannot
 * implemented by this class but are left for implementation in
 * derived classes.
 *
 * Furthermore, this class assumes that shape functions of the
 * FiniteElement under consideration do <em>not</em> depend on the
 * actual shape of the cells in real space, i.e. update_once()
 * includes <tt>update_values</tt>. For FiniteElements whose shape
 * functions depend on the cells in real space, the update_once() and
 * update_each() functions must be overloaded.
 *
 * @todo Since nearly all functions for spacedim != dim are
 * specialized, this class needs cleaning up.
 *
 * @author Ralf Hartmann 2004, Guido Kanschat, 2009
 **/
template <class POLY, int dim=POLY::dimension, int spacedim=dim>
class FE_Poly : public FiniteElement<dim,spacedim>
{
public:
  /**
   * Constructor.
   */
  FE_Poly (const POLY &poly_space,
           const FiniteElementData<dim> &fe_data,
           const std::vector<bool> &restriction_is_additive_flags,
           const std::vector<ComponentMask> &nonzero_components);

  /**
   * Return the polynomial degree of this finite element, i.e. the value
   * passed to the constructor.
   */
  unsigned int get_degree () const;

  /**
   * Return the numbering of the underlying polynomial space compared to
   * lexicographic ordering of the basis functions. Returns
   * POLY::get_numbering().
   */
  std::vector<unsigned int> get_poly_space_numbering() const;

  /**
   * Return the inverse numbering of the underlying polynomial space. Returns
   * POLY::get_numbering_inverse().
   */
  std::vector<unsigned int> get_poly_space_numbering_inverse() const;

  /**
   * Return the value of the <tt>i</tt>th shape function at the point
   * <tt>p</tt>. See the FiniteElement base class for more information about
   * the semantics of this function.
   */
  virtual double shape_value (const unsigned int i,
                              const Point<dim> &p) const;

  /**
   * Return the value of the <tt>component</tt>th vector component of the
   * <tt>i</tt>th shape function at the point <tt>p</tt>. See the
   * FiniteElement base class for more information about the semantics of this
   * function.
   *
   * Since this element is scalar, the returned value is the same as if the
   * function without the <tt>_component</tt> suffix were called, provided
   * that the specified component is zero.
   */
  virtual double shape_value_component (const unsigned int i,
                                        const Point<dim> &p,
                                        const unsigned int component) const;

  /**
   * Return the gradient of the <tt>i</tt>th shape function at the point
   * <tt>p</tt>. See the FiniteElement base class for more information about
   * the semantics of this function.
   */
  virtual Tensor<1,dim> shape_grad (const unsigned int  i,
                                    const Point<dim>   &p) const;

  /**
   * Return the gradient of the <tt>component</tt>th vector component of the
   * <tt>i</tt>th shape function at the point <tt>p</tt>. See the
   * FiniteElement base class for more information about the semantics of this
   * function.
   *
   * Since this element is scalar, the returned value is the same as if the
   * function without the <tt>_component</tt> suffix were called, provided
   * that the specified component is zero.
   */
  virtual Tensor<1,dim> shape_grad_component (const unsigned int i,
                                              const Point<dim> &p,
                                              const unsigned int component) const;

  /**
   * Return the tensor of second derivatives of the <tt>i</tt>th shape
   * function at point <tt>p</tt> on the unit cell. See the FiniteElement base
   * class for more information about the semantics of this function.
   */
  virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
                                         const Point<dim> &p) const;

  /**
   * Return the second derivative of the <tt>component</tt>th vector component
   * of the <tt>i</tt>th shape function at the point <tt>p</tt>. See the
   * FiniteElement base class for more information about the semantics of this
   * function.
   *
   * Since this element is scalar, the returned value is the same as if the
   * function without the <tt>_component</tt> suffix were called, provided
   * that the specified component is zero.
   */
  virtual Tensor<2,dim> shape_grad_grad_component (const unsigned int i,
                                                   const Point<dim> &p,
                                                   const unsigned int component) const;

protected:

  virtual
  typename Mapping<dim,spacedim>::InternalDataBase *
  get_data (const UpdateFlags,
            const Mapping<dim,spacedim> &mapping,
            const Quadrature<dim> &quadrature) const ;

  virtual void
  fill_fe_values (const Mapping<dim,spacedim>                           &mapping,
                  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                  const Quadrature<dim>                                 &quadrature,
                  typename Mapping<dim,spacedim>::InternalDataBase      &mapping_internal,
                  typename Mapping<dim,spacedim>::InternalDataBase      &fe_internal,
                  FEValuesData<dim,spacedim>                            &data,
                  CellSimilarity::Similarity                       &cell_similarity) const;

  virtual void
  fill_fe_face_values (const Mapping<dim,spacedim> &mapping,
                       const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                       const unsigned int                    face_no,
                       const Quadrature<dim-1>                &quadrature,
                       typename Mapping<dim,spacedim>::InternalDataBase      &mapping_internal,
                       typename Mapping<dim,spacedim>::InternalDataBase      &fe_internal,
                       FEValuesData<dim,spacedim> &data) const ;

  virtual void
  fill_fe_subface_values (const Mapping<dim,spacedim> &mapping,
                          const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                          const unsigned int                    face_no,
                          const unsigned int                    sub_no,
                          const Quadrature<dim-1>                &quadrature,
                          typename Mapping<dim,spacedim>::InternalDataBase      &mapping_internal,
                          typename Mapping<dim,spacedim>::InternalDataBase      &fe_internal,
                          FEValuesData<dim,spacedim> &data) const ;


  /**
   * Determine the values that need to be computed on the unit cell to be able
   * to compute all values required by <tt>flags</tt>.
   *
   * For the purpuse of this function, refer to the documentation in
   * FiniteElement.
   *
   * This class assumes that shape functions of this FiniteElement do
   * <em>not</em> depend on the actual shape of the cells in real
   * space. Therefore, the effect in this element is as follows: if
   * <tt>update_values</tt> is set in <tt>flags</tt>, copy it to the
   * result. All other flags of the result are cleared, since everything else
   * must be computed for each cell.
   */
  virtual UpdateFlags update_once (const UpdateFlags flags) const;

  /**
   * Determine the values that need to be computed on every cell to be able to
   * compute all values required by <tt>flags</tt>.
   *
   * For the purpuse of this function, refer to the documentation in
   * FiniteElement.
   *
   * This class assumes that shape functions of this FiniteElement do
   * <em>not</em> depend on the actual shape of the cells in real space.
   *
   * The effect in this element is as follows:
   * <ul>
   *
   * <li> if <tt>update_gradients</tt> is set, the result will contain
   * <tt>update_gradients</tt> and <tt>update_covariant_transformation</tt>.
   * The latter is required to transform the gradient on the unit cell to the
   * real cell. Remark, that the action required by
   * <tt>update_covariant_transformation</tt> is actually performed by the
   * Mapping object used in conjunction with this finite element.
   *
   * <li> if <tt>update_hessians</tt> is set, the result will contain
   * <tt>update_hessians</tt> and <tt>update_covariant_transformation</tt>.
   * The rationale is the same as above and no higher derivatives of the
   * transformation are required, since we use difference quotients for the
   * actual computation.
   *
   * </ul>
   */
  virtual UpdateFlags update_each (const UpdateFlags flags) const;


  /**
   * Fields of cell-independent data.
   *
   * For information about the general purpose of this class, see the
   * documentation of the base class.
   */
  class InternalData : public FiniteElement<dim,spacedim>::InternalDataBase
  {
  public:
    /**
     * Array with shape function values in quadrature points. There is one row
     * for each shape function, containing values for each quadrature point.
     *
     * In this array, we store the values of the shape function in the
     * quadrature points on the unit cell. Since these values do not change
     * under transformation to the real cell, we only need to copy them over
     * when visiting a concrete cell.
     */
    std::vector<std::vector<double> > shape_values;

    /**
     * Array with shape function gradients in quadrature points. There is one
     * row for each shape function, containing values for each quadrature
     * point.
     *
     * We store the gradients in the quadrature points on the unit cell. We
     * then only have to apply the transformation (which is a matrix-vector
     * multiplication) when visiting an actual cell.
     */
    std::vector<std::vector<Tensor<1,dim> > > shape_gradients;
  };

  /**
   * The polynomial space. Its type is given by the template parameter POLY.
   */
  POLY poly_space;
};

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif

// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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

#ifndef __deal2__fe_poly_face_h
#define __deal2__fe_poly_face_h


#include <deal.II/fe/fe.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup febase */
/*@{*/

/**
 * @warning This class is not sufficiently tested yet!
 *
 * This class gives a unified framework for the implementation of
 * FiniteElement classes only located on faces of the mesh. They are
 * based on polynomial spaces like the TensorProductPolynomials or a
 * PolynomialSpace classes.
 *
 * Every class that implements the following functions can be used as
 * template parameter POLY.
 *
 * @code
 * double compute_value (const unsigned int i,
 *                       const Point<dim> &p) const;
 * @endcode
 * Example classes are TensorProductPolynomials, PolynomialSpace or
 * PolynomialsP.
 *
 * This class is not a fully implemented FiniteElement class. Instead there
 * are several pure virtual functions declared in the FiniteElement class
 * which cannot be implemented by this class but are left for implementation
 * in derived classes.
 *
 * Furthermore, this class assumes that shape functions of the
 * FiniteElement under consideration do <em>not</em> depend on the
 * actual shape of the cells in real space, i.e. update_once()
 * includes <tt>update_values</tt>. For FiniteElements whose shape
 * functions depend on the cells in real space, the update_once() and
 * update_each() functions must be overloaded.
 *
 * @author Guido Kanschat, 2009
 **/
template <class POLY, int dim=POLY::dimension+1, int spacedim=dim>
class FE_PolyFace : public FiniteElement<dim,spacedim>
{
public:
  /**
   * Constructor.
   */
  FE_PolyFace (const POLY &poly_space,
               const FiniteElementData<dim> &fe_data,
               const std::vector<bool> &restriction_is_additive_flags);

  /**
   * Return the polynomial degree
   * of this finite element,
   * i.e. the value passed to the
   * constructor.
   */
  unsigned int get_degree () const;

  /**
   * Return the value of the
   * <tt>i</tt>th shape function at
   * the point <tt>p</tt>. See the
   * FiniteElement base class
   * for more information about the
   * semantics of this function.
   */
//    virtual double shape_value (const unsigned int i,
//                              const Point<dim> &p) const;

  /**
   * Return the value of the
   * <tt>component</tt>th vector
   * component of the <tt>i</tt>th
   * shape function at the point
   * <tt>p</tt>. See the
   * FiniteElement base class
   * for more information about the
   * semantics of this function.
   *
   * Since this element is scalar,
   * the returned value is the same
   * as if the function without the
   * <tt>_component</tt> suffix
   * were called, provided that the
   * specified component is zero.
   */
//    virtual double shape_value_component (const unsigned int i,
//                                        const Point<dim> &p,
//                                        const unsigned int component) const;

  /**
   * Return the gradient of the
   * <tt>i</tt>th shape function at
   * the point <tt>p</tt>. See the
   * FiniteElement base class
   * for more information about the
   * semantics of this function.
   */
//    virtual Tensor<1,dim> shape_grad (const unsigned int  i,
//                                    const Point<dim>   &p) const;

  /**
   * Return the gradient of the
   * <tt>component</tt>th vector
   * component of the <tt>i</tt>th
   * shape function at the point
   * <tt>p</tt>. See the
   * FiniteElement base class
   * for more information about the
   * semantics of this function.
   *
   * Since this element is scalar,
   * the returned value is the same
   * as if the function without the
   * <tt>_component</tt> suffix
   * were called, provided that the
   * specified component is zero.
   */
//    virtual Tensor<1,dim> shape_grad_component (const unsigned int i,
//                                              const Point<dim> &p,
//                                              const unsigned int component) const;

  /**
   * Return the tensor of second
   * derivatives of the
   * <tt>i</tt>th shape function at
   * point <tt>p</tt> on the unit
   * cell. See the
   * FiniteElement base class
   * for more information about the
   * semantics of this function.
   */
//    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
//                                         const Point<dim> &p) const;

  /**
   * Return the second derivative
   * of the <tt>component</tt>th
   * vector component of the
   * <tt>i</tt>th shape function at
   * the point <tt>p</tt>. See the
   * FiniteElement base class
   * for more information about the
   * semantics of this function.
   *
   * Since this element is scalar,
   * the returned value is the same
   * as if the function without the
   * <tt>_component</tt> suffix
   * were called, provided that the
   * specified component is zero.
   */
//    virtual Tensor<2,dim> shape_grad_grad_component (const unsigned int i,
//                                                   const Point<dim> &p,
//                                                   const unsigned int component) const;

protected:

  virtual
  typename Mapping<dim,spacedim>::InternalDataBase *
  get_data (const UpdateFlags,
            const Mapping<dim,spacedim> &mapping,
            const Quadrature<dim> &quadrature) const ;

  typename Mapping<dim,spacedim>::InternalDataBase *
  get_face_data (const UpdateFlags,
                 const Mapping<dim,spacedim> &mapping,
                 const Quadrature<dim-1>& quadrature) const ;

  typename Mapping<dim,spacedim>::InternalDataBase *
  get_subface_data (const UpdateFlags,
                    const Mapping<dim,spacedim> &mapping,
                    const Quadrature<dim-1>& quadrature) const ;

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
   * Determine the values that need
   * to be computed on the unit
   * cell to be able to compute all
   * values required by
   * <tt>flags</tt>.
   *
   * For the purpuse of this
   * function, refer to the
   * documentation in
   * FiniteElement.
   *
   * This class assumes that shape
   * functions of this
   * FiniteElement do <em>not</em>
   * depend on the actual shape of
   * the cells in real
   * space. Therefore, the effect
   * in this element is as follows:
   * if <tt>update_values</tt> is
   * set in <tt>flags</tt>, copy it
   * to the result. All other flags
   * of the result are cleared,
   * since everything else must be
   * computed for each cell.
   */
  virtual UpdateFlags update_once (const UpdateFlags flags) const;

  /**
   * Determine the values that need
   * to be computed on every cell
   * to be able to compute all
   * values required by
   * <tt>flags</tt>.
   *
   * For the purpuse of this
   * function, refer to the
   * documentation in
   * FiniteElement.
   *
   * This class assumes that shape
   * functions of this
   * FiniteElement do <em>not</em>
   * depend on the actual shape of
   * the cells in real
   * space.
   *
   * The effect in this element is
   * as follows:
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
   * For information about the
   * general purpose of this class,
   * see the documentation of the
   * base class.
   */
  class InternalData : public FiniteElement<dim,spacedim>::InternalDataBase
  {
  public:
    /**
     * Array with shape function
     * values in quadrature
     * points on one face. There is one
     * row for each shape
     * function, containing
     * values for each quadrature
     * point.
     *
     * In this array, we store
     * the values of the shape
     * function in the quadrature
     * points on one face of the unit
     * cell. Since these values
     * do not change under
     * transformation to the real
     * cell, we only need to copy
     * them over when visiting a
     * concrete cell.
     *
     * In particular, we can
     * simply copy the same set
     * of values to each of the
     * faces.
     */
    std::vector<std::vector<double> > shape_values;
  };

  /**
   * The polynomial space. Its type
   * is given by the template
   * parameter POLY.
   */
  POLY poly_space;
};

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif

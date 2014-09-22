// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2014 by the deal.II authors
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

#ifndef __deal2__fe_abf_h
#define __deal2__fe_abf_h

#include <deal.II/base/config.h>
#include <deal.II/base/table.h>
#include <deal.II/base/polynomials_abf.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_poly_tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim> class MappingQ;


/*!@addtogroup fe */
/*@{*/

/**
 * Implementation of Arnold-Boffi-Falk (ABF) elements, conforming with the
 * space H<sup>div</sup>. These elements generate vector fields with
 * normal components continuous between mesh cells.
 *
 * These elements are based on an article from Arnold, Boffi and Falk:
 * Quadrilateral H(div) finite elements, SIAM J. Numer. Anal.
 * Vol.42, No.6, pp.2429-2451
 *
 * In this article, the authors demonstrate that the usual RT elements
 * and also BDM and other proposed finite dimensional subspaces of
 * H(div) do not work properly on arbitrary FE grids. I.e. the
 * convergence rates deteriorate on these meshes. As a solution the
 * authors propose the ABF elements, which are implemented in this
 * module.
 *
 * This class is not implemented for the codimension one case
 * (<tt>spacedim != dim</tt>).
 *
 * @todo Even if this element is implemented for two and three space
 * dimensions, the definition of the node values relies on
 * consistently oriented faces in 3D. Therefore, care should be taken
 * on complicated meshes.
 *
 * <h3>Interpolation</h3>
 *
 * The @ref GlossInterpolation "interpolation" operators associated
 * with the RT element are constructed such that interpolation and
 * computing the divergence are commuting operations. We require this
 * from interpolating arbitrary functions as well as the #restriction
 * matrices.  It can be achieved by two interpolation schemes, the
 * simplified one in FE_RaviartThomasNodal and the original one here:
 *
 * <h4>Node values on edges/faces</h4>
 *
 * On edges or faces, the @ref GlossNodes "node values" are the moments of
 * the normal component of the interpolated function with respect to
 * the traces of the RT polynomials. Since the normal trace of the RT
 * space of degree <i>k</i> on an edge/face is the space
 * <i>Q<sub>k</sub></i>, the moments are taken with respect to this
 * space.
 *
 * <h4>Interior node values</h4>
 *
 * Higher order RT spaces have interior nodes. These are moments taken
 * with respect to the gradient of functions in <i>Q<sub>k</sub></i>
 * on the cell (this space is the matching space for RT<sub>k</sub> in
 * a mixed formulation).
 *
 * <h4>Generalized support points</h4>
 *
 * The node values above rely on integrals, which will be computed by
 * quadrature rules themselves. The generalized support points are a
 * set of points such that this quadrature can be performed with
 * sufficient accuracy. The points needed are those of
 * QGauss<sub>k+1</sub> on each face as well as QGauss<sub>k</sub> in
 * the interior of the cell (or none for RT<sub>0</sub>).
 * See the @ref GlossGeneralizedSupport "glossary entry on generalized support points"
 * for more information.
 *
 *
 * @author Oliver Kayser-Herold, 2006, based on previous work
 * by Guido Kanschat and Wolfgang Bangerth
 */
template <int dim>
class FE_ABF : public FE_PolyTensor<PolynomialsABF<dim>, dim>
{
public:
  /**
   * Constructor for the ABF
   * element of degree @p p.
   */
  FE_ABF (const unsigned int p);

  /**
   * Return a string that uniquely
   * identifies a finite
   * element. This class returns
   * <tt>FE_ABF<dim>(degree)</tt>, with
   * @p dim and @p degree
   * replaced by appropriate
   * values.
   */
  virtual std::string get_name () const;

  /**
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero function values somewhere on the face @p face_index.
   *
   * Right now, this is only implemented for RT0 in 1D. Otherwise, returns
   * always @p true.
   */
  virtual bool has_support_on_face (const unsigned int shape_index,
                                    const unsigned int face_index) const;

  virtual void interpolate(std::vector<double>                &local_dofs,
                           const std::vector<double> &values) const;
  virtual void interpolate(std::vector<double>                &local_dofs,
                           const std::vector<Vector<double> > &values,
                           unsigned int offset = 0) const;
  virtual void interpolate(
    std::vector<double> &local_dofs,
    const VectorSlice<const std::vector<std::vector<double> > > &values) const;
  virtual std::size_t memory_consumption () const;
  virtual FiniteElement<dim> *clone() const;

private:
  /**
   * The order of the
   * ABF element. The
   * lowest order elements are
   * usually referred to as RT0,
   * even though their shape
   * functions are piecewise
   * quadratics.
   */
  const unsigned int rt_order;

  /**
   * Only for internal use. Its
   * full name is
   * @p get_dofs_per_object_vector
   * function and it creates the
   * @p dofs_per_object vector that is
   * needed within the constructor to
   * be passed to the constructor of
   * @p FiniteElementData.
   */
  static std::vector<unsigned int>
  get_dpo_vector (const unsigned int degree);

  /**
   * Initialize the @p
   * generalized_support_points
   * field of the FiniteElement
   * class and fill the tables with
   * interpolation weights
   * (#boundary_weights and
   * #interior_weights). Called
   * from the constructor.
  *
  * See the @ref GlossGeneralizedSupport "glossary entry on generalized support points"
  * for more information.
   */
  void initialize_support_points (const unsigned int rt_degree);

  /**
   * Initialize the interpolation
   * from functions on refined mesh
   * cells onto the father
   * cell. According to the
   * philosophy of the
   * Raviart-Thomas element, this
   * restriction operator preserves
   * the divergence of a function
   * weakly.
   */
  void initialize_restriction ();

  /**
   * Given a set of flags indicating
   * what quantities are requested
   * from a @p FEValues object,
   * return which of these can be
   * precomputed once and for
   * all. Often, the values of
   * shape function at quadrature
   * points can be precomputed, for
   * example, in which case the
   * return value of this function
   * would be the logical and of
   * the input @p flags and
   * @p update_values.
   *
   * For the present kind of finite
   * element, this is exactly the
   * case.
   */
  virtual UpdateFlags update_once (const UpdateFlags flags) const;

  /**
   * This is the opposite to the
   * above function: given a set of
   * flags indicating what we want
   * to know, return which of these
   * need to be computed each time
   * we visit a new cell.
   *
   * If for the computation of one
   * quantity something else is
   * also required (for example, we
   * often need the covariant
   * transformation when gradients
   * need to be computed), include
   * this in the result as well.
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
  class InternalData : public FiniteElement<dim>::InternalDataBase
  {
  public:
    /**
     * Array with shape function
     * values in quadrature
     * points. There is one row
     * for each shape function,
     * containing values for each
     * quadrature point. Since
     * the shape functions are
     * vector-valued (with as
     * many components as there
     * are space dimensions), the
     * value is a tensor.
     *
     * In this array, we store
     * the values of the shape
     * function in the quadrature
     * points on the unit
     * cell. The transformation
     * to the real space cell is
     * then simply done by
     * multiplication with the
     * Jacobian of the mapping.
     */
    std::vector<std::vector<Tensor<1,dim> > > shape_values;

    /**
     * Array with shape function
     * gradients in quadrature
     * points. There is one
     * row for each shape
     * function, containing
     * values for each quadrature
     * point.
     *
     * We store the gradients in
     * the quadrature points on
     * the unit cell. We then
     * only have to apply the
     * transformation (which is a
     * matrix-vector
     * multiplication) when
     * visiting an actual cell.
     */
    std::vector<std::vector<Tensor<2,dim> > > shape_gradients;
  };

  /**
   * These are the factors
   * multiplied to a function in
   * the
   * #generalized_face_support_points
   * when computing the
   * integration. They are
   * organized such that there is
   * one row for each generalized
   * face support point and one
   * column for each degree of
   * freedom on the face.
   */
  Table<2, double> boundary_weights;
  /**
   * Precomputed factors for
   * interpolation of interior
   * degrees of freedom. The
   * rationale for this Table is
   * the same as for
   * #boundary_weights. Only, this
   * table has a third coordinate
   * for the space direction of the
   * component evaluated.
   */
  Table<3, double> interior_weights;



  /**
   * These are the factors
   * multiplied to a function in
   * the
   * #generalized_face_support_points
   * when computing the
   * integration. They are
   * organized such that there is
   * one row for each generalized
   * face support point and one
   * column for each degree of
   * freedom on the face.
   */
  Table<2, double> boundary_weights_abf;
  /**
   * Precomputed factors for
   * interpolation of interior
   * degrees of freedom. The
   * rationale for this Table is
   * the same as for
   * #boundary_weights. Only, this
   * table has a third coordinate
   * for the space direction of the
   * component evaluated.
   */
  Table<3, double> interior_weights_abf;


  /**
   * Allow access from other
   * dimensions.
   */
  template <int dim1> friend class FE_ABF;
};



/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif

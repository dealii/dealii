// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2015 by the deal.II authors
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

#ifndef dealii__fe_dgp_monomial_h
#define dealii__fe_dgp_monomial_h

#include <deal.II/base/config.h>
#include <deal.II/base/polynomials_p.h>
#include <deal.II/fe/fe_poly.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup fe */
/*@{*/

/**
 * Discontinuous finite elements based on monomials.
 *
 * This finite element implements complete polynomial spaces, that is, dim-
 * dimensional polynomials of degree p. For example, in 2d the element
 * FE_DGP(1) would represent the span of the functions $\{1,\hat x,\hat y\}$,
 * which is in contrast to the element FE_DGQ(1) that is formed by the span of
 * $\{1,\hat x,\hat y,\hat x\hat y\}$. Since the DGP space has only three
 * unknowns for each quadrilateral, it is immediately clear that this element
 * can not be continuous.
 *
 * The basis functions for this element are chosen to be the monomials listed
 * above. Note that this is the main difference to the FE_DGP class that uses
 * a set of polynomials of complete degree <code>p</code> that form a Legendre
 * basis on the unit square. Thus, there, the mass matrix is diagonal, if the
 * grid cells are parallelograms. The basis here does not have this property;
 * however, it is simpler to compute. On the other hand, this element has the
 * additional disadvantage that the local cell matrices usually have a worse
 * condition number than the ones originating from the FE_DGP element.
 *
 * This class is not implemented for the codimension one case (<tt>spacedim !=
 * dim</tt>).
 *
 * <h3>Transformation properties</h3>
 *
 * It is worth noting that under a (bi-, tri-)linear mapping, the space
 * described by this element does not contain $P(k)$, even if we use a basis
 * of polynomials of degree $k$. Consequently, for example, on meshes with
 * non-affine cells, a linear function can not be exactly represented by
 * elements of type FE_DGP(1) or FE_DGPMonomial(1).
 *
 * This can be understood by the following 2-d example: consider the cell with
 * vertices at $(0,0),(1,0),(0,1),(s,s)$:
 * @image html dgp_doesnt_contain_p.png
 *
 * For this cell, a bilinear transformation $F$ produces the relations $x=\hat
 * x+\hat x\hat y$ and $y=\hat y+\hat x\hat y$ that correlate reference
 * coordinates $\hat x,\hat y$ and coordinates in real space $x,y$. Under this
 * mapping, the constant function is clearly mapped onto itself, but the two
 * other shape functions of the $P_1$ space, namely $\phi_1(\hat x,\hat
 * y)=\hat x$ and $\phi_2(\hat x,\hat y)=\hat y$ are mapped onto
 * $\phi_1(x,y)=\frac{x-t}{t(s-1)},\phi_2(x,y)=t$ where
 * $t=\frac{y}{s-x+sx+y-sy}$.
 *
 * For the simple case that $s=1$, i.e. if the real cell is the unit square,
 * the expressions can be simplified to $t=y$ and
 * $\phi_1(x,y)=x,\phi_2(x,y)=y$. However, for all other cases, the functions
 * $\phi_1(x,y),\phi_2(x,y)$ are not linear any more, and neither is any
 * linear combination of them. Consequently, the linear functions are not
 * within the range of the mapped $P_1$ polynomials.
 *
 *
 * <h3>Visualization of shape functions</h3> In 2d, the shape functions of
 * this element look as follows.
 *
 * <h4>$P_0$ element</h4>
 *
 * <table> <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P1/P1_DGPMonomial_shape0000.png
 * </td>
 *
 * <td align="center"> </td> </tr> <tr> <td align="center"> $P_0$ element,
 * shape function 0 </td>
 *
 * <td align="center"></tr> </table>
 *
 * <h4>$P_1$ element</h4>
 *
 * <table> <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P1/P1_DGPMonomial_shape0000.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P1/P1_DGPMonomial_shape0001.png
 * </td> </tr> <tr> <td align="center"> $P_1$ element, shape function 0 </td>
 *
 * <td align="center"> $P_1$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P1/P1_DGPMonomial_shape0002.png
 * </td>
 *
 * <td align="center"> </td> </tr> <tr> <td align="center"> $P_1$ element,
 * shape function 2 </td>
 *
 * <td align="center"></td> </tr> </table>
 *
 *
 * <h4>$P_2$ element</h4>
 *
 * <table> <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P2/P2_DGPMonomial_shape0000.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P2/P2_DGPMonomial_shape0001.png
 * </td> </tr> <tr> <td align="center"> $P_2$ element, shape function 0 </td>
 *
 * <td align="center"> $P_2$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P2/P2_DGPMonomial_shape0002.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P2/P2_DGPMonomial_shape0003.png
 * </td> </tr> <tr> <td align="center"> $P_2$ element, shape function 2 </td>
 *
 * <td align="center"> $P_2$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P2/P2_DGPMonomial_shape0004.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P2/P2_DGPMonomial_shape0005.png
 * </td> </tr> <tr> <td align="center"> $P_2$ element, shape function 4 </td>
 *
 * <td align="center"> $P_2$ element, shape function 5 </td> </tr> </table>
 *
 *
 * <h4>$P_3$ element</h4>
 *
 * <table> <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P3/P3_DGPMonomial_shape0000.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P3/P3_DGPMonomial_shape0001.png
 * </td> </tr> <tr> <td align="center"> $P_3$ element, shape function 0 </td>
 *
 * <td align="center"> $P_3$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P3/P3_DGPMonomial_shape0002.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P3/P3_DGPMonomial_shape0003.png
 * </td> </tr> <tr> <td align="center"> $P_3$ element, shape function 2 </td>
 *
 * <td align="center"> $P_3$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P3/P3_DGPMonomial_shape0004.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P3/P3_DGPMonomial_shape0005.png
 * </td> </tr> <tr> <td align="center"> $P_3$ element, shape function 4 </td>
 *
 * <td align="center"> $P_3$ element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P3/P3_DGPMonomial_shape0006.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P3/P3_DGPMonomial_shape0007.png
 * </td> </tr> <tr> <td align="center"> $P_3$ element, shape function 6 </td>
 *
 * <td align="center"> $P_3$ element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P3/P3_DGPMonomial_shape0008.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P3/P3_DGPMonomial_shape0009.png
 * </td> </tr> <tr> <td align="center"> $P_3$ element, shape function 8 </td>
 *
 * <td align="center"> $P_3$ element, shape function 9 </td> </tr> </table>
 *
 *
 * <h4>$P_4$ element</h4> <table> <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0000.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0001.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 0 </td>
 *
 * <td align="center"> $P_4$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0002.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0003.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 2 </td>
 *
 * <td align="center"> $P_4$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0004.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0005.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 4 </td>
 *
 * <td align="center"> $P_4$ element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0006.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0007.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 6 </td>
 *
 * <td align="center"> $P_4$ element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0008.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0009.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 8 </td>
 *
 * <td align="center"> $P_4$ element, shape function 9 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0010.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0011.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 10 </td>
 *
 * <td align="center"> $P_4$ element, shape function 11 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0012.png
 * </td>
 *
 * <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0013.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 12 </td>
 *
 * <td align="center"> $P_4$ element, shape function 13 </td> </tr>
 *
 * <tr> <td align="center">
 * @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0014.png
 * </td>
 *
 * <td align="center"> </td> </tr> <tr> <td align="center"> $P_4$ element,
 * shape function 14 </td>
 *
 * <td align="center"></td> </tr> </table>
 *
 * @author Ralf Hartmann, 2004
 */
template <int dim>
class FE_DGPMonomial : public FE_Poly<PolynomialsP<dim>,dim>
{
public:
  /**
   * Constructor for the polynomial space of degree <tt>p</tt>.
   */
  FE_DGPMonomial (const unsigned int p);

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_DGPMonomial<dim>(degree)</tt>, with <tt>dim</tt> and
   * <tt>p</tt> replaced by appropriate values.
   */
  virtual std::string get_name () const;

  /**
   * @name Functions to support hp
   * @{
   */

  /**
   * If, on a vertex, several finite elements are active, the hp code first
   * assigns the degrees of freedom of each of these FEs different global
   * indices. It then calls this function to find out which of them should get
   * identical values, and consequently can receive the same global DoF index.
   * This function therefore returns a list of identities between DoFs of the
   * present finite element object with the DoFs of @p fe_other, which is a
   * reference to a finite element object representing one of the other finite
   * elements active on this particular vertex. The function computes which of
   * the degrees of freedom of the two finite element objects are equivalent,
   * both numbered between zero and the corresponding value of dofs_per_vertex
   * of the two finite elements. The first index of each pair denotes one of
   * the vertex dofs of the present element, whereas the second is the
   * corresponding index of the other finite element.
   *
   * This being a discontinuous element, the set of such constraints is of
   * course empty.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_vertex_dof_identities (const FiniteElement<dim> &fe_other) const;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on lines.
   *
   * This being a discontinuous element, the set of such constraints is of
   * course empty.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_line_dof_identities (const FiniteElement<dim> &fe_other) const;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on quads.
   *
   * This being a discontinuous element, the set of such constraints is of
   * course empty.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_quad_dof_identities (const FiniteElement<dim> &fe_other) const;

  /**
   * Return whether this element implements its hanging node constraints in
   * the new way, which has to be used to make elements "hp compatible".
   *
   * For the FE_DGPMonomial class the result is always true (independent of
   * the degree of the element), as it has no hanging nodes (being a
   * discontinuous element).
   */
  virtual bool hp_constraints_are_implemented () const;

  /**
   * Return whether this element dominates the one given as argument when they
   * meet at a common face, whether it is the other way around, whether
   * neither dominates, or if either could dominate.
   *
   * For a definition of domination, see FiniteElementBase::Domination and in
   * particular the
   * @ref hp_paper "hp paper".
   */
  virtual
  FiniteElementDomination::Domination
  compare_for_face_domination (const FiniteElement<dim> &fe_other) const;

  /**
   * @}
   */

  /**
   * Return the matrix interpolating from the given finite element to the
   * present one. The size of the matrix is then @p dofs_per_cell times
   * <tt>source.dofs_per_cell</tt>.
   *
   * These matrices are only available if the source element is also a @p FE_Q
   * element. Otherwise, an exception of type
   * FiniteElement<dim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_interpolation_matrix (const FiniteElement<dim> &source,
                            FullMatrix<double>           &matrix) const;

  /**
   * Return the matrix interpolating from a face of of one element to the face
   * of the neighboring element. The size of the matrix is then @p
   * dofs_per_face times <tt>source.dofs_per_face</tt>.
   *
   * Derived elements will have to implement this function. They may only
   * provide interpolation matrices for certain source finite elements, for
   * example those from the same family. If they don't implement interpolation
   * from a given element, then they must throw an exception of type
   * FiniteElement<dim>::ExcInterpolationNotImplemented.
   */
  virtual void
  get_face_interpolation_matrix (const FiniteElement<dim> &source,
                                 FullMatrix<double>       &matrix) const;

  /**
   * Return the matrix interpolating from a face of of one element to the face
   * of the neighboring element. The size of the matrix is then @p
   * dofs_per_face times <tt>source.dofs_per_face</tt>.
   *
   * Derived elements will have to implement this function. They may only
   * provide interpolation matrices for certain source finite elements, for
   * example those from the same family. If they don't implement interpolation
   * from a given element, then they must throw an exception of type
   * FiniteElement<dim>::ExcInterpolationNotImplemented.
   */
  virtual void
  get_subface_interpolation_matrix (const FiniteElement<dim> &source,
                                    const unsigned int        subface,
                                    FullMatrix<double>       &matrix) const;

  /**
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero function values somewhere on the face @p face_index.
   */
  virtual bool has_support_on_face (const unsigned int shape_index,
                                    const unsigned int face_index) const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   *
   * This function is made virtual, since finite element objects are usually
   * accessed through pointers to their base class, rather than the class
   * itself.
   */
  virtual std::size_t memory_consumption () const;

protected:

  /**
   * @p clone function instead of a copy constructor.
   *
   * This function is needed by the constructors of @p FESystem.
   */
  virtual FiniteElement<dim> *clone() const;

private:

  /**
   * Only for internal use. Its full name is @p get_dofs_per_object_vector
   * function and it creates the @p dofs_per_object vector that is needed
   * within the constructor to be passed to the constructor of @p
   * FiniteElementData.
   */
  static std::vector<unsigned int> get_dpo_vector (const unsigned int degree);

  /**
   * Initialize the embedding matrices. Called from the constructor.
   */
  void initialize_embedding ();

  /**
   * Initialize the restriction matrices. Called from the constructor.
   */
  void initialize_restriction ();
};

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif

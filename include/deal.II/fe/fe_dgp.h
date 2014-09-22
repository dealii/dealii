// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2014 by the deal.II authors
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

#ifndef __deal2__fe_dgp_h
#define __deal2__fe_dgp_h

#include <deal.II/base/config.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/fe/fe_poly.h>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim> class MappingQ;

/*!@addtogroup fe */
/*@{*/

/**
 * Discontinuous finite elements based on Legendre polynomials.
 *
 * This finite element implements complete polynomial spaces, that is,
 * dim-dimensional polynomials of degree p. For example, in 2d the
 * element FE_DGP(1) would represent the span of the functions
 * $\{1,\hat x,\hat y\}$, which is in contrast to the element FE_DGQ(1)
 * that is formed by the span of $\{1,\hat x,\hat y,\hat x\hat y\}$. Since the
 * DGP space has only three unknowns for each quadrilateral, it is
 * immediately clear that this element can not be continuous.
 *
 * The basis functions used in this element for the space described above
 * are chosen to form a Legendre basis on the unit square. As a consequence,
 * the first basis function of this element is always the function that
 * is constant and equal to one.  As a result of the orthogonality of
 * the basis functions, the mass matrix is diagonal if the
 * grid cells are parallelograms. Note that this is in contrast to the
 * FE_DGPMonomial class that actually uses the monomial basis listed
 * above as basis functions.
 *
 * The shape functions are defined in the class PolynomialSpace. The
 * polynomials used inside PolynomialSpace are Polynomials::Legendre
 * up to degree <tt>p</tt> given in FE_DGP. For the ordering of the
 * basis functions, refer to PolynomialSpace, remembering that the
 * Legendre polynomials are ordered by ascending degree.
 *
 * @note This element is not defined by finding shape functions within
 * the given function space that interpolate a particular set of points.
 * Consequently, there are no support points to which a given function
 * could be interpolated; finding a finite element function that approximates
 * a given function is therefore only possible through projection, rather
 * than interpolation. Secondly, the shape functions of this element do not
 * jointly add up to one. As a consequence of this, adding or subtracting
 * a constant value -- such as one would do to make a function have mean
 * value zero -- can not be done by simply subtracting the constant value
 * from each degree of freedom. Rather, one needs to use the fact that the
 * first basis function is constant equal to one and simply subtract the
 * constant from the value of the degree of freedom corresponding to this
 * first shape function on each cell.
 *
 *
 * @note This class is only partially implemented for the codimension one case
 * (<tt>spacedim != dim </tt>), since no passage of information
 * between meshes of different refinement level is possible because
 * the embedding and projection matrices are not computed in the class
 * constructor.
 *
 * <h3>Transformation properties</h3>
 *
 * It is worth noting that under a (bi-, tri-)linear mapping, the
 * space described by this element does not contain $P(k)$, even if we
 * use a basis of polynomials of degree $k$. Consequently, for
 * example, on meshes with non-affine cells, a linear function can not
 * be exactly represented by elements of type FE_DGP(1) or
 * FE_DGPMonomial(1).
 *
 * This can be understood by the following 2-d example: consider the
 * cell with vertices at $(0,0),(1,0),(0,1),(s,s)$:
 * @image html dgp_doesnt_contain_p.png
 *
 * For this cell, a bilinear transformation $F$ produces the relations
 * $x=\hat x+\hat x\hat y$ and $y=\hat y+\hat x\hat y$ that correlate
 * reference coordinates $\hat x,\hat y$ and coordinates in real space
 * $x,y$. Under this mapping, the constant function is clearly mapped
 * onto itself, but the two other shape functions of the $P_1$ space,
 * namely $\phi_1(\hat x,\hat y)=\hat x$ and $\phi_2(\hat x,\hat
 * y)=\hat y$ are mapped onto
 * $\phi_1(x,y)=\frac{x-t}{t(s-1)},\phi_2(x,y)=t$ where
 * $t=\frac{y}{s-x+sx+y-sy}$.
 *
 * For the simple case that $s=1$, i.e. if the real cell is the unit
 * square, the expressions can be simplified to $t=y$ and
 * $\phi_1(x,y)=x,\phi_2(x,y)=y$. However, for all other cases, the
 * functions $\phi_1(x,y),\phi_2(x,y)$ are not linear any more, and
 * neither is any linear combincation of them. Consequently, the
 * linear functions are not within the range of the mapped $P_1$
 * polynomials.
 *
 *
 * @author Guido Kanschat, 2001, 2002, Ralf Hartmann 2004
 */
template <int dim, int spacedim=dim>
class FE_DGP : public FE_Poly<PolynomialSpace<dim>,dim,spacedim>
{
public:
  /**
   * Constructor for tensor product
   * polynomials of degree @p p.
   */
  FE_DGP (const unsigned int p);

  /**
   * Return a string that uniquely
   * identifies a finite
   * element. This class returns
   * <tt>FE_DGP<dim>(degree)</tt>, with
   * @p dim and @p degree
   * replaced by appropriate
   * values.
   */
  virtual std::string get_name () const;

  /**
   * @name Functions to support hp
   * @{
   */

  /**
   * If, on a vertex, several finite elements are active, the hp code
   * first assigns the degrees of freedom of each of these FEs
   * different global indices. It then calls this function to find out
   * which of them should get identical values, and consequently can
   * receive the same global DoF index. This function therefore
   * returns a list of identities between DoFs of the present finite
   * element object with the DoFs of @p fe_other, which is a reference
   * to a finite element object representing one of the other finite
   * elements active on this particular vertex. The function computes
   * which of the degrees of freedom of the two finite element objects
   * are equivalent, both numbered between zero and the corresponding
   * value of dofs_per_vertex of the two finite elements. The first
   * index of each pair denotes one of the vertex dofs of the present
   * element, whereas the second is the corresponding index of the
   * other finite element.
   *
   * This being a discontinuous element, the set of such constraints
   * is of course empty.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_vertex_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats
   * degrees of freedom on lines.
   *
   * This being a discontinuous element, the set of such constraints
   * is of course empty.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_line_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const;

  /**
   * Same as
   * hp_vertex_dof_indices(),
   * except that the function
   * treats degrees of freedom on
   * quads.
   *
   * This being a discontinuous element,
   * the set of such constraints is of
   * course empty.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_quad_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const;

  /**
   * Return whether this element
   * implements its hanging node
   * constraints in the new way,
   * which has to be used to make
   * elements "hp compatible".
   *
   * For the FE_DGP class the result is
   * always true (independent of the degree
   * of the element), as it has no hanging
   * nodes (being a discontinuous element).
   */
  virtual bool hp_constraints_are_implemented () const;

  /**
   * Return whether this element dominates
   * the one given as argument when they
   * meet at a common face,
   * whether it is the other way around,
   * whether neither dominates, or if
   * either could dominate.
   *
   * For a definition of domination, see
   * FiniteElementBase::Domination and in
   * particular the @ref hp_paper "hp paper".
   */
  virtual
  FiniteElementDomination::Domination
  compare_for_face_domination (const FiniteElement<dim,spacedim> &fe_other) const;

  /**
   * @}
   */

  /**
   * Return the matrix
   * interpolating from a face of
   * of one element to the face of
   * the neighboring element.
   * The size of the matrix is
   * then <tt>source.dofs_per_face</tt> times
   * <tt>this->dofs_per_face</tt>.
   *
   * Derived elements will have to
   * implement this function. They
   * may only provide interpolation
   * matrices for certain source
   * finite elements, for example
   * those from the same family. If
   * they don't implement
   * interpolation from a given
   * element, then they must throw
   * an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented.
   */
  virtual void
  get_face_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
                                 FullMatrix<double>       &matrix) const;

  /**
   * Return the matrix
   * interpolating from a face of
   * of one element to the face of
   * the neighboring element.
   * The size of the matrix is
   * then <tt>source.dofs_per_face</tt> times
   * <tt>this->dofs_per_face</tt>.
   *
   * Derived elements will have to
   * implement this function. They
   * may only provide interpolation
   * matrices for certain source
   * finite elements, for example
   * those from the same family. If
   * they don't implement
   * interpolation from a given
   * element, then they must throw
   * an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented.
   */
  virtual void
  get_subface_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
                                    const unsigned int        subface,
                                    FullMatrix<double>       &matrix) const;

  /**
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero function values somewhere on the face @p face_index.
   */
  virtual bool has_support_on_face (const unsigned int shape_index,
                                    const unsigned int face_index) const;

  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   *
   * This function is made virtual,
   * since finite element objects
   * are usually accessed through
   * pointers to their base class,
   * rather than the class itself.
   */
  virtual std::size_t memory_consumption () const;


  /**
   * Declare a nested class which
   * will hold static definitions of
   * various matrices such as
   * constraint and embedding
   * matrices. The definition of
   * the various static fields are
   * in the files <tt>fe_dgp_[123]d.cc</tt>
   * in the source directory.
   */
  struct Matrices
  {
    /**
     * As @p embedding but for
     * projection matrices.
     */
    static const double *const projection_matrices[][GeometryInfo<dim>::max_children_per_cell];

    /**
     * As
     * @p n_embedding_matrices
     * but for projection
     * matrices.
     */
    static const unsigned int n_projection_matrices;
  };

  /**
   * Returns a list of constant modes of the element. For this element, the
   * first entry is true, all other are false.
   */
  virtual std::pair<Table<2,bool>, std::vector<unsigned int> >
  get_constant_modes () const;

protected:

  /**
   * @p clone function instead of
   * a copy constructor.
   *
   * This function is needed by the
   * constructors of @p FESystem.
   */
  virtual FiniteElement<dim,spacedim> *clone() const;

private:

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
  static std::vector<unsigned int> get_dpo_vector (const unsigned int degree);
};

/* @} */
#ifndef DOXYGEN


// declaration of explicit specializations of member variables, if the
// compiler allows us to do that (the standard says we must)
#ifndef DEAL_II_MEMBER_VAR_SPECIALIZATION_BUG
template <>
const double *const FE_DGP<1>::Matrices::projection_matrices[][GeometryInfo<1>::max_children_per_cell];

template <>
const unsigned int FE_DGP<1>::Matrices::n_projection_matrices;

template <>
const double *const FE_DGP<2>::Matrices::projection_matrices[][GeometryInfo<2>::max_children_per_cell];

template <>
const unsigned int FE_DGP<2>::Matrices::n_projection_matrices;

template <>
const double *const FE_DGP<3>::Matrices::projection_matrices[][GeometryInfo<3>::max_children_per_cell];

template <>
const unsigned int FE_DGP<3>::Matrices::n_projection_matrices;

//codimension 1
template <>
const double *const FE_DGP<1,2>::Matrices::projection_matrices[][GeometryInfo<1>::max_children_per_cell];

template <>
const unsigned int FE_DGP<1,2>::Matrices::n_projection_matrices;

template <>
const double *const FE_DGP<2,3>::Matrices::projection_matrices[][GeometryInfo<2>::max_children_per_cell];

template <>
const unsigned int FE_DGP<2,3>::Matrices::n_projection_matrices;

#endif

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

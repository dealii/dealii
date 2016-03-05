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

#ifndef dealii__fe_base_h
#define dealii__fe_base_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/table.h>
#include <deal.II/base/vector_slice.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_indices.h>
#include <deal.II/fe/fe_update_flags.h>

#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * A namespace solely for the purpose of defining the Domination enum as well
 * as associated operators.
 */
namespace FiniteElementDomination
{
  /**
   * An enum that describes the outcome of comparing two elements for mutual
   * domination. If one element dominates another, then the restriction of the
   * space described by the dominated element to a face of the cell is
   * strictly larger than that of the dominating element. For example, in 2-d
   * Q(2) elements dominate Q(4) elements, because the traces of Q(4) elements
   * are quartic polynomials which is a space strictly larger than the
   * quadratic polynomials (the restriction of the Q(2) element). In general,
   * Q(k) dominates Q(k') if $k\le k'$.
   *
   * This enum is used in the FiniteElement::compare_for_face_domination()
   * function that is used in the context of hp finite element methods when
   * determining what to do at faces where two different finite elements meet
   * (see the
   * @ref hp_paper "hp paper"
   * for a more detailed description of the following). In that case, the
   * degrees of freedom of one side need to be constrained to those on the
   * other side. The determination which side is which is based on the outcome
   * of a comparison for mutual domination: the dominated side is constrained
   * to the dominating one.
   *
   * A similar situation happens in 3d, where we have to consider different
   * elements meeting at only an edge, not an entire face. Such comparisons
   * are then implemented in the FiniteElement::compare_for_line_domination()
   * function.
   *
   * Note that there are situations where neither side dominates. The
   * @ref hp_paper "hp paper"
   * lists two case, with the simpler one being that a $Q_2\times Q_1$ vector-
   * valued element (i.e. a <code>FESystem(FE_Q(2),1,FE_Q(1),1)</code>) meets
   * a $Q_1\times Q_2$ element: here, for each of the two vector-components,
   * we can define a domination relationship, but it is different for the two
   * components.
   *
   * It is clear that the concept of domination doesn't matter for
   * discontinuous elements. However, discontinuous elements may be part of
   * vector-valued elements and may therefore be compared against each other
   * for domination. They should return
   * <code>either_element_can_dominate</code> in that case. Likewise, when
   * comparing two identical finite elements, they should return this code;
   * the reason is that we can not decide which element will dominate at the
   * time we look at the first component of, for example, two $Q_2\times Q_1$
   * and $Q_2\times Q_2$ elements, and have to keep our options open until we
   * get to the second base element.
   *
   * Finally, the code no_requirements exists for cases where elements impose
   * no continuity requirements. The case is primarily meant for FE_Nothing
   * which is an element that has no degrees of freedom in a subdomain. It
   * could also be used by discontinuous elements, for example.
   *
   * More details on domination can be found in the
   * @ref hp_paper "hp paper".
   */
  enum Domination
  {
    this_element_dominates,
    other_element_dominates,
    neither_element_dominates,
    either_element_can_dominate,
    no_requirements
  };


  /**
   * A generalization of the binary <code>and</code> operator to a comparison
   * relationship. The way this works is pretty much as when you would want to
   * define a comparison relationship for vectors: either all elements of the
   * first vector are smaller, equal, or larger than those of the second
   * vector, or some are and some are not.
   *
   * This operator is pretty much the same: if both arguments are
   * <code>this_element_dominates</code> or
   * <code>other_element_dominates</code>, then the returned value is that
   * value. On the other hand, if one of the values is
   * <code>either_element_can_dominate</code>, then the returned value is that
   * of the other argument. If either argument is
   * <code>neither_element_dominates</code>, or if the two arguments are
   * <code>this_element_dominates</code> and
   * <code>other_element_dominates</code>, then the returned value is
   * <code>neither_element_dominates</code>.
   */
  inline Domination operator & (const Domination d1,
                                const Domination d2);
}


/**
 * A class that declares a number of scalar constant variables that describe
 * basic properties of a finite element implementation. This includes, for
 * example, the number of degrees of freedom per vertex, line, or cell; the
 * number of vector components; etc.
 *
 * The kind of information stored here is computed during initialization of a
 * finite element object and is passed down to this class via its constructor.
 * The data stored by this class is part of the public interface of the
 * FiniteElement class (which derives from the current class). See there for
 * more information.
 *
 * @ingroup febase
 * @author Wolfgang Bangerth, Guido Kanschat, 1998, 1999, 2000, 2001, 2003,
 * 2005
 */
template <int dim>
class FiniteElementData
{
public:
  /**
   * Enumerator for the different types of continuity a finite element may
   * have. Continuity is measured by the Sobolev space containing the
   * constructed finite element space and is also called this way.
   *
   * Note that certain continuities may imply others. For instance, a function
   * in <i>H<sup>1</sup></i> is in <i>H<sup>curl</sup></i> and
   * <i>H<sup>div</sup></i> as well.
   *
   * If you are interested in continuity in the classical sense, then the
   * following relations hold:
   *
   * <ol>
   *
   * <li> <i>H<sup>1</sup></i> implies that the function is continuous over
   * cell boundaries.
   *
   * <li> <i>H<sup>2</sup></i> implies that the function is continuously
   * differentiable over cell boundaries.
   *
   * <li> <i>L<sup>2</sup></i> indicates that the element is discontinuous.
   * Since discontinuous elements have no topological couplings between grid
   * cells and code may actually depend on this property, <i>L<sup>2</sup></i>
   * conformity is handled in a special way in the sense that it is <b>not</b>
   * implied by any higher conformity.
   * </ol>
   *
   * In order to test if a finite element conforms to a certain space, use
   * FiniteElementData<dim>::conforms().
   */
  enum Conformity
  {
    /**
     * Indicates incompatible continuities of a system.
     */
    unknown = 0x00,

    /**
     * Discontinuous elements. See above!
     */
    L2 = 0x01,

    /**
     * Conformity with the space <i>H<sup>curl</sup></i> (continuous
     * tangential component of a vector field)
     */
    Hcurl = 0x02,

    /**
     * Conformity with the space <i>H<sup>div</sup></i> (continuous normal
     * component of a vector field)
     */
    Hdiv = 0x04,

    /**
     * Conformity with the space <i>H<sup>1</sup></i> (continuous)
     */
    H1 = Hcurl | Hdiv,

    /**
     * Conformity with the space <i>H<sup>2</sup></i> (continuously
     * differentiable)
     */
    H2 = 0x0e
  };

  /**
   * The dimension of the finite element, which is the template parameter
   * <tt>dim</tt>
   */
  static const unsigned int dimension = dim;

  /**
   * Number of degrees of freedom on a vertex.
   */
  const unsigned int dofs_per_vertex;

  /**
   * Number of degrees of freedom in a line; not including the degrees of
   * freedom on the vertices of the line.
   */
  const unsigned int dofs_per_line;

  /**
   * Number of degrees of freedom in a quadrilateral; not including the
   * degrees of freedom on the lines and vertices of the quadrilateral.
   */
  const unsigned int dofs_per_quad;

  /**
   * Number of degrees of freedom in a hexahedron; not including the degrees
   * of freedom on the quadrilaterals, lines and vertices of the hexahedron.
   */
  const unsigned int dofs_per_hex;

  /**
   * First index of dof on a line.
   */
  const unsigned int first_line_index;

  /**
   * First index of dof on a quad.
   */
  const unsigned int first_quad_index;

  /**
   * First index of dof on a hexahedron.
   */
  const unsigned int first_hex_index;

  /**
   * First index of dof on a line for face data.
   */
  const unsigned int first_face_line_index;

  /**
   * First index of dof on a quad for face data.
   */
  const unsigned int first_face_quad_index;

  /**
   * Number of degrees of freedom on a face. This is the accumulated number of
   * degrees of freedom on all the objects of dimension up to <tt>dim-1</tt>
   * constituting a face.
   */
  const unsigned int dofs_per_face;

  /**
   * Total number of degrees of freedom on a cell. This is the accumulated
   * number of degrees of freedom on all the objects of dimension up to
   * <tt>dim</tt> constituting a cell.
   */
  const unsigned int dofs_per_cell;

  /**
   * Number of vector components of this finite element, and dimension of the
   * image space. For vector-valued finite elements (i.e. when this number is
   * greater than one), the number of vector components is in many cases equal
   * to the number of base elements glued together with the help of the
   * FESystem class. However, for elements like the Nedelec element, the
   * number is greater than one even though we only have one base element.
   */
  const unsigned int components;

  /**
   * Maximal polynomial degree of a shape function in a single coordinate
   * direction.
   */
  const unsigned int degree;

  /**
   * Indicate the space this element conforms to.
   */
  const Conformity conforming_space;

  /**
   * Storage for an object describing the sizes of each block of a compound
   * element. For an element which is not an FESystem, this contains only a
   * single block with length #dofs_per_cell.
   */
  const BlockIndices block_indices_data;

  /**
   * Constructor, computing all necessary values from the distribution of dofs
   * to geometrical objects.
   *
   * @param[in] dofs_per_object A vector that describes the number of degrees
   * of freedom on geometrical objects for each dimension. This vector must
   * have size dim+1, and entry 0 describes the number of degrees of freedom
   * per vertex, entry 1 the number of degrees of freedom per line, etc. As an
   * example, for the common $Q_1$ Lagrange element in 2d, this vector would
   * have elements <code>(1,0,0)</code>. On the other hand, for a $Q_3$
   * element in 3d, it would have entries <code>(1,2,4,8)</code>.
   *
   * @param[in] n_components Number of vector components of the element.
   *
   * @param[in] degree The maximal polynomial degree of any of the shape
   * functions of this element in any variable on the reference element. For
   * example, for the $Q_1$ element (in any space dimension), this would be
   * one; this is so despite the fact that the element has a shape function of
   * the form $\hat x\hat y$ (in 2d) and $\hat x\hat y\hat z$ (in 3d), which,
   * although quadratic and cubic polynomials, are still only linear in each
   * reference variable separately. The information provided by this variable
   * is typically used in determining what an appropriate quadrature formula
   * is.
   *
   * @param[in] conformity A variable describing which Sobolev space this
   * element conforms to. For example, the $Q_p$ Lagrange elements
   * (implemented by the FE_Q class) are $H^1$ conforming, whereas the
   * Raviart-Thomas element (implemented by the FE_RaviartThomas class) is
   * $H_\text{div}$ conforming; finally, completely discontinuous elements
   * (implemented by the FE_DGQ class) are only $L_2$ conforming.
   *
   * @param[in] block_indices An argument that describes how the base elements
   * of a finite element are grouped. The default value constructs a single
   * block that consists of all @p dofs_per_cell degrees of freedom. This is
   * appropriate for all "atomic" elements (including non-primitive ones) and
   * these can therefore omit this argument. On the other hand, composed
   * elements such as FESystem will want to pass a different value here.
   */
  FiniteElementData (const std::vector<unsigned int> &dofs_per_object,
                     const unsigned int               n_components,
                     const unsigned int               degree,
                     const Conformity                 conformity = unknown,
                     const BlockIndices              &block_indices = BlockIndices());

  /**
   * Number of dofs per vertex.
   */
  unsigned int n_dofs_per_vertex () const;

  /**
   * Number of dofs per line. Not including dofs on lower dimensional objects.
   */
  unsigned int n_dofs_per_line () const;

  /**
   * Number of dofs per quad. Not including dofs on lower dimensional objects.
   */
  unsigned int n_dofs_per_quad () const;

  /**
   * Number of dofs per hex. Not including dofs on lower dimensional objects.
   */
  unsigned int n_dofs_per_hex () const;

  /**
   * Number of dofs per face, accumulating degrees of freedom of all lower
   * dimensional objects.
   */
  unsigned int n_dofs_per_face () const;

  /**
   * Number of dofs per cell, accumulating degrees of freedom of all lower
   * dimensional objects.
   */
  unsigned int n_dofs_per_cell () const;

  /**
   * Return the number of degrees per structdim-dimensional object. For
   * structdim==0, the function therefore returns dofs_per_vertex, for
   * structdim==1 dofs_per_line, etc. This function is mostly used to allow
   * some template trickery for functions that should work on all sorts of
   * objects without wanting to use the different names (vertex, line, ...)
   * associated with these objects.
   */
  template <int structdim>
  unsigned int n_dofs_per_object () const;

  /**
   * Number of components. See
   * @ref GlossComponent "the glossary"
   * for more information.
   */
  unsigned int n_components () const;

  /**
   * Number of blocks. See
   * @ref GlossBlock "the glossary"
   * for more information.
   */
  unsigned int n_blocks () const;

  /**
   * Detailed information on block sizes.
   */
  const BlockIndices &block_indices() const;

  /**
   * Return whether the entire finite element is primitive, in the sense that
   * all its shape functions are primitive. If the finite element is scalar,
   * then this is always the case.
   *
   * Since this is an extremely common operation, the result is cached in the
   * #cached_primitivity variable which is computed in the constructor.
   */
  bool is_primitive () const;

  /**
   * Maximal polynomial degree of a shape function in a single coordinate
   * direction.
   *
   * This function can be used to determine the optimal quadrature rule.
   */
  unsigned int tensor_degree () const;

  /**
   * Test whether a finite element space conforms to a certain Sobolev space.
   *
   * @note This function will return a true value even if the finite element
   * space has higher regularity than asked for.
   */
  bool conforms (const Conformity) const;

  /**
   * Comparison operator.
   */
  bool operator == (const FiniteElementData &) const;

protected:

  /**
   * Set the primitivity of the element. This is usually done by the
   * constructor of a derived class.  See
   * @ref GlossPrimitive "primitive"
   * for details.
   */
  void set_primitivity(const bool value);

private:
  /**
   * Store whether all shape functions are primitive. Since finding this out
   * is a very common operation, we cache the result, i.e. compute the value
   * in the constructor for simpler access.
   */
  bool cached_primitivity;
};



// --------- inline and template functions ---------------


#ifndef DOXYGEN

namespace FiniteElementDomination
{
  inline
  Domination operator & (const Domination d1,
                         const Domination d2)
  {
    // go through the entire list of possibilities. note that if we were into
    // speed, obfuscation and cared enough, we could implement this operator
    // by doing a bitwise & (and) if we gave these values to the enum values:
    // neither_element_dominates=0, this_element_dominates=1,
    // other_element_dominates=2, either_element_can_dominate=3
    // =this_element_dominates|other_element_dominates
    switch (d1)
      {
      case this_element_dominates:
        if ((d2 == this_element_dominates) ||
            (d2 == either_element_can_dominate) ||
            (d2 == no_requirements))
          return this_element_dominates;
        else
          return neither_element_dominates;

      case other_element_dominates:
        if ((d2 == other_element_dominates) ||
            (d2 == either_element_can_dominate) ||
            (d2 == no_requirements))
          return other_element_dominates;
        else
          return neither_element_dominates;

      case neither_element_dominates:
        return neither_element_dominates;

      case either_element_can_dominate:
        if (d2 == no_requirements)
          return either_element_can_dominate;
        else
          return d2;

      case no_requirements:
        return d2;

      default:
        // shouldn't get here
        Assert (false, ExcInternalError());
      }

    return neither_element_dominates;
  }
}


template <int dim>
inline
unsigned int
FiniteElementData<dim>::n_dofs_per_vertex () const
{
  return dofs_per_vertex;
}



template <int dim>
inline
unsigned int
FiniteElementData<dim>::n_dofs_per_line () const
{
  return dofs_per_line;
}



template <int dim>
inline
unsigned int
FiniteElementData<dim>::n_dofs_per_quad () const
{
  return dofs_per_quad;
}



template <int dim>
inline
unsigned int
FiniteElementData<dim>::n_dofs_per_hex () const
{
  return dofs_per_hex;
}



template <int dim>
inline
unsigned int
FiniteElementData<dim>::n_dofs_per_face () const
{
  return dofs_per_face;
}



template <int dim>
inline
unsigned int
FiniteElementData<dim>::n_dofs_per_cell () const
{
  return dofs_per_cell;
}



template <int dim>
template <int structdim>
inline
unsigned int
FiniteElementData<dim>::n_dofs_per_object () const
{
  switch (structdim)
    {
    case 0:
      return dofs_per_vertex;
    case 1:
      return dofs_per_line;
    case 2:
      return dofs_per_quad;
    case 3:
      return dofs_per_hex;
    default:
      Assert (false, ExcInternalError());
    }
  return numbers::invalid_unsigned_int;
}



template <int dim>
inline
unsigned int
FiniteElementData<dim>::n_components () const
{
  return components;
}



template <int dim>
inline
bool
FiniteElementData<dim>::is_primitive () const
{
  return cached_primitivity;
}


template <int dim>
inline
void
FiniteElementData<dim>::set_primitivity (const bool value)
{
  cached_primitivity = value;
}


template <int dim>
inline
const BlockIndices &
FiniteElementData<dim>::block_indices () const
{
  return block_indices_data;
}



template <int dim>
inline
unsigned int
FiniteElementData<dim>::n_blocks () const
{
  return block_indices_data.size();
}



template <int dim>
inline
unsigned int
FiniteElementData<dim>::tensor_degree () const
{
  return degree;
}


template <int dim>
inline
bool
FiniteElementData<dim>::conforms (const Conformity space) const
{
  return ((space & conforming_space) == space);
}


#endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif

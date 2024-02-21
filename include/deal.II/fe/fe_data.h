// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_data_h
#define dealii_fe_data_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/grid/reference_cell.h>

#include <deal.II/lac/block_indices.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations:
#ifndef DOXYGEN
template <int dim>
class FiniteElementData;
#endif

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
   * quadratic polynomials (the restriction of the Q(2) element). Similar
   * reasonings apply for vertices and cells as well. In general, Q(k) dominates
   * Q(k') if $k\le k'$.
   *
   * This enum is used in the FiniteElement::compare_for_domination() function
   * that is used in the context of hp-finite element methods when determining
   * what to do at faces where two different finite elements meet (see the
   * @ref hp_paper "hp-paper"
   * for a more detailed description of the following). In that case, the
   * degrees of freedom of one side need to be constrained to those on the
   * other side. The determination which side is which is based on the outcome
   * of a comparison for mutual domination: the dominated side is constrained
   * to the dominating one.
   *
   * Note that there are situations where neither side dominates. The
   * @ref hp_paper "hp-paper"
   * lists two case, with the simpler one being that a $Q_2\times Q_1$
   * vector-valued element (i.e. a <code>FESystem(FE_Q(2),1,FE_Q(1),1)</code>)
   * meets a $Q_1\times Q_2$ element: here, for each of the two
   * vector-components, we can define a domination relationship, but it is
   * different for the two components.
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
   * @ref hp_paper "hp-paper".
   */
  enum Domination
  {
    /**
     * The current element dominates.
     */
    this_element_dominates,
    /**
     * The other element dominates.
     */
    other_element_dominates,
    /**
     * Neither element dominates.
     */
    neither_element_dominates,
    /**
     * Either element may dominate.
     */
    either_element_can_dominate,
    /**
     * There are no requirements.
     */
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
  inline Domination
  operator&(const Domination d1, const Domination d2);
} // namespace FiniteElementDomination

namespace internal
{
  /**
   * Internal data structure for setting up FiniteElementData. It stores for
   * each object the (inclusive/exclusive) number of degrees of freedoms, as
   * well as the index of its first degree of freedom within a cell and the
   * index of the first d-dimensional object within each face. Here, inclusive
   * means "the number of DoFs located on this object as well as the
   * lower-dimensional objects that bound it", and exclusive is then the
   * number not including the lower-dimensional objects.
   *
   * The information is saved as a vector of vectors. One can query the
   * inclusive number of dofs of the i-th d-dimensional object via:
   * dofs_per_object_inclusive[d][i].
   *
   * As an example, the data is shown for a quadratic wedge. Which consists of
   * 6 vertices, 9 lines, and 5 faces (two triangles and three quadrilaterals).
   * @code
   *            vertices            lines                  faces      cell
   * dpo_excl  1 1 1 1 1 1 | 1 1 1 1  1  1  1  1  1 |  0  0  1  1  1 |  0
   * dpo_incl  1 1 1 1 1 1 | 3 3 3 3  3  3  3  3  3 |  6  6  9  9  9 | 18
   * obj_index 0 1 2 3 4 5 | 6 7 8 9 10 11 12 13 14 | 15 15 15 16 17 | 18
   * @endcode
   *
   * The above table has these numbers because of the following considerations:
   *
   * - For each triangular face:
   * @code
   * dpo_excl  1  1  1 | 1  1  1 |  0
   * obj_index 0  1  2 | 3  4  5 |  6
   * @endcode
   *
   * - For each quadrilateral face:
   * @code
   * dpo_excl  1  1  1  1 | 1  1  1  1 |  1
   * obj_index 0  1  2  3 | 4  5  6  7 |  8
   * @endcode
   *
   * The index of the first d-dimensional object within each face results as:
   * @code
   *                         vertices      lines       face
   * first_obj_index_on_face 0 0 0 0 0 | 3 3 4 4 4 | 6 6 8 8 8
   * @endcode
   *
   */
  struct GenericDoFsPerObject
  {
    /**
     * Exclusive number of degrees of freedom per object.
     */
    std::vector<std::vector<unsigned int>> dofs_per_object_exclusive;

    /**
     * Inclusive number of degrees of freedom per object.
     */
    std::vector<std::vector<unsigned int>> dofs_per_object_inclusive;

    /**
     * First index of an object.
     */
    std::vector<std::vector<unsigned int>> object_index;

    /**
     * First index of an object within a face.
     */
    std::vector<std::vector<unsigned int>> first_object_index_on_face;

    /**
     * Function that fills the fields based on a provided finite element.
     */
    template <int dim>
    static GenericDoFsPerObject
    generate(const FiniteElementData<dim> &fe);
  };
} // namespace internal

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
  static constexpr unsigned int dimension = dim;

private:
  /**
   * Reference cell type.
   */
  const ReferenceCell reference_cell_kind;

  /**
   * Number of unique two-dimensional sub-objects. If all
   * two-dimensional sub-objects have the same type, the value is one;
   * else it equals the number of quads.
   */
  const unsigned int number_of_unique_2d_subobjects;

  /**
   * Number of unique faces. If all faces have the same type, the value is
   * one; else it equals the number of faces.
   */
  const unsigned int number_unique_faces;

public:
  /**
   * Number of degrees of freedom on a vertex.
   */
  const unsigned int dofs_per_vertex;

  /**
   * Number of degrees of freedom in a line; not including the degrees of
   * freedom on the vertices of the line.
   */
  const unsigned int dofs_per_line;

private:
  /**
   * Number of degrees of freedom on quads. If all quads have the same
   * number of degrees freedoms the values equal dofs_per_quad.
   */
  const std::vector<unsigned int> n_dofs_on_quad;

public:
  /**
   * Number of degrees of freedom in a quadrilateral; not including the
   * degrees of freedom on the lines and vertices of the quadrilateral.
   */
  const unsigned int dofs_per_quad;

private:
  /**
   * Maximum number of degrees of freedom on any quad.
   */
  const unsigned int dofs_per_quad_max;

public:
  /**
   * Number of degrees of freedom in a hexahedron; not including the degrees
   * of freedom on the quadrilaterals, lines and vertices of the hexahedron.
   */
  const unsigned int dofs_per_hex;

  /**
   * First index of dof on a line.
   */
  const unsigned int first_line_index;

private:
  /**
   * First index of a quad. If all quads have the same number of degrees of
   * freedom, only the first index of the first quad is stored since the
   * indices of the others can be simply recomputed.
   */
  const std::vector<unsigned int> first_index_of_quads;

public:
  /**
   * First index of dof on a quad.
   */
  const unsigned int first_quad_index;

  /**
   * First index of dof on a hexahedron.
   */
  const unsigned int first_hex_index;

private:
  /**
   * Index of the first line of all faces.
   */
  const std::vector<unsigned int> first_line_index_of_faces;

public:
  /**
   * First index of dof on a line for face data.
   */
  const unsigned int first_face_line_index;

private:
  /**
   * Index of the first quad of all faces.
   */
  const std::vector<unsigned int> first_quad_index_of_faces;

public:
  /**
   * First index of dof on a quad for face data.
   */
  const unsigned int first_face_quad_index;

private:
  /**
   * Number of degrees of freedom on faces. If all faces have the same
   * number of degrees freedoms the values equal dofs_per_quad.
   */
  const std::vector<unsigned int> n_dofs_on_face;

public:
  /**
   * Number of degrees of freedom on a face. This is the accumulated number of
   * degrees of freedom on all the objects of dimension up to <tt>dim-1</tt>
   * constituting a face.
   */
  const unsigned int dofs_per_face;

private:
  /**
   * Maximum number of degrees of freedom on any face.
   */
  const unsigned int dofs_per_face_max;

public:
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
  FiniteElementData(const std::vector<unsigned int> &dofs_per_object,
                    const unsigned int               n_components,
                    const unsigned int               degree,
                    const Conformity                 conformity = unknown,
                    const BlockIndices &block_indices = BlockIndices());

  /**
   * The same as above but with the difference that also the type of the
   * underlying geometric entity can be specified.
   */
  FiniteElementData(const std::vector<unsigned int> &dofs_per_object,
                    const ReferenceCell              reference_cell,
                    const unsigned int               n_components,
                    const unsigned int               degree,
                    const Conformity                 conformity = unknown,
                    const BlockIndices &block_indices = BlockIndices());

  /**
   * The same as above but instead of passing a vector containing the degrees
   * of freedoms per object a struct of type GenericDoFsPerObject. This allows
   * that 2d objects might have different number of degrees of freedoms, which
   * is particular useful for cells with triangles and quadrilaterals as faces.
   */
  FiniteElementData(const internal::GenericDoFsPerObject &data,
                    const ReferenceCell                   reference_cell,
                    const unsigned int                    n_components,
                    const unsigned int                    degree,
                    const Conformity                      conformity = unknown,
                    const BlockIndices &block_indices = BlockIndices());

  /**
   * Return the kind of reference cell this element is defined on: For
   * example, whether the element's reference cell is a square or
   * triangle, or similar choices in higher dimensions.
   */
  ReferenceCell
  reference_cell() const;

  /**
   * Number of unique 2d subobjects. If all two-dimensional subobjects
   * are of the same kind, the value is one; else it equals the number
   * of two-dimensional subobjects. For example, for hex reference
   * cells with the usual finite elements, each face has the same
   * geometric shape (a square) and each face has the same number and
   * kind of shape functions associated with it; as a consequence, the
   * returned value is one. On the other hand, for a wedge element,
   * the two-dimensional subobjects can be both quadrilaterals and
   * triangles, and so the returned value equals the number of faces
   * of these objects (i.e., five).
   */
  unsigned int
  n_unique_2d_subobjects() const;

  /**
   * Number of unique faces. If all faces have the same type (i.e.,
   * have the same shape and also have the same kind and number of
   * DoFs associated with them), the value is one; else it equals the
   * number of faces.
   */
  unsigned int
  n_unique_faces() const;

  /**
   * Number of dofs per vertex.
   */
  unsigned int
  n_dofs_per_vertex() const;

  /**
   * Number of dofs per line. Not including dofs on lower dimensional objects.
   */
  unsigned int
  n_dofs_per_line() const;

  /**
   * Number of dofs per quad. Not including dofs on lower dimensional objects.
   */
  unsigned int
  n_dofs_per_quad(unsigned int face_no = 0) const;

  /**
   * Maximum number of dofs per quad. Not including dofs on lower dimensional
   * objects.
   */
  unsigned int
  max_dofs_per_quad() const;

  /**
   * Number of dofs per hex. Not including dofs on lower dimensional objects.
   */
  unsigned int
  n_dofs_per_hex() const;

  /**
   * Number of dofs per face, accumulating degrees of freedom of all lower
   * dimensional objects.
   */
  unsigned int
  n_dofs_per_face(unsigned int face_no = 0, unsigned int child = 0) const;

  /**
   * Maximum number of dofs per face, accumulating degrees of freedom of all
   * lower dimensional objects.
   */
  unsigned int
  max_dofs_per_face() const;

  /**
   * Number of dofs per cell, accumulating degrees of freedom of all lower
   * dimensional objects.
   */
  unsigned int
  n_dofs_per_cell() const;

  /**
   * Return the number of degrees per structdim-dimensional object. For
   * structdim==0, the function therefore returns dofs_per_vertex, for
   * structdim==1 dofs_per_line, etc. This function is mostly used to allow
   * some template trickery for functions that should work on all sorts of
   * objects without wanting to use the different names (vertex, line, ...)
   * associated with these objects.
   */
  template <int structdim>
  unsigned int
  n_dofs_per_object(const unsigned int i = 0) const;

  /**
   * Number of components. See
   * @ref GlossComponent "the glossary"
   * for more information.
   */
  unsigned int
  n_components() const;

  /**
   * Number of blocks. See
   * @ref GlossBlock "the glossary"
   * for more information.
   */
  unsigned int
  n_blocks() const;

  /**
   * Detailed information on block sizes.
   */
  const BlockIndices &
  block_indices() const;

  /**
   * Maximal polynomial degree of a shape function in a single coordinate
   * direction.
   *
   * This function can be used to determine the optimal quadrature rule.
   */
  unsigned int
  tensor_degree() const;

  /**
   * Test whether a finite element space conforms to a certain Sobolev space.
   *
   * @note This function will return a true value even if the finite element
   * space has higher regularity than asked for.
   */
  bool
  conforms(const Conformity) const;

  /**
   * Comparison operator.
   */
  bool
  operator==(const FiniteElementData &) const;

  /**
   * Return first index of dof on a line.
   */
  unsigned int
  get_first_line_index() const;

  /**
   * Return first index of dof on a quad.
   */
  unsigned int
  get_first_quad_index(const unsigned int quad_no = 0) const;

  /**
   * Return first index of dof on a hexahedron.
   */
  unsigned int
  get_first_hex_index() const;

  /**
   * Return first index of dof on a line for face data.
   */
  unsigned int
  get_first_face_line_index(const unsigned int face_no = 0) const;

  /**
   * Return first index of dof on a quad for face data.
   */
  unsigned int
  get_first_face_quad_index(const unsigned int face_no = 0) const;
};

namespace internal
{
  /**
   * Utility function to convert "dofs per object" information
   * of a @p dim dimensional reference cell @p reference_cell.
   */
  internal::GenericDoFsPerObject
  expand(const unsigned int               dim,
         const std::vector<unsigned int> &dofs_per_object,
         const ReferenceCell              reference_cell);
} // namespace internal



// --------- inline and template functions ---------------


#ifndef DOXYGEN

namespace FiniteElementDomination
{
  inline Domination
  operator&(const Domination d1, const Domination d2)
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
              (d2 == either_element_can_dominate) || (d2 == no_requirements))
            return this_element_dominates;
          else
            return neither_element_dominates;

        case other_element_dominates:
          if ((d2 == other_element_dominates) ||
              (d2 == either_element_can_dominate) || (d2 == no_requirements))
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
          DEAL_II_ASSERT_UNREACHABLE();
      }

    return neither_element_dominates;
  }
} // namespace FiniteElementDomination


template <int dim>
inline ReferenceCell
FiniteElementData<dim>::reference_cell() const
{
  return reference_cell_kind;
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::n_unique_2d_subobjects() const
{
  return number_of_unique_2d_subobjects;
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::n_unique_faces() const
{
  return number_unique_faces;
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::n_dofs_per_vertex() const
{
  return dofs_per_vertex;
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::n_dofs_per_line() const
{
  return dofs_per_line;
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::n_dofs_per_quad(unsigned int face_no) const
{
  return n_dofs_on_quad[n_dofs_on_quad.size() == 1 ? 0 : face_no];
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::max_dofs_per_quad() const
{
  return dofs_per_quad_max;
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::n_dofs_per_hex() const
{
  return dofs_per_hex;
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::n_dofs_per_face(unsigned int face_no,
                                        unsigned int child_no) const
{
  (void)child_no;

  return n_dofs_on_face[n_dofs_on_face.size() == 1 ? 0 : face_no];
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::max_dofs_per_face() const
{
  return dofs_per_face_max;
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::n_dofs_per_cell() const
{
  return dofs_per_cell;
}



template <int dim>
template <int structdim>
inline unsigned int
FiniteElementData<dim>::n_dofs_per_object(const unsigned int i) const
{
  switch (structdim)
    {
      case 0:
        return n_dofs_per_vertex();
      case 1:
        return n_dofs_per_line();
      case 2:
        return n_dofs_per_quad((structdim == 2 && dim == 3) ? i : 0);
      case 3:
        return n_dofs_per_hex();
      default:
        DEAL_II_ASSERT_UNREACHABLE();
    }
  return numbers::invalid_unsigned_int;
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::n_components() const
{
  return components;
}



template <int dim>
inline const BlockIndices &
FiniteElementData<dim>::block_indices() const
{
  return block_indices_data;
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::n_blocks() const
{
  return block_indices_data.size();
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::tensor_degree() const
{
  return degree;
}


template <int dim>
inline bool
FiniteElementData<dim>::conforms(const Conformity space) const
{
  return ((space & conforming_space) == space);
}



template <int dim>
unsigned int
FiniteElementData<dim>::get_first_line_index() const
{
  return first_line_index;
}

template <int dim>
unsigned int
FiniteElementData<dim>::get_first_quad_index(const unsigned int quad_no) const
{
  if (first_index_of_quads.size() == 1)
    return first_index_of_quads[0] + quad_no * n_dofs_per_quad(0);
  else
    return first_index_of_quads[quad_no];
}

template <int dim>
unsigned int
FiniteElementData<dim>::get_first_hex_index() const
{
  return first_hex_index;
}

template <int dim>
unsigned int
FiniteElementData<dim>::get_first_face_line_index(
  const unsigned int face_no) const
{
  return first_line_index_of_faces[first_line_index_of_faces.size() == 1 ?
                                     0 :
                                     face_no];
}

template <int dim>
unsigned int
FiniteElementData<dim>::get_first_face_quad_index(
  const unsigned int face_no) const
{
  return first_quad_index_of_faces[first_quad_index_of_faces.size() == 1 ?
                                     0 :
                                     face_no];
}

template <int dim>
internal::GenericDoFsPerObject
internal::GenericDoFsPerObject::generate(const FiniteElementData<dim> &fe)
{
  const auto reference_cell = fe.reference_cell();

  internal::GenericDoFsPerObject result;

  result.dofs_per_object_exclusive.resize(4);
  result.dofs_per_object_inclusive.resize(4);
  result.object_index.resize(4);

  unsigned int counter = 0;

  for (const unsigned int v : reference_cell.vertex_indices())
    {
      const auto c = fe.template n_dofs_per_object<0>(v);

      result.dofs_per_object_exclusive[0].emplace_back(c);
      result.dofs_per_object_inclusive[0].emplace_back(c);
      result.object_index[0].emplace_back(counter);

      counter += c;
    }

  if (dim >= 2)
    for (const unsigned int l : reference_cell.line_indices())
      {
        const auto c = fe.template n_dofs_per_object<1>(l);

        result.dofs_per_object_exclusive[1].emplace_back(c);
        result.dofs_per_object_inclusive[1].emplace_back(
          c + 2 * fe.template n_dofs_per_object<0>());
        result.object_index[1].emplace_back(counter);

        counter += c;
      }

  if (dim == 3)
    for (const unsigned int f : reference_cell.face_indices())
      {
        const auto c = fe.template n_dofs_per_object<2>(f);

        result.dofs_per_object_exclusive[2].emplace_back(c);
        result.dofs_per_object_inclusive[2].emplace_back(fe.n_dofs_per_face(f));
        result.object_index[2].emplace_back(counter);

        counter += c;
      }

  {
    const auto c = fe.template n_dofs_per_object<dim>();

    result.dofs_per_object_exclusive[dim].emplace_back(c);
    result.dofs_per_object_inclusive[dim].emplace_back(fe.n_dofs_per_cell());
    result.object_index[dim].emplace_back(counter);

    counter += c;
  }

  for (unsigned int d = dim + 1; d <= 3; ++d)
    {
      result.dofs_per_object_exclusive[d].emplace_back(0);
      result.dofs_per_object_inclusive[d].emplace_back(0);
      result.object_index[d].emplace_back(counter);
    }

  result.first_object_index_on_face.resize(3);
  for (const unsigned int face_no : reference_cell.face_indices())
    {
      result.first_object_index_on_face[0].emplace_back(0);

      result.first_object_index_on_face[1].emplace_back(
        fe.get_first_face_line_index(face_no));

      result.first_object_index_on_face[2].emplace_back(
        fe.get_first_face_quad_index(face_no));
    }

  return result;
}


#endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif

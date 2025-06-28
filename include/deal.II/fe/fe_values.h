// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_values_h
#define dealii_fe_values_h


#include <deal.II/base/config.h>

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/std_cxx20/iota_view.h>
#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values_base.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_related_data.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/read_vector.h>

#include <algorithm>
#include <memory>
#include <type_traits>

DEAL_II_NAMESPACE_OPEN

/**
 * Finite element evaluated in quadrature points of a cell.
 *
 * This function implements the initialization routines for FEValuesBase, if
 * values in quadrature points of a cell are needed. For further documentation
 * see this class.
 *
 * @ingroup feaccess
 */
template <int dim, int spacedim = dim>
class FEValues : public FEValuesBase<dim, spacedim>
{
public:
  /**
   * Dimension of the object over which we integrate. For the present class,
   * this is equal to <code>dim</code>.
   */
  static constexpr unsigned int integral_dimension = dim;

  /**
   * Constructor. Gets cell independent data from mapping and finite element
   * objects, matching the quadrature rule and update flags.
   */
  FEValues(const Mapping<dim, spacedim>       &mapping,
           const FiniteElement<dim, spacedim> &fe,
           const Quadrature<dim>              &quadrature,
           const UpdateFlags                   update_flags);

  /**
   * Like the function above, but taking a collection of quadrature rules.
   *
   * @note We require, in contrast to FEFaceValues, that the number of quadrature
   *   rules in the collection is one.
   */
  FEValues(const Mapping<dim, spacedim>       &mapping,
           const FiniteElement<dim, spacedim> &fe,
           const hp::QCollection<dim>         &quadrature,
           const UpdateFlags                   update_flags);

  /**
   * Constructor. This constructor is equivalent to the other one except that
   * it makes the object use a $Q_1$ mapping (i.e., an object of type
   * MappingQ(1)) implicitly.
   */
  FEValues(const FiniteElement<dim, spacedim> &fe,
           const Quadrature<dim>              &quadrature,
           const UpdateFlags                   update_flags);

  /**
   * Like the function above, but taking a collection of quadrature rules.
   *
   * @note We require, in contrast to FEFaceValues, that the number of quadrature
   *   rules in the collection is one.
   */
  FEValues(const FiniteElement<dim, spacedim> &fe,
           const hp::QCollection<dim>         &quadrature,
           const UpdateFlags                   update_flags);

  /**
   * Reinitialize the gradients, Jacobi determinants, etc for the given cell
   * of type "iterator into a DoFHandler object", and the finite element
   * associated with this object. It is assumed that the finite element used
   * by the given cell is also the one used by this FEValues object.
   */
  template <bool level_dof_access>
  void
  reinit(
    const TriaIterator<DoFCellAccessor<dim, spacedim, level_dof_access>> &cell);

  /**
   * Reinitialize the gradients, Jacobi determinants, etc for the given cell
   * of type "iterator into a Triangulation object", and the given finite
   * element. Since iterators into triangulation alone only convey information
   * about the geometry of a cell, but not about degrees of freedom possibly
   * associated with this cell, you will not be able to call some functions of
   * this class if they need information about degrees of freedom. These
   * functions are, above all, the
   * <tt>get_function_value/gradients/hessians/laplacians/third_derivatives</tt>
   * functions. If you want to call these functions, you have to call the @p
   * reinit variants that take iterators into DoFHandler or other DoF handler
   * type objects.
   */
  void
  reinit(const typename Triangulation<dim, spacedim>::cell_iterator &cell);

  /**
   * Return a reference to the copy of the quadrature formula stored by this
   * object.
   */
  const Quadrature<dim> &
  get_quadrature() const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  std::size_t
  memory_consumption() const;

  /**
   * Return a reference to this very object.
   *
   * Though it seems that it is not very useful, this function is there to
   * provide capability to the hp::FEValues class, in which case it provides
   * the FEValues object for the present cell (remember that for hp-finite
   * elements, the actual FE object used may change from cell to cell, so we
   * also need different FEValues objects for different cells; once you
   * reinitialize the hp::FEValues object for a specific cell, it retrieves
   * the FEValues object for the FE on that cell and returns it through a
   * function of the same name as this one; this function here therefore only
   * provides the same interface so that one can templatize on FEValues and
   * hp::FEValues).
   */
  const FEValues<dim, spacedim> &
  get_present_fe_values() const;

private:
  /**
   * Store a copy of the quadrature formula here.
   */
  const Quadrature<dim> quadrature;

  /**
   * Do work common to the two constructors.
   */
  void
  initialize(const UpdateFlags update_flags);

  /**
   * The reinit() functions do only that part of the work that requires
   * knowledge of the type of iterator. After setting present_cell(), they
   * pass on to this function, which does the real work, and which is
   * independent of the actual type of the cell iterator.
   */
  void
  do_reinit();
};


/**
 * Extend the interface of FEValuesBase to values that only make sense when
 * evaluating something on the surface of a cell. All the data that is
 * available in the interior of cells is also available here.
 *
 * See FEValuesBase
 *
 * @ingroup feaccess
 */
template <int dim, int spacedim = dim>
class FEFaceValuesBase : public FEValuesBase<dim, spacedim>
{
public:
  /**
   * Dimension of the object over which we integrate. For the present class,
   * this is equal to <code>dim-1</code>.
   */
  static constexpr unsigned int integral_dimension = dim - 1;

  /**
   * Constructor. Call the constructor of the base class and set up the arrays
   * of this class with the right sizes.  Actually filling these arrays is a
   * duty of the derived class's constructors.
   *
   * @p n_faces_or_subfaces is the number of faces or subfaces that this
   * object is to store. The actual number depends on the derived class, for
   * FEFaceValues it is <tt>2*dim</tt>, while for the FESubfaceValues class it
   * is <tt>2*dim*(1<<(dim-1))</tt>, i.e. the number of faces times the number
   * of subfaces per face.
   */
  FEFaceValuesBase(const unsigned int                  dofs_per_cell,
                   const UpdateFlags                   update_flags,
                   const Mapping<dim, spacedim>       &mapping,
                   const FiniteElement<dim, spacedim> &fe,
                   const Quadrature<dim - 1>          &quadrature);

  /**
   * Like the function above, but taking a collection of quadrature rules. This
   * allows to assign each face a different quadrature rule. In the case that
   * the collection only contains a single face quadrature, this quadrature
   * rule is use on all faces.
   */
  FEFaceValuesBase(const unsigned int                  dofs_per_cell,
                   const UpdateFlags                   update_flags,
                   const Mapping<dim, spacedim>       &mapping,
                   const FiniteElement<dim, spacedim> &fe,
                   const hp::QCollection<dim - 1>     &quadrature);

  /**
   * Boundary form of the transformation of the cell at the <tt>q_point</tt>th
   * quadrature point.  See
   * @ref GlossBoundaryForm.
   *
   * @dealiiRequiresUpdateFlags{update_boundary_forms}
   */
  const Tensor<1, spacedim> &
  boundary_form(const unsigned int q_point) const;

  /**
   * Return the list of outward normal vectors times the Jacobian of the
   * surface mapping.
   *
   * @dealiiRequiresUpdateFlags{update_boundary_forms}
   */
  const std::vector<Tensor<1, spacedim>> &
  get_boundary_forms() const;

  /**
   * Return the number of the face selected the last time the reinit() function
   * was called.
   */
  unsigned int
  get_face_number() const;

  /**
   * Return the index of the face selected the last time the reinit() function
   * was called.
   */
  unsigned int
  get_face_index() const;

  /**
   * Return a reference to the copy of the quadrature formula stored by this
   * object.
   */
  const Quadrature<dim - 1> &
  get_quadrature() const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  std::size_t
  memory_consumption() const;

protected:
  /**
   * Number of the face selected the last time the reinit() function was
   * called.
   */
  unsigned int present_face_no;

  /**
   * Index of the face selected the last time the reinit() function was
   * called.
   */
  unsigned int present_face_index;

  /**
   * Store a copy of the quadrature formula here.
   */
  const hp::QCollection<dim - 1> quadrature;
};



/**
 * Finite element evaluated in quadrature points on a face.
 *
 * This class adds the functionality of FEFaceValuesBase to FEValues; see
 * there for more documentation.
 *
 * Since finite element functions and their derivatives may be discontinuous
 * at cell boundaries, there is no restriction of this function to a mesh
 * face. But, there are limits of these values approaching the face from
 * either of the neighboring cells.
 *
 * @ingroup feaccess
 */
template <int dim, int spacedim = dim>
class FEFaceValues : public FEFaceValuesBase<dim, spacedim>
{
public:
  /**
   * Dimension in which this object operates.
   */

  static constexpr unsigned int dimension = dim;

  static constexpr unsigned int space_dimension = spacedim;

  /**
   * Dimension of the object over which we integrate. For the present class,
   * this is equal to <code>dim-1</code>.
   */
  static constexpr unsigned int integral_dimension = dim - 1;

  /**
   * Constructor. Gets cell independent data from mapping and finite element
   * objects, matching the quadrature rule and update flags.
   */
  FEFaceValues(const Mapping<dim, spacedim>       &mapping,
               const FiniteElement<dim, spacedim> &fe,
               const Quadrature<dim - 1>          &quadrature,
               const UpdateFlags                   update_flags);

  /**
   * Like the function above, but taking a collection of quadrature rules. This
   * allows to assign each face a different quadrature rule. In the case that
   * the collection only contains a single face quadrature, this quadrature
   * rule is use on all faces.
   */
  FEFaceValues(const Mapping<dim, spacedim>       &mapping,
               const FiniteElement<dim, spacedim> &fe,
               const hp::QCollection<dim - 1>     &quadrature,
               const UpdateFlags                   update_flags);

  /**
   * Constructor. This constructor is equivalent to the other one except that
   * it makes the object use a $Q_1$ mapping (i.e., an object of type
   * MappingQ(1)) implicitly.
   */
  FEFaceValues(const FiniteElement<dim, spacedim> &fe,
               const Quadrature<dim - 1>          &quadrature,
               const UpdateFlags                   update_flags);

  /**
   * Like the function above, but taking a collection of quadrature rules. This
   * allows to assign each face a different quadrature rule. In the case that
   * the collection only contains a single face quadrature, this quadrature
   * rule is use on all faces.
   */
  FEFaceValues(const FiniteElement<dim, spacedim> &fe,
               const hp::QCollection<dim - 1>     &quadrature,
               const UpdateFlags                   update_flags);

  /**
   * Reinitialize the gradients, Jacobi determinants, etc for the face with
   * number @p face_no of @p cell and the given finite element.
   */
  template <bool level_dof_access>
  void
  reinit(
    const TriaIterator<DoFCellAccessor<dim, spacedim, level_dof_access>> &cell,
    const unsigned int face_no);

  /**
   * Reinitialize the gradients, Jacobi determinants, etc for face @p face
   * and cell @p cell.
   *
   * @note @p face must be one of @p cell's face iterators.
   */
  template <bool level_dof_access>
  void
  reinit(
    const TriaIterator<DoFCellAccessor<dim, spacedim, level_dof_access>> &cell,
    const typename Triangulation<dim, spacedim>::face_iterator           &face);

  /**
   * Reinitialize the gradients, Jacobi determinants, etc for the given face
   * on a given cell of type "iterator into a Triangulation object", and the
   * given finite element. Since iterators into a triangulation alone only
   * convey information about the geometry of a cell, but not about degrees of
   * freedom possibly associated with this cell, you will not be able to call
   * some functions of this class if they need information about degrees of
   * freedom. These functions are, above all, the
   * <tt>get_function_value/gradients/hessians/third_derivatives</tt>
   * functions. If you want to call these functions, you have to call the @p
   * reinit variants that take iterators into DoFHandler or other DoF handler
   * type objects.
   */
  void
  reinit(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
         const unsigned int                                          face_no);

  /*
   * Reinitialize the gradients, Jacobi determinants, etc for the given face
   * on a given cell of type "iterator into a Triangulation object", and the
   * given finite element. Since iterators into a triangulation alone only
   * convey information about the geometry of a cell, but not about degrees of
   * freedom possibly associated with this cell, you will not be able to call
   * some functions of this class if they need information about degrees of
   * freedom. These functions are, above all, the
   * <tt>get_function_value/gradients/hessians/third_derivatives</tt>
   * functions. If you want to call these functions, you have to call the @p
   * reinit variants that take iterators into DoFHandler or other DoF handler
   * type objects.
   *
   * @note @p face must be one of @p cell's face iterators.
   */
  void
  reinit(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
         const typename Triangulation<dim, spacedim>::face_iterator &face);

  /**
   * Return a reference to this very object.
   *
   * Though it seems that it is not very useful, this function is there to
   * provide capability to the hp::FEValues class, in which case it provides
   * the FEValues object for the present cell (remember that for hp-finite
   * elements, the actual FE object used may change from cell to cell, so we
   * also need different FEValues objects for different cells; once you
   * reinitialize the hp::FEValues object for a specific cell, it retrieves
   * the FEValues object for the FE on that cell and returns it through a
   * function of the same name as this one; this function here therefore only
   * provides the same interface so that one can templatize on FEValues and
   * hp::FEValues).
   */
  const FEFaceValues<dim, spacedim> &
  get_present_fe_values() const;

private:
  /**
   * Do work common to the two constructors.
   */
  void
  initialize(const UpdateFlags update_flags);

  /**
   * The reinit() functions do only that part of the work that requires
   * knowledge of the type of iterator. After setting present_cell(), they
   * pass on to this function, which does the real work, and which is
   * independent of the actual type of the cell iterator.
   */
  void
  do_reinit(const unsigned int face_no);
};


/**
 * Finite element evaluated in quadrature points on a face.
 *
 * This class adds the functionality of FEFaceValuesBase to FEValues; see
 * there for more documentation.
 *
 * This class is used for faces lying on a refinement edge. In this case, the
 * neighboring cell is refined. To be able to compute differences between
 * interior and exterior function values, the refinement of the neighboring
 * cell must be simulated on this cell. This is achieved by applying a
 * quadrature rule that simulates the refinement. The resulting data fields
 * are split up to reflect the refinement structure of the neighbor: a subface
 * number corresponds to the number of the child of the neighboring face.
 *
 * @ingroup feaccess
 */
template <int dim, int spacedim = dim>
class FESubfaceValues : public FEFaceValuesBase<dim, spacedim>
{
public:
  /**
   * Dimension in which this object operates.
   */
  static constexpr unsigned int dimension = dim;

  /**
   * Dimension of the space in which this object operates.
   */
  static constexpr unsigned int space_dimension = spacedim;

  /**
   * Dimension of the object over which we integrate. For the present class,
   * this is equal to <code>dim-1</code>.
   */
  static constexpr unsigned int integral_dimension = dim - 1;

  /**
   * Constructor. Gets cell independent data from mapping and finite element
   * objects, matching the quadrature rule and update flags.
   */
  FESubfaceValues(const Mapping<dim, spacedim>       &mapping,
                  const FiniteElement<dim, spacedim> &fe,
                  const Quadrature<dim - 1>          &face_quadrature,
                  const UpdateFlags                   update_flags);

  /**
   * Like the function above, but taking a collection of quadrature rules.
   *
   * @note We require, in contrast to FEFaceValues, that the number of quadrature
   *   rules in the collection is one.
   */
  FESubfaceValues(const Mapping<dim, spacedim>       &mapping,
                  const FiniteElement<dim, spacedim> &fe,
                  const hp::QCollection<dim - 1>     &face_quadrature,
                  const UpdateFlags                   update_flags);

  /**
   * Constructor. This constructor is equivalent to the other one except that
   * it makes the object use a $Q_1$ mapping (i.e., an object of type
   * MappingQ(1)) implicitly.
   */
  FESubfaceValues(const FiniteElement<dim, spacedim> &fe,
                  const Quadrature<dim - 1>          &face_quadrature,
                  const UpdateFlags                   update_flags);

  /**
   * Like the function above, but taking a collection of quadrature rules.
   *
   * @note We require, in contrast to FEFaceValues, that the number of quadrature
   *   rules in the collection is one.
   */
  FESubfaceValues(const FiniteElement<dim, spacedim> &fe,
                  const hp::QCollection<dim - 1>     &face_quadrature,
                  const UpdateFlags                   update_flags);

  /**
   * Reinitialize the gradients, Jacobi determinants, etc for the given cell
   * of type "iterator into a DoFHandler object", and the finite element
   * associated with this object. It is assumed that the finite element used
   * by the given cell is also the one used by this FESubfaceValues object.
   */
  template <bool level_dof_access>
  void
  reinit(
    const TriaIterator<DoFCellAccessor<dim, spacedim, level_dof_access>> &cell,
    const unsigned int face_no,
    const unsigned int subface_no);

  /**
   * Alternative reinitialization function that takes, as arguments, iterators
   * to the face and subface instead of their numbers.
   */
  template <bool level_dof_access>
  void
  reinit(
    const TriaIterator<DoFCellAccessor<dim, spacedim, level_dof_access>> &cell,
    const typename Triangulation<dim, spacedim>::face_iterator           &face,
    const typename Triangulation<dim, spacedim>::face_iterator &subface);

  /**
   * Reinitialize the gradients, Jacobi determinants, etc for the given
   * subface on a given cell of type "iterator into a Triangulation object", and
   * the given finite element. Since iterators into a triangulation alone only
   * convey information about the geometry of a cell, but not about degrees of
   * freedom possibly associated with this cell, you will not be able to call
   * some functions of this class if they need information about degrees of
   * freedom. These functions are, above all, the
   * <tt>get_function_value/gradients/hessians/third_derivatives</tt>
   * functions. If you want to call these functions, you have to call the @p
   * reinit variants that take iterators into DoFHandler or other DoF handler
   * type objects.
   */
  void
  reinit(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
         const unsigned int                                          face_no,
         const unsigned int subface_no);

  /**
   * Reinitialize the gradients, Jacobi determinants, etc for the given
   * subface on a given cell of type "iterator into a Triangulation object", and
   * the given finite element. Since iterators into a triangulation alone only
   * convey information about the geometry of a cell, but not about degrees of
   * freedom possibly associated with this cell, you will not be able to call
   * some functions of this class if they need information about degrees of
   * freedom. These functions are, above all, the
   * <tt>get_function_value/gradients/hessians/third_derivatives</tt>
   * functions. If you want to call these functions, you have to call the @p
   * reinit variants that take iterators into DoFHandler or other DoF handler
   * type objects.
   *
   * This does the same thing as the previous function but takes iterators
   * instead of numbers as arguments.
   *
   * @note @p face and @p subface must correspond to a face (and a subface of
   * that face) of @p cell.
   */
  void
  reinit(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
         const typename Triangulation<dim, spacedim>::face_iterator &face,
         const typename Triangulation<dim, spacedim>::face_iterator &subface);

  /**
   * Return a reference to this very object.
   *
   * Though it seems that it is not very useful, this function is there to
   * provide capability to the hp::FEValues class, in which case it provides
   * the FEValues object for the present cell (remember that for hp-finite
   * elements, the actual FE object used may change from cell to cell, so we
   * also need different FEValues objects for different cells; once you
   * reinitialize the hp::FEValues object for a specific cell, it retrieves
   * the FEValues object for the FE on that cell and returns it through a
   * function of the same name as this one; this function here therefore only
   * provides the same interface so that one can templatize on FEValues and
   * hp::FEValues).
   */
  const FESubfaceValues<dim, spacedim> &
  get_present_fe_values() const;

  /**
   * @todo Document this
   *
   * @ingroup Exceptions
   */
  DeclException0(ExcReinitCalledWithBoundaryFace);

  /**
   * @todo Document this
   *
   * @ingroup Exceptions
   */
  DeclException0(ExcFaceHasNoSubfaces);

private:
  /**
   * Do work common to the two constructors.
   */
  void
  initialize(const UpdateFlags update_flags);

  /**
   * The reinit() functions do only that part of the work that requires
   * knowledge of the type of iterator. After setting present_cell(), they
   * pass on to this function, which does the real work, and which is
   * independent of the actual type of the cell iterator.
   */
  void
  do_reinit(const unsigned int face_no, const unsigned int subface_no);
};


#ifndef DOXYGEN


/*--------------------- Inline functions: FEValues --------------------------*/


template <int dim, int spacedim>
inline const Quadrature<dim> &
FEValues<dim, spacedim>::get_quadrature() const
{
  return quadrature;
}



template <int dim, int spacedim>
inline const FEValues<dim, spacedim> &
FEValues<dim, spacedim>::get_present_fe_values() const
{
  return *this;
}


/*---------------------- Inline functions: FEFaceValuesBase -----------------*/


template <int dim, int spacedim>
inline unsigned int
FEFaceValuesBase<dim, spacedim>::get_face_number() const
{
  return present_face_no;
}


template <int dim, int spacedim>
inline unsigned int
FEFaceValuesBase<dim, spacedim>::get_face_index() const
{
  return present_face_index;
}


/*----------------------- Inline functions: FE*FaceValues -------------------*/

template <int dim, int spacedim>
inline const Quadrature<dim - 1> &
FEFaceValuesBase<dim, spacedim>::get_quadrature() const
{
  return quadrature[quadrature.size() == 1 ? 0 : present_face_no];
}



template <int dim, int spacedim>
inline const FEFaceValues<dim, spacedim> &
FEFaceValues<dim, spacedim>::get_present_fe_values() const
{
  return *this;
}



template <int dim, int spacedim>
inline const FESubfaceValues<dim, spacedim> &
FESubfaceValues<dim, spacedim>::get_present_fe_values() const
{
  return *this;
}



template <int dim, int spacedim>
inline const Tensor<1, spacedim> &
FEFaceValuesBase<dim, spacedim>::boundary_form(const unsigned int q_point) const
{
  AssertIndexRange(q_point, this->mapping_output.boundary_forms.size());
  Assert(this->update_flags & update_boundary_forms,
         (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
           "update_boundary_forms")));

  return this->mapping_output.boundary_forms[q_point];
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif

// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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

#ifndef __deal2__fe_update_flags_h
#define __deal2__fe_update_flags_h


#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup feaccess */
/*@{*/

/**
 * The enum type given to the constructors of FEValues, FEFaceValues
 * and FESubfaceValues, telling those objects which data will be
 * needed on each mesh cell.
 *
 * Selecting these flags in a restrictive way is crucial for the
 * efficiency of FEValues::reinit(), FEFaceValues::reinit() and
 * FESubfaceValues::reinit(). Therefore, only the flags actually
 * needed should be selected. It is the responsibility of the involved
 * Mapping and FiniteElement to add additional flags according to
 * their own requirements. For instance, most finite elements will
 * add #update_covariant_transformation if #update_gradients is
 * selected.

 * By default, all flags are off, i.e. no reinitialization will be
 * done.
 *
 * You can select more than one flag by concatenation
 * using the bitwise or operator|(UpdateFlags,UpdateFlags).
 *
 * <h3>Generating the actual flags</h3>
 *
 * When given a set of UpdateFlags @p flags, the FEValues object must
 * determine, which values will have to be computed once only for the
 * reference cell and which values will have to be updated for each
 * cell. Here, it is important to note that in many cases, the
 * FiniteElement will require additional updates from the Mapping. To
 * this end, several auxiliary functions have been implemented:
 *
 * FiniteElement::update_once(flags) and
 * FiniteElement::update_each(flags) determine the values required by
 * the FiniteElement once or on each cell. The same functions exist in Mapping.
 *
 * Since the FiniteElement does not know if a value required from
 * Mapping should be computed once or for each cell,
 * FEValuesBase::compute_update_flags() is used to compute the union
 * of all values to be computed ever. It does this by first adding to
 * the flags set by the user all flags (once and each) added by the
 * FiniteElement. This new set of flags is then given to the Mapping
 * and all flags required there are added, again once and each.
 *
 * This union of all flags is given to Mapping::fill_fe_values() and
 * FiniteElement::fill_fe_values, where it is split again into the
 * information generated only once and the information that must be
 * updated on each cell.
 *
 * The flags finally stored in FEValues then are the union of all the
 * flags required by the user, by FiniteElement and by Mapping, for
 * computation once or on each cell. Subsequent calls to the functions
 * @p update_once and @p update_each should just select among these
 * flags, but should not add new flags.
 *
 * The mechanism by which all this is accomplished is also discussed
 * on the page on @ref UpdateFlagsEssay.
 */
enum UpdateFlags
{
  //! No update
  update_default = 0,
  //! Shape function values
  /**
   * Compute the values of the shape functions at the quadrature points on
   * the real space cell. For the usual Lagrange elements, these values are
   * equal to the values of the shape functions at the quadrature points on
   * the unit cell, but they are different for more complicated elements,
   * such as FE_RaviartThomas elements.
   */
  update_values = 0x0001,
  //! Shape function gradients
  /**
   * Compute the gradients of the shape functions in coordinates of the
   * real cell.
   */
  update_gradients = 0x0002,
  //! Second derivatives of shape functions
  /**
   * Compute the second derivatives of the shape functions in coordinates
   * of the real cell.
   */
  update_hessians = 0x0004,
  //! Outer normal vector, not normalized
  /**
   * Vector product of tangential vectors, yielding a normal vector with
   * a length corresponding to the surface element; may be more efficient
   * than computing both.
   */
  update_boundary_forms = 0x0008,
  //! Transformed quadrature points
  /**
   * Compute the quadrature points transformed into real cell coordinates.
   */
  update_quadrature_points = 0x0010,
  //! Transformed quadrature weights
  /**
   * Compute the quadrature weights on the real cell, i.e. the weights
   * of the quadrature rule multiplied with the determinant of the Jacoian
   * of the transformation from reference to realcell.
   */
  update_JxW_values = 0x0020,
  //! Normal vectors
  /**
   * Compute the normal vectors, either for a face or for a cell of
   * codimension one. Setting this flag for any other object will raise
   * an error.
   */
  update_normal_vectors = 0x0040,
  /**
   * @deprecated Use #update_normal_vectors
   */
  update_face_normal_vectors = update_normal_vectors,
  /**
   * @deprecated Use #update_normal_vectors
   */
  update_cell_normal_vectors = update_normal_vectors,
  //! Volume element
  /**
   * Compute the Jacobian of the transformation from the reference cell
   * to the real cell.
   */
  update_jacobians = 0x0080,
  //! Gradient of volume element
  /**
   * Compute the dervatives of the Jacobian of the transformation.
   */
  update_jacobian_grads = 0x0100,
  //! Volume element
  /**
   * Compute the inverse Jacobian of the transformation from the reference
   * cell to the real cell.
   */
  update_inverse_jacobians = 0x0200,
  //! Covariant transformation
  /**
   * Compute all values the Mapping needs to perform a contravariant
   * transformation of vectors. For special mappings like MappingCartesian
   * this may be simpler than #update_inverse_jacobians.
   */
  update_covariant_transformation = 0x0400,
  //! Contravariant transformation
  /**
   * Compute all values the Mapping needs to perform a contravariant
   * transformation of vectors. For special mappings like MappingCartesian
   * this may be simpler than #update_jacobians.
   */
  update_contravariant_transformation = 0x0800,
  //! Shape function values of transformation
  /**
   * Compute the shape function values of the transformation defined by
   * the Mapping.
   */
  update_transformation_values = 0x1000,
  //! Shape function gradients of transformation
  /**
   * Compute the shape function gradients of the transformation defined
   * by the Mapping.
   */
  update_transformation_gradients = 0x2000,
  //! Determinant of the Jacobian
  /**
   * Compute the volume element in each quadrature point.
   */
  update_volume_elements = 0x4000,
  /**
   * @deprecated: this flag has no effect.
   */
  update_support_points = 0x10000,
  /**
   * Update the Jacobian of the mapping in generalized support points.
   */
  update_support_jacobians = 0x20000,
  /**
   * Update the inverse Jacobian of the mapping in generalized support
   * points.
   */
  update_support_inverse_jacobians = 0x40000,
  /**
   * @deprecated Update quadrature points
   */
  update_q_points = update_quadrature_points,
  /**
   * @deprecated Update second derivatives.
   */
  update_second_derivatives = update_hessians,
  //! Values needed for Piola transform
  /**
   * Combination of the flags needed for Piola transform of Hdiv elements.
   */
  update_piola = update_volume_elements | update_contravariant_transformation
};


/**
 * Output operator which outputs update flags as a set of or'd text values.
 *
 * @ref UpdateFlags
 */
template <class STREAM>
inline
STREAM &operator << (STREAM &s, UpdateFlags u)
{
  s << " UpdateFlags|";
  if (u & update_values)                       s << "values|";
  if (u & update_gradients)                    s << "gradients|";
  if (u & update_hessians)                     s << "hessians|";
  if (u & update_quadrature_points)            s << "quadrature_points|";
  if (u & update_JxW_values)                   s << "JxW_values|";
  if (u & update_normal_vectors)               s << "normal_vectors|";
  if (u & update_jacobians)                    s << "jacobians|";
  if (u & update_inverse_jacobians)            s << "inverse_jacobians|";
  if (u & update_jacobian_grads)               s << "jacobian_grads|";
  if (u & update_covariant_transformation)     s << "covariant_transformation|";
  if (u & update_contravariant_transformation) s << "contravariant_transformation|";
  if (u & update_transformation_values)        s << "transformation_values|";
  if (u & update_transformation_gradients)     s << "transformation_gradients|";
  if (u & update_support_points)               s << "support_points|";
  if (u & update_support_jacobians)            s << "support_jacobians|";
  if (u & update_support_inverse_jacobians)    s << "support_inverse_jacobians|";

//TODO: check that 'u' really only has the flags set that are handled above
  return s;
}


/**
 * Global operator which returns an object in which all bits are set
 * which are either set in the first or the second argument. This
 * operator exists since if it did not then the result of the bit-or
 * <tt>operator |</tt> would be an integer which would in turn trigger
 * a compiler warning when we tried to assign it to an object of type
 * UpdateFlags.
 *
 * @ref UpdateFlags
 */
inline
UpdateFlags
operator | (UpdateFlags f1, UpdateFlags f2)
{
  return static_cast<UpdateFlags> (
           static_cast<unsigned int> (f1) |
           static_cast<unsigned int> (f2));
}




/**
 * Global operator which sets the bits from the second argument also
 * in the first one.
 *
 * @ref UpdateFlags
 */
inline
UpdateFlags &
operator |= (UpdateFlags &f1, UpdateFlags f2)
{
  f1 = f1 | f2;
  return f1;
}


/**
 * Global operator which returns an object in which all bits are set
 * which are set in the first as well as the second argument. This
 * operator exists since if it did not then the result of the bit-and
 * <tt>operator &</tt> would be an integer which would in turn trigger
 * a compiler warning when we tried to assign it to an object of type
 * UpdateFlags.
 *
 * @ref UpdateFlags
 */
inline
UpdateFlags
operator & (UpdateFlags f1, UpdateFlags f2)
{
  return static_cast<UpdateFlags> (
           static_cast<unsigned int> (f1) &
           static_cast<unsigned int> (f2));
}


/**
 * Global operator which clears all the bits in the first argument if
 * they are not also set in the second argument.
 *
 * @ref UpdateFlags
 */
inline
UpdateFlags &
operator &= (UpdateFlags &f1, UpdateFlags f2)
{
  f1 = f1 & f2;
  return f1;
}



/**
 * This enum definition is used for storing similarities of the current cell
 * to the previously visited cell. This information is used for reusing data
 * when calling the method FEValues::reinit() (like derivatives, which do
 * not change if one cell is just a translation of the previous). Currently,
 * this variable does only recognize a translation and an inverted translation
 * (if dim<spacedim). However, this concept
 * makes it easy to add additional staties to be detected in
 * FEValues/FEFaceValues for making use of these similarities as well.
 */
namespace CellSimilarity
{
  enum Similarity
  {
    none,
    translation,
    inverted_translation,
    invalid_next_cell
  };
}


/*@}*/



DEAL_II_NAMESPACE_CLOSE

#endif

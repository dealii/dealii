//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__fe_update_flags_h
#define __deal2__fe_update_flags_h


#include <base/config.h>


/*!@addtogroup febase */
/*@{*/

/**
 * Provide a set of flags which tells the FEValues::reinit
 * function, which fields are to be updated for each cell. E.g. if you
 * do not need the gradients since you want to assemble the mass
 * matrix, you can switch that off. By default, all flags are off,
 * i.e. no reinitialization will be done.
 *
 * A variable of this type has to be passed to the constructor of the
 * FEValues object. You can select more than one flag by concatenation
 * using the <tt>|</tt> (bitwise <tt>or</tt>) operator.
 *
 *
 * @sect2{Description of Flags}
 *
 * The following flags are declared:
 * <ul>
 * <li> <tt>update_default</tt>: Default: update nothing.
 * <li> <tt>update_values</tt>: Compute the values of the shape
 *     functions at the quadrature points on the real space cell. For
 *     the usual Lagrange elements, these values are equal to the
 *     values of the shape functions at the quadrature points on the
 *     unit cell, but they are different for more complicated
 *     elements, such as BDM or Raviart-Thomas elements.
 * <li> <tt>update_gradients</tt>: Transform gradients on unit cell
 *     to gradients on real cell.
 * <li> <tt>update_second_derivatives</tt>: Update the second
 *     derivatives of the shape functions on the real cell.
 * <li> <tt>update_boundary_forms</tt>: Update boundary forms on the
 *     face.  This flag is only evaluated by the FEFaceValues class.
 *     Giving this flag to the FEValues class will result in an error,
 *     since boundary forms only exist on the boundary.
 * <li> <tt>update_q_points</tt>: Compute quadrature points in real
 *     space (not on unit cell).
 * <li> <tt>update_JxW_values</tt>: Compute the JxW values (Jacobian
 *     determinant at the quadrature point times the weight of this
 *     point).
 * <li> <tt>update_normal_vectors</tt>: Update the outward normal
 *     vectors to the face relative to this cell.  This flag is only
 *     evaluated by the FEFaceValues class.  Giving this flag to the
 *     FEValues class will result in an error, since normal vectors
 *     are not useful in that case.
 * <li> <tt>update_jacobians</tt>: Compute jacobian matrices of the
 *     transform between unit and real cell in the evaluation points.
 * <li> <tt>update_jacobian_grads</tt>: Update gradients of the
 *     jacobian. These are used to compute second derivatives.
 * <li> <tt>update_covariant_transformation</tt>: Update co-variant
 *     transformation.  This flag is used internally to tell Mapping
 *     objects to compute the transformation matrices for co-variant
 *     vectors.
 * <li> <tt>update_contravariant_transformation</tt>: Update
 *     contra-variant transformation.  This flag is used internally to
 *     tell Mapping objects to compute the transformation matrices for
 *     contra-variant vectors.
 * <li> <tt>update_transformation_values</tt>: Update the shape
 *     function values of the transformation.
 * <li> <tt>update_transformation_gradients</tt>: Update the
 *     gradients of the shape functions of the transformation.
 * </ul>
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1998, 1999, 2000, 2001, Ralf Hartmann 2004
 */
enum UpdateFlags
{
				       //! No update
      update_default                      = 0,
				       //! Shape function values
      update_values                       = 0x0001,
				       //! Shape function gradients
      update_gradients                    = 0x0002,
				       //! Second derivatives of shape functions
      update_second_derivatives           = 0x0004,
				       /*! Normal vector times surface
					* element; may be more
					* efficient than computing both.
					*/
      update_boundary_forms               = 0x0008,
				       //! Transformed quadrature points
      update_q_points                     = 0x0010,
				       //! Transformed quadrature weights
      update_JxW_values                   = 0x0020,
				       //! Transformed normal vectors
      update_normal_vectors               = 0x0040,
				       //! Volume element
      update_jacobians                    = 0x0080,
				       //! Gradient of volume element
      update_jacobian_grads               = 0x0100,
				       //! Covariant transformation
      update_covariant_transformation     = 0x0200,
				       //! Contravariant transformation
      update_contravariant_transformation = 0x0400,
				       //! Shape function values of transformation
      update_transformation_values        = 0x0800,
				       //! Shape function gradients of transformation
      update_transformation_gradients     = 0x1000
};





/**
 * Global operator which returns an object in which all bits are set
 * which are either set in the first or the second argument. This
 * operator exists since if it did not then the result of the bit-or
 * <tt>operator |</tt> would be an integer which would in turn trigger
 * a compiler warning when we tried to assign it to an object of type
 * UpdateFlags.
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
 */
inline
UpdateFlags &
operator &= (UpdateFlags &f1, UpdateFlags f2)
{
  f1 = f1 & f2;
  return f1;
}



/*@}*/

#endif

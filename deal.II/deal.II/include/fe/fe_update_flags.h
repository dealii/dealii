//----------------------------  fe_update_flags.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_update_flags.h  ---------------------------
#ifndef __deal2__fe_update_flags_h
#define __deal2__fe_update_flags_h


/**
 * Provide a set of flags which tells the @p{FEValues<>::reinit} function, which
 * fields are to be updated for each cell. E.g. if you do not need the
 * gradients since you want to assemble the mass matrix, you can switch that
 * off. By default, all flags are off, i.e. no reinitialization will be done.
 *
 * A variable of this type has to be passed to the constructor of the
 * @p{FEValues} object. You can select more than one flag by concatenation
 * using the @p{|} (bitwise @p{or}) operator.
 *
 *
 * @sect2{Description of Flags}
 *
 * The following flags are declared:
 * @begin{itemize}
 *   @item @p{update_default}: Default: update nothing.
 *   @item @p{update_values}: Compute the values of the shape
 *     functions at the quadrature points on the real space cell. For the
 *     usual Lagrange elements, these values are equal to the values of
 *     the shape functions at the quadrature points on the unit cell, but
 *     they are different for more complicated elements, such as BDM or
 *     Raviart-Thomas elements.
 *   @item @p{update_q_points}: Compute quadrature points in real
 *     space (not on unit cell).
 *   @item @p{update_gradients}: Transform gradients on unit cell to
 *     gradients on real cell.
 *   @item @p{update_jacobians}: Compute jacobian matrices of the
 *     transform between unit and real cell
 *     in the evaluation points.
 *   @item @p{update_JxW_values}: Compute the JxW values (Jacobian
 *     determinant at the quadrature point
 *     times the weight of this point).
 *   @item @p{update_normal_vectors}: Update the outward normal vectors
 *     to the face relative to this cell.
 *     This flag is only evaluated by
 *     the @p{FEFaceValues} class.
 *     Giving this flag to the
 *     @p{FEValues} class will result in
 *     an error, since normal vectors are
 *     not useful in that case.
 *   @item @p{update_second_derivatives}: Update the second derivatives of the
 *     shape functions on the real cell.
 * @end{itemize}
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1998, 1999, 2000, 2001
 */
enum UpdateFlags {
				       /**
					* Default: update nothing.
					*/
      update_default  = 0,
				       /**
					* Compute quadrature points in
					* real space (not on unit
					* cell).
					*/
      update_q_points = 0x0001,

				       /**
					* Compute the JxW values
					* (Jacobian determinant at the
					* quadrature point times the
					* weight of this point).
					*/
      update_JxW_values = 0x0002,

				       /**
					* Update boundary forms on the
				        * face.  This flag is only
				        * evaluated by the
				        * @ref{FEFaceValues} class.
					*
					* Giving this flag to the
					* @ref{FEValues} class will
					* result in an error, since
					* boundary forms only exist on
					* the boundary.
					*/
      update_boundary_forms = 0x0004,

				       /**
					* Compute the values of the
					* shape functions at the
					* quadrature points on the
					* real space cell. For the
					* usual Lagrange elements,
					* these values are equal to
					* the values of the shape
					* functions at the quadrature
					* points on the unit cell, but
					* they are different for more
					* complicated elements, such
					* as BDM or Raviart-Thomas
					* elements.
					*/
      update_values = 0x0010,
				       /**
					* Transform gradients on unit
					* cell to gradients on real
					* cell.
					*/
      update_gradients = 0x0020,

				       /**
					* Update the second
					* derivatives of the shape
					* functions on the real cell.
				        */
      update_second_derivatives = 0x0040,

				       /**
					* Compute jacobian matrices of
					* the transform between unit
					* and real cell in the
					* evaluation points.
					*/
      update_jacobians = 0x0080,

				       /**
					* Update co-variant
					* transformation.  This flag
					* is used internally to tell
					* Mapping objects to compute
					* the transformation matrices
					* for co-variant vectors.
					*/
      update_covariant_transformation = 0x0100,

				       /**
					* Update contra-variant
					* transformation.  This flag
					* is used internally to tell
					* Mapping objects to compute
					* the transformation matrices
					* for contra-variant vectors.
					*/
      update_contravariant_transformation = 0x0200,

				       /**
					* Update gradients of the
					* jacobian.  These are used to
					* compute second derivatives.
					*/
      update_jacobian_grads = 0x0400,

				       /**
					* Update shape function values
					* of the transformation.
					*/

      update_transformation_values = 0x0800,

				       /**
					* Update gradients of
					* transformation shape
					* functions.
					*/
      update_transformation_gradients = 0x1000,

				       /**
					* Compute the points on the
					* real cell on which the trial
					* functions are located.
					*
					* Giving this flag to the
					* @ref{FESubfaceValues} class
					* will result in an error,
					* since support points are not
					* useful in that case.
					*/
      update_support_points = 0x2000,

				       /**
				        * Update the outward normal
				        * vectors to the face relative
				        * to this cell.  This flag is
				        * only evaluated by the
				        * @ref{FEFaceValues} class.
					*
					* Giving this flag to the
					* @ref{FEValues} class will
					* result in an error, since
					* normal vectors are not
					* useful in that case.
				        */
      update_normal_vectors = 0x4000
};


inline
UpdateFlags&
operator |= (UpdateFlags& f1, const UpdateFlags& f2)
{
  return ((UpdateFlags) (((int)f1) |= f2));
}


inline
UpdateFlags
operator | (const UpdateFlags& f1, const UpdateFlags& f2)
{
  UpdateFlags result = f1;
  result |= f2;
  return result;
}

inline
UpdateFlags&
operator &= (UpdateFlags& f1, const UpdateFlags& f2)
{
  return ((UpdateFlags) (((int)f1) &= f2));
}


inline
UpdateFlags
operator & (const UpdateFlags& f1, const UpdateFlags& f2)
{
  UpdateFlags result = f1;
  result &= f2;
  return result;
}

#endif

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
 *   @item @p{update_default  = 0}: Default: update nothing.
 *   @item @p{update_values = 1}: Compute the values of the shape
 *     functions at the quadrature points on the real space cell. For the
 *     usual Lagrange elements, these values are equal to the values of
 *     the shape functions at the quadrature points on the unit cell, but
 *     they are different for more complicated elements, such as BDM or
 *     Raviart-Thomas elements.
 *   @item @p{update_q_points = 2}: Compute quadrature points in real
 *     space (not on unit cell).
 *   @item @p{update_gradients = 4}: Transform gradients on unit cell to
 *     gradients on real cell.
 *   @item @p{update_jacobians = 8}: Compute jacobian matrices of the
 *     transform between unit and real cell
 *     in the evaluation points.
 *   @item @p{update_JxW_values = 16}: Compute the JxW values (Jacobian
 *     determinant at the quadrature point
 *     times the weight of this point).
 *   @item @p{update_support_points = 32}: Compute the points on the real cell
 *     on which the trial functions are located.
 *     Giving this flag to the
 *     @p{FESubfaceValues} class will result
 *     in an error, since support points are
 *     not useful in that case.
 *   @item @p{update_normal_vectors = 64}: Update the outward normal vectors
 *     to the face relative to this cell.
 *     This flag is only evaluated by
 *     the @p{FEFaceValues} class.
 *     Giving this flag to the
 *     @p{FEValues} class will result in
 *     an error, since normal vectors are
 *     not useful in that case.
 *   @item @p{update_second_derivatives = 128}: Update the second derivatives of the
 *     shape functions on the real cell.
 * @end{itemize}
 *
 * @author Wolfgang Bangerth, 1998, 1999
 */
enum UpdateFlags
{
      update_default            = 0,
      update_values             = 1,
      update_q_points           = 2,
      update_gradients          = 4,
      update_jacobians          = 8,
      update_JxW_values         = 16,
      update_support_points     = 32,
      update_normal_vectors     = 64,
      update_second_derivatives = 128
};


#endif

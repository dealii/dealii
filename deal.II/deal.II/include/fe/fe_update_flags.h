/*----------------------------   fe_update_flags.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __fe_update_flags_H
#define __fe_update_flags_H
/*----------------------------   fe_update_flags.h     ---------------------------*/


/**
 * Provide a set of flags which tells the #FEValues<>::reinit# function, which
 * fields are to be updated for each cell. E.g. if you do not need the
 * gradients since you want to assemble the mass matrix, you can switch that
 * off. By default, all flags are off, i.e. no reinitialization will be done.
 *
 * A variable of this type has to be passed to the constructor of the
 * #FEValues# object. You can select more than one flag by concatenation
 * using the #|# (bitwise #or#) operator.
 */
enum UpdateFlags {
				       /**
					* Default: update nothing.
					*/
      update_default  = 0,
				       /**
					* Compute quadrature points in real
					* space (not on unit cell).
					*/
      update_q_points = 1,
				       /**
					* Transform gradients on unit cell to
					* gradients on real cell.
					*/
      update_gradients = 2,
				       /**
					* Compute jacobian matrices of the
					* transform between unit and real cell
					* in the evaluation points.
					*/
      update_jacobians = 4,
				       /**
					* Compute the JxW values (Jacobian
					* determinant at the quadrature point
					* times the weight of this point).
					*/
      update_JxW_values = 8,
				       /**
					* Compute the points on the real cell
					* on which the trial functions are
					* located.
					*
					* Giving this flag to the
					* #FESubfaceValues# class will result
					* in an error, since support points are
					* not useful in that case.
					*/
      update_support_points = 16,
				       /**
				        * Update the outward normal vectors
				        * to the face relative to this cell.
				        * This flag is only evaluated by
				        * the #FEFaceValues# class.
					*
					* Giving this flag to the
					* #FEValues# class will result in
					* an error, since normal vectors are
					* not useful in that case.
				        */
      update_normal_vectors = 32,

				       /**
					* Update the second derivatives of the
					* shape functions on the real cell.
				        */
      update_second_derivatives
};



/*----------------------------   fe_update_flags.h     ---------------------------*/
/* end of #ifndef __fe_update_flags_H */
#endif
/*----------------------------   fe_update_flags.h     ---------------------------*/

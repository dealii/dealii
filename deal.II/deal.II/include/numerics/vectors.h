/*----------------------------   vectors.h     ---------------------------*/
/*      $Id$                 */
#ifndef __vectors_H
#define __vectors_H
/*----------------------------   vectors.h     ---------------------------*/


template <int dim> class DoFHandler;
template <int dim> class Function;
template <int dim> class Quadrature;
template <int dim> class Boundary;
template <int dim> class FiniteElement;
class ConstraintMatrix;
class dVector;



/**
 * Provide a class which assembles some standard vectors. Among these are
 * interpolations and projections of continuous functions to the finite
 * element space and other operations.
 */
template <int dim>
class VectorCreator {
  public:
				     /**
				      * Compute the interpolation of
				      * #function# at the ansatz points to
				      * the finite element space.
				      */
    static void interpolate (const DoFHandler<dim>    &dof,
			     const FiniteElement<dim> &fe,
			     const Function<dim>      &function,
			     dVector                  &vec);

				     /**
				      * Compute the projection of
				      * #function# to the finite element space.
				      */
    static void project (const DoFHandler<dim>    &dof,
			 const ConstraintMatrix   &constraints,
			 const FiniteElement<dim> &fe,
			 const Quadrature<dim>    &q,
			 const Boundary<dim>      &boundary,
			 const Function<dim>      &function,
			 dVector                  &vec);
};



/*----------------------------   vectors.h     ---------------------------*/
/* end of #ifndef __vectors_H */
#endif
/*----------------------------   vectors.h     ---------------------------*/

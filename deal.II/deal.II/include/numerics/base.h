/*----------------------------   problem_base.h     ---------------------------*/
/*      $Id$                 */
#ifndef __problem_base_H
#define __problem_base_H
/*----------------------------   problem_base.h     ---------------------------*/


#include <lac/dsmatrix.h>


// forward declaration
template <int dim> class Triangulation;
template <int dim> class DoFHandler;
template <int dim> class FiniteElement;
template <int dim> class Quadrature;





template <int dim>
class ProblemBase {
  public:
    ProblemBase (Triangulation<dim> *tria,
		 DoFHandler<dim>    *dof_handler);

    virtual void assemble (const Equation<dim> &equation,
			   const Quadrature<dim> &q);

  protected:
    Triangulation<dim> *tria;
    DoFHandler<dim>    *dof_handler;

    dSMatrix            system_matrix;
    dSMatrixStruct      system_sparsity;

    vector<dVector*>    right_hand_sides;
    dVector             solution;

  friend class Assembler<dim>;
};

    


/*----------------------------   problem_base.h     ---------------------------*/
/* end of #ifndef __problem_base_H */
#endif
/*----------------------------   problem_base.h     ---------------------------*/

/* $Id$ */


#include <basic/function.h>
#include <grid/dof.h>
#include <grid/dof_constraints.h>
#include <fe/fe.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <lac/dvector.h>
#include <lac/dsmatrix.h>

#include "../../../mia/control.h"
#include "../../../mia/vectormemory.h"
#include "../../../mia/cg.h"


template <int dim>
void VectorCreator<dim>::interpolate (const DoFHandler<dim>    &dof,
				      const FiniteElement<dim> &fe,
				      const Function<dim>      &function,
				      dVector                  &vec) {
  ;
};




template <int dim>
void VectorCreator<dim>::project (const DoFHandler<dim>    &dof,
				  const ConstraintMatrix   &constraints,
				  const FiniteElement<dim> &fe,
				  const Quadrature<dim>    &q,
				  const Boundary<dim>      &boundary,
				  const Function<dim>      &function,
				  dVector                  &vec) {
  vec.reinit (dof.n_dofs());
  
  dSMatrixStruct sparsity(dof.n_dofs(),
			  dof.n_dofs(),
			  dof.max_couplings_between_dofs());
  dof.make_sparsity_pattern (sparsity);

  dSMatrix mass_matrix (sparsity);
  dVector tmp (mass_matrix.n());
  MatrixCreator<dim>::create_mass_matrix (dof, fe, q, boundary,
					  mass_matrix, function, vec);

  int    max_iter  = 4000;
  double tolerance = 1.e-16;
  Control                          control1(max_iter,tolerance);
  PrimitiveVectorMemory<dVector>   memory(tmp.size());
  CG<dSMatrix,dVector>             cg(control1,memory);

				   // solve
  cg (mass_matrix, vec, tmp);
				   // distribute solution
  constraints.distribute (vec);
};




template VectorCreator<1>;
template VectorCreator<2>;

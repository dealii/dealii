/* $Id$ */


#include <basic/function.h>
#include <grid/dof.h>
#include <grid/dof_accessor.h>
#include <grid/tria_iterator.h>
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
				      const Boundary<dim>      &boundary,
				      const Function<dim>      &function,
				      dVector                  &vec) {
  DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
					endc = dof.end();
  vector<int>         dofs_on_cell (fe.total_dofs);
  vector<double>      dof_values_on_cell (fe.total_dofs);
  vector<Point<dim> > ansatz_points (fe.total_dofs);
  for (; cell!=endc; ++cell) 
    {
				       // for each cell:
				       // get location of finite element
				       // off-points
      fe.get_ansatz_points (cell, boundary, ansatz_points);
				       // get function values at these points
      function.value_list (ansatz_points, dof_values_on_cell);
				       // get indices of the dofs on this cell
      cell->get_dof_indices (dofs_on_cell);
				       // distribute function values to the
				       // whole vector
      for (unsigned int i=0; i<fe.total_dofs; ++i)
	vec(dofs_on_cell[i]) = dof_values_on_cell[i];
    };
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
  constraints.condense (sparsity);
  
  dSMatrix mass_matrix (sparsity);
  dVector tmp (mass_matrix.n());
  MatrixCreator<dim>::create_mass_matrix (dof, fe, q, boundary,
					  mass_matrix, function, tmp);

  constraints.condense (mass_matrix);
  constraints.condense (tmp);

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

//----------------------------  base.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  base.cc  ---------------------------


#include <numerics/assembler.h>
#include <numerics/base.h>
#include <numerics/matrices.h>
#include <numerics/vectors.h>
#include <dofs/dof_constraints.h>
#include <grid/tria_iterator.h>
#include <numerics/data_out.h>
#include <dofs/dof_tools.h>
#include <base/function.h>
#include <fe/fe.h>
#include <base/quadrature.h>
#include <lac/vector.h>
#include <lac/precondition.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>

#include <map>
#include <numeric>
#include <algorithm>
#include <cmath>


template <int dim>
ProblemBase<dim>::ProblemBase () :
		tria(0),
		dof_handler(0),
		system_sparsity(),        // dummy initialisation, is later reinit'd
		system_matrix()           // dummy initialisation, is later reinit'd
{};


template <int dim>
void ProblemBase<dim>::set_tria_and_dof (Triangulation<dim> *t,
					 DoFHandler<dim>    *d)
{
  tria        = t;
  dof_handler = d;

				   // allow a user to reset both tria and
				   // dof to null pointers, but don't allow
				   // something other
  if ((tria != 0) || (dof_handler != 0))
    {
      Assert ((tria!=0) && (dof_handler!=0), ExcNoTriaSelected());
      Assert (tria == &dof_handler->get_tria(), ExcDofAndTriaDontMatch());
    };
};


template <int dim>
void ProblemBase<dim>::clear () {
  if (tria)        { delete tria;         tria        = 0; };
  if (dof_handler) { delete dof_handler;  dof_handler = 0; };
  system_sparsity.reinit (0,0,1);
  system_matrix.clear ();
  right_hand_side.reinit (0);
  solution.reinit (0);
  constraints.clear ();
};


template <int dim>
ProblemBase<dim>::~ProblemBase () {};


template <int dim>
void ProblemBase<dim>::assemble (const Equation<dim>      &equation,
				 const Quadrature<dim>    &quadrature,
				 const UpdateFlags         update_flags,
				 const FunctionMap        &dirichlet_bc)
{
  Assert ((tria!=0) && (dof_handler!=0), ExcNoTriaSelected());
  
  system_sparsity.reinit (dof_handler->n_dofs(),
			  dof_handler->n_dofs(),
			  dof_handler->max_couplings_between_dofs());
  right_hand_side.reinit (dof_handler->n_dofs());
  
				   // make up sparsity pattern and
				   // compress with constraints
  constraints.clear ();
  DoFTools::make_hanging_node_constraints (*dof_handler, constraints);
  constraints.close ();
  DoFTools::make_sparsity_pattern (*dof_handler, system_sparsity);
  constraints.condense (system_sparsity);

				   // reinite system matrix
  system_matrix.reinit (system_sparsity);
				   // reinit solution vector, preset
				   // with zeroes.
  solution.reinit (dof_handler->n_dofs());
  
				   // create assembler
  Assembler<dim>::AssemblerData data (*dof_handler,
				      true, true, //assemble matrix and rhs
				      system_matrix,
				      right_hand_side,
				      quadrature,
				      update_flags);
  active_assemble_iterator assembler (tria,
				      tria->begin_active()->level(),
				      tria->begin_active()->index(),
				      &data);
				   // loop over all cells, fill matrix and rhs
  do 
    {
      assembler->assemble (equation);
    }
  while ((++assembler).state() == valid);

				   // condense system matrix in-place
  constraints.condense (system_matrix);

				   // condense right hand side in-place
  constraints.condense (right_hand_side);

				   // apply Dirichlet bc as described
				   // in the docs
  map<unsigned int, double> boundary_value_list;

  for (typename FunctionMap::const_iterator dirichlet = dirichlet_bc.begin() ;
       dirichlet != dirichlet_bc.end() ; ++dirichlet)
    VectorTools::interpolate_boundary_values (*dof_handler,
					      dirichlet->first,
					      *(dirichlet->second), 
					      boundary_value_list);
  MatrixTools<dim>::apply_boundary_values (boundary_value_list,
					   system_matrix, solution,
					   right_hand_side,
					   true);
};


template <int dim>
void ProblemBase<dim>::solve ()
{
  Assert ((tria!=0) && (dof_handler!=0), ExcNoTriaSelected());
  
  SolverControl           control(4000, 1e-16);
  PrimitiveVectorMemory<> memory;
  SolverCG<>              cg(control,memory);

  PreconditionRelaxation<>
    prec(system_matrix,
	 &SparseMatrix<double>::template precondition_SSOR<double>,
	 1.2);

				   // solve
  cg.solve (system_matrix, solution, right_hand_side, prec);
				   // distribute solution
  constraints.distribute (solution);
};


template <int dim>
void ProblemBase<dim>::fill_data (DataOut<dim> &out) const {
  Assert ((tria!=0) && (dof_handler!=0), ExcNoTriaSelected());
  
  out.clear ();
  out.attach_dof_handler (*dof_handler);

  out.add_data_vector (solution, get_solution_name());
};


template <int dim>
string ProblemBase<dim>::get_solution_name () const {
  return "solution";
};


// explicit instantiations
template class ProblemBase<deal_II_dimension>;

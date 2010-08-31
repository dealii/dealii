#include "../include/linear_elastic.h"

template <int dim>
LinearElastic<dim>::LinearElastic()
  :
  dbc(dim),
  nbc(dim),
  bf(dim)
{}


template <int dim>
LinearElastic<dim>::~LinearElastic()
{
  A.clear();
}

template <int dim>
void LinearElastic<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Dirichlet Data");
  Functions::ParsedFunction<dim>::declare_parameters(prm, dim);
  prm.leave_subsection();

  prm.enter_subsection("f - Body Source");
  Functions::ParsedFunction<dim>::declare_parameters(prm, dim);
  prm.leave_subsection();

  prm.enter_subsection("Elastic Moduli");
  ParsedSymmetricTensorFunction<4, dim>::declare_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Neumann Data");
  Functions::ParsedFunction<dim>::declare_parameters(prm, dim);
  prm.leave_subsection();

}

template <int dim>
void LinearElastic<dim>::parse_parameters(ParameterHandler &prm)
{ 
  prm.enter_subsection("f - Body Source");
  bf.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Neumann Data");
  nbc.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Elastic Moduli");
  C.parse_parameters(prm);
  prm.leave_subsection();

  prm.enter_subsection("Dirichlet Data");
  dbc.parse_parameters(prm);
  prm.leave_subsection();

}

template <int dim>
void LinearElastic<dim>::reinit(VectorSpace<dim> &vspace)
{

  sp_A.reinit(vspace.n_dofs(),
	      vspace.n_dofs(),
	      vspace.get_dh().max_couplings_between_dofs());
  
  DoFTools::make_sparsity_pattern (static_cast<DoFHandler<dim> &> (vspace.get_dh()), 
				   sp_A);
  vspace.get_hang().condense (sp_A);
  sp_A.compress();
  
  A.reinit(sp_A);

  sol_total.reinit(vspace.n_dofs());
  sol_increment.reinit(vspace.n_dofs());
  rhs.reinit(vspace.n_dofs());

}

template <int dim>
void LinearElastic<dim>::build_matrix(ParameterHandler &prm,
				      VectorSpace<dim> &vspace)
{

  LocalAssembleElasticMatrix<dim> local_elastic_matrix;
  local_elastic_matrix.reinit(vspace.get_fe());
  local_elastic_matrix.parameters(prm);

  MyTools::assemble(vspace.get_dh(), vspace.get_hang(),
		    -1, A, local_elastic_matrix);

}

template <int dim>
void LinearElastic<dim>::reinit_step(double &time)
{
  
  rhs = 0;
  //sol_increment = 0;
  bf.set_time(time);
  nbc.set_time(time);
  dbc.set_time(time);
}

template <int dim>
void LinearElastic<dim>::build_rhs(VectorSpace<dim> &vspace)
{

  LocalAssembleElasticRHS<dim> elastic_local_rhs;

  elastic_local_rhs.reinit(vspace.get_fe(),
			   bf, nbc,
			   vspace.neumann_bc);
  
  MyTools::assemble_rhs(vspace.get_dh(), vspace.get_hang(),
			-1, rhs, elastic_local_rhs);
  
}

template <int dim>
void LinearElastic<dim>::solve(VectorSpace<dim> &vspace,
			       double &tolerance)
{

  //get the DBC map
  std::map<unsigned int, double>  bv;
  vspace.interpolate_dirichlet_bc(dbc, bv);
  
  //we want to leave the matrix untouched for future projections
  //so we use a filtered matrix object
  FilteredMatrix<Vector<double> > filtered_A(A);
  filtered_A.add_constraints(bv);
  filtered_A.apply_constraints(rhs, true);

  deallog << "RHS L2 Norm After DBC: " << rhs.l2_norm() << std::endl;

  //make the preconditioner
  PreconditionJacobi<SparseMatrix<double> > precon;
  precon.initialize(A, 0.8);
  FilteredMatrix<Vector<double> > filtered_precon(precon);

  SolverControl control (vspace.n_dofs(),
			 tolerance*rhs.linfty_norm(),
			 false, true);

  GrowingVectorMemory<Vector<double> > mem;
  SolverCG<Vector<double> > solver (control, mem);

  //SparseDirectUMFPACK direct_solver;
  //direct_solver.initialize(filtered_A);

  solver.solve(filtered_A,sol_total,rhs,precon);

  vspace.get_hang().distribute(sol_total);

  //sol_total += sol_increment;

}


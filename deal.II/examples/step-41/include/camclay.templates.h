#include "../include/camclay.h"

template <int dim>
CamClay<dim>::CamClay()
 :
yield_stress(1)
{}

template <int dim>
CamClay<dim>::~CamClay()
{
  MM.clear();
}

template <int dim>
void CamClay<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Yield Stress");
  Functions::ParsedFunction<dim>::declare_parameters(prm, 1);
  prm.leave_subsection();

}

template <int dim>
void CamClay<dim>::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Yield Stress");
  yield_stress.parse_parameters(prm);
  prm.leave_subsection();
}

template <int dim>
void CamClay<dim>::reinit(VectorSpace<dim> &vspace)
{

  //the total number of quadrature points
  double total_qp = std::pow(double(2 * vspace.get_fe().degree + 1), dim);

  //Resize the tables to hold the plastic strains and other internal variables.
  plastic_strain.reinit(vspace.get_tria().n_active_cells(), int(total_qp));
  hardening.reinit(vspace.get_tria().n_active_cells(), int(total_qp));
  iter_table.reinit(vspace.get_tria().n_active_cells(), int(total_qp));

  solution.reinit(vspace.n_dofs());
  sol_hard_iter.reinit(vspace.n_dofs());

  iterations.reinit(vspace.get_tria().n_active_cells());
  iterations = 0; 

  //setup for projection mass matrix
  sp_MM.reinit(vspace.n_dofs(), 
	       vspace.n_dofs(), 
	       vspace.get_dh().max_couplings_between_dofs());

  DoFTools::make_sparsity_pattern (static_cast<DoFHandler<dim> &> (vspace.get_dh()), sp_MM);
  sp_MM.compress();

  MM.reinit(sp_MM);
}

template <int dim>
void CamClay<dim>::reinit_step()
{

  solution = 0;
  iterations = 0;
  sol_hard_iter = 0;

}

template <int dim>
void CamClay<dim>::initial_conditions(VectorSpace<dim> &vspace)
{

  QGauss<dim>   qf_v((2*vspace.get_fe().degree) + 1);
  
  FEValues<dim> fe_v (vspace.get_fe(), qf_v, 
		      update_values   | update_gradients |
		      update_quadrature_points | update_JxW_values);
  
  const unsigned int n_qp_v = qf_v.n_quadrature_points;
  

  typename MGDoFHandler<dim>::active_cell_iterator cell = vspace.get_dh().begin_active(),
    endc = vspace.get_dh().end();
  for (; cell!=endc; ++cell) {
    fe_v.reinit(cell);

    for (unsigned int qp = 0; qp<n_qp_v; ++qp) {

      hardening(cell->index(), qp) = 2e8;

      //plastic_strain(cell->index(), qp)[0][1] = 0;


    }
  }



}

template <int dim>
void CamClay<dim>::update_internal_variables(VectorSpace<dim> &vspace,
					     Vector<double> &elastic_solution,
                                             ParsedSymmetricTensorFunction<4, dim> &C)
{

  //first thing, we need to loop over the cells 
  //and then the quadrature points
  QGauss<dim>   qf_v((2*vspace.get_fe().degree) + 1);
  
  FEValues<dim> fe_v (vspace.get_fe(), qf_v, 
		      update_values   | update_gradients |
		      update_quadrature_points | update_JxW_values);
  
  const unsigned int n_qp_v = qf_v.n_quadrature_points;

  std::vector<Point<dim> > points(n_qp_v);
  
  SymmetricTensor<2,dim> trial_strain;
  SymmetricTensor<2,dim> trial_stress;
  SymmetricTensor<2,dim> total_strain;

  std::vector< std::vector<Tensor<1,dim> > > total_grads(n_qp_v, std::vector<Tensor<1,dim> >(dim));

  //average number of NR Iterations per cell 
  double nr_iterations_ave = 0;

  //This is the material parameter M, for the yield function/surface
  //right now, I just consider it constant
  M = 1.0;
  unsigned int counter = 0;
  bool marker = false;

  //loop only over the number of cells on this process
  //we will thread this later
  typename MGDoFHandler<dim>::active_cell_iterator cell = vspace.get_dh().begin_active(),
    endc = vspace.get_dh().end();
  for (; cell!=endc; ++cell) {
    //reinit with the correct cell
    fe_v.reinit(cell);
    int cell_index = cell->index();
    nr_iterations_ave = 0;
    marker = false;

    //get the gradients of the PURELY elastic solution
    fe_v.get_function_gradients(elastic_solution, total_grads);

    //get the quadrature points for the cell
    points = fe_v.get_quadrature_points();

    //loop over the quadrature points---------------------------------
    for (unsigned int qp = 0; qp<n_qp_v; ++qp) {

      //get the elastic moduli at this point
      C_qp = C(points[qp]);
     
      //calculate the infinitesimal strain
      for(unsigned int m=0; m<dim; ++m) {
	for(unsigned int n=0; n<dim; ++n) {
	  total_strain[m][n] = 0.5*(total_grads[qp][m][n]+total_grads[qp][n][m]);
	}
      }
  
      //get the trial strain and hardening.  Here we subtract the existing plastic
      //strain from the new total strain, giving the trial elastic strain
      //from this step.  The trial hardening is simply the exisiting hardening
      //parameter
      trial_strain = (total_strain - plastic_strain(cell_index, qp));
      //trial_strain *= -1.0;

      //calculate the trial stress.  This is of course just the contraction
      //of the elasticity tensor and the trial strain
      trial_stress = C_qp * trial_strain;
  
      //if the yield function is greater than zero, solve 
      //for the new values
      if ( yield_function(p(trial_stress), q(xi(trial_stress)), 
			  hardening(cell_index, qp) ) > 0 ) 
	{
	  marker = true;
	  //deallog << "Cell: " << cell_index
	  //	  << " stress01 value: " << trial_stress[0][1] << std::endl; 
	  
	   nr_iterations_ave += solve(cell_index,
	  			     qp,
	  			     total_strain);

       
	}
      
    } //quad points
    
    nr_iterations_ave /= n_qp_v;
    
    iterations(cell_index) = nr_iterations_ave;
    
    if (marker == true) ++counter;

  } //cells

  deallog << "Plastic Deformation in: "
	  << counter << " of " << vspace.get_tria().n_active_cells() << " cells." << std::endl;
  
}

template <int dim>
inline double CamClay<dim>::p(const SymmetricTensor<2, dim> &stress)
{

  //expicitly assumes specific 2D symmetry - NO SHEAR 
  double pval = (2*stress[0][0] + stress[1][1]);

  pval *= -1.0/3.0;

  return pval;

}

template <int dim>
inline SymmetricTensor<2,3> CamClay<dim>::xi(const SymmetricTensor<2,dim> &stress)
{

  SymmetricTensor<2,3> xival;

  xival[0][0] = stress[0][0];
  xival[1][1] = stress[0][0];
  xival[2][2] = stress[1][1];

  for(unsigned int i=0; i<3; ++i) xival[i][i] += p(stress);

  return xival;
}

template <int dim>
inline double CamClay<dim>::q(const SymmetricTensor<2, 3> &xi)
{

  double qval = std::sqrt(3.0/2.0) * xi.norm();

  return qval;

}

template <int dim>
double CamClay<dim>::solve(const int &index,
                           const unsigned int &qp,
                           const SymmetricTensor<2,dim> &total_strain)
{
  /*What we have is a system of non-linear ODE that need to
    be solved for the proper strains, plastic multiplier, and
    hardening coefficient.  The plastic multiplier is relevant only to
    this function, so we will only see it here*/
  
  /*This version is simpler - it is NOT using the invariant based
    methods of Claudio Tamagnini*/
  
  //these numbers are the 
  //stress components+hardening+plastic multiplier
  unsigned int size =( (dim*dim) + dim + 4 )/2;
  
  //we need vectors and a matrix
  Vector<double> cell_res(size);
  Vector<double> cell_sol(size);
  Vector<double> cell_prev(size);
  FullMatrix<double> cell_jac(size,size);
  
  //these will be needed later
  double first_norm = 0;
  unsigned int iters = 0;
  
  //create the guess, we start with the previously converged values of
  //the plastic strain and hardening, and zero for the plastic mult
  if (dim == 1) {
    cell_sol(0) = plastic_strain(index, qp)[0][0]; //e_p
    cell_sol(1) = hardening(index, qp); //k
    cell_sol(2) = 0; //delta gamma
  }
  if (dim == 2) { 
    cell_sol(0) = plastic_strain(index, qp)[0][0]; //e_p 11
    cell_sol(1) = plastic_strain(index, qp)[1][1]; //e_p 22
    cell_sol(2) = plastic_strain(index, qp)[0][1]; //e_p 12/21
    cell_sol(3) = hardening(index, qp); //k
    cell_sol(4) = 0; //delta gamma
  }
  if (dim == 3) {
    cell_sol(0) = plastic_strain(index, qp)[0][0]; //e_p xx
    cell_sol(1) = plastic_strain(index, qp)[1][1]; //e_p yy
    cell_sol(2) = plastic_strain(index, qp)[2][2]; //e_p zz
    cell_sol(3) = plastic_strain(index, qp)[1][2]; //e_p yz
    cell_sol(4) = plastic_strain(index, qp)[0][2]; //e_p xz
    cell_sol(5) = plastic_strain(index, qp)[0][1]; //e_p xy
    cell_sol(6) = hardening(index, qp); //k
    cell_sol(7) = 0; //delta gamma 
  }
  //since the guess of the solution is the previously converged values
  cell_prev = cell_sol;
  first_norm = 0;
  //loop over a newton raphson scheme.
  //for right now, max nr iterations is set to 10
  for(unsigned int n=0; n<10; ++n) {
    cell_res = 0;
    cell_jac = 0; 

    //compute the residual using the guess
    compute_res_jac(cell_res, cell_jac, 
		    cell_sol, cell_prev, total_strain); 
    
    cell_res *= -1.0;

    if ( (index == 20) && (qp == 0)) {
      deallog << "Cell Index: " << index
	     << " qp: " << qp << " Iters: " << iters 
	     << " Rnorm: " << cell_res.linfty_norm() << std::endl;
      
      // for (unsigned int i=0; i<size; ++i) 
      //deallog << "     Residual Comp: " << cell_res(i)  << std::endl;
    }

    
    //check for convergence, for right now just use 1e-10
    if(n == 0) first_norm = cell_res.linfty_norm();
    if((cell_res.linfty_norm()/first_norm) < 1e-13) break;
    
    //if we didn't converge, we are doing an iteration
    ++iters;


     //deallog << "Solution: " << cell_sol(size-2) << std::endl;

     // std::ofstream out("matrix.txt");    

    //cell_jac.print_formatted(out,3,true,0,"0",1,0);

    //invert the 5x5 full matrix
    cell_jac.gauss_jordan();
    
    //vmult the inverted matrix and add it to the solution
    cell_jac.vmult_add(cell_sol, cell_res);
    
  }

  //put the values into the proper places
  for(unsigned int i=0; i<dim; ++i) {
    for(unsigned int j=0; j<dim; ++j) {
      if (i >= j) plastic_strain(index, qp)[i][j] = cell_sol(sym2voigt(i,j));
    }
  }
  hardening(index,qp) = cell_sol(size-2);
  iter_table(index,qp) = iters;  

  return iters;
}

template <int dim>
void CamClay<dim>::compute_res_jac(Vector<double> &res,
				   FullMatrix<double> &jac,
				   const Vector<double> &sol,
				   const Vector<double> &prev_sol,
				   const SymmetricTensor<2,dim> &total_strain)
{
  //the critical function, computing the residual and the jacobian
  //basically, if we can assemble the residual using the sacado
  //doubles, then we have done all of the work.  We just need to be careful about
  //NOT exploiting the fact that the total strain is symmetric.

  typedef Sacado::Fad::DFad<double> fad_double;
  unsigned int size = sol.size();

  std::vector< std::vector<fad_double> > estrain(dim, std::vector<fad_double>(dim));

  //vectors of the unknowns and independent vars
  std::vector<fad_double> x(size); //this is ep11,ep22,ep12,k,dgamma
  std::vector<fad_double> off_diag(((dim*dim)-dim)/2);
  std::vector<fad_double> R(size); //this is ep11,ep22,ep12,k,dgamma

  for (unsigned int i=0; i<(2*(size-1))-dim; ++i) {
    //assign the values of the independent variables
    if( i<size) {
      x[i] = sol(i);
      x[i].diff(i,(2*(size-1))-dim);
    }

    if(i>=size) {
      off_diag[i-size] = sol(i-3);
      off_diag[i-size].diff(i,(2*(size-1))-dim);
    }

  }

  for(unsigned int i=0; i<dim; ++i) {
    for (unsigned int j=0; j<dim; ++j) {

      if(i >= j) estrain[i][j] =  (total_strain[i][j] - x[sym2voigt(i,j)]);

      if(i < j) estrain[i][j] = (total_strain[i][j] - off_diag[sym2voigt(i,j)]);

    }
  }

  //now we are set, we have estrain, with the right values, 
  //but with different sacado variables, and useable notation
  std::vector< std::vector<fad_double> > stress(dim, std::vector<fad_double>(dim));
  fad_double p;
  p=0;
  for(unsigned int i=0; i<dim; ++i) {
    for(unsigned int j=0; j<dim; ++j) {
      for (unsigned int m=0; m<dim; ++m) {
	for (unsigned int n=0; n<dim; ++n) {
	  
	  stress[i][j] += C_qp[i][j][m][n] * estrain[m][n];
	  
	}
      }

      if (i==j) p += stress[i][j];

    }
  }
  
  p *= (-1.0/double(dim));

  //deallog << "P value: " << p.val() <<std::endl;

  std::vector< std::vector<fad_double> > s(dim, std::vector<fad_double>(dim));
  fad_double q;
  fad_double trs;
  trs = 0;

  for (unsigned int i=0; i<dim; ++i) {
    for (unsigned int j=0; j<dim; ++j) {

      s[i][j] = stress[i][j];

      if (i == j) {

	s[i][j] += p;
	trs += s[i][j];
	
      }

       q += s[i][j]*s[i][j];

    }
  }

  q = std::sqrt(3.0*q/2.0);
  //deallog << " s00 value: " << estrain[0][0].val() << std::endl;
  //deallog << " s01 value: " << estrain[1][1].val() << std::endl;
  //deallog << " q value: " << q.val() << std::endl;
  //deallog << " dg value: " << x[size-1].val() << std::endl;

  for(unsigned int i=0; i<(size-2); ++i) {

    R[i] = x[i] - prev_sol(i) - ( x[size-1]*3.0*s[voigt2sym(i)[0]][voigt2sym(i)[1]]/M/M );

    //deallog << R[i] << std::endl;

    if(i<dim) R[i] -= x[size-1]*( ((x[size-2] - 2.0*p)/double(dim)) + (3.0*trs/M/M/double(dim)) ); 

  }

  //if (R[0].fastAccessDx(0) - R[1].fastAccessDx(1) > 1e-8) deallog << "Shit is Broke" <<std::endl;

  R[size-2] = x[size-2] - prev_sol(size-2) - (16.0*x[size-2]*x[size-1]*(2.0*p - x[size-2]));

  R[size-1] = (q*q/M/M) + (p*(p-x[size-2]));

  for (unsigned int i=0; i<size; ++i) {

    res(i) = R[i].val();
    //deallog << " R value: " << R[i] << std::endl;

    for (unsigned int j=0; j<size; ++j) {
      
      jac(i,j) = R[i].fastAccessDx(j);

      
    }//j
    
  }//i


}


template <int dim>
inline std::vector<unsigned int> CamClay<dim>::voigt2sym(const unsigned int i)
{

  std::vector<unsigned int> return_indices(2);

  if (dim == 2) {
    if (i == 0) { return_indices[0]=0; return_indices[1]=0; }
    if (i == 1) { return_indices[0]=1; return_indices[1]=1; }
    if (i == 2) { return_indices[0]=0; return_indices[1]=1; }
  }

  if (dim == 3) {
    //finish this at some time
  }

  return return_indices;

}

template <int dim>
inline double CamClay<dim>::delta(const unsigned int i,
				  const unsigned int j)
{
  double delta_val=0;

  if(i == j) delta_val = 1;

  return delta_val;

}

template <int dim>
inline double CamClay<dim>::delta(const unsigned int i)
{
  double delta_val=1;

  if(i >= dim) delta_val = 0;

  return delta_val;

}


template <int dim>
inline unsigned int CamClay<dim>::sym2voigt(const unsigned int i,
					    const unsigned int j)
{
  unsigned int voigt_value = 0;

  if (dim == 2) {

    if ( (i == 0) && (j == 0) ) voigt_value = 0;
    if ( (i == 1) && (j == 1) ) voigt_value = 1;
    if ( (i == 0) && (j == 1) ) voigt_value = 2;
    if ( (i == 1) && (j == 0) ) voigt_value = 2;


  }

  if (dim == 3) {
    //needs finished

  }

  return voigt_value;

}

template <int dim>
double CamClay<dim>::h(const double p,
                       const double k)
{
  
  double hval = k * ( 2*p - k ); 
  
  return hval;
  
}

template <int dim>
Vector<double> CamClay<dim>::dF_dstress(const SymmetricTensor<2, dim> &stress,
                                        const double k)
{
  SymmetricTensor<2,dim> tmp;
  Vector<double> tmpv(((dim*dim)+dim)/2);
  
  double p = (1.0/3.0) * first_invariant(stress);
  
  tmp = stress;
  for(unsigned int i=0; i<dim; ++i) tmp[i][i] -= p;
  
  double norm_xi = tmp.norm();
  
  double q = (std::sqrt(2.0/3.0)) * norm_xi;
  
  tmp = stress;
  
  tmp *= (q/M/M) * std::sqrt(8.0/3.0) * (1.0/norm_xi);
  
  for(unsigned int i=0; i<dim; ++i )
    tmp[i][i] += (2.0/3.0)*p - (1.0/3.0)*k - (q*p/M/M)*std::sqrt(8.0/3.0)*(1.0/norm_xi);
  
  tmpv(0) = tmp[0][0];
  tmpv(1) = tmp[1][1];
  tmpv(2) = tmp[0][1];
  
  return tmpv;
  
}

template <int dim>
inline double CamClay<dim>::yield_function(const double p, const double q,
					   const double k)
{
  double f_value = p*(p-k) + (q*q/M/M);

  return f_value;

}


template <int dim>
void CamClay<dim>::project_strain(ParameterHandler &prm,
				  VectorSpace<dim> &vspace,
				  SparseMatrix<double> &A)
{

  LocalAssemblePlasticProject<dim> local_plastic_project;
  local_plastic_project.reinit(vspace.get_fe(), plastic_strain);
  local_plastic_project.parameters(prm);

  MyTools::assemble_rhs(vspace.get_dh(), vspace.get_hang(),
			-1, solution, local_plastic_project);

  SparseDirectUMFPACK  direct_solver;

  direct_solver.initialize(A);

  direct_solver.solve(solution);

  vspace.get_hang().distribute(solution);

}

template <int dim>
void CamClay<dim>::build_matrix(VectorSpace<dim> &vspace)
{
  QGauss<dim> quad((2*vspace.get_fe().degree) + 1);

  MatrixCreator::create_mass_matrix(vspace.get_dh(), quad, MM);

}

template <int dim>
void CamClay<dim>::project_hardening(VectorSpace<dim> &vspace)
{
  //we need a mass matrix for a scalar variable over the mesh
  //so, lets build one from scratch

  LocalAssembleScalarProject<dim>  local_ass;
  local_ass.reinit(vspace.get_fe(), hardening, iter_table);

  MyTools::assemble_rhs(vspace.get_dh(), vspace.get_hang(),
			-1, sol_hard_iter, local_ass);

  SparseDirectUMFPACK  direct_solver;

  direct_solver.initialize(MM);

  direct_solver.solve(sol_hard_iter);

  vspace.get_hang().distribute(sol_hard_iter);

  //project the iterations while we are at it
  //local_ass.reinit(vspace_cc, iter_table);
  //MyTools::assemble_rhs(vspace_cc.get_dh(), vspace_cc.get_hang(),
  //			-1, sol_iterations, local_ass);



}


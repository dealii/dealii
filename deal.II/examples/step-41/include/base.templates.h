#include "../include/base.h"

template <int dim>
Base<dim>::Base ()
  :
  exact_solution(dim)
{}

template <int dim>
Base<dim>::~Base () 
{}

template <int dim>
void Base<dim>::run()
{
  //setup the depth of screen reporting
  deallog.depth_console (10);

  //get the parameters
  parameters ();

  //run as many convergence cycles as we need
  for(unsigned int cc=0; cc<num_cc; ++cc) {
    run_cc (cc);
  }

  //at the end of the problem, output the error table
  error_handler.output_table(0);
  
}

template <int dim>
void Base<dim>::run_cc(const unsigned int &cc)
{
  deallog.push("CC-"+Utilities::int_to_string(cc));
  deallog << "###############  Convergence Cycle: " 
	  << Utilities::int_to_string(cc) 
	  << " ###############" << std::endl;
  
  //for each convergence cycle, refine
  //for now, we just leave things global
  if(cc > 0) vspace.get_tria().refine_global(1);

  vspace.redistribute_dofs();
  
  elastic.reinit(vspace);
  
  plastic.reinit(vspace);

  plastic.initial_conditions(vspace);

  deallog << "Building Matrix..." << std::endl;
  elastic.build_matrix(prm, vspace);
  plastic.build_matrix(vspace);

  write_files (cc, 0,
	       elastic.sol_total,
	       plastic.solution,
	       plastic.sol_hard_iter);

  write_plot_values(cc, 0, elastic.sol_total, plastic.sol_hard_iter,
		    plastic.solution, elastic.C);

  //just loop over the loading steps
  for(unsigned int step=1; step<num_steps; ++step) {
    
    run_step(cc, step);
    
  }

  //we are interested in the errors, after the whole run
  exact_solution.set_time(end_time);
  
  error_handler.error_from_exact(vspace.get_dh(),
				 elastic.sol_total,
				 exact_solution,
				 0, 0);
  
  deallog.pop();

}

template <int dim>
void Base<dim>::run_step(const unsigned int &cc,
			 const unsigned int &step)
{
  deallog << std::endl;
  deallog.push("Step-" + Utilities::int_to_string(step,2));

  double time = end_time * (double(step)/double(num_steps));

  deallog << "Time: " << time <<std::endl;

  elastic_predictor (time);

  plastic_corrector ();

  write_files (cc, step,
	       elastic.sol_total,
	       plastic.solution,
	       plastic.sol_hard_iter);

  write_stress_strain (cc, step, elastic.sol_total);

  write_plot_values(cc, step, elastic.sol_total, plastic.sol_hard_iter,
		    plastic.solution, elastic.C);

  deallog.pop();
}

template <int dim>
void Base<dim>::elastic_predictor (double &time)
{
  deallog.push("ELASTIC");

  elastic.reinit_step(time);

  elastic.build_rhs(vspace);

  elastic.solve(vspace, lin_red_tol);

  deallog.pop();
}

template <int dim>
void Base<dim>::plastic_corrector ()
{
  deallog.push("PLASTIC");
  
  plastic.reinit_step();
  
  plastic.update_internal_variables(vspace,
				    elastic.sol_total,
                                    elastic.C);
  
  Timer timer;
  timer.start();
  plastic.project_strain(prm, vspace, elastic.A);
  plastic.project_hardening(vspace);
  timer.stop();
  deallog << "Time for plastic projection: " << timer() << std::endl;
  
  deallog.pop();
}


template<int dim>
void Base<dim>::write_files (const unsigned int &cc,
			     const unsigned int &step,
			     const Vector<double> &my_elas,
			     const Vector<double> &my_plas,
			     const Vector<double> &my_other)
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler (vspace.get_dh());

  std::vector<std::string> elastic_names;
  std::vector<std::string> plastic_names;
  std::vector<std::string> other_names;
  if(dim == 1) {
    elastic_names.push_back("u1_elastic");
    plastic_names.push_back("u1_plastic");
  } 

  if(dim == 2) {
    elastic_names.push_back("u1_elastic");
    elastic_names.push_back("u2_elastic");
    plastic_names.push_back("u1_plastic");
    plastic_names.push_back("u2_plastic");
    other_names.push_back("hardening");
    other_names.push_back("Newton_Iterations");
  } 

  if(dim == 3) {
    elastic_names.push_back("u1_elastic");
    elastic_names.push_back("u2_elastic");
    elastic_names.push_back("u3_elastic");
    plastic_names.push_back("u1_plastic");
    plastic_names.push_back("u2_plastic");
    plastic_names.push_back("u3_plastic");
    other_names.push_back("hardening");
    other_names.push_back("Newton_Iterations");
    other_names.push_back("Nothing");
  } 

  data_out.add_data_vector(my_elas, elastic_names);
  data_out.add_data_vector(my_plas, plastic_names);
  data_out.add_data_vector(my_other, other_names);
  data_out.add_data_vector(plastic.iterations, "Average_Newton_Iterations");

  std::vector<unsigned int> partition_int (vspace.get_tria().n_active_cells());
  GridTools::get_subdomain_association (vspace.get_tria(), partition_int);
  const Vector<double> partitioning(partition_int.begin(),
				    partition_int.end());
  data_out.add_data_vector (partitioning, "partitioning");

  data_out.build_patches (vspace.get_fe().degree);

  std::ostringstream filename;
  filename << "out/solution-";
  //first the convergence cycle
  filename << std::setfill('0');
  filename << std::setw(2) << cc;
  filename << "-";
  //then the step
  filename << std::setfill('0');
  filename << std::setw(5) << step;
  filename << ".vtk";

  std::ofstream output (filename.str().c_str());
  data_out.write_vtk (output);

}

template<int dim>
void Base<dim>::write_stress_strain (const unsigned int &cc,
				     const unsigned int &step,
				     const Vector<double> &my_elas)
{
  //reopen the just written file
  std::ostringstream filename;
  filename << "out/solution-";
  //first the convergence cycle
  filename << std::setfill('0');
  filename << std::setw(2) << cc;
  filename << "-";
  //then the step
  filename << std::setfill('0');
  filename << std::setw(5) << step;
  filename << ".vtk";
  std::ofstream output (filename.str().c_str(),std::ios::app);

  //WARNING - THE FOLLOWING IS NOT DIM INDEPENDENT!!!

  //create an fe object to get the gradient values
  QIterated<dim> qf(QTrapez<1>(),2);
  FEValues<dim> fe_v(vspace.get_fe(), qf,
		     update_gradients | update_quadrature_points);
  std::vector< std::vector< Tensor<1,dim> > > total_grads(qf.n_quadrature_points,
							  std::vector<Tensor<1,dim> >(dim));
  std::vector< Point<dim> > points(qf.n_quadrature_points);

  output << "SCALARS uyy double 1" << std::endl;
  output << "LOOKUP_TABLE default" << std::endl;

  typename MGDoFHandler<dim>::active_cell_iterator cell = vspace.get_dh().begin_active(),
    endc = vspace.get_dh().end();
  for (; cell!=endc; ++cell) {

    fe_v.reinit(cell);
    fe_v.get_function_gradients(my_elas, total_grads);

    points = fe_v.get_quadrature_points();

    //loop over the quadrature points---------------------------------
    for (unsigned int qp = 0; qp<qf.n_quadrature_points; ++qp) {

      output << total_grads[qp][1][1] << " ";

    }
  }


  //be responsible - close your files
  output.close();

}

template<int dim>
void Base<dim>::write_plot_values(const unsigned int &cc,
				  const unsigned int &step,
				  const Vector<double> &my_elas,
				  const Vector<double> &my_hardening,
				  const Vector<double> &my_plas,
				  const ParsedSymmetricTensorFunction<4,dim> &C)
{
  double p=0;
  double q=0;
  double deve=0;
  double k = 0;
  //std::vector< std::vector< double> > s(dim, std::vector<double>(dim));

  //reopen the just written file
  std::ostringstream filename;
  filename << "plot_data/data-";
  //first the convergence cycle
  filename << std::setfill('0');
  filename << std::setw(2) << cc;
  filename << ".txt";
  std::ofstream output;

  if(step == 0) {
   output.open(filename.str().c_str(),std::ios::trunc);
  } else {
   output.open(filename.str().c_str(),std::ios::app);
  }

  //create an fe object to get the gradient values at the nodes
  //QIterated<dim> qf(QTrapez<1>(),1);
  QTrapez<dim> qf;
  FEValues<dim> fe_v(vspace.get_fe(), qf,
		     update_values | update_gradients | update_quadrature_points);
  std::vector< std::vector< Tensor<1,dim> > > elastic_grads(qf.n_quadrature_points,
							    std::vector<Tensor<1,dim> >(dim));
  std::vector< std::vector< Tensor<1,dim> > > plastic_grads(qf.n_quadrature_points,
							    std::vector<Tensor<1,dim> >(dim));
  std::vector< Vector<double> > hard_values(qf.n_quadrature_points, Vector<double>(dim));
  std::vector< Point<dim> > points(qf.n_quadrature_points);

  //std::vector< std::vector<double> > strain(dim, std::vector<double>(dim) );

  //calculate p,q,deve
  typename MGDoFHandler<dim>::active_cell_iterator cell = vspace.get_dh().begin_active(),
    endc = vspace.get_dh().end();
  for (; cell!=endc; ++cell) {

    fe_v.reinit(cell);

    if (cell->index() == 20 ) {
      fe_v.get_function_gradients(my_elas, elastic_grads);
      fe_v.get_function_gradients(my_plas, plastic_grads);
      fe_v.get_function_values(my_hardening, hard_values);
      
      points = fe_v.get_quadrature_points();

      SymmetricTensor<2,dim> strain;
      SymmetricTensor<2,dim> stress;
      SymmetricTensor<2,3> s;

      for(unsigned int a=0; a<dim; ++a) {
	for(unsigned int b=0; b<dim; ++b) {
	  strain[a][b] = ( elastic_grads[0][a][b] - plastic_grads[0][a][b]);
	}
      }

      for(unsigned int a=0; a<dim; ++a) {
	for(unsigned int b=0; b<dim; ++b) {
	  for(unsigned int m=0; m<dim; ++m) {
	    for(unsigned int n=0; n<dim; ++n) {
	      stress[a][b] = C(points[0])[a][b][m][n] * strain[m][n];
	    }
	  }
	}
      }
      //explicitly assuming the loading/symmetry condition

      p = (-1.0/3.0) * (2*stress[0][0] + stress[1][1]);

      s[0][0] = stress[0][0];
      s[1][1] = stress[0][0];
      s[2][2] = stress[1][1];

      for(unsigned int kk=0; kk<3; ++kk) s[kk][kk] += p;

      q = std::sqrt(3.0/2.0) * s.norm();

      k = hard_values[0](0);

    }
  }

  //write what we want
  output << step << " "
	 << p << " "
	 << q << " "
	 << deve << " " 
	 << k << " " << std::endl;


  //be responsible - close your files
  output.close();

}

template<int dim>
void Base<dim>::parameters ()
{
  deallog.push("PARAMETERS");

  domain.declare_parameters(prm);

  //add settings for the elastic vector space to the prm file
  std::string space_name = "Vector Space Parameters";
  vspace.declare_parameters(prm, space_name);

  error_handler.declare_parameters(prm);

  elastic.declare_parameters(prm);

  plastic.declare_parameters(prm);
  
  prm.enter_subsection("General Parameters");
  prm.declare_entry ("Linear Solver Reducation Tolerance",
		     "1.0E-16",
		     Patterns::Double(),
		     "Linear Solver Reduction Tolerance");
  prm.declare_entry ("Number of Convergence Cycles",
		     "1",
		     Patterns::Integer(),
		     "Number of Convergence Cycles Past Initial Refinement");
  prm.declare_entry ("Console Depth",
		     "10",
		     Patterns::Integer(),
		     "Determines level of screen output");
  prm.leave_subsection();

  prm.enter_subsection("Loading Parameters");
  prm.declare_entry ("End Time for Simulation",
		     "1.0",
		     Patterns::Double(),
		     "End time for the simulation");
  prm.declare_entry ("Number of Loading Steps",
		     "10",
		     Patterns::Integer(),
		     "Number of Loading Cycles to split time interval");
  prm.leave_subsection();

  prm.enter_subsection("Exact Solution");
  Functions::ParsedFunction<dim>::declare_parameters(prm, dim);
  prm.leave_subsection();
  
  // ==============================
  
  prm.read_input("CamClay.prm");
  
  // ==============================
  //Initialize the domain    
  domain.reinit(prm);
  /*
  domain.get_tria().clear();

  Point<dim> lowerleft(0.0,0.0);
  Point<dim> upperright(0.3,0.1);
  std::vector<unsigned int> subs(2);
  subs[0] = 30; subs[1]=10;
  GridGenerator::subdivided_hyper_rectangle(domain.get_tria(),
						 subs,
						 lowerleft,
						 upperright, true);
  */  
  vspace.reinit(prm, domain.get_tria());
  
  error_handler.parse_parameters(prm);

  elastic.parse_parameters(prm);

  plastic.parse_parameters(prm);
  
  prm.enter_subsection("Exact Solution");
  exact_solution.parse_parameters(prm);
  prm.leave_subsection();
  
  prm.enter_subsection("General Parameters");
  lin_red_tol = prm.get_double("Linear Solver Reducation Tolerance");
  num_cc = prm.get_integer("Number of Convergence Cycles");
  console_depth = prm.get_integer("Console Depth");
  prm.leave_subsection();

  prm.enter_subsection("Loading Parameters");
  num_steps = prm.get_integer("Number of Loading Steps");
  end_time = prm.get_double("End Time for Simulation");
  prm.leave_subsection();

  deallog.pop();

}



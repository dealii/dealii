/*
    Implementation of the Navier_Stokes_projection class
    by Abner Salgado.
*/
#include "../include/NavierStokes.h"



#include <sys/time.h>

/// debug
bool show;
double Solve_time;
///----


// Constructor
template<int dim> Navier_Stokes_Projection<dim>::Navier_Stokes_Projection( const Data_Storage &data):
                    type( data.form ), deg( data.pressure_degree ), dt( data.initial_dt ), t_0( data.initial_time ),
                    T( data.final_time ), Re( data.Reynolds ), rhs( data.initial_time ), vel_exact( data.initial_time ),
                    dof_handler_velocity(triangulation), dof_handler_pressure(triangulation),
                    fe_velocity(deg+1), fe_pressure(deg), quadrature_pressure(deg+1), quadrature_velocity(deg+2),
                    vel_max_its( data.vel_max_iterations ), vel_Krylov_size( data.vel_Krylov_size ),
                    vel_off_diagonals( data.vel_off_diagonals ),
                    vel_update_prec( data.vel_update_prec ), vel_eps( data.vel_eps ), vel_diag_strength( data.vel_diag_strength),
                    proj_max_its( data.proj_max_iterations ), proj_off_diagonals( data.proj_off_diagonals ),
                    proj_eps( data.proj_eps ), proj_diag_strength( data.proj_diag_strength ),
                    pres_max_its( data.pres_max_iterations), pres_off_diagonals( data.pres_off_diagonals ),
                    pres_eps( data.pres_eps ), pres_diag_strength( data.pres_diag_strength )
{
  // After having initialized the bunch of data we do nothing
  // NOTE TO SELF: Do I need to do this check?
  if(deg < 1)
  std::cout<<" WARNING: The chosen pair of finite element spaces is not stable."<<std::endl
           <<" The obtained results will be nonsense"<<std::endl;
  // We check that the time step is within the allowed limits
  AssertThrow( not ( ( dt <= 0. ) or ( dt > .5*T ) ), ExcInvalidTimeStep( dt, .5*T ) );
}



// Destructor
template<int dim> Navier_Stokes_Projection<dim>::~Navier_Stokes_Projection(){
  dof_handler_velocity.clear();
  dof_handler_pressure.clear();
}



// Set time step
template<int dim> void Navier_Stokes_Projection<dim>::set_dt( const double ddt ){
  // We just check that it is within the permitted limits
  AssertThrow( not ( ( ddt <= 0. ) or ( ddt > .5*T ) ), ExcInvalidTimeStep( ddt, .5*T ) );
  dt = ddt;
}



// Initialization of the velocity matrices and assembly of those that do not depend on dt
template<int dim> void Navier_Stokes_Projection<dim>::init_velocity_matrices(){
  //// Init the sparsity pattern for the velocity
  spar_pattern_velocity.reinit( dof_handler_velocity.n_dofs(), dof_handler_velocity.n_dofs(),
                                dof_handler_velocity.max_couplings_between_dofs() );
  DoFTools::make_sparsity_pattern( dof_handler_velocity, spar_pattern_velocity );
  spar_pattern_velocity.compress();

  //// Init the matrices for the velocity
  vel_Laplace_plus_Mass.reinit( spar_pattern_velocity );
  for( unsigned int d=0; d<dim; ++d )
    vel_it_matrix[d].reinit( spar_pattern_velocity );
  vel_Mass.reinit( spar_pattern_velocity );
  vel_Laplace.reinit( spar_pattern_velocity );

  // We assemble the mass matrix
  MatrixCreator::create_mass_matrix( dof_handler_velocity, quadrature_velocity, vel_Mass );

  // We assemble the Laplace matrix
  MatrixCreator::create_laplace_matrix( dof_handler_velocity, quadrature_velocity, vel_Laplace );
}



template<int dim> void Navier_Stokes_Projection<dim>::init_pressure_matrices(){
  //// Init the sparsity pattern for the pressure
  spar_pattern_pressure.reinit( dof_handler_pressure.n_dofs(), dof_handler_pressure.n_dofs(),
                                dof_handler_pressure.max_couplings_between_dofs() );
  DoFTools::make_sparsity_pattern( dof_handler_pressure, spar_pattern_pressure );

  // Before we close the sparsity pattern, we need to
  // init the constraints for the Laplace operator on the pressure space
  pres_regularization.clear();
//  DoFTools::make_hanging_node_constraints( dof_handler_pressure, pres_regularization ); //??

  // Add the only constraint that we have
  pres_regularization.add_line(0);

  // close it
  pres_regularization.close();

  // condense the sparsity pattern
  pres_regularization.condense( spar_pattern_pressure );

  // Compress the sparsity pattern for the pressure
  spar_pattern_pressure.compress();

  //// Init the matrices for the pressure
  pres_Laplace.reinit( spar_pattern_pressure );
  pres_Mass.reinit( spar_pattern_pressure );

  // Now we assemble the matrices
  // The Laplace operator is the projection matrix
  MatrixCreator::create_laplace_matrix( dof_handler_pressure, quadrature_pressure, pres_Laplace );
  // The pressure mass matrix to do the pressure update
  MatrixCreator::create_mass_matrix(  dof_handler_pressure, quadrature_pressure, pres_Mass );

  // Finally we condense the Laplace operator
  pres_regularization.condense( pres_Laplace );
}



template<int dim> void Navier_Stokes_Projection<dim>::init_gradient_operator(){
  //// Init the sparsity pattern for the gradient operator
  spar_pattern_pres_vel.reinit( dof_handler_velocity.n_dofs(), dof_handler_pressure.n_dofs(),
                                dof_handler_velocity.max_couplings_between_dofs() );
  DoFTools::make_sparsity_pattern( dof_handler_velocity, dof_handler_pressure, spar_pattern_pres_vel );
  spar_pattern_pres_vel.compress();

  /*
      To assemble each component of the gradient operator
      we need to make a loop over all cells and compute the products of
          velocity * d_#(pressure)
      where # is the current space component. For this reason we need two cell
      iterators, one for the velocity and one for the pressure space.
  */
  typename DoFHandler<dim>::active_cell_iterator cell_init       = dof_handler_velocity.begin_active(),
                                                   cell_end        = dof_handler_velocity.end(),
                                                   cell,
                                                   bogus_cell_init = dof_handler_pressure.begin_active(),
                                                   bogus_cell;
  /*
    The FEValues extractors.
    For the velocity we need values,
    For the pressure we only need gradients.
  */
  FEValues<dim> fe_values_velocity( fe_velocity, quadrature_velocity, update_values | update_JxW_values ),
                fe_values_pressure( fe_pressure, quadrature_velocity, update_gradients );

  // Usual, and useful, abreviations
  const unsigned int vel_dofs_per_cell  = fe_velocity.dofs_per_cell,
                     pres_dofs_per_cell = fe_pressure.dofs_per_cell,
                     n_q_points         = quadrature_velocity.size();

  // The local gradient operator
  FullMatrix<double> local_grad( vel_dofs_per_cell, pres_dofs_per_cell );

  // Local to global DoF's map
  std::vector<unsigned int> vel_local_dof_indices( vel_dofs_per_cell ), pres_local_dof_indices( pres_dofs_per_cell );

  for( unsigned int d=0; d<dim; ++d ){
    // we reinit the given component of the gradient operator
    pres_Diff[d].reinit( spar_pattern_pres_vel );

    // The usual cell loop
    for( cell = cell_init, bogus_cell = bogus_cell_init; cell not_eq cell_end; ++cell, ++bogus_cell ){
      fe_values_velocity.reinit( cell );
      fe_values_pressure.reinit( bogus_cell );

      cell->get_dof_indices( vel_local_dof_indices );
      bogus_cell->get_dof_indices( pres_local_dof_indices );

      // local contributions
      local_grad = 0.;
      for( unsigned int q=0; q<n_q_points; ++q ){
        for( unsigned int i=0; i<vel_dofs_per_cell; ++i )
          for( unsigned int j=0; j<pres_dofs_per_cell; ++j )
            local_grad( i, j ) += fe_values_velocity.JxW(q)*fe_values_velocity.shape_value( i, q )
                                  *fe_values_pressure.shape_grad( j, q)[d];
      }

      // Add it to the global contributions
      for( unsigned int i=0; i<vel_dofs_per_cell; ++i )
        for( unsigned int j=0; j<pres_dofs_per_cell; ++j)
          pres_Diff[d].add( vel_local_dof_indices[i], pres_local_dof_indices[j], local_grad( i, j)  );
    }
  }
}



// Create Triangulation method
template<int dim> void Navier_Stokes_Projection<dim>::Create_Triangulation( const unsigned int n_of_refines ){
  // A disk
  GridGenerator::hyper_ball( triangulation );
  static const HyperBallBoundary<dim> boundary;
  triangulation.set_boundary( 0, boundary );

/*  // Our domain is a unit square.
  GridGenerator::hyper_cube(triangulation);*/
  triangulation.refine_global( n_of_refines );
  std::cout<<" Number of active cells: "<<triangulation.n_active_cells()<<std::endl;

  // Having created the mesh, we create the mesh-dependent data

  //// Distribute the degrees of freedom and renumber them
  dof_handler_velocity.distribute_dofs( fe_velocity );
  DoFRenumbering::boost::Cuthill_McKee( dof_handler_velocity );
  dof_handler_pressure.distribute_dofs( fe_pressure );
  DoFRenumbering::boost::Cuthill_McKee( dof_handler_pressure );

  // Initialize the matrices for the velocity
  init_velocity_matrices();

  // Initialize the matrices for the pressure
  init_pressure_matrices();

  // Initialize the gradient operator
  init_gradient_operator();


  //// Set the correct size for the various vectors that we need
  pres_n.reinit( dof_handler_pressure.n_dofs() );
  pres_n_minus_1.reinit( dof_handler_pressure.n_dofs() );

  phi_n.reinit( dof_handler_pressure.n_dofs() );
  phi_n_minus_1.reinit( dof_handler_pressure.n_dofs() );

  pres_tmp.reinit( dof_handler_pressure.n_dofs() );

  for(unsigned int d=0; d<dim; ++d){
    u_n[d].reinit( dof_handler_velocity.n_dofs() );
    u_n_minus_1[d].reinit( dof_handler_velocity.n_dofs() );
    u_star[d].reinit( dof_handler_velocity.n_dofs() );
    force[d].reinit( dof_handler_velocity.n_dofs() );
  }
  v_tmp.reinit( dof_handler_velocity.n_dofs() );

  std::cout<<" dim( X_h ) = "<<( dof_handler_velocity.n_dofs()*dim )<<std::endl
           <<" dim( M_h ) = "<<dof_handler_pressure.n_dofs()<<std::endl;
}



// Initialization of the time dependent data
template<int dim> void Navier_Stokes_Projection<dim>::Initialize(){

  // The constant part of the matrix for the velocity is 1.5/dt*Mass + 1/Re*Delta_h
  vel_Laplace_plus_Mass = 0.;
  vel_Laplace_plus_Mass.add( 1./Re, vel_Laplace );
  vel_Laplace_plus_Mass.add( 1.5/dt, vel_Mass );

  // Now we start the vectors and load them with the initial data
  ////pressure
  Pressure<dim> pres( t_0 );
  VectorTools::interpolate( dof_handler_pressure, pres, pres_n_minus_1 );
  pres.advance_time( dt );
  VectorTools::interpolate( dof_handler_pressure, pres, pres_n );

  ////phi
  phi_n = 0.;
  phi_n_minus_1 = 0.;

  ////velocity
  Velocity<dim> v;
  for(unsigned int d=0; d<dim; ++d){
    v.set_time( t_0 );
    v.set_component(d);

    VectorTools::interpolate( dof_handler_velocity, v, u_n_minus_1[d] );

    v.advance_time( dt );
    VectorTools::interpolate( dof_handler_velocity, v, u_n[d] );
  }
  plot_solution(1);
}



/*
    A time marching step, the output messages and name of the functions are self
    explanatory, so no comments here
*/
template<int dim> void Navier_Stokes_Projection<dim>::run( const bool verbose, const unsigned int n_of_plots ){
  // We set the number of steps
  unsigned int n_steps = ( T - t_0 )/dt;
  // Set the source term to the correct time
  rhs.set_time( 2.*dt );
  vel_exact.set_time( 2.*dt );
  // do as many steps as necessary
  timeval init_time, end_time;

  /// debug
  show = false;
  Solve_time = 0.;
  ///-----

  gettimeofday( &init_time, 0 );
  for( unsigned int n = 2; n<=n_steps; ++n ){

    /// debug
    if( n == n_steps )
      show = true;
    ///----

    if( verbose )
      std::cout<<" Step = "<<n<<" Time = "<<(n*dt)<<std::endl;
    if( verbose )
      std::cout<<"  Interpolating the velocity "<<std::endl;
    interpolate_velocity();
    if( verbose )
      std::cout<<"  Diffusion Step"<<std::endl;
    if( ( n%vel_update_prec==0 ) and ( verbose ) )
      std::cout<<"   With reinitialization of the preconditioner"<<std::endl;
    diffusion_step( (n%vel_update_prec == 0 ) or ( n == 2) );
    if( verbose )
      std::cout<<"  Projection Step"<<std::endl;
    projection_step( ( n == 2 ) );
    if( verbose )
      std::cout<<"  Updating the Pressure"<<std::endl;
    update_pressure( ( n == 2 ) );
    if( n%n_of_plots ==0 ){
      if( verbose )
        std::cout<<" Plotting Solution"<<std::endl;
      plot_solution(n);
    }
    //  advance in time
    rhs.advance_time(dt);
    vel_exact.advance_time(dt);
  }
  gettimeofday( &end_time, 0 );
  double run_time = double( end_time.tv_sec - init_time.tv_sec ) + 1e-6*double( end_time.tv_usec - init_time.tv_usec );
  run_time /= n_steps;
  run_time /= double( dim*dof_handler_velocity.n_dofs() + dof_handler_pressure.n_dofs() );
//   if( verbose )
    std::cout<<" Run time = "<<run_time<<std::endl;
}



// Simple procedure that interpolates the velocity
template<int dim> void Navier_Stokes_Projection<dim>::interpolate_velocity(){
  for( unsigned int d=0; d<dim; ++d ){
    u_star[d] = 0.;
    u_star[d].add( 2., u_n[d], -1., u_n_minus_1[d] );
  }
}



template<int dim> void Navier_Stokes_Projection<dim>::assemble_advection_term( const unsigned int d ){
  // Cell iterator
  typename DoFHandler<dim>::active_cell_iterator cell     = dof_handler_velocity.begin_active(),
                                                   cell_end = dof_handler_velocity.end();

  /*
    The FEValues extractors.
    We need values because of the terms involving u^*
    also we need the gradient to assemble the advection part
  */
  FEValues<dim> fe_values( fe_velocity, quadrature_velocity, update_values | update_JxW_values | update_gradients );

  // Usual, and useful, abreviations
  const unsigned int dofs_per_cell = fe_velocity.dofs_per_cell,
                     n_q_points    = quadrature_velocity.size();

  /*
    We need to compute the values of several finite element functions (and some gradients)
    at quadrature points. For this we use all these arrays.
    The names are self-explanatory
  */
  std::vector<double> u_old_local( n_q_points );
  std::vector< Point<dim> > u_star_local( n_q_points );
  std::vector< Tensor<1,dim> > grad_u_star(n_q_points);

  // Local to global DoF's map
  std::vector<unsigned int> local_dof_indices( dofs_per_cell );

  // Local contribution of the advection operator
  Vector<double> local_rhs( dofs_per_cell );

  // loop over cells
  for( ; cell not_eq cell_end; ++cell ){
    // reinit the stuff
    fe_values.reinit( cell );
    cell->get_dof_indices( local_dof_indices );

    // create the local u^*
    for( unsigned int s=0; s<dim; ++s ){
      fe_values.get_function_values( u_star[s], u_old_local );
      for( unsigned int q=0; q<n_q_points; ++q )
        u_star_local[q](s) = u_old_local[q];
    }

    fe_values.get_function_gradients( u_star[d], grad_u_star );

    local_rhs = 0.;
    for( unsigned int q =0; q<n_q_points; ++q )
      for( unsigned int i=0; i<dofs_per_cell; ++i )
        // The local contribution of the advection term
        local_rhs(i) -= ( ( u_star_local[q]*grad_u_star[q] )*fe_values.shape_value( i, q ) )*fe_values.JxW(q) ;

    // Add local contributions to the global matrix
    for( unsigned int i=0; i<dofs_per_cell; ++i )
      force[d]( local_dof_indices[i] ) += local_rhs(i);
  }
}




template<int dim> void Navier_Stokes_Projection<dim>::threaded_solve( const unsigned int d){

  // Create the solver that we are going to use
  SolverControl solver_control( vel_max_its, vel_eps*force[d].l2_norm() );
  SolverCG<> cg( solver_control );
  // Solve
  cg.solve( vel_it_matrix[d], u_n[d], force[d], prec_velocity[d] );
}



// A diffusion step
template<int dim> void Navier_Stokes_Projection<dim>::diffusion_step( const bool reinit_prec ){
  // Used to solve in parallel the components of the velocity
  Threads::ThreadGroup<> threads;

  /// debug
  static double assembly_time = 0.;
  timeval init_time_assembly, end_time_assembly, init_time_solve, end_time_solve;
  gettimeofday( &init_time_assembly, 0 );
  ///----

  /*
    Define the pressure interpolant. It is defined by
        p^* = p^{n} + 4/3 phi^{n} - 1/3 phi^{n-1}
  */
  pres_tmp = pres_n;
  pres_tmp.add(4./3., phi_n, -1./3., phi_n_minus_1);
  pres_tmp *= -1.;

  for( unsigned int d=0; d<dim; ++d ){
    // First we need to construct the right hand sides
    rhs.set_component(d);
    //// compute the source term
    VectorTools::create_right_hand_side( dof_handler_velocity, quadrature_velocity, rhs, force[d] );

    // We assemble the advection term
    threads += Threads::spawn( *this, &Navier_Stokes_Projection<dim>::assemble_advection_term )(d);
    //assemble_advection_term(d);

    /*
      Compute the value of the contribution of the old velocities.
      Recall that BDF2 reads
      ( 0.5/dt)*( 3u^{n+1} - 4u^{n} + u^{n-1} )
      the terms involving u^{n} and u^{n-1} are known, so they go into the
      right hand side
    */
    v_tmp = 0.;
    v_tmp.add( 2./dt,u_n[d],-.5/dt,u_n_minus_1[d] );
    vel_Mass.vmult_add( force[d], v_tmp );

    // And add the contribution of the pressure part
    pres_Diff[d].vmult_add( force[d], pres_tmp );

    // Copy the old value
    u_n_minus_1[d] = u_n[d];


    // Create the matrix
    vel_it_matrix[d].copy_from( vel_Laplace_plus_Mass );
    //// Create the boundary values
    vel_exact.set_component(d);
    VectorTools::interpolate_boundary_values( dof_handler_velocity, 0, vel_exact, boundary_values );

    // apply BCs
    MatrixTools::apply_boundary_values( boundary_values, vel_it_matrix[d], u_n[d], force[d] );
  }
  threads.join_all();

  ///  debug
  gettimeofday( &end_time_assembly, 0 );
  assembly_time += double( end_time_assembly.tv_sec - init_time_assembly.tv_sec )
                            + 1e-6*double( end_time_assembly.tv_usec - init_time_assembly.tv_usec );
  gettimeofday( &init_time_solve, 0 );
  ///-----

  for(unsigned int d=0; d<dim; ++d ){
    // Check if the preconditioner needs reinitialization
  //  if( reinit_prec )
      prec_velocity[d].initialize( vel_it_matrix[d]/*, SparseILU<double>::AdditionalData( vel_diag_strength, vel_off_diagonals )*/ );

    //prec_velocity[d].factorize( vel_it_matrix[d] );
//    threaded_solve(d);

//     And solve
//    threads += Threads::spawn( *this, &Navier_Stokes_Projection<dim>::threaded_solve )(d);

//     prec_velocity[d].solve( force[d] );
//     u_n[d] = force[d];
    prec_velocity[d].vmult( force[d], u_n[d] );
  }
//   threads.join_all();


  ///debug
  gettimeofday( &end_time_solve, 0 );
  Solve_time += double( end_time_solve.tv_sec - init_time_solve.tv_sec )
                            + 1e-6*double( end_time_solve.tv_usec - init_time_solve.tv_usec );
  if(show)
    std::cout<<" Assembly time = "<<assembly_time<<std::endl
             <<" Solve time = "<<Solve_time<<std::endl;
  ///----
}



// A projection step
template<int dim> void Navier_Stokes_Projection<dim>::projection_step( const bool reinit_prec ){

  /// debug
  timeval init_time_assemble, init_time_solve, end_time_assemble, end_time_solve, init_time_factorize, end_time_factorize;
  static double assemble_time=0., solve_time=0., factorize_time = 0.;
  ///----


  // Check if the preconditioner needs reinitialization
  if( reinit_prec ){

    /// debug
    gettimeofday( &init_time_factorize, 0 );
    ///-----
// prec_pressure.initialize( pres_Laplace, SparseILU<double>::AdditionalData( proj_diag_strength, proj_off_diagonals ) );
    prec_pressure.initialize( pres_Laplace );

    /// debug
    gettimeofday( &end_time_factorize, 0 );
    factorize_time += double( end_time_factorize.tv_sec - init_time_factorize.tv_sec )
                  + 1e-6*double( end_time_factorize.tv_usec - init_time_factorize.tv_usec );
    ///--------
  }

  /// debug
  gettimeofday( &init_time_assemble, 0 );
  ///-----

  // To make things faster we assembled the discrete gradient operator
  // so we just need to multiply by its transpose
  pres_tmp = 0.;
  for( unsigned d=0; d<dim; ++d )
    pres_Diff[d].Tvmult_add( pres_tmp, u_n[d] );

  /// debug
  gettimeofday( &end_time_assemble, 0 );
  assemble_time += double( end_time_assemble.tv_sec - init_time_assemble.tv_sec )
                 + 1e-6*double( end_time_assemble.tv_usec - init_time_assemble.tv_usec );
  ///-----

  // Up to this point, pres_tmp has the term < u_n, grad q >
  // we will assemble this only once, so we
  // Copy the old value of phi
  phi_n_minus_1 = phi_n;

  phi_n = pres_tmp;
  phi_n *= 1.5/dt;

  /// debug
  gettimeofday( &init_time_solve, 0 );
  ///-------

  // condense because of the constraint
  pres_regularization.condense( phi_n );


  prec_pressure.solve( phi_n );
  // The method that we are going to use
//   SolverControl solver_control( proj_max_its, proj_eps*pres_tmp.l2_norm() );
//   SolverCG<> cg( solver_control );


  // And solve
//   cg.solve( pres_Laplace, phi_n, pres_tmp, prec_pressure );

  // After solving we need to distribute the constrained stuff
  pres_regularization.distribute( phi_n );

  /// debug
  gettimeofday( &end_time_solve, 0 );
  solve_time += double( end_time_solve.tv_sec - init_time_solve.tv_sec )
                 + 1e-6*double( end_time_solve.tv_usec - init_time_solve.tv_usec );
  if(show)
    std::cout<<" DIV assembly time = "<<assemble_time<<std::endl
             <<" LAP solve time = "<<solve_time<<std::endl
             <<" LAP factorize time = "<<factorize_time<<std::endl;
  ///------
}



// Pressure update step
template<int dim> void Navier_Stokes_Projection<dim>::update_pressure( const bool reinit_prec ){

  /// debug
  static double solve_time = 0.;
  timeval init_time_solve, end_time_solve;
  ///---------

  // Copy the old pressure
  pres_n_minus_1 = pres_n;
  switch( type ){
    case METHOD_STANDARD:
      /*
        Simple standard case
            p^{n+1} = p^n + \phi^{n+1}
      */
      pres_n += phi_n;
      break;
    case METHOD_ROTATIONAL:
    {
      /*
        Rotational form
            p^{n+1} = p^n + \phi^{n+1} - 1/Re div u^{n+1}
      */

      /// debug
      gettimeofday( &init_time_solve, 0 );
      ///------

      // The preconditioner for the mass matrix problem
      if( reinit_prec )
        prec_mass.initialize( pres_Mass/*, SparseILU<double>::AdditionalData( pres_diag_strength, pres_off_diagonals )*/ );

      // The solver for  the mass matrix
//       SolverControl solver_control( pres_max_its, pres_eps*pres_tmp.l2_norm() );
//       SolverCG<> cg( solver_control );
//
//       // Solve the mass matrix problem
//       cg.solve( pres_Mass, pres_n, pres_tmp, prec_mass );

      pres_n = pres_tmp;
      prec_mass.solve( pres_n );

      pres_n.sadd(1./Re, 1., pres_n_minus_1, 1., phi_n );

      /// debug
      gettimeofday( &end_time_solve, 0 );
      solve_time += double( end_time_solve.tv_sec - init_time_solve.tv_sec )
                 + 1e-6*double( end_time_solve.tv_usec - init_time_solve.tv_usec );
      ///----------

      break;
    }
    default:
      Assert( false, ExcNotImplemented() );
  };

  /// debug
  if(show)
    std::cout<<" PRES_UPD solve time = "<<solve_time<<std::endl;
  ///---

}



// Post Processing
template<int dim> void Navier_Stokes_Projection<dim>::Post_Process(){
  double tmp, vel_err_L2=0., vel_err_H1=0., pres_err_L2;

  Vector<double> differences( triangulation.n_active_cells() );

  Velocity<dim> u_exact(T);
  for( unsigned int d=0; d<dim; ++d ){
    u_exact.set_component(d);

    // L2 error for this component of the velocity
    differences = 0.;
    VectorTools::integrate_difference( dof_handler_velocity, u_n[d], u_exact, differences, quadrature_velocity, VectorTools::L2_norm );
    tmp = differences.l2_norm();
    vel_err_L2 += tmp*tmp;

    // H1 error for this component of the velocity
    differences = 0.;
    VectorTools::integrate_difference( dof_handler_velocity, u_n[d], u_exact, differences, quadrature_velocity, VectorTools::H1_seminorm );
    tmp = differences.l2_norm();
    vel_err_H1 += tmp*tmp;
  }
  vel_err_L2 = std::sqrt( vel_err_L2 );
  vel_err_H1 = std::sqrt( vel_err_H1 );


  /*
    The actual pressure needs to have mean value zero.
    We compute the mean value and subtract it from the computed pressure
  */
  double pres_mean_value = VectorTools::compute_mean_value( dof_handler_pressure, quadrature_pressure, pres_n, 0 );
  for( unsigned int i=0; i<pres_n.size(); ++i )
    pres_n(i) -= pres_mean_value;

  Pressure<dim> pres_exact(T);

  // L2 error of the pressure
  differences = 0.;
  VectorTools::integrate_difference( dof_handler_pressure, pres_n, pres_exact, differences, quadrature_pressure, VectorTools::L2_norm );
  pres_err_L2 = differences.l2_norm();

  // Add the obtained results to the convergence table
  convergence_table.add_value( "dt"     , dt );
  convergence_table.add_value( "u_L2"   , vel_err_L2 );
  convergence_table.add_value( "u_H1"   , vel_err_H1 );
  convergence_table.add_value( "pres_L2", pres_err_L2 );

  convergence_table.set_precision( "dt"     , 5 );
  convergence_table.set_precision( "u_L2"   , 5 );
  convergence_table.set_precision( "u_H1"   , 5 );
  convergence_table.set_precision( "pres_L2", 5 );

  convergence_table.set_scientific( "u_L2"   , true );
  convergence_table.set_scientific( "u_H1"   , true );
  convergence_table.set_scientific( "pres_L2", true );

  // And show it
  convergence_table.write_text(std::cout);
}



// Output of the solution
template<int dim> void Navier_Stokes_Projection<dim>::plot_solution( const unsigned int step ){
  // This only works in 2d
  Assert( dim==2, ExcNotImplemented() );

  // We need to assemble the vorticy
  FEValues<dim> fe_values( fe_velocity, quadrature_velocity, update_values | update_gradients | update_JxW_values );
  const unsigned int dofs_per_cell = fe_velocity.n_dofs_per_cell(),
                     n_q_points    = quadrature_velocity.size();
  std::vector< Tensor<1,dim> > grad_vel_1( n_q_points ), grad_vel_2( n_q_points );
  std::vector< unsigned int> local_dof_indices( dofs_per_cell );
  Vector<double> local_rhs( dofs_per_cell );
  double vorticity;

  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_velocity.begin_active(),
                                                   cend = dof_handler_velocity.end();

  force[0] = 0.;
  //We start the usual loop
  for( ; cell not_eq cend; ++cell ){
    // reinit all the needed stuff
    local_rhs = 0.;
    fe_values.reinit( cell );
    cell->get_dof_indices( local_dof_indices );

    // get the gradients of each function
    fe_values.get_function_gradients( u_n[0], grad_vel_1 );
    fe_values.get_function_gradients( u_n[1], grad_vel_2 );

    // usual loop over quad points and local dofs
    for( unsigned int q=0; q<n_q_points; ++q ){
      vorticity = grad_vel_2[q][1] - grad_vel_1[q][0];
      for( unsigned int i=0; i<dofs_per_cell; ++i )
        local_rhs(i) += fe_values.JxW(q)*vorticity*fe_values.shape_value( i, q );
    }

    // And copy it to the local one
    for( unsigned int i=0; i<dofs_per_cell; ++i )
      force[0]( local_dof_indices[i] ) += local_rhs(i);
  }

  // We solve a mass matrix problem
  SolverControl solver_control( 500, 1e-6*force[0].l2_norm() );
  SolverCG<> cg( solver_control );
  static bool is_prec_initted = false;
  static SparseILU<double> prec;
  if( not is_prec_initted ){
    prec.initialize( vel_Mass, SparseILU<double>::AdditionalData( 1e-5, 70 ) );
    is_prec_initted = true;
  }
  cg.solve( vel_Mass, force[1], force[0], prec );

  // Once we have the vorticity we can output it
  DataOut<dim> data_out;
  data_out.attach_dof_handler( dof_handler_velocity );
  data_out.add_data_vector( force[1], "vorticity" );
  data_out.build_patches();
  std::ostringstream filename;
  filename<<"vorticity"<<step<<".vtk";
  std::ofstream output( filename.str().c_str() );
  data_out.write_vtk( output );
}


// explicit template instantiation
template class Navier_Stokes_Projection<deal_II_dimension>;

#include "../include/FileReader.h"



// Constructor. Here we set up the
// Data Format we are going to use
Data_Storage::Data_Storage(){
  prm.declare_entry( "Method_Form", "rotational", Patterns::Selection( "rotational|standard" ),
                    " Used to select the type of method that we are going to use. " );

  // Physical data of the problem
  prm.enter_subsection( "Physical data" );
    prm.declare_entry( "initial_time", "0.", Patterns::Double( 0. ), " The initial time of the simulation. " );
    prm.declare_entry( "final_time", "1.", Patterns::Double( 0. ), " The final time of the simulation. " );
    prm.declare_entry( "Reynolds", "1.", Patterns::Double( 0. ), " The Reynolds number. " );
  prm.leave_subsection();

  // Time stepping data of the problem
  prm.enter_subsection( "Time step data" );
    prm.declare_entry( "initial_dt", "0.1", Patterns::Double( 0. ), " The initial time step size. " );
    prm.declare_entry( "final_dt", "5e-4",  Patterns::Double( 0. ), " The final time step size. " );
    prm.declare_entry( "dt_decrement", "2.", Patterns::Double( 1.5 ), " The factor by which the time step will be divided. " );
  prm.leave_subsection();

  // Space discretization data
  prm.enter_subsection( "Space discretization" );
    prm.declare_entry( "n_of_refines", "5", Patterns::Integer( 1, 15), " The number of global refines we do on the mesh. " );
    prm.declare_entry( "pressure_fe_degree", "1", Patterns::Integer( 1, 5 ), " The polynomial degree for the pressure space. " );
  prm.leave_subsection();

  // Velocity solution data
  prm.enter_subsection( "Data solve velocity" );
    prm.declare_entry( "max_iterations", "1000", Patterns::Integer( 1, 1000 ), " The maximal number of iterations GMRES must make. " );
    prm.declare_entry( "eps", "1e-12", Patterns::Double( 0. ), " The stopping criterion. " );
    prm.declare_entry( "Krylov_size", "30", Patterns::Integer(1), " The size of the Krylov subspace to be used. " );
    prm.declare_entry( "off_diagonals", "60", Patterns::Integer(0), " The number of off-diagonal elements ILU must compute. " );
    prm.declare_entry( "diag_strength", "0.01", Patterns::Double( 0. ), " Diagonal strengthening coefficient. " );
    prm.declare_entry( "update_prec", "15", Patterns::Integer(1), " This number indicates how often we need to update the preconditioner" );
  prm.leave_subsection();

  // Projection step data
  prm.enter_subsection( "Data solve projection" );
    prm.declare_entry( "max_iterations", "1000", Patterns::Integer( 1, 1000 ), " The maximal number of iterations CG must make. " );
    prm.declare_entry( "eps", "1e-12", Patterns::Double( 0. ), " The stopping criterion. " );
    prm.declare_entry( "off_diagonals", "100", Patterns::Integer(1), " The number of off-diagonal elements ILU must compute" );
    prm.declare_entry( "diag_strength", "0.1", Patterns::Double( 0. ), " Diagonal strengthening coefficient. " );
  prm.leave_subsection();

  // Pressure update data
  prm.enter_subsection( "Data solve pressure update" );
    prm.declare_entry( "max_iterations", "1000", Patterns::Integer( 1, 1000 ), " The maximal number of iterations CG must make. " );
    prm.declare_entry( "eps", "1e-12", Patterns::Double( 0. ), " The stopping criterion. " );
    prm.declare_entry( "off_diagonals", "10", Patterns::Integer(0), " The number of off-diagonal elements that ILU must compute" );
    prm.declare_entry( "diag_strength", "0.", Patterns::Double(0), " Diagonal strengthening coefficient" );
  prm.leave_subsection();

  // Verbosity of output
  prm.declare_entry( "verbose", "true", Patterns::Bool(), " This indicates whether the output of the solution process should be verbose. " );

  // How often we want the data to be outputted
  prm.declare_entry( "output", "10", Patterns::Integer(1), " This indicates between how many time steps we print the solution. " );
}



// Destructor. Does nothing
Data_Storage::~Data_Storage(){
}



// Here is where all happens. We read all the data from the indicated file
void Data_Storage::read_data( char *filename ){
  std::ifstream file;
  file.open( filename );
  if( not file )
    throw ExcFileNotOpen( filename );
  prm.read_input( file );

  std::string token = prm.get( "Method_Form" );
  if( token == std::string("rotational") )
    form = METHOD_ROTATIONAL;
  else
    form = METHOD_STANDARD;

  // Physical data of the problem
  prm.enter_subsection( "Physical data" );
    initial_time = prm.get_double( "initial_time" );
    final_time   = prm.get_double( "final_time" );
    Reynolds     = prm.get_double( "Reynolds" );
  prm.leave_subsection();

  // Time stepping data of the problem
  prm.enter_subsection( "Time step data" );
    initial_dt   = prm.get_double( "initial_dt" );
    final_dt     = prm.get_double( "final_dt" );
    dt_decrement = prm.get_double( "dt_decrement" );
  prm.leave_subsection();

  // Space discretization data
  prm.enter_subsection( "Space discretization" );
    n_of_global_refines = prm.get_integer( "n_of_refines" );
    pressure_degree     = prm.get_integer( "pressure_fe_degree" );
  prm.leave_subsection();

  // Velocity solution data
  prm.enter_subsection( "Data solve velocity" );
    vel_max_iterations = prm.get_double( "max_iterations" );
    vel_eps            = prm.get_double( "eps" );
    vel_Krylov_size    = prm.get_integer( "Krylov_size" );
    vel_off_diagonals  = prm.get_integer( "off_diagonals" );
    vel_diag_strength  = prm.get_double( "diag_strength" );
    vel_update_prec    = prm.get_integer( "update_prec" );
  prm.leave_subsection();

  // Projection step data
  prm.enter_subsection( "Data solve projection" );
    proj_max_iterations = prm.get_integer( "max_iterations" );
    proj_eps            = prm.get_double( "eps" );
    proj_off_diagonals  = prm.get_integer( "off_diagonals" );
    proj_diag_strength  = prm.get_double( "diag_strength" );
  prm.leave_subsection();

  // Pressure update data
  prm.enter_subsection( "Data solve pressure update" );
    pres_max_iterations = prm.get_integer( "max_iterations" );
    pres_eps            = prm.get_double( "eps" );
    pres_off_diagonals  = prm.get_integer( "off_diagonals" );
    pres_diag_strength  = prm.get_double( "diag_strength" );
  prm.leave_subsection();

  // Verbosity
  verbose = prm.get_bool( "verbose" );

  // Output frequency
  output = prm.get_integer( "output" );

  file.close();
}



// Prints the current values of the data
// Mostly (if not only) used for debugging purposes
void Data_Storage::print_status() const{
  std::cout<<"Method Form = "<<( (form == METHOD_ROTATIONAL)?"rotational":"standard" )<<std::endl
           <<"Physical Data:"<<std::endl
           <<"  initial_time = "<<initial_time<<std::endl
           <<"  final_time = "<<final_time<<std::endl
           <<"  Re = "<<Reynolds<<std::endl<<std::endl
           <<"Time step data:"<<std::endl
           <<"  initial_dt = "<<initial_dt<<std::endl
           <<"  final_dt = "<<final_dt<<std::endl
           <<"  dt_decrement = "<<dt_decrement<<std::endl<<std::endl
           <<"Space discretization data:"<<std::endl
           <<"  n_of_global_refines = "<<n_of_global_refines<<std::endl
           <<"  Pressure space degree = "<<pressure_degree<<std::endl<<std::endl
           <<"Data for the solution of the diffusion:"<<std::endl
           <<"  max iterations = "<<vel_max_iterations<<std::endl
           <<"  eps = "<<vel_eps<<std::endl
           <<"  Krylov subspace size = "<<vel_Krylov_size<<std::endl
           <<"  off diagonals = "<<vel_off_diagonals<<std::endl
           <<"  diagonal strengthening = "<<vel_diag_strength<<std::endl
           <<"  update prec every "<<vel_update_prec<<" steps"<<std::endl<<std::endl
           <<"Data for the projection step:"<<std::endl
           <<"  max iterations = "<<proj_max_iterations<<std::endl
           <<"  eps = "<<proj_eps<<std::endl
           <<"  off diagonals = "<<proj_off_diagonals<<std::endl
           <<"  diagonal strengthening = "<<proj_diag_strength<<std::endl<<std::endl
           <<"Data for the pressure update step"<<std::endl
           <<"  max iterations = "<<pres_max_iterations<<std::endl
           <<"  eps = "<<pres_eps<<std::endl
           <<"  off_diagonals = "<<pres_off_diagonals<<std::endl
           <<"  diagonal strengthening = "<<pres_diag_strength<<std::endl<<std::endl
           <<"Verbose = "<<(verbose?"true":"false")<<std::endl<<std::endl
           <<"Output ="<<output<<std::endl<<std::endl;
}



// This method just prints the format the
// input data file should have
void Data_Storage::print_usage(){
  static const char *message = "";
  std::cout<<message<<std::endl;
  prm.print_parameters( std::cout, ParameterHandler::Text );
}

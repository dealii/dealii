/* $Id$ */
/* Version: $Name:  $ */
/*      */
/*    Copyright (C) 2007, 2008, 2009 by the deal.II authors */
/*    Author: Abner Salgado, Texas A&M University 2009 */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */


				 // @sect3{Include Files}

				 // We start by including all the necessary
				 // dealii header files and some C++ related
				 // ones. Each one of them has been discussed
				 // in previous tutorial programs, so we will
				 // not get into details here.
#include <base/parameter_handler.h>
#include <base/point.h>
#include <base/function.h>
#include <base/quadrature_lib.h>
#include <base/multithread_info.h>
#include <base/thread_management.h>
#include <base/work_stream.h>
#include <base/parallel.h>

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_in.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_constraints.h>

#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <fe/fe_tools.h>
#include <fe/fe_system.h>

#include <lac/vector.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <lac/solver_gmres.h>
#include <lac/sparse_ilu.h>
#include <lac/sparse_direct.h>

#include <numerics/matrices.h>
#include <numerics/vectors.h>
#include <numerics/data_out.h>

#include <fstream>
#include <iostream>
#include <cmath>
#include <iostream>

// Finally we import all the dealii names to the global namespace
using namespace dealii;



// @sect3{Run Time parameters}
// Since our method has several options that can be fine-tuned we decided to group all these into
// an external file, so that these can be determined at run-time.<br>
// First, the formulation of the method, which we set as a <code>enum</code>.
// Next, we declare a class that is going to read and store all the parameters that our program
// needs to run.
namespace RunTimeParameters{
  enum Method_Formulation{
    METHOD_STANDARD,
    METHOD_ROTATIONAL
  };

  class Data_Storage{
    public:
      Data_Storage();
      ~Data_Storage();
      void read_data( const char *filename );
      Method_Formulation form;
      double initial_time,
             final_time,
             Reynolds;
      double dt;
      unsigned int n_of_global_refines,
                   pressure_degree;
      unsigned int vel_max_iterations,
                   vel_Krylov_size,
                   vel_off_diagonals,
                   vel_update_prec;
      double vel_eps,
             vel_diag_strength;
      bool verbose;
      unsigned int output;
    protected:
      ParameterHandler prm;
  };

// In the constructor of this class we declare all the parameters.
// The details of how this works have been discussed somewhere else ***
// so let's not elaborate on that
  Data_Storage::Data_Storage(){
    prm.declare_entry( "Method_Form", "rotational", Patterns::Selection( "rotational|standard" ),
                      " Used to select the type of method that we are going to use. " );
    prm.enter_subsection( "Physical data" );
      prm.declare_entry( "initial_time", "0.", Patterns::Double( 0. ), " The initial time of the simulation. " );
      prm.declare_entry( "final_time", "1.", Patterns::Double( 0. ), " The final time of the simulation. " );
      prm.declare_entry( "Reynolds", "1.", Patterns::Double( 0. ), " The Reynolds number. " );
    prm.leave_subsection();

    prm.enter_subsection( "Time step data" );
      prm.declare_entry( "dt", "5e-4",  Patterns::Double( 0. ), " The time step size. " );
    prm.leave_subsection();

    prm.enter_subsection( "Space discretization" );
      prm.declare_entry( "n_of_refines", "0", Patterns::Integer( 0, 15),
                         " The number of global refines we do on the mesh. " );
      prm.declare_entry( "pressure_fe_degree", "1", Patterns::Integer( 1, 5 ),
                         " The polynomial degree for the pressure space. " );
    prm.leave_subsection();

    prm.enter_subsection( "Data solve velocity" );
      prm.declare_entry( "max_iterations", "1000", Patterns::Integer( 1, 1000 ),
                         " The maximal number of iterations GMRES must make. " );
      prm.declare_entry( "eps", "1e-12", Patterns::Double( 0. ), " The stopping criterion. " );
      prm.declare_entry( "Krylov_size", "30", Patterns::Integer(1), " The size of the Krylov subspace to be used. " );
      prm.declare_entry( "off_diagonals", "60", Patterns::Integer(0),
                         " The number of off-diagonal elements ILU must compute. " );
      prm.declare_entry( "diag_strength", "0.01", Patterns::Double( 0. ),
                          " Diagonal strengthening coefficient. " );
      prm.declare_entry( "update_prec", "15", Patterns::Integer(1),
                         " This number indicates how often we need to update the preconditioner" );
    prm.leave_subsection();

    prm.declare_entry( "verbose", "true", Patterns::Bool(),
                       " This indicates whether the output of the solution process should be verbose. " );

    prm.declare_entry( "output", "1", Patterns::Integer(1),
                       " This indicates between how many time steps we print the solution. " );
  }

  Data_Storage::~Data_Storage(){}

  void Data_Storage::read_data( const char *filename ){
    std::ifstream file( filename );
    if( not file )
      throw ExcFileNotOpen( filename );
    prm.read_input( file );

    std::string token = prm.get( "Method_Form" );
    if( token == std::string( "rotational" ) )
      form = METHOD_ROTATIONAL;
    else
      form = METHOD_STANDARD;

    prm.enter_subsection( "Physical data" );
      initial_time = prm.get_double( "initial_time" );
      final_time   = prm.get_double( "final_time" );
      Reynolds     = prm.get_double( "Reynolds" );
      token        = prm.get( "filename" );
    prm.leave_subsection();

    prm.enter_subsection( "Time step data" );
      dt = prm.get_double( "dt" );
    prm.leave_subsection();

    prm.enter_subsection( "Space discretization" );
      n_of_global_refines = prm.get_integer( "n_of_refines" );
      pressure_degree     = prm.get_integer( "pressure_fe_degree" );
    prm.leave_subsection();

    prm.enter_subsection( "Data solve velocity" );
      vel_max_iterations = prm.get_double( "max_iterations" );
      vel_eps            = prm.get_double( "eps" );
      vel_Krylov_size    = prm.get_integer( "Krylov_size" );
      vel_off_diagonals  = prm.get_integer( "off_diagonals" );
      vel_diag_strength  = prm.get_double( "diag_strength" );
      vel_update_prec    = prm.get_integer( "update_prec" );
    prm.leave_subsection();

    verbose = prm.get_bool( "verbose" );

    output = prm.get_integer( "output" );

    file.close();
  }
}



// @sect3{The Equation Data}
// Here we declare the initial and boundary conditions.
namespace EquationData{
  // Because of our implementation, we do not take advantage of the capabilities of the library
  // to handle vector valued problems. Whether this is a good or bad idea is another issue. The
  // point here is that we want to be able to write an interface for the equation data that is
  // somehow dimension indenpendent. To be able to do that, our functions should be able to know
  // on which space component we are currently working, and we should be able to have a
  // common interface to do that. The following class is an attempt in that direction.
  template<int dim> class MultiComponentFunction: public Function<dim>{
    public:
      MultiComponentFunction( const double initial_time = 0. );
      void set_component( const unsigned int d );
    protected:
      unsigned int comp;
  };

  template<int dim> MultiComponentFunction<dim>::MultiComponentFunction( const double initial_time ):
                                                  Function<dim>( 1, initial_time ), comp(0){}

  template<int dim> void MultiComponentFunction<dim>::set_component(const unsigned int d ){
    Assert( d<dim, ExcIndexRange( d, 0, dim ) );
    comp = d;
  }

  template<int dim> class Velocity: public MultiComponentFunction<dim>{
    public:
      Velocity( const double initial_time = 0.0 );
      virtual double value( const Point<dim> &p, const unsigned int component = 0 ) const;
      virtual void value_list( const std::vector< Point<dim> > &points, std::vector<double> &values,
                                const unsigned int component = 0 ) const;
  };

  template<int dim> Velocity<dim>::Velocity( const double initial_time ):
                                    MultiComponentFunction<dim>( initial_time ){}

  template<int dim> void  Velocity<dim>::value_list( const std::vector<Point<dim> > &points,
                                                      std::vector<double> &values, const unsigned int ) const{
    const unsigned int n_points = points.size();
    Assert( values.size() == n_points, ExcDimensionMismatch( values.size(), n_points ) );
    for (unsigned int i=0; i<n_points; ++i)
      values[i] = Velocity<dim>::value( points[i] );
  }

  template<int dim> inline double Velocity<dim>::value( const Point<dim> &p, const unsigned int ) const{
    double return_value = 0.;
    static const double Um = 1.5, H = 4.1;
    if( MultiComponentFunction<dim>::comp == 0 )
      return_value = 4.*Um*p(1)*( H - p(1) )/(H*H);
    return return_value;
  }

  template<int dim> class Pressure: public Function<dim>{
    public:
      Pressure( const double initial_time = 0.0 );
      virtual double value( const Point<dim> &p, const unsigned int component = 0 ) const;
      virtual void value_list( const std::vector< Point<dim> > &points, std::vector<double> &values,
                                const unsigned int component = 0 ) const;
  };

  template<int dim> Pressure<dim>::Pressure( const double initial_time ): Function<dim>( 1, initial_time ){}

  template<int dim> inline double Pressure<dim>::value( const Point<dim> &p, const unsigned int ) const{
    return 0.;
  }

  template<int dim> void Pressure<dim>::value_list( const std::vector<Point<dim> > &points,
                                                      std::vector<double> &values,
                                                      const unsigned int ) const{
    const unsigned int n_points = points.size();
    Assert( values.size() == n_points, ExcDimensionMismatch( values.size(), n_points ) );
    for (unsigned int i=0; i<n_points; ++i)
      values[i] = Pressure<dim>::value( points[i] );
  }
}



// @sect3{The <code>Navier_Stokes_Projection</code> class}
// This is the main class of the program. It implements the various avatars of the projection
// methods for Navier-Stokes equations.
// The names for all the methods and attributes are self-explanatory.
template<int dim> class Navier_Stokes_Projection{
  public:
    Navier_Stokes_Projection( const RunTimeParameters::Data_Storage &data );
    ~Navier_Stokes_Projection();
    void run( const bool verbose = false, const unsigned int n_of_plots = 10 );
  protected:
    RunTimeParameters::Method_Formulation type;

    unsigned int deg;
    double dt;

    double t_0, T, Re;

    EquationData::Velocity<dim> vel_exact;
    std::map<unsigned int, double> boundary_values;
    std::vector<unsigned char> boundary_indicators;

    Triangulation<dim> triangulation;
    DoFHandler<dim> dof_handler_velocity, dof_handler_pressure;
    FE_Q<dim> fe_velocity, fe_pressure;
    QGauss<dim> quadrature_pressure, quadrature_velocity;

    SparsityPattern spar_pattern_velocity, spar_pattern_pressure, spar_pattern_pres_vel;
    SparseMatrix<double> vel_Laplace_plus_Mass, vel_it_matrix[dim], vel_Mass, vel_Laplace,
                         vel_Advection,
                         pres_Laplace, pres_Mass, pres_Diff[dim];
    ConstraintMatrix pres_regularization;
    Vector<double> pres_n, pres_n_minus_1, phi_n, phi_n_minus_1, u_n[dim], u_n_minus_1[dim],
                  u_star[dim],
                  force[dim],
                  v_tmp, pres_tmp,
                  rot_u;

    SparseILU<double> prec_velocity[dim];
    SparseDirectUMFPACK prec_mass, prec_pressure, prec_vel_mass;

    DeclException2( ExcInvalidTimeStep, double, double, <<" The time step "<<arg1<<" is out of range."<<std::endl
                                                                <<" The permitted range is (0,"<<arg2<<"]");

    void Create_Triangulation( const unsigned int n_of_refines );
    void Initialize();

    inline void interpolate_velocity();
    inline void diffusion_step( const bool reinit_prec );
    inline void projection_step( const bool reinit_prec );
    inline void update_pressure( const bool reinit_prec );
  private:
    unsigned int vel_max_its, vel_Krylov_size, vel_off_diagonals, vel_update_prec;
    double vel_eps, vel_diag_strength;

    inline void init_velocity_matrices();
    inline void init_pressure_matrices();
    inline void init_gradient_operator();

    typedef std_cxx1x::tuple< typename DoFHandler<dim>::active_cell_iterator,
                                typename DoFHandler<dim>::active_cell_iterator
                             > IteratorTuple;
    typedef parallel::internal::SynchronousIterators<IteratorTuple> SIterators;
    struct InitGradPerTaskData{
      unsigned int d, vel_dpc, pres_dpc;
      FullMatrix<double> local_grad;
      std::vector<unsigned int> vel_local_dof_indices, pres_local_dof_indices;
      InitGradPerTaskData( const unsigned int dd, const unsigned int vdpc, const unsigned int pdpc ):
                                             d(dd), vel_dpc( vdpc ), pres_dpc( pdpc ),
                                             local_grad( vdpc, pdpc ), vel_local_dof_indices( vdpc ),
                                             pres_local_dof_indices( pdpc ){}
    };
    struct InitGradScratchData{
      unsigned int nqp;
      FEValues<dim> fe_val_vel, fe_val_pres;
      InitGradScratchData( const FE_Q<dim> &fe_v, const FE_Q<dim> &fe_p, const QGauss<dim> &quad,
                           const UpdateFlags flags_v, const UpdateFlags flags_p ) :
                              nqp( quad.size() ), fe_val_vel( fe_v, quad, flags_v ),
                              fe_val_pres( fe_p, quad, flags_p ){}
      InitGradScratchData( const InitGradScratchData &data ): nqp( data.nqp ),
                                        fe_val_vel( data.fe_val_vel.get_fe(), data.fe_val_vel.get_quadrature(),
                                                    data.fe_val_vel.get_update_flags() ),
                                        fe_val_pres( data.fe_val_pres.get_fe(), data.fe_val_pres.get_quadrature(),
                                                     data.fe_val_pres.get_update_flags() ) {}
    };
    void assemble_one_cell_of_gradient( const SIterators &SI, InitGradScratchData &scratch,
                                        InitGradPerTaskData &data );
    void copy_gradient_local_to_global( const InitGradPerTaskData &data );

    inline void assemble_advection_term();
    struct AdvectionPerTaskData{
      FullMatrix<double> local_advection;
      std::vector<unsigned int> local_dof_indices;
      AdvectionPerTaskData( const unsigned int dpc ): local_advection( dpc, dpc ), local_dof_indices( dpc ) {}
    };
    struct AdvectionScratchData{
      unsigned int nqp, dpc;
      std::vector< Point<dim> > u_star_local;
      std::vector< Tensor<1,dim> > grad_u_star;
      std::vector<double> u_star_tmp;
      FEValues<dim> fe_val;
      AdvectionScratchData( const FE_Q<dim> &fe, const QGauss<dim> &quad, const UpdateFlags flags ):
                            nqp( quad.size() ), dpc( fe.dofs_per_cell ),
                            u_star_local( nqp ), grad_u_star( nqp ), u_star_tmp( nqp ),
                            fe_val( fe, quad, flags ){}
      AdvectionScratchData( const AdvectionScratchData &data ): nqp( data.nqp ), dpc( data.dpc ),
                            u_star_local( nqp ), grad_u_star( nqp ), u_star_tmp( nqp ),
                            fe_val( data.fe_val.get_fe(), data.fe_val.get_quadrature(),
                                    data.fe_val.get_update_flags() ) {}
    };

    void assemble_one_cell_of_advection( const typename DoFHandler<dim>::active_cell_iterator &cell,
                                         AdvectionScratchData &scratch, AdvectionPerTaskData &data );
    void copy_advection_local_to_global( const AdvectionPerTaskData &data );
    inline void diffusion_component_solve( const unsigned int d );

    inline void plot_solution( const unsigned int step );
    inline void assemble_vorticity( const bool reinit_prec );
};


template<int dim> Navier_Stokes_Projection<dim>::~Navier_Stokes_Projection(){
  dof_handler_velocity.clear();
  dof_handler_pressure.clear();
}



// @sect4{ <code>Navier_Stokes_Projection::Navier_Stokes_Projection</code> }
// In the constructor, we just read all the data from the <code>Data_Storage</code>
// object that is passed as an argument, verify that the data we read is reasonable
// and, finally, create the triangulation and load the initial data.
template<int dim> Navier_Stokes_Projection<dim>::Navier_Stokes_Projection(
                                        const RunTimeParameters::Data_Storage &data ):
                    type( data.form ), deg( data.pressure_degree ), dt( data.dt ), t_0( data.initial_time ),
                    T( data.final_time ), Re( data.Reynolds ), vel_exact( data.initial_time ),
                    dof_handler_velocity( triangulation ), dof_handler_pressure( triangulation ),
                    fe_velocity( deg+1 ), fe_pressure( deg ),
                    quadrature_pressure( deg+1 ), quadrature_velocity( deg+2 ),
                    vel_max_its( data.vel_max_iterations ), vel_Krylov_size( data.vel_Krylov_size ),
                    vel_off_diagonals( data.vel_off_diagonals ),
                    vel_update_prec( data.vel_update_prec ), vel_eps( data.vel_eps ),
                    vel_diag_strength( data.vel_diag_strength)
{
  if(deg < 1)
  std::cout<<" WARNING: The chosen pair of finite element spaces is not stable."<<std::endl
           <<" The obtained results will be nonsense"<<std::endl;

  AssertThrow( not ( ( dt <= 0. ) or ( dt > .5*T ) ), ExcInvalidTimeStep( dt, .5*T ) );

  Create_Triangulation( data.n_of_global_refines );
  Initialize();
}


// @sect4{ <code>Navier_Stokes_Projection::Create_Triangulation</code> }
// The method that creates the triangulation and refines it the needed number of times.
// After creating the triangulation, it creates the mesh dependent data, i.e. it distributes
// degrees of freedom and renumbers them, and initializes the matrices and vectors
// that we will use.
template<int dim> void Navier_Stokes_Projection<dim>::Create_Triangulation( const unsigned int n_of_refines ){
  GridIn<dim> grid_in;
  grid_in.attach_triangulation( triangulation );

  std::string filename = "nsbench2.inp";
  std::ifstream file( filename.c_str() );
  Assert( file, ExcFileNotOpen( filename.c_str() ) );
  grid_in.read_ucd( file );
  file.close();

  std::cout<<" Number of refines = "<<n_of_refines<<std::endl;
  triangulation.refine_global( n_of_refines );
  std::cout<<" Number of active cells: "<<triangulation.n_active_cells()<<std::endl;

  boundary_indicators = triangulation.get_boundary_indicators();

  dof_handler_velocity.distribute_dofs( fe_velocity );
  DoFRenumbering::boost::Cuthill_McKee( dof_handler_velocity );
  dof_handler_pressure.distribute_dofs( fe_pressure );
  DoFRenumbering::boost::Cuthill_McKee( dof_handler_pressure );

  init_velocity_matrices();
  init_pressure_matrices();
  init_gradient_operator();

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
  rot_u.reinit( dof_handler_velocity.n_dofs() );

  std::cout<<" dim( X_h ) = "<<( dof_handler_velocity.n_dofs()*dim )<<std::endl
           <<" dim( M_h ) = "<<dof_handler_pressure.n_dofs()<<std::endl
           <<" Re         = "<<Re<<std::endl;
}


// @sect4{ <code>Navier_Stokes_Projection::Initialize</code> }
// This method creates the constant matrices and loads the initial data
template<int dim> void Navier_Stokes_Projection<dim>::Initialize(){
  vel_Laplace_plus_Mass = 0.;
  vel_Laplace_plus_Mass.add( 1./Re, vel_Laplace );
  vel_Laplace_plus_Mass.add( 1.5/dt, vel_Mass );

  EquationData::Pressure<dim> pres( t_0 );
  VectorTools::interpolate( dof_handler_pressure, pres, pres_n_minus_1 );
  pres.advance_time( dt );
  VectorTools::interpolate( dof_handler_pressure, pres, pres_n );
  phi_n = 0.;
  phi_n_minus_1 = 0.;
  for(unsigned int d=0; d<dim; ++d){
    vel_exact.set_time( t_0 );
    vel_exact.set_component(d);
    VectorTools::interpolate( dof_handler_velocity, vel_exact, u_n_minus_1[d] );
    vel_exact.advance_time( dt );
    VectorTools::interpolate( dof_handler_velocity, vel_exact, u_n[d] );
  }
}


// @sect4{ The <code>Navier_Stokes_Projection::init_*_matrices</code> methods }
// In this set of methods we initialize the sparsity patterns, the constraints (if any)
// and assemble the matrices that do not depend on the timestep $\Delta t$.
template<int dim> void Navier_Stokes_Projection<dim>::init_velocity_matrices(){
  spar_pattern_velocity.reinit( dof_handler_velocity.n_dofs(), dof_handler_velocity.n_dofs(),
                                dof_handler_velocity.max_couplings_between_dofs() );
  DoFTools::make_sparsity_pattern( dof_handler_velocity, spar_pattern_velocity );
  spar_pattern_velocity.compress();

  vel_Laplace_plus_Mass.reinit( spar_pattern_velocity );
  for( unsigned int d=0; d<dim; ++d )
    vel_it_matrix[d].reinit( spar_pattern_velocity );
  vel_Mass.reinit( spar_pattern_velocity );
  vel_Laplace.reinit( spar_pattern_velocity );
  vel_Advection.reinit( spar_pattern_velocity );

  MatrixCreator::create_mass_matrix( dof_handler_velocity, quadrature_velocity, vel_Mass );
  MatrixCreator::create_laplace_matrix( dof_handler_velocity, quadrature_velocity, vel_Laplace );
}

// For the initialization of the matrices that act on the pressure space it is worth noticing
// one small detail. Since the projection step involves the solution
// of a Poisson equation with homogeneous Neumann boundary conditions, we need somehow to
// regularize this problem, that is to pick a solution.
// The way we do it is by setting the value of the solution at the first node (wherever it
// is) to zero. This regularizes the problem and does not increase the sparsity pattern we
// use.
template<int dim> void Navier_Stokes_Projection<dim>::init_pressure_matrices(){
  spar_pattern_pressure.reinit( dof_handler_pressure.n_dofs(), dof_handler_pressure.n_dofs(),
                                dof_handler_pressure.max_couplings_between_dofs() );
  DoFTools::make_sparsity_pattern( dof_handler_pressure, spar_pattern_pressure );

  pres_regularization.clear();
  pres_regularization.add_line(0);
  pres_regularization.close();
  pres_regularization.condense( spar_pattern_pressure );

  spar_pattern_pressure.compress();

  pres_Laplace.reinit( spar_pattern_pressure );
  pres_Mass.reinit( spar_pattern_pressure );

  MatrixCreator::create_laplace_matrix( dof_handler_pressure, quadrature_pressure, pres_Laplace );
  MatrixCreator::create_mass_matrix(  dof_handler_pressure, quadrature_pressure, pres_Mass );

  pres_regularization.condense( pres_Laplace );
}


// For the gradient operator, we start by initializing the sparsity pattern and compressing it.
// It is important to notice here that the gradient operator acts from the pressure space
// into the velocity space, so we have to deal with two different finite element spaces. To keep
// the loops synchronized, we use the <code>typedef</code>'s that we have defined before, namely
// <code>PairedIterators</code> and <code>SIterators</code>.
template<int dim> void Navier_Stokes_Projection<dim>::init_gradient_operator(){
  spar_pattern_pres_vel.reinit( dof_handler_velocity.n_dofs(), dof_handler_pressure.n_dofs(),
                                dof_handler_velocity.max_couplings_between_dofs() );
  DoFTools::make_sparsity_pattern( dof_handler_velocity, dof_handler_pressure, spar_pattern_pres_vel );
  spar_pattern_pres_vel.compress();

  InitGradPerTaskData per_task_data( 0, fe_velocity.dofs_per_cell, fe_pressure.dofs_per_cell );
  InitGradScratchData scratch_data( fe_velocity, fe_pressure, quadrature_velocity,
                                    update_values | update_JxW_values, update_gradients );

  for( unsigned int d=0; d<dim; ++d ){
    pres_Diff[d].reinit( spar_pattern_pres_vel );
    per_task_data.d = d;
    WorkStream::run( SIterators( IteratorTuple( dof_handler_velocity.begin_active(),
                                                dof_handler_pressure.begin_active()
                                              )
                               ),
                     SIterators( IteratorTuple( dof_handler_velocity.end(),
                                                dof_handler_pressure.end()
                                              )
                               ),
                     *this,
                     &Navier_Stokes_Projection<dim>::assemble_one_cell_of_gradient,
                     &Navier_Stokes_Projection<dim>::copy_gradient_local_to_global,
                     scratch_data,
                     per_task_data
                   );
  }
}

template<int dim> void Navier_Stokes_Projection<dim>::assemble_one_cell_of_gradient( const SIterators &SI,
                                                              InitGradScratchData &scratch,
                                                              InitGradPerTaskData &data ){
  scratch.fe_val_vel.reinit( std_cxx1x::get<0>( SI.iterators ) );
  scratch.fe_val_pres.reinit( std_cxx1x::get<1>( SI.iterators ) );

  std_cxx1x::get<0>( SI.iterators )->get_dof_indices( data.vel_local_dof_indices );
  std_cxx1x::get<1>( SI.iterators )->get_dof_indices( data.pres_local_dof_indices );

  data.local_grad = 0.;
  for( unsigned int q=0; q<scratch.nqp; ++q ){
    for( unsigned int i=0; i<data.vel_dpc; ++i )
      for( unsigned int j=0; j<data.pres_dpc; ++j )
        data.local_grad( i, j ) += scratch.fe_val_vel.JxW(q)*scratch.fe_val_vel.shape_value( i, q )
                                   *scratch.fe_val_pres.shape_grad( j, q )[data.d];
  }
}

template<int dim> void Navier_Stokes_Projection<dim>::copy_gradient_local_to_global(
                                                                  const InitGradPerTaskData &data ){
  for( unsigned int i=0; i<data.vel_dpc; ++i )
    for( unsigned int j=0; j<data.pres_dpc; ++j)
      pres_Diff[data.d].add( data.vel_local_dof_indices[i], data.pres_local_dof_indices[j],
                             data.local_grad( i, j)  );
}


// @sect4{ <code>Navier_Stokes_Projection::run</code> }
// This is the time marching function, which starting at <code>t_0</code> advances in time
// using the projection method with time step <code>dt</code> until <code>T</code>. <br>
// The boolean parameter, <code>verbose</code>, that it takes is to enable
// information about what the method is doing at the given moment, i.e. diffusion, projection
// substep; updating preconditioners etc. This is useful mostly for debugging purposes
// and so it is by default set to false
template<int dim> void Navier_Stokes_Projection<dim>::run( const bool verbose, const unsigned int n_of_plots ){
  unsigned int n_steps = ( T - t_0 )/dt;
  vel_exact.set_time( 2.*dt );
  plot_solution(1);
  for( unsigned int n = 2; n<=n_steps; ++n ){
    if( n%n_of_plots == 0 ){
      if( verbose )
        std::cout<<" Plotting Solution"<<std::endl;
      plot_solution(n);
    }
    if( verbose )
      std::cout<<" Step = "<<n<<" Time = "<<(n*dt)<<std::endl;
    if( verbose )
      std::cout<<"  Interpolating the velocity "<<std::endl;
    interpolate_velocity();
    if( verbose )
      std::cout<<"  Diffusion Step"<<std::endl;
    if( ( n%vel_update_prec == 0 ) and ( verbose ) )
      std::cout<<"   With reinitialization of the preconditioner"<<std::endl;
    diffusion_step( (n%vel_update_prec == 0 ) or ( n == 2) );
    if( verbose )
      std::cout<<"  Projection Step"<<std::endl;
    projection_step( ( n == 2 ) );
    if( verbose )
      std::cout<<"  Updating the Pressure"<<std::endl;
    update_pressure( ( n == 2 ) );
    vel_exact.advance_time(dt);
  }
  plot_solution( n_steps );
}

template<int dim> void Navier_Stokes_Projection<dim>::interpolate_velocity(){
  for( unsigned int d=0; d<dim; ++d )
    u_star[d].equ( 2., u_n[d], -1, u_n_minus_1[d] );
}


// @sect4{<code>Navier_Stokes_Projection::diffusion_step</code>}
// The implementation of a diffusion step.
template<int dim> void Navier_Stokes_Projection<dim>::diffusion_step( const bool reinit_prec ){
  pres_tmp.equ( -1., pres_n, -4./3., phi_n, 1./3., phi_n_minus_1 );

  assemble_advection_term();

  for( unsigned int d=0; d<dim; ++d ){
    force[d] = 0.;
    v_tmp.equ( 2./dt,u_n[d],-.5/dt,u_n_minus_1[d] );
    vel_Mass.vmult_add( force[d], v_tmp );

    pres_Diff[d].vmult_add( force[d], pres_tmp );
    u_n_minus_1[d] = u_n[d];

    vel_it_matrix[d].copy_from( vel_Laplace_plus_Mass );
    vel_it_matrix[d].add( 1., vel_Advection );

    vel_exact.set_component(d);
    std::vector<unsigned char>::const_iterator boundaries = boundary_indicators.begin(),
                                               b_end = boundary_indicators.end();
    boundary_values.clear();
    for( ; boundaries not_eq b_end; ++boundaries ){
      switch( *boundaries ){
        case 1:
          VectorTools::interpolate_boundary_values( dof_handler_velocity, *boundaries,
                                                    ZeroFunction<dim>(), boundary_values );
          break;
        case 2:
          VectorTools::interpolate_boundary_values( dof_handler_velocity, *boundaries,
                                                    vel_exact, boundary_values );
          break;
        case 3:
          if( d not_eq 0 )
            VectorTools::interpolate_boundary_values( dof_handler_velocity, *boundaries,
                                                      ZeroFunction<dim>(), boundary_values );
          break;
        case 4:
          VectorTools::interpolate_boundary_values( dof_handler_velocity, *boundaries,
                                                    ZeroFunction<dim>(), boundary_values );
          break;
        default:
          Assert( false, ExcNotImplemented() );
      }
    }
    MatrixTools::apply_boundary_values( boundary_values, vel_it_matrix[d], u_n[d], force[d] );
  }


  Threads::TaskGroup<void> tasks;
  for(unsigned int d=0; d<dim; ++d ){
    if( reinit_prec )
      prec_velocity[d].initialize( vel_it_matrix[d],
                                  SparseILU<double>::AdditionalData( vel_diag_strength, vel_off_diagonals ) );
    tasks += Threads::new_task( &Navier_Stokes_Projection<dim>::diffusion_component_solve, *this, d );
  }
  tasks.join_all();
}

template<int dim> void Navier_Stokes_Projection<dim>::diffusion_component_solve( const unsigned int d ){
  SolverControl solver_control( vel_max_its, vel_eps*force[d].l2_norm() );
  SolverGMRES<> gmres( solver_control, SolverGMRES<>::AdditionalData( vel_Krylov_size ) );
  gmres.solve( vel_it_matrix[d], u_n[d], force[d], prec_velocity[d] );
}


// @sect4{ The <code>Navier_Stokes_Projection::assemble_advection_term</code>  method and related}
template<int dim> void Navier_Stokes_Projection<dim>::assemble_advection_term(){
  vel_Advection = 0.;
  AdvectionPerTaskData data( fe_velocity.dofs_per_cell );
  AdvectionScratchData scratch( fe_velocity, quadrature_velocity,
                                update_values | update_JxW_values | update_gradients );
  WorkStream::run( dof_handler_velocity.begin_active(), dof_handler_velocity.end(), *this,
                   &Navier_Stokes_Projection<dim>::assemble_one_cell_of_advection,
                   &Navier_Stokes_Projection<dim>::copy_advection_local_to_global, scratch, data);
}

template<int dim> void Navier_Stokes_Projection<dim>::assemble_one_cell_of_advection(
                                                  const typename DoFHandler<dim>::active_cell_iterator &cell,
                                                  AdvectionScratchData &scratch, AdvectionPerTaskData &data ){
  scratch.fe_val.reinit(cell);
  cell->get_dof_indices( data.local_dof_indices );
  for( unsigned int d=0; d<dim; ++d ){
    scratch.fe_val.get_function_values( u_star[d], scratch.u_star_tmp );
    for( unsigned int q=0; q<scratch.nqp; ++q )
      scratch.u_star_local[q](d) = scratch.u_star_tmp[q];
  }

  for( unsigned int d=0; d<dim; ++d ){
    scratch.fe_val.get_function_gradients( u_star[d], scratch.grad_u_star );
    for( unsigned int q=0; q<scratch.nqp; ++q ){
      if( d==0 )
        scratch.u_star_tmp[q] = 0.;
      scratch.u_star_tmp[q] += scratch.grad_u_star[q][d];
    }
  }

  data.local_advection = 0.;
  for( unsigned int q=0; q<scratch.nqp; ++q )
    for( unsigned int i=0; i<scratch.dpc; ++i )
      for( unsigned int j=0; j<scratch.dpc; ++j )
        data.local_advection(i,j) += (
                scratch.u_star_local[q]*scratch.fe_val.shape_grad( j, q )*scratch.fe_val.shape_value( i, q )
                + 0.5*scratch.u_star_tmp[q]*scratch.fe_val.shape_value( i, q )*scratch.fe_val.shape_value( j, q )
                                    )*scratch.fe_val.JxW(q) ;
}

template<int dim> void Navier_Stokes_Projection<dim>::copy_advection_local_to_global(
                                                        const AdvectionPerTaskData &data ){
  for( unsigned int i=0; i<fe_velocity.dofs_per_cell; ++i )
    for( unsigned int j=0; j<fe_velocity.dofs_per_cell; ++j )
      vel_Advection.add( data.local_dof_indices[i] , data.local_dof_indices[j], data.local_advection(i,j) );
}



// @sect4{<code>Navier_Stokes_Projection::projection_step</code>}
// This implements the projection step.
template<int dim> void Navier_Stokes_Projection<dim>::projection_step( const bool reinit_prec ){
  if( reinit_prec )
    prec_pressure.initialize( pres_Laplace );

  pres_tmp = 0.;
  for( unsigned d=0; d<dim; ++d )
    pres_Diff[d].Tvmult_add( pres_tmp, u_n[d] );

  phi_n_minus_1 = phi_n;
  phi_n = pres_tmp;
  phi_n *= 1.5/dt;

  pres_regularization.condense( phi_n );
  prec_pressure.solve( phi_n );
  pres_regularization.distribute( phi_n );
}


// @sect4{ <code>Navier_Stokes_Projection::update_pressure</code> }
// This is the pressure update step of the projection method. It implements the
// standard formulation of the method, that is
// @f[
//      p^{n+1} = p^n + \phi^{n+1},
// @f]
// or the rotational form, which is
// @f[
//      p^{n+1} = p^n + \phi^{n+1} - \frac{1}{Re} \nabla\cdot u^{n+1}.
// @f]
template<int dim> void Navier_Stokes_Projection<dim>::update_pressure( const bool reinit_prec ){
  pres_n_minus_1 = pres_n;
  switch( type ){
    case RunTimeParameters::METHOD_STANDARD:
      pres_n += phi_n;
      break;
    case RunTimeParameters::METHOD_ROTATIONAL:
      if( reinit_prec )
        prec_mass.initialize( pres_Mass );
      pres_n = pres_tmp;
      prec_mass.solve( pres_n );
      pres_n.sadd(1./Re, 1., pres_n_minus_1, 1., phi_n );
      break;
    default:
      Assert( false, ExcNotImplemented() );
  };
}


// @sect4{ <code>Navier_Stokes_Projection::plot_solution</code> }
// This method plots the current solution. It is an adaptation of
// step-31 **** and so I will not elaborate on it.
// There is one small detail here. It is often interested to see
// the vorticity of the flow. But, since we are using it here only for plotting
// purposes, we are not going to compute it at every time step, but only when
// we are going to plot it.
template<int dim> void Navier_Stokes_Projection<dim>::plot_solution( const unsigned int step ){
  assemble_vorticity( ( step == 1 ) );
  const FESystem<dim> joint_fe( fe_velocity, dim, fe_pressure, 1, fe_velocity, 1 );
  DoFHandler<dim> joint_dof_handler( triangulation );
  joint_dof_handler.distribute_dofs( joint_fe );
  Assert( joint_dof_handler.n_dofs() == (dim + 1)*dof_handler_velocity.n_dofs() + dof_handler_pressure.n_dofs(),
          ExcInternalError() );
  static Vector<double> joint_solution( joint_dof_handler.n_dofs() );
  std::vector<unsigned int> loc_joint_dof_indices( joint_fe.dofs_per_cell ),
                            loc_vel_dof_indices( fe_velocity.dofs_per_cell ),
                            loc_pres_dof_indices( fe_pressure.dofs_per_cell );
  typename DoFHandler<dim>::active_cell_iterator
                              joint_cell = joint_dof_handler.begin_active(),
                              joint_endc = joint_dof_handler.end(),
                              vel_cell   = dof_handler_velocity.begin_active(),
                              pres_cell  = dof_handler_pressure.begin_active();
  for( ; joint_cell not_eq joint_endc; ++joint_cell, ++vel_cell, ++pres_cell ){
    joint_cell->get_dof_indices( loc_joint_dof_indices );
    vel_cell->get_dof_indices( loc_vel_dof_indices ),
    pres_cell->get_dof_indices( loc_pres_dof_indices );
    for( unsigned int i=0; i<joint_fe.dofs_per_cell; ++i )
      switch( joint_fe.system_to_base_index(i).first.first ){
        case 0:
          Assert( joint_fe.system_to_base_index(i).first.second < dim, ExcInternalError() );
          joint_solution( loc_joint_dof_indices[i] ) =
                              u_n[ joint_fe.system_to_base_index(i).first.second ]
                                      ( loc_vel_dof_indices[ joint_fe.system_to_base_index(i).second ] );
          break;
        case 1:
          Assert( joint_fe.system_to_base_index(i).first.second == 0, ExcInternalError() );
          joint_solution( loc_joint_dof_indices[i] ) =
                                          pres_n( loc_pres_dof_indices[ joint_fe.system_to_base_index(i).second ] );
          break;
        case 2:
          Assert( joint_fe.system_to_base_index(i).first.second == 0, ExcInternalError() );
          joint_solution( loc_joint_dof_indices[i] ) =
                          rot_u( loc_vel_dof_indices[ joint_fe.system_to_base_index(i).second ] );
          break;
        default:
          Assert( false, ExcInternalError() );
      }
  }
  std::vector<std::string> joint_solution_names( dim, "v" );
  joint_solution_names.push_back( "p" );
  joint_solution_names.push_back( "rot_u" );
  DataOut<dim> data_out;
  data_out.attach_dof_handler (joint_dof_handler);
  std::vector< DataComponentInterpretation::DataComponentInterpretation >
            component_interpretation( dim+2, DataComponentInterpretation::component_is_part_of_vector );
  component_interpretation[dim]   = DataComponentInterpretation::component_is_scalar;
  component_interpretation[dim+1] = DataComponentInterpretation::component_is_scalar;
  data_out.add_data_vector( joint_solution, joint_solution_names, DataOut<dim>::type_dof_data,
                            component_interpretation );
  data_out.build_patches( deg + 1 );
  std::ostringstream filename;
  filename<<"solution-"<<step<<".vtk";
  std::ofstream output( filename.str().c_str() );
  data_out.write_vtk( output );
}



// Since this function is supposed to be called only when the plot is going to be made,
// which should not be every time step, we do not parallelize it. Of course, if needed,
// this can be done as in the other cases. Moreover, the implementation that we have
// here only works for 2d, so we bail if that is not the case.
template<int dim> void Navier_Stokes_Projection<dim>::assemble_vorticity( const bool reinit_prec ){
  Assert( dim == 2, ExcNotImplemented() );
  if( reinit_prec )
    prec_vel_mass.initialize( vel_Mass );

  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_velocity.begin_active(),
                                                   end = dof_handler_velocity.end();
  FEValues<dim> fe_val_vel( fe_velocity, quadrature_velocity,
                            update_gradients | update_JxW_values | update_values );
  const unsigned int dpc = fe_velocity.dofs_per_cell,
                     nqp = quadrature_velocity.size();
  std::vector<unsigned int> ldi( dpc );
  Vector<double> loc_rot( dpc );

  std::vector< Tensor<1,dim> > grad_u1( nqp ), grad_u2( nqp );
  rot_u = 0.;
  for( ; cell not_eq end; ++cell ){
    fe_val_vel.reinit( cell );
    cell->get_dof_indices( ldi );
    fe_val_vel.get_function_gradients( u_n[0], grad_u1 );
    fe_val_vel.get_function_gradients( u_n[1], grad_u2 );
    loc_rot = 0.;
    for( unsigned int q=0; q<nqp; ++q )
      for( unsigned int i=0; i<dpc; ++i )
        loc_rot(i) += fe_val_vel.JxW(q)*( grad_u2[q][0] - grad_u1[q][1] )*fe_val_vel.shape_value( i, q );

    for( unsigned int i=0; i<dpc; ++i )
      rot_u( ldi[i] ) += loc_rot(i);
  }

  prec_vel_mass.solve( rot_u );
}


// @sect3{ The main function }
// The main function looks very much like in all the other tutorial programs.
int main(){
  try{
    RunTimeParameters::Data_Storage data;
    data.read_data( "parameter-file.prm" );
    deallog.depth_console( data.verbose?2:0 );
    Navier_Stokes_Projection<2> test( data );
    test.run( data.verbose, data.output );
  }
  catch (std::exception &exc){
    std::cerr << std::endl << std::endl
            << "----------------------------------------------------"
            << std::endl;
    std::cerr << "Exception on processing: " << std::endl
            << exc.what() << std::endl
            << "Aborting!" << std::endl
            << "----------------------------------------------------"
            << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  std::cout<<"----------------------------------------------------"
           <<std::endl
           <<"Apparently everything went fine!"
           <<std::endl
           <<"Don't forget to brush your teeth :-)"
           <<std::endl<<std::endl;
  return 0;
}

/*
    Definition of the Projection Method for constant density
    Navier Stokes
    by Abner Salgado.
*/
#ifndef _NAVIER_STOKES_H_
#define _NAVIER_STOKES_H_




/*
    First we include the local files, where we have the definitions of
    1.- The exact solution
    2.- The class that reads runtime parameters from a file
*/
#include "EqData.h"
#include "FileReader.h"



/*
    These includes are to get the quadrature
    and to handle the convergence table
*/
#include <base/quadrature_lib.h>
#include <base/convergence_table.h>



/*
    The includes needed to manage multithreading
*/
#include <base/multithread_info.h>
#include <base/thread_management.h>



/*
    Grid management includes.
    Triangulation handler
    Grid generator & refinement, accessor and iterators
*/
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_in.h>



/*
    Degrees of freedom management includes
    Handling, accessing, computing useful stuff from and for them and renumbering.
*/
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_constraints.h>



/*
    Finite element management includes
    Definition of finite elements and extraction of info from them
    and from finite element functions.
*/
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <fe/fe_tools.h>



/*
    Linear algebra includes
    Sparse matrices, vectors, preconditioners and solvers.
*/
#include <lac/vector.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <lac/solver_gmres.h>
#include <lac/sparse_ilu.h>
#include <lac/sparse_direct.h>



/*
    Numerical algorithms includes
    Assembly of right hand sides, computation of errors, output of the solution, etc
*/
#include <numerics/matrices.h>
#include <numerics/vectors.h>
#include <numerics/data_out.h>


/*
    This class is the core of all the program
    It is the implementation of the projection method for the Navier-Stokes equations.
*/
template<int dim> class Navier_Stokes_Projection{
  public:
    /*
      Constructor. It takes as a parameter a reference to a Data_Storage class, which is
      basically a container of all the data this class needs to be defined (and some other stuff)
      Here we basically just copy the needed data and/or assign the values that are needed.
    */
    Navier_Stokes_Projection( const Data_Storage &data );
    /*
      Destructor. Completely trivial.
    */
    ~Navier_Stokes_Projection();
    /*
      This function creates the triangulation and refines it the required number of times
      It also initializes all the quantities that are mesh dependent but not time stepping dependent
      i.e. the size of matrices and vectors.
    */
    void Create_Triangulation( const unsigned int n_of_refines );
    /*
      Having created a mesh and once the time step is set, the next step is to initialize the method
      i.e. compute the matrices that are never going to change, load the initial data, etc.
      This is what this function does.
    */
    void Initialize();
    /*
      This is the time marching function, which starting at t_0 advances in time using the projection method
      with time step dt until T.
      The boolean parameter that it takes is to enable information about what the function is doing at the present
      moment, i.e. diffusion, projection substep; updating preconditioners etc.
      This is useful mostly for debuggin purposes and so it is by default set to false
    */
    void run( const bool verbose = false, const unsigned int n_of_plots = 10 );
    /*
      Having reached the final time T, we want to measure the error that we have made.
      This method is responsible for that. Saves the results in a ConvergenceTable object
      which later we can print or compute things with it
    */
    void Post_Process();
    /*
      The whole reason for this class, and all this code after all, is to make convergence tests
      with respect to the time discretization. For this reason we need to be able to vary the time step
      we are going to work with. This method sets up the time step that is to be used in the given run
    */
    void set_dt( const double ddt );

  protected:
    // The type of method
    Method_Formulation type;

    // Discretization data
    //// The polynomial degree
    unsigned int deg;
    //// The time step
    double dt;

    // Physical data:
    //// initial time, final time, Reynolds number
    double t_0, T, Re;
    //// The external driving force
    Force<dim> rhs;
    //// The exact velocity (to apply boundary values)
    Velocity<dim> vel_exact;
    //// The boundary conditions
    std::map<unsigned int, double> boundary_values;

    // Finite Element Spaces data
    //// The mesh
    Triangulation<dim> triangulation;
    //// The DoF handlers
    DoFHandler<dim> dof_handler_velocity, dof_handler_pressure;
    //// The polynomial spaces
    FE_Q<dim> fe_velocity, fe_pressure;
    //// The quadrature formulae
    QGauss<dim> quadrature_pressure, quadrature_velocity;

    /*
      Linear Algebra Data
        The sparsity patterns where the matrices will live
    */
    SparsityPattern spar_pattern_velocity, spar_pattern_pressure, spar_pattern_pres_vel;
    /*
      The actual matrices. The projection matrix never changes.
      For the velocity a part of the matrix never changes (neither in time nor for each component),
      namely the part related with the time derivative and the diffusion
      but, the advection part changes and so we need dim+1 matrices. One to store
      the constant part, and dim for each matrix at each iteration
    */
    SparseMatrix<double> vel_Laplace_plus_Mass, vel_it_matrix[dim], vel_Mass, vel_Laplace,
                         pres_Laplace, pres_Mass, pres_Diff[dim];
    /*
      We need to regularize the Laplace operator for the pressure
      for that we use a constraint
    */
    ConstraintMatrix pres_regularization;
    //// The solutions at times n and (n-1)
    Vector<double> pres_n, pres_n_minus_1, phi_n, phi_n_minus_1, u_n[dim], u_n_minus_1[dim],
    //// The time extrapolation of the velocity, used to make the semi-implicit time discretization
    //// of the advection term
                  u_star[dim],
    //// Right hand side for the component of the momentum equation
                  force[dim],
    //// Temporary arrays
    v_tmp, pres_tmp;
    //// The preconditioners
//     SparseILU<double> prec_velocity[dim];
    SparseDirectUMFPACK prec_velocity[dim], prec_mass, prec_pressure;

    // The convergence table
    ConvergenceTable convergence_table;

    // Exception we will throw if the time step is invalid
    DeclException2( ExcInvalidTimeStep, double, double, <<" The time step "<<arg1<<" is out of range."<<std::endl
                                                                <<" The permitted range is (0,"<<arg2<<"]");

    /*
      A diffusion step. The boolean parameter that it takes indicates if we want the preconditioners
      for each component of the velocity to be updated. Recall that the operator that these matrices
      represent is of the form
          1.5/dt Mass + 1/Re Laplace_h + ( u^*.D )
      so the dominating part is
          1.5/dt Mass + 1/Re Laplace_h
      and we can recycle the preconditioner for several steps to save computational time
    */
    inline void diffusion_step( const bool reinit_prec );
    /*
        A projection step. The parameter that it takes indicates if we want the preconditioner to be
        updated. Since this matrix never changes, ideally this should be done just once.
    */
    inline void projection_step( const bool reinit_prec );
    /*
      A pressure update step.
      Capable of recognizing which formulation of the method we are using.
      If it is the standard formulation, then the pressure update is done by
          p^{n+1} = p^{n} + phi^{n+1}
      which involves just adding numerical values.
      On the other hand, for the rotational formulation of the method the pressure update
      is done by the formula
          p^{n+1} = p^{n} + phi^{n+1} +1/Re div u^{n+1}
      which in any case involves the solution of a mass matrix problem.
      NOTE: Here we assume that the mesh is uniformly refined, and thus we do not use
      a preconditioner for this mass matrix problem
    */
    inline void update_pressure( const bool reinit_prec );
  private:

    // A bunch of data for the solvers and preconditioners
    // data solve velocity
    unsigned int vel_max_its, vel_Krylov_size, vel_off_diagonals, vel_update_prec;
    double vel_eps, vel_diag_strength;
    // data solve projection
    unsigned int proj_max_its, proj_off_diagonals;
    double proj_eps, proj_diag_strength;
    // data solve pressure update
    unsigned int pres_max_its, pres_off_diagonals;
    double pres_eps, pres_diag_strength;

    // Extremely simple function that computes the velocity extrapolant. Added for readability
    inline void interpolate_velocity();

    //////// initialization of constant wrt dt matrices
    inline void init_velocity_matrices();
    inline void init_pressure_matrices();
    inline void init_gradient_operator();

    // Assembly of the advection term
    inline void assemble_advection_term( const unsigned int d );

    // Solution of each component of the velocity on the diffusion step
    inline void threaded_solve( const unsigned int d );
    inline void plot_solution( const unsigned int step );
};



#endif

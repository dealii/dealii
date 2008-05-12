/* $Id$ */
/* Author: David Neckels, Boulder, Colorado, 2007, 2008 */

/*    $Id$       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 2007, 2008 by the deal.II authors and David Neckels */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

				 // @sect3{Include files}

                                 // First a standard set of deal.II
                                 // includes. Nothing special to comment on
                                 // here:
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/parameter_handler.h>
#include <base/function_parser.h>

#include <lac/vector.h>
#include <lac/sparse_matrix.h>
#include <lac/vector_memory.h>

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>
#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_in.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>

#include <fe/fe_values.h>
#include <fe/fe_system.h>
#include <fe/mapping_q1.h>
#include <fe/fe_q.h>

#include <numerics/data_out.h>
#include <numerics/vectors.h>
#include <numerics/solution_transfer.h>

                                 // Then, as mentioned in the introduction, we
                                 // use various Trilinos packages as linear
                                 // solvers as well as for automatic
                                 // differentiation. These are in the
                                 // following include files.
                                 //
                                 // In particular, Epetra is the basic
                                 // trilinos vector/matrix library and comes
                                 // with several header files pertaining to
                                 // individual aspects of it that will become
                                 // clear later on:
#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
                                 // Next, Teuchos is a Trilinos utility
                                 // library that is used to set parameters
                                 // within the Aztec solver library:
#include <Teuchos_ParameterList.hpp>

                                 // Aztec itself is the iterative solver
                                 // library:
#include <AztecOO.h>
#include <AztecOO_Operator.h>

                                 // Amesos is a direct solver package within
                                 // Trilinos:
#include <Amesos.h>

                                 // Finally, Sacado is the automatic
                                 // differentiation package, which is used to
                                 // find the Jacobian for a fully implicit
                                 // Newton iteration:
#include <Sacado.hpp>


				 // And this again is C++:
#include <iostream>
#include <fstream>
#include <vector>

				 // To end this section, introduce everythin
				 // in the dealii library into the current
				 // namespace:
using namespace dealii;


				 // @sect3{Flux function definition}

				 // Here we define the flux function for this
				 // particular system of conservation laws,
				 // the Euler equations for gas dynamics. We
				 // group all this into a structure that
				 // defines everything that has to do with the
				 // flux. All members of this structures are
				 // static, i.e. the structure has no actual
				 // state specified by instance member
				 // variables. The better way to do this,
				 // rather than a structure with all static
				 // members would be to use a namespace -- but
				 // namespaces can't be templatized and we
				 // want some of the member variables of the
				 // structure to depend on the space
				 // dimension, which we in our usual way
				 // introduce using a template parameter:
template <int dim>
struct EulerEquations
{
				   // First a few variables that
				   // describe the various components of our
				   // solution vector in a generic way. This
				   // includes the number of components in the
				   // system (Euler's equations have one entry
				   // for momenta in each spatial direction,
				   // plus the energy and density components,
				   // for a total of <code>dim+2</code>
				   // components), as well as functions that
				   // describe the index within the solution
				   // vector of the first momentum component,
				   // the density component, and the energy
				   // density component. Note that all these
				   // %numbers depend on the space dimension;
				   // defining them in a generic way (rather
				   // than by implicit convention) makes our
				   // code more flexible and makes it easier
				   // to later extend it, for example by
				   // adding more components to the equations.
    static const unsigned int n_components             = dim + 2;
    static const unsigned int first_momentum_component = 0;
    static const unsigned int density_component        = dim;
    static const unsigned int energy_component         = dim+1;

				   // Next, we define the gas
				   // constant. We will set it to 1.4
				   // in its definition immediately
				   // following the declaration of
				   // this class (unlike integer
				   // variables, like the ones above,
				   // static const floating point
				   // member variables cannot be
				   // initialized within the class
				   // declaration in C++). This value
				   // of 1.4 is representative of a
				   // gas that consists of molecules
				   // composed of two atoms, such as
				   // air which consists up to small
				   // traces almost entirely of $N_2$
				   // and $O_2$.
    static const double gas_gamma;

				     // We define the flux function
				     // $F(W)$ as one large matrix.
				     // Each row of this matrix
				     // represents a scalar
				     // conservation law for the
				     // component in that row.  The
				     // exact form of this matrix is
				     // given in the
				     // introduction. Note that we
				     // know the size of the matrix:
				     // it has as many rows as the
				     // system has components, and
				     // <code>dim</code> columns;
				     // rather than using a FullMatrix
				     // object for such a matrix
				     // (which has a variable number
				     // of rows and columns and must
				     // therefore allocate memory on
				     // the heap each time such a
				     // matrix is created), we use a
				     // rectangular array of numbers
				     // right away.
				     //
				     // We templatize the numerical
				     // type of the flux function so
				     // that we may use the automatic
				     // differentiation type here.
				     // The flux functions are defined
				     // in terms of the conserved
				     // variables $\rho w_0, \dots,
				     // \rho w_{d-1}, \rho, E$, so
				     // they do not look exactly like
				     // the Euler equations one is
				     // used to seeing.
    template <typename number>
    static
    void flux_matrix (const std::vector<number> &W,
		      number (&flux)[n_components][dim])
      {
					 // First compute the pressure
					 // that appears in the flux
					 // matrix, based on the
					 // energy density and the
					 // kinetic energy $\frac 12
					 // \rho |\mathbf v|^2 =
					 // \frac{|\rho \mathbf
					 // v|^2}{2\rho}$ (note that
					 // the independent variables
					 // contain the momentum
					 // components $\rho v_i$, not
					 // the velocities $v_i$):
	number kinetic_energy = 0;
	for (unsigned int d=0; d<dim; ++d)
	  kinetic_energy += W[first_momentum_component+d] *
			    W[first_momentum_component+d];
	kinetic_energy *= 1./(2 * W[density_component]);
	
	const number pressure = (gas_gamma-1.0)*(W[energy_component] - kinetic_energy);

					 // Then compute the first
					 // <code>dim</code> columns
					 // of the matrix that
					 // correspond to the momentum
					 // terms:
	for (unsigned int d=0; d<dim; ++d)
	  {
	    for (unsigned int e=0; e<dim; ++e)
	      flux[first_momentum_component+d][e] = W[first_momentum_component+d] *
						    W[first_momentum_component+e] /
						    W[density_component];
	  
	    flux[first_momentum_component+d][d] += pressure;
	  }
	
					 // Then the terms for the
					 // density (i.e. mass
					 // conservation), and,
					 // lastly, conservation of
					 // energy:
	for (unsigned int d=0; d<dim; ++d)
	  flux[density_component][d] = W[first_momentum_component+d]; 

	for (unsigned int d=0; d<dim; ++d)
	  flux[energy_component][d] = W[first_momentum_component+d] /
				      W[density_component] *
				      (W[energy_component] + pressure);
      }


				     // On the boundaries of the
				     // domain and across hanging
				     // nodes we use a numerical flux
				     // function to enforce boundary
				     // conditions.  This routine is
				     // the basic Lax-Friedrich's flux
				     // with a stabilization parameter
				     // $\alpha$. It's form has also
				     // been given already in the
				     // introduction:
    template <typename number>
    static
    void numerical_normal_flux(const Point<dim>          &normal,
			       const std::vector<number> &Wplus,
			       const std::vector<number> &Wminus,
			       const double alpha,
			       Sacado::Fad::DFad<double> (&normal_flux)[n_components])
      {
	Sacado::Fad::DFad<double> iflux[n_components][dim];
	Sacado::Fad::DFad<double> oflux[n_components][dim];
	  
	flux_matrix(Wplus, iflux);
	flux_matrix(Wminus, oflux);
	  
	for (unsigned int di=0; di<n_components; ++di)
	  {
	    normal_flux[di] = 0;
	    for (unsigned int d=0; d<dim; ++d) 
	      normal_flux[di] += 0.5*(iflux[di][d] + oflux[di][d]) * normal(d);
	      
	    normal_flux[di] += 0.5*alpha*(Wplus[di] - Wminus[di]);
	  }
    }
};
    

template <int dim>
const double EulerEquations<dim>::gas_gamma = 1.4;


				 // @sect3{Conservation Law class}

				 // Here we define a Conservation Law
				 // class that helps group operations
				 // and data for our Euler equations
				 // into a manageable entity.  Member
				 // functions will be described as
				 // their definitions appear.
template <int dim>
class ConsLaw
{
  public:
    ConsLaw ();
    ~ConsLaw ();

    void run ();
    void declare_parameters();
    void load_parameters(const char *);
    
  private:
    void build_fe();
    void setup_system ();
    void initialize_system ();
    void assemble_system (double &res_norm);
    void solve (Vector<double> &solution, int &, double &);
    void refine_grid ();
    void output_results (const unsigned int cycle) const;
    void initialize();
    void zero_matrix();
    void estimate();
    void postprocess();
    void compute_predictor();
    
    Triangulation<dim>   triangulation;
    const MappingQ1<dim> mapping;
    
    
    FESystem<dim>        *fe_ptr;

    DoFHandler<dim>      dof_handler;

    SparsityPattern      sparsity_pattern;
    const QGauss<dim>   quadrature;
    const QGauss<dim-1> face_quadrature;
    
                                     // The actual solution to the Euler equation
    Vector<double>       solution;
                                     // The current value of the solution during the Newton iteration
    Vector<double>       nlsolution;
                                     // An estimate of the next time value; used for adaptivity and as a
                                     // guess for the next Newton iteration.
    Vector<double>       predictor;
                                     // Values after post-processing (used to output the physical variables).
    Vector<double>       ppsolution;
                                     // The solution to the linear problem during the Newton iteration
    Vector<double>       dsolution;
    Vector<double>       right_hand_side;

    Epetra_SerialComm    communicator;
    
  public:

    void assemble_cell_term(const FEValues<dim>& fe_v,
                            std::vector<unsigned int> &dofs,
                            unsigned int cell_no
    );
    
    void assemble_face_term(
      int face_no,
      const FEFaceValuesBase<dim>& fe_v,
      const FEFaceValuesBase<dim>& fe_v_neighbor,
      std::vector<unsigned int> &dofs,
      std::vector<unsigned int> &dofs_neighbor,
      int boundary = -1
    );

  private:
				     // T = current time, dT = time step, TF = final time.
    double T, dT, TF;
    double face_diameter;
    double cell_diameter;
				     // An object to handle parsing the input deck.
    ParameterHandler prm;
				     // Name of the mesh to read in.
    string mesh;
    FunctionParser<dim> initial_conditions;

				     // Enums for the various supported boundary conditions.
    typedef enum {INFLOW_BC = 1, OUTFLOW_BC=2, NO_PENETRATION_BC=3, PRESSURE_BC=4} bc_type;

				     // For each boundary we store a map from boundary # to the type
				     // of boundary condition.  If the boundary condition is prescribed,
				     // we store a pointer to a function object that will hold the expression
				     // for that boundary condition.
    typedef typename std::map<unsigned int, std::pair<std::vector<bc_type>, Function<dim>*> > bdry_map_type;
    bdry_map_type bdry_map;

    void add_boundary(unsigned int bd, std::vector<bc_type>& flags, Function<dim> *bf);

				     // An object to store parameter information about the Aztec solver.
    typedef struct {
	int LIN_OUTPUT;
	typedef enum { GMRES = 0, DIRECT = 1} solver_type;
	solver_type SOLVER;
	typedef enum { QUIET = 0, VERBOSE = 1 } output_type;
	output_type OUTPUT;
					 // Linear residual tolerance.
	double RES;
	int MAX_ITERS;
					 // We use the ILUT preconditioner.  This is similar
					 // to the ILU.  FILL is the number of extra entries
					 // to add when forming the ILU decomposition.
	double ILUT_FILL;
					 // When forming the preconditioner, for certain problems
					 // bad conditioning (or just bad luck) can cause the
					 // preconditioner to be very poorly conditioned.  Hence
					 // it can help to add diagonal perturbations to the
					 // original matrix and form the preconditioner for this
					 // slightly better matrix.  ATOL is an absolute perturbation
					 // that is added to the diagonal before forming the
					 // prec, and RTOL is a scaling factor $rtol >= 1$.
	double ILUT_ATOL;
	double ILUT_RTOL;
					 // The ILUT will drop any values that have magnitude less
					 // than this value.  This is a way to
					 // manage the amount of memory used by this preconditioner.
	double ILUT_DROP;
    } solver_params_type;

    solver_params_type solver_params;

				     // Some refinement parameters.
    typedef struct {
	typedef enum { NONE = 0, FIXED_NUMBER = 1, SHOCK = 2} refine_type;
	double high_frac;
	double low_frac;
	refine_type refine;
	double high_frac_sav;
	double max_cells;
	double shock_val;
	double shock_levels;
    } refinement_params_type;

    refinement_params_type refinement_params;

				     // The user can set the stabilization parameter $\alpha$ 
				     // in the Lax-Friedrich's flux.  
    typedef struct {
	typedef enum {CONSTANT=1,MESH=2} LF_stab_type;
	LF_stab_type LF_stab;
	double LF_stab_value;
    } flux_params_type;

    flux_params_type flux_params;

    bool is_stationary;

				     // Power for the mesh stabilization term.
    double diffusion_power;
    double gravity;
				     // If true, we output the squared gradient of the
				     // density instead of density.  Using this one can
				     // create shock plots.
    bool schlieren_plot;
				     // How often to create an output file.
    double output_step;

    Epetra_Map         *Map;
    Epetra_CrsMatrix   *Matrix;
    Vector<double>      indicator;
 
				     // Crank-Nicolson value
    const double        theta; 

};


				 // Asign a row of the conservation law a specified
				 // boundary type and (possibly) function.
template <int dim>
void ConsLaw<dim>::add_boundary(unsigned int bd,
				std::vector<bc_type> &flags, Function<dim> *bf) {

  std::pair<std::vector<bc_type>, Function<dim> *> entry(flags, bf);
  bdry_map[bd] = entry;
}


				 // Apply the initialial condition.  Simultaneously
				 // initialize the non-linear solution.
template <int dim>
void ConsLaw<dim>::initialize() {
  VectorTools::interpolate(dof_handler,
                           initial_conditions, solution);
  nlsolution = solution;
}

				 // @sect3{Assembly}
				 // @sect4{%Function: assemble_cell_term}
				 //
                                 // Assembles the cell term, adding minus the residual
                                 // to the right hand side, and adding in the Jacobian
                                 // contributions.
template <int dim>
void ConsLaw<dim>::assemble_cell_term(
  const FEValues<dim> &fe_v,
  std::vector<unsigned int> &dofs,
  unsigned int /*cell_no*/
) 
{
  unsigned int dofs_per_cell = fe_v.dofs_per_cell;
  unsigned int n_q_points = fe_v.n_quadrature_points;

				   // We will define the dofs on this cell in these fad variables.
  std::vector<Sacado::Fad::DFad<double> > DOF(dofs_per_cell);

				   // Values of the conservative variables at the quadrature points.
  std::vector<std::vector<Sacado::Fad::DFad<double> > > W (n_q_points,
					    std::vector<Sacado::Fad::DFad<double> >(EulerEquations<dim>::n_components));

				   // Values at the last time step of the conservative variables.
				   // Note that these do not use fad variables, since they do
				   // not depend on the 'variables to be sought'=DOFS.
  std::vector<std::vector<double > > Wl (n_q_points,
					 std::vector<double >(EulerEquations<dim>::n_components));

				   // Here we will hold the averaged values of the conservative
				   // variables that we will linearize around (cn=Crank Nicholson).
  std::vector<std::vector<Sacado::Fad::DFad<double> > > Wcn (n_q_points,
					      std::vector<Sacado::Fad::DFad<double> >(EulerEquations<dim>::n_components));

				   // Gradients of the current variables.  It is a
				   // bit of a shame that we have to compute these; we almost don't.
				   // The nice thing about a simple conservation law is that the
				   // the flux doesn't generally involve any gradients.  We do
				   // need these, however, for the diffusion stabilization. 
  std::vector<std::vector<std::vector<Sacado::Fad::DFad<double> > > > Wgrads (n_q_points,
							      std::vector<std::vector<Sacado::Fad::DFad<double> > >(EulerEquations<dim>::n_components,
												    std::vector<Sacado::Fad::DFad<double> >(dim)));


				   // Here is the magical point where we declare a subset
				   // of the fad variables as degrees of freedom.  All 
				   // calculations that reference these variables (either
				   // directly or indirectly) will accumulate sensitivies
				   // with respect to these dofs.
  for (unsigned int in = 0; in < dofs_per_cell; in++) {
    DOF[in] = nlsolution(dofs[in]);
    DOF[in].diff(in, dofs_per_cell);
  }

				   // Here we compute the shape function values and gradients
				   // at the quadrature points.  Ideally, we could call into 
				   // something like get_function_values, get_function_grads,
				   // but since we don't want to make the entire solution vector
				   // fad types, only the local cell variables, we explicitly
				   // code this loop;
  for (unsigned int q = 0; q < n_q_points; q++) {
    for (unsigned int di = 0; di < EulerEquations<dim>::n_components; di++) {
      W[q][di] = 0;
      Wl[q][di] = 0;
      Wcn[q][di] = 0;
      for (unsigned int d = 0; d < dim; d++) {
        Wgrads[q][di][d] = 0;
      }
    }
    for (unsigned int sf = 0; sf < dofs_per_cell; sf++) {
      int di = fe_v.get_fe().system_to_component_index(sf).first;
      W[q][di] +=
	DOF[sf]*fe_v.shape_value_component(sf, q, di);
      Wl[q][di] +=
	solution(dofs[sf])*fe_v.shape_value_component(sf, q, di);
      Wcn[q][di] +=
	(theta*DOF[sf]+(1-theta)*solution(dofs[sf]))*fe_v.shape_value_component(sf, q, di);

      for (unsigned int d = 0; d < dim; d++) {
	Wgrads[q][di][d] += DOF[sf]*
			    fe_v.shape_grad_component(sf, q, di)[d];
      } // for d

    }

  } // for q

                                   // Gather the flux values for all components at
                                   // all of the quadrature points.  This also
                                   // computes the matrix of sensitivities.  Perhaps
                                   // this could be done in a better way, since this
                                   // could be a rather large object, but for now it 
                                   // seems to work just fine.
  typedef Sacado::Fad::DFad<double> FluxMatrix[EulerEquations<dim>::n_components][dim];
  FluxMatrix *flux = new FluxMatrix[n_q_points];
  
  for (unsigned int q=0; q < n_q_points; ++q)
    EulerEquations<dim>::flux_matrix(Wcn[q], flux[q]);
  

				   // We now have all of the function values/grads/fluxes,
				   // so perform the assembly.  We have an outer loop
				   // through the components of the system, and an
				   // inner loop over the quadrature points, where we
				   // accumulate contributions to the ith residual.
				   //
				   // We initialy sum all contributions of the residual
				   // in the positive sense, so that we don't need to
				   // negative the Jacobian entries.  Then, when we sum
				   // into the <code> right_hand_side </code> vector,
				   // we negate this residual.
  for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i) 
    {
				       // Find which component this dof contributes to.
      const unsigned int
	component_i = fe_v.get_fe().system_to_component_index(i).first;

				       // The residual for each row (i) will be accumulating 
				       // into this fad variable.  At the end of the assembly
				       // for this row, we will query for the sensitivities
				       // to this variable and add them into the Jacobian.
      Sacado::Fad::DFad<double> F_i;

      for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
	{
					   // Integrate the flux times gradient of the test function
	  for (unsigned int d=0; d<dim; d++) 
	    F_i -= flux[point][component_i][d] *
		   fe_v.shape_grad_component(i, point, component_i)[d] *
		   fe_v.JxW(point);

					   // The mass term (if the simulation is non-stationary).
	  if (!is_stationary)
	    F_i += 1.0/dT*(W[point][component_i] - Wl[point][component_i]) *
		   fe_v.shape_value_component(i, point, component_i) *
		   fe_v.JxW(point);

					   // Stabilization (cell wise diffusion)
	  for (unsigned int d = 0; d < dim; d++)
	    F_i += 1.0*std::pow(cell_diameter, diffusion_power) *
		   fe_v.shape_grad_component(i, point, component_i)[d] *
		   Wgrads[point][component_i][d] *
		   fe_v.JxW(point);
          
					   // The gravity component only enters into the energy 
					   // equation and into the vertical component of the 
					   // velocity.
	  if (component_i == dim - 1)
	    F_i += gravity *
		   Wcn[point][EulerEquations<dim>::density_component] *
		   fe_v.shape_value_component(i,point, component_i) *
		   fe_v.JxW(point);
	  else if (component_i == EulerEquations<dim>::energy_component)
	    F_i += gravity *
		   Wcn[point][EulerEquations<dim>::density_component] *
		   Wcn[point][dim-1] *
		   fe_v.shape_value_component(i,point, component_i) *
		   fe_v.JxW(point);
	}

				       // Here we gain access to the array of sensitivities
				       // of the residual.  We then sum these into the
				       // Epetra matrix.
      double *values = &(F_i.fastAccessDx(0));
      Matrix->SumIntoGlobalValues(dofs[i],
				  dofs_per_cell,
				  values,
				  reinterpret_cast<int*>(&dofs[0]));
 
				       // Add minus the residual to the right hand side.
      right_hand_side(dofs[i]) -= F_i.val();
    }

  delete[] flux;
}
				 // @sect4{%Function: assemble_face_term}
				 // These are either
				 // boundary terms or terms across differing 
				 // levels of refinement.  In the first case,
				 // fe_v==fe_v_neighbor and dofs==dofs_neighbor.
				 // The int boundary < 0 if not at a boundary,
				 // otherwise it is the boundary indicator.
template <int dim>
void ConsLaw<dim>::assemble_face_term(
  int face_no,
  const FEFaceValuesBase<dim>& fe_v,
  const FEFaceValuesBase<dim>& fe_v_neighbor,      
  std::vector<unsigned int> &dofs,
  std::vector<unsigned int> &dofs_neighbor,
  int boundary
) 
{
  Sacado::Fad::DFad<double> F_i;
  const unsigned int n_q_points = fe_v.n_quadrature_points;
  const unsigned int dofs_per_cell = fe_v.get_fe().dofs_per_cell;
  const unsigned int ndofs_per_cell = fe_v_neighbor.get_fe().dofs_per_cell;
  Assert(dofs_per_cell == ndofs_per_cell,
	 ExcDimensionMismatch(dofs_per_cell, ndofs_per_cell));

				   // As above, the fad degrees of freedom
  std::vector<Sacado::Fad::DFad<double> > DOF(dofs_per_cell+ndofs_per_cell);

				   // The conservative variables for this cell,
				   // and for 
  std::vector<std::vector<Sacado::Fad::DFad<double> > > Wplus (n_q_points,
						std::vector<Sacado::Fad::DFad<double> >(EulerEquations<dim>::n_components));
  std::vector<std::vector<Sacado::Fad::DFad<double> > > Wminus (n_q_points,
						 std::vector<Sacado::Fad::DFad<double> >(EulerEquations<dim>::n_components));


  const std::vector<Point<dim> > &normals = fe_v.get_normal_vectors ();


				   // If we are at a boundary, then dofs_neighbor are
				   // the same as dofs, so we do not want to duplicate them.
				   // If there is a neighbor cell, then we want to include 
				   // them.
  int ndofs = (boundary < 0 ? dofs_per_cell + ndofs_per_cell : dofs_per_cell);
				   // Set the local DOFS.
  for (unsigned int in = 0; in < dofs_per_cell; in++) {
    DOF[in] = nlsolution(dofs[in]);
    DOF[in].diff(in, ndofs);
  }
				   // If present, set the neighbor dofs.
  if (boundary < 0)
    for (unsigned int in = 0; in < ndofs_per_cell; in++) {
      DOF[in+dofs_per_cell] = nlsolution(dofs_neighbor[in]);
      DOF[in+dofs_per_cell].diff(in+dofs_per_cell, ndofs);
    }

				   // Set the values of the local conservative variables.
				   // Initialize all variables to zero.
  for (unsigned int q = 0; q < n_q_points; q++) {
    for (unsigned int di = 0; di < EulerEquations<dim>::n_components; di++) {
      Wplus[q][di] = 0;
      Wminus[q][di] = 0;
    }
    for (unsigned int sf = 0; sf < dofs_per_cell; sf++) {
      int di = fe_v.get_fe().system_to_component_index(sf).first;
      Wplus[q][di] +=
	(theta*DOF[sf]+(1.0-theta)*solution(dofs[sf]))*fe_v.shape_value_component(sf, q, di);
    }


				     // If there is a cell across, then initialize
				     // the exterior trace as a function of the other
				     // cell degrees of freedom.
    if (boundary < 0) {
      for (unsigned int sf = 0; sf < ndofs_per_cell; sf++) {
	int di = fe_v_neighbor.get_fe().system_to_component_index(sf).first;
	Wminus[q][di] +=
	  (theta*DOF[sf+dofs_per_cell]+(1.0-theta)*solution(dofs_neighbor[sf]))*
	  fe_v_neighbor.shape_value_component(sf, q, di);
      }
    } 
  } // for q

				   // If this is a boundary, then the values of $W^-$ will
				   // be either functions of $W^+$, or they will be prescribed.
				   // This switch sets them appropriately.  Since we are
				   // using fad variables here, sensitivities will be updated 
				   // appropriately.  These sensitivities would be tremendously
				   // difficult to manage without fad!!!
  if (boundary >= 0) {
				     // Get the boundary descriptor.
    typename bdry_map_type::iterator bme = bdry_map.find(boundary);
    assert(bme != bdry_map.end());

				     // Evaluate the function object.  This is a bit
				     // tricky; a given boundary might have both prescribed
				     // and implicit values.  If a particular component is not
				     // prescribed, the values evaluate to zero and are
				     // ignored, below.
    std::vector<Vector<double> > bvals(n_q_points, Vector<double>(EulerEquations<dim>::n_components));
    bme->second.second->vector_value_list(fe_v.get_quadrature_points(), bvals);

				     // We loop the quadrature points, and we treat each
				     // component individualy.
    for (unsigned int q = 0; q < n_q_points; q++) {
      for (unsigned int di = 0; di < EulerEquations<dim>::n_components; di++) {

					 // An inflow/dirichlet type of boundary condition
        if (bme->second.first[di] == INFLOW_BC) {
          Wminus[q][di] = bvals[q](di);
        } else if (bme->second.first[di] == PRESSURE_BC) {
					   // A prescribed pressure boundary condition.  This boundary
					   // condition is complicated by the fact that even though
					   // the pressure is prescribed, we really are setting
					   // the energy index here, which will depend on velocity
					   // and pressure. So even though this seems like a dirichlet
					   // type boundary condition, we get sensitivities of
					   // energy to velocity and density (unless these
					   // are also prescribed.
          Sacado::Fad::DFad<double> rho_vel_sqr = 0;
          Sacado::Fad::DFad<double> dens;
          
          dens = bme->second.first[EulerEquations<dim>::density_component] == INFLOW_BC ? bvals[q](EulerEquations<dim>::density_component) :
                 Wplus[q][EulerEquations<dim>::density_component];

          for (unsigned int d=0; d < dim; d++) {
            if (bme->second.first[d] == INFLOW_BC)
              rho_vel_sqr += bvals[q](d)*bvals[q](d);
            else
              rho_vel_sqr += Wplus[q][d]*Wplus[q][d];
          }
          rho_vel_sqr /= dens;
					   // Finally set the energy value as determined by the
					   // prescribed pressure and the other variables.
          Wminus[q][di] = bvals[q](di)/(EulerEquations<dim>::gas_gamma-1.0) +
			  0.5*rho_vel_sqr;

        } else if (bme->second.first[di] == OUTFLOW_BC) {
					   // A free/outflow boundary, very simple.
          Wminus[q][di] = Wplus[q][di];

        } else { 
					   // We must be at a no-penetration boundary.  We
					   // prescribe the velocity (we are dealing with a
					   // particular component here so that the average
					   // of the velocities is orthogonal to the surface
					   // normal.  This creates sensitivies of across
					   // the velocity components.
          Sacado::Fad::DFad<double> vdotn = 0;
          for (unsigned int d = 0; d < dim; d++) {
            vdotn += Wplus[q][d]*normals[q](d);
          }

          Wminus[q][di] = Wplus[q][di] - 2.0*vdotn*normals[q](di);
        }
      }
    } // for q
  } // b>= 0
   
				   // Determine the Lax-Friedrich's stability parameter,
				   // and evaluate the numerical flux function at the quadrature points
  typedef Sacado::Fad::DFad<double> NormalFlux[EulerEquations<dim>::n_components];
  NormalFlux *normal_fluxes = new NormalFlux[n_q_points];

  double alpha = 1;

  switch(flux_params.LF_stab) {
    case flux_params_type::CONSTANT:
	  alpha = flux_params.LF_stab_value;
	  break;
    case flux_params_type::MESH:
	  alpha = face_diameter/(2.0*dT);
	  break;
  }

  for (unsigned int q=0; q<n_q_points; ++q)
    EulerEquations<dim>::numerical_normal_flux(normals[q], Wplus[q], Wminus[q], alpha,
					       normal_fluxes[q]);

				   // Now assemble the face term
  for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
    {
      if (!fe_v.get_fe().has_support_on_face(i, face_no))
	continue;
      
      F_i = 0;
      for (unsigned int point=0; point<n_q_points; ++point)
	{
	  const unsigned int
	    component_i = fe_v.get_fe().system_to_component_index(i).first;
	  
	  F_i += normal_fluxes[point][component_i] *
		 fe_v.shape_value_component(i, point, component_i) *
		 fe_v.JxW(point);
	} 

				       // Retrieve a pointer to the jacobian.
      double *values = &(F_i.fastAccessDx(0));
      Assert (values != 0, ExcInternalError());

				       // Update the matrix.  Depending on whether there
				       // is/isn't a neighboring cell, we add more/less
				       // entries.
      Matrix->SumIntoGlobalValues(dofs[i],
				  dofs_per_cell, &values[0], reinterpret_cast<int*>(&dofs[0]));
      if (boundary < 0) {
	Matrix->SumIntoGlobalValues(dofs[i],
				    dofs_per_cell, &values[dofs_per_cell], reinterpret_cast<int*>(&dofs_neighbor[0]));
      }

				       // And add into the residual
      right_hand_side(dofs[i]) -= F_i.val();
    }

  delete[] normal_fluxes;
}
                                 // @sect4{Assembling the whole system}
                                 // Now we put all of the assembly pieces together
                                 // in a routine that dispatches the correct
                                 // piece for each cell/face.  We keep track of
                                 // the norm of the resdual for the Newton iteration.
template <int dim>
void ConsLaw<dim>::assemble_system (double &res_norm) 
{
  FESystem<dim> &fe = *fe_ptr;
  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

				   // We track the dofs on this cell and (if necessary)
				   // the adjacent cell.
  std::vector<unsigned int> dofs (dofs_per_cell);
  std::vector<unsigned int> dofs_neighbor (dofs_per_cell);

				   // First we create the
				   // ``UpdateFlags'' for the
				   // ``FEValues'' and the
				   // ``FEFaceValues'' objects.
  UpdateFlags update_flags = update_values
			     | update_gradients
			     | update_q_points
			     | update_JxW_values;

				   // Note, that on faces we do not
				   // need gradients but we need
				   // normal vectors.
  UpdateFlags face_update_flags = update_values
				  | update_q_points
				  | update_JxW_values
				  | update_normal_vectors;
  
				   // On the neighboring cell we only
				   // need the shape values. Given a
				   // specific face, the quadrature
				   // points and `JxW values' are the
				   // same as for the current cells,
				   // the normal vectors are known to
				   // be the negative of the normal
				   // vectors of the current cell.
  UpdateFlags neighbor_face_update_flags = update_values;
   
				   // Then we create the ``FEValues''
				   // object. Note, that since version
				   // 3.2.0 of deal.II the constructor
				   // of this class takes a
				   // ``Mapping'' object as first
				   // argument. Although the
				   // constructor without ``Mapping''
				   // argument is still supported it
				   // is recommended to use the new
				   // constructor. This reduces the
				   // effect of `hidden magic' (the
				   // old constructor implicitely
				   // assumes a ``MappingQ1'' mapping)
				   // and makes it easier to change
				   // the mapping object later.
  FEValues<dim> fe_v (
    mapping, fe, quadrature, update_flags);
  
				   // Similarly we create the
				   // ``FEFaceValues'' and
				   // ``FESubfaceValues'' objects for
				   // both, the current and the
				   // neighboring cell. Within the
				   // following nested loop over all
				   // cells and all faces of the cell
				   // they will be reinited to the
				   // current cell and the face (and
				   // subface) number.
  FEFaceValues<dim> fe_v_face (
    mapping, fe, face_quadrature, face_update_flags);
  FESubfaceValues<dim> fe_v_subface (
    mapping, fe, face_quadrature, face_update_flags);
  FEFaceValues<dim> fe_v_face_neighbor (
    mapping, fe, face_quadrature, neighbor_face_update_flags);
  FESubfaceValues<dim> fe_v_subface_neighbor (
    mapping, fe, face_quadrature, neighbor_face_update_flags);

				   // Furthermore we need some cell
				   // iterators.
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

				   // Now we start the loop over all
				   // active cells.
  unsigned int cell_no = 0;
  for (;cell!=endc; ++cell, ++cell_no) 
    {
      
				       // Now we reinit the ``FEValues''
				       // object for the current cell
      fe_v.reinit (cell);

                                       // Collect the local dofs and
                                       // asssemble the cell term.
      cell->get_dof_indices (dofs);

      cell_diameter = cell->diameter();

      assemble_cell_term(fe_v,
                         dofs,
                         cell_no);

                                       // We use the DG style loop through faces
                                       // to determine if we need to apply a
                                       // 'hanging node' flux calculation or a boundary
                                       // computation.
      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
	{
					   // First we set the face
					   // iterator
	  typename DoFHandler<dim>::face_iterator face=cell->face(face_no);
          face_diameter = face->diameter();
	  
	  if (face->at_boundary())
	    {
					       // We reinit the
					       // ``FEFaceValues''
					       // object to the
					       // current face
	      fe_v_face.reinit (cell, face_no);

					       // and assemble the
					       // corresponding face
					       // terms.  We send the same
                                               // fe_v and dofs as described
                                               // in the assembly routine.
	      assemble_face_term(
		face_no, fe_v_face,
		fe_v_face,
		dofs,
		dofs,
		face->boundary_indicator());
	    }
	  else
	    {
					       // Now we are not on
					       // the boundary of the
					       // domain, therefore
					       // there must exist a
					       // neighboring cell.
	      typename DoFHandler<dim>::cell_iterator neighbor=
		cell->neighbor(face_no);;

	      if (face->has_children())
		{
						   // case I: This cell refined compared to neighbor

		  const unsigned int neighbor2=
		    cell->neighbor_of_neighbor(face_no);
		  
		  
						   // We loop over
						   // subfaces
		  for (unsigned int subface_no=0;
		       subface_no<GeometryInfo<dim>::subfaces_per_face;
		       ++subface_no)
		    {
		      typename DoFHandler<dim>::active_cell_iterator
                        neighbor_child
                        = cell->neighbor_child_on_subface (face_no, subface_no);

                      face_diameter = neighbor_child->diameter();  // working on subface
		      
		      Assert (neighbor_child->face(neighbor2) == face->child(subface_no),
			      ExcInternalError());
		      Assert (!neighbor_child->has_children(), ExcInternalError());

		      fe_v_subface.reinit (cell, face_no, subface_no);
		      fe_v_face_neighbor.reinit (neighbor_child, neighbor2);
		      neighbor_child->get_dof_indices (dofs_neighbor);

						       // Assemble as if we are working with
						       // a DG element.
		      assemble_face_term(
			face_no, fe_v_subface,
			fe_v_face_neighbor,
			dofs,
			dofs_neighbor);
		      
		    }
						   // End of ``if
						   // (face->has_children())''
		}
	      else
		{
						   // We have no children, but 
						   // the neighbor cell may be refine
						   // compared to use
		  neighbor->get_dof_indices (dofs_neighbor);
		  if (neighbor->level() != cell->level()) 
		    {
						       // case II: This is refined compared to neighbor
		      Assert(neighbor->level() < cell->level(), ExcInternalError());
		      const std::pair<unsigned int, unsigned int> faceno_subfaceno=
			cell->neighbor_of_coarser_neighbor(face_no);
		      const unsigned int neighbor_face_no=faceno_subfaceno.first,
				      neighbor_subface_no=faceno_subfaceno.second;

		      Assert (neighbor->neighbor_child_on_subface (neighbor_face_no,
                                                                   neighbor_subface_no)
                              == cell,
                              ExcInternalError());

						       // Reinit the
						       // appropriate
						       // ``FEFaceValues''
						       // and assemble
						       // the face
						       // terms.
		      fe_v_face.reinit (cell, face_no);
		      fe_v_subface_neighbor.reinit (neighbor, neighbor_face_no,
						    neighbor_subface_no);
		      
		      assemble_face_term(
			face_no, fe_v_face,
			fe_v_subface_neighbor,
			dofs,
			dofs_neighbor);

		    }

		} 
					       // End of ``face not at boundary'':
	    }
					   // End of loop over all faces:
	} 
      
				       // End iteration through cells.
    } 

				   // Notify Epetra that the matrix is done.
  Matrix->FillComplete();

				   // Compute the nonlinear residual.
  res_norm = right_hand_side.l2_norm();
    
}

				 // Create a conservation law with some defaults.
template <int dim>
ConsLaw<dim>::ConsLaw ()
		:
		mapping (),
                fe_ptr(NULL),
		dof_handler (triangulation),
		quadrature (2),
		face_quadrature (2),
                T(0),
                dT(0.05),
                TF(10),
		initial_conditions (EulerEquations<dim>::n_components),
                is_stationary(false),
                Map(NULL),
                Matrix(NULL),
                theta(0.5) 
{}

				 // At one time this example could work for both DG and
				 // continuous finite elements.  The choice was made here.
template <int dim>
void ConsLaw<dim>::build_fe() {
  fe_ptr = new FESystem<dim>(FE_Q<dim>(1), EulerEquations<dim>::n_components);
}

				 // Bye bye Conservation law.
template <int dim>
ConsLaw<dim>::~ConsLaw () 
{
  dof_handler.clear ();
  delete fe_ptr;
}

				 // @sect3{Initialize System}
				 // Sizes all of the vectors and sets up the
				 // sparsity patter.  This function is called at
				 // the very beginning of a simulation.  The function
				 // <code> setup_system </code> repeats some of these
				 // chores and is called after adaptivity in leiu
				 // of this function.
template <int dim>
void ConsLaw<dim>::initialize_system ()
{
				   // First we need to distribute the
				   // DoFs.
  dof_handler.clear();
  dof_handler.distribute_dofs (*fe_ptr);
  
                                   // Size all of the fields.
  solution.reinit (dof_handler.n_dofs());
  nlsolution.reinit (dof_handler.n_dofs());
  predictor.reinit (dof_handler.n_dofs());
  ppsolution.reinit (dof_handler.n_dofs());
  dsolution.reinit (dof_handler.n_dofs());
  right_hand_side.reinit (dof_handler.n_dofs());
  indicator.reinit(triangulation.n_active_cells());
}

				 // @sect3{Setup System}
				 // We call this function to build the sparsity
				 // and the matrix.
template <int dim>
void ConsLaw<dim>::setup_system ()
{

				   // The DoFs of a cell are coupled
				   // with all DoFs of all neighboring
				   // cells.  Therefore the maximum
				   // number of matrix entries per row
				   // is needed when all neighbors of
				   // a cell are once more refined
				   // than the cell under
				   // consideration.
  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   (GeometryInfo<dim>::faces_per_cell
			    *GeometryInfo<dim>::subfaces_per_face+1)*fe_ptr->dofs_per_cell);
  
                                   // Since the continuous sparsity pattern is
                                   // a subset of the DG one, and since we need
                                   // the DG terms for handling hanging nodes, we use
                                   // the flux pattern.
  DoFTools::make_flux_sparsity_pattern (dof_handler, sparsity_pattern);
  
  sparsity_pattern.compress();
  
                                   // Rebuild the map.  In serial this doesn't do much,
                                   // but is needed.  In parallel, this would desribe
                                   // the parallel dof layout.
  if (Map) delete Map;
  Map  = new Epetra_Map(dof_handler.n_dofs(), 0, communicator);

                                   // Epetra can build a more efficient matrix if
                                   // one knows ahead of time the maximum number of
                                   // columns in any row entry
  std::vector<int> row_lengths (dof_handler.n_dofs());
  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    row_lengths[i] = sparsity_pattern.row_length (i);

				   // Now we build the matrix, using
				   // the constructor that optimizes
				   // with the existing lengths per row
				   // variable.
  if (Matrix) delete Matrix;
  Matrix = new Epetra_CrsMatrix(Copy, *Map, &row_lengths[0], true);

				   // We add the sparsity pattern to the matrix by
				   // inserting zeros.
  const unsigned int max_nonzero_entries = *std::max_element (row_lengths.begin(),
							      row_lengths.end());
  std::vector<double> vals(max_nonzero_entries, 0);
  std::vector<int> row_indices(max_nonzero_entries);
 
  unsigned int cur_row = 0;
  unsigned int cur_col = 0;
  for (SparsityPattern::iterator s_i = sparsity_pattern.begin(); 
       s_i != sparsity_pattern.end(); s_i++) {
    if (s_i->row() != cur_row) {
      Matrix->InsertGlobalValues(cur_row, cur_col, &vals[0], &row_indices[0]);
      cur_col = 0;
      cur_row = s_i->row();
    }
    row_indices[cur_col++] = s_i->column();
  }
				   // The last row.
  Matrix->InsertGlobalValues(cur_row, cur_col, &vals[0], &row_indices[0]);

				   // Epetra requires this function after building or
				   // filling a matrix.  It typically does some parallel
				   // bookeeping; perhaps more.
  Matrix->FillComplete();

}

                                 // @sect3{Solving the linear system}
                                 // Actually solve the linear system, using either
                                 // Aztec or Amesos.
template <int dim>
void ConsLaw<dim>::solve (Vector<double> &dsolution, int &niter, double &lin_residual) 
{

				   // We must hand the solvers Epetra vectors.
				   // Luckily, they support the concept of a 
				   // 'view', so we just send in a pointer to our
				   // dealii vectors.
  Epetra_Vector x(View, *Map, dsolution.begin());
  Epetra_Vector b(View, *Map, right_hand_side.begin());

				   // The Direct option selects the Amesos solver.
  if (solver_params.SOLVER == solver_params_type::DIRECT) {
   
				     // Setup for solving with
				     // Amesos. Other solvers are
				     // available and may be selected by
				     // changing th string given to the
				     // <code>Create</code> function.
    Epetra_LinearProblem prob;
    prob.SetOperator(Matrix);
    Amesos_BaseSolver *solver = Amesos().Create ("Amesos_Klu", prob);

    Assert (solver != NULL, ExcInternalError());

				     // There are two parts to the direct solve.
				     // As I understand, the symbolic part figures
				     // out the sparsity patterns, and then the
				     // numerical part actually performs Gaussian
				     // elimination or whatever the approach is.
    if (solver_params.OUTPUT == solver_params_type::VERBOSE)
      std::cout << "Starting Symbolic fact\n" << std::flush;

    solver->SymbolicFactorization();

    if (solver_params.OUTPUT == solver_params_type::VERBOSE)
      std::cout << "Starting Numeric fact\n" << std::flush;

    solver->NumericFactorization();

    
				     // Define the linear problem by setting the
				     // right hand and left hand sides.
    prob.SetRHS(&b);
    prob.SetLHS(&x);
				     // And finally solve the problem.
    if (solver_params.OUTPUT == solver_params_type::VERBOSE)
      std::cout << "Starting solve\n" << std::flush;
    solver->Solve();
    niter = 0;
    lin_residual = 0;

				     // We must free the solver that was created
				     // for us.
    delete solver;

  } else if (solver_params.SOLVER == solver_params_type::GMRES) {

				     // For the iterative solvers, we use Aztec.
    AztecOO Solver;

				     // Select the appropriate level of verbosity.
    if (solver_params.OUTPUT == solver_params_type::QUIET)
      Solver.SetAztecOption(AZ_output, AZ_none);

    if (solver_params.OUTPUT == solver_params_type::VERBOSE)
      Solver.SetAztecOption(AZ_output, AZ_all);

				     // Select gmres.  Other solvers are available.
    Solver.SetAztecOption(AZ_solver, AZ_gmres);
    Solver.SetRHS(&b);
    Solver.SetLHS(&x);

				     // Set up the ILUT preconditioner.  I do not know
				     // why, but we must pretend like we are in parallel
				     // using domain decomposition or the preconditioner
				     // refuses to activate.
    Solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
    Solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
    Solver.SetAztecOption(AZ_overlap, 0);
    Solver.SetAztecOption(AZ_reorder, 0);

				     // ILUT parameters as described above.
    Solver.SetAztecParam(AZ_drop, solver_params.ILUT_DROP);
    Solver.SetAztecParam(AZ_ilut_fill, solver_params.ILUT_FILL);
    Solver.SetAztecParam(AZ_athresh, solver_params.ILUT_ATOL);
    Solver.SetAztecParam(AZ_rthresh, solver_params.ILUT_RTOL);
    Solver.SetUserMatrix(Matrix);

				     // Run the solver iteration.  Collect the number
				     // of iterations and the residual.
    Solver.Iterate(solver_params.MAX_ITERS, solver_params.RES);
    niter = Solver.NumIters();
    lin_residual = Solver.TrueResidual();
  }
}

                                 // @sect3{Postprocessing and Output}
                                 // Recover the physical variables from the conservative
                                 // variables so that output will be (perhaps) more
                                 // meaningfull.
template <int dim>
void ConsLaw<dim>::postprocess() {
  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
  std::vector<unsigned int> dofs (dofs_per_cell);
  UpdateFlags update_flags = update_values
			     | update_gradients
			     | update_q_points
			     | update_JxW_values;
  UpdateFlags update_flags1 = update_values
			      | update_gradients
			      | update_q_points
			      | update_JxW_values;

  QGauss<dim>  quadrature_formula(4);

  const std::vector<Point<dim> > &us = fe_ptr->base_element(0).get_unit_support_points();


  Quadrature<dim>  unit_support(us);

  int n_q_points = quadrature_formula.n_quadrature_points;
  int n_uq_points = unit_support.n_quadrature_points;

  FEValues<dim> fe_v (
    mapping, *fe_ptr, quadrature_formula, update_flags);

  FEValues<dim> fe_v_unit (
    mapping, *fe_ptr, unit_support, update_flags1);

  std::vector<Vector<double> > U(n_uq_points,
                                 Vector<double>(EulerEquations<dim>::n_components));
  std::vector<Vector<double> > UU(n_q_points,
				  Vector<double>(EulerEquations<dim>::n_components));
  std::vector<std::vector<Tensor<1,dim> > > dU(n_uq_points,
					       std::vector<Tensor<1,dim> >(EulerEquations<dim>::n_components));
  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

				   // Loop the cells
  for (unsigned int cell_no=0; cell!=endc; ++cell, ++cell_no) {
    cell->get_dof_indices (dofs);
    fe_v_unit.reinit(cell);
    fe_v.reinit(cell);

    fe_v_unit.get_function_values(solution, U);
    fe_v_unit.get_function_grads(solution, dU);
    fe_v.get_function_values(solution, UU);

    for (unsigned int q = 0; q < fe_v.get_fe().base_element(0).n_dofs_per_cell(); q++) {
      unsigned int didx = fe_v.get_fe().component_to_system_index(EulerEquations<dim>::density_component, q);
      unsigned int eidx = fe_v.get_fe().component_to_system_index(EulerEquations<dim>::energy_component, q);
      double rho_normVsqr = 0;
      for (unsigned int d = 0; d < dim; d++) {
        unsigned int vidx = fe_v.get_fe().component_to_system_index(d, q);
        ppsolution(dofs[vidx]) = solution(dofs[vidx])/solution(dofs[didx]);
        rho_normVsqr += solution(dofs[vidx])*solution(dofs[vidx]);
      }
      rho_normVsqr /= solution(dofs[didx]);
				       // Pressure
      ppsolution(dofs[eidx]) = (EulerEquations<dim>::gas_gamma-1.0)*(solution(dofs[eidx]) - 0.5*rho_normVsqr);

				       // Either output density or gradient squared of density,
				       // depending on what the user wants.
      if (!schlieren_plot) {
        ppsolution(dofs[didx]) = solution(dofs[didx]);
      } else {
        double ng = 0;
        for (unsigned int i = 0; i < dim; i++) ng += dU[q][EulerEquations<dim>::density_component][i]*dU[q][EulerEquations<dim>::density_component][i];
        ng = std::sqrt(ng);
        ppsolution(dofs[didx]) = ng;
      }
    }

  } // cell

}

				 // Loop and assign a value for refinement.  We
				 // simply use the density squared, which selects
				 // shocks with some success.
template <int dim>
void ConsLaw<dim>::estimate() {
  
  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
  std::vector<unsigned int> dofs (dofs_per_cell);
  UpdateFlags update_flags = update_values
			     | update_gradients
			     | update_q_points
			     | update_JxW_values;

  QGauss<dim>  quadrature_formula(1);
  unsigned int n_q_points = quadrature_formula.n_quadrature_points;


  FEValues<dim> fe_v (
    mapping, *fe_ptr, quadrature_formula, update_flags);

  std::vector<Vector<double> > U(n_q_points,
                                 Vector<double>(EulerEquations<dim>::n_components));
  std::vector<std::vector<Tensor<1,dim> > > dU(n_q_points,
					       std::vector<Tensor<1,dim> >(EulerEquations<dim>::n_components));
  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (unsigned int cell_no=0; cell!=endc; ++cell, ++cell_no) {
    fe_v.reinit(cell);

    fe_v.get_function_values(predictor, U);
    fe_v.get_function_grads(predictor, dU);

    indicator(cell_no) = 0;
    for (unsigned int q = 0; q < n_q_points; q++) {
      double ng = 0;
      for (unsigned int d = 0; d < dim; d++) ng += dU[q][EulerEquations<dim>::density_component][d]*dU[q][EulerEquations<dim>::density_component][d];

      indicator(cell_no) += std::log(1+std::sqrt(ng));
      
    } 
    indicator(cell_no) /= n_q_points;

  } 
}

template <int dim>
void ConsLaw<dim>::refine_grid ()
{

  SolutionTransfer<dim, double> soltrans(dof_handler);

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

				   // Loop cells.  If the indicator
				   // for the cell matches the refinement criterion,
				   // refine, else unrefine.  The unrefinement has
				   // a slight hysterisis to avoid 'flashing' from refined
				   // to unrefined.
  for (unsigned int cell_no=0; cell!=endc; ++cell, ++cell_no) {
    cell->clear_coarsen_flag();
    cell->clear_refine_flag();
    if (cell->level() < refinement_params.shock_levels &&
        std::fabs(indicator(cell_no)) > refinement_params.shock_val ) {
      cell->set_refine_flag();
    } else {
      if (cell->level() > 0 &&
	  std::fabs(indicator(cell_no)) < 0.75*refinement_params.shock_val)
	cell->set_coarsen_flag();
    }
  }

				   // The following code prolongs the solution
				   // to the new grid and carries out the refinement.
  std::vector<Vector<double> > interp_in;
  std::vector<Vector<double> > interp_out;

  interp_in.push_back(solution);
  interp_in.push_back(predictor);

  triangulation.prepare_coarsening_and_refinement();
  soltrans.prepare_for_coarsening_and_refinement(interp_in);

  triangulation.execute_coarsening_and_refinement ();

  dof_handler.clear();
  dof_handler.distribute_dofs (*fe_ptr);

  {
    Vector<double> new_solution(1);
    Vector<double> new_predictor(1);

    interp_out.push_back(new_solution);
    interp_out.push_back(new_predictor);
    interp_out[0].reinit(dof_handler.n_dofs());
    interp_out[1].reinit(dof_handler.n_dofs());
  }

  soltrans.interpolate(interp_in, interp_out);
  
				   // Let the vector delete a very small vector
  solution.reinit(1);
  predictor.reinit(1);
  solution.swap(interp_out[0]);
  predictor.swap(interp_out[1]);

				   // resize these vectors for the new grid.
  nlsolution.reinit(dof_handler.n_dofs());
  ppsolution.reinit(dof_handler.n_dofs());
  nlsolution = solution;
  dsolution.reinit (dof_handler.n_dofs());
  right_hand_side.reinit (dof_handler.n_dofs());

  indicator.reinit(triangulation.n_active_cells());

}

template <int dim>
void ConsLaw<dim>::output_results (const unsigned int cycle) const
{
  char filename[512];
  std::sprintf(filename, "solution-%03d.vtk", cycle);
  std::ofstream output (filename);

  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  std::vector<std::string> solution_names (dim, "velocity");
  solution_names.push_back ("density");
  solution_names.push_back ("pressure");

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation
    .push_back (DataComponentInterpretation::component_is_scalar);
  data_component_interpretation
    .push_back (DataComponentInterpretation::component_is_scalar);
  
  data_out.add_data_vector (ppsolution, solution_names,
			    DataOut<dim>::type_dof_data,
			    data_component_interpretation);

  data_out.add_data_vector (indicator, "error");
  data_out.build_patches ();
  data_out.write_vtk (output);

  output.close();
}

				 // @sect3{Parsing the Input Deck}
				 // Declare the parameters for the
				 // input deck.  We assume a certain
				 // maximum number of boundaries and process
				 // any boundary the user supplies up to
				 // that maximum number.  We
				 // leave a detailed explanation of these 
				 // parameters to our description of the input
				 // sample file.
const unsigned int MAX_BD = 10;
template <int dim>
void ConsLaw<dim>::declare_parameters() {

                                   // Global scope parameters/
  prm.declare_entry("mesh", "grid.inp",
                    Patterns::Anything(),
                    "intput file");

  prm.declare_entry("diffusion power", "2.0",
		    Patterns::Double(),
		    "power of mesh size for diffusion");

  prm.declare_entry("gravity", "0.0",
		    Patterns::Double(),
		    "gravity forcing");

                                   // Time stepping block
  prm.enter_subsection("time stepping");
  prm.declare_entry("time step", "0.1",
		    Patterns::Double(),
		    "simulation time step");
  prm.declare_entry("final time", "10.0",
		    Patterns::Double(),
		    "simulation end time");
  prm.leave_subsection();


				   // Declare the boundary parameters
  for (unsigned int b = 0; b < MAX_BD; b++) {
    char bd[512];
    std::sprintf(bd, "boundary_%d", b);
    prm.enter_subsection(bd);
    prm.declare_entry("no penetration", "false",
		      Patterns::Selection("true|false"),
		      "<true|false>");
				     // declare a slot for each of the conservative
				     // variables.
    for (unsigned int di = 0; di < EulerEquations<dim>::n_components; di++) {
      char var[512];
      std::sprintf(var, "w_%d", di);
      prm.declare_entry(var, "outflow",
			Patterns::Selection(
			  "inflow|outflow|pressure"),
			"<inflow|outflow|pressure>");
      
				       // for dirichlet, a function in x,y,z
      std::sprintf(var, "w_%d value", di);
      prm.declare_entry(var, "0.0",
			Patterns::Anything(),
			"expression in x,y,z");
    }

    prm.leave_subsection();
  }

				   // Initial condition block.
  prm.enter_subsection("initial condition");
  for (unsigned int di = 0; di < EulerEquations<dim>::n_components; di++) {
    char var[512];
    std::sprintf(var, "w_%d", di);
      
				     // for dirichlet, a function in x,y,z
    std::sprintf(var, "w_%d value", di);
    prm.declare_entry(var, "0.0",
		      Patterns::Anything(),
		      "expression in x,y,z");
  }
  prm.leave_subsection();

				   // The linear solver block.
  prm.enter_subsection("linear solver");
  prm.declare_entry("output", "quiet",
		    Patterns::Selection(
		      "quiet|verbose"),
		    "<quiet|verbose>");
  prm.declare_entry("method", "gmres",
		    Patterns::Selection(
		      "gmres|direct"),
		    "<gmres|direct>");
  prm.declare_entry("residual", "1e-10",
		    Patterns::Double(),
		    "linear solver residual");
  prm.declare_entry("max iters", "300",
		    Patterns::Double(),
		    "maximum solver iterations");
  prm.declare_entry("ilut fill", "2",
		    Patterns::Double(),
		    "ilut preconditioner fill");
  prm.declare_entry("ilut absolute tolerance", "1e-9",
		    Patterns::Double(),
		    "ilut preconditioner tolerance");
  prm.declare_entry("ilut relative tolerance", "1.1",
		    Patterns::Double(),
		    "rel tol");
  prm.declare_entry("ilut drop tolerance", "1e-10",
		    Patterns::Double(),
		    "ilut drop tol");
  prm.leave_subsection();


				   // A refinement controller block.
  prm.enter_subsection("refinement");
  prm.declare_entry("refinement", "none",
		    Patterns::Selection(
		      "none|fixed number|shock"),
		    "<on|off>");
  prm.declare_entry("refinement fraction", "0.1",
		    Patterns::Double(),
		    "Fraction of high refinement");
  prm.declare_entry("unrefinement fraction", "0.1",
		    Patterns::Double(),
		    "Fraction of low unrefinement");
  prm.declare_entry("max elements", "1000000",
		    Patterns::Double(),
		    "maximum number of elements");
  prm.declare_entry("shock value", "4.0",
		    Patterns::Double(),
		    "value for shock indicator");
  prm.declare_entry("shock levels", "3.0",
		    Patterns::Double(),
		    "number of shock refinement levels");
  prm.leave_subsection();

				   // Output control.
  prm.enter_subsection("output");
  prm.declare_entry("density", "standard",
		    Patterns::Selection(
		      "standard|schlieren"),
		    "<standard|schlieren>");
  prm.declare_entry("step", "-1",
		    Patterns::Double(),
		    "output once per this period");
  prm.leave_subsection();

				   // Flux control
  prm.enter_subsection("flux");
  prm.declare_entry("stab", "alpha",
		    Patterns::Selection(
		      "alpha|constant|mesh"),
		    "<alpha|constant|mesh>");
  prm.declare_entry("stab value", "1",
		    Patterns::Double(),
		    "alpha stabilization");
  prm.leave_subsection();


}

				 // Code to actually parse an input file.  This function
				 // matches the declarations above.
template <int dim>
void ConsLaw<dim>::load_parameters(const char *infile){

  prm.read_input(infile);

				   // The global parameters.
  mesh = prm.get("mesh");

  diffusion_power = prm.get_double("diffusion power");

  gravity = prm.get_double("gravity");

				   // The time stepping.
  prm.enter_subsection("time stepping");
  dT = prm.get_double("time step");
  std::cout << "dT=" << dT << std::endl;
  if (dT == 0) {
    is_stationary = true;
    dT = 1.0;
    TF = 1.0;
    std::cout << "Stationary mode" << std::endl;
  }
  TF = prm.get_double("final time");
  std::cout << "TF=" << TF << std::endl;
  prm.leave_subsection();

				   // The boundary info
  for (unsigned int b = 0; b < MAX_BD; b++) {
    std::vector<bc_type> flags(EulerEquations<dim>::n_components, OUTFLOW_BC);

				     // Define a parser for every boundary, though it may be
				     // unused.
    FunctionParser<dim> *sd = new FunctionParser<dim>(EulerEquations<dim>::n_components);

    std::vector<std::string> expressions(EulerEquations<dim>::n_components, "0.0");
    
    char bd[512];
    std::sprintf(bd, "boundary_%d", b);
    prm.enter_subsection(bd);

    const std::string nopen = prm.get("no penetration");

				     // Determine how each component is handled.
    for (unsigned int di = 0; di < EulerEquations<dim>::n_components; di++) {
      char var[512];
      std::sprintf(var, "w_%d", di);
      std::string btype = prm.get(var);
      std::sprintf(var, "w_%d value", di);
      std::string var_value = prm.get(var);

      if (di < dim && nopen == "true") {
        flags[di] = NO_PENETRATION_BC;
      } else if (btype == "inflow") {
        flags[di] = INFLOW_BC;
        expressions[di] = var_value;
      } else if (btype == "pressure") {
        flags[di] = PRESSURE_BC;
        expressions[di] = var_value;
      }
    } 
    prm.leave_subsection();

				     // Add the boundary condition to the law.
    sd->initialize (FunctionParser<dim>::default_variable_names(),
		    expressions,
		    std::map<std::string, double>());
    add_boundary(b, flags, sd);
  }

				   // Initial conditions.
  prm.enter_subsection("initial condition");
  {
    std::vector<std::string> expressions (EulerEquations<dim>::n_components,
					  "0.0");
    for (unsigned int di = 0; di < EulerEquations<dim>::n_components; di++)
      {
	char var[512];
	std::sprintf(var, "w_%d value", di);
	expressions[di] = prm.get(var);
      }
    initial_conditions.initialize (FunctionParser<dim>::default_variable_names(),
				   expressions,
				   std::map<std::string, double>());
  }
  prm.leave_subsection();

				   // The linear solver.
  prm.enter_subsection("linear solver");
  const std::string &op = prm.get("output");
  if (op == "verbose") solver_params.OUTPUT = solver_params_type::VERBOSE;
  if (op == "quiet") solver_params.OUTPUT = solver_params_type::QUIET;
  const std::string &sv = prm.get("method");
  if (sv == "direct") {
    solver_params.SOLVER = solver_params_type::DIRECT;
  } else if (sv == "gmres") {
    solver_params.SOLVER = solver_params_type::GMRES;
  } 

  solver_params.RES = prm.get_double("residual");
  solver_params.MAX_ITERS = (int) prm.get_double("max iters");
  solver_params.ILUT_FILL = prm.get_double("ilut fill");
  solver_params.ILUT_ATOL = prm.get_double("ilut absolute tolerance");
  solver_params.ILUT_RTOL = prm.get_double("ilut relative tolerance");
  solver_params.ILUT_DROP = prm.get_double("ilut drop tolerance");
  solver_params.RES = prm.get_double("residual");
  prm.leave_subsection();


				   // And refiement.
  prm.enter_subsection("refinement");
  const std::string &ref = prm.get("refinement");
  if (ref == "none") {
    refinement_params.refine = refinement_params_type::NONE;
  } else if (ref == "fixed number") {
    refinement_params.refine = refinement_params_type::FIXED_NUMBER;
  } else if (ref == "shock") {
    refinement_params.refine = refinement_params_type::SHOCK;
  } else
    refinement_params.high_frac = prm.get_double("refinement fraction");
  refinement_params.high_frac_sav = refinement_params.high_frac;
  refinement_params.low_frac = prm.get_double("unrefinement fraction");
  refinement_params.max_cells = prm.get_double("max elements");
  refinement_params.shock_val = prm.get_double("shock value");
  refinement_params.shock_levels = prm.get_double("shock levels");
  prm.leave_subsection();
    
				   // Output control.
  prm.enter_subsection("output");
  const std::string &dens = prm.get("density");
  schlieren_plot = dens == "schlieren" ? true : false;
  output_step = prm.get_double("step");
  prm.leave_subsection();

				   // Flux control.
  prm.enter_subsection("flux");
  const std::string &stab = prm.get("stab");
  if (stab == "constant") {
    flux_params.LF_stab = flux_params_type::CONSTANT;
  } else if (stab == "mesh ") {
    flux_params.LF_stab = flux_params_type::MESH;
  }
  flux_params.LF_stab_value = prm.get_double("stab value");
  prm.leave_subsection();


}

template<int dim>
void ConsLaw<dim>::zero_matrix() {
  Matrix->PutScalar(0); Matrix->FillComplete();
}

				 // We use a predictor to try and make adaptivity
				 // work better.  The idea is to try and refine ahead
				 // of a front, rather than stepping into a coarse
				 // set of elements and smearing the solution.  This
				 // simple time extrapolator does the job.
template<int dim>
void ConsLaw<dim>::compute_predictor() {
  predictor = nlsolution;
  predictor.sadd(3/2.0, -1/2.0, solution);
}

				 // @sect3{Run the simulation}
				 // Contains the initialization
				 // the time loop, and the inner Newton iteration.
template <int dim>
void ConsLaw<dim>::run () 
{

				   // Open and load the mesh.
  GridIn<dim> grid_in;
  grid_in.attach_triangulation(triangulation);
  std::cout << "Opening mesh <" << mesh << ">" << std::endl;
  std::ifstream input_file(mesh.c_str());

  Assert (input_file, ExcFileNotOpen(mesh.c_str()));

  grid_in.read_ucd(input_file);   
  input_file.close();
  
  build_fe();

  unsigned int nstep = 0;
  
				   // Initialize fields and matrices.
  initialize_system (); 
  setup_system();
  initialize(); 
  predictor = solution;

				   // Initial refinement.  We apply the ic,
				   // estimate, refine, and repeat until
				   // happy.
  if (refinement_params.refine != refinement_params_type::NONE)
    for (unsigned int i = 0; i < refinement_params.shock_levels; i++) {
      estimate();
      refine_grid();
      setup_system();
      initialize(); 
      predictor = solution;
    }
  postprocess();
  output_results (nstep);

				   // Determine when we will output next.
  double next_output = T + output_step;

				   // @sect4{Main time stepping loop}
  predictor = solution;
  while(T < TF)
    {
      std::cout << "T=" << T << ", ";


      std::cout << "   Number of active cells:       "
		<< triangulation.n_active_cells()
		<< std::endl;


      std::cout << "   Number of degrees of freedom: "
		<< dof_handler.n_dofs()
		<< std::endl;

      bool nonlin_done = false;
      double res_norm;
      int lin_iter;

				       // Print some relevant information during the
				       // Newton iteration.
      std::cout << "NonLin Res:       Lin Iter     Lin Res" << std::endl;
      std::cout << "______________________________________" << std::endl;

      const unsigned int max_nonlin = 7;
      unsigned int nonlin_iter = 0;
      double lin_res;

				       // @sect5{Newton iteration}
      nlsolution = predictor;
      while (!nonlin_done) {
        lin_iter = 0;
        zero_matrix();
        right_hand_side = 0;
        assemble_system (res_norm);
					 // Flash a star to the screen so one can
					 // know when the assembly has stopped and the linear
					 // solution is starting.
        std::cout << "* " << std::flush;

					 // Test against a (hardcoded) nonlinear tolderance.
					 // Do not solve the linear system at the last step 
					 // (since it would be a waste).
                      
        if (fabs(res_norm) < 1e-10) {
          nonlin_done = true;
        } else {
					   // Solve the linear system and update with the
					   // delta.
	  dsolution = 0;
	  solve (dsolution, lin_iter, lin_res);
	  nlsolution.add(1.0, dsolution);
        }

					 // Print the residuals.
        std::printf("%-16.3e %04d        %-5.2e\n",
		    res_norm, lin_iter, lin_res);

        ++nonlin_iter;

	AssertThrow (nonlin_iter <= max_nonlin,
		     ExcMessage ("No convergence in nonlinear solver"));
      } 

				       // Various post convergence tasks.
      compute_predictor();

      solution = nlsolution;


      estimate();

      postprocess();

      T += dT;

				       // Output if it is time.
      if (output_step < 0) {
        output_results (++nstep);
      } else if (T >= next_output) {
        output_results (++nstep);
        next_output += output_step;
      }

				       // Refine, if refinement is selected.
      if (refinement_params.refine != refinement_params_type::NONE) {
        refine_grid();
        setup_system();
      }
    }
}

				 // The following ``main'' function is
				 // similar to previous examples and
				 // need not to be commented on. Note
				 // that the program aborts if no
				 // input file name is given on the
				 // command line.
int main (int argc, char *argv[]) 
{
  if (argc != 2)
    {
      std::cout << "Usage:" << argv[0] << " infile" << std::endl;
      std::exit(1);
    }
  
  try
    {
      ConsLaw<2> cons;
      cons.declare_parameters();
      cons.load_parameters(argv[1]);
      cons.run ();
    }
  catch (std::exception &exc)
    {
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
  catch (...) 
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };
  
  return 0;
}



/* Authors: Jean-Paul Pelteret, University of Cape Town,            */
/*          Andrew McBride, University of Erlangen-Nuremberg, 2010  */
/*                                                                  */
/*    Copyright (C) 2010, 2011 by the deal.II authors                     */
/*                        & Jean-Paul Pelteret and Andrew McBride   */
/*                                                                  */
/*    This file is subject to QPL and may not be  distributed       */
/*    without copyright and license information. Please refer       */
/*    to the file deal.II/doc/license.html for the  text  and       */
/*    further information on this license.                          */

#include <base/function.h>
#include <base/parameter_handler.h>
#include <base/point.h>
#include <base/quadrature_lib.h>
#include <base/symmetric_tensor.h>
#include <base/tensor.h>
#include <base/timer.h>
#include <base/work_stream.h>

#include <dofs/dof_constraints.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_tools.h>

#include <grid/grid_generator.h>
#include <grid/grid_tools.h>
#include <grid/grid_in.h>
#include <grid/tria.h>
#include <grid/tria_boundary_lib.h>

#include <fe/fe_dgp_monomial.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <fe/fe_tools.h>
#include <fe/fe_values.h>

#include <fe/mapping_q_eulerian.h>

#include <lac/block_sparse_matrix.h>
#include <lac/block_vector.h>
#include <lac/compressed_sparsity_pattern.h>
#include <lac/full_matrix.h>
#include <lac/precondition.h>
#include <lac/solver_cg.h>
#include <lac/sparse_direct.h>

#include <numerics/data_out.h>
#include <numerics/vectors.h>

#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace dealii;

// @sect3{Run-time parameters}
namespace Parameters
{
// Finite Element system
struct FESystem
{
    int poly_degree;
    int quad_order;

    static void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
};

void FESystem::declare_parameters (ParameterHandler &prm)
{
    prm.enter_subsection("Finite element system");
    {
	prm.declare_entry("Polynomial degree",
			  "1",
			  Patterns::Integer(),
			  "Displacement system polynomial order");

	prm.declare_entry("Quadrature order",
			  "2",
			  Patterns::Integer(),
			  "Gauss quadrature order");
    }
    prm.leave_subsection();
}

void FESystem::parse_parameters (ParameterHandler &prm)
{
    prm.enter_subsection("Finite element system");
    {
	poly_degree  = prm.get_integer("Polynomial degree");
	quad_order  = prm.get_integer("Quadrature order");
    }
    prm.leave_subsection();
}

// Geometry
struct Geometry
{
    int global_refinement;
    double scale;
    double p_p0;

    static void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
};

void Geometry::declare_parameters (ParameterHandler &prm)
{
    prm.enter_subsection("Geometry");
    {
	prm.declare_entry("Global refinement",
			  "2",
			  Patterns::Integer(),
			  "Global refinement level");

	prm.declare_entry("Grid scale",
			  "1.0",
			  Patterns::Double(),
			  "Global grid scaling factor");

	prm.declare_entry("Pressure ratio p/p0",
			  "40",
			  Patterns::Selection("20|40|60|80|100"),
			  "Ratio of applied pressure to reference pressure");
    }
    prm.leave_subsection();
}

void Geometry::parse_parameters (ParameterHandler &prm)
{
    prm.enter_subsection("Geometry");
    {
	global_refinement = prm.get_integer("Global refinement");
	scale  = prm.get_double("Grid scale");
	p_p0= prm.get_double("Pressure ratio p/p0");
    }
    prm.leave_subsection();
}

// Materials
struct Materials
{
    double nu;
    double mu;

    static void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
};

void Materials::declare_parameters (ParameterHandler &prm)
{
    prm.enter_subsection("Material properties");
    {
	prm.declare_entry("Poisson's ratio",
			  "0.49",
			  Patterns::Double(),
			  "Poisson's ratio");

	prm.declare_entry("Shear modulus",
			  "1.0e6",
			  Patterns::Double(),
			  "Shear modulus");
    }
    prm.leave_subsection();
}

void Materials::parse_parameters (ParameterHandler &prm)
{
    prm.enter_subsection("Material properties");
    {
	nu  = prm.get_double("Poisson's ratio");
	mu  = prm.get_double("Shear modulus");
    }
    prm.leave_subsection();
}

// Linear solver
struct LinearSolver
{
    std::string type_lin;
    double tol_lin;
    double max_iterations_lin;
    double ssor_relaxation;

    static void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
};

void LinearSolver::declare_parameters (ParameterHandler &prm)
{
    prm.enter_subsection("Linear solver");
    {
	prm.declare_entry("Solver type",
			  "CG",
			  Patterns::Selection("CG|Direct"),
			  "Type of solver used to solve the linear system");

	prm.declare_entry("Residual",
			  "1e-6",
			  Patterns::Double(),
			  "Linear solver residual (scaled by residual norm)");

	prm.declare_entry("Max iteration multiplier",
			  "2",
			  Patterns::Double(),
			  "Linear solver iterations (multiples of the system matrix size)");

	prm.declare_entry("SSOR Relaxation",
			  "0.6",
			  Patterns::Double(),
			  "SSOR preconditioner relaxation value");
    }
    prm.leave_subsection();
}

void LinearSolver::parse_parameters (ParameterHandler &prm)
{
    prm.enter_subsection("Linear solver");
    {
	type_lin = prm.get("Solver type");
	tol_lin = prm.get_double("Residual");
	max_iterations_lin  = prm.get_double("Max iteration multiplier");
	ssor_relaxation = prm.get_double("SSOR Relaxation");
    }
    prm.leave_subsection();
}

// Nonlinear solver
struct NonlinearSolver
{
    unsigned int max_iterations_NR;
    double tol_f;
    double tol_u;

    static void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
};

void NonlinearSolver::declare_parameters (ParameterHandler &prm)
{
    prm.enter_subsection("Nonlinear solver");
    {
	prm.declare_entry("Max iterations Newton-Raphson",
			  "10",
			  Patterns::Integer(),
			  "Number of Newton-Raphson iterations allowed");

	prm.declare_entry("Tolerance force",
			  "1.0e-9",
			  Patterns::Double(),
			  "Force residual tolerance");

	prm.declare_entry("Tolerance displacement",
			  "1.0e-3",
			  Patterns::Double(),
			  "Displacement error tolerance");
    }
    prm.leave_subsection();
}

void NonlinearSolver::parse_parameters (ParameterHandler &prm)
{
    prm.enter_subsection("Nonlinear solver");
    {
	max_iterations_NR  = prm.get_integer("Max iterations Newton-Raphson");
	tol_f = prm.get_double("Tolerance force");
	tol_u = prm.get_double("Tolerance displacement");
    }
    prm.leave_subsection();
}

// Time
struct Time
{
    double end_time;
    double delta_t;

    static void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
};

void Time::declare_parameters (ParameterHandler &prm)
{
    prm.enter_subsection("Time");
    {
	prm.declare_entry("End time",
			  "1",
			  Patterns::Double(),
			  "End time");

	prm.declare_entry("Time step size",
			  "0.1",
			  Patterns::Double(),
			  "Time step size");
    }
    prm.leave_subsection();
}

void Time::parse_parameters (ParameterHandler &prm)
{
    prm.enter_subsection("Time");
    {
	end_time  = prm.get_double("End time");
	delta_t  = prm.get_double("Time step size");
    }
    prm.leave_subsection();
}

// All parameters
struct AllParameters
	:
	public FESystem,
	public Geometry,
	public Materials,
	public LinearSolver,
	public NonlinearSolver,
	public Time

{
    AllParameters (const std::string & input_file);

    static void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
};

AllParameters::AllParameters (const std::string & input_file)
{
    ParameterHandler prm;
    declare_parameters(prm);
    prm.read_input (input_file);
    parse_parameters(prm);
}

void AllParameters::declare_parameters (ParameterHandler &prm)
{
    FESystem::declare_parameters(prm);
    Geometry::declare_parameters(prm);
    Materials::declare_parameters(prm);
    LinearSolver::declare_parameters(prm);
    NonlinearSolver::declare_parameters(prm);
    Time::declare_parameters(prm);
}

void AllParameters::parse_parameters (ParameterHandler &prm)
{
    FESystem::parse_parameters(prm);
    Geometry::parse_parameters(prm);
    Materials::parse_parameters(prm);
    LinearSolver::parse_parameters(prm);
    NonlinearSolver::parse_parameters(prm);
    Time::parse_parameters(prm);
}
}

// @sect3{General tools}
namespace AdditionalTools
{
template <typename MatrixType>
void extract_submatrix(const std::vector< unsigned int > &row_index_set,
		       const std::vector< unsigned int > &column_index_set,
		       const MatrixType &matrix,
		       FullMatrix< double > &sub_matrix )
{

    const unsigned int n_rows_submatrix = row_index_set.size();
    const unsigned int n_cols_submatrix = column_index_set.size();

    sub_matrix.reinit(n_rows_submatrix, n_cols_submatrix);

    for (unsigned int sub_row = 0; sub_row < n_rows_submatrix; ++sub_row) {
	const unsigned int row = row_index_set[sub_row];
	Assert (row<=matrix.m(), ExcIndexRange(row, 0, matrix.m()));

	for (unsigned int sub_col = 0; sub_col < n_cols_submatrix; ++sub_col) {
	    const unsigned int col = column_index_set[sub_col];
	    Assert (col<=matrix.n(), ExcIndexRange(col, 0, matrix.n()));

	    sub_matrix(sub_row,sub_col) = matrix(row, col);
	}
    }
}

template <typename MatrixType>
void replace_submatrix(const std::vector< unsigned int > &row_index_set,
		       const std::vector< unsigned int > &column_index_set,
		       const MatrixType &sub_matrix,
		       FullMatrix< double >  &matrix)
{
    const unsigned int n_rows_submatrix = row_index_set.size();
    Assert (n_rows_submatrix<=sub_matrix.m(), ExcIndexRange(n_rows_submatrix, 0, sub_matrix.m()));
    const unsigned int n_cols_submatrix = column_index_set.size();
    Assert (n_cols_submatrix<=sub_matrix.n(), ExcIndexRange(n_cols_submatrix, 0, sub_matrix.n()));

    for (unsigned int sub_row = 0; sub_row < n_rows_submatrix; ++sub_row) {
	const unsigned int row = row_index_set[sub_row];
	Assert (row<=matrix.m(), ExcIndexRange(row, 0, matrix.m()));

	for (unsigned int sub_col = 0; sub_col < n_cols_submatrix; ++sub_col) {
	    const unsigned int col = column_index_set[sub_col];
	    Assert (col<=matrix.n(), ExcIndexRange(col, 0, matrix.n()));

	    matrix(row, col) = sub_matrix(sub_row, sub_col);

	}
    }
}

}

// @sect3{Time class}
class Time {
public:
    Time (const double & time_end,
          const double & delta_t)
	:
	  timestep (0),
	  time_current (0.0),
	  time_end (time_end),
	  delta_t (delta_t)
    {}
    virtual ~Time (void) {}

    const double & current (void) const {return time_current;}
    const double & end (void) const {return time_end;}
    const double & get_delta_t (void) const {return delta_t;}
    const unsigned int & get_timestep (void) const {return timestep;}
    void increment (void) {time_current += delta_t; ++timestep;}

private:
    unsigned int timestep;
    double time_current;
    const double time_end;
    const double delta_t;
};

// @sect3{Neo-Hookean material}
template <int dim>
class Material_NH
{
public:
    /// \brief Class constructor
    Material_NH (const double & lambda,
		 const double & mu)
	:
	  lambda_0 (lambda),
	  mu_0 (mu),
	  kappa_0 (lambda + 2.0/3.0*mu)
    { }
    virtual ~Material_NH (void) {};

    // Stress and constitutive tensors
    virtual SymmetricTensor<2, dim> get_T (const double & J,
					   const SymmetricTensor <2, dim> & B)
    {
	const double dW_dJ  = get_dU_dtheta (J);
	return mu_0*B + dW_dJ*J*I;
    }

    virtual SymmetricTensor<4, dim> get_JC (const double & J,
					    const SymmetricTensor <2, dim> & B)
    {
	const double dW_dJ   = get_dU_dtheta (J);
	const double d2W_dJ2 = get_d2U_dtheta2 (J);
	return  J*(  (dW_dJ + J*d2W_dJ2)*IxI - (2.0*dW_dJ)*II  );
    }

    // Volumetric quantities methods
    double get_dU_dtheta    (const double & d) {return kappa_0*(d - 1.0/d);}
    double get_d2U_dtheta2  (const double & d) {return kappa_0*(1.0 + 1.0/(d*d));}

protected:
    // Material properties
    const double lambda_0; // Lame modulus
    const double mu_0;     // Shear modulus
    const double kappa_0;  // Bulk modulus

    static SymmetricTensor<2, dim> const I;
    static SymmetricTensor<4, dim> const IxI;
    static SymmetricTensor<4, dim> const II;
};

template <int dim> SymmetricTensor<2, dim> const Material_NH<dim>::I   = SymmetricTensor<2, dim> (unit_symmetric_tensor <dim> ());
template <int dim> SymmetricTensor<4, dim> const Material_NH<dim>::IxI = SymmetricTensor<4, dim> (outer_product (I, I));
template <int dim> SymmetricTensor<4, dim> const Material_NH<dim>::II  = SymmetricTensor<4, dim> (identity_tensor <dim> ());

// @sect3{Quadrature point history}
template <int dim>
class PointHistory
{
public:
    PointHistory (void)
	:
	  material (NULL),
	  dilatation_n (1.0),
	  pressure_n (0.0)
    { }
    virtual ~PointHistory (void) {delete material;}

    void setup_lqp ( Parameters::AllParameters & parameters )
    {
	const double lambda = 2.0*parameters.mu*parameters.nu / (1.0-2.0*parameters.nu);
	material = new Material_NH<dim> (lambda,
					 parameters.mu);

        // Initialise all tensors correctly
	update_values (Tensor <2,dim> (), 0.0, 1.0);
    }

    // Total Variables
    void update_values (const Tensor<2, dim> & grad_u_n,
			const double & pressure,
			const double & dilatation)
    {
	// Calculated variables from displacement, displacement gradients
	const Tensor <2,dim>  F = static_cast <Tensor < 2, dim> > (unit_symmetric_tensor <dim> ()) + grad_u_n;
	J     = determinant(F);
	F_inv = invert(F);
	B_bar = std::pow(get_J(), -2.0/3.0) * symmetrize ( F* transpose (F) );

        // Precalculated pressure, dilatation
	pressure_n = pressure;
	dilatation_n = dilatation;

        // Now that all the necessary variables are set, we can update the stress tensors
        // Stress update can only update the stresses once the
        // dilatation has been set as p = p(d)
        T_bar = material->get_T (get_J(), get_B_bar());
        T_iso = dev_P*get_T_bar(); // Note: T_iso depends on T_bar
        T_vol = get_pressure()*get_J()*I;
    }

    // Displacement and strain
    const double & get_dilatation(void) const {return dilatation_n;}
    const double & get_J (void) const {return J;}
    const Tensor <2,dim> & get_F_inv (void) const {return F_inv;}
    const SymmetricTensor <2,dim> & get_B_bar (void) const {return B_bar;}

    // Volumetric terms
    double get_dU_dtheta (void) {
	return material->get_dU_dtheta(get_dilatation());
    }

    double get_d2U_dtheta2 (void) {
	return material->get_d2U_dtheta2(get_dilatation());
    }

    // Stress
    double get_pressure(void) {return pressure_n;}
    const SymmetricTensor<2, dim> & get_T_iso (void) const {return T_iso;}
    const SymmetricTensor<2, dim> & get_T_vol (void) const {return T_vol;};

    // Tangent matrices
    SymmetricTensor <4,dim> get_C_iso(void)
    {
        const double & J = get_J();
        const SymmetricTensor<2, dim> & B_bar = get_B_bar();
        const SymmetricTensor<2, dim> & T_iso = get_T_iso();

        const SymmetricTensor <4,dim> T_iso_x_I = outer_product(T_iso, I);
        const SymmetricTensor <4,dim> I_x_T_iso = outer_product(I, T_iso);
	const SymmetricTensor <4,dim> CC_bar = material->get_JC (J, B_bar);

	return     2.0/3.0*trace(get_T_bar())*dev_P
		-  2.0/3.0*(T_iso_x_I + I_x_T_iso)
		+  dev_P*CC_bar*dev_P;
    }

    SymmetricTensor <4,dim> get_C_vol(void)
    {
	const double & p = get_pressure();
	const double & J = get_J();
	return p*J*(IxI - 2.0*II);
    }

private:
    // === MATERIAL ===
    Material_NH <dim>* material;

    // ==== VOLUME, DISPLACEMENT AND STRAIN VARIABLES ====
    double                  dilatation_n;   // Current dilatation
    double                  J;
    Tensor <2,dim>	    F_inv;
    SymmetricTensor <2,dim> B_bar;
    SymmetricTensor <2,dim> E;

    // ==== STRESS VARIABLES ====
    double                  pressure_n; // Current pressure
    SymmetricTensor<2, dim> T_bar;
    SymmetricTensor<2, dim> T_iso;
    SymmetricTensor<2, dim> T_vol;
    const SymmetricTensor<2, dim> & get_T_bar (void) const {return T_bar;}

    // Basis tensors
    static SymmetricTensor<2, dim> const I;
    static SymmetricTensor<4, dim> const IxI;
    static SymmetricTensor<4, dim> const II;
    static SymmetricTensor<4, dim> const dev_P;
};

template <int dim> SymmetricTensor<2,dim> const PointHistory<dim>::I
= SymmetricTensor<2,dim> (unit_symmetric_tensor <dim> ());
template <int dim> SymmetricTensor<4,dim> const PointHistory<dim>::IxI
= SymmetricTensor<4,dim> (outer_product (I, I));
template <int dim> SymmetricTensor<4,dim> const PointHistory<dim>::II
= SymmetricTensor<4,dim> (identity_tensor <dim> ());
template <int dim> SymmetricTensor<4,dim> const PointHistory<dim>::dev_P
= SymmetricTensor<4,dim> (II - 1.0/3.0*IxI);


// @sect3{Quasi-static quasi-incompressible finite-strain solid}
template <int dim>
class Solid
{
public:
    Solid (const std::string & input_file);
    virtual ~Solid (void);
    void run (void);

private:

    // === DATA STRUCTS ===

    struct PerTaskData_K
    {
	FullMatrix<double>          cell_matrix;
	std::vector<unsigned int>   local_dof_indices;

	PerTaskData_K (const unsigned int dofs_per_cell)
	    :
	      cell_matrix        (dofs_per_cell,
		  dofs_per_cell),
	      local_dof_indices  (dofs_per_cell)
	{ }

        void reset (void) {
            cell_matrix = 0.0;
        }
    };

    struct ScratchData_K
    {
	FEValues <dim> fe_values_ref;

	std::vector < std::vector< double > >                  Nx;
	std::vector < std::vector< Tensor<2, dim> > >          grad_Nx;
	std::vector < std::vector< SymmetricTensor<2, dim> > > symm_grad_Nx;

	ScratchData_K ( const FiniteElement <dim> & fe_cell,
		       const QGauss <dim> & qf_cell,
		       const UpdateFlags uf_cell)
	    :
	      fe_values_ref   (fe_cell,
		  qf_cell,
		  uf_cell),
	      Nx              (qf_cell.size(),
		  std::vector< double >(fe_cell.dofs_per_cell)),
	      grad_Nx         (qf_cell.size(),
		  std::vector< Tensor<2, dim> >(fe_cell.dofs_per_cell)),
	      symm_grad_Nx    (qf_cell.size(),
		  std::vector< SymmetricTensor<2, dim> >(fe_cell.dofs_per_cell))
	{  }

	ScratchData_K ( const ScratchData_K & rhs ) :
	    fe_values_ref ( rhs.fe_values_ref.get_fe(),
		rhs.fe_values_ref.get_quadrature(),
		rhs.fe_values_ref.get_update_flags() ),
	    Nx (rhs.Nx),
	    grad_Nx (rhs.grad_Nx),
	    symm_grad_Nx (rhs.symm_grad_Nx)
	{  }

	void reset (void) {
            for (unsigned int q_point=0; q_point < grad_Nx.size(); ++q_point) {
                for (unsigned int k=0; k < Nx.size(); ++k) {
		    Nx[q_point][k] = 0.0;
		    grad_Nx[q_point][k] = 0.0;
		    symm_grad_Nx[q_point][k] = 0.0;
		}
	    }
	}

    };

    struct PerTaskData_F
    {
	Vector<double>              cell_rhs;
	std::vector<unsigned int>   local_dof_indices;

	PerTaskData_F (const unsigned int dofs_per_cell)
	    :
	      cell_rhs           (dofs_per_cell),
	      local_dof_indices  (dofs_per_cell)
	{ }

	void reset (void) { cell_rhs = 0.0; }
    };

    struct ScratchData_F
    {
	FEValues <dim>     fe_values_ref;
	FEFaceValues <dim> fe_face_values_ref;

	std::vector < std::vector< double > > Nx;
	std::vector < std::vector< SymmetricTensor<2, dim> > > symm_grad_Nx;
	std::vector< Vector<double> > rhs_values;

	// Solution data
	std::vector< std::vector<Tensor <1,dim> > > solution_grads;

	ScratchData_F ( const FiniteElement <dim> & fe_cell,
		       const QGauss <dim> & qf_cell,
		       const UpdateFlags uf_cell,
		       const QGauss <dim-1> & qf_face,
		       const UpdateFlags uf_face)
	    :
	      fe_values_ref   (fe_cell,
		  qf_cell,
		  uf_cell),
	      fe_face_values_ref   (fe_cell,
		  qf_face,
		  uf_face),
	      Nx              (qf_cell.size(),
		  std::vector< double >(fe_cell.dofs_per_cell)),
	      symm_grad_Nx    (qf_cell.size(),
		  std::vector< SymmetricTensor<2, dim> >(fe_cell.dofs_per_cell)),
	      rhs_values   (qf_cell.size(),
		  Vector<double>(dim))
	{  }

	ScratchData_F ( const ScratchData_F & rhs )
	    :
	      fe_values_ref ( rhs.fe_values_ref.get_fe(),
		  rhs.fe_values_ref.get_quadrature(),
		  rhs.fe_values_ref.get_update_flags() ),
	      fe_face_values_ref ( rhs.fe_face_values_ref.get_fe(),
		  rhs.fe_face_values_ref.get_quadrature(),
		  rhs.fe_face_values_ref.get_update_flags() ),
	      Nx (rhs.Nx),
	      symm_grad_Nx (rhs.symm_grad_Nx),
	      rhs_values (rhs.rhs_values)
	{  }

	void reset (void) {
	    for (unsigned int q_point=0; q_point < symm_grad_Nx.size(); ++q_point) {
		for (unsigned int k=0; k < symm_grad_Nx[q_point].size(); ++k) {
		    Nx[q_point][k] = 0.0;
		    symm_grad_Nx[q_point][k] = 0.0;
		    rhs_values[q_point] = 0.0;
		}
	    }
	}

    };

    struct PerTaskData_SC
    {
        FullMatrix<double>          cell_matrix;
        std::vector<unsigned int>   local_dof_indices;

	// Calculation matrices (auto resized)
	FullMatrix<double> K_orig;
	FullMatrix<double> K_pu;
	FullMatrix<double> K_pt;
	FullMatrix<double> K_tt;
	// Calculation matrices (manual resized)
	FullMatrix<double> K_pt_inv;
	FullMatrix<double> K_tt_inv;
	FullMatrix<double> K_con;
	FullMatrix<double> A;
	FullMatrix<double> B;
	FullMatrix<double> C;

	PerTaskData_SC (const unsigned int & dofs_per_cell,
			const unsigned int & n_u,
			const unsigned int & n_p,
			const unsigned int & n_t)
            :
	      cell_matrix        (dofs_per_cell,
		  dofs_per_cell),
	      local_dof_indices  (dofs_per_cell),
	      K_pt_inv (n_t, n_p),
	      K_tt_inv (n_t, n_t),
	      K_con (n_u, n_u),
	      A (n_t, n_u),
	      B (n_t, n_u),
	      C (n_p, n_u)
	{  }

	// Choose not to reset any data
	// The matrix extraction and replacement tools will take care of this
	void reset(void) { }
    };

    // Dummy struct for TBB
    struct ScratchData_SC
    {
	ScratchData_SC (void) { }
	ScratchData_SC (const ScratchData_SC & rhs) { }
	void reset (void) { }
    };

    // Dummy struct for TBB
    struct PerTaskData_UQPH
    {
	PerTaskData_UQPH (void) { }
	void reset(void) { }
    };

    struct ScratchData_UQPH
    {
	FEValues<dim> fe_values_ref;
	std::vector< Tensor< 2, dim> > solution_grads_u_total;
	std::vector <double> solution_values_p_total;
	std::vector <double> solution_values_t_total;
	const BlockVector <double> & solution_total;

	ScratchData_UQPH (const FiniteElement <dim> & fe_cell,
			  const QGauss <dim> & qf_cell,
			  const UpdateFlags uf_cell,
			  const BlockVector <double> & solution_total)
	    :
	      fe_values_ref (fe_cell,
		  qf_cell,
		  uf_cell),
	      solution_grads_u_total (qf_cell.size()),
	      solution_values_p_total (qf_cell.size()),
	      solution_values_t_total (qf_cell.size()),
	      solution_total (solution_total)
	{ }

	ScratchData_UQPH (const ScratchData_UQPH & rhs)
	    :
	      fe_values_ref (rhs.fe_values_ref.get_fe(),
		  rhs.fe_values_ref.get_quadrature(),
		  rhs.fe_values_ref.get_update_flags()),
	      solution_grads_u_total (rhs.solution_grads_u_total),
	      solution_values_p_total (rhs.solution_values_p_total),
	      solution_values_t_total (rhs.solution_values_t_total),
	      solution_total (rhs.solution_total)
	{ }

	void reset (void)
	{
	    // Is this necessary? Won't the call to fe_values.get_gradient overwrite this data?
	    for (unsigned int q=0; q < qf_cell.size(); ++q)
	    {
		solution_grads_u_total[q] = 0.0;
		solution_values_p_total[q] = 0.0;
		solution_values_t_total[q] = 0.0;
	    }
	}
    };

    // === METHODS ===

    /// \brief Print out a greeting for the user
    void make_grid (void);
    /// \brief Setup the Finite Element system to be solved
    void system_setup (void);
    void determine_component_extractors(void);

    /// \brief Assemble the system and right hand side matrices using multi-threading
    void assemble_system_K          (void);
    void assemble_system_K_one_cell (const typename DoFHandler<dim>::active_cell_iterator & cell,
				     ScratchData_K & scratch,
				     PerTaskData_K & data);
    void copy_local_to_global_K     (const PerTaskData_K & data);
    void assemble_system_F          (void);
    void assemble_system_F_one_cell (const typename DoFHandler<dim>::active_cell_iterator & cell,
				     ScratchData_F & scratch,
				     PerTaskData_F & data);
    void copy_local_to_global_F     (const PerTaskData_F & data);
    void assemble_SC                (void);
    void assemble_SC_one_cell       (const typename DoFHandler<dim>::active_cell_iterator & cell,
                                     ScratchData_SC & scratch,
                                     PerTaskData_SC & data);
    void copy_local_to_global_SC    (const PerTaskData_SC & data);
    /// \brief Apply Dirichlet boundary values
    void make_constraints (const int & it_nr,
			   ConstraintMatrix & constraints);

    //    /// \brief Setup the quadrature point history for each cell
    void setup_qph(void);
    //    /// \brief Update the quadrature points stress and strain values, and fibre directions
    void update_qph_incremental ( const BlockVector <double> & solution_delta );
    void update_qph_incremental_one_cell (const typename DoFHandler<dim>::active_cell_iterator & cell,
					  ScratchData_UQPH & scratch,
					  PerTaskData_UQPH & data);
    void copy_local_to_global_UQPH    (const PerTaskData_UQPH & data) {}
    /// \brief Solve for the displacement using a Newton-Rhapson method
    void solve_nonlinear_timestep (BlockVector <double> & solution_delta);
    void solve_linear_system (BlockVector <double> & newton_update);

    /// \brief Error measurement
    void get_error_res (const BlockVector <double> & residual, BlockVector <double> & error_res);
    void get_error_update (const BlockVector <double> & newton_update, BlockVector <double> & error_update);
    double get_error_dil (void);

    // Solution
    BlockVector <double> get_solution_total (const BlockVector <double> & solution_delta);

    // Postprocessing
    void output_results(void);

    // === ATTRIBUTES ===
    // Parameters
    Parameters::AllParameters parameters;

    // Geometry
    Triangulation<dim> triangulation; // Describes the triangulation

    // Time
    Time time;
    TimerOutput timer;

    // === Quadrature points ===
    std::vector< PointHistory <dim> > quadrature_point_history; // Quadrature point history

    // === Finite element system ===
    DoFHandler<dim>     dof_handler_ref; // Describes the degrees of freedom
    const unsigned int  degree;
    const FESystem<dim> fe; // Describes the global FE system

    unsigned int dofs_per_cell; // Number of degrees of freedom on each cell
    const FEValuesExtractors::Vector u_fe;
    const FEValuesExtractors::Scalar p_fe;
    const FEValuesExtractors::Scalar t_fe;

    // Block description
    static const unsigned int n_blocks  = 3;
    static const unsigned int n_components = dim + 2;
    static const unsigned int first_u_component = 0;
    static const unsigned int p_component = dim;
    static const unsigned int t_component = dim + 1;

    enum {u_dof=0 , p_dof, t_dof};
    std::vector<unsigned int> dofs_per_block;
    std::vector<unsigned int> element_indices_u;
    std::vector<unsigned int> element_indices_p;
    std::vector<unsigned int> element_indices_t;

    // === Quadrature ===
    QGauss<dim> qf_cell; // Cell quadrature formula
    QGauss<dim-1> qf_face; // Face quadrature formula
    unsigned int n_q_points; // Number of quadrature points in a cell
    unsigned int n_q_points_f; // Number of quadrature points in a face

    // === Stiffness matrix setup ====
    ConstraintMatrix constraints; // Matrix to keep track of all constraints
    BlockSparsityPattern sparsity_pattern; // Sparsity pattern for the stiffness matrix
    BlockSparseMatrix <double> tangent_matrix; // Global stiffness matrix
    BlockVector <double> residual; // Holds the residual vector
    BlockVector <double> solution_n; // Holds the solution vector: Total displacement over all time-steps
};

// @sect3{Implementation of the <code>Solid</code> class}

// @sect4{Public interface}
template <int dim>
Solid<dim>::Solid (const std::string & input_file)
    :
      parameters (input_file),
      triangulation (Triangulation<dim>::maximum_smoothing),
      time (parameters.end_time, parameters.delta_t),
      timer (std::cout,
	  TimerOutput::summary,
	  TimerOutput::wall_times),
      dof_handler_ref (triangulation),
      degree (parameters.poly_degree),
      fe (FE_Q<dim>(parameters.poly_degree), dim,    // displacement
	  FE_DGPMonomial<dim>(parameters.poly_degree-1), 1,  // pressure
	  FE_DGPMonomial<dim>(parameters.poly_degree-1), 1), // dilatation
      u_fe (first_u_component),
      p_fe (p_component),
      t_fe (t_component),
      dofs_per_block (n_blocks),
      qf_cell (parameters.quad_order),
      qf_face (parameters.quad_order)
{
    n_q_points = qf_cell.size();
    n_q_points_f = qf_face.size();
    dofs_per_cell = fe.dofs_per_cell;
    determine_component_extractors();
}

template <int dim>
Solid<dim>::~Solid (void)
{
    dof_handler_ref.clear ();
}

template <int dim>
void Solid<dim>::run (void)
{
    // Pre-processing
    make_grid ();
    system_setup ();
    output_results (); // Output initial grid position
    time.increment();

    BlockVector <double> solution_delta (dofs_per_block);
    solution_delta.collect_sizes ();

    while (time.current() <= time.end()) {
	solution_delta = 0.0;

	// Solve step and update total solution vector
	solve_nonlinear_timestep (solution_delta);
	solution_n += solution_delta;

	output_results ();
	time.increment();
    }
}

// @sect4{Solid::make_grid}
template <int dim>
void Solid<dim>::make_grid (void)
{
    GridGenerator::hyper_rectangle ( triangulation,
				    Point<dim> (0.0, 0.0, 0.0),
				    Point<dim> (1.0, 1.0, 1.0),
				    true );
    GridTools::scale (parameters.scale, triangulation);

    // Need to refine at least once for the indentation problem
    if (parameters.global_refinement == 0) triangulation.refine_global (1);
    else triangulation.refine_global (parameters.global_refinement);

    // Apply different BC's to a patch on the top surface
    typename Triangulation<dim>::active_cell_iterator
	    cell = triangulation.begin_active(),
	    endc = triangulation.end();
    for (; cell!=endc; ++cell)
    {
        if (cell->at_boundary() == true) {
	    for (unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; ++face) {
		// Find faces on the +y surface
		if (   cell->face(face)->at_boundary() == true
			&& cell->face(face)->center()[2] == 1.0*parameters.scale)
		{
		    if (   cell->face(face)->center()[0] < 0.5*parameters.scale
			    && cell->face(face)->center()[1] < 0.5*parameters.scale)
		    {
			cell->face(face)->set_boundary_indicator (6); // Set a new boundary id on a patch
		    }
		}
	    }
	}
    }
}

// @sect4{Solid::system_setup}
template <int dim>
void Solid<dim>::system_setup (void)
{
    timer.enter_subsection ("Setup system");

    // Number of components per block
    std::vector<unsigned int> block_component (n_components, u_dof); // Displacement
    block_component[p_component] = p_dof; // Pressure
    block_component[t_component] = t_dof; // Dilatation

    // Setup DOF handler
    dof_handler_ref.distribute_dofs (fe);
    DoFRenumbering::Cuthill_McKee (dof_handler_ref);
    DoFRenumbering::component_wise (dof_handler_ref, block_component);
    // Count number of dofs per block
    DoFTools::count_dofs_per_block (dof_handler_ref, dofs_per_block, block_component);

    std::cout
	    << "Triangulation:"
	    << "\n\t Number of active cells: " << triangulation.n_active_cells()
	    << "\n\t Number of degrees of freedom: " << dof_handler_ref.n_dofs()
	    << std::endl;

    // the global system matrix will have the following structure
    //      | K'_uu |   K_up    |     0     |         | dU_u |         | dR_u |
    // K =  | K_pu  |   K_tt^-1 |   K_pt^-1 | , dU =  | dU_p | , dR =  | dR_p |
    //      |   0   |   K_tp    |   K_tt    |         | dU_t |         | dR_t |
    // reflect this structure in the sparsity pattern
    Table<2,DoFTools::Coupling> coupling (n_components, n_components);
    for (unsigned int ii = 0; ii < n_components; ++ii) {
        for (unsigned int jj = ii; jj < n_components; ++jj) {
            if ((ii < p_component) && (jj == t_component)) {
            	coupling[jj][ii] = DoFTools::none;
            	coupling[ii][jj] = DoFTools::none;
            }
            else {
            	coupling[ii][jj] = DoFTools::always;
            	coupling[jj][ii] = DoFTools::always;
            }
        }
    }

    // Setup system matrix
    tangent_matrix.clear ();
    {
	const unsigned int n_dofs_u = dofs_per_block[u_dof];
	const unsigned int n_dofs_p = dofs_per_block[p_dof];
	const unsigned int n_dofs_t = dofs_per_block[t_dof];

        BlockCompressedSimpleSparsityPattern csp (n_blocks, n_blocks);

        csp.block(u_dof,u_dof).reinit (n_dofs_u, n_dofs_u);
        csp.block(u_dof,p_dof).reinit (n_dofs_u, n_dofs_p);
        csp.block(u_dof,t_dof).reinit (n_dofs_u, n_dofs_t);

        csp.block(p_dof,u_dof).reinit (n_dofs_p, n_dofs_u);
        csp.block(p_dof,p_dof).reinit (n_dofs_p, n_dofs_p);
        csp.block(p_dof,t_dof).reinit (n_dofs_p, n_dofs_t);

        csp.block(t_dof,u_dof).reinit (n_dofs_t, n_dofs_u);
        csp.block(t_dof,p_dof).reinit (n_dofs_t, n_dofs_p);
        csp.block(t_dof,t_dof).reinit (n_dofs_t, n_dofs_t);
        csp.collect_sizes();

	DoFTools::make_sparsity_pattern (dof_handler_ref, csp);
	//        DoFTools::make_sparsity_pattern (dof_handler_ref, csp, constraints, false);
	//        DoFTools::make_sparsity_pattern (dof_handler_ref, coupling, csp, constraints, false);
        sparsity_pattern.copy_from (csp);
    }


    tangent_matrix.reinit (sparsity_pattern);

    // Setup storage vectors
    residual.reinit (dofs_per_block);
    residual.collect_sizes ();

    solution_n.reinit (dofs_per_block);
    solution_n.collect_sizes ();
    solution_n.block(t_dof) = 1.0; // Dilatation is 1 in the initial configuration

    // Set up the quadrature point history
    setup_qph ();

    timer.leave_subsection();
}

// A way to extract subblocks from the matrix
template <int dim>
void Solid<dim>::determine_component_extractors(void)
{
    element_indices_u.clear();
    element_indices_p.clear();
    element_indices_t.clear();

    for (unsigned int k=0; k < fe.dofs_per_cell; ++k) {
	// 0 = u, 1 = p, 2 = dilatation interpolation fields
	const unsigned int k_group = fe.system_to_base_index(k).first.first;
	if (k_group == u_dof) {
	    element_indices_u.push_back(k);
	}
	else if (k_group == p_dof) {
	    element_indices_p.push_back(k);
	}
	else if (k_group == t_dof) {
	    element_indices_t.push_back(k);
	}
	else {
	    Assert (k_group <= t_dof, ExcInternalError());
	}
    }
}

// @sect4{Solid::setup_qph}
template <int dim>
void Solid<dim>::setup_qph (void)
{
    std::cout << "    Setting up quadrature point data..." << std::endl;

    {
    	typename Triangulation<dim>::active_cell_iterator
		cell = triangulation.begin_active(),
		endc = triangulation.end();

    	unsigned int our_cells = 0;
    	for (; cell != endc; ++cell) {
    	    cell->clear_user_pointer();
    	    ++our_cells;
    	}

    	{
    	    std::vector<PointHistory <dim> > tmp;
    	    tmp.swap(quadrature_point_history);
    	}

    	quadrature_point_history.resize(our_cells * n_q_points);

    	unsigned int history_index = 0;
    	for (cell = triangulation.begin_active(); cell != endc; ++cell) {
    	    cell->set_user_pointer(&quadrature_point_history[history_index]);
    	    history_index += n_q_points;
    	}

    	Assert(history_index == quadrature_point_history.size(), ExcInternalError());
    }

    // Setup initial data
    typename DoFHandler<dim>::active_cell_iterator
    	    cell = dof_handler_ref.begin_active(),
    	    endc = dof_handler_ref.end();
    for (; cell != endc; ++cell) {
    	PointHistory<dim>* lqph = reinterpret_cast<PointHistory<dim>*> (cell->user_pointer());
    	Assert(lqph >= &quadrature_point_history.front(), ExcInternalError());
    	Assert(lqph < &quadrature_point_history.back(), ExcInternalError());

    	// Setup any initial information at displacement gauss points
    	for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
	    lqph[q_point].setup_lqp( parameters );
    	}
    }
}

// @sect4{Solid::update_qph_incremental}
template <int dim>
void Solid<dim>::update_qph_incremental (const BlockVector <double> & solution_delta)
{
    timer.enter_subsection("Update QPH data");
    std::cout << "Update QPH data..."<< std::endl;

    // Get total solution as it stands at this update increment
    const BlockVector <double> solution_total = get_solution_total(solution_delta);
    const UpdateFlags uf_UQPH ( update_values | update_gradients );
    PerTaskData_UQPH per_task_data_UQPH;
    ScratchData_UQPH scratch_data_UQPH (fe,
					qf_cell,
					uf_UQPH,
					solution_total);

    WorkStream::run (  dof_handler_ref.begin_active(),
		     dof_handler_ref.end(),
		     *this,
		     &Solid::update_qph_incremental_one_cell,
		     &Solid::copy_local_to_global_UQPH,
		     scratch_data_UQPH,
		     per_task_data_UQPH);

    timer.leave_subsection();
}

template <int dim>
void Solid<dim>::update_qph_incremental_one_cell (const typename DoFHandler<dim>::active_cell_iterator & cell,
						  ScratchData_UQPH & scratch,
						  PerTaskData_UQPH & data)
{
    PointHistory<dim>* lqph = reinterpret_cast<PointHistory<dim>*> (cell->user_pointer());
    Assert(lqph >= &quadrature_point_history.front(), ExcInternalError());
    Assert(lqph < &quadrature_point_history.back(), ExcInternalError());

    Assert(scratch.solution_grads_u_total.size()  == n_q_points, ExcInternalError());
    Assert(scratch.solution_values_p_total.size() == n_q_points, ExcInternalError());
    Assert(scratch.solution_values_t_total.size() == n_q_points, ExcInternalError());

    // Find the values and gradients at quadrature points inside the current cell
    scratch.fe_values_ref.reinit(cell);
    scratch.fe_values_ref[u_fe].get_function_gradients (scratch.solution_total, scratch.solution_grads_u_total);
    scratch.fe_values_ref[p_fe].get_function_values (scratch.solution_total, scratch.solution_values_p_total);
    scratch.fe_values_ref[t_fe].get_function_values (scratch.solution_total,scratch. solution_values_t_total);

    // === UPDATE DATA AT EACH GAUSS POINT ===
    // Update displacement and deformation gradient at all quadrature points
    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
	lqph[q_point].update_values (scratch.solution_grads_u_total [q_point],
				     scratch.solution_values_p_total[q_point],
				     scratch.solution_values_t_total[q_point]);
    }
}

// @sect4{Solid::solve_nonlinear_timestep}
template <int dim>
void Solid<dim>::solve_nonlinear_timestep (BlockVector <double> & solution_delta)
{
    //    timer.enter_subsection("Nonlinear solver");
    std::cout
	    << "Timestep " << time.get_timestep()
	    << std::endl;

    // Newton update vector
    BlockVector <double> newton_update (dofs_per_block);
    newton_update.collect_sizes ();

    // Solution error vectors
    BlockVector <double> soln_error_res (dofs_per_block); // Holds the true residual vector
    BlockVector <double> soln_error_update (dofs_per_block); // Holds the update error vector
    soln_error_res.collect_sizes ();
    soln_error_update .collect_sizes ();

    double res_u = 0.0, res_f = 0.0;
    double res_u_0 = 1.0, res_f_0 = 1.0;
    for (unsigned int it_nr=0; it_nr < parameters.max_iterations_NR; ++ it_nr)
    {
	std::cout
		<< std::endl
		<< "Newton iteration: " << it_nr
		<< std::endl;

	tangent_matrix = 0.0;
	residual = 0.0;

	// Check residual
	make_constraints (it_nr, constraints); // Make boundary conditions
	assemble_system_F (); // Assemble RHS
	get_error_res(residual, soln_error_res);
	// Residual scaling factors
	res_f = soln_error_res.block(u_dof).l2_norm();
	if (it_nr == 0) res_f_0 = res_f;

	// Check for solution convergence
	if (   it_nr > 0
		&& res_u/res_u_0 <= parameters.tol_u
		&& res_f/res_f_0 <= parameters.tol_f)
	{
	    std::cout
		    << std::endl
		    << "Solution for timestep " << time.get_timestep()
		    << " converged on Newton iteration " << it_nr-1 << "."
		    << std::endl
		    << "Relative displacement error: " << res_u/res_u_0
		    << "\t Relative force error: " << res_f/res_f_0
		    << "\t Dilatation error: " << get_error_dil()
		    << std::endl << std::endl;

	    //	    timer.leave_subsection();
	    return;
	}

	// No convergence -> continue with calculations
	// Assemble stiffness matrix
	assemble_system_K ();

	// Do the static condensation to make K'_uu, and put K_pt^{-1}
	// in the K_pt block and K_tt^{-1} in the K_pp block
        assemble_SC();


        // Do the static condensation to make K'_uu, and put K_pt^{-1}
        // in the K_pt block and K_tt^{-1} in the K_pp block
        assemble_SC();

	constraints.condense (tangent_matrix, residual); // Apply BC's
	solve_linear_system (newton_update);
	constraints.distribute(newton_update); // Populate the constrained DOF's with their values

	// Newton update error
	get_error_update(newton_update, soln_error_update);
	res_u = soln_error_update.block(u_dof).l2_norm();

	// Residual scaling factors
	if (it_nr == 0) res_u_0 = res_u;
	std::cout
		<< "Nonlinear system error: "
		<< std::endl << std::scientific
		<< " Solution update \t ||dU||: " << soln_error_update.l2_norm()
                << "\t ||dU_u||: " << soln_error_update.block(u_dof).l2_norm()
                << "\t ||dU_p||: " << soln_error_update.block(p_dof).l2_norm()
                << "\t ||dU_t||: " << soln_error_update.block(t_dof).l2_norm()
		<< std::endl;
	std::cout << std::scientific
		  << " Residual     \t ||dF||: " << soln_error_res.l2_norm()
		  << "\t ||dR_u||: " << soln_error_res.block(u_dof).l2_norm()
		  << "\t ||dR_p||: " << soln_error_res.block(p_dof).l2_norm()
		  << "\t ||dR_t||: " << soln_error_res.block(t_dof).l2_norm()
		  << std::endl;
	std::cout << std::scientific
		  << " Relative displacement error: " << res_u/res_u_0
		  << "\t Relative force error: " << res_f/res_f_0
		  << "\t Dilatation error: " << get_error_dil()
		  << std::endl;

	// Update and continue iterating
	solution_delta += newton_update; // Update current solution
	update_qph_incremental (solution_delta); // Update quadrature point information
    }

    throw(ExcMessage("No convergence in nonlinear solver!"));
}

template <int dim>
void Solid<dim>::get_error_res (const BlockVector <double> & residual, BlockVector <double> & error_res)
{
    for (unsigned int i=0; i < dof_handler_ref.n_dofs(); ++i)
    	if (!constraints.is_constrained(i))
	    error_res(i) = residual(i);
}

template <int dim>
void Solid<dim>::get_error_update (const BlockVector <double> & newton_update, BlockVector <double> & error_update)
{
    for (unsigned int i=0; i < dof_handler_ref.n_dofs(); ++i)
    	if (!constraints.is_constrained(i))
	    error_update(i) = newton_update(i);
}

template <int dim>
double Solid<dim>::get_error_dil (void)
{
    double v_e = 0.0; // Volume in current configuration
    double V_e = 0.0; // Volume in reference configuration

    FEValues<dim> fe_values_ref (fe, qf_cell, update_JxW_values);

    typename DoFHandler<dim>::active_cell_iterator
            cell = dof_handler_ref.begin_active(),
            endc = dof_handler_ref.end();
    for (; cell != endc; ++cell) {
        fe_values_ref.reinit (cell);
        PointHistory<dim>* lqph = reinterpret_cast<PointHistory<dim>*> (cell->user_pointer());
    	Assert(lqph >= &quadrature_point_history.front(), ExcInternalError());
    	Assert(lqph < &quadrature_point_history.back(), ExcInternalError());

        for (unsigned int q_point=0; q_point < n_q_points; ++q_point) {
            v_e += lqph[q_point].get_dilatation() * fe_values_ref.JxW(q_point);
            V_e += fe_values_ref.JxW(q_point);
        }
    }

    return std::abs((v_e - V_e)/V_e); // Difference between initial and current volume
}

// Solution (valid at any Newton step)
template <int dim>
BlockVector <double> Solid<dim>::get_solution_total (const BlockVector <double> & solution_delta)
{
    BlockVector <double> solution_total (solution_n);
    solution_total += solution_delta;

    return solution_total;
}

// @sect4{Solid::solve_linear_system}
template <int dim>
void Solid<dim>::solve_linear_system (BlockVector <double> & newton_update)
{
    std::cout << "Solve linear system..." << std::endl;

    BlockVector <double> A (dofs_per_block);
    BlockVector <double> B (dofs_per_block);
    A.collect_sizes ();
    B.collect_sizes ();

    //      | K'_uu |   K_up    |     0     |         | dU_u |         | dR_u |
    // K =  | K_pu  |   K_tt^-1 |   K_pt^-1 | , dU =  | dU_p | , dR =  | dR_p |
    //      |   0   |   K_tp    |   K_tt    |         | dU_t |         | dR_t |

    // Solve for du
    {

	// K'uu du = Ru − Kup Ktp^-1 (Rt − Ktt Kpt^{-1} Rp)
	tangent_matrix.block(p_dof, t_dof).vmult(A.block(t_dof), residual.block(p_dof));
	tangent_matrix.block(t_dof, t_dof).vmult (B.block(t_dof), A.block(t_dof));
	A.block(t_dof).equ(1.0, residual.block(t_dof), -1.0, B.block(t_dof));
	tangent_matrix.block(p_dof, t_dof).Tvmult(A.block(p_dof), A.block(t_dof));
	tangent_matrix.block(u_dof, p_dof).vmult(A.block(u_dof), A.block(p_dof));
	residual.block(u_dof) -= A.block(u_dof);

	timer.enter_subsection("Linear solver");
	if (parameters.type_lin == "CG")
	{
	    const int solver_its = tangent_matrix.block(u_dof, u_dof).m() * parameters.max_iterations_lin;
	    const double tol_sol = parameters.tol_lin * residual.block(u_dof).l2_norm();

	    SolverControl solver_control (solver_its , tol_sol);

	    GrowingVectorMemory < Vector<double> > GVM;
	    SolverCG < Vector<double> >  solver_CG (solver_control, GVM);

	    // SSOR -> much better than Jacobi for symmetric systems
	    PreconditionSSOR <SparseMatrix<double> > preconditioner;
	    preconditioner.initialize (tangent_matrix.block(u_dof, u_dof), parameters.ssor_relaxation);

	    solver_CG.solve (tangent_matrix.block(u_dof, u_dof),
			     newton_update.block(u_dof),
			     residual.block(u_dof),
			     preconditioner);

	    std::cout
		    << "\t Iterations: " << solver_control.last_step()
		    << "\n\t Residual: " << solver_control.last_value()
		    << std::endl;
	}
	else if (parameters.type_lin == "Direct")
	{
	    SparseDirectUMFPACK  A_direct;
	    A_direct.initialize(tangent_matrix.block(u_dof, u_dof));
	    A_direct.vmult (newton_update.block(u_dof),
			    residual.block(u_dof));
	}
	else throw (ExcMessage("Linear solver type not implemented"));
	timer.leave_subsection();
    }

    timer.enter_subsection("Linear solver postprocessing");
    // Postprocess for dp
    {
	// dp = Ktp^{-1} ( Rt − Ktt Kpt^{-1} (Rp − Kpu du) )
	tangent_matrix.block(p_dof, u_dof).vmult (A.block(p_dof), newton_update.block(u_dof));
	B.block(p_dof).equ(1.0, residual.block(p_dof), -1.0, A.block(p_dof));
	tangent_matrix.block(p_dof, t_dof).vmult(A.block(t_dof), B.block(p_dof));
	tangent_matrix.block(t_dof, t_dof).vmult(B.block(t_dof), A.block(t_dof));
	A.block(t_dof).equ (1.0, residual.block(t_dof), -1.0, B.block(t_dof));
	tangent_matrix.block(p_dof, t_dof).Tvmult (newton_update.block(p_dof), A.block(t_dof));
    }

    // Postprocess for dt
    {
	// dt = Ktt^{-1} (Rt − Ktp dp)
	tangent_matrix.block(t_dof, p_dof).vmult (A.block(t_dof), newton_update.block(p_dof));
	residual.block(t_dof) -= A.block(t_dof);
	tangent_matrix.block(p_dof, p_dof).vmult (newton_update.block(t_dof), residual.block(t_dof));
    }
    timer.leave_subsection();
}

// @sect4{Solid::assemble_system_K}
template <int dim>
void Solid<dim>::assemble_system_K (void)
{
    timer.enter_subsection("Assemble system matrix");
    std::cout << "Assemble system matrix..."<< std::endl;

    tangent_matrix = 0.0; // Clear the matrix

    const UpdateFlags uf_cell ( update_values | update_gradients | update_JxW_values  );

    PerTaskData_K per_task_data (dofs_per_cell); // Initialise members of per_task_data to the correct sizes.
    ScratchData_K scratch_data (fe, qf_cell, uf_cell);

    WorkStream::run (  dof_handler_ref.begin_active(),
		     dof_handler_ref.end(),
		     *this,
		     &Solid::assemble_system_K_one_cell,
		     &Solid::copy_local_to_global_K,
		     scratch_data,
		     per_task_data);

    timer.leave_subsection();
}

template <int dim>
void Solid<dim>::copy_local_to_global_K (const PerTaskData_K & data)
{
    // Add the local contribution to the system matrix
    for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	    tangent_matrix.add (data.local_dof_indices[i],
				data.local_dof_indices[j],
				data.cell_matrix(i,j));
}

template <int dim>
void Solid<dim>::assemble_system_K_one_cell (const typename DoFHandler<dim>::active_cell_iterator & cell,
					     ScratchData_K & scratch,
					     PerTaskData_K & data)
{
    data.reset(); // Reset data in the PerTaskData_K storage unit
    scratch.reset(); // Reset data in the Scratch storage unit
    scratch.fe_values_ref.reinit (cell);
    cell->get_dof_indices (data.local_dof_indices); // Find out which global numbers the degrees of freedom on this cell have
    PointHistory<dim> *lqph = reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

    // Set up cell shape function gradients
    static const SymmetricTensor<2, dim> I = unit_symmetric_tensor <dim> ();
    for (unsigned int q_point=0; q_point < n_q_points; ++q_point) {
	const Tensor<2, dim> F_inv = lqph[q_point].get_F_inv();

	for (unsigned int k=0; k< dofs_per_cell; ++k) {
	    const unsigned int k_group = fe.system_to_base_index(k).first.first;

	    if (k_group == u_dof) {
		scratch.grad_Nx[q_point][k] = scratch.fe_values_ref[u_fe].gradient(k, q_point) * F_inv;
		scratch.symm_grad_Nx[q_point][k] = symmetrize(scratch.grad_Nx[q_point][k]);
	    }
	    else if (k_group == p_dof) {
		scratch.Nx[q_point][k] = scratch.fe_values_ref[p_fe].value(k, q_point);
	    }
	    else if (k_group == t_dof) {
		scratch.Nx[q_point][k] = scratch.fe_values_ref[t_fe].value(k, q_point);
	    }
	    else {
		Assert (k_group <= t_dof, ExcInternalError());
	    }
	}
    }

    // Build cell stiffness matrix
    // Global and local system matrices are symmetric
    //  => Take advantage of this:  Build only the lower half of the local matrix
    // Only assemble 1/2 of the K_uu, K_pp = 0, K_tt blocks and the whole K_pt, K_ut, K_up blocks
    for (unsigned int q_point=0; q_point < n_q_points; ++q_point) {
	const Tensor <2,dim>          T   = static_cast < Tensor<2, dim> > (lqph[q_point].get_T_iso() + lqph[q_point].get_T_vol());
	const SymmetricTensor <4,dim> C   = lqph[q_point].get_C_iso() + lqph[q_point].get_C_vol();
	const double                  C_v = lqph[q_point].get_d2U_dtheta2();
	const double                  J   = lqph[q_point].get_J();

	const std::vector<double> & N = scratch.Nx[q_point];
	const std::vector< SymmetricTensor <2,dim> > & symm_B = scratch.symm_grad_Nx[q_point];
	const std::vector< Tensor <2,dim> > & B = scratch.grad_Nx[q_point];
	const double & JxW = scratch.fe_values_ref.JxW(q_point);

	for (unsigned int i=0; i < dofs_per_cell; ++i) {

	    const unsigned int component_i = fe.system_to_component_index(i).first;
	    const unsigned int i_group = fe.system_to_base_index(i).first.first;

	    // Only assemble the lower diagonal part of the local matrix
	    for (unsigned int j=0; j <= i; ++j) {

		const unsigned int component_j = fe.system_to_component_index(j).first;
		const unsigned int j_group = fe.system_to_base_index(j).first.first;

		if (   (i_group == j_group) && (i_group == u_dof ) ) {
		    data.cell_matrix(i,j)
			    += ( symm_B[i] * C * symm_B[j]   // Material stiffness
				+  ( component_i == component_j ?
					B[i][component_i] * T * B[j][component_j]  :
					0.0 ) // Geometric stiffness. Only add this along local diagonals
				) * JxW;  // K_uu
		}
		else if ( (i_group == p_dof) && (j_group == u_dof) ) {
		    data.cell_matrix(i,j) += N[i]*J*(symm_B[j]*I)*JxW; // K_pu
		}
		else if ( (i_group == t_dof) && (j_group == p_dof) ) {
		    data.cell_matrix(i,j) -= N[i]*N[j]*JxW; // K_tp
		}
		else if ( (i_group == j_group) && (i_group == t_dof)  ) {
		    data.cell_matrix(i,j) +=  N[i]*C_v*N[j]*JxW; // K_tt
		}
		else Assert ((i_group <= t_dof) && (j_group <= t_dof), ExcInternalError());
	    } // END j LOOP
	} // END i LOOP

    } // END q_point LOOP

    // Global and local system matrices are symmetric
    // => Copy the upper half of the local matrix in the bottom  half of the local matrix
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
	for (unsigned int j=i+1; j<dofs_per_cell; ++j) {
	    data.cell_matrix(i,j) = data.cell_matrix(j,i);
	}
    }
}

// @sect4{Solid::assemble_system_F}
template <int dim>
void Solid<dim>::assemble_system_F (void)
{
    timer.enter_subsection("Assemble system RHS");
    std::cout << "Assemble system RHS..."<< std::endl;

    residual  = 0.0; // Clear the vector

    const UpdateFlags uf_cell ( update_values | update_gradients | update_JxW_values );
    const UpdateFlags uf_face ( update_values | update_normal_vectors | update_JxW_values);

    PerTaskData_F per_task_data (dofs_per_cell); // Initialise members of per_task_data to the correct sizes.
    ScratchData_F scratch_data (fe,
				qf_cell,
				uf_cell,
				qf_face,
				uf_face);

    WorkStream::run ( dof_handler_ref.begin_active(),
		     dof_handler_ref.end(),
		     *this,
		     &Solid::assemble_system_F_one_cell,
		     &Solid::copy_local_to_global_F,
		     scratch_data,
		     per_task_data );

    timer.leave_subsection();
}

template <int dim>
void Solid<dim>::copy_local_to_global_F (const PerTaskData_F & data)
{
    // Add the local contribution to the system RHS vector
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
	residual(data.local_dof_indices[i]) += data.cell_rhs(i);
    }
}

template <int dim>
void Solid<dim>::assemble_system_F_one_cell (const typename DoFHandler<dim>::active_cell_iterator & cell,
					     ScratchData_F & scratch,
					     PerTaskData_F & data)
{
    data.reset(); // Reset data in the PerTaskData_K storage unit
    scratch.reset(); // Reset data in the ScratchData_F storage unit
    scratch.fe_values_ref.reinit (cell);
    cell->get_dof_indices (data.local_dof_indices); // Find out which global numbers the degrees of freedom on this cell have
    PointHistory<dim> *lqph = reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

    // Precompute some data
    for (unsigned int q_point=0; q_point < n_q_points; ++q_point) {
	const Tensor<2, dim> F_inv = lqph[q_point].get_F_inv();

	for (unsigned int k=0; k<dofs_per_cell; ++k) {
	    const unsigned int k_group = fe.system_to_base_index(k).first.first;

	    if (k_group == u_dof) {
		scratch.symm_grad_Nx[q_point][k] = symmetrize(scratch.fe_values_ref[u_fe].gradient(k, q_point) * F_inv);
	    }
	    else if (k_group == p_dof) {
		scratch.Nx[q_point][k] = scratch.fe_values_ref[p_fe].value(k, q_point);
	    }
	    else if (k_group == t_dof) {
		scratch.Nx[q_point][k] = scratch.fe_values_ref[t_fe].value(k, q_point);
	    }
	    else Assert (k_group <= t_dof, ExcInternalError());
	}
    }

    // Assembly for residual contribution
    for (unsigned int q_point=0; q_point < n_q_points; ++q_point) {
	const SymmetricTensor <2,dim>  T = lqph[q_point].get_T_iso() + lqph[q_point].get_T_vol();
	const double  J = lqph[q_point].get_J();
	const double  D = lqph[q_point].get_dilatation();
	const double  p = lqph[q_point].get_pressure();
	const double  p_star = lqph[q_point].get_dU_dtheta();

	const std::vector< double > & N = scratch.Nx[q_point];
	const std::vector< SymmetricTensor <2,dim> > & symm_B = scratch.symm_grad_Nx[q_point];
	const double  JxW = scratch.fe_values_ref.JxW(q_point);

	for (unsigned int i=0; i<dofs_per_cell; ++i) {
	    const unsigned int i_group = fe.system_to_base_index(i).first.first;

	    if (i_group == u_dof) {
		data.cell_rhs(i) -= ( symm_B[i]*T )*JxW; // R_u
	    }
	    else if (i_group == p_dof ) {
		data.cell_rhs(i) -= N[i]*(J - D)*JxW; // R_p
	    }
	    else if ( i_group == t_dof) {
		data.cell_rhs(i) -= N[i]*(p_star-p)*JxW; // R_t
	    }
	    else Assert (i_group <= t_dof, ExcInternalError());
	} // END i LOOP
    } // END q_point LOOP

    // Assembly for Neumann RHS contribution
    if (cell->at_boundary() == true)
    {
	static const Tensor <2, dim> I = static_cast < Tensor <2, dim> > ( unit_symmetric_tensor <dim> () );

	for (unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; ++face)
	{
	    if (    cell->face(face)->at_boundary() == true
		    &&  cell->face(face)->boundary_indicator() == 6 )
	    {
		scratch.fe_face_values_ref.reinit (cell, face);

		for (unsigned int f_q_point=0; f_q_point < n_q_points_f; ++f_q_point)
		{
		    const Tensor <1, dim> & N = scratch.fe_face_values_ref.normal_vector(f_q_point);

		    // Traction in reference configuration
		    // t_0 = p*N
		    static const double p0 = -4.0/(parameters.scale*parameters.scale); // Reference pressure of 4 Pa
		    const double time_ramp = (time.current() / time.end()); // Linearly ramp up the pressure with time
		    const double pressure = p0 * parameters.p_p0 * time_ramp;
		    const Tensor <1,dim> traction = pressure * N;

		    for (unsigned int i=0; i < dofs_per_cell; ++i) {
			// Determine the dimensional component that matches the dof component (i.e. i % dim)
			const unsigned int i_group = fe.system_to_base_index(i).first.first;

			if (i_group == u_dof) {
			    const unsigned int component_i = fe.system_to_component_index(i).first;
			    const double & Ni = scratch.fe_face_values_ref.shape_value(i,f_q_point);
			    const double & JxW = scratch.fe_face_values_ref.JxW(f_q_point);

			    // Add traction vector contribution to the local RHS vector (displacement dofs only)
			    data.cell_rhs(i) += (Ni * traction[component_i])  // Contribution from external forces
				    * JxW;
			}
		    } // END i LOOP
		} // END face q_point LOOP
	    } // END at boundary check LOOP

	} // END face LOOP
    }
}

// @sect4{Solid::assemble_system_SC}
template <int dim>
void Solid<dim>::assemble_SC  (void)
{
    timer.enter_subsection("Perform static condensation");

    PerTaskData_SC per_task_data (dofs_per_cell,
				  element_indices_u.size(),
				  element_indices_p.size(),
				  element_indices_t.size()); // Initialise members of per_task_data to the correct sizes.
    ScratchData_SC scratch_data;

    WorkStream::run (  dof_handler_ref.begin_active(),
		     dof_handler_ref.end(),
		     *this,
		     &Solid::assemble_SC_one_cell,
		     &Solid::copy_local_to_global_SC,
		     scratch_data,
		     per_task_data  );

    timer.leave_subsection();
}

template <int dim>
void Solid<dim>::copy_local_to_global_SC (const PerTaskData_SC & data)
{
    // Add the local contribution to the system matrix
    for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	    tangent_matrix.add (data.local_dof_indices[i],
				data.local_dof_indices[j],
				data.cell_matrix(i,j));
}

template <int dim>
void Solid<dim>::assemble_SC_one_cell (const typename DoFHandler<dim>::active_cell_iterator & cell,
				       ScratchData_SC & scratch,
				       PerTaskData_SC & data)
{
    data.reset();
    scratch.reset();
    cell->get_dof_indices (data.local_dof_indices); // Find out which global numbers the degrees of freedom on this cell have

    // The local stifness matrix K_e is:
    //  | K_uu  |   K_up   |   0   |
    //  | K_pu  |     0    |  K_pt |
    //  |   0   |   K_tp   |  K_tt |
    //
    // We are going to exploit the zeros for post-processing as:
    //  | K'_uu |   K_up    |     0     |
    //  | K_pu  |   K_tt^-1 |   K_pt^-1 |
    //  |   0   |   K_tp    |   K_tt    |
    // with K'_uu = K_uu + Kup Ktp^{-1} Ktt Kpt^{-1} Kpu

    // NOTE:
    // GLOBAL Data already exists in the K_uu, K_pt, K_tp subblocks
    //
    // For the K_uu block in particular, this means that contributions have been
    // added from the surrounding cells, so we need to be careful when we manipulate this block.
    // We can't just erase the subblocks and
    // Additionally the copy_local_to_global operation is a "+=" operation -> need to take this
    // into account
    //
    // So the intermediate matrix that we need to get from what we have in K_uu and what we
    // are actually wanting is:
    //  | K'_uu - K_uu |     0    |        0        |
    //  |       0      |  K_tt^-1 |  K_pt^-1 - K_pt |
    //  |       0      |     0    |        0        |
    //
    // Strategy to get the subblocks we want:
    // K'_uu: Since we don't have access to K_uu^h, but we know its contribution is added to the global
    //        K_uu matrix, we just want to add the element wise static-condensation
    //        K'_uu^h = K_uu^h + K_up^h K_tp^{-1}^h K_tt^h K_pt^{-1}^h K_pu^h
    //        Since we already have K_uu^h in the system matrix, we just need to do the following
    //        K'_uu^h == (K_uu^h += K_up^h K_tp^{-1}^h K_tt^h K_pt^{-1}^h K_pu^h)
    // K_pt^-1: Similarly, K_pt exists in the subblock. Since the copy operation is a += operation, we need
    //          to subtract the existing K_pt submatrix in addition to "adding" that which we wish to
    //          replace it with.
    // K_tp^-1: Same as above
    // K_tt^-1: Nothing exists in the original K_pp subblock, so we can just add this contribution as is.

    // Extract element data from the system matrix

    AdditionalTools::extract_submatrix(data.local_dof_indices,
				       data.local_dof_indices,
				       tangent_matrix,
				       data.K_orig);
    AdditionalTools::extract_submatrix(element_indices_p,
				       element_indices_u,
				       data.K_orig,
				       data.K_pu);
    AdditionalTools::extract_submatrix(element_indices_p,
				       element_indices_t,
				       data.K_orig,
				       data.K_pt);
    AdditionalTools::extract_submatrix(element_indices_t,
				       element_indices_t,
				       data.K_orig,
				       data.K_tt);

    // Place K_pt^-1 in the K_pt block
    data.K_pt_inv.invert(data.K_pt);
    data.K_pt_inv.add (-1.0, data.K_pt);
    AdditionalTools::replace_submatrix(element_indices_p,
				       element_indices_t,
				       data.K_pt_inv,
				       data.cell_matrix);

    // Place K_tt^-1 in the K_pp block
    data.K_tt_inv.invert(data.K_tt);
    AdditionalTools::replace_submatrix(element_indices_p,
				       element_indices_p,
				       data.K_tt_inv,
				       data.cell_matrix);

    // Make condensation terms to add to the K_uu block
    data.K_pt_inv.mmult(data.A, data.K_pu);
    data.K_tt.mmult(data.B, data.A);
    data.K_pt_inv.Tmmult(data.C, data.B); // Symmetric matrix
    data.K_pu.Tmmult(data.K_con, data.C); // Symmetric matrix
    AdditionalTools::replace_submatrix(element_indices_u,
				       element_indices_u,
				       data.K_con,
				       data.cell_matrix);
}

// @sect4{Solid::make_constraints}
template <int dim>
void Solid<dim>::make_constraints (const int & it_nr,
				   ConstraintMatrix & constraints)
{
    std::cout << "Make constraints..."<< std::endl;

    constraints.clear();
    const bool apply_dirichlet_bc = (it_nr == 0);

    // Boundary conditions:
    // b_id 0: -x face: Zero x-component of displacement : Symmetry plane
    // b_id 2: -y face: Zero y-component of displacement : Symmetry plane
    // b_id 4: -z face: Zero z-component of displacement : Symmetry plane

    // b_id 5: +z face: Zero x-component and Zero y-component
    // b_id 6: Applied pressure face: Zero x-component and Zero y-component
    // b_id 1: +x face: Traction free
    // b_id 3: +y face: Traction free
    {
	const int boundary_id = 0;

	std::vector< bool > components (n_components, false);
	components[0] = true;

	if (apply_dirichlet_bc == true) {
	    VectorTools::interpolate_boundary_values ( dof_handler_ref, boundary_id, ZeroFunction<dim>(n_components), constraints, components );
	}
	else {
	    VectorTools::interpolate_boundary_values ( dof_handler_ref, boundary_id, ZeroFunction<dim>(n_components), constraints, components );
	}
    }
    {
	const int boundary_id = 2;

	std::vector< bool > components (n_components, false);
	components[1] = true;

	if (apply_dirichlet_bc == true) {
	    VectorTools::interpolate_boundary_values ( dof_handler_ref, boundary_id, ZeroFunction<dim>(n_components), constraints, components );
	}
	else {
	    VectorTools::interpolate_boundary_values ( dof_handler_ref, boundary_id, ZeroFunction<dim>(n_components), constraints, components );
	}
    }
    {
	const int boundary_id = 4;
	std::vector< bool > components (n_components, false);
	components[2] = true;

	if (apply_dirichlet_bc == true) {
	    VectorTools::interpolate_boundary_values ( dof_handler_ref, boundary_id, ZeroFunction<dim>(n_components), constraints, components );
	}
	else {
	    VectorTools::interpolate_boundary_values ( dof_handler_ref, boundary_id, ZeroFunction<dim>(n_components), constraints, components );
	}
    }
    {
	const int boundary_id = 5;
	std::vector< bool > components (n_components, true);
	components[2] = false;

	if (apply_dirichlet_bc == true) {
	    VectorTools::interpolate_boundary_values ( dof_handler_ref, boundary_id, ZeroFunction<dim>(n_components), constraints, components );
	}
	else {
	    VectorTools::interpolate_boundary_values ( dof_handler_ref, boundary_id, ZeroFunction<dim>(n_components), constraints, components );
	}
    }
    {
	const int boundary_id = 6;
	std::vector< bool > components (n_components, true);
	components[2] = false;

	if (apply_dirichlet_bc == true) {
	    VectorTools::interpolate_boundary_values ( dof_handler_ref, boundary_id, ZeroFunction<dim>(n_components), constraints, components );
	}
	else {
	    VectorTools::interpolate_boundary_values ( dof_handler_ref, boundary_id, ZeroFunction<dim>(n_components), constraints, components );
	}
    }

    constraints.close();
}

// @sect4{Solid::output_results}
template <int dim>
void Solid<dim>::output_results(void)
{
    DataOut<dim> data_out;

    std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation (dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back (DataComponentInterpretation::component_is_scalar);
    data_component_interpretation.push_back (DataComponentInterpretation::component_is_scalar);

    std::vector<std::string> solution_name (dim, "displacement");
    solution_name.push_back ("pressure");
    solution_name.push_back ("dilatation");

    data_out.attach_dof_handler (dof_handler_ref);
    data_out.add_data_vector (solution_n,
			      solution_name,
			      DataOut<dim>::type_dof_data, data_component_interpretation);
    //    MappingQEulerian<dim> q_mapping (degree, solution_n.block(u_dof), dof_handler_ref);
    //    MappingQEulerian<dim> q_mapping (degree, solution_n, dof_handler_ref);
    Vector<double> soln;
    soln.reinit(solution_n.size());
    for (unsigned int i=0; i < soln.size(); ++i) soln(i) = solution_n(i);
    MappingQEulerian<dim> q_mapping (degree, soln, dof_handler_ref);
    data_out.build_patches (q_mapping,degree);

    std::ostringstream filename;
    filename << "solution-"
	     << time.get_timestep()
	     << ".vtk";

    std::ofstream output (filename.str().c_str());
    data_out.write_vtk (output);
}

// @sect3{Main function}
int main ()
{
    try
    {
	deallog.depth_console (0);

	Solid<3> solid_3d ("parameters.prm");
	solid_3d.run();
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
    }

    return 0;
}



/* Authors: Jean-Paul Pelteret, University of Cape Town,            */
/*          Andrew McBride, University of Erlangen-Nuremberg, 2010  */
/*                                                                  */
/*    Copyright (C) 2010 by the deal.II authors                     */
/*                        & Jean-Paul Pelteret and Andrew McBride   */
/*                                                                  */
/*    This file is subject to QPL and may not be  distributed       */
/*    without copyright and license information. Please refer       */
/*    to the file deal.II/doc/license.html for the  text  and       */
/*    further information on this license.                          */

// @sect3{Include files}
// We start by including all the necessary
// deal.II header files and some C++ related
// ones. They have been discussed in detail
// in previous tutorial programs, so you need
// only refer to past tutorials for details.

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

#include <iostream>
#include <fstream>

// Next we import all the deal.II
// function and class names to the global namespace
using namespace dealii;

// @sect3{Run-time parameters}
//
// There are several parameters that can be set
// so we choose to set up a parameter
// handler object so that we can read in choices
// at run-time.
namespace Parameters
{
// @sect4{Finite Element system}
// Change the polynomial order used to approximate the solution.
// The quadrature should be adjusted accordingly.
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
        poly_degree = prm.get_integer("Polynomial degree");
        quad_order = prm.get_integer("Quadrature order");
    }
    prm.leave_subsection();
}

// @sect4{Geometry}
// Make adjustments to the problem geometry and the applied load.
// Since the problem modelled here is quite specific, the load
// scale can be altered to specific values to attain results given
// in the literature.
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
        scale = prm.get_double("Grid scale");
        p_p0 = prm.get_double("Pressure ratio p/p0");
    }
    prm.leave_subsection();
}

// @sect{Materials}
// Store the shear modulus and Lame constant
// for the Neo-Hookean material
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
        nu = prm.get_double("Poisson's ratio");
        mu = prm.get_double("Shear modulus");
    }
    prm.leave_subsection();
}

// @sect4{Linear solver}
// Choose both CG solver and SSOR preconditioner settings.
// The default values are optimal for this particular problem.
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
        max_iterations_lin = prm.get_double("Max iteration multiplier");
        ssor_relaxation = prm.get_double("SSOR Relaxation");
    }
    prm.leave_subsection();
}

// Nonlinear solver
// Define the tolerances and maximum number of iterations for the
// Newton-Raphson nono-linear solver.
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
        max_iterations_NR = prm.get_integer("Max iterations Newton-Raphson");
        tol_f = prm.get_double("Tolerance force");
        tol_u = prm.get_double("Tolerance displacement");
    }
    prm.leave_subsection();
}

// @sect4{Time}
// Set the timestep size and the simulation end-time.
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

// sect4{All parameters}
// Finally we consolidate all of the above structures into
// a single container that holds all of our run-time selections.
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

}  // End Parameters namespace

// @sect3{General tools}
// We need to perform some specific operations that are not defined
// in the deal.II library yet. We place these common operations
// in a seperate namespace for convenience.
namespace AdditionalTools
{
// Define an operation that takes two tensors \f$ \mathbf{A} \f$ and
// \f$ \mathbf{B} \f$ such that their outer-product
// \f$ \mathbf{A} \bar{\otimes} \mathbf{B} \Rightarrow C_{ijkl} = A_{ik} B_{jl} \f$
template <int dim>
SymmetricTensor<4,dim> outer_product_T23 (const SymmetricTensor<2,dim> & A,
                                          const SymmetricTensor<2,dim> & B)
{
    SymmetricTensor<4,dim> A_ik_B_jl;

    for (unsigned int i=0; i<dim; ++i) {
        for (unsigned int j=i; j<dim; ++j) {
            for (unsigned int k=0; k<dim; ++k) {
                for (unsigned int l=k; k<dim; ++k) {
                    A_ik_B_jl[i][j][k][l] += A[i][k] * B[j][l];
                }
            }
        }
    }

    return A_ik_B_jl;
}

// The \a extract_submatrix function takes specific entries from a \a matrix,
// and copies them to a \a sub_matrix. The copied entries are defined by the
// first two parameters which hold the row and column entries to be extracted.
// The \a matrix is automatically resized to size \f$ r \times c \f$.
template <typename MatrixType>
void extract_submatrix(const std::vector< unsigned int > &row_index_set,
                       const std::vector< unsigned int > &column_index_set,
                       const MatrixType &matrix,
                       FullMatrix< double > &sub_matrix)
{

    const unsigned int n_rows_submatrix = row_index_set.size();
    const unsigned int n_cols_submatrix = column_index_set.size();

    sub_matrix.reinit(n_rows_submatrix, n_cols_submatrix);

    for (unsigned int sub_row = 0; sub_row < n_rows_submatrix; ++sub_row) {
        const unsigned int row = row_index_set[sub_row];
        Assert (row<=matrix.m(), ExcInternalError());

        for (unsigned int sub_col = 0; sub_col < n_cols_submatrix; ++sub_col) {
            const unsigned int col = column_index_set[sub_col];
            Assert (col<=matrix.n(), ExcInternalError());

            sub_matrix(sub_row,sub_col) = matrix(row, col);
        }
    }
}

template <>
void extract_submatrix < dealii::BlockSparseMatrix <double> >(const std::vector< unsigned int > &row_index_set,
                                                             const std::vector< unsigned int > &column_index_set,
                                                             const dealii::BlockSparseMatrix <double> &matrix,
                                                             FullMatrix< double > &sub_matrix)
{

    const unsigned int n_rows_submatrix = row_index_set.size();
    const unsigned int n_cols_submatrix = column_index_set.size();

    sub_matrix.reinit(n_rows_submatrix, n_cols_submatrix);

    for (unsigned int sub_row = 0; sub_row < n_rows_submatrix; ++sub_row) {
        const unsigned int row = row_index_set[sub_row];
        Assert (row<=matrix.m(), ExcInternalError());

        for (unsigned int sub_col = 0; sub_col < n_cols_submatrix; ++sub_col) {
            const unsigned int col = column_index_set[sub_col];
            Assert (col<=matrix.n(), ExcInternalError());
            if (matrix.get_sparsity_pattern().exists(row, col) == false) continue;

            sub_matrix(sub_row,sub_col) = matrix(row, col);
        }
    }
}

// The \a replace_submatrix function takes specific entries from a \a matrix,
// and copies them to a \a sub_matrix. The copied entries are defined by the
// first two parameters which hold the row and column entries to be replaced.
// The \a matrix expected to be of the correct size.
template <typename MatrixType>
void replace_submatrix(const std::vector< unsigned int > &row_index_set,
                       const std::vector< unsigned int > &column_index_set,
                       const MatrixType &sub_matrix,
                       FullMatrix< double >  &matrix)
{
    const unsigned int n_rows_submatrix = row_index_set.size();
    Assert (n_rows_submatrix<=sub_matrix.m(), ExcInternalError());
    const unsigned int n_cols_submatrix = column_index_set.size();
    Assert (n_cols_submatrix<=sub_matrix.n(), ExcInternalError());

    for (unsigned int sub_row = 0; sub_row < n_rows_submatrix; ++sub_row) {
        const unsigned int row = row_index_set[sub_row];
        Assert (row<=matrix.m(), ExcInternalError());

        for (unsigned int sub_col = 0; sub_col < n_cols_submatrix; ++sub_col) {
            const unsigned int col = column_index_set[sub_col];
            Assert (col<=matrix.n(), ExcInternalError());

            matrix(row, col) = sub_matrix(sub_row, sub_col);

        }
    }
}

}

// @sect3{Time class}
// A simple class to store time data is created. Its
// functioning is transparent so no discussion is
// necessary.
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
// The entire domain is to made of a Neo-Hookean material
// with constant properties throughout. This class defines
// the behaviour of this material. Neo-Hookean materials
// can be described by a strain-energy function (SEF)
// \f$ \phi = \phi_{B} + \phi_{V}  \f$
// where the bulk deformation is given by
// \f$ \phi_{B} = C_{1} \left( I_{1} - 3 \right)  \f$
// where \f$ C_{1} - \frac{\mu}{2} \f$ and $I_{1}$ is the first
// invariant of the left- or right- Cauchy deformation tensors.
// In this example the SEF that governs the volumetric
// response is defined as
// \f$ \phi_{V} = \kappa \left( \frac{1}{2} \left( \theta^{2} - 1 \right) - ln \left( \theta \right) \right)  \f$
// where $\kappa$ is the bulk modulus.
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
    ~Material_NH (void) {}

    // The Kirchhoff stress tensor is required in the formulation
    // used in this work. This is obtained from the SEF by
    // \f$ \mathbf{T} = \mathbf{F} \frac{\partial \hat{\phi}}{\partial \mathbf{C}} \mathbf{F}^{T} = \frac{\partial \phi}{\partial \mathbf{B}} \mathbf{B} \f$
    SymmetricTensor<2, dim> get_T (const double & J,
                                   const SymmetricTensor <2, dim> & B)
    {
	const double dW_dJ  = get_dU_dtheta (J);
	return mu_0*B + dW_dJ*J*I;
    }

    // The tangent matrix for this material is also calculated from the SEF by
    // \f$ JC_{ijkl} = F_{iA} F_{jB} H_{ABCD} F_{kC} F_{lD}\f$
    // with
    // \f$ \mathbf{H} = \frac{\partial^{2} \hat{\phi}}{\partial \mathbf{C} \partial \mathbf{C}}
    SymmetricTensor<4, dim> get_JC (const double & J,
                                    const SymmetricTensor <2, dim> & B)
    {
	const double dW_dJ   = get_dU_dtheta (J);
	const double d2W_dJ2 = get_d2U_dtheta2 (J);
	return  J*(  (dW_dJ + J*d2W_dJ2)*IxI - (2.0*dW_dJ)*II  );
    }

    // From the volumetric strain-energy function we calculate the
    // first and second derivatives with respect to the dilatation
    double get_dU_dtheta    (const double & d) {return kappa_0*(d - 1.0/d);}
    double get_d2U_dtheta2  (const double & d) {return kappa_0*(1.0 + 1.0/(d*d));}

protected:
    // Material properties
    const double lambda_0; // Lame constant
    const double mu_0;     // Shear modulus
    const double kappa_0;  // Bulk modulus

    // We also choose to precalculate and store some frequently used
    // second and fourth-order tensors.
    static SymmetricTensor<2, dim> const I;
    static SymmetricTensor<4, dim> const IxI;
    static SymmetricTensor<4, dim> const II;
};

template <int dim> SymmetricTensor<2, dim> const Material_NH<dim>::I   = SymmetricTensor<2, dim> (unit_symmetric_tensor <dim> ());
template <int dim> SymmetricTensor<4, dim> const Material_NH<dim>::IxI = SymmetricTensor<4, dim> (outer_product (I, I));
template <int dim> SymmetricTensor<4, dim> const Material_NH<dim>::II  = SymmetricTensor<4, dim> (identity_tensor <dim> ());

// @sect3{Quadrature point history}
// As seen in step-18, the point history class offers
// a method of storing data defined at the quadrature points.
// As this method requires the nonlinear stress and
// material tangents to be evaluated at these points,
// we used this class to perform these operations.
//
// We introduce the multiplicative decomposition of the
// deformation gradient into a volume-preserving and volume
// changing component:
// \f$ \mathbf{F} = \hat{\mathbf{F}} \bar{\mathbf{F}} \f$
// where the volumetric part is
// \f$ \hat{\mathbf{F}} = J^{\frac{1}{3}} \mathbf{I} \f$
// and the isochoric part is given by
// \f$ \bar{\mathbf{F}} = J^{-\frac{1}{3}} \mathbf{F} \f$
// . From this, the deviatoric left Cauchy-Green deformation
// tensor can be defined as
// \f$ \bar{\mathbf{B}} = \bar{\mathbf{F}} \bar{\mathbf{F}}^{T} = J^{-\frac{2}{3}} \mathbf{F} \mathbf{F}^{T} \f$
//
// Here we also introduce an additive volumetric-deviatoric split
// in the material reponse. We can express the governing SEF as
// \f$ \phi = \phi_{V} + \phi_{I} \f$
// with the result that the Kirchhoff stress is additively
// decomposed into
// \f$ \mathbf{\tau} = \mathbf{\tau}_{V} + \mathbf{\tau}_{I} \f$
// as is the tangent matrix
// \f$ J\mathbf{C} = J\mathbf{C}_{V} + J\mathbf{C}_{I} \f$.
//
// These quantities are calculated as
// \f$  \mathbf{\tau}_{I} = pJ\mathbf{I} \f$
// \f$  \mathbf{\tau}_{V} = \mathcal{P}:\bar{\mathbf{\tau}} \f$
// with \f$ \bar{\mathbf{\tau}} = \mathbf{\tau} \vert_{\mathbf{B} = \bar{\mathbf{B}}} \f$
// and the deviatoric tensor \f$ \mathcal{P} = \mathcal{I} - \mathbf{I} \otimes \mathbf{I} \f$
// \f$  J\mathbf{C}_{I} = pJ(\mathbf{I} \otimes \mathbf{I} - 2 \mathcal{I}) \f$
// \f$  J\mathbf{C}_{V} = \frac{2}{3} tr\left(\bar{\mathbf{\tau}}\right) \mathcal{P} - \frac{2}{3} \left(\mathbf{\tau}_{I}\otimes\mathbf{I} + \mathbf{I}\otimes\mathbf{\tau}_{I} \right) +  \mathcal{P}:\bar{\mathcal{C}}:\mathcal{P} \f$
// with \f$ \bar{\mathcal{C}} = \mathcal{C} \vert_{\mathbf{B} = \bar{\mathbf{B}}} \f$
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

    // We first create a material object based on the data sent in.
    // This object could potentially be shared amoung QPH objects
    // but this could cause data-race issues when assembling the system matrix.
    void setup_lqp ( Parameters::AllParameters & parameters )
    {
	const double lambda = 2.0*parameters.mu*parameters.nu / (1.0-2.0*parameters.nu);
	material = new Material_NH<dim> (lambda,
					 parameters.mu);

        // Initialise all tensors correctly
        update_values (Tensor <2,dim> (),
                       0.0,
                       1.0);
    }

    // We can update the stored values and stresses based on the current
    // deformation configuration and pressure and dilation field values
    void update_values (const Tensor<2, dim> & grad_u_n,
			const double & pressure,
			const double & dilatation)
    {
        // Deformation variables calculated from displacement, displacement gradients
        static const Tensor < 2, dim> I =  static_cast <Tensor < 2, dim> > (unit_symmetric_tensor <dim> ());
        const Tensor <2,dim>  F = I + grad_u_n;
	J     = determinant(F);
	F_inv = invert(F);
	B_bar = std::pow(get_J(), -2.0/3.0) * symmetrize ( F* transpose (F) );

	// Store the precalculated pressure and dilatation
	pressure_n = pressure;
	dilatation_n = dilatation;

        // Now that all the necessary variables are set, we can update the stress tensors
        // Stress update can only update the stresses once the
        // dilatation has been set as p = p(d).
        // Note that T_iso depends on T_bar so it must be calculated afterwards.
        T_bar = material->get_T (get_J(), get_B_bar());
        T_iso = dev_P*get_T_bar();
        T_vol =-get_pressure()*get_J()*I;
    }

    // We offer and interface to retrieve certain data.
    // Here are the displacement and strain variables
    const double & get_dilatation(void) const {return dilatation_n;}
    const double & get_J (void) const {return J;}
    const Tensor <2,dim> & get_F_inv (void) const {return F_inv;}

    //, the volumetric SEF quantities
    double get_dU_dtheta (void) { return material->get_dU_dtheta(get_dilatation()); }
    double get_d2U_dtheta2 (void) { return material->get_d2U_dtheta2(get_dilatation()); }

    // and stress-based variables. These are used in the material and global
    // tangent matrix and residual assembly operations so we compute these and
    // store them.
    double get_pressure(void) {return pressure_n;}
    const SymmetricTensor<2, dim> & get_T_iso (void) const {return T_iso;}
    const SymmetricTensor<2, dim> & get_T_vol (void) const {return T_vol;}

    // Here we provide the local material tangent matrix contribution.
    // Since they are only used in the tangent matrix assembly process
    // we compute them as required.
    // This is the isochoric contribution
    SymmetricTensor <4,dim> get_C_iso(void)
    {
        const double & J = get_J();
        const SymmetricTensor<2, dim> & B_bar = get_B_bar();
        const SymmetricTensor<2, dim> & T_iso = get_T_iso();

        const SymmetricTensor <4,dim> T_iso_x_I = outer_product(T_iso, I);
        const SymmetricTensor <4,dim> I_x_T_iso = outer_product(I, T_iso);
        const SymmetricTensor <4,dim> C_bar = material->get_JC (J, B_bar);

	return     2.0/3.0*trace(get_T_bar())*dev_P
		-  2.0/3.0*(T_iso_x_I + I_x_T_iso)
		+  dev_P*C_bar*dev_P;
    }
    // and the volumetric contribution
    SymmetricTensor <4,dim> get_C_vol(void)
    {
        const double & p = get_pressure();
	const double & J = get_J();
	return -p*J*(IxI - 2.0*II);
    }

private:
    // We specify that each QP has a copy of a material
    // type in case different materials are used
    // in different regions of the domain. This also
    // deals with the issue of preventing data-races during
    // multi-threading operations when using shared objects.
    Material_NH <dim>* material;

    // These are all the volume, displacement and strain variables
    double                  dilatation_n;
    double                  J;
    Tensor <2,dim>	    F_inv;
    SymmetricTensor <2,dim> B_bar;
    SymmetricTensor <2,dim> E;
    const SymmetricTensor <2,dim> & get_B_bar (void) const {return B_bar;}

    // and the stress-type variables
    double                  pressure_n;
    SymmetricTensor<2, dim> T_bar;
    SymmetricTensor<2, dim> T_iso;
    SymmetricTensor<2, dim> T_vol;
    const SymmetricTensor<2, dim> & get_T_bar (void) const {return T_bar;}

    // Some higher-order tensors are frequently used but
    // remain unchanged. We calculate these once-off
    // and store them such that they are shared between
    // all QPH objects.
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

    // Threaded building-blocks data structures
    struct PerTaskData_K;
    struct ScratchData_K;
    struct PerTaskData_F;
    struct ScratchData_F;
    struct PerTaskData_SC;
    struct ScratchData_SC;
    struct PerTaskData_UQPH;
    struct ScratchData_UQPH;

    // Build the grid
    void make_grid (void);

    // Setup the Finite Element system to be solved
    void system_setup (void);
    void determine_component_extractors(void);

    // Assemble the system and right hand side matrices using multi-threading
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

    // Create and update the quadrature points stress and strain values
    void setup_qph(void);
    void update_qph_incremental ( const BlockVector <double> & solution_delta );
    void update_qph_incremental_one_cell (const typename DoFHandler<dim>::active_cell_iterator & cell,
					  ScratchData_UQPH & scratch,
					  PerTaskData_UQPH & data);
    void copy_local_to_global_UQPH    (const PerTaskData_UQPH & data) {}

    // Solve for the displacement using a Newton-Rhapson method
    void solve_nonlinear_timestep (BlockVector <double> & solution_delta);
    std::pair <unsigned int, double> solve_linear_system (BlockVector <double> & newton_update);

    // Solution retrieval
    BlockVector <double> get_solution_total (const BlockVector <double> & solution_delta);

    // Postprocessing and writing data to file
    void output_results(void);

    // A collection of the parameters used to describe the problem setup
    Parameters::AllParameters parameters;

    // Description of the geometry on which the problem is solved
    Triangulation<dim> triangulation;

    // Keep track of the current time and the time spent evaluating certain functions
    Time time;
    TimerOutput timer;

    // A storage object for quadrature point information
    std::vector< PointHistory <dim> > quadrature_point_history;

    // A desciption of the finite-element system including the displacement polynomial degree,
    // the degree-of-freedom handler, number of dof's per cell and the extractor objects used
    // to retrieve information from the solution vectors
    const unsigned int  degree;
    const FESystem<dim> fe;
    DoFHandler<dim>     dof_handler_ref;
    unsigned int dofs_per_cell;
    const FEValuesExtractors::Vector u_fe;
    const FEValuesExtractors::Scalar p_fe;
    const FEValuesExtractors::Scalar t_fe;

    // Description of how the block-system is arranged
    // There are 3 blocks, the first contains a vector DOF
    // while the other two describe scalar DOFs.
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

    // Rules for gauss-quadrature on both the cell and faces. The
    // number of quadrature points on both cells and faces is
    // recorded.
    QGauss<dim> qf_cell;
    QGauss<dim-1> qf_face;
    unsigned int n_q_points;
    unsigned int n_q_points_f;

    // Objects that store the converged solution and residual vectors,
    // as well as the tangent matrix. There is a ConstraintMatrix object
    // used to keep track of constraints for the nonlinear problem.
    ConstraintMatrix constraints;
    BlockSparsityPattern sparsity_pattern;
    BlockSparseMatrix <double> tangent_matrix;
    BlockVector <double> residual;
    BlockVector <double> solution_n;

    // Then define a number of variables to store residual and update
    // norms and normalisation factors.
    struct Errors
    {
        Errors (void) : norm(1.0), u (1.0), p(1.0), t(1.0) {}
        double norm,u,p,t;
        void reset (void) {norm = 1.0; u = 1.0; p = 1.0; t = 1.0;}
        void normalise (const Errors & rhs)
        {
            if (rhs.norm != 0.0) norm /= rhs.norm;
            if (rhs.u != 0.0) u /= rhs.u;
            if (rhs.p != 0.0) p /= rhs.p;
            if (rhs.t != 0.0) t /= rhs.t;
        }
    }
    error_residual, error_residual_0, error_residual_norm,
    error_update, error_update_0, error_update_norm;

    // Methods to calculate error measures
    void get_error_residual (Errors & error_residual);
    void get_error_update (const BlockVector <double> & newton_update,
                           Errors & error_update);
    double get_error_dil (void);

    // Print information to screen
    void print_conv_header (void);
    void print_conv_footer (void);
};

// @sect3{Implementation of the <code>Solid</code> class}

// @sect4{Public interface}
// We initialise the the solid class using data extracted
// from the parameter file.
template <int dim>
Solid<dim>::Solid (const std::string & input_file)
    :
      parameters (input_file),
      triangulation (Triangulation<dim>::maximum_smoothing),
      time (parameters.end_time,
            parameters.delta_t),
      timer (std::cout,
          TimerOutput::summary,
          TimerOutput::wall_times),
      degree (parameters.poly_degree),
      // The Finite Element System is composed of dim continuous
      // displacment DOFs and linear discontinuous pressure and
      // dilatation DOFs. In an attempt to satisfy the LBB conditions,
      // we setup a Q(n)-P(n-1)-P(n-1) system. Q2-P1 element satisfy
      // this condition, while Q1-P0 elements do not. However, it
      // has been shown that they demonstrate good convergence
      // characteristics nonetheless.
      fe (FE_Q<dim>(parameters.poly_degree), dim,
          FE_DGPMonomial<dim>(parameters.poly_degree-1), 1,
          FE_DGPMonomial<dim>(parameters.poly_degree-1), 1),
      dof_handler_ref (triangulation),
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

// The class destructor simply needs to clear the data held by the DOFHandler
template <int dim>
Solid<dim>::~Solid (void)
{
    dof_handler_ref.clear ();
}

// In solving the quasti-static problem, the time
// becomes a loading parameter. We choose to increment
// time linearly using a constant timestep size.
template <int dim>
void Solid<dim>::run (void)
{
    // After preprocessing, we output the initial grid
    // before starting the simulation proper.
    make_grid ();
    system_setup ();
    output_results ();
    time.increment();

    BlockVector <double> solution_delta (dofs_per_block);
    solution_delta.collect_sizes ();

    while (time.current() < time.end()) {
        // We need to reset the solution update
        // for this timestep
	solution_delta = 0.0;

	// Solve the current timestep and update total
	// solution vector
	solve_nonlinear_timestep (solution_delta);
	solution_n += solution_delta;
	output_results ();

        time.increment();
    }
}

// @sect3{Private interface}

// @sect4{Threaded-building-blocks structures}
// We choose to use TBB to perform as many computationally intensive
// distributed tasks as possible. In particular, we assemble the
// tangent matrix and residual vector, assemble the static
// condensation contributions and update data stored
// at the quadrature points.

// Firstly we deal with the tangent matrix assembly structures.
// The PerTaskData object stores local contributions.
template <int dim>
struct Solid<dim>::PerTaskData_K
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
// while the ScratchData object stores the larger objects
// such as the shape-function values object and a shape function
// values and gradient vector which we will precompute later.
template <int dim>
struct Solid<dim>::ScratchData_K
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

// Next are the same data structures used for the residual assembly.
// The PerTaskData object again stores local contributions
template <int dim>
struct Solid<dim>::PerTaskData_F
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
// and the ScratchData object the shape function object
// and precomputed values vector
template <int dim>
struct Solid<dim>::ScratchData_F
{
    FEValues <dim>     fe_values_ref;
    FEFaceValues <dim> fe_face_values_ref;

    std::vector < std::vector< double > > Nx;
    std::vector < std::vector< SymmetricTensor<2, dim> > > symm_grad_Nx;

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
              std::vector< SymmetricTensor<2, dim> >(fe_cell.dofs_per_cell))
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
          symm_grad_Nx (rhs.symm_grad_Nx)
    {  }

    void reset (void) {
        for (unsigned int q_point=0; q_point < symm_grad_Nx.size(); ++q_point) {
            for (unsigned int k=0; k < symm_grad_Nx[q_point].size(); ++k) {
                Nx[q_point][k] = 0.0;
                symm_grad_Nx[q_point][k] = 0.0;
            }
        }
    }

};

// Here we define structures to assemble the static condensation contributions.
// As the operations are matrix-based, we need to setup a number of matrices
// to store the local contributions from a number of the tangent matrix subblocks.
// We place these in the PerTaskData struct.
template <int dim>
struct Solid<dim>::PerTaskData_SC
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
          cell_matrix (dofs_per_cell,
                       dofs_per_cell),
          local_dof_indices  (dofs_per_cell),
          K_pt_inv (n_t, n_p),
          K_tt_inv (n_t, n_t),
          K_con (n_u, n_u),
          A (n_t, n_u),
          B (n_t, n_u),
          C (n_p, n_u)
    {  }

    // Choose not to reset any data as the matrix extraction and
    // replacement tools will take care of this
    void reset(void) { }
};
// The ScratchData object is not strictly necessary for the
// operations we wish to perform, but it still needs to be defined for the
// current implementation of TBB in deal.II.So we creatre a dummy struct for this purpose.
template <int dim>
struct Solid<dim>::ScratchData_SC
{
    ScratchData_SC (void) { }
    ScratchData_SC (const ScratchData_SC & rhs) { }
    void reset (void) { }
};

// And finally we define the structures to assist with updating the quadrature
// point information. Similar to the SC assembly process, we choose not to use
// the PerTaskData object to store any information but must define one nonetheless.
template <int dim>
struct Solid<dim>::PerTaskData_UQPH
{
    PerTaskData_UQPH (void) { }
    void reset(void) { }
};
// The ScratchData object will be used to store a alias fort the solution vector
// so that we don't have to copy this large data structure. We then define
// a number of vectors to extract the solution values and gradients at the
// quadrature points.
template <int dim>
struct Solid<dim>::ScratchData_UQPH
{
    const BlockVector <double> & solution_total;

    std::vector< Tensor< 2, dim> > solution_grads_u_total;
    std::vector <double> solution_values_p_total;
    std::vector <double> solution_values_t_total;

    FEValues<dim> fe_values_ref;

    ScratchData_UQPH (const FiniteElement <dim> & fe_cell,
                      const QGauss <dim> & qf_cell,
                      const UpdateFlags uf_cell,
                      const BlockVector <double> & solution_total)
        :
          solution_total (solution_total),
          solution_grads_u_total (qf_cell.size()),
          solution_values_p_total (qf_cell.size()),
          solution_values_t_total (qf_cell.size()),
          fe_values_ref (fe_cell,
              qf_cell,
              uf_cell)
    { }

    ScratchData_UQPH (const ScratchData_UQPH & rhs)
        :
          solution_total (rhs.solution_total),
          solution_grads_u_total (rhs.solution_grads_u_total),
          solution_values_p_total (rhs.solution_values_p_total),
          solution_values_t_total (rhs.solution_values_t_total),
          fe_values_ref (rhs.fe_values_ref.get_fe(),
                            rhs.fe_values_ref.get_quadrature(),
                            rhs.fe_values_ref.get_update_flags())
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

// @sect4{Solid::make_grid}
// Here we create the grid on which the minimisation problem is to be solved.
template <int dim>
void Solid<dim>::make_grid (void)
{
    // Create a unit cube with each face given a boundary ID number
    GridGenerator::hyper_rectangle ( triangulation,
				     Point<dim> (0.0, 0.0, 0.0),
				     Point<dim> (1.0, 1.0, 1.0),
				     true );
    GridTools::scale (parameters.scale,
                      triangulation);

    // The grid must be refined at least once for the indentation problem
    if (parameters.global_refinement == 0)
        triangulation.refine_global (1);
    else
        triangulation.refine_global (parameters.global_refinement);

    // Since we wish to apply a Neumann BC to a patch on the top surface,
    // we must find the cell faces in this part of the domain and
    // mark them with a distinct boundary ID number
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
// Next we describe how the FE system is setup.
template <int dim>
void Solid<dim>::system_setup (void)
{
    timer.enter_subsection ("Setup system");

    // We first describe the number of components per block. Since the
    // displacement is a vector component, the first dim components
    // belong to it, while the next two describe scalar pressure and
    // dilatation DOFs.
    std::vector<unsigned int> block_component (n_components, u_dof); // Displacement
    block_component[p_component] = p_dof; // Pressure
    block_component[t_component] = t_dof; // Dilatation

    // DOF handler is then initialised and we renumber the grid in an
    // efficient manner. We also record the number of DOF's per block.
    dof_handler_ref.distribute_dofs (fe);
    DoFRenumbering::Cuthill_McKee (dof_handler_ref);
    DoFRenumbering::component_wise (dof_handler_ref,
                                    block_component);
    DoFTools::count_dofs_per_block (dof_handler_ref,
                                    dofs_per_block,
                                    block_component);

    std::cout
	    << "Triangulation:"
	    << "\n\t Number of active cells: " << triangulation.n_active_cells()
	    << "\n\t Number of degrees of freedom: " << dof_handler_ref.n_dofs()
	    << std::endl;

    // Setup the sparsity pattern and tangent matrix
    tangent_matrix.clear ();
    {
	const unsigned int n_dofs_u = dofs_per_block[u_dof];
	const unsigned int n_dofs_p = dofs_per_block[p_dof];
	const unsigned int n_dofs_t = dofs_per_block[t_dof];

        BlockCompressedSimpleSparsityPattern csp (n_blocks,
                                                  n_blocks);

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

        // The global system matrix will have the following structure
        //      | K'_uu |  K_up  |     0     |         | dU_u |         | dR_u |
        // K =  | K_pu  |    0   |   K_pt^-1 | , dU =  | dU_p | , dR =  | dR_p |
        //      |   0   |  K_tp  |   K_tt    |         | dU_t |         | dR_t |
        // We optimise the sparsity pattern to reflect this structure
        // and prevent unnecessary data creation for the right-diagonal
        // block components.
        Table<2,DoFTools::Coupling> coupling (n_components, n_components);
        for (unsigned int ii = 0; ii < n_components; ++ii) {
            for (unsigned int jj = 0; jj < n_components; ++jj) {

                if (    ( (ii <  p_component) && (jj == t_component) )
                     || ( (ii == t_component) && (jj <  p_component) )
                     || ( (ii == p_component) && (jj == p_component) ) )
                {
                    coupling[ii][jj] = DoFTools::none;
                }
                else {
                    coupling[ii][jj] = DoFTools::always;
                }
            }
        }
        DoFTools::make_sparsity_pattern (dof_handler_ref, coupling, csp, constraints, false);
        sparsity_pattern.copy_from (csp);
    }
    
    tangent_matrix.reinit (sparsity_pattern);

    // Setup storage vectors noting that the dilatation is unity
    // in the reference configuration
    residual.reinit (dofs_per_block);
    residual.collect_sizes ();

    solution_n.reinit (dofs_per_block);
    solution_n.collect_sizes ();
    solution_n.block(t_dof) = 1.0;

    // and finally set up the quadrature point history
    setup_qph ();

    timer.leave_subsection();
}

// We next get information from the FE system
// that describes which local element DOFs are
// attached to which block component.
// This is used later to extract subblocks from the global matrix.
template <int dim>
void Solid<dim>::determine_component_extractors(void)
{
    element_indices_u.clear();
    element_indices_p.clear();
    element_indices_t.clear();

    for (unsigned int k=0; k < fe.dofs_per_cell; ++k) {
        // The next call has the FE System indicate to which block component
        // the current DOF is attached to.
        // Currently, the interpotation fields are setup such that
        // 0 indicates a displacement DOF, 1 a pressure DOF and 2 a dilatation DOF.
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
// The method used to store quadrature information is already described in
// tutorial 18. Here we implement a similar setup for a SMP machine.
template <int dim>
void Solid<dim>::setup_qph (void)
{
    std::cout << "    Setting up quadrature point data..." << std::endl;

    // Firstly the actual QPH data objects are created. This must be done
    // only once the grid is refined to its finest level.
    {
        quadrature_point_history = std::vector< PointHistory <dim> > (triangulation.n_active_cells() * n_q_points);

    	unsigned int history_index = 0;
        typename Triangulation<dim>::active_cell_iterator
                        cell = triangulation.begin_active(),
                        endc = triangulation.end();
    	for (cell = triangulation.begin_active(); cell != endc; ++cell) {
    	    cell->set_user_pointer(&quadrature_point_history[history_index]);
    	    history_index += n_q_points;
    	}

    	Assert(history_index == quadrature_point_history.size(), ExcInternalError());
    }

    // Next we setup the initial QP data
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
// As the update of QP information occurs frequently and involves a number of
// expensive operations, we define a multi-threaded approach to distributing
// the task across a number of CPU cores.
template <int dim>
void Solid<dim>::update_qph_incremental (const BlockVector <double> & solution_delta)
{
    timer.enter_subsection("Update QPH data");
    std::cout << " UQPH "<< std::flush;

    // Firstly we need to attain the total solution as it stands
    // at this Newton increment
    const BlockVector <double> solution_total = get_solution_total(solution_delta);

    // Next we create the initial copy of TBB objects
    const UpdateFlags uf_UQPH ( update_values | update_gradients );
    PerTaskData_UQPH per_task_data_UQPH;
    ScratchData_UQPH scratch_data_UQPH (fe,
					qf_cell,
					uf_UQPH,
					solution_total);

    // and pass them and the one-cell update function to the workstream to be processed
    WorkStream::run (  dof_handler_ref.begin_active(),
		       dof_handler_ref.end(),
		       *this,
		       &Solid::update_qph_incremental_one_cell,
		       &Solid::copy_local_to_global_UQPH,
		       scratch_data_UQPH,
		       per_task_data_UQPH);

    timer.leave_subsection();
}

// Now we describe how we extract data from the solution vector and pass it
// along to each QP storage object for processing.
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

    // Firstly we need to find the values and gradients at quadrature points
    // inside the current cell
    scratch.fe_values_ref.reinit(cell);
    scratch.fe_values_ref[u_fe].get_function_gradients (scratch.solution_total, scratch.solution_grads_u_total);
    scratch.fe_values_ref[p_fe].get_function_values (scratch.solution_total, scratch.solution_values_p_total);
    scratch.fe_values_ref[t_fe].get_function_values (scratch.solution_total,scratch. solution_values_t_total);

    // and then we update the each local QP using the displacment deformation gradient
    // and total pressure and dilatation solution values.
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
            << std::endl
	    << "Timestep " << time.get_timestep()
	    << " @ " << time.current() << "s"
	    << std::endl;

    // We create a new vector to store the current Newton update step
    BlockVector <double> newton_update (dofs_per_block);
    newton_update.collect_sizes ();

    // Reset the error storage objects
    error_residual.reset();
    error_residual_0.reset();
    error_residual_norm.reset();
    error_update.reset();
    error_update_0.reset();
    error_update_norm.reset();

    // Print solver header
    print_conv_header();

    // We now perform a number of Newton iterations to iteratively solve
    // the nonlinear problem.
    for (unsigned int it_nr=0; it_nr < parameters.max_iterations_NR; ++ it_nr)
    {
        // Print Newton iteration
	std::cout
		<< " "
		<< std::setw(2)
		<< it_nr
		<< " "
		<< std::flush;

	// Since the problem is fully nonlinear and we are using a
	// full Newton method, the data stored in the tangent matrix
	// and residual vector is not reusable and must be cleared
	// at each Newton step.
	tangent_matrix = 0.0;
	residual = 0.0;

	// We initially build the residual vector to check for convergence.
	// The unconstrained DOF's of the residual vector hold the out-of-balance
	// forces. This is done before assembling the system matrix as the latter
	// is an expensive operation and we can potentially avoid an extra
	// assembly process by not assembling the tangent matrix when convergence
	// is attained.
	assemble_system_F (); // Assemble RHS
	get_error_residual(error_residual);

	// We store the residual errors after the first iteration
	// in order to normalise by their value
	if (it_nr == 0) error_residual_0 = error_residual;

	// We can now determine the normalised residual error
	error_residual_norm = error_residual;
	error_residual_norm.normalise(error_residual_0);

	// Check for solution convergence
	if (   it_nr > 0
	       && error_update_norm.u <= parameters.tol_u
	       && error_residual_norm.u <= parameters.tol_f)
	{
	    std::cout
		    << " CONVERGED! "
		    << std::endl;

	    print_conv_footer();

	    //	    timer.leave_subsection();
	    return;
	}


	assemble_system_K (); // Assemble stiffness matrix
	make_constraints (it_nr, constraints); // Make boundary conditions
	constraints.condense (tangent_matrix,
			      residual); // Apply BC's

	const std::pair <unsigned int, double> lin_solver_output = solve_linear_system (newton_update);
	constraints.distribute(newton_update); // Populate the constrained DOF's with their values

	get_error_update(newton_update,
			 error_update);
	if (it_nr == 0) error_update_0 = error_update;
	// We can now determine the normalised newton update error
	error_update_norm = error_update;
	error_update_norm.normalise(error_update_0);

	// The current solution state unacceptable, so we need to update
	// the solution increment for this timestep, update all quadrature
	// point inforation pertaining to this new displacment and stress state
	// and continue iterating.
	solution_delta += newton_update;
	update_qph_incremental (solution_delta);

	std::cout
		<< " | "
		<< std::fixed
		<< std::setprecision(3)
		<< std::setw(7)
		<< std::scientific
		<< lin_solver_output.first << "  "
		<< lin_solver_output.second << "  "
		<< error_residual_norm.norm << "  "
		<< error_residual_norm.u << "  "
		<< error_residual_norm.p << "  "
		<< error_residual_norm.t << "  "
		<< error_update_norm.norm << "  "
		<< error_update_norm.u << "  "
		<< error_update_norm.p << "  "
		<< error_update_norm.t << "  "
		<< std::endl;
    }

    throw(ExcMessage("No convergence in nonlinear solver!"));
}

// We print out data in a nice table that is updated
// on a per-iteration basis. Here we set up the table
// header
template <int dim>
void Solid<dim>::print_conv_header (void)
{
    static const unsigned int l_width = 155;

    for (unsigned int i=0; i < l_width; ++i)
        std::cout << "_";
    std::cout << std::endl;

    std::cout
            << "                 "
            << "SOLVER STEP"
            << "                  "
            << " | "
            << " LIN_IT  "
            << " LIN_RES   "
            << " RES_NORM    "
            << " RES_U    "
            << " RES_P     "
            << " RES_T    "
            << " NU_NORM     "
            << " NU_U      "
            << " NU_P      "
            << " NU_T "
            << std::endl;

    for (unsigned int i=0; i < l_width; ++i)
        std::cout << "_";
    std::cout << std::endl;
}
// and here the footer
template <int dim>
void Solid<dim>::print_conv_footer (void)
{
    static const unsigned int l_width = 155;

    for (unsigned int i=0; i < l_width; ++i)
        std::cout << "_";
    std::cout << std::endl;


    std::cout
            << "Relative errors:" << std::endl
            << "Displacement:\t" << error_update.u/error_update_0.u << std::endl
            << "Force: \t\t" << error_residual.u/error_residual_0.u << std::endl
            << "Dilatation:\t" << get_error_dil()
            << std::endl;
}

// Calculate the ratio of the volume of the domain in the
// current configuration and the reference configuration
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

// Determine the true residual error for the problem
template <int dim>
void Solid<dim>::get_error_residual (Errors & error_residual)
{
    BlockVector <double> error_res (dofs_per_block);
    error_res.collect_sizes ();

    // Need to ignore constrained DOFs
    for (unsigned int i=0; i < dof_handler_ref.n_dofs(); ++i)
        if (!constraints.is_constrained(i))
            error_res(i) = residual(i);

    error_residual.norm = error_res.l2_norm();
    error_residual.u = error_res.block(u_dof).l2_norm();
    error_residual.p = error_res.block(p_dof).l2_norm();
    error_residual.t = error_res.block(t_dof).l2_norm();
}

// Determine the true Newton update error for the problem
template <int dim>
void Solid<dim>::get_error_update (const BlockVector <double> & newton_update,
                                   Errors & error_update)
{
    BlockVector <double> error_ud (dofs_per_block);
    error_ud.collect_sizes ();

    // Need to ignore constrained DOFs as they have a prescribed
    // value
    for (unsigned int i=0; i < dof_handler_ref.n_dofs(); ++i)
        if (!constraints.is_constrained(i))
            error_ud(i) = newton_update(i);

    error_update.norm = error_ud.l2_norm();
    error_update.u = error_ud.block(u_dof).l2_norm();
    error_update.p = error_ud.block(p_dof).l2_norm();
    error_update.t = error_ud.block(t_dof).l2_norm();
}

// This function provides the total solution, which is valid at any Newton step.
// This is required as, to reduce computational error, the total solution is
// only updated at the end of the timestep.
template <int dim>
BlockVector <double> Solid<dim>::get_solution_total (const BlockVector <double> & solution_delta)
{
    BlockVector <double> solution_total (solution_n);
    solution_total += solution_delta;

    return solution_total;
}

// @sect4{Solid::assemble_system_K}
// Since we use TBB for assembly, we simply setup a copy of the
// data structures required for the process and pass them, along
// with the memory addresses of the assembly functions to the
// WorkStream object for processing. Note that we must ensure that
// the matrix is reset before any assembly operations can occur.
template <int dim>
void Solid<dim>::assemble_system_K (void)
{
    timer.enter_subsection("Assemble tangent matrix");
    std::cout << " ASM_K " << std::flush;

    tangent_matrix = 0.0;

    const UpdateFlags uf_cell (update_values | update_gradients | update_JxW_values);

    PerTaskData_K per_task_data (dofs_per_cell);
    ScratchData_K scratch_data (fe,
                                qf_cell,
                                uf_cell);

    WorkStream::run (  dof_handler_ref.begin_active(),
                       dof_handler_ref.end(),
                       *this,
                       &Solid::assemble_system_K_one_cell,
                       &Solid::copy_local_to_global_K,
                       scratch_data,
                       per_task_data);

    timer.leave_subsection();
}

// This function adds the local contribution to the system matrix.
// Note that we choose not to use the constraint matrix to do the
// job for us because the tangent matrix and residual processes have
// been split up into two seperate functions.
template <int dim>
void Solid<dim>::copy_local_to_global_K (const PerTaskData_K & data)
{
    for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
            tangent_matrix.add (data.local_dof_indices[i],
                               data.local_dof_indices[j],
                               data.cell_matrix(i,j));
}

// Here we define how we assemble the tangent matrix contribution for a
// single cell.
template <int dim>
void Solid<dim>::assemble_system_K_one_cell (const typename DoFHandler<dim>::active_cell_iterator & cell,
                                             ScratchData_K & scratch,
                                             PerTaskData_K & data)
{
    // We first need to reset and initialise some of the data structures and retrieve some
    // basic information regarding the DOF numbering on this cell
    data.reset();
    scratch.reset();
    scratch.fe_values_ref.reinit (cell);
    cell->get_dof_indices (data.local_dof_indices);
    PointHistory<dim> *lqph = reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

    // We can precalculate the cell shape function values and gradients. Note that the
    // shape function gradients are defined in the current configuration for this problem.
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

    // Now we build the local cell stiffness matrix. Since the global and local system
    // matrices are symmetric, we can exploit this property by building only the lower
    // half of the local matrix and copying those values to the upper half.
    // So we only assemble half of the K_uu, K_pp (= 0), K_tt blocks, while the whole
    // K_pt, K_ut, K_up blocks are built.
    for (unsigned int q_point=0; q_point < n_q_points; ++q_point) {
        // We first extract some configuration dependent variables from our
        // QPH history objects that remain constant at each QP.
        const Tensor <2,dim>          T   = static_cast < Tensor<2, dim> > (lqph[q_point].get_T_iso() + lqph[q_point].get_T_vol());
        const SymmetricTensor <4,dim> C   = lqph[q_point].get_C_iso() + lqph[q_point].get_C_vol();
        const double                  C_v = lqph[q_point].get_d2U_dtheta2();
        const double                  J   = lqph[q_point].get_J();

	// Next we define some aliases to make the assembly process easier to follow
	const std::vector<double> & N = scratch.Nx[q_point];
	const std::vector< SymmetricTensor <2,dim> > & symm_B = scratch.symm_grad_Nx[q_point];
	const std::vector< Tensor <2,dim> > & B = scratch.grad_Nx[q_point];
	const double & JxW = scratch.fe_values_ref.JxW(q_point);

	for (unsigned int i=0; i < dofs_per_cell; ++i) {
	    const unsigned int component_i = fe.system_to_component_index(i).first;
	    // Determine the dimensional component that matches the dof component (i.e. i % dim)
	    const unsigned int i_group = fe.system_to_base_index(i).first.first;

	    for (unsigned int j=0; j <= i; ++j) {
		const unsigned int component_j = fe.system_to_component_index(j).first;
		const unsigned int j_group = fe.system_to_base_index(j).first.first;

		// This is the K_{uu} contribution. It comprises of a material stiffness
		// contribution and a geometric stiffness contribution which is only
		// added along the local matrix diagonals
		if (   (i_group == j_group) && (i_group == u_dof ) ) {
		    data.cell_matrix(i,j) += symm_B[i] * C * symm_B[j] * JxW;
		    if (component_i == component_j)
			data.cell_matrix(i,j) += B[i][component_i] * T * B[j][component_j] * JxW;
		}
		// Next is the K_{pu} contibution
		else if ( (i_group == p_dof) && (j_group == u_dof) ) {
		    data.cell_matrix(i,j) -= N[i]*J*(symm_B[j]*I)*JxW;
		}
		// and the K_{tp} contibution
		else if ( (i_group == t_dof) && (j_group == p_dof) ) {
		    data.cell_matrix(i,j) += N[i]*N[j]*JxW;
		}
		// and lastly the K_{tt} contibution
		else if ( (i_group == j_group) && (i_group == t_dof)  ) {
		    data.cell_matrix(i,j) -=  N[i]*C_v*N[j]*JxW;
		}
		else Assert ((i_group <= t_dof) && (j_group <= t_dof), ExcInternalError());
	    }
	}
    }

    // Here we copy the lower half of the local matrix in the upper
    // half of the local matrix
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
        for (unsigned int j=i+1; j<dofs_per_cell; ++j) {
            data.cell_matrix(i,j) = data.cell_matrix(j,i);
        }
    }
}

// @sect4{Solid::assemble_system_F}
// The setup of the residual assembly process is similar to the
// tangent matrix, so we will not describe it in too much detail.
// Note that since we are describing a problem with Neumann BCs,
// we will need the face normals and so must specify this in the
// update flags.
template <int dim>
void Solid<dim>::assemble_system_F (void)
{
    timer.enter_subsection("Assemble residual");
    std::cout << " ASM_R "<< std::flush;

    residual  = 0.0;

    const UpdateFlags uf_cell (update_values | update_gradients | update_JxW_values);
    const UpdateFlags uf_face (update_values | update_normal_vectors | update_JxW_values);

    PerTaskData_F per_task_data (dofs_per_cell);
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
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
        residual(data.local_dof_indices[i]) += data.cell_rhs(i);
    }
}

template <int dim>
void Solid<dim>::assemble_system_F_one_cell (const typename DoFHandler<dim>::active_cell_iterator & cell,
                                             ScratchData_F & scratch,
                                             PerTaskData_F & data)
{
    // Again we reset the data structures
    data.reset();
    scratch.reset();
    scratch.fe_values_ref.reinit (cell);
    cell->get_dof_indices (data.local_dof_indices);
    PointHistory<dim> *lqph = reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

    // and then precompute some shape function data
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

    // and can now assemble the residual contribution
    for (unsigned int q_point=0; q_point < n_q_points; ++q_point) {
        // We fist retrieve data that remains constant a QP
        const SymmetricTensor <2,dim>  T = lqph[q_point].get_T_iso() + lqph[q_point].get_T_vol();
        const double  J = lqph[q_point].get_J();
        const double  D = lqph[q_point].get_dilatation();
        const double  p = lqph[q_point].get_pressure();
        const double  p_star = lqph[q_point].get_dU_dtheta();

	// define some shortcuts
	const std::vector< double > & N = scratch.Nx[q_point];
	const std::vector< SymmetricTensor <2,dim> > & symm_B = scratch.symm_grad_Nx[q_point];
	const double  JxW = scratch.fe_values_ref.JxW(q_point);

	for (unsigned int i=0; i<dofs_per_cell; ++i) {
	    const unsigned int i_group = fe.system_to_base_index(i).first.first;
	    // Add the contribution to the R_{u} block
	    if (i_group == u_dof) {
		data.cell_rhs(i) -= ( symm_B[i]*T )*JxW;
	    }
	    // the R_{p} block
	    else if (i_group == p_dof ) {
		data.cell_rhs(i) += N[i]*(J - D)*JxW;
	    }
	    // and finally the R_{t} block
	    else if ( i_group == t_dof) {
		data.cell_rhs(i) += N[i]*(p_star-p)*JxW;
	    }
	    else Assert (i_group <= t_dof, ExcInternalError());
	}
    }

    // Next we assemble the Neumann contribution. We first check to see
    // it the cell face exists on a boundary on which a traction is
    // applied and add the contribution if this is the case.
    if (cell->at_boundary() == true) {
        for (unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; ++face) {
            if (    cell->face(face)->at_boundary() == true
                    &&  cell->face(face)->boundary_indicator() == 6 ) {
                scratch.fe_face_values_ref.reinit (cell, face);

		for (unsigned int f_q_point=0; f_q_point < n_q_points_f; ++f_q_point) {
		    // We retrieve the face normal at this QP
		    const Tensor <1, dim> & N = scratch.fe_face_values_ref.normal_vector(f_q_point);

		    // and specify the traction in reference configuration. For this problem,
		    // a defined pressure is applied in the reference configuration. so the
		    // traction defined using the first Piola-Kirchhoff stress is simply
		    // t_0 = P*N = (pI)*N = p*N
		    // We choose to use the time variable to linearly ramp up the pressure
		    // load.
		    static const double p0 = -4.0/(parameters.scale*parameters.scale);
		    const double time_ramp = (time.current() / time.end());
		    const double pressure = p0 * parameters.p_p0 * time_ramp;
		    const Tensor <1,dim> traction = pressure * N;

		    for (unsigned int i=0; i < dofs_per_cell; ++i) {
			const unsigned int i_group = fe.system_to_base_index(i).first.first;

			if (i_group == u_dof) {
			    // More shortcuts being assigned
			    const unsigned int component_i = fe.system_to_component_index(i).first;
			    const double & Ni = scratch.fe_face_values_ref.shape_value(i,f_q_point);
			    const double & JxW = scratch.fe_face_values_ref.JxW(f_q_point);

			    // And finally we can add the traction vector contribution to
			    // the local RHS vector. Note that this contribution is present
			    // on displacement DOFs only.
			    data.cell_rhs(i) += (Ni * traction[component_i]) * JxW;
			}
		    }
		}
	    }
	}
    }
}

// @sect4{Solid::make_constraints}
// The constraints for this problem are simple to describe.
// However, since we are dealing with an iterative Newton method,
// it should be noted that any displacement constraints should only
// be specified at the zeroth iteration and subsequently no
// additional contributions are to be made since the constraints
// are already exactly satisfied. So we describe this process for
// completeness although for this problem the constraints are
// trivial and it would not have made a difference if this had
// not been accounted for in this problem.
template <int dim>
void Solid<dim>::make_constraints (const int & it_nr,
                                   ConstraintMatrix & constraints)
{
    std::cout << " CST "<< std::flush;

    // Since the constraints are different at Newton iterations,
    // we need to clear the constraints matrix and completely
    // rebuild it. However, after the first iteration, the
    // constraints remain the same and we can simply skip the
    // rebuilding step if we do not clear it.
    if (it_nr > 1) return;
    constraints.clear();
    const bool apply_dirichlet_bc = (it_nr == 0);

    // The boundary conditions for the indentation problem are as follows:
    // On the -x, -y and -z faces (ID's 0,2,4) we set up a symmetry condition
    // to allow only planar movement while the +x and +y faces (ID's 1,3) are
    // traction free. In this contrived problem, part of the +z face (ID 5) is
    // set to have no motion in the x- and y-component. Finally, as described
    // earlier, the other part of the +z face has an the applied pressure but
    // is also constrained in the x- and y-directions.
    {
        const int boundary_id = 0;

	std::vector< bool > components (n_components, false);
	components[0] = true;

	if (apply_dirichlet_bc == true) {
	    VectorTools::interpolate_boundary_values ( dof_handler_ref,
						       boundary_id,
						       ZeroFunction<dim>(n_components),
						       constraints,
						       components );
	}
	else {
	    VectorTools::interpolate_boundary_values ( dof_handler_ref,
						       boundary_id,
						       ZeroFunction<dim>(n_components),
						       constraints,
						       components );
	}
    }
    {
        const int boundary_id = 2;

	std::vector< bool > components (n_components, false);
	components[1] = true;

	if (apply_dirichlet_bc == true) {
	    VectorTools::interpolate_boundary_values ( dof_handler_ref,
						       boundary_id,
						       ZeroFunction<dim>(n_components),
						       constraints,
						       components );
	}
	else {
	    VectorTools::interpolate_boundary_values ( dof_handler_ref,
						       boundary_id,
						       ZeroFunction<dim>(n_components),
						       constraints,
						       components );
	}
    }
    {
        const int boundary_id = 4;
        std::vector< bool > components (n_components, false);
        components[2] = true;

	if (apply_dirichlet_bc == true) {
	    VectorTools::interpolate_boundary_values ( dof_handler_ref,
						       boundary_id,
						       ZeroFunction<dim>(n_components),
						       constraints,
						       components );
	}
	else {
	    VectorTools::interpolate_boundary_values ( dof_handler_ref,
						       boundary_id,
						       ZeroFunction<dim>(n_components),
						       constraints,
						       components );
	}
    }
    {
        const int boundary_id = 5;
        std::vector< bool > components (n_components, true);
        components[2] = false;

	if (apply_dirichlet_bc == true) {
	    VectorTools::interpolate_boundary_values ( dof_handler_ref,
						       boundary_id,
						       ZeroFunction<dim>(n_components),
						       constraints,
						       components );
	}
	else {
	    VectorTools::interpolate_boundary_values ( dof_handler_ref,
						       boundary_id,
						       ZeroFunction<dim>(n_components),
						       constraints,
						       components );
	}
    }
    {
        const int boundary_id = 6;
        std::vector< bool > components (n_components, true);
        components[2] = false;

	if (apply_dirichlet_bc == true) {
	    VectorTools::interpolate_boundary_values ( dof_handler_ref,
						       boundary_id,
						       ZeroFunction<dim>(n_components),
						       constraints,
						       components );
	}
	else {
	    VectorTools::interpolate_boundary_values ( dof_handler_ref,
						       boundary_id,
						       ZeroFunction<dim>(n_components),
						       constraints,
						       components );
	}
    }

    constraints.close();
}

// @sect4{Solid::solve_linear_system}
// Solving the entire block system is a bit problematic as there are no
// contributions to the K_{pp} block, rendering it non-invertable.
// Since the pressure and dilatation variables DOFs are discontinuous, we can
// condense them out to form a smaller displacement-only system which
// we will then solve and subsequently post-process to retrieve the
// pressure and dilatation solutions.
template <int dim>
std::pair <unsigned int, double> Solid<dim>::solve_linear_system (BlockVector <double> & newton_update)
{
    // Need to create two temporary vectors so that the static condensation operation can be performed
    BlockVector <double> A (dofs_per_block);
    BlockVector <double> B (dofs_per_block);
    A.collect_sizes ();
    B.collect_sizes ();

    // Store the number of linear solver iterations and residual
    unsigned int lin_it = 0;
    double lin_res = 0.0;

    //      | K'_uu |   K_up  |     0     |         | dU_u |         | dR_u |
    // K =  | K_pu  |     0   |   K_pt^-1 | , dU =  | dU_p | , dR =  | dR_p |
    //      |   0   |   K_tp  |   K_tt    |         | dU_t |         | dR_t |

    // Solve for du
    {
        // Do the static condensation to make K'_uu,
        // and put K_pt^{-1} in the K_pt block
        assemble_SC();

	// K'uu du = Ru'
	// with Ru' = Ru − Kup Ktp^-1 (Rt − Ktt Kpt^{-1} Rp)
	// Assemble the RHS vector to solve for du
	tangent_matrix.block(p_dof, t_dof).vmult(A.block(t_dof), residual.block(p_dof));
	tangent_matrix.block(t_dof, t_dof).vmult (B.block(t_dof), A.block(t_dof));
	A.block(t_dof).equ(1.0, residual.block(t_dof), -1.0, B.block(t_dof));
	tangent_matrix.block(p_dof, t_dof).Tvmult(A.block(p_dof), A.block(t_dof));
	tangent_matrix.block(u_dof, p_dof).vmult(A.block(u_dof), A.block(p_dof));
	residual.block(u_dof) -= A.block(u_dof);

	timer.enter_subsection("Linear solver");
	std::cout << " SLV " << std::flush;
	if (parameters.type_lin == "CG")
	{
	    const int solver_its = tangent_matrix.block(u_dof, u_dof).m() * parameters.max_iterations_lin;
	    const double tol_sol = parameters.tol_lin * residual.block(u_dof).l2_norm();

	    SolverControl solver_control (solver_its , tol_sol);

	    GrowingVectorMemory < Vector<double> > GVM;
	    SolverCG < Vector<double> >  solver_CG (solver_control, GVM);

	    // We've chosen a SSOR preconditioner as it appears to provide
	    // the fastest solver convergence characteristics for this problem.
	    PreconditionSSOR <SparseMatrix<double> > preconditioner;
	    preconditioner.initialize (tangent_matrix.block(u_dof, u_dof), parameters.ssor_relaxation);

	    solver_CG.solve (tangent_matrix.block(u_dof, u_dof),
			     newton_update.block(u_dof),
			     residual.block(u_dof),
			     preconditioner);

	    lin_it = solver_control.last_step();
	    lin_res = solver_control.last_value();
	}
	else if (parameters.type_lin == "Direct")
	{
	    // Otherwise if the problem is small enough, a direct solver
	    // can be utilised.
	    SparseDirectUMFPACK  A_direct;
	    A_direct.initialize(tangent_matrix.block(u_dof, u_dof));
	    A_direct.vmult (newton_update.block(u_dof),
			    residual.block(u_dof));

	    lin_it = 1;
	    lin_res = 0.0;
	}
	else throw (ExcMessage("Linear solver type not implemented"));
	timer.leave_subsection();
    }

    timer.enter_subsection("Linear solver postprocessing");
    std::cout << " PP " << std::flush;
    // Now that we've solved the displacement problem, we can post-process
    // to get the dilatation solution from the substitution
    // dt = Kpt^{-1} ( Rp - Kpu du )
    {
        tangent_matrix.block(p_dof, u_dof).vmult (A.block(p_dof), newton_update.block(u_dof));
        A.block(p_dof) *= -1.0;
        A.block(p_dof) += residual.block(p_dof);
        tangent_matrix.block(p_dof, t_dof).Tvmult (newton_update.block(t_dof), A.block(p_dof));
    }
    // and finally we solve for the pressure update with the substitution
    // dp = Ktp^{-1} ( Rt - Ktt dt )
    {
        tangent_matrix.block(t_dof, t_dof).vmult (A.block(t_dof), newton_update.block(t_dof));
    	A.block(t_dof) *= -1.0;
        A.block(t_dof) += residual.block(t_dof);
        tangent_matrix.block(p_dof, t_dof).vmult (newton_update.block(p_dof), A.block(t_dof));
    }
    timer.leave_subsection();

    return std::make_pair(lin_it, lin_res);
}

// @sect4{Solid::assemble_system_SC}
// The static condensation process could be performed at a global level
// but we need the inverse of one of the blocks. However, since the
// pressure and dilatation variables are discontinous, the SC operation
// can be done on a per-cell basis and we can produce the inverse of the
// block-diagonal K_{pt} block by inverting the local blocks. We can
// again use TBB to do this since each operation will be independent of
// one another.
template <int dim>
void Solid<dim>::assemble_SC  (void)
{
    timer.enter_subsection("Perform static condensation");
    std::cout << " ASM_SC " << std::flush;

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

// We need to describe how to add the local contribution to the tangent matrix.
template <int dim>
void Solid<dim>::copy_local_to_global_SC (const PerTaskData_SC & data)
{
    for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
            tangent_matrix.add (data.local_dof_indices[i],
                                data.local_dof_indices[j],
                                data.cell_matrix(i,j));
}

// Now we describe the static condensation process.
template <int dim>
void Solid<dim>::assemble_SC_one_cell (const typename DoFHandler<dim>::active_cell_iterator & cell,
                                       ScratchData_SC & scratch,
                                       PerTaskData_SC & data)
{
    // As per usual, we must first find out which global numbers the
    // degrees of freedom on this cell have and reset some data structures
    data.reset();
    scratch.reset();
    cell->get_dof_indices (data.local_dof_indices);

    // Currently the the local stifness matrix K_e is of the form
    //  | K_uu  |   K_up   |   0   |
    //  | K_pu  |     0    |  K_pt |
    //  |   0   |   K_tp   |  K_tt |
    //
    // We now need to modify it such that it appear as
    //  | K'_uu |   K_up   |     0     |
    //  | K_pu  |     0    |   K_pt^-1 |
    //  |   0   |   K_tp   |   K_tt    |
    // with K'_uu = K_uu + Kup Ktp^{-1} Ktt Kpt^{-1} Kpu
    //
    // At this point, we need to take note of the fact that
    // global data already exists in the K_uu, K_pt, K_tp subblocks.
    // So if we are to modify them, we must account for the data that is
    // already there (i.e. simply add to it or remove it if necessary).
    // Since the copy_local_to_global operation is a "+=" operation,
    // we need to take this into account
    //
    // For the K_uu block in particular, this means that contributions have been
    // added from the surrounding cells, so we need to be careful when we manipulate this block.
    // We can't just erase the subblocks.
    //
    // So the intermediate matrix that we need to get from what we have in K_uu and what we
    // are actually wanting is:
    //  | K'_uu - K_uu |   0   |        0        |
    //  |       0      |   0   |  K_pt^-1 - K_pt |
    //  |       0      |   0   |        0        |
    //
    // This is the strategy we will employ to get the subblocks we want:
    // K'_{uu}: Since we don't have access to K_{uu}^h, but we know its contribution is added to the global
    //        K_{uu} matrix, we just want to add the element wise static-condensation
    //        K'_{uu}^h = K_{uu}^h + K_{up}^h K_{tp}^{-1, h} K_{tt}^h K_{pt}^{-1, h} K_{pu}^h
    //        Since we already have K_uu^h in the system matrix, we just need to do the following
    //        K'_{uu}^h == (K_{uu}^h += K_{up}^h K_{tp}^{-1}^h K_{tt}^h K_{pt}^{-1, h} K_{pu}^h)
    // K_{pt}^-1: Similarly, K_pt exists in the subblock. Since the copy operation is a += operation, we need
    //          to subtract the existing K_pt submatrix in addition to "adding" that which we wish to
    //          replace it with.
    // K_{tp}^-1: Since the global matrix is symmetric, this block is the same as the one above
    //          and we can simply use K_pt^-1 as a substitute for this one

    // We first extract element data from the system matrix. So first
    // we get the entire subblock for the cell
    AdditionalTools::extract_submatrix(data.local_dof_indices,
                                       data.local_dof_indices,
                                       tangent_matrix,
                                       data.K_orig);
    // and next the local matrices for K_{pu}, K_{pt} and K_{tt}
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

    // To get the inverse of K_{pt}, we invert it directly.
    // This operation is relatively inexpensive since
    // K_{pt} is block-diagonal.
    data.K_pt_inv.invert(data.K_pt);

    // Now we can make condensation terms to add to the
    // K_{uu} block and put them in the cell local matrix
    data.K_pt_inv.mmult(data.A, data.K_pu);
    data.K_tt.mmult(data.B, data.A);
    data.K_pt_inv.Tmmult(data.C, data.B);
    data.K_pu.Tmmult(data.K_con, data.C);
    AdditionalTools::replace_submatrix(element_indices_u,
                                       element_indices_u,
                                       data.K_con,
                                       data.cell_matrix);

    // Next we place K_{pt}^-1 in the K_{pt} block for post-processing
    // Note again that we need to remove the K_pt contribution that
    // already exists there.
    data.K_pt_inv.add (-1.0, data.K_pt);
    AdditionalTools::replace_submatrix(element_indices_p,
                                       element_indices_t,
                                       data.K_pt_inv,
                                       data.cell_matrix);
}

// @sect4{Solid::output_results}
// Here we present how the results are written to file to be viewed
// using Paraview. The method is similar to that shown in previous
// tutorials so will not be discussed in detail.
template <int dim>
void Solid<dim>::output_results(void)
{
    DataOut<dim> data_out;
    std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation (dim,
                                                                                                         DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back (DataComponentInterpretation::component_is_scalar);
    data_component_interpretation.push_back (DataComponentInterpretation::component_is_scalar);

    std::vector<std::string> solution_name (dim,
                                            "displacement");
    solution_name.push_back ("pressure");
    solution_name.push_back ("dilatation");

    data_out.attach_dof_handler (dof_handler_ref);
    data_out.add_data_vector (solution_n,
			      solution_name,
			      DataOut<dim>::type_dof_data,
			      data_component_interpretation);

    // Since we are dealing with a large deformation problem, it would be nice
    // to display the result on a displaced grid! The MappingQEulerian class
    // linked with the DataOut class provides an interface through which this
    // can be achieved without physically moving the grid points ourselves.
    // We first need to copy the solution to a temporary vector and then
    // create the Eularian mapping. We also specify the polynomial degree
    // to the DataOut object in order to produce a more refined output dataset
    // when higher order polynomials are used.
    Vector<double> soln (solution_n.size());
    for (unsigned int i=0; i < soln.size(); ++i)
        soln(i) = solution_n(i);
    MappingQEulerian<dim> q_mapping (degree,
                                     soln,
                                     dof_handler_ref);
    data_out.build_patches (q_mapping,
                            degree);

    std::ostringstream filename;
    filename << "solution-"
             << time.get_timestep()
             << ".vtk";

    std::ofstream output (filename.str().c_str());
    data_out.write_vtk (output);
}

// @sect3{Main function}
// Lastly we provide the main driver function which appears
// no different to the other tutorials.
int main (void)
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

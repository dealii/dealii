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

// We start by including all the necessary
// deal.II header files and some C++ related
// ones. They have been discussed in detail
// in previous tutorial programs, so you need
// only refer to past tutorials for details.
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_constraints.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/fe/fe_dgp_monomial.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vectors.h>

#include <iostream>
#include <fstream>

// Next we import all the deal.II
// function and class names to the global namespace
using namespace dealii;

// @sect3{Run-time parameters}
//
// There are several parameters that can be set
// in the code so we set up a ParameterHandler
// object to read in the choices at run-time.
namespace Parameters {
// @sect4{Finite Element system}
// As mentioned in the introduction, a different order
// interpolation should be used for the displacement
// $\mathbf{u}$ than for the pressure $p$ and
// the dilatation $\widetilde{J}$.
// Choosing $p$ and $\widetilde{J}$ as discontinuous (constant)
// functions at the element level leads to the
// mean-dilatation method. The discontinuous approximation
// allows $p$ and $\widetilde{J}$ to be condensed out
// and a classical displacement based method is recovered.
// Here we specify the polynomial order used to
// approximate the solution.
// The quadrature order should be adjusted accordingly.
struct FESystem {
	int poly_degree;
	int quad_order;

	static void
	declare_parameters(ParameterHandler &prm);
	void
	parse_parameters(ParameterHandler &prm);
};

void FESystem::declare_parameters(ParameterHandler &prm) {
	prm.enter_subsection("Finite element system");
	{
		prm.declare_entry("Polynomial degree", "2", Patterns::Integer(),
				"Displacement system polynomial order");

		prm.declare_entry("Quadrature order", "3", Patterns::Integer(),
				"Gauss quadrature order");
	}
	prm.leave_subsection();
}

void FESystem::parse_parameters(ParameterHandler &prm) {
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
struct Geometry {
	int global_refinement;
	double scale;
	double p_p0;

	static void
	declare_parameters(ParameterHandler &prm);
	void
	parse_parameters(ParameterHandler &prm);
};

void Geometry::declare_parameters(ParameterHandler &prm) {
	prm.enter_subsection("Geometry");
	{
		prm.declare_entry("Global refinement", "2", Patterns::Integer(),
				"Global refinement level");

		prm.declare_entry("Grid scale", "1e-3", Patterns::Double(),
				"Global grid scaling factor");

		prm.declare_entry("Pressure ratio p/p0", "100",
				Patterns::Selection("20|40|60|80|100"),
				"Ratio of applied pressure to reference pressure");
	}
	prm.leave_subsection();
}

void Geometry::parse_parameters(ParameterHandler &prm) {
	prm.enter_subsection("Geometry");
	{
		global_refinement = prm.get_integer("Global refinement");
		scale = prm.get_double("Grid scale");
		p_p0 = prm.get_double("Pressure ratio p/p0");
	}
	prm.leave_subsection();
}

// @sect4{Materials}
// Need the shear modulus $ \mu $
// and Poisson ration $ \nu $
// for the neo-Hookean material.
struct Materials {
	double nu;
	double mu;

	static void
	declare_parameters(ParameterHandler &prm);
	void
	parse_parameters(ParameterHandler &prm);
};

// ToDo: add a range check
void Materials::declare_parameters(ParameterHandler &prm) {
	prm.enter_subsection("Material properties");
	{
		prm.declare_entry("Poisson's ratio", "0.4999", Patterns::Double(),
				"Poisson's ratio");

		prm.declare_entry("Shear modulus", "80.194e6", Patterns::Double(),
				"Shear modulus");
	}
	prm.leave_subsection();
}

void Materials::parse_parameters(ParameterHandler &prm) {
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
struct LinearSolver {
	std::string type_lin;
	double tol_lin;
	double max_iterations_lin;
	double ssor_relaxation;

	static void
	declare_parameters(ParameterHandler &prm);
	void
	parse_parameters(ParameterHandler &prm);
};

void LinearSolver::declare_parameters(ParameterHandler &prm) {
	prm.enter_subsection("Linear solver");
	{
		prm.declare_entry("Solver type", "CG", Patterns::Selection("CG|Direct"),
				"Type of solver used to solve the linear system");

		prm.declare_entry("Residual", "1e-6", Patterns::Double(),
				"Linear solver residual (scaled by residual norm)");

		prm.declare_entry(
				"Max iteration multiplier",
				"1",
				Patterns::Double(),
				"Linear solver iterations (multiples of the system matrix size)");

		prm.declare_entry("SSOR Relaxation", "0.65", Patterns::Double(),
				"SSOR preconditioner relaxation value");
	}
	prm.leave_subsection();
}

void LinearSolver::parse_parameters(ParameterHandler &prm) {
	prm.enter_subsection("Linear solver");
	{
		type_lin = prm.get("Solver type");
		tol_lin = prm.get_double("Residual");
		max_iterations_lin = prm.get_double("Max iteration multiplier");
		ssor_relaxation = prm.get_double("SSOR Relaxation");
	}
	prm.leave_subsection();
}

// @sect4{Nonlinear solver}
// A Newton-Raphson scheme is used to
// solve the nonlinear system of governing equations.
// Define the tolerances and the maximum number of
// iterations for the Newton-Raphson nonlinear solver.
struct NonlinearSolver {
	unsigned int max_iterations_NR;
	double tol_f;
	double tol_u;

	static void
	declare_parameters(ParameterHandler &prm);
	void
	parse_parameters(ParameterHandler &prm);
};

void NonlinearSolver::declare_parameters(ParameterHandler &prm) {
	prm.enter_subsection("Nonlinear solver");
	{
		prm.declare_entry("Max iterations Newton-Raphson", "10",
				Patterns::Integer(),
				"Number of Newton-Raphson iterations allowed");

		prm.declare_entry("Tolerance force", "1.0e-9", Patterns::Double(),
				"Force residual tolerance");

		prm.declare_entry("Tolerance displacement", "1.0e-6",
				Patterns::Double(), "Displacement error tolerance");
	}
	prm.leave_subsection();
}

void NonlinearSolver::parse_parameters(ParameterHandler &prm) {
	prm.enter_subsection("Nonlinear solver");
	{
		max_iterations_NR = prm.get_integer("Max iterations Newton-Raphson");
		tol_f = prm.get_double("Tolerance force");
		tol_u = prm.get_double("Tolerance displacement");
	}
	prm.leave_subsection();
}

// @sect4{Time}
// Set the timestep size $ \varDelta t $
// and the simulation end-time.
struct Time {
	double delta_t;
	double end_time;

	static void
	declare_parameters(ParameterHandler &prm);
	void
	parse_parameters(ParameterHandler &prm);
};

void Time::declare_parameters(ParameterHandler &prm) {
	prm.enter_subsection("Time");
	{
		prm.declare_entry("End time", "1", Patterns::Double(), "End time");

		prm.declare_entry("Time step size", "0.1", Patterns::Double(),
				"Time step size");
	}
	prm.leave_subsection();
}

void Time::parse_parameters(ParameterHandler &prm) {
	prm.enter_subsection("Time");
	{
		end_time = prm.get_double("End time");
		delta_t = prm.get_double("Time step size");
	}
	prm.leave_subsection();
}

// @sect4{All parameters}
// Finally we consolidate all of the above structures into
// a single container that holds all of our run-time selections.
struct AllParameters: public FESystem,
public Geometry,
public Materials,
public LinearSolver,
public NonlinearSolver,
public Time

{
	AllParameters(const std::string & input_file);

	static void
	declare_parameters(ParameterHandler &prm);
	void
	parse_parameters(ParameterHandler &prm);
};

AllParameters::AllParameters(const std::string & input_file) {
	ParameterHandler prm;
	declare_parameters(prm);
	prm.read_input(input_file);
	parse_parameters(prm);
}

void AllParameters::declare_parameters(ParameterHandler &prm) {
	FESystem::declare_parameters(prm);
	Geometry::declare_parameters(prm);
	Materials::declare_parameters(prm);
	LinearSolver::declare_parameters(prm);
	NonlinearSolver::declare_parameters(prm);
	Time::declare_parameters(prm);
}

void AllParameters::parse_parameters(ParameterHandler &prm) {
	FESystem::parse_parameters(prm);
	Geometry::parse_parameters(prm);
	Materials::parse_parameters(prm);
	LinearSolver::parse_parameters(prm);
	NonlinearSolver::parse_parameters(prm);
	Time::parse_parameters(prm);
}
}

// @sect3{General tools}
// We need to perform some specific operations that are not defined
// in the deal.II library yet.
// We place these common operations
// in a separate namespace for convenience.
// We also include some widely used operators
namespace AdditionalTools {
// Define an operation that takes two
// symmetric second-order tensors
// $\mathbf{A}$ and $\mathbf{B}$
// such that their outer-product
// $ \mathbf{A} \overline{\otimes} \mathbf{B} \Rightarrow C_{ijkl} = A_{ik} B_{jl} $
template<int dim>
SymmetricTensor<4, dim> outer_product_T23(const SymmetricTensor<2, dim> & A,
		const SymmetricTensor<2, dim> & B) {
	SymmetricTensor<4, dim> A_ik_B_jl;

	for (unsigned int i = 0; i < dim; ++i) {
		for (unsigned int j = i; j < dim; ++j) {
			for (unsigned int k = 0; k < dim; ++k) {
				for (unsigned int l = k; k < dim; ++k) {
					A_ik_B_jl[i][j][k][l] += A[i][k] * B[j][l];
				}
			}
		}
	}

	return A_ik_B_jl;
}

// The  extract_submatrix function
// takes specific entries from a matrix,
// and copies them to a  sub_matrix.
// The copied entries are defined by the
// first two parameters which hold the
// row and columns to be extracted.
// The  matrix is automatically resized
// to size $ r \times c $.
template<typename MatrixType>
void extract_submatrix(const std::vector<unsigned int> &row_index_set,
		const std::vector<unsigned int> &column_index_set,
		const MatrixType &matrix, FullMatrix<double> &sub_matrix) {

	const unsigned int n_rows_submatrix = row_index_set.size();
	const unsigned int n_cols_submatrix = column_index_set.size();
	// check the size of the input vectors
	Assert(n_rows_submatrix > 0, ExcInternalError());
	Assert(n_cols_submatrix > 0, ExcInternalError());

	sub_matrix.reinit(n_rows_submatrix, n_cols_submatrix);

	for (unsigned int sub_row = 0; sub_row < n_rows_submatrix; ++sub_row) {
		const unsigned int row = row_index_set[sub_row];
		Assert(row<=matrix.m(), ExcInternalError());

		for (unsigned int sub_col = 0; sub_col < n_cols_submatrix; ++sub_col) {
			const unsigned int col = column_index_set[sub_col];
			Assert(col<=matrix.n(), ExcInternalError());

			sub_matrix(sub_row, sub_col) = matrix(row, col);
		}
	}
}

// As above, but to extract entries from
// a <code> BlockSparseMatrix </code>.
template<>
void extract_submatrix<dealii::BlockSparseMatrix<double> >(
		const std::vector<unsigned int> &row_index_set,
		const std::vector<unsigned int> &column_index_set,
		const dealii::BlockSparseMatrix<double> &matrix,
		FullMatrix<double> &sub_matrix) {

	const unsigned int n_rows_submatrix = row_index_set.size();
	const unsigned int n_cols_submatrix = column_index_set.size();

	// check the size of the input vectors
	Assert(n_rows_submatrix > 0, ExcInternalError());
	Assert(n_cols_submatrix > 0, ExcInternalError());

	sub_matrix.reinit(n_rows_submatrix, n_cols_submatrix);

	for (unsigned int sub_row = 0; sub_row < n_rows_submatrix; ++sub_row) {
		const unsigned int row = row_index_set[sub_row];
		Assert(row<=matrix.m(), ExcInternalError());

		for (unsigned int sub_col = 0; sub_col < n_cols_submatrix; ++sub_col) {
			const unsigned int col = column_index_set[sub_col];
			Assert(col<=matrix.n(), ExcInternalError());
			if (matrix.get_sparsity_pattern().exists(row, col) == false)
				continue;

			sub_matrix(sub_row, sub_col) = matrix(row, col);
		}
	}
}

// The replace_submatrix function takes
// specific entries from a  sub_matrix,
// and copies them into a  matrix.
// The copied entries are defined by the
// first two parameters which hold the
// row and column entries to be replaced.
// The matrix expected to be of the correct size.
template<typename MatrixType>
void replace_submatrix(const std::vector<unsigned int> &row_index_set,
		const std::vector<unsigned int> &column_index_set,
		const MatrixType &sub_matrix, FullMatrix<double> &matrix) {
	const unsigned int n_rows_submatrix = row_index_set.size();
	Assert(n_rows_submatrix<=sub_matrix.m(), ExcInternalError());
	const unsigned int n_cols_submatrix = column_index_set.size();
	Assert(n_cols_submatrix<=sub_matrix.n(), ExcInternalError());

	for (unsigned int sub_row = 0; sub_row < n_rows_submatrix; ++sub_row) {
		const unsigned int row = row_index_set[sub_row];
		Assert(row<=matrix.m(), ExcInternalError());

		for (unsigned int sub_col = 0; sub_col < n_cols_submatrix; ++sub_col) {
			const unsigned int col = column_index_set[sub_col];
			Assert(col<=matrix.n(), ExcInternalError());

			matrix(row, col) = sub_matrix(sub_row, sub_col);

		}
	}
}

// Define some frequently used
// second and fourth-order tensors:
template<int dim>
class StandardTensors {
public:

	// $\mathbf{I}$
	static SymmetricTensor<2, dim> const I;
	// $\mathbf{I} \otimes \mathbf{I}$
	static SymmetricTensor<4, dim> const IxI;
	// $\mathcal{S}$, note that as we only use
	// this fourth-order unit tensor to operate
	// on symmetric second-order tensors.
	// To maintain notation consistent with Holzapfel (2001)
	// we name the tensor $\mathcal{I}$
	static SymmetricTensor<4, dim> const II;
	// Fourth-order deviatoric such that
	// $\textrm{dev}(\bullet) = (\bullet) - (1/3)[(\bullet):\mathbf{I}]\mathbf{I}$
	static SymmetricTensor<4, dim> const dev_P;
};

template<int dim>
SymmetricTensor<2, dim> const StandardTensors<dim>::I = SymmetricTensor<2, dim>(
		unit_symmetric_tensor<dim>());
template<int dim>
SymmetricTensor<4, dim> const StandardTensors<dim>::IxI =
		SymmetricTensor<4, dim>(outer_product(I, I));
template<int dim>
SymmetricTensor<4, dim> const StandardTensors<dim>::II =
		SymmetricTensor<4, dim>(identity_tensor<dim>());
template<int dim>
SymmetricTensor<4, dim> const StandardTensors<dim>::dev_P = (II
		- (1.0 / dim) * IxI);
}

// @sect3{Time class}
// A simple class to store time data. Its
// functioning is transparent so no discussion is
// necessary. For simplicity we assume a constant
// time step size.
class Time {
public:
	Time(const double time_end, const double delta_t) :
		timestep(0),
		time_current(0.0),
		time_end(time_end),
		delta_t(delta_t) {
	}
	virtual ~Time(void) {
	}

	double current(void) const {
		return time_current;
	}
	double end(void) const {
		return time_end;
	}
	double get_delta_t(void) const {
		return delta_t;
	}
	unsigned int get_timestep(void) const {
		return timestep;
	}
	void increment(void) {
		time_current += delta_t;
		++timestep;
	}

private:
	unsigned int timestep;
	double time_current;
	const double time_end;
	const double delta_t;
};

// @sect3{Compressible neo-Hookean material}

// As discussed in the Introduction,
// Neo-Hookean materials are a
// type of hyperelastic materials.
// The entire domain is assumed 
// to be composed of a compressible neo-Hookean material. 
// This class defines
// the behaviour of this material. 
// Compressible neo-Hookean materials
// can be described by a strain-energy function (SEF)
// $ \Psi = \Psi_{\text{iso}}(\overline{\mathbf{b}}) + \Psi_{\text{vol}}(J) $.
//
// The isochoric response is given by
// $ \Psi_{\text{iso}}(\mathbf{b}) = c_{1} [\overline{I}_{1} - 3]  $
// where $ c_{1} = \frac{\mu}{2} $ and $\overline{I}_{1}$ is the first
// invariant of the left- or right- isochoric Cauchy-Green deformation tensors.
// That is $\overline{I}_1 :=\textrm{tr}(\overline{\mathbf{b}})$.
// In this example the SEF that governs the volumetric
// response is defined as
// $ \Psi_{\text{vol}}(\widetilde{J})  = \kappa \bigl[ \frac{1}{2} [ \widetilde{J}^{2} - 1 ] - \textrm{ln}( \widetilde{J}) ] \bigr]  $
// where $\kappa:= \lambda + 2/3 \mu$ is the bulk modulus and
// $\lambda$ is a Lame moduli.
template<int dim>
class Material_Compressible_Neo_Hook_Three_Field {
public:
	Material_Compressible_Neo_Hook_Three_Field(const double mu, const double nu) :
		kappa((2.0 * mu * (1.0 + nu)) / (3.0 * (1.0 - 2.0 * nu))),
		c_1(mu / 2.0),
		det_F(1.0),
		p_tilde(0.0),
		J_tilde(1.0),
		b_bar(AdditionalTools::StandardTensors<dim>::I) {
		Assert(kappa > 0, ExcInternalError());
	}

	~Material_Compressible_Neo_Hook_Three_Field(void) {
	}

	// The Kirchhoff stress tensor $\boldsymbol{\tau}$ is
	// the chosen stress measure.
	// Recall that
	// $\boldsymbol{\tau} = \chi_{*}(\mathbf{S})$, i.e.
	// $\boldsymbol{\tau} = \mathbf{F} \mathbf{S} \mathbf{F}^{T}$.
	// Furthermore,
	// $\boldsymbol{\tau} = 2 \mathbf{F} \frac{\partial \Psi(\mathbf{C})}{\partial \mathbf{C}}  \mathbf{F}^{T} = 2 \mathbf{b} \frac{\partial \Psi(\mathbf{b})}{\partial \mathbf{b}}$.
	// Therefore,
	// $\boldsymbol{\tau} = 2 \mathbf{b} \bigl[ \frac{\partial \Psi_{\text{iso}}(\mathbf{b})}{\partial \mathbf{b}} + \frac{\partial \Psi_{\text{vol}}(J)}{\partial J}\frac{\partial J}{\partial \mathbf{b}} \bigr] =  2 \mathbf{b} \frac{\partial \Psi_{\text{iso}}(\mathbf{b})}{\partial \mathbf{b}} + J\frac{\partial \Psi_{\text{vol}}(J)}{\partial J}\mathbf{I} $

	// We update the material model with various deformation
	// dependent data based on F
	void update_material_data(const Tensor<2, dim> & F,
			const double p_tilde_in,
			const double J_tilde_in
			) {
		det_F = determinant(F);
		b_bar = std::pow(det_F, -2.0 / 3.0) * symmetrize(F * transpose(F));
		p_tilde = p_tilde_in;
		J_tilde = J_tilde_in;

		// include a coupled of checks on the input data
		Assert(det_F > 0, ExcInternalError());
	}

	// Determine the Kirchhoff stress
	// $\boldsymbol{\tau} = \boldsymbol{\tau}_{\textrm{iso}} + \boldsymbol{\tau}_{\textrm{vol}}$
	SymmetricTensor<2, dim> get_tau(void) {
		return get_tau_iso() + get_tau_vol();
	}

	// The fourth-order elasticity tensor in the spatial setting
	// $\mathfrak{c}$ is calculated from the SEF $\Psi$ as
	// $ J \mathfrak{c}_{ijkl} = F_{iA} F_{jB} \mathfrak{C}_{ABCD} F_{kC} F_{lD}$
	// where
	// $ \mathfrak{C} = 4 \frac{\partial^2 \Psi(\mathbf{C})}{\partial \mathbf{C} \partial \mathbf{C}}$
	SymmetricTensor<4, dim> get_Jc(void) const {
		return get_Jc_vol() + get_Jc_iso();
	}

	// Derivative of the volumetric free energy wrt $\widetilde{J}$
	// return $\frac{\partial \Psi_{\text{vol}}(\widetilde{J})}{\partial \widetilde{J}}$
	double get_dPsi_vol_dJ(void) const {
		return kappa * (J_tilde - 1.0 / J_tilde);
	}

	// Second derivative of the volumetric free energy wrt $\widetilde{J}$
	// We need the following computation explicitly in the tangent
	// so we make it public.
	// calculate
	// $\frac{\partial^2 \Psi_{\textrm{vol}}(\widetilde{J})}{\partial \widetilde{J} \partial \widetilde{J}}$
	double get_d2Psi_vol_dJ2(void) const {
		return kappa * (1.0 + 1.0 / (J_tilde * J_tilde));
	}


	double get_det_F(void) const {
		return det_F;
	}

	double get_p_tilde(void) const {
		return p_tilde;
	}

	double get_J_tilde(void) const {
		return J_tilde;
	}

protected:
	// Model properties $\kappa$ and $c_1$
	const double kappa; // Bulk modulus
	const double c_1; // neo-Hookean model parameter

	// Model specific data that is convenient to store with the material
	double det_F;
	double p_tilde;
	double J_tilde;

	SymmetricTensor<2, dim> b_bar;

	// Determine the volumetric Kirchhoff stress
	// $\boldsymbol{\tau}_{\textrm{vol}}$
	SymmetricTensor<2, dim> get_tau_vol(void) const {
		return p_tilde * det_F * AdditionalTools::StandardTensors<dim>::I;
	}

	// Determine the isochoric Kirchhoff stress
	// $\boldsymbol{\tau}_{\textrm{iso}} = \mathcal{P}:\overline{\boldsymbol{\tau}}$
	SymmetricTensor<2, dim> get_tau_iso(void) const {
		return AdditionalTools::StandardTensors<dim>::dev_P * get_tau_bar();
	}

	// Determine the fictitious Kirchhoff stress
	SymmetricTensor<2, dim> get_tau_bar(void) const {
		return 2.0 * c_1 * b_bar;
	}

	// Calculate the volumetric part of the tangent $J \mathfrak{c}_\textrm{vol}$
	SymmetricTensor<4, dim> get_Jc_vol(void) const {

		return p_tilde * det_F
				* ( AdditionalTools::StandardTensors<dim>::IxI
						- (2.0 * AdditionalTools::StandardTensors<dim>::II) );
	}

	// Calculate the isochoric part of the tangent $J \mathfrak{c}_\textrm{iso}$
	SymmetricTensor<4, dim> get_Jc_iso(void) const {
		const SymmetricTensor<2, dim> tau_bar = get_tau_bar();
		const SymmetricTensor<2, dim> tau_iso = get_tau_iso();
		const SymmetricTensor<4, dim> tau_iso_x_I = outer_product(tau_iso,
				AdditionalTools::StandardTensors<dim>::I);
		const SymmetricTensor<4, dim> I_x_tau_iso = outer_product(
				AdditionalTools::StandardTensors<dim>::I, tau_iso);
		const SymmetricTensor<4, dim> c_bar = get_c_bar();

		return (2.0 / 3.0) * trace(tau_bar)
				* AdditionalTools::StandardTensors<dim>::dev_P
				- (2.0 / 3.0) * (tau_iso_x_I + I_x_tau_iso)
				+ AdditionalTools::StandardTensors<dim>::dev_P * c_bar
				* AdditionalTools::StandardTensors<dim>::dev_P;
	}

	// Calculate the fictitious elasticity tensor $\overline{\mathfrak{c}}$
	SymmetricTensor<4, dim> get_c_bar() const {
		SymmetricTensor<4, dim> c_bar;
		c_bar = 0.0;
		return c_bar;
	}
};

// @sect3{Quadrature point history}
// As seen in step-18, the <code> PointHistory </code> class offers
// a method for storing data at the quadrature points.
// We need to evaluate the Kirchhoff stress $\boldsymbol{\tau}$ and
// the tangent $J\mathfrak{c}$ at the quadrature points.

template<int dim>
class PointHistory {
public:
	PointHistory(void) :
		material(NULL),
		F_inv(AdditionalTools::StandardTensors<dim>::I),
		tau(SymmetricTensor<2, dim>()),
		d2Psi_vol_dJ2(0.0),
		dPsi_vol_dJ(0.0),
		Jc(SymmetricTensor<4, dim>()) {
	}
	virtual ~PointHistory(void) {
		delete material;
		material = NULL;
	}

	// We first create a material object.
	void setup_lqp(Parameters::AllParameters & parameters) {

		// Create an instance of a neo-Hookean material
		material = new Material_Compressible_Neo_Hook_Three_Field<dim>(
				parameters.mu, parameters.nu);

		// Initialise all tensors correctly
		update_values(Tensor<2, dim>(), 0.0, 1.0);
	}

	// Update the stored values and stresses based on the current
	// deformation configuration, pressure $p$ and
	// dilation $\widetilde{J}$ field values.
	// The input is the material gradient of the displacement
	// $\textrm{Grad}\mathbf{u}_{\textrm{n}}$
	void update_values(const Tensor<2, dim> & Grad_u_n,
			const double p_tilde,
			const double J_tilde) {
		// Store the calculated pressure $p$
		// and dilatation $\widetilde{J}$

		// Various deformation gradient $\mathbf{F}$ from the
		// displacement gradient $\textrm{Grad}\mathbf{u}$, i.e.
		// $\mathbf{F}(\mathbf{u}) = \mathbf{I} + \textrm{Grad} \mathbf{u}$
		static const Tensor<2, dim> I =
				static_cast<Tensor<2, dim> >(AdditionalTools::StandardTensors<
						dim>::I);
		const Tensor<2, dim> F = I + Grad_u_n;

		// We use the inverse of $\mathbf{F}$ frequently so we store it
		F_inv = invert(F);

		// Now we update the material model with the new deformation measures
		material->update_material_data(F, p_tilde, J_tilde);

		// The material has been updated so we now calculate the
		// Kirchhoff stress $\mathbf{\tau}$ and the tangent $J\mathfrak{c}$
		tau = material->get_tau();

		Jc = material->get_Jc();
		dPsi_vol_dJ = material->get_dPsi_vol_dJ();
		d2Psi_vol_dJ2 = material->get_d2Psi_vol_dJ2();

	}

	// We offer an interface to retrieve certain data.
	// Here are the kinematic variables
	double get_J_tilde(void) const {
		return material->get_J_tilde();
	}
	double get_det_F(void) const {
		return material->get_det_F();
	}
	Tensor<2, dim> get_F_inv(void) const {
		return F_inv;
	}

	// and the kinetic variables.
	// These are used in the material and global
	// tangent matrix and residual assembly operations
	// so we compute these and store them.
	double get_p_tilde(void) const {
		return material->get_p_tilde();
	}
	SymmetricTensor<2, dim> get_tau(void) const {
		return tau;
	}

	double get_dPsi_vol_dJ(void) const {
		return dPsi_vol_dJ;
	}

	double get_d2Psi_vol_dJ2(void) const {
		return d2Psi_vol_dJ2;
	}

	// and finally the tangent
	SymmetricTensor<4, dim> get_Jc(void) const {
		return Jc;
	}

private:
	// We specify that each QP has a copy of a material
	// type in case different materials are used
	// in different regions of the domain.
	// This also
	// deals with the issue of preventing data-races during
	// multi-threading operations when using shared objects.
	Material_Compressible_Neo_Hook_Three_Field<dim>* material;

	// These are all the volume, displacement and strain variables
	Tensor<2, dim> F_inv;

	// and the stress-type variables
	SymmetricTensor<2, dim> tau;
	double d2Psi_vol_dJ2;
	double dPsi_vol_dJ;

	// and the tangent
	SymmetricTensor<4, dim> Jc;
};

// @sect3{Quasi-static quasi-incompressible finite-strain solid}
template<int dim>
class Solid {
public:
	Solid(const std::string & input_file);
	virtual
	~Solid(void);
	void
	run(void);

private:

	// Threaded building-blocks data structures:
	// for the tangent matrix
	struct PerTaskData_K;
	struct ScratchData_K;
	// for the right-hand side
	struct PerTaskData_RHS;
	struct ScratchData_RHS;
	// for the static-condensation
	struct PerTaskData_SC;
	struct ScratchData_SC;
	// for the updating of the quadrature points
	struct PerTaskData_UQPH;
	struct ScratchData_UQPH;

	// Build the grid
	void
	make_grid(void);

	// Setup the Finite Element system to be solved
	void
	system_setup(void);
	void
	determine_component_extractors(void);

	// Assemble the system and right hand side matrices using multi-threading
	void
	assemble_system_tangent(void);
	void
	assemble_system_tangent_one_cell(
			const typename DoFHandler<dim>::active_cell_iterator & cell,
			ScratchData_K & scratch, PerTaskData_K & data);
	void
	copy_local_to_global_K(const PerTaskData_K & data);
	void
	assemble_system_rhs(void);
	void
	assemble_system_rhs_one_cell(
			const typename DoFHandler<dim>::active_cell_iterator & cell,
			ScratchData_RHS & scratch, PerTaskData_RHS & data);
	void
	copy_local_to_global_rhs(const PerTaskData_RHS & data);
	void
	assemble_sc(void);
	void
	assemble_sc_one_cell(
			const typename DoFHandler<dim>::active_cell_iterator & cell,
			ScratchData_SC & scratch, PerTaskData_SC & data);
	void
	copy_local_to_global_sc(const PerTaskData_SC & data);
	// Apply Dirichlet boundary values
	void
	make_constraints(const int & it_nr, ConstraintMatrix & constraints);

	// Create and update the quadrature points stress and strain values
	void
	setup_qph(void);
	void
	update_qph_incremental(const BlockVector<double> & solution_delta);
	void
	update_qph_incremental_one_cell(
			const typename DoFHandler<dim>::active_cell_iterator & cell,
			ScratchData_UQPH & scratch, PerTaskData_UQPH & data);
	void copy_local_to_global_UQPH(const PerTaskData_UQPH & data) {
	}

	// Solve for the displacement using a Newton-Rhapson method
	void
	solve_nonlinear_timestep(BlockVector<double> & solution_delta);
	std::pair<unsigned int, double>
	solve_linear_system(BlockVector<double> & newton_update);

	// Solution retrieval
	BlockVector<double>
	get_solution_total(const BlockVector<double> & solution_delta) const;

	// Post-processing and writing data to file
	void
	output_results(void) const;

	// A collection of the parameters used to describe the problem setup
	Parameters::AllParameters parameters;

	// The volume of the reference and current configurations
	double vol_reference;
	double vol_current;

	// Description of the geometry on which the problem is solved
	Triangulation<dim> triangulation;

	// Keep track of the current time and the time spent evaluating certain functions
	Time time;
	TimerOutput timer;

	// A storage object for quadrature point information
	std::vector<PointHistory<dim> > quadrature_point_history;

	// A description of the finite-element system including the displacement polynomial degree,
	// the degree-of-freedom handler, number of dof's per cell and the extractor objects used
	// to retrieve information from the solution vectors
	const unsigned int degree;
	const FESystem<dim> fe;
	DoFHandler<dim> dof_handler_ref;
	unsigned int dofs_per_cell;
	const FEValuesExtractors::Vector u_fe;
	const FEValuesExtractors::Scalar p_fe;
	const FEValuesExtractors::Scalar J_fe;

	// Description of how the block-system is arranged
	// There are 3 blocks, the first contains a vector DOF $\mathbf{u}$
	// while the other two describe scalar DOFs, $p$ and $\widetilde{J}$.
	static const unsigned int n_blocks = 3;
	static const unsigned int n_components = dim + 2;
	static const unsigned int first_u_component = 0;
	static const unsigned int p_component = dim;
	static const unsigned int J_component = dim + 1;

	enum {
		u_dof = 0, p_dof, J_dof
	};
	std::vector<unsigned int> dofs_per_block;
	std::vector<unsigned int> element_indices_u;
	std::vector<unsigned int> element_indices_p;
	std::vector<unsigned int> element_indices_J;

	// Rules for Gauss-quadrature on both the cell and faces. The
	// number of quadrature points on both cells and faces is
	// recorded.
	QGauss<dim> qf_cell;
	QGauss<dim - 1> qf_face;
	unsigned int n_q_points;
	unsigned int n_q_points_f;

	// Objects that store the converged solution and right-hand side vectors,
	// as well as the tangent matrix. There is a ConstraintMatrix object
	// used to keep track of constraints.
	ConstraintMatrix constraints;
	BlockSparsityPattern sparsity_pattern;
	BlockSparseMatrix<double> tangent_matrix;
	BlockVector<double> system_rhs;
	BlockVector<double> solution_n;

	// Then define a number of variables to store norms and update
	// norms and normalisation factors.
	struct Errors {
		Errors(void) :
			norm(1.0), u(1.0), p(1.0), J(1.0) {
		}
		double norm, u, p, J;
		void reset(void) {
			norm = 1.0;
			u = 1.0;
			p = 1.0;
			J = 1.0;
		}
		void normalise(const Errors & rhs) {
			if (rhs.norm != 0.0)
				norm /= rhs.norm;
			if (rhs.u != 0.0)
				u /= rhs.u;
			if (rhs.p != 0.0)
				p /= rhs.p;
			if (rhs.J != 0.0)
				J /= rhs.J;
		}
	} error_residual, error_residual_0, error_residual_norm, error_update,
	error_update_0, error_update_norm;

	// Methods to calculate error measures
	void
	get_error_residual(Errors & error_residual);
	void
	get_error_update(const BlockVector<double> & newton_update,
			Errors & error_update);
	double
	get_error_dil(void);

	// Print information to screen
	void
	print_conv_header(void);
	void
	print_conv_footer(void);
};

// @sect3{Implementation of the <code>Solid</code> class}

// @sect4{Public interface}
// We initialise the Solid class using data extracted
// from the parameter file.
template<int dim>
Solid<dim>::Solid(const std::string & input_file) :
parameters(input_file), triangulation(
		Triangulation<dim>::maximum_smoothing), time(
				parameters.end_time, parameters.delta_t), timer(std::cout,
						TimerOutput::summary, TimerOutput::wall_times), degree(
								parameters.poly_degree),
								// The Finite Element System is composed of dim continuous
								// displacement DOFs, and discontinuous pressure and
								// dilatation DOFs. In an attempt to satisfy the LBB conditions,
								// we setup a Q(n)-P(n-1)-P(n-1) system. Q2-P1-P1 elements satisfy
								// this condition, while Q1-P0-P0 elements do not. However, it
								// has been shown that the latter demonstrate good convergence
								// characteristics nonetheless.
								fe(FE_Q<dim>(parameters.poly_degree), dim, // displacement
										FE_DGPMonomial<dim>(parameters.poly_degree - 1), 1, // pressure
										FE_DGPMonomial<dim>(parameters.poly_degree - 1), 1), // dilatation
										dof_handler_ref(triangulation), u_fe(first_u_component), p_fe(
												p_component), J_fe(J_component), dofs_per_block(n_blocks), qf_cell(
														parameters.quad_order), qf_face(parameters.quad_order) {
	n_q_points = qf_cell.size();
	n_q_points_f = qf_face.size();
	dofs_per_cell = fe.dofs_per_cell;
	determine_component_extractors();
}

// The class destructor simply clears the data held by the DOFHandler
template<int dim>
Solid<dim>::~Solid(void) {
	dof_handler_ref.clear();
}

// In solving the quasi-static problem, the time
// becomes a loading parameter. We choose to increment
// time linearly using a constant time step size.
template<int dim>
void Solid<dim>::run(void) {
	// After preprocessing, we output the initial grid
	// before starting the simulation proper.
	make_grid();
	system_setup();
	output_results();
	time.increment();

	// Here we define
	// $\varDelta \mathbf{\Xi}:= \{\varDelta \mathbf{u},\varDelta p, \varDelta \widetilde{J} \}$.
	BlockVector<double> solution_delta(dofs_per_block);
	solution_delta.collect_sizes();

	// Now we loop over the time domain
	while (time.current() < time.end()) {
		// We need to reset the solution update
		// for this time step
		solution_delta = 0.0;

		// Solve the current time step and update total
		// solution vector
		solve_nonlinear_timestep(solution_delta);
		// $\varDelta \mathbf{\Xi}_{\textrm{n}} = \varDelta \mathbf{\Xi}_{\textrm{n-1}} + \varDelta \mathbf{\Xi}$
		solution_n += solution_delta;
		output_results();

		time.increment();
	}
}

// @sect3{Private interface}

// @sect4{Threaded-building-blocks structures}
// We use TBB to perform as many computationally intensive
// distributed tasks as possible. In particular, we assemble the
// tangent matrix and residual vector, the static
// condensation contributions, and update data stored
// at the quadrature points using TBB.

// Firstly we deal with the tangent matrix assembly structures.
// The PerTaskData object stores local contributions. 
template<int dim>
struct Solid<dim>::PerTaskData_K {
	FullMatrix<double> cell_matrix;
	std::vector<unsigned int> local_dof_indices;

	PerTaskData_K(const unsigned int dofs_per_cell) :
		cell_matrix(dofs_per_cell, dofs_per_cell), local_dof_indices(
				dofs_per_cell) {
	}

	void reset(void) {
		cell_matrix = 0.0;
	}
};
// while the ScratchData object stores the larger objects
// such as the shape-function values object and a shape function
// gradient and symmetric gradient vector which we will precompute later.
template<int dim>
struct Solid<dim>::ScratchData_K {
	FEValues<dim> fe_values_ref;

	std::vector<std::vector<double> > Nx;
	std::vector<std::vector<Tensor<2, dim> > > grad_Nx;
	std::vector<std::vector<SymmetricTensor<2, dim> > > symm_grad_Nx;

	ScratchData_K(const FiniteElement<dim> & fe_cell,
			const QGauss<dim> & qf_cell, const UpdateFlags uf_cell) :
				fe_values_ref(fe_cell, qf_cell, uf_cell), Nx(qf_cell.size(),
						std::vector<double>(fe_cell.dofs_per_cell)), grad_Nx(
								qf_cell.size(),
								std::vector<Tensor<2, dim> >(fe_cell.dofs_per_cell)), symm_grad_Nx(
										qf_cell.size(),
										std::vector<SymmetricTensor<2, dim> >(
												fe_cell.dofs_per_cell)) {
	}

	ScratchData_K(const ScratchData_K & rhs) :
		fe_values_ref(rhs.fe_values_ref.get_fe(),
				rhs.fe_values_ref.get_quadrature(),
				rhs.fe_values_ref.get_update_flags()), Nx(rhs.Nx), grad_Nx(
						rhs.grad_Nx), symm_grad_Nx(rhs.symm_grad_Nx) {
	}

	void reset(void) {
		const unsigned int n_q_points = Nx.size();
		const unsigned int n_dofs_per_cell = Nx[0].size();
		for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
			Assert( Nx[q_point].size() == n_dofs_per_cell, ExcInternalError());
			Assert( grad_Nx[q_point].size() == n_dofs_per_cell,
					ExcInternalError());
			Assert( symm_grad_Nx[q_point].size() == n_dofs_per_cell,
					ExcInternalError());
			for (unsigned int k = 0; k < n_dofs_per_cell; ++k) {
				Nx[q_point][k] = 0.0;
				grad_Nx[q_point][k] = 0.0;
				symm_grad_Nx[q_point][k] = 0.0;
			}
		}
	}

};

// Next are the same data structures used for the
// right-hand side assembly.
// The PerTaskData object again stores local contributions
template<int dim>
struct Solid<dim>::PerTaskData_RHS {
	Vector<double> cell_rhs;
	std::vector<unsigned int> local_dof_indices;

	PerTaskData_RHS(const unsigned int dofs_per_cell) :
		cell_rhs(dofs_per_cell), local_dof_indices(dofs_per_cell) {
	}

	void reset(void) {
		cell_rhs = 0.0;
	}
};
// and the ScratchData object the shape function object
// and precomputed values vector
template<int dim>
struct Solid<dim>::ScratchData_RHS {
	FEValues<dim> fe_values_ref;
	FEFaceValues<dim> fe_face_values_ref;

	std::vector<std::vector<double> > Nx;
	std::vector<std::vector<SymmetricTensor<2, dim> > > symm_grad_Nx;

	ScratchData_RHS(const FiniteElement<dim> & fe_cell,
			const QGauss<dim> & qf_cell, const UpdateFlags uf_cell,
			const QGauss<dim - 1> & qf_face, const UpdateFlags uf_face) :
				fe_values_ref(fe_cell, qf_cell, uf_cell), fe_face_values_ref(
						fe_cell, qf_face, uf_face), Nx(qf_cell.size(),
								std::vector<double>(fe_cell.dofs_per_cell)), symm_grad_Nx(
										qf_cell.size(),
										std::vector<SymmetricTensor<2, dim> >(
												fe_cell.dofs_per_cell)) {
	}

	ScratchData_RHS(const ScratchData_RHS & rhs) :
		fe_values_ref(rhs.fe_values_ref.get_fe(),
				rhs.fe_values_ref.get_quadrature(),
				rhs.fe_values_ref.get_update_flags()), fe_face_values_ref(
						rhs.fe_face_values_ref.get_fe(),
						rhs.fe_face_values_ref.get_quadrature(),
						rhs.fe_face_values_ref.get_update_flags()), Nx(rhs.Nx), symm_grad_Nx(
								rhs.symm_grad_Nx) {
	}

	void reset(void) {
		const unsigned int n_q_points = Nx.size();
		const unsigned int n_dofs_per_cell = Nx[0].size();
		for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
			Assert( Nx[q_point].size() == n_dofs_per_cell, ExcInternalError());
			Assert( symm_grad_Nx[q_point].size() == n_dofs_per_cell,
					ExcInternalError());
			for (unsigned int k = 0; k < n_dofs_per_cell; ++k) {
				Nx[q_point][k] = 0.0;
				symm_grad_Nx[q_point][k] = 0.0;
			}
		}
	}

};

// Here we define structures to assemble the statically
// condensed tangent matrix. Recall that we wish to solve 
// for a displacement-based formulation. 
// We do the condensation at the element
// level as the $p$ and $\widetilde{J}$
// fields are element-wise discontinuous. 
// As these operations are matrix-based, 
// we need to setup a number of matrices
// to store the local contributions from 
// a number of the tangent matrix sub-blocks.
// We place these in the PerTaskData struct.
template<int dim>
struct Solid<dim>::PerTaskData_SC {
	FullMatrix<double> cell_matrix;
	std::vector<unsigned int> local_dof_indices;

	// Calculation matrices (auto resized)
	FullMatrix<double> k_orig;
	FullMatrix<double> k_pu;
	FullMatrix<double> k_pJ;
	FullMatrix<double> k_JJ;
	// Calculation matrices (manual resized)
	FullMatrix<double> k_pJ_inv;
	FullMatrix<double> k_bbar;
	FullMatrix<double> A;
	FullMatrix<double> B;
	FullMatrix<double> C;

	PerTaskData_SC(const unsigned int dofs_per_cell, const unsigned int n_u,
			const unsigned int n_p, const unsigned int n_J) :
				cell_matrix(dofs_per_cell, dofs_per_cell),
				local_dof_indices(dofs_per_cell),
				k_orig(dofs_per_cell, dofs_per_cell),
				k_pu(n_p, n_u),
				k_pJ(n_p, n_J),
				k_JJ(n_J, n_J),
				k_pJ_inv(n_p, n_J),
				k_bbar(n_u, n_u),
				A(n_J,n_u),
				B(n_J, n_u),
				C(n_p, n_u) {
	}

	// Choose not to reset any data as the matrix extraction and
	// replacement tools will take care of this
	void reset(void) {
	}
};
// The ScratchData object is not strictly necessary for the
// operations we wish to perform, but it still needs to be defined for the
// current implementation of TBB in deal.II. 
// So we create a dummy struct for this purpose.
template<int dim>
struct Solid<dim>::ScratchData_SC {
	ScratchData_SC(void) {
	}
	ScratchData_SC(const ScratchData_SC & rhs) {
	}
	void reset(void) {
	}
};

// And finally we define the structures to assist with updating the quadrature
// point information. Similar to the SC assembly process, we choose not to use
// the PerTaskData object to store any information but must define one nonetheless.
template<int dim>
struct Solid<dim>::PerTaskData_UQPH {
	PerTaskData_UQPH(void) {
	}
	void reset(void) {
	}
};
// The ScratchData object will be used to store an alias for the solution vector
// so that we don't have to copy this large data structure. We then define
// a number of vectors to extract the solution values and gradients at the
// quadrature points.
template<int dim>
struct Solid<dim>::ScratchData_UQPH {

	 const BlockVector<double> & solution_total;

	std::vector<Tensor<2, dim> > solution_grads_u_total;
	std::vector<double> solution_values_p_total;
	std::vector<double> solution_values_J_total;

	FEValues<dim> fe_values_ref;

	ScratchData_UQPH(const FiniteElement<dim> & fe_cell,
			const QGauss<dim> & qf_cell,
			const UpdateFlags uf_cell,
			const BlockVector<double> & solution_total) :
				solution_total(solution_total),
				solution_grads_u_total(qf_cell.size()),
				solution_values_p_total(qf_cell.size()),
				solution_values_J_total(qf_cell.size()),
				fe_values_ref(fe_cell, qf_cell, uf_cell) {
	}

	ScratchData_UQPH(const ScratchData_UQPH & rhs) :
		solution_total(rhs.solution_total), solution_grads_u_total(
				rhs.solution_grads_u_total), solution_values_p_total(
						rhs.solution_values_p_total), solution_values_J_total(
								rhs.solution_values_J_total), fe_values_ref(
										rhs.fe_values_ref.get_fe(),
										rhs.fe_values_ref.get_quadrature(),
										rhs.fe_values_ref.get_update_flags()) {
	}

	void reset(void) {
		const unsigned int n_q_points = solution_grads_u_total.size();
		for (unsigned int q = 0; q < n_q_points; ++q) {
			solution_grads_u_total[q] = 0.0;
			solution_values_p_total[q] = 0.0;
			solution_values_J_total[q] = 0.0;
		}
	}
};

// @sect4{Solid::make_grid}
// Here we create the triangulation of the domain
template<int dim>
void Solid<dim>::make_grid(void) {
	// Create a unit cube with each face given a boundary ID number
	GridGenerator::hyper_rectangle(triangulation, Point<dim>(0.0, 0.0, 0.0),
			Point<dim>(1.0, 1.0, 1.0), true);
	GridTools::scale(parameters.scale, triangulation);

	// The grid must be refined at least once for the indentation problem
	if (parameters.global_refinement == 0)
		triangulation.refine_global(1);
	else
		triangulation.refine_global(parameters.global_refinement);

	// determine the volume of the reference configuration
	vol_reference = GridTools::volume(triangulation);
	vol_current = vol_reference;
	std::cout << "Grid:\n\t Reference volume: " << vol_reference << std::endl;

	// Since we wish to apply a Neumann BC to a patch on the top surface,
	// we must find the cell faces in this part of the domain and
	// mark them with a distinct boundary ID number
	typename Triangulation<dim>::active_cell_iterator cell =
			triangulation.begin_active(), endc = triangulation.end();
	for (; cell != endc; ++cell) {
		if (cell->at_boundary() == true) {
			for (unsigned int face = 0;
					face < GeometryInfo<dim>::faces_per_cell; ++face) {
				// Find faces on the +y surface
				if (cell->face(face)->at_boundary() == true
						&& cell->face(face)->center()[2]
						                              == 1.0 * parameters.scale) {
					if (cell->face(face)->center()[0] < 0.5 * parameters.scale
							&& cell->face(face)->center()[1]
							                              < 0.5 * parameters.scale) {
						cell->face(face)->set_boundary_indicator(6); // Set a new boundary id on a patch
					}
				}
			}
		}
	}
}

// @sect4{Solid::system_setup}
// Next we describe how the FE system is setup.
template<int dim>
void Solid<dim>::system_setup(void) {
	timer.enter_subsection("Setup system");

	// We first describe the number of components per block. Since the
	// displacement is a vector component, the first dim components
	// belong to it, while the next two describe scalar pressure and
	// dilatation DOFs.
	std::vector<unsigned int> block_component(n_components, u_dof); // Displacement
	block_component[p_component] = p_dof; // Pressure
	block_component[J_component] = J_dof; // Dilatation

	// DOF handler is then initialised and we renumber the grid in an
	// efficient manner. We also record the number of DOF's per block.
	dof_handler_ref.distribute_dofs(fe);
	DoFRenumbering::Cuthill_McKee(dof_handler_ref);
	DoFRenumbering::component_wise(dof_handler_ref, block_component);
	DoFTools::count_dofs_per_block(dof_handler_ref, dofs_per_block,
			block_component);

	std::cout << "Triangulation:"
			<< "\n\t Number of active cells: " << triangulation.n_active_cells()
			<< "\n\t Number of degrees of freedom: " << dof_handler_ref.n_dofs()
			<< std::endl;

	// Setup the sparsity pattern and tangent matrix
	tangent_matrix.clear();
	{
		const unsigned int n_dofs_u = dofs_per_block[u_dof];
		const unsigned int n_dofs_p = dofs_per_block[p_dof];
		const unsigned int n_dofs_J = dofs_per_block[J_dof];

		BlockCompressedSimpleSparsityPattern csp(n_blocks, n_blocks);

		csp.block(u_dof, u_dof).reinit(n_dofs_u, n_dofs_u);
		csp.block(u_dof, p_dof).reinit(n_dofs_u, n_dofs_p);
		csp.block(u_dof, J_dof).reinit(n_dofs_u, n_dofs_J);

		csp.block(p_dof, u_dof).reinit(n_dofs_p, n_dofs_u);
		csp.block(p_dof, p_dof).reinit(n_dofs_p, n_dofs_p);
		csp.block(p_dof, J_dof).reinit(n_dofs_p, n_dofs_J);

		csp.block(J_dof, u_dof).reinit(n_dofs_J, n_dofs_u);
		csp.block(J_dof, p_dof).reinit(n_dofs_J, n_dofs_p);
		csp.block(J_dof, J_dof).reinit(n_dofs_J, n_dofs_J);
		csp.collect_sizes();

		// In order to perform the static condensation efficiently,
		// we choose to exploit the symmetry of the the system matrix.
		// The global system matrix has the following structure
		//      | K_con |  K_up  |     0     |         | dU_u |         | R_u |
		// K =  | K_pu  |    0   |   K_pJ^-1 | , dU =  | dU_p | , R =   | R_p |
		//      |   0   |  K_Jp  |   K_JJ    |         | dU_J |         | R_J |
		// We optimise the sparsity pattern to reflect this structure
		// and prevent unnecessary data creation for the right-diagonal
		// block components.
		Table<2, DoFTools::Coupling> coupling(n_components, n_components);
		for (unsigned int ii = 0; ii < n_components; ++ii) {
			for (unsigned int jj = 0; jj < n_components; ++jj) {
				if (((ii < p_component) && (jj == J_component))
						|| ((ii == J_component) && (jj < p_component))
						|| ((ii == p_component) && (jj == p_component))) {
					coupling[ii][jj] = DoFTools::none;
				} else {
					coupling[ii][jj] = DoFTools::always;
				}
			}
		}
		DoFTools::make_sparsity_pattern(dof_handler_ref, coupling, csp,
				constraints, false);
		sparsity_pattern.copy_from(csp);
	}

	tangent_matrix.reinit(sparsity_pattern);

	// Setup storage vectors noting that the dilatation is unity
	// (i.e. $\widetilde{J} = 1$)
	// in the undeformed configuration
	system_rhs.reinit(dofs_per_block);
	system_rhs.collect_sizes();

	solution_n.reinit(dofs_per_block);
	solution_n.collect_sizes();
	solution_n.block(J_dof) = 1.0;

	// and finally set up the quadrature point history
	setup_qph();

	timer.leave_subsection();
}

// We next get information from the FE system
// that describes which local element DOFs are
// attached to which block component.
// This is used later to extract sub-blocks from the global matrix.
template<int dim>
void Solid<dim>::determine_component_extractors(void) {
	element_indices_u.clear();
	element_indices_p.clear();
	element_indices_J.clear();

	for (unsigned int k = 0; k < fe.dofs_per_cell; ++k) {
		// The next call has the FE System indicate to which block component
		// the current DOF is attached to.
		// Currently, the interpolation fields are setup such that
		// 0 indicates a displacement DOF, 1 a pressure DOF and 2 a dilatation DOF.
		const unsigned int k_group = fe.system_to_base_index(k).first.first;
		if (k_group == u_dof) {
			element_indices_u.push_back(k);
		} else if (k_group == p_dof) {
			element_indices_p.push_back(k);
		} else if (k_group == J_dof) {
			element_indices_J.push_back(k);
		} else {
			Assert(k_group <= J_dof, ExcInternalError());
		}
	}
}

// @sect4{Solid::setup_qph}
// The method used to store quadrature information is already described in
// step-18. Here we implement a similar setup for a SMP machine.
template<int dim>
void Solid<dim>::setup_qph(void) {
	std::cout << "    Setting up quadrature point data..." << std::endl;

	// Firstly the actual QPH data objects are created. This must be done
	// only once the grid is refined to its finest level.
	{
		triangulation.clear_user_data();

		{
			std::vector<PointHistory<dim> > tmp;
			tmp.swap(quadrature_point_history);
		}

		quadrature_point_history.resize(
				triangulation.n_active_cells() * n_q_points);

		unsigned int history_index = 0;
		for (typename Triangulation<dim>::active_cell_iterator cell =
				triangulation.begin_active(); cell != triangulation.end();
				++cell) {
			cell->set_user_pointer(&quadrature_point_history[history_index]);
			history_index += n_q_points;
		}

		Assert(history_index == quadrature_point_history.size(),
				ExcInternalError());
	}

	// Next we setup the initial QP data
	for (typename Triangulation<dim>::active_cell_iterator cell =
			triangulation.begin_active(); cell != triangulation.end(); ++cell) {
		PointHistory<dim>* lqph =
				reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

		Assert(lqph >= &quadrature_point_history.front(), ExcInternalError());
		Assert(lqph <= &quadrature_point_history.back(), ExcInternalError());

		// Setup any initial information at Gauss points
		for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
			lqph[q_point].setup_lqp(parameters);
		}
	}
}

// @sect4{Solid::update_qph_incremental}
// As the update of QP information occurs frequently and involves a number of
// expensive operations, we define a multi-threaded approach to distributing
// the task across a number of CPU cores.
template<int dim>
void Solid<dim>::update_qph_incremental(
		const BlockVector<double> & solution_delta) {
	timer.enter_subsection("Update QPH data");
	std::cout << " UQPH " << std::flush;

	// Firstly we need to obtain the total solution as it stands
	// at this Newton increment
	const BlockVector<double> solution_total(
			get_solution_total(solution_delta));

	// Next we create the initial copy of TBB objects
	const UpdateFlags uf_UQPH(update_values | update_gradients);
	PerTaskData_UQPH per_task_data_UQPH;
	ScratchData_UQPH scratch_data_UQPH(fe, qf_cell, uf_UQPH, solution_total);

	// and pass them and the one-cell update function to the workstream to be processed
	WorkStream::run(dof_handler_ref.begin_active(), dof_handler_ref.end(),
			*this, &Solid::update_qph_incremental_one_cell,
			&Solid::copy_local_to_global_UQPH, scratch_data_UQPH,
			per_task_data_UQPH);

	timer.leave_subsection();
}

// Now we describe how we extract data from the solution vector and pass it
// along to each QP storage object for processing.
template<int dim>
void Solid<dim>::update_qph_incremental_one_cell(
		const typename DoFHandler<dim>::active_cell_iterator & cell,
		ScratchData_UQPH & scratch, PerTaskData_UQPH & data) {
	PointHistory<dim>* lqph =
			reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

	Assert(lqph >= &quadrature_point_history.front(), ExcInternalError());
	Assert(lqph <= &quadrature_point_history.back(), ExcInternalError());

	Assert(scratch.solution_grads_u_total.size() == n_q_points,
			ExcInternalError());
	Assert(scratch.solution_values_p_total.size() == n_q_points,
			ExcInternalError());
	Assert(scratch.solution_values_J_total.size() == n_q_points,
			ExcInternalError());

	scratch.reset();

	// Firstly we need to find the values and gradients at quadrature points
	// inside the current cell
	scratch.fe_values_ref.reinit(cell);

	scratch.fe_values_ref[u_fe].get_function_gradients(scratch.solution_total,
			scratch.solution_grads_u_total);
	scratch.fe_values_ref[p_fe].get_function_values(scratch.solution_total,
			scratch.solution_values_p_total);
	scratch.fe_values_ref[J_fe].get_function_values(scratch.solution_total,
			scratch.solution_values_J_total);

	// and then we update each local QP
	// using the displacement gradient
	// and total pressure and dilatation solution values.
	for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
		lqph[q_point].update_values(scratch.solution_grads_u_total[q_point],
				scratch.solution_values_p_total[q_point],
				scratch.solution_values_J_total[q_point]);
	}
}

// @sect4{Solid::solve_nonlinear_timestep}
template<int dim>
void Solid<dim>::solve_nonlinear_timestep(
		BlockVector<double> & solution_delta) {
	//    timer.enter_subsection("Nonlinear solver");
	std::cout << std::endl << "Timestep " << time.get_timestep() << " @ "
			<< time.current() << "s" << std::endl;

	// We create a new vector to store the current Newton update step
	BlockVector<double> newton_update(dofs_per_block);
	newton_update.collect_sizes();

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
	for (unsigned int it_nr = 0; it_nr < parameters.max_iterations_NR;
			++it_nr) {
		// Print Newton iteration
		std::cout << " " << std::setw(2) << it_nr << " " << std::flush;

		// Since the problem is fully nonlinear and we are using a
		// full Newton method, the data stored in the tangent matrix
		// and right-hand side vector is not reusable and must be cleared
		// at each Newton step.
		tangent_matrix = 0.0;
		system_rhs = 0.0;

		// We initially build the right-hand side vector to check for convergence.
		// The unconstrained DOF's of the rhs vector hold the out-of-balance
		// forces. The building is done before assembling the system matrix as the latter
		// is an expensive operation and we can potentially avoid an extra
		// assembly process by not assembling the tangent matrix when convergence
		// is attained.
		assemble_system_rhs(); // Assemble RHS
		get_error_residual(error_residual);

		// We store the residual errors after the first iteration
		// in order to normalise by their value
		if (it_nr == 0)
			error_residual_0 = error_residual;

		// We can now determine the normalised residual error
		error_residual_norm = error_residual;
		error_residual_norm.normalise(error_residual_0);

		// Check for solution convergence
		if (it_nr > 0 && error_update_norm.u <= parameters.tol_u
				&& error_residual_norm.u <= parameters.tol_f) {
			std::cout << " CONVERGED! " << std::endl;
			print_conv_footer();
			return;
		}

		assemble_system_tangent(); // Assemble stiffness matrix
		make_constraints(it_nr, constraints); // Make boundary conditions
		constraints.condense(tangent_matrix, system_rhs); // Apply BC's

		const std::pair<unsigned int, double> lin_solver_output =
				solve_linear_system(newton_update);

		get_error_update(newton_update, error_update);
		if (it_nr == 0)
			error_update_0 = error_update;

		// We can now determine the normalised Newton update error
		error_update_norm = error_update;
		error_update_norm.normalise(error_update_0);

		// The current solution state is unacceptable, so we need to update
		// the solution increment for this time step, update all quadrature
		// point information pertaining to this new displacement and stress state
		// and continue iterating.
		solution_delta += newton_update;
		update_qph_incremental(solution_delta);

		std::cout << " | " << std::fixed << std::setprecision(3) << std::setw(7)
		<< std::scientific << lin_solver_output.first << "  "
		<< lin_solver_output.second << "  " << error_residual_norm.norm
		<< "  " << error_residual_norm.u << "  "
		<< error_residual_norm.p << "  " << error_residual_norm.J
		<< "  " << error_update_norm.norm << "  " << error_update_norm.u
		<< "  " << error_update_norm.p << "  " << error_update_norm.J
		<< "  " << std::endl;
	}

	throw(ExcMessage("No convergence in nonlinear solver!"));
}

// We print out data in a nice table that is updated
// on a per-iteration basis. Here we set up the table
// header
template<int dim>
void Solid<dim>::print_conv_header(void) {
	static const unsigned int l_width = 155;

	for (unsigned int i = 0; i < l_width; ++i)
		std::cout << "_";
	std::cout << std::endl;

	std::cout << "                 " << "SOLVER STEP" << "                  "
			<< " | " << " LIN_IT  " << " LIN_RES   " << " RES_NORM    "
			<< " RES_U    " << " RES_P     " << " RES_J    " << " NU_NORM     "
			<< " NU_U      " << " NU_P      " << " NU_J " << std::endl;

	for (unsigned int i = 0; i < l_width; ++i)
		std::cout << "_";
	std::cout << std::endl;
}
// and here the footer
template<int dim>
void Solid<dim>::print_conv_footer(void) {
	static const unsigned int l_width = 155;

	for (unsigned int i = 0; i < l_width; ++i)
		std::cout << "_";
	std::cout << std::endl;

	std::cout << "Relative errors:" << std::endl
			<< "Displacement:\t" << error_update.u / error_update_0.u << std::endl
			<< "Force: \t\t" << error_residual.u / error_residual_0.u << std::endl
			<< "Dilatation:\t" << get_error_dil() << std::endl
			<< "v / V_0:\t" << vol_current << " / " << vol_reference << " = " << vol_current / vol_reference << std::endl;

}

// Calculate how well the dilatation $\widetilde{J}$ 
// agrees with $J := \textrm{det}\mathbf{F}$
// from the $L^2$ error
// $ \bigl[ \int_{\Omega_0} {[ J - \widetilde{J}]}^{2}\textrm{d}V \bigr]^{1/2}$ 
// which is then normalised by the current volume
// $\int_{\Omega_0}  J ~\textrm{d}V = \int_\Omega  ~\textrm{d}v$.
template<int dim>
// ToDO: return the ratio of the reference and current volumes
double Solid<dim>::get_error_dil(void) {

	double dil_L2_error = 0.0;
	vol_current = 0.0;

	FEValues<dim> fe_values_ref(fe, qf_cell, update_JxW_values);

	for (typename Triangulation<dim>::active_cell_iterator cell =
				triangulation.begin_active(); cell != triangulation.end(); ++cell) {
		fe_values_ref.reinit(cell);

		PointHistory<dim>* lqph =
				reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

		Assert(lqph >= &quadrature_point_history.front(), ExcInternalError());
		Assert(lqph <= &quadrature_point_history.back(), ExcInternalError());

		for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {

			const double det_F_qp = lqph[q_point].get_det_F();
			const double J_tilde_qp = lqph[q_point].get_J_tilde();
			const double the_error_qp_squared = std::pow(
					(det_F_qp - J_tilde_qp), 2);
			const double JxW = fe_values_ref.JxW(q_point);

			dil_L2_error += the_error_qp_squared * JxW;
			vol_current += det_F_qp * JxW;
		}Assert(vol_current > 0, ExcInternalError());
	}

	return (std::sqrt(dil_L2_error));
}

// Determine the true residual error for the problem. 
// That is, determine the error in the residual for
// unconstrained dof.
template<int dim>
void Solid<dim>::get_error_residual(Errors & error_residual) {
	BlockVector<double> error_res(dofs_per_block);
	error_res.collect_sizes();

	// Need to ignore constrained DOFs
	for (unsigned int i = 0; i < dof_handler_ref.n_dofs(); ++i)
		if (!constraints.is_constrained(i))
			error_res(i) = system_rhs(i);

	error_residual.norm = error_res.l2_norm();
	error_residual.u = error_res.block(u_dof).l2_norm();
	error_residual.p = error_res.block(p_dof).l2_norm();
	error_residual.J = error_res.block(J_dof).l2_norm();
}

// Determine the true Newton update error for the problem
template<int dim>
void Solid<dim>::get_error_update(const BlockVector<double> & newton_update,
		Errors & error_update) {
	BlockVector<double> error_ud(dofs_per_block);
	error_ud.collect_sizes();

	// Need to ignore constrained DOFs as they have a prescribed
	// value
	for (unsigned int i = 0; i < dof_handler_ref.n_dofs(); ++i)
		if (!constraints.is_constrained(i))
			error_ud(i) = newton_update(i);

	error_update.norm = error_ud.l2_norm();
	error_update.u = error_ud.block(u_dof).l2_norm();
	error_update.p = error_ud.block(p_dof).l2_norm();
	error_update.J = error_ud.block(J_dof).l2_norm();
}

// This function provides the total solution, which is valid at any Newton step.
// This is required as, to reduce computational error, the total solution is
// only updated at the end of the timestep.
template<int dim>
BlockVector<double> Solid<dim>::get_solution_total(
		const BlockVector<double> & solution_delta) const {
	BlockVector<double> solution_total(solution_n);
	solution_total += solution_delta;
	return solution_total;

}

// @sect4{Solid::assemble_system_tangent}
// Since we use TBB for assembly, we simply setup a copy of the
// data structures required for the process and pass them, along
// with the memory addresses of the assembly functions to the
// WorkStream object for processing. Note that we must ensure that
// the matrix is reset before any assembly operations can occur.
template<int dim>
void Solid<dim>::assemble_system_tangent(void) {
	timer.enter_subsection("Assemble tangent matrix");
	std::cout << " ASM_K " << std::flush;

	tangent_matrix = 0.0;

	const UpdateFlags uf_cell(
			update_values | update_gradients | update_JxW_values);

	PerTaskData_K per_task_data(dofs_per_cell);
	ScratchData_K scratch_data(fe, qf_cell, uf_cell);

	WorkStream::run(dof_handler_ref.begin_active(), dof_handler_ref.end(),
			*this, &Solid::assemble_system_tangent_one_cell,
			&Solid::copy_local_to_global_K, scratch_data, per_task_data);

	timer.leave_subsection();
}

// This function adds the local contribution to the system matrix.
// Note that we choose not to use the constraint matrix to do the
// job for us because the tangent matrix and residual processes have
// been split up into two separate functions.
template<int dim>
void Solid<dim>::copy_local_to_global_K(const PerTaskData_K & data) {
	for (unsigned int i = 0; i < dofs_per_cell; ++i)
		for (unsigned int j = 0; j < dofs_per_cell; ++j)
			tangent_matrix.add(data.local_dof_indices[i],
					data.local_dof_indices[j], data.cell_matrix(i, j));
}

// Here we define how we assemble the tangent matrix contribution for a
// single cell.
template<int dim>
void Solid<dim>::assemble_system_tangent_one_cell(
		const typename DoFHandler<dim>::active_cell_iterator & cell,
		ScratchData_K & scratch, PerTaskData_K & data) {
	// We first need to reset and initialise some
	// of the data structures and retrieve some
	// basic information regarding the DOF numbering on this cell
	data.reset();
	scratch.reset();
	scratch.fe_values_ref.reinit(cell);
	cell->get_dof_indices(data.local_dof_indices);
	PointHistory<dim> *lqph =
			reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

	// We can precalculate the cell shape function values and gradients. Note that the
	// shape function gradients are defined wrt the current configuration.
	// That is
	// $\textrm{grad}\boldsymbol{\varphi} = \textrm{Grad}\boldsymbol{\varphi} \mathbf{F}^{-1}$
	static const SymmetricTensor<2, dim> I = AdditionalTools::StandardTensors<
			dim>::I;
	for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
		const Tensor<2, dim> F_inv = lqph[q_point].get_F_inv();
		for (unsigned int k = 0; k < dofs_per_cell; ++k) {
			const unsigned int k_group = fe.system_to_base_index(k).first.first;

			if (k_group == u_dof) {
				scratch.grad_Nx[q_point][k] =
						scratch.fe_values_ref[u_fe].gradient(k, q_point)
						* F_inv;
				scratch.symm_grad_Nx[q_point][k] = symmetrize(
						scratch.grad_Nx[q_point][k]);
			} else if (k_group == p_dof) {
				scratch.Nx[q_point][k] = scratch.fe_values_ref[p_fe].value(k,
						q_point);
			} else if (k_group == J_dof) {
				scratch.Nx[q_point][k] = scratch.fe_values_ref[J_fe].value(k,
						q_point);
			} else {
				Assert(k_group <= J_dof, ExcInternalError());
			}
		}
	}

	// Now we build the local cell stiffness matrix. Since the global and local system
	// matrices are symmetric, we can exploit this property by building only the lower
	// half of the local matrix and copying the values to the upper half.
	// So we only assemble half of the K_uu, K_pp (= 0), K_JJ blocks, while the whole
	// K_pJ, K_uJ (=0), K_up blocks are built.
	for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
		// We first extract some configuration dependent variables from our
		// QPH history objects that for the current q_point.
		// Get the current stress state $\boldsymbol{\tau}$
		const Tensor<2, dim> tau =
				static_cast<Tensor<2, dim> >(lqph[q_point].get_tau());
		const SymmetricTensor<4, dim> Jc = lqph[q_point].get_Jc();
		const double d2Psi_vol_dJ2 = lqph[q_point].get_d2Psi_vol_dJ2();
		const double det_F = lqph[q_point].get_det_F();

		// Next we define some aliases to make the assembly process easier to follow
		const std::vector<double> & N = scratch.Nx[q_point];
		const std::vector<SymmetricTensor<2, dim> > & symm_grad_Nx =
				scratch.symm_grad_Nx[q_point];
		const std::vector<Tensor<2, dim> > & grad_Nx = scratch.grad_Nx[q_point];
		const double JxW = scratch.fe_values_ref.JxW(q_point);

		for (unsigned int i = 0; i < dofs_per_cell; ++i) {
			const unsigned int component_i =
					fe.system_to_component_index(i).first;
			// Determine the dimensional component that matches the dof component (i.e. i % dim)
			const unsigned int i_group = fe.system_to_base_index(i).first.first;

			for (unsigned int j = 0; j <= i; ++j) {
				const unsigned int component_j =
						fe.system_to_component_index(j).first;
				const unsigned int j_group =
						fe.system_to_base_index(j).first.first;

				// This is the K_{uu} contribution. It comprises of a material
				// contribution and a geometrical stress contribution which is only
				// added along the local matrix diagonals
				if ((i_group == j_group) && (i_group == u_dof)) {
					// The material contribution:
					data.cell_matrix(i, j) += symm_grad_Nx[i] * Jc
							* symm_grad_Nx[j] * JxW;
					if (component_i == component_j) // geometrical stress contribution
						data.cell_matrix(i, j) += grad_Nx[i][component_i] * tau
						* grad_Nx[j][component_j] * JxW;
				}
				// Next is the K_{pu} contribution
				else if ((i_group == p_dof) && (j_group == u_dof)) {
					data.cell_matrix(i, j) += N[i] * det_F
							* (symm_grad_Nx[j]
							                * AdditionalTools::StandardTensors<dim>::I)
							                * JxW;
				}
				// and the K_{Jp} contribution
				else if ((i_group == J_dof) && (j_group == p_dof)) {
					data.cell_matrix(i, j) -= N[i] * N[j] * JxW;
				}
				// and lastly the K_{JJ} contribution
				else if ((i_group == j_group) && (i_group == J_dof)) {
					data.cell_matrix(i, j) += N[i] * d2Psi_vol_dJ2 * N[j] * JxW;
				} else
					Assert((i_group <= J_dof) && (j_group <= J_dof),
							ExcInternalError());
			}
		}
	}

	// Here we copy the lower half of the local matrix in the upper
	// half of the local matrix
	for (unsigned int i = 0; i < dofs_per_cell; ++i) {
		for (unsigned int j = i + 1; j < dofs_per_cell; ++j) {
			data.cell_matrix(i, j) = data.cell_matrix(j, i);
		}
	}
}

// @sect4{Solid::assemble_system_rhs}
// The assembly of the right-hand side process is similar to the
// tangent matrix, so we will not describe it in too much detail.
// Note that since we are describing a problem with Neumann BCs,
// we will need the face normals and so must specify this in the
// update flags.
template<int dim>
void Solid<dim>::assemble_system_rhs(void) {
	timer.enter_subsection("Assemble system right-hand side");
	std::cout << " ASM_R " << std::flush;

	system_rhs = 0.0;

	const UpdateFlags uf_cell(
			update_values | update_gradients | update_JxW_values);
	const UpdateFlags uf_face(
			update_values | update_normal_vectors | update_JxW_values);

	PerTaskData_RHS per_task_data(dofs_per_cell);
	ScratchData_RHS scratch_data(fe, qf_cell, uf_cell, qf_face, uf_face);

	WorkStream::run(dof_handler_ref.begin_active(), dof_handler_ref.end(),
			*this, &Solid::assemble_system_rhs_one_cell,
			&Solid::copy_local_to_global_rhs, scratch_data, per_task_data);

	timer.leave_subsection();
}

template<int dim>
void Solid<dim>::copy_local_to_global_rhs(const PerTaskData_RHS & data) {
	for (unsigned int i = 0; i < dofs_per_cell; ++i) {
		system_rhs(data.local_dof_indices[i]) += data.cell_rhs(i);
	}
}

template<int dim>
void Solid<dim>::assemble_system_rhs_one_cell(
		const typename DoFHandler<dim>::active_cell_iterator & cell,
		ScratchData_RHS & scratch, PerTaskData_RHS & data) {
	// Again we reset the data structures
	data.reset();
	scratch.reset();
	scratch.fe_values_ref.reinit(cell);
	cell->get_dof_indices(data.local_dof_indices);
	PointHistory<dim> *lqph =
			reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

	// and then precompute some shape function data
	for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
		const Tensor<2, dim> F_inv = lqph[q_point].get_F_inv();

		for (unsigned int k = 0; k < dofs_per_cell; ++k) {
			const unsigned int k_group = fe.system_to_base_index(k).first.first;

			if (k_group == u_dof) {
				scratch.symm_grad_Nx[q_point][k] = symmetrize(
						scratch.fe_values_ref[u_fe].gradient(k, q_point)
						* F_inv);
			} else if (k_group == p_dof) {
				scratch.Nx[q_point][k] = scratch.fe_values_ref[p_fe].value(k,
						q_point);
			} else if (k_group == J_dof) {
				scratch.Nx[q_point][k] = scratch.fe_values_ref[J_fe].value(k,
						q_point);
			} else
				Assert(k_group <= J_dof, ExcInternalError());
		}
	}

	// and can now assemble the right-hand side contribution
	for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
		// We fist retrieve data stored at the qp
		const SymmetricTensor<2, dim> tau = lqph[q_point].get_tau();
		const double det_F = lqph[q_point].get_det_F();
		const double J_tilde = lqph[q_point].get_J_tilde();
		const double p_tilde = lqph[q_point].get_p_tilde();
		const double dPsi_vol_dJ = lqph[q_point].get_dPsi_vol_dJ();

		// define some shortcuts
		const std::vector<double> & N = scratch.Nx[q_point];
		const std::vector<SymmetricTensor<2, dim> > & symm_grad_Nx =
				scratch.symm_grad_Nx[q_point];
		const double JxW = scratch.fe_values_ref.JxW(q_point);

		// We first compute the contributions from the internal forces.
		// Note, by definition of the rhs as the negative of the residual,
		// these contributions are subtracted.
		for (unsigned int i = 0; i < dofs_per_cell; ++i) {
			const unsigned int i_group = fe.system_to_base_index(i).first.first;
			// Add the contribution to the F_u block
			if (i_group == u_dof) {
				data.cell_rhs(i) -= (symm_grad_Nx[i] * tau) * JxW;
			}
			// the F_p block
			else if (i_group == p_dof) {
				data.cell_rhs(i) -= N[i] * (det_F - J_tilde) * JxW;
			}
			// and finally the F_J block
			else if (i_group == J_dof) {
				data.cell_rhs(i) -= N[i] * (dPsi_vol_dJ - p_tilde) * JxW;
			} else
				Assert(i_group <= J_dof, ExcInternalError());
		}
	}

	// Next we assemble the Neumann contribution. We first check to see
	// it the cell face exists on a boundary on which a traction is
	// applied and add the contribution if this is the case.
	if (cell->at_boundary() == true) {
		for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
				++face) {
			if (cell->face(face)->at_boundary() == true
					&& cell->face(face)->boundary_indicator() == 6) {
				scratch.fe_face_values_ref.reinit(cell, face);

				for (unsigned int f_q_point = 0; f_q_point < n_q_points_f;
						++f_q_point) {
					// We retrieve the face normal at this QP
					const Tensor<1, dim> & N =
							scratch.fe_face_values_ref.normal_vector(f_q_point);

					// and specify the traction in reference configuration. For this problem,
					// a defined pressure is applied in the reference configuration.
					// The direction of the applied traction is assumed
					// not to evolve with the deformation of the domain. The
					// traction is defined using the first Piola-Kirchhoff stress is simply
					// t_0 = P*N = (pI)*N = p*N
					// We choose to use the time variable to linearly ramp up the pressure
					// load.
					static const double p0 = -4.0
							/ (parameters.scale * parameters.scale);
					const double time_ramp = (time.current() / time.end());
					const double pressure = p0 * parameters.p_p0 * time_ramp;
					const Tensor<1, dim> traction = pressure * N;

					for (unsigned int i = 0; i < dofs_per_cell; ++i) {
						const unsigned int i_group =
								fe.system_to_base_index(i).first.first;

						if (i_group == u_dof) {
							// More shortcuts being assigned
							const unsigned int component_i =
									fe.system_to_component_index(i).first;
							const double Ni =
									scratch.fe_face_values_ref.shape_value(i,
											f_q_point);
							const double JxW = scratch.fe_face_values_ref.JxW(
									f_q_point);

							// And finally we can add the traction vector contribution to
							// the local RHS vector. Note that this contribution is present
							// on displacement DOFs only.
							data.cell_rhs(i) += (Ni * traction[component_i])
											* JxW;
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
template<int dim>
void Solid<dim>::make_constraints(const int & it_nr,
		ConstraintMatrix & constraints) {
	std::cout << " CST " << std::flush;

	// Since the constraints are different at Newton iterations,
	// we need to clear the constraints matrix and completely
	// rebuild it. However, after the first iteration, the
	// constraints remain the same and we can simply skip the
	// rebuilding step if we do not clear it.
	if (it_nr > 1)
		return;
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

		std::vector<bool> components(n_components, false);
		components[0] = true;

		if (apply_dirichlet_bc == true) {
			VectorTools::interpolate_boundary_values(dof_handler_ref,
					boundary_id, ZeroFunction<dim>(n_components), constraints,
					components);
		} else {
			VectorTools::interpolate_boundary_values(dof_handler_ref,
					boundary_id, ZeroFunction<dim>(n_components), constraints,
					components);
		}
	}
	{
		const int boundary_id = 2;

		std::vector<bool> components(n_components, false);
		components[1] = true;

		if (apply_dirichlet_bc == true) {
			VectorTools::interpolate_boundary_values(dof_handler_ref,
					boundary_id, ZeroFunction<dim>(n_components), constraints,
					components);
		} else {
			VectorTools::interpolate_boundary_values(dof_handler_ref,
					boundary_id, ZeroFunction<dim>(n_components), constraints,
					components);
		}
	}
	{
		const int boundary_id = 4;
		std::vector<bool> components(n_components, false);
		components[2] = true;

		if (apply_dirichlet_bc == true) {
			VectorTools::interpolate_boundary_values(dof_handler_ref,
					boundary_id, ZeroFunction<dim>(n_components), constraints,
					components);
		} else {
			VectorTools::interpolate_boundary_values(dof_handler_ref,
					boundary_id, ZeroFunction<dim>(n_components), constraints,
					components);
		}
	}
	{
		const int boundary_id = 5;
		std::vector<bool> components(n_components, true);
		components[2] = false;

		if (apply_dirichlet_bc == true) {
			VectorTools::interpolate_boundary_values(dof_handler_ref,
					boundary_id, ZeroFunction<dim>(n_components), constraints,
					components);
		} else {
			VectorTools::interpolate_boundary_values(dof_handler_ref,
					boundary_id, ZeroFunction<dim>(n_components), constraints,
					components);
		}
	}
	{
		const int boundary_id = 6;
		std::vector<bool> components(n_components, true);
		components[2] = false;

		if (apply_dirichlet_bc == true) {
			VectorTools::interpolate_boundary_values(dof_handler_ref,
					boundary_id, ZeroFunction<dim>(n_components), constraints,
					components);
		} else {
			VectorTools::interpolate_boundary_values(dof_handler_ref,
					boundary_id, ZeroFunction<dim>(n_components), constraints,
					components);
		}
	}

	constraints.close();
}

// @sect4{Solid::solve_linear_system}
// Solving the entire block system is a bit problematic as there are no
// contributions to the K_{JJ} block, rendering it non-invertible.
// Since the pressure and dilatation variables DOFs are discontinuous, we can
// condense them out to form a smaller displacement-only system which
// we will then solve and subsequently post-process to retrieve the
// pressure and dilatation solutions.
template<int dim>
std::pair<unsigned int, double> Solid<dim>::solve_linear_system(
		BlockVector<double> & newton_update) {
	// Need two temporary vectors to help
	// with the static condensation.
	BlockVector<double> A(dofs_per_block);
	BlockVector<double> B(dofs_per_block);
	A.collect_sizes();
	B.collect_sizes();

	// Store the number of linear solver iterations
	// the (hopefully converged) residual
	unsigned int lin_it = 0;
	double lin_res = 0.0;

	//      		| K_con |   K_up  |     0     |         | du |         | F_u |
	// K_store =  	| K_pu  |     0   |   K_pJ^-1 | , dXi = | dp | , R =   | F_p |
	//      		|   0   |   K_Jp  |   K_JJ    |         | dJ |         | F_J |

	// Solve for the incremental displacement du
	{
		// Perform static condensation to make
		// K_con = K_uu + K_bbar,
		// and put K_pJ^{-1} in the original K_pJ block.
		// That is, we make K_store.
		assemble_sc();

		// K_con du = F_con
		// with F_con = F_u + K_up [- K_Jp^-1 F_j + K_bar F_p]
		// Assemble the RHS vector to solve for du
		// A_J = K_pJ^-1 F_p
		tangent_matrix.block(p_dof, J_dof).vmult(A.block(J_dof),
				system_rhs.block(p_dof));
		// B_J = K_JJ  K_pJ^-1  F_p
		tangent_matrix.block(J_dof, J_dof).vmult(B.block(J_dof),
				A.block(J_dof));
		// A_J = F_J - K_JJ  K_pJ^-1  F_p
		A.block(J_dof).equ(1.0, system_rhs.block(J_dof), -1.0, B.block(J_dof));
		// A_p = K_Jp^-1 [  F_J - K_JJ  K_pJ^-1  F_p ]
		tangent_matrix.block(p_dof, J_dof).Tvmult(A.block(p_dof),
				A.block(J_dof));
		// A_u = K_up  K_Jp^-1 [  F_J - K_JJ  K_pJ^-1  F_p ]
		tangent_matrix.block(u_dof, p_dof).vmult(A.block(u_dof),
				A.block(p_dof));
		// F_con = F_u -  K_up  K_Jp^-1 [  F_J - K_JJ  K_pJ^-1  F_p ]
		system_rhs.block(u_dof) -= A.block(u_dof);

		timer.enter_subsection("Linear solver");
		std::cout << " SLV " << std::flush;
		if (parameters.type_lin == "CG") {
			const int solver_its = tangent_matrix.block(u_dof, u_dof).m()
							* parameters.max_iterations_lin;
			const double tol_sol = parameters.tol_lin
					* system_rhs.block(u_dof).l2_norm();

			SolverControl solver_control(solver_its, tol_sol);

			GrowingVectorMemory<Vector<double> > GVM;
			SolverCG<Vector<double> > solver_CG(solver_control, GVM);

			// We've chosen a SSOR preconditioner as it appears to provide
			// the fastest solver convergence characteristics for this problem.
			PreconditionSSOR<SparseMatrix<double> > preconditioner;
			preconditioner.initialize(tangent_matrix.block(u_dof, u_dof),
					parameters.ssor_relaxation);

			solver_CG.solve(tangent_matrix.block(u_dof, u_dof),
					newton_update.block(u_dof), system_rhs.block(u_dof),
					preconditioner);

			lin_it = solver_control.last_step();
			lin_res = solver_control.last_value();
		} else if (parameters.type_lin == "Direct") {
			// Otherwise if the problem is small enough, a direct solver
			// can be utilised.
			SparseDirectUMFPACK A_direct;
			A_direct.initialize(tangent_matrix.block(u_dof, u_dof));
			A_direct.vmult(newton_update.block(u_dof), system_rhs.block(u_dof));

			lin_it = 1;
			lin_res = 0.0;
		} else
			throw(ExcMessage("Linear solver type not implemented"));
		timer.leave_subsection();
	}

	// distribute the constrained dof back to the Newton update
	constraints.distribute(newton_update);

	timer.enter_subsection("Linear solver postprocessing");
	std::cout << " PP " << std::flush;

	// Now that we've solved the displacement problem, we can post-process
	// to get the dilatation solution from the substitution
	// dJ = KpJ^{-1} (F_p - K_pu du )
	{
		// A_p  = K_pu du
		tangent_matrix.block(p_dof, u_dof).vmult(A.block(p_dof),
				newton_update.block(u_dof));
		// A_p  = -K_pu du
		A.block(p_dof) *= -1.0;
		// A_p  = F_p - K_pu du
		A.block(p_dof) += system_rhs.block(p_dof);
		// d_J = K_pJ^{-1} [ F_p - K_pu du ]
		tangent_matrix.block(p_dof, J_dof).vmult(newton_update.block(J_dof),
				A.block(p_dof));
	}

	constraints.distribute(newton_update);

	// and finally we solve for the pressure update with the substitution
	// dp = KJp^{-1} [ R_J - K_JJ dJ ]
	{
		// A_J = K_JJ dJ
		tangent_matrix.block(J_dof, J_dof).vmult(A.block(J_dof),
				newton_update.block(J_dof));
		// A_J = -K_JJ dJ
		A.block(J_dof) *= -1.0;
		// A_J = F_J - K_JJ dJ
		A.block(J_dof) += system_rhs.block(J_dof);
		// dp = K_Jp^{-1}   [F_J - K_JJ dJ]
		tangent_matrix.block(p_dof, J_dof).Tvmult(newton_update.block(p_dof),
				A.block(J_dof));
	}

	// distribute the constrained dof back to the Newton update
	constraints.distribute(newton_update);

	timer.leave_subsection();

	return std::make_pair(lin_it, lin_res);
}

// @sect4{Solid::assemble_system_SC}
// The static condensation process could be performed at a global level
// but we need the inverse of one of the blocks. However, since the
// pressure and dilatation variables are discontinuous, the SC operation
// can be done on a per-cell basis and we can produce the inverse of the
// block-diagonal K_{pt} block by inverting the local blocks. We can
// again use TBB to do this since each operation will be independent of
// one another.
template<int dim>
void Solid<dim>::assemble_sc(void) {
	timer.enter_subsection("Perform static condensation");
	std::cout << " ASM_SC " << std::flush;

	PerTaskData_SC per_task_data(dofs_per_cell, element_indices_u.size(),
			element_indices_p.size(), element_indices_J.size()); // Initialise members of per_task_data to the correct sizes.
	ScratchData_SC scratch_data;

	// Using TBB, we assemble the contributions to add to
	// K_uu to form K_con from each elements contributions.
	// These contributions are then added to the glabal stiffness
	// matrix.
	WorkStream::run(dof_handler_ref.begin_active(), dof_handler_ref.end(),
			*this, &Solid::assemble_sc_one_cell,
			&Solid::copy_local_to_global_sc, scratch_data, per_task_data);

	timer.leave_subsection();
}

// We need to describe how to add the local contributions
// to K to form K_store
template<int dim>
void Solid<dim>::copy_local_to_global_sc(const PerTaskData_SC & data) {
	for (unsigned int i = 0; i < dofs_per_cell; ++i)
		for (unsigned int j = 0; j < dofs_per_cell; ++j)
			tangent_matrix.add(data.local_dof_indices[i],
					data.local_dof_indices[j], data.cell_matrix(i, j));
}

// Now we describe the static condensation process.
template<int dim>
void Solid<dim>::assemble_sc_one_cell(
		const typename DoFHandler<dim>::active_cell_iterator & cell,
		ScratchData_SC & scratch, PerTaskData_SC & data) {
	// As per usual, we must first find out which global numbers the
	// degrees of freedom on this cell have and reset some data structures
	data.reset();
	scratch.reset();
	cell->get_dof_indices(data.local_dof_indices);

	// We now extract the contribution of
	// the  dof associated with the current cell
	// to the global stiffness matrix.
	// The discontinuous nature of the p and J
	// interpolations mean that their is no
	// coupling of the local contributions at the
	// global level. This is not the case with the u dof.
	// In other words, k_Jp, k_pJ and k_JJ, when extracted
	// from the global stiffness matrix are the element
	// contributions. This is not the case for k_uu.

	// Currently the matrix corresponding to
	// the dof associated with the current element
	// (denoted somewhat loosely as k) is of the form
	//  | k_uu  |   k_up   |   0   |
	//  | k_pu  |     0    |  k_pJ |
	//  |   0   |   k_Jp   |  k_JJ |
	//
	// We now need to modify it such that it appear as
	//  | k_con |   k_up   |     0     |
	//  | k_pu  |     0    |   k_pJ^-1 |
	//  |   0   |   k_Jp   |   k_JJ    |
	// with k_con = k_uu + k_bbar
	// where
	// k_bbar = k_up k_bar k_pu
	// and
	// k_bar = k_Jp^{-1} k_JJ kpJ^{-1}
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

	// extract k for the dof associated with the current element
	AdditionalTools::extract_submatrix(data.local_dof_indices,
			data.local_dof_indices, tangent_matrix, data.k_orig);
	// and next the local matrices for k_pu, k_pJ and k_JJ
	AdditionalTools::extract_submatrix(element_indices_p, element_indices_u,
			data.k_orig, data.k_pu);
	AdditionalTools::extract_submatrix(element_indices_p, element_indices_J,
			data.k_orig, data.k_pJ);
	AdditionalTools::extract_submatrix(element_indices_J, element_indices_J,
			data.k_orig, data.k_JJ);

	// To get the inverse of k_pJ, we invert it directly.
	// This operation is relatively inexpensive since
	// k_pJ is block-diagonal.
	data.k_pJ_inv.invert(data.k_pJ);

	// Now we can make condensation terms to add to the
	// k_uu block and put them in the cell local matrix
	// A = k_pJ^-1 k_pu
	data.k_pJ_inv.mmult(data.A, data.k_pu);
	// B = k_JJ k_pJ^-1 k_pu
	data.k_JJ.mmult(data.B, data.A);
	// C = k_Jp^-1 k_JJ k_pJ^-1 k_pu
	data.k_pJ_inv.Tmmult(data.C, data.B);
	// k_bbar = k_up k_Jp^-1 k_JJ k_pJ^-1 k_pu
	data.k_pu.Tmmult(data.k_bbar, data.C);
	AdditionalTools::replace_submatrix(element_indices_u, element_indices_u,
			data.k_bbar, data.cell_matrix);

	// Next we place k_{pJ}^-1 in the k_{pJ} block for post-processing.
	// Note again that we need to remove the k_pJ contribution that
	// already exists there.
	data.k_pJ_inv.add(-1.0, data.k_pJ);
	AdditionalTools::replace_submatrix(element_indices_p, element_indices_J,
			data.k_pJ_inv, data.cell_matrix);
}

// @sect4{Solid::output_results}
// Here we present how the results are written to file to be viewed
// using ParaView. The method is similar to that shown in previous
// tutorials so will not be discussed in detail.
template<int dim>
void Solid<dim>::output_results(void) const {
	DataOut<dim> data_out;
	std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation(
			dim, DataComponentInterpretation::component_is_part_of_vector);
	data_component_interpretation.push_back(
			DataComponentInterpretation::component_is_scalar);
	data_component_interpretation.push_back(
			DataComponentInterpretation::component_is_scalar);

	std::vector<std::string> solution_name(dim, "displacement");
	solution_name.push_back("pressure");
	solution_name.push_back("dilatation");

	data_out.attach_dof_handler(dof_handler_ref);
	data_out.add_data_vector(solution_n, solution_name,
			DataOut<dim>::type_dof_data, data_component_interpretation);

	// Since we are dealing with a large deformation problem, it would be nice
	// to display the result on a displaced grid! The MappingQEulerian class
	// linked with the DataOut class provides an interface through which this
	// can be achieved without physically moving the grid points ourselves.
	// We first need to copy the solution to a temporary vector and then
	// create the Eulerian mapping. We also specify the polynomial degree
	// to the DataOut object in order to produce a more refined output dataset
	// when higher order polynomials are used.
	Vector<double> soln(solution_n.size());
	for (unsigned int i = 0; i < soln.size(); ++i)
		soln(i) = solution_n(i);
	MappingQEulerian<dim> q_mapping(degree, soln, dof_handler_ref);
	data_out.build_patches(q_mapping, degree);

	std::ostringstream filename;
	filename << "solution-" << time.get_timestep() << ".vtk";

	std::ofstream output(filename.str().c_str());
	data_out.write_vtk(output);
}

// @sect3{Main function}
// Lastly we provide the main driver function which appears
// no different to the other tutorials.
int main(void) {
	try {
		deallog.depth_console(0);

		Solid<3> solid_3d("parameters.prm");
		solid_3d.run();
	} catch (std::exception &exc) {
		std::cerr << std::endl << std::endl
				<< "----------------------------------------------------"
				<< std::endl;
		std::cerr << "Exception on processing: " << std::endl << exc.what()
						<< std::endl << "Aborting!" << std::endl
						<< "----------------------------------------------------"
						<< std::endl;

		return 1;
	} catch (...) {
		std::cerr << std::endl << std::endl
				<< "----------------------------------------------------"
				<< std::endl;
		std::cerr << "Unknown exception!" << std::endl << "Aborting!"
				<< std::endl
				<< "----------------------------------------------------"
				<< std::endl;
		return 1;
	}

	return 0;
}

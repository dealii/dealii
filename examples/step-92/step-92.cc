//#define NO_FULLY_DISTRIBUTED_TRIA	// do not use parallel::fullydistributed (instead: parallel::distributed)

/* ---------------------------------------------------------------------
 *
 * This work is partially based in deal.II tutorial step-81 (as on 2023-12-20)
 * Authors of deal.II tutorial step-81 are:
 *     Manaswinee Bezbaruah, Matthias Maier, Texas A&M University, 2021.
 *
 *
 * modifications/enhancements: Copyright (C) 2023 - 2024 by Dr.-Ing. Stephan Voss,
 *     Neunkirchen am Brand, Germany, stvoss@gmx.de
 *
 * ---------------------------------------------------------------------
 *
 * This is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found at
 * https://www.gnu.org/licenses/old-licenses/lgpl-2.1
 * https://www.gnu.org/licenses/licenses/lgpl-3.0
 *
 * ---------------------------------------------------------------------
 *
 * S. Voss' solver modifications / enhancements:
 * - complex valued (frequency domain) solvers for both, E- and H- field (see S1)...S3))
 * - parallelization based on MPI
 * - adopted parameter file
 * - mesh-file import (.msh / usually generated using gmsh - see TO DO 2 for different model/mesh generation)
 * - data post-processors for output on cells, faces and also
 *
 *
 * S1) 'ELECTRIC' solver (maintained/running):
 *
 * "SV_ELECTRIC" delivers an implementation of Poisson's PDE
 *    -div (kappa + i omega ) grad phi_e = div (J + i omega rho_f)
 * for
 *  - solved for complex-valued 3D electric potential phi_e, hence:
 *  - electric field (steady state) problems,
 *  - low-frequency, time-harmonic applications,
 *  - importing 3D mesh-files (*.msh4), e.g. meshed/created using gmsh.
 *
 * This implementation (S1) is actually NEGLECTING:
 *  - skin effect inside conductors,
 *  - wave propagation:
 *  - effects related to relativity
 *
 *  see also TO DO s at documentation of compiler flags and throughout the code
 *
 * ---------------------------------------------------------------------*/

/* ---------------------------------------------------------------------
 *
 * using metric units!
 *
 * symbols:
 *
 * i, j		: []			imaginary number i^2 = -1 = j^2
 *
 * t		: [s]			time
 * f,
 * freq,
 * frequency: [Hz]			frequency
 * omega	: [1/s]			angular frequency (omega := 2 pi f)
 *
 * phi_e	: [V] 			scalar electric potential (real or complex)
 * E		: [V / m]		electric field strength (real or complex)
 * D		: [A s / m^2]	electric flux density (real or complex)
 *
 * rho_f	: [A s / m^3] 	electric charge density (in a volume) (real or complex)
 * J		: [A / m^2] 	current density (real or complex)
 *
 * mu		: [V s / (A m)] magnetic permeability (index '_0' or '_vacuum' means vacuum permeability)
 * epsilon  : [A s / (V m)] electric permittivity (index '_0' or '_vacuum' means vacuum permittivity)
 * kappa	: [A / (V m)] 	electric conductivity
 *
 * alpha	: [1/m]			skin constant, alpha^2 := i omega kappa mu ; [1/m^2]
 * beta		: [1/m]			wave number, beta^2 := omega^2 epsilon mu = -(i omega)^2 epsilon mu ; [1/m^2]
 * gamma	: [1/m]			complex wave number: gamma^2 := alpha^2 - beta^2 ; [1/m^2]
 *
 *
 * x,y,z	: [m], [m], [m]			cartesian coordinates
 * rho,phi,z: [m], [radians], [m]	cylinder coordinates
 *
 *
 * usage of differential operators (div, grad, rot, curl, nabla) on scalar or vector fields :
 *
 * 		rot vector_F = nabla x vector_F = curl vector_F  (rot operator is same as curl operator)
 * 		div vector_F = nabla * vector_F
 * 		grad scalar_f = nabla * scalar_f
 *
 * ---------------------------------------------------------------------*/

/* ---------------------------------------------------------------------
 *
 * deal.II: Copyright (C) 2021 - 2023 by the deal.II authors
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 */

// @sect3{Include files}

// The fe/fe_nedelec_sz.h file allows to use the FE_NedelecSZ elements.
// It is an implementation of the $H^{curl}$ conforming Nédélec Elements
// (for details refer to the documentation of the FE_NedelecSZ element).

#include <deal.II/base/function.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/generic_linear_algebra.h>


// This program can use either PETSc or Trilinos for its parallel
// algebra needs. By default, if deal.II has been configured with
// PETSc, it will use PETSc. Otherwise, the following few lines will
// check that deal.II has been configured with Trilinos and take that.
//
// But there may be cases where you want to use Trilinos, even though
// deal.II has *also* been configured with PETSc, for example to
// compare the performance of these two libraries. To do this,
// add the following \#define to the source code:
// @code
// #define FORCE_USE_OF_TRILINOS
// @endcode
//
// Using this logic, the following lines will then import either the
// PETSc or Trilinos wrappers into the namespace `LA` (for linear
// algebra). In the former case, we are also defining the macro
// `USE_PETSC_LA` so that we can detect if we are using PETSc (see
// solve() for an example where this is necessary).
namespace LA
{
#if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \
  !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
  using namespace dealii::LinearAlgebraPETSc;
#  define USE_PETSC_LA
#elif defined(DEAL_II_WITH_TRILINOS)
  using namespace dealii::LinearAlgebraTrilinos;
#else
#  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required
#endif
} // namespace LA

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>


#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>


#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_faces.h>
#include <deal.II/numerics/error_estimator.h>


// Utilities::System namespace will be used to query things like the
// number of processors associated with the current MPI universe, or the
// number within this universe the processor this job runs on is:
#include <deal.II/base/utilities.h>

// class ConditionOStream allows to write
// code that would output things to a stream (such as <code>std::cout</code>)
// on every processor but throws the text away on all but one of them. 
#include <deal.II/base/conditional_ostream.h>

// Indicate which elements a particular
// processor owns or need to know of. This is the realm of the IndexSet class:
// if there are a total of $N$ cells, degrees of freedom, or vector elements,
// associated with (non-negative) integral indices $[0,N)$, then both the set
// of elements the current processor owns as well as the (possibly larger) set
// of indices it needs to know about are subsets of the set $[0,N)$. IndexSet
// is a class that stores subsets of this set in an efficient format:
#include <deal.II/base/index_set.h>

// The next header file is necessary for a single function,
// SparsityTools::distribute_sparsity_pattern:
#include <deal.II/lac/sparsity_tools.h>

// parallel::distributed::Triangulation provide meshes distributed
// across a potentially very large number of processors, while the second
// provides the namespace parallel::distributed::GridRefinement that offers
// functions that can adaptively refine such distributed meshes:
#ifndef NO_FULLY_DISTRIBUTED_TRIA
#include <deal.II/distributed/fully_distributed_tria.h>
#else
#include <deal.II/distributed/tria.h>
#endif
#include <deal.II/distributed/grid_refinement.h>


#include <fstream>
#include <iostream>
#include <memory>


// @sect3{Class Template Declarations}
//
// Declaration of all classes with their data structures and methods

namespace SV_Maxwell
{
  template <typename T> T sum (const T &, const T &);

  template <typename T> double sum (const double &a, const double &b) { return a+b; };

  using namespace dealii;
  using namespace std::complex_literals;

  // @sect4{Parameters Class}

  // The Parameters class inherits ParameterAcceptor, and instantiates all the
  // coefficients in the variational equations.
  // These coefficients are passed through ParameterAcceptor and are editable
  // through a .prm file.



  template <int dim>
  class ProblemParameters : public ParameterAcceptor
  {
  public:
    ProblemParameters();
    ~ProblemParameters();
    void finalize();

    using absphase_type = Tensor<1, 2, double>;

    using rank0_type = std::complex<double>;

    using rank1_type = Tensor<1, dim, std::complex<double>>;

    using phys_ID = unsigned int;

    using rank2_type = Tensor<2, dim, rank0_type>;

    using curl_type = Tensor<1, dim == 2 ? 1 : dim, rank0_type>;

    using Material_tuple = std::tuple<std::string,unsigned int,std::complex<double>,std::complex<double>,double>;
    using FaceOperationalData_tuple = std::tuple<std::string,phys_ID,std::complex<double> >;


#define SV_NAN	1e40

    class SV_physical {
        public:
        types::material_id   material_id;
        std::string			 name;
        std::complex<double> epsilon;		// absolute epsilon (not relative !! --> eps_r * eps_0)
        std::complex<double> mu;			// absolute mu (not relative !! --> mu_r * mu_0)
        double               kappa;

        std::complex<double> alpha_sq, beta_sq, gamma_sq;	// absolute values for skin (alpha) and wave constants (beta, gamma resp.)!
    };

    class SV_operational_data {
        public:
        types::material_id   material_id;
        std::string			 name;
        std::complex<double> potential; // [V] electric potential
    };


  public:

    double get_omega();
    double get_omega_sq();

    double get_sigma(const Point<dim> & x,
                     types::material_id left,
                     types::material_id right);

    //rank1_type get_J_a(const Point<dim> &point, types::material_id material);

    SV_physical *get_p_physical(types::material_id material);
    SV_physical get_physical(types::material_id material);

    SV_operational_data *get_p_operational_data(types::material_id material);
    SV_operational_data get_operational_data(types::material_id material);

    std::string get_mesh_filename();
    std::string get_cells_solution_filename();

    double         mu_vacuum;		// vacuum permeability: 4 pi 10e-7 V s / (A m)
    double         epsilon_vacuum;  // vacuum permittivity: 8.8541878188(14)×10^−12 A s / (V m)

    std::string mesh_filename,
				cells_solution_filename;

  private:
    double         frequency;		// [Hz] frequency
    double         omega;			// [1/s] angular frequency omega = 2 pi frequency
    double         omega_sq;		// = omega^2

    int				default_physical_ix=-1;

    std::list<Material_tuple> material_list;
    std::list<FaceOperationalData_tuple> operational_data_list;

    /* TO DO: scaling of material properties and domains has not been maintained and may be faulty/buggy !!
     * --> better remove this !!
     */
    #define SV_SCALE_OMEGA		1.0
	#define SV_SCALE_KAPPA		1.0 //1.0e7
	#define SV_SCALE_MU			1.0 //mu_vacuum
	#define SV_SCALE_EPSILON	1.0 //epsilon_vacuum


  public:
    rank1_type r1null;

    unsigned int	N_components;	// 1 for real (omega==0), 2 for complex (omega!=0)
    bool			is_complex;
    unsigned long N_physicals=0;
    #define SV_N_max_physicals 20
    SV_physical *pphysicals = NULL;
    unsigned long N_operational_data=0;
    SV_operational_data *p_operational_data = NULL;

  };


  template <int dim>
  ProblemParameters<dim>::ProblemParameters()
    : ParameterAcceptor("Problem")
  {
	mu_vacuum = 4.0e-7*dealii::numbers::PI; // [ V s / (A m) ]
	epsilon_vacuum = 8.8541878176204199e-12; // [ A s / (V m) ]


	frequency = 50.0;			// [Hz]
	add_parameter("frequency",
				  frequency,
				  "frequency");

	mesh_filename = "";
	add_parameter("mesh filename", mesh_filename, "mesh file name");
	cells_solution_filename = "";
	add_parameter("cells solution filename", cells_solution_filename, "filename for exporting cells solution");

    add_parameter("materials", material_list, "list of material definitions: \"name : ID : eps_r_r,eps_r_i : mu_r_r,mu_r_i : kappa_abs\"");
	add_parameter("potentials", operational_data_list, "list of potential definitions: \"name : ID : epot_r, epot_i [V, abs.]\"");

  }

  template <int dim>
  void ProblemParameters<dim>::finalize()
  {
	const std::complex<double> cnull(0,0);
	for (unsigned int ix=0;ix<dim;++ix)
		r1null[ix] = cnull;

    omega = 2.0*dealii::numbers::PI*frequency;
    omega_sq = omega*omega;
    N_components = ( omega==0 ? 1 : 2 );
    is_complex = (N_components>1);

    constexpr std::complex<double> real(1,0),imag(0,1);

    N_physicals = material_list.size()+1;

    pphysicals = new SV_physical[N_physicals];

    int ix = 0;
    default_physical_ix=-1;
    for (auto const& item : material_list) {
	    std::string name =  std::get<0>(item);
	    unsigned int id =  std::get<1>(item);
	    std::complex<double> eps_r =  std::get<2>(item);
	    std::complex<double> mu_r =  std::get<3>(item);
	    double kappa =  std::get<4>(item);

    	pphysicals[ix].name 		  = name;

    	pphysicals[ix].material_id    = id;
    	if (id==0)
    		default_physical_ix = ix;

    	pphysicals[ix].epsilon        = eps_r * epsilon_vacuum;
        pphysicals[ix].mu             = mu_r * mu_vacuum;
        pphysicals[ix].kappa          = kappa;

        pphysicals[ix].alpha_sq       = imag*omega*pphysicals[ix].mu*pphysicals[ix].kappa;
        pphysicals[ix].beta_sq        = real*omega_sq*pphysicals[ix].mu*pphysicals[ix].epsilon;
        pphysicals[ix].gamma_sq        = pphysicals[ix].alpha_sq - pphysicals[ix].beta_sq;

	    ix ++;
    }

    N_physicals = ix;

    if (default_physical_ix<0) {
    	N_physicals++;
    	pphysicals[ix].name 		  = "default";

    	default_physical_ix = ix;
    	pphysicals[ix].material_id    = 0;

    	pphysicals[ix].epsilon        = epsilon_vacuum;
        pphysicals[ix].mu             = mu_vacuum;
        pphysicals[ix].kappa          = 0.0;

        pphysicals[ix].alpha_sq       = imag*omega*pphysicals[ix].mu*pphysicals[ix].kappa;
        pphysicals[ix].beta_sq        = real*omega_sq*pphysicals[ix].mu*pphysicals[ix].epsilon;
        pphysicals[ix].gamma_sq       = pphysicals[ix].alpha_sq - pphysicals[ix].beta_sq;
    }

    N_operational_data = operational_data_list.size();
    p_operational_data = new SV_operational_data[N_operational_data];

    unsigned int ipot = 0;
    SV_operational_data *t_p_opdata;
    for (auto const& item : operational_data_list) {
	    std::string name =  std::get<0>(item);
	    auto id =  std::get<1>(item);
	    std::complex<double> potential =  std::get<2>(item);

	    t_p_opdata = &p_operational_data[ipot];
	    t_p_opdata->material_id = id;
	    t_p_opdata->name = name;
	    t_p_opdata->potential = potential;

	    ipot++;
    }

  }

  template <int dim>
  ProblemParameters<dim>::~ProblemParameters()
  {
  }

  template <int dim>
  std::string ProblemParameters<dim>::get_mesh_filename()
  {
    return mesh_filename;
  }

  template <int dim>
  std::string ProblemParameters<dim>::get_cells_solution_filename()
  {
    return cells_solution_filename;
  }

  template <int dim>
  double ProblemParameters<dim>::get_omega()
  {
    return omega;
  }

  template <int dim>
  double ProblemParameters<dim>::get_omega_sq()
  {
    return omega_sq;
  }

  template <int dim>
  typename ProblemParameters<dim>::SV_physical *
  ProblemParameters<dim>::get_p_physical(types::material_id material)
  {
    for (unsigned int ix=0;ix<N_physicals;ix++) {
        if (material==(pphysicals[ix].material_id)) {
            return &pphysicals[ix];
        }
    }

    return NULL;
  }

  template <int dim>
  typename ProblemParameters<dim>::SV_physical
  ProblemParameters<dim>::get_physical(types::material_id material)
  {
    SV_physical *retval_physical=get_p_physical(material);

    if (retval_physical==NULL)
    	retval_physical = &pphysicals[default_physical_ix];

    return *retval_physical;
  }

  template <int dim>
  typename ProblemParameters<dim>::SV_operational_data *
  ProblemParameters<dim>::get_p_operational_data(types::material_id material)
  {
    for (unsigned int ix=0;ix<N_operational_data;ix++) {
        if (material==(p_operational_data[ix].material_id)) {
            return &p_operational_data[ix];
        }
    }

    return NULL;
  }

  template <int dim>
  typename ProblemParameters<dim>::SV_operational_data
  ProblemParameters<dim>::get_operational_data(types::material_id material)
  {
	  SV_operational_data *retval_operational_data=get_p_operational_data(material);

    if (retval_operational_data==NULL)
    	retval_operational_data = &p_operational_data[0];

    return *retval_operational_data;
  }

  // @sect4{DirichletBoundaryValues Class}

  // This class implements functions to set boundary conditions

  template <int dim>
  class DirichletBoundaryValues : public Function<dim>
  {
  public:
	unsigned int N_components;
	double *m_values;


	DirichletBoundaryValues(unsigned int _N_components, // e.g.: 1 for real numbers, 2 for complex numbers
							double *pvalues=NULL)		// pointer to array of _N_components values to set all components
      : Function<dim>(_N_components)
	  , N_components(_N_components)
    {
    	m_values = new double[_N_components];

        for (unsigned int ix=0;ix<N_components;ix++)
        	m_values[ix] = pvalues!=NULL ? pvalues[ix] : 0.0;
    }


    ~DirichletBoundaryValues()
    {
    	delete m_values;
    }


    // return (component of) function value at a point:
    virtual double 	value (const Point< dim > &, const unsigned int _component) const override
    {
    	return(m_values[_component]);
    }

    // return gradient of function value at a point:
    virtual Tensor< 1, dim, double> gradient (const Point< dim > &, const unsigned int) const override
    {
    	return(Tensor<1,dim,double>({0,0,0}));
    }

    // set function vector at point
    virtual void vector_value(const Point<dim> & /*p*/,
                              Vector<double> &values) const override
    {
      AssertDimension(values.size(), N_components);

      for (unsigned int ix=0;ix<N_components;ix++) {
    	  values(ix) = m_values[ix];
      }
    }

    // set entire value list
    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  value_list) const override
    {
      AssertDimension(value_list.size(), points.size());

      for (unsigned int p = 0; p < points.size(); ++p)
        DirichletBoundaryValues<dim>::vector_value(points[p], value_list[p]);
    }
  };




  // @sect4{Maxwell_Poisson Class}
  // Declaration of all major building blocks of
  // finite element program which consists of the usual setup and
  // assembly routines.

  template <int dim>
  class Maxwell_Poisson : public ParameterAcceptor
  {
  public:
    Maxwell_Poisson();
    ~Maxwell_Poisson();
    void run();

    using rank1_type = Tensor<1, dim, std::complex<double>>;

    ProblemParameters<dim>      problem_parameters;

  private:
    MPI_Comm mpi_communicator;

    /* run time parameters */
    unsigned int N_refinements,N_materials;
    unsigned int fe_E_order;
    unsigned int quadrature_order;
    unsigned int N_BC_func_components;

    void parse_parameters_callback();
    void make_grid();
    void setup_system();
    void assemble_system();
    void solve();
    void refine_grid();
    void output_results();

    ConditionalOStream pcout;
    TimerOutput        computing_timer;

#ifndef NO_FULLY_DISTRIBUTED_TRIA
    parallel::fullydistributed::Triangulation<dim> triangulation;
#else
    parallel::distributed::Triangulation<dim> triangulation;
#endif

    std::unique_ptr<FiniteElement<dim>> fe;
    DoFHandler<dim> dof_handler;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    AffineConstraints<double> constraints;

    SparsityPattern           sparsity_pattern;
    LA::MPI::SparseMatrix     system_matrix;
    LA::MPI::Vector           locally_relevant_solution;
    LA::MPI::Vector           system_rhs;

    class CellsPostprocessor;
    class FacesPostprocessor;
  };

  // @sect3{Class Template Definitions and Implementation}
  //
  // @sect4{The Constructor}
  // The Constructor simply consists of default initialization a number of
  // discretization parameters (such as the domain size, mesh refinement,
  // and the order of finite elements and quadrature) and declaring a
  // corresponding entry via ParameterAcceptor::add_parameter(). All of
  // these can be modified by editing the .prm file. Absorbing boundary
  // conditions can be controlled with the absorbing_boundary boolean. If
  // absorbing boundary conditions are disabled, simply
  // homogeneous Dirichlet conditions on the tangential component of the
  // magnetic field are enforced.

  template <int dim>
  Maxwell_Poisson<dim>::~Maxwell_Poisson()
  {
      //if (locally_relevant_pots!=NULL) delete locally_relevant_pots;
      //if (locally_relevant_curls!=NULL) delete locally_relevant_curls;
  }


  template <int dim>
  Maxwell_Poisson<dim>::Maxwell_Poisson()
    : ParameterAcceptor(
    		"Maxwell_E_solver"
    		)
    , mpi_communicator(MPI_COMM_WORLD)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    , computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::never,
                      TimerOutput::wall_times)
    , triangulation(mpi_communicator)
    , dof_handler(triangulation)
  {

    ParameterAcceptor::parse_parameters_call_back.connect(
      [&]() { parse_parameters_callback(); });

    N_refinements = 0;
    add_parameter("refinements",
                  N_refinements,
                  "number of refinements of the geometry");

    N_materials = 0;
    add_parameter("materials",
                  N_materials,
                  "number of materials");

    fe_E_order = 0;
    add_parameter("fe E order", fe_E_order, "order of the finite element space for scalar electric potential");


    quadrature_order = 1;
    add_parameter("quadrature order",
                  quadrature_order,
                  "order of the quadrature");
  }


  template <int dim>
  void Maxwell_Poisson<dim>::parse_parameters_callback()
  {
	  problem_parameters.finalize();

	  N_BC_func_components = 0;

	  pcout << "fe_E_order = " << fe_E_order << std::endl;
	  N_BC_func_components += problem_parameters.N_components;
	  pcout << "quadrature order = " << quadrature_order << std::endl;

	  fe = std::make_unique<FESystem<dim>>(FE_Q<dim>(fe_E_order) ^ problem_parameters.N_components);

  }

  // The Maxwell::make_grid() routine
  // reads the mesh-file for the computational domain.
  // A block decomposition into real and imaginary matrices
  // for the solution matrices is used.

  template <int dim>
  void Maxwell_Poisson<dim>::make_grid()
  {
    TimerOutput::Scope t(computing_timer, "01. make grid");

    std::string infilename = problem_parameters.get_mesh_filename();
    infilename += ".msh";

    GridIn<dim> gi;

#ifndef NO_FULLY_DISTRIBUTED_TRIA
    const unsigned int mpi_size =
      Utilities::MPI::n_mpi_processes(mpi_communicator);
    auto construction_data = TriangulationDescription::Utilities::
      create_description_from_triangulation_in_groups<dim, dim>(
        [&](Triangulation<dim> &tria) {
    	    pcout
              << "MPI-comm: " << Utilities::MPI::this_mpi_process(mpi_communicator)
              << ", file: \"" << infilename << "\"" << std::endl;
            gi.attach_triangulation(tria);
            gi.read_msh(infilename);
        },
        [&](Triangulation<dim> &tria_serial,
            const MPI_Comm /*mpi_comm*/,
            const unsigned int /*group_size*/) {
          GridTools::partition_triangulation(mpi_size, tria_serial);
        },
        mpi_communicator,
        1);
    triangulation.create_triangulation(construction_data);
#else //  #ifndef NO_FULLY_DISTRIBUTED_TRIA
    Triangulation<dim> tria;
	gi.attach_triangulation(tria);
	gi.read_msh(infilename);

	triangulation.copy_triangulation(tria);
#endif
  }

  // The Maxwell::setup_system() routine follows the usual routine of
  // enumerating all the degrees of freedom and setting up the matrix and
  // vector objects to hold the system data. Enumerating is done by using
  // DoFHandler::distribute_dofs().

  template <int dim>
  void Maxwell_Poisson<dim>::setup_system()
  {
    TimerOutput::Scope t(computing_timer, "02. setup");

    dof_handler.distribute_dofs(*fe);
    pcout << "  number of fe_systems: " << fe->n_components() << std::endl
        	<< "  setup: number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;


    // Two index sets that provide information about which degrees of
    // freedom are owned by the current processor and an index set that
    // indicates which degrees of freedom are locally relevant.
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);

    // Initialize the solution and right hand side vectors.
    // Solution vector we seek does not only store
    // elements locally owned, but also ghost entries; 
    // Right hand side vector only needs to have the entries the current processor
    // owns (It will only be written locally- never read.) 
    // (Linear solvers will read from it, but they do not care about 
    // the geometric location of degrees of freedom).
    locally_relevant_solution.reinit(locally_owned_dofs,
                                     locally_relevant_dofs,
                                     mpi_communicator);
    system_rhs.reinit(locally_owned_dofs, mpi_communicator);


    // Compute hanging node and boundary value
    // constraints, which are combined into a single object storing all
    // constraints.
    //
    // As with all other things in %parallel, the mantra must be that no
    // processor can store all information about the entire universe. As a
    // consequence, the AffineConstraints object needs information for which
    // degrees of freedom it can store constraints and for which it may not
    // expect any information to store. 
    // The degrees of freedom it needs to care about on
    // each processor are the locally relevant ones, so this is passed to the
    // AffineConstraints::reinit function. 
    
    // As a side note, if it's forgotten to
    // pass this argument, the AffineConstraints class will allocate an array
    // with length equal to the largest DoF index it has seen so far. For
    // processors with high MPI process number, this may be very large --
    // maybe on the order of billions. The program would then allocate more
    // memory than for likely all other operations combined for this single
    // array.
    constraints.clear();
    // seen in tutorial step-41: constraints.reinit(locally_relevant_dofs);

    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

#ifndef kdshfkahfs// not defined(DEBUG)

    const FEValuesExtractors::Scalar 	phi_E_re(                         0     ),
    									phi_E_im((problem_parameters.is_complex ? 1 : 0));

    for (unsigned int ix=0;ix<problem_parameters.N_operational_data;ix++)
	{
		pcout << "set physical " << problem_parameters.p_operational_data[ix].material_id;

		const auto potential = problem_parameters.p_operational_data[ix].potential;
		const auto pot_real = potential.real(), pot_imag = potential.imag();
		if (pot_real<SV_NAN and -pot_real<SV_NAN and pot_imag<SV_NAN and -pot_imag<SV_NAN)
		{
			pcout << " potential " << problem_parameters.p_operational_data[ix].potential << " V";

			VectorTools::interpolate_boundary_values(
			  dof_handler, problem_parameters.p_operational_data[ix].material_id,
			  Functions::ConstantFunction<dim,double>(pot_real,N_BC_func_components),
			  constraints,
			  fe->component_mask(phi_E_re));

			if (problem_parameters.is_complex)
			{
				VectorTools::interpolate_boundary_values(
				  dof_handler, problem_parameters.p_operational_data[ix].material_id,
				  Functions::ConstantFunction<dim,double>(pot_imag,N_BC_func_components),
				  constraints,
				  fe->component_mask(phi_E_im));
			}
		}


		pcout << std::endl;

	}
#endif // _NOT_IN_USE_


    constraints.close();

    DynamicSparsityPattern dsp(locally_relevant_dofs);

    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern(dsp,
                                               dof_handler.locally_owned_dofs(),
                                               mpi_communicator,
                                               locally_relevant_dofs);

    system_matrix.reinit(locally_owned_dofs,
                         locally_owned_dofs,
                         dsp,
                         mpi_communicator);

    pcout << "  system set up." << std::endl;
  }

  // This is a helper function that takes the tangential component of a tensor.
  template <int dim>
  DEAL_II_ALWAYS_INLINE inline Tensor<1, dim, std::complex<double>>
  tangential_part(const Tensor<1, dim, std::complex<double>> &tensor,
                  const Tensor<1, dim> &                      normal)
  {
    auto result = tensor;
    result[0]   = normal[1] * (tensor[0] * normal[1] - tensor[1] * normal[0]);
    result[1]   = -normal[0] * (tensor[0] * normal[1] - tensor[1] * normal[0]);
    return result;
  }

  // Assemble the stiffness matrix and the right-hand side:
  template <int dim>
  void Maxwell_Poisson<dim>::assemble_system()
  {
    TimerOutput::Scope t(computing_timer, "03. assemble");

    QGauss<dim>     quadrature_formula(quadrature_order);
    QGauss<dim - 1> face_quadrature_formula(quadrature_order);

    FEValues<dim>     fe_values(*fe,
                                 quadrature_formula,
                                 update_values | update_gradients |
                                 update_quadrature_points |
                                 update_JxW_values);
                                 
    FEFaceValues<dim> fe_face_values(*fe,
                                     face_quadrature_formula,
                                     update_values | update_gradients |
                                     update_quadrature_points |
                                     update_normal_vectors |
                                     update_JxW_values);

    const auto dofs_per_cell = fe->dofs_per_cell;

    const unsigned int n_q_points      = quadrature_formula.size();
    [[maybe_unused]] const unsigned int n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // This is assembling the interior of the domain on the left hand side.
    // In doing so, test functions $\varphi_i$ and $\varphi_j$ are needed, and the
    // curl of these test variables.

    const auto omega = problem_parameters.get_omega();
    const std::complex<double> real{1.0,0.0},imag{0., problem_parameters.is_complex ? 1.0 : 0.0};
    constexpr std::complex<double> cnull{0., 0.};

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
          if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);

            const FEValuesExtractors::Scalar phi_E_real(0),
            								 phi_E_imag((problem_parameters.is_complex ? 1 : 0));

            cell_matrix = 0.;
            cell_rhs    = 0.;

            cell->get_dof_indices(local_dof_indices);
            const auto material_id = cell->material_id();

            const auto   material = problem_parameters.get_physical(material_id);

            std::complex<double> psi_E{	/* real */ omega==0 ? (material.kappa==0 ? material.epsilon.real() : material.kappa) : material.kappa-omega*material.epsilon.imag(),
            							/* imag */ omega==0 ? 0 : omega*material.epsilon.real()};

            for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
              {

                for(const auto i : fe_values.dof_indices())
                  {
                    const auto grad_phi_E_i_r =        fe_values[phi_E_real].gradient(i, q_point),
                               grad_phi_E_i_i = imag * fe_values[phi_E_imag].gradient(i, q_point);

                    const auto grad_phi_E_i_conj      = grad_phi_E_i_r - grad_phi_E_i_i;

                    for (const auto j : fe_values.dof_indices())
                      {
                        auto temp = cnull;
                        const auto grad_phi_E_j_r =        fe_values[phi_E_real].gradient(j, q_point);
						const auto grad_phi_E_j_i = imag * fe_values[phi_E_imag].gradient(j, q_point);
                        const auto grad_phi_E_j   = grad_phi_E_j_r + grad_phi_E_j_i;

                        // el. PDE: psi_E * nabla phi_E_j * nabla phi_E_i'
                        /* E.2 */
						temp  += psi_E * scalar_product( grad_phi_E_j, grad_phi_E_i_conj ) * fe_values.JxW(q_point);

                        cell_matrix(i, j) += temp.real();
                      } // for j
                  } // for i

              } // for qpoint



/*------------------------------------------------------------------------------------------
 *
 * FACES integral
 *
 *------------------------------------------------------------------------------------------*/



            // Assemble the face and the boundary.
            /*
             * ....this section is not deleted for partial reusage while implementing following TO DO:
             *
             * TO DO: for magnetic field an integral boundary condition will be needed for connecting current sources to windings
             * This connection will be done on boundaries (with boundary IDs) and integrating vector potential: A = J / (i omega kappa)
             */


            if constexpr(dim>2)
            {

                for (const auto &face : cell->face_iterators())
				  {
					if (face->at_boundary())
					  {

						fe_face_values.reinit(cell, face);

			            std::complex<double> currdens_J(0,0);

			            const FEValuesExtractors::Scalar phi_E_real(                         0     ),
			            								 phi_E_imag((problem_parameters.is_complex ? 1 : 0));
						for (unsigned int q_point = 0; q_point < n_face_q_points; q_point++)
						  {
							const auto normal = fe_face_values.normal_vector(q_point);

							for (const auto i : fe_face_values.dof_indices())
							  {

							    const auto 	phi_E_i_r =        fe_face_values[phi_E_real].value(i, q_point),
											phi_E_i_i = imag * fe_face_values[phi_E_imag].value(i, q_point);

								const auto phi_E_i_conj = phi_E_i_r - phi_E_i_i;

								for (const auto j : fe_face_values.dof_indices())
								  {
									std::complex<double> temp(0,0);

									{
										const auto 	grad_phi_E_j_r =        fe_face_values[phi_E_real].gradient(j, q_point),
													grad_phi_E_j_i = imag * fe_face_values[phi_E_imag].gradient(j, q_point);

										const auto grad_phi_E_j = grad_phi_E_j_r + grad_phi_E_j_i;

										// el. PDE: -psi_E (nabla phi_E_j) * phi_E_i' * e_n
										/* E.3 */
										temp -= (psi_E /*-material.kappa*/) * phi_E_i_conj * scalar_product(grad_phi_E_j,normal) * fe_face_values.JxW(q_point);
									}

									cell_matrix(i, j) += temp.real();
								  } /* for(j) */

							  }     /* for(i) */

						  }         /* for(q_point) */
					  } // if (face->at_boundary())
				  } // END: for (const auto &face : cell->face_iterators())
            } // END: if (dim>2)




            constraints.distribute_local_to_global(
            cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);


          } // END: if (cell->is_locally_owned())

    } // END: for (const auto &cell : dof_handler.active_cell_iterators())


    // In the operations above, specifically the call to
    // `distribute_local_to_global()` in the last line, every MPI
    // process was only working on its local data. If the operation
    // required adding something to a matrix or vector entry that is
    // not actually stored on the current process, then the matrix or
    // vector object keeps track of this for a later data exchange,
    // but for efficiency reasons, this part of the operation is only
    // queued up, rather than executed right away. But now that we got
    // here, it is time to send these queued-up additions to those
    // processes that actually own these matrix or vector entries. 
    
    // In other words: This is "finalization" of the global data
    // structures. This is done by invoking the function `compress()`
    // on both the matrix and vector objects. See
    // @ref GlossCompress "Compressing distributed objects"
    // for more information on what `compress()` actually does.
    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
    
    pcout << "  system assembled." << std::endl;
  }

  // Use a direct solver to solve the system:
  // SparseDirectMUMPS for MPI parallel implementation
  template <int dim>
  void Maxwell_Poisson<dim>::solve()
  {
    TimerOutput::Scope t(computing_timer, "04. solve");

    LA::MPI::Vector completely_distributed_solution(locally_owned_dofs,
                                                    mpi_communicator);

    SolverControl solver_control;//(dof_handler.n_dofs(), 1e-12);
    dealii::PETScWrappers::SparseDirectMUMPS solver(solver_control);//, mpi_communicator);
    solver.set_symmetric_mode(true);
    solver.solve(system_matrix, completely_distributed_solution, system_rhs);
        
    auto nsteps = solver_control.last_step();
    pcout << "  solved in " << nsteps << " iteration" << (nsteps>1 ? "s." : ".")
          << std::endl;

    constraints.distribute(completely_distributed_solution);

    locally_relevant_solution = completely_distributed_solution;

  }


  // @sect4{Maxwell::refine_grid}

  // The function that flags cells to be refined is in namespace
  // parallel::distributed::GridRefinement -- a namespace that has functions
  // that can communicate between all involved processors and determine global
  // thresholds to use in deciding which cells to refine and which to coarsen.
  //
  // Note that nothing special to be done about the
  // KellyErrorEstimator class: just give it a vector with as many elements
  // as the local triangulation has cells (locally owned cells, ghost cells,
  // and artificial ones), but it only fills those entries that correspond to
  // cells that are locally owned.
  template <int dim>
  void Maxwell_Poisson<dim>::refine_grid()
  {
    TimerOutput::Scope t(computing_timer, "refine");

    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
    KellyErrorEstimator<dim>::estimate(
      dof_handler,
      QGauss<dim - 1>(fe->degree + 1),
      std::map<types::boundary_id, const Function<dim> *>(),
      locally_relevant_solution,
      estimated_error_per_cell);
    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
      triangulation, estimated_error_per_cell, 0.3, 0.03);
    triangulation.execute_coarsening_and_refinement();

  }




  // @sect4{Maxwell::output_results}

  // Generate output - using a class PostProcessor that
  // inherits from the class DataPostprocessor, which can be attached to
  // DataOut. This allows to output derived quantities from the solution,
  // It overloads the
  // virtual function DataPostprocessor::evaluate_vector_field(),
  // which is then internally called from DataOut::build_patches(). 
  // It's given values of the numerical solution, its derivatives, normals to the
  // cell, the actual evaluation points and any additional quantities.

  template <int dim>
  class Maxwell_Poisson<dim>::CellsPostprocessor
    : public DataPostprocessor<dim>
  {
  public:
    CellsPostprocessor(Maxwell_Poisson<dim> * problem);

	  virtual void evaluate_field(
			const DataPostprocessorInputs::Scalar<dim> *psinputs,
			const DataPostprocessorInputs::Vector<dim> *pvinputs,
		    std::vector<Vector<double>> &computed_quantities) const;


    virtual void evaluate_scalar_field(
      const DataPostprocessorInputs::Scalar<dim> &inputs,
      std::vector<Vector<double>> &computed_quantities) const override;

    virtual void evaluate_vector_field(
      const DataPostprocessorInputs::Vector<dim> &inputs,
      std::vector<Vector<double>> &computed_quantities) const override;

    virtual std::vector<std::string> get_names() const override;

    virtual std::vector<
      DataComponentInterpretation::DataComponentInterpretation>
    get_data_component_interpretation() const override;

    virtual UpdateFlags get_needed_update_flags() const override;

  public:
    Maxwell_Poisson<dim> *problem;
  };


  template <int dim>
  Maxwell_Poisson<dim>::CellsPostprocessor::CellsPostprocessor(
    Maxwell_Poisson<dim> * problem)
    : problem(problem)
  {
  }


  // Define the names for the variables in the output.
  template <int dim>
  std::vector<std::string>
  Maxwell_Poisson<dim>::CellsPostprocessor::get_names() const
  {
    std::vector<std::string> solution_names = { };

    solution_names.emplace_back("Phi_real");
    solution_names.emplace_back("Phi_real2_unused");
    solution_names.emplace_back("Phi_real3_unused");
    solution_names.emplace_back("Phi_imag");
    solution_names.emplace_back("Phi_imag2_unused");
    solution_names.emplace_back("Phi_imag3_unused");

    solution_names.emplace_back("Phi");

    solution_names.emplace_back("Ex_real");
    solution_names.emplace_back("Ey_real");
    solution_names.emplace_back("Ez_real");
    solution_names.emplace_back("Ex_imag");
    solution_names.emplace_back("Ey_imag");
    solution_names.emplace_back("Ez_imag");

    solution_names.emplace_back("E_real");
    solution_names.emplace_back("E_imag");
    solution_names.emplace_back("E");

    solution_names.emplace_back("Dx_real");
    solution_names.emplace_back("Dy_real");
    solution_names.emplace_back("Dz_real");
    solution_names.emplace_back("Dx_imag");
    solution_names.emplace_back("Dy_imag");
    solution_names.emplace_back("Dz_imag");

    solution_names.emplace_back("D_real");
    solution_names.emplace_back("D_imag");
    solution_names.emplace_back("D");

    solution_names.emplace_back("Jx_real");
    solution_names.emplace_back("Jy_real");
    solution_names.emplace_back("Jz_real");
    solution_names.emplace_back("Jx_imag");
    solution_names.emplace_back("Jy_imag");
    solution_names.emplace_back("Jz_imag");

    solution_names.emplace_back("J_real");
    solution_names.emplace_back("J_imag");
    solution_names.emplace_back("J");

    solution_names.emplace_back("material_ID");
    solution_names.emplace_back("mur_real");
    solution_names.emplace_back("mur_imag");
    solution_names.emplace_back("epsilonr_real");
    solution_names.emplace_back("epsilonr_imag");
    solution_names.emplace_back("kappa");


    return solution_names;
  }


  // TO DO: only "SV_ELECTRIC" considered !!
  template <int dim>
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  Maxwell_Poisson<dim>::CellsPostprocessor::get_data_component_interpretation()
    const
  {

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      interpretation;

    interpretation.push_back(DataComponentInterpretation::component_is_scalar); // Phi real
    interpretation.push_back(DataComponentInterpretation::component_is_scalar); // (not in use)
    interpretation.push_back(DataComponentInterpretation::component_is_scalar); // (not in use)
    interpretation.push_back(DataComponentInterpretation::component_is_scalar); // Phi imag
    interpretation.push_back(DataComponentInterpretation::component_is_scalar); // (not in use)
    interpretation.push_back(DataComponentInterpretation::component_is_scalar); // (not in use)

    interpretation.push_back(DataComponentInterpretation::component_is_scalar); // Phi abs

    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector); // E x,y,z, real
    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector); // E x,y,z, imag
    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);

    interpretation.push_back(DataComponentInterpretation::component_is_scalar); // E real
    interpretation.push_back(DataComponentInterpretation::component_is_scalar); // E imag
    interpretation.push_back(DataComponentInterpretation::component_is_scalar); // E abs

    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector); // D x,y,z, real
    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector); // D x,y,z, imag
    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);

    interpretation.push_back(DataComponentInterpretation::component_is_scalar); // D real
    interpretation.push_back(DataComponentInterpretation::component_is_scalar); // D imag
    interpretation.push_back(DataComponentInterpretation::component_is_scalar); // D abs

    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector); // J x,y,z, imag
    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector); // J x,y,z, imag
    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
    interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);

    interpretation.push_back(DataComponentInterpretation::component_is_scalar); // J real
    interpretation.push_back(DataComponentInterpretation::component_is_scalar); // J imag
    interpretation.push_back(DataComponentInterpretation::component_is_scalar); // J abs

    interpretation.push_back(DataComponentInterpretation::component_is_scalar); // mat_ID

    interpretation.push_back(DataComponentInterpretation::component_is_scalar); // mu_r
    interpretation.push_back(DataComponentInterpretation::component_is_scalar);
    interpretation.push_back(DataComponentInterpretation::component_is_scalar); // epsilon_r
    interpretation.push_back(DataComponentInterpretation::component_is_scalar);
    interpretation.push_back(DataComponentInterpretation::component_is_scalar); // kappa

    return interpretation;
  }



  template <int dim>
  UpdateFlags
  Maxwell_Poisson<dim>::CellsPostprocessor::get_needed_update_flags() const
  {
    return update_values | update_gradients | update_quadrature_points
    		;
  }


  // This is the function that computes the derived quantities. 
  template <int dim>
  void Maxwell_Poisson<dim>::CellsPostprocessor::evaluate_scalar_field(
    const DataPostprocessorInputs::Scalar<dim> &inputs,
    std::vector<Vector<double>> &               computed_quantities) const
  {
	  evaluate_field(&inputs,NULL,computed_quantities);
  }

  // This is the function that computes the derived quantities.
  template <int dim>
  void Maxwell_Poisson<dim>::CellsPostprocessor::evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &inputs,
    std::vector<Vector<double>> &               computed_quantities) const
  {
	  evaluate_field(NULL,&inputs,computed_quantities);
  }



  // This is the function that computes the derived quantities.
  template <int dim>
  void Maxwell_Poisson<dim>::CellsPostprocessor::evaluate_field(
	const DataPostprocessorInputs::Scalar<dim> *psinputs,
	const DataPostprocessorInputs::Vector<dim> *pvinputs,
    std::vector<Vector<double>> &computed_quantities) const
  {
	unsigned int n_evaluation_points;
	if (psinputs==NULL)
	{
		n_evaluation_points= pvinputs->solution_values.size();
		Assert(pvinputs->solution_gradients.size() == n_evaluation_points,
			   ExcInternalError());
	}
	else
	{
		n_evaluation_points= psinputs->solution_values.size();
		Assert(psinputs->solution_gradients.size() == n_evaluation_points,
			   ExcInternalError());
	}
    Assert(computed_quantities.size() == n_evaluation_points,
           ExcInternalError());
    //Assert(inputs.solution_values[0].size() == 2*dim, ExcInternalError());

    const bool is_complex = problem->problem_parameters.is_complex;

    const std::complex<double> real(1,0),imag(0,is_complex ? 1 : 0);

    const typename DoFHandler<dim>::cell_iterator current_cell =  psinputs==NULL ? pvinputs->template get_cell<dim>() : psinputs->template get_cell<dim>();

    const auto material_id = current_cell->material_id();

    const auto material = problem->problem_parameters.get_physical(material_id);

    for (unsigned int p = 0; p < n_evaluation_points; p++)
      {
        Maxwell_Poisson::rank1_type sJ;
        double J_real=0,J_imag=0;
        double  E_real=0,E_imag=0,D_real=0,D_imag=0;

        Maxwell_Poisson::rank1_type sE,sD;
        std::complex<double> Phi;


    	if (psinputs==NULL)
    	{
    		const auto t_sol_val = pvinputs->solution_values[p];
    		const auto t_sol_grads = pvinputs->solution_gradients[p];

			Phi = { std::complex<double>(			  t_sol_val(              0          ),
										 is_complex ? t_sol_val((is_complex ? 1 : 0)) : 0) };
			for (unsigned int d = 0; d < dim; d++)
	        {
	          sE[d] = { /* real */  			-t_sol_grads[0                   ][d]  ,
	                    /* imag */ is_complex ? -t_sol_grads[(is_complex ? 1 : 0)][d]  : 0};
	        }
    	}
    	else // if (psinputs==NULL)
    	{
    		const auto t_sol_val = psinputs->solution_values[p];
    		const auto t_sol_grads = psinputs->solution_gradients[p];

			Phi = {std::complex<double>(t_sol_val,0) };

			for (unsigned int d = 0; d < dim; d++)
	        {
	          sE[d] = { /* real */ -t_sol_grads[d] ,
	                    /* imag */  0 };
	        } // for (unsigned int d = 0; d < dim; d++)
    	} // if (psinputs==NULL)

        sD = sE * material.epsilon;

        sJ = sE * material.kappa;

        #define SV_OFF_E     (3+2*dim)
        const unsigned int 	off_E_ePot = 0,
        					off_E_vE   = 7,
							off_E_vD   = off_E_vE  + SV_OFF_E,
							off_E_vJe  = off_E_vD  + SV_OFF_E,
							off_E_mat  = off_E_vJe + SV_OFF_E;
        

        computed_quantities[p](off_E_ePot + 0 ) = Phi.real();
        computed_quantities[p](off_E_ePot + 1 ) = 0;
        computed_quantities[p](off_E_ePot + 2 ) = 0;
        computed_quantities[p](off_E_ePot + 3 ) = Phi.imag();
        computed_quantities[p](off_E_ePot + 4 ) = 0;
        computed_quantities[p](off_E_ePot + 5 ) = 0;
        computed_quantities[p](off_E_ePot + 6 ) = std::sqrt(Phi.real()*Phi.real()+Phi.imag()*Phi.imag());

        for (unsigned int d = 0; d < dim; d++)
        {
          const auto d2 = dim+d;

          computed_quantities[p](off_E_vE + d)  = sE[d].real();
          computed_quantities[p](off_E_vE + d2) = sE[d].imag();

          computed_quantities[p](off_E_vD + d)  = sD[d].real();
          computed_quantities[p](off_E_vD + d2) = sD[d].imag();

          computed_quantities[p](off_E_vJe + d)  = sJ[d].real();
          computed_quantities[p](off_E_vJe + d2) = sJ[d].imag();

          E_real += sE[d].real()*sE[d].real();
          E_imag += sE[d].imag()*sE[d].imag();

          D_real += sD[d].real()*sD[d].real();
          D_imag += sD[d].imag()*sD[d].imag();

          J_real += sJ[d].real()*sJ[d].real();
          J_imag += sJ[d].imag()*sJ[d].imag();

        }

        const auto E_total = std::sqrt(E_real+E_imag);
        E_real = std::sqrt(E_real);
        E_imag = std::sqrt(E_imag);
        const auto D_total = std::sqrt(D_real+D_imag);
        D_real = std::sqrt(D_real);
        D_imag = std::sqrt(D_imag);
        const auto J_total = std::sqrt(J_real+J_imag);
        J_real = std::sqrt(J_real);
        J_imag = std::sqrt(J_imag);

        computed_quantities[p](off_E_vD -3 + 0) = E_real;
        computed_quantities[p](off_E_vD -3 + 1) = E_imag;
        computed_quantities[p](off_E_vD -3 + 2) = E_total;

        computed_quantities[p](off_E_vJe -3 + 0) = D_real;
        computed_quantities[p](off_E_vJe -3 + 1) = D_imag;
        computed_quantities[p](off_E_vJe -3 + 2) = D_total;

        computed_quantities[p](off_E_mat -3 + 0) = J_real;
        computed_quantities[p](off_E_mat -3 + 1) = J_imag;
        computed_quantities[p](off_E_mat -3 + 2) = J_total;
        
        computed_quantities[p](off_E_mat    + 0) = material_id;
        computed_quantities[p](off_E_mat    + 1) = material.mu.real()/problem->problem_parameters.mu_vacuum;
        computed_quantities[p](off_E_mat    + 2) = material.mu.imag()/problem->problem_parameters.mu_vacuum;
        computed_quantities[p](off_E_mat    + 3) = material.epsilon.real()/problem->problem_parameters.epsilon_vacuum;
        computed_quantities[p](off_E_mat    + 4) = material.epsilon.imag()/problem->problem_parameters.epsilon_vacuum;
        computed_quantities[p](off_E_mat    + 5) = material.kappa;

      } // for (unsigned int p = 0; p < n_evaluation_points; p++)

  } // Maxwell_Poisson<dim>::CellsPostprocessor::evaluate_field


  // @sect4{Maxwell::output_results}

  // Write field outputs into files (for cells and faces)
  template <int dim>
  void Maxwell_Poisson<dim>::output_results()
  {
    TimerOutput::Scope t(computing_timer, "05. post-process");

    std::string cells_solution_filename = problem_parameters.get_cells_solution_filename();

    if (cells_solution_filename.size()>0)
    {
        CellsPostprocessor cell_postprocessor(this);//Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));

        cells_solution_filename += ".vtu";

		DataOut<dim> data_out;
		data_out.attach_dof_handler(dof_handler);
		data_out.add_data_vector(locally_relevant_solution, cell_postprocessor);

		data_out.build_patches();
		data_out.write_vtu_in_parallel(cells_solution_filename, mpi_communicator);

		pcout << "exported cells solution into file \"" << cells_solution_filename << "\"" << std::endl;
    }
  }



  // @sect4{Maxwell::run}

  // the central starting-point to run all simulation steps
  template <int dim>
  void Maxwell_Poisson<dim>::run()
  {
    pcout << "Running with "
#ifdef USE_PETSC_LA
          << "PETSc"
#else
          << "Trilinos"
#endif
          << " on " << Utilities::MPI::n_mpi_processes(mpi_communicator)
          << " MPI rank(s)..." << std::endl;

    make_grid();

    for (unsigned int cycle=0;cycle<=N_refinements;cycle++)
    {
    	pcout << "cycle: " << cycle+1 << " / " << N_refinements+1 << std::endl;
    	if (cycle>0)
    	{
    		refine_grid();
    	}

		setup_system();
        assemble_system();
		solve();
    }

    output_results();

    computing_timer.print_summary();
    computing_timer.reset();

    pcout << std::endl;
  }

} // namespace SV_Maxwell


// main function call: set up MPI, Maxwell classes,
// ParameterAcceptor, and call the run() function.

int main(
        int argc, char *argv[]
)
{

  try
    {
      using namespace dealii;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      SV_Maxwell::Maxwell_Poisson<3> maxwell_solver;


      ParameterAcceptor::initialize(argc<2 ? "parameters.prm": argv[1]);

      maxwell_solver.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
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
      std::cerr << std::endl
                << std::endl
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

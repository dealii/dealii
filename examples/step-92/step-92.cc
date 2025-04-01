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

/* ---------------------------------------------------------------------
 *
 * deal.II tutorial step-92
 * 
 * 2024-11, Stephan Voss, Neunkirchen am Brand, Germany, stvoss@gmx.de
 *
 * This solver is derived from several deal.II examples
 * and tutorial steps 1, 3, 49 and 81.
 * Thank You to all deal.II contributors.
 * 
 * ---------------------------------------------------------------------
 */


// @sect4{distributed vs. fullydistributed triangulation}
//
// A small mesh import performance test will be done at the
// end of this tutorial for comparing duration for
// importing meshes when using parallel::distributed
// vs. parallel::fullydistributed triangulation.
// 
// If the option USE_FULLY_DISTRIBUTED_TRIA is defined,
// the solver will be compiled for using
// parallel::fullydistributed triangulation
// otherwise, parallel::distributed triangulation will be used.
#define USE_FULLY_DISTRIBUTED_TRIA


// @sect3{Include files}

#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#ifdef USE_FULLY_DISTRIBUTED_TRIA
#  include <deal.II/distributed/fully_distributed_tria.h>
#else
#  include <deal.II/distributed/tria.h>
#endif

#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>



// @sect3{Class Template Declarations}
//
// Declaration of all classes with their data structures and methods

namespace Step92
{

  using namespace dealii;

  constexpr double epsilon_vacuum = 8.8541878176204199e-12; // vacuum permittivity  [ A s / (V m) ]

  // @sect4{Parameters Class}

  // The Parameters class inherits ParameterAcceptor, and instantiates all the
  // coefficients in the variational equations.
  // These coefficients are passed through ParameterAcceptor and are editable
  // through a .prm file.

  class ProblemParameters : public ParameterAcceptor
  {
  public:
    ProblemParameters();
    ~ProblemParameters() = default;
    void finalize();

    using physical_ID = unsigned int;

    using DomainProperties_tuple = std::tuple<std::string, physical_ID, double, double>;

    class domain_properties
    {
    public:
      double get_domain_parameter_a();
      
      types::material_id    domain_id;      // domain ID number
      std::string           name;           // clear-text name for this domain
      double                epsilon;        // absolute electric permittivity epsilon [ A s / (V m) ]
      double                kappa;          // absolute electric conductivity
    };

    using BoundaryPotential_tuple = std::tuple<std::string, physical_ID, double>;

    class boundary_potential
    {
    public:
      types::material_id boundary_id;   // boundary ID number
      std::string        name;          // clear-text name for this boundary
      double             potential;     // [V] electric potential used as dirichlet boundary value
    };


  public:

    domain_properties get_domain_properties(types::material_id domain_id);

    boundary_potential get_boundary_value(types::material_id boundary_id);

    std::string mesh_filename, cells_output_data_filename;

  private:

    std::list<DomainProperties_tuple>          domains_list;
    std::list<BoundaryPotential_tuple> boundary_values_list;

  public:
    int default_domain_index;
    std::vector<domain_properties> domains_vector;

    const types::material_id ID_NO_BOUNDARY_VALUE = -1;
    const boundary_potential no_boundary_value = { ID_NO_BOUNDARY_VALUE, "none", 0 };
    std::vector<boundary_potential> boundaries_vector;
  };


  double ProblemParameters::domain_properties::get_domain_parameter_a()
  {
    return (kappa>0 ? kappa : epsilon);
  }


  // @sect5{ProblemParameters class constructor}
  // In the constructor ProblemParameters::ProblemParameters(), all parameters are defined,
  // which shall be read when parsing the input parameter file:
  // <ul>
  //   <li>
  //     "mesh filename" is a string specifying a gmsh file (without suffix ".msh") containing the triangulation
  //   </li>
  //   <li>
  //     "cells solution filename" is sepcifying filename (without suffix) for exporting output from cells solution data.
  //     This filename will automatically be extended with a suffix ".vtu"
  //   </li>
  //   <li>
  //     "domains" is a list of physical domains.
  // 
  //     Each comma-separated entry is formatted as "name : ID : epsilon_r : kappa", with
  //     <ol>
  //       <li> "name" : clear-text name of a domain </li>
  //       <li> "ID" is a unique positive ID number (or zero for default domain) </li>
  //       <li> "epsilon_r" is the relative electric permittivity </li>
  //       <li> "kappa" is the electric conductivity [A/Vm] </li>
  //     </ol>
  //   </li>
  //   <li>
  //     "potentials" is a list of boundary values.
  //
  //     Each comma-separated entry is formatted as "name : ID : el. potential", with
  //     <ol>
  //       <li> "name" is a clear-text name of a boundary </li>
  //       <li> "ID" is a unique positive ID number (or zero for default boundary) </li>
  //       <li> "el. potential" is the electrostatic potential [V] used for dirichlet boundary condition </li>
  //     </ol>
  //   </li>
  // </ul>
  
  ProblemParameters::ProblemParameters()
    : ParameterAcceptor("Problem")
  {
    ParameterAcceptor::parse_parameters_call_back.connect(
      [&]() { finalize(); });


    add_parameter("mesh filename", mesh_filename, "mesh file name");

    Assert(mesh_filename.size() > 0, ExcInternalError());

    add_parameter("cells output data filename",
                  cells_output_data_filename,
                  "filename for exporting cells solution");

    Assert(cells_output_data_filename.size() > 0, ExcInternalError());

    add_parameter(
      "domains",
      domains_list,
      "list of domain definitions and its material properties: \"name : ID : epsilon_r : kappa [A/Vm]\"");

    add_parameter("potentials",
                  boundary_values_list,
                  "list of potential definitions for dirichlet boundaries: \"name : ID : phi_e [V]\"");
  }


  // @sect5{ProblemParameters::finalize()}
  // After parsing all parameters, some processing of the domains data and the boundary values is done:
  void ProblemParameters::finalize()
  {
    // 1. copy the domain properties entries into a std::vector
    // Most domain properties are directly copied into memory.
    // Only the electric permittivity is converted into an absolute value [As/Vm]
    // by multiplying the relative permittivity number from parameter file with vacuum permittivity
    // If a domain ID 0 is assigned, its properties will be defined as default values
    
    domain_properties tmp_domain;
    
    int domain_index  = 0;
    default_domain_index = -1;
    
    for (const auto &item : domains_list)
      {
        tmp_domain.name      = std::get<0>(item);
        tmp_domain.domain_id = std::get<1>(item);
        tmp_domain.epsilon   = epsilon_vacuum * std::get<2>(item);;
        tmp_domain.kappa     = std::get<3>(item);
        
        domains_vector.push_back(tmp_domain);

        if (tmp_domain.domain_id == 0) default_domain_index = domain_index;

        domain_index++;
      }

    // if no default domain properties have been defined in the parameter file,
    // following values will be taken as default properties
    if (default_domain_index < 0)
      {
        domain_index++;
        default_domain_index                = domain_index;

        tmp_domain.name        = "default";
        tmp_domain.domain_id   = 0;
        tmp_domain.epsilon     = epsilon_vacuum;
        tmp_domain.kappa       = 0.0;

        domains_vector.push_back(tmp_domain);
      }

    // 2. copy the boundary potentials entries into a std::vector
    boundary_potential tmp_boundary;

    for (const auto &item : boundary_values_list)
      {
        tmp_boundary.name        = std::get<0>(item);
        tmp_boundary.boundary_id = std::get<1>(item);
        tmp_boundary.potential   = std::get<2>(item);
        
        boundaries_vector.push_back(tmp_boundary);
      }
  }


  // Function <code>ProblemParameters::get_domain_properties</code> 
  // will be used to retrieve the domain properties by providing a domain ID number.
  typename ProblemParameters::domain_properties
  ProblemParameters::get_domain_properties(types::material_id domain_id)
  {
    for (const auto & domain : domains_vector)
      {
        if (domain_id == domain.domain_id)
          {
            return domain;
          }
      }
      
    return domains_vector[default_domain_index];
  }
  

  // Function <code>ProblemParameters::get_boundary_value</code> 
  // will be used to retrieve the boundary value by providing a boundary ID number.
  typename ProblemParameters::boundary_potential
  ProblemParameters::get_boundary_value(types::material_id boundary_id)
  {
    for (const auto & boundary : boundaries_vector)
      {
        if (boundary_id == boundary.boundary_id)
          {
            return boundary;
          }
      }
      
    return no_boundary_value;
  }


  // @sect4{Laplace Class}
  // Declaration of all major building blocks for the finite element program
  // using a direct solver on a Laplace problem for electristatic field computation
  // consisting of conductive and dielectric domains.

  template <int dim>
  class Laplace
  {
  public:
    Laplace();
    void run();

    using fieldvector_type = Tensor<1, dim, double>;

    ProblemParameters problem_parameters;

  private:
    MPI_Comm mpi_communicator;

    void make_grid();
    void setup_system();
    void assemble_system();
    void solve();
    void output_results();

    ConditionalOStream pcout;
    TimerOutput        computing_timer;

#ifdef USE_FULLY_DISTRIBUTED_TRIA
    parallel::fullydistributed::Triangulation<dim> triangulation;
#else
    parallel::distributed::Triangulation<dim> triangulation;
#endif

    FE_Q<dim>           fe;
    DoFHandler<dim>     dof_handler;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    AffineConstraints<double> constraints;

    SparsityPattern       sparsity_pattern;
    dealii::LinearAlgebraPETSc::MPI::SparseMatrix system_matrix;
    dealii::LinearAlgebraPETSc::MPI::Vector       locally_relevant_solution;
    dealii::LinearAlgebraPETSc::MPI::Vector       system_rhs;

    class CellsPostprocessor;
  };


  // @sect3{Class Template Definitions and Implementation}
  //
  // @sect4{The Constructor}
  // The Constructor simply consists of default initialization of
  // - the MPI communicator
  // - a timer the performance of several routines
  // - the triangulation
  // - the finite element space
  // - the handler of degrees of freedom

  template <int dim>
  Laplace<dim>::Laplace()
    : mpi_communicator(MPI_COMM_WORLD)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    , computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::never,
                      TimerOutput::wall_times)
    , triangulation(mpi_communicator)
    , fe(/*polinomial degree = */ 1)
    , dof_handler(triangulation)
  {
  }


  // The Laplace::make_grid() routine
  // imports the mesh from a file as specified in the parameter-file.

  template <int dim>
  void Laplace<dim>::make_grid()
  {
    TimerOutput::Scope t(computing_timer, "01. make grid");

    std::string infilename = problem_parameters.mesh_filename + ".msh";

    GridIn<dim> gi;

#ifdef USE_FULLY_DISTRIBUTED_TRIA
    const unsigned int mpi_size =
      Utilities::MPI::n_mpi_processes(mpi_communicator);
      
    auto construction_data = TriangulationDescription::Utilities::
      create_description_from_triangulation_in_groups<dim, dim>(
        [&](Triangulation<dim> &tria) {
          pcout << "mesh-file: \"" << infilename << "\"" << std::endl;
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
#else //  #ifdef USE_FULLY_DISTRIBUTED_TRIA
    Triangulation<dim>                        tria;
    pcout << "mesh-file: \"" << infilename << "\"" << std::endl;
    gi.attach_triangulation(tria);
    gi.read_msh(infilename);

    triangulation.copy_triangulation(tria);
#endif
  }

  // The Laplace::setup_system() routine follows the usual routine of
  // enumerating all the degrees of freedom and setting up the matrix and
  // vector objects to hold the system data. Enumerating is done by using
  // DoFHandler::distribute_dofs().

  template <int dim>
  void Laplace<dim>::setup_system()
  {
    TimerOutput::Scope t(computing_timer, "02. setup");

    dof_handler.distribute_dofs(fe);
    pcout << "setup: number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl;


    // Two index sets that provide information about which degrees of
    // freedom are owned by the current processor and an index set that
    // indicates which degrees of freedom are locally relevant.
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);

    // Initialize the solution and right hand side vectors.
    locally_relevant_solution.reinit(locally_owned_dofs,
                                     locally_relevant_dofs,
                                     mpi_communicator);
    system_rhs.reinit(locally_owned_dofs, mpi_communicator);


    // Compute hanging node and boundary value
    // constraints, which are combined into a single object storing all
    // constraints.
    //
    constraints.clear();

    DoFTools::make_hanging_node_constraints(dof_handler, constraints);


    // Assign electrostatic potential values on the Dirichlet boundaries.
    const FEValuesExtractors::Scalar phi_e(0);

    for (const auto & boundary : problem_parameters.boundaries_vector)
      {
        pcout << "set " << boundary.name
              << " boundary potential: phi_e_"
              << boundary.boundary_id
              << " = "
              << boundary.potential << " V";

        VectorTools::interpolate_boundary_values(
          dof_handler,
          boundary.boundary_id,
          Functions::ConstantFunction<dim, double>(
            boundary.potential),
          constraints,
          fe.component_mask(phi_e));

        pcout << std::endl;
      }


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


  // Assemble the stiffness matrix and the right-hand side:
  template <int dim>
  void Laplace<dim>::assemble_system()
  {
    TimerOutput::Scope t(computing_timer, "03. assemble");

    QGauss<dim>     quadrature_formula(fe.degree + 1);
    QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_gradients |
                            update_quadrature_points | 
                            update_JxW_values);


    FEFaceValues<dim> fe_face_values(fe,
                                     face_quadrature_formula,
                                     update_values | update_gradients |
                                       update_quadrature_points |
                                       update_normal_vectors |
                                       update_JxW_values);


    const auto dofs_per_cell = fe.dofs_per_cell;

    const unsigned int n_q_points      = quadrature_formula.size();
    const unsigned int n_face_q_points      = face_quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // This is assembling the interior of the domain on the left hand side.
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);

            const FEValuesExtractors::Scalar phi_e(0);

            cell_matrix = 0.;
            cell_rhs    = 0.;

            cell->get_dof_indices(local_dof_indices);
            const auto domain_id = cell->material_id();

            auto domain = problem_parameters.get_domain_properties(domain_id);

            const double a = domain.get_domain_parameter_a();

            for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
              {
                for (const auto i : fe_values.dof_indices())
                  {
                    const auto grad_phi_e_i = fe_values[phi_e].gradient(i, q_point);

                    for (const auto j : fe_values.dof_indices())
                      {
                        const auto grad_phi_e_j = fe_values[phi_e].gradient(j, q_point);

                        cell_matrix(i, j) += a * grad_phi_e_j * grad_phi_e_i * fe_values.JxW(q_point);
                      } /* END:  for j */
                  } /* END:  for i */

              } /* END:  for q_point */

            
            // Assemble only the boundary faces without Dirichlet boundary condition
            for (const auto &face : cell->face_iterators())
              {
                if (face->at_boundary())
                  {
                    const auto boundary_id = face->boundary_id();

                    auto boundary_value = problem_parameters.get_boundary_value(boundary_id);
                    if (boundary_value.boundary_id!=problem_parameters.ID_NO_BOUNDARY_VALUE) continue;
                    
                    fe_face_values.reinit(cell, face);

                    for (unsigned int q_point = 0; q_point < n_face_q_points; q_point++)
                      {
                        for (const auto i : fe_face_values.dof_indices())
                          {
                            const auto phi_e_i = fe_face_values[phi_e].value(i, q_point);

                            for (const auto j : fe_face_values.dof_indices())
                              {
                                const auto grad_phi_e_j = fe_face_values[phi_e].gradient(j, q_point);

                                cell_matrix(i, j) += -a * phi_e_i * grad_phi_e_j * fe_face_values.normal_vector(q_point) * fe_face_values.JxW(q_point);
                              } /* END:  for j */

                          } /* END:  for i */

                      } /* END:  for q_point */
                  } /* END:  if (face->at_boundary()) */
              } /* END:  for (const auto &face : cell->face_iterators()) */

            constraints.distribute_local_to_global(cell_matrix,
                                                   cell_rhs,
                                                   local_dof_indices,
                                                   system_matrix,
                                                   system_rhs);


          } /* END: if (cell->is_locally_owned()) */

      } /* END: for (const auto &cell : dof_handler.active_cell_iterators()) */

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);

    pcout << "  system assembled." << std::endl;
  }


  // Use a direct solver to solve the system:
  // SparseDirectMUMPS for MPI parallel implementation
  template <int dim>
  void Laplace<dim>::solve()
  {
    TimerOutput::Scope t(computing_timer, "04. solve");

    dealii::LinearAlgebraPETSc::MPI::Vector completely_distributed_solution(locally_owned_dofs,
                                                    mpi_communicator);

    SolverControl                            solver_control;
    dealii::PETScWrappers::SparseDirectMUMPS solver(solver_control);
    solver.set_symmetric_mode(true);
    solver.solve(system_matrix, completely_distributed_solution, system_rhs);

    pcout << "  solved."  << std::endl;

    constraints.distribute(completely_distributed_solution);

    locally_relevant_solution = completely_distributed_solution;
  }


  // @sect4{Data Output}
  // After solving the problem, some output data shall be written into a vtk-file.
  // The output-filename has been specified in the parameter file and shall contain
  // not only the solution of the Laplace PDE, but also electric field strengths,
  // current densities and domain properties for each evaluation point.
  // For this purpose a Cells Postprocessor is used.
  
  // @sect5{Cells Postprocessor for Data Output}

  // Generate output - using a class PostProcessor that
  // inherits from the class DataPostprocessor, which can be attached to
  // DataOut. This allows to output derived quantities from the solution,
  // For the solved Laplace-PDE, it overloads the
  // virtual function DataPostprocessor::evaluate_scalar_field(),
  // which is then internally called from DataOut::build_patches().

  template <int dim>
  class Laplace<dim>::CellsPostprocessor : public DataPostprocessor<dim>
  {
  public:
    CellsPostprocessor(Laplace<dim> *problem);

    virtual void evaluate_scalar_field(
      const DataPostprocessorInputs::Scalar<dim> &inputs,
      std::vector<Vector<double>> &computed_quantities) const override;

    virtual std::vector<std::string> get_names() const override;

    virtual std::vector<
      DataComponentInterpretation::DataComponentInterpretation>
    get_data_component_interpretation() const override;

    virtual UpdateFlags get_needed_update_flags() const override;

  public:
    Laplace<dim> *problem;
  };


  template <int dim>
  Laplace<dim>::CellsPostprocessor::CellsPostprocessor(Laplace<dim> *problem)
    : problem(problem)
  {}


  // Define the names for the variables in the output.
  // The output data file shall contain following values at the evaluation point locations:
  // - "phi_e",        the scalar electrostatic potential [V]
  // - "Ex",           the electric field strength, cartesian x-component [V/m]
  // - "Ey",           the electric field strength, cartesian y-component [V/m]
  // - "Ez",           the electric field strength, cartesian z-component [V/m]
  // - "Jx",           the electric current density, cartesian x-component [A/m^2]
  // - "Jy",           the electric current density, cartesian y-component [A/m^2]
  // - "Jz",           the electric current density, cartesian z-component [A/m^2]
  // - "domain_ID",    the domain ID number 
  // - "epsilon_r",    the relative electric permittivity number
  // - "kappa"         the electric conductivity [A/Vm]
  template <int dim>
  std::vector<std::string> Laplace<dim>::CellsPostprocessor::get_names() const
  {
    return { 
            "phi_e",        // scalar electrostatic potential
            "Ex",           // electric field strength, cartesian x-component
            "Ey",           // electric field strength, cartesian y-component
            "Ez",           // electric field strength, cartesian z-component
            "Jx",           // electric current density, cartesian x-component
            "Jy",           // electric current density, cartesian y-component
            "Jz",           // electric current density, cartesian z-component
            "domain_ID",    // domain ID number 
            "epsilon_r",    // relative electric permittivity
            "kappa"         // electric conductivity
            };
  }


  // Define the interpretation of the variables in the output
  // - The electrostatic potential is a scalar value
  // - Electric field strength and current density are vectors
  // - ID numbers and material parameters ar scalar values
  //   (such as relative electric permittivity or electric conductivity)
  template <int dim>
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  Laplace<dim>::CellsPostprocessor::get_data_component_interpretation() const
  {
    return { 
            DataComponentInterpretation::component_is_scalar,           // phi_e [V]
            DataComponentInterpretation::component_is_part_of_vector,   // Ex [V/m]
            DataComponentInterpretation::component_is_part_of_vector,   // Ey [V/m]
            DataComponentInterpretation::component_is_part_of_vector,   // Ez [V/m]
            DataComponentInterpretation::component_is_part_of_vector,   // Jx [A/m^2]
            DataComponentInterpretation::component_is_part_of_vector,   // Jy [A/m^2]
            DataComponentInterpretation::component_is_part_of_vector,   // Jz [A/m^2]
            DataComponentInterpretation::component_is_scalar,           // domain_ID  []
            DataComponentInterpretation::component_is_scalar,           // epsilon_r  []
            DataComponentInterpretation::component_is_scalar            // kappa      [A/(V m)]
            };
  }


  // The solution itself is the electrostatic potential.
  // In order to export it, the <code>update_values</code> flag is needed.
  // For plotting electric streamlines or currents, the gradients of the solution 
  // will be neded and hence, the <code>update_gradients</code> has to be set.
  template <int dim>
  UpdateFlags Laplace<dim>::CellsPostprocessor::get_needed_update_flags() const
  {
    return update_values | update_gradients;
  }


  // <code>Laplace<dim>::CellsPostprocessor::evaluate_scalar_field</code>
  // is the function that takes the solution values and its gradients from <code>inputs</code>,
  // calculates electrostatic potential, electric field strengths, electric current density,
  // and domain properties (such as ID number, relative electric permittivity and electric conductivity)
  // and arranges it into <code>computed_quantities</code> for each evaluation point
  template <int dim>
  void Laplace<dim>::CellsPostprocessor::evaluate_scalar_field(
    const DataPostprocessorInputs::Scalar<dim> &inputs,
    std::vector<Vector<double>>                &computed_quantities) const
  {
    unsigned int n_evaluation_points;
    n_evaluation_points = inputs.solution_values.size();
    Assert(inputs.solution_gradients.size() == n_evaluation_points,
           ExcInternalError());

    Assert(computed_quantities.size() == n_evaluation_points,
           ExcInternalError());

    const typename DoFHandler<dim>::cell_iterator current_cell =
      inputs.template get_cell<dim>();

    const auto domain_id = current_cell->material_id();

    const auto domain = problem->problem_parameters.get_domain_properties(domain_id);

    for (unsigned int p = 0; p < n_evaluation_points; p++)
      {
        Laplace::fieldvector_type solution_E_field, solution_current_density_J;

        // The <code>solution_values</code> can be directly taken as the electrostatic potential $\varphi_e$ 
        // (<code>solution_phi_e</code> resp.)
        // as the electric field strength $\mathbf{E}$ is calculated using the <code>solution_gradients</code>.
        // The current density $\mathbf{J}$ is derived from electric field strength under terms of Ohm's law $\mathbf{J}=\kappa\mathbf{E}$.
        const auto solution_phi_e    = inputs.solution_values[p];
        const auto solution_gradient = inputs.solution_gradients[p];

        for (unsigned int d = 0; d < dim; d++)
          {
            solution_E_field[d] = -solution_gradient[d];
          }

        solution_current_density_J = solution_E_field * domain.kappa; // Ohm's law


        // For arranging the output data in the sequence as defined in function 
        // <code>Laplace<dim>::CellsPostprocessor::get_names()</code>
        // the following offset-values are used to indicate the correct positions
        constexpr unsigned int  offset_phi_e                = 0,    // position index for electrostatic potential
                                offset_E_field              = 1,    // position index for first electric field component (Ex)
                                offset_current_density_J    = offset_E_field           + dim, // position index for first current density component (Jx)
                                offset_material_properties  = offset_current_density_J + dim; // position index for domain ID number


        computed_quantities[p](offset_phi_e) = solution_phi_e;

        for (unsigned int d = 0; d < dim; d++)
          {
            computed_quantities[p](offset_E_field           + d) = solution_E_field          [d];
            computed_quantities[p](offset_current_density_J + d) = solution_current_density_J[d];
          }

        computed_quantities[p](offset_material_properties + 0) = domain_id;
        computed_quantities[p](offset_material_properties + 1) = domain.epsilon / epsilon_vacuum;
        computed_quantities[p](offset_material_properties + 2) = domain.kappa;

      } /* END: for (unsigned int p = 0; p < n_evaluation_points; p++) */

  } /* END: Laplace<dim>::CellsPostprocessor::evaluate_field */



  // Initialize the Cells Postprocessor 
  // and write field output data for cells into a file.
  template <int dim>
  void Laplace<dim>::output_results()
  {
    TimerOutput::Scope t(computing_timer, "05. post-process");

    std::string cells_output_data_filename = problem_parameters.cells_output_data_filename + ".vtu";

    CellsPostprocessor cell_postprocessor(this);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(locally_relevant_solution, cell_postprocessor);

    data_out.build_patches();
    data_out.write_vtu_in_parallel(cells_output_data_filename,
                                   mpi_communicator);

    pcout << "exported cells solution and output data into file \""
          << cells_output_data_filename << "\"" << std::endl;
  }



  // @sect4{Main algorithm}

  // This is the central starting-point to run all simulation steps
  template <int dim>
  void Laplace<dim>::run()
  {
    pcout << "Running with PETSc on "
          << Utilities::MPI::n_mpi_processes(mpi_communicator)
          << " MPI rank(s)..."
          << std::endl;

    make_grid();

    setup_system();
    assemble_system();
    solve();

    output_results();

    computing_timer.print_summary();
    computing_timer.reset();

    pcout << std::endl;
  }

} // namespace Step92



// @sect4{Program entry function}
// An optional parameter can be passed via the command line call: a parameter filename
// If this optional parameter filename is not provided, parameter file "example-1.prm" will be taken per default.
// In the main function
// - MPI is initialized
// - Laplace class object for dim=3 is instantiated,
// - the ParameterAcceptor is initialized and
// finally the main algorithm <code>run()</code> is called.


int main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      Step92::Laplace<3> laplace_solver;


      ParameterAcceptor::initialize(argc < 2 ? "example-1.prm" : argv[1]);

      laplace_solver.run();
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

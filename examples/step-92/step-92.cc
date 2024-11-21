/* ---------------------------------------------------------------------
 *
 * This example can be compiled for usage with parallel::distributed or
 * parallel::fullydistributed triangulations
 *
 * In order to force usage of parallel::fullydistributed triangulations
 * define
 * #define USE_FULLY_DISTRIBUTED_TRIA
 *
 * otherwise (if not defined), parallel::distributed triangulation
 * will be used
 * ---------------------------------------------------------------------
 */
#define USE_FULLY_DISTRIBUTED_TRIA

/* ---------------------------------------------------------------------
 *
 * 2024-11, Stephan Voss, Neunkirchen am Brand, Germany, stvoss@gmx.de
 *
 * This solver is derived from several deal.II examples and tutorials.
 * Thank You to all deal.II contributors.
 *
 * ---------------------------------------------------------------------
 */

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
// To compare the performance of PETSc and Trilinos
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
#ifdef USE_FULLY_DISTRIBUTED_TRIA
#  include <deal.II/distributed/fully_distributed_tria.h>
#else
#  include <deal.II/distributed/tria.h>
#endif
#include <deal.II/distributed/grid_refinement.h>


#include <fstream>
#include <iostream>
#include <memory>


// @sect3{Class Template Declarations}
//
// Declaration of all classes with their data structures and methods

namespace Step92
{
  template <typename T>
  T sum(const T &, const T &);

  template <typename T>
  double sum(const double &a, const double &b)
  {
    return a + b;
  };

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

    using fieldvector_type = Tensor<1, dim, double>;

    using physical_ID = unsigned int;

    using Material_tuple = std::tuple<std::string, physical_ID, double, double>;

    using BoundaryPotential_tuple =
      std::tuple<std::string, physical_ID, double>;

    class material_definition
    {
    public:
      types::material_id material_id;
      std::string        name;
      double epsilon; // absolute electric permittivity epsilon (not relative !!
                      // --> eps_r * eps_0)
      double kappa;   // absolute electric conductivity
    };

    class boundary_potential
    {
    public:
      types::material_id material_id;
      std::string        name;
      double             potential; // [V] electric potential
    };


  public:
    material_definition get_physical(types::material_id material);

    boundary_potential get_boundary_value(types::material_id material);

    std::string get_mesh_filename();
    std::string get_cells_solution_filename();

    // double mu_vacuum;      // vacuum permeability: 4 pi 10e-7 V s / (A m)
    double epsilon_vacuum; // vacuum permittivity: 8.8541878188(14)×10^−12 A s /
                           // (V m)

    std::string mesh_filename, cells_solution_filename;

  private:
    int default_physical_ix = -1;

    std::list<Material_tuple>          materials_list;
    std::list<BoundaryPotential_tuple> boundary_values_list;

  public:
    unsigned long        N_physicals = 0;
    material_definition *p_physicals = NULL;

    unsigned long       N_boundary_values = 0;
    boundary_potential *p_boundary_values = NULL;
  };


  template <int dim>
  ProblemParameters<dim>::ProblemParameters()
    : ParameterAcceptor("Problem")
  {
    // mu_vacuum      = 4.0e-7 * dealii::numbers::PI; // [ V s / (A m) ]
    epsilon_vacuum = 8.8541878176204199e-12; // [ A s / (V m) ]


    mesh_filename = "";
    add_parameter("mesh filename", mesh_filename, "mesh file name");
    cells_solution_filename = "";
    add_parameter("cells solution filename",
                  cells_solution_filename,
                  "filename for exporting cells solution");

    add_parameter(
      "materials",
      materials_list,
      "list of material definitions: \"name : ID : eps_r : kappa\"");
    add_parameter("potentials",
                  boundary_values_list,
                  "list of potential definitions: \"name : ID : epot [V]\"");
  }

  template <int dim>
  void ProblemParameters<dim>::finalize()
  {
    N_physicals = materials_list.size() + 1;

    p_physicals = new material_definition[N_physicals];

    int ix              = 0;
    default_physical_ix = -1;
    for (const auto &item : materials_list)
      {
        std::string  name  = std::get<0>(item);
        unsigned int id    = std::get<1>(item);
        double       eps_r = std::get<2>(item);
        double       kappa = std::get<3>(item);

        p_physicals[ix].name = name;

        p_physicals[ix].material_id = id;
        if (id == 0)
          default_physical_ix = ix;

        p_physicals[ix].epsilon = eps_r * epsilon_vacuum;
        p_physicals[ix].kappa   = kappa;

        ix++;
      }

    N_physicals = ix;

    if (default_physical_ix < 0)
      {
        N_physicals++;
        p_physicals[ix].name = "default";

        default_physical_ix         = ix;
        p_physicals[ix].material_id = 0;

        p_physicals[ix].epsilon = epsilon_vacuum;
        p_physicals[ix].kappa   = 0.0;
      }

    N_boundary_values = boundary_values_list.size();
    p_boundary_values = new boundary_potential[N_boundary_values];

    unsigned int        ipot = 0;
    boundary_potential *t_p_opdata;
    for (const auto &item : boundary_values_list)
      {
        std::string name      = std::get<0>(item);
        auto        id        = std::get<1>(item);
        double      potential = std::get<2>(item);

        t_p_opdata              = &p_boundary_values[ipot];
        t_p_opdata->material_id = id;
        t_p_opdata->name        = name;
        t_p_opdata->potential   = potential;

        ipot++;
      }
  }

  template <int dim>
  ProblemParameters<dim>::~ProblemParameters()
  {}

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
  typename ProblemParameters<dim>::material_definition
  ProblemParameters<dim>::get_physical(types::material_id material)
  {
    for (unsigned int ix = 0; ix < N_physicals; ix++)
      {
        if (material == (p_physicals[ix].material_id))
          {
            return p_physicals[ix];
          }
      }

    return p_physicals[default_physical_ix];
  }

  template <int dim>
  typename ProblemParameters<dim>::boundary_potential
  ProblemParameters<dim>::get_boundary_value(types::material_id material)
  {
    for (unsigned int ix = 0; ix < N_boundary_values; ix++)
      {
        if (material == (p_boundary_values[ix].material_id))
          {
            return p_boundary_values[ix];
          }
      }

    return *p_boundary_values[0];
  }

  // @sect4{Laplace Class}
  // Declaration of all major building blocks of
  // finite element program which consists of the usual setup and
  // assembly routines.

  template <int dim>
  class Laplace : public ParameterAcceptor
  {
  public:
    Laplace();
    void run();

    using fieldvector_type = Tensor<1, dim, double>;

    ProblemParameters<dim> problem_parameters;

  private:
    MPI_Comm mpi_communicator;

    /* run time parameters */
    unsigned int N_materials;
    unsigned int fe_order;
    unsigned int quadrature_order;

    void parse_parameters_callback();
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

    std::unique_ptr<FiniteElement<dim>> fe;
    DoFHandler<dim>                     dof_handler;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    AffineConstraints<double> constraints;

    SparsityPattern       sparsity_pattern;
    LA::MPI::SparseMatrix system_matrix;
    LA::MPI::Vector       locally_relevant_solution;
    LA::MPI::Vector       system_rhs;

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
  // these can be modified by editing the .prm file.

  template <int dim>
  Laplace<dim>::Laplace()
    : ParameterAcceptor("Solver")
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

    N_materials = 0;
    add_parameter("materials", N_materials, "number of materials");

    fe_order = 0;
    add_parameter(
      "fe order",
      fe_order,
      "order of the finite element space for scalar electric potential");

    quadrature_order = 1;
    add_parameter("quadrature order",
                  quadrature_order,
                  "order of the quadrature");
  }


  template <int dim>
  void Laplace<dim>::parse_parameters_callback()
  {
    problem_parameters.finalize();

    fe = std::make_unique<FESystem<dim>>(FE_Q<dim>(fe_order));
  }

  // The Laplace::make_grid() routine
  // reads the mesh-file for the computational domain.
  // A block decomposition into real and imaginary matrices
  // for the solution matrices is used.

  template <int dim>
  void Laplace<dim>::make_grid()
  {
    TimerOutput::Scope t(computing_timer, "01. make grid");

    std::string infilename = problem_parameters.get_mesh_filename();
    infilename += ".msh";

    GridIn<dim> gi;

#ifdef USE_FULLY_DISTRIBUTED_TRIA
    const unsigned int mpi_size =
      Utilities::MPI::n_mpi_processes(mpi_communicator);
    auto construction_data = TriangulationDescription::Utilities::
      create_description_from_triangulation_in_groups<dim, dim>(
        [&](Triangulation<dim> &tria) {
          pcout << "MPI-comm: "
                << Utilities::MPI::this_mpi_process(mpi_communicator)
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
#else //  #ifdef USE_FULLY_DISTRIBUTED_TRIA
    Triangulation<dim>                        tria;
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

    dof_handler.distribute_dofs(*fe);
    pcout << "  number of fe_systems: " << fe->n_components() << std::endl
          << "  setup: number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl;


    // Two index sets that provide information about which degrees of
    // freedom are owned by the current processor and an index set that
    // indicates which degrees of freedom are locally relevant.
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);

    // Initialize the solution and right hand side vectors.
    // Solution vector we seek does not only store
    // elements locally owned, but also ghost entries;
    // Right hand side vector only needs to have the entries the current
    // processor owns (It will only be written locally- never read.) (Linear
    // solvers will read from it, but they do not care about the geometric
    // location of degrees of freedom).
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

    const FEValuesExtractors::Scalar phi_E(0);

    for (unsigned int ix = 0; ix < problem_parameters.N_boundary_values; ix++)
      {
        pcout << "set " << problem_parameters.p_boundary_values[ix].name
              << " boundary (ID "
              << problem_parameters.p_boundary_values[ix].material_id
              << ") potential: "
              << problem_parameters.p_boundary_values[ix].potential << " V";

        VectorTools::interpolate_boundary_values(
          dof_handler,
          problem_parameters.p_boundary_values[ix].material_id,
          Functions::ConstantFunction<dim, double>(
            problem_parameters.p_boundary_values[ix].potential),
          constraints,
          fe->component_mask(phi_E));

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

    QGauss<dim>     quadrature_formula(quadrature_order);
    QGauss<dim - 1> face_quadrature_formula(quadrature_order);

    FEValues<dim> fe_values(*fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    FEFaceValues<dim> fe_face_values(*fe,
                                     face_quadrature_formula,
                                     update_values | update_gradients |
                                       update_quadrature_points |
                                       update_normal_vectors |
                                       update_JxW_values);

    const auto dofs_per_cell = fe->dofs_per_cell;

    const unsigned int                  n_q_points = quadrature_formula.size();
    [[maybe_unused]] const unsigned int n_face_q_points =
      face_quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // This is assembling the interior of the domain on the left hand side.
    // In doing so, test functions $\varphi_i$ and $\varphi_j$ are needed, and
    // the curl of these test variables.

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);

            const FEValuesExtractors::Scalar phi_E(0);

            cell_matrix = 0.;
            cell_rhs    = 0.;

            cell->get_dof_indices(local_dof_indices);
            const auto material_id = cell->material_id();

            const auto material = problem_parameters.get_physical(material_id);

            double psi_E =
              (material.kappa <= 0 ? material.epsilon : material.kappa);

            for (unsigned int q_point = 0; q_point < n_q_points; q_point++)
              {
                for (const auto i : fe_values.dof_indices())
                  {
                    const auto grad_phi_E_i =
                      fe_values[phi_E].gradient(i, q_point);

                    for (const auto j : fe_values.dof_indices())
                      {
                        const auto grad_phi_E_j =
                          fe_values[phi_E].gradient(j, q_point);

                        const auto temp =
                          psi_E * scalar_product(grad_phi_E_j, grad_phi_E_i) *
                          fe_values.JxW(q_point);

                        cell_matrix(i, j) += temp;
                      } // for j
                  }     // for i

              } // for qpoint



            /*------------------------------------------------------------------------------------------
             *
             * FACES integral
             *
             *------------------------------------------------------------------------------------------*/

            // Assemble the face and the boundary.

            for (const auto &face : cell->face_iterators())
              {
                if (face->at_boundary())
                  {
                    fe_face_values.reinit(cell, face);

                    const FEValuesExtractors::Scalar phi_E(0);
                    for (unsigned int q_point = 0; q_point < n_face_q_points;
                         q_point++)
                      {
                        const auto normal =
                          fe_face_values.normal_vector(q_point);

                        for (const auto i : fe_face_values.dof_indices())
                          {
                            const auto phi_E_i =
                              fe_face_values[phi_E].value(i, q_point);

                            for (const auto j : fe_face_values.dof_indices())
                              {
                                const auto grad_phi_E_j =
                                  fe_face_values[phi_E].gradient(j, q_point);

                                const auto temp =
                                  -psi_E * phi_E_i *
                                  scalar_product(grad_phi_E_j, normal) *
                                  fe_face_values.JxW(q_point);

                                cell_matrix(i, j) += temp;
                              } /* for(j) */

                          } /* for(i) */

                      } /* for(q_point) */
                  }     // if (face->at_boundary())
              }         // END: for (const auto &face : cell->face_iterators())



            constraints.distribute_local_to_global(cell_matrix,
                                                   cell_rhs,
                                                   local_dof_indices,
                                                   system_matrix,
                                                   system_rhs);


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
  void Laplace<dim>::solve()
  {
    TimerOutput::Scope t(computing_timer, "04. solve");

    LA::MPI::Vector completely_distributed_solution(locally_owned_dofs,
                                                    mpi_communicator);

    SolverControl                            solver_control;
    dealii::PETScWrappers::SparseDirectMUMPS solver(solver_control);
    solver.set_symmetric_mode(true);
    solver.solve(system_matrix, completely_distributed_solution, system_rhs);

    auto nsteps = solver_control.last_step();
    pcout << "  solved in " << nsteps << " iteration"
          << (nsteps > 1 ? "s." : ".") << std::endl;

    constraints.distribute(completely_distributed_solution);

    locally_relevant_solution = completely_distributed_solution;
  }


  // @sect4{Cells Postprocessor for Data Output}

  // Generate output - using a class PostProcessor that
  // inherits from the class DataPostprocessor, which can be attached to
  // DataOut. This allows to output derived quantities from the solution,
  // It overloads the
  // virtual function DataPostprocessor::evaluate_vector_field(),
  // which is then internally called from DataOut::build_patches().
  // It's given values of the numerical solution, its derivatives, normals to
  // the cell, the actual evaluation points and any additional quantities.

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
  template <int dim>
  std::vector<std::string> Laplace<dim>::CellsPostprocessor::get_names() const
  {
    std::vector<std::string> solution_names = {};

    solution_names.emplace_back("phi_e"); // electric potential [V]

    solution_names.emplace_back("Ex"); // electric field strength [V/m]
    solution_names.emplace_back("Ey");
    solution_names.emplace_back("Ez");

    solution_names.emplace_back("Jx"); // electric current density [A/m^2]
    solution_names.emplace_back("Jy");
    solution_names.emplace_back("Jz");

    solution_names.emplace_back("material_ID");
    solution_names.emplace_back("epsilon_r");
    solution_names.emplace_back("kappa");

    return solution_names;
  }


  template <int dim>
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  Laplace<dim>::CellsPostprocessor::get_data_component_interpretation() const
  {
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      interpretation;

    interpretation.push_back(
      DataComponentInterpretation::component_is_scalar); // phi_e

    interpretation.push_back(
      DataComponentInterpretation::component_is_part_of_vector); // E x,y,z
    interpretation.push_back(
      DataComponentInterpretation::component_is_part_of_vector);
    interpretation.push_back(
      DataComponentInterpretation::component_is_part_of_vector);

    interpretation.push_back(
      DataComponentInterpretation::component_is_part_of_vector); // J x,y,z
    interpretation.push_back(
      DataComponentInterpretation::component_is_part_of_vector);
    interpretation.push_back(
      DataComponentInterpretation::component_is_part_of_vector);


    interpretation.push_back(
      DataComponentInterpretation::component_is_scalar); // mat_ID
    interpretation.push_back(
      DataComponentInterpretation::component_is_scalar); // epsilon_r
    interpretation.push_back(
      DataComponentInterpretation::component_is_scalar); // kappa

    return interpretation;
  }



  template <int dim>
  UpdateFlags Laplace<dim>::CellsPostprocessor::get_needed_update_flags() const
  {
    return update_values | update_gradients | update_quadrature_points;
  }


  // This is the function that computes the derived quantities.
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

    const auto material_id = current_cell->material_id();

    const auto material = problem->problem_parameters.get_physical(material_id);

    for (unsigned int p = 0; p < n_evaluation_points; p++)
      {
        Laplace::fieldvector_type sE, sJ;
        double                    Phi;


        const auto t_sol_val   = inputs.solution_values[p];
        const auto t_sol_grads = inputs.solution_gradients[p];

        Phi = t_sol_val;

        for (unsigned int d = 0; d < dim; d++)
          {
            sE[d] = -t_sol_grads[d];
          } // for (unsigned int d = 0; d < dim; d++)

        sJ = sE * material.kappa;

        const unsigned int off_E_ePot = 0, off_E_vE = 1,
                           off_E_vJ  = off_E_vE + dim,
                           off_E_mat = off_E_vJ + dim;


        computed_quantities[p](off_E_ePot) = Phi;

        for (unsigned int d = 0; d < dim; d++)
          {
            computed_quantities[p](off_E_vE + d) = sE[d];
            computed_quantities[p](off_E_vJ + d) = sJ[d];
          }

        computed_quantities[p](off_E_mat + 0) = material_id;
        computed_quantities[p](off_E_mat + 1) =
          material.epsilon / problem->problem_parameters.epsilon_vacuum;
        computed_quantities[p](off_E_mat + 2) = material.kappa;

      } // for (unsigned int p = 0; p < n_evaluation_points; p++)

  } // Laplace<dim>::CellsPostprocessor::evaluate_field


  // @sect4{Data Output}

  // Write field outputs into files (for cells and faces)
  template <int dim>
  void Laplace<dim>::output_results()
  {
    TimerOutput::Scope t(computing_timer, "05. post-process");

    std::string cells_solution_filename =
      problem_parameters.get_cells_solution_filename();

    if (cells_solution_filename.size() > 0)
      {
        CellsPostprocessor cell_postprocessor(this);

        cells_solution_filename += ".vtu";

        DataOut<dim> data_out;
        data_out.attach_dof_handler(dof_handler);
        data_out.add_data_vector(locally_relevant_solution, cell_postprocessor);

        data_out.build_patches();
        data_out.write_vtu_in_parallel(cells_solution_filename,
                                       mpi_communicator);

        pcout << "exported cells solution into file \""
              << cells_solution_filename << "\"" << std::endl;
      }
  }



  // @sect4{Main algorithm}

  // the central starting-point to run all simulation steps
  template <int dim>
  void Laplace<dim>::run()
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

    setup_system();
    assemble_system();
    solve();

    output_results();

    computing_timer.print_summary();
    computing_timer.reset();

    pcout << std::endl;
  }

} // namespace Step92


// main function call: set up MPI, Laplace classes,
// ParameterAcceptor, and call the run() function.

int main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      Step92::Laplace<3> laplace_solver;


      ParameterAcceptor::initialize(argc < 2 ? "parameters.prm" : argv[1]);

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

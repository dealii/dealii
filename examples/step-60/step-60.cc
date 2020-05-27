/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2018 - 2020 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 *
 * Authors: Luca Heltai, Giovanni Alzetta,
 * International School for Advanced Studies, Trieste, 2018
 */

// @sect3{Include files}
// Most of these have been introduced elsewhere, we'll comment only on the new
// ones.

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>

// The parameter acceptor class is the first novelty of this tutorial program:
// in general parameter files are used to steer the execution of a program at
// run time. While even a simple approach saves compile time, as the same
// executable can be run with different parameter settings, it can become
// difficult to handle hundreds of parameters simultaneously while maintaining
// compatibility between different programs. This is where the class
// ParameterAcceptor proves useful.
//
// This class is used to define a public interface for classes that want to use
// a single global ParameterHandler to handle parameters. The class provides a
// static ParameterHandler member, namely ParameterAcceptor::prm, and
// implements the "Command design pattern" (see, for example, E. Gamma, R. Helm,
// R. Johnson, J. Vlissides, Design Patterns: Elements of Reusable
// Object-Oriented Software, Addison-Wesley Professional, 1994.
// https://goo.gl/FNYByc).
//
// ParameterAcceptor provides a global subscription mechanism. Whenever an
// object of a class derived from ParameterAcceptor is constructed, a pointer
// to that object-of-derived-type is registered, together with a section entry
// in the parameter file. Such registry is traversed upon invocation of the
// single function ParameterAcceptor::initialize("file.prm") which in turn makes
// sure that all classes stored in the global registry declare the parameters
// they will be using, and after having declared them, it reads the content of
// `file.prm` to parse the actual parameters.
//
// If you call the method ParameterHandler::add_parameter for each of the
// parameters you want to use in your code, there is nothing else you need to
// do. If you are using an already existing class that provides the two
// functions `declare_parameters` and `parse_parameters`, you can still use
// ParameterAcceptor, by encapsulating the existing class into a
// ParameterAcceptorProxy class.
//
// In this example, we'll use both strategies, using ParameterAcceptorProxy for
// deal.II classes, and deriving our own parameter classes directly from
// ParameterAcceptor.
#include <deal.II/base/parameter_acceptor.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

// The other new include file is the one that contains the GridTools::Cache
// class. The structure of deal.II, as many modern numerical libraries, is
// organized following a Directed Acyclic Graph (DAG). A DAG is a directed graph
// with topological ordering: each node structurally represents an object, and
// is connected to child nodes by one (or more) oriented edges, from the parent
// to the child. The most significant example of this structure is the
// Triangulation and its Triangulation::cell_iterator structure. From a
// Triangulation (the main node), we can access each cell (children nodes of the
// triangulation). From the cells themselves we can access over all vertices of
// the cell. In this simple example, the DAG structure can be represented as
// three node types (the triangulation, the cell iterator, and the vertex)
// connected by oriented edges from the triangulation to the cell iterators, and
// from the cell iterator to the vertices. This has several advantages, but it
// intrinsically creates “asymmetries”, making certain operations fast and their
// inverse very slow: finding the vertices of a cell has low computational cost,
// and can be done by simply traversing the DAG, while finding all the cells
// that share a vertex requires a non-trivial computation unless a new DAG data
// structure is added that represents the inverse search.
//
// Since inverse operations are usually not needed in a finite element code,
// these are implemented in GridTools without the use of extra data structures
// related to the Triangulation which would make them much faster. One such data
// structure, for example, is a map from the vertices of a Triangulation to all
// cells that share those vertices, which would reduce the computations needed
// to answer to the previous question.
//
// Some methods, for example GridTools::find_active_cell_around_point, make
// heavy usage of these non-standard operations. If you need to call these
// methods more than once, it becomes convenient to store those data structures
// somewhere. GridTools::Cache does exactly this, giving you access to
// previously computed objects, or computing them on the fly (and then storing
// them inside the class for later use), and making sure that whenever the
// Triangulation is updated, also the relevant data structures are recomputed.
#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

// In this example, we will be using a reference domain to describe an embedded
// Triangulation, deformed through a finite element vector field.
//
// The next two include files contain the definition of two classes that can be
// used in these cases. MappingQEulerian allows one to describe a domain through
// a *displacement* field, based on a FESystem[FE_Q(p)^spacedim] finite element
// space. The second is a little more generic, and allows you to use arbitrary
// vector FiniteElement spaces, as long as they provide a *continuous*
// description of your domain. In this case, the description is done through the
// actual *deformation* field, rather than a *displacement* field.
//
// Which one is used depends on how the user wants to specify the reference
// domain, and/or the actual configuration. We'll provide both options, and
// experiment a little in the results section of this tutorial program.
#include <deal.II/fe/mapping_q_eulerian.h>
#include <deal.II/fe/mapping_fe_field.h>

#include <deal.II/dofs/dof_tools.h>

// The parsed function class is another new entry. It allows one to create a
// Function object, starting from a string in a parameter file which is parsed
// into an object that you can use anywhere deal.II accepts a Function (for
// example, for interpolation, boundary conditions, etc.).
#include <deal.II/base/parsed_function.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

// This is the last new entry for this tutorial program. The namespace
// NonMatching contains a few methods that are useful when performing
// computations on non-matching grids, or on curves that are not aligned with
// the underlying mesh.
//
// We'll discuss its use in detail later on in the `setup_coupling` method.
#include <deal.II/non_matching/coupling.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/linear_operator_tools.h>

#include <iostream>
#include <fstream>

namespace Step60
{
  using namespace dealii;

  // @sect3{DistributedLagrangeProblem}
  //
  // In the DistributedLagrangeProblem, we need two parameters describing the
  // dimensions of the domain $\Gamma$ (`dim`) and of the domain $\Omega$
  // (`spacedim`).
  //
  // These will be used to initialize a Triangulation<dim,spacedim> (for
  // $\Gamma$) and a Triangulation<spacedim,spacedim> (for $\Omega$).
  //
  // A novelty with respect to other tutorial programs is the heavy use of
  // std::unique_ptr. These behave like classical pointers, with the advantage
  // of doing automatic house-keeping: the contained object is automatically
  // destroyed as soon as the unique_ptr goes out of scope, even if it is inside
  // a container or there's an exception. Moreover it does not allow for
  // duplicate pointers, which prevents ownership problems. We do this, because
  // we want to be able to i) construct the problem, ii) read the parameters,
  // and iii) initialize all objects according to what is specified in a
  // parameter file.
  //
  // We construct the parameters of our problem in the internal class
  // `Parameters`, derived from ParameterAcceptor. The
  // `DistributedLagrangeProblem` class takes a const reference to a
  // `Parameters` object, so that it is not possible
  // to modify the parameters from within the DistributedLagrangeProblem class
  // itself.
  //
  // We could have initialized the parameters first, and then pass the
  // parameters to the DistributedLagrangeProblem assuming all entries are set
  // to the desired values, but this has two disadvantages:
  //
  // - We should not make assumptions on how the user initializes a class that
  // is not under our direct control. If the user fails to initialize the
  // class, we should notice and throw an exception;
  //
  // - Not all objects that need to read parameters from a parameter file may
  // be available when we construct the Parameters;
  // this is often the case for complex programs, with multiple physics, or
  // where we reuse existing code in some external classes. We simulate this by
  // keeping some "complex" objects, like ParsedFunction objects, inside the
  // `DistributedLagrangeProblem` instead of inside the
  // `Parameters`.
  //
  // Here we assume that upon construction, the classes that build up our
  // problem are not usable yet. Parsing the parameter file is what ensures we
  // have all ingredients to build up our classes, and we design them so that if
  // parsing fails, or is not executed, the run is aborted.

  template <int dim, int spacedim = dim>
  class DistributedLagrangeProblem
  {
  public:
    // The `Parameters` class is derived from ParameterAcceptor. This allows us
    // to use the ParameterAcceptor::add_parameter() method in its constructor.
    //
    // The members of this function are all non-const, but the
    // `DistributedLagrangeProblem` class takes a const reference to a
    // `Parameters` object: this ensures that
    // parameters are not modified from within the `DistributedLagrangeProblem`
    // class.
    class Parameters : public ParameterAcceptor
    {
    public:
      Parameters();

      // The parameters now described can all be set externally using a
      // parameter file: if no parameter file is present when running the
      // executable, the program will create a "parameters.prm" file with the
      // default values defined here, and then abort to give the user a chance
      // to modify the parameters.prm file.

      // Initial refinement for the embedding grid, corresponding to the domain
      // $\Omega$.
      unsigned int initial_refinement = 4;

      // The interaction between the embedded grid $\Omega$ and the embedding
      // grid $\Gamma$ is handled through the computation of $C$, which
      // involves all cells of $\Omega$ overlapping with parts of $\Gamma$:
      // a higher refinement of such cells might improve quality of our
      // computations.
      // For this reason we define `delta_refinement`: if it is greater
      // than zero, then we mark each cell of the space grid that contains
      // a vertex of the embedded grid and its neighbors, execute the
      // refinement, and repeat this process `delta_refinement` times.
      unsigned int delta_refinement = 3;

      // Starting refinement of the embedded grid, corresponding to the domain
      // $\Gamma$.
      unsigned int initial_embedded_refinement = 8;

      // The list of boundary ids where we impose homogeneous Dirichlet boundary
      // conditions. On the remaining boundary ids (if any), we impose
      // homogeneous Neumann boundary conditions.
      // As a default problem we have zero Dirichlet boundary conditions on
      // $\partial \Omega$
      std::list<types::boundary_id> homogeneous_dirichlet_ids{0, 1, 2, 3};

      // FiniteElement degree of the embedding space: $V_h(\Omega)$
      unsigned int embedding_space_finite_element_degree = 1;

      // FiniteElement degree of the embedded space: $Q_h(\Gamma)$
      unsigned int embedded_space_finite_element_degree = 1;

      // FiniteElement degree of the space used to describe the deformation
      // of the embedded domain
      unsigned int embedded_configuration_finite_element_degree = 1;

      // Order of the quadrature formula used to integrate the coupling
      unsigned int coupling_quadrature_order = 3;

      // If set to true, then the embedded configuration function is
      // interpreted as a displacement function
      bool use_displacement = false;

      // Level of verbosity to use in the output
      unsigned int verbosity_level = 10;

      // A flag to keep track if we were initialized or not
      bool initialized = false;
    };

    DistributedLagrangeProblem(const Parameters &parameters);

    // Entry point for the DistributedLagrangeProblem
    void run();

  private:
    // Object containing the actual parameters
    const Parameters &parameters;

    // The following functions are similar to all other tutorial programs, with
    // the exception that we now need to set up things for two different
    // families of objects, namely the ones related to the *embedding* grids,
    // and the ones related to the *embedded* one.

    void setup_grids_and_dofs();

    void setup_embedding_dofs();

    void setup_embedded_dofs();

    // The only unconventional function we have here is the `setup_coupling()`
    // method, used to generate the sparsity patter for the coupling matrix $C$.

    void setup_coupling();

    void assemble_system();

    void solve();

    void output_results();


    // first we gather all the objects related to the embedding space geometry

    std::unique_ptr<Triangulation<spacedim>> space_grid;
    std::unique_ptr<GridTools::Cache<spacedim, spacedim>>
                                             space_grid_tools_cache;
    std::unique_ptr<FiniteElement<spacedim>> space_fe;
    std::unique_ptr<DoFHandler<spacedim>>    space_dh;

    // Then the ones related to the embedded grid, with the DoFHandler
    // associated to the Lagrange multiplier `lambda`

    std::unique_ptr<Triangulation<dim, spacedim>> embedded_grid;
    std::unique_ptr<FiniteElement<dim, spacedim>> embedded_fe;
    std::unique_ptr<DoFHandler<dim, spacedim>>    embedded_dh;

    // And finally, everything that is needed to *deform* the embedded
    // triangulation
    std::unique_ptr<FiniteElement<dim, spacedim>> embedded_configuration_fe;
    std::unique_ptr<DoFHandler<dim, spacedim>>    embedded_configuration_dh;
    Vector<double>                                embedded_configuration;

    // The ParameterAcceptorProxy class is a "transparent" wrapper derived
    // from both ParameterAcceptor and the type passed as its template
    // parameter. At construction, the arguments are split into two parts: the
    // first argument is an std::string, forwarded to the ParameterAcceptor
    // class, and containing the name of the section that should be used for
    // this class, while all the remaining arguments are forwarded to the
    // constructor of the templated type, in this case, to the
    // Functions::ParsedFunction constructor.
    //
    // This class allows you to use existing classes in conjunction with the
    // ParameterAcceptor registration mechanism, provided that those classes
    // have the members `declare_parameters()` and `parse_parameters()`.
    //
    // This is the case here, making it fairly easy to exploit the
    // Functions::ParsedFunction class: instead of requiring users to create new
    // Function objects in their code for the RHS, boundary functions, etc.,
    // (like it is done in most of the other tutorials), here we allow the user
    // to use deal.II interface to muParser (http://muparser.beltoforion.de),
    // where the specification of the function is not done at compile time, but
    // at run time, using a string that is parsed into an actual Function
    // object.
    //
    // In this case, the `embedded_configuration_function` is a vector valued
    // Function that can be interpreted as either a *deformation* or a
    // *displacement* according to the boolean value of
    // `parameters.use_displacement`. The number of components is specified
    // later on in the construction.

    ParameterAcceptorProxy<Functions::ParsedFunction<spacedim>>
      embedded_configuration_function;

    std::unique_ptr<Mapping<dim, spacedim>> embedded_mapping;

    // We do the same thing to specify the value of the function $g$,
    // which is what we want our solution to be in the embedded space.
    // In this case the Function is a scalar one.
    ParameterAcceptorProxy<Functions::ParsedFunction<spacedim>>
      embedded_value_function;

    // Similarly to what we have done with the Functions::ParsedFunction class,
    // we repeat the same for the ReductionControl class, allowing us to
    // specify all possible stopping criteria for the Schur complement
    // iterative solver we'll use later on.
    ParameterAcceptorProxy<ReductionControl> schur_solver_control;

    // Next we gather all SparsityPattern, SparseMatrix, and Vector objects
    // we'll need
    SparsityPattern stiffness_sparsity;
    SparsityPattern coupling_sparsity;

    SparseMatrix<double> stiffness_matrix;
    SparseMatrix<double> coupling_matrix;

    AffineConstraints<double> constraints;

    Vector<double> solution;
    Vector<double> rhs;

    Vector<double> lambda;
    Vector<double> embedded_rhs;
    Vector<double> embedded_value;

    // The TimerOutput class is used to provide some statistics on
    // the performance of our program.
    TimerOutput monitor;
  };

  // @sect3{DistributedLagrangeProblem::Parameters}
  //
  // At construction time, we initialize also the ParameterAcceptor class, with
  // the section name we want our problem to use when parsing the parameter
  // file.
  //
  // Parameter files can be organized into section/subsection/etc.:
  // this has the advantage that defined objects share parameters when
  // sharing the same section/subsection/etc. ParameterAcceptor allows
  // to specify the section name using Unix conventions on paths.
  // If the section name starts with a slash ("/"), then the section is
  // interpreted as an *absolute path*, ParameterAcceptor enters a subsection
  // for each directory in the path, using the last name it encountered as
  // the landing subsection for the current class.
  //
  // For example, if you construct your class using
  // `ParameterAcceptor("/first/second/third/My Class")`, the parameters will be
  // organized as follows:
  //
  // @code
  // # Example parameter file
  // subsection first
  //   subsection second
  //     subsection third
  //       subsection My Class
  //        ... # all the parameters
  //       end
  //     end
  //   end
  // end
  // @endcode
  //
  // Internally, the *current path* stored in ParameterAcceptor is now
  // considered to be "/first/second/third/", i.e. when you specify an
  // absolute path, ParameterAcceptor *changes* the current section to the
  // current path, i.e. to the path of the section name until the *last* "/".
  //
  // You can now construct another class derived from ParameterAcceptor using a
  // relative path (e.g., `ParameterAcceptor("My Other Class")`) instead of the
  // absolute one (e.g. `ParameterAcceptor("/first/second/third/My Other
  // Class")`), obtaining:
  // @code
  // # Example parameter file
  // subsection first
  //   subsection second
  //     subsection third
  //       subsection My Class
  //         ... # all the parameters
  //       end
  //       subsection My Other Class
  //         ... # all the parameters of MyOtherClass
  //       end
  //     end
  //   end
  // end
  // @endcode
  //
  // If the section name *ends* with a slash then subsequent classes will
  // interpret this as a full path: for example, similar to the one above, if
  // we have two classes, one initialized with
  // `ParameterAcceptor("/first/second/third/My Class/")`
  // and the other with `ParameterAcceptor("My Other Class")`, then the
  // resulting parameter file will look like:
  //
  // @code
  // # Example parameter file
  // subsection first
  //   subsection second
  //     subsection third
  //       subsection My Class
  //         ... # all the parameters of MyClass
  //         ... # notice My Class subsection does not end here
  //         subsection My Other Class
  //           ... # all the parameters of MyOtherClass
  //         end # of subsection My Other Class
  //       end # of subsection My Class
  //     end
  //   end
  // end
  // @endcode
  //
  // We are going to exploit this, by making our
  // `Parameters` the *parent* of all subsequently
  // constructed classes. Since most of the other classes are members of
  // `DistributedLagrangeProblem` this allows, for example, to construct two
  // `DistributedLagrangeProblem` for two different dimensions, without having
  // conflicts in the parameters for the two problems.
  template <int dim, int spacedim>
  DistributedLagrangeProblem<dim, spacedim>::Parameters::Parameters()
    : ParameterAcceptor("/Distributed Lagrange<" +
                        Utilities::int_to_string(dim) + "," +
                        Utilities::int_to_string(spacedim) + ">/")
  {
    // The ParameterAcceptor::add_parameter() function does a few things:
    //
    // - enters the subsection specified at construction time to
    // ParameterAcceptor
    //
    // - calls the ParameterAcceptor::prm.add_parameter() function
    //
    // - calls any signal you may have attached to
    // ParameterAcceptor::declare_parameters_call_back
    //
    // - leaves the subsection
    //
    // In turn, ParameterAcceptor::prm.add_parameter
    //
    // - declares an entry in the parameter handler for the given variable;
    //
    // - takes the current value of the variable
    //
    // - transforms it to a string, used as the default value for the parameter
    // file
    //
    // - attaches an *action* to ParameterAcceptor::prm that monitors when a
    // file is parsed, or when an entry is set, and when this happens, it
    // updates the value of the variable passed to `add_parameter()` by setting
    // it to whatever was specified in the input file (of course, after the
    // input file has been parsed and the text representation converted to the
    // type of the variable).
    add_parameter("Initial embedding space refinement", initial_refinement);

    add_parameter("Initial embedded space refinement",
                  initial_embedded_refinement);

    add_parameter("Local refinements steps near embedded domain",
                  delta_refinement);

    add_parameter("Homogeneous Dirichlet boundary ids",
                  homogeneous_dirichlet_ids);

    add_parameter("Use displacement in embedded interface", use_displacement);

    add_parameter("Embedding space finite element degree",
                  embedding_space_finite_element_degree);

    add_parameter("Embedded space finite element degree",
                  embedded_space_finite_element_degree);

    add_parameter("Embedded configuration finite element degree",
                  embedded_configuration_finite_element_degree);

    add_parameter("Coupling quadrature order", coupling_quadrature_order);

    add_parameter("Verbosity level", verbosity_level);

    // Once the parameter file has been parsed, then the parameters are good to
    // go. Set the internal variable `initialized` to true.
    parse_parameters_call_back.connect([&]() -> void { initialized = true; });
  }

  // The constructor is pretty standard, with the exception of the
  // `ParameterAcceptorProxy` objects, as explained earlier.
  template <int dim, int spacedim>
  DistributedLagrangeProblem<dim, spacedim>::DistributedLagrangeProblem(
    const Parameters &parameters)
    : parameters(parameters)
    , embedded_configuration_function("Embedded configuration", spacedim)
    , embedded_value_function("Embedded value")
    , schur_solver_control("Schur solver control")
    , monitor(std::cout, TimerOutput::summary, TimerOutput::cpu_and_wall_times)
  {
    // Here is a way to set default values for a ParameterAcceptor class
    // that was constructed using ParameterAcceptorProxy.
    //
    // In this case, we set the default deformation of the embedded grid to be a
    // circle with radius $R$ and center $(Cx, Cy)$, we set the default value
    // for the embedded_value_function to be the constant one, and specify some
    // sensible values for the SolverControl object.
    //
    // It is fundamental for $\Gamma$ to be embedded: from the definition of
    // $C_{\alpha j}$ is clear that, if $\Gamma \not\subseteq \Omega$, certain
    // rows of the matrix $C$ will be zero. This would be a problem, as the
    // Schur complement method requires $C$ to have full column rank.
    embedded_configuration_function.declare_parameters_call_back.connect(
      []() -> void {
        ParameterAcceptor::prm.set("Function constants", "R=.3, Cx=.4, Cy=.4");


        ParameterAcceptor::prm.set("Function expression",
                                   "R*cos(2*pi*x)+Cx; R*sin(2*pi*x)+Cy");
      });

    embedded_value_function.declare_parameters_call_back.connect(
      []() -> void { ParameterAcceptor::prm.set("Function expression", "1"); });

    schur_solver_control.declare_parameters_call_back.connect([]() -> void {
      ParameterAcceptor::prm.set("Max steps", "1000");
      ParameterAcceptor::prm.set("Reduction", "1.e-12");
      ParameterAcceptor::prm.set("Tolerance", "1.e-12");
    });
  }

  // @sect3{Set up}
  //
  // The function `DistributedLagrangeProblem::setup_grids_and_dofs()` is used
  // to set up the finite element spaces. Notice how `std::make_unique` is
  // used to create objects wrapped inside `std::unique_ptr` objects.
  template <int dim, int spacedim>
  void DistributedLagrangeProblem<dim, spacedim>::setup_grids_and_dofs()
  {
    TimerOutput::Scope timer_section(monitor, "Setup grids and dofs");

    // Initializing $\Omega$: constructing the Triangulation and wrapping it
    // into a `std::unique_ptr` object
    space_grid = std::make_unique<Triangulation<spacedim>>();

    // Next, we actually create the triangulation using
    // GridGenerator::hyper_cube(). The last argument is set to true: this
    // activates colorization (i.e., assigning different boundary indicators to
    // different parts of the boundary), which we use to assign the Dirichlet
    // and Neumann conditions.
    GridGenerator::hyper_cube(*space_grid, 0, 1, true);

    // Once we constructed a Triangulation, we refine it globally according to
    // the specifications in the parameter file, and construct a
    // GridTools::Cache with it.
    space_grid->refine_global(parameters.initial_refinement);
    space_grid_tools_cache =
      std::make_unique<GridTools::Cache<spacedim, spacedim>>(*space_grid);

    // The same is done with the embedded grid. Since the embedded grid is
    // deformed, we first need to setup the deformation mapping. We do so in the
    // following few lines:
    embedded_grid = std::make_unique<Triangulation<dim, spacedim>>();
    GridGenerator::hyper_cube(*embedded_grid);
    embedded_grid->refine_global(parameters.initial_embedded_refinement);

    embedded_configuration_fe = std::make_unique<FESystem<dim, spacedim>>(
      FE_Q<dim, spacedim>(
        parameters.embedded_configuration_finite_element_degree),
      spacedim);

    embedded_configuration_dh =
      std::make_unique<DoFHandler<dim, spacedim>>(*embedded_grid);

    embedded_configuration_dh->distribute_dofs(*embedded_configuration_fe);
    embedded_configuration.reinit(embedded_configuration_dh->n_dofs());

    // Once we have defined a finite dimensional space for the deformation, we
    // interpolate the `embedded_configuration_function` defined in the
    // parameter file:
    VectorTools::interpolate(*embedded_configuration_dh,
                             embedded_configuration_function,
                             embedded_configuration);

    // Now we can interpret it according to what the user has specified in the
    // parameter file: as a displacement, in which case we construct a mapping
    // that *displaces* the position of each support point of our configuration
    // finite element space by the specified amount on the corresponding
    // configuration vector, or as an absolution position.
    //
    // In the first case, the class MappingQEulerian offers its services, while
    // in the second one, we'll use the class MappingFEField. They are in fact
    // very similar. MappingQEulerian will only work for systems of FE_Q finite
    // element spaces, where the displacement vector is stored in the first
    // `spacedim` components of the FESystem, and the degree given as a
    // parameter at construction time, must match the degree of the first
    // `spacedim` components.
    //
    // The class MappingFEField is slightly more general, in that it allows you
    // to select arbitrary FiniteElement types when constructing your
    // approximation. Naturally some choices may (or may not) make sense,
    // according to the type of FiniteElement you choose. MappingFEField
    // implements the pure iso-parametric concept, and can be used, for example,
    // to implement iso-geometric analysis codes in deal.II, by combining it
    // with the FE_Bernstein finite element class. In this example, we'll use
    // the two interchangeably, by taking into account the fact that one
    // configuration will be a `displacement`, while the other will be an
    // absolute `deformation` field.

    if (parameters.use_displacement == true)
      embedded_mapping =
        std::make_unique<MappingQEulerian<dim, Vector<double>, spacedim>>(
          parameters.embedded_configuration_finite_element_degree,
          *embedded_configuration_dh,
          embedded_configuration);
    else
      embedded_mapping =
        std::make_unique<MappingFEField<dim,
                                        spacedim,
                                        Vector<double>,
                                        DoFHandler<dim, spacedim>>>(
          *embedded_configuration_dh, embedded_configuration);

    setup_embedded_dofs();

    // In this tutorial program we not only refine $\Omega$ globally,
    // but also allow a local refinement depending on the position of $\Gamma$,
    // according to the value of `parameters.delta_refinement`, that we use to
    // decide how many rounds of local refinement we should do on $\Omega$,
    // corresponding to the position of $\Gamma$.
    //
    // With the mapping in place, it is now possible to query what is the
    // location of all support points associated with the `embedded_dh`, by
    // calling the method DoFTools::map_dofs_to_support_points.
    //
    // This method has two variants. One that does *not* take a Mapping, and
    // one that takes a Mapping. If you use the second type, like we are doing
    // in this case, the support points are computed through the specified
    // mapping, which can manipulate them accordingly.
    //
    // This is precisely what the `embedded_mapping` is there for.
    std::vector<Point<spacedim>> support_points(embedded_dh->n_dofs());
    if (parameters.delta_refinement != 0)
      DoFTools::map_dofs_to_support_points(*embedded_mapping,
                                           *embedded_dh,
                                           support_points);

    // Once we have the support points of the embedded finite element space, we
    // would like to identify what cells of the embedding space contain what
    // support point, to get a chance at refining the embedding grid where it is
    // necessary, i.e., where the embedded grid is. This can be done manually,
    // by looping over each support point, and then calling the method
    // Mapping::transform_real_to_unit_cell for each cell of the embedding
    // space, until we find one that returns points in the unit reference cell,
    // or it can be done in a more intelligent way.
    //
    // The GridTools::find_active_cell_around_point is a possible option that
    // performs the above task in a cheaper way, by first identifying the
    // closest vertex of the embedding Triangulation to the target point, and
    // then by calling Mapping::transform_real_to_unit_cell only for those cells
    // that share the found vertex.
    //
    // In fact, there are algorithms in the GridTools namespace that exploit a
    // GridTools::Cache object, and possibly a KDTree object to speed up these
    // operations as much as possible.
    //
    // The simplest way to exploit the maximum speed is by calling a
    // specialized method, GridTools::compute_point_locations, that will store a
    // lot of useful information and data structures during the first point
    // search, and then reuse all of this for subsequent points.
    //
    // GridTools::compute_point_locations returns a tuple where the first
    // element is a vector of cells containing the input points, in this
    // case support_points. For refinement, this is the only information we
    // need, and this is exactly what happens now.
    //
    // When we need to assemble a coupling matrix, however, we'll also need the
    // reference location of each point to evaluate the basis functions of the
    // embedding space. The other elements of the tuple returned by
    // GridTools::compute_point_locations allow you to reconstruct, for each
    // point, what cell contains it, and what is the location in the reference
    // cell of the given point. Since this information is better grouped into
    // cells, then this is what the algorithm returns: a tuple, containing a
    // vector of all cells that have at least one point in them, together with a
    // list of all reference points and their corresponding index in the
    // original vector.
    //
    // In the following loop, we will be ignoring all returned objects except
    // the first, identifying all cells contain at least one support point of
    // the embedded space. This allows for a simple adaptive refinement
    // strategy: refining these cells and their neighbors.
    //
    // Notice that we need to do some sanity checks, in the sense that we want
    // to have an embedding grid which is well refined around the embedded grid,
    // but where two consecutive support points lie either in the same cell, or
    // in neighbor embedding cells.
    //
    // This is only possible if we ensure that the smallest cell size of the
    // embedding grid is nonetheless bigger than the largest cell size of the
    // embedded grid. Since users can modify both levels of refinements, as well
    // as the amount of local refinement they want around the embedded grid, we
    // make sure that the resulting meshes satisfy our requirements, and if this
    // is not the case, we bail out with an exception.
    for (unsigned int i = 0; i < parameters.delta_refinement; ++i)
      {
        const auto point_locations =
          GridTools::compute_point_locations(*space_grid_tools_cache,
                                             support_points);
        const auto &cells = std::get<0>(point_locations);
        for (auto &cell : cells)
          {
            cell->set_refine_flag();
            for (unsigned int face_no : GeometryInfo<spacedim>::face_indices())
              if (!cell->at_boundary(face_no))
                cell->neighbor(face_no)->set_refine_flag();
          }
        space_grid->execute_coarsening_and_refinement();
      }

    // In order to construct a well posed coupling interpolation operator $C$,
    // there are some constraints on the relative dimension of the grids between
    // the embedding and the embedded domains. The coupling operator $C$ and the
    // spaces $V$ and $Q$ have to satisfy an inf-sup condition in order for the
    // problem to have a solution. It turns out that the non-matching $L^2$
    // projection satisfies such inf-sup, provided that the spaces $V$ and $Q$
    // are compatible between each other (for example, provided that they are
    // chosen to be the ones described in the introduction).
    //
    // However, the *discrete* inf-sup condition must also hold. No
    // complications arise here, but it turns out that the discrete inf-sup
    // constant deteriorates when the non-matching grids have local diameters
    // that are too far away from each other. In particular, it turns out that
    // if you choose an embedding grid which is *finer* with respect to the
    // embedded grid, the inf-sup constant deteriorates much more than if you
    // let the embedded grid be finer.
    //
    // In order to avoid issues, in this tutorial we will throw an exception if
    // the parameters chosen by the user are such that the maximal diameter of
    // the embedded grid is greater than the minimal diameter of the embedding
    // grid.
    //
    // This choice guarantees that almost every cell of the embedded grid spans
    // no more than two cells of the embedding grid, with some rare exceptions,
    // that are negligible in terms of the resulting inf-sup.
    const double embedded_space_maximal_diameter =
      GridTools::maximal_cell_diameter(*embedded_grid, *embedded_mapping);
    double embedding_space_minimal_diameter =
      GridTools::minimal_cell_diameter(*space_grid);

    deallog << "Embedding minimal diameter: "
            << embedding_space_minimal_diameter
            << ", embedded maximal diameter: "
            << embedded_space_maximal_diameter << ", ratio: "
            << embedded_space_maximal_diameter /
                 embedding_space_minimal_diameter
            << std::endl;

    AssertThrow(embedded_space_maximal_diameter <
                  embedding_space_minimal_diameter,
                ExcMessage(
                  "The embedding grid is too refined (or the embedded grid "
                  "is too coarse). Adjust the parameters so that the minimal "
                  "grid size of the embedding grid is larger "
                  "than the maximal grid size of the embedded grid."));

    // $\Omega$ has been refined and we can now set up its DoFs
    setup_embedding_dofs();
  }

  // We now set up the DoFs of $\Omega$ and $\Gamma$: since they are
  // fundamentally independent (except for the fact that $\Omega$'s mesh is more
  // refined "around"
  // $\Gamma$) the procedure is standard.
  template <int dim, int spacedim>
  void DistributedLagrangeProblem<dim, spacedim>::setup_embedding_dofs()
  {
    space_dh = std::make_unique<DoFHandler<spacedim>>(*space_grid);
    space_fe = std::make_unique<FE_Q<spacedim>>(
      parameters.embedding_space_finite_element_degree);
    space_dh->distribute_dofs(*space_fe);

    DoFTools::make_hanging_node_constraints(*space_dh, constraints);
    for (auto id : parameters.homogeneous_dirichlet_ids)
      {
        VectorTools::interpolate_boundary_values(
          *space_dh, id, Functions::ZeroFunction<spacedim>(), constraints);
      }
    constraints.close();

    // By definition the stiffness matrix involves only $\Omega$'s DoFs
    DynamicSparsityPattern dsp(space_dh->n_dofs(), space_dh->n_dofs());
    DoFTools::make_sparsity_pattern(*space_dh, dsp, constraints);
    stiffness_sparsity.copy_from(dsp);
    stiffness_matrix.reinit(stiffness_sparsity);
    solution.reinit(space_dh->n_dofs());
    rhs.reinit(space_dh->n_dofs());

    deallog << "Embedding dofs: " << space_dh->n_dofs() << std::endl;
  }

  template <int dim, int spacedim>
  void DistributedLagrangeProblem<dim, spacedim>::setup_embedded_dofs()
  {
    embedded_dh = std::make_unique<DoFHandler<dim, spacedim>>(*embedded_grid);
    embedded_fe = std::make_unique<FE_Q<dim, spacedim>>(
      parameters.embedded_space_finite_element_degree);
    embedded_dh->distribute_dofs(*embedded_fe);

    // By definition the rhs of the system we're solving involves only a zero
    // vector and $G$, which is computed using only $\Gamma$'s DoFs
    lambda.reinit(embedded_dh->n_dofs());
    embedded_rhs.reinit(embedded_dh->n_dofs());
    embedded_value.reinit(embedded_dh->n_dofs());

    deallog << "Embedded dofs: " << embedded_dh->n_dofs() << std::endl;
  }

  // Creating the coupling sparsity pattern is a complex operation,
  // but it can be easily done using the
  // NonMatching::create_coupling_sparsity_pattern, which requires the
  // two DoFHandler objects, the quadrature points for the coupling,
  // a DynamicSparsityPattern (which then needs to be copied into the
  // sparsity one, as usual), the component mask for the embedding and
  // embedded Triangulation (which we leave empty) and the mappings
  // for both the embedding and the embedded Triangulation.
  template <int dim, int spacedim>
  void DistributedLagrangeProblem<dim, spacedim>::setup_coupling()
  {
    TimerOutput::Scope timer_section(monitor, "Setup coupling");

    QGauss<dim> quad(parameters.coupling_quadrature_order);

    DynamicSparsityPattern dsp(space_dh->n_dofs(), embedded_dh->n_dofs());

    NonMatching::create_coupling_sparsity_pattern(*space_grid_tools_cache,
                                                  *space_dh,
                                                  *embedded_dh,
                                                  quad,
                                                  dsp,
                                                  AffineConstraints<double>(),
                                                  ComponentMask(),
                                                  ComponentMask(),
                                                  *embedded_mapping);
    coupling_sparsity.copy_from(dsp);
    coupling_matrix.reinit(coupling_sparsity);
  }

  // @sect3{Assembly}
  //
  // The following function creates the matrices: as noted before computing the
  // stiffness matrix and the rhs is a standard procedure.
  template <int dim, int spacedim>
  void DistributedLagrangeProblem<dim, spacedim>::assemble_system()
  {
    {
      TimerOutput::Scope timer_section(monitor, "Assemble system");

      // Embedding stiffness matrix $K$, and the right hand side $G$.
      MatrixTools::create_laplace_matrix(
        *space_dh,
        QGauss<spacedim>(2 * space_fe->degree + 1),
        stiffness_matrix,
        static_cast<const Function<spacedim> *>(nullptr),
        constraints);

      VectorTools::create_right_hand_side(*embedded_mapping,
                                          *embedded_dh,
                                          QGauss<dim>(2 * embedded_fe->degree +
                                                      1),
                                          embedded_value_function,
                                          embedded_rhs);
    }
    {
      TimerOutput::Scope timer_section(monitor, "Assemble coupling system");

      // To compute the coupling matrix we use the
      // NonMatching::create_coupling_mass_matrix tool, which works similarly to
      // NonMatching::create_coupling_sparsity_pattern.
      QGauss<dim> quad(parameters.coupling_quadrature_order);
      NonMatching::create_coupling_mass_matrix(*space_grid_tools_cache,
                                               *space_dh,
                                               *embedded_dh,
                                               quad,
                                               coupling_matrix,
                                               AffineConstraints<double>(),
                                               ComponentMask(),
                                               ComponentMask(),
                                               *embedded_mapping);

      VectorTools::interpolate(*embedded_mapping,
                               *embedded_dh,
                               embedded_value_function,
                               embedded_value);
    }
  }

  // @sect3{Solve}
  //
  // All parts have been assembled: we solve the system
  // using the Schur complement method
  template <int dim, int spacedim>
  void DistributedLagrangeProblem<dim, spacedim>::solve()
  {
    TimerOutput::Scope timer_section(monitor, "Solve system");

    // Start by creating the inverse stiffness matrix
    SparseDirectUMFPACK K_inv_umfpack;
    K_inv_umfpack.initialize(stiffness_matrix);

    // Initializing the operators, as described in the introduction
    auto K  = linear_operator(stiffness_matrix);
    auto Ct = linear_operator(coupling_matrix);
    auto C  = transpose_operator(Ct);

    auto K_inv = linear_operator(K, K_inv_umfpack);

    // Using the Schur complement method
    auto                     S = C * K_inv * Ct;
    SolverCG<Vector<double>> solver_cg(schur_solver_control);
    auto S_inv = inverse_operator(S, solver_cg, PreconditionIdentity());

    lambda = S_inv * embedded_rhs;

    solution = K_inv * Ct * lambda;

    constraints.distribute(solution);
  }

  // The following function simply generates standard result output on two
  // separate files, one for each mesh.
  template <int dim, int spacedim>
  void DistributedLagrangeProblem<dim, spacedim>::output_results()
  {
    TimerOutput::Scope timer_section(monitor, "Output results");

    DataOut<spacedim> embedding_out;

    std::ofstream embedding_out_file("embedding.vtu");

    embedding_out.attach_dof_handler(*space_dh);
    embedding_out.add_data_vector(solution, "solution");
    embedding_out.build_patches(
      parameters.embedding_space_finite_element_degree);
    embedding_out.write_vtu(embedding_out_file);

    // The only difference between the two output routines is that in the
    // second case, we want to output the data on the current configuration, and
    // not on the reference one. This is possible by passing the actual
    // embedded_mapping to the DataOut::build_patches function. The mapping will
    // take care of outputting the result on the actual deformed configuration.

    DataOut<dim, DoFHandler<dim, spacedim>> embedded_out;

    std::ofstream embedded_out_file("embedded.vtu");

    embedded_out.attach_dof_handler(*embedded_dh);
    embedded_out.add_data_vector(lambda, "lambda");
    embedded_out.add_data_vector(embedded_value, "g");
    embedded_out.build_patches(*embedded_mapping,
                               parameters.embedded_space_finite_element_degree);
    embedded_out.write_vtu(embedded_out_file);
  }

  // Similar to all other tutorial programs, the `run()` function simply calls
  // all other methods in the correct order. Nothing special to note, except
  // that we check if parsing was done before we actually attempt to run our
  // program.
  template <int dim, int spacedim>
  void DistributedLagrangeProblem<dim, spacedim>::run()
  {
    AssertThrow(parameters.initialized, ExcNotInitialized());
    deallog.depth_console(parameters.verbosity_level);

    setup_grids_and_dofs();
    setup_coupling();
    assemble_system();
    solve();
    output_results();
  }
} // namespace Step60



int main(int argc, char **argv)
{
  try
    {
      using namespace dealii;
      using namespace Step60;

      const unsigned int dim = 1, spacedim = 2;

      // Differently to what happens in other tutorial programs, here we use
      // ParameterAcceptor style of initialization, i.e., all objects are first
      // constructed, and then a single call to the static method
      // ParameterAcceptor::initialize is issued to fill all parameters of the
      // classes that are derived from ParameterAcceptor.
      //
      // We check if the user has specified a parameter file name to use when
      // the program was launched. If so, try to read that parameter file,
      // otherwise, try to read the file "parameters.prm".
      //
      // If the parameter file that was specified (implicitly or explicitly)
      // does not exist, ParameterAcceptor::initialize will create one for you,
      // and exit the program.

      DistributedLagrangeProblem<dim, spacedim>::Parameters parameters;
      DistributedLagrangeProblem<dim, spacedim>             problem(parameters);

      std::string parameter_file;
      if (argc > 1)
        parameter_file = argv[1];
      else
        parameter_file = "parameters.prm";

      ParameterAcceptor::initialize(parameter_file, "used_parameters.prm");
      problem.run();
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

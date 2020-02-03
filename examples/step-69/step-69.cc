/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the deal.II authors
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
 * Authors: Matthias Maier, Texas A&M University;
 *          Ignacio Tomas, Texas A&M University, Sandia National Laboratories
 */

// @sect3{Include files}

// The set of include files is quite standard. The most intriguing part is
// the fact that we will rely solely on deal.II data structures for MPI
// parallelization, in particular distributed::Triangulation and
// LinearAlgebra::distributed::Vector included through
// <code>distributed/tria.h</code> and
// <code>lac/la_parallel_vector.h</code>. Instead of a Trilinos, or PETSc
// specific matrix class, we will use a non-distributed
// dealii::SparseMatrix (<code>lac/sparse_matrix.h</code>) to store the local
// part of the $c_{ij}$, $n_{ij}$ and $d_{ij}$ matrices.
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_matrix.templates.h>
#include <deal.II/lac/vector.h>

#include <deal.II/meshworker/scratch_data.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

// In addition to above deal.II specific includes, we also include four
// boost headers. The first two are for binary archives that we will use
// for implementing a check-pointing and restart mechanism.
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

// The last two boost header files are for creating custom iterator ranges
// over integer intervals.
#include <boost/range/irange.hpp>
#include <boost/range/iterator_range.hpp>

// @sect3{Class template declarations}
//
// We begin our actual implementation by declaring all classes with its
// data structures and methods upfront. In contrast to previous example
// steps we use a more fine-grained encapsulation of concepts, data
// structures, and parameters into individual classes. A single class thus
// usually centers around either a  single data structure (such as the
// Triangulation) in the <code>Discretization</code> class, or a single
// method (such as the <code>step()</code> function of the
// <code>TimeStep</code> class). We typically declare parameter variables
// and scratch data object private and make methods and data structures
// used by other classes public.
//
// @note: A cleaner approach would be to guard access to all data
// structures by <a
// href="https://en.wikipedia.org/wiki/Mutator_method">getter/setter
// functions</a>. For the sake of brevity, we refrain from that approach,
// though.

namespace Step69
{
  using namespace dealii;

  // We start with an enum describing all possible boundary conditions
  // encountered in this tutorial step. Such an enum allows us to refer to
  // boundary types by a mnemonic (such as
  // <code>Boundary::do_nothing</code>) rather than a numerical value.
  enum Boundary : types::boundary_id
  {
    do_nothing = 0,
    slip       = 1,
    dirichlet  = 2,
  };

  // @sect4{The <code>Discretization</code> class}
  //
  // The class <code>Discretization</code> contains all data structures
  // concerning the mesh (triangulation) and discretization (mapping,
  // finite element, quadrature) of the problem. We use the
  // ParameterAcceptor class to automatically populate problem-specific
  // parameters, such as the geometry information
  // (<code>length</code>, etc.) or the refinement level
  // (<code>refinement</code>) from a parameter file. This requires us to
  // split the initialization of data structures into two functions: We
  // initialize everything that does not depend on parameters in the
  // constructor, and defer the creation of the mesh to the
  // <code>setup()</code> method that can be called once all parameters are
  // read-in via ParameterAcceptor::initialize().
  //
  template <int dim>
  class Discretization : public ParameterAcceptor
  {
  public:
    Discretization(const MPI_Comm &   mpi_communicator,
                   TimerOutput &      computing_timer,
                   const std::string &subsection = "Discretization");

    void setup();

    const MPI_Comm &mpi_communicator;

    parallel::distributed::Triangulation<dim> triangulation;

    const MappingQ<dim>   mapping;
    const FE_Q<dim>       finite_element;
    const QGauss<dim>     quadrature;
    const QGauss<dim - 1> face_quadrature;

  private:
    TimerOutput &computing_timer;

    double length;
    double height;
    double disc_position;
    double disc_diameter;

    unsigned int refinement;
  };

  // @sect4{The <code>OfflineData</code> class}
  //
  // The class <code>OfflineData</code> contains pretty much all components
  // of the discretization that do not evolve in time, in particular, the
  // DoFHandler, SparsityPattern, boundary maps, the lumped mass, $c_{ij}$,
  // and $n_{ij}$ matrices.
  //
  // Here, the term <i>offline</i> refers to the fact that all the class
  // members of <code>OfflineData</code> have well-defined values
  // independent of the current time step. This means that they can be
  // initialized ahead of time (at <i>time step zero</i>) and are not meant
  // to be modified at any other later time step. For instance, the
  // sparsity pattern should not change as we advance in time (we are not
  // doing any form of adaptivity in space). Similarly, the entries of the
  // lumped mass matrix should not be modified as we advance in time
  // either.
  //
  // We also compute and store a <code>boundary_normal_map</code> that
  // contains a map from a global index of type `types:global_dof_index` of
  // a boundary degree of freedom to a tuple consisting of a normal vector,
  // the boundary id, and the position associated with the degree of
  // freedom. We actually have to compute and store this geometric
  // information in this class because we won't have access to geometric
  // (or cell-based) information later on in the algebraic loops over the
  // sparsity pattern.
  //
  // @note Even though this class currently does not have any parameters
  // that could be read in from a parameter file we nevertheless derive
  // from ParameterAcceptor and follow the same idiom of providing a
  // <code>setup()</code> (and <code>assemble()</code>) method as for the
  // class Discretization.

  template <int dim>
  class OfflineData : public ParameterAcceptor
  {
  public:
    using BoundaryNormalMap =
      std::map<types::global_dof_index,
               std::tuple<Tensor<1, dim>, types::boundary_id, Point<dim>>>;

    OfflineData(const MPI_Comm &           mpi_communicator,
                TimerOutput &              computing_timer,
                const Discretization<dim> &discretization,
                const std::string &        subsection = "OfflineData");

    void setup();
    void assemble();

    DoFHandler<dim> dof_handler;

    std::shared_ptr<const Utilities::MPI::Partitioner> partitioner;

    unsigned int n_locally_owned;
    unsigned int n_locally_relevant;

    SparsityPattern sparsity_pattern;

    BoundaryNormalMap boundary_normal_map;

    SparseMatrix<double>                  lumped_mass_matrix;
    std::array<SparseMatrix<double>, dim> cij_matrix;
    std::array<SparseMatrix<double>, dim> nij_matrix;
    SparseMatrix<double>                  norm_matrix;

  private:
    const MPI_Comm &mpi_communicator;
    TimerOutput &   computing_timer;

    SmartPointer<const Discretization<dim>> discretization;
  };

  // @sect4{The <code>ProblemDescription</code> class}
  //
  // The member functions of this class are utility functions specific to
  // Euler's equations:
  // - The type alias <code>rank1_type</code> is used for the states
  //   $\mathbf{U}_i^n$
  // - The type alias <code>rank2_type</code> is used for the fluxes
  //   $\mathbb{f}(\mathbf{U}_j^n)$.
  // - The <code>momentum</code> function extracts $\textbf{m}$
  //   out of the state vector $[\rho,\textbf{m},E]$ and stores it in a
  //   <code>Tensor<1, dim></code>.
  // - The <code>internal_energy</code> function computes $E -
  //   \frac{|\textbf{m}|^2}{2\rho}$ from a given state vector
  //   $[\rho,\textbf{m},E]$.
  //
  // The purpose of the class members <code>component_names</code>,
  // <code>pressure</code>, and <code>speed_of_sound</code>, is evident
  // from their names. We also provide a function
  // <code>compute_lambda_max</code>, that computes the wave speed estimate
  // mentioned above, $\lambda_{max}(\mathbf{U},\mathbf{V},\mathbf{n})$,
  // which is used in the computation of the $d_{ij}$ matrix.
  //
  // @note The <code>DEAL_II_ALWAYS_INLINE</code> macro expands to a
  // (compiler specific) pragma that ensures that the corresponding
  // function defined in this class is always inlined, i.e., the function
  // body is put in place for every invocation of the function, and no call
  // (and code indirection) is generated. This is stronger than the
  // <code>inline</code> keyword, which is more or less a (mild) suggestion
  // to the compiler that the programmer things it would be beneficial to
  // inline the function. <code>DEAL_II_ALWAYS_INLINE</code> should only be
  // used rarely and with caution in situations such as this one, where we
  // actually know (due to benchmarking) that inlining the function in
  // question actually improves performance.
  // 
  // Finally we note that:
  //  - This is the only class in this tutorial step that is tied to a
  //    particular "physics" or "hyperbolic conservation law" (in this
  //    case Euler's equations). All the other classes are primarily
  //    "discretization" classes, very much agnostic of the particular physics
  //    being solved.
  //  - This is a "pure static" class (the antithesis of a
  //    "pure virtual" class). It's just a convenient way to wrap-up a
  //    collection of related methods into a single object. Note that we will 
  //    be able to invoke such methods without without creating an instance of 
  //    the class. Similarly, we will not have to provide a constructor 
  //    for this class.

  template <int dim>
  class ProblemDescription
  {
  public:

    /* constexpr tells the compiler to evaluate "2 + dim" just once at compile 
       time rather than everytime problem_dimension is invoked. */
    static constexpr unsigned int problem_dimension = 2 + dim;

    using rank1_type = Tensor<1, problem_dimension>;
    using rank2_type = Tensor<1, problem_dimension, Tensor<1, dim>>;

    const static std::array<std::string, dim + 2> component_names;

    static constexpr double gamma = 7. / 5.;

    static DEAL_II_ALWAYS_INLINE inline Tensor<1, dim>
    momentum(const rank1_type &U);

    static DEAL_II_ALWAYS_INLINE inline double
    internal_energy(const rank1_type &U);

    static DEAL_II_ALWAYS_INLINE inline double pressure(const rank1_type &U);

    static DEAL_II_ALWAYS_INLINE inline double
    speed_of_sound(const rank1_type &U);

    static DEAL_II_ALWAYS_INLINE inline rank2_type f(const rank1_type &U);

    static DEAL_II_ALWAYS_INLINE inline double
    compute_lambda_max(const rank1_type &    U_i,
                       const rank1_type &    U_j,
                       const Tensor<1, dim> &n_ij);
  };

  // @sect4{The <code>InitialValues</code> class}
  //
  // The class <code>InitialValues</code>'s only public data attribute is a
  // std::function <code>initial_state</code> that computes the initial
  // state of a given point and time. For the purpose of this example step
  // we simply implement a homogeneous uniform flow field for which the
  // direction and a 1D primitive state (density, velocity, pressure) are
  // read from the parameter file.
  //
  // It would be desirable to initialize the class in a single shot:
  // initialize/set the parameters and define the class members that 
  // depend on these default parameters. However, since we do not know the 
  // actual final values for the parameters, this would be sort of 
  // meaningless an unsafe in general (we would like to have mechanisms to 
  // check the consistency of the input parameters). Instead of defining 
  // another <code>setup()</code> method to be called (by-hand) after the 
  // call to <code> ParameterAcceptor::initialize() </code> we provide an 
  // "implementation" for the class member 
  // <code>parse_parameters_call_back</code> which is automatically called when
  // invoking <code> ParameterAcceptor::initialize() </code> for every class 
  // that inherits from ParameterAceptor.

  template <int dim>
  class InitialValues : public ParameterAcceptor
  {
  public:
    using rank1_type = typename ProblemDescription<dim>::rank1_type;

    InitialValues(const std::string &subsection = "InitialValues");

    std::function<rank1_type(const Point<dim> &point, double t)> initial_state;

  private:

    /* Auxiliary void function to be hooked to the inherited class member
       ParameterAcceptor::parse_parameters_call_back. */
    void parse_parameters_callback();

    Tensor<1, dim> initial_direction;
    Tensor<1, 3>   initial_1d_state;
  };

  // @sect4{The <code>TimeStep</code> class}
  //
  // With the <code>OfflineData</code> and <code>ProblemDescription</code>
  // classes at hand we can now implement the explicit time-stepping scheme
  // that was introduced in the discussion above. The main method of the
  // <code>TimeStep</code> class is <code>step(vector_type &U, double
  // t)</code>. That takes a reference to a state vector <code>U</code> and
  // a time point <code>t</code> as arguments, computes the updated solution, 
  // stores it in the vector <code>temp</code>, swaps its contents with the
  // vector <code>U</code>, and returns the chosen step-size $\tau$.
  //
  // The other important method is <code>prepare()</code> which primarily sets
  // the proper partition and sparsity pattern for the auxiliary vector 
  // <code>temp</code> and the matrix <code>dij_matrix</code>.
  //

  template <int dim>
  class TimeStep : public ParameterAcceptor
  {
  public:
    static constexpr unsigned int problem_dimension =
      ProblemDescription<dim>::problem_dimension;

    using rank1_type = typename ProblemDescription<dim>::rank1_type;
    using rank2_type = typename ProblemDescription<dim>::rank2_type;

    typedef std::array<LinearAlgebra::distributed::Vector<double>,
                       problem_dimension>
      vector_type;

    TimeStep(const MPI_Comm &          mpi_communicator,
             TimerOutput &             computing_timer,
             const OfflineData<dim> &  offline_data,
             const InitialValues<dim> &initial_values,
             const std::string &       subsection = "TimeStep");

    void prepare();

    double step(vector_type &U, double t);

  private:
    const MPI_Comm &mpi_communicator;
    TimerOutput &   computing_timer;

    SmartPointer<const OfflineData<dim>>   offline_data;
    SmartPointer<const InitialValues<dim>> initial_values;

    SparseMatrix<double> dij_matrix;

    vector_type temp;

    double cfl_update;
  };

  // @sect4{The <code>SchlierenPostprocessor</code> class}
  //
  // At its core, the Schlieren class implements the class member
  // <code>compute_schlieren</code>. The main purpose of this class member
  // is to compute an auxiliary finite element field
  // <code>schlieren</code>, that is defined at each node by
  // \f[ \text{schlieren}[i] = e^{\beta \frac{ |\nabla r_i|
  // - \min_j |\nabla r_j| }{\max_j |\nabla r_j| - \min_j |\nabla r_j| } }, \f]
  // where $r$ can in principle be any scalar quantitiy, in practice
  // though, the density is a natural candidate, viz. $r := \rho$.
  // Schlieren postprocessing is a standard method for enhancing the
  // contrast of a visualization inspired by actual experimental X-ray and
  // shadowgraphy techniques of visualization.

  template <int dim>
  class SchlierenPostprocessor : public ParameterAcceptor
  {
  public:
    static constexpr unsigned int problem_dimension =
      ProblemDescription<dim>::problem_dimension;

    using rank1_type = typename ProblemDescription<dim>::rank1_type;

    using vector_type =
      std::array<LinearAlgebra::distributed::Vector<double>, problem_dimension>;

    SchlierenPostprocessor(
      const MPI_Comm &        mpi_communicator,
      TimerOutput &           computing_timer,
      const OfflineData<dim> &offline_data,
      const std::string &     subsection = "SchlierenPostprocessor");

    void prepare();

    void compute_schlieren(const vector_type &U);

    LinearAlgebra::distributed::Vector<double> schlieren;

  private:
    const MPI_Comm &mpi_communicator;
    TimerOutput &   computing_timer;

    SmartPointer<const OfflineData<dim>> offline_data;

    Vector<double> r;

    unsigned int schlieren_index;
    double       schlieren_beta;
  };

  // @sect4{The <code>TimeLoop</code> class}
  //
  // Now, all that is left to do is to chain the methods implemented in the
  // <code>TimeStep</code>, <code>InitialValues</code>, and
  // <code>SchlierenPostprocessor</code> classes together. We do this in a
  // separate class <code>TimeLoop</code> that contains an object of every
  // class and again reads in a number of parameters with the help of the
  // ParameterAcceptor class.

  template <int dim>
  class TimeLoop : public ParameterAcceptor
  {
  public:
    using vector_type = typename TimeStep<dim>::vector_type;

    TimeLoop(const MPI_Comm &mpi_comm);

    void run();

  private:
    vector_type interpolate_initial_values(double t = 0);

    void output(const vector_type &U,
                const std::string &name,
                double             t,
                unsigned int       cycle,
                bool               checkpoint = false);

    const MPI_Comm &   mpi_communicator;
    std::ostringstream timer_output;
    TimerOutput        computing_timer;

    ConditionalOStream pcout;

    std::string base_name;
    double      t_final;
    double      output_granularity;
    bool        enable_compute_error;

    bool resume;

    Discretization<dim>         discretization;
    OfflineData<dim>            offline_data;
    InitialValues<dim>          initial_values;
    TimeStep<dim>               time_step;
    SchlierenPostprocessor<dim> schlieren_postprocessor;

    std::unique_ptr<std::ofstream> filestream;

    std::thread output_thread;
    vector_type output_vector;
  };

  // @sect3{Implementation}

  // @sect4{Grid generation, setup of data structures}

  // The first major task at hand is the typical triplet of grid
  // generation, setup of data structures, and assembly. A notable novelty
  // in this example step is the use of the ParameterAcceptor class that we
  // use to populate parameter values: We first initialize the
  // ParameterAcceptor class by calling its constructor with a string
  // <code>subsection</code> denoting the correct subsection in the
  // parameter file. Then, in the constructor body every parameter value is
  // initialized to a sensible default value and registered with the
  // ParameterAcceptor class with a call to
  // ParameterAcceptor::add_parameter.

  template <int dim>
  Discretization<dim>::Discretization(const MPI_Comm &   mpi_communicator,
                                      TimerOutput &      computing_timer,
                                      const std::string &subsection)
    : ParameterAcceptor(subsection)
    , mpi_communicator(mpi_communicator)
    , triangulation(mpi_communicator)
    , mapping(1)
    , finite_element(1)
    , quadrature(3)
    , face_quadrature(3)
    , computing_timer(computing_timer)
  {
    length = 4.;
    add_parameter("immersed disc - length",
                  length,
                  "Immersed disc: length of computational domain");

    height = 2.;
    add_parameter("immersed disc - height",
                  height,
                  "Immersed disc: height of computational domain");

    disc_position = 0.6;
    add_parameter("immersed disc - object position",
                  disc_position,
                  "Immersed disc: x position of immersed disc center point");

    disc_diameter = 0.5;
    add_parameter("immersed disc - object diameter",
                  disc_diameter,
                  "Immersed disc: diameter of immersed disc");

    refinement = 5;
    add_parameter("initial refinement",
                  refinement,
                  "Initial refinement of the geometry");
  }

  // Note that in the previous constructor we only passed the MPI
  // communicator to the <code>triangulation</code>but we still have not
  // initialized the underlying geometry/mesh. As mentioned earlier, we
  // have to postpone this task to the <code>setup()</code> function that
  // gets called after the ParameterAcceptor::initialize() function has
  // populated all parameter variables with the final values read from the
  // parameter file.
  //
  // The <code>setup()</code> function is the last class member that has to
  // be implemented. It creates the actual triangulation that is a
  // benchmark configuration consisting of a channel with a disc obstacle, see
  // @cite GuermondEtAl2018. We construct the geometry by modifying the
  // mesh generated by GridGenerator::hyper_cube_with_cylindrical_hole().
  // We refer to Step-49, Step-53, and Step-54 for an overview how to
  // create advanced meshes.

  template <int dim>
  void Discretization<dim>::setup()
  {
    TimerOutput::Scope t(computing_timer, "discretization - setup");

    triangulation.clear();

    // We first create 4 temporary (non distributed) coarse triangulations
    // that we stitch together with the
    // GridGenerator::merge_triangulation() function. We center the disc at
    // $(0,0)$ with a diameter of <code>disc_diameter</code>. The lower
    // left corner of the channel has coordinates
    // (<code>-disc_position</code>, <code>-height/2</code>) and the upper
    // right corner has (<code>length-disc_position</code>,
    // <code>height/2</code>).

    Triangulation<dim> tria1, tria2, tria3, tria4;

    GridGenerator::hyper_cube_with_cylindrical_hole(
      tria1, disc_diameter / 2., disc_diameter, 0.5, 1, false);

    GridGenerator::subdivided_hyper_rectangle(
      tria2,
      {2, 1},
      Point<2>(-disc_diameter, disc_diameter),
      Point<2>(disc_diameter, height / 2.));

    GridGenerator::subdivided_hyper_rectangle(
      tria3,
      {2, 1},
      Point<2>(-disc_diameter, -disc_diameter),
      Point<2>(disc_diameter, -height / 2.));

    GridGenerator::subdivided_hyper_rectangle(
      tria4,
      {6, 4},
      Point<2>(disc_diameter, -height / 2.),
      Point<2>(length - disc_position, height / 2.));

    GridGenerator::merge_triangulations({&tria1, &tria2, &tria3, &tria4},
                                        triangulation,
                                        1.e-12,
                                        true);

    triangulation.set_manifold(0, PolarManifold<2>(Point<2>()));

    // We have to fix up the left edge that is currently located at
    // $x=-$<code>disc_diameter</code> and has to be shifted to
    // $x=-$<code>disc_position</code>:

    for (auto cell : triangulation.active_cell_iterators())
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
        {
          auto &vertex = cell->vertex(v);
          if (vertex[0] <= -disc_diameter + 1.e-6)
            vertex[0] = -disc_position;
        }

    // As a last step the boundary has to be colorized with
    // <code>Boundary::do_nothing</code> on the right,
    // <code>Boundary::dirichlet</code> on the left and
    // <code>Boundary::slip</code> on the upper and lower outer boundaries
    // and the obstacle:

    for (auto cell : triangulation.active_cell_iterators())
      {
        for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f)
          {
            const auto face = cell->face(f);

            if (!face->at_boundary())
              continue;

            const auto center = face->center();

            if (center[0] > length - disc_position - 1.e-6)
              face->set_boundary_id(Boundary::do_nothing);
            else if (center[0] < -disc_position + 1.e-6)
              face->set_boundary_id(Boundary::dirichlet);
            else
              face->set_boundary_id(Boundary::slip);
          }
      }

    triangulation.refine_global(refinement);
  }

  // @sect4{Assembly of offline matrices}

  // Not much is done in the constructor of <code>OfflineData</code> other
  // than initializing the corresponding class members in the
  // initialization list.

  template <int dim>
  OfflineData<dim>::OfflineData(const MPI_Comm &           mpi_communicator,
                                TimerOutput &              computing_timer,
                                const Discretization<dim> &discretization,
                                const std::string &        subsection)
    : ParameterAcceptor(subsection)
    , mpi_communicator(mpi_communicator)
    , computing_timer(computing_timer)
    , discretization(&discretization)
  {}

  // Now we can initialize the DoFHandler, extract the IndexSet objects for
  // locally owned and locally relevant DOFs, and initialize a
  // Utilities::MPI::Partitioner object that is needed for distributed
  // vectors.

  template <int dim>
  void OfflineData<dim>::setup()
  {
    IndexSet locally_owned;
    IndexSet locally_relevant;

    {
      TimerOutput::Scope t(computing_timer, "offline_data - distribute dofs");

      dof_handler.initialize(discretization->triangulation,
                             discretization->finite_element);

      DoFRenumbering::Cuthill_McKee(dof_handler);

      locally_owned   = dof_handler.locally_owned_dofs();
      n_locally_owned = locally_owned.n_elements();

      DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant);
      n_locally_relevant = locally_relevant.n_elements();

      partitioner.reset(new Utilities::MPI::Partitioner(locally_owned,
                                                        locally_relevant,
                                                        mpi_communicator));
    }

    const auto dofs_per_cell = discretization->finite_element.dofs_per_cell;

    // @sect4{Translation to local index ranges}

    // We are now in a position to create the sparsity pattern for our
    // matrices. There are quite a few peculiarities that need a detailed
    // explanation. We avoid using a distributed matrix class (as for
    // example provided by Trilinos or PETSc) and instead rely on deal.II's
    // own SparseMatrix object to store the local part of all matrices.
    // This design decision is motivated by the fact that we actually never
    // perform a matrix-vector multiplication. Instead, we will have to
    // perform nonlinear updates while iterating over (the local part) of a
    // connectivity stencil; a task for which deal.II own SparsityPattern
    // is specificially optimized for.
    //
    // This design consideration has a caveat, though. What makes the
    // deal.II SparseMatrix class fast is the <a
    // href="https://en.wikipedia.org/wiki/Sparse_matrix">compressed row
    // storage (CSR)</a> used in the SparsityPattern (see @ref Sparsity).
    // This, unfortunately, does not play nicely with a global distributed
    // index range because a sparsity pattern with CSR cannot contain
    // "holes" in the index range. The distributed matrices offered by
    // deal.II avoid this by translating from a global index range into a
    // contiguous local index range. But this is the precisely the type of
    // index manipulation we want to avoid.
    //
    // Lucky enough, the Utilities::MPI::Partitioner used for distributed
    // vectors provides exactly what we need: It manages a translation from
    // a global index range to a contiguous local (per MPI rank) index
    // range. We therefore simply create a "local" sparsity pattern for the
    // contiguous index range $[0,$<code>n_locally_relevant</code>$)$ and
    // translate between global dof indices and the above local range with
    // the help of the Utilities::MPI::Partitioner::global_to_local()
    // function. All that is left to do is to ensure that we always access
    // elements of a distributed vector by a call to
    // LinearAlgebra::distributed::Vector::local_element(). That way we
    // avoid index translations altogether.

    {
      TimerOutput::Scope t(
        computing_timer,
        "offline_data - create sparsity pattern and set up matrices");

      // We have to create the "local" sparsity pattern by hand. We
      // therefore loop over all locally owned and ghosted cells (see @ref
      // GlossArtificialCell) and extract the (global)
      // <code>dof_indices</code> associated to the cell DOFs and renumber
      // them using <code>partitioner->global_to_local(index)</code>.
      //
      // @note In the case of a locally owned dof, such renumbering consist
      // of applying a shift (i.e. we subtract an offset) such that now they
      // will become a number in the integer interval
      // $[0,$<code>n_locally_owned</code>$)$. However, in the case of a
      // ghosted dof (i.e. not locally owned) the situation is quite
      // different, since the global indices associated to ghosted DOFs will
      // not be (in general) a contiguous set of integers.

      DynamicSparsityPattern dsp(n_locally_relevant, n_locally_relevant);

      std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

      for (auto cell : dof_handler.active_cell_iterators())
        {
          if (cell->is_artificial())
            continue;

          cell->get_dof_indices(dof_indices);
          std::transform(dof_indices.begin(),
                         dof_indices.end(),
                         dof_indices.begin(),
                         [&](auto index) {
                           return partitioner->global_to_local(index);
                         });

          for (const auto dof : dof_indices)
            dsp.add_entries(dof, dof_indices.begin(), dof_indices.end());
        }

      sparsity_pattern.copy_from(dsp);

      lumped_mass_matrix.reinit(sparsity_pattern);
      norm_matrix.reinit(sparsity_pattern);
      for (auto &matrix : cij_matrix)
        matrix.reinit(sparsity_pattern);
      for (auto &matrix : nij_matrix)
        matrix.reinit(sparsity_pattern);
    }
  }

  // This concludes the setup of the DoFHandler and SparseMatrix objects
  // Next, we have to assemble various matrices. We next define a number of
  // helper functions and data structures in an anonymous namespace.

  namespace
  {
    // <code>CopyData</code> class that will be used to assemble the
    // offline data matrices using WorkStream. It acts as a container: it
    // is just a struct where WorkStream stores the local cell
    // contributions. Note that it also contains a class member
    // <code>local_boundary_normal_map</code> used to store the local
    // contributions required to compute the normals at the boundary.

    template <int dim>
    struct CopyData
    {
      bool                                         is_artificial;
      std::vector<types::global_dof_index>         local_dof_indices;
      typename OfflineData<dim>::BoundaryNormalMap local_boundary_normal_map;
      FullMatrix<double>                           cell_lumped_mass_matrix;
      std::array<FullMatrix<double>, dim>          cell_cij_matrix;
    };

    // Next we introduce a number of helper functions that are all
    // concerned about reading and writing matrix and vector entries. They
    // are mainly motivated by providing slightly more efficient code and
    // <a href="https://en.wikipedia.org/wiki/Syntactic_sugar"> syntactic
    // sugar</a> for otherwise somewhat tedious code.

    // The first function we introduce, <code>get_entry</code>, will be
    // used to read the value stored at the entry pointed by a
    // SparsityPattern iterator <code>it</code> of <code>matrix</code>. The
    // function works around a small deficiency in the SparseMatrix
    // interface: The SparsityPattern is concerned with all index
    // operations of the sparse matrix stored in CRS format. As such the
    // iterator already knows the global index of the corresponding matrix
    // entry in the low-level vector stored in the SparseMatrix object. Due
    // to the lack of an interface in the SparseMatrix for accessing the
    // element directly with a SparsityPattern iterator, we unfortunately
    // have to create a temporary SparseMatrix iterator. We simply hide
    // this in the <code>get_entry</code> function.

    template <typename Matrix, typename Iterator>
    DEAL_II_ALWAYS_INLINE inline typename Matrix::value_type
    get_entry(const Matrix &matrix, const Iterator &it)
    {
      const auto                            global_index = it->global_index();
      const typename Matrix::const_iterator matrix_iterator(&matrix,
                                                            global_index);
      return matrix_iterator->value();
    }

    // The <code>set_entry</code> helper is the inverse operation of
    // <code>get_value</code>: Given an iterator and a value, it sets the
    // entry pointed to by the iterator in the matrix.

    template <typename Matrix, typename Iterator>
    DEAL_II_ALWAYS_INLINE inline void
    set_entry(Matrix &                    matrix,
              const Iterator &            it,
              typename Matrix::value_type value)
    {
      const auto                global_index = it->global_index();
      typename Matrix::iterator matrix_iterator(&matrix, global_index);
      matrix_iterator->value() = value;
    }

    // <code>gather_get_entry</code>: we note that $\mathbf{c}_{ij} \in
    // \mathbb{R}^d$. If $d=2$ then $\mathbf{c}_{ij} =
    // [\mathbf{c}_{ij}^1,\mathbf{c}_{ij}^2]^\top$. Which basically implies
    // that we need one matrix per space dimension to store the
    // $\mathbf{c}_{ij}$ vectors. Similar observation follows for the
    // matrix $\mathbf{n}_{ij}$. The purpose of
    // <code>gather_get_entry</code> is to retrieve those entries a store
    // them into a <code>Tensor<1, dim></code> for our convenience.

    template <typename T1, std::size_t k, typename T2>
    DEAL_II_ALWAYS_INLINE inline Tensor<1, k>
    gather_get_entry(const std::array<T1, k> &U, const T2 it)
    {
      Tensor<1, k> result;
      for (unsigned int j = 0; j < k; ++j)
        result[j] = get_entry(U[j], it);
      return result;
    }

    // <code>gather</code> (first interface): this first function
    // signature, having three input arguments, will be used to retrieve
    // the individual components <code>(i,l)</code> of a matrix. The
    // functionality of <code>gather_get_entry</code> and
    // <code>gather</code> is very much the same, but their context is
    // different: the function <code>gather</code> is meant to be used in
    // exceptional/limited number of cases. The reader should be aware that
    // accessing an arbitrary <code>(i,l)</code> entry of a matrix (say for
    // instance Trilinos or PETSc  matrices) is very expensive. Here is
    // where we might want to keep an eye on complexity: we want this
    // operation to have constant complexity (and that's the case of this
    // implementation using deal.ii matrices).

    template <typename T1, std::size_t k, typename T2, typename T3>
    DEAL_II_ALWAYS_INLINE inline Tensor<1, k>
    gather(const std::array<T1, k> &U, const T2 i, const T3 l)
    {
      Tensor<1, k> result;
      for (unsigned int j = 0; j < k; ++j)
        result[j] = U[j](i, l);
      return result;
    }

    // <code>gather</code> (second interface): this second function
    // signature having two input arguments will be used to gather the
    // state at a node <code>i</code> and return <code>Tensor<1,
    // problem_dimension></code> for our convenience.

    template <typename T1, std::size_t k, typename T2>
    DEAL_II_ALWAYS_INLINE inline Tensor<1, k> gather(const std::array<T1, k> &U,
                                                     const T2                 i)
    {
      Tensor<1, k> result;
      for (unsigned int j = 0; j < k; ++j)
        result[j] = U[j].local_element(i);
      return result;
    }

    // <code>scatter</code>: this function has three input arguments, the
    // first one is meant to be a global object (say a locally owned
    // vector), the second argument which could be a
    // <code>Tensor<1,problem_dimension></code>, and the last argument
    // which represents a index of the global object. This function will be
    // primarily used to write the updated nodal values, stored as
    // <code>Tensor<1,problem_dimension></code>, into the globally owned
    // vector.

    template <typename T1, std::size_t k1, typename T2, typename T3>
    DEAL_II_ALWAYS_INLINE inline void
    scatter(std::array<T1, k1> &U, const T2 &result, const T3 i)
    {
      for (unsigned int j = 0; j < k1; ++j)
        U[j].local_element(i) = result[j];
    }
  } // namespace

  // We are now in a position to assemble all matrices stored in
  // <code>OfflineData</code>: the lumped mass entries $m_i$, the
  // vector-valued matrices $\mathbf{c}_{ij}$ and $\mathbf{n}_{ij} =
  // \frac{\mathbf{c}_{ij}}{|\mathbf{c}_{ij}|}$, and the boundary normals
  // $\boldsymbol{\nu}_i$.
  //
  // In order to exploit thread parallelization we use WorkStream approach
  // detailed in the @ref threads "Parallel computing with multiple processors
  // accessing shared memory". As customary this requires
  // definition of
  //  - Scratch data (i.e. input info required to carry out computations): in 
  //    this case it is <code>scratch_data</code>.
  //  - The worker: in the case it is <code>local_assemble_system</code> that
  //    actually computes the local (i.e. current cell) contributions from the 
  //    scratch data.
  //  - A copy data: a struct that contains all the local assembly
  //    contributions, in this case <code>CopyData<dim>()</code>.
  //  - A copy data routine: in this case it is
  //    <code>copy_local_to_global</code> in charge of actually coping these
  //    local contributions into the global objects (matrices and/or vectors)
  //
  // Most the following lines are spent in the definition of the worker
  // <code>local_assemble_system</code> and the copy data routine
  // <code>copy_local_to_global</code>. There is not much to say about the
  // WorkStream framework since the vast majority of ideas are reasonably
  // well-documented in Step-9, Step-13 and Step-32 among others.
  //
  // Finally, assuming that $\mathbf{x}_i$ is a support point at the boundary,
  // the normals are defined as
  //
  // $\widehat{\boldsymbol{\nu}}_i :=
  // \frac{\boldsymbol{\nu}_i}{|\boldsymbol{\nu}_i|}$ where
  // $\boldsymbol{\nu}_i := \sum_{T \in \text{supp}(\phi_i)}
  // \sum_{F \subset \partial T \cap \partial \Omega}
  // \sum_{\mathbf{x}_{q,F}} \nu(\mathbf{x}_{q,F})
  // \phi_i(\mathbf{x}_{q,F})$
  //
  // here $T$ denotes elements,
  // $\text{supp}(\phi_i)$ the support of the shape function $\phi_i$,
  // $F$ are faces of the element $T$, and $\mathbf{x}_{q,F}$
  // are quadrature points on such face.
  // Other more sophisticated definitions for $\nu_i$ are
  // possible but none of them have much influence in theory or practice.

  template <int dim>
  void OfflineData<dim>::assemble()
  {
    lumped_mass_matrix = 0.;
    norm_matrix        = 0.;
    for (auto &matrix : cij_matrix)
      matrix = 0.;
    for (auto &matrix : nij_matrix)
      matrix = 0.;

    const unsigned int dofs_per_cell =
      discretization->finite_element.dofs_per_cell;
    const unsigned int n_q_points = discretization->quadrature.size();

    /* This is the implementation of the scratch data required by WorkStream */
    MeshWorker::ScratchData<dim> scratch_data(
      discretization->mapping,
      discretization->finite_element,
      discretization->quadrature,
      update_values | update_gradients | update_quadrature_points |
        update_JxW_values,
      discretization->face_quadrature,
      update_normal_vectors | update_values | update_JxW_values);

    {
      TimerOutput::Scope t(
        computing_timer,
        "offline_data - assemble lumped mass matrix, and c_ij");

      /* This is the implementation of the "worker" required by WorkStream */
      const auto local_assemble_system = [&](const auto &cell,
                                             auto &      scratch,
                                             auto &      copy) {
        auto &is_artificial             = copy.is_artificial;
        auto &local_dof_indices         = copy.local_dof_indices;
        auto &local_boundary_normal_map = copy.local_boundary_normal_map;
        auto &cell_lumped_mass_matrix   = copy.cell_lumped_mass_matrix;
        auto &cell_cij_matrix           = copy.cell_cij_matrix;

        is_artificial = cell->is_artificial();
        if (is_artificial)
          return;

        local_boundary_normal_map.clear();
        cell_lumped_mass_matrix.reinit(dofs_per_cell, dofs_per_cell);
        for (auto &matrix : cell_cij_matrix)
          matrix.reinit(dofs_per_cell, dofs_per_cell);

        const auto &fe_values = scratch.reinit(cell);

        local_dof_indices.resize(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);

        std::transform(local_dof_indices.begin(),
                       local_dof_indices.end(),
                       local_dof_indices.begin(),
                       [&](auto index) {
                         return partitioner->global_to_local(index);
                       });

        /* We compute the local contributions for the lumped mass
         matrix entries m_i and and vectors c_ij */
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {
            const auto JxW = fe_values.JxW(q_point);

            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                const auto value_JxW = fe_values.shape_value(j, q_point) * JxW;
                const auto grad_JxW  = fe_values.shape_grad(j, q_point) * JxW;

                cell_lumped_mass_matrix(j, j) += value_JxW;

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  {
                    const auto value = fe_values.shape_value(i, q_point);
                    for (unsigned int d = 0; d < dim; ++d)
                      cell_cij_matrix[d](i, j) += (value * grad_JxW)[d];

                  } /* for i */
              }     /* for j */
          }         /* for q */

        /* Now we have to compute the boundary normals. Note that the
           following loop does not actually do much unless the the element
           has faces on the boundary of the domain */
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
          {
            const auto face = cell->face(f);
            const auto id   = face->boundary_id();

            if (!face->at_boundary())
              continue;

            const auto &fe_face_values = scratch.reinit(cell, f);

            const unsigned int n_face_q_points =
              fe_face_values.get_quadrature().size();

            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                if (!discretization->finite_element.has_support_on_face(j, f))
                  continue;

                /* Note that "normal" will only represent the contributions
                   from one of the faces in the support of the shape
                   function phi_j. So we cannot normalize this local
                   contribution right here, we have to take it "as is", store
                   it and pass it to the copy data routine. The proper
                   normalization requires an additional loop on nodes.*/
                Tensor<1, dim> normal;
                if (id == Boundary::slip)
                  {
                    for (unsigned int q = 0; q < n_face_q_points; ++q)
                      normal += fe_face_values.normal_vector(q) *
                                fe_face_values.shape_value(j, q);
                  }

                const auto index = local_dof_indices[j];

                Point<dim> position;
                const auto global_index = partitioner->local_to_global(index);
                for (unsigned int v = 0;
                     v < GeometryInfo<dim>::vertices_per_cell;
                     ++v)
                  if (cell->vertex_dof_index(v, 0) == global_index)
                    position = cell->vertex(v);

                const auto old_id =
                  std::get<1>(local_boundary_normal_map[index]);
                local_boundary_normal_map[index] =
                  std::make_tuple(normal, std::max(old_id, id), position);
              } /* done with the loop on shape functions */
          }     /* done with the loop on faces */
      };        /* done with the definition of the worker */

      /* This is the copy data routine for WorkStream */
      const auto copy_local_to_global = [&](const auto &copy) {
        const auto &is_artificial             = copy.is_artificial;
        const auto &local_dof_indices         = copy.local_dof_indices;
        const auto &local_boundary_normal_map = copy.local_boundary_normal_map;
        const auto &cell_lumped_mass_matrix   = copy.cell_lumped_mass_matrix;
        const auto &cell_cij_matrix           = copy.cell_cij_matrix;

        if (is_artificial)
          return;

        for (const auto &it : local_boundary_normal_map)
          {
            auto &normal   = std::get<0>(boundary_normal_map[it.first]);
            auto &id       = std::get<1>(boundary_normal_map[it.first]);
            auto &position = std::get<2>(boundary_normal_map[it.first]);

            const auto &new_normal   = std::get<0>(it.second);
            const auto &new_id       = std::get<1>(it.second);
            const auto &new_position = std::get<2>(it.second);

            normal += new_normal;
            id       = std::max(id, new_id);
            position = new_position;
          }

        lumped_mass_matrix.add(local_dof_indices, cell_lumped_mass_matrix);

        for (int k = 0; k < dim; ++k)
          {
            cij_matrix[k].add(local_dof_indices, cell_cij_matrix[k]);
            nij_matrix[k].add(local_dof_indices, cell_cij_matrix[k]);
          }
      }; /* end of the copy data routine */

      WorkStream::run(dof_handler.begin_active(),
                      dof_handler.end(),
                      local_assemble_system,
                      copy_local_to_global,
                      scratch_data,
                      CopyData<dim>());
    } /* We are done with m_i and c_{ij} */

    // At this point in time we are done with the computation of $m_i$ and
    // $\mathbf{c}_{ij}$, but so far the matrix <code>nij_matrix</code>
    // contains a just copy of the matrix <code>cij_matrix</code>.
    // That's not what we really
    // want: we have to normalize its entries. In addition, we have not filled
    // the entries of the matrix <code>norm_matrix</code>  and the
    // vectors stored in the map
    // <code>OfflineData<dim>::BoundaryNormalMap</code> are not normalized.
    //
    // In principle, this is just offline data, it doesn't make much sense
    // to over-optimize their computation, since their cost will get amortized
    // over the many time steps that we are going to use. However,
    // computing/storing the entries of the matrix
    // <code>norm_matrix</code> and the normalization of <code>nij_matrix</code>
    // are perfect to illustrate thread-parallel node-loops:
    // - we want to visit every node $i$ in the mesh/sparsity graph,
    // - and for every such node we want to visit to every $j$ such that
    // $\mathbf{c}_{ij} \not \equiv 0$.
    //
    // From an algebraic point of view, this is equivalent to: visiting
    // every row in the matrix and for each one of these rows execute a loop on 
    // the columns. Node-loops is a core theme of this tutorial step (see 
    // the pseudo-code in the introduction) that will repeat over and over 
    // again. That's why this is the right time to introduce them.
    //
    // We have the thread paralellization capability
    // parallel::apply_to_subranges that is somehow more general than the
    // WorkStream framework. In particular, parallel::apply_to_subranges can
    // be used for our node-loops. This functionality requires four input 
    // arguments which we explain in detail (for the specific case of our 
    // thread-parallel node loops):
    // - The iterator <code>indices.begin()</code> points to
    //   to a row index.
    // - The iterator <code>indices.end()</code> points to a numerically higher 
    //   row index.
    // - The function <code>on_subranges(i1,i2)</code> (where <code>i1</code>
    //   and <code>i2</code> define sub-range within the range spanned by
    //   the end and begin iterators defined in the two previous bullets)
    //   applies operation for every iterator in such subrange. We may as well 
    //   call <code>on_subranges</code> the worker.
    // - Grainsize: minimum number of iterators (in this case representing 
    //   rows) processed by each thread. We decided for a minimum of 4096 
    //   rows.
    //
    // A minor caveat here is that the iterators <code>indices.begin()</code> 
    // and <code>indices.end()</code> supplied to
    // parallel::apply_to_subranges have to be random access iterators:
    // internally, apply_to_subranges will break the range defined by the
    // <code>indices.begin()</code> and <code>indices.end()</code> iterators
    // into subranges (we want to be able to read any entry in those
    // subranges with constant complexity). In order to provide such
    // iterators we resort to boost::irange.
    //
    // The bulk of the following piece of code is spent defining
    // the "worker" <code>on_subranges</code>: i.e. the  operation applied at 
    // each row of the sub-range. Given a fixed <code>row_index</code> 
    // we want to visit every column/entry in such row. In order to execute 
    // such columns-loops we use
    // <a href="http://www.cplusplus.com/reference/algorithm/for_each/">
    // std::for_each</a>
    // from the standard library, where:
    //  - <code>sparsity_pattern.begin(row_index)</code>
    //    gives us an iterator starting at the first column of the row,
    //  - <code>sparsity_pattern.end(row_index)</code> is an iterator pointing
    //    at the last column of the row,
    //  - the last argument required by std::for_each is the operation
    //    applied at each column (a lambda expression in this case) of
    //    such row.
    //
    // We note that, parallel::apply_to_subranges will operate on disjoint sets
    // of rows (the subranges) and our goal is to write into these rows.
    // Because of the simple nature of the operations we want to carry out
    // (computation and storage of normals, and normalization of the
    // $\mathbf{c}_{ij}$ of entries) threads cannot conflict attempting to
    // write the same entry (we do not need a scheduler).

    {
      TimerOutput::Scope t(computing_timer,
                           "offline_data - compute |c_ij|, and n_ij");

      /* Here [i1,i2] represent a subrange of rows */
      const auto on_subranges = [&](auto i1, const auto i2) {
        for (; i1 < i2; ++i1)
          {
            const auto row_index = *i1;

            /* First column-loop: we compute and store the entries of the matrix
            norm_matrix */
            std::for_each(sparsity_pattern.begin(row_index),
                          sparsity_pattern.end(row_index),
                          [&](const auto &jt) {
                            const auto value =
                              gather_get_entry(cij_matrix, &jt);
                            const double norm = value.norm();
                            set_entry(norm_matrix, &jt, norm);
                          });

            /* Second column-loop: we normalize the entries of the matrix
            nij_matrix */
            for (auto &matrix : nij_matrix)
              {
                auto nij_entry = matrix.begin(row_index);
                std::for_each(norm_matrix.begin(row_index),
                              norm_matrix.end(row_index),
                              [&](const auto &it) {
                                const auto norm = it.value();
                                nij_entry->value() /= norm;
                                ++nij_entry;
                              });
              }

          } /* row_index */
      };    /* done with the definition of "on_subranges" */

      const auto indices = boost::irange<unsigned int>(0, n_locally_relevant);
      parallel::apply_to_subranges(indices.begin(),
                                   indices.end(),
                                   on_subranges,
                                   4096);

    // Finally, we normalize the vector stored in
    // <code>OfflineData<dim>::BoundaryNormalMap</code>. This operation has
    // not been thread paralellized as it would neither illustrate any important
    // concept nor lead to any noticeable speed gain.

      for (auto &it : boundary_normal_map)
        {
          auto &normal = std::get<0>(it.second);
          normal /= (normal.norm() + std::numeric_limits<double>::epsilon());
        }
    }

    // In order to implement reflecting boundary conditions
    // $\mathbf{m} \cdot \boldsymbol{\nu}_i =0$ (or equivalently $\mathbf{v}
    // \cdot \boldsymbol{\nu}_i =0$ ) the vectors $\mathbf{c}_{ij}$ at the
    // boundary have to be modified as:
    //
    // $\mathbf{c}_{ij} \, +\!\!= \int_{\partial \Omega}
    // (\boldsymbol{\nu}_j - \boldsymbol{\nu}(s)) \phi_j \, \mathrm{d}s$
    //
    // Otherwise we will not be able to claim conservation. The ideas repeat
    // themselves: we use Workstream in order to compute this correction, most
    // of the following code is about the definition of the worker
    // <code>local_assemble_system</code>.

    {
      TimerOutput::Scope t(computing_timer,
                           "offline_data - fix slip boundary c_ij");

      const auto local_assemble_system = [&](const auto &cell,
                                             auto &      scratch,
                                             auto &      copy) {
        auto &is_artificial     = copy.is_artificial;
        auto &local_dof_indices = copy.local_dof_indices;

        auto &cell_cij_matrix = copy.cell_cij_matrix;

        is_artificial = cell->is_artificial();
        if (is_artificial)
          return;

        for (auto &matrix : cell_cij_matrix)
          matrix.reinit(dofs_per_cell, dofs_per_cell);

        local_dof_indices.resize(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);
        std::transform(local_dof_indices.begin(),
                       local_dof_indices.end(),
                       local_dof_indices.begin(),
                       [&](auto index) {
                         return partitioner->global_to_local(index);
                       });

        for (auto &matrix : cell_cij_matrix)
          matrix = 0.;

        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
          {
            const auto face = cell->face(f);
            const auto id   = face->boundary_id();

            if (!face->at_boundary())
              continue;

            if (id != Boundary::slip)
              continue;

            const auto &fe_face_values = scratch.reinit(cell, f);

            const unsigned int n_face_q_points =
              fe_face_values.get_quadrature().size();

            for (unsigned int q = 0; q < n_face_q_points; ++q)
              {
                const auto JxW      = fe_face_values.JxW(q);
                const auto normal_q = fe_face_values.normal_vector(q);

                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  {
                    if (!discretization->finite_element.has_support_on_face(j,
                                                                            f))
                      continue;

                    const auto &[normal_j, _1, _2] =
                      boundary_normal_map[local_dof_indices[j]];

                    const auto value_JxW =
                      fe_face_values.shape_value(j, q) * JxW;

                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                      {
                        const auto value = fe_face_values.shape_value(i, q);

                        /* This is the correction of the boundary c_ij */
                        for (unsigned int d = 0; d < dim; ++d)
                          cell_cij_matrix[d](i, j) +=
                            (normal_j[d] - normal_q[d]) * (value * value_JxW);
                      } /* i */
                  }     /* j */
              }         /* q */
          }             /* f */
      };                /* Done with the definition of the worker */

      const auto copy_local_to_global = [&](const auto &copy) {
        const auto &is_artificial     = copy.is_artificial;
        const auto &local_dof_indices = copy.local_dof_indices;
        const auto &cell_cij_matrix   = copy.cell_cij_matrix;

        if (is_artificial)
          return;

        for (int k = 0; k < dim; ++k)
          cij_matrix[k].add(local_dof_indices, cell_cij_matrix[k]);
      };

      WorkStream::run(dof_handler.begin_active(),
                      dof_handler.end(),
                      local_assemble_system,
                      copy_local_to_global,
                      scratch_data,
                      CopyData<dim>());
    }
  } /* assemble() */

  // At this point we are very much done with anything related to offline data.

  // @sect4{The class <code>ProblemDescription</code> implementation.}

  // In this section we describe the implementation of the class members of
  // <code>ProblemDescription</code>. All these class member only have meaning
  // in the context of Euler's equations using with ideal gas law. If we wanted
  // to re-purpose Step-69 for a different conservation law (say for instance
  // shallow water equations) the implementation of this entire class would
  // have to change (or wiped out in its entirety). But most of the other 
  // classes, in particular those defining loop structures, would remain 
  // unchanged.
  //
  // Now we define the implementation of the utility
  // functions <code>momentum</code>,
  // <code>internal_energy</code>, <code>pressure</code>,
  // <code>speed_of_sound</code>, and <code>f</code> (the flux of the system).
  // The functionality of each one of these functions is self-explanatory from
  // their names.

  template <int dim>
  DEAL_II_ALWAYS_INLINE inline Tensor<1, dim>
  ProblemDescription<dim>::momentum(const rank1_type &U)
  {
    Tensor<1, dim> result;
    std::copy(&U[1], &U[1 + dim], &result[0]);
    return result;
  }

  template <int dim>
  DEAL_II_ALWAYS_INLINE inline double
  ProblemDescription<dim>::internal_energy(const rank1_type &U)
  {
    const double &rho = U[0];
    const auto    m   = momentum(U);
    const double &E   = U[dim + 1];
    return E - 0.5 * m.norm_square() / rho;
  }

  template <int dim>
  DEAL_II_ALWAYS_INLINE inline double
  ProblemDescription<dim>::pressure(const rank1_type &U)
  {
    return (gamma - 1.) * internal_energy(U);
  }

  template <int dim>
  DEAL_II_ALWAYS_INLINE inline double
  ProblemDescription<dim>::speed_of_sound(const rank1_type &U)
  {
    const double &rho = U[0];
    const double  p   = pressure(U);

    return std::sqrt(gamma * p / rho);
  }

  template <int dim>
  DEAL_II_ALWAYS_INLINE inline typename ProblemDescription<dim>::rank2_type
  ProblemDescription<dim>::f(const rank1_type &U)
  {
    const double &rho = U[0];
    const auto    m   = momentum(U);
    const auto    p   = pressure(U);
    const double &E   = U[dim + 1];

    rank2_type result;

    result[0] = m;
    for (unsigned int i = 0; i < dim; ++i)
      {
        result[1 + i] = m * m[i] / rho;
        result[1 + i][i] += p;
      }
    result[dim + 1] = m / rho * (E + p);

    return result;
  }

  // Now we discuss the computation of $\lambda_{\text{max}}
  // (\mathbf{U}_i^{n},\mathbf{U}_j^{n}, \textbf{n}_{ij})$. The analysis
  // and derivation of sharp upper-bounds of maximum wavespeeds of Riemann
  // problems is a very technical endeavor and we cannot include an
  // advanced discussion about it in this tutorial. In this portion
  // of the documentation we will limit ourselves to sketch the main
  // functionality of these auxiliary functions and point to specific
  // academic references in order to help (the interested) reader trace the
  // source (and proper mathematical justification) of these ideas.
  //
  // In general, obtaining a sharp guaranteed upper-bound on the maximum
  // wavespeed requires solving a quite expensive scalar nonlinear problem.
  // In order to simplify the presentation we decided not to include such
  // iterative scheme. Here we have taken the following shortcut: formulas
  // (2.11) (3.7), (3.8) and (4.3) from
  //
  // - J-L Guermond, B. Popov, Fast estimation of the maximum wave speed in
  //  the Riemann problem for the Euler equations, JCP, 2016,
  //
  // are enough to define a guaranteed upper bound on the maximum
  // wavespeed. This estimate is returned by the a call to the function
  // <code>lambda_max_two_rarefaction</code>. At its core the construction
  // of such upper bound uses the so-called two-rarefaction approximation
  // for the intermediate pressure $p^*$, see for instance
  //
  // - Formula (4.46), page 128 in: E.Toro, Riemann Solvers and Numerical
  //   Methods for Fluid Dynamics, 2009.
  //
  // The estimate <code>lambda_max_two_rarefaction</code>
  // is in general very sharp and it would be enough for the
  // purposes of this code. However, for some specific situations (in
  // particular when one of states is close to vacuum conditions) such
  // estimate will be overly pessimistic. That's why we used a second
  // estimate to avoid this degeneracy that will be invoked by a call to
  // the function <code>lambda_max_expansion</code>. The most important 
  // function here is <code>compute_lambda_max</code> which takes the minimum 
  // between the estimates
  // - <code>lambda_max_two_rarefaction</code>
  // - <code>lambda_max_expansion</code>
  //
  // The remaining functions
  // - <code>riemann_data_from_state</code>
  // - <code>positive_part</code>
  // - <code>negative_part</code>
  // - <code>lambda1_minus</code>
  // - <code>lambda2_minus</code>
  //
  // are just auxiliary functions required in order to compute both estimates.

  namespace
  {
    template <int dim>
    DEAL_II_ALWAYS_INLINE inline std::array<double, 4> riemann_data_from_state(
      const typename ProblemDescription<dim>::rank1_type U,
      const Tensor<1, dim> &                             n_ij)
    {
      Tensor<1, 3> projected_U;
      projected_U[0] = U[0];

      const auto m   = ProblemDescription<dim>::momentum(U);
      projected_U[1] = n_ij * m;

      const auto perpendicular_m = m - projected_U[1] * n_ij;
      projected_U[2] = U[1 + dim] - 0.5 * perpendicular_m.norm_square() / U[0];

      std::array<double, 4> result;
      result[0] = projected_U[0];
      result[1] = projected_U[1] / projected_U[0];
      result[2] = ProblemDescription<1>::pressure(projected_U);
      result[3] = ProblemDescription<1>::speed_of_sound(projected_U);

      return result;
    }


    DEAL_II_ALWAYS_INLINE inline double positive_part(const double number)
    {
      return (std::abs(number) + number) / 2.0;
    }


    DEAL_II_ALWAYS_INLINE inline double negative_part(const double number)
    {
      return (std::fabs(number) - number) / 2.0;
    }


    /* Implements formula (3.7) in Guermond-Popov-2016 */
    DEAL_II_ALWAYS_INLINE inline double
    lambda1_minus(const std::array<double, 4> &riemann_data,
                  const double                 p_star)
    {
      constexpr double gamma             = ProblemDescription<1>::gamma;
      const auto &[rho_Z, u_Z, p_Z, a_Z] = riemann_data;

      const double factor = (gamma + 1.0) / 2.0 / gamma;
      const double tmp    = positive_part((p_star - p_Z) / p_Z);
      return u_Z - a_Z * std::sqrt(1.0 + factor * tmp);
    }


    /* Implements formula (3.8) in Guermond-Popov-2016 */
    DEAL_II_ALWAYS_INLINE inline double
    lambda3_plus(const std::array<double, 4> &riemann_data, const double p_star)
    {
      constexpr double gamma             = ProblemDescription<1>::gamma;
      const auto &[rho_Z, u_Z, p_Z, a_Z] = riemann_data;

      const double factor = (gamma + 1.0) / 2.0 / gamma;
      const double tmp    = positive_part((p_star - p_Z) / p_Z);
      return u_Z + a_Z * std::sqrt(1.0 + factor * tmp);
    }


    /* Implements formula (2.11) in Guermond-Popov-2016*/
    DEAL_II_ALWAYS_INLINE inline double
    lambda_max_two_rarefaction(const std::array<double, 4> &riemann_data_i,
                               const std::array<double, 4> &riemann_data_j)
    {
      constexpr double gamma             = ProblemDescription<1>::gamma;
      const auto &[rho_i, u_i, p_i, a_i] = riemann_data_i;
      const auto &[rho_j, u_j, p_j, a_j] = riemann_data_j;

      const double numerator = a_i + a_j - (gamma - 1.) / 2. * (u_j - u_i);

      const double denominator =
        a_i * std::pow(p_i / p_j, -1. * (gamma - 1.) / 2. / gamma) + a_j * 1.;

      /* Formula (4.3) in Guermond-Popov-2016 */
      const double p_star =
        p_j * std::pow(numerator / denominator, 2. * gamma / (gamma - 1));

      const double lambda1 = lambda1_minus(riemann_data_i, p_star);
      const double lambda3 = lambda3_plus(riemann_data_j, p_star);

      /* Returns formula (2.11) in Guermond-Popov-2016 */
      return std::max(positive_part(lambda3), negative_part(lambda1));
    };


    /* This estimate is, in general, not as sharp as the two-rarefaction
       estimate. But it will save the day in the context of near vacuum
       conditions when the two-rarefaction approximation will tend to
       exaggerate the maximum wave speed. */
    DEAL_II_ALWAYS_INLINE inline double
    lambda_max_expansion(const std::array<double, 4> &riemann_data_i,
                         const std::array<double, 4> &riemann_data_j)
    {
      const auto &[rho_i, u_i, p_i, a_i] = riemann_data_i;
      const auto &[rho_j, u_j, p_j, a_j] = riemann_data_j;

      /* Here the constant 5.0 multiplying the soundspeeds is NOT
         an ad-hoc constant or tuning parameter. It defines a upper bound
         for any $\gamma \in [0,5/3]$. Do not play with it! */
      return std::max(std::abs(u_i), std::abs(u_j)) + 5. * std::max(a_i, a_j);
    }
  } // namespace

  // The is the main function that we are going to call in order to compute
  // $\lambda_{\text{max}}
  // (\mathbf{U}_i^{n},\mathbf{U}_j^{n}, \textbf{n}_{ij})$.
  template <int dim>
  DEAL_II_ALWAYS_INLINE inline double
  ProblemDescription<dim>::compute_lambda_max(const rank1_type &    U_i,
                                              const rank1_type &    U_j,
                                              const Tensor<1, dim> &n_ij)
  {
    const auto riemann_data_i = riemann_data_from_state(U_i, n_ij);
    const auto riemann_data_j = riemann_data_from_state(U_j, n_ij);

    const double lambda_1 =
      lambda_max_two_rarefaction(riemann_data_i, riemann_data_j);

    const double lambda_2 =
      lambda_max_expansion(riemann_data_i, riemann_data_j);

    return std::min(lambda_1, lambda_2);
  }

  // Here <code>component_names</code> are just tags
  // that we will use for the output. We consider the template specializations
  // for dimensions dimensions one, two and three.

  template <>
  const std::array<std::string, 3> ProblemDescription<1>::component_names{"rho",
                                                                          "m",
                                                                          "E"};

  template <>
  const std::array<std::string, 4> ProblemDescription<2>::component_names{"rho",
                                                                          "m_1",
                                                                          "m_2",
                                                                          "E"};

  template <>
  const std::array<std::string, 5> ProblemDescription<3>::component_names{"rho",
                                                                          "m_1",
                                                                          "m_2",
                                                                          "m_3",
                                                                          "E"};

  // @sect4{Class <code>InitialValues</code> implementation}

  // Constructor for the class InitialValues. We add some parameters with 
  // some default values. We also provide a non-empty an implementation 
  // for the class member <code>parse_parameters_call_back</code>.
  //
  // The class member <code>parse_parameters_call_back</code> (inherited 
  // ParameterAcceptor) has an empty implementation by default.
  // This function will only be invoked for every class that is derived 
  // from ParameterAceptor after the call to ParameterAcceptor::initialize. In 
  // that regard, its use is appropriate for situations where the parameters 
  // have to be postprocessed (in some sense) or some consistency 
  // condition between the parameters has to be checked.

  template <int dim>
  InitialValues<dim>::InitialValues(const std::string &subsection)
    : ParameterAcceptor(subsection)
  {
    /* We wire-up InitialValues<dim>::parse_parameters_callback (declared 
       a few lines below) to ParameterAcceptor::parse_parameters_call_back */
    ParameterAcceptor::parse_parameters_call_back.connect(
      std::bind(&InitialValues<dim>::parse_parameters_callback, this));

    initial_direction[0] = 1.;
    add_parameter("initial direction",
                  initial_direction,
                  "Initial direction of the uniform flow field");

    static constexpr auto gamma = ProblemDescription<dim>::gamma;
    initial_1d_state[0]         = gamma;
    initial_1d_state[1]         = 3.;
    initial_1d_state[2]         = 1.;
    add_parameter("initial 1d state",
                  initial_1d_state,
                  "Initial 1d state (rho, u, p) of the uniform flow field");
  }

  // So far the constructor of <code>InitialValues</code> has defined 
  // default values for the two private members <code>initial_direction</code> 
  // and <code>initial_1d_state</code> and added them to the parameter list. 
  // But we have not defined an implementation for the only public member that 
  // we really care about, which is <code>initial_state</code> (the 
  // function that we are going to call to actually evaluate the initial 
  // solution at the mesh nodes).
  //
  // As commented, we could have avoided using the method 
  // <code>parse_parameters_call_back </code> and define a class member 
  // <code>setup()</code> in order to define the implementation of 
  // <code>initial_state</code>. But this illustrates a different way to use 
  // inheritance of ParameterAceptor to our benefit.

  template <int dim>
  void InitialValues<dim>::parse_parameters_callback()
  {
    AssertThrow(initial_direction.norm() != 0.,
                ExcMessage(
                  "Initial shock front direction is set to the zero vector."));
    initial_direction /= initial_direction.norm();

    static constexpr auto gamma = ProblemDescription<dim>::gamma;

    /* Function that translates primitive 1d states in to conserved 2d states.
       Note that we have some room for freedom to change the direction of the 
       flow. */
    const auto from_1d_state =
      [=](const Tensor<1, 3, double> &state_1d) -> rank1_type {
      const auto &rho = state_1d[0];
      const auto &u   = state_1d[1];
      const auto &p   = state_1d[2];

      rank1_type state;

      state[0] = rho;
      for (unsigned int i = 0; i < dim; ++i)
        state[1 + i] = rho * u * initial_direction[i];
      state[dim + 1] = p / (gamma - 1.) + 0.5 * rho * u * u;

      return state;
    };

    initial_state = [=](const Point<dim> & /*point*/, double /*t*/) {
      return from_1d_state(initial_1d_state);
    };
  }

  // @sect4{Class <code>TimeStep</code> implementation}

  template <int dim>
  TimeStep<dim>::TimeStep(const MPI_Comm &          mpi_communicator,
                          TimerOutput &             computing_timer,
                          const OfflineData<dim> &  offline_data,
                          const InitialValues<dim> &initial_values,
                          const std::string &       subsection /*= "TimeStep"*/)
    : ParameterAcceptor(subsection)
    , mpi_communicator(mpi_communicator)
    , computing_timer(computing_timer)
    , offline_data(&offline_data)
    , initial_values(&initial_values)
  {
    cfl_update = 0.80;
    add_parameter("cfl update",
                  cfl_update,
                  "relative CFL constant used for update");
  }

  // In the class member <code>prepare()</code> we set the partition of the 
  // auxiliary vector <code>temp</code> (locally owned + ghosted layer) and 
  // set the sparsity pattern for <code>dij_matrix</code> (borrowed from 
  // offline_data, a pointer to the unique OfflineData instance).
  // The vector <code>temp</code> will be used to store temporarily the 
  // solution update, to later swap its contents with the old vector.

  template <int dim>
  void TimeStep<dim>::prepare()
  {
    TimerOutput::Scope time(computing_timer,
                            "time_step - prepare scratch space");

    const auto &partitioner = offline_data->partitioner;
    for (auto &it : temp)
      it.reinit(partitioner);

    const auto &sparsity = offline_data->sparsity_pattern;
    dij_matrix.reinit(sparsity);
  }

  // An efficient implementation of the class member
  // <code>TimeStep<dim>::step</code>
  // should only compute the quantities that evolve for
  // every time-step (the fluxes  $\mathbb{f}(\mathbf{U}_j^{n})$ and
  // the viscosities $d_{ij}$) and assemble the new solution 
  // $\mathbf{U}_i^{n+1}$:
  // - We execute thread-parallel node-loops using
  //   <code>parallel::apply_to_subranges</code> for all the necessary tasks. 
  //   Pretty much all the ideas used to compute/store the entries of the
  //   matrix <code>norm_matrix</code> and the normalization of 
  //   <code>nij_matrix</code> (described a few hundreds of lines above) 
  //   are used here again. Most of the code intricacies lie around the 
  //   definition of the new "workers" <code>on_subranges</code> required for 
  //   the new tasks.
  // - The first step is computing the matrix the viscosities of $d_{ij}$. 
  //   It is important to highlight that viscosities are bound to the 
  //   constraint $d_{ij} = d_{ji}$ and our algorithm should reflect that. 
  //   In this regard we note here that
  //   $\int_{\Omega} \nabla \phi_j \phi_i \, \mathrm{d}\mathbf{x}= -
  //   \int_{\Omega} \nabla \phi_i \phi_j \, \mathrm{d}\mathbf{x}$
  //   (or equivanlently $\mathbf{c}_{ij} = - \mathbf{c}_{ji}$) provided 
  //   either $\mathbf{x}_i$ or $\mathbf{x}_j$ is a support point at the 
  //   boundary. In such case we can check that 
  //   $\lambda_{\text{max}} (\mathbf{U}_i^{n}, \mathbf{U}_j^{n},
  //   \textbf{n}_{ij}) = \lambda_{\text{max}} (\mathbf{U}_j^{n},
  //   \mathbf{U}_i^{n},\textbf{n}_{ji})$
  //   by construction, which guarantees the property $d_{ij} = d_{ji}$.
  //   However, if both support points $\mathbf{x}_i$ or $\mathbf{x}_j$ happen 
  //   to lie on the boundary then the equalities $\mathbf{c}_{ij} = - 
  //   \mathbf{c}_{ji}$ and $\lambda_{\text{max}}
  //   (\mathbf{U}_i^{n}, \mathbf{U}_j^{n},
  //   \textbf{n}_{ij}) = \lambda_{\text{max}} (\mathbf{U}_j^{n},
  //   \mathbf{U}_i^{n},
  //   \textbf{n}_{ji})$ are not necessarily true. The only mathematically
  //   safe solution for this dilemma is to compute both of them and take the
  //   largest one.
  //
  // In order to increase the efficiency we only compute the
  // upper-triangular entries of $d_{ij}$ and copy the corresponding
  // entries to the lower-triangular part. Note that this strategy
  // intrinsically makes the assumption that memory access to the lower
  // triangular entries is inexpensive (they are cached, or somehow local 
  // memorywise).
  //
  // *** IT: Clarify, why is this the case? I don't think CRS has anything to 
  // do with it. Is the Cuthill_McKee inducing/creating data locality 
  // here? ***
  //

  template <int dim>
  double TimeStep<dim>::step(vector_type &U, double t)
  {
    const auto &n_locally_owned    = offline_data->n_locally_owned;
    const auto &n_locally_relevant = offline_data->n_locally_relevant;

    const auto indices_owned = boost::irange<unsigned int>(0, n_locally_owned);
    const auto indices_relevant =
      boost::irange<unsigned int>(0, n_locally_relevant);

    const auto &sparsity = offline_data->sparsity_pattern;

    const auto &lumped_mass_matrix = offline_data->lumped_mass_matrix;
    const auto &norm_matrix        = offline_data->norm_matrix;
    const auto &nij_matrix         = offline_data->nij_matrix;
    const auto &cij_matrix         = offline_data->cij_matrix;

    const auto &boundary_normal_map = offline_data->boundary_normal_map;

    {
      TimerOutput::Scope time(computing_timer, "time_step - 1 compute d_ij");

      /* Definition of the "worker" that computes the viscosity d_{ij} */
      const auto on_subranges = [&](auto i1, const auto i2) {
        for (const auto i : boost::make_iterator_range(i1, i2))
          {
            const auto U_i = gather(U, i);

            /* Column-loop */
            for (auto jt = sparsity.begin(i); jt != sparsity.end(i); ++jt)
              {
                const auto j = jt->column();

                /* We compute only dij if i < j (upper triangular entries) and 
                   later we copy this entry into dji. */
                if (j >= i)
                  continue;

                const auto U_j = gather(U, j);

                const auto   n_ij = gather_get_entry(nij_matrix, jt);
                const double norm = get_entry(norm_matrix, jt);

                const auto lambda_max =
                  ProblemDescription<dim>::compute_lambda_max(U_i, U_j, n_ij);

                double d = norm * lambda_max;

                /* If both support points happen to be at the boundary
                   we have to compute dji too and then take max(dij,dji) */
                if (boundary_normal_map.count(i) != 0 &&
                    boundary_normal_map.count(j) != 0)
                  {
                    const auto n_ji = gather(nij_matrix, j, i);
                    const auto lambda_max_2 =
                      ProblemDescription<dim>::compute_lambda_max(U_j,
                                                                  U_i,
                                                                  n_ji);
                    const double norm_2 = norm_matrix(j, i);

                    d = std::max(d, norm_2 * lambda_max_2);
                  }

                /* We set the upper triangular entry */
                set_entry(dij_matrix, jt, d);
                /* We set the lower triangular entry */
                dij_matrix(j, i) = d;
              } /* End of column-loop */
          }     /* End of row-loop */
      };        /* End of definition of on_subranges */

      parallel::apply_to_subranges(indices_relevant.begin(),
                                   indices_relevant.end(),
                                   on_subranges,
                                   4096);
    } /* End of the computation of the off-diagonal entries of dij_matrix */

    // So far the matrix <code>dij_matrix</code> contains the off-diagonal
    // components. We still have to fill its diagonal entries defined as
    // $d_{ii}^n = - \sum_{j \in \mathcal{I}(i)\backslash \{i\}} d_{ij}^n$. We
    // use again <code>parallel::apply_to_subranges</code> for this purpose. 
    // While in the process of computing the $d_{ii}$'s we also record the 
    // largest admissible time-step, which is defined as
    //
    // \f[ \tau_n := c_{\text{cfl}}\,\min_{
    // i\in\mathcal{V}}\left(\frac{m_i}{-2\,d_{ii}^{n}}\right) \, . \f]
    //
    // Note that the operation $\min_{i \in \mathcal{V}}$ is intrinsically
    // global, it operates on all nodes: first we would have to first take the
    // $\min$ among all threads and finally take the $\min$ among all MPI
    // processes. In the current implementation:
    // - We do not take the $\min$ among threads: we simply define
    //  <code>tau_max</code> as <a
    //  href="http://www.cplusplus.com/reference/atomic/atomic/">
    //  std::atomic<double> </a>. The internal implementation of std::atomic
    //  will take care of resolving any possible conflict when more than
    //  one thread attempts read or write tau_max at the same time.
    // - In order to take the min among all MPI process we use the utility
    //   <code>Utilities::MPI::min</code>.

    /* We define tau_max as an atomic double in order to avoid any read/write 
       conflicts between threads and initialize it as the largest possible 
       number that can be represented by the float-type double. */
    std::atomic<double> tau_max{std::numeric_limits<double>::infinity()};

    {
      TimerOutput::Scope time(computing_timer,
                              "time_step - 2 compute d_ii, and tau_max");

      const auto on_subranges = [&](auto i1, const auto i2) {
        double tau_max_on_subrange = std::numeric_limits<double>::infinity();

        for (const auto i : boost::make_iterator_range(i1, i2))
          {
            double d_sum = 0.;

            /* See the definition of dii in the introduction. */
            for (auto jt = sparsity.begin(i); jt != sparsity.end(i); ++jt)
              {
                const auto j = jt->column();

                if (j == i)
                  continue;

                d_sum -= get_entry(dij_matrix, jt);
              }

            dij_matrix.diag_element(i) = d_sum;

            const double mass = lumped_mass_matrix.diag_element(i);
            /* See the definition of time-step constraint (CFL) */
            const double tau    = cfl_update * mass / (-2. * d_sum);
            tau_max_on_subrange = std::min(tau_max_on_subrange, tau);
          }

        double current_tau_max = tau_max.load();
        while (
          current_tau_max > tau_max_on_subrange &&
          !tau_max.compare_exchange_weak(current_tau_max, tau_max_on_subrange))
          ;
      }; /* End of definition of the worker on_subranges */

      /* Thread-parallel loop on locally owned rows */
      parallel::apply_to_subranges(indices_relevant.begin(),
                                   indices_relevant.end(),
                                   on_subranges,
                                   4096);

      /* We find the tau_max min among all MPI processes */
      tau_max.store(Utilities::MPI::min(tau_max.load(), mpi_communicator));

      AssertThrow(!std::isnan(tau_max) && !std::isinf(tau_max) && tau_max > 0.,
                  ExcMessage("I'm sorry, Dave. I'm afraid I can't "
                             "do that. - We crashed."));
    } /* End of the computation of the diagonal entries of dij_matrix */

    // At this point, we have computed all viscosity coefficients $d_{ij}$ and
    // we know what is the maximum time-step size we can use (which is,
    // strictly speaking, a consequence of the size of the viscosity
    // coefficients). So we compute the update as:
    //
    // \f[\mathbf{U}_i^{n+1} = \mathbf{U}_i^{n} - \frac{\tau_{\text{max}} }{m_i}
    //  \sum_{j \in \mathcal{I}(i)} (\mathbb{f}(\mathbf{U}_j^{n}) -
    //  \mathbb{f}(\mathbf{U}_i^{n})) \cdot \mathbf{c}_{ij} - d_{ij}
    //  (\mathbf{U}_j^{n} - \mathbf{U}_i^{n})\f]
    //
    // This update formula is different from that one used in the
    // pseudo-code. However, it can be shown that it is algebraically
    // equivalent (it will produce the same numerical values). We favor
    // this second formula since it has natural cancellation properties
    // that might help avoid numerical artifacts.

    {
      TimerOutput::Scope time(computing_timer, "time_step - 3 perform update");

      /* We define the "worker" for the subranges of rows */
      const auto on_subranges = [&](auto i1, const auto i2) {
        for (const auto i : boost::make_iterator_range(i1, i2))
          {
            Assert(i < n_locally_owned, ExcInternalError());

            const auto U_i = gather(U, i);

            const auto   f_i = ProblemDescription<dim>::f(U_i);
            const double m_i = lumped_mass_matrix.diag_element(i);

            auto U_i_new = U_i;

            /* This is the loop on the columns */
            for (auto jt = sparsity.begin(i); jt != sparsity.end(i); ++jt)
              {
                const auto j = jt->column();

                const auto U_j = gather(U, j);
                const auto f_j = ProblemDescription<dim>::f(U_j);

                const auto c_ij = gather_get_entry(cij_matrix, jt);
                const auto d_ij = get_entry(dij_matrix, jt);

                /* We define use the update formula here */
                for (unsigned int k = 0; k < problem_dimension; ++k)
                  {
                    U_i_new[k] +=
                      tau_max / m_i *
                      (-(f_j[k] - f_i[k]) * c_ij + d_ij * (U_j[k] - U_i[k]));
                  }
              }

            scatter(temp, U_i_new, i);
          }
      };

      /* Thread-parallel loop on locally owned rows */
      parallel::apply_to_subranges(indices_owned.begin(),
                                   indices_owned.end(),
                                   on_subranges,
                                   4096);
    } /* End of the computation of the new solution */

    // The vast majority of the updated values is right, except those at the
    // boundary which have to be corrected. This is known as
    // explicit treatment of the boundary conditions:
    // - You advance in time satisfying no boundary condition at all,
    // - At the end of the time step you enforce them (you post process
    //   your solution).
    //
    // When solving parabolic and/or elliptic equations, we know that: in order
    // to enforce essential boundary conditions we should make them part
    // of the approximation space, while natural boundary conditions
    // should become part of the variational formulation. We also know
    // that explicit treatment of the boundary conditions (in the context of
    // parabolic PDE) almost surely leads to catastrophic consequences.
    // However, in the context of nonlinear hyperbolic equations there is enough
    // numerical evidence suggesting that explicit treatment of essential
    // boundary conditions is stable (at least in the eye-ball norm) and does
    // not introduce any loss in accuracy (convergence rates). In addition,
    // it is relatively straightforward to prove that (for the case of
    // reflecting boundary conditions) explicit treatment of boundary
    // conditions is not only conservative but also guarantees preservation of
    // the invariant set. We are not aware of any theoretical result showing
    // that it is possible to provide such invariant-set guarantees when
    // using either direct enforcement of boundary conditions into the
    // approximation space and/or weak enforcement using Nitsche penalty
    // method (e.g. widely used in dG schemes).
    //
    // Here the worker <code>on_subranges</code> executes the correction
    //
    // $\mathbf{m}_i := \mathbf{m}_i - (\boldsymbol{\nu}_i \cdot \mathbf{m}_i)
    //   \boldsymbol{\nu}_i$
    //
    // which removes the normal component of $\mathbf{m}$. We note that
    // conservation is not just a consequence of this correction but also a
    // consequence of modification of the $\mathbf{c}_{ij}$ coefficients at the
    // boundary (see the third thread-parallel loop on nodes in
    // <code>OfflineData<dim>::assemble()</code>).

    {
      TimerOutput::Scope time(computing_timer,
                              "time_step - 4 fix boundary states");

      const auto on_subranges = [&](const auto it1, const auto it2) {
        for (auto it = it1; it != it2; ++it)
          {
            const auto i = it->first;

            /* Only iterate over locally owned subset */
            if (i >= n_locally_owned)
              continue;

            const auto &normal   = std::get<0>(it->second);
            const auto &id       = std::get<1>(it->second);
            const auto &position = std::get<2>(it->second);

            /* Skip constrained degrees of freedom */
            if (++sparsity.begin(i) == sparsity.end(i))
              continue;

            auto U_i = gather(temp, i);

            /* On boundary 1 remove the normal component of the momentum: */

            if (id == Boundary::slip)
              {
                auto m = ProblemDescription<dim>::momentum(U_i);
                m -= 1. * (m * normal) * normal;
                for (unsigned int k = 0; k < dim; ++k)
                  U_i[k + 1] = m[k];
              }

            /* On boundary 2 enforce initial conditions: */
            else if (id == Boundary::dirichlet)
              {
                U_i = initial_values->initial_state(position, t + tau_max);
              }

            scatter(temp, U_i, i);
          }
      };

      on_subranges(boundary_normal_map.begin(), boundary_normal_map.end());
    }

    for (auto &it : temp)
      it.update_ghost_values();

    U.swap(temp);

    return tau_max;
  } /* End of TimeStep<dim>::step */

  // @sect4{Class <code>SchlierenPostprocessor</code> implementation}

  // Here
  // - schlieren_beta: is an ad-hoc positive amplification factor in order to
  //   enhance/exaggerate contrast in the visualization. Its actual value is a
  //   matter of taste.
  // - schlieren_index: is a integer indicates which component of the 
  //   state $[\rho, \mathbf{m},E]$ are we going to use in order generate 
  //   the visualization.

  template <int dim>
  SchlierenPostprocessor<dim>::SchlierenPostprocessor(
    const MPI_Comm &        mpi_communicator,
    TimerOutput &           computing_timer,
    const OfflineData<dim> &offline_data,
    const std::string &     subsection /*= "SchlierenPostprocessor"*/)
    : ParameterAcceptor(subsection)
    , mpi_communicator(mpi_communicator)
    , computing_timer(computing_timer)
    , offline_data(&offline_data)
  {
    schlieren_beta = 10.;
    add_parameter("schlieren beta",
                  schlieren_beta,
                  "Beta factor used in Schlieren-type postprocessor");

    schlieren_index = 0;
    add_parameter("schlieren index",
                  schlieren_index,
                  "Use the corresponding component of the state vector for the "
                  "schlieren plot");
  }

  // Here <code>prepare()</code> initializes the vector <code>r</code>
  // and <code>schlieren</code> with proper sizes.

  template <int dim>
  void SchlierenPostprocessor<dim>::prepare()
  {
    TimerOutput::Scope t(computing_timer,
                         "schlieren_postprocessor - prepare scratch space");

    const auto &n_locally_relevant = offline_data->n_locally_relevant;
    const auto &partitioner        = offline_data->partitioner;

    r.reinit(n_locally_relevant);
    schlieren.reinit(partitioner);
  }

  // We now discuss the implementation of the class member 
  // <code>SchlierenPostprocessor<dim>::compute_schlieren</code>, which 
  // basically takes a component of the state vector <code>U</code> and 
  // computes the Schlieren indicator for such component (the formula of the
  // Schlieren indicator can be found just before the declaration of the class
  // <code>SchlierenPostprocessor</code>). We start by noting 
  // that this formula requires the "nodal gradients" $\nabla r_j$. 
  // However, nodal values of gradients are not defined for $\mathcal{C}^0$ 
  // finite element functions. More generally, pointwise values of gradients
  // are not defined for $W^{1,p}(\Omega)$ functions (though weak 
  // derivatives are). The simplest technique we can use to recover gradients 
  // at nodes is weighted-averaging i.e.
  //
  // \f[ \nabla r_j := \frac{1}{\int_{S_i} \omega_i(\mathbf{x}) \, 
  // \mathrm{d}\mathbf{x}}
  //  \int_{S_i} r_h(\mathbf{x}) \omega_i(\mathbf{x}) \, \mathrm{d}\mathbf{x}
  // \ \ \ \ \ \mathbf{(*)} \f]
  //
  // where $S_i$ is the support of the shape function $\phi_i$, and 
  // $\omega_i(\mathbf{x})$ is the weight. The weight could be any 
  // positive function such as 
  // $\omega_i(\mathbf{x}) \equiv 1$ (that would allow us to recover the usual 
  // notion of mean value). But as usual, the goal is to reuse the off-line 
  // data as much as it could be possible. In sense this, the most natural  
  // choice of weight is $\omega_i = \phi_i$. Inserting this choice of 
  // weight and the expansion $r_h(\mathbf{x}) = \sum_{j \in \mathcal{V}} 
  // r_j \phi_j(\mathbf{x})$ into $\mathbf{(*)}$ we get :
  //
  // \f[ \nabla r_j := \frac{1}{m_i} \sum_{j \in \mathcal{I}(i)} r_j
  //  \mathbf{c}_{ij} \ \ \ \ \ \mathbf{(**)} \, . \f] 
  //
  // Using this last formula we can recover averaged nodal gradients without
  // resorting to any form of quadrature. This idea aligns quite well with
  // the whole spirit of edge-based schemes (or algebraic schemes) where
  // we want to operate on matrices and vectors as directly as
  // it could be possible avoiding by all means assembly of bilinear
  // forms, cell-loops, quadrature, or any other
  // intermediate construct/operation between the input arguments (the state
  // from the previous time-step) and the actual matrices and vectors
  // required to compute the update.
  //
  // The second thing to note is that we have to compute global minimum and
  // maximums $\max_j |\nabla r_j|$ and $\min_j |\nabla r_j|$. Following the
  // same ideas used to compute the time step size in the class member
  // <code>TimeStep<dim>::step</code> we define $\max_j |\nabla r_j|$ and
  // $\min_j |\nabla r_j|$ as atomic doubles in order to
  // resolve any conflicts between threads. As usual, we use
  // <code>Utilities::MPI::max</code> and <code>Utilities::MPI::min</code> to
  // find the global maximum/minimum among all MPI processes.
  //
  // Finally, it is not possible to compute the Schlieren indicator in a single
  // loop over all nodes. The entire operation requires two loops over nodes:
  //
  // - The first loop computes $|\nabla r_i|$ for all $i \in \mathcal{V}$ in
  //   the mesh, and the bounds $\max_j |\nabla r_j|$ and
  //   $\min_j |\nabla r_j|$.
  // - The second loop finally computes the Schlieren indicator using the 
  //   formula
  //
  // \f[ \text{schlieren}[i] = e^{\beta \frac{ |\nabla r_i|
  // - \min_j |\nabla r_j| }{\max_j |\nabla r_j| - \min_j |\nabla r_j| } }
  // \, . \f]
  //
  // This means that we will have to define two workers
  // <code>on_subranges</code> for each one of these stages.

  template <int dim>
  void SchlierenPostprocessor<dim>::compute_schlieren(const vector_type &U)
  {
    TimerOutput::Scope t(computing_timer,
                         "schlieren_postprocessor - compute schlieren plot");

    const auto &sparsity            = offline_data->sparsity_pattern;
    const auto &lumped_mass_matrix  = offline_data->lumped_mass_matrix;
    const auto &cij_matrix          = offline_data->cij_matrix;
    const auto &boundary_normal_map = offline_data->boundary_normal_map;

    const auto &n_locally_owned = offline_data->n_locally_owned;
    const auto  indices = boost::irange<unsigned int>(0, n_locally_owned);

    /* We define the r_i_max and r_i_min in the current MPI process as
       atomic doubles in order to resolve conflicts among threads. */
    std::atomic<double> r_i_max{0.};
    std::atomic<double> r_i_min{std::numeric_limits<double>::infinity()};

    /* Implementation of the first worker: computes the averaged gradient 
       at each node and the global max and mins of such gradients. */
    {
      const auto on_subranges = [&](auto i1, const auto i2) {
        double r_i_max_on_subrange = 0.;
        double r_i_min_on_subrange = std::numeric_limits<double>::infinity();

        for (; i1 < i2; ++i1)
          {
            const auto i = *i1;

            Assert(i < n_locally_owned, ExcInternalError());

            Tensor<1, dim> r_i;

            /* This is the loop on the columns */
            /* We compute the numerator of expression (**) */
            for (auto jt = sparsity.begin(i); jt != sparsity.end(i); ++jt)
              {
                const auto j = jt->column();

                if (i == j)
                  continue;
 
                /* Usual practice is that schlieren_index = 0 (density of the 
                  system). In this tutorial step schlieren_index is set by the 
                  constructor. */
                const auto U_js = U[schlieren_index].local_element(j);
                const auto c_ij = gather_get_entry(cij_matrix, jt);

                r_i += c_ij * U_js;
              }

            const auto bnm_it = boundary_normal_map.find(i);
            if (bnm_it != boundary_normal_map.end())
              {
                const auto &normal = std::get<0>(bnm_it->second);
                const auto &id     = std::get<1>(bnm_it->second);

                if (id == Boundary::slip)
                  r_i -= 1. * (r_i * normal) * normal;
                else
                  r_i = 0.;
              }

            /* Here we remind the reader that we are not interested in the 
               nodal gradients per se. We want their norms in order to 
               compute the Schlieren indicator. Finally, we have to 
               divide r[i] by m_i. */
            const double m_i = lumped_mass_matrix.diag_element(i);
            r[i]             = r_i.norm() / m_i;

            r_i_max_on_subrange = std::max(r_i_max_on_subrange, r[i]);
            r_i_min_on_subrange = std::min(r_i_min_on_subrange, r[i]);
          }

        /* We compare the current_r_i_max and current_r_i_min (in the current 
           subrange) with r_i_max and r_i_min (for the current MPI process) 
           and update them if necessary */
        double current_r_i_max = r_i_max.load();
        while (
          current_r_i_max < r_i_max_on_subrange &&
          !r_i_max.compare_exchange_weak(current_r_i_max, r_i_max_on_subrange))
          ;

        double current_r_i_min = r_i_min.load();
        while (
          current_r_i_min > r_i_min_on_subrange &&
          !r_i_min.compare_exchange_weak(current_r_i_min, r_i_min_on_subrange))
          ;
      };

      parallel::apply_to_subranges(indices.begin(),
                                   indices.end(),
                                   on_subranges,
                                   4096);
    }

    r_i_max.store(Utilities::MPI::max(r_i_max.load(), mpi_communicator));
    r_i_min.store(Utilities::MPI::min(r_i_min.load(), mpi_communicator));

    /* Implementation of the second worker: we have the vector r_i and the 
       scalars r_i_max and r_i_min at our disposal. Now we are in position of 
       actually computing the Schlieren indicator. */

    {
      const auto on_subranges = [&](auto i1, const auto i2) {
        for (; i1 < i2; ++i1)
          {
            const auto i = *i1;

            Assert(i < n_locally_owned, ExcInternalError());

            /* It's just the Schlieren formula */
            /* There is no loop on columns for this case, we don't need it */
            schlieren.local_element(i) =
              1. - std::exp(-schlieren_beta * (r[i] - r_i_min) /
                            (r_i_max - r_i_min));
          }
      };

      parallel::apply_to_subranges(indices.begin(),
                                   indices.end(),
                                   on_subranges,
                                   4096);
    }

    schlieren.update_ghost_values();
  }

  // @sect4{The Timeloop class implementation.}

  // Constructor of the class <code>Timeloop</code>. Note that this class wraps 
  // up pretty much all the other classes that we have discussed so far. 
  // More precisely the constructor has to initialize an instance of 
  //   - <code>Discretization<dim> </code>
  //   - <code>OfflineData<dim> </code>
  //   - <code>InitialValues<dim> </code>
  //   - <code>TimeStep<dim> </code>
  //   - <code>SchlierenPostprocessor<dim> </code>
  //
  // Most of the functionality of the class
  // <code>Timeloop</code> comes from the methods of those five classes. In
  // itself, the class <code>TimeLoop<dim></code> only requires the 
  // implementation of three new class members/methods:
  //  - <code>TimeLoop<dim>::run </code>.
  //  - <code>TimeLoop<dim>::interpolate_initial_values </code>
  //  - <code>TimeLoop<dim>::output </code>
  //
  // Note that in the construction we also add the boolean parameter 
  // "resume" which will be used to restart interrupted computations.

  template <int dim>
  TimeLoop<dim>::TimeLoop(const MPI_Comm &mpi_comm)
    : ParameterAcceptor("A - TimeLoop")
    , mpi_communicator(mpi_comm)
    , computing_timer(mpi_communicator,
                      timer_output,
                      TimerOutput::never,
                      TimerOutput::cpu_and_wall_times)
    , pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    , discretization(mpi_communicator, computing_timer, "B - Discretization")
    , offline_data(mpi_communicator,
                   computing_timer,
                   discretization,
                   "C - OfflineData")
    , initial_values("D - InitialValues")
    , time_step(mpi_communicator,
                computing_timer,
                offline_data,
                initial_values,
                "E - TimeStep")
    , schlieren_postprocessor(mpi_communicator,
                              computing_timer,
                              offline_data,
                              "F - SchlierenPostprocessor")
  {
    base_name = "test";
    add_parameter("basename", base_name, "Base name for all output files");

    t_final = 4.;
    add_parameter("final time", t_final, "Final time");

    output_granularity = 0.02;
    add_parameter("output granularity",
                  output_granularity,
                  "time interval for output");

    resume = false;
    add_parameter("resume", resume, "Resume an interrupted computation.");
  }

  // We define an auxiliary namespace to be used in the implementation of
  // the class member <code>TimeLoop<dim>::run()</code>. It's only content
  // is the void function <code>print_head</code> used to output
  // messages in the terminal with a "nice" format.

  namespace
  {
    void print_head(ConditionalOStream &pcout,
                    std::string         header,
                    std::string         secondary = "")
    {
      const auto header_size   = header.size();
      const auto padded_header = std::string((34 - header_size) / 2, ' ') +
                                 header +
                                 std::string((35 - header_size) / 2, ' ');

      const auto secondary_size = secondary.size();
      const auto padded_secondary =
        std::string((34 - secondary_size) / 2, ' ') + secondary +
        std::string((35 - secondary_size) / 2, ' ');

      /* clang-format off */
      pcout << std::endl;
      pcout << "    ####################################################" << std::endl;
      pcout << "    #########                                  #########" << std::endl;
      pcout << "    #########"     <<  padded_header   <<     "#########" << std::endl;
      pcout << "    #########"     << padded_secondary <<     "#########" << std::endl;
      pcout << "    #########                                  #########" << std::endl;
      pcout << "    ####################################################" << std::endl;
      pcout << std::endl;
      /* clang-format on */
    }
  } // namespace

  // The class member <code>TimeLoop<dim>::run()</code> is one of only three 
  // class member we actually have to implement. We initialize the 
  // (global) parameter list, setup all the accessory classes (discretization, 
  // offline_data, time_step, and schlieren_postprocessor), interpolate the 
  // initial data, and run a forward-Euler time loop.
  //
  // We note here that the (unique) call to ParameterAcceptor::initialize
  // initializes the global ParameterHandler with the
  // parameters contained in the classes derived from ParameterAceptor. 
  // This function enters the subsection returned by get_section_name() for 
  // each derived class, and declares all parameters that were added using 
  // add_parameter()

  template <int dim>
  void TimeLoop<dim>::run()
  {
    pcout << "Reading parameters and allocating objects... " << std::flush;

    /* Initialization of the global ParameterHandler. */
    ParameterAcceptor::initialize("step-69.prm");
    pcout << "done" << std::endl;

    print_head(pcout, "create triangulation");
    discretization.setup();

    print_head(pcout, "compute offline data");
    offline_data.setup();
    offline_data.assemble();

    print_head(pcout, "set up time step");
    time_step.prepare();
    schlieren_postprocessor.prepare();

    double       t            = 0.;
    unsigned int output_cycle = 0;

    print_head(pcout, "interpolate initial values");
    /* The vector U and time_step.temp are the only ones in the entire code
       storing the old and/or new state of the system. */
    auto U = interpolate_initial_values();

    /* By default resume is false, but that could have changed after reading 
       the input file when calling ParameterAcceptor::initialize */
    if (resume)
      {
        print_head(pcout, "restore interrupted computation");

        const auto &       triangulation = discretization.triangulation;
        const unsigned int i    = triangulation.locally_owned_subdomain();
        std::string        name = base_name + "-checkpoint-" +
                           Utilities::int_to_string(i, 4) + ".archive";
        std::ifstream file(name, std::ios::binary);

        boost::archive::binary_iarchive ia(file);
        ia >> t >> output_cycle;

        for (auto &it1 : U)
          {
            for (auto &it2 : it1)
              ia >> it2;
            it1.update_ghost_values();
          }
      }

    output(U, base_name + "-solution", t, output_cycle++);

    print_head(pcout, "enter main loop");

    for (unsigned int cycle = 1; t < t_final; ++cycle)
      {
        std::ostringstream head;
        head << "Cycle  " << Utilities::int_to_string(cycle, 6) << "  ("
             << std::fixed << std::setprecision(1) << t / t_final * 100 << "%)";
        std::ostringstream secondary;
        secondary << "at time t = " << std::setprecision(8) << std::fixed << t;
        print_head(pcout, head.str(), secondary.str());

        t += time_step.step(U, t);

        if (t > output_cycle * output_granularity)
          output(U, base_name + "-solution", t, output_cycle++, true);

      } /* End of time loop */

    if (output_thread.joinable())
      output_thread.join();

    computing_timer.print_summary();
    pcout << timer_output.str() << std::endl;
  }

  // Implementation of the class member <code>interpolate_initial_values</code>.
  // This function takes an initial time "t" as input argument in order to
  // evaluate an analytic expression (a function of space and time)
  // and returns a <code>vector_type</code> containing the initial values.

  template <int dim>
  typename TimeLoop<dim>::vector_type
  TimeLoop<dim>::interpolate_initial_values(double t)
  {
    pcout << "TimeLoop<dim>::interpolate_initial_values(t = " << t << ")"
          << std::endl;
    TimerOutput::Scope timer(computing_timer,
                             "time_loop - setup scratch space");

    vector_type U;

    const auto &partitioner = offline_data.partitioner;
    for (auto &it : U)
      it.reinit(partitioner);

    constexpr auto problem_dimension =
      ProblemDescription<dim>::problem_dimension;

    for (unsigned int i = 0; i < problem_dimension; ++i)
      VectorTools::interpolate(offline_data.dof_handler,
                               ScalarFunctionFromFunctionObject<dim, double>(
                                 [&](const auto &p) {
                                   return initial_values.initial_state(p, t)[i];
                                 }),
                               U[i]);

    for (auto &it : U)
      it.update_ghost_values();

    return U;
  }

  // Implementation of the class member <code>output</code>. Most of the 
  // following lines of code are invested in the implementation of the 
  // <code>output_worker</code> in order to write the output. We note that:
  //  - Before calling the <code>output_worker</code>, we create a copy of 
  //    <code>U[i]</code> (the vector we want to output). This copy is stored in
  //    <code>output_vector</code>.
  //  - the task <code>output_worker</code> is assigned to a thread 
  //  - this task is later moved to the thread <code>output_thread</code>.
  //
  // Since <code>output_vector</code> and <code>output_thread</code> are class 
  // members of <code>TimeLoop</code>, their scope extends beyond that one of 
  // anything defined inside <code>output_worker</code>. This allows the 
  // output task to continue its execution even when we 
  // <code>TimeLoop<dim>::output</code> releases its control to the function 
  // that called it. This is how (ideally) writing to disk becomes a 
  // background process and not a locking method.
  //
  // The only penalty is the copy of the vector we want to output. This 
  // penalty could be minimized by defining a class member 
  // TimeLoop<dim>::prepare() in order to allocate a priori the space for 
  // <code>output_vector</code> as we did with the vector <code>temp</code> in 
  // TimeStep<dim>::prepare().

  template <int dim>
  void TimeLoop<dim>::output(const typename TimeLoop<dim>::vector_type &U,
                             const std::string &                        name,
                             double                                     t,
                             unsigned int                               cycle,
                             bool checkpoint)
  {
    pcout << "TimeLoop<dim>::output(t = " << t
          << ", checkpoint = " << checkpoint << ")" << std::endl;

    /* We check if the thread is still running */
    /* If so, we wait to for it to join. */
    if (output_thread.joinable())
      {
        TimerOutput::Scope timer(computing_timer, "time_loop - stalled output");
        output_thread.join();
      }

    constexpr auto problem_dimension =
      ProblemDescription<dim>::problem_dimension;
    const auto &component_names = ProblemDescription<dim>::component_names;

    /* We make a copy the vector we want to output */
    for (unsigned int i = 0; i < problem_dimension; ++i)
      {
        output_vector[i] = U[i];
        output_vector[i].update_ghost_values();
      }

    schlieren_postprocessor.compute_schlieren(output_vector);

    /* We define the lambda function "output_worker" */
    const auto output_worker = [this, name, t, cycle, checkpoint]() {
      constexpr auto problem_dimension =
        ProblemDescription<dim>::problem_dimension;
      const auto &dof_handler   = offline_data.dof_handler;
      const auto &triangulation = discretization.triangulation;
      const auto &mapping       = discretization.mapping;

      if (checkpoint)
        {
          const unsigned int i    = triangulation.locally_owned_subdomain();
          std::string        name = base_name + "-checkpoint-" +
                             Utilities::int_to_string(i, 4) + ".archive";

          // FIXME: Refactor to Boost (this is C++17)
          // if (std::filesystem::exists(name))
          //   std::filesystem::rename(name, name + "~");

          std::ofstream file(name, std::ios::binary | std::ios::trunc);

          boost::archive::binary_oarchive oa(file);
          oa << t << cycle;
          for (const auto &it1 : output_vector)
            for (const auto &it2 : it1)
              oa << it2;
        }

      DataOut<dim> data_out;
      data_out.attach_dof_handler(dof_handler);

      for (unsigned int i = 0; i < problem_dimension; ++i)
        data_out.add_data_vector(output_vector[i], component_names[i]);

      data_out.add_data_vector(schlieren_postprocessor.schlieren,
                               "schlieren_plot");

      data_out.build_patches(mapping, discretization.finite_element.degree - 1);

      DataOutBase::VtkFlags flags(t,
                                  cycle,
                                  true,
                                  DataOutBase::VtkFlags::best_speed);
      data_out.set_flags(flags);

      data_out.write_vtu_with_pvtu_record("", name, cycle, 6, mpi_communicator);

       /* There is no return statement, we don't need it this is a void-like 
          lambda expression */
    }; 

    /* We launch the thread that executing the output and abandon the 
       function TimeLoop<dim>::output (returning the control to the 
       function that called it). */
    output_thread = std::move(std::thread(output_worker));
  } 

} /* End of namespace Step69 */

// @sect4{The main()}

int main(int argc, char *argv[])
{
  constexpr int dim = 2;

  using namespace dealii;
  using namespace Step69;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

  MPI_Comm      mpi_communicator(MPI_COMM_WORLD);
  TimeLoop<dim> time_loop(mpi_communicator);

  time_loop.run();
}

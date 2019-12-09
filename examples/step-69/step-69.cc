/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2012 - 2019 by the deal.II authors
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
// The set of include files is quite standard. The most intriguing part 
// is that: either though this code is a "thread and mpi parallel"
// we are using neither Trilinos nor PETSC vectors. Actually we are using dealii
// distributed vectors <code>la_parallel_vector.h</code> and the regular dealii
// sparse matrices <code>sparse_matrix.h</code>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/graph_coloring.h>
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

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/core/demangle.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/iterator_range.hpp>

#include <filesystem>


// @sect3{Declaration/s of the namespace Step69}

namespace Step69
{
  enum Boundary : dealii::types::boundary_id
  {
    do_nothing = 0,
    slip       = 2,
    dirichlet  = 3,
  };


  // @sect4{Declaration of <code>Discretization</code> class template}

  // The main goal of this class is to digest the input file and act as a
  // "container" of members that may be changed/decided at run time (through the
  // input file). It was natural to derive this class from
  // <code>dealii::ParameterAcceptor</code>. This class is in charge of
  // initializing mpi comunicator, geometry dimensions, triangulation, mapping,
  // finite element, mapping, and quadratures. If we think of the class
  // <code>Discretization</code> as a "container": it doesn't contain any
  // memmory demanding class member such a dof_handlers, vectors or matrices.
  // The most memmory thirsty class member is the <code>
  // dealii::parallel::distributed::Triangulation<dim></code>.

  template <int dim>
  class Discretization : public dealii::ParameterAcceptor
  {
  public:
    Discretization(const MPI_Comm &     mpi_communicator,
                   dealii::TimerOutput &computing_timer,
                   const std::string &  subsection = "Discretization");

    void setup();

    const MPI_Comm &mpi_communicator;

    dealii::parallel::distributed::Triangulation<dim> triangulation;

    const dealii::MappingQ<dim>   mapping;
    const dealii::FE_Q<dim>       finite_element;
    const dealii::QGauss<dim>     quadrature;
    const dealii::QGauss<dim - 1> face_quadrature;

  private:
    dealii::TimerOutput &computing_timer;

    double immersed_disc_length;
    double immersed_disc_height;
    double immersed_disc_object_position;
    double immersed_disc_object_diameter;

    unsigned int refinement;
  };

  // @sect4{Declaration of <code>OfflineData</code> class template}

  // The class OfflineData is initializes (and "owns")
  // pretty much all the components of the discretization that
  // do not evolve in time. In particular: dof_handler, sparsity
  // patterns, boundary maps, lumped mass matrix, and other matrices
  // and vectors that do not change in time are members of this class.
  // The term "offline" here refers to the idea that all the class members
  // of <code>OfflineData</code> are initialized and assigned values
  // a "time step zero" and are not meant to be modified at any other later
  // time step. For instance, the sparsity pattern should not
  // change as we advance in time (we are not doing any form of adaptivity in
  // space). Similarly, the entries of the vector
  // <code>lumped_mass_matrix</code> should not be modified as we advance in
  // time either.
  //
  // Placeholder: Say something about BoundaryNormalMap.

  template <int dim>
  class OfflineData : public dealii::ParameterAcceptor
  {
  public:
    using BoundaryNormalMap =
      std::map<dealii::types::global_dof_index,
               std::tuple<dealii::Tensor<1, dim> /* normal */,
                          dealii::types::boundary_id /* id */,
                          dealii::Point<dim>> /* position */>;

    OfflineData(const MPI_Comm &           mpi_communicator,
                dealii::TimerOutput &      computing_timer,
                const Discretization<dim> &discretization,
                const std::string &        subsection = "OfflineData");

    void setup();
    void assemble();

    dealii::DoFHandler<dim> dof_handler;

    std::shared_ptr<const dealii::Utilities::MPI::Partitioner> partitioner;

    unsigned int n_locally_owned;
    unsigned int n_locally_relevant;

    dealii::SparsityPattern sparsity_pattern;

    BoundaryNormalMap boundary_normal_map;

    dealii::SparseMatrix<double>                  lumped_mass_matrix;
    std::array<dealii::SparseMatrix<double>, dim> cij_matrix;
    std::array<dealii::SparseMatrix<double>, dim> nij_matrix;
    dealii::SparseMatrix<double>                  norm_matrix;

  private:
    const MPI_Comm &     mpi_communicator;
    dealii::TimerOutput &computing_timer;

    dealii::SmartPointer<const Discretization<dim>> discretization;
  };

  // @sect4{Declaration of <code>ProblemDescription</code> class template}

  // Most of the implementations of the members of this class will be utility
  // classes/functions specific for Euler's equations:
  // - The type alias <code>rank1_type</code> will be used for the states
  // $\mathbf{U}_i^n$
  // - The type alias <code>rank2_type</code> will be used for the fluxes
  // $\mathbb{f}(\mathbf{U}_j^n)$.
  // - The implementation of <code>momentum</code> will extract $\textbf{m}$
  //   (out of the state vector $[\rho,\textbf{m},E]$) and store it in a
  //   <code>Tensor<1, dim></code> for our convenience.
  // - The implementation of <code>internal_energy</code> will compute
  //   $E - \frac{|\textbf{m}|^2}{2\rho}$ from the state vector
  //   $[\rho,\textbf{m},E]$.
  //
  // The purpose of the remaining class members <code>component_names</code>,
  // <code>pressure</code>, and <code>speed_of_sound</code>,
  // is evident from their names. Most notably, the last
  // one <code>compute_lambda_max</code> is in charge of computing
  // $\lambda_{max}(\mathbf{U},\mathbf{V},\mathbf{n})$ which is required
  // to compute the first order viscosity $d_{ij}$ as detailed in the section
  // <b>Description of the scheme</b>.

  template <int dim>
  class ProblemDescription
  {
  public:
    static constexpr unsigned int problem_dimension = 2 + dim;

    using rank1_type = dealii::Tensor<1, problem_dimension>;
    using rank2_type =
      dealii::Tensor<1, problem_dimension, dealii::Tensor<1, dim>>;

    const static std::array<std::string, dim + 2> component_names;

    static constexpr double gamma = 7. / 5.;

    static DEAL_II_ALWAYS_INLINE inline dealii::Tensor<1, dim>
    momentum(const rank1_type U);

    static DEAL_II_ALWAYS_INLINE inline double
    internal_energy(const rank1_type U);

    static DEAL_II_ALWAYS_INLINE inline double pressure(const rank1_type U);

    static DEAL_II_ALWAYS_INLINE inline double
    speed_of_sound(const rank1_type U);

    static DEAL_II_ALWAYS_INLINE inline rank2_type f(const rank1_type U);

    static DEAL_II_ALWAYS_INLINE inline double
    compute_lambda_max(const rank1_type              U_i,
                       const rank1_type              U_j,
                       const dealii::Tensor<1, dim> &n_ij);
  };

  // @sect4{Declaration of <code>InitialValues</code> class template}

  // Placeholder here

  template <int dim>
  class InitialValues : public dealii::ParameterAcceptor
  {
  public:
    using rank1_type = typename ProblemDescription<dim>::rank1_type;

    InitialValues(const std::string &subsection = "InitialValues");

    std::function<rank1_type(const dealii::Point<dim> &point, double t)>
      initial_state;

  private:
    void parse_parameters_callback();

    dealii::Tensor<1, dim> initial_direction;
    dealii::Tensor<1, 3>   initial_1d_state;
  };

  // @sect4{Declaration of <code>TimeStep</code> class template}

  // Placeholder here

  template <int dim>
  class TimeStep : public dealii::ParameterAcceptor
  {
  public:
    static constexpr unsigned int problem_dimension =
      ProblemDescription<dim>::problem_dimension;

    using rank1_type = typename ProblemDescription<dim>::rank1_type;
    using rank2_type = typename ProblemDescription<dim>::rank2_type;

    typedef std::array<dealii::LinearAlgebra::distributed::Vector<double>,
                       problem_dimension>
      vector_type;

    TimeStep(const MPI_Comm &          mpi_communicator,
             dealii::TimerOutput &     computing_timer,
             const OfflineData<dim> &  offline_data,
             const InitialValues<dim> &initial_values,
             const std::string &       subsection = "TimeStep");

    void prepare();

    double step(vector_type &U, double t);

  private:
    const MPI_Comm &     mpi_communicator;
    dealii::TimerOutput &computing_timer;

    dealii::SmartPointer<const OfflineData<dim>>   offline_data;
    dealii::SmartPointer<const InitialValues<dim>> initial_values;

    dealii::SparseMatrix<double> dij_matrix;

    vector_type temp;

    double cfl_update;
  };

  // @sect4{Declaration of <code>SchlierenPostprocessor</code> class template}

  // At its core, the Schilieren class implements the class member
  // <code>compute_schlieren</code>. The main purpose of this class member
  // is to compute auxiliary finite element field <code>schlieren</code>
  // at each node, defined as
  // \f[ \text{schlieren}[i] = e^{\beta \frac{ |\nabla r_i|
  // - \min_j |\nabla r_j| }{\max_j |\nabla r_j| - \min_j |\nabla r_j| } } \f]
  // where $r$ in principle could be any scalar finite element field.
  // The natural candidate is choosing $r := \rho$. Schlieren postprocessing
  // is a standard methodology to enhance the contrast of the visualization
  // inspired in actual X-ray and shadowgraphy experimental techniques of
  // visualization.

  template <int dim>
  class SchlierenPostprocessor : public dealii::ParameterAcceptor
  {
  public:
    static constexpr unsigned int problem_dimension =
      ProblemDescription<dim>::problem_dimension;

    using rank1_type = typename ProblemDescription<dim>::rank1_type;

    using vector_type =
      std::array<dealii::LinearAlgebra::distributed::Vector<double>,
                 problem_dimension>;

    SchlierenPostprocessor(
      const MPI_Comm &        mpi_communicator,
      dealii::TimerOutput &   computing_timer,
      const OfflineData<dim> &offline_data,
      const std::string &     subsection = "SchlierenPostprocessor");

    void prepare();

    void compute_schlieren(const vector_type &U);

    dealii::LinearAlgebra::distributed::Vector<double> schlieren;

  private:
    const MPI_Comm &     mpi_communicator;
    dealii::TimerOutput &computing_timer;

    dealii::SmartPointer<const OfflineData<dim>> offline_data;

    dealii::Vector<double> r;

    unsigned int schlieren_index;
    double       schlieren_beta;
  };

  // @sect4{Declaration of <code>TimeLoop</code> class template}

  // Placeholder here

  template <int dim>
  class TimeLoop : public dealii::ParameterAcceptor
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

    const MPI_Comm &    mpi_communicator;
    std::ostringstream  timer_output;
    dealii::TimerOutput computing_timer;

    dealii::ConditionalOStream pcout;

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

} // namespace Step69


// @sect3{Implementation of the classes in namespace <code>Step69</code>}

namespace Step69
{
  using namespace dealii;

  // @sect4{Implementation of the members of the class <code>Discretization</code>}

  // Not much is done here other that initializing the corresponding
  // class members in the initialization list.

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
    immersed_disc_length = 4.;
    add_parameter("immersed disc - length",
                  immersed_disc_length,
                  "Immersed disc: length of computational domain");

    immersed_disc_height = 2.;
    add_parameter("immersed disc - height",
                  immersed_disc_height,
                  "Immersed disc: height of computational domain");

    immersed_disc_object_position = 0.6;
    add_parameter("immersed disc - object position",
                  immersed_disc_object_position,
                  "Immersed disc: x position of immersed disc center point");

    immersed_disc_object_diameter = 0.5;
    add_parameter("immersed disc - object diameter",
                  immersed_disc_object_diameter,
                  "Immersed disc: diameter of immersed disc");

    refinement = 5;
    add_parameter("initial refinement",
                  refinement,
                  "Initial refinement of the geometry");
  }

  // Note that in the previous constructor we only passed the MPI
  // communicator to the <code>triangulation</code>but we still have not
  // initialized the underlying geometry/mesh. In order to define the geometry
  // we will use the class <code>create_immersed_disc_geometry</code>
  // that uses the tools in GridGenerator in order to create a
  // rectangular domain with a whole.

  // The following is just a dummy implementation/placeholder that does
  // nothing other than throwing an exception if we want to run this program
  // with a space dimension that is not 2.

  template <int dim>
  void
  create_immersed_disc_geometry(parallel::distributed::Triangulation<dim> &,
                                const double /*length*/,
                                const double /*height*/,
                                const double /*step_position*/,
                                const double /*step_height*/)
  {
    AssertThrow(false, ExcNotImplemented());
  }

  // For the two-dimensional case we have the following template
  // specialization that creates the geometry.

  template <>
  void create_immersed_disc_geometry<2>(
    parallel::distributed::Triangulation<2> &triangulation,
    const double                             length,
    const double                             height,
    const double                             disc_position,
    const double                             disc_diameter)
  {
    constexpr int dim = 2;

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

    for (auto cell : triangulation.active_cell_iterators())
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
        {
          auto &vertex = cell->vertex(v);
          if (vertex[0] <= -disc_diameter + 1.e-6)
            vertex[0] = -disc_position;
        }

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
  }

  // This the the last class member to be implemented of the class.
  // <code>Discretization</code>: it initializes the actual mesh
  // of the triangulation by calling <code>create_immersed_disc_geometry</code>.

  template <int dim>
  void Discretization<dim>::setup()
  {
    TimerOutput::Scope t(computing_timer, "discretization - setup");

    triangulation.clear();

    create_immersed_disc_geometry(triangulation,
                                  immersed_disc_length,
                                  immersed_disc_height,
                                  immersed_disc_object_position,
                                  immersed_disc_object_diameter);

    triangulation.refine_global(refinement);
  }


  // @sect4{Implementation of the members of the class <code>OfflineData</code>}

  // Not much is done here other that initializing the corresponding
  // class members in the initialization list.
  // Constructor of the class <code>OfflineData</code>.

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

  // Now the class member <code>OfflineData<dim>::setup()</code> will take care
  // of initializating
  //   - The <code>dof_handler</code>.
  //   - The IndexSets corresponding to locally owned and locally relevant DOFs.
  //   - The partitioner.

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
    }

    {
      TimerOutput::Scope t(
        computing_timer,
        "offline_data - create partitioner and affine constraints");

      locally_relevant.clear();
      DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant);
      n_locally_relevant = locally_relevant.n_elements();

      partitioner.reset(new Utilities::MPI::Partitioner(locally_owned,
                                                        locally_relevant,
                                                        mpi_communicator));
    }

    const auto dofs_per_cell = discretization->finite_element.dofs_per_cell;

    // Here we create the sparsity patterns for the off-line data. There are
    // quite a few peculiarities that deserve our attention:
    // - Our "local" dynamic sparsity pattern (<code>dsp</code>)
    //   will be of dimensions
    //   <code>n_locally_relevant</code> $\times$
    //   <code>n_locally_relevant</code> (this choice is definitely unusual).
    //   The goal behind such choice is to reduce communication: we do
    //   not want to request (to another mpi process) ghosted offline data (such
    //   as the vectors $\mathbf{c}_{ij}$ when $j$ is not locally owned) for
    //   every time step. It is more efficient to simply take more memory by
    //   storing (locally) all relevant off-line data.
    // - We loop on all locally owned and ghosted cells (see @ref
    //   GlossArtificialCell "this glossary entry") order to:
    //   <ul>
    //   <li> Extract the <code>dof_indices</code> associated to the cell DOFs
    //   (having global numbering) and renumber them using
    //   <code>partitioner->global_to_local(index)</code>. For the case
    //   of locally owned DOFs: such renumbering consist in applying a
    //   shift (i.e. we subtract a number) such that now they will
    //   become a number in the integer interval
    //   $[0,$<code>n_locally_owned</code>$)$. However, for the case of
    //   "ghosted DOFs" (i.e. not locally owned) the situation is quite
    //   different, since the global indices associated to ghosted DOFs
    //   will not be (in general) a contiguous set of integers.
    //   </li>
    //   <li> Once, we are done with that, we add the corresponding entries to
    //   the rows of the dynamic sparsity pattern with
    //   <code>dsp.add_entries</code></li>
    //   </ul>
    // Finally we use <code>dsp</code> to initialize the actual sparsity
    // pattern <code>sparsity_pattern</code>.

    {
      TimerOutput::Scope t(computing_timer,
                           "offline_data - create sparsity pattern");

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
    }

    // We initialize the off-line data matrices. Note that these matrices do
    // not require an mpi communicator (that's the idea).

    {
      TimerOutput::Scope t(computing_timer, "offline_data - set up matrices");

      lumped_mass_matrix.reinit(sparsity_pattern);
      norm_matrix.reinit(sparsity_pattern);
      for (auto &matrix : cij_matrix)
        matrix.reinit(sparsity_pattern);
      for (auto &matrix : nij_matrix)
        matrix.reinit(sparsity_pattern);
    }
  }

  // In this last brace we finished with the implementation of the
  // <code>OfflineData<dim>::setup()</code>.
  //
  // Now we define a collection of assembly utilities:
  // - <code>CopyData</code>: This will only be used to compute the off-line
  //   data using WorkStream. It acts as a container: it is just a
  //  struct where WorkStream stores the local cell contributions.
  // - <code>get_entry</code>: it reads the value stored at the entry
  //  pointed by the iterator <code>it</code> of <code>matrix</code>. Here is
  //  where we might want to keep an eye on complexity: we want this operation
  //  to have constant complexity (that's the case of this implementation).
  //  Note also that the return argument (<code>Matrix::value_type</code>) is
  //  going to be (in general) a double.
  // - <code>set_entry</code>: it sets <code>value</code> at the entry
  //  pointed by the iterator <code>it</code> of <code>matrix</code>.
  // - <code>gather_get_entry</code>: we note that
  //   $\mathbf{c}_{ij} \in \mathbb{R}^d$. If $d=2$ then
  //   $\mathbf{c}_{ij} = [\mathbf{c}_{ij}^1,\mathbf{c}_{ij}^2]^\top$.
  //   Which basically implies
  //   that we need one matrix per space dimension to store the
  //  $\mathbf{c}_{ij}$ vectors. Similar observation follows for the matrix
  //  $\mathbf{n}_{ij}$. The purpose of <code>gather_get_entry</code>
  //  is to retrieve those entries a store them into a
  //  <code>Tensor<1, dim></code> for our convenience.
  // - <code>gather</code> (first interface):
  // - <code>gather</code> (second interface):
  // - <code>scatter</code>:

  namespace
  {
    template <int dim>
    struct CopyData
    {
      bool                                         is_artificial;
      std::vector<dealii::types::global_dof_index> local_dof_indices;
      typename OfflineData<dim>::BoundaryNormalMap local_boundary_normal_map;
      dealii::FullMatrix<double>                   cell_lumped_mass_matrix;
      std::array<dealii::FullMatrix<double>, dim>  cell_cij_matrix;
    };


    template <typename Matrix, typename Iterator>
    DEAL_II_ALWAYS_INLINE inline typename Matrix::value_type
    get_entry(const Matrix &matrix, const Iterator &it)
    {
      const auto                            global_index = it->global_index();
      const typename Matrix::const_iterator matrix_iterator(&matrix,
                                                            global_index);
      return matrix_iterator->value();
    }


    template <typename Matrix, typename Iterator>
    inline DEAL_II_ALWAYS_INLINE void
    set_entry(Matrix &                    matrix,
              const Iterator &            it,
              typename Matrix::value_type value)
    {
      const auto                global_index = it->global_index();
      typename Matrix::iterator matrix_iterator(&matrix, global_index);
      matrix_iterator->value() = value;
    }


    template <typename T1, std::size_t k, typename T2>
    DEAL_II_ALWAYS_INLINE inline dealii::Tensor<1, k>
    gather_get_entry(const std::array<T1, k> &U, const T2 it)
    {
      dealii::Tensor<1, k> result;
      for (unsigned int j = 0; j < k; ++j)
        result[j] = get_entry(U[j], it);
      return result;
    }


    template <typename T1, std::size_t k, typename T2, typename T3>
    DEAL_II_ALWAYS_INLINE inline dealii::Tensor<1, k>
    gather(const std::array<T1, k> &U, const T2 i, const T3 l)
    {
      dealii::Tensor<1, k> result;
      for (unsigned int j = 0; j < k; ++j)
        result[j] = U[j](i, l);
      return result;
    }


    template <typename T1, std::size_t k, typename T2>
    DEAL_II_ALWAYS_INLINE inline dealii::Tensor<1, k>
    gather(const std::array<T1, k> &U, const T2 i)
    {
      dealii::Tensor<1, k> result;
      for (unsigned int j = 0; j < k; ++j)
        result[j] = U[j].local_element(i);
      return result;
    }


    template <typename T1, std::size_t k1, typename T2, typename T3>
    DEAL_II_ALWAYS_INLINE inline void
    scatter(std::array<T1, k1> &U, const T2 &result, const T3 i)
    {
      for (unsigned int j = 0; j < k1; ++j)
        U[j].local_element(i) = result[j];
    }
  } // end of namespace.

  // The following piece of code implements the class member
  // <code>OfflineData<dim>::assemble()</code> which (in short)
  // computes the lumped mass entries $m_i$, the vectors $\mathbf{c}_{ij}$,
  // the vector $\mathbf{n}_{ij} = \frac{\mathbf{c}_{ij}}{|\mathbf{c}_{ij}|}$,
  // and the boundary normals $\boldsymbol{\nu}_i$.
  //
  // In order to exploit thread parallelization we use WorkStream approach
  // detailed in the @ref threads "Parallel computing with multiple processors
  // accessing shared memory". As customary this requires
  // definition of
  //  - Scratch data: in this case it is <code>scratch_data</code>.
  //  - The worker: in the case it is <code>local_assemble_system</code> that
  //    actually computes the local (i.e. current cell) contributions.
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
  // Finally the boundary normals are defined as
  // $\widehat{\boldsymbol{\nu}}_i =
  // \frac{\boldsymbol{\nu}_i}{|\boldsymbol{\nu}_i|}$ where
  // $\boldsymbol{\nu}_i = \sum_{F \subset \text{supp}(\phi_i)}
  // \sum_{\mathbf{x}_{q,F}} \nu(\mathbf{x}_{q,F})
  // \phi_i(\mathbf{x}_{q,F})$, here: $F \subset \partial \Omega$ denotes 
  // faces of elements at the boundary of the domain, and $\mathbf{x}_{q,F}$ 
  // are quadrature points on such face.
  // Other more sophisticated definitions for $\nu_i$ are 
  // possible but none of them have much influence in theory or practice.
  // We remind the reader that <code>CopyData</code> includes the class member
  // <code>local_boundary_normal_map</code> in order to store these local
  // contributions for the boundary map.

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

        /* We compute the local contributions for the  lumped mass
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
           following loop does not actually do much unless the faces of the
           cell are actually faces on the boundary of the domain */
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
                   function \phi_j. So we cannot normalize this local 
                   contribution right here, we have to take it "as is" and pass 
                   it to the copy data routine. */
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
            auto &[normal, id, position] = boundary_normal_map[it.first];
            auto &[new_normal, new_id, new_position] = it.second;

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
    // want: we have to normalize its entries. In addition, we have not even
    // touched the entries of the matrix <code>norm_matrix</code> yet, and the 
    // vectors stored in the map
    // <code>OfflineData<dim>::BoundaryNormalMap</code> are not normalized.
    //
    // In principle, this is just offline data, it doesn't make much sense
    // to over-optimize their computation, since their cost will get amortized
    // over the many time steps that we are going to use. However, 
    // computing/storing the entries of the matrix
    // <code>norm_matrix</code> and the normalization of <code>nij_matrix</code>
    // are perfect to illustrate thread-parallel node-loops: 
    // - We want to visit every node $i$ in the mesh/sparsity graph, 
    // - and for every such node we want to visit to every $j$ such that
    // $\mathbf{c}_{ij} \not \equiv 0$.
    //
    // From an algebraic point of view, this is equivalent to: visiting 
    // every row in the matrix (equivalently sparsity
    // pattern) and for each one of these rows execute a loop on the columns.
    // Node-loops is a core theme of this tutorial step (see the pseudo-code
    // in the introduction) that will repeat over and over again. That's why 
    // this is the right time to introduce them.
    //
    // We have the thread paralellization capability
    // parallel::apply_to_subranges that is somehow more general than the
    // WorkStream framework. In particular, it can be used for our 
    // node-loops.
    // This functionality requires four input arguments:
    // - A begin iterator: <code>indices.begin()</code>
    // - A end iterator: <code>indices.end()</code>
    // - A function f(i1,i2), where <code>i1</code> and <code>i2</code> define a
    //   sub-range with the range spanned by the the end and begin iterators
    //   of the previous two bullets. The function <code>f(i1,i2)</code> is
    //   called <code>on_subranges</code> in this example. It applies an
    //   operation for every "abstract element" in the subrange. In this case
    //   each "element" is a row of the sparsity pattern.
    // - Grainsize: minimum number of "elements" (in this case rows) processed
    // by
    //   each thread. We decided for a minimum of 4096 rows.
    //
    // We start by defining the operation <code>on_subranges</code> to be
    // applied at each row in the sub-range. Given a fixed
    // <code>row_index</code> we want to visit every entry in such row. In order
    // to execute such columns-loops we use <a
    // href="http://www.cplusplus.com/reference/algorithm/for_each/">
    // std::for_each</a>
    // from the standard library, where:
    // <code>sparsity_pattern.begin(row_index)</code>
    // gives us an iterator starting at the first column,
    // <code>sparsity_pattern.end(row_index)</code> is an iterator pointing at
    // the last column of the row. The last
    // argument required by std::for_each is the operation applied at each
    // column (a lambda expression in this case) of such row. We note that
    // because of the nature of the data that we want to modify (we want to
    // modify entries of a entire row at a time) threads cannot collide
    // attempting to write the same entry (we do not need a scheduler). This
    // advantage appears to be a particular characteristic of edge-based finite
    // element schemes when they are properly implemented.
    //
    // Finally, we normalize the vector stored in
    // <code>OfflineData<dim>::BoundaryNormalMap</code>. This operation has 
    // not been thread paralellized as it would not illustrate any important 
    // concept.

    {
      TimerOutput::Scope t(computing_timer,
                           "offline_data - compute |c_ij|, and n_ij");

      /* Here [i1,i2] represent a subrange of rows */
      const auto on_subranges = [&](auto i1, const auto i2) {
        for (; i1 < i2; ++i1)
          {
            const auto row_index = *i1;

            /* First column-loop: we compute/store the entries of the matrix
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

      /* We normalize the normals at the boundary. */
      /* This is not thread parallelized, too bad! */
      for (auto &it : boundary_normal_map)
        {
          auto &[normal, id, _] = it.second;
          normal /= (normal.norm() + std::numeric_limits<double>::epsilon());
        }
    }

    // In order to implement reflecting boundary conditions
    // $\mathbf{m} \cdot \boldsymbol{\nu}_i =0$ (or equivalently $\mathbf{v}
    // \cdot \boldsymbol{\nu}_i =0$ ) the vectors $\mathbf{c}_{ij}$ at the 
    // boundary have to be modified as:
    //
    // $\mathbf{c}_{ij} += \int_{\partial \Omega}
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
  //
  // Now we define the implementation of <code>momentum</code>,
  // <code>internal_energy</code>, <code>pressure</code>,
  // <code>speed_of_sound</code>, and <code>f</code> (the flux of the system).
  // The functionality of each one of these functions is self-explanatory from 
  // their names.

  template <int dim>
  DEAL_II_ALWAYS_INLINE inline dealii::Tensor<1, dim>
  ProblemDescription<dim>::momentum(const rank1_type U)
  {
    dealii::Tensor<1, dim> result;
    std::copy(&U[1], &U[1 + dim], &result[0]);
    return result;
  }

  template <int dim>
  DEAL_II_ALWAYS_INLINE inline double
  ProblemDescription<dim>::internal_energy(const rank1_type U)
  {
    const double &rho = U[0];
    const auto    m   = momentum(U);
    const double &E   = U[dim + 1];
    return E - 0.5 * m.norm_square() / rho;
  }

  template <int dim>
  DEAL_II_ALWAYS_INLINE inline double
  ProblemDescription<dim>::pressure(const rank1_type U)
  {
    return (gamma - 1.) * internal_energy(U);
  }

  template <int dim>
  DEAL_II_ALWAYS_INLINE inline double
  ProblemDescription<dim>::speed_of_sound(const rank1_type U)
  {
    const double &rho = U[0];
    const double  p   = pressure(U);

    return std::sqrt(gamma * p / rho);
  }

  template <int dim>
  DEAL_II_ALWAYS_INLINE inline typename ProblemDescription<dim>::rank2_type
  ProblemDescription<dim>::f(const rank1_type U)
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

  // The following function, <code>riemann_data_from_state</code>, takes the 
  // full state $\mathbf{u} = [\rho,\mathbf{m},E]^\top$ defines a new 
  // "projected state" defined as
  //
  // $\widetilde{\mathbf{u}} = [\rho,
  // \mathbf{m} - (\mathbf{m}\cdot \mathbf{n}_{ij})\mathbf{n}_{ij},
  //  E - \tfrac{(\mathbf{m}\cdot \mathbf{n}_{ij})^2}{2\rho} ]^\top$
  //
  // Projected states appear naturally when attempting to compute a maximum 
  // wavespeed appearing in Riemann problems.

  namespace
  {
    template <int dim>
    DEAL_II_ALWAYS_INLINE inline std::array<double, 4> riemann_data_from_state(
      const typename ProblemDescription<dim>::rank1_type U,
      const dealii::Tensor<1, dim> &                     n_ij)
    {
      dealii::Tensor<1, 3> projected_U;
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


    DEAL_II_ALWAYS_INLINE inline double
    lambda3_plus(const std::array<double, 4> &riemann_data, const double p_star)
    {
      constexpr double gamma             = ProblemDescription<1>::gamma;
      const auto &[rho_Z, u_Z, p_Z, a_Z] = riemann_data;

      const double factor = (gamma + 1.0) / 2.0 / gamma;
      const double tmp    = positive_part((p_star - p_Z) / p_Z);
      return u_Z + a_Z * std::sqrt(1.0 + factor * tmp);
    }


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

      const double p_star =
        p_j * std::pow(numerator / denominator, 2. * gamma / (gamma - 1));

      const double lambda1 = lambda1_minus(riemann_data_i, p_star);
      const double lambda3 = lambda3_plus(riemann_data_j, p_star);

      return std::max(positive_part(lambda3), negative_part(lambda1));
    };


    DEAL_II_ALWAYS_INLINE inline double
    lambda_max_expansion(const std::array<double, 4> &riemann_data_i,
                         const std::array<double, 4> &riemann_data_j)
    {
      const auto &[rho_i, u_i, p_i, a_i] = riemann_data_i;
      const auto &[rho_j, u_j, p_j, a_j] = riemann_data_j;

      return std::max(std::abs(u_i), std::abs(u_j)) + 5. * std::max(a_i, a_j);
    }
  } /* End of namespace dedicated to the computation of the maximum wavespeed */

  // Placeholder here.

  template <int dim>
  DEAL_II_ALWAYS_INLINE inline double
  ProblemDescription<dim>::compute_lambda_max(
    const rank1_type              U_i,
    const rank1_type              U_j,
    const dealii::Tensor<1, dim> &n_ij)
  {
    const auto riemann_data_i = riemann_data_from_state(U_i, n_ij);
    const auto riemann_data_j = riemann_data_from_state(U_j, n_ij);

    const double lambda_1 =
      lambda_max_two_rarefaction(riemann_data_i, riemann_data_j);

    const double lambda_2 =
      lambda_max_expansion(riemann_data_i, riemann_data_j);

    return std::min(lambda_1, lambda_2);
  }

  // Placeholder here.

  template <>
  const std::array<std::string, 3> //
    ProblemDescription<1>::component_names{"rho", "m", "E"};

  template <>
  const std::array<std::string, 4> //
    ProblemDescription<2>::component_names{"rho", "m_1", "m_2", "E"};

  template <>
  const std::array<std::string, 5> //
    ProblemDescription<3>::component_names{"rho", "m_1", "m_2", "m_3", "E"};

  // Placeholder here.

  template <int dim>
  InitialValues<dim>::InitialValues(const std::string &subsection)
    : ParameterAcceptor(subsection)
  {
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

  // Placeholder here.

  template <int dim>
  void InitialValues<dim>::parse_parameters_callback()
  {
    AssertThrow(initial_direction.norm() != 0.,
                ExcMessage(
                  "Initial shock front direction is set to the zero vector."));
    initial_direction /= initial_direction.norm();

    static constexpr auto gamma = ProblemDescription<dim>::gamma;

    const auto from_1d_state =
      [=](const dealii::Tensor<1, 3, double> &state_1d) -> rank1_type {
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

    initial_state = [=](const dealii::Point<dim> & /*point*/, double /*t*/) {
      return from_1d_state(initial_1d_state);
    };
  }

  // Placeholder here.

  template <int dim>
  TimeStep<dim>::TimeStep(const MPI_Comm &          mpi_communicator,
                          dealii::TimerOutput &     computing_timer,
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

  // Placeholder here.

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

  // Placeholder here.

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

      const auto on_subranges = [&](auto i1, const auto i2) {
        for (const auto i : boost::make_iterator_range(i1, i2))
          {
            const auto U_i = gather(U, i);

            /* Column-loop */
            for (auto jt = sparsity.begin(i); jt != sparsity.end(i); ++jt)
              {
                const auto j = jt->column();

                if (j >= i)
                  continue;

                const auto U_j = gather(U, j);

                const auto   n_ij = gather_get_entry(nij_matrix, jt);
                const double norm = get_entry(norm_matrix, jt);

                const auto lambda_max =
                  ProblemDescription<dim>::compute_lambda_max(U_i, U_j, n_ij);

                double d = norm * lambda_max;

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

                set_entry(dij_matrix, jt, d);
                dij_matrix(j, i) = d;
              } /* End of column-loop */
          }     /* End of row-loop */
      };        /* End of definition of on_subranges */

      parallel::apply_to_subranges(indices_relevant.begin(),
                                   indices_relevant.end(),
                                   on_subranges,
                                   4096);
    } /* End of the computation of the off-diagonal entries of dij_matrix */

    std::atomic<double> tau_max{std::numeric_limits<double>::infinity()};

    {
      TimerOutput::Scope time(computing_timer,
                              "time_step - 2 compute d_ii, and tau_max");

      const auto on_subranges = [&](auto i1, const auto i2) {
        double tau_max_on_subrange = std::numeric_limits<double>::infinity();

        for (const auto i : boost::make_iterator_range(i1, i2))
          {
            double d_sum = 0.;

            for (auto jt = sparsity.begin(i); jt != sparsity.end(i); ++jt)
              {
                const auto j = jt->column();

                if (j == i)
                  continue;

                d_sum -= get_entry(dij_matrix, jt);
              }

            dij_matrix.diag_element(i) = d_sum;

            const double mass   = lumped_mass_matrix.diag_element(i);
            const double tau    = cfl_update * mass / (-2. * d_sum);
            tau_max_on_subrange = std::min(tau_max_on_subrange, tau);
          }

        double current_tau_max = tau_max.load();
        while (
          current_tau_max > tau_max_on_subrange &&
          !tau_max.compare_exchange_weak(current_tau_max, tau_max_on_subrange))
          ;
      };

      parallel::apply_to_subranges(indices_relevant.begin(),
                                   indices_relevant.end(),
                                   on_subranges,
                                   4096);

      tau_max.store(Utilities::MPI::min(tau_max.load(), mpi_communicator));

      AssertThrow(!std::isnan(tau_max) && !std::isinf(tau_max) && tau_max > 0.,
                  ExcMessage("I'm sorry, Dave. I'm afraid I can't "
                             "do that. - We crashed."));
    } /* End of the computation of the diagonal entries of dij_matrix */

    {
      TimerOutput::Scope time(computing_timer, "time_step - 3 perform update");

      const auto on_subranges = [&](auto i1, const auto i2) {
        for (const auto i : boost::make_iterator_range(i1, i2))
          {
            Assert(i < n_locally_owned, ExcInternalError());

            const auto U_i = gather(U, i);

            const auto   f_i = ProblemDescription<dim>::f(U_i);
            const double m_i = lumped_mass_matrix.diag_element(i);

            auto U_i_new = U_i;

            for (auto jt = sparsity.begin(i); jt != sparsity.end(i); ++jt)
              {
                const auto j = jt->column();

                const auto U_j = gather(U, j);
                const auto f_j = ProblemDescription<dim>::f(U_j);

                const auto c_ij = gather_get_entry(cij_matrix, jt);
                const auto d_ij = get_entry(dij_matrix, jt);

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

      /* Only iterate over locally owned subset! */
      parallel::apply_to_subranges(indices_owned.begin(),
                                   indices_owned.end(),
                                   on_subranges,
                                   4096);
    } /* End of the computation of the new solution */

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

            const auto &[normal, id, position] = it->second;

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

            if (id == Boundary::dirichlet)
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

  // Placeholder here.

  template <int dim>
  SchlierenPostprocessor<dim>::SchlierenPostprocessor(
    const MPI_Comm &        mpi_communicator,
    dealii::TimerOutput &   computing_timer,
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

    std::atomic<double> r_i_max{0.};
    std::atomic<double> r_i_min{std::numeric_limits<double>::infinity()};

    {
      const auto on_subranges = [&](auto i1, const auto i2) {
        double r_i_max_on_subrange = 0.;
        double r_i_min_on_subrange = std::numeric_limits<double>::infinity();

        for (; i1 < i2; ++i1)
          {
            const auto i = *i1;

            Assert(i < n_locally_owned, ExcInternalError());

            Tensor<1, dim> r_i;

            for (auto jt = sparsity.begin(i); jt != sparsity.end(i); ++jt)
              {
                const auto j = jt->column();

                if (i == j)
                  continue;

                const auto U_js = U[schlieren_index].local_element(j);
                const auto c_ij = gather_get_entry(cij_matrix, jt);

                r_i += c_ij * U_js;
              }

            const auto bnm_it = boundary_normal_map.find(i);
            if (bnm_it != boundary_normal_map.end())
              {
                const auto [normal, id, _] = bnm_it->second;
                if (id == Boundary::slip)
                  {
                    r_i -= 1. * (r_i * normal) * normal;
                  }
                else
                  {
                    r_i = 0.;
                  }
              }

            const double m_i = lumped_mass_matrix.diag_element(i);
            r[i]             = r_i.norm() / m_i;

            r_i_max_on_subrange = std::max(r_i_max_on_subrange, r[i]);
            r_i_min_on_subrange = std::min(r_i_min_on_subrange, r[i]);
          }

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

    {
      const auto on_subranges = [&](auto i1, const auto i2) {
        for (; i1 < i2; ++i1)
          {
            const auto i = *i1;

            Assert(i < n_locally_owned, ExcInternalError());

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

  // Placeholder here.

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

  // Placeholder here.

  namespace
  {
    void print_head(dealii::ConditionalOStream &pcout,
                    std::string                 header,
                    std::string                 secondary = "")
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


  // Implementation of the class member <code >interpolate_initial_values
  // </code>.

  template <int dim>
  void TimeLoop<dim>::run()
  {
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        pcout << "Reading parameters and allocating objects... " << std::flush;
        ParameterAcceptor::initialize("step-69.prm");
        pcout << "done" << std::endl;
      }
    else
      {
        ParameterAcceptor::initialize("step-69.prm");
      }

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
    auto U = interpolate_initial_values();

    if (resume)
      {
        print_head(pcout, "restore interrupted computation");

        const auto &       triangulation = discretization.triangulation;
        const unsigned int i    = triangulation.locally_owned_subdomain();
        std::string        name = base_name + "-checkpoint-" +
                           dealii::Utilities::int_to_string(i, 4) + ".archive";
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

      } /* end of loop */

    if (output_thread.joinable())
      output_thread.join();

    computing_timer.print_summary();
    pcout << timer_output.str() << std::endl;
  }

  // Implementation of the class member <code >interpolate_initial_values
  // </code>.

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

  // Implementation of the class member <code >TimeLoop </code>.

  template <int dim>
  void TimeLoop<dim>::output(const typename TimeLoop<dim>::vector_type &U,
                             const std::string &                        name,
                             double                                     t,
                             unsigned int                               cycle,
                             bool checkpoint)
  {
    pcout << "TimeLoop<dim>::output(t = " << t
          << ", checkpoint = " << checkpoint << ")" << std::endl;

    if (output_thread.joinable())
      {
        TimerOutput::Scope timer(computing_timer, "time_loop - stalled output");
        output_thread.join();
      }

    constexpr auto problem_dimension =
      ProblemDescription<dim>::problem_dimension;
    const auto &component_names = ProblemDescription<dim>::component_names;

    for (unsigned int i = 0; i < problem_dimension; ++i)
      {
        output_vector[i] = U[i];
        output_vector[i].update_ghost_values();
      }

    schlieren_postprocessor.compute_schlieren(output_vector);

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
                             dealii::Utilities::int_to_string(i, 4) +
                             ".archive";

          if (std::filesystem::exists(name))
            std::filesystem::rename(name, name + "~");

          std::ofstream file(name, std::ios::binary | std::ios::trunc);

          boost::archive::binary_oarchive oa(file);
          oa << t << cycle;
          for (const auto &it1 : output_vector)
            for (const auto &it2 : it1)
              oa << it2;
        }

      dealii::DataOut<dim> data_out;
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

      const auto filename = [&](const unsigned int i) -> std::string {
        const auto seq = dealii::Utilities::int_to_string(i, 4);
        return name + "-" + Utilities::int_to_string(cycle, 6) + "-" + seq +
               ".vtu";
      };

      const unsigned int i = triangulation.locally_owned_subdomain();
      std::ofstream      output(filename(i));
      data_out.write_vtu(output);

      if (dealii::Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        {
          const unsigned int n_mpi_processes =
            dealii::Utilities::MPI::n_mpi_processes(mpi_communicator);
          std::vector<std::string> filenames;
          for (unsigned int i = 0; i < n_mpi_processes; ++i)
            filenames.push_back(filename(i));

          std::ofstream output(name + "-" + Utilities::int_to_string(cycle, 6) +
                               ".pvtu");
          data_out.write_pvtu_record(output, filenames);
        }
    };

    output_thread = std::move(std::thread(output_worker));
  }

} // namespace Step69


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

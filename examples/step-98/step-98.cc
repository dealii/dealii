/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * This program was contributed by Siarhei Uzunbajakau, www.cembooks.nl, 2025.
 */

#include <deal.II/base/tensor.h>
#include <deal.II/base/function.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/types.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_data.h>
#include <deal.II/fe/fe_update_flags.h>

#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_component_interpretation.h>
#include <deal.II/numerics/data_out.h>

#include <string>
#include <ostream>

using namespace dealii;

// @sect3{Settings}

// The following is the control panel of the program. The scaling of the
// program can be changed by setting `mu_0 = 1.0`. Then the computed magnetic
// field, $B$, will have to be multiplied by a factor of
// $\mu_0 = 1.25664 \cdot 10^{-6}$. The free-current density, $\vec{J}_f$, and
// the current vector potential, $T$,  do not depend on scaling. The parameter
// `fe_degree` encodes the degree of the FE_Nedelec, FE_RaviartThomas, and
// FE_DGQ finite elements. The degree of the FE_Q finite elements is computed
// as `fe_degree + 1`. The reason for this is that the lowermost degree of the
// FE_Q finite elements equals 1 while the lowermost degree of the FE_Nedelec,
// FE_RaviartThomas, and FE_DGQ finite elements equals 0 by convention.
// The boundary, material, and manifold IDs are set in the circle.geo file.
// The IDs listed below must match the corresponding IDs set in the circle.geo
// file. If `project_exact_solution = true`, the program projects the exact
// solutions for $T$, $\vec{J}_f$, $\vec{A}$, and $B$ onto the corresponding
// function spaces and saves the results into corresponding vtu files next to
// the numerical solutions. If there are no bugs in the program, the projected
// exact solution and the numerical solution for $T$, $\vec{J}_f$, and $B$ will
// look alike. The situation with $\vec{A}$ is a bit more complicated, see
// section
// [Possibilities for extensions](@ref Step98_PossibilitiesForExtensions). The
// projected exact solution can be used for debugging. In case of problems with
// the convergence of the CG solver the parameter $\eta^2$ must be increased
// just a bit.

namespace Settings
{
  const double permeability_fs = 1.2566370614359172954e-6;
  const double mu_0            = permeability_fs; // Permeability of free space.
  const double mu_r = 4; // Relative permeability of the magnetic material.
  const double mu_1 = mu_0 * mu_r; // Permeability of the magnetic material.
  const double b1   = 0.3;         // Radius of the magnetic core.
  const double a2   = 0.5;         // Inner radius of the free-current region.
  const double b2   = 0.7;         // Outer radius of the free-current region.

  const double K0 = 1.0; // Magnitude of the free-current density.

  /* All IDs are set in the circle.geo file. See the comments there. */
  const types::material_id material_id_free_space =
    1; // Material ID of the free space.
  const types::material_id material_id_core =
    2; // Material ID of the magnetic core.
  const types::material_id material_id_free_current =
    3; // Material ID of the free-current region.
  const types::boundary_id outer_boundary_id = 1; // Boundary ID.
  const types::manifold_id spherical_manifold_id =
    1; /* The ID of the manifold assigned to all of the mesh except for the
        magnetic core. */
  const types::manifold_id flat_manifold_id =
    2; /* The ID of the manifold assigned to the square region in the middle
          of the mesh. */
  const types::manifold_id transfinite_interpolation_manifold_id =
    3; /* The ID of the manifold assigned to the cells in between the square
          region in the middle of the mesh and the innermost circular interface
          (the boundary of the magnetic core). */

  const unsigned int mapping_degree = 2; // Mapping degree used in all solvers.
  const unsigned int fe_degree      = 0; // Degree of the finite elements.

  const double eta_squared = 0.0; // eta^2 when solving for A$.

  const unsigned int n_threads_max = 0; // If >0 limits the number of threads.
  const double eps = 1e-12; // Two doubles are equal if their difference < eps.

  /*The following two parameters control the behavior of the conjugate gradient
   * algorithm*/
  const unsigned int CG_n_mult   = 1; // Defines maximal amount of iterations.
  const double       CG_tol_mult = 1e-8; // Defines the tolerance.

  const bool CG_log_convergence = false; // Save the CG convergence data.
  const bool print_time_tables  = false; // Print the time tables on the screen.
  const bool project_exact_solution = false; // Save the exact solution.
} // namespace Settings

// @sect3{Convergence table}
// This class describes a convergence table. The convergence tables are saved on
// disk in TeX format.
class MainOutputTable : public ConvergenceTable
{
public:
  MainOutputTable() = delete;

  MainOutputTable(const unsigned int dimensions)
    : ConvergenceTable()
    , dimensions(dimensions)
  {}

  void set_new_order(const std::vector<std::string> &new_order_in)
  {
    new_order = new_order_in;
  }

  void append_new_order(const std::string &new_column)
  {
    new_order.push_back(new_column);
  }

  void format()
  {
    set_precision("L2", 2);

    set_scientific("L2", true);

    evaluate_convergence_rates("L2",
                               "ncells",
                               ConvergenceTable::reduction_rate_log2,
                               dimensions);

    set_tex_caption("p", "p");
    set_tex_caption("r", "r");
    set_tex_caption("ncells", "nr. cells");
    set_tex_caption("ndofs", "nr. dofs");
    set_tex_caption("L2", "L2 norm");

    set_column_order(new_order);
  }

  void save(const std::string &fname)
  {
    format();

    std::ofstream ofs(fname + ".tex");
    write_tex(ofs);
  }

private:
  const unsigned int       dimensions;
  std::vector<std::string> new_order = {"p", "r", "ncells", "ndofs", "L2"};
};

// @sect3{Equations}
// This name space contains all closed-form analytical expressions mentioned in
// the introduction to this tutorial.
namespace ExactSolutions
{

  // This function describes the free-current density, $\vec{J}_f$, inside the
  // current region. Two classes implement the current density in this tutorial,
  // `ExactSolutions::FreeCurrentDensity`  and
  // `BaseClasses::FreeCurrentDensity`. There is a subtle difference in how
  // these two classes compute the free-current density. Both classes, however,
  // utilize the same expression for the free-current density. This function
  // describes the expression.
  inline Tensor<1, 2> volume_free_current_density(const Point<2> &p,
                                                  const double    K0)
  {
    return Tensor<1, 2>({-K0 * p[1], K0 * p[0]});
  }

  // This class implements the closed-form analytical expression for the
  // free-current density, $\vec{J}_f$, in the entire domain. The free-current
  // density is computed purely on the basis of the spatial coordinates of the
  // field point. In come coordinates, out comes the free-current density. The
  // information on the material ID of the mesh cells and any other information
  // on the mesh is ignored. This function is used for computing $L_2$ error
  // norms and for computing the projected exact solution. The $\vec{J}_f$ on
  // the right-hand side of the div-grad equation is implemented by the function
  // `BaseSolver::FreeCurrentDensity`.
  class FreeCurrentDensity : public Function<2>
  {
  public:
    FreeCurrentDensity()
      : Function<2>(2)
    {}

    virtual void
    vector_value_list(const std::vector<Point<2>> &p,
                      std::vector<Vector<double>> &values) const override final
    {
      Assert(values.size() == p.size(),
             ExcDimensionMismatch(values.size(), p.size()));

      double r;

      for (unsigned int i = 0; i < values.size(); i++)
        {
          r = p[i].norm();

          if ((r >= Settings::a2) && (r <= Settings::b2))
            {
              const Tensor<1, 2> Jf =
                volume_free_current_density(p[i], Settings::K0);

              values[i][0] = Jf[0];
              values[i][1] = Jf[1];
            }
          else
            {
              values[i][0] = 0.0;
              values[i][1] = 0.0;
            }
        }
    }
  };

  // This class implements the closed-form analytical expression for the
  // magnetic vector potential, $\vec{A}$.
  class MagneticVectorPotential : public Function<2>
  {
  public:
    MagneticVectorPotential()
      : Function<2>(2)
    {}

    virtual void
    vector_value_list(const std::vector<Point<2>> &p,
                      std::vector<Vector<double>> &values) const override final
    {
      Assert(values.size() == p.size(),
             ExcDimensionMismatch(values.size(), p.size()));

      double r;
      double A;

      for (unsigned int i = 0; i < values.size(); i++)
        {
          r = p[i].norm();

          using namespace Settings;

          if (r < b1)
            {
              A = (mu_1 * K0 / 4.0) * (b2 * b2 - a2 * a2);
            }
          else if (r < a2)
            {
              A = (mu_0 * K0 / 4.0) * (b2 * b2 - a2 * a2) *
                  (r + b1 * b1 * (mu_r - 1.0) / r) / r;
            }
          else if (r < b2)
            {
              A = (mu_0 * K0 / 2.0) *
                  (b2 * b2 * r / 2.0 - std::pow(r, 3) / 4.0 +
                   (a2 * (b2 * b2 - a2 * a2) *
                      (a2 + b1 * b1 * (mu_r - 1.0) / a2) / 2.0 -
                    std::pow(b2 * a2, 2) / 2.0 + std::pow(a2, 4) / 4.0) /
                     r) /
                  r;
            }
          else
            {
              A = (mu_0 * K0 / 2.0) * b2 *
                  (std::pow(b2, 3) / 4.0 +
                   (a2 * (b2 * b2 - a2 * a2) *
                      (a2 + b1 * b1 * (mu_r - 1.0) / a2) / 2.0 -
                    std::pow(b2 * a2, 2) / 2.0 + std::pow(a2, 4) / 4.0) /
                     b2) /
                  (r * r);
            }

          values[i][0] = -A * p[i][1];
          values[i][1] = A * p[i][0];
        }
    }
  };

  // This class implements the closed-form analytical expression for the current
  // vector potential, $T$.
  class CurrentVectorPotential : public Function<2>
  {
  public:
    CurrentVectorPotential()
      : Function<2>()
    {}

    virtual void
    value_list(const std::vector<Point<2>> &p,
               std::vector<double>         &values,
               const unsigned int           component = 0) const override final
    {
      Assert(values.size() == p.size(),
             ExcDimensionMismatch(values.size(), p.size()));

      Assert(component == 0,
             ExcMessage("This line is to avoid compiler warnings."));

      double r;

      for (unsigned int i = 0; i < values.size(); i++)
        {
          r = p[i].norm();
          using namespace Settings;

          if (r < a2)
            {
              values[i] = K0 * (b2 * b2 - a2 * a2) / 2.0;
            }
          else if (r < b2)
            {
              values[i] = -K0 * (r * r - b2 * b2) / 2.0;
            }
          else
            {
              values[i] = 0.0;
            }
        }
    }
  };

  // This class implements the closed-form analytical expression for the
  // magnetic field, $B$.
  class MagneticField : public Function<2>
  {
  public:
    MagneticField()
      : Function<2>()
    {}

    virtual void
    value_list(const std::vector<Point<2>> &p,
               std::vector<double>         &values,
               const unsigned int           component = 0) const override final
    {
      Assert(values.size() == p.size(),
             ExcDimensionMismatch(values.size(), p.size()));

      Assert(component == 0,
             ExcMessage("This line is to avoid compiler warnings."));

      double r;

      for (unsigned int i = 0; i < values.size(); i++)
        {
          r = p[i].norm();
          using namespace Settings;

          if (r < b1)
            {
              values[i] = mu_1 * K0 * (b2 * b2 - a2 * a2) / 2.0;
            }
          else if (r < a2)
            {
              values[i] = mu_0 * K0 * (b2 * b2 - a2 * a2) / 2.0;
            }
          else if (r < b2)
            {
              values[i] = mu_0 * K0 * (b2 * b2 - r * r) / 2.0;
            }
          else
            {
              values[i] = 0.0;
            }
        }
    }
  };

} // namespace ExactSolutions

// As discussed above, this name space aggregates the code common to all
// four solvers used in the tutorial. All four solvers are derived from the
// `BaseSolver` class.
namespace BaseClasses
{
  // In general, the `BaseSolver` class needs to know the type of the solver
  // being implemented by the derived class. This information is communicated to
  // `BaseSolver` by passing an argument of the type `SolverType` to the
  // constructor.
  enum SolverType
  {
    ScalarSolver, // Solves a boundary value problem based on div-grad PDE.
    VectorSolver, // Solves a boundary value problem based on curl-curl PDE.
    Projector     // Computes a derivative of a potential.
  };

  // The computation of the integrals of the functionals is delegated to the
  // derived classes. The function that computes the integrals,
  // `BaseSolver::system_matrix_local`, is virtual and must be overridden by
  // the derived classes. However, the objects of the types FEValues are
  // initialized at the level of the `BaseSolver` class together with
  // `AssemblyScratchData`. For the initialization to work properly, the
  // `BaseSolver` class must know which cell data to compute for a particular
  // implementation of the solver down the hierarchy. This information is
  // communicated to the `BaseSolver` by passing an argument of the type
  // `UpdateFlagsCollection` to the constructor.
  struct UpdateFlagsCollection
  {
    UpdateFlags fe_values;
    UpdateFlags fe_values_rhs;
  };

  // The following class implements permeability, $\mu$, in the entire problem
  // domain.
  class Permeability
  {
  public:
    void value_list(const types::material_id mid,
                    std::vector<double>     &values) const
    {
      if (mid == Settings::material_id_core)
        std::fill(values.begin(), values.end(), Settings::mu_1);
      else
        std::fill(values.begin(), values.end(), Settings::mu_0);
    }
  };

  // The following class implements the free-current density, $\vec{J}_f$, on
  // the right-hand side on the div-grad equation when computing the current
  // vector potential, $T$. Unlike `ExactSolutions::FreeCurrentDensity`, this
  // function takes into account the material ID of the cell. If the field point
  // is located outside the cells that constitute the current region,
  // $\vec{J}_f$ is set to zero.
  class FreeCurrentDensity
  {
  public:
    void vector_list(const std::vector<Point<2>> &p,
                     const types::material_id     mid,
                     std::vector<Tensor<1, 2>>   &values) const
    {
      Assert(p.size() == values.size(),
             ExcDimensionMismatch(p.size(), values.size()));

      if ((mid == Settings::material_id_free_space) ||
          (mid == Settings::material_id_core))
        std::fill(values.begin(), values.end(), Tensor<1, 2>());

      if (mid == Settings::material_id_free_current)
        for (unsigned int i = 0; i < values.size(); i++)
          values[i] =
            ExactSolutions::volume_free_current_density(p[i], Settings::K0);
    }
  };

  // The derived classes can control the behavior of the conjugate gradient
  // solver by passing an argument of the type `CgSettings` to the constructor
  // of the `BaseSolver` class.
  struct CgSettings
  {
    unsigned int n_mult;   // Defines maximal amount of iterations.
    double       tol_mult; // Defines the tolerance.
  };

  // All four solvers used in this tutorial are derived from the following
  // class.
  class BaseSolver
  {
  public:
    // As discussed in the introduction, each iteration of the program consists
    // of four stages. At the first stage the mesh is uploaded and refined if
    // necessary. The following constructor is used for implementing the solver
    // of the first stage.
    BaseSolver(const SolverType            solver_type,
               const unsigned int          mapping_degree,
               const UpdateFlagsCollection update_flags_collection,
               const double                eta_squared,
               const std::string          &fname,
               const Function<2>          *exact_solution,
               const CgSettings            cg_settings);

    // The following constructor is used for implementing the solvers of the
    // second, third, and fourth stages. These solvers reuse the mesh loaded at
    // the first stage. They also require a potential computed at one of the
    // preceding stages to be supplied in a form of a field function. To this
    // end, the following constructor contains three extra arguments,
    // `triangulation_rhs`, `dof_handler_rhs`, and `solution_rhs`. The first
    // argument supplies a reference to the triangulation, the second and the
    // third arguments supply the potential computed at one of the preceding
    // stages.
    BaseSolver(const Triangulation<2>     &triangulation_rhs,
               const DoFHandler<2>        &dof_handler_rhs,
               const Vector<double>       &solution_rhs,
               const SolverType            solver_type,
               const unsigned int          mapping_degree,
               const UpdateFlagsCollection update_flags_collection,
               const double                eta_squared,
               const std::string          &fname,
               const Function<2>          *exact_solution,
               const CgSettings            cg_settings);

    // The following eight functions implement various steps of a single
    // simulation. The function `run` aggregates these steps. This arrangement
    // of functions is quite standard in deal.II, see step-3, for instance.
    // Normally, the `setup` function begins by distributing the dofs. The
    // problem is: the type of the finite elements is not known in the
    // `BaseSolver` class. It is specified in the derived classes. Therefore, it
    // could be reasonable to make the setup function virtual as well. Instead,
    // we move the dof distribution code to the end of the `make_mesh` function
    // which is, in fact, virtual and must be overridden in the derived classes
    // anyway.
    virtual void make_mesh() = 0; /* If used in the first stage, loads the mesh
                                     and distributes the dofs. In other stages
                                     just distributes the dofs.*/
    void setup();      /* Applies Dirichlet boundary condition, setups the
                          dof pattern, and initializes vectors, matrices. */
    void assemble();   // Assembles the system of linear equations.
    void solve();      // Solves the system of linear equations.
    void save() const; // Saves the result into a vtu file.
    void compute_error_norms();        // Computes L^2 error norm.
    void project_exact_solution_fcn(); // Projects exact solution.
    void clear();                      // Clears the memory for the next solver.
    void run(); /* Executes the last eight functions in the proper order
                   and measures the execution time for each function. */

    // The following three functions provide access to protected data. These
    // data is used to fill the convergence tables in the main loop below.
    double                  get_L2_norm();
    unsigned int            get_n_cells() const;
    types::global_dof_index get_n_dofs() const;

    // The following three functions provide references to the triangulation,
    // dof handler, and the dofs that describe the solution. These data is fed
    // to one of the next solvers.
    const Triangulation<2> &get_tria() const;
    const DoFHandler<2>    &get_dof_handler() const;
    const Vector<double>   &get_solution() const;

  protected:
    // The following three data members describe the solver. The first data
    // member is supplied as an argument to the constructor. The other two are
    // deduced in  the `setup` function.
    const SolverType solver_type;

    bool fe_is_vector;
    bool fe_is_FE_Q;

    // The following data members store the input from one of the preceding
    // solvers. If there is no preceding solver and these data is not provided,
    // i.e., first of the two constructors above was used, the three data
    // members below point to `triangulation`, `dof_handler`, and `solution`,
    // respectively. See the implementation of the first constructor below.
    const Triangulation<2> &triangulation_rhs;
    const DoFHandler<2>    &dof_handler_rhs;
    const Vector<double>   &solution_rhs;

    // If the derived class implements the solver of the first stage, the
    // following data member stores the uploaded and refined mesh. This data
    // member is not used if the derived class implements solvers of the second,
    // third, and fourth stages.
    Triangulation<2> triangulation;

    // The following two data members are the dof handler and the dofs that
    // describe the solution.
    DoFHandler<2>  dof_handler;
    Vector<double> solution;

    // These two data members describe the system of linear equations to be
    // solved by the solver. The components of the system matrix and that of the
    // right-hand side are computed by the `assemble` function.
    SparseMatrix<double> system_matrix;
    Vector<double>       system_rhs;

    // The data member `exact_solution` points to the closed-form analytical
    // solution the solver attempts to compute. If
    // `Settings::project_exact_solution=true`, the exact solution
    // is projected onto a proper function space and the data member
    // `projected_exact_solution` is populated by the dofs of the projected
    // exact solution. The corresponding dof handler is `dof_handler`. Together
    // `dof_handler` and `projected_exact_solution` constitute the field
    // function which describes the exact solution. It is saved into the vtu
    // file next to the solution and the $L_2$ error norm.
    const Function<2> *exact_solution;
    Vector<double>     projected_exact_solution;

    // The constraints are used to apply the Dirichlet boundary conditions and
    // distribute the local (cell specific) system matrix and right-hand side to
    // their global counterparts. The are no hanging nodes and hanging node
    // constraints in this program.
    AffineConstraints<double> constraints;

    // The following stores the dynamic sparsity pattern. The tutorial step-2
    // discusses the rationale behind the dynamic sparsity pattern.
    SparsityPattern sparsity_pattern;

    // The following five data members simply store the data supplied as
    // arguments to the constructor.
    const unsigned int          mapping_degree;
    const UpdateFlagsCollection update_flags_collection;
    const double                eta_squared;
    const std::string           fname;
    const CgSettings            cg_settings;

    // The homogeneous Dirichlet boundary condition to be used with the static
    // solver when solving for $T$.
    const Functions::ZeroFunction<2> dirichlet_boundary_condition;

    // The following two data members are computed by the function
    // `compute_error_norms`. `L2_per_cell` is saved into the vtu file next to
    // the solution. `L2_norm` is reported in the convergence table.
    Vector<double> L2_per_cell;
    double         L2_norm;

    // The timer is used in the `run` function to measure the execution time of
    // various sections of the code. The execution time is reported in the time
    // table. To print the time table on the screen, set
    // `Settings::print_time_table = true`.
    TimerOutput timer;

    // The solvers at the second, third, and fourth stages use two dof handlers,
    // `dof_handler` and `dof_handler_rhs`. The WorkStream needs to walk through
    // the two dof handlers synchronously. For this purpose we will pair two
    // active cell iterators (one from `dof_handler`, another from
    // `dof_handler_rhs`). For that we need the `IteratorPair` type. The solver
    // at the first stage (that is, a  solver constructed by invoking the first
    // constructor above) uses only one dof handler. In this case, however, the
    // constructor makes `dof_handler_rhs` to reference `dof_handler`. In
    // effect, both iterators of the tuple will iterate the same dof handler. In
    // the case of the first-stage solver we will use only the first iterator.
    using IteratorTuple =
      std::tuple<typename DoFHandler<2>::active_cell_iterator,
                 typename DoFHandler<2>::active_cell_iterator>;

    using IteratorPair = SynchronousIterators<IteratorTuple>;

    // The program utilizes the WorkStream technology. The Step-9 tutorial
    // does a much better job of explaining the workings of WorkStream.
    // Reading the "WorkStream paper", see the glossary, is recommended.
    // The integrals of the functionals are computed by the function
    // `system_matrix_local`. This function is specific to every solver and must
    // be overridden in the derived classes.
    struct AssemblyScratchData
    {
      AssemblyScratchData(const DoFHandler<2>        &dof_handler,
                          const DoFHandler<2>        &dof_handler_rhs,
                          const Vector<double>       &dofs_rhs,
                          const unsigned int          mapping_degree,
                          const UpdateFlagsCollection update_flags_collection,
                          const bool                  fe_is_FE_Q,
                          const double                eta_squared);

      AssemblyScratchData(const AssemblyScratchData &scratch_data);

      const Permeability permeability;

      const FreeCurrentDensity free_current_density;

      MappingQ<2> mapping;

      FEValues<2> fe_values;
      FEValues<2> fe_values_rhs;

      const unsigned int dofs_per_cell;
      const unsigned int n_q_points;

      std::vector<double>                    permeability_list;
      std::vector<double>                    values_list_rhs;
      std::vector<Tensor<1, 2>>              vectors_list_rhs;
      std::vector<std::vector<Tensor<1, 2>>> vectors_vectors_list_rhs;

      const FEValuesExtractors::Scalar se;
      const FEValuesExtractors::Vector ve;

      const DoFHandler<2>  &dof_handler_rhs;
      const Vector<double> &dofs_rhs;

      const double eta_squared;
    };

    struct AssemblyCopyData
    {
      FullMatrix<double>                   cell_matrix;
      Vector<double>                       cell_rhs;
      std::vector<types::global_dof_index> local_dof_indices;
    };

    virtual void system_matrix_local(const IteratorPair  &IP,
                                     AssemblyScratchData &scratch_data,
                                     AssemblyCopyData    &copy_data) = 0;

    void copy_local_to_global(const AssemblyCopyData &copy_data);
  }; // class BaseSolver

  // The implementation of the first constructor. It is used only for
  // constructing the first-stage solver. Note that `triangulation_rhs`,
  // `dof_handler_rhs`, and `solution_rhs` reference `triangulation`,
  // `dof_handler`, and `solution`, respectively.
  BaseSolver::BaseSolver(const SolverType            solver_type,
                         const unsigned int          mapping_degree,
                         const UpdateFlagsCollection update_flags_collection,
                         const double                eta_squared,
                         const std::string          &fname,
                         const Function<2>          *exact_solution,
                         const CgSettings            cg_settings)
    : solver_type(solver_type)
    , triangulation_rhs(triangulation)
    , dof_handler_rhs(dof_handler)
    , solution_rhs(solution)
    , exact_solution(exact_solution)
    , mapping_degree(mapping_degree)
    , update_flags_collection(update_flags_collection)
    , eta_squared(eta_squared)
    , fname(fname)
    , cg_settings(cg_settings)
    , timer(std::cout,
            (Settings::print_time_tables) ? TimerOutput::summary :
                                            TimerOutput::never,
            TimerOutput::cpu_and_wall_times_grouped)
  {}

  // The implementation of the second constructor. It is used for constructing
  // solvers at the second, third, and the fourth stages.
  BaseSolver::BaseSolver(const Triangulation<2>     &triangulation_rhs,
                         const DoFHandler<2>        &dof_handler_rhs,
                         const Vector<double>       &solution_rhs,
                         const SolverType            solver_type,
                         const unsigned int          mapping_degree,
                         const UpdateFlagsCollection update_flags_collection,
                         const double                eta_squared,
                         const std::string          &fname,
                         const Function<2>          *exact_solution,
                         const CgSettings            cg_settings)
    : solver_type(solver_type)
    , triangulation_rhs(triangulation_rhs)
    , dof_handler_rhs(dof_handler_rhs)
    , solution_rhs(solution_rhs)
    , exact_solution(exact_solution)
    , mapping_degree(mapping_degree)
    , update_flags_collection(update_flags_collection)
    , eta_squared(eta_squared)
    , fname(fname)
    , cg_settings(cg_settings)
    , timer(std::cout,
            (Settings::print_time_tables) ? TimerOutput::summary :
                                            TimerOutput::never,
            TimerOutput::cpu_and_wall_times_grouped)
  {}

  // This function applies the Dirichlet boundary condition, setups a sparsity
  // pattern, and initializes the vectors and matrices. It is common for the
  // setup function to distribute the dofs. The type of the finite elements,
  // however, is not known at the level of `BaseSolver` class. The type of the
  // finite elements is chosen in the derived classes. For this reason, the
  // task of distributing the dofs is shifted to the end of the `make_mesh`
  // function which is a virtual function.
  void BaseSolver::setup()
  {
    // First we check if the type of the finite elements is chosen properly.
    Assert(
      (dof_handler.get_fe(0).get_name().compare(0, 4, "FE_Q") == 0) ||
        (dof_handler.get_fe(0).get_name().compare(0, 6, "FE_DGQ") == 0) ||
        (dof_handler.get_fe(0).get_name().compare(0, 10, "FE_Nedelec") == 0) ||
        (dof_handler.get_fe(0).get_name().compare(0, 16, "FE_RaviartThomas") ==
         0),
      ExcMessage(
        "The program assumes one of the following types of the finite elements:\
       FE_Q, FE_Nedelec, FE_RaviartThomas, FE_DGQ. These finite elements\
       discretize the four function spaces of the De Rham complex. Other\
       finite elements are not allowed."));

    if (solver_type == ScalarSolver)
      Assert((dof_handler.get_fe(0).get_name().compare(0, 4, "FE_Q") == 0),
             ExcMessage("The scalar solver (div-grad) must use FE_Q finite\
         elements."));

    if (solver_type == VectorSolver)
      Assert(
        (dof_handler.get_fe(0).get_name().compare(0, 10, "FE_Nedelec") == 0),
        ExcMessage(
          "The vector solver (curl-curl) must use FE_Nedelec finite elements."));

    if (solver_type == Projector)
      Assert(!(dof_handler.get_fe(0).get_name().compare(0, 4, "FE_Q") == 0),
             ExcMessage(
               "Projector cannot use FE_Q finite elements. It is impossible to\
               project a derivative of a potential onto H(grad) function space,\
               see the Bossavit's diagrams or De Rham complex."));

    // Second, we deduce if the finite elements are vector elements and if the
    // finite elements are of the FE_Q type.
    fe_is_vector =
      ((dof_handler.get_fe(0).get_name().compare(0, 10, "FE_Nedelec") == 0) ||
       (dof_handler.get_fe(0).get_name().compare(0, 16, "FE_RaviartThomas") ==
        0));

    fe_is_FE_Q = (dof_handler.get_fe(0).get_name().compare(0, 4, "FE_Q") == 0);

    // Third, we apply the Dirichlet boundary condition. As discussed in
    // step-97, the Dirichlet boundary condition is an essential condition in
    // the case of the curl-curl partial differential equation. The same method
    // (interrogating the functional) can be used to show that the Dirichlet
    // boundary condition is essential in the case of the div-grad partial
    // differential equation as well. Essential boundary conditions must be
    // enforced by constraining the system matrix. This segment of the code does
    // the constraining. In this tutorial the Dirichlet boundary condition is
    // used only in the case of the static solver when solving for $T$. In all
    // other cases we leave the constraints empty.
    constraints.clear();

    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    if (solver_type == ScalarSolver)
      VectorTools::interpolate_boundary_values(MappingQ<2>(mapping_degree),
                                               dof_handler,
                                               Settings::outer_boundary_id,
                                               dirichlet_boundary_condition,
                                               constraints);

    constraints.close();

    // Fourth, we setup the sparsity pattern.
    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    sparsity_pattern.copy_from(dsp);

    // Fifth, we initialize the matrices and arrays.
    system_matrix.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    if (Settings::project_exact_solution && exact_solution)
      projected_exact_solution.reinit(dof_handler.n_dofs());

    if (exact_solution)
      L2_per_cell.reinit(triangulation.n_active_cells());
  }

  // Formally, the following function assembles the system of linear equations.
  // In reality, however, it just spells all the magic words to get the
  // WorkStream going. The interesting part, i.e., computing the components of
  // the system matrix and the right-hand side, happens in the function
  // `system_matrix_local`. The last is a virtual function. That is to say, the
  // recipe for the functional is implemented in the derived class by overriding
  // the virtual function `system_matrix_local`.
  void BaseSolver::assemble()
  {
    WorkStream::run(IteratorPair(IteratorTuple(dof_handler.begin_active(),
                                               dof_handler_rhs.begin_active())),
                    IteratorPair(
                      IteratorTuple(dof_handler.end(), dof_handler_rhs.end())),
                    *this,
                    &BaseSolver::system_matrix_local,
                    &BaseSolver::copy_local_to_global,
                    AssemblyScratchData(dof_handler,
                                        dof_handler_rhs,
                                        solution_rhs,
                                        mapping_degree,
                                        update_flags_collection,
                                        fe_is_FE_Q,
                                        eta_squared),
                    AssemblyCopyData());
  }

  // The following two constructors initialize scratch data from the input
  // parameters and from another object of the same type, i.e., a copy
  // constructor.
  BaseSolver::AssemblyScratchData::AssemblyScratchData(
    const DoFHandler<2>        &dof_handler,
    const DoFHandler<2>        &dof_handler_rhs,
    const Vector<double>       &dofs_rhs,
    const unsigned int          mapping_degree,
    const UpdateFlagsCollection update_flags_collection,
    const bool                  fe_is_FE_Q,
    const double                eta_squared)
    : permeability()
    , free_current_density()
    , mapping(mapping_degree)
    , fe_values(mapping,
                dof_handler.get_fe(),
                QGauss<2>((fe_is_FE_Q) ? (dof_handler.get_fe().degree + 1) :
                                         (dof_handler.get_fe().degree + 2)),
                update_flags_collection.fe_values)
    , fe_values_rhs(mapping,
                    dof_handler_rhs.get_fe(),
                    QGauss<2>((fe_is_FE_Q) ? (dof_handler.get_fe().degree + 1) :
                                             (dof_handler.get_fe().degree + 2)),
                    update_flags_collection.fe_values_rhs)
    , dofs_per_cell(fe_values.dofs_per_cell)
    , n_q_points(fe_values.get_quadrature().size())
    , permeability_list(n_q_points)
    , values_list_rhs(n_q_points)
    , vectors_list_rhs(n_q_points)
    , vectors_vectors_list_rhs(n_q_points, std::vector<Tensor<1, 2>>(2))
    , se(0)
    , ve(0)
    , dof_handler_rhs(dof_handler_rhs)
    , dofs_rhs(dofs_rhs)
    , eta_squared(eta_squared)
  {}

  BaseSolver::AssemblyScratchData::AssemblyScratchData(
    const AssemblyScratchData &scratch_data)
    : permeability()
    , free_current_density()
    , mapping(scratch_data.mapping.get_degree())
    , fe_values(mapping,
                scratch_data.fe_values.get_fe(),
                scratch_data.fe_values.get_quadrature(),
                scratch_data.fe_values.get_update_flags())
    , fe_values_rhs(mapping,
                    scratch_data.fe_values_rhs.get_fe(),
                    scratch_data.fe_values_rhs.get_quadrature(),
                    scratch_data.fe_values_rhs.get_update_flags())
    , dofs_per_cell(fe_values.dofs_per_cell)
    , n_q_points(fe_values.get_quadrature().size())
    , permeability_list(n_q_points)
    , values_list_rhs(n_q_points)
    , vectors_list_rhs(n_q_points)
    , vectors_vectors_list_rhs(n_q_points, std::vector<Tensor<1, 2>>(2))
    , se(0)
    , ve(0)
    , dof_handler_rhs(scratch_data.dof_handler_rhs)
    , dofs_rhs(scratch_data.dofs_rhs)
    , eta_squared(scratch_data.eta_squared)
  {}

  // This function copies the components of a cell matrix and a cell right-hand
  // side into the system matrix, $A_{ij}$, and the system right-hand side,
  // $b_i$.
  void BaseSolver::copy_local_to_global(const AssemblyCopyData &copy_data)
  {
    constraints.distribute_local_to_global(copy_data.cell_matrix,
                                           copy_data.cell_rhs,
                                           copy_data.local_dof_indices,
                                           system_matrix,
                                           system_rhs);
  }

  // This function solves the system of linear equations.
  // If `Settings::log_cg_convergence == true`, the convergence data
  // is saved into a file.  The maximal number of iteration steps is defined as
  // `cg_settings.n_mult * system_rhs.size()`. The stopping condition is
  // $|\boldsymbol{b}-\boldsymbol{A}\boldsymbol{c}|<\alpha|\boldsymbol{b}|$,
  // where $\alpha$ is defined by `cg_settings.tol_mult`. As soon as we use
  // constraints, we must not forget to distribute them.
  void BaseSolver::solve()
  {
    SolverControl control(cg_settings.n_mult * system_rhs.size(),
                          cg_settings.tol_mult * system_rhs.l2_norm(),
                          false,
                          false);

    if (Settings::CG_log_convergence)
      control.enable_history_data();

    GrowingVectorMemory<Vector<double>> memory;
    SolverCG<Vector<double>>            cg(control, memory);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    constraints.distribute(solution);

    if (Settings::CG_log_convergence)
      {
        const std::vector<double> &history_data = control.get_history_data();

        std::ofstream ofs(fname + "_cg_convergence.csv");

        for (unsigned int i = 1; i < history_data.size(); i++)
          ofs << i << ", " << history_data[i] << "\n";
      }
  }

  // The following two functions compute the error norms and project the exact
  // solution.
  void BaseSolver::compute_error_norms()
  {
    if (exact_solution)
      {
        VectorTools::integrate_difference(MappingQ<2>(mapping_degree),
                                          dof_handler,
                                          solution,
                                          *exact_solution,
                                          L2_per_cell,
                                          QGauss<2>(
                                            dof_handler.get_fe(0).degree + 4),
                                          VectorTools::L2_norm);

        L2_norm = VectorTools::compute_global_error(triangulation_rhs,
                                                    L2_per_cell,
                                                    VectorTools::L2_norm);
      }
  }

  void BaseSolver::project_exact_solution_fcn()
  {
    if (Settings::project_exact_solution && exact_solution)
      {
        AffineConstraints<double> constraints_empty;
        constraints_empty.close();

        VectorTools::project(MappingQ<2>(mapping_degree),
                             dof_handler,
                             constraints_empty,
                             QGauss<2>((fe_is_FE_Q) ?
                                         (dof_handler.get_fe().degree + 1) :
                                         (dof_handler.get_fe().degree + 2)),
                             *exact_solution,
                             projected_exact_solution);
      }
  }

  // This function saves the solution and the $L_2$ error norm into a vtu file.
  // If `Settings::project_exact_solution = true`, the projected
  // exact solution is saved as well.
  void BaseSolver::save() const
  {
    const std::string  name = (fe_is_vector) ? "VectorField" : "ScalarField";
    const unsigned int components = (fe_is_vector) ? 2 : 1;
    const DataComponentInterpretation::DataComponentInterpretation
      component_interpretation =
        (fe_is_vector) ?
          DataComponentInterpretation::component_is_part_of_vector :
          DataComponentInterpretation::component_is_scalar;

    std::vector<std::string> solution_names(components, name);
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(components, component_interpretation);

    DataOut<2> data_out;

    data_out.add_data_vector(dof_handler,
                             solution,
                             solution_names,
                             data_component_interpretation);

    if (exact_solution)
      {
        data_out.add_data_vector(L2_per_cell, "L2norm");

        if (Settings::project_exact_solution)
          {
            std::vector<std::string> solution_names_exact(components,
                                                          name + "Exact");
            data_out.add_data_vector(dof_handler,
                                     projected_exact_solution,
                                     solution_names_exact,
                                     data_component_interpretation);
          }
      }

    DataOutBase::VtkFlags flags;
    flags.write_higher_order_cells = true;
    data_out.set_flags(flags);

    const MappingQ<2> mapping(mapping_degree);

    data_out.build_patches(mapping,
                           dof_handler.get_fe(0).degree + 2,
                           DataOut<2>::CurvedCellRegion::curved_inner_cells);

    std::ofstream ofs(fname + ".vtu");
    data_out.write_vtu(ofs);
  }

  // This function clears the memory for the next solver.
  void BaseSolver::clear()
  {
    system_matrix.clear();
    system_rhs.reinit(0);
  }

  // This functions calls all the constituent functions in the right order and
  // measures the execution time.
  void BaseSolver::run()
  {
    {
      TimerOutput::Scope timer_section(timer, "Make mesh");
      make_mesh();
    }
    {
      TimerOutput::Scope timer_section(timer, "Setup");
      setup();
    }
    {
      TimerOutput::Scope timer_section(timer, "Assemble");
      assemble();
    }
    {
      TimerOutput::Scope timer_section(timer, "Solve");
      solve();
    }
    {
      TimerOutput::Scope timer_section(timer, "Compute error norms");
      compute_error_norms();
    }
    {
      TimerOutput::Scope timer_section(timer, "Project exact solution");
      project_exact_solution_fcn();
    }
    {
      TimerOutput::Scope timer_section(timer, "Save");
      save();
    }
    {
      TimerOutput::Scope timer_section(timer, "Clear");
      clear();
    }
  }

  // The following is the straightforward implementation of the get functions.
  double BaseSolver::get_L2_norm()
  {
    return L2_norm;
  }

  unsigned int BaseSolver::get_n_cells() const
  {
    return triangulation_rhs.n_active_cells();
  }

  types::global_dof_index BaseSolver::get_n_dofs() const
  {
    return dof_handler.n_dofs();
  }

  const Triangulation<2> &BaseSolver::get_tria() const
  {
    return triangulation_rhs;
  }

  const DoFHandler<2> &BaseSolver::get_dof_handler() const
  {
    return dof_handler;
  }

  const Vector<double> &BaseSolver::get_solution() const
  {
    return solution;
  }

} // namespace BaseClasses

// This variable facilitates derivation of the solver classes. If a derived
// solver class can be constructed using the global settings, just pass
// `cg_settings_global` as the last parameter to the constructor. If a
// solver requires tailored settings, construct a new structure of the type
// `BaseClasses::CgSettings`.
BaseClasses::CgSettings cg_settings_global = {Settings::CG_n_mult,
                                              Settings::CG_tol_mult};

// @sect3{Solver - T}

// This name space contains the code related to the computation of the current
// vector potential, $T$.
namespace SolverT
{
  using namespace BaseClasses;

  // We derive the solver from the `BaseSolver` class. What is left to do is
  // to initialize the `BaseSolver` properly and override two virtual functions,
  // `make_mesh` and `system_matrix_local`.
  class Solver : public BaseSolver
  {
  public:
    Solver() = delete;
    Solver(const unsigned int p, // Degree of the FE_Q finite elements.
           const unsigned int r, // The mesh refinement parameter.
           const unsigned int mapping_degree,
           const std::string &fname          = "data",
           const Function<2> *exact_solution = nullptr);

    ~Solver() = default;

  private:
    virtual void make_mesh() override final;
    virtual void
    system_matrix_local(const IteratorPair  &IP,
                        AssemblyScratchData &scratch_data,
                        AssemblyCopyData    &copy_data) override final;

    const FE_Q<2>      fe;
    const unsigned int refinement_parameter;
  };

  // This solver is used at the first stage, see the introduction. For this
  // reason, it needs to upload the mesh. Consequently, we use the first
  // constructor of the `BaseSolver` class.
  Solver::Solver(const unsigned int p,
                 const unsigned int r,
                 const unsigned int mapping_degree,
                 const std::string &fname,
                 const Function<2> *exact_solution)
    : BaseSolver(SolverType::ScalarSolver,
                 mapping_degree,
                 {update_quadrature_points | update_gradients |
                    update_JxW_values,
                  update_default},
                 0.0,
                 fname,
                 exact_solution,
                 cg_settings_global)

    , fe(p)
    , refinement_parameter(r)
  {}

  // This function uploads the mesh, creates manifolds, bounds the manifolds
  // to the manifold IDs, refines the mesh, and distributes the dofs. Note that
  // we are allowed to create the manifolds locally in the function as the
  // triangulation object keeps copies of the manifolds, see step-65. Also
  // recall that we have shifted the task of distributing the dofs from the
  // `setup` function to the `make_mesh` function to evade the necessity to make
  // `setup` function virtual. This allows us to keep one `setup` function in
  // the `BaseSolver` class that serves the needs of all derived classes.
  void Solver::make_mesh()
  {
    GridIn<2> gridin;

    gridin.attach_triangulation(triangulation);
    gridin.read_msh("circle.msh");

    using namespace Settings;

    triangulation.set_manifold(spherical_manifold_id, SphericalManifold<2>());
    triangulation.set_manifold(flat_manifold_id, FlatManifold<2>());

    TransfiniteInterpolationManifold<2> transfinite_manifold;
    transfinite_manifold.initialize(triangulation);
    triangulation.set_manifold(transfinite_interpolation_manifold_id,
                               transfinite_manifold);

    triangulation.refine_global(refinement_parameter);

    dof_handler.reinit(triangulation);
    dof_handler.distribute_dofs(fe);
  }

  // This function assembles a fraction of the system matrix and the system
  // right-hand side related to a single cell. These fractions are
  // `copy_data.cell_matrix` and `copy_data.cell_rhs`. They are copied into
  // the system matrix, $A_{ij}$, and the right-hand side, $b_i$, by the
  // function `Solver::copy_local_to_global()`.
  void Solver::system_matrix_local(const IteratorPair  &IP,
                                   AssemblyScratchData &scratch_data,
                                   AssemblyCopyData    &copy_data)
  {
    copy_data.cell_matrix.reinit(scratch_data.dofs_per_cell,
                                 scratch_data.dofs_per_cell);

    copy_data.cell_rhs.reinit(scratch_data.dofs_per_cell);

    copy_data.local_dof_indices.resize(scratch_data.dofs_per_cell);

    auto cell = std::get<0>(*IP);

    scratch_data.fe_values.reinit(cell);

    scratch_data.free_current_density.vector_list(
      scratch_data.fe_values.get_quadrature_points(),
      cell->material_id(),
      scratch_data.vectors_list_rhs);

    for (unsigned int q_index = 0; q_index < scratch_data.n_q_points; ++q_index)
      {
        for (unsigned int i = 0; i < scratch_data.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < scratch_data.dofs_per_cell; ++j)
              {
                copy_data.cell_matrix(i, j) += // Integral I_a1.
                  (scratch_data.fe_values[scratch_data.se].gradient(
                     i, q_index) * // grad N_i
                   scratch_data.fe_values[scratch_data.se].gradient(
                     j, q_index) // grad N_j
                   ) *
                  scratch_data.fe_values.JxW(q_index); // dS
              }
            copy_data.cell_rhs(i) += // Integral I_b3-1.
              (scratch_data.vectors_list_rhs[q_index][0] *
                 scratch_data.fe_values[scratch_data.se].gradient(i,
                                                                  q_index)[1] -
               scratch_data.vectors_list_rhs[q_index][1] *
                 scratch_data.fe_values[scratch_data.se].gradient(i,
                                                                  q_index)[0]) *
              scratch_data.fe_values.JxW(q_index); // J_f.(curlv N_i)dS.
          }
      }

    cell->get_dof_indices(copy_data.local_dof_indices);
  }
} // namespace SolverT

// @sect3{Solver - A}

// This name space contains all the code related to the computation of the
// magnetic vector potential, $\vec{A}$.
namespace SolverA
{
  using namespace BaseClasses;

  // We derive the solver from the `BaseSolver` class. What is left to do is
  // to initialize the `BaseSolver` properly and override two virtual functions,
  // `make_mesh` and `system_matrix_local`.
  class Solver : public BaseSolver
  {
  public:
    Solver() = delete;
    Solver(const unsigned int p, // Degree of the FE_Nedelec finite elements.
           const unsigned int mapping_degree,
           const Triangulation<2> &triangulation_ext,
           const DoFHandler<2>    &dof_handler_ext,
           const Vector<double>   &solution_ext,
           const double            eta_squared    = 0.0,
           const std::string      &fname          = "data",
           const Function<2>      *exact_solution = nullptr);

    ~Solver() = default;

  private:
    virtual void make_mesh() override final;
    virtual void
    system_matrix_local(const IteratorPair  &IP,
                        AssemblyScratchData &scratch_data,
                        AssemblyCopyData    &copy_data) override final;

    FE_Nedelec<2> fe;
  };

  // This solver is used at the third stage, see the introduction. Consequently,
  // we use the second constructor of the `BaseSolver` class.
  Solver::Solver(const unsigned int      p,
                 const unsigned int      mapping_degree,
                 const Triangulation<2> &triangulation_rhs,
                 const DoFHandler<2>    &dof_handler_rhs,
                 const Vector<double>   &solution_rhs,
                 const double            eta_squared,
                 const std::string      &fname,
                 const Function<2>      *exact_solution)
    : BaseSolver(triangulation_rhs,
                 dof_handler_rhs,
                 solution_rhs,
                 SolverType::VectorSolver,
                 mapping_degree,
                 {update_values | update_gradients | update_JxW_values,
                  update_values},
                 eta_squared,
                 fname,
                 exact_solution,
                 {100, Settings::CG_tol_mult})
    , fe(p)
  {}

  // At the third stage we do not upload the mesh. We reuse the mesh uploaded at
  // the first stage. Consequently, we just need to distribute the dofs.
  void Solver::make_mesh()
  {
    dof_handler.reinit(triangulation_rhs);
    dof_handler.distribute_dofs(fe);
  }

  // This function assembles a fraction of the system matrix and the system
  // right-hand side related to a single cell. These fractions are
  // `copy_data.cell_matrix` and `copy_data.cell_rhs`. They are copied into
  // the system matrix, $A_{ij}$, and the right-hand side, $b_i$, by the
  // function `Solver::copy_local_to_global()`.
  void Solver::system_matrix_local(const IteratorPair  &IP,
                                   AssemblyScratchData &scratch_data,
                                   AssemblyCopyData    &copy_data)
  {
    copy_data.cell_matrix.reinit(scratch_data.dofs_per_cell,
                                 scratch_data.dofs_per_cell);

    copy_data.cell_rhs.reinit(scratch_data.dofs_per_cell);

    copy_data.local_dof_indices.resize(scratch_data.dofs_per_cell);

    auto cell     = std::get<0>(*IP);
    auto cell_ext = std::get<1>(*IP);

    scratch_data.fe_values.reinit(cell);
    scratch_data.fe_values_rhs.reinit(cell_ext);

    scratch_data.permeability.value_list(cell->material_id(),
                                         scratch_data.permeability_list);

    scratch_data.fe_values_rhs.get_function_values(
      scratch_data.dofs_rhs, scratch_data.values_list_rhs);

    for (unsigned int q_index = 0; q_index < scratch_data.n_q_points; ++q_index)
      {
        for (unsigned int i = 0; i < scratch_data.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < scratch_data.dofs_per_cell; ++j)
              {
                copy_data.cell_matrix(i, j) += // Integral I_a1+I_a3.
                  (1.0 / scratch_data.permeability_list[q_index]) * // 1 / mu
                  (scratch_data.fe_values[scratch_data.ve].curl(
                     i, q_index) * // curls N_i
                     scratch_data.fe_values[scratch_data.ve].curl(
                       j, q_index)              // curls N_j
                   + scratch_data.eta_squared * // eta^2
                       scratch_data.fe_values[scratch_data.ve].value(
                         i, q_index) * // N_i
                       scratch_data.fe_values[scratch_data.ve].value(
                         j, q_index) // N_j
                   ) *
                  scratch_data.fe_values.JxW(q_index); // dS
              }
            copy_data.cell_rhs(i) += // Integral I_b3-1.
              (scratch_data.values_list_rhs[q_index] *
               scratch_data.fe_values[scratch_data.ve].curl(i, q_index))[0] *
              scratch_data.fe_values.JxW(q_index); // T(curls N_i)dS
          }
      }

    cell->get_dof_indices(copy_data.local_dof_indices);
  }
} // namespace SolverA

// @sect3{Projector from H(grad) to H(div)}

// This name space contains all the code related to the computation of the
// free-current density, $\vec{J}_f$.
namespace ProjectorHgradToHdiv
{
  using namespace BaseClasses;

  // We derive the solver from the `BaseSolver` class. What is left to do is
  // to initialize the `BaseSolver` properly and override two virtual functions,
  // `make_mesh` and `system_matrix_local`.
  class Solver : public BaseSolver
  {
  public:
    Solver() = delete;
    Solver(const unsigned int p, /* Degree of the FE_RaviartThomas finite
                                    elements. */
           const unsigned int      mapping_degree,
           const Triangulation<2> &triangulation_rhs,
           const DoFHandler<2>    &dof_handler_rhs,
           const Vector<double>   &solution_rhs,
           const std::string      &fname          = "data",
           const Function<2>      *exact_solution = nullptr);

    ~Solver() = default;

  private:
    virtual void make_mesh() override final;
    virtual void
    system_matrix_local(const IteratorPair  &IP,
                        AssemblyScratchData &scratch_data,
                        AssemblyCopyData    &copy_data) override final;

    FE_RaviartThomas<2> fe;
  };

  // This solver is used at the second stage, see the introduction.
  // Consequently, we use the second constructor of the `BaseSolver` class.
  Solver::Solver(const unsigned int      p,
                 const unsigned int      mapping_degree,
                 const Triangulation<2> &triangulation_rhs,
                 const DoFHandler<2>    &dof_handler_rhs,
                 const Vector<double>   &solution_rhs,
                 const std::string      &fname,
                 const Function<2>      *exact_solution)
    : BaseSolver(triangulation_rhs,
                 dof_handler_rhs,
                 solution_rhs,
                 SolverType::Projector,
                 mapping_degree,
                 {update_values | update_JxW_values, update_gradients},
                 0.0,
                 fname,
                 exact_solution,
                 cg_settings_global)
    , fe(p)
  {}

  // At the second stage we do not upload the mesh. We reuse the mesh uploaded
  // at the first stage. Consequently, we just need to distribute the dofs.
  void Solver::make_mesh()
  {
    dof_handler.reinit(triangulation_rhs);
    dof_handler.distribute_dofs(fe);
  }

  // This function assembles a fraction of the system matrix and the system
  // right-hand side related to a single cell. These fractions are
  // `copy_data.cell_matrix` and `copy_data.cell_rhs`. They are copied into
  // the system matrix, $A_{ij}$, and the right-hand side, $b_i$, by the
  // function `Solver::copy_local_to_global()`.
  void Solver::system_matrix_local(const IteratorPair  &IP,
                                   AssemblyScratchData &scratch_data,
                                   AssemblyCopyData    &copy_data)
  {
    copy_data.cell_matrix.reinit(scratch_data.dofs_per_cell,
                                 scratch_data.dofs_per_cell);

    copy_data.cell_rhs.reinit(scratch_data.dofs_per_cell);

    copy_data.local_dof_indices.resize(scratch_data.dofs_per_cell);

    scratch_data.fe_values.reinit(std::get<0>(*IP));
    scratch_data.fe_values_rhs.reinit(std::get<1>(*IP));

    scratch_data.fe_values_rhs.get_function_gradients(
      scratch_data.dofs_rhs, scratch_data.vectors_list_rhs);

    for (unsigned int q_index = 0; q_index < scratch_data.n_q_points; ++q_index)
      {
        for (unsigned int i = 0; i < scratch_data.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < scratch_data.dofs_per_cell; ++j)
              {
                copy_data.cell_matrix(i, j) += // Integral I_a
                  scratch_data.fe_values[scratch_data.ve].value(i, q_index) *
                  scratch_data.fe_values[scratch_data.ve].value(j,
                                                                q_index) *
                  scratch_data.fe_values.JxW(q_index); // Ni.Nj dS
              }

            copy_data.cell_rhs(i) += // Integral I_b
              (scratch_data.vectors_list_rhs[q_index][1] *
                 scratch_data.fe_values[scratch_data.ve].value(i, q_index)[0] -
               scratch_data.vectors_list_rhs[q_index][0] *
                 scratch_data.fe_values[scratch_data.ve].value(i, q_index)[1]) *
              scratch_data.fe_values.JxW(q_index); // (curlv T).Ni
          }
      }

    std::get<0>(*IP)->get_dof_indices(copy_data.local_dof_indices);
  }
} // namespace ProjectorHgradToHdiv

// @sect3{Projector from H(curl) to L2}

// This name space contains all the code related to the computation of the
// magnetic field, $B$.
namespace ProjectorHcurlToL2
{
  using namespace BaseClasses;

  // We derive the solver from the `BaseSolver` class. What is left to do is
  // to initialize the `BaseSolver` properly and override two virtual functions,
  // `make_mesh` and `system_matrix_local`.
  class Solver : public BaseSolver
  {
  public:
    Solver() = delete;
    Solver(const unsigned int      p, // Degree of the FE_DGQ finite elements.
           const unsigned int      mapping_degree,
           const Triangulation<2> &triangulation_rhs,
           const DoFHandler<2>    &dof_handler_rhs,
           const Vector<double>   &solution_rhs,
           const std::string      &fname          = "data",
           const Function<2>      *exact_solution = nullptr);

    ~Solver() = default;

  private:
    virtual void make_mesh() override final;
    virtual void
    system_matrix_local(const IteratorPair  &IP,
                        AssemblyScratchData &scratch_data,
                        AssemblyCopyData    &copy_data) override final;

    FE_DGQ<2> fe;
  };

  // This solver is used at the fourth stage, see the introduction.
  // Consequently, we use the second constructor of the `BaseSolver` class.
  Solver::Solver(const unsigned int      p,
                 const unsigned int      mapping_degree,
                 const Triangulation<2> &triangulation_ext,
                 const DoFHandler<2>    &dof_handler_ext,
                 const Vector<double>   &solution_ext,
                 const std::string      &fname,
                 const Function<2>      *exact_solution)
    : BaseSolver(triangulation_ext,
                 dof_handler_ext,
                 solution_ext,
                 SolverType::Projector,
                 mapping_degree,
                 {update_values | update_JxW_values, update_gradients},
                 0.0,
                 fname,
                 exact_solution,
                 cg_settings_global)
    , fe(p)
  {}

  // At the fourth stage we do not upload the mesh. We reuse the mesh uploaded
  // at the first stage. Consequently, we just need to distribute the dofs.
  void Solver::make_mesh()
  {
    dof_handler.reinit(triangulation_rhs);
    dof_handler.distribute_dofs(fe);
  }

  // This function assembles a fraction of the system matrix and the system
  // right-hand side related to a single cell. These fractions are
  // `copy_data.cell_matrix` and `copy_data.cell_rhs`. They are copied into
  // the system matrix, $A_{ij}$, and the right-hand side, $b_i$, by the
  // function `Solver::copy_local_to_global()`.
  void Solver::system_matrix_local(const IteratorPair  &IP,
                                   AssemblyScratchData &scratch_data,
                                   AssemblyCopyData    &copy_data)
  {
    copy_data.cell_matrix.reinit(scratch_data.dofs_per_cell,

                                 scratch_data.dofs_per_cell);

    copy_data.cell_rhs.reinit(scratch_data.dofs_per_cell);

    copy_data.local_dof_indices.resize(scratch_data.dofs_per_cell);

    scratch_data.fe_values.reinit(std::get<0>(*IP));
    scratch_data.fe_values_rhs.reinit(std::get<1>(*IP));

    scratch_data.fe_values_rhs.get_function_gradients(
      scratch_data.dofs_rhs, scratch_data.vectors_vectors_list_rhs);

    for (unsigned int q_index = 0; q_index < scratch_data.n_q_points; ++q_index)
      {
        for (unsigned int i = 0; i < scratch_data.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < scratch_data.dofs_per_cell; ++j)
              {
                copy_data.cell_matrix(i, j) += // Integral I_a
                  scratch_data.fe_values[scratch_data.se].value(i, q_index) *
                  scratch_data.fe_values[scratch_data.se].value(j, q_index) *
                  scratch_data.fe_values.JxW(q_index); // Ni Nj dS
              }

            copy_data.cell_rhs(i) += // Integral I_b
              (scratch_data.vectors_vectors_list_rhs[q_index][1][0] -
               scratch_data.vectors_vectors_list_rhs[q_index][0][1]) *
              scratch_data.fe_values[scratch_data.se].value(i, q_index) *
              scratch_data.fe_values.JxW(q_index); // (curls A) Ni dS
          }
      }

    std::get<0>(*IP)->get_dof_indices(copy_data.local_dof_indices);
  }
} // namespace ProjectorHcurlToL2


class MagneticProblem
{
public:
  void run()
  {
    if (Settings::n_threads_max)
      MultithreadInfo::set_thread_limit(Settings::n_threads_max);

    MainOutputTable table_T(2);
    MainOutputTable table_Jf(2);
    MainOutputTable table_A(2);
    MainOutputTable table_B(2);

    table_T.clear();
    table_Jf.clear();
    table_A.clear();
    table_B.clear();

    std::cout << "Solving for (p' = " << Settings::fe_degree + 1
              << "; p = " << Settings::fe_degree << "): " << std::flush;

    for (unsigned int r = 1; r < 5; r++) // Mesh refinement parameter.
      {
        table_T.add_value("r", r);
        table_T.add_value("p", Settings::fe_degree + 1);

        table_Jf.add_value("r", r);
        table_Jf.add_value("p", Settings::fe_degree);

        table_A.add_value("r", r);
        table_A.add_value("p", Settings::fe_degree);

        table_B.add_value("r", r);
        table_B.add_value("p", Settings::fe_degree);

        // Stage 1. Computing $T$.

        std::cout << "T " << std::flush;

        ExactSolutions::CurrentVectorPotential T_exact;

        SolverT::Solver T(Settings::fe_degree + 1,
                          r,
                          Settings::mapping_degree,
                          "T_p" + std::to_string(Settings::fe_degree + 1) +
                            "_r" + std::to_string(r),
                          &T_exact);

        T.run();
        table_T.add_value("ndofs", T.get_n_dofs());
        table_T.add_value("ncells", T.get_n_cells());
        table_T.add_value("L2", T.get_L2_norm());

        // Stage 2. Computing $\vec{J}_f$.

        std::cout << "Jf " << std::flush;

        ExactSolutions::FreeCurrentDensity Jf_exact;

        ProjectorHgradToHdiv::Solver Jf(Settings::fe_degree,
                                        Settings::mapping_degree,
                                        T.get_tria(),
                                        T.get_dof_handler(),
                                        T.get_solution(),
                                        "Jf_p" +
                                          std::to_string(Settings::fe_degree) +
                                          "_r" + std::to_string(r),
                                        &Jf_exact);

        Jf.run();

        table_Jf.add_value("ndofs", Jf.get_n_dofs());
        table_Jf.add_value("ncells", Jf.get_n_cells());
        table_Jf.add_value("L2", Jf.get_L2_norm());

        // Stage 3. Computing $\vec{A}$.

        std::cout << "A " << std::flush;

        ExactSolutions::MagneticVectorPotential A_exact;

        SolverA::Solver A(Settings::fe_degree,
                          Settings::mapping_degree,
                          T.get_tria(),
                          T.get_dof_handler(),
                          T.get_solution(),
                          Settings::eta_squared,
                          "A_p" + std::to_string(Settings::fe_degree) + "_r" +
                            std::to_string(r),
                          &A_exact);

        A.run();

        table_A.add_value("ndofs", A.get_n_dofs());
        table_A.add_value("ncells", A.get_n_cells());
        table_A.add_value("L2", A.get_L2_norm());

        // Stage 4. Computing $B$.

        std::cout << "B " << std::flush;

        ExactSolutions::MagneticField B_exact;

        ProjectorHcurlToL2::Solver B(Settings::fe_degree,
                                     Settings::mapping_degree,
                                     T.get_tria(),
                                     A.get_dof_handler(),
                                     A.get_solution(),
                                     "B_p" +
                                       std::to_string(Settings::fe_degree) +
                                       "_r" + std::to_string(r),
                                     &B_exact);
        B.run();

        table_B.add_value("ndofs", B.get_n_dofs());
        table_B.add_value("ncells", B.get_n_cells());
        table_B.add_value("L2", B.get_L2_norm());
        // End stage 4.
      }

    table_T.save("table_T_p" + std::to_string(Settings::fe_degree + 1));
    table_Jf.save("table_Jf_p" + std::to_string(Settings::fe_degree));
    table_A.save("table_A_p" + std::to_string(Settings::fe_degree));
    table_B.save("table_B_p" + std::to_string(Settings::fe_degree));
    std::cout << std::endl;
  }
};

int main()
{
  try
    {
      ParameterHandler parameters;



      MagneticProblem problem;
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

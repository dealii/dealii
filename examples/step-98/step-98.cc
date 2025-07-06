/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2024-2026 by the deal.II authors
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
 * This program was contributed by Siarhei Uzunbajakau, www.cembooks.nl, 2026.
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
// The boundary, material, and manifold IDs are set in the `circle.geo` file.
// The IDs listed below must match the corresponding IDs set in the `circle.geo`
// file. If `project_exact_solution = true`, the program projects the exact
// solutions for $T$, $\vec{J}_f$, $\vec{A}$, and $B$ onto the corresponding
// function spaces and saves the results into the corresponding `.vtu` files
// next to the numerical solutions. By default, this feature is switched off.
// It is used only in the
// [Possibilities for extensions](@ref Step98_PossibilitiesForExtensions)
// section.
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

  const double eta_squared = 0.0; // eta^2 when solving for A.

  const unsigned int n_threads_max = 0; // If >0 limits the number of threads.

  const bool project_exact_solution = false; // Save the exact solution.
} // namespace Settings

// @sect3{Convergence table}
// The following class describes a convergence table. The convergence tables are
// saved on disk in TeX format.
class MainOutputTable : public ConvergenceTable
{
public:
  MainOutputTable() = delete;

  MainOutputTable(const unsigned int dim)
    : ConvergenceTable()
    , dim(dim)
  {}

  void save(const std::string &file_name)
  {
    set_precision("L2", 2);

    set_scientific("L2", true);

    evaluate_convergence_rates("L2",
                               "ncells",
                               ConvergenceTable::reduction_rate_log2,
                               dim);

    set_tex_caption("p", "p");
    set_tex_caption("r", "r");
    set_tex_caption("ncells", "nr. cells");
    set_tex_caption("ndofs", "nr. dofs");
    set_tex_caption("L2", "L2 norm");

    set_column_order({"p", "r", "ncells", "ndofs", "L2"});

    std::ofstream ofs(file_name + ".tex");
    write_tex(ofs);
  }

private:
  const unsigned int dim;
};

// @sect3{Equations}
// The following namespace contains closed-form analytical expressions
// for $T$, $\vec{J}_f$, $\vec{A}$, $B$, mentioned in the introduction to this
// tutorial.
namespace ExactSolutions
{

  // The following function describes the free-current density, $\vec{J}_f$,
  // inside the current region. The current density in this tutorial is
  // implemented by `ExactSolutions::FreeCurrentDensity` class and by
  // `SolverT::Solver::free_current_density` member function . There is a
  // subtle difference in how these two implementation compute the free-current
  // density. Both classes, however, utilize the same expression for the
  // free-current density. This function describes the expression.
  inline Tensor<1, 2> volume_free_current_density(const Point<2> &p,
                                                  const double    K0)
  {
    return Tensor<1, 2>({-K0 * p[1], K0 * p[0]});
  }

  // The following class implements the closed-form analytical
  // [expression](@ref Step98_Equation_1_Jf) for the free-current density,
  // $\vec{J}_f$, in the entire domain. The free-current density is computed
  // purely on the basis of the spatial coordinates of the field point. In go
  // coordinates, out comes the free-current density. The information on the
  // material ID of the mesh cells and any other information on the mesh is
  // ignored. This function is used for computing $L_2$ error norms and for
  // computing the projected exact solution. The $\vec{J}_f$ on the right-hand
  // side of the div-grad equation is implemented by the member function
  // `SolverT::Solver::free_current_density`
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

      for (unsigned int i = 0; i < values.size(); i++)
        {
          const double r = p[i].norm();

          if ((r >= Settings::a2) && (r <= Settings::b2))
            {
              const Tensor<1, 2> Jf =
                volume_free_current_density(p[i], Settings::K0);

              values[i][0] = Jf[0];
              values[i][1] = Jf[1];
            }
          else
            values[i] = 0;
        }
    }
  };

  // The following class implements the closed-form analytical
  // [expression](@ref Step98_Equation_4_A)
  // for the magnetic vector potential, $\vec{A}$.
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

      for (unsigned int i = 0; i < values.size(); i++)
        {
          const double r = p[i].norm();
          double       A;

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

  // The following class implements the closed-form analytical
  // [expression](@ref Step98_Equation_3_T) for
  // the current vector potential, $T$.
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

      for (unsigned int i = 0; i < values.size(); i++)
        {
          const double r = p[i].norm();
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

  // The following class implements the closed-form analytical
  // [expression](@ref Step98_Equation_2_B)
  // for the magnetic field, $B$.
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

      for (unsigned int i = 0; i < values.size(); i++)
        {
          const double r = p[i].norm();
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

// @sect3{Base Solver}

// As discussed above, the following namespace aggregates the code common to all
// four solvers used in the tutorial. All four solvers are derived from the
// `BaseSolver` class.
namespace BaseClasses
{
  // The computation of the integrals of the functionals is delegated to the
  // derived classes. The function that computes the integrals,
  // `BaseSolver::system_matrix_local`, is virtual and must be overridden by
  // the derived classes. However, the objects of the types FEValues are
  // initialized at the level of the `BaseSolver` class together with
  // `AssemblyScratchData`. For the initialization to work properly, the
  // `BaseSolver` class must know which cell data to compute for a particular
  // implementation of the solver down the hierarchy. This information is
  // communicated to the `BaseSolver` by passing an argument of the type
  // `UpdateFlagsCollection` to the constructor. In all solvers, with exception
  // of `SolverT`, we have two types of finite elements. One type of finite
  // elements models the solution to the partial differential equation. Another
  // models the physical quantity on the right-hand side of the partial
  // differential equation. The `solution_update_flags` data member below
  // contains the flags for updating the values of the finite elements that
  // model the solution. The `rhs_update_flags` data member contains the flags
  // for updating the values of the finite elements that model the physical
  // quantity on the right-hand side of the partial differential equation.
  struct UpdateFlagsCollection
  {
    UpdateFlags solution_update_flags;
    UpdateFlags rhs_update_flags;
  };

  // Each iteration of the program consists of
  // [four stages](@ref Step98_FourStages). Each stage utilizes one solver. All
  // four solvers used in this tutorial are derived from the following class.
  // The solver used in the first stage loads the mesh and refines it if
  // necessary. The solvers in the other three stages reuse the mesh prepared at
  // the first stage. Furthermore, the solver in the first stage expects a
  // closed-form analytical expression on the right-hand side of the partial
  // differential equation. Each solver in the other three stages expects a
  // potential computed at one of the preceding stages, i.e., a field in a form
  // of linear superposition of the shape functions. At the top of the following
  // class, we declare two constructors. The first constructor must be used for
  // constructing solver for the first stage. The second constructor must be
  // used for constructing the solvers for the second, third, and fourth stages.
  // The second constructor has three extra arguments, `triangulation_rhs`,
  // `dof_handler_rhs`, and `solution_rhs` for accommodating the mesh and the
  // potential computed at one of the preceding stages.
  //
  // Following the constructors are the eight functions that implement various
  // steps typical for every solver. The function `run` aggregates these steps.
  // This arrangement of functions is quite standard in deal.II, see step-3,
  // for instance. Normally, the `setup` function begins by distributing the
  // dofs. The problem is: the type of the finite elements is not known in the
  // `BaseSolver` class. It is specified in the derived classes. Therefore, it
  // could be reasonable to make the setup function virtual as well. Instead,
  // we move the dof distribution code to the end of the `make_mesh` function
  // which is, in fact, virtual and must be overridden in the derived classes
  // anyway.
  //
  // After that, we declare six get functions that provide access to protected
  // data. These functions are called from outside the solver, see
  // `MagneticProblem::run` function. The data provided by these get-functions
  // is used to fill the convergence tables and to pass the references to the
  // triangulation, the dof handler, and the dofs between solvers.
  class BaseSolver
  {
  public:
    BaseSolver(const unsigned int          mapping_degree,
               const UpdateFlagsCollection update_flags_collection,
               const std::string          &file_name,
               const Function<2>          *exact_solution);

    BaseSolver(const Triangulation<2>     &triangulation_rhs,
               const DoFHandler<2>        &dof_handler_rhs,
               const Vector<double>       &solution_rhs,
               const unsigned int          stage,
               const unsigned int          mapping_degree,
               const UpdateFlagsCollection update_flags_collection,
               const std::string          &file_name,
               const Function<2>          *exact_solution);

    virtual void make_mesh() = 0; /* If used in the first stage, loads the mesh
                                     and distributes the dofs. In other stages
                                     just distributes the dofs.*/
    void setup();    /* Applies Dirichlet boundary condition, setups the
                        dof pattern, and initializes vectors, matrices. */
    void assemble(); // Assembles the system of linear equations.
    void solve();    // Solves the system of linear equations.
    void output_results() const;       // Saves the result into a `.vtu` file.
    void compute_error_norms();        // Computes L^2 error norm.
    void project_exact_solution_fcn(); // Projects exact solution.
    void clear();                      // Clears the memory for the next solver.
    void run(); /* Executes the last eight functions in the proper order
                   and measures the execution time for each function. */

    double                  get_L2_norm() const;
    unsigned int            get_n_cells() const;
    types::global_dof_index get_n_dofs() const;
    const Triangulation<2> &get_tria() const;
    const DoFHandler<2>    &get_dof_handler() const;
    const Vector<double>   &get_solution() const;

    // We begin the `protected` section of the `BaseSolver` class by declaring
    // three data members that store the input from one of the preceding
    // solvers. If there is no preceding solver and these data is not provided,
    // i.e., first of the two constructors above has been used, these three data
    // members point to `triangulation`, `dof_handler`, and `solution`, of the
    // current solver.
    //
    // Following are the three data members that store the triangulation, dof
    // handler, and dofs vector of the current solver. In the case of the
    // first-stage solver, the `triangulation` data member stores the loaded and
    // refined mesh. This data member is not used in the case of the solvers of
    // the second, third, and fourth stages. The `dof_handler` and `solution`
    // represent the result of the solver, i.e., the computed potential or
    // field.
    //
    // Next, we declare four data members that describe the system of linear
    // equations to be solved by the linear solver. The components of the system
    // matrix and that of the right-hand side are computed by the `assemble`
    // function. The affine constraints are used to apply the Dirichlet boundary
    // conditions and to distribute the local (cell specific) system matrix and
    // right-hand side to `system_matrix` and `system_rhs`. The are no hanging
    // nodes and hanging node constraints in this program. The last data member
    // in this block describes the dynamic sparsity pattern. The tutorial step-2
    // discusses the rationale behind the dynamic sparsity pattern.
    //
    // The next block contains two data members related to the exact solution.
    // The data member `exact_solution` points to the closed-form analytical
    // solution the solver attempts to compute.
    // If `Settings::project_exact_solution=true`, the exact solution is
    // projected onto a proper function space and the data member
    // `projected_exact_solution` is populated by the dofs of the projected
    // exact solution. The corresponding dof handler is `dof_handler`. Together
    // `dof_handler` and `projected_exact_solution` constitute the field
    // function which describes the exact solution. It is saved into the `.vtu`
    // file next to the solution and the $L_2$ error norm.
    //
    // The next four data members simply store the data supplied as arguments
    // to the constructor. The `stage` data member stores the number of the
    // current stage. The `mapping_degree` data member contains the degree of
    // mapping from the reference cell to a mesh cell and back. The
    // `update_flags_collection` contains the information on which finite
    // element values must be computed for each cell. The names of the output
    // files are derived by appending strings to `file_name`.
    //
    // The next block contains two data members computed by the function
    // `BaseSolver::compute_error_norms`. The `L2_per_cell` data member
    // contains one value of the $L^2$ error norm per mesh cell. It is saved
    // into the `.vtu` file next to the solution. The `L2_norm` data member
    // contains one value of the $L^2$ error norm per mesh. It is reported in
    // the convergence table.
  protected:
    const Triangulation<2> &triangulation_rhs;
    const DoFHandler<2>    &dof_handler_rhs;
    const Vector<double>   &solution_rhs;

    Triangulation<2> triangulation;
    DoFHandler<2>    dof_handler;
    Vector<double>   solution;

    SparseMatrix<double>      system_matrix;
    Vector<double>            system_rhs;
    AffineConstraints<double> constraints;
    SparsityPattern           sparsity_pattern;

    const Function<2> *exact_solution;
    Vector<double>     projected_exact_solution;

    const unsigned int          stage;
    const unsigned int          mapping_degree;
    const UpdateFlagsCollection update_flags_collection;
    const std::string           file_name;

    Vector<double> L2_per_cell;
    double         L2_norm;

    // The program utilizes the WorkStream technology. The step-9 tutorial
    // does a much better job of explaining the workings of WorkStream.
    // Reading the @ref workstream_paper "WorkStream paper" is recommended.
    // In very simple terms, the workings of the WorkStream can be envisioned as
    // the following. Let us assume we have a task of computing components of
    // the system matrix, $A_{ij}$, and the right-hand side vector, $b_i$.
    // Simply put, we need to fill in the matrix `system_matrix` and vector
    // `system_rhs`. This is a big task as the number of dofs is large. The
    // idea is to split the big task on a number of small tasks and feed them
    // to multiple threads to speed up the calculation process. Each small
    // task consists of computing the contributions of a single mesh cell to
    // `system_matrix` and `system_rhs`. These contributions are stored
    // temporary in `cell_matrix` and `cell_rhs` for each cell. These
    // contributions are then copied to `system_matrix` and `system_rhs`.
    // WorkStream creates and schedules the small tasks and takes care
    // of copying `cell_matrix` and `cell_rhs` to `system_matrix` and
    // `system_rhs`. The rest of the declarations in the `protected` section
    // of the `BaseSolver` class help to communicate to WorkStream information
    // which is necessary for its operation.
    //
    // First, the `CellIteratorPair` type is declared in the block of code
    // below. The solvers in the second, third, and fourth stages use two dof
    // handlers, `dof_handler` and `dof_handler_rhs`. The WorkStream needs to
    // walk through the two dof handlers synchronously. For this purpose we
    // pair two active cell iterators (one from `dof_handler`, another from
    // `dof_handler_rhs`). For that we need the `CellIteratorPair` type. The
    // solver at the first stage (a  solver constructed by invoking the first
    // constructor above) uses only one dof handler. In this case the
    // constructor makes `dof_handler_rhs` to reference `dof_handler`. In
    // effect, both iterators of the tuple will iterate the same dof handler.
    // In the case of the first-stage solver we will use only the first
    // iterator.
    //
    // Next, we declare the `AssemblyScratchData` type. An object of this type
    // contains the relevant information on the current mesh cell which is used
    // as an input for computing the components of the system matrix and the
    // right-hand side. WorkStream creates an object of this type and passes
    // it to the function `system_matrix_local` which, in turn, computes
    // components of `cell_matrix` and `cell_rhs`.
    //
    // Next, the type `AssemblyCopyData` is declared. An objects of this type
    // contains the cell specific contributions to the system matrix and the
    // right-hand side. The WorkStream creates object of this type and passes
    // it to function `system_matrix_local` along with an object of the type
    // `AssemblyScratchData`. The function `system_matrix_local`, in turn,
    // takes input data form the object of the type `AssemblyScratchData`,
    // computes the relevant integrals and places the result into the object
    // of the type `AssemblyCopyData`. The WorkStream copies the content of
    // the `AssemblyCopyData` object, `cell_matrix` and `cell_rhs`, into
    // `system_matrix` and `system_rhs`. The data member
    // `AssemblyCopyData::local_dof_indices` contains the indices of global
    // components to which cell-specific data must be copied.
    //
    // WorkStream calls the function `system_matrix_local` to compute the
    // cell-specific components. Likewise, WorkStream uses calls to
    // `copy_local_to_global` to copy the cell-specific data into the system
    // matrix and right-hand side. The functions `system_matrix_local` and
    // `copy_local_to_global` are declared last.
    using IteratorTuple =
      std::tuple<typename DoFHandler<2>::active_cell_iterator,
                 typename DoFHandler<2>::active_cell_iterator>;

    using CellIteratorPair = SynchronousIterators<IteratorTuple>;

    struct AssemblyScratchData
    {
      AssemblyScratchData(const DoFHandler<2>        &dof_handler,
                          const DoFHandler<2>        &dof_handler_rhs,
                          const Vector<double>       &dofs_rhs,
                          const unsigned int          mapping_degree,
                          const UpdateFlagsCollection update_flags_collection,
                          const unsigned int          stage);

      AssemblyScratchData(const AssemblyScratchData &scratch_data);

      MappingQ<2> mapping;

      FEValues<2> fe_values;
      FEValues<2> fe_values_rhs;

      const unsigned int dofs_per_cell;
      const unsigned int n_q_points;

      std::vector<double>                    permeability_list;
      std::vector<double>                    values_list_rhs;
      std::vector<Tensor<1, 2>>              vectors_list_rhs;
      std::vector<std::vector<Tensor<1, 2>>> vectors_vectors_list_rhs;

      const DoFHandler<2>  &dof_handler_rhs;
      const Vector<double> &dofs_rhs;
    };

    struct AssemblyCopyData
    {
      FullMatrix<double>                   cell_matrix;
      Vector<double>                       cell_rhs;
      std::vector<types::global_dof_index> local_dof_indices;
    };

    virtual void system_matrix_local(const CellIteratorPair &IP,
                                     AssemblyScratchData    &scratch_data,
                                     AssemblyCopyData       &copy_data) = 0;

    void copy_local_to_global(const AssemblyCopyData &copy_data);
  }; // class BaseSolver

  // The following are the implementations of the two constructors
  // of the `BaseSolver` class.
  BaseSolver::BaseSolver(const unsigned int          mapping_degree,
                         const UpdateFlagsCollection update_flags_collection,
                         const std::string          &file_name,
                         const Function<2>          *exact_solution)
    : triangulation_rhs(triangulation)
    , dof_handler_rhs(dof_handler)
    , solution_rhs(solution)
    , exact_solution(exact_solution)
    , stage(1)
    , mapping_degree(mapping_degree)
    , update_flags_collection(update_flags_collection)
    , file_name(file_name)
  {}

  BaseSolver::BaseSolver(const Triangulation<2>     &triangulation_rhs,
                         const DoFHandler<2>        &dof_handler_rhs,
                         const Vector<double>       &solution_rhs,
                         const unsigned int          stage,
                         const unsigned int          mapping_degree,
                         const UpdateFlagsCollection update_flags_collection,
                         const std::string          &file_name,
                         const Function<2>          *exact_solution)
    : triangulation_rhs(triangulation_rhs)
    , dof_handler_rhs(dof_handler_rhs)
    , solution_rhs(solution_rhs)
    , exact_solution(exact_solution)
    , stage(stage)
    , mapping_degree(mapping_degree)
    , update_flags_collection(update_flags_collection)
    , file_name(file_name)
  {}

  // The following function applies the Dirichlet boundary condition, sets
  // up a sparsity pattern, and initializes the vectors and matrices. It is
  // common for the setup function to distribute the dofs. The type of the
  // finite elements, however, is not known at the level of `BaseSolver` class.
  // The type of the finite elements is chosen in the derived classes. For
  // this reason, the task of distributing the dofs is shifted to the end of
  // the `make_mesh` function which is a virtual function.
  //
  // The Dirichlet boundary condition must be enforced only in the first
  // [stage](@ref Step98_FourStages)
  // where the div-grad equation is solved for the current vector potential, T.
  // For this reason we have `if (stage == 1)` filter in the beginning of the
  // function. The boundary value problem for the curl-curl equation utilizes
  // the Neumann boundary condition. It is a natural boundary condition. It is
  // enforced by minimization of the functional. The two projectors,
  // $T \rightarrow \vec{J}_f$ and $\vec{A} \rightarrow B$, use no boundary
  // conditions.
  void BaseSolver::setup()
  {
    constraints.clear();

    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    if (stage == 1)
      VectorTools::interpolate_boundary_values(MappingQ<2>(mapping_degree),
                                               dof_handler,
                                               Settings::outer_boundary_id,
                                               Functions::ZeroFunction<2>(),
                                               constraints);

    constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    sparsity_pattern.copy_from(dsp);

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
  // `system_matrix_local` which is a virtual function. That is to say, the
  // recipe for the functional is implemented in the derived class by overriding
  // the virtual function `system_matrix_local`.
  void BaseSolver::assemble()
  {
    WorkStream::run(CellIteratorPair({dof_handler.begin_active(),
                                      dof_handler_rhs.begin_active()}),
                    CellIteratorPair(
                      {dof_handler.end(), dof_handler_rhs.end()}),
                    *this,
                    &BaseSolver::system_matrix_local,
                    &BaseSolver::copy_local_to_global,
                    AssemblyScratchData(dof_handler,
                                        dof_handler_rhs,
                                        solution_rhs,
                                        mapping_degree,
                                        update_flags_collection,
                                        stage),
                    AssemblyCopyData());
  }

  // The following are the implementation of the constructors of the
  // `AssemblyScratchData`. The first constructor initializes the scratch
  // data from the input parameters. The second - from another object of
  // the same type, i.e., a copy constructor.
  BaseSolver::AssemblyScratchData::AssemblyScratchData(
    const DoFHandler<2>        &dof_handler,
    const DoFHandler<2>        &dof_handler_rhs,
    const Vector<double>       &dofs_rhs,
    const unsigned int          mapping_degree,
    const UpdateFlagsCollection update_flags_collection,
    const unsigned int          stage)
    : mapping(mapping_degree)
    , fe_values(mapping,
                dof_handler.get_fe(),
                QGauss<2>((stage == 1) ? (dof_handler.get_fe().degree + 1) :
                                         (dof_handler.get_fe().degree + 2)),
                update_flags_collection.solution_update_flags)
    , fe_values_rhs(mapping,
                    dof_handler_rhs.get_fe(),
                    QGauss<2>((stage == 1) ? (dof_handler.get_fe().degree + 1) :
                                             (dof_handler.get_fe().degree + 2)),
                    update_flags_collection.rhs_update_flags)
    , dofs_per_cell(fe_values.dofs_per_cell)
    , n_q_points(fe_values.get_quadrature().size())
    , permeability_list(n_q_points)
    , values_list_rhs(n_q_points)
    , vectors_list_rhs(n_q_points)
    , vectors_vectors_list_rhs(n_q_points, std::vector<Tensor<1, 2>>(2))
    , dof_handler_rhs(dof_handler_rhs)
    , dofs_rhs(dofs_rhs)
  {}

  BaseSolver::AssemblyScratchData::AssemblyScratchData(
    const AssemblyScratchData &scratch_data)
    : mapping(scratch_data.mapping.get_degree())
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
    , dof_handler_rhs(scratch_data.dof_handler_rhs)
    , dofs_rhs(scratch_data.dofs_rhs)
  {}

  // The following function copies the components of a cell matrix and a cell
  // right-hand side into the system matrix, $A_{ij}$, and the system right-hand
  // side, $b_i$.
  void BaseSolver::copy_local_to_global(const AssemblyCopyData &copy_data)
  {
    constraints.distribute_local_to_global(copy_data.cell_matrix,
                                           copy_data.cell_rhs,
                                           copy_data.local_dof_indices,
                                           system_matrix,
                                           system_rhs);
  }

  // The following function solves the system of linear equations.
  // The stopping condition for the iteration algorithm is
  // $\|\boldsymbol{b}-\boldsymbol{A}\boldsymbol{c}\|<10^{-8}\|\boldsymbol{b}\|$.
  // The maximum number of iteration steps is set to `system_rhs.size()` as
  // the conjugate gradient algorithm is supposed to find the solution in at
  // most $m$ steps for an $m \times m$ system matrix. This function also
  // distributes constraints. The constraints are only used to enforce the
  // Dirichlet boundary condition.
  void BaseSolver::solve()
  {
    SolverControl control(system_rhs.size(), 1e-8 * system_rhs.l2_norm());

    GrowingVectorMemory<Vector<double>> memory;
    SolverCG<Vector<double>>            cg(control, memory);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    constraints.distribute(solution);
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

        constraints_empty.clear();

        DoFTools::make_hanging_node_constraints(dof_handler, constraints_empty);

        constraints_empty.close();

        VectorTools::project(MappingQ<2>(mapping_degree),
                             dof_handler,
                             constraints_empty,
                             QGauss<2>((stage == 1) ?
                                         (dof_handler.get_fe().degree + 1) :
                                         (dof_handler.get_fe().degree + 2)),
                             *exact_solution,
                             projected_exact_solution);
      }
  }

  // The following function saves the solution and the $L_2$ error norm into a
  // `.vtu` file. If `Settings::project_exact_solution = true`, the projected
  // exact solution is saved as well.
  void BaseSolver::output_results() const
  {
    const bool fe_is_vector = ((stage == 2) || (stage == 3));

    const std::string  name = (fe_is_vector ? "VectorField" : "ScalarField");
    const unsigned int components = (fe_is_vector ? 2 : 1);
    const DataComponentInterpretation::DataComponentInterpretation
      component_interpretation =
        (fe_is_vector ?
           DataComponentInterpretation::component_is_part_of_vector :
           DataComponentInterpretation::component_is_scalar);

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

    std::ofstream ofs(file_name + ".vtu");
    data_out.write_vtu(ofs);
  }

  // The following function clears the memory for the next solver.
  void BaseSolver::clear()
  {
    system_matrix.clear();
    system_rhs.reinit(0);
  }

  // The following functions calls all the constituent functions in the right
  // order.
  void BaseSolver::run()
  {
    make_mesh();
    setup();
    assemble();
    solve();
    compute_error_norms();
    project_exact_solution_fcn();
    output_results();
    clear();
  }

  // The following is the straightforward implementation of the get functions.
  double BaseSolver::get_L2_norm() const
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

// @sect3{Solver - T}

// The following namespace contains the code related to the computation of the
// current vector potential, $T$.
namespace SolverT
{
  using namespace BaseClasses;

  // We derive the solver from the `BaseSolver` class. What is left to do is
  // to initialize the `BaseSolver`, override the two virtual functions
  // (`make_mesh` and `system_matrix_local`), and implement the function
  // `free_current_density`. The function `free_current_density` implements
  // the closed-form analytical expression for $\vec{J}_f$ on the right-hand
  // side of the [div-grad equation](@ref Step98_FourStages).
  class Solver : public BaseSolver
  {
  public:
    Solver() = delete;
    Solver(const unsigned int p, // Degree of the FE_Q finite elements.
           const unsigned int r, // The mesh refinement parameter.
           const unsigned int mapping_degree,
           const std::string &file_name      = "data",
           const Function<2> *exact_solution = nullptr);

  private:
    virtual void make_mesh() override final;
    virtual void
    system_matrix_local(const CellIteratorPair &IP,
                        AssemblyScratchData    &scratch_data,
                        AssemblyCopyData       &copy_data) override final;

    const FE_Q<2>      fe;
    const unsigned int refinement_parameter;

    void free_current_density(const std::vector<Point<2>> &p,
                              const types::material_id     material_id,
                              std::vector<Tensor<1, 2>>   &values) const;
  };

  // Following is the implementation of the constructor.
  // We use the first constructor of the `BaseSolver` class as
  // the solver is used at the
  // [first stage](@ref Step98_FourStages). By looking at the
  // [expressions](@ref Step98_Numerical_Recipe_T)
  // for $A_{ij}$ and $b_i$ we can conclude that to compute them we need
  // gradients of the shape functions, quadrature points, and the quadrature
  // weights multiplied by the Jacobian determinant(`JxW`). The quadrature
  // points are needed to sample the closed-form analytical expression for
  // $\vec{J}_f$. Accordingly, we use the update flags
  // `update_quadrature_points`, `update_gradients`, and `update_JxW_values`
  // for the FE_Q finite elements. This time we do not use the finite elements
  // that model the potential on the right-hand side of the equation, so we
  // set `update_default` for the right-hand side finite elements.
  Solver::Solver(const unsigned int p,
                 const unsigned int r,
                 const unsigned int mapping_degree,
                 const std::string &file_name,
                 const Function<2> *exact_solution)
    : BaseSolver(mapping_degree,
                 {update_quadrature_points | update_gradients |
                    update_JxW_values,
                  update_default},
                 file_name,
                 exact_solution)

    , fe(p)
    , refinement_parameter(r)
  {}

  // The following function loads the mesh, creates manifolds, bounds the
  // manifolds to the manifold IDs, refines the mesh, and distributes the dofs.
  // Note that we are allowed to create the manifolds locally in the function as
  // the triangulation object keeps copies of the manifolds, see step-65. Also
  // recall that we have shifted the task of distributing the dofs from the
  // `setup` function to the `make_mesh` function to evade the necessity to make
  // the `setup` function virtual. This allows us to keep one `setup` function
  // in the `BaseSolver` class that serves the needs of all derived classes.
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

  // The following function assembles a fraction of
  // [the system matrix and the system right-hand side](@ref Step98_Numerical_Recipe_T)
  // related to a single cell. These fractions are
  // `copy_data.cell_matrix` and `copy_data.cell_rhs`. They are copied to
  // `system_matrix` and `system_rhs` by WorkStream.
  void Solver::system_matrix_local(const CellIteratorPair &IP,
                                   AssemblyScratchData    &scratch_data,
                                   AssemblyCopyData       &copy_data)
  {
    const FEValuesExtractors::Scalar se(0);

    copy_data.cell_matrix.reinit(scratch_data.dofs_per_cell,
                                 scratch_data.dofs_per_cell);

    copy_data.cell_rhs.reinit(scratch_data.dofs_per_cell);

    copy_data.local_dof_indices.resize(scratch_data.dofs_per_cell);

    const typename DoFHandler<2>::active_cell_iterator cell = std::get<0>(*IP);

    scratch_data.fe_values.reinit(cell);

    Solver::free_current_density(scratch_data.fe_values.get_quadrature_points(),
                                 cell->material_id(),
                                 scratch_data.vectors_list_rhs);

    for (unsigned int q_index = 0; q_index < scratch_data.n_q_points; ++q_index)
      {
        for (unsigned int i = 0; i < scratch_data.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < scratch_data.dofs_per_cell; ++j)
              {
                copy_data.cell_matrix(i, j) += // Integral I_a1.
                  (scratch_data.fe_values[se].gradient(
                     i, q_index) * // grad phi_i(x_q)
                   scratch_data.fe_values[se].gradient(
                     j, q_index) // grad phi_j(x_q)
                   ) *
                  scratch_data.fe_values.JxW(q_index); // dx
              }
            copy_data.cell_rhs(i) += // Integral I_b3-1.
              (scratch_data.vectors_list_rhs[q_index][0] *
                 scratch_data.fe_values[se].gradient(i, q_index)[1] -
               scratch_data.vectors_list_rhs[q_index][1] *
                 scratch_data.fe_values[se].gradient(i,
                                                     q_index)[0]) *
              scratch_data.fe_values.JxW(
                q_index); // J_f(x_q).(curlv phi_i(x_q))dx.
          }
      }

    cell->get_dof_indices(copy_data.local_dof_indices);
  }

  // The following function implements the closed form analytical
  // [expression](@ref Step98_Equation_1_Jf)
  // for $\vec{J}_f$ on the right-hand side of the
  // [div-grad equation](@ref Step98_PDE_T).
  void Solver::free_current_density(const std::vector<Point<2>> &p,
                                    const types::material_id     material_id,
                                    std::vector<Tensor<1, 2>>   &values) const
  {
    Assert(p.size() == values.size(),
           ExcDimensionMismatch(p.size(), values.size()));

    if ((material_id == Settings::material_id_free_space) ||
        (material_id == Settings::material_id_core))
      std::fill(values.begin(), values.end(), Tensor<1, 2>());

    if (material_id == Settings::material_id_free_current)
      for (unsigned int i = 0; i < values.size(); i++)
        values[i] =
          ExactSolutions::volume_free_current_density(p[i], Settings::K0);
  }
} // namespace SolverT

// @sect3{Projector from H(grad) to H(div)}

// The following namespace contains all the code related to the computation of
// the free-current density, $\vec{J}_f$.
namespace ProjectorHgradToHdiv
{
  using namespace BaseClasses;

  // We derive the solver from the `BaseSolver` class. What is left to do is
  // to initialize the `BaseSolver` and override the two virtual functions
  // (`make_mesh` and `system_matrix_local`).
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
           const std::string      &file_name      = "data",
           const Function<2>      *exact_solution = nullptr);

  private:
    virtual void make_mesh() override final;
    virtual void
    system_matrix_local(const CellIteratorPair &IP,
                        AssemblyScratchData    &scratch_data,
                        AssemblyCopyData       &copy_data) override final;

    FE_RaviartThomas<2> fe;
  };

  // Following is the implementation of the constructor.
  // We use the second constructor of the `BaseSolver` class as
  // the solver is used at the
  // [second stage](@ref Step98_FourStages). By looking at the
  // [expressions](@ref Step98_Numerical_Recipe_Jf)
  // for $A_{ij}$ and $b_i$ we can conclude that to compute them we need
  // values of the shape functions and the quadrature weights multiplied by the
  // Jacobian determinant(`JxW`) from the FE_RaviartThomas finite elements.
  // Accordingly, we use the update flags `update_values`, and
  // `update_JxW_values` for the FE_RaviartThomas finite elements. This time
  // there is a numerically computed potential, $T$, on the right-hand side of
  // the [equation](@ref Step98_PDE_Jf). It is modeled by the FE_Q finite elements.
  // To compute the right-hand side, we need gradients of the shape functions.
  // Accordingly, we use the update flag `update_gradients` for the FE_Q finite
  // elements.
  Solver::Solver(const unsigned int      p,
                 const unsigned int      mapping_degree,
                 const Triangulation<2> &triangulation_rhs,
                 const DoFHandler<2>    &dof_handler_rhs,
                 const Vector<double>   &solution_rhs,
                 const std::string      &file_name,
                 const Function<2>      *exact_solution)
    : BaseSolver(triangulation_rhs,
                 dof_handler_rhs,
                 solution_rhs,
                 2,
                 mapping_degree,
                 {update_values | update_JxW_values, update_gradients},
                 file_name,
                 exact_solution)
    , fe(p)
  {}

  // At the second stage we do not load the mesh. We reuse the mesh loaded
  // at the first stage. Consequently, we just need to distribute the dofs.
  void Solver::make_mesh()
  {
    dof_handler.reinit(triangulation_rhs);
    dof_handler.distribute_dofs(fe);
  }

  // The following function assembles a fraction of
  // [the system matrix and the system right-hand side](@ref Step98_Numerical_Recipe_Jf)
  // related to a single cell. These fractions are
  // `copy_data.cell_matrix` and `copy_data.cell_rhs`. They are copied to
  // `system_matrix` and `system_rhs` by WorkStream.
  void Solver::system_matrix_local(const CellIteratorPair &IP,
                                   AssemblyScratchData    &scratch_data,
                                   AssemblyCopyData       &copy_data)
  {
    const FEValuesExtractors::Vector ve(0);

    copy_data.cell_matrix.reinit(scratch_data.dofs_per_cell,
                                 scratch_data.dofs_per_cell);

    copy_data.cell_rhs.reinit(scratch_data.dofs_per_cell);

    copy_data.local_dof_indices.resize(scratch_data.dofs_per_cell);

    const typename DoFHandler<2>::active_cell_iterator cell = std::get<0>(*IP);
    const typename DoFHandler<2>::active_cell_iterator cell_rhs =
      std::get<1>(*IP);

    scratch_data.fe_values.reinit(cell);
    scratch_data.fe_values_rhs.reinit(cell_rhs);

    scratch_data.fe_values_rhs.get_function_gradients(
      scratch_data.dofs_rhs, scratch_data.vectors_list_rhs);

    for (unsigned int q_index = 0; q_index < scratch_data.n_q_points; ++q_index)
      {
        for (unsigned int i = 0; i < scratch_data.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < scratch_data.dofs_per_cell; ++j)
              {
                copy_data.cell_matrix(i, j) +=                   // Integral I_a
                  scratch_data.fe_values[ve].value(i, q_index) * // phi_i(x_q)
                  scratch_data.fe_values[ve].value(j,
                                                   q_index) * // phi_j(x_q)
                  scratch_data.fe_values.JxW(q_index);        // dx
              }

            copy_data.cell_rhs(i) += // Integral I_b
              (scratch_data.vectors_list_rhs[q_index][1] *
                 scratch_data.fe_values[ve].value(i, q_index)[0] -
               scratch_data.vectors_list_rhs[q_index][0] *
                 scratch_data.fe_values[ve].value(i, q_index)[1]) *
              scratch_data.fe_values.JxW(
                q_index); // [curlv T(x_q)].phi_i(x_q)dx
          }
      }

    cell->get_dof_indices(copy_data.local_dof_indices);
  }
} // namespace ProjectorHgradToHdiv



// @sect3{Solver - A}

// The following namespace contains all the code related to the computation of
// the magnetic vector potential, $\vec{A}$.
namespace SolverA
{
  using namespace BaseClasses;

  // We derive the solver from the `BaseSolver` class. What is left to do is
  // to initialize the `BaseSolver`, override the two virtual functions
  // (`make_mesh` and `system_matrix_local`), and implement the function
  // `permeability`.
  class Solver : public BaseSolver
  {
  public:
    Solver() = delete;
    Solver(const unsigned int p, // Degree of the FE_Nedelec finite elements.
           const unsigned int mapping_degree,
           const Triangulation<2> &triangulation_ext,
           const DoFHandler<2>    &dof_handler_ext,
           const Vector<double>   &solution_ext,
           const std::string      &file_name      = "data",
           const Function<2>      *exact_solution = nullptr);

  private:
    virtual void make_mesh() override final;
    virtual void
    system_matrix_local(const CellIteratorPair &IP,
                        AssemblyScratchData    &scratch_data,
                        AssemblyCopyData       &copy_data) override final;

    FE_Nedelec<2> fe;

    void permeability(const types::material_id material_id,
                      std::vector<double>     &values) const;
  };

  // Following is the implementation of the constructor.
  // We use the second constructor of the `BaseSolver` class as
  // the solver is used at the
  // [third stage](@ref Step98_FourStages). By looking at the
  // [expressions](@ref Step98_Numerical_Recipe_A)
  // for $A_{ij}$ and $b_i$ we can conclude that to compute them we need
  // values of the shape functions, their gradients, and the quadrature
  // weights multiplied by the Jacobian determinant(`JxW`) from the FE_Nedelec
  // finite elements. Accordingly, we use the update flags `update_values`,
  // `update_gradients`, and `update_JxW_values` for the FE_Nedelec finite
  // elements. This time there is a numerically computed potential, $T$, on
  // the right-hand side of the
  // [equation](@ref Step98_PDE_A). It is modeled by the FE_Q finite elements.
  // To compute it, we need values of the shape functions. Accordingly, we use
  // the update flag `update_values` for the FE_Q finite elements.
  Solver::Solver(const unsigned int      p,
                 const unsigned int      mapping_degree,
                 const Triangulation<2> &triangulation_rhs,
                 const DoFHandler<2>    &dof_handler_rhs,
                 const Vector<double>   &solution_rhs,
                 const std::string      &file_name,
                 const Function<2>      *exact_solution)
    : BaseSolver(triangulation_rhs,
                 dof_handler_rhs,
                 solution_rhs,
                 3,
                 mapping_degree,
                 {update_values | update_gradients | update_JxW_values,
                  update_values},
                 file_name,
                 exact_solution)
    , fe(p)
  {}

  // At the third stage we do not load the mesh. We reuse the mesh loaded at
  // the first stage. Consequently, we just need to distribute the dofs.
  void Solver::make_mesh()
  {
    dof_handler.reinit(triangulation_rhs);
    dof_handler.distribute_dofs(fe);
  }

  // The following function assembles a fraction of
  // [the system matrix and the system right-hand side](@ref Step98_Numerical_Recipe_A)
  // related to a single cell. These fractions are
  // `copy_data.cell_matrix` and `copy_data.cell_rhs`. They are copied to
  // `system_matrix` and `system_rhs` by WorkStream.
  void Solver::system_matrix_local(const CellIteratorPair &IP,
                                   AssemblyScratchData    &scratch_data,
                                   AssemblyCopyData       &copy_data)
  {
    const FEValuesExtractors::Vector ve(0);

    copy_data.cell_matrix.reinit(scratch_data.dofs_per_cell,
                                 scratch_data.dofs_per_cell);

    copy_data.cell_rhs.reinit(scratch_data.dofs_per_cell);

    copy_data.local_dof_indices.resize(scratch_data.dofs_per_cell);

    const typename DoFHandler<2>::active_cell_iterator cell = std::get<0>(*IP);
    const typename DoFHandler<2>::active_cell_iterator cell_rhs =
      std::get<1>(*IP);

    scratch_data.fe_values.reinit(cell);
    scratch_data.fe_values_rhs.reinit(cell_rhs);

    Solver::permeability(cell->material_id(), scratch_data.permeability_list);

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
                  (scratch_data.fe_values[ve].curl(
                     i, q_index) * // curls phi_i(x_q)
                     scratch_data.fe_values[ve].curl(
                       j, q_index) // curls phi_j(x_q)
                   +
                   Settings::eta_squared * // eta^2
                     scratch_data.fe_values[ve].value(i,
                                                      q_index) *  // phi_i(x_q)
                     scratch_data.fe_values[ve].value(j, q_index) // phi_j(x_q)
                   ) *
                  scratch_data.fe_values.JxW(q_index); // dx
              }
            copy_data.cell_rhs(i) += // Integral I_b3-1.
              (scratch_data.values_list_rhs[q_index] *
               scratch_data.fe_values[ve].curl(i, q_index)) *
              scratch_data.fe_values.JxW(q_index); // T(x_q)(curls phi_i(x_q))dx
          }
      }

    cell->get_dof_indices(copy_data.local_dof_indices);
  }

  // The following function implements the
  // [equation](@ref Step98_Equation_MU)
  // for permeability.
  void Solver::permeability(const types::material_id material_id,
                            std::vector<double>     &values) const
  {
    if (material_id == Settings::material_id_core)
      std::fill(values.begin(), values.end(), Settings::mu_1);
    else
      std::fill(values.begin(), values.end(), Settings::mu_0);
  }
} // namespace SolverA

// @sect3{Projector from H(curl) to L2}

// The following namespace contains all the code related to the computation of
// the magnetic field, $B$.
namespace ProjectorHcurlToL2
{
  using namespace BaseClasses;

  // We derive the solver from the `BaseSolver` class. What is left to do is
  // to initialize the `BaseSolver` and override the two virtual functions
  // (`make_mesh` and `system_matrix_local`).
  class Solver : public BaseSolver
  {
  public:
    Solver() = delete;
    Solver(const unsigned int      p, // Degree of the FE_DGQ finite elements.
           const unsigned int      mapping_degree,
           const Triangulation<2> &triangulation_rhs,
           const DoFHandler<2>    &dof_handler_rhs,
           const Vector<double>   &solution_rhs,
           const std::string      &file_name      = "data",
           const Function<2>      *exact_solution = nullptr);

  private:
    virtual void make_mesh() override final;
    virtual void
    system_matrix_local(const CellIteratorPair &IP,
                        AssemblyScratchData    &scratch_data,
                        AssemblyCopyData       &copy_data) override final;

    FE_DGQ<2> fe;
  };

  // Following is the implementation of the constructor.
  // We use the second constructor of the `BaseSolver` class as
  // the solver is used at the
  // [fourth stage](@ref Step98_FourStages). By looking at the
  // [expressions](@ref Step98_Numerical_Recipe_B)
  // for $A_{ij}$ and $b_i$ we can conclude that to compute them we need
  // values of the shape functions and the quadrature weights multiplied by the
  // Jacobian determinant(`JxW`) from the FE_DGQ finite elements.
  // Accordingly, we use the update flags `update_values`, and
  // `update_JxW_values` for the FE_DGQ finite elements. This time there is a
  // numerically computed potential, $\vec{A}$, on the right-hand side of the
  // [equation](@ref Step98_PDE_B). It is modeled by the FE_Nedelec finite
  // elements. To compute the right-hand side, we need gradients of the shape
  // functions. Accordingly, we use the update flag `update_gradients` for the
  // FE_Nedelec finite elements.
  Solver::Solver(const unsigned int      p,
                 const unsigned int      mapping_degree,
                 const Triangulation<2> &triangulation_ext,
                 const DoFHandler<2>    &dof_handler_ext,
                 const Vector<double>   &solution_ext,
                 const std::string      &file_name,
                 const Function<2>      *exact_solution)
    : BaseSolver(triangulation_ext,
                 dof_handler_ext,
                 solution_ext,
                 4,
                 mapping_degree,
                 {update_values | update_JxW_values, update_gradients},
                 file_name,
                 exact_solution)
    , fe(p)
  {}

  // At the fourth stage we do not load the mesh. We reuse the mesh loaded
  // at the first stage. Consequently, we just need to distribute the dofs.
  void Solver::make_mesh()
  {
    dof_handler.reinit(triangulation_rhs);
    dof_handler.distribute_dofs(fe);
  }

  // The following function assembles a fraction of
  // [the system matrix and the system right-hand side](@ref Step98_Numerical_Recipe_B)
  // related to a single cell. These fractions are
  // `copy_data.cell_matrix` and `copy_data.cell_rhs`. They are copied to
  // `system_matrix` and `system_rhs` by WorkStream.
  void Solver::system_matrix_local(const CellIteratorPair &IP,
                                   AssemblyScratchData    &scratch_data,
                                   AssemblyCopyData       &copy_data)
  {
    const FEValuesExtractors::Scalar se(0);

    copy_data.cell_matrix.reinit(scratch_data.dofs_per_cell,

                                 scratch_data.dofs_per_cell);

    copy_data.cell_rhs.reinit(scratch_data.dofs_per_cell);

    copy_data.local_dof_indices.resize(scratch_data.dofs_per_cell);

    const typename DoFHandler<2>::active_cell_iterator cell = std::get<0>(*IP);
    const typename DoFHandler<2>::active_cell_iterator cell_rhs =
      std::get<1>(*IP);

    scratch_data.fe_values.reinit(cell);
    scratch_data.fe_values_rhs.reinit(cell_rhs);

    scratch_data.fe_values_rhs.get_function_gradients(
      scratch_data.dofs_rhs, scratch_data.vectors_vectors_list_rhs);

    for (unsigned int q_index = 0; q_index < scratch_data.n_q_points; ++q_index)
      {
        for (unsigned int i = 0; i < scratch_data.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < scratch_data.dofs_per_cell; ++j)
              {
                copy_data.cell_matrix(i, j) +=                   // Integral I_a
                  scratch_data.fe_values[se].value(i, q_index) * // phi_i(x_q)
                  scratch_data.fe_values[se].value(j, q_index) * // phi_j(x_q)
                  scratch_data.fe_values.JxW(q_index);           // dx
              }

            copy_data.cell_rhs(i) += // Integral I_b
              (scratch_data.vectors_vectors_list_rhs[q_index][1][0] -
               scratch_data.vectors_vectors_list_rhs[q_index][0][1]) *
              scratch_data.fe_values[se].value(i, q_index) *
              scratch_data.fe_values.JxW(
                q_index); // (curls A(x_q)) phi_i(x_q) dx
          }
      }

    cell->get_dof_indices(copy_data.local_dof_indices);
  }
} // namespace ProjectorHcurlToL2

// @sect3{The main loop}

// The `MagneticProblem` class hosts the main loop inside the `run`
// function. The implementation of the loop is straightforward -
// we create and run the four solvers one-by-one and copy the relevant
// data into convergence tables.
class MagneticProblem
{
public:
  void run()
  {
    if (Settings::n_threads_max != 0)
      MultithreadInfo::set_thread_limit(Settings::n_threads_max);

    MainOutputTable table_T(2);
    MainOutputTable table_Jf(2);
    MainOutputTable table_A(2);
    MainOutputTable table_B(2);

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

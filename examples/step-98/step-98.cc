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
// fe_degree encodes the degree of the FE_Nedelec, FE_RaviartThomas, and
// FE_DGQ finite elements. The degree of the FE_Q finite elements is computed
// as `fe_degree + 1`. The reason for this is that the lowermost degree of the
// FE_Q finite elements equals 1 while the lowermost degree of the FE_Nedelec,
// FE_RaviartThomas, and FE_DGQ finite elements equals 0 by convention.
// The boundary ID is set in the circle.geo file. This is done by specifying the
// physical line, i.e.,
// @code
// Physical Line(1) = {21, 22, 40, 60, 79, 99, 119, 139};
// @endcode
// The boundary ID in the geo file must match the boundary ID setting below.
// If `project_exact_solution = true`, the program projects the exact solutions
// for $T$, $\vec{J}_f$, $\vec{A}$, and $B$ onto the corresponding function
// spaces and saves the results into corresponding vtu files next to the
// numerical solutions. If there are no bugs in the program, the projected exact
// solution and the numerical solution for $T$, $\vec{J}_f$, and $B$ will look
// alike. The situation with $\vec{A}$ is a bit more complicated, see section
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

  const double d1 = 0.1; // Half-side of the square in the center of the mesh.
  const double b1 = 0.3; // Radius of the magnetic core.
  const double a2 = 0.5; // Inner radius of the free-current region.
  const double b2 = 0.7; // Outer radius of the free-current region.

  const double K0 = 1.0; // Magnitude of the free-current density.

  const types::material_id material_id_free_space =
    1; // Material ID of the free space.
  const types::material_id material_id_core =
    2; // Material ID of the magnetic core.
  const types::material_id material_id_free_current =
    3; // Material ID of the free-current region.
  const types::boundary_id outer_boundary_id =
    1; // Boundary ID. Set in the geo file.

  const unsigned int mapping_degree = 2; // Mapping degree used in all solvers.
  const unsigned int fe_degree      = 0; // Degree of the finite elements.

  const double eta_squared = 0.0; // eta^2 when solving for A$.

  const unsigned int n_threads_max = 0; // If >0 limits the number of threads.
  const double eps = 1e-12; // Two doubles are equal if their difference < eps.
  const bool   log_cg_convergence = false; // Save the CG convergence data.
  const bool print_time_tables = false; // Print the time tables on the screen.
  const bool project_exact_solution = false; // Save the exact solution.
} // namespace Settings

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
  // current region.
  inline Tensor<1, 2> volume_free_current_density(const Point<2> &p,
                                                  const double    K0)
  {
    return Tensor<1, 2>({-K0 * p[1], K0 * p[0]});
  }

  // This class implements the closed-form analytical expression for the
  // free-current density, $\vec{J}_f$, in the entire domain.
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

// @sect3{Solver - T}

// This name space contains all the code related to the computation of the
// current vector potential, $T$.
namespace SolverT
{
  // This class describes the free-current density, $\vec{J}_f$, on the
  // right-hand side of the div-grad equation, see equation (i) in the
  // boundary value problem for $T$. The free-current density is given as
  // a closed-form analytical expression by the definition of the problem.
  class FreeCurrentDensity
  {
  public:
    void value_list(const std::vector<Point<2>> &p, // Quadrature points.
                    const types::material_id   mid, // Material ID of the cell.
                    std::vector<Tensor<1, 2>> &values) const
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

  // This class implements the solver that minimizes the functional $F(T)$, see
  // the introduction. The mesh is loaded and refined in this class. All other
  // solvers use a reference to this mesh.
  class Solver
  {
  public:
    Solver() = delete;
    Solver(const unsigned int p, // Degree of the FE_Q finite elements.
           const unsigned int r, // The mesh refinement parameter.
           const unsigned int mapping_degree,
           const std::string &fname          = "data",
           const Function<2> *exact_solution = nullptr);

    void make_mesh();  /* Loads and refines the mesh. Assigns material IDs.
                          Attaches spherical manifold.*/
    void setup();      // Initializes dofs, vectors, matrices.
    void assemble();   // Assembles the system of linear equations.
    void solve();      // Solves the system of linear equations.
    void save() const; // Saves computed T into a vtu file.
    void compute_error_norms();        // Computes L^2 error norm.
    void project_exact_solution_fcn(); // Projects exact solution.
    void clear()
    { // Clears the memory for the next solver.
      system_matrix.clear();
      system_rhs.reinit(0);
    }
    void run(); /* Executes the last eight functions in the proper order
                   and measures the execution time for each function. */

    double get_L2_norm()
    {
      return L2_norm;
    };

    unsigned int get_n_cells() const
    {
      return triangulation.n_active_cells();
    }

    types::global_dof_index get_n_dofs() const
    {
      return dof_handler.n_dofs();
    }

    // These three get functions are used to channel the mesh and the solution
    // to the next solver.
    const Triangulation<2> &get_tria() const
    {
      return triangulation;
    }
    const DoFHandler<2> &get_dof_handler() const
    {
      return dof_handler;
    }
    const Vector<double> &get_solution() const
    {
      return solution;
    }

  private:
    // The following data members are typical for all deal.II simulations:
    // triangulation, finite elements, dof handlers, etc. The constraints
    // are used to enforce the Dirichlet boundary condition. The names of the
    // data members are self-explanatory.

    Triangulation<2> triangulation;

    const FE_Q<2>  fe;
    DoFHandler<2>  dof_handler;
    Vector<double> solution;

    SparseMatrix<double> system_matrix;
    Vector<double>       system_rhs;

    const Function<2> *exact_solution;

    Vector<double> projected_exact_solution;

    AffineConstraints<double> constraints;
    SparsityPattern           sparsity_pattern;

    SphericalManifold<2> sphere;

    const unsigned int refinement_parameter;
    const unsigned int mapping_degree;

    Vector<double> L2_per_cell;
    double         L2_norm;

    const std::string fname;
    TimerOutput       timer;

    // The program utilizes the WorkStream technology. The Step-9 tutorial
    // does a much better job of explaining the workings of WorkStream.
    // Reading the @ref workstream_paper "WorkStream paper" is recommended.
    // The following structures and functions are related to WorkStream.
    struct AssemblyScratchData
    {
      AssemblyScratchData(const FiniteElement<2> &fe,
                          const unsigned int      mapping_degree);

      AssemblyScratchData(const AssemblyScratchData &scratch_data);

      const FreeCurrentDensity Jf;

      MappingQ<2> mapping;
      FEValues<2> fe_values;

      const unsigned int dofs_per_cell;
      const unsigned int n_q_points;

      std::vector<Tensor<1, 2>> Jf_list;

      const FEValuesExtractors::Scalar se;
    };

    struct AssemblyCopyData
    {
      FullMatrix<double>                   cell_matrix;
      Vector<double>                       cell_rhs;
      std::vector<types::global_dof_index> local_dof_indices;
    };

    void system_matrix_local(
      const typename DoFHandler<2>::active_cell_iterator &cell,
      AssemblyScratchData                                &scratch_data,
      AssemblyCopyData                                   &copy_data);

    void copy_local_to_global(const AssemblyCopyData &copy_data);
  };

  Solver::Solver(const unsigned int p,
                 const unsigned int r,
                 const unsigned int mapping_degree,
                 const std::string &fname,
                 const Function<2> *exact_solution)
    : fe(p + 1)
    , exact_solution(exact_solution)
    , refinement_parameter(r)
    , mapping_degree(mapping_degree)
    , fname(fname)
    , timer(std::cout,
            (Settings::print_time_tables) ? TimerOutput::summary :
                                            TimerOutput::never,
            TimerOutput::cpu_and_wall_times_grouped)
  {
    Assert(exact_solution != nullptr,
           ExcMessage("The exact solution is missing."));
  }

  // The following function loads and refines the mesh, assigns material
  // IDs to all cells, and attaches the spherical manifold to the mesh.
  // The material IDs are assigned on the basis of the distance from the
  // center of a cell to the origin. The spherical manifold is attached
  // to a face if all vertices of the face are at the same distance from
  // the origin provided the cell is outside the square in the center of
  // the mesh, see mesh description in the introduction. Note that the
  // material ID and the manifold ID are inherited by child cells, see
  // the glossary entries on @ref GlossMaterialId "Material id" and
  // @ref GlossManifoldIndicator "Manifold indicator". For this reason,
  // we can assign the IDs first and then refine the mesh.
  void Solver::make_mesh()
  {
    GridIn<2> gridin;

    gridin.attach_triangulation(triangulation);
    std::ifstream ifs("circle.msh");
    gridin.read_msh(ifs);

    triangulation.reset_all_manifolds();

    for (auto cell : triangulation.active_cell_iterators())
      {
        cell->set_material_id(
          Settings::material_id_free_space); // The cell is in free space.

        if (cell->center().norm() < Settings::b1)
          cell->set_material_id(
            Settings::material_id_core); // The cell is inside the core.

        if ((cell->center().norm() > Settings::a2) &&
            (cell->center().norm() < Settings::b2))
          cell->set_material_id(
            Settings::material_id_free_current); /* The cell is inside the Jf
                                                    region. */

        for (unsigned int f = 0; f < cell->n_faces(); f++)
          {
            double dif_norm = 0.0;
            for (unsigned int v = 1; v < cell->face(f)->n_vertices(); v++)
              dif_norm += std::abs(cell->face(f)->vertex(0).norm() -
                                   cell->face(f)->vertex(v).norm());

            if ((dif_norm < Settings::eps) &&
                (cell->center().norm() > Settings::d1))
              cell->face(f)->set_all_manifold_ids(1);
          }
      }

    triangulation.set_manifold(1, sphere);

    if (refinement_parameter > 0)
      triangulation.refine_global(refinement_parameter);
  }

  // This function initializes the dofs, applies the Dirichlet boundary
  // condition, and initializes the vectors and matrices.
  void Solver::setup()
  {
    dof_handler.reinit(triangulation);
    dof_handler.distribute_dofs(fe);

    // The following segment of the code applies the homogeneous Dirichlet
    // boundary condition. As discussed in the introduction, the Dirichlet
    // boundary condition is an essential condition and must be enforced
    // by constraining the system matrix. This segment of the code does the
    // constraining.
    constraints.clear();

    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    VectorTools::interpolate_boundary_values(MappingQ<2>(mapping_degree),
                                             dof_handler,
                                             Settings::outer_boundary_id,
                                             Functions::ZeroFunction<2>(),
                                             constraints);

    constraints.close();

    // The rest of the function arranges the dofs in a sparsity pattern and
    // initializes the system matrices and vectors.
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
  // WorkStream going. The interesting part, i.e., the actual assembling of the
  // system matrix and the right-hand side, happens below in the function
  // `Solver::system_matrix_local`.
  void Solver::assemble()
  {
    WorkStream::run(dof_handler.begin_active(),
                    dof_handler.end(),
                    *this,
                    &Solver::system_matrix_local,
                    &Solver::copy_local_to_global,
                    AssemblyScratchData(fe, mapping_degree),
                    AssemblyCopyData());
  }

  // The following two constructors initialize scratch data from the input
  // parameters and from another object of the same type, i.e., a copy
  // constructor.
  Solver::AssemblyScratchData::AssemblyScratchData(
    const FiniteElement<2> &fe,
    const unsigned int      mapping_degree)
    : Jf()
    , mapping(mapping_degree)
    , fe_values(mapping,
                fe,
                QGauss<2>(fe.degree + 1),
                update_gradients | update_values | update_quadrature_points |
                  update_JxW_values)
    , dofs_per_cell(fe_values.dofs_per_cell)
    , n_q_points(fe_values.get_quadrature().size())
    , Jf_list(n_q_points, Tensor<1, 2>())
    , se(0)
  {}

  Solver::AssemblyScratchData::AssemblyScratchData(
    const AssemblyScratchData &scratch_data)
    : Jf()
    , mapping(scratch_data.mapping.get_degree())
    , fe_values(mapping,
                scratch_data.fe_values.get_fe(),
                scratch_data.fe_values.get_quadrature(),
                update_gradients | update_values | update_quadrature_points |
                  update_JxW_values)
    , dofs_per_cell(fe_values.dofs_per_cell)
    , n_q_points(fe_values.get_quadrature().size())
    , Jf_list(n_q_points, Tensor<1, 2>())
    , se(0)
  {}

  // This function assembles a fraction of the system matrix and the system
  // right-hand side related to a single cell. These fractions are
  // `copy_data.cell_matrix` and `copy_data.cell_rhs`. They are copied into
  // the system matrix, $A_{ij}$, and the right-hand side, $b_i$, by the
  // function `Solver::copy_local_to_global()`.
  void Solver::system_matrix_local(
    const typename DoFHandler<2>::active_cell_iterator &cell,
    AssemblyScratchData                                &scratch_data,
    AssemblyCopyData                                   &copy_data)
  {
    // First, we reinitialize the matrices and vectors related to the current
    // cell and compute the FE values.
    copy_data.cell_matrix.reinit(scratch_data.dofs_per_cell,
                                 scratch_data.dofs_per_cell);

    copy_data.cell_rhs.reinit(scratch_data.dofs_per_cell);

    copy_data.local_dof_indices.resize(scratch_data.dofs_per_cell);

    scratch_data.fe_values.reinit(cell);

    // Second, we compute the free-current density, $\vec{J}_f$, at the
    // quadrature points.
    scratch_data.Jf.value_list(scratch_data.fe_values.get_quadrature_points(),
                               cell->material_id(),
                               scratch_data.Jf_list);

    // Third, we compute the components of the cell matrix and cell right-hand
    // side. The labels of the integrals are the same as in the introduction to
    // this tutorial.
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
              (scratch_data.Jf_list[q_index][0] *
                 scratch_data.fe_values[scratch_data.se].gradient(i,
                                                                  q_index)[1] -
               scratch_data.Jf_list[q_index][1] *
                 scratch_data.fe_values[scratch_data.se].gradient(i,
                                                                  q_index)[0]) *
              scratch_data.fe_values.JxW(q_index); // J_f.(curlv N_i)dS.
          }
      }

    // Finally, we query the dof indices on the current cell and store them
    // in the copy data structure, so we know to which locations of the system
    // matrix and right-hand side the components of the cell matrix and
    // cell right-hand side must be copied.
    cell->get_dof_indices(copy_data.local_dof_indices);
  }

  // This function copies the components of a cell matrix and a cell right-hand
  // side into the system matrix, $A_{ij}$, and the system right-hand side,
  // $b_i$.
  void Solver::copy_local_to_global(const AssemblyCopyData &copy_data)
  {
    constraints.distribute_local_to_global(copy_data.cell_matrix,
                                           copy_data.cell_rhs,
                                           copy_data.local_dof_indices,
                                           system_matrix,
                                           system_rhs);
  }

  // The following two functions compute the error norms and project the exact
  // solution.
  void Solver::compute_error_norms()
  {
    VectorTools::integrate_difference(MappingQ<2>(mapping_degree),
                                      dof_handler,
                                      solution,
                                      *exact_solution,
                                      L2_per_cell,
                                      QGauss<2>(fe.degree + 4),
                                      VectorTools::L2_norm);

    L2_norm = VectorTools::compute_global_error(triangulation,
                                                L2_per_cell,
                                                VectorTools::L2_norm);
  }

  void Solver::project_exact_solution_fcn()
  {
    AffineConstraints<double> constraints_empty;
    constraints_empty.close();

    VectorTools::project(MappingQ<2>(mapping_degree),
                         dof_handler,
                         constraints_empty,
                         QGauss<2>(fe.degree + 1),
                         *exact_solution,
                         projected_exact_solution);
  }

  // This function solves the system of linear equations.
  // If `Settings::log_cg_convergence == true`, the convergence data is saved
  // into a file.
  // The stopping condition is
  // \f[ |\boldsymbol{b} - \boldsymbol{A}\boldsymbol{c}|
  // < 10^{-8} |\boldsymbol{b}|.
  // \f]
  // As soon as we use constraints, we must not forget to distribute them.
  void Solver::solve()
  {
    SolverControl control(system_rhs.size(),
                          1.0e-8 * system_rhs.l2_norm(),
                          false,
                          false);

    if (Settings::log_cg_convergence)
      control.enable_history_data();

    GrowingVectorMemory<Vector<double>> memory;
    SolverCG<Vector<double>>            cg(control, memory);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    constraints.distribute(solution);

    if (Settings::log_cg_convergence)
      {
        const std::vector<double> history_data = control.get_history_data();

        std::ofstream ofs(fname + "_cg_convergence.csv");

        for (unsigned int i = 1; i < history_data.size(); i++)
          ofs << i << ", " << history_data[i] << "\n";
      }
  }

  // This function saves the computed current vector potential into a vtu file.
  void Solver::save() const
  {
    DataOut<2> data_out;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "ScalarField");

    if (exact_solution)
      {
        data_out.add_data_vector(L2_per_cell, "L2norm");

        if (Settings::project_exact_solution)
          data_out.add_data_vector(projected_exact_solution,
                                   "ScalarFieldExact");
      }

    DataOutBase::VtkFlags flags;
    flags.write_higher_order_cells = true;
    data_out.set_flags(flags);

    const MappingQ<2> mapping(mapping_degree);

    data_out.build_patches(mapping,
                           fe.degree + 2,
                           DataOut<2>::CurvedCellRegion::curved_inner_cells);

    std::ofstream ofs(fname + ".vtu");
    data_out.write_vtu(ofs);
  }

  void Solver::run()
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
    if (exact_solution)
      {
        {
          TimerOutput::Scope timer_section(timer, "Compute error norms");
          compute_error_norms();
        }

        if (Settings::project_exact_solution)
          {
            {
              TimerOutput::Scope timer_section(timer, "Project exact solution");
              project_exact_solution_fcn();
            }
          }
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
} // namespace SolverT

// @sect3{Solver - A}

// This name space contains all the code related to the computation of the
// magnetic vector potential, $\vec{A}$. The main difference between this
// solver and the solver for the current vector potential, $T$, is in how
// the information on the source is fed to respective solvers. The solver for
// $T$ is fed data sampled from the analytical closed-form expression for
// $\vec{J}_f$. The solver for $\vec{A}$ is fed a field function, i.e., a
// numerically computed current vector potential, $T$.
namespace SolverA
{
  // This class describes the permeability in the entire problem domain. The
  // permeability is given by the definition of the problem, see the
  // introduction.
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

  // This class implements the solver that minimizes the functional
  // $F(\vec{A})$. The numerically computed current vector potential, $T$,
  // is fed to this solver by means of the input parameters `dof_handler_T` and
  // `solution_T`. Moreover, this solver reuses the mesh on which $T$ has
  // been computed. The reference to the mesh is passed via the input parameter
  // `triangulation_T`.
  class Solver
  {
  public:
    Solver() = delete;
    Solver(const unsigned int      p, // Degree of the Nedelec finite elements.
           const unsigned int      mapping_degree,
           const Triangulation<2> &triangulation_T,
           const DoFHandler<2>    &dof_handler_T,
           const Vector<double>   &solution_T,
           const double            eta_squared    = 0.0,
           const std::string      &fname          = "data",
           const Function<2>      *exact_solution = nullptr);

    double get_L2_norm()
    {
      return L2_norm;
    };

    unsigned int get_n_cells() const
    {
      return triangulation_T.n_active_cells();
    }

    types::global_dof_index get_n_dofs() const
    {
      return dof_handler.n_dofs();
    }

    void setup();               // Initializes dofs, vectors, matrices.
    void assemble();            // Assembles the system of linear equations.
    void solve();               // Solves the system of linear equations.
    void save() const;          // Saves computed A into a vtu file.
    void compute_error_norms(); // Computes L^2 error norm.
    void project_exact_solution_fcn(); // Projects exact solution.
    void clear()                       // Clears the memory for the next solver.
    {
      system_matrix.clear();
      system_rhs.reinit(0);
    }
    void run(); /* Executes the last seven functions in the proper order
                   and measures the execution time for each function. */

    const DoFHandler<2> &get_dof_handler() const
    {
      return dof_handler;
    }
    const Vector<double> &get_solution() const
    {
      return solution;
    }

  private:
    const Triangulation<2> &triangulation_T;
    const DoFHandler<2>    &dof_handler_T;
    const Vector<double>   &solution_T;

    // The following data members are typical for all deal.II simulations:
    // triangulation, finite elements, dof handlers, etc. The names of the
    // data members are self-explanatory.
    const FE_Nedelec<2> fe;
    DoFHandler<2>       dof_handler;
    Vector<double>      solution;

    SparseMatrix<double> system_matrix;
    Vector<double>       system_rhs;

    Vector<double> projected_exact_solution;

    AffineConstraints<double> constraints;
    SparsityPattern           sparsity_pattern;

    const Function<2> *exact_solution;

    const unsigned int mapping_degree;
    const double       eta_squared;

    Vector<double> L2_per_cell;
    double         L2_norm;

    const std::string fname;
    TimerOutput       timer;

    // This time we have two dof handlers, `dof_handler_T` for $T$ and
    // `dof_handler` for $\vec{A}$. The WorkStream needs to walk through
    // the two dof handlers synchronously. For this purpose we will pair two
    // active cell iterators (one from `dof_handler_T`, another from
    // `dof_handler`). For that we need the `IteratorPair` type.
    using IteratorTuple =
      std::tuple<typename DoFHandler<2>::active_cell_iterator,
                 typename DoFHandler<2>::active_cell_iterator>;

    using IteratorPair = SynchronousIterators<IteratorTuple>;

    // The program utilizes the WorkStream technology. The Step-9 tutorial
    // does a much better job of explaining the workings of WorkStream.
    // Reading the @ref workstream_paper "WorkStream paper" is recommended.
    // The following structures and functions are related to WorkStream.
    struct AssemblyScratchData
    {
      AssemblyScratchData(const FiniteElement<2> &fe,
                          const DoFHandler<2>    &dof_hand_T,
                          const Vector<double>   &dofs_T,
                          const unsigned int      mapping_degree,
                          const double            eta_squared);

      AssemblyScratchData(const AssemblyScratchData &scratch_data);

      const Permeability permeability;

      MappingQ<2> mapping;
      FEValues<2> fe_values;

      FEValues<2> fe_values_T;

      const unsigned int dofs_per_cell;
      const unsigned int n_q_points;

      std::vector<double> permeability_list;
      std::vector<double> T_values;

      const FEValuesExtractors::Vector ve;

      const DoFHandler<2>  &dof_hand_T;
      const Vector<double> &dofs_T;

      const double eta_squared;
    };

    struct AssemblyCopyData
    {
      FullMatrix<double>                   cell_matrix;
      Vector<double>                       cell_rhs;
      std::vector<types::global_dof_index> local_dof_indices;
    };

    void system_matrix_local(const IteratorPair  &IP,
                             AssemblyScratchData &scratch_data,
                             AssemblyCopyData    &copy_data);

    void copy_local_to_global(const AssemblyCopyData &copy_data);
  };

  Solver::Solver(const unsigned int      p,
                 const unsigned int      mapping_degree,
                 const Triangulation<2> &triangulation_T,
                 const DoFHandler<2>    &dof_handler_T,
                 const Vector<double>   &solution_T,
                 const double            eta_squared,
                 const std::string      &fname,
                 const Function<2>      *exact_solution)
    : triangulation_T(triangulation_T)
    , dof_handler_T(dof_handler_T)
    , solution_T(solution_T)
    , fe(p)
    , exact_solution(exact_solution)
    , mapping_degree(mapping_degree)
    , eta_squared(eta_squared)
    , fname(fname)
    , timer(std::cout,
            (Settings::print_time_tables) ? TimerOutput::summary :
                                            TimerOutput::never,
            TimerOutput::cpu_and_wall_times_grouped)
  {}

  // This function initializes the dofs, vectors and matrices. We use the
  // homogeneous Neumann boundary condition. It is a natural boundary
  // condition. It is enforced by minimization of the functional. We do not
  // apply any essential boundary conditions this time around. For this reason
  // we use empty constraints.
  void Solver::setup()
  {
    constraints.close();

    dof_handler.reinit(triangulation_T);
    dof_handler.distribute_dofs(fe);

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);

    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    if (Settings::project_exact_solution && exact_solution)
      projected_exact_solution.reinit(dof_handler.n_dofs());

    if (exact_solution)
      L2_per_cell.reinit(triangulation_T.n_active_cells());
  }

  // Formally, this function assembles the system of linear equations. In
  // reality, however, it just spells all the magic words to get the WorkStream
  // going. The interesting part, i.e., the actual assembling of the system
  // matrix and the right-hand side happens below in the
  // `Solver::system_matrix_local` function. Note that this time the first two
  // input parameters to `WorkStream::run` are pairs of iterators, not
  // iterators themselves as per usual. Note also the order in which we package
  // the iterators: first the iterator of `dof_handler`, then the iterator of
  // the `dof_handler_T`. We will extract them in the same order.
  void Solver::assemble()
  {
    WorkStream::run(
      IteratorPair(IteratorTuple(dof_handler.begin_active(),
                                 dof_handler_T.begin_active())),
      IteratorPair(IteratorTuple(dof_handler.end(), dof_handler_T.end())),
      *this,
      &Solver::system_matrix_local,
      &Solver::copy_local_to_global,
      AssemblyScratchData(
        fe, dof_handler_T, solution_T, mapping_degree, eta_squared),
      AssemblyCopyData());
  }

  // The following two constructors initialize scratch data from the input
  // parameters and from another object of the same type, i.e., a copy
  // constructor.
  Solver::AssemblyScratchData::AssemblyScratchData(
    const FiniteElement<2> &fe,
    const DoFHandler<2>    &dof_hand_T,
    const Vector<double>   &dofs_T,
    const unsigned int      mapping_degree,
    const double            eta_squared)
    : permeability()
    , mapping(mapping_degree)
    , fe_values(mapping,
                fe,
                QGauss<2>(fe.degree + 2),
                update_gradients | update_values | update_JxW_values)
    , fe_values_T(mapping,
                  dof_hand_T.get_fe(),
                  QGauss<2>(fe.degree + 2),
                  update_values | update_gradients)
    , dofs_per_cell(fe_values.dofs_per_cell)
    , n_q_points(fe_values.get_quadrature().size())
    , permeability_list(n_q_points)
    , T_values(n_q_points)
    , ve(0)
    , dof_hand_T(dof_hand_T)
    , dofs_T(dofs_T)
    , eta_squared(eta_squared)
  {}

  Solver::AssemblyScratchData::AssemblyScratchData(
    const AssemblyScratchData &scratch_data)
    : permeability()
    , mapping(scratch_data.mapping.get_degree())
    , fe_values(mapping,
                scratch_data.fe_values.get_fe(),
                scratch_data.fe_values.get_quadrature(),
                update_gradients | update_values | update_JxW_values)
    , fe_values_T(mapping,
                  scratch_data.fe_values_T.get_fe(),
                  scratch_data.fe_values_T.get_quadrature(),
                  update_values | update_gradients)
    , dofs_per_cell(fe_values.dofs_per_cell)
    , n_q_points(fe_values.get_quadrature().size())
    , permeability_list(n_q_points)
    , T_values(n_q_points)
    , ve(0)
    , dof_hand_T(scratch_data.dof_hand_T)
    , dofs_T(scratch_data.dofs_T)
    , eta_squared(scratch_data.eta_squared)
  {}

  // This function assembles a fraction of the system matrix and the system
  // right-hand side related to a single cell. These fractions are
  // `copy_data.cell_matrix` and `copy_data.cell_rhs`. They are copied into
  // to the system matrix, $A_{ij}$, and the right-hand side, $b_i$, by the
  // function `Solver::copy_local_to_global()`.
  void Solver::system_matrix_local(const IteratorPair  &IP,
                                   AssemblyScratchData &scratch_data,
                                   AssemblyCopyData    &copy_data)
  {
    // First we reinitialize the matrices and vectors related to the current
    // cell.
    copy_data.cell_matrix.reinit(scratch_data.dofs_per_cell,
                                 scratch_data.dofs_per_cell);

    copy_data.cell_rhs.reinit(scratch_data.dofs_per_cell);

    copy_data.local_dof_indices.resize(scratch_data.dofs_per_cell);

    // Second, we extract the cells from the pair. We extract them in the
    // correct order, see above.
    auto cell   = std::get<0>(*IP);
    auto cell_T = std::get<1>(*IP);

    // Third, we compute the ordered FE values, the permeability, and the values
    // of the current vector potential, $T$, on the cell.
    scratch_data.fe_values.reinit(cell);
    scratch_data.fe_values_T.reinit(cell_T);

    scratch_data.permeability.value_list(cell->material_id(),
                                         scratch_data.permeability_list);

    scratch_data.fe_values_T.get_function_values(scratch_data.dofs_T,
                                                 scratch_data.T_values);

    // Fourth, we compute the components of the cell matrix and cell right-hand
    // side. The labels of the integrals are the same as in the introduction to
    // this tutorial.
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
              (scratch_data.T_values[q_index] *
               scratch_data.fe_values[scratch_data.ve].curl(i, q_index))[0] *
              scratch_data.fe_values.JxW(q_index); // T(curls N_i)dS
          }
      }

    // Finally, we query the dof indices on the current cell and store them
    // in the copy data structure, so we know to which locations of the system
    // matrix and right-hand side the components of the cell matrix and
    // cell right-hand side must be copied.
    cell->get_dof_indices(copy_data.local_dof_indices);
  }

  // This function copies the components of a cell matrix and a cell right-hand
  // side into the system matrix, $A_{ij}$, and the system right-hand side,
  // $b_i$.
  void Solver::copy_local_to_global(const AssemblyCopyData &copy_data)
  {
    constraints.distribute_local_to_global(copy_data.cell_matrix,
                                           copy_data.cell_rhs,
                                           copy_data.local_dof_indices,
                                           system_matrix,
                                           system_rhs);
  }

  // The following two functions compute the error norms and project the exact
  // solution.
  void Solver::compute_error_norms()
  {
    VectorTools::integrate_difference(MappingQ<2>(mapping_degree),
                                      dof_handler,
                                      solution,
                                      *exact_solution,
                                      L2_per_cell,
                                      QGauss<2>(fe.degree + 4),
                                      VectorTools::L2_norm);

    L2_norm = VectorTools::compute_global_error(triangulation_T,
                                                L2_per_cell,
                                                VectorTools::L2_norm);
  }

  void Solver::project_exact_solution_fcn()
  {
    AffineConstraints<double> constraints_empty;
    constraints_empty.close();

    VectorTools::project(MappingQ<2>(mapping_degree),
                         dof_handler,
                         constraints_empty,
                         QGauss<2>(fe.degree + 2),
                         *exact_solution,
                         projected_exact_solution);
  }

  // This function solves the system of linear equations.
  // If `Settings::log_cg_convergence == true`, the convergence data is saved
  // into a file. In theory, a CG solver can solve an $m \times m$ system of
  // linear equations in at most $m$ steps. In practice, it can take more steps
  // to converge. The convergence of the algorithm depends on the spectral
  // properties of the system matrix. The best case is if the eigenvalues form a
  // compact cluster away from zero. In our case, however, the eigenvalues are
  // spread in between zero and the maximal eigenvalue. Consequently, we expect
  // a poor convergence and increase the maximal number of iteration steps by a
  // factor of 10, i.e., `10*system_rhs.size()`. The stopping condition is \f[
  // |\boldsymbol{b} - \boldsymbol{A}\boldsymbol{c}|
  // < 10^{-8} |\boldsymbol{b}|.
  // \f]
  // As soon as we use constraints (they are empty this time, thought), we must
  // not forget to distribute them.
  void Solver::solve()
  {
    SolverControl control(10 * system_rhs.size(),
                          1.0e-8 * system_rhs.l2_norm(),
                          false,
                          false);

    if (Settings::log_cg_convergence)
      control.enable_history_data();

    GrowingVectorMemory<Vector<double>> memory;
    SolverCG<Vector<double>>            cg(control, memory);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    constraints.distribute(solution);

    if (Settings::log_cg_convergence)
      {
        const std::vector<double> history_data = control.get_history_data();

        std::ofstream ofs(fname + "_cg_convergence.csv");

        for (unsigned int i = 1; i < history_data.size(); i++)
          ofs << i << ", " << history_data[i] << "\n";
      }
  }

  // This function saves the computed magnetic vector potential into a vtu file.
  void Solver::save() const
  {
    std::vector<std::string> solution_names(2, "VectorField");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      interpretation(2,
                     DataComponentInterpretation::component_is_part_of_vector);

    DataOut<2> data_out;

    data_out.add_data_vector(dof_handler,
                             solution,
                             solution_names,
                             interpretation);
    if (Settings::project_exact_solution)
      {
        std::vector<std::string> solution_names_ex(2, "VectorFieldExact");

        data_out.add_data_vector(dof_handler,
                                 projected_exact_solution,
                                 solution_names_ex,
                                 interpretation);
      }

    if (exact_solution)
      {
        data_out.add_data_vector(L2_per_cell, "L2norm");
      }

    DataOutBase::VtkFlags flags;
    flags.write_higher_order_cells = true;
    data_out.set_flags(flags);

    const MappingQ<2> mapping(mapping_degree);

    data_out.build_patches(mapping,
                           fe.degree + 2,
                           DataOut<2>::CurvedCellRegion::curved_inner_cells);

    std::ofstream ofs(fname + ".vtu");
    data_out.write_vtu(ofs);
  }

  void Solver::run()
  {
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
    if (exact_solution)
      {
        {
          TimerOutput::Scope timer_section(timer, "Compute error norms");
          compute_error_norms();
        }

        if (Settings::project_exact_solution)
          {
            {
              TimerOutput::Scope timer_section(timer, "Project exact solution");
              project_exact_solution_fcn();
            }
          }
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
} // namespace SolverA

// @sect3{Projector from H(grad) to H(div)}
// This name space contains all the code related to the conversion of the
// current vector potential, $T$, into free-current density, $\vec{J}_f$.
// The current vector potential is modeled by the FE_Q finite elements,
// while the free-current density is modeled by the FE_RaviartThomas finite
// elements.
namespace ProjectorHgradToHdiv
{
  // This class implements the solver that minimizes the functional
  // $F(\vec{J}_f)$, see the introduction. The input out-of-plane vector field,
  // $T$, is fed to the solver by means of the input parameters
  // `dof_handler_Hgrad` and `solution_Hgrad`. Moreover, this solver reuses the
  // mesh on which $T$ has been computed. The reference to the mesh is passed
  // via the input parameter `triangulation_Hgrad`.
  class Solver
  {
  public:
    Solver() = delete;
    Solver(const unsigned int p, // Degree of the Raviart-Thomas finite elements
           const unsigned int mapping_degree,
           const Triangulation<2> &triangulation_Hgrad,
           const DoFHandler<2>    &dof_handler_Hgrad,
           const Vector<double>   &solution_Hgrad,
           const std::string      &fname          = "data",
           const Function<2>      *exact_solution = nullptr);

    double get_L2_norm()
    {
      return L2_norm;
    };

    unsigned int get_n_cells() const
    {
      return triangulation_Hgrad.n_active_cells();
    }

    types::global_dof_index get_n_dofs() const
    {
      return dof_handler_Hdiv.n_dofs();
    }

    void setup();               // Initializes dofs, vectors, matrices.
    void assemble();            // Assembles the system of linear equations.
    void solve();               // Solves the system of linear equations.
    void save() const;          // Saves computed Jf into a vtu file.
    void compute_error_norms(); // Computes L^2 error norm.
    void project_exact_solution_fcn(); // Projects exact solution.
    void clear()
    {
      system_matrix.clear();
      system_rhs.reinit(0);
    }
    void run(); /* Executes the last seven functions in the proper order
                   and measures the execution time for each function. */

  private:
    const Triangulation<2> &triangulation_Hgrad;
    const DoFHandler<2>    &dof_handler_Hgrad;
    const Vector<double>   &solution_Hgrad;

    // The following data members are typical for all deal.II simulations:
    // triangulation, finite elements, dof handlers, etc. The constraints
    // are used to enforce the Dirichlet boundary conditions. The names of the
    // data members are self-explanatory.
    const FE_RaviartThomas<2> fe_Hdiv;
    DoFHandler<2>             dof_handler_Hdiv;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution_Hdiv;
    Vector<double> system_rhs;

    Vector<double> projected_exact_solution;

    AffineConstraints<double> constraints;

    const Function<2> *exact_solution;

    const unsigned int mapping_degree;

    Vector<double> L2_per_cell;
    double         L2_norm;

    const std::string fname;
    TimerOutput       timer;

    // This time we have two dof handlers, `dof_handler_Hgrad` for $T$ and
    // `dof_handler_Hdiv` for $\vec{J}_f$. The WorkStream needs to walk through
    // the two dof handlers synchronously. For this purpose we will pair two
    // active cells iterators (one from `dof_handler_Hgrad`, another from
    // `dof_handler_Hdiv`) to be walked through synchronously. For that we
    // need the `IteratorPair` type.
    using IteratorTuple =
      std::tuple<typename DoFHandler<2>::active_cell_iterator,
                 typename DoFHandler<2>::active_cell_iterator>;

    using IteratorPair = SynchronousIterators<IteratorTuple>;

    // The program utilizes the WorkStream technology. The Step-9 tutorial
    // does a much better job of explaining the workings of WorkStream.
    // Reading the @ref workstream_paper "WorkStream paper" is recommended.
    // The following structures and functions are related to WorkStream.
    struct AssemblyScratchData
    {
      AssemblyScratchData(const FiniteElement<2> &fe,
                          const DoFHandler<2>    &dof_handr_Hgrad,
                          const Vector<double>   &dofs_Hgrad,
                          const unsigned int      mapping_degree);

      AssemblyScratchData(const AssemblyScratchData &scratch_data);

      MappingQ<2> mapping;
      FEValues<2> fe_values_Hdiv;
      FEValues<2> fe_values_Hgrad;

      const unsigned int dofs_per_cell;
      const unsigned int n_q_points;

      std::vector<Tensor<1, 2>> grad_T; // We will convert grad T into curlv T

      const FEValuesExtractors::Vector ve;

      const DoFHandler<2>  &dof_hand_Hgrad;
      const Vector<double> &dofs_Hgrad;
    };

    struct AssemblyCopyData
    {
      FullMatrix<double>                   cell_matrix;
      Vector<double>                       cell_rhs;
      std::vector<types::global_dof_index> local_dof_indices;
    };

    void system_matrix_local(const IteratorPair  &IP,
                             AssemblyScratchData &scratch_data,
                             AssemblyCopyData    &copy_data);

    void copy_local_to_global(const AssemblyCopyData &copy_data);
  };

  Solver::Solver(const unsigned int      p,
                 const unsigned int      mapping_degree,
                 const Triangulation<2> &triangulation_Hgrad,
                 const DoFHandler<2>    &dof_handler_Hgrad,
                 const Vector<double>   &solution_Hgrad,
                 const std::string      &fname,
                 const Function<2>      *exact_solution)
    : triangulation_Hgrad(triangulation_Hgrad)
    , dof_handler_Hgrad(dof_handler_Hgrad)
    , solution_Hgrad(solution_Hgrad)
    , fe_Hdiv(p)
    , exact_solution(exact_solution)
    , mapping_degree(mapping_degree)
    , fname(fname)
    , timer(std::cout,
            (Settings::print_time_tables) ? TimerOutput::summary :
                                            TimerOutput::never,
            TimerOutput::cpu_and_wall_times_grouped)
  {
    Assert(exact_solution != nullptr,
           ExcMessage("The exact solution is missing."));
  }

  // This function initializes the dofs, vectors and matrices. This time the
  // constraints are empty as we do not apply Dirichlet boundary condition.
  void Solver::setup()
  {
    constraints.close();

    dof_handler_Hdiv.reinit(triangulation_Hgrad);
    dof_handler_Hdiv.distribute_dofs(fe_Hdiv);

    DynamicSparsityPattern dsp(dof_handler_Hdiv.n_dofs(),
                               dof_handler_Hdiv.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler_Hdiv, dsp, constraints, false);

    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
    solution_Hdiv.reinit(dof_handler_Hdiv.n_dofs());
    system_rhs.reinit(dof_handler_Hdiv.n_dofs());

    if (Settings::project_exact_solution && exact_solution)
      projected_exact_solution.reinit(dof_handler_Hdiv.n_dofs());

    if (exact_solution)
      L2_per_cell.reinit(triangulation_Hgrad.n_active_cells());
  }

  // Formally, this function assembles the system of linear equations. In
  // reality, however, it just spells all the magic words to get the WorkStream
  // going. The interesting part, i.e., the actual assembling of the system
  // matrix and the right-hand side happens below in the
  // `Solver::system_matrix_local` function.
  void Solver::assemble()
  {
    WorkStream::run(IteratorPair(
                      IteratorTuple(dof_handler_Hdiv.begin_active(),
                                    dof_handler_Hgrad.begin_active())),
                    IteratorPair(IteratorTuple(dof_handler_Hdiv.end(),
                                               dof_handler_Hgrad.end())),
                    *this,
                    &Solver::system_matrix_local,
                    &Solver::copy_local_to_global,
                    AssemblyScratchData(fe_Hdiv,
                                        dof_handler_Hgrad,
                                        solution_Hgrad,
                                        mapping_degree),
                    AssemblyCopyData());
  }

  // The following two constructors initialize scratch data from the input
  // parameters and from another object of the same type, i.e., a copy
  // constructor.
  Solver::AssemblyScratchData::AssemblyScratchData(
    const FiniteElement<2> &fe,
    const DoFHandler<2>    &dof_hand_Hgrad,
    const Vector<double>   &dofs_Hgrad,
    const unsigned int      mapping_degree)
    : mapping(mapping_degree)
    , fe_values_Hdiv(mapping,
                     fe,
                     QGauss<2>(fe.degree + 2),
                     update_values | update_JxW_values)
    , fe_values_Hgrad(mapping,
                      dof_hand_Hgrad.get_fe(),
                      QGauss<2>(fe.degree + 2),
                      update_gradients)
    , dofs_per_cell(fe_values_Hdiv.dofs_per_cell)
    , n_q_points(fe_values_Hdiv.get_quadrature().size())
    , grad_T(n_q_points)
    , ve(0)
    , dof_hand_Hgrad(dof_hand_Hgrad)
    , dofs_Hgrad(dofs_Hgrad)
  {}

  Solver::AssemblyScratchData::AssemblyScratchData(
    const AssemblyScratchData &scratch_data)
    : mapping(scratch_data.mapping.get_degree())
    , fe_values_Hdiv(mapping,
                     scratch_data.fe_values_Hdiv.get_fe(),
                     scratch_data.fe_values_Hdiv.get_quadrature(),
                     update_values | update_JxW_values)
    , fe_values_Hgrad(mapping,
                      scratch_data.fe_values_Hgrad.get_fe(),
                      scratch_data.fe_values_Hgrad.get_quadrature(),
                      update_gradients)
    , dofs_per_cell(fe_values_Hdiv.dofs_per_cell)
    , n_q_points(fe_values_Hdiv.get_quadrature().size())
    , grad_T(scratch_data.n_q_points)
    , ve(0)
    , dof_hand_Hgrad(scratch_data.dof_hand_Hgrad)
    , dofs_Hgrad(scratch_data.dofs_Hgrad)
  {}

  // This function assembles a fraction of the system matrix and the system
  // right-hand side related to a single cell. These fractions are
  // `copy_data.cell_matrix` and `copy_data.cell_rhs`. They are copied into
  // to the system matrix, $A_{ij}$, and the right-hand side, $b_i$, by the
  // function `Solver::copy_local_to_global`.
  void Solver::system_matrix_local(const IteratorPair  &IP,
                                   AssemblyScratchData &scratch_data,
                                   AssemblyCopyData    &copy_data)
  {
    // First we reinitialize the matrices and vectors related to the current
    // cell, update the FE values, and compute $\vec{\nabla}T$. We will convert
    // $\vec{\nabla}T$ into
    // $\big(\vec{\nabla}\overset{V}{\times}T\big)\cdot\vec{N}_i$ in the code
    // below.
    copy_data.cell_matrix.reinit(scratch_data.dofs_per_cell,
                                 scratch_data.dofs_per_cell);

    copy_data.cell_rhs.reinit(scratch_data.dofs_per_cell);

    copy_data.local_dof_indices.resize(scratch_data.dofs_per_cell);

    scratch_data.fe_values_Hdiv.reinit(std::get<0>(*IP));
    scratch_data.fe_values_Hgrad.reinit(std::get<1>(*IP));

    scratch_data.fe_values_Hgrad.get_function_gradients(scratch_data.dofs_Hgrad,
                                                        scratch_data.grad_T);

    // Second, we compute the components of the cell matrix and cell right-hand
    // side. The labels of the integrals are the same as in the introduction to
    // this tutorial.
    for (unsigned int q_index = 0; q_index < scratch_data.n_q_points; ++q_index)
      {
        for (unsigned int i = 0; i < scratch_data.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < scratch_data.dofs_per_cell; ++j)
              {
                copy_data.cell_matrix(i, j) += // Integral I_a
                  scratch_data.fe_values_Hdiv[scratch_data.ve].value(i,
                                                                     q_index) *
                  scratch_data.fe_values_Hdiv[scratch_data.ve].value(j,
                                                                     q_index) *
                  scratch_data.fe_values_Hdiv.JxW(q_index); // Ni.Nj dS
              }

            copy_data.cell_rhs(i) += // Integral I_b
              (scratch_data.grad_T[q_index][1] *
                 scratch_data.fe_values_Hdiv[scratch_data.ve].value(
                   i, q_index)[0] -
               scratch_data.grad_T[q_index][0] *
                 scratch_data.fe_values_Hdiv[scratch_data.ve].value(
                   i, q_index)[1]) *
              scratch_data.fe_values_Hdiv.JxW(q_index); // (curlv T).Ni
          }
      }

    // Finally, we query the dof indices on the current cell and store them
    // in the copy data structure, so we know to which locations of the system
    // matrix and right-hand side the components of the cell matrix and
    // cell right-hand side must be copied.
    std::get<0>(*IP)->get_dof_indices(copy_data.local_dof_indices);
  }

  // This function copies the components of a cell matrix and a cell right-hand
  // side into the system matrix, $A_{ij}$, and the system right-hand side,
  // $b_i$.
  void Solver::copy_local_to_global(const AssemblyCopyData &copy_data)
  {
    constraints.distribute_local_to_global(copy_data.cell_matrix,
                                           copy_data.cell_rhs,
                                           copy_data.local_dof_indices,
                                           system_matrix,
                                           system_rhs);
  }

  // The following two functions compute the error norms and project the exact
  // solution.
  void Solver::compute_error_norms()
  {
    VectorTools::integrate_difference(MappingQ<2>(mapping_degree),
                                      dof_handler_Hdiv,
                                      solution_Hdiv,
                                      *exact_solution,
                                      L2_per_cell,
                                      QGauss<2>(fe_Hdiv.degree + 4),
                                      VectorTools::L2_norm);

    L2_norm = VectorTools::compute_global_error(triangulation_Hgrad,
                                                L2_per_cell,
                                                VectorTools::L2_norm);
  }

  void Solver::project_exact_solution_fcn()
  {
    AffineConstraints<double> constraints_empty;
    constraints_empty.close();

    VectorTools::project(MappingQ<2>(mapping_degree),
                         dof_handler_Hdiv,
                         constraints_empty,
                         QGauss<2>(fe_Hdiv.degree + 2),
                         *exact_solution,
                         projected_exact_solution);
  }

  // This function solves the system of linear equations. This time we are
  // dealing with a mass matrix. It has good spectral properties. Consequently,
  // we do not use the factor of 10 as in `SolverA`. The stopping condition is
  // \f[
  // |\boldsymbol{b} - \boldsymbol{A}\boldsymbol{c}|
  // < 10^{-8} |\boldsymbol{b}|.
  // \f]
  void Solver::solve()
  {
    SolverControl control(system_rhs.size(),
                          1.0e-8 * system_rhs.l2_norm(),
                          false,
                          false);

    if (Settings::log_cg_convergence)
      control.enable_history_data();

    GrowingVectorMemory<Vector<double>> memory;
    SolverCG<Vector<double>>            cg(control, memory);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve(system_matrix, solution_Hdiv, system_rhs, preconditioner);

    if (Settings::log_cg_convergence)
      {
        const std::vector<double> history_data = control.get_history_data();

        std::ofstream ofs(fname + "_cg_convergence.csv");

        for (unsigned int i = 1; i < history_data.size(); i++)
          ofs << i << ", " << history_data[i] << "\n";
      }
  }

  // This function saves the computed $\vec{J}_f$, $L^2$ error norm, and
  // the projected exact solution into a vtu file. The exact solution is
  // only saved if `Settings::project_exact_solution = true`
  void Solver::save() const
  {
    std::vector<std::string> solution_names(2, "VectorField");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      interpretation(2,
                     DataComponentInterpretation::component_is_part_of_vector);

    DataOut<2> data_out;

    data_out.add_data_vector(dof_handler_Hdiv,
                             solution_Hdiv,
                             solution_names,
                             interpretation);

    if (Settings::project_exact_solution)
      {
        std::vector<std::string> solution_names_ex(2, "VectorFieldExact");

        data_out.add_data_vector(dof_handler_Hdiv,
                                 projected_exact_solution,
                                 solution_names_ex,
                                 interpretation);
      }

    if (exact_solution)
      {
        data_out.add_data_vector(L2_per_cell, "L2norm");
      }

    DataOutBase::VtkFlags flags;
    flags.write_higher_order_cells = true;
    data_out.set_flags(flags);

    const MappingQ<2> mapping(mapping_degree);

    data_out.build_patches(mapping,
                           fe_Hdiv.degree + 2,
                           DataOut<2>::CurvedCellRegion::curved_inner_cells);

    std::ofstream ofs(fname + ".vtu");
    data_out.write_vtu(ofs);
  }

  void Solver::run()
  {
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

    if (exact_solution)
      {
        {
          TimerOutput::Scope timer_section(timer, "Compute error norms");
          compute_error_norms();
        }

        if (Settings::project_exact_solution)
          {
            {
              TimerOutput::Scope timer_section(timer, "Project exact solution");
              project_exact_solution_fcn();
            }
          }
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
} // namespace ProjectorHgradToHdiv

// @sect3{Projector from H(curl) to L2}
// This name space contains all the code related to the conversion of the
// magnetic vector potential, $\vec{A}$, into magnetic field, $B$.
// The magnetic vector potential is modeled by the FE_Nedelec finite elements,
// while the magnetic field is modeled by the FE_DGQ finite elements.
namespace ProjectorHcurlToL2
{
  // This class implements the solver that minimizes the functional $F(B)$,
  // see the introduction. The input vector field, $\vec{A}$, is fed to the
  // solver by means of the input parameters `dof_handler_Hcurl` and
  // `solution_Hcurl`. Moreover, this solver reuses the mesh on which $\vec{A}$
  // has been computed. The reference to the mesh is passed via the input
  // parameter `triangulation_Hcurl`. The constraints are empty this time
  // around as we are not going to enforce any essential conditions.
  class Solver
  {
  public:
    Solver() = delete;
    Solver(const unsigned int      p, // Degree of the FE_DGQ finite elements
           const unsigned int      mapping_degree,
           const Triangulation<2> &triangulation_Hcurl,
           const DoFHandler<2>    &dof_handler_Hcurl,
           const Vector<double>   &solution_Hcurl,
           const std::string      &fname          = "data",
           const Function<2>      *exact_solution = nullptr);

    double get_L2_norm()
    {
      return L2_norm;
    };

    unsigned int get_n_cells() const
    {
      return triangulation_Hcurl.n_active_cells();
    }

    types::global_dof_index get_n_dofs() const
    {
      return dof_handler_L2.n_dofs();
    }

    void setup();               // Initializes dofs, vectors, matrices.
    void assemble();            // Assembles the system of linear equations.
    void solve();               // Solves the system of linear equations.
    void save() const;          // Saves computed B into a vtu file.
    void compute_error_norms(); // Computes L^2 error norm.
    void project_exact_solution_fcn(); // Projects exact solution.
    void clear()
    {
      system_matrix.clear();
      system_rhs.reinit(0);
    }
    void run(); /* Executes the last seven functions in the proper order
                   and measures the execution time for each function. */

  private:
    const Triangulation<2> &triangulation_Hcurl;
    const DoFHandler<2>    &dof_handler_Hcurl;
    const Vector<double>   &solution_Hcurl;

    // The following data members are typical for all deal.II simulations:
    // triangulation, finite elements, dof handlers, etc. The constraints
    // are used to enforce the Dirichlet boundary conditions. The names of the
    // data members are self-explanatory.
    const FE_DGQ<2> fe_L2;
    DoFHandler<2>   dof_handler_L2;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution_L2;
    Vector<double> system_rhs;

    Vector<double> projected_exact_solution;

    AffineConstraints<double> constraints;

    const Function<2> *exact_solution;

    const unsigned int mapping_degree;

    Vector<double> L2_per_cell;
    double         L2_norm;

    const std::string fname;
    TimerOutput       timer;

    // This time we have two dof handlers, `dof_handler_Hcurl` for $\vec{A}$
    // and `dof_handler_L2` for $B$. The WorkStream needs to walk through the
    // two dof handlers synchronously. For this purpose we will pair two active
    // cells iterators (one from `dof_handler_Hcurl`, another from
    // `dof_handler_L2`) to be walked through synchronously. For that we need
    // the `IteratorPair` type.
    using IteratorTuple =
      std::tuple<typename DoFHandler<2>::active_cell_iterator,
                 typename DoFHandler<2>::active_cell_iterator>;

    using IteratorPair = SynchronousIterators<IteratorTuple>;

    // The program utilizes the WorkStream technology. The Step-9 tutorial
    // does a much better job of explaining the workings of WorkStream.
    // Reading the @ref workstream_paper "WorkStream paper" is recommended.
    // The following structures and functions are related to WorkStream.
    struct AssemblyScratchData
    {
      AssemblyScratchData(const FiniteElement<2> &fe,
                          const DoFHandler<2>    &dof_handr_Hcurl,
                          const Vector<double>   &dofs_Hcurl,
                          const unsigned int      mapping_degree);

      AssemblyScratchData(const AssemblyScratchData &scratch_data);

      MappingQ<2> mapping;
      FEValues<2> fe_values_L2;
      FEValues<2> fe_values_Hcurl;

      const unsigned int dofs_per_cell;
      const unsigned int n_q_points;

      std::vector<std::vector<Tensor<1, 2>>>
        grad_A; // We will convert grad A into curls A

      const FEValuesExtractors::Scalar se;

      const DoFHandler<2>  &dof_hand_Hcurl;
      const Vector<double> &dofs_Hcurl;
    };

    struct AssemblyCopyData
    {
      FullMatrix<double>                   cell_matrix;
      Vector<double>                       cell_rhs;
      std::vector<types::global_dof_index> local_dof_indices;
    };

    void system_matrix_local(const IteratorPair  &IP,
                             AssemblyScratchData &scratch_data,
                             AssemblyCopyData    &copy_data);

    void copy_local_to_global(const AssemblyCopyData &copy_data);
  };

  Solver::Solver(const unsigned int      p,
                 const unsigned int      mapping_degree,
                 const Triangulation<2> &triangulation_Hcurl,
                 const DoFHandler<2>    &dof_handler_Hcurl,
                 const Vector<double>   &solution_Hcurl,
                 const std::string      &fname,
                 const Function<2>      *exact_solution)
    : triangulation_Hcurl(triangulation_Hcurl)
    , dof_handler_Hcurl(dof_handler_Hcurl)
    , solution_Hcurl(solution_Hcurl)
    , fe_L2(p)
    , exact_solution(exact_solution)
    , mapping_degree(mapping_degree)
    , fname(fname)
    , timer(std::cout,
            (Settings::print_time_tables) ? TimerOutput::summary :
                                            TimerOutput::never,
            TimerOutput::cpu_and_wall_times_grouped)
  {
    Assert(exact_solution != nullptr,
           ExcMessage("The exact solution is missing."));
  }

  // This function initializes the dofs, vectors and matrices. This time the
  // constraints are empty as we do not have essential conditions to enforce.
  void Solver::setup()
  {
    constraints.close();

    dof_handler_L2.reinit(triangulation_Hcurl);
    dof_handler_L2.distribute_dofs(fe_L2);

    DynamicSparsityPattern dsp(dof_handler_L2.n_dofs(),
                               dof_handler_L2.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler_L2, dsp, constraints, false);

    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
    solution_L2.reinit(dof_handler_L2.n_dofs());
    system_rhs.reinit(dof_handler_L2.n_dofs());

    if (Settings::project_exact_solution && exact_solution)
      projected_exact_solution.reinit(dof_handler_L2.n_dofs());

    if (exact_solution)
      L2_per_cell.reinit(triangulation_Hcurl.n_active_cells());
  }

  // Formally, this function assembles the system of linear equations. In
  // reality, however, it just spells all the magic words to get the WorkStream
  // going. The interesting part, i.e., the actual assembling of the system
  // matrix and the right-hand side happens below in the
  // `Solver::system_matrix_local` function.
  void Solver::assemble()
  {
    WorkStream::run(IteratorPair(
                      IteratorTuple(dof_handler_L2.begin_active(),
                                    dof_handler_Hcurl.begin_active())),
                    IteratorPair(IteratorTuple(dof_handler_L2.end(),
                                               dof_handler_Hcurl.end())),
                    *this,
                    &Solver::system_matrix_local,
                    &Solver::copy_local_to_global,
                    AssemblyScratchData(
                      fe_L2, dof_handler_Hcurl, solution_Hcurl, mapping_degree),
                    AssemblyCopyData());
  }

  // The following two constructors initialize scratch data from the input
  // parameters and from another object of the same type, i.e., a copy
  // constructor.
  Solver::AssemblyScratchData::AssemblyScratchData(
    const FiniteElement<2> &fe,
    const DoFHandler<2>    &dof_hand_Hcurl,
    const Vector<double>   &dofs_Hcurl,
    const unsigned int      mapping_degree)
    : mapping(mapping_degree)
    , fe_values_L2(mapping,
                   fe,
                   QGauss<2>(fe.degree + 2),
                   update_values | update_JxW_values)
    , fe_values_Hcurl(mapping,
                      dof_hand_Hcurl.get_fe(),
                      QGauss<2>(fe.degree + 2),
                      update_gradients)
    , dofs_per_cell(fe_values_L2.dofs_per_cell)
    , n_q_points(fe_values_L2.get_quadrature().size())
    , grad_A(n_q_points, std::vector<Tensor<1, 2>>(2))
    , se(0)
    , dof_hand_Hcurl(dof_hand_Hcurl)
    , dofs_Hcurl(dofs_Hcurl)
  {}

  Solver::AssemblyScratchData::AssemblyScratchData(
    const AssemblyScratchData &scratch_data)
    : mapping(scratch_data.mapping.get_degree())
    , fe_values_L2(mapping,
                   scratch_data.fe_values_L2.get_fe(),
                   scratch_data.fe_values_L2.get_quadrature(),
                   update_values | update_JxW_values)
    , fe_values_Hcurl(mapping,
                      scratch_data.fe_values_Hcurl.get_fe(),
                      scratch_data.fe_values_Hcurl.get_quadrature(),
                      update_gradients)
    , dofs_per_cell(fe_values_L2.dofs_per_cell)
    , n_q_points(fe_values_L2.get_quadrature().size())
    , grad_A(scratch_data.n_q_points, std::vector<Tensor<1, 2>>(2))
    , se(0)
    , dof_hand_Hcurl(scratch_data.dof_hand_Hcurl)
    , dofs_Hcurl(scratch_data.dofs_Hcurl)
  {}

  // This function assembles a fraction of the system matrix and the system
  // right-hand side related to a single cell. These fractions are
  // `copy_data.cell_matrix` and `copy_data.cell_rhs`. They are copied into
  // to the system matrix, $A_{ij}$, and the right-hand side, $b_i$, by the
  // function `Solver::copy_local_to_global()`.
  void Solver::system_matrix_local(const IteratorPair  &IP,
                                   AssemblyScratchData &scratch_data,
                                   AssemblyCopyData    &copy_data)
  {
    // First we reinitialize the matrices and vectors related to the current
    // cell, update the FE values, and compute gradients of the components of
    // $\vec{A}$. We will convert the gradients into
    // $\big(\vec{\nabla}\overset{S}{\times}\vec{A}\big)N_i$ in the code
    // below.
    copy_data.cell_matrix.reinit(scratch_data.dofs_per_cell,
                                 scratch_data.dofs_per_cell);

    copy_data.cell_rhs.reinit(scratch_data.dofs_per_cell);

    copy_data.local_dof_indices.resize(scratch_data.dofs_per_cell);

    scratch_data.fe_values_L2.reinit(std::get<0>(*IP));
    scratch_data.fe_values_Hcurl.reinit(std::get<1>(*IP));

    scratch_data.fe_values_Hcurl.get_function_gradients(scratch_data.dofs_Hcurl,
                                                        scratch_data.grad_A);

    // Second, we compute the components of the cell matrix and cell right-hand
    // side. The labels of the integrals are the same as in the introduction to
    // this tutorial.
    for (unsigned int q_index = 0; q_index < scratch_data.n_q_points; ++q_index)
      {
        for (unsigned int i = 0; i < scratch_data.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < scratch_data.dofs_per_cell; ++j)
              {
                copy_data.cell_matrix(i, j) += // Integral I_a
                  scratch_data.fe_values_L2[scratch_data.se].value(i, q_index) *
                  scratch_data.fe_values_L2[scratch_data.se].value(j, q_index) *
                  scratch_data.fe_values_L2.JxW(q_index); // Ni Nj dS
              }

            copy_data.cell_rhs(i) += // Integral I_b
              (scratch_data.grad_A[q_index][1][0] -
               scratch_data.grad_A[q_index][0][1]) *
              scratch_data.fe_values_L2[scratch_data.se].value(i, q_index) *
              scratch_data.fe_values_L2.JxW(q_index); // (curls A) Ni dS
          }
      }

    // Finally, we query the dof indices on the current cell and store them
    // in the copy data structure, so we know to which locations of the system
    // matrix and right-hand side the components of the cell matrix and
    // cell right-hand side must be copied.
    std::get<0>(*IP)->get_dof_indices(copy_data.local_dof_indices);
  }

  // This function copies the components of a cell matrix and a cell right-hand
  // side into the system matrix, $A_{ij}$, and the system right-hand side,
  // $b_i$.
  void Solver::copy_local_to_global(const AssemblyCopyData &copy_data)
  {
    constraints.distribute_local_to_global(copy_data.cell_matrix,
                                           copy_data.cell_rhs,
                                           copy_data.local_dof_indices,
                                           system_matrix,
                                           system_rhs);
  }

  // The following two functions compute the error norms and project the exact
  // solution.
  void Solver::compute_error_norms()
  {
    VectorTools::integrate_difference(MappingQ<2>(mapping_degree),
                                      dof_handler_L2,
                                      solution_L2,
                                      *exact_solution,
                                      L2_per_cell,
                                      QGauss<2>(fe_L2.degree + 4),
                                      VectorTools::L2_norm);

    L2_norm = VectorTools::compute_global_error(triangulation_Hcurl,
                                                L2_per_cell,
                                                VectorTools::L2_norm);
  }

  void Solver::project_exact_solution_fcn()
  {
    AffineConstraints<double> constraints_empty;
    constraints_empty.close();

    VectorTools::project(MappingQ<2>(mapping_degree),
                         dof_handler_L2,
                         constraints_empty,
                         QGauss<2>(fe_L2.degree + 2),
                         *exact_solution,
                         projected_exact_solution);
  }

  // This function solves the system of linear equations. This time we are
  // dealing with a mass matrix. It has good spectral properties. Consequently,
  // we do not use the factor of 10 as in `SolverA`. The stopping condition is
  // \f[
  // |\boldsymbol{b} - \boldsymbol{A}\boldsymbol{c}|
  // < 10^{-8} |\boldsymbol{b}|.
  // \f]
  void Solver::solve()
  {
    SolverControl control(system_rhs.size(),
                          1.0e-8 * system_rhs.l2_norm(),
                          false,
                          false);

    if (Settings::log_cg_convergence)
      control.enable_history_data();

    GrowingVectorMemory<Vector<double>> memory;
    SolverCG<Vector<double>>            cg(control, memory);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve(system_matrix, solution_L2, system_rhs, preconditioner);

    if (Settings::log_cg_convergence)
      {
        const std::vector<double> history_data = control.get_history_data();

        std::ofstream ofs(fname + "_cg_convergence.csv");

        for (unsigned int i = 1; i < history_data.size(); i++)
          ofs << i << ", " << history_data[i] << "\n";
      }
  }

  // This function saves the computed $B$, $L^2$ error norm, and the projected
  // exact solution into a vtu file. The exact solution is only saved if
  // `Settings::project_exact_solution = true`.
  void Solver::save() const
  {
    std::vector<std::string> solution_names(1, "ScalarField");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        1, DataComponentInterpretation::component_is_scalar);

    DataOut<2> data_out;
    data_out.attach_dof_handler(dof_handler_L2);
    data_out.add_data_vector(solution_L2,
                             solution_names,
                             DataOut<2>::type_dof_data,
                             data_component_interpretation);

    if (exact_solution)
      {
        solution_names[0] = "L2norm";

        data_out.add_data_vector(L2_per_cell,
                                 solution_names,
                                 DataOut<2>::type_cell_data,
                                 data_component_interpretation);

        if (Settings::project_exact_solution)
          {
            solution_names[0] = "ScalarFieldExact";

            data_out.add_data_vector(projected_exact_solution,
                                     solution_names,
                                     DataOut<2>::type_dof_data,
                                     data_component_interpretation);
          }
      }

    DataOutBase::VtkFlags flags;
    flags.write_higher_order_cells = true;
    data_out.set_flags(flags);

    const MappingQ<2> mapping(mapping_degree);

    data_out.build_patches(mapping,
                           fe_L2.degree + 2,
                           DataOut<2>::CurvedCellRegion::curved_inner_cells);

    std::ofstream ofs(fname + ".vtu");
    data_out.write_vtu(ofs);
  }

  void Solver::run()
  {
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

    if (exact_solution)
      {
        {
          TimerOutput::Scope timer_section(timer, "Compute error norms");
          compute_error_norms();
        }

        if (Settings::project_exact_solution)
          {
            {
              TimerOutput::Scope timer_section(timer, "Project exact solution");
              project_exact_solution_fcn();
            }
          }
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
} // namespace ProjectorHcurlToL2


// @sect3{The main loop}

// This class contains the main loop of the program.
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

    for (unsigned int r = 0; r < 4; r++) // Mesh refinement parameter.
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

        SolverT::Solver T(Settings::fe_degree,
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

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

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_raviart_thomas.h>
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

#define VE scratch_data.ve

#define TMR(__name) TimerOutput::Scope timer_section(timer, __name)

using namespace dealii;

// This enumeration is used to switch between boundary conditions when solving
// for the magnetic vector potential, $\vec{A}$. The boundary condition is fixed
// (Dirichlet) when solving for the current vector potential, $\vec{T}$.
enum BoundaryConditionType
{
  Dirichlet,
  Neumann,
  Robin
};

// @sect3{Settings}

// The following is the control panel of the program. The scaling of the
// program can be changed by setting `mu_0 = 1.0`. Then the computed magnetic
// field, $\vec{B}$, will have to be multiplied by a factor of
// $\mu_0 = 1.25664 \cdot 10^{-6}$. The free-current density, $\vec{J}_f$, does
// not depend on scaling. The setting `boundary_condition_type_A` can be used
// to switch between Dirichlet, Neumann, and Robin (ABC) boundary conditions.
// This has an effect only on the solver that solves for the magnetic vector
// potential, $\vec{A}$. The solver that solves for the current vector
// potential,
// $\vec{T}$, applies the Dirichlet boundary condition. This cannot be changed.
// Recall, that forcing to zero the tangential component of $\vec{T}$ on the
// boundary allows us to skip the boundary integral $I_{b3-2}$. The boundary ID
// is set in the geo files that describe the mesh geometry. This is done by
// specifying the physical surface, i.e.,
// @code
// Physical Surface(2) = {115, 479, 665, 297, 847, 1029};
// @endcode
// The boundary ID in the geo file must match the boundary ID setting below.
// If `project_exact_solution = true`, the program projects the exact
// solutions for $\vec{J}_f$ and $\vec{B}$ onto $H(\text{div})$ function space
// and saves the results into corresponding vtu files next to the numerical
// solutions. If the projected exact solution and the numerical solution look
// alike, it is a good sign. In case of problems with the convergence of the
// CG solver the parameter $\eta^2$ must be increased just a bit.
class Settings
{
public:
  Settings(){};

  const double permeability_fs = 1.2566370614359172954e-6;
  const double mu_0            = permeability_fs; // Permeability of free space.

  const double mu_r = 4; // Relative permeability of the magnetic material.
  const double mu_1 = mu_0 * mu_r; // Permeability of the magnetic material.
  const double d1   = 0.1; // Half-side of the cube in the center of the mesh.
  const double a1   = 0.3; // Inner radius of the magnetic core.
  const double b1   = 0.6; // Outer radius of the magnetic core.
  const double a2   = 0.9; // Inner radius of the free-current region.
  const double b2   = 1.2; // Outer radius of the free-current region.
  const double d2   = 2.0; // Error norms are computed if r < d2.
  const double d3 = 3.0; // Radius of the outer boundary of the problem domain.
  const double K0 = 1.0; // Magnitude of the free-current density.
  const double H0 =      // Magnitude of H-field inside the coil, core removed.
    (1.0 / 3.0) * K0 * (pow(b2, 2) - pow(a2, 2));
  const types::material_id mid_1 = 1; // Material ID of the free space.
  const types::material_id mid_2 = 2; // Material ID of the magnetic core.
  const types::material_id mid_3 = 3; // Material ID of the free-current region.
  const types::boundary_id bid   = 2; // Boundary ID. Set in the geo file.

  const unsigned int          fe_degree = 0; // Degree of the finite elements.
  const BoundaryConditionType boundary_condition_type_A = Robin;

  const double eta_squared_T = 0.0; // eta^2 when solving for T$.
  const double eta_squared_A = 0.0; // eta^2 when solving for A$.

  const unsigned int nr_threads_max = 0; // If >0 limits the number of threads.
  const double eps = 1e-12; // Two doubles are equal if their difference < eps.
  const bool   log_cg_convergence = false; // Save the CG convergence data.
  const bool print_time_tables = false; // Print the time tables on the screen.
  const bool project_exact_solution = false; // Save the exact solution.
};

// Next comes the weight function used to limit the region in which the $L^2$
// error norms are computed. The error norms are computed if $r < d2$.
class Weight : public Function<3>, private Settings
{
public:
  virtual double value(const Point<3> &p,
                       const unsigned int) const override final
  {
    if (p.norm() > d2)
      return 0.0;

    return 1.0;
  }
};

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
  void append_new_order(const std::string &new_clmn)
  {
    new_order.push_back(new_clmn);
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
  inline Tensor<1, 3>
  volume_free_current_density(double x, double y, double z, double K0)
  {
    Tensor<1, 3> Jf;

    Jf[0] = -K0 * y;
    Jf[1] = K0 * x;
    Jf[2] = 0.0 * z;

    return Jf;
  }

  // This function computes the magnetic field induced by the free current in
  // absence of the magnetic core, see the expression for $\vec{B}_J$ in the
  // introduction.
  inline Tensor<1, 3> magnetic_field_coil(double x,
                                          double y,
                                          double z,
                                          double K0,
                                          double mu_0,
                                          double a2,
                                          double b2)
  {
    const double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

    const double cos_theta = z / r;
    const double sin_theta = sqrt(pow(x, 2) + pow(y, 2)) / r;

    const double cos_phi = x / sqrt(pow(x, 2) + pow(y, 2));
    const double sin_phi = y / sqrt(pow(x, 2) + pow(y, 2));

    const Tensor<1, 3> r_hat     = {x / r, y / r, z / r};
    const Tensor<1, 3> theta_hat = {cos_theta * cos_phi,
                                    cos_theta * sin_phi,
                                    -sin_theta};

    F1 = (2.0 / 3.0) * mu_0 * K0 * (cos_theta * r_hat - sin_theta * theta_hat);
    F2 = (2.0 / 3.0) * mu_0 * K0 *
         (cos_theta * r_hat + 0.5 * sin_theta * theta_hat) / pow(r, 3);

    if (r <= a2)
      {
        Bj = 0.5 * (pow(b2, 2) - pow(a2, 2)) * F1;
      }
    else if (r >= b2)
      {
        Bj = 0.2 * (pow(b2, 5) - pow(a2, 5)) * F2;
      }
    else
      {
        Bj = 0.5 * (pow(b2, 2) - pow(r, 2)) * F1 +
             0.2 * (pow(r, 5) - pow(a2, 5)) * F2;
      }

    return Bj;
  }

  // This function computes the magnetic field induced by the magnetic core, see
  // the expression for $\vec{B}_{\mu}$ in the introduction.
  inline Tensor<1, 3> magnetic_field_core(double x,
                                          double y,
                                          double z,
                                          double H0,
                                          double mur,
                                          double mu0,
                                          double a1,
                                          double b1)
  {
    double a3 = std::pow(a1, 3);
    double b3 = std::pow(b1, 3);

    double OMEGA   = ((mur - 1.0) / (mur + 2.0)) * (a3 / b3);
    double gamma_1 = (-3.0 * b3 * H0 * OMEGA) /
                     ((2.0 * mur + 1.0) - 2.0 * (mur - 1.0) * OMEGA);
    double beta_1  = ((2.0 * mur + 1.0) * gamma_1) / ((mur - 1.0) * a3);
    double alpha_1 = (-b3 * H0 + 2.0 * mur * gamma_1 - mur * b3 * beta_1) / 2.0;
    double delta_1 = (mur * a3 * beta_1 - 2.0 * mur * gamma_1) / a3;

    double r  = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    double r3 = std::pow(r, 3);
    double r5 = std::pow(r, 5);

    double zz = -3.0 * z * z / r5;
    double xz = -3.0 * x * z / r5;
    double yz = -3.0 * y * z / r5;

    double mu;

    Tensor<1, 3> grad_psi;
    Tensor<1, 3> Bmu;
    Tensor<1, 3> Happl;

    if (r <= a1)
      {
        grad_psi[0] = 0.0;
        grad_psi[1] = 0.0;
        grad_psi[2] = delta_1;
        mu          = mu0;
      }
    else if (r >= b1)
      {
        grad_psi[0] = alpha_1 * xz;
        grad_psi[1] = alpha_1 * yz;
        grad_psi[2] = -H0 + alpha_1 / r3 + alpha_1 * zz;
        mu          = mu0;
      }
    else
      {
        grad_psi[0] = gamma_1 * xz;
        grad_psi[1] = gamma_1 * yz;
        grad_psi[2] = beta_1 + gamma_1 / r3 + gamma_1 * zz;
        mu          = mur * mu0;
      }

    Happl[0] = 0.0;
    Happl[1] = 0.0;
    Happl[2] = H0;

    Bmu = -mu * grad_psi - mu0 * Happl;

    return Bmu;
  }

  // This class implements the closed-form analytical expression for the
  // free-current density, $\vec{J}_f$, in the entire domain.
  class FreeCurrentDensity : public Function<3>, public Settings
  {
  public:
    FreeCurrentDensity()
      : Function<3>(3)
    {}

    virtual void
    vector_value_list(const std::vector<Point<3>> &p,
                      std::vector<Vector<double>> &values) const override final
    {
      Assert(values.size() == p.size(),
             ExcDimensionMismatch(values.size(), p.size()));

      Tensor<1, 3> Jf;
      double       r;

      for (unsigned int i = 0; i < values.size(); i++)
        {
          r = p[i].norm();

          if ((r >= a2) && (r <= b2))
            {
              Jf = volume_free_current_density(p[i][0], p[i][1], p[i][2], K0);

              values[i][0] = Jf[0];
              values[i][1] = Jf[1];
              values[i][2] = Jf[2];
            }
          else
            {
              values[i][0] = 0.0;
              values[i][1] = 0.0;
              values[i][2] = 0.0;
            }
        }
    }
  };

  // This class implements the closed-form analytical expression for the
  // magnetic field induced by the coil, $\vec{B} = \vec{B}_J + \vec{B}_{\mu}$,
  // see the introduction.
  class MagneticField : public Function<3>, public Settings
  {
  public:
    MagneticField()
      : Function<3>(3)
    {}

    virtual void
    vector_value_list(const std::vector<Point<3>> &p,
                      std::vector<Vector<double>> &values) const override final
    {
      Assert(values.size() == p.size(),
             ExcDimensionMismatch(values.size(), p.size()));

      Tensor<1, 3> B;

      for (unsigned int i = 0; i < values.size(); i++)
        {
          B = magnetic_field_coil(p[i][0], p[i][1], p[i][2], K0, mu_0, a2, b2) +
              magnetic_field_core(
                p[i][0], p[i][1], p[i][2], H0, mu_r, mu_0, a1, b1);

          values[i][0] = B[0];
          values[i][1] = B[1];
          values[i][2] = B[2];
        }
    }
  };

} // namespace ExactSolutions

// @sect3{Solver - T}

// This name space contains all the code related to the computation of the
// current vector potential, $\vec{T}$.
namespace SolverT
{
  // This class describes the free-current density, $\vec{J}_f$, on the
  // right-hand side of the curl-curl equation, see equation (i) in the
  // boundary value problem for $\vec{T}$. The free-current density is given as
  // a closed-form analytical expression by the definition of the problem.
  class FreeCurrentDensity : public Settings
  {
  public:
    void value_list(const std::vector<Point<3>> &p, // Quadrature points.
                    types::material_id         mid, // Material ID of the cell.
                    std::vector<Tensor<1, 3>> &values) const
    {
      Assert(p.size() == values.size(),
             ExcDimensionMismatch(p.size(), values.size()));

      if ((mid == mid_1) || (mid == mid_2))
        for (unsigned int i = 0; i < values.size(); i++)
          {
            values[i][0] = 0.0;
            values[i][1] = 0.0;
            values[i][2] = 0.0;
          }

      if (mid == mid_3)
        {
          Tensor<1, 3> Jf;

          for (unsigned int i = 0; i < values.size(); i++)
            {
              Jf = ExactSolutions::volume_free_current_density(p[i][0],
                                                               p[i][1],
                                                               p[i][2],
                                                               K0);

              values[i][0] = Jf[0];
              values[i][1] = Jf[1];
              values[i][2] = Jf[2];
            }
        }
    }
  };

  // This class implements the solver that minimizes the functional
  // $F(\vec{T})$, see the introduction. The mesh is loaded in this class.
  // All other solvers use a reference to this mesh.
  class Solver : public Settings
  {
  public:
    Solver() = delete;
    Solver(unsigned int p, // Degree of the Nedelec finite elements.
           unsigned int r, // The mesh refinement parameter.
           unsigned int mapping_degree,
           double       eta_squared = 0.0,
           std::string  fname       = "data");

    void make_mesh();  // Loads the mesh. Assigns mat. IDs. Attaches manifold.
    void setup();      // Initializes dofs, vectors, matrices.
    void assemble();   // Assembles the system of linear equations.
    void solve();      // Solves the system of linear equations.
    void save() const; // Saves computed T into a vtu file.
    void clear()
    { // Clears the memory for the next solver.
      system_matrix.clear();
      system_rhs.reinit(0);
    }

    // These three get-functions are used to channel the mesh and the solution
    // to the next solver.
    const Triangulation<3> &get_tria() const
    {
      return triangulation;
    }
    const DoFHandler<3> &get_dof_handler() const
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
    // are used to enforce the Dirichlet boundary conditions. The names of the
    // data members are self-explanatory.

    Triangulation<3> triangulation;

    const FE_Nedelec<3> fe;
    DoFHandler<3>       dof_handler;
    Vector<double>      solution;

    SparseMatrix<double> system_matrix;
    Vector<double>       system_rhs;

    AffineConstraints<double> constraints;
    SparsityPattern           sparsity_pattern;

    SphericalManifold<3> sphere;

    const unsigned int refinement_parameter;
    const unsigned int mapping_degree;
    const double       eta_squared;
    const std::string  fname;

    // The program utilizes the WorkStream technology. The Step-9 tutorial
    // does a much better job of explaining the workings of WarkStream.
    // Reading the "WorkStream paper", see the glossary, is recommended.
    // The following structures and functions are related to WorkStream.
    struct AssemblyScratchData
    {
      AssemblyScratchData(const FiniteElement<3> &fe,
                          double                  eta_squared,
                          unsigned int            mapping_degree);

      AssemblyScratchData(const AssemblyScratchData &scratch_data);

      FreeCurrentDensity Jf;

      MappingQ<3> mapping;
      FEValues<3> fe_values;

      const unsigned int dofs_per_cell;
      const unsigned int n_q_points;

      std::vector<Tensor<1, 3>> Jf_list;

      const double                     eta_squared;
      const FEValuesExtractors::Vector ve;
    };

    struct AssemblyCopyData
    {
      FullMatrix<double>                   cell_matrix;
      Vector<double>                       cell_rhs;
      std::vector<types::global_dof_index> local_dof_indices;
    };

    void system_matrix_local(
      const typename DoFHandler<3>::active_cell_iterator &cell,
      AssemblyScratchData                                &scratch_data,
      AssemblyCopyData                                   &copy_data);

    void copy_local_to_global(const AssemblyCopyData &copy_data);
  };

  // The following is the implementation of the solver's constructor.
  // After initializing the constant data members the constructor adjusts the
  // timer which is used to measure execution time. Next, it runs all the
  // functions that implement the functionality of the solver in the correct
  // order. The brackets around the functions such as `{ TMR("Save"); ...}` are
  // needed for a proper operation of the timer.
  Solver::Solver(unsigned int p,
                 unsigned int r,
                 unsigned int mapping_degree,
                 double       eta_squared,
                 std::string  fname)
    : fe(p)
    , refinement_parameter(r)
    , mapping_degree(mapping_degree)
    , eta_squared(eta_squared)
    , fname(fname)
  {
    TimerOutput::OutputFrequency tf =
      (print_time_tables) ? TimerOutput::summary : TimerOutput::never;

    TimerOutput timer(std::cout, tf, TimerOutput::cpu_and_wall_times_grouped);

    {
      TMR("Make mesh");
      make_mesh();
    }
    {
      TMR("Setup");
      setup();
    }
    {
      TMR("Assemble");
      assemble();
    }
    {
      TMR("Solve");
      solve();
    }
    {
      TMR("Save");
      save();
    }
  }

  // The following function loads the mesh, assigns material IDs to all cells,
  // and attaches the spherical manifold to the mesh. The material IDs are
  // assigned on the basis of the distance from the center of a cell to the
  // origin. The spherical manifold is attached to a face if all vertices of
  // the face are at the same distance from the origin provided the cell is
  // outside the cube in the center of the mesh, see mesh description in the
  // introduction.
  void Solver::make_mesh()
  {
    GridIn<3> gridin;

    gridin.attach_triangulation(Solver::triangulation);
    std::ifstream ifs("sphere_r" + std::to_string(refinement_parameter) +
                      ".msh");
    gridin.read_msh(ifs);

    Solver::triangulation.reset_all_manifolds();

    for (auto cell : Solver::triangulation.active_cell_iterators())
      {
        cell->set_material_id(mid_1); // The cell is in free space.

        if ((cell->center().norm() > a1) && (cell->center().norm() < b1))
          cell->set_material_id(mid_2); // The cell is inside the core.

        if ((cell->center().norm() > a2) && (cell->center().norm() < b2))
          cell->set_material_id(mid_3); // The cell is inside the Jf region.

        for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; f++)
          {
            double dif_norm = 0.0;
            for (unsigned int v = 1; v < GeometryInfo<3>::vertices_per_face;
                 v++)
              dif_norm += std::abs(cell->face(f)->vertex(0).norm() -
                                   cell->face(f)->vertex(v).norm());

            if ((dif_norm < eps) && (cell->center().norm() > d1))
              cell->face(f)->set_all_manifold_ids(1);
          }
      }

    Solver::triangulation.set_manifold(1, sphere);
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

    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler,
      0,
      Functions::ZeroFunction<3>(3),
      bid,
      constraints,
      MappingQ<3>(mapping_degree));

    constraints.close();

    // The rest of the function arranges the dofs in a sparsity pattern and
    // initializes the system matrices and vectors.
    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }

  // Formally, this function assembles the system of linear equations. In
  // reality, however, it just spells all the magic words to get the WorkStream
  // going. The interesting part, i.e., the actual assembling of the system
  // matrix and the right-hand side, happens below in the function
  // Solver::system_matrix_local.
  void Solver::assemble()
  {
    WorkStream::run(dof_handler.begin_active(),
                    dof_handler.end(),
                    *this,
                    &Solver::system_matrix_local,
                    &Solver::copy_local_to_global,
                    AssemblyScratchData(fe, eta_squared, mapping_degree),
                    AssemblyCopyData());
  }

  // The following two constructors initialize scratch data from the input
  // parameters and from another object of the same type, i.e., a copy
  // constructor.
  Solver::AssemblyScratchData::AssemblyScratchData(const FiniteElement<3> &fe,
                                                   double       eta_squared,
                                                   unsigned int mapping_degree)
    : mapping(mapping_degree)
    , fe_values(mapping,
                fe,
                QGauss<3>(fe.degree + 2),
                update_gradients | update_values | update_quadrature_points |
                  update_JxW_values)
    , dofs_per_cell(fe_values.dofs_per_cell)
    , n_q_points(fe_values.get_quadrature().size())
    , Jf_list(n_q_points, Tensor<1, 3>())
    , eta_squared(eta_squared)
    , ve(0)
  {}

  Solver::AssemblyScratchData::AssemblyScratchData(
    const AssemblyScratchData &scratch_data)
    : mapping(scratch_data.mapping.get_degree())
    , fe_values(mapping,
                scratch_data.fe_values.get_fe(),
                scratch_data.fe_values.get_quadrature(),
                update_gradients | update_values | update_quadrature_points |
                  update_JxW_values)
    , dofs_per_cell(fe_values.dofs_per_cell)
    , n_q_points(fe_values.get_quadrature().size())
    , Jf_list(n_q_points, Tensor<1, 3>())
    , eta_squared(scratch_data.eta_squared)
    , ve(0)
  {}

  // This function assembles a fraction of the system matrix and the system
  // right-hand side related to a single cell. These fractions are
  // `copy_data.cell_matrix` and `copy_data.cell_rhs`. They are copied into
  // to the system matrix, $A_{ij}$, and the right-hand side, $b_i$, by the
  // function `Solver::copy_local_to_global()`.
  void Solver::system_matrix_local(
    const typename DoFHandler<3>::active_cell_iterator &cell,
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
                copy_data.cell_matrix(i, j) += // Integral I_a1+I_a3.
                  (scratch_data.fe_values[VE].curl(i, q_index) * // curl N_i
                     scratch_data.fe_values[VE].curl(j, q_index) // curl N_j
                   + scratch_data.eta_squared *                  // eta^2
                       scratch_data.fe_values[VE].value(i, q_index) * // N_i
                       scratch_data.fe_values[VE].value(j, q_index)   // N_j
                   ) *
                  scratch_data.fe_values.JxW(q_index); // dV
              }
            copy_data.cell_rhs(i) += // Integral I_b3-1.
              (scratch_data.Jf_list[q_index][0] *
                 scratch_data.fe_values[VE].curl(i, q_index)[0] +
               scratch_data.Jf_list[q_index][1] *
                 scratch_data.fe_values[VE].curl(i, q_index)[1] +
               scratch_data.Jf_list[q_index][2] *
                 scratch_data.fe_values[VE].curl(i, q_index)[2]) *
              scratch_data.fe_values.JxW(q_index); // J_f.(curl N_i)dV.
          }
      }

    // Finally, we query the dof indices on the current cell and store them
    // in the copy data structure, so we know to which locations of the system
    // matrix and right-hand side the components of the cell matrix and
    // cell right-hand side must be copied.
    cell->get_dof_indices(copy_data.local_dof_indices);
  }

  // This function copies the components of a cell matrix and a cell right-hand
  // side into the system matrix, $A_{i,j}$, and the system right-hand side,
  // $b_i$.
  void Solver::copy_local_to_global(const AssemblyCopyData &copy_data)
  {
    constraints.distribute_local_to_global(copy_data.cell_matrix,
                                           copy_data.cell_rhs,
                                           copy_data.local_dof_indices,
                                           system_matrix,
                                           system_rhs);
  }

  // This function solves the system of linear equations.
  // If `Settings::log_cg_convergence == true`, the convergence data is saved
  // into a file. In theory, a CG solver can solve an $m \times m$ system of
  // linear equations in at most $m$ steps. In practice, it can take more steps
  // to converge. The convergence of the algorithm depends on the spectral
  // properties of the system matrix. The best if the eigenvalues form a compact
  // cluster away from zero. In our case, however, the eigenvalues are spread in
  // between zero and the maximal eigenvalue. Consequently, we expect a poor
  // convergence and increase the maximal number of iteration steps by a factor
  // of 1000, i.e., `1000*system_rhs.size()`. The stopping condition is
  // \f[
  // |\boldsymbol{b} - \boldsymbol{A}\boldsymbol{c}|
  // < 10^{-6} |\boldsymbol{b}|.
  // \f]
  // As soon as we use constraints, we must not forget to distribute them.
  void Solver::solve()
  {
    SolverControl control(1000 * system_rhs.size(),
                          1.0e-6 * system_rhs.l2_norm(),
                          false,
                          false);

    if (log_cg_convergence)
      control.enable_history_data();

    GrowingVectorMemory<Vector<double>> memory;
    SolverCG<Vector<double>>            cg(control, memory);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    constraints.distribute(solution);

    if (log_cg_convergence)
      {
        const std::vector<double> history_data = control.get_history_data();

        std::ofstream ofs(fname + "_cg_convergence.csv");

        unsigned int i = 1;
        for (auto item : history_data)
          {
            ofs << i << ", " << item << "\n";
            i++;
          }

        ofs.close();
      }
  }

  // This function saves the computed current vector potential into a vtu file.
  void Solver::save() const
  {
    std::vector<std::string> solution_names(3, "VectorField");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      interpretation(3,
                     DataComponentInterpretation::component_is_part_of_vector);

    DataOut<3> data_out;

    data_out.add_data_vector(dof_handler,
                             solution,
                             solution_names,
                             interpretation);

    std::ofstream ofs;

    DataOutBase::VtkFlags flags;
    flags.write_higher_order_cells = true;
    data_out.set_flags(flags);

    const MappingQ<3> mapping(mapping_degree);

    data_out.build_patches(mapping,
                           fe.degree + 2,
                           DataOut<3>::CurvedCellRegion::curved_inner_cells);

    ofs.open(fname + ".vtu");
    data_out.write_vtu(ofs);

    ofs.close();
  }
} // namespace SolverT

// @sect3{Solver - A}

// This name space contains all the code related to the computation of the
// magnetic vector potential, $\vec{A}$. The main difference between this
// solver and the solver for the current vector potential, $\vec{T}$, is in how
// the information on the source is fed to respective solvers. The solver for
// $\vec{T}$ is fed data sampled from the analytical closed-form expression for
// $\vec{J}_f$. The solver for $\vec{A}$ is fed a field function, i.e., a
// numerically computed current vector potential, $\vec{T}$.
namespace SolverA
{
  // This class describes the permeability in the entire problem domain. The
  // permeability is given by the definition of the problem, see the
  // introduction.
  class Permeability : private Settings
  {
  public:
    void value_list(types::material_id mid, std::vector<double> &values) const
    {
      if ((mid == mid_1) || (mid == mid_3))
        for (unsigned int i = 0; i < values.size(); i++)
          values[i] = mu_0;

      if (mid == mid_2)
        for (unsigned int i = 0; i < values.size(); i++)
          values[i] = mu_1;
    }
  };

  // This class describes the parameter $\gamma$ in the Robin boundary
  // condition. As soon as it is evaluated on the boundary, the permeability
  // equals to that of free space. Therefore, we evaluate the parameter gamma as
  // \f[
  // \gamma = \dfrac{1}{\mu_0 r}.
  // \f]
  class Gamma : private Settings
  {
  public:
    void value_list(const std::vector<Point<3>> &r,
                    std::vector<double>         &values) const
    {
      Assert(r.size() == values.size(),
             ExcDimensionMismatch(r.size(), values.size()));

      for (unsigned int i = 0; i < values.size(); i++)
        values[i] = 1.0 / (mu_0 * r[i].norm());
    }
  };

  // This class implements the solver that minimizes the functional
  // $F(\vec{A})$. The numerically computed current vector potential, $\vec{T}$,
  // is fed to this solver by means of the input parameters `dof_handler_T` and
  // `solution_T`. Moreover, this solver reuses the mesh on which $\vec{T}$ has
  // been computed. The reference to the mesh is passed via the input parameter
  // `triangulation_T`.
  class Solver : Settings
  {
  public:
    Solver() = delete;
    Solver(unsigned int            p, // Degree of the Nedelec finite elements.
           unsigned int            mapping_degree,
           const Triangulation<3> &triangulation_T,
           const DoFHandler<3>    &dof_handler_T,
           const Vector<double>   &solution_T,
           double                  eta_squared = 0.0,
           std::string             fname       = "data");

    void setup();      // Initializes dofs, vectors, matrices.
    void assemble();   // Assembles the system of linear equations.
    void solve();      // Solves the system of linear equations.
    void save() const; // Saves computed T into a vtu file.
    void clear()       // Clears the memory for the next solver.
    {
      system_matrix.clear();
      system_rhs.reinit(0);
    }

    const DoFHandler<3> &get_dof_handler() const
    {
      return dof_handler;
    }
    const Vector<double> &get_solution() const
    {
      return solution;
    }

  private:
    const Triangulation<3> &triangulation_T;
    const DoFHandler<3>    &dof_handler_T;
    const Vector<double>   &solution_T;

    // The following data members are typical for all deal.II simulations:
    // triangulation, finite elements, dof handlers, etc. The constraints
    // are used to enforce the Dirichlet boundary conditions. The names of the
    // data members are self-explanatory.
    const FE_Nedelec<3> fe;
    DoFHandler<3>       dof_handler;
    Vector<double>      solution;

    SparseMatrix<double> system_matrix;
    Vector<double>       system_rhs;

    AffineConstraints<double> constraints;
    SparsityPattern           sparsity_pattern;

    const unsigned int mapping_degree;
    const double       eta_squared;
    const std::string  fname;

    // This time we have two dof handlers, `dof_handler_T` for $\vec{T}$ and
    // `dof_handler` for $\vec{A}$. The WorkStream needs to walk through
    // the two dof handlers synchronously. For this purpose we will pair two
    // active cell iterators (one from `dof_handler_T`, another from
    // `dof_handler`). For that we need the `IteratorPair` type.
    using IteratorTuple =
      std::tuple<typename DoFHandler<3>::active_cell_iterator,
                 typename DoFHandler<3>::active_cell_iterator>;

    using IteratorPair = SynchronousIterators<IteratorTuple>;

    // The program utilizes the WorkStream technology. The Step-9 tutorial
    // does a much better job of explaining the workings of WorkStream.
    // Reading the "WorkStream paper", see the glossary, is recommended.
    // The following structures and functions are related to WorkStream.
    struct AssemblyScratchData
    {
      AssemblyScratchData(const FiniteElement<3> &fe,
                          const DoFHandler<3>    &dof_hand_T,
                          const Vector<double>   &dofs_T,
                          unsigned int            mapping_degree,
                          double                  eta_squared,
                          BoundaryConditionType   boundary_condition_type);

      AssemblyScratchData(const AssemblyScratchData &scratch_data);

      Permeability permeability;
      Gamma        gamma;

      MappingQ<3>     mapping;
      FEValues<3>     fe_values;
      FEFaceValues<3> fe_face_values;

      FEValues<3> fe_values_T;

      const unsigned int dofs_per_cell;
      const unsigned int n_q_points;
      const unsigned int n_q_points_face;

      std::vector<double>       permeability_list;
      std::vector<double>       gamma_list;
      std::vector<Tensor<1, 3>> T_values;

      const FEValuesExtractors::Vector ve;

      const DoFHandler<3>  &dof_hand_T;
      const Vector<double> &dofs_T;

      const double                eta_squared;
      const BoundaryConditionType boundary_condition_type;
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

  // The following is the implementation of the solver's constructor.
  // After initializing the constant data members the constructor adjusts the
  // timer which is used to measure execution time. Next, it runs all the
  // functions that implement the functionality of the solver in the correct
  // order. The brackets around the functions such as `{ TMR("Save"); ...}` are
  // needed for a proper operation of the timer.
  Solver::Solver(unsigned int            p,
                 unsigned int            mapping_degree,
                 const Triangulation<3> &triangulation_T,
                 const DoFHandler<3>    &dof_handler_T,
                 const Vector<double>   &solution_T,
                 double                  eta_squared,
                 std::string             fname)
    : triangulation_T(triangulation_T)
    , dof_handler_T(dof_handler_T)
    , solution_T(solution_T)
    , fe(p)
    , mapping_degree(mapping_degree)
    , eta_squared(eta_squared)
    , fname(fname)
  {
    TimerOutput::OutputFrequency tf =
      (print_time_tables) ? TimerOutput::summary : TimerOutput::never;

    TimerOutput timer(std::cout, tf, TimerOutput::cpu_and_wall_times_grouped);

    {
      TMR("Setup");
      setup();
    }
    {
      TMR("Assemble");
      assemble();
    }
    {
      TMR("Solve");
      solve();
    }
    {
      TMR("Save");
      save();
    }
  }

  // This function initializes the dofs, applies the Dirichlet boundary
  // condition, and initializes the vectors and matrices.
  void Solver::setup()
  {
    dof_handler.reinit(triangulation_T);
    dof_handler.distribute_dofs(fe);

    // The following segment of the code applies the homogeneous Dirichlet
    // boundary condition. As discussed in the introduction, the Dirichlet
    // boundary condition is an essential condition and must be enforced
    // by constraining the system matrix. This segment of code does the
    // constraining.
    constraints.clear();

    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    if (boundary_condition_type_A == Dirichlet)
      VectorTools::project_boundary_values_curl_conforming_l2(
        dof_handler,
        0,
        Functions::ZeroFunction<3>(3),
        bid,
        constraints,
        MappingQ<3>(mapping_degree));

    constraints.close();

    // The rest of the function arranges the dofs in a sparsity pattern and
    // initializes the system matrix and the system vectors.
    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);

    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }

  // Formally, this function assembles the system of linear equations. In
  // reality, however, it just spells all the magic words to get the WorkStream
  // going. The interesting part, i.e., the actual assembling of the system
  // matrix and the right-hand side happens below in the
  // Solver::system_matrix_local function. Note that this time the first two
  // input parameters to `WorkStream::run` are pairs of iterators, not
  // iterators themselves as per usual. Note also the order in which we package
  // the iterators: first the iterator of `dof_handler`, then the iterator of
  // the `dof_handler_T`. We will extract them in the same order.
  void Solver::assemble()
  {
    WorkStream::run(IteratorPair(IteratorTuple(dof_handler.begin_active(),
                                               dof_handler_T.begin_active())),
                    IteratorPair(
                      IteratorTuple(dof_handler.end(), dof_handler_T.end())),
                    *this,
                    &Solver::system_matrix_local,
                    &Solver::copy_local_to_global,
                    AssemblyScratchData(fe,
                                        dof_handler_T,
                                        solution_T,
                                        mapping_degree,
                                        eta_squared,
                                        boundary_condition_type_A),
                    AssemblyCopyData());
  }

  // The following two constructors initialize scratch data from the input
  // parameters and from another object of the same type, i.e., a copy
  // constructor.
  Solver::AssemblyScratchData::AssemblyScratchData(
    const FiniteElement<3> &fe,
    const DoFHandler<3>    &dof_hand_T,
    const Vector<double>   &dofs_T,
    unsigned int            mapping_degree,
    double                  eta_squared,
    BoundaryConditionType   boundary_condition_type)
    : mapping(mapping_degree)
    , fe_values(mapping,
                fe,
                QGauss<3>(fe.degree + 2),
                update_gradients | update_values | update_JxW_values)
    , fe_face_values(mapping,
                     fe,
                     QGauss<2>(fe.degree + 2),
                     update_values | update_normal_vectors |
                       update_quadrature_points | update_JxW_values)
    , fe_values_T(mapping,
                  dof_hand_T.get_fe(),
                  QGauss<3>(fe.degree + 2),
                  update_values | update_gradients)
    , dofs_per_cell(fe_values.dofs_per_cell)
    , n_q_points(fe_values.get_quadrature().size())
    , n_q_points_face(fe_face_values.get_quadrature().size())
    , permeability_list(n_q_points)
    , gamma_list(n_q_points_face)
    , T_values(n_q_points)
    , ve(0)
    , dof_hand_T(dof_hand_T)
    , dofs_T(dofs_T)
    , eta_squared(eta_squared)
    , boundary_condition_type(boundary_condition_type)
  {}

  Solver::AssemblyScratchData::AssemblyScratchData(
    const AssemblyScratchData &scratch_data)
    : mapping(scratch_data.mapping.get_degree())
    , fe_values(mapping,
                scratch_data.fe_values.get_fe(),
                scratch_data.fe_values.get_quadrature(),
                update_gradients | update_values | update_JxW_values)
    , fe_face_values(mapping,
                     scratch_data.fe_face_values.get_fe(),
                     scratch_data.fe_face_values.get_quadrature(),
                     update_values | update_normal_vectors |
                       update_quadrature_points | update_JxW_values)
    , fe_values_T(mapping,
                  scratch_data.fe_values_T.get_fe(),
                  scratch_data.fe_values_T.get_quadrature(),
                  update_values | update_gradients)
    , dofs_per_cell(fe_values.dofs_per_cell)
    , n_q_points(fe_values.get_quadrature().size())
    , n_q_points_face(fe_face_values.get_quadrature().size())
    , permeability_list(n_q_points)
    , gamma_list(n_q_points_face)
    , T_values(n_q_points)
    , ve(0)
    , dof_hand_T(scratch_data.dof_hand_T)
    , dofs_T(scratch_data.dofs_T)
    , eta_squared(scratch_data.eta_squared)
    , boundary_condition_type(scratch_data.boundary_condition_type)
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
    // of the current vector potential, $\vec{T}$, on the cell.
    scratch_data.fe_values.reinit(cell);
    scratch_data.fe_values_T.reinit(cell_T);

    scratch_data.permeability.value_list(cell->material_id(),
                                         scratch_data.permeability_list);

    scratch_data.fe_values_T[VE].get_function_values(scratch_data.dofs_T,
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
                  (scratch_data.fe_values[VE].curl(i, q_index) *    // curl N_i
                     scratch_data.fe_values[VE].curl(j, q_index)    // curl N_j
                   + scratch_data.eta_squared *                     // eta^2
                       scratch_data.fe_values[VE].value(i, q_index) * // N_i
                       scratch_data.fe_values[VE].value(j, q_index)   // N_j
                   ) *
                  scratch_data.fe_values.JxW(q_index); // dV
              }
            copy_data.cell_rhs(i) += // Integral I_b3-1.
              (scratch_data.T_values[q_index][0] *
                 scratch_data.fe_values[VE].curl(i, q_index)[0] +
               scratch_data.T_values[q_index][1] *
                 scratch_data.fe_values[VE].curl(i, q_index)[1] +
               scratch_data.T_values[q_index][2] *
                 scratch_data.fe_values[VE].curl(
                   i, q_index)[2]) *               // T . (curl N_i)
              scratch_data.fe_values.JxW(q_index); // dV
          }
      }

    // If the Robin boundary condition (first-order ABC) is ordered,
    // we compute an extra integral over the boundary.
    if (scratch_data.boundary_condition_type == BoundaryConditionType::Robin)
      {
        for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; ++f)
          {
            if (cell->face(f)->at_boundary())
              {
                scratch_data.fe_face_values.reinit(cell, f);

                for (unsigned int q_index_face = 0;
                     q_index_face < scratch_data.n_q_points_face;
                     ++q_index_face)
                  {
                    for (unsigned int i = 0; i < scratch_data.dofs_per_cell;
                         ++i)
                      {
                        scratch_data.gamma.value_list(
                          scratch_data.fe_face_values.get_quadrature_points(),
                          scratch_data.gamma_list);

                        for (unsigned int j = 0; j < scratch_data.dofs_per_cell;
                             ++j)
                          {
                            copy_data.cell_matrix(i, j) += // Integral I_a2.
                              scratch_data.gamma_list[q_index_face] * // gamma
                              ((scratch_data.fe_face_values.normal_vector(
                                  q_index_face)[1] *
                                  scratch_data.fe_face_values[VE].value(
                                    i, q_index_face)[2] -
                                scratch_data.fe_face_values.normal_vector(
                                  q_index_face)[2] *
                                  scratch_data.fe_face_values[VE].value(
                                    i, q_index_face)[1]) *
                                 (scratch_data.fe_face_values.normal_vector(
                                    q_index_face)[1] *
                                    scratch_data.fe_face_values[VE].value(
                                      j, q_index_face)[2] -
                                  scratch_data.fe_face_values.normal_vector(
                                    q_index_face)[2] *
                                    scratch_data.fe_face_values[VE].value(
                                      j, q_index_face)[1]) +
                               (scratch_data.fe_face_values.normal_vector(
                                  q_index_face)[0] *
                                  scratch_data.fe_face_values[VE].value(
                                    i, q_index_face)[2] -
                                scratch_data.fe_face_values.normal_vector(
                                  q_index_face)[2] *
                                  scratch_data.fe_face_values[VE].value(
                                    i, q_index_face)[0]) *
                                 (scratch_data.fe_face_values.normal_vector(
                                    q_index_face)[0] *
                                    scratch_data.fe_face_values[VE].value(
                                      j, q_index_face)[2] -
                                  scratch_data.fe_face_values.normal_vector(
                                    q_index_face)[2] *
                                    scratch_data.fe_face_values[VE].value(
                                      j, q_index_face)[0]) +
                               (scratch_data.fe_face_values.normal_vector(
                                  q_index_face)[0] *
                                  scratch_data.fe_face_values[VE].value(
                                    i, q_index_face)[1] -
                                scratch_data.fe_face_values.normal_vector(
                                  q_index_face)[1] *
                                  scratch_data.fe_face_values[VE].value(
                                    i, q_index_face)[0]) *
                                 (scratch_data.fe_face_values.normal_vector(
                                    q_index_face)[0] *
                                    scratch_data.fe_face_values[VE].value(
                                      j, q_index_face)[1] -
                                  scratch_data.fe_face_values.normal_vector(
                                    q_index_face)[1] *
                                    scratch_data.fe_face_values[VE].value(
                                      j, q_index_face)[0])) // (n x N_i) . (n x
                                                            // N_j)
                              * scratch_data.fe_face_values.JxW(
                                  q_index_face); // dS
                          }                      // for j
                      }                          // for i
                  }                              // for q_index_face
              } // if (cell->face(f)->at_boundary())
          }     // for f
      }         // if (scratch_data.boundary_condition_type ==
                // BoundaryConditionType::Robin)

    // Finally, we query the dof indices on the current cell and store them
    // in the copy data structure, so we know to which locations of the system
    // matrix and right-hand side the components of the cell matrix and
    // cell right-hand side must be copied.
    cell->get_dof_indices(copy_data.local_dof_indices);
  }

  // This function copies the components of a cell matrix and a cell right-hand
  // side into the system matrix, $A_{i,j}$, and the system right-hand side,
  // $b_i$.
  void Solver::copy_local_to_global(const AssemblyCopyData &copy_data)
  {
    constraints.distribute_local_to_global(copy_data.cell_matrix,
                                           copy_data.cell_rhs,
                                           copy_data.local_dof_indices,
                                           system_matrix,
                                           system_rhs);
  }

  // This function solves the system of linear equations.
  // If `Settings::log_cg_convergence == true`, the convergence data is saved
  // into a file. In theory, a CG solver can solve an $m \times m$ system of
  // linear equations in at most $m$ steps. In practice, it can take more steps
  // to converge. The convergence of the algorithm depends on the spectral
  // properties of the system matrix. The best if the eigenvalues form a compact
  // cluster away from zero. In our case, however, the eigenvalues are spread in
  // between zero and the maximal eigenvalue. Consequently, we expect a poor
  // convergence and increase the maximal number of iteration steps by a factor
  // of 1000, i.e., `1000*system_rhs.size()`. The stopping condition is
  // \f[
  // |\boldsymbol{b} - \boldsymbol{A}\boldsymbol{c}|
  // < 10^{-6} |\boldsymbol{b}|.
  // \f]
  // As soon as we use constraints, we must not forget to distribute them.
  void Solver::solve()
  {
    SolverControl control(1000 * system_rhs.size(),
                          1.0e-6 * system_rhs.l2_norm(),
                          false,
                          false);

    if (log_cg_convergence)
      control.enable_history_data();

    GrowingVectorMemory<Vector<double>> memory;
    SolverCG<Vector<double>>            cg(control, memory);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    constraints.distribute(solution);

    if (log_cg_convergence)
      {
        const std::vector<double> history_data = control.get_history_data();

        std::ofstream ofs(fname + "_cg_convergence.csv");

        unsigned int i = 1;
        for (auto item : history_data)
          {
            ofs << i << ", " << item << "\n";
            i++;
          }
        ofs.close();
      }
  }

  // This function saves the computed magnetic vector potential into a vtu file.
  void Solver::save() const
  {
    std::vector<std::string> solution_names(3, "VectorField");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      interpretation(3,
                     DataComponentInterpretation::component_is_part_of_vector);

    DataOut<3> data_out;

    data_out.add_data_vector(dof_handler,
                             solution,
                             solution_names,
                             interpretation);

    std::ofstream ofs;

    DataOutBase::VtkFlags flags;
    flags.write_higher_order_cells = true;
    data_out.set_flags(flags);

    const MappingQ<3> mapping(mapping_degree);

    data_out.build_patches(mapping,
                           fe.degree + 2,
                           DataOut<3>::CurvedCellRegion::curved_inner_cells);

    ofs.open(fname + ".vtu");
    data_out.write_vtu(ofs);

    ofs.close();
  }
} // namespace SolverA


// @sect3{Projector from H(curl) to H(div)}
// This name space contains all the code related to the conversion of the
// magnetic vector potential, $\vec{A}$, into magnetic field, $\vec{B}$.
// The magnetic vector potential is modeled by the FE_Nedelec finite elements,
// while the magnetic field is modeled by the FE_RaviartThomas finite elements.
// This code is also used for converting the current vector potential,
// $\vec{T}$ into the free-current density, $\vec{J}_f$.
namespace ProjectorHcurlToHdiv
{

  // This class implements the solver that minimizes the functional $F(\vec{B})$
  // or $F(\vec{J}_f)$, see the introduction. The input vector field,
  // $\vec{A}$ or $\vec{T}$, is fed to the solver by means of the input
  // parameters `dof_handler_Hcurl` and `solution_Hcurl`. Moreover, this solver
  // reuses the mesh on which the input vector field has been computed. The
  // reference to the mesh is passed via the input parameter
  // `triangulation_Hcurl`. There is no constraints this time around as we are
  // not going to apply the Dirichlet boundary condition.
  class Solver : public Settings
  {
  public:
    Solver() = delete;
    Solver(unsigned int p, // Degree of the Raviart-Thomas finite elements
           unsigned int mapping_degree,
           const Triangulation<3> &triangulation_Hcurl,
           const DoFHandler<3>    &dof_handler_Hcurl,
           const Vector<double>   &solution_Hcurl,
           std::string             fname          = "data",
           const Function<3>      *exact_solution = nullptr);

    double get_L2_norm()
    {
      return L2_norm;
    };

    unsigned int get_n_cells() const
    {
      return static_cast<unsigned int>(triangulation_Hcurl.n_active_cells());
    }

    unsigned int get_n_dofs() const
    {
      return static_cast<unsigned int>(dof_handler_Hdiv.n_dofs());
    }

    void clear()
    {
      system_matrix.clear();
      system_rhs.reinit(0);
    }

  private:
    void setup();               // Initializes dofs, vectors, matrices.
    void assemble();            // Assembles the system of linear equations.
    void solve();               // Solves the system of linear equations.
    void save() const;          // Saves computed T into a vtu file.
    void compute_error_norms(); // Computes L^2 error norm.
    void project_exact_solution_fcn(); // Projects exact solution.

    const std::string fname;

    const Triangulation<3> &triangulation_Hcurl;
    const DoFHandler<3>    &dof_handler_Hcurl;
    const Vector<double>   &solution_Hcurl;

    // The following data members are typical for all deal.II simulations:
    // triangulation, finite elements, dof handlers, etc. The constraints
    // are used to enforce the Dirichlet boundary conditions. The names of the
    // data members are self-explanatory.
    const FE_RaviartThomas<3> fe_Hdiv;
    DoFHandler<3>             dof_handler_Hdiv;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution_Hdiv;
    Vector<double> system_rhs;

    Vector<double> projected_exact_solution;

    AffineConstraints<double> constraints;

    const Function<3> *exact_solution;

    const unsigned int mapping_degree;

    Vector<double> L2_per_cell;
    double         L2_norm;

    // This time we have two dof handlers, `dof_handler_Hcurl` for the input
    // vector field and `dof_handler_Hdiv` for the output vector field. The
    // WorkStream needs to walk through the two dof handlers synchronously.
    // For this purpose we will pair two active cells iterators (one from
    // `dof_handler_Hcurl`, another from `dof_handler_Hdiv`) to be walked
    // through synchronously. For that we need the `IteratorPair` type.
    using IteratorTuple =
      std::tuple<typename DoFHandler<3>::active_cell_iterator,
                 typename DoFHandler<3>::active_cell_iterator>;

    using IteratorPair = SynchronousIterators<IteratorTuple>;

    // The program utilizes the WorkStream technology. The Step-9 tutorial
    // does a much better job of explaining the workings of WarkStream.
    // Reading the "WorkStream paper", see the glossary, is recommended.
    // The following structures and functions are related to WorkStream.
    struct AssemblyScratchData
    {
      AssemblyScratchData(const FiniteElement<3> &fe,
                          const DoFHandler<3>    &dof_handr_Hcurl,
                          const Vector<double>   &dofs_Hcurl,
                          unsigned int            mapping_degree);

      AssemblyScratchData(const AssemblyScratchData &scratch_data);

      MappingQ<3> mapping;
      FEValues<3> fe_values_Hdiv;
      FEValues<3> fe_values_Hcurl;

      const unsigned int dofs_per_cell;
      const unsigned int n_q_points;

      std::vector<std::vector<Tensor<1, 3>>> vector_gradients;
      Tensor<1, 3>                           curl_vec_in_Hcurl;

      const FEValuesExtractors::Vector ve;

      const DoFHandler<3>  &dof_hand_Hcurl;
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

  // The following is the implementation of the solver's constructor.
  // After initializing the constant data members the constructor adjusts the
  // timer which is used to measure execution time. Next, it runs all the
  // functions that implement the functionality of the solver in the correct
  // order. The brackets around the functions such as `{ TMR("Save"); ...}` are
  // needed for a proper operation of the timer.
  Solver::Solver(unsigned int            p,
                 unsigned int            mapping_degree,
                 const Triangulation<3> &triangulation_Hcurl,
                 const DoFHandler<3>    &dof_handler_Hcurl,
                 const Vector<double>   &solution_Hcurl,
                 std::string             fname,
                 const Function<3>      *exact_solution)
    : fname(fname)
    , triangulation_Hcurl(triangulation_Hcurl)
    , dof_handler_Hcurl(dof_handler_Hcurl)
    , solution_Hcurl(solution_Hcurl)
    , fe_Hdiv(p)
    , exact_solution(exact_solution)
    , mapping_degree(mapping_degree)
  {
    Assert(exact_solution != nullptr,
           ExcMessage("The exact solution is missing."));

    TimerOutput::OutputFrequency tf =
      (print_time_tables) ? TimerOutput::summary : TimerOutput::never;

    TimerOutput timer(std::cout, tf, TimerOutput::cpu_and_wall_times_grouped);

    {
      TMR("Setup");
      setup();
    }
    {
      TMR("Assemble");
      assemble();
    }
    {
      TMR("Solve");
      solve();
    }

    if (exact_solution)
      {
        {
          TMR("Compute error norms");
          compute_error_norms();
        }

        if (project_exact_solution)
          {
            {
              TMR("Project exact solution");
              project_exact_solution_fcn();
            }
          }
      }
    {
      TMR("Save");
      save();
    }
  }

  // This function initializes the dofs, vectors and matrices. This time there
  // are no constraints as we do not apply Dirichlet boundary condition.
  void Solver::setup()
  {
    constraints.close();

    dof_handler_Hdiv.reinit(triangulation_Hcurl);
    dof_handler_Hdiv.distribute_dofs(fe_Hdiv);

    DynamicSparsityPattern dsp(dof_handler_Hdiv.n_dofs(),
                               dof_handler_Hdiv.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler_Hdiv, dsp, constraints, false);

    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
    solution_Hdiv.reinit(dof_handler_Hdiv.n_dofs());
    system_rhs.reinit(dof_handler_Hdiv.n_dofs());

    if ((project_exact_solution) && (exact_solution))
      projected_exact_solution.reinit(dof_handler_Hdiv.n_dofs());

    if (exact_solution)
      L2_per_cell.reinit(triangulation_Hcurl.n_active_cells());
  }



  // Formally, this function assembles the system of linear equations. In
  // reality, however, it just spells all the magic words to get the WorkStream
  // going. The interesting part, i.e., the actual assembling of the system
  // matrix and the right-hand side happens below in the
  // Solver::system_matrix_local function.
  void Solver::assemble()
  {
    WorkStream::run(IteratorPair(
                      IteratorTuple(dof_handler_Hdiv.begin_active(),
                                    dof_handler_Hcurl.begin_active())),
                    IteratorPair(IteratorTuple(dof_handler_Hdiv.end(),
                                               dof_handler_Hcurl.end())),
                    *this,
                    &Solver::system_matrix_local,
                    &Solver::copy_local_to_global,
                    AssemblyScratchData(fe_Hdiv,
                                        dof_handler_Hcurl,
                                        solution_Hcurl,
                                        mapping_degree),
                    AssemblyCopyData());
  }

  // The following two constructors initialize scratch data from the input
  // parameters and from another object of the same type, i.e., a copy
  // constructor.
  Solver::AssemblyScratchData::AssemblyScratchData(
    const FiniteElement<3> &fe,
    const DoFHandler<3>    &dof_hand_Hcurl,
    const Vector<double>   &dofs_Hcurl,
    unsigned int            mapping_degree)
    : mapping(mapping_degree)
    , fe_values_Hdiv(mapping,
                     fe,
                     QGauss<3>(fe.degree + 2),
                     update_values | update_JxW_values)
    , fe_values_Hcurl(mapping,
                      dof_hand_Hcurl.get_fe(),
                      QGauss<3>(fe.degree + 2),
                      update_gradients)
    , dofs_per_cell(fe_values_Hdiv.dofs_per_cell)
    , n_q_points(fe_values_Hdiv.get_quadrature().size())
    , vector_gradients(n_q_points, std::vector<Tensor<1, 3>>(3))
    , ve(0)
    , dof_hand_Hcurl(dof_hand_Hcurl)
    , dofs_Hcurl(dofs_Hcurl)
  {}

  Solver::AssemblyScratchData::AssemblyScratchData(
    const AssemblyScratchData &scratch_data)
    : mapping(scratch_data.mapping.get_degree())
    , fe_values_Hdiv(mapping,
                     scratch_data.fe_values_Hdiv.get_fe(),
                     scratch_data.fe_values_Hdiv.get_quadrature(),
                     update_values | update_JxW_values)
    , fe_values_Hcurl(mapping,
                      scratch_data.fe_values_Hcurl.get_fe(),
                      scratch_data.fe_values_Hcurl.get_quadrature(),
                      update_gradients)
    , dofs_per_cell(fe_values_Hdiv.dofs_per_cell)
    , n_q_points(fe_values_Hdiv.get_quadrature().size())
    , vector_gradients(n_q_points, std::vector<Tensor<1, 3>>(3))
    , ve(0)
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
    // cell, update the FE values, and compute the gradients of the input vector
    // field.
    copy_data.cell_matrix.reinit(scratch_data.dofs_per_cell,
                                 scratch_data.dofs_per_cell);

    copy_data.cell_rhs.reinit(scratch_data.dofs_per_cell);

    copy_data.local_dof_indices.resize(scratch_data.dofs_per_cell);

    scratch_data.fe_values_Hdiv.reinit(std::get<0>(*IP));
    scratch_data.fe_values_Hcurl.reinit(std::get<1>(*IP));

    scratch_data.fe_values_Hcurl.get_function_gradients(
      scratch_data.dofs_Hcurl, scratch_data.vector_gradients);

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
                  scratch_data.fe_values_Hdiv[VE].value(i, q_index) *
                  scratch_data.fe_values_Hdiv[VE].value(j, q_index) *
                  scratch_data.fe_values_Hdiv.JxW(q_index);
              }

            // The variable `curl_vec_in_Hcurl` denotes the curl of the input
            // vector field, $\vec{\nabla} \times \vec{T}$ or
            // $\vec{\nabla} \times \vec{A}$, depending on the context.
            scratch_data.curl_vec_in_Hcurl[0] =
              scratch_data.vector_gradients.at(q_index).at(2)[1] -
              scratch_data.vector_gradients.at(q_index).at(1)[2];

            scratch_data.curl_vec_in_Hcurl[1] =
              scratch_data.vector_gradients.at(q_index).at(0)[2] -
              scratch_data.vector_gradients.at(q_index).at(2)[0];

            scratch_data.curl_vec_in_Hcurl[2] =
              scratch_data.vector_gradients.at(q_index).at(1)[0] -
              scratch_data.vector_gradients.at(q_index).at(0)[1];

            copy_data.cell_rhs(i) += // Integral I_b
              scratch_data.curl_vec_in_Hcurl *
              scratch_data.fe_values_Hdiv[VE].value(i, q_index) *
              scratch_data.fe_values_Hdiv.JxW(q_index);
          }
      }

    // Finally, we query the dof indices on the current cell and store them
    // in the copy data structure, so we know to which locations of the system
    // matrix and right-hand side the components of the cell matrix and
    // cell right-hand side must be copied.
    std::get<0>(*IP)->get_dof_indices(copy_data.local_dof_indices);
  }

  // This function copies the components of a cell matrix and a cell right-hand
  // side into the system matrix, $A_{i,j}$, and the system right-hand side,
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
    Weight                     weight;
    const Function<3, double> *mask = &weight;

    VectorTools::integrate_difference(MappingQ<3>(mapping_degree),
                                      dof_handler_Hdiv,
                                      solution_Hdiv,
                                      *exact_solution,
                                      L2_per_cell,
                                      QGauss<3>(fe_Hdiv.degree + 4),
                                      VectorTools::L2_norm,
                                      mask);

    L2_norm = VectorTools::compute_global_error(triangulation_Hcurl,
                                                L2_per_cell,
                                                VectorTools::L2_norm);
  }

  void Solver::project_exact_solution_fcn()
  {
    AffineConstraints<double> constraints_empty;
    constraints_empty.close();

    VectorTools::project(MappingQ<3>(mapping_degree),
                         dof_handler_Hdiv,
                         constraints_empty,
                         QGauss<3>(fe_Hdiv.degree + 2),
                         *exact_solution,
                         projected_exact_solution);
  }

  // This function solves the system of linear equations. This time we are
  // dealing with a mass matrix. It has good spectral properties. Consequently,
  // we do not use a factor of 1000 as in preceding two solvers.
  // The stopping condition is
  // \f[
  // |\boldsymbol{b} - \boldsymbol{A}\boldsymbol{c}|
  // < 10^{-6} |\boldsymbol{b}|.
  // \f]
  void Solver::solve()
  {
    SolverControl control(system_rhs.size(),
                          1.0e-6 * system_rhs.l2_norm(),
                          false,
                          false);

    if (log_cg_convergence)
      control.enable_history_data();

    GrowingVectorMemory<Vector<double>> memory;
    SolverCG<Vector<double>>            cg(control, memory);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve(system_matrix, solution_Hdiv, system_rhs, preconditioner);

    if (log_cg_convergence)
      {
        const std::vector<double> history_data = control.get_history_data();

        std::ofstream ofs(fname + "_cg_convergence.csv");

        unsigned int i = 1;
        for (auto item : history_data)
          {
            ofs << i << ", " << item << "\n";
            i++;
          }
        ofs.close();
      }
  }

  // This function saves the computed fields into a vtu file. This time we also
  // save the projected exact solution and the $L^2$ error norm. The exact
  // solution is only saved if the `Settings::project_exact_solution = true`
  void Solver::save() const
  {
    std::vector<std::string> solution_names(3, "VectorField");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      interpretation(3,
                     DataComponentInterpretation::component_is_part_of_vector);

    DataOut<3> data_out;

    data_out.add_data_vector(dof_handler_Hdiv,
                             solution_Hdiv,
                             solution_names,
                             interpretation);

    if (project_exact_solution)
      {
        std::vector<std::string> solution_names_ex(3, "VectorFieldExact");

        data_out.add_data_vector(dof_handler_Hdiv,
                                 projected_exact_solution,
                                 solution_names_ex,
                                 interpretation);
      }

    if (exact_solution)
      {
        data_out.add_data_vector(L2_per_cell, "L2norm");
      }

    std::ofstream ofs;

    DataOutBase::VtkFlags flags;
    flags.write_higher_order_cells = true;
    data_out.set_flags(flags);

    const MappingQ<3> mapping(mapping_degree);

    data_out.build_patches(mapping,
                           fe_Hdiv.degree + 2,
                           DataOut<3>::CurvedCellRegion::curved_inner_cells);

    ofs.open(fname + ".vtu");
    data_out.write_vtu(ofs);

    ofs.close();
  }
} // namespace ProjectorHcurlToHdiv


// @sect3{The main loop}

// This class contains the main loop of the program.
class Batch : public Settings
{
public:
  void run()
  {
    if (Settings::nr_threads_max)
      MultithreadInfo::set_thread_limit(Settings::nr_threads_max);

    const unsigned int m = 2; // Mapping degree.
    const unsigned int p = Settings::fe_degree;

    MainOutputTable table_Jf(3);
    MainOutputTable table_B(3);

    table_Jf.clear();
    table_B.clear();

    std::cout << "Solving for (p = " << Settings::fe_degree
              << "): " << std::flush;

    for (unsigned int r = 6; r < 10; r++) // Mesh refinement parameter.
      {
        table_Jf.add_value("r", r);
        table_Jf.add_value("p", p);

        table_B.add_value("r", r);
        table_B.add_value("p", p);

        // Stage 1. Computing $\vec{T}$.

        std::cout << "T " << std::flush;

        SolverT::Solver T(p,
                          r,
                          m,
                          Settings::eta_squared_T,
                          "T_p" + std::to_string(p) + "_r" + std::to_string(r));

        T.clear();

        // Stage 2. Computing $\vec{J}_f$.

        std::cout << "Jf " << std::flush;

        ExactSolutions::FreeCurrentDensity Jf_exact;

        ProjectorHcurlToHdiv::Solver Jf(p,
                                        m,
                                        T.get_tria(),
                                        T.get_dof_handler(),
                                        T.get_solution(),
                                        "Jf_p" + std::to_string(p) + "_r" +
                                          std::to_string(r),
                                        &Jf_exact);

        Jf.clear();

        table_Jf.add_value("ndofs", Jf.get_n_dofs());
        table_Jf.add_value("ncells", Jf.get_n_cells());
        table_Jf.add_value("L2", Jf.get_L2_norm());

        // Stage 3. Computing $\vec{A}$.

        std::cout << "A " << std::flush;

        SolverA::Solver A(p,
                          m,
                          T.get_tria(),
                          T.get_dof_handler(),
                          T.get_solution(),
                          Settings::eta_squared_A,
                          "A_p" + std::to_string(p) + "_r" + std::to_string(r));

        A.clear();

        // Stage 4. Computing $\vec{B}$.

        std::cout << "B " << std::flush;

        ExactSolutions::MagneticField B_exact;

        ProjectorHcurlToHdiv::Solver B(p,
                                       m,
                                       T.get_tria(),
                                       A.get_dof_handler(),
                                       A.get_solution(),
                                       "B_p" + std::to_string(p) + "_r" +
                                         std::to_string(r),
                                       &B_exact);

        table_B.add_value("ndofs", B.get_n_dofs());
        table_B.add_value("ncells", B.get_n_cells());
        table_B.add_value("L2", B.get_L2_norm());
        // End stage 4.

        table_Jf.save("table_Jf_p" + std::to_string(p));
        table_B.save("table_B_p" + std::to_string(p));
      }
    std::cout << std::endl;
  }
};

int main()
{
  try
    {
      Batch b;
      b.run();
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

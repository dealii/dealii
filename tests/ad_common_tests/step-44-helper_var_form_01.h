// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Header file:
// This is a modified version of step-44, which tests the implementation of
// cell-level auto-differentiation (via a variational formulation).

#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_point_data.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/differentiation/ad.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/packaged_operation.h>
#include <deal.II/lac/precondition_selector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_selector.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <functional>
#include <iostream>

#include "../tests.h"
namespace Step44
{
  namespace AD = dealii::Differentiation::AD;
  namespace Parameters
  {
    struct FESystem
    {
      unsigned int poly_degree;
      unsigned int quad_order;
      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };
    void
    FESystem::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Finite element system");
      {
        prm.declare_entry("Polynomial degree",
                          "2",
                          Patterns::Integer(0),
                          "Displacement system polynomial order");
        prm.declare_entry("Quadrature order",
                          "3",
                          Patterns::Integer(0),
                          "Gauss quadrature order");
      }
      prm.leave_subsection();
    }
    void
    FESystem::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Finite element system");
      {
        poly_degree = prm.get_integer("Polynomial degree");
        quad_order  = prm.get_integer("Quadrature order");
      }
      prm.leave_subsection();
    }
    struct Geometry
    {
      unsigned int global_refinement;
      double       scale;
      double       p_p0;
      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };
    void
    Geometry::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry");
      {
        prm.declare_entry("Global refinement",
                          "2",
                          Patterns::Integer(0),
                          "Global refinement level");
        prm.declare_entry("Grid scale",
                          "1e-3",
                          Patterns::Double(0.0),
                          "Global grid scaling factor");
        prm.declare_entry("Pressure ratio p/p0",
                          "100",
                          Patterns::Selection("20|40|60|80|100"),
                          "Ratio of applied pressure to reference pressure");
      }
      prm.leave_subsection();
    }
    void
    Geometry::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry");
      {
        global_refinement = prm.get_integer("Global refinement");
        scale             = prm.get_double("Grid scale");
        p_p0              = prm.get_double("Pressure ratio p/p0");
      }
      prm.leave_subsection();
    }
    struct Materials
    {
      double nu;
      double mu;
      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };
    void
    Materials::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Material properties");
      {
        prm.declare_entry("Poisson's ratio",
                          "0.4999",
                          Patterns::Double(-1.0, 0.5),
                          "Poisson's ratio");
        prm.declare_entry("Shear modulus",
                          "80.194e6",
                          Patterns::Double(),
                          "Shear modulus");
      }
      prm.leave_subsection();
    }
    void
    Materials::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Material properties");
      {
        nu = prm.get_double("Poisson's ratio");
        mu = prm.get_double("Shear modulus");
      }
      prm.leave_subsection();
    }
    struct LinearSolver
    {
      std::string type_lin;
      double      tol_lin;
      bool        use_static_condensation;
      std::string preconditioner_type;
      double      preconditioner_relaxation;
      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };
    void
    LinearSolver::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Linear solver");
      {
        prm.declare_entry("Solver type",
                          "CG",
                          Patterns::Selection("CG|Direct"),
                          "Type of solver used to solve the linear system");
        prm.declare_entry("Residual",
                          "1e-6",
                          Patterns::Double(0.0),
                          "Linear solver residual (scaled by residual norm)");
        prm.declare_entry("Use static condensation",
                          "true",
                          Patterns::Bool(),
                          "Solve the full block system or a reduced problem");
        prm.declare_entry("Preconditioner type",
                          "ssor",
                          Patterns::Selection("jacobi|ssor"),
                          "Type of preconditioner");
        prm.declare_entry("Preconditioner relaxation",
                          "0.65",
                          Patterns::Double(0.0),
                          "Preconditioner relaxation value");
      }
      prm.leave_subsection();
    }
    void
    LinearSolver::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Linear solver");
      {
        type_lin                  = prm.get("Solver type");
        tol_lin                   = prm.get_double("Residual");
        use_static_condensation   = prm.get_bool("Use static condensation");
        preconditioner_type       = prm.get("Preconditioner type");
        preconditioner_relaxation = prm.get_double("Preconditioner relaxation");
      }
      prm.leave_subsection();
    }
    struct NonlinearSolver
    {
      unsigned int max_iterations_NR;
      double       tol_f;
      double       tol_u;
      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };
    void
    NonlinearSolver::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Nonlinear solver");
      {
        prm.declare_entry("Max iterations Newton-Raphson",
                          "10",
                          Patterns::Integer(0),
                          "Number of Newton-Raphson iterations allowed");
        prm.declare_entry("Tolerance force",
                          "1.0e-9",
                          Patterns::Double(0.0),
                          "Force residual tolerance");
        prm.declare_entry("Tolerance displacement",
                          "1.0e-6",
                          Patterns::Double(0.0),
                          "Displacement error tolerance");
      }
      prm.leave_subsection();
    }
    void
    NonlinearSolver::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Nonlinear solver");
      {
        max_iterations_NR = prm.get_integer("Max iterations Newton-Raphson");
        tol_f             = prm.get_double("Tolerance force");
        tol_u             = prm.get_double("Tolerance displacement");
      }
      prm.leave_subsection();
    }
    struct Time
    {
      double delta_t;
      double end_time;
      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };
    void
    Time::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Time");
      {
        prm.declare_entry("End time", "1", Patterns::Double(), "End time");
        prm.declare_entry("Time step size",
                          "0.1",
                          Patterns::Double(),
                          "Time step size");
      }
      prm.leave_subsection();
    }
    void
    Time::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Time");
      {
        end_time = prm.get_double("End time");
        delta_t  = prm.get_double("Time step size");
      }
      prm.leave_subsection();
    }
    struct AllParameters : public FESystem,
                           public Geometry,
                           public Materials,
                           public LinearSolver,
                           public NonlinearSolver,
                           public Time
    {
      AllParameters(const std::string &input_file);
      static void
      declare_parameters(ParameterHandler &prm);
      void
      parse_parameters(ParameterHandler &prm);
    };
    AllParameters::AllParameters(const std::string &input_file)
    {
      ParameterHandler prm;
      declare_parameters(prm);
      prm.parse_input(input_file);
      parse_parameters(prm);
    }
    void
    AllParameters::declare_parameters(ParameterHandler &prm)
    {
      FESystem::declare_parameters(prm);
      Geometry::declare_parameters(prm);
      Materials::declare_parameters(prm);
      LinearSolver::declare_parameters(prm);
      NonlinearSolver::declare_parameters(prm);
      Time::declare_parameters(prm);
    }
    void
    AllParameters::parse_parameters(ParameterHandler &prm)
    {
      FESystem::parse_parameters(prm);
      Geometry::parse_parameters(prm);
      Materials::parse_parameters(prm);
      LinearSolver::parse_parameters(prm);
      NonlinearSolver::parse_parameters(prm);
      Time::parse_parameters(prm);
    }
  } // namespace Parameters
  template <int dim>
  class StandardTensors
  {
  public:
    static const SymmetricTensor<2, dim> I;
    static const SymmetricTensor<4, dim> IxI;
    static const SymmetricTensor<4, dim> II;
    static const SymmetricTensor<4, dim> dev_P;
  };
  template <int dim>
  const SymmetricTensor<2, dim> StandardTensors<dim>::I =
    unit_symmetric_tensor<dim>();
  template <int dim>
  const SymmetricTensor<4, dim> StandardTensors<dim>::IxI = outer_product(I, I);
  template <int dim>
  const SymmetricTensor<4, dim> StandardTensors<dim>::II =
    identity_tensor<dim>();
  template <int dim>
  const SymmetricTensor<4, dim> StandardTensors<dim>::dev_P =
    deviator_tensor<dim>();
  class Time
  {
  public:
    Time(const double time_end, const double delta_t)
      : timestep(0)
      , time_current(0.0)
      , time_end(time_end)
      , delta_t(delta_t)
    {}
    virtual ~Time()
    {}
    double
    current() const
    {
      return time_current;
    }
    double
    end() const
    {
      return time_end;
    }
    double
    get_delta_t() const
    {
      return delta_t;
    }
    unsigned int
    get_timestep() const
    {
      return timestep;
    }
    void
    increment()
    {
      time_current += delta_t;
      ++timestep;
    }

  private:
    unsigned int timestep;
    double       time_current;
    const double time_end;
    const double delta_t;
  };
  template <int dim>
  class Material_Compressible_Neo_Hook_Three_Field
  {
  public:
    Material_Compressible_Neo_Hook_Three_Field(const double mu, const double nu)
      : kappa((2.0 * mu * (1.0 + nu)) / (3.0 * (1.0 - 2.0 * nu)))
      , c_1(mu / 2.0)
    {
      Assert(kappa > 0, ExcInternalError());
    }
    ~Material_Compressible_Neo_Hook_Three_Field()
    {}

    template <typename NumberType>
    NumberType
    get_Psi_iso(const SymmetricTensor<2, dim, NumberType> &C_bar)
    {
      return c_1 * (trace(C_bar) - dim);
    }
    template <typename NumberType>
    NumberType
    get_Psi_vol(const NumberType &J_tilde)
    {
      return (kappa / 4.0) * (J_tilde * J_tilde - 1.0 - 2.0 * log(J_tilde));
    }

  protected:
    const double kappa;
    const double c_1;
  };
  template <int dim>
  class PointHistory
  {
  public:
    PointHistory()
    {}
    virtual ~PointHistory()
    {}
    void
    setup_lqp(const Parameters::AllParameters &parameters)
    {
      material.reset(
        new Material_Compressible_Neo_Hook_Three_Field<dim>(parameters.mu,
                                                            parameters.nu));
    }

    template <typename NumberType>
    NumberType
    get_Psi(const Tensor<2, dim, NumberType> &F,
            const NumberType                 &p_tilde,
            const NumberType                 &J_tilde) const
    {
      const SymmetricTensor<2, dim, NumberType> C =
        symmetrize(transpose(F) * F);
      const NumberType                    det_F = determinant(F);
      SymmetricTensor<2, dim, NumberType> C_bar(C);
      C_bar *= pow(det_F, -2.0 / dim);

      NumberType psi_CpJ = material->get_Psi_iso(C_bar);
      psi_CpJ += material->get_Psi_vol(J_tilde);
      psi_CpJ += p_tilde * (det_F - J_tilde);

      return psi_CpJ;
    }

  private:
    std::shared_ptr<Material_Compressible_Neo_Hook_Three_Field<dim>> material;
  };
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  class Solid
  {
  public:
    Solid(const std::string &input_file);
    virtual ~Solid();
    void
    run();

  private:
    struct PerTaskData_ASM;
    struct ScratchData_ASM;
    struct PerTaskData_SC;
    struct ScratchData_SC;
    void
    make_grid();
    void
    system_setup();
    void
    determine_component_extractors();
    void
    assemble_system(const BlockVector<double> &solution_delta);
    void
    assemble_system_one_cell(
      const typename DoFHandler<dim>::active_cell_iterator &cell,
      ScratchData_ASM                                      &scratch,
      PerTaskData_ASM                                      &data) const;
    void
    copy_local_to_global_system(const PerTaskData_ASM &data);
    void
    assemble_sc();
    void
    assemble_sc_one_cell(
      const typename DoFHandler<dim>::active_cell_iterator &cell,
      ScratchData_SC                                       &scratch,
      PerTaskData_SC                                       &data);
    void
    copy_local_to_global_sc(const PerTaskData_SC &data);
    void
    make_constraints(const int &it_nr);
    void
    setup_qph();
    void
    solve_nonlinear_timestep(BlockVector<double> &solution_delta);
    std::pair<unsigned int, double>
    solve_linear_system(BlockVector<double> &newton_update);
    BlockVector<double>
    get_total_solution(const BlockVector<double> &solution_delta) const;
    void
                              output_results() const;
    Parameters::AllParameters parameters;
    double                    vol_reference;
    Triangulation<dim>        triangulation;
    Time                      time;
    mutable TimerOutput       timer;
    CellDataStorage<typename Triangulation<dim>::cell_iterator,
                    PointHistory<dim>>
                                     quadrature_point_history;
    const unsigned int               degree;
    const FESystem<dim>              fe;
    DoFHandler<dim>                  dof_handler_ref;
    const unsigned int               dofs_per_cell;
    const FEValuesExtractors::Vector u_fe;
    const FEValuesExtractors::Scalar p_fe;
    const FEValuesExtractors::Scalar J_fe;
    static const unsigned int        n_blocks          = 3;
    static const unsigned int        n_components      = dim + 2;
    static const unsigned int        first_u_component = 0;
    static const unsigned int        p_component       = dim;
    static const unsigned int        J_component       = dim + 1;
    enum
    {
      u_dof = 0,
      p_dof = 1,
      J_dof = 2
    };
    std::vector<types::global_dof_index> dofs_per_block;
    std::vector<types::global_dof_index> element_indices_u;
    std::vector<types::global_dof_index> element_indices_p;
    std::vector<types::global_dof_index> element_indices_J;
    const QGauss<dim>                    qf_cell;
    const QGauss<dim - 1>                qf_face;
    const unsigned int                   n_q_points;
    const unsigned int                   n_q_points_f;
    AffineConstraints<double>            constraints;
    BlockSparsityPattern                 sparsity_pattern;
    BlockSparseMatrix<double>            tangent_matrix;
    BlockVector<double>                  system_rhs;
    BlockVector<double>                  solution_n;
    struct Errors
    {
      Errors()
        : norm(1.0)
        , u(1.0)
        , p(1.0)
        , J(1.0)
      {}
      void
      reset()
      {
        norm = 1.0;
        u    = 1.0;
        p    = 1.0;
        J    = 1.0;
      }
      void
      normalise(const Errors &rhs)
      {
        if (rhs.norm != 0.0)
          norm /= rhs.norm;
        if (rhs.u != 0.0)
          u /= rhs.u;
        if (rhs.p != 0.0)
          p /= rhs.p;
        if (rhs.J != 0.0)
          J /= rhs.J;
      }
      double norm, u, p, J;
    };
    Errors error_residual, error_residual_0, error_residual_norm, error_update,
      error_update_0, error_update_norm;
    void
    get_error_residual(Errors &error_residual);
    void
    get_error_update(const BlockVector<double> &newton_update,
                     Errors                    &error_update);
    std::pair<double, double>
    get_error_dilation(const BlockVector<double> &solution_total) const;
    void
    print_conv_header();
    void
    print_conv_footer(const BlockVector<double> &solution_delta);
  };
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  Solid<dim, number_t, ad_type_code>::Solid(const std::string &input_file)
    : parameters(input_file)
    , triangulation(Triangulation<dim>::maximum_smoothing)
    , time(parameters.end_time, parameters.delta_t)
    , timer(std::cout, TimerOutput::never, TimerOutput::wall_times)
    , degree(parameters.poly_degree)
    , fe(FE_Q<dim>(parameters.poly_degree),
         dim, // displacement
         FE_DGP<dim>(parameters.poly_degree - 1),
         1, // pressure
         FE_DGP<dim>(parameters.poly_degree - 1),
         1)
    , // dilatation
    dof_handler_ref(triangulation)
    , dofs_per_cell(fe.dofs_per_cell)
    , u_fe(first_u_component)
    , p_fe(p_component)
    , J_fe(J_component)
    , dofs_per_block(n_blocks)
    , qf_cell(parameters.quad_order)
    , qf_face(parameters.quad_order)
    , n_q_points(qf_cell.size())
    , n_q_points_f(qf_face.size())
  {
    Assert(dim == 2 || dim == 3,
           ExcMessage("This problem only works in 2 or 3 space dimensions."));
    determine_component_extractors();
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  Solid<dim, number_t, ad_type_code>::~Solid()
  {
    dof_handler_ref.clear();
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  void
  Solid<dim, number_t, ad_type_code>::run()
  {
    make_grid();
    system_setup();
    {
      AffineConstraints<double> constraints;
      constraints.close();
      const ComponentSelectFunction<dim> J_mask(J_component, n_components);
      VectorTools::project(dof_handler_ref,
                           constraints,
                           QGauss<dim>(degree + 2),
                           J_mask,
                           solution_n);
    }
    output_results();
    time.increment();
    BlockVector<double> solution_delta(dofs_per_block);
    while (time.current() < time.end())
      {
        solution_delta = 0.0;
        solve_nonlinear_timestep(solution_delta);
        solution_n += solution_delta;
        output_results();
        // Output displacement at centre of traction surface
        {
          const Point<dim> soln_pt(
            dim == 3 ? Point<dim>(0.0, 1.0, 0.0) * parameters.scale :
                       Point<dim>(0.0, 1.0) * parameters.scale);
          typename DoFHandler<dim>::active_cell_iterator
            cell = dof_handler_ref.begin_active(),
            endc = dof_handler_ref.end();
          for (; cell != endc; ++cell)
            for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
              if (cell->vertex(v).distance(soln_pt) < 1e-6 * parameters.scale)
                {
                  Tensor<1, dim> soln;
                  for (unsigned int d = 0; d < dim; ++d)
                    soln[d] = solution_n(cell->vertex_dof_index(v, u_dof + d));
                  deallog << "Timestep " << time.get_timestep() << ": " << soln
                          << std::endl;
                }
        }
        time.increment();
      }
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  struct Solid<dim, number_t, ad_type_code>::PerTaskData_ASM
  {
    FullMatrix<double>                   cell_matrix;
    Vector<double>                       cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;
    PerTaskData_ASM(const unsigned int dofs_per_cell)
      : cell_matrix(dofs_per_cell, dofs_per_cell)
      , cell_rhs(dofs_per_cell)
      , local_dof_indices(dofs_per_cell)
    {}
    void
    reset()
    {
      cell_matrix = 0.0;
      cell_rhs    = 0.0;
    }
  };
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  struct Solid<dim, number_t, ad_type_code>::ScratchData_ASM
  {
    const BlockVector<double> &solution_total;
    FEValues<dim>              fe_values_ref;
    FEFaceValues<dim>          fe_face_values_ref;
    ScratchData_ASM(const FiniteElement<dim>  &fe_cell,
                    const QGauss<dim>         &qf_cell,
                    const UpdateFlags          uf_cell,
                    const QGauss<dim - 1>     &qf_face,
                    const UpdateFlags          uf_face,
                    const BlockVector<double> &solution_total)
      : solution_total(solution_total)
      , fe_values_ref(fe_cell, qf_cell, uf_cell)
      , fe_face_values_ref(fe_cell, qf_face, uf_face)
    {}
    ScratchData_ASM(const ScratchData_ASM &rhs)
      : solution_total(rhs.solution_total)
      , fe_values_ref(rhs.fe_values_ref.get_fe(),
                      rhs.fe_values_ref.get_quadrature(),
                      rhs.fe_values_ref.get_update_flags())
      , fe_face_values_ref(rhs.fe_face_values_ref.get_fe(),
                           rhs.fe_face_values_ref.get_quadrature(),
                           rhs.fe_face_values_ref.get_update_flags())
    {}
    void
    reset()
    {}
  };
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  struct Solid<dim, number_t, ad_type_code>::PerTaskData_SC
  {
    FullMatrix<double>                   cell_matrix;
    std::vector<types::global_dof_index> local_dof_indices;
    FullMatrix<double>                   k_orig;
    FullMatrix<double>                   k_pu;
    FullMatrix<double>                   k_pJ;
    FullMatrix<double>                   k_JJ;
    FullMatrix<double>                   k_pJ_inv;
    FullMatrix<double>                   k_bbar;
    FullMatrix<double>                   A;
    FullMatrix<double>                   B;
    FullMatrix<double>                   C;
    PerTaskData_SC(const unsigned int dofs_per_cell,
                   const unsigned int n_u,
                   const unsigned int n_p,
                   const unsigned int n_J)
      : cell_matrix(dofs_per_cell, dofs_per_cell)
      , local_dof_indices(dofs_per_cell)
      , k_orig(dofs_per_cell, dofs_per_cell)
      , k_pu(n_p, n_u)
      , k_pJ(n_p, n_J)
      , k_JJ(n_J, n_J)
      , k_pJ_inv(n_p, n_J)
      , k_bbar(n_u, n_u)
      , A(n_J, n_u)
      , B(n_J, n_u)
      , C(n_p, n_u)
    {}
    void
    reset()
    {}
  };
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  struct Solid<dim, number_t, ad_type_code>::ScratchData_SC
  {
    void
    reset()
    {}
  };
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  void
  Solid<dim, number_t, ad_type_code>::make_grid()
  {
    GridGenerator::hyper_rectangle(
      triangulation,
      (dim == 3 ? Point<dim>(0.0, 0.0, 0.0) : Point<dim>(0.0, 0.0)),
      (dim == 3 ? Point<dim>(1.0, 1.0, 1.0) : Point<dim>(1.0, 1.0)),
      true);
    GridTools::scale(parameters.scale, triangulation);
    triangulation.refine_global(std::max(1U, parameters.global_refinement));
    vol_reference = GridTools::volume(triangulation);
    std::cout << "Grid:\n\t Reference volume: " << vol_reference << std::endl;
    typename Triangulation<dim>::active_cell_iterator cell = triangulation
                                                               .begin_active(),
                                                      endc =
                                                        triangulation.end();
    for (; cell != endc; ++cell)
      for (const unsigned int face : GeometryInfo<dim>::face_indices())
        {
          if (cell->face(face)->at_boundary() == true &&
              cell->face(face)->center()[1] == 1.0 * parameters.scale)
            {
              if (dim == 3)
                {
                  if (cell->face(face)->center()[0] < 0.5 * parameters.scale &&
                      cell->face(face)->center()[2] < 0.5 * parameters.scale)
                    cell->face(face)->set_boundary_id(6);
                }
              else
                {
                  if (cell->face(face)->center()[0] < 0.5 * parameters.scale)
                    cell->face(face)->set_boundary_id(6);
                }
            }
        }
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  void
  Solid<dim, number_t, ad_type_code>::system_setup()
  {
    timer.enter_subsection("Setup system");
    std::vector<unsigned int> block_component(n_components,
                                              u_dof); // Displacement
    block_component[p_component] = p_dof;             // Pressure
    block_component[J_component] = J_dof;             // Dilatation
    dof_handler_ref.distribute_dofs(fe);
    DoFRenumbering::Cuthill_McKee(dof_handler_ref);
    DoFRenumbering::component_wise(dof_handler_ref, block_component);
    dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof_handler_ref, block_component);
    std::cout << "Triangulation:"
              << "\n\t Number of active cells: "
              << triangulation.n_active_cells()
              << "\n\t Number of degrees of freedom: "
              << dof_handler_ref.n_dofs() << std::endl;
    tangent_matrix.clear();
    {
      const types::global_dof_index n_dofs_u = dofs_per_block[u_dof];
      const types::global_dof_index n_dofs_p = dofs_per_block[p_dof];
      const types::global_dof_index n_dofs_J = dofs_per_block[J_dof];
      BlockDynamicSparsityPattern   dsp(n_blocks, n_blocks);
      dsp.block(u_dof, u_dof).reinit(n_dofs_u, n_dofs_u);
      dsp.block(u_dof, p_dof).reinit(n_dofs_u, n_dofs_p);
      dsp.block(u_dof, J_dof).reinit(n_dofs_u, n_dofs_J);
      dsp.block(p_dof, u_dof).reinit(n_dofs_p, n_dofs_u);
      dsp.block(p_dof, p_dof).reinit(n_dofs_p, n_dofs_p);
      dsp.block(p_dof, J_dof).reinit(n_dofs_p, n_dofs_J);
      dsp.block(J_dof, u_dof).reinit(n_dofs_J, n_dofs_u);
      dsp.block(J_dof, p_dof).reinit(n_dofs_J, n_dofs_p);
      dsp.block(J_dof, J_dof).reinit(n_dofs_J, n_dofs_J);
      dsp.collect_sizes();
      Table<2, DoFTools::Coupling> coupling(n_components, n_components);
      for (unsigned int ii = 0; ii < n_components; ++ii)
        for (unsigned int jj = 0; jj < n_components; ++jj)
          if (((ii < p_component) && (jj == J_component)) ||
              ((ii == J_component) && (jj < p_component)) ||
              ((ii == p_component) && (jj == p_component)))
            coupling[ii][jj] = DoFTools::none;
          else
            coupling[ii][jj] = DoFTools::always;
      DoFTools::make_sparsity_pattern(
        dof_handler_ref, coupling, dsp, constraints, false);
      sparsity_pattern.copy_from(dsp);
    }
    tangent_matrix.reinit(sparsity_pattern);
    system_rhs.reinit(dofs_per_block);
    system_rhs.collect_sizes();
    solution_n.reinit(dofs_per_block);
    solution_n.collect_sizes();
    setup_qph();
    timer.leave_subsection();
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  void
  Solid<dim, number_t, ad_type_code>::determine_component_extractors()
  {
    element_indices_u.clear();
    element_indices_p.clear();
    element_indices_J.clear();
    for (unsigned int k = 0; k < fe.dofs_per_cell; ++k)
      {
        const unsigned int k_group = fe.system_to_base_index(k).first.first;
        if (k_group == u_dof)
          element_indices_u.push_back(k);
        else if (k_group == p_dof)
          element_indices_p.push_back(k);
        else if (k_group == J_dof)
          element_indices_J.push_back(k);
        else
          {
            Assert(k_group <= J_dof, ExcInternalError());
          }
      }
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  void
  Solid<dim, number_t, ad_type_code>::setup_qph()
  {
    std::cout << "    Setting up quadrature point data..." << std::endl;
    quadrature_point_history.initialize(triangulation.begin_active(),
                                        triangulation.end(),
                                        n_q_points);
    for (typename Triangulation<dim>::active_cell_iterator cell =
           triangulation.begin_active();
         cell != triangulation.end();
         ++cell)
      {
        const std::vector<std::shared_ptr<PointHistory<dim>>> lqph =
          quadrature_point_history.get_data(cell);
        Assert(lqph.size() == n_q_points, ExcInternalError());
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          lqph[q_point]->setup_lqp(parameters);
      }
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  void
  Solid<dim, number_t, ad_type_code>::solve_nonlinear_timestep(
    BlockVector<double> &solution_delta)
  {
    std::cout << std::endl
              << "Timestep " << time.get_timestep() << " @ " << time.current()
              << "s" << std::endl;
    BlockVector<double> newton_update(dofs_per_block);
    error_residual.reset();
    error_residual_0.reset();
    error_residual_norm.reset();
    error_update.reset();
    error_update_0.reset();
    error_update_norm.reset();
    print_conv_header();
    unsigned int newton_iteration = 0;
    for (; newton_iteration < parameters.max_iterations_NR; ++newton_iteration)
      {
        std::cout << " " << std::setw(2) << newton_iteration << " "
                  << std::flush;
        make_constraints(newton_iteration);
        assemble_system(solution_delta);
        get_error_residual(error_residual);
        if (newton_iteration == 0)
          error_residual_0 = error_residual;
        error_residual_norm = error_residual;
        error_residual_norm.normalise(error_residual_0);
        if (newton_iteration > 0 && error_update_norm.u <= parameters.tol_u &&
            error_residual_norm.u <= parameters.tol_f)
          {
            std::cout << " CONVERGED! " << std::endl;
            print_conv_footer(solution_delta);
            break;
          }
        const std::pair<unsigned int, double> lin_solver_output =
          solve_linear_system(newton_update);
        get_error_update(newton_update, error_update);
        if (newton_iteration == 0)
          error_update_0 = error_update;
        error_update_norm = error_update;
        error_update_norm.normalise(error_update_0);
        solution_delta += newton_update;
        std::cout << " | " << std::fixed << std::setprecision(3) << std::setw(7)
                  << std::scientific << lin_solver_output.first << "  "
                  << lin_solver_output.second << "  "
                  << error_residual_norm.norm << "  " << error_residual_norm.u
                  << "  " << error_residual_norm.p << "  "
                  << error_residual_norm.J << "  " << error_update_norm.norm
                  << "  " << error_update_norm.u << "  " << error_update_norm.p
                  << "  " << error_update_norm.J << "  " << std::endl;
      }
    AssertThrow(newton_iteration <= parameters.max_iterations_NR,
                ExcMessage("No convergence in nonlinear solver!"));
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  void
  Solid<dim, number_t, ad_type_code>::print_conv_header()
  {
    static const unsigned int l_width = 144;
    for (unsigned int i = 0; i < l_width; ++i)
      std::cout << "_";
    std::cout << std::endl;
    std::cout << "           SOLVER STEP             "
              << " |  LIN_IT   LIN_RES    RES_NORM    "
              << " RES_U     RES_P      RES_J     NU_NORM     "
              << " NU_U       NU_P       NU_J " << std::endl;
    for (unsigned int i = 0; i < l_width; ++i)
      std::cout << "_";
    std::cout << std::endl;
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  void
  Solid<dim, number_t, ad_type_code>::print_conv_footer(
    const BlockVector<double> &solution_delta)
  {
    static const unsigned int l_width = 144;
    for (unsigned int i = 0; i < l_width; ++i)
      std::cout << "_";
    std::cout << std::endl;
    const std::pair<double, double> error_dil =
      get_error_dilation(get_total_solution(solution_delta));
    std::cout << "Relative errors:" << std::endl
              << "Displacement:\t" << error_update.u / error_update_0.u
              << std::endl
              << "Force: \t\t" << error_residual.u / error_residual_0.u
              << std::endl
              << "Dilatation:\t" << error_dil.first << std::endl
              << "v / V_0:\t" << error_dil.second * vol_reference << " / "
              << vol_reference << " = " << error_dil.second << std::endl;
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  std::pair<double, double>
  Solid<dim, number_t, ad_type_code>::get_error_dilation(
    const BlockVector<double> &solution_total) const
  {
    double        vol_current  = 0.0;
    double        dil_L2_error = 0.0;
    FEValues<dim> fe_values_ref(
      fe, qf_cell, update_values | update_gradients | update_JxW_values);
    std::vector<Tensor<2, dim>> solution_grads_u_total(qf_cell.size());
    std::vector<double>         solution_values_J_total(qf_cell.size());
    for (typename DoFHandler<dim>::active_cell_iterator cell =
           dof_handler_ref.begin_active();
         cell != dof_handler_ref.end();
         ++cell)
      {
        fe_values_ref.reinit(cell);
        fe_values_ref[u_fe].get_function_gradients(solution_total,
                                                   solution_grads_u_total);
        fe_values_ref[J_fe].get_function_values(solution_total,
                                                solution_values_J_total);
        const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
          quadrature_point_history.get_data(cell);
        Assert(lqph.size() == n_q_points, ExcInternalError());
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {
            const double det_F_qp = determinant(
              StandardTensors<dim>::I + solution_grads_u_total[q_point]);
            const double J_tilde_qp = solution_values_J_total[q_point];
            const double the_error_qp_squared = pow((det_F_qp - J_tilde_qp), 2);
            const double JxW                  = fe_values_ref.JxW(q_point);
            dil_L2_error += the_error_qp_squared * JxW;
            vol_current += det_F_qp * JxW;
          }
      }
    Assert(vol_current > 0.0, ExcInternalError());

    return std::make_pair(std::sqrt(dil_L2_error), vol_current / vol_reference);
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  void
  Solid<dim, number_t, ad_type_code>::get_error_residual(Errors &error_residual)
  {
    BlockVector<double> error_res(dofs_per_block);
    for (unsigned int i = 0; i < dof_handler_ref.n_dofs(); ++i)
      if (!constraints.is_constrained(i))
        error_res(i) = system_rhs(i);
    error_residual.norm = error_res.l2_norm();
    error_residual.u    = error_res.block(u_dof).l2_norm();
    error_residual.p    = error_res.block(p_dof).l2_norm();
    error_residual.J    = error_res.block(J_dof).l2_norm();
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  void
  Solid<dim, number_t, ad_type_code>::get_error_update(
    const BlockVector<double> &newton_update,
    Errors                    &error_update)
  {
    BlockVector<double> error_ud(dofs_per_block);
    for (unsigned int i = 0; i < dof_handler_ref.n_dofs(); ++i)
      if (!constraints.is_constrained(i))
        error_ud(i) = newton_update(i);
    error_update.norm = error_ud.l2_norm();
    error_update.u    = error_ud.block(u_dof).l2_norm();
    error_update.p    = error_ud.block(p_dof).l2_norm();
    error_update.J    = error_ud.block(J_dof).l2_norm();
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  BlockVector<double>
  Solid<dim, number_t, ad_type_code>::get_total_solution(
    const BlockVector<double> &solution_delta) const
  {
    BlockVector<double> solution_total(solution_n);
    solution_total += solution_delta;
    return solution_total;
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  void
  Solid<dim, number_t, ad_type_code>::assemble_system(
    const BlockVector<double> &solution_delta)
  {
    timer.enter_subsection("Assemble system");
    std::cout << " ASM_SYS " << std::flush;
    tangent_matrix = 0.0;
    system_rhs     = 0.0;
    const BlockVector<double> solution_total(
      get_total_solution(solution_delta));
    const UpdateFlags uf_cell(update_values | update_gradients |
                              update_JxW_values);
    const UpdateFlags uf_face(update_values | update_normal_vectors |
                              update_JxW_values);
    PerTaskData_ASM   per_task_data(dofs_per_cell);
    ScratchData_ASM   scratch_data(
      fe, qf_cell, uf_cell, qf_face, uf_face, solution_total);
    // ADOL-C is incompatible with TBB
    // WorkStream::run(dof_handler_ref.begin_active(),
    //                 dof_handler_ref.end(),
    //                 std::bind(&Solid<dim,number_t,ad_type_code>::assemble_system_one_cell,
    //                                 this,
    //                                 std::placeholders::_1,
    //                                 std::placeholders::_2,
    //                                 std::placeholders::_3),
    //                 std::bind(&Solid<dim,number_t,ad_type_code>::copy_local_to_global_system,
    //                                 this,
    //                                 std::placeholders::_1),
    //                 scratch_data,
    //                 per_task_data);
    for (typename DoFHandler<dim>::active_cell_iterator cell =
           dof_handler_ref.begin_active();
         cell != dof_handler_ref.end();
         ++cell)
      {
        assemble_system_one_cell(cell, scratch_data, per_task_data);
        copy_local_to_global_system(per_task_data);
      }
    timer.leave_subsection();
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  void
  Solid<dim, number_t, ad_type_code>::copy_local_to_global_system(
    const PerTaskData_ASM &data)
  {
    if (data.cell_matrix.frobenius_norm() > 1e-12)
      constraints.distribute_local_to_global(data.cell_matrix,
                                             data.cell_rhs,
                                             data.local_dof_indices,
                                             tangent_matrix,
                                             system_rhs);
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  void
  Solid<dim, number_t, ad_type_code>::assemble_system_one_cell(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    ScratchData_ASM                                      &scratch,
    PerTaskData_ASM                                      &data) const
  {
    data.reset();
    scratch.reset();
    scratch.fe_values_ref.reinit(cell);
    cell->get_dof_indices(data.local_dof_indices);

    const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
      quadrature_point_history.get_data(cell);
    Assert(lqph.size() == n_q_points, ExcInternalError());

    const unsigned int n_independent_variables = data.local_dof_indices.size();

    using ADHelper     = AD::EnergyFunctional<ad_type_code, number_t>;
    using ADNumberType = typename ADHelper::ad_type;
    ADHelper ad_helper(n_independent_variables);
    ad_helper.set_tape_buffer_sizes(); // Increase the buffer size from the
                                       // default values

    const int  tape_no = 1;
    const bool is_recording =
      ad_helper.start_recording_operations(tape_no, //  material_id
                                           true,    // overwrite_tape
                                           true);   // keep

    if (is_recording == true)
      {
        // Set the values for all DoFs
        ad_helper.register_dof_values(scratch.solution_total,
                                      data.local_dof_indices);
        const std::vector<ADNumberType> dof_values_ad =
          ad_helper.get_sensitive_dof_values();

        // Compute all values, gradients etc. based on sensitive AD DoF values
        std::vector<Tensor<2, dim, ADNumberType>> Grad_u(
          n_q_points, Tensor<2, dim, ADNumberType>());
        std::vector<ADNumberType> p_tilde(n_q_points, ADNumberType(0.0));
        std::vector<ADNumberType> J_tilde(n_q_points, ADNumberType(0.0));
        scratch.fe_values_ref[u_fe]
          .get_function_gradients_from_local_dof_values(dof_values_ad, Grad_u);
        scratch.fe_values_ref[p_fe].get_function_values_from_local_dof_values(
          dof_values_ad, p_tilde);
        scratch.fe_values_ref[J_fe].get_function_values_from_local_dof_values(
          dof_values_ad, J_tilde);

        // Compute the total energy = (internal - external) energies
        ADNumberType energy_ad = ADNumberType(0.0);
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {
            const Tensor<2, dim, ADNumberType> F =
              unit_symmetric_tensor<dim>() + Grad_u[q_point];
            const Tensor<2, dim, ADNumberType> F_inv = invert(F);
            const ADNumberType                 det_F = determinant(F);
            Assert(numbers::value_is_greater_than(det_F, 0.0),
                   ExcMessage("Negative jacobian detected!"));

            const double JxW = scratch.fe_values_ref.JxW(q_point);

            energy_ad +=
              lqph[q_point]->get_Psi(F, p_tilde[q_point], J_tilde[q_point]) *
              JxW;
          }
        for (const unsigned int face : GeometryInfo<dim>::face_indices())
          if (cell->face(face)->at_boundary() == true &&
              cell->face(face)->boundary_id() == 6)
            {
              scratch.fe_face_values_ref.reinit(cell, face);

              std::vector<Tensor<1, dim, ADNumberType>> u_face(
                n_q_points_f, Tensor<1, dim, ADNumberType>());
              scratch.fe_face_values_ref[u_fe]
                .get_function_values_from_local_dof_values(dof_values_ad,
                                                           u_face);

              for (unsigned int f_q_point = 0; f_q_point < n_q_points_f;
                   ++f_q_point)
                {
                  const Tensor<1, dim> &N =
                    scratch.fe_face_values_ref.normal_vector(f_q_point);
                  static const double p0 =
                    -4.0 / (parameters.scale * parameters.scale);
                  const double time_ramp = (time.current() / time.end());
                  const double pressure  = p0 * parameters.p_p0 * time_ramp;
                  const Tensor<1, dim> traction = pressure * N;

                  const double JxW = scratch.fe_face_values_ref.JxW(f_q_point);

                  energy_ad -= (u_face[f_q_point] * traction) * JxW;
                }
            }

        ad_helper.register_energy_functional(energy_ad);
        ad_helper.stop_recording_operations(false /*write_tapes_to_file*/);
      }
    else
      {
        Assert(is_recording == true, ExcInternalError());
      }

    // Unnecessary when keep == true
    // ad_helper.activate_recorded_tape(tape_no);
    // ad_helper.set_dof_values(scratch.solution_total,
    //                          data.local_dof_indices);

    // Compute the residual values and their jacobian for the new evaluation
    // point
    ad_helper.compute_residual(data.cell_rhs);
    data.cell_rhs *= -1.0; // RHS = - residual
    ad_helper.compute_linearization(data.cell_matrix);
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  void
  Solid<dim, number_t, ad_type_code>::make_constraints(const int &it_nr)
  {
    std::cout << " CST " << std::flush;
    if (it_nr > 1)
      return;
    constraints.clear();
    const bool                       apply_dirichlet_bc = (it_nr == 0);
    const FEValuesExtractors::Scalar x_displacement(0);
    const FEValuesExtractors::Scalar y_displacement(1);
    {
      const int boundary_id = 0;
      if (apply_dirichlet_bc == true)
        VectorTools::interpolate_boundary_values(
          dof_handler_ref,
          boundary_id,
          Functions::ZeroFunction<dim>(n_components),
          constraints,
          fe.component_mask(x_displacement));
      else
        VectorTools::interpolate_boundary_values(
          dof_handler_ref,
          boundary_id,
          Functions::ZeroFunction<dim>(n_components),
          constraints,
          fe.component_mask(x_displacement));
    }
    {
      const int boundary_id = 2;
      if (apply_dirichlet_bc == true)
        VectorTools::interpolate_boundary_values(
          dof_handler_ref,
          boundary_id,
          Functions::ZeroFunction<dim>(n_components),
          constraints,
          fe.component_mask(y_displacement));
      else
        VectorTools::interpolate_boundary_values(
          dof_handler_ref,
          boundary_id,
          Functions::ZeroFunction<dim>(n_components),
          constraints,
          fe.component_mask(y_displacement));
    }
    if (dim == 3)
      {
        const FEValuesExtractors::Scalar z_displacement(2);
        {
          const int boundary_id = 3;
          if (apply_dirichlet_bc == true)
            VectorTools::interpolate_boundary_values(
              dof_handler_ref,
              boundary_id,
              Functions::ZeroFunction<dim>(n_components),
              constraints,
              (fe.component_mask(x_displacement) |
               fe.component_mask(z_displacement)));
          else
            VectorTools::interpolate_boundary_values(
              dof_handler_ref,
              boundary_id,
              Functions::ZeroFunction<dim>(n_components),
              constraints,
              (fe.component_mask(x_displacement) |
               fe.component_mask(z_displacement)));
        }
        {
          const int boundary_id = 4;
          if (apply_dirichlet_bc == true)
            VectorTools::interpolate_boundary_values(
              dof_handler_ref,
              boundary_id,
              Functions::ZeroFunction<dim>(n_components),
              constraints,
              fe.component_mask(z_displacement));
          else
            VectorTools::interpolate_boundary_values(
              dof_handler_ref,
              boundary_id,
              Functions::ZeroFunction<dim>(n_components),
              constraints,
              fe.component_mask(z_displacement));
        }
        {
          const int boundary_id = 6;
          if (apply_dirichlet_bc == true)
            VectorTools::interpolate_boundary_values(
              dof_handler_ref,
              boundary_id,
              Functions::ZeroFunction<dim>(n_components),
              constraints,
              (fe.component_mask(x_displacement) |
               fe.component_mask(z_displacement)));
          else
            VectorTools::interpolate_boundary_values(
              dof_handler_ref,
              boundary_id,
              Functions::ZeroFunction<dim>(n_components),
              constraints,
              (fe.component_mask(x_displacement) |
               fe.component_mask(z_displacement)));
        }
      }
    else
      {
        {
          const int boundary_id = 3;
          if (apply_dirichlet_bc == true)
            VectorTools::interpolate_boundary_values(
              dof_handler_ref,
              boundary_id,
              Functions::ZeroFunction<dim>(n_components),
              constraints,
              (fe.component_mask(x_displacement)));
          else
            VectorTools::interpolate_boundary_values(
              dof_handler_ref,
              boundary_id,
              Functions::ZeroFunction<dim>(n_components),
              constraints,
              (fe.component_mask(x_displacement)));
        }
        {
          const int boundary_id = 6;
          if (apply_dirichlet_bc == true)
            VectorTools::interpolate_boundary_values(
              dof_handler_ref,
              boundary_id,
              Functions::ZeroFunction<dim>(n_components),
              constraints,
              (fe.component_mask(x_displacement)));
          else
            VectorTools::interpolate_boundary_values(
              dof_handler_ref,
              boundary_id,
              Functions::ZeroFunction<dim>(n_components),
              constraints,
              (fe.component_mask(x_displacement)));
        }
      }
    constraints.close();
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  void
  Solid<dim, number_t, ad_type_code>::assemble_sc()
  {
    timer.enter_subsection("Perform static condensation");
    std::cout << " ASM_SC " << std::flush;
    PerTaskData_SC per_task_data(dofs_per_cell,
                                 element_indices_u.size(),
                                 element_indices_p.size(),
                                 element_indices_J.size());
    ScratchData_SC scratch_data;
    WorkStream::run(dof_handler_ref.begin_active(),
                    dof_handler_ref.end(),
                    *this,
                    &Solid::assemble_sc_one_cell,
                    &Solid::copy_local_to_global_sc,
                    scratch_data,
                    per_task_data);
    timer.leave_subsection();
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  void
  Solid<dim, number_t, ad_type_code>::copy_local_to_global_sc(
    const PerTaskData_SC &data)
  {
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
        tangent_matrix.add(data.local_dof_indices[i],
                           data.local_dof_indices[j],
                           data.cell_matrix(i, j));
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  void
  Solid<dim, number_t, ad_type_code>::assemble_sc_one_cell(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    ScratchData_SC                                       &scratch,
    PerTaskData_SC                                       &data)
  {
    data.reset();
    scratch.reset();
    cell->get_dof_indices(data.local_dof_indices);
    data.k_orig.extract_submatrix_from(tangent_matrix,
                                       data.local_dof_indices,
                                       data.local_dof_indices);
    data.k_pu.extract_submatrix_from(data.k_orig,
                                     element_indices_p,
                                     element_indices_u);
    data.k_pJ.extract_submatrix_from(data.k_orig,
                                     element_indices_p,
                                     element_indices_J);
    data.k_JJ.extract_submatrix_from(data.k_orig,
                                     element_indices_J,
                                     element_indices_J);
    data.k_pJ_inv.invert(data.k_pJ);
    data.k_pJ_inv.mmult(data.A, data.k_pu);
    data.k_JJ.mmult(data.B, data.A);
    data.k_pJ_inv.Tmmult(data.C, data.B);
    data.k_pu.Tmmult(data.k_bbar, data.C);
    data.k_bbar.scatter_matrix_to(element_indices_u,
                                  element_indices_u,
                                  data.cell_matrix);
    data.k_pJ_inv.add(-1.0, data.k_pJ);
    data.k_pJ_inv.scatter_matrix_to(element_indices_p,
                                    element_indices_J,
                                    data.cell_matrix);
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  std::pair<unsigned int, double>
  Solid<dim, number_t, ad_type_code>::solve_linear_system(
    BlockVector<double> &newton_update)
  {
    unsigned int lin_it  = 0;
    double       lin_res = 0.0;
    if (parameters.use_static_condensation == true)
      {
        BlockVector<double> A(dofs_per_block);
        BlockVector<double> B(dofs_per_block);
        {
          assemble_sc();
          tangent_matrix.block(p_dof, J_dof)
            .vmult(A.block(J_dof), system_rhs.block(p_dof));
          tangent_matrix.block(J_dof, J_dof)
            .vmult(B.block(J_dof), A.block(J_dof));
          A.block(J_dof) = system_rhs.block(J_dof);
          A.block(J_dof) -= B.block(J_dof);
          tangent_matrix.block(p_dof, J_dof)
            .Tvmult(A.block(p_dof), A.block(J_dof));
          tangent_matrix.block(u_dof, p_dof)
            .vmult(A.block(u_dof), A.block(p_dof));
          system_rhs.block(u_dof) -= A.block(u_dof);
          timer.enter_subsection("Linear solver");
          std::cout << " SLV " << std::flush;
          if (parameters.type_lin == "CG")
            {
              const unsigned int solver_its =
                tangent_matrix.block(u_dof, u_dof).m();
              const double tol_sol =
                parameters.tol_lin * system_rhs.block(u_dof).l2_norm();
              SolverControl solver_control(solver_its, tol_sol, false, false);
              GrowingVectorMemory<Vector<double>> GVM;
              SolverCG<Vector<double>> solver_CG(solver_control, GVM);
              PreconditionSelector<SparseMatrix<double>, Vector<double>>
                preconditioner(parameters.preconditioner_type,
                               parameters.preconditioner_relaxation);
              preconditioner.use_matrix(tangent_matrix.block(u_dof, u_dof));
              solver_CG.solve(tangent_matrix.block(u_dof, u_dof),
                              newton_update.block(u_dof),
                              system_rhs.block(u_dof),
                              preconditioner);
              lin_it  = solver_control.last_step();
              lin_res = solver_control.last_value();
            }
          else if (parameters.type_lin == "Direct")
            {
              SparseDirectUMFPACK A_direct;
              A_direct.initialize(tangent_matrix.block(u_dof, u_dof));
              A_direct.vmult(newton_update.block(u_dof),
                             system_rhs.block(u_dof));
              lin_it  = 1;
              lin_res = 0.0;
            }
          else
            Assert(false, ExcMessage("Linear solver type not implemented"));
          timer.leave_subsection();
        }
        constraints.distribute(newton_update);
        timer.enter_subsection("Linear solver postprocessing");
        std::cout << " PP " << std::flush;
        {
          tangent_matrix.block(p_dof, u_dof)
            .vmult(A.block(p_dof), newton_update.block(u_dof));
          A.block(p_dof) *= -1.0;
          A.block(p_dof) += system_rhs.block(p_dof);
          tangent_matrix.block(p_dof, J_dof)
            .vmult(newton_update.block(J_dof), A.block(p_dof));
        }
        constraints.distribute(newton_update);
        {
          tangent_matrix.block(J_dof, J_dof)
            .vmult(A.block(J_dof), newton_update.block(J_dof));
          A.block(J_dof) *= -1.0;
          A.block(J_dof) += system_rhs.block(J_dof);
          tangent_matrix.block(p_dof, J_dof)
            .Tvmult(newton_update.block(p_dof), A.block(J_dof));
        }
        constraints.distribute(newton_update);
        timer.leave_subsection();
      }
    else
      {
        std::cout << " ------ " << std::flush;
        timer.enter_subsection("Linear solver");
        std::cout << " SLV " << std::flush;
        if (parameters.type_lin == "CG")
          {
            const Vector<double> &f_u = system_rhs.block(u_dof);
            const Vector<double> &f_p = system_rhs.block(p_dof);
            const Vector<double> &f_J = system_rhs.block(J_dof);
            Vector<double>       &d_u = newton_update.block(u_dof);
            Vector<double>       &d_p = newton_update.block(p_dof);
            Vector<double>       &d_J = newton_update.block(J_dof);
            const auto            K_uu =
              linear_operator(tangent_matrix.block(u_dof, u_dof));
            const auto K_up =
              linear_operator(tangent_matrix.block(u_dof, p_dof));
            const auto K_pu =
              linear_operator(tangent_matrix.block(p_dof, u_dof));
            const auto K_Jp =
              linear_operator(tangent_matrix.block(J_dof, p_dof));
            const auto K_JJ =
              linear_operator(tangent_matrix.block(J_dof, J_dof));
            PreconditionSelector<SparseMatrix<double>, Vector<double>>
              preconditioner_K_Jp_inv("jacobi");
            preconditioner_K_Jp_inv.use_matrix(
              tangent_matrix.block(J_dof, p_dof));
            ReductionControl solver_control_K_Jp_inv(
              tangent_matrix.block(J_dof, p_dof).m(),
              1.0e-30,
              parameters.tol_lin);
            SolverSelector<Vector<double>> solver_K_Jp_inv;
            solver_K_Jp_inv.select("cg");
            solver_K_Jp_inv.set_control(solver_control_K_Jp_inv);
            const auto K_Jp_inv =
              inverse_operator(K_Jp, solver_K_Jp_inv, preconditioner_K_Jp_inv);
            const auto K_pJ_inv     = transpose_operator(K_Jp_inv);
            const auto K_pp_bar     = K_Jp_inv * K_JJ * K_pJ_inv;
            const auto K_uu_bar_bar = K_up * K_pp_bar * K_pu;
            const auto K_uu_con     = K_uu + K_uu_bar_bar;
            PreconditionSelector<SparseMatrix<double>, Vector<double>>
              preconditioner_K_con_inv(parameters.preconditioner_type,
                                       parameters.preconditioner_relaxation);
            preconditioner_K_con_inv.use_matrix(
              tangent_matrix.block(u_dof, u_dof));
            ReductionControl solver_control_K_con_inv(
              tangent_matrix.block(u_dof, u_dof).m(),
              1.0e-30,
              parameters.tol_lin);
            SolverSelector<Vector<double>> solver_K_con_inv;
            solver_K_con_inv.select("cg");
            solver_K_con_inv.set_control(solver_control_K_con_inv);
            const auto K_uu_con_inv =
              inverse_operator(K_uu_con,
                               solver_K_con_inv,
                               preconditioner_K_con_inv);
            d_u =
              K_uu_con_inv * (f_u - K_up * (K_Jp_inv * f_J - K_pp_bar * f_p));
            timer.leave_subsection();
            timer.enter_subsection("Linear solver postprocessing");
            std::cout << " PP " << std::flush;
            d_J     = K_pJ_inv * (f_p - K_pu * d_u);
            d_p     = K_Jp_inv * (f_J - K_JJ * d_J);
            lin_it  = solver_control_K_con_inv.last_step();
            lin_res = solver_control_K_con_inv.last_value();
          }
        else if (parameters.type_lin == "Direct")
          {
            SparseDirectUMFPACK A_direct;
            A_direct.initialize(tangent_matrix);
            A_direct.vmult(newton_update, system_rhs);
            lin_it  = 1;
            lin_res = 0.0;
            std::cout << " -- " << std::flush;
          }
        else
          Assert(false, ExcMessage("Linear solver type not implemented"));
        timer.leave_subsection();
        constraints.distribute(newton_update);
      }
    return std::make_pair(lin_it, lin_res);
  }
  template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
  void
  Solid<dim, number_t, ad_type_code>::output_results() const
  {
    DataOut<dim> data_out;
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);
    std::vector<std::string> solution_name(dim, "displacement");
    solution_name.push_back("pressure");
    solution_name.push_back("dilatation");
    data_out.attach_dof_handler(dof_handler_ref);
    data_out.add_data_vector(solution_n,
                             solution_name,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation);
    Vector<double> soln(solution_n.size());
    for (unsigned int i = 0; i < soln.size(); ++i)
      soln(i) = solution_n(i);
    MappingQEulerian<dim> q_mapping(degree, dof_handler_ref, soln);
    data_out.build_patches(q_mapping, degree);
    std::ostringstream filename;
    filename << "solution-" << dim << "d-" << time.get_timestep() << ".vtk";
    std::ofstream output(filename.str());
    data_out.write_vtk(output);
  }
} // namespace Step44

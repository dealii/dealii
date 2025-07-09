/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 1999 - 2025 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

// @sect3{Include files}

// The DPG method requires a large breadth of elements type which are included
// below:

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_trace.h>

// The rest of the includes are some well-known files:

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_faces.h>
#include <deal.II/numerics/vector_tools.h>

#include <chrono>
#include <fstream>
#include <iostream>

namespace Step100
{
  using namespace dealii;

  // In the following function declaration, we create the analytical
  // solutions of the velocity field ($\mathbf{u}$) and the pressure field
  // ($p^*$) as well as the associated boundary values. However, in this
  // tutorial, we will avoid the use of deal.ii complex arithmetic capabilities
  // and only use the complex functions that are defined in the C++ standard
  // library. Consequently, in what follows, we will separate the real and
  // imaginary component of our spaces. Therefore, we will also define two
  // implementations of each function, one for the real component and one for
  // the imaginary one.

  // We start and create the analytical solution class for kinematic pressure
  // ($p^*$). The analytical solution depends on the wavenumber $k$ and the
  // angle
  // $\theta$ which are passed to the constructor.
  template <int dim>
  class AnalyticalSolution_p_real : public Function<dim>
  {
  public:
    AnalyticalSolution_p_real(const double wavenumber, const double theta)
      : Function<dim>()
      , wavenumber(wavenumber)
      , theta(theta)
    {}
    virtual double value(const Point<dim>  &p,
                         const unsigned int component) const override;

  private:
    const double wavenumber;
    const double theta;
  };

  template <int dim>
  double AnalyticalSolution_p_real<dim>::value(
    const Point<dim>                   &p,
    [[maybe_unused]] const unsigned int component) const
  {
    constexpr std::complex<double> imag(0., 1.);

    return std::exp(-imag * wavenumber *
                    (p[0] * std::cos(theta) + p[1] * std::sin(theta)))
      .real();
  }

  // The same goes for the imaginary part of the analytical solution.
  template <int dim>
  class AnalyticalSolution_p_imag : public Function<dim>
  {
  public:
    AnalyticalSolution_p_imag(const double wavenumber, const double theta)
      : Function<dim>()
      , wavenumber(wavenumber)
      , theta(theta)
    {}
    virtual double value(const Point<dim>  &p,
                         const unsigned int component) const override;

  private:
    const double wavenumber;
    const double theta;
  };

  template <int dim>
  double AnalyticalSolution_p_imag<dim>::value(
    const Point<dim>                   &p,
    [[maybe_unused]] const unsigned int component) const
  {
    constexpr std::complex<double> imag(0., 1.);

    return std::exp(-imag * wavenumber *
                    (p[0] * std::cos(theta) + p[1] * std::sin(theta)))
      .imag();
  }


  // Now we will do the same for the velocity field ($\mathbf{u}$) and create
  // two class that will return the real or the imaginary part of our
  // analytical solution. This class is similar to the previous ones but
  // since the velocity field is a vector, we will need the function to return
  // $dim$ component. For our problem of interest, $dim = 2$.
  template <int dim>
  class AnalyticalSolution_u_real : public Function<dim>
  {
  public:
    AnalyticalSolution_u_real(const double wavenumber, const double theta)
      : Function<dim>(2)
      , wavenumber(wavenumber)
      , theta(theta)
    {}
    virtual double value(const Point<dim>  &p,
                         const unsigned int component) const override;

  private:
    const double wavenumber;
    const double theta;
  };

  template <int dim>
  double
  AnalyticalSolution_u_real<dim>::value(const Point<dim>  &p,
                                        const unsigned int component) const
  {
    AssertIndexRange(component, 2);

    constexpr std::complex<double> imag(0., 1.);

    if (component == 0)
      return (std::cos(theta) *
              std::exp(-imag * wavenumber *
                       (p[0] * std::cos(theta) + p[1] * std::sin(theta))))
        .real();


    else
      return (std::sin(theta) *
              std::exp(-imag * wavenumber *
                       (p[0] * std::cos(theta) + p[1] * std::sin(theta))))
        .real();
  }

  template <int dim>
  class AnalyticalSolution_u_imag : public Function<dim>
  {
  public:
    AnalyticalSolution_u_imag(double wavenumber, double theta)
      : Function<dim>(2)
      , wavenumber(wavenumber)
      , theta(theta)
    {}
    virtual double value(const Point<dim>  &p,
                         const unsigned int component) const override;

  private:
    const double wavenumber;
    const double theta;
  };

  template <int dim>
  double
  AnalyticalSolution_u_imag<dim>::value(const Point<dim>  &p,
                                        const unsigned int component) const
  {
    AssertIndexRange(component, 2);

    constexpr std::complex<double> imag(0., 1.);

    if (component == 0)
      return (std::cos(theta) *
              std::exp(-imag * wavenumber *
                       (p[0] * std::cos(theta) + p[1] * std::sin(theta))))
        .imag();

    else
      return (std::sin(theta) *
              std::exp(-imag * wavenumber *
                       (p[0] * std::cos(theta) + p[1] * std::sin(theta))))
        .imag();
  }

  // Similar classes are required for the boundary values functions. The main
  // difference is that the number of components will now be 4 because those
  // functions will be applied to our space of skeletons unknowns via
  // VectorTools::interpolate_boundary_values. This space has 4 components,
  // because the skeleton unknowns on faces for the velocity field are scalars
  // from the definition $\hat{u}_n = \mathbf{u} \cdot n$ and there are the real
  // and imaginary part of both fields.

  template <int dim>
  class BoundaryValues_p_real : public Function<dim>
  {
  public:
    BoundaryValues_p_real(const double wavenumber, const double theta)
      : Function<dim>(4)
      , wavenumber(wavenumber)
      , theta(theta)
    {}
    virtual double value(const Point<dim>  &p,
                         const unsigned int component) const override;

  private:
    double wavenumber;
    double theta;
  };

  template <int dim>
  double BoundaryValues_p_real<dim>::value(
    const Point<dim>                   &p,
    [[maybe_unused]] const unsigned int component) const
  {
    constexpr std::complex<double> imag(0., 1.);

    return std::exp(-imag * wavenumber * p[1] * std::sin(theta)).real();
  }

  template <int dim>
  class BoundaryValues_p_imag : public Function<dim>
  {
  public:
    BoundaryValues_p_imag(const double wavenumber, const double theta)
      : Function<dim>(4)
      , wavenumber(wavenumber)
      , theta(theta)
    {}
    virtual double value(const Point<dim>  &p,
                         const unsigned int component) const override;

  private:
    double wavenumber;
    double theta;
  };

  template <int dim>
  double BoundaryValues_p_imag<dim>::value(
    const Point<dim>                   &p,
    [[maybe_unused]] const unsigned int component) const
  {
    constexpr std::complex<double> imag(0., 1.);

    return std::exp(-imag * wavenumber * p[1] * std::sin(theta)).imag();
  }
  template <int dim>
  class BoundaryValues_u_real : public Function<dim>
  {
  public:
    BoundaryValues_u_real(const double wavenumber, const double theta)
      : Function<dim>(4)
      , wavenumber(wavenumber)
      , theta(theta)
    {}
    virtual double value(const Point<dim>  &p,
                         const unsigned int component) const override;

  private:
    double wavenumber;
    double theta;
  };

  template <int dim>
  double BoundaryValues_u_real<dim>::value(
    const Point<dim>                   &p,
    [[maybe_unused]] const unsigned int component) const
  {
    constexpr std::complex<double> imag(0., 1.);
    return -1 * (std::sin(theta) *
                 std::exp(-imag * wavenumber * p[0] * std::cos(theta)))
                  .real();
  }
  template <int dim>
  class BoundaryValues_u_imag : public Function<dim>
  {
  public:
    BoundaryValues_u_imag(const double wavenumber, const double theta)
      : Function<dim>(4)
      , wavenumber(wavenumber)
      , theta(theta)
    {}
    virtual double value(const Point<dim>  &p,
                         const unsigned int component) const override;

  private:
    double wavenumber;
    double theta;
  };

  template <int dim>
  double BoundaryValues_u_imag<dim>::value(
    const Point<dim>                   &p,
    [[maybe_unused]] const unsigned int component) const
  {
    constexpr std::complex<double> imag(0., 1.);

    return -1 * (std::sin(theta) *
                 std::exp(-imag * wavenumber * p[0] * std::cos(theta)))
                  .imag();
  }

  // @sect3{The <code>DPGHelmholtz</code> class declaration}

  // Next let's declare the main class of this program. The main difference from
  // other examples lies in the fact that we rely on multiple DoFHandler and
  // FESystem. The DoFHandlers that we rely on are the following:
  // - The <code>dof_handler_trial_interior</code> is for the unknowns in the
  // interior of the cells;
  // - The <code>dof_handler_trial_skeleton</code> is for the unknowns in the
  // skeleton;
  // - The <code>dof_handler_test</code> is for the test functions. Although we
  // do not use the unknowns associated with this DoFHandler, it enables us to
  // evaluate the test function we will use in DPG.

  // The same applies for the three FESystem:
  // <code>fe_system_trial_interior</code>,
  // <code>fe_system_trial_skeleton</code> and <code>fe_system_test</code>. In
  // each one of these, we will store the relevant finite element space in the
  // same order to avoid confusion. The first component will therefore always be
  // related to the real part of the velocity, the second component to the its
  // imaginary part, the third component to the real part of the pressure and
  // the fourth component to its imaginary part.

  template <int dim>
  class DPGHelmholtz
  {
  public:
    // The constructor takes as an argument the <code>degree</code> of the
    // trial space as well as the  <code>delta_degree</code> between the trial
    // space and the test space which is necessary to construct the DPG problem.
    // The <code>delta_degree</code> must be at least 1 to ensure that the DPG
    // method is well posed. The parameter <code>theta</code> determines the
    // angle of the incident plane wave. That angle must be in the closed
    // interval $0$ and
    // $\frac{\pi}{2}$ for the boundary conditions to make sense. Those
    // restrictions are asserted in the constructor.

    DPGHelmholtz(unsigned int degree,
                 unsigned int delta_degree,
                 double       wavenumber,
                 double       theta);
    void run();

  private:
    // The <code>setup_system</code> function initializes the three DoFHandlers,
    // the system matrix and right-hand side and establishes the boundary
    // conditions that rely on AffineConstraints (e.g. the Dirichlet and Neumann
    // boundary conditions).

    void setup_system();

    // The <code>assemble_system</code> assembles both the right-hand side and
    // the system matrix. This function is used twice per problem and it has
    // two functions.
    // - When <code>solve_interior = false</code> the system is assembled and is
    // locally condensed such that the resulting system only contains the
    // skeleton unknowns. This is achieved by local condensation.
    // - When <code>solve_interior = true</code> the system is assembled and the
    // skeleton degrees of freedom solution is used to
    // reconstruct the interior solution.

    void assemble_system(bool solve_interior = false);

    // <code>solve_linear_system_skeleton</code> solves the linear system of
    // equations for the skeleton degree of freedom.
    void solve_linear_system_skeleton();

    // <code>refine_grid</code> refines the mesh uniformly.
    void refine_grid(unsigned int cycle);

    // <code>output_results</code> write the skeleton and the interior unknowns
    // into two different VTK file.
    void output_results(unsigned int cycle);

    // <code> calculate_L2_error</code> calculates the $L^2$ norm of the error
    // using the analytical solution.
    void calculate_L2_error();

    // We define multiple variables that are used throughout the class
    // starting with the triangulation.
    Triangulation<dim> triangulation;

    // The variables for the interior unknowns.
    const FESystem<dim> fe_trial_interior;
    DoFHandler<dim>     dof_handler_trial_interior;
    Vector<double>      solution_interior;

    // The variables for the skeleton and, consequently, the linear system.
    const FESystem<dim>       fe_trial_skeleton;
    DoFHandler<dim>           dof_handler_trial_skeleton;
    Vector<double>            solution_skeleton;
    Vector<double>            system_rhs;
    SparsityPattern           sparsity_pattern;
    SparseMatrix<double>      system_matrix;
    AffineConstraints<double> constraints;

    // The variables for the test spaces.
    const FESystem<dim> fe_test;
    DoFHandler<dim>     dof_handler_test;

    // A container for the $L^2$ error and other related quantities.
    ConvergenceTable error_table;

    // The coefficients which are used to define the plane wave problem.
    const double wavenumber;
    const double theta;

    // Some extractors that will be used at multiple places in the class to
    // select the relevant components for the calculation. Those are created
    // here to avoid repetition throughout the class.
    const FEValuesExtractors::Vector extractor_u_real;
    const FEValuesExtractors::Vector extractor_u_imag;
    const FEValuesExtractors::Scalar extractor_p_real;
    const FEValuesExtractors::Scalar extractor_p_imag;

    // However, the skeleton space does not have the same number of components
    // as the interior or the test space because the space $H^{-1/2}$ related to
    // the velocity field is a scalar field. Consequently, we need to define
    // additional extractor specifically for the skeleton spaces.
    const FEValuesExtractors::Scalar extractor_u_hat_real;
    const FEValuesExtractors::Scalar extractor_u_hat_imag;
    const FEValuesExtractors::Scalar extractor_p_hat_real;
    const FEValuesExtractors::Scalar extractor_p_hat_imag;
  };

  // @sect3{DPGHelmholtz Constructor}
  // In the constructor, we assign the relevant finite element to each FESystem
  // following the nomenclature described above in the class description:
  // - <code>fe_system_trial_interior</code> contains $\Re(\mathbf{u})$,
  // $\Im(\mathbf{u})$, $\Re(p^*)$, $\Im(p^*)$ ;
  // - <code>fe_system_trial_skeleton</code> contains $\Re(\hat{u}_n)$,
  // $\Im(\hat{u}_n)$, $\Re(\hat{p}^*)$, $\Im(\hat{p}^*)$ ;
  // - <code>fe_system_test</code> contains $\Re(\mathbf{v})$,
  // $\Im(\mathbf{v})$, $\Re(q)$, $\Im(q)$.

  // Note that the Q elements have a higher degree than the others because their
  // numbering start at 1 instead of 0.

  template <int dim>
  DPGHelmholtz<dim>::DPGHelmholtz(const unsigned int degree,
                                  const unsigned int delta_degree,
                                  double             wavenumber,
                                  double             theta)
    : fe_trial_interior(FE_DGQ<dim>(degree) ^ dim,
                        FE_DGQ<dim>(degree) ^ dim,
                        FE_DGQ<dim>(degree),
                        FE_DGQ<dim>(degree))
    , dof_handler_trial_interior(triangulation)
    , fe_trial_skeleton(FE_FaceQ<dim>(degree),
                        FE_FaceQ<dim>(degree),
                        FE_TraceQ<dim>(degree + 1),
                        FE_TraceQ<dim>(degree + 1))
    , dof_handler_trial_skeleton(triangulation)
    , fe_test(FE_RaviartThomas<dim>(degree + delta_degree),
              FE_RaviartThomas<dim>(degree + delta_degree),
              FE_Q<dim>(degree + delta_degree + 1),
              FE_Q<dim>(degree + delta_degree + 1))
    , dof_handler_test(triangulation)
    , wavenumber(wavenumber)
    , theta(theta)

    // We also initialize the FEValuesExtractors that will be used according to
    // our FESystems.
    , extractor_u_real(0)
    , extractor_u_imag(dim)
    , extractor_p_real(2 * dim)
    , extractor_p_imag(2 * dim + 1)

    , extractor_u_hat_real(0)
    , extractor_u_hat_imag(1)
    , extractor_p_hat_real(2)
    , extractor_p_hat_imag(3)

  // In the constructor we also put some assertions to  check if everything is
  // correctly defined for our problem. The first assertion checks that the
  // dimension is 2 because the problem is not defined in 3D. We also verify
  // that the delta_degree variable is at least 1 since the degree of the test
  // space must be at least one degree higher than the trial space. Finally, we
  // check that the wavenumber is positive since it is the magnitude of the wave
  // vector and that the angle theta is in the interval [0, $\pi$/2] because, as
  // stated above,  // other angles would not be compatible with the current
  // boundary definitions.
  {
    AssertThrow(dim == 2,
                ExcMessage("The step-100 example only works for dim==2"));

    AssertThrow(delta_degree >= 1,
                ExcMessage("The delta_degree needs to be at least 1."));

    AssertThrow(wavenumber > 0, ExcMessage("The wavenumber must be positive."));

    AssertThrow(theta >= 0 && theta <= M_PI / 2,
                ExcMessage(
                  "The angle theta must be in the interval [0, pi/2]."));
  }

  // @sect3{DPG::setup_system}
  // This function is similar to the other examples. The main difference lies in
  // the fact that we need to set up multiple DOFHandlers for the interior, the
  // skeleton and the test space.
  template <int dim>
  void DPGHelmholtz<dim>::setup_system()
  {
    dof_handler_trial_skeleton.distribute_dofs(fe_trial_skeleton);
    dof_handler_trial_interior.distribute_dofs(fe_trial_interior);
    dof_handler_test.distribute_dofs(fe_test);

    // We print the number of degrees of freedoms for each of the DoFHandler as
    // well as adding this information to the ConvergenceTable.
    std::cout << std::endl
              << "Number of dofs for the interior: "
              << dof_handler_trial_interior.n_dofs() << std::endl;

    error_table.add_value("dofs_interior", dof_handler_trial_interior.n_dofs());

    std::cout << "Number of dofs for the skeleton: "
              << dof_handler_trial_skeleton.n_dofs() << std::endl;

    error_table.add_value("dofs_skeleton", dof_handler_trial_skeleton.n_dofs());

    std::cout << "Number of dofs for the test space: "
              << dof_handler_test.n_dofs() << std::endl;

    error_table.add_value("dofs_test", dof_handler_test.n_dofs());


    constraints.clear();

    DoFTools::make_hanging_node_constraints(dof_handler_trial_skeleton,
                                            constraints);


    // We need to specify different boundary conditions for the four unknowns on
    // the faces. Therefore, we instantiate the functions that are used to
    // establish these boundary conditions.
    BoundaryValues_p_real<dim> p_real(wavenumber, theta);
    BoundaryValues_p_imag<dim> p_imag(wavenumber, theta);
    BoundaryValues_u_real<dim> u_real(wavenumber, theta);
    BoundaryValues_u_imag<dim> u_imag(wavenumber, theta);

    // Using the functions and the FEValuesExtractors, we impose the four
    // different constraints. As stated in the problem description, we first
    // impose a Dirichlet boundary condition on the pressure field for the left
    // boundary (id=0).
    VectorTools::interpolate_boundary_values(dof_handler_trial_skeleton,
                                             0,
                                             p_real,
                                             constraints,
                                             fe_trial_skeleton.component_mask(
                                               extractor_p_hat_real));
    VectorTools::interpolate_boundary_values(dof_handler_trial_skeleton,
                                             0,
                                             p_imag,
                                             constraints,
                                             fe_trial_skeleton.component_mask(
                                               extractor_p_hat_imag));

    // Then we impose a Neumann boundary condition on pressure by applying a
    // Dirichlet on the pressure "flux", which is the normal velocity field on
    // the bottom boundary (id=2).
    VectorTools::interpolate_boundary_values(dof_handler_trial_skeleton,
                                             2,
                                             u_real,
                                             constraints,
                                             fe_trial_skeleton.component_mask(
                                               extractor_u_hat_real));
    VectorTools::interpolate_boundary_values(dof_handler_trial_skeleton,
                                             2,
                                             u_imag,
                                             constraints,
                                             fe_trial_skeleton.component_mask(
                                               extractor_u_hat_imag));

    // The Robin boundary conditions do not require constraints and will be
    // dealt with in the assembly so we can close our constraints.
    constraints.close();

    // The linear system that we form is only related to the skeleton degrees of
    // freedom. We initialize all the necessary variables.
    solution_skeleton.reinit(dof_handler_trial_skeleton.n_dofs());
    system_rhs.reinit(dof_handler_trial_skeleton.n_dofs());
    solution_interior.reinit(dof_handler_trial_interior.n_dofs());

    DynamicSparsityPattern dsp(dof_handler_trial_skeleton.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler_trial_skeleton,
                                    dsp,
                                    constraints,
                                    false);
    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
  };

  // @sect3{DPG::assemble_system}
  // Now comes the most significant part of the program where we assemble the
  // DPG system.
  template <int dim>
  void DPGHelmholtz<dim>::assemble_system(const bool solve_interior)
  {
    // We first define quadrature rules and related variables. Since the
    // quadrature formula should be the same for both trial and test FE and that
    // the test space have higher polynomial degree than the others by
    // construction, we use it to define the quadrature formula.
    const QGauss<dim>     quadrature_formula(fe_test.degree + 1);
    const QGauss<dim - 1> face_quadrature_formula(fe_test.degree + 1);
    const unsigned int    n_q_points      = quadrature_formula.size();
    const unsigned int    n_face_q_points = face_quadrature_formula.size();

    // We then create the corresponding FEValues and FEFaceValues objects. Note
    // that only the test space needs gradients because of the ultraweak
    // formulation. Similarly, because everything is on the same triangulation,
    // we only need to update the quadrature points and JxW values in one of the
    // spaces. Here we choose the trial space.
    FEValues<dim> fe_values_trial_interior(fe_trial_interior,
                                           quadrature_formula,
                                           update_values |
                                             update_quadrature_points |
                                             update_JxW_values);

    FEValues<dim> fe_values_test(fe_test,
                                 quadrature_formula,
                                 update_values | update_gradients);

    FEFaceValues<dim> fe_values_trial_skeleton(fe_trial_skeleton,
                                               face_quadrature_formula,
                                               update_values |
                                                 update_quadrature_points |
                                                 update_normal_vectors |
                                                 update_JxW_values);

    FEFaceValues<dim> fe_face_values_test(fe_test,
                                          face_quadrature_formula,
                                          update_values);


    // We also create all the relevant matrices and vector to build the DPG
    // system. To do so we first need the number of dofs per cell for each of
    // the finite element spaces.
    const unsigned int dofs_per_cell_test = fe_test.n_dofs_per_cell();
    const unsigned int dofs_per_cell_trial_interior =
      fe_trial_interior.n_dofs_per_cell();
    const unsigned int dofs_per_cell_trial_skeleton =
      fe_trial_skeleton.n_dofs_per_cell();

    // We will need the shapes
    // functions values, gradient and divergence at the quadrature points for
    // the trial and test spaces. To avoid the query of the FEValues at each
    // quadrature point, we will use containers that will store the desired
    // quantities beforehand.

    std::vector<Tensor<1, dim, std::complex<double>>> v(dofs_per_cell_test);
    std::vector<Tensor<1, dim, std::complex<double>>> v_conj(
      dofs_per_cell_test);
    std::vector<std::complex<double>> div_v(dofs_per_cell_test);
    std::vector<std::complex<double>> div_v_conj(dofs_per_cell_test);
    std::vector<std::complex<double>> q(dofs_per_cell_test);
    std::vector<std::complex<double>> q_conj(dofs_per_cell_test);
    std::vector<Tensor<1, dim, std::complex<double>>> grad_q(
      dofs_per_cell_test);
    std::vector<Tensor<1, dim, std::complex<double>>> grad_q_conj(
      dofs_per_cell_test);
    std::vector<std::complex<double>> v_face_n(dofs_per_cell_test);
    std::vector<std::complex<double>> v_face_n_conj(dofs_per_cell_test);
    std::vector<std::complex<double>> q_face(dofs_per_cell_test);
    std::vector<std::complex<double>> q_face_conj(dofs_per_cell_test);

    std::vector<Tensor<1, dim, std::complex<double>>> u(
      dofs_per_cell_trial_interior);
    std::vector<std::complex<double>> p(dofs_per_cell_trial_interior);

    std::vector<std::complex<double>> u_hat_n(dofs_per_cell_trial_skeleton);
    std::vector<std::complex<double>> u_hat_n_conj(
      dofs_per_cell_trial_skeleton);
    std::vector<std::complex<double>> p_hat(dofs_per_cell_trial_skeleton);
    std::vector<std::complex<double>> p_hat_conj(dofs_per_cell_trial_skeleton);

    // We create the DPG local matrices before condensation:
    LAPACKFullMatrix<double> G_matrix(dofs_per_cell_test, dofs_per_cell_test);
    LAPACKFullMatrix<double> B_matrix(dofs_per_cell_test,
                                      dofs_per_cell_trial_interior);
    LAPACKFullMatrix<double> B_hat_matrix(dofs_per_cell_test,
                                          dofs_per_cell_trial_skeleton);
    LAPACKFullMatrix<double> D_matrix(dofs_per_cell_trial_skeleton,
                                      dofs_per_cell_trial_skeleton);

    // We create the DPG local vectors:
    Vector<double> g_vector(dofs_per_cell_trial_skeleton);
    Vector<double> l_vector(dofs_per_cell_test);

    // We create the condensation matrices:
    LAPACKFullMatrix<double> M1_matrix(dofs_per_cell_trial_interior,
                                       dofs_per_cell_trial_interior);
    LAPACKFullMatrix<double> M2_matrix(dofs_per_cell_trial_interior,
                                       dofs_per_cell_trial_skeleton);
    LAPACKFullMatrix<double> M3_matrix(dofs_per_cell_trial_skeleton,
                                       dofs_per_cell_trial_skeleton);
    LAPACKFullMatrix<double> M4_matrix(dofs_per_cell_trial_interior,
                                       dofs_per_cell_test);
    LAPACKFullMatrix<double> M5_matrix(dofs_per_cell_trial_skeleton,
                                       dofs_per_cell_test);

    // During the calculation of matrix vector product, we require intermediary
    // matrices and vector that we also allocate here.
    LAPACKFullMatrix<double> tmp_matrix(dofs_per_cell_trial_skeleton,
                                        dofs_per_cell_trial_interior);

    LAPACKFullMatrix<double> tmp_matrix2(dofs_per_cell_trial_skeleton,
                                         dofs_per_cell_trial_skeleton);

    LAPACKFullMatrix<double> tmp_matrix3(dofs_per_cell_trial_skeleton,
                                         dofs_per_cell_test);

    Vector<double> tmp_vector(dofs_per_cell_trial_interior);

    // We create the cell matrix and the RHS that will be distributed in the
    // full system after the assembly along with the index’s mapping.
    FullMatrix<double> cell_matrix(dofs_per_cell_trial_skeleton,
                                   dofs_per_cell_trial_skeleton);
    Vector<double>     cell_skeleton_rhs(dofs_per_cell_trial_skeleton);
    std::vector<types::global_dof_index> local_dof_indices(
      dofs_per_cell_trial_skeleton);

    // Finally, when reconstructing the interior solution from the skeleton, we
    // require additional vectors that we allocate here.
    Vector<double> cell_interior_rhs(dofs_per_cell_trial_interior);
    Vector<double> cell_interior_solution(dofs_per_cell_trial_interior);
    Vector<double> cell_skeleton_solution(dofs_per_cell_trial_skeleton);

    // We also define the imaginary unit and
    // two complex constant that will be used during the following assembly.
    // Note that even if the system and matrix that we build are real, we still
    // make use of the standard library complex operation to facilitate some
    // computations as it is done in step-81. Also, because the phase velocity
    // is 1, the wavenumber $k$ is equivalent to the angular
    // frequency $omega$.
    constexpr std::complex<double> imag(0., 1.);
    const std::complex<double>     iomega      = imag * wavenumber;
    const std::complex<double>     iomega_conj = conj(iomega);

    // As it is standard we first loop over the cells of the triangulation. Here
    // we have the choice of the DoFHandler to perform this loop. We use the
    // DoFHandler associated with the trial space.
    for (const auto &cell : dof_handler_trial_interior.active_cell_iterators())
      {
        // We reinitialize the FEValues objects to the current cell.
        fe_values_trial_interior.reinit(cell);

        // We will also need to reinitialize the FEValues for the test
        // space and make sure that is the same cell as the one used for the
        // trial space.
        const typename DoFHandler<dim>::active_cell_iterator cell_test =
          cell->as_dof_handler_iterator(dof_handler_test);
        fe_values_test.reinit(cell_test);

        // Similarly, we reinitialize the FEValues for the trial space on the
        // skeleton, but this will not be used before we also loop on the cells
        // faces.
        const typename DoFHandler<dim>::active_cell_iterator cell_skeleton =
          cell->as_dof_handler_iterator(dof_handler_trial_skeleton);

        // We then reinitialize all the matrices that we are aggregating
        // information for each cell.
        G_matrix     = 0;
        B_matrix     = 0;
        B_hat_matrix = 0;
        D_matrix     = 0;
        g_vector     = 0;
        l_vector     = 0;

        // We also need to reinitialize the $M_1$ condensation matrix between
        // each iteration on cell to get rid of its inverse status.
        M1_matrix = 0;

        // We loop over all quadrature points of
        // our cell.
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {
            // To avoid unnecessary computation, we fill the shape values
            // containers for the real and imaginary parts of the velocity and
            // pressure fields.
            const double &JxW = fe_values_trial_interior.JxW(q_point);

            for (unsigned int k : fe_values_test.dof_indices())
              {
                v[k] =
                  fe_values_test[extractor_u_real].value(k, q_point) +
                  imag * fe_values_test[extractor_u_imag].value(k, q_point);
                v_conj[k] =
                  fe_values_test[extractor_u_real].value(k, q_point) -
                  imag * fe_values_test[extractor_u_imag].value(k, q_point);

                div_v[k] =
                  fe_values_test[extractor_u_real].divergence(k, q_point) +
                  imag *
                    fe_values_test[extractor_u_imag].divergence(k, q_point);
                div_v_conj[k] =
                  fe_values_test[extractor_u_real].divergence(k, q_point) -
                  imag *
                    fe_values_test[extractor_u_imag].divergence(k, q_point);

                q[k] =
                  fe_values_test[extractor_p_real].value(k, q_point) +
                  imag * fe_values_test[extractor_p_imag].value(k, q_point);
                q_conj[k] =
                  fe_values_test[extractor_p_real].value(k, q_point) -
                  imag * fe_values_test[extractor_p_imag].value(k, q_point);

                grad_q[k] =
                  fe_values_test[extractor_p_real].gradient(k, q_point) +
                  imag * fe_values_test[extractor_p_imag].gradient(k, q_point);
                grad_q_conj[k] =
                  fe_values_test[extractor_p_real].gradient(k, q_point) -
                  imag * fe_values_test[extractor_p_imag].gradient(k, q_point);
              }

            for (unsigned int k : fe_values_trial_interior.dof_indices())
              {
                u[k] =
                  fe_values_trial_interior[extractor_u_real].value(k, q_point) +
                  imag *
                    fe_values_trial_interior[extractor_u_imag].value(k,
                                                                     q_point);

                p[k] =
                  fe_values_trial_interior[extractor_p_real].value(k, q_point) +
                  imag *
                    fe_values_trial_interior[extractor_p_imag].value(k,
                                                                     q_point);
              }

            // We are now ready to loop on our test space indices "i".
            for (const auto i : fe_values_test.dof_indices())
              {
                // When building the different part of the matrices, we need to
                // know in what finite element the current dof is associated
                // with. The following command does exactly that for the test
                // space.
                const unsigned int current_element_test_i =
                  fe_test.system_to_base_index(i).first.first;

                // If we are in the test space for pressure, we build the load
                // vector. For our plane wave problem this source term is null.
                if ((current_element_test_i == 2) ||
                    (current_element_test_i == 3))
                  {
                    l_vector(i) += 0;
                  }

                // The only two matrices that have some value in the interior of
                // the cell is the $G$ and $B$ matrices. We
                // will first construct the Gram matrix, so we will loop over
                // test space dofs a second time.
                for (const auto j : fe_values_test.dof_indices())
                  {
                    const unsigned int current_element_test_j =
                      fe_test.system_to_base_index(j).first.first;

                    // Now, using the information of which finite element is
                    // owning the dof <code>i</code> and  <code>j</code>, we can
                    // assemble the desired terms.

                    // For example, if both <code>i</code> and <code>j</code>
                    // are in test space associated to the test functions
                    // $\mathbf{v}$, we build the terms $(\mathbf{v},
                    // \overline{\mathbf{v}})_{\Omega_h}
                    // +
                    // (\nabla \cdot \mathbf{v}, \overline{\nabla \cdot
                    // \mathbf{v}})_{\Omega_h} + (\overline{i\omega}\mathbf{v},
                    // i\omega \overline{ \mathbf{v}})_{\Omega_h}$ :
                    if (((current_element_test_i == 0) ||
                         (current_element_test_i == 1)) &&
                        ((current_element_test_j == 0) ||
                         (current_element_test_j == 1)))
                      {
                        G_matrix(i, j) +=
                          (((v[j] * v_conj[i]) + (div_v[j] * div_v_conj[i]) +
                            (iomega_conj * v[j] * iomega * v_conj[i])) *
                           JxW)
                            .real();
                      }
                    // If the dof <code>i</code> is in test function $mathbf{v}$
                    // and dof <code>j</code> in test function $q$ we build the
                    // terms $-(\nabla q, i\omega \overline{
                    // \mathbf{v}})_{\Omega_h} - (\overline{ i\omega} q,
                    // \overline{\nabla \cdot \mathbf{v}})_{\Omega_h}$ :
                    else if (((current_element_test_i == 0) ||
                              (current_element_test_i == 1)) &&
                             ((current_element_test_j == 2) ||
                              (current_element_test_j == 3)))
                      {
                        G_matrix(i, j) -=
                          (((grad_q[j] * iomega * v_conj[i]) +
                            (iomega_conj * q[j] * div_v_conj[i])) *
                           JxW)
                            .real();
                      }
                    // If the dof <code>i</code> is in test function $q$ and the
                    // dof <code>j</code> is in the test function $\mathbf{v}$,
                    // we build the terms $-(\overline{i\omega} \mathbf{v},
                    // \overline{\nabla q})_{\Omega_h} - (\nabla \cdot
                    // \mathbf{v}, i\omega \overline{ q})_{\Omega_h}$ :
                    else if (((current_element_test_i == 2) ||
                              (current_element_test_i == 3)) &&
                             ((current_element_test_j == 0) ||
                              (current_element_test_j == 1)))
                      {
                        G_matrix(i, j) -=
                          (((iomega_conj * v[j] * grad_q_conj[i]) +
                            (div_v[j] * iomega * q_conj[i])) *
                           JxW)
                            .real();
                      }
                    // Finally, if both in test functions are in $q$, we build
                    // the terms $(q, \overline{q})_{\Omega_h} + (\nabla q,
                    // \overline{\nabla q})_{\Omega_h} + (\overline{ i\omega} q,
                    // i\omega\overline{ q})_{\Omega_h}$:
                    else if (((current_element_test_i == 2) ||
                              (current_element_test_i == 3)) &&
                             ((current_element_test_j == 2) ||
                              (current_element_test_j == 3)))
                      {
                        G_matrix(i, j) +=
                          (((q[j] * q_conj[i]) + (grad_q[j] * grad_q_conj[i]) +
                            (iomega_conj * q[j] * iomega * q_conj[i])) *
                           JxW)
                            .real();
                      }
                  }

                // Now we will build the matrix $B$ associated
                // with the operator of our problem on the interior element. We
                // loop over trial space dofs on <code>j</code> and preform a
                // similar procedure as for the Gram matrix.
                for (const auto j : fe_values_trial_interior.dof_indices())
                  {
                    // Again, we get the information to map the dof to the
                    // right shape function.
                    const unsigned int current_element_trial_j =
                      fe_trial_interior.system_to_base_index(j).first.first;

                    // If dof <code>i</code> in test function $\mathbf{v}$ and
                    // dof <code>j</code> in trial function $\mathbf{u}$ we
                    // build the term $(i\omega \mathbf{u},
                    // \overline{\mathbf{v}})_{\Omega_h}$:
                    if (((current_element_test_i == 0) ||
                         (current_element_test_i == 1)) &&
                        ((current_element_trial_j == 0) ||
                         (current_element_trial_j == 1)))
                      {
                        B_matrix(i, j) +=
                          ((iomega * u[j] * v_conj[i]) * JxW).real();
                      }
                    // If dof <code>i</code> in test function $\mathbf{v}$ and
                    // dof <code>j</code> in trial function $p$ we build the
                    // term $ -(p^*, \overline{\nabla \cdot
                    // \mathbf{v}})_{\Omega_h}$:
                    else if (((current_element_test_i == 0) ||
                              (current_element_test_i == 1)) &&
                             ((current_element_trial_j == 2) ||
                              (current_element_trial_j == 3)))
                      {
                        B_matrix(i, j) -= ((p[j] * div_v_conj[i]) * JxW).real();
                      }

                    // If dof <code>i</code> in test function $q$ and dof
                    // <code>j</code> in trial function $\mathbf{u}$ we build
                    // the term $-(\mathbf{u}, \overline{\nabla q})_{\Omega_h}$:
                    else if (((current_element_test_i == 2) ||
                              (current_element_test_i == 3)) &&
                             ((current_element_trial_j == 0) ||
                              (current_element_trial_j == 1)))
                      {
                        B_matrix(i, j) -=
                          ((u[j] * grad_q_conj[i]) * JxW).real();
                      }
                    // If dof <code>i</code> in test function $q$ and dof
                    // <code>j</code> in trial function $p$ we build the term
                    // $(i\omega p^*, \overline{q})_{\Omega_h}$:
                    else if (((current_element_test_i == 2) ||
                              (current_element_test_i == 3)) &&
                             ((current_element_trial_j == 2) ||
                              (current_element_trial_j == 3)))
                      {
                        B_matrix(i, j) +=
                          ((iomega * p[j] * q_conj[i]) * JxW).real();
                      }
                  }
              }
          }
        // We now build the skeleton terms. Similarly, we choose to loop on the
        // skeleton trial space faces.
        for (const auto &face : cell_skeleton->face_iterators())
          {
            // We reinitialize the FEFaceValues objects to the current faces.
            fe_face_values_test.reinit(cell_test, face);
            fe_values_trial_skeleton.reinit(cell_skeleton, face);

            // In what follows we will need to add terms that are only relevant
            // for the Robin boundary conditions. We also get the face number
            // that will be used later on to determine the correct normal
            // orientation of the fluxes.
            const auto face_no             = cell->face_iterator_to_index(face);
            const auto current_boundary_id = face->boundary_id();

            // The Robin boundary conditions in our plane wave problem have a
            // factor $\frac{\mathbf{k} \cdot \mathbf{n}}{\omega}$. However,
            // $\omega = k c_s$ with $c_s=1$ and $\mathbf{k} \cdot \mathbf{n}$
            // is either $k\cos{\theta}$ for the right boundary or
            // $k\sin{\theta}$ for the top boundary because of our domain
            // geometry. The wavenumber simplifies and we are left with the
            // cosine and sine for the value of this factor, which we compute in
            // advance according to the boundary id. Note that, in our case,
            // the term is always real and  $\overline{\frac{k_n}{\omega}} =
            // \frac{k_n}{\omega}$ in what follows.
            const double kn_omega = (current_boundary_id == 1) ? cos(theta) :
                                    (current_boundary_id == 3) ? sin(theta) :
                                                                 1.;

            // As for the cell loop, we go over the quadrature points, but for
            // the face this time.
            for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point)
              {
                // Again, we assign the constant values and fill the shape
                // values containers but for the face variables this time.
                const auto &normal =
                  fe_values_trial_skeleton.normal_vector(q_point);
                const double &JxW_face = fe_values_trial_skeleton.JxW(q_point);

                for (unsigned int k : fe_face_values_test.dof_indices())
                  {
                    v_face_n[k] =
                      normal *
                      (fe_face_values_test[extractor_u_real].value(k, q_point) +
                       imag *
                         fe_face_values_test[extractor_u_imag].value(k,
                                                                     q_point));
                    v_face_n_conj[k] =
                      normal *
                      (fe_face_values_test[extractor_u_real].value(k, q_point) -
                       imag *
                         fe_face_values_test[extractor_u_imag].value(k,
                                                                     q_point));

                    q_face[k] =
                      fe_face_values_test[extractor_p_real].value(k, q_point) +
                      imag *
                        fe_face_values_test[extractor_p_imag].value(k, q_point);
                    q_face_conj[k] =
                      fe_face_values_test[extractor_p_real].value(k, q_point) -
                      imag *
                        fe_face_values_test[extractor_p_imag].value(k, q_point);
                  }
                for (unsigned int k : fe_values_trial_skeleton.dof_indices())
                  {
                    u_hat_n[k] =
                      fe_values_trial_skeleton[extractor_u_hat_real].value(
                        k, q_point) +
                      imag * fe_values_trial_skeleton[extractor_u_hat_imag]
                               .value(k, q_point);
                    u_hat_n_conj[k] =
                      fe_values_trial_skeleton[extractor_u_hat_real].value(
                        k, q_point) -
                      imag * fe_values_trial_skeleton[extractor_u_hat_imag]
                               .value(k, q_point);

                    p_hat[k] =
                      fe_values_trial_skeleton[extractor_p_hat_real].value(
                        k, q_point) +
                      imag * fe_values_trial_skeleton[extractor_p_hat_imag]
                               .value(k, q_point);
                    p_hat_conj[k] =
                      fe_values_trial_skeleton[extractor_p_hat_real].value(
                        k, q_point) -
                      imag * fe_values_trial_skeleton[extractor_p_hat_imag]
                               .value(k, q_point);
                  }

                // We loop over the test space face dofs with indices
                // <code>i</code>, create the face test basis functions and get
                // the information to map the index to the right shape function.
                for (const auto i : fe_face_values_test.dof_indices())
                  {
                    const unsigned int current_element_test_i =
                      fe_test.system_to_base_index(i).first.first;

                    // We then verify if we are at one of the two Robin
                    // boundaries. If so, we add additional terms to the Gram
                    // matrix using the same procedure as for the interior
                    // terms.
                    if (current_boundary_id == 1 || current_boundary_id == 3)
                      {
                        for (const auto j : fe_face_values_test.dof_indices())
                          {
                            const unsigned int current_element_test_j =
                              fe_test.system_to_base_index(j).first.first;


                            // We build $ \langle \mathbf{v} \cdot \mathbf{n},
                            //  \overline{\mathbf{v} \cdot \mathbf{n}}
                            //  \rangle_{\Gamma_1 \cup \Gamma_3}$ :
                            if (((current_element_test_i == 0) ||
                                 (current_element_test_i == 1)) &&
                                ((current_element_test_j == 0) ||
                                 (current_element_test_j == 1)))
                              {
                                G_matrix(i, j) +=
                                  (v_face_n[j] * v_face_n_conj[i] * JxW_face)
                                    .real();
                              }
                            //  We build $\langle
                            //  \overline{\frac{k_n}{\omega}}q,
                            //  \overline{\mathbf{v} \cdot \mathbf{n}}
                            //  \rangle_{\Gamma_1 \cup \Gamma_3}$ :
                            else if (((current_element_test_i == 0) ||
                                      (current_element_test_i == 1)) &&
                                     ((current_element_test_j == 2) ||
                                      (current_element_test_j == 3)))
                              {
                                G_matrix(i, j) += (kn_omega * q_face[j] *
                                                   v_face_n_conj[i] * JxW_face)
                                                    .real();
                              }
                            // We build $\langle \mathbf{v} \cdot \mathbf{n},
                            // \frac{k_n}{\omega}\overline{q} \rangle_{\Gamma_1
                            // \cup \Gamma_3}$ :
                            else if (((current_element_test_i == 2) ||
                                      (current_element_test_i == 3)) &&
                                     ((current_element_test_j == 0) ||
                                      (current_element_test_j == 1)))
                              {
                                G_matrix(i, j) += (v_face_n[j] * kn_omega *
                                                   q_face_conj[i] * JxW_face)
                                                    .real();
                              }
                            // We build $\langle \overline{\frac{k_n}{\omega}}q,
                            // \frac{k_n}{\omega}\overline{q} \rangle_{\Gamma_1
                            // \cup \Gamma_3}$ :
                            else if (((current_element_test_i == 2) ||
                                      (current_element_test_i == 3)) &&
                                     ((current_element_test_j == 2) ||
                                      (current_element_test_j == 3)))
                              {
                                G_matrix(i, j) +=
                                  (kn_omega * q_face[j] * kn_omega *
                                   q_face_conj[i] * JxW_face)
                                    .real();
                              }
                          }
                      }


                    // Again, same procedure, we construct the $\hat{B}$
                    // matrix and to do so, we loop over trial space dofs on
                    // <code>j</code>.
                    for (const auto j : fe_values_trial_skeleton.dof_indices())
                      {
                        const unsigned int current_element_trial_j =
                          fe_trial_skeleton.system_to_base_index(j).first.first;

                        // We build the term $\left\langle \hat{p}^*,
                        // \overline{\mathbf{v} \cdot \mathbf{n}}
                        // \right\rangle_{\partial \Omega_h}$:
                        if (((current_element_test_i == 0) ||
                             (current_element_test_i == 1)) &&
                            ((current_element_trial_j == 2) ||
                             (current_element_trial_j == 3)))
                          {
                            B_hat_matrix(i, j) +=
                              ((p_hat[j] * v_face_n_conj[i]) * JxW_face).real();
                          }

                        // We build the term $\left\langle \hat{u}_n,
                        // \overline{q} \right\rangle_{\partial \Omega_h}$ :
                        else if (((current_element_test_i == 2) ||
                                  (current_element_test_i == 3)) &&
                                 ((current_element_trial_j == 0) ||
                                  (current_element_trial_j == 1)))
                          {
                            // Now, the FE_FaceQ elements describe the trace of
                            // H(div) conforming elements. However, they do not
                            // have an orientation embedded in them. Therefore,
                            // to make sure that the flux that cross a face in a
                            // given cell is equal to the flux that crosses the
                            // same face in the adjacent cell (for which the
                            // normal is opposite), we add a factor that will be
                            // either 1 or -1. Since this type of element is
                            // only continuous across faces, we can create an
                            // orientation rule that only depends on a given
                            // face without the need for a global cell
                            // orientation. The rule used here is based on the
                            // cell index. The flux is always oriented from the
                            // lowest cell index to the highest cell index. At
                            // boundaries we follow the standard convention
                            // that the normal is oriented outward from the
                            // domain and the flux is align with the normal
                            // (i.e., it equals 1).
                            const int neighbor_cell_id =
                              face->at_boundary() ?
                                INT_MAX :
                                cell->neighbor(face_no)->index();
                            const double flux_orientation =
                              neighbor_cell_id > cell->index() ? 1. : -1.;

                            B_hat_matrix(i, j) +=
                              (flux_orientation * u_hat_n[j] * q_face_conj[i] *
                               JxW_face)
                                .real();
                          }
                      }
                  }

                // Finally, we will build the $\mathbf{D}$ matrix and the source
                // term vector $G$ for the Robin boundaries. We first
                // verify if the face is at one of them.
                if (current_boundary_id == 1 || current_boundary_id == 3)
                  {
                    // We adjust the flux orientation, but because we are at a
                    // boundary, it is always in the same direction as the
                    // normal.
                    const double flux_orientation = 1.;

                    // We loop over the trial space dofs on <code>i</code> and
                    // create the associated trial basis functions.
                    for (const auto i : fe_values_trial_skeleton.dof_indices())
                      {
                        const unsigned int current_element_trial_i =
                          fe_trial_skeleton.system_to_base_index(i).first.first;

                        // As for the load vector, we assemble the source terms
                        // for the Robin boundary even if it is null to
                        // illustrate the procedure.
                        if ((current_element_trial_i == 0) ||
                            (current_element_trial_i == 1))
                          {
                            g_vector(i) -=
                              ((0) * u_hat_n_conj[i]).real() * JxW_face;
                          }
                        else if ((current_element_trial_i == 2) ||
                                 (current_element_trial_i == 3))
                          {
                            g_vector(i) +=
                              ((0.) * kn_omega * p_hat_conj[i]).real() *
                              JxW_face;
                          }

                        // The matrix $\mathbf{D}$ put in relation the trial
                        // skeleton space with itself so we loop again on the
                        // trial space dofs on <code>j</code>.
                        for (const auto j :
                             fe_values_trial_skeleton.dof_indices())
                          {
                            const unsigned int current_element_trial_j =
                              fe_trial_skeleton.system_to_base_index(j)
                                .first.first;

                            // We build the term $- \langle \hat{u}_n,
                            // \overline{\hat{u}_n} \rangle_{\Gamma_1 \cup
                            // \Gamma_3}$:
                            if (((current_element_trial_i == 0) ||
                                 (current_element_trial_i == 1)) &&
                                ((current_element_trial_j == 0) ||
                                 (current_element_trial_j == 1)))
                              {
                                D_matrix(i, j) -=
                                  (flux_orientation * u_hat_n[j] *
                                   flux_orientation * u_hat_n_conj[i] *
                                   JxW_face)
                                    .real();
                              }
                            // We build the term $ \langle \hat{u}_n,
                            // \overline{\frac{k_n}{\omega} \hat{p}^*}
                            // \rangle_{\Gamma_1 \cup \Gamma_3}$:
                            else if (((current_element_trial_i == 0) ||
                                      (current_element_trial_i == 1)) &&
                                     ((current_element_trial_j == 2) ||
                                      (current_element_trial_j == 3)))
                              {
                                D_matrix(i, j) +=
                                  (kn_omega * p_hat[j] * flux_orientation *
                                   u_hat_n_conj[i] * JxW_face)
                                    .real();
                              }
                            // We build the term $\langle \frac{k_n}{\omega}
                            // \hat{p}^*, \overline{\hat{u}_n} \rangle_{\Gamma_1
                            // \cup \Gamma_3}$:
                            else if (((current_element_trial_i == 2) ||
                                      (current_element_trial_i == 3)) &&
                                     ((current_element_trial_j == 0) ||
                                      (current_element_trial_j == 1)))
                              {
                                D_matrix(i, j) +=
                                  (flux_orientation * u_hat_n[j] * kn_omega *
                                   p_hat_conj[i] * JxW_face)
                                    .real();
                              }
                            // We build the term $- \langle \frac{k_n}{\omega}
                            // \hat{p}^*, \overline{\frac{k_n}{\omega}
                            // \hat{p}^*} \rangle_{\Gamma_1 \cup \Gamma_3}$:
                            else if (((current_element_trial_i == 2) ||
                                      (current_element_trial_i == 3)) &&
                                     ((current_element_trial_j == 2) ||
                                      (current_element_trial_j == 3)))
                              {
                                D_matrix(i, j) -=
                                  (kn_omega * p_hat[j] * kn_omega *
                                   p_hat_conj[i] * JxW_face)
                                    .real();
                              }
                          }
                      }
                  }
              }
          }
        // Finally, after having assembled all the matrices and vectors, we
        // build the condensed version of the system.

        // We only need the inverse of the Gram matrix $G$, so we
        // invert it.
        G_matrix.invert();

        // We construct $M_4 = B^\dagger G^{-1}$ and $M_5 = \hat{B}^\dagger
        // G^{-1}$ with it:
        B_matrix.Tmmult(M4_matrix, G_matrix);
        B_hat_matrix.Tmmult(M5_matrix, G_matrix);

        // Then using $M_4$ we compute the condensed matrix $M_1 = B^\dagger
        // G^{-1} B$ and $M_2 = B^\dagger G^{-1} \hat{B}$:
        M4_matrix.mmult(M1_matrix, B_matrix);
        M4_matrix.mmult(M2_matrix, B_hat_matrix);

        // We also compute the matrix $M_3 = \hat{B}^\dagger G^{-1} \hat{B}$
        // using $M_5$ and then subtract $D$:
        M5_matrix.mmult(M3_matrix, B_hat_matrix);
        M3_matrix.add(-1.0, D_matrix);

        // Finally, as for the $G$ matrix, we invert the $M_1$
        // matrix:
        M1_matrix.invert();

        // If the flag solves interior is set to true, we have already the
        // solution on the skeleton and only need to perform $u_h = M_1^{-1}
        // (M_4 l - M_2 \hat{u}_h)$ on each cell.
        if (solve_interior)
          {
            // We first get the solution vector for this cell.
            cell_skeleton->get_dof_values(solution_skeleton,
                                          cell_skeleton_solution);

            // Then we do the matrix-vector product to obtain the interior
            // unknowns.
            M2_matrix.vmult(tmp_vector, cell_skeleton_solution);
            M4_matrix.vmult(cell_interior_rhs, l_vector);
            cell_interior_rhs -= tmp_vector;
            M1_matrix.vmult(cell_interior_solution, cell_interior_rhs);

            // Finally, we map the cell interior solution to the global interior
            // solution.
            cell->distribute_local_to_global(cell_interior_solution,
                                             solution_interior);
          }
        // If the flag solves interior is set to false, we have to compute the
        // local matrix and the local RHS for the condensed system.
        else
          {
            // So the cell matrix is obtained with the formula $(M_3 -
            // M_2^\dagger M_1^{-1} M_2)$:
            M2_matrix.Tmmult(tmp_matrix, M1_matrix);
            tmp_matrix.mmult(tmp_matrix2, M2_matrix);
            tmp_matrix2.add(-1.0, M3_matrix);
            tmp_matrix2 *= -1.0;
            // This line is used to convert the LAPACK matrix to a full matrix
            // so we can perform the distribution to the global system below.
            cell_matrix = tmp_matrix2;

            // Then we compute the cell RHS using $(M_5 -
            // M_2^\dagger M_1^{-1} M_4)l -
            // G$.
            tmp_matrix.mmult(tmp_matrix3, M4_matrix);
            M5_matrix.add(-1.0, tmp_matrix3);
            M5_matrix.vmult(cell_skeleton_rhs, l_vector);
            cell_skeleton_rhs -= g_vector;

            // Map to global matrix
            cell_skeleton->get_dof_indices(local_dof_indices);
            constraints.distribute_local_to_global(cell_matrix,
                                                   cell_skeleton_rhs,
                                                   local_dof_indices,
                                                   system_matrix,
                                                   system_rhs);
          }
      }
  }

  // @sect3{DPG::solve}
  // This function is in charge of solving the linear system assembled and has
  // nothing specific to DPG per se. Nonetheless, the method allows us to use
  // the Conjugate Gradient iterative solver . Note that because we do not have
  // any preconditioner, the number of iterations can be quite high. For
  // simplicity, we put a high upper limit on the number of iterations, but in
  // practice one would want to change this function to have a more robust
  // solver. The tolerance for the convergence here is defined proportional to
  // the $L^2$ norm of the RHS vector so the stopping criterion is scaled aware.
  template <int dim>
  void DPGHelmholtz<dim>::solve_linear_system_skeleton()
  {
    std::cout << std::endl << "Solving the DPG system..." << std::endl;

    SolverControl solver_control(100000, 1e-10 * system_rhs.l2_norm());
    SolverCG<Vector<double>> solver(solver_control);
    solver.solve(system_matrix,
                 solution_skeleton,
                 system_rhs,
                 PreconditionIdentity());
    constraints.distribute(solution_skeleton);

    std::cout << "   " << solver_control.last_step()
              << " CG iterations needed to obtain convergence. \n"
              << std::endl;

    error_table.add_value("n_iter", solver_control.last_step());
  }

  // @sect3{DPG::output_results}
  // This function is also quite standard and is in charge of outputting the
  // solution to a file that can be visualized with Paraview or Visit (VTK
  // format). However, because the skeleton solution lives only on the mesh
  // faces, we made use of the DataOutFaces class to output it properly.
  template <int dim>
  void DPGHelmholtz<dim>::output_results(const unsigned int cycle)
  {
    // We first attach the DoFHandler to the DataOut object for the interior
    // solution, define the names of the solution components, and specify how to
    // interpret each component (as part of a vector for the velocity real and
    // imaginary parts, and as scalars for the pressure real and imaginary
    // parts).
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler_trial_interior);

    std::vector<std::string> solution_interior_names(dim, "velocity_real");
    for (unsigned int i = 0; i < dim; ++i)
      {
        solution_interior_names.emplace_back("velocity_imag");
      }
    solution_interior_names.emplace_back("pressure_real");
    solution_interior_names.emplace_back("pressure_imag");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);

    // Once the setup is done, we add the solution vector to the DataOut object,
    // build the patches for visualization, and write the output.
    data_out.add_data_vector(solution_interior,
                             solution_interior_names,
                             DataOut<dim>::type_automatic,
                             data_component_interpretation);

    data_out.build_patches(fe_trial_interior.degree);

    std::ofstream output("solution_planewave_square-" + std::to_string(cycle) +
                         ".vtk");
    data_out.write_vtk(output);

    // Then we do the same for the skeleton solution. The main difference here,
    // beside the use of DataOutFaces, is that every component are scalars.
    DataOutFaces<dim> data_out_faces(false);
    data_out_faces.attach_dof_handler(dof_handler_trial_skeleton);

    std::vector<std::string> solution_skeleton_names(1, "velocity_hat_real");
    solution_skeleton_names.emplace_back("velocity_hat_imag");
    solution_skeleton_names.emplace_back("pressure_hat_real");
    solution_skeleton_names.emplace_back("pressure_hat_imag");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation_skeleton(
        4, DataComponentInterpretation::component_is_scalar);

    data_out_faces.add_data_vector(solution_skeleton,
                                   solution_skeleton_names,
                                   DataOutFaces<dim>::type_automatic,
                                   data_component_interpretation_skeleton);

    data_out_faces.build_patches(fe_trial_skeleton.degree);
    std::ofstream output_face("solution_face_planewave_square-" +
                              std::to_string(cycle) + ".vtk");
    data_out_faces.write_vtk(output_face);
  }

  // @sect3{DPG::calculate_error}
  // In this function, we compute the $L^2$ error of each component of our
  // solution, i.e., for the real and imaginary parts of the velocity and
  // pressure for both the interior and skeleton solutions. Because we want to
  // have all the error for all the different components of our solution vectors
  // separately, we cannot use the VectorTools::integrate_difference function
  // directly. Instead, we will perform the computation "by hand" by looping
  // over all the cells and faces, interpolating the solution at the quadrature
  // points, and computing the error with respect to the analytical solution at
  // those points.
  template <int dim>
  void DPGHelmholtz<dim>::calculate_L2_error()
  {
    // We first need to instantiate the FEValues and FEFaceValues objects as it
    // has been done during the assembly to manage the solution when we loop
    // over the cells and faces.
    QGauss<dim>           quadrature_formula(fe_test.degree + 1);
    FEValues<dim>         fe_values_trial_interior(fe_trial_interior,
                                           quadrature_formula,
                                           update_values |
                                             update_quadrature_points |
                                             update_JxW_values);
    const QGauss<dim - 1> face_quadrature_formula(fe_test.degree + 1);
    FEFaceValues<dim>     fe_values_trial_skeleton(fe_trial_skeleton,
                                               face_quadrature_formula,
                                               update_values |
                                                 update_quadrature_points |
                                                 update_normal_vectors |
                                                 update_JxW_values);

    const unsigned int n_q_points      = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    // We create variables that will store all the integration result we are
    // interested in.
    double L2_error_p_real     = 0;
    double L2_error_p_imag     = 0;
    double L2_error_p_hat_real = 0;
    double L2_error_p_hat_imag = 0;
    double L2_error_u_real     = 0;
    double L2_error_u_imag     = 0;
    double L2_error_u_hat_real = 0;
    double L2_error_u_hat_imag = 0;

    // To compute the error for the velocity field on the faces we will need to
    // compute the normal component of our analytical solution on those faces.
    // Here we create a variable to store the result of those scalar product.
    double u_hat_n_analytical_real = 0.;
    double u_hat_n_analytical_imag = 0.;

    // When looping on each cell or face, we will extract the different field
    // solution obtain numerically. The containers used to store the
    // interpolated solution at the quadrature points are declared below.
    std::vector<Tensor<1, dim>> local_u_real(n_q_points);
    std::vector<Tensor<1, dim>> local_u_imag(n_q_points);
    std::vector<double>         local_p_real(n_q_points);
    std::vector<double>         local_p_imag(n_q_points);
    std::vector<double>         local_u_hat_real(n_face_q_points);
    std::vector<double>         local_u_hat_imag(n_face_q_points);
    std::vector<double>         local_p_hat_real(n_face_q_points);
    std::vector<double>         local_p_hat_imag(n_face_q_points);

    // Here we initialize analytical solution object from the previously defined
    // classes.
    AnalyticalSolution_p_real<dim> analytical_solution_p_real(wavenumber,
                                                              theta);
    AnalyticalSolution_p_imag<dim> analytical_solution_p_imag(wavenumber,
                                                              theta);
    AnalyticalSolution_u_real<dim> analytical_solution_u_real(wavenumber,
                                                              theta);
    AnalyticalSolution_u_imag<dim> analytical_solution_u_imag(wavenumber,
                                                              theta);


    // We first loop over all the cells of our mesh and extract the solution in
    // the interior.
    for (const auto &cell : dof_handler_trial_interior.active_cell_iterators())
      {
        fe_values_trial_interior.reinit(cell);
        const typename DoFHandler<dim>::active_cell_iterator cell_skeleton =
          cell->as_dof_handler_iterator(dof_handler_trial_skeleton);

        fe_values_trial_interior[extractor_u_real].get_function_values(
          solution_interior, local_u_real);
        fe_values_trial_interior[extractor_u_imag].get_function_values(
          solution_interior, local_u_imag);
        fe_values_trial_interior[extractor_p_real].get_function_values(
          solution_interior, local_p_real);
        fe_values_trial_interior[extractor_p_imag].get_function_values(
          solution_interior, local_p_imag);

        // Then we loop on the cell quadrature points because it is where we
        // will compute the error.
        const auto &quadrature_points =
          fe_values_trial_interior.get_quadrature_points();

        for (const unsigned int q_index :
             fe_values_trial_interior.quadrature_point_indices())
          {
            const double JxW      = fe_values_trial_interior.JxW(q_index);
            const auto  &position = quadrature_points[q_index];


            for (unsigned int i = 0; i < dim; ++i)
              {
                L2_error_u_real +=
                  pow((local_u_real[q_index][i] -
                       analytical_solution_u_real.value(position, i)),
                      2) *
                  JxW;
                L2_error_u_imag +=
                  pow((local_u_imag[q_index][i] -
                       analytical_solution_u_imag.value(position, i)),
                      2) *
                  JxW;
              }

            L2_error_p_real +=
              pow((local_p_real[q_index] -
                   analytical_solution_p_real.value(position, 0)),
                  2) *
              JxW;

            L2_error_p_imag +=
              pow((local_p_imag[q_index] -
                   analytical_solution_p_imag.value(position, 0)),
                  2) *
              JxW;
          }

        // Then we loop on the face for the skeleton solution error and use the
        // same approach.
        for (const auto &face : cell->face_iterators())
          {
            fe_values_trial_skeleton.reinit(cell_skeleton, face);
            const auto face_no = cell_skeleton->face_iterator_to_index(face);

            fe_values_trial_skeleton[extractor_u_hat_real].get_function_values(
              solution_skeleton, local_u_hat_real);
            fe_values_trial_skeleton[extractor_u_hat_imag].get_function_values(
              solution_skeleton, local_u_hat_imag);
            fe_values_trial_skeleton[extractor_p_hat_real].get_function_values(
              solution_skeleton, local_p_hat_real);
            fe_values_trial_skeleton[extractor_p_hat_imag].get_function_values(
              solution_skeleton, local_p_hat_imag);

            const auto &quadrature_points =
              fe_values_trial_skeleton.get_quadrature_points();

            for (const unsigned int &q_index :
                 fe_values_trial_skeleton.quadrature_point_indices())
              {
                const double JxW      = fe_values_trial_skeleton.JxW(q_index);
                const auto  &position = quadrature_points[q_index];
                const auto  &normal =
                  fe_values_trial_skeleton.normal_vector(q_index);

                // Because we are looping over cells and cells share faces in
                // the interior we add a condition to only integrate once per
                // face. We use a similar idea than to define the flux
                // orientation during the assembly. Each face is only integrated
                // when we are at the cell with the lowest index among the two
                // cells sharing the face. Boundary faces are always integrated.
                int neighbor_cell_id = face->at_boundary() ?
                                         INT_MAX :
                                         cell->neighbor(face_no)->index();
                if (neighbor_cell_id < cell->index())
                  {
                    continue;
                  }

                // Has mentioned earlier, to compute the error for the
                // $\hat{u}_n$ we need to do the scalar product of the
                // analytical solution with the normal vector at the face
                // quadrature points. We also take the absolute of the values
                // for both the numerical and analytical solutions because the
                // sign only depends on the orientation of flux from the
                // convention with respect to the normal vector and therefore
                // change from one face to another.
                u_hat_n_analytical_real = 0.;
                u_hat_n_analytical_imag = 0.;
                for (unsigned int i = 0; i < dim; ++i)
                  {
                    u_hat_n_analytical_real +=
                      normal[i] * analytical_solution_u_real.value(position, i);
                    u_hat_n_analytical_imag +=
                      normal[i] * analytical_solution_u_imag.value(position, i);
                  }

                L2_error_u_hat_real += pow(abs(local_u_hat_real[q_index]) -
                                             abs(u_hat_n_analytical_real),
                                           2) *
                                       JxW;
                L2_error_u_hat_imag += pow(abs(local_u_hat_imag[q_index]) -
                                             abs(u_hat_n_analytical_imag),
                                           2) *
                                       JxW;

                L2_error_p_hat_real +=
                  pow((local_p_hat_real[q_index] -
                       analytical_solution_p_real.value(position, 0)),
                      2) *
                  JxW;
                L2_error_p_hat_imag +=
                  pow((local_p_hat_imag[q_index] -
                       analytical_solution_p_imag.value(position, 0)),
                      2) *
                  JxW;
              }
          }
      }

    // Finally, we output the results to the terminal and store the errors in
    // the <code>error_table</code>.
    std::cout << "Velocity real part L2 error is : "
              << std::sqrt(L2_error_u_real) << std::endl;
    std::cout << "Velocity imag part L2 error is : "
              << std::sqrt(L2_error_u_imag) << std::endl;
    std::cout << "Pressure real part L2 error is : "
              << std::sqrt(L2_error_p_real) << std::endl;
    std::cout << "Pressure imag part L2 error is : "
              << std::sqrt(L2_error_p_imag) << std::endl;
    std::cout << "Velocity skeleton real part L2 error is : "
              << std::sqrt(L2_error_u_hat_real) << std::endl;
    std::cout << "Velocity skeleton imag part L2 error is : "
              << std::sqrt(L2_error_u_hat_imag) << std::endl;
    std::cout << "Pressure skeleton real part L2 error is : "
              << std::sqrt(L2_error_p_hat_real) << std::endl;
    std::cout << "Pressure skeleton imag part L2 error is : "
              << std::sqrt(L2_error_p_hat_imag) << std::endl;

    error_table.add_value("eL2_u_r", std::sqrt(L2_error_u_real));
    error_table.add_value("eL2_u_i", std::sqrt(L2_error_u_imag));
    error_table.add_value("eL2_p_r", std::sqrt(L2_error_p_real));
    error_table.add_value("eL2_p_i", std::sqrt(L2_error_p_imag));
    error_table.add_value("eL2_u_hat_r", std::sqrt(L2_error_u_hat_real));
    error_table.add_value("eL2_u_hat_i", std::sqrt(L2_error_u_hat_imag));
    error_table.add_value("eL2_p_hat_r", std::sqrt(L2_error_p_hat_real));
    error_table.add_value("eL2_p_hat_i", std::sqrt(L2_error_p_hat_imag));
  }

  // @sect3{DPG::refine_grid}
  // This function creates the mesh for the first cycle and then refines it
  // uniformly for subsequent cycles. It also records the number of cells and
  // the maximum cell diameter in the error table for convergence analysis.
  template <int dim>
  void DPGHelmholtz<dim>::refine_grid(const unsigned int cycle)
  {
    if (cycle == 0)
      {
        const Point<dim> p1{0., 0.};
        const Point<dim> p2{1., 1.};

        std::vector<unsigned int> repetitions({2, 2});
        GridGenerator::subdivided_hyper_rectangle(
          triangulation, repetitions, p1, p2, true);
        triangulation.refine_global(0);
      }
    else
      {
        triangulation.refine_global();
      }

    std::cout << "Number of active cells: " << triangulation.n_active_cells()
              << std::endl;

    error_table.add_value("cycle", cycle);
    error_table.add_value("n_cells", triangulation.n_active_cells());
    error_table.add_value("cell_size",
                          GridTools::maximal_cell_diameter<dim>(triangulation));
  }

  // @sect3{DPG::run}
  // This function is the main loop of the program using all the previously
  // defined functions. It is also where the convergence rates are obtained
  // after all the refinement cycles. Note again, after solving the skeleton
  // system, we call the assembly function another time to solve for the
  // interior.
  template <int dim>
  void DPGHelmholtz<dim>::run()
  {
    for (unsigned int cycle = 0; cycle < 8; ++cycle)
      {
        std::cout << "===========================================" << std::endl
                  << "Cycle " << cycle << ':' << std::endl;

        refine_grid(cycle);
        setup_system();
        assemble_system(false);
        solve_linear_system_skeleton();
        assemble_system(true);
        calculate_L2_error();
        output_results(cycle);
      }


    error_table.evaluate_convergence_rates(
      "eL2_u_r", "n_cells", ConvergenceTable::reduction_rate_log2);
    error_table.evaluate_convergence_rates(
      "eL2_u_i", "n_cells", ConvergenceTable::reduction_rate_log2);
    error_table.evaluate_convergence_rates(
      "eL2_p_r", "n_cells", ConvergenceTable::reduction_rate_log2);
    error_table.evaluate_convergence_rates(
      "eL2_p_i", "n_cells", ConvergenceTable::reduction_rate_log2);
    error_table.evaluate_convergence_rates(
      "eL2_u_hat_r", "n_cells", ConvergenceTable::reduction_rate_log2);
    error_table.evaluate_convergence_rates(
      "eL2_u_hat_i", "n_cells", ConvergenceTable::reduction_rate_log2);
    error_table.evaluate_convergence_rates(
      "eL2_p_hat_r", "n_cells", ConvergenceTable::reduction_rate_log2);
    error_table.evaluate_convergence_rates(
      "eL2_p_hat_i", "n_cells", ConvergenceTable::reduction_rate_log2);

    std::cout << "===========================================" << std::endl;
    std::cout << "Convergence table:" << std::endl;
    error_table.write_text(std::cout);
  }
} // end of namespace Step100

// @sect3{The <code>main</code> function}

// This is the main function of the program. It creates an instance of the
// <code>DPGHelmholtz</code> class and calls its run method.
int main()
{
  const unsigned int dim = 2;

  try
    {
      // Here we create the necessary variables for our 2D DPG Helmholtz, i.e.,
      // the degree $p$ of the trial space <code>degree</code>, the degree
      // difference $\Delta p$ between the test and trial spaces
      // <code>delta_degree</code>, the wavenumber $k$ <code>wavenumber</code>,
      // and the angle of incidence of the plane wave $\theta$
      // <code>theta</code> in radians.
      int    degree       = 2;
      int    delta_degree = 1;
      double wavenumber   = 20 * M_PI;
      double theta        = M_PI / 4.;

      std::cout << "===========================================" << std::endl
                << "Trial order: " << degree << std::endl
                << "Test order: " << delta_degree + degree << std::endl
                << "===========================================" << std::endl
                << std::endl;

      Step100::DPGHelmholtz<dim> dpg_helmholtz(degree,
                                               delta_degree,
                                               wavenumber,
                                               theta);

      dpg_helmholtz.run();

      std::cout << std::endl;
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

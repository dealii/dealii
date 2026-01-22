/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2025 by the deal.II authors
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

// The DPG method requires a large breadth of element types which are included
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
#include <deal.II/base/tensor_function.h>

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

#include <fstream>
#include <iostream>

const double pi = dealii::numbers::PI;

namespace Step100
{
  using namespace dealii;

  // In the following function declaration, we create the analytical
  // solutions of the velocity field ($\mathbf{u}$) and the pressure field
  // ($p^*$) as well as the associated boundary values. However, in this
  // tutorial, we will avoid the use of deal.II complex arithmetic capabilities
  // and only use the complex functions that are defined in the C++ standard
  // library. Consequently, in what follows, we will separate the real and
  // imaginary component of our spaces. Therefore, we will also define two
  // implementations of each function, one for the real component and one for
  // the imaginary one.

  // We start and create the analytical solution class for kinematic pressure
  // ($p^*$). The analytical solution depends on the wavenumber $k$ and the
  // angle $\theta$ which are passed to the constructor. Since we are splitting
  // the solution in real and imaginary part, we can directly take $\Re(p^*) =
  // \Re(e^{-i k (x \cos(\theta) + y \sin(\theta))}) = \cos(k (x \cos(\theta) +
  // y \sin(\theta)))$ and $\Im(p^*) = \Im(e^{-i k (x \cos(\theta) + y
  // \sin(\theta))}) = -\sin(k (x \cos(\theta) + y \sin(\theta)))$. The same
  // goes for the velocity field real and imaginary parts. The only difference
  // is that the velocity field is a vector field, so it will be derived from
  // the TensorFunction class and return a Tensor<1,dim> in the value function.
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
  double
  AnalyticalSolution_p_real<dim>::value(const Point<dim> &p,
                                        const unsigned int /*component*/) const
  {
    return std::cos(wavenumber *
                    (p[0] * std::cos(theta) + p[1] * std::sin(theta)));
  }

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
  double
  AnalyticalSolution_p_imag<dim>::value(const Point<dim> &p,
                                        const unsigned int /*component*/) const
  {
    return -std::sin(wavenumber *
                     (p[0] * std::cos(theta) + p[1] * std::sin(theta)));
  }

  template <int dim>
  class AnalyticalSolution_u_real : public TensorFunction<1, dim>
  {
  public:
    AnalyticalSolution_u_real(const double wavenumber, const double theta)
      : TensorFunction<1, dim>()
      , wavenumber(wavenumber)
      , theta(theta)
    {}
    virtual Tensor<1, dim> value(const Point<dim> &p) const override;

  private:
    const double wavenumber;
    const double theta;
  };

  template <int dim>
  Tensor<1, dim>
  AnalyticalSolution_u_real<dim>::value(const Point<dim> &p) const
  {
    AssertDimension(dim, 2);

    Tensor<1, dim> return_value;
    return_value[0] =
      std::cos(theta) *
      std::cos(wavenumber * (p[0] * std::cos(theta) + p[1] * std::sin(theta)));
    return_value[1] =
      std::sin(theta) *
      std::cos(wavenumber * (p[0] * std::cos(theta) + p[1] * std::sin(theta)));
    return return_value;
  }

  template <int dim>
  class AnalyticalSolution_u_imag : public TensorFunction<1, dim>
  {
  public:
    AnalyticalSolution_u_imag(const double wavenumber, const double theta)
      : TensorFunction<1, dim>()
      , wavenumber(wavenumber)
      , theta(theta)
    {}
    virtual Tensor<1, dim> value(const Point<dim> &p) const override;

  private:
    const double wavenumber;
    const double theta;
  };

  template <int dim>
  Tensor<1, dim>
  AnalyticalSolution_u_imag<dim>::value(const Point<dim> &p) const
  {
    AssertDimension(dim, 2);

    Tensor<1, dim> return_value;
    return_value[0] =
      std::cos(theta) *
      -std::sin(wavenumber * (p[0] * std::cos(theta) + p[1] * std::sin(theta)));
    return_value[1] =
      std::sin(theta) *
      -std::sin(wavenumber * (p[0] * std::cos(theta) + p[1] * std::sin(theta)));
    return return_value;
  }

  // Similar classes are required for the boundary values functions. The main
  // difference is that the number of components will now be 4 because those
  // functions will be applied to our space of skeleton unknowns via
  // VectorTools::interpolate_boundary_values. This space has 4 components,
  // because the skeleton unknowns on faces for the velocity field are scalars
  // from the definition $\hat{u}_n = \mathbf{u} \cdot n$ and there are the real
  // and imaginary part of both fields. When the function is not called with the
  // right component index, it will simply return 0.0.

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
    if (component == 2)
      {
        return std::cos(wavenumber * p[1] * std::sin(theta));
      }

    return 0.0;
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
    if (component == 3)
      {
        return -std::sin(wavenumber * p[1] * std::sin(theta));
      }
    return 0.0;
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
    if (component == 0)
      {
        return -1 * (std::sin(theta) *
                     std::cos(wavenumber * p[0] * std::cos(theta)));
      }
    return 0.0;
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
    if (component == 1)
      {
        return std::sin(theta) * std::sin(wavenumber * p[0] * std::cos(theta));
      }
    return 0.0;
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
  // evaluate the test functions we will use in DPG.

  // The same applies for the three FESystem:
  // <code>fe_system_trial_interior</code>,
  // <code>fe_system_trial_skeleton</code> and <code>fe_system_test</code>. In
  // each one of these, we will store the relevant finite element space in the
  // same order to avoid confusion. The first component will therefore always be
  // related to the real part of the velocity, the second component to the its
  // imaginary part, the third component to the real part of the pressure and
  // the fourth component to its imaginary part.

  // The constructor of the class takes four arguments to define the problem.
  // The first two are related to the finite element spaces degree, i.e.,
  // <code>degree</code> defines the polynomial degree of the trial space and
  // <code>delta_degree</code> defines the difference in degree between the
  // trial and the test space. Since the test space needs to be enriched
  // compared to the trial space in DPG, <code>delta_degree</code> must be at
  // least 1 to ensure that the method is well posed. The last two arguments are
  // related to the plane wave parameters. The parameter <code>wavenumber</code>
  // defines the wavenumber $k$ of the plane wave problem while
  // <code>theta</code> defines the incident angle. That angle must be in the
  // closed interval $0$ and $\frac{\pi}{2}$ for the boundary conditions to make
  // sense. All those restrictions are asserted in the constructor.

  // The class also provides a number of member functions that are responsible
  // for setting up, solving, and postprocessing the DPG formulation. The
  // <code>setup_system()</code> function initializes the three DoFHandler
  // objects associated with the interior trial space, the skeleton trial space,
  // and the test space. In addition, it sets up the sparsity pattern, the
  // system matrix, and the right-hand side vector, and establishes the boundary
  // conditions using AffineConstraints, including both Dirichlet and Neumann
  // conditions.

  // The assembly of the linear system is handled by
  // <code>assemble_system()</code>, which is called twice for each problem.
  // When <code>solve_interior = false</code>, the bilinear and linear forms are
  // assembled and the system is locally condensed so that the resulting global
  // system only involves the skeleton unknowns. When
  // <code>solve_interior = true</code>, the system is assembled again and the
  // previously computed skeleton solution is used to reconstruct the interior
  // solution variables. This two-step approach is interesting to reduces the
  // size of the global system that needs to be solve which helps for memory
  // consumption, but requires assembling the system twice. To avoid duplication
  // of the assembly process, we introduce a boolean flag that indicates whether
  // we are assembling for the skeleton solve or for the interior
  // reconstruction.

  // The function <code>solve_linear_system_skeleton()</code> solves the
  // resulting linear system for the skeleton degrees of freedom. Mesh
  // refinement is performed by <code>refine_grid()</code>, which applies
  // uniform refinement to the triangulation. The function
  // <code>output_results()</code> writes both the skeleton and interior
  // solutions to separate VTU files for visualization, while
  // <code>calculate_L2_error()</code> computes the $L^2$ norm of the error
  // using the known analytical solution.

  // In addition to these member functions, the class defines a number of member
  // variables that are used throughout the implementation. These include the
  // triangulation, finite element systems, DoFHandler objects, solution
  // vectors, linear system data structures, and a ConvergenceTable used to
  // store the $L^2$ error and related quantities. The coefficients defining the
  // incident plane wave, namely the wavenumber and the angle of incidence, are
  // also stored as class members.

  // Finally, the class defines several FEValuesExtractors that are reused at
  // multiple points in the implementation to select the appropriate components
  // of the solution fields. These extractors provide access to the real and
  // imaginary parts of the velocity and pressure variables. Since the skeleton
  // space does not have the same number of components as the interior or test
  // spaces (because the $H^{-1/2}$ space associated with the velocity field is
  // scalar) additional extractors are defined specifically for the skeleton
  // variables.
  template <int dim>
  class DPGHelmholtz
  {
  public:
    DPGHelmholtz(const unsigned int degree,
                 const unsigned int delta_degree,
                 const double       wavenumber,
                 const double       theta);
    void run();

  private:
    void setup_system();
    void assemble_system(bool solve_interior);
    void solve_linear_system_skeleton();
    void refine_grid(unsigned int cycle);
    void output_results(unsigned int cycle);
    void calculate_L2_error();

    Triangulation<dim> triangulation;

    const FESystem<dim> fe_trial_interior;
    DoFHandler<dim>     dof_handler_trial_interior;
    Vector<double>      solution_interior;

    const FESystem<dim>       fe_trial_skeleton;
    DoFHandler<dim>           dof_handler_trial_skeleton;
    Vector<double>            solution_skeleton;
    Vector<double>            system_rhs;
    SparsityPattern           sparsity_pattern;
    SparseMatrix<double>      system_matrix;
    AffineConstraints<double> constraints;

    const FESystem<dim> fe_test;
    DoFHandler<dim>     dof_handler_test;

    ConvergenceTable error_table;

    const double wavenumber;
    const double theta;

    const FEValuesExtractors::Vector extractor_u_real;
    const FEValuesExtractors::Vector extractor_u_imag;
    const FEValuesExtractors::Scalar extractor_p_real;
    const FEValuesExtractors::Scalar extractor_p_imag;

    const FEValuesExtractors::Scalar extractor_u_hat_real;
    const FEValuesExtractors::Scalar extractor_u_hat_imag;
    const FEValuesExtractors::Scalar extractor_p_hat_real;
    const FEValuesExtractors::Scalar extractor_p_hat_imag;
  };

  // @sect3{<code>DPGHelmholtz</code> Constructor}
  // In the constructor, we assign the relevant finite element to each FESystem
  // following the nomenclature described above in the class description:
  // - <code>fe_system_trial_interior</code> contains $\Re(\mathbf{u})$,
  // $\Im(\mathbf{u})$, $\Re(p^*)$, $\Im(p^*)$ ;
  // - <code>fe_system_trial_skeleton</code> contains $\Re(\hat{u}_n)$,
  // $\Im(\hat{u}_n)$, $\Re(\hat{p}^*)$, $\Im(\hat{p}^*)$ ;
  // - <code>fe_system_test</code> contains $\Re(\mathbf{v})$,
  // $\Im(\mathbf{v})$, $\Re(q)$, $\Im(q)$.

  // Note that the Q and FE_TraceQ elements have a higher degree than the others
  // because their numbering start at 1 instead of 0. This is to ensure that the
  // spaces chosen follow the exact sequence of energy spaces $Q_{k+1}
  // \rightarrow Nédélec_k \rightarrow Raviart-Thomas_k \rightarrow DGQ_k$.

  // We also initialize the FEValuesExtractors that will be used according
  // to our FESystems nomenclature put the assertions to check if everything is
  // correctly defined for our problem. The first assertion checks that the
  // dimension is 2 because the problem is not implemented in 3D. We also verify
  // that the `delta_degree` variable is at least 1 since the degree of the test
  // space must be at least one degree higher than the trial space. Finally, we
  // check that the wavenumber is positive since it is the magnitude of the wave
  // vector and that the angle theta is in the interval $[0, \pi/2]$ because, as
  // stated above, other angles would not be compatible with the current
  // boundary definitions.

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
    , extractor_u_real(0)
    , extractor_u_imag(dim)
    , extractor_p_real(2 * dim)
    , extractor_p_imag(2 * dim + 1)
    , extractor_u_hat_real(0)
    , extractor_u_hat_imag(1)
    , extractor_p_hat_real(2)
    , extractor_p_hat_imag(3)

  {
    AssertThrow(dim == 2,
                ExcMessage("This tutorial example only works for dim==2"));

    AssertThrow(delta_degree >= 1,
                ExcMessage("The delta_degree needs to be at least 1."));

    AssertThrow(wavenumber > 0, ExcMessage("The wavenumber must be positive."));

    AssertThrow(theta >= 0 && theta <= pi / 2,
                ExcMessage(
                  "The angle theta must be in the interval [0, pi/2]."));
  }

  // @sect3{DPGHelmholtz::setup_system}
  // This function is similar to the other examples. The main difference lies in
  // the fact that we need to set up multiple DOFHandlers for the interior, the
  // skeleton and the test space. The corresponding degrees of freedom are
  // distributed first, and the number of DoFs associated with each space is
  // printed and recorded in the ConvergenceTable for later reference.

  // Since the global linear system is posed exclusively in terms of the
  // skeleton unknowns, constraints are only built for the corresponding
  // DoFHandler. These include hanging-node constraints as well as boundary
  // conditions imposed weakly or strongly depending on their type. In
  // particular, Dirichlet boundary conditions are enforced on selected
  // components of the skeleton variables by interpolating analytical boundary
  // data onto the appropriate trace spaces using component masks and
  // FEValuesExtractors. A Dirichlet condition is first applied to the pressure
  // trace on the left boundary (id=0), while a Neumann condition on the
  // pressure is enforced by prescribing the normal component of the velocity
  // trace on the bottom boundary (id=2). Robin boundary conditions are not
  // enforced through constraints and are instead incorporated later during the
  // assembly of the bilinear and linear forms.

  // Once all constraints have been specified and closed, the vectors and
  // matrices associated with the global linear system are initialized. Because
  // the system only involves skeleton degrees of freedom, the sparsity pattern,
  // system matrix, and right-hand side are constructed accordingly. The
  // solution vectors for both the skeleton and interior unknowns are also
  // initialized at this stage, preparing the class for the subsequent assembly
  // and solution steps.
  template <int dim>
  void DPGHelmholtz<dim>::setup_system()
  {
    dof_handler_trial_skeleton.distribute_dofs(fe_trial_skeleton);
    dof_handler_trial_interior.distribute_dofs(fe_trial_interior);
    dof_handler_test.distribute_dofs(fe_test);

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

    const BoundaryValues_p_real<dim> p_real(wavenumber, theta);
    const BoundaryValues_p_imag<dim> p_imag(wavenumber, theta);
    const BoundaryValues_u_real<dim> u_real(wavenumber, theta);
    const BoundaryValues_u_imag<dim> u_imag(wavenumber, theta);

    VectorTools::interpolate_boundary_values(dof_handler_trial_skeleton,
                                             types::boundary_id(0),
                                             p_real,
                                             constraints,
                                             fe_trial_skeleton.component_mask(
                                               extractor_p_hat_real));
    VectorTools::interpolate_boundary_values(dof_handler_trial_skeleton,
                                             types::boundary_id(0),
                                             p_imag,
                                             constraints,
                                             fe_trial_skeleton.component_mask(
                                               extractor_p_hat_imag));

    VectorTools::interpolate_boundary_values(dof_handler_trial_skeleton,
                                             types::boundary_id(2),
                                             u_real,
                                             constraints,
                                             fe_trial_skeleton.component_mask(
                                               extractor_u_hat_real));
    VectorTools::interpolate_boundary_values(dof_handler_trial_skeleton,
                                             types::boundary_id(2),
                                             u_imag,
                                             constraints,
                                             fe_trial_skeleton.component_mask(
                                               extractor_u_hat_imag));

    constraints.close();

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

  // @sect3{DPGHelmholtz::assemble_system}

  // This function implements the core of the DPG method by assembling the local
  // contributions of the bilinear and linear forms and, depending on the value
  // of <code>solve_interior</code>, either performing local static condensation
  // or reconstructing the interior solution from the skeleton unknowns.

  // We begin by defining volume and face quadrature rules. Since the test space
  // has a higher polynomial degree than the trial spaces by construction, the
  // quadrature order is chosen based on the test finite element to ensure
  // sufficient accuracy for all integrals. The number of quadrature points for
  // both cell and face integration is stored for later use.

  // Next, we create FEValues and FEFaceValues objects for the interior trial,
  // skeleton trial, and test spaces. In the ultraweak formulation used here,
  // gradients are only required for the test functions, while values are needed
  // for all spaces. Because all spaces are defined on the same triangulation,
  // the update of quadrature points and JxW values is only required for one of
  // the FEValues objects, which we choose to be the interior trial space.

  // We then query and store the number of degrees of freedom per cell
  // associated with each finite element space. These values determine the sizes
  // of all local matrices and vectors used during assembly.

  // To avoid repeated queries to FEValues objects at each quadrature point, we
  // allocate containers to store shape function values, gradients, divergences,
  // and their complex conjugates. The first group of containers corresponds to
  // the test space quantities, including vector-valued test functions, their
  // divergence, scalar test functions, and their gradients, both in the cell
  // interior and on faces. The second group stores the interior trial
  // variables, namely the velocity and pressure fields. The third group
  // contains the skeleton trial variables, which represent the normal velocity
  // and pressure traces and their complex conjugates.

  // The local DPG matrices are then allocated. These include the Gram matrix
  // $G$ of the test space, the coupling matrix between test and interior trial
  // spaces $B$, the coupling matrix between test and skeleton trial spaces
  // $\hat{B}$, and a matrix associated with skeleton-only terms $D$. Together,
  // these matrices define the uncondensed local DPG system. In addition, local
  // vectors corresponding to the linear functional in the test space $l$ and to
  // the skeleton trial space $g$ are defined. These vectors are later combined
  // with the local matrices during condensation.

  // To perform local static condensation, we allocate a set of auxiliary
  // matrices that represent intermediate block operators arising in the
  // elimination of interior degrees of freedom ($M_1$, $M_2$, $M_3$, $M_4$ and
  // $M_5$). Further temporary matrices and vectors are also created to store
  // intermediate results during matrix–matrix and matrix–vector products. These
  // temporary objects are labeled with a "temp" prefix for clarity.

  // We also define local cell matrix and right-hand side vector associated with
  // the skeleton degrees of freedom that will be used to solve our system,
  // together with a local-to-global DoF index map used for distribution into
  // the global system. In addition, we define vectors to store the interior
  // solution that will be use when <code>solve_interior</code> is set to
  // <code>true</code>.
  //
  // Finally, since the Helmholtz problem is complex-valued, we also define the
  // imaginary unit and several complex constants that appear in the bilinear
  // and linear forms. Although the global linear system is real-valued, complex
  // arithmetic from the C++ standard library is used locally to simplify the
  // formulation, in the same spirit as in step-81.

  template <int dim>
  void DPGHelmholtz<dim>::assemble_system(const bool solve_interior)
  {
    const QGauss<dim>     quadrature_formula(fe_test.degree + 1);
    const QGauss<dim - 1> face_quadrature_formula(fe_test.degree + 1);
    const unsigned int    n_q_points      = quadrature_formula.size();
    const unsigned int    n_face_q_points = face_quadrature_formula.size();

    FEValues<dim>     fe_values_trial_interior(fe_trial_interior,
                                           quadrature_formula,
                                           update_values |
                                             update_quadrature_points |
                                             update_JxW_values);
    FEValues<dim>     fe_values_test(fe_test,
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

    const unsigned int dofs_per_cell_test = fe_test.n_dofs_per_cell();
    const unsigned int dofs_per_cell_trial_interior =
      fe_trial_interior.n_dofs_per_cell();
    const unsigned int dofs_per_cell_trial_skeleton =
      fe_trial_skeleton.n_dofs_per_cell();

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

    LAPACKFullMatrix<double> G_matrix(dofs_per_cell_test, dofs_per_cell_test);
    LAPACKFullMatrix<double> B_matrix(dofs_per_cell_test,
                                      dofs_per_cell_trial_interior);
    LAPACKFullMatrix<double> B_hat_matrix(dofs_per_cell_test,
                                          dofs_per_cell_trial_skeleton);
    LAPACKFullMatrix<double> D_matrix(dofs_per_cell_trial_skeleton,
                                      dofs_per_cell_trial_skeleton);

    Vector<double> g_vector(dofs_per_cell_trial_skeleton);
    Vector<double> l_vector(dofs_per_cell_test);

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

    LAPACKFullMatrix<double> tmp_matrix(dofs_per_cell_trial_skeleton,
                                        dofs_per_cell_trial_interior);
    LAPACKFullMatrix<double> tmp_matrix2(dofs_per_cell_trial_skeleton,
                                         dofs_per_cell_trial_skeleton);
    LAPACKFullMatrix<double> tmp_matrix3(dofs_per_cell_trial_skeleton,
                                         dofs_per_cell_test);
    Vector<double>           tmp_vector(dofs_per_cell_trial_interior);

    FullMatrix<double> cell_matrix(dofs_per_cell_trial_skeleton,
                                   dofs_per_cell_trial_skeleton);
    Vector<double>     cell_skeleton_rhs(dofs_per_cell_trial_skeleton);
    std::vector<types::global_dof_index> local_dof_indices(
      dofs_per_cell_trial_skeleton);

    Vector<double> cell_interior_rhs(dofs_per_cell_trial_interior);
    Vector<double> cell_interior_solution(dofs_per_cell_trial_interior);
    Vector<double> cell_skeleton_solution(dofs_per_cell_trial_skeleton);

    constexpr std::complex<double> imag(0., 1.);
    const std::complex<double>     iomega      = imag * wavenumber;
    const std::complex<double>     iomega_conj = std::conj(iomega);

    // We now assemble the local contributions of the DPG formulation. As usual,
    // we loop over all active cells of the triangulation. We choose the
    // DoFHandler associated with the interior trial space as the primary
    // iterator, since the cell-wise assembly is naturally tied to the interior
    // unknowns. For each such cell, we explicitly obtain the corresponding
    // iterators for the test space and for the skeleton trial space to ensure
    // that all FEValues objects are reinitialized on the same physical cell.

    // The loop on cells is used to assemble all the local matrices and vectors
    // entering the DPG static condensation procedure ($G$, $B$, $\hat{B}$, $D$,
    // $g$, $l$). It follows that all these objects are reinitialized to zero at
    // the beginning of each cell loop. In addition, we need to reset the local
    // condensation matrix $M_1$ because LAPACKFullMatrix keeps track of its
    // inverse status between iterations and forbids to invert it again if it
    // has the inverted status.

    // To avoid unnecessary computations, at each quadrature point, we evaluate
    // and cache the values, gradients, and divergences of the test functions
    // ($\mathbf{v}$ and $q$), as well as the values of the trial functions
    // ($\mathbf{u}$ and $p$). These quantities are stored as complex-valued
    // expressions, together with their complex conjugates, in order to directly
    // form the sesquilinear forms appearing in the time-harmonic formulation.

    // For each quadrature point, we then loop over the test space degrees of
    // freedom. In a first nested loop over test indices, we assemble the Gram
    // matrix $G$. Depending on whether the test basis functions correspond to
    // the velocity or pressure components, we add the appropriate
    // contributions:
    //  - If both <code>i</code> and <code>j</code> are in test space associated
    //  to the test functions $\mathbf{v}$, we build the terms $(\mathbf{v},
    //  \mathbf{v})_{\Omega_h} + (\nabla \cdot \mathbf{v}, \nabla \cdot
    //  \mathbf{v})_{\Omega_h} + (i\omega\mathbf{v}, i\omega
    //  \mathbf{v})_{\Omega_h}$;
    //  - If the dof <code>i</code> is in test function $\mathbf{v}$ and dof
    //  <code>j</code> in test function $q$ we build the terms $(i\omega
    //  \mathbf{v}, \nabla q)_{\Omega_h} + (\nabla \cdot \mathbf{v}, i\omega
    //  q)_{\Omega_h}$;
    //  - If the dof <code>i</code> is in test function $q$ and the dof
    //  <code>j</code> is in the test function $\mathbf{v}$, we build the terms
    //  $(\nabla q, i\omega \mathbf{v})_{\Omega_h} + (i\omega q, \nabla \cdot
    //  \mathbf{v})_{\Omega_h}$;
    //  - Finally, if both in test functions are in $q$, we build the terms $(q,
    //  q)_{\Omega_h} + (\nabla q,\nabla q)_{\Omega_h} + (i\omega q, i\omega
    //  q)_{\Omega_h}$.

    // In a second nested loop over interior trial space degrees of freedom, we
    // assemble the operator matrix B. Here again, the contributions depend on
    // the pairing of test and trial components:
    //  - If dof <code>i</code> in test function $\mathbf{v}$ and dof
    //  <code>j</code> in trial function $\mathbf{u}$ we build the term
    //  $(\mathbf{v}, i\omega \mathbf{u})_{\Omega_h}$;
    //  - If dof <code>i</code> in test function $\mathbf{v}$ and dof
    //  <code>j</code> in trial function $p$ we build the term $ -( \nabla \cdot
    //  \mathbf{v}, p^*)_{\Omega_h}$;
    //  - If dof <code>i</code> in test function $q$ and dof <code>j</code> in
    //  trial function $\mathbf{u}$ we build the term $-(\nabla q,
    //  \mathbf{u})_{\Omega_h}$;
    //  -  If dof <code>i</code> in test function $q$ and dof <code>j</code> in
    //  trial function $p$ we build the term $(q, i\omega p^*)_{\Omega_h}$.
    //
    // Finally, we assemble the load vector l. In the present plane-wave
    // configuration, the volumetric source term is zero, but we nevertheless
    // include the corresponding code path for completeness by assembling
    //  $(q, l)_{\Omega_h}$ over the cell.
    for (const auto &cell : dof_handler_trial_interior.active_cell_iterators())
      {
        fe_values_trial_interior.reinit(cell);

        const typename DoFHandler<dim>::active_cell_iterator cell_test =
          cell->as_dof_handler_iterator(dof_handler_test);
        fe_values_test.reinit(cell_test);

        const typename DoFHandler<dim>::active_cell_iterator cell_skeleton =
          cell->as_dof_handler_iterator(dof_handler_trial_skeleton);

        G_matrix     = 0;
        B_matrix     = 0;
        B_hat_matrix = 0;
        D_matrix     = 0;
        g_vector     = 0;
        l_vector     = 0;
        M1_matrix    = 0;

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {
            const double JxW = fe_values_trial_interior.JxW(q_point);

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

            for (const auto i : fe_values_test.dof_indices())
              {
                for (const auto j : fe_values_test.dof_indices())
                  {
                    if ((fe_test.shape_function_belongs_to(i,
                                                           extractor_u_real) ||
                         fe_test.shape_function_belongs_to(i,
                                                           extractor_u_imag)) &&
                        (fe_test.shape_function_belongs_to(j,
                                                           extractor_u_real) ||
                         fe_test.shape_function_belongs_to(j,
                                                           extractor_u_imag)))
                      {
                        G_matrix(i, j) +=
                          (((v_conj[i] * v[j]) + (div_v_conj[i] * div_v[j]) +
                            (iomega_conj * v_conj[i] * iomega * v[j])) *
                           JxW)
                            .real();
                      }

                    else if ((fe_test.shape_function_belongs_to(
                                i, extractor_u_real) ||
                              fe_test.shape_function_belongs_to(
                                i, extractor_u_imag)) &&
                             (fe_test.shape_function_belongs_to(
                                j, extractor_p_real) ||
                              fe_test.shape_function_belongs_to(
                                j, extractor_p_imag)))
                      {
                        G_matrix(i, j) +=
                          (((iomega_conj * v_conj[i] * grad_q[j]) +
                            (div_v_conj[i] * iomega * q[j])) *
                           JxW)
                            .real();
                      }

                    else if ((fe_test.shape_function_belongs_to(
                                i, extractor_p_real) ||
                              fe_test.shape_function_belongs_to(
                                i, extractor_p_imag)) &&
                             (fe_test.shape_function_belongs_to(
                                j, extractor_u_real) ||
                              fe_test.shape_function_belongs_to(
                                j, extractor_u_imag)))
                      {
                        G_matrix(i, j) +=
                          (((grad_q_conj[i] * iomega * v[j]) +
                            (iomega_conj * q_conj[i] * div_v[j])) *
                           JxW)
                            .real();
                      }

                    else if ((fe_test.shape_function_belongs_to(
                                i, extractor_p_real) ||
                              fe_test.shape_function_belongs_to(
                                i, extractor_p_imag)) &&
                             (fe_test.shape_function_belongs_to(
                                j, extractor_p_real) ||
                              fe_test.shape_function_belongs_to(
                                j, extractor_p_imag)))
                      {
                        G_matrix(i, j) +=
                          (((q_conj[i] * q[j]) + (grad_q[j] * grad_q_conj[i]) +
                            (iomega_conj * q_conj[i] * iomega * q[j])) *
                           JxW)
                            .real();
                      }
                  }

                for (const auto j : fe_values_trial_interior.dof_indices())
                  {
                    if ((fe_test.shape_function_belongs_to(i,
                                                           extractor_u_real) ||
                         fe_test.shape_function_belongs_to(i,
                                                           extractor_u_imag)) &&
                        (fe_trial_interior.shape_function_belongs_to(
                           j, extractor_u_real) ||
                         fe_trial_interior.shape_function_belongs_to(
                           j, extractor_u_imag)))
                      {
                        B_matrix(i, j) +=
                          ((v_conj[i] * iomega * u[j]) * JxW).real();
                      }

                    else if ((fe_test.shape_function_belongs_to(
                                i, extractor_u_real) ||
                              fe_test.shape_function_belongs_to(
                                i, extractor_u_imag)) &&
                             (fe_trial_interior.shape_function_belongs_to(
                                j, extractor_p_real) ||
                              fe_trial_interior.shape_function_belongs_to(
                                j, extractor_p_imag)))
                      {
                        B_matrix(i, j) -= ((div_v_conj[i] * p[j]) * JxW).real();
                      }

                    else if ((fe_test.shape_function_belongs_to(
                                i, extractor_p_real) ||
                              fe_test.shape_function_belongs_to(
                                i, extractor_p_imag)) &&
                             (fe_trial_interior.shape_function_belongs_to(
                                j, extractor_u_real) ||
                              fe_trial_interior.shape_function_belongs_to(
                                j, extractor_u_imag)))
                      {
                        B_matrix(i, j) -=
                          ((grad_q_conj[i] * u[j]) * JxW).real();
                      }

                    else if ((fe_test.shape_function_belongs_to(
                                i, extractor_p_real) ||
                              fe_test.shape_function_belongs_to(
                                i, extractor_p_imag)) &&
                             (fe_trial_interior.shape_function_belongs_to(
                                j, extractor_p_real) ||
                              fe_trial_interior.shape_function_belongs_to(
                                j, extractor_p_imag)))
                      {
                        B_matrix(i, j) +=
                          ((q_conj[i] * iomega * p[j]) * JxW).real();
                      }
                  }

                if (fe_test.shape_function_belongs_to(i, extractor_p_real) ||
                    fe_test.shape_function_belongs_to(i, extractor_p_imag))
                  {
                    double source_term = 0.0;
                    l_vector(i) += (q_conj[i] * source_term * JxW).real();
                  }
              }
          }

        // We now assemble the skeleton (face) contributions of the DPG
        // formulation. For this purpose, we loop over all faces of the current
        // cell using the DoFHandler associated with the skeleton trial space.
        // On each face, we reinitialize the FEFaceValues objects for both the
        // test space and the skeleton trial space, ensuring that all quantities
        // are evaluated on the same geometric entity.

        // In addition to the standard face integrals, this loop also accounts
        // for Robin boundary conditions, which in the present plane-wave
        // configuration are imposed on two boundaries of the domain (id=1 and
        // id=3). We therefore record both the boundary id. We also record the
        // face number that will be used to consistently determine the
        // orientation of normal fluxes.

        // The Robin terms involve the factor $\frac{\mathbf{k} \cdot
        // \mathbf{n}}{\omega}$. In our configuration, $\omega = k c_s$ with
        // $c_s=1$, and the geometry of the domain implies that $\mathbf{k}
        // \cdot \mathbf{n}$ reduces to either $k\cos{\theta}$ for the right
        // boundary (id=1) or $k\sin{\theta}$ for the top boundary (id=3). The
        // wavenumber cancels out,  and we are left with the cosine and sine of
        // the propagation direction. Note that, in our case study, the term is
        // always real and
        // $\overline{\frac{k_n}{\omega}} = \frac{k_n}{\omega}$ in what follows.

        // Now, as for the cell-wise assembly, we loop over the face quadrature
        // points. We evaluate and cache the relevant quantities to assemble
        // both the Gram matrix face contributions and the face operator matrix
        // \hat{B}. And then, we loop over the test space face degrees of
        // freedom and trial space face degrees of freedom to assemble the
        // corresponding contributions. The face contributions to the Gram
        // matrix only arise when the face lies on a Robin boundary. In that
        // case, we add the appropriate terms based on the shape function the
        // dofs belong to (as for the cell loop):
        //  - If both <code>i</code> and <code>j</code> are in test space
        //  associated to the test functions $\mathbf{v}$ we build, $\langle
        //  \mathbf{v} \cdot \mathbf{n}, \mathbf{v} \cdot \mathbf{n}
        //  \rangle_{\Gamma_1 \cup \Gamma_3}$;
        //  - If the dof <code>i</code> is in test function $\mathbf{v}$ and dof
        //  <code>j</code> in test function $q$, we build $\langle
        //  \mathbf{v} \cdot \mathbf{n},\frac{k_n}{\omega}q
        //  \rangle_{\Gamma_1 \cup \Gamma_3}$
        //  - If the dof <code>i</code> is in test function $q$ and the dof
        //  <code>j</code> is in the test function $\mathbf{v}$, we build
        //  $\langle \frac{k_n}{\omega}q, \mathbf{v}
        // \cdot \mathbf{n} \rangle_{\Gamma_1 \cup
        // \Gamma_3}$;
        //  - Finally, if both in test functions are in $q$, we build t$\langle
        //  \frac{k_n}{\omega}q,
        // \frac{k_n}{\omega}q \rangle_{\Gamma_1
        // \cup \Gamma_3}$.
        //
        // Independently of the boundary type, we also assemble the face
        // operator matrix \hat{B} by looping over the skeleton trial space
        // degrees of freedom. The two terms are :
        //  - If dof <code>i</code> in test function $\mathbf{v}$ and dof
        //  <code>j</code> in trial function $\hat{p}^*$ we build the
        //  term$\left\langle
        // \mathbf{v} \cdot \mathbf{n}, \hat{p}^*
        // \right\rangle_{\partial \Omega_h}$;
        //  - If dof <code>i</code> in test function $q$ and dof <code>j</code>
        //  in trial function $\hat{u}_n$ we build the term $\left\langle q,
        //  \hat{u}_n
        // \right\rangle_{\partial \Omega_h}$.
        //
        // Since the FE_FaceQ elements used to represent the trace of H(div)
        // conforming fields do not encode an intrinsic orientation, special
        // care must be taken to ensure consistency of the numerical flux across
        // shared faces (i.e., the flux that cross a face in a given cell is
        // equal to the flux that crosses the same face in the adjacent cell for
        // which the normal is opposite). To this end, we introduce a sign
        // factor that enforces a unique orientation rule: the flux is always
        // oriented from the cell with the smaller active cell index toward the
        // cell with the larger one. On boundary faces, we recover the standard
        // convention in which the flux is aligned with the outward normal by
        // defining the neighbor cell index as std::numeric_limits<unsigned
        // int>::max().This local rule guarantees that the flux contributions
        // are consistent across neighboring cells without requiring a global
        // orientation of the mesh.

        for (const auto &face : cell_skeleton->face_iterators())
          {
            fe_face_values_test.reinit(cell_test, face);
            fe_values_trial_skeleton.reinit(cell_skeleton, face);

            const auto face_no             = cell->face_iterator_to_index(face);
            const auto current_boundary_id = face->boundary_id();

            const double kn_omega =
              (current_boundary_id == 1) ?
                cos(theta) :
                ((current_boundary_id == 3) ? sin(theta) : 1.);

            for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point)
              {
                const Tensor<1, dim> normal =
                  fe_values_trial_skeleton.normal_vector(q_point);
                const double JxW_face = fe_values_trial_skeleton.JxW(q_point);

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

                for (const auto i : fe_face_values_test.dof_indices())
                  {
                    if (current_boundary_id == 1 || current_boundary_id == 3)
                      {
                        for (const auto j : fe_face_values_test.dof_indices())
                          {
                            if ((fe_test.shape_function_belongs_to(
                                   i, extractor_u_real) ||
                                 fe_test.shape_function_belongs_to(
                                   i, extractor_u_imag)) &&
                                (fe_test.shape_function_belongs_to(
                                   j, extractor_u_real) ||
                                 fe_test.shape_function_belongs_to(
                                   j, extractor_u_imag)))
                              {
                                G_matrix(i, j) +=
                                  (v_face_n_conj[i] * v_face_n[j] * JxW_face)
                                    .real();
                              }

                            else if ((fe_test.shape_function_belongs_to(
                                        i, extractor_u_real) ||
                                      fe_test.shape_function_belongs_to(
                                        i, extractor_u_imag)) &&
                                     ((fe_test.shape_function_belongs_to(
                                         j, extractor_p_real) ||
                                       fe_test.shape_function_belongs_to(
                                         j, extractor_p_imag))))
                              {
                                G_matrix(i, j) += (v_face_n_conj[i] * kn_omega *
                                                   q_face[j] * JxW_face)
                                                    .real();
                              }

                            else if ((fe_test.shape_function_belongs_to(
                                        i, extractor_p_real) ||
                                      fe_test.shape_function_belongs_to(
                                        i, extractor_p_imag)) &&
                                     (fe_test.shape_function_belongs_to(
                                        j, extractor_u_real) ||
                                      fe_test.shape_function_belongs_to(
                                        j, extractor_u_imag)))
                              {
                                G_matrix(i, j) += (kn_omega * q_face_conj[i] *
                                                   v_face_n[j] * JxW_face)
                                                    .real();
                              }

                            else if ((fe_test.shape_function_belongs_to(
                                        i, extractor_p_real) ||
                                      fe_test.shape_function_belongs_to(
                                        i, extractor_p_imag)) &&
                                     ((fe_test.shape_function_belongs_to(
                                         j, extractor_p_real) ||
                                       fe_test.shape_function_belongs_to(
                                         j, extractor_p_imag))))
                              {
                                G_matrix(i, j) +=
                                  (kn_omega * q_face_conj[i] * kn_omega *
                                   q_face[j] * JxW_face)
                                    .real();
                              }
                          }
                      }

                    for (const auto j : fe_values_trial_skeleton.dof_indices())
                      {
                        if ((fe_test.shape_function_belongs_to(
                               i, extractor_u_real) ||
                             fe_test.shape_function_belongs_to(
                               i, extractor_u_imag)) &&
                            (fe_trial_skeleton.shape_function_belongs_to(
                               j, extractor_p_hat_real) ||
                             fe_trial_skeleton.shape_function_belongs_to(
                               j, extractor_p_hat_imag)))
                          {
                            B_hat_matrix(i, j) +=
                              ((v_face_n_conj[i] * p_hat[j]) * JxW_face).real();
                          }

                        else if ((fe_test.shape_function_belongs_to(
                                    i, extractor_p_real) ||
                                  fe_test.shape_function_belongs_to(
                                    i, extractor_p_imag)) &&
                                 (fe_trial_skeleton.shape_function_belongs_to(
                                    j, extractor_u_hat_real) ||
                                  fe_trial_skeleton.shape_function_belongs_to(
                                    j, extractor_u_hat_imag)))
                          {
                            const unsigned int neighbor_cell_id =
                              face->at_boundary() ?
                                std::numeric_limits<unsigned int>::max() :
                                cell->neighbor(face_no)->active_cell_index();
                            const double flux_orientation =
                              neighbor_cell_id > cell->active_cell_index() ?
                                1. :
                                -1.;

                            B_hat_matrix(i, j) +=
                              (q_face_conj[i] * flux_orientation * u_hat_n[j] *
                               JxW_face)
                                .real();
                          }
                      }
                  }

                // Finally, we assemble the matrix $D$ and the corresponding
                // source-term vector $g$ associated with the Robin boundary
                // conditions. We first verify whether the current face belongs
                // to one of the Robin boundaries, which in this configuration
                // correspond to boundary ids 1 and 3.

                // On these boundaries, the numerical flux orientation is
                // unambiguous: since the face lies on the exterior boundary of
                // the domain, the flux is always aligned with the outward
                // normal. Consequently, the flux orientation factor is set to
                // +1.

                // We then loop over the skeleton trial-space degrees of freedom
                // and assemble both the boundary load vector $g$ and the $D$
                // matrix.

                // The matrix $D$ couples the skeleton trial space with itself,
                // and its entries are therefore obtained by looping again over
                // the skeleton trial-space degrees of freedom. Depending on
                // whether the basis functions correspond to velocity or
                // pressure traces, we assemble the following Robin boundary
                // contributions:
                //  //  - If both <code>i</code> and <code>j</code> are in the
                //  skeleton trace associated to $\hat{u}_n$, we build the term
                //  $- \langle \hat{u}_n, \hat{u}_n \rangle_{\Gamma_1 \cup
                //  \Gamma_3}$;
                //  - If <code>i</code> is in the skeleton trace associated to
                //  $\hat{u}_n$ and <code>j</code> in the skeleton trace
                //  associated to
                //  $\hat{p}$, we build the term $\langle \hat{u}_n,
                //  \frac{k_n}{\omega} \hat{p}^* \rangle_{\Gamma_1 \cup
                //  \Gamma_ 3}$;
                //  - If <code>i</code> is in the skeleton trace associated to
                //  $\hat{p}$ and <code>j</code> in the skeleton trace
                //  associated to
                //  $\hat{u}_n$, we build the term $\langle \frac{k_n}{\omega}
                //  \hat{p}^*, \hat{u}_n \rangle_{\Gamma_1 \cup \Gamma_3}$;
                //  - If both <code>i</code> and <code>j</code> are in the
                //  skeleton trace associated to $\hat{p}$, we build the term $-
                //  \langle \frac{k_n}{\omega} \hat{p}^*, \frac{k_n}{\omega}
                //  \hat{p}^* \rangle_{\Gamma_1 \cup \Gamma_3}$.

                // Although the source term associated with the Robin condition
                // $g$ is zero for the present plane-wave problem, we still
                // assemble it explicitly in order to illustrate the general
                // procedure:
                //
                //  -  If both <code>i</code> are in the skeleton trace
                //  associated to
                //  $\hat{u}_n$, we assemble the term $- \langle \hat{u}_n, g_R
                //  \rangle_{\Gamma_1 \cup \Gamma_3}$;
                //  -  If both <code>i</code> are in the skeleton trace
                //  associated to
                //  $\hat{p}$, we assemble the term $\langle (k_n / \omega)
                //  \hat{p}^*, g_R \rangle_{\Gamma_1 \cup \Gamma_3}$.
                //

                if (current_boundary_id == 1 || current_boundary_id == 3)
                  {
                    const double flux_orientation = 1.;

                    for (const auto i : fe_values_trial_skeleton.dof_indices())
                      {
                        double source_term = 0.;
                        if (fe_trial_skeleton.shape_function_belongs_to(
                              i, extractor_u_hat_real) ||
                            fe_trial_skeleton.shape_function_belongs_to(
                              i, extractor_u_hat_imag))
                          {
                            g_vector(i) -=
                              (u_hat_n_conj[i] * source_term).real() * JxW_face;
                          }
                        else if (fe_trial_skeleton.shape_function_belongs_to(
                                   i, extractor_p_hat_real) ||
                                 fe_trial_skeleton.shape_function_belongs_to(
                                   i, extractor_p_hat_imag))
                          {
                            g_vector(i) +=
                              (kn_omega * p_hat_conj[i] * source_term).real() *
                              JxW_face;
                          }

                        for (const auto j :
                             fe_values_trial_skeleton.dof_indices())
                          {
                            if ((fe_trial_skeleton.shape_function_belongs_to(
                                   i, extractor_u_hat_real) ||
                                 fe_trial_skeleton.shape_function_belongs_to(
                                   i, extractor_u_hat_imag)) &&
                                (fe_trial_skeleton.shape_function_belongs_to(
                                   j, extractor_u_hat_real) ||
                                 fe_trial_skeleton.shape_function_belongs_to(
                                   j, extractor_u_hat_imag)))
                              {
                                D_matrix(i, j) -=
                                  (flux_orientation * u_hat_n_conj[i] *
                                   flux_orientation * u_hat_n[j] * JxW_face)
                                    .real();
                              }

                            else if ((fe_trial_skeleton
                                        .shape_function_belongs_to(
                                          i, extractor_u_hat_real) ||
                                      fe_trial_skeleton
                                        .shape_function_belongs_to(
                                          i, extractor_u_hat_imag)) &&
                                     (fe_trial_skeleton
                                        .shape_function_belongs_to(
                                          j, extractor_p_hat_real) ||
                                      fe_trial_skeleton
                                        .shape_function_belongs_to(
                                          j, extractor_p_hat_imag)))
                              {
                                D_matrix(i, j) +=
                                  (flux_orientation * u_hat_n_conj[i] *
                                   kn_omega * p_hat[j] * JxW_face)
                                    .real();
                              }

                            else if ((fe_trial_skeleton
                                        .shape_function_belongs_to(
                                          i, extractor_p_hat_real) ||
                                      fe_trial_skeleton
                                        .shape_function_belongs_to(
                                          i, extractor_p_hat_imag)) &&
                                     (fe_trial_skeleton
                                        .shape_function_belongs_to(
                                          j, extractor_u_hat_real) ||
                                      fe_trial_skeleton
                                        .shape_function_belongs_to(
                                          j, extractor_u_hat_imag)))
                              {
                                D_matrix(i, j) +=
                                  (kn_omega * p_hat_conj[i] * flux_orientation *
                                   u_hat_n[j] * JxW_face)
                                    .real();
                              }

                            else if ((fe_trial_skeleton
                                        .shape_function_belongs_to(
                                          i, extractor_p_hat_real) ||
                                      fe_trial_skeleton
                                        .shape_function_belongs_to(
                                          i, extractor_p_hat_imag)) &&
                                     (fe_trial_skeleton
                                        .shape_function_belongs_to(
                                          j, extractor_p_hat_real) ||
                                      fe_trial_skeleton
                                        .shape_function_belongs_to(
                                          j, extractor_p_hat_imag)))
                              {
                                D_matrix(i, j) -=
                                  (kn_omega * p_hat_conj[i] * kn_omega *
                                   p_hat[j] * JxW_face)
                                    .real();
                              }
                          }
                      }
                  }
              }
          }
        // After assembling all local matrices and vectors, we perform the
        // cell-wise static condensation associated with the DPG formulation. We
        // first invert the test-space Gram matrix G and use it to form the
        // auxiliary operators $M_4 = B^\dagger G^{-1}$ and $M_5 =
        // \hat{B}^\dagger G^{-1}$ . These are then used to construct the
        // condensed blocks $M_1 = B^\dagger G^{-1} B$, $M_2 = B^\dagger G^{-1}
        // \hat{B}$ and $M_3 = \hat{B}^\dagger G^{-1} \hat{B} - D$. Then, if
        // <code>solve_interior</code> is <code>true</code>, the skeleton
        // solution $\hat{u}_h$ is assumed known and we recover the interior
        // unknowns on each cell by solving $u_h = M1^{-1} (M4 l - M2
        // \hat{u}_h)$, followed by distribution to the global interior solution
        // vector. Otherwise, we assemble the fully condensed local system for
        // the skeleton unknowns by forming the Schur complement $(M3 - M2^†
        // M1^{-1} M2)$, together with the corresponding right-hand side $(M5 -
        // M2^† M1^{-1} M4) l - g$, and distribute the resulting local matrix
        // and vector to the global skeleton system while enforcing constraints.
        G_matrix.invert();

        B_matrix.Tmmult(M4_matrix, G_matrix);
        B_hat_matrix.Tmmult(M5_matrix, G_matrix);

        M4_matrix.mmult(M1_matrix, B_matrix);
        M4_matrix.mmult(M2_matrix, B_hat_matrix);

        M5_matrix.mmult(M3_matrix, B_hat_matrix);
        M3_matrix.add(-1.0, D_matrix);

        M1_matrix.invert();

        if (solve_interior)
          {
            cell_skeleton->get_dof_values(solution_skeleton,
                                          cell_skeleton_solution);

            M2_matrix.vmult(tmp_vector, cell_skeleton_solution);
            M4_matrix.vmult(cell_interior_rhs, l_vector);
            cell_interior_rhs -= tmp_vector;
            M1_matrix.vmult(cell_interior_solution, cell_interior_rhs);

            cell->distribute_local_to_global(cell_interior_solution,
                                             solution_interior);
          }

        else
          {
            M2_matrix.Tmmult(tmp_matrix, M1_matrix);
            tmp_matrix.mmult(tmp_matrix2, M2_matrix);
            tmp_matrix2.add(-1.0, M3_matrix);
            tmp_matrix2 *= -1.0;
            cell_matrix = tmp_matrix2;

            tmp_matrix.mmult(tmp_matrix3, M4_matrix);
            M5_matrix.add(-1.0, tmp_matrix3);
            M5_matrix.vmult(cell_skeleton_rhs, l_vector);
            cell_skeleton_rhs -= g_vector;

            cell_skeleton->get_dof_indices(local_dof_indices);
            constraints.distribute_local_to_global(cell_matrix,
                                                   cell_skeleton_rhs,
                                                   local_dof_indices,
                                                   system_matrix,
                                                   system_rhs);
          }
      }
  }

  // @sect3{DPGHelmholtz::solve_linear_system_skeleton}
  // This function is in charge of solving the linear system assembled and has
  // nothing specific to DPG per se. Even though the original PDE was
  // indefinite, the way we solve it (using a least-squares approach) means that
  // the linear system is symmetric and positive definite. As a consequence, the
  // method allows us to use the Conjugate Gradient iterative solver. Note that
  // because we do not have any preconditioner, the number of iterations can be
  // quite high. For simplicity, we put a high upper limit on the number of
  // iterations, but in practice one would want to change this function to have
  // a more robust solver. The tolerance for the convergence here is defined
  // proportional to the $l^2$ norm of the RHS vector so the stopping criterion
  // independent of whatever scaling we apply to the equation. Note that the
  // tolerance is quite small, but this is needed to recover the accurate order
  // of convergence for the smallest mesh sizes.
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

  // @sect3{DPGHelmholtz::output_results}
  // This function outputs both the interior and skeleton solutions in VTU
  // format for visualization in ParaView or VisIt. The interior solution is
  // written using the standard DataOut class by attaching the interior
  // DoFHandler and providing appropriate component names and interpretations:
  // the real and imaginary parts of the velocity are treated as vector-valued
  // fields, while the real and imaginary parts of the pressure are treated as
  // scalar fields. The skeleton solution, which is defined only on mesh faces,
  // is handled separately using the DataOutFaces class like presented in
  // step-51 for the HDG method. Note also that in this case all components are
  // scalar-valued. For both outputs, visualization patches are built using the
  // corresponding polynomial degree, and the results are written to
  // cycle-dependent VTU files.
  template <int dim>
  void DPGHelmholtz<dim>::output_results(const unsigned int cycle)
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler_trial_interior);

    std::vector<std::string> solution_interior_names;
    for (unsigned int i = 0; i < dim; ++i)
      {
        solution_interior_names.emplace_back("velocity_real");
      }
    for (unsigned int i = 0; i < dim; ++i)
      {
        solution_interior_names.emplace_back("velocity_imag");
      }
    solution_interior_names.emplace_back("pressure_real");
    solution_interior_names.emplace_back("pressure_imag");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation;
    for (unsigned int i = 0; i < dim; ++i)
      {
        data_component_interpretation.push_back(
          DataComponentInterpretation::component_is_part_of_vector);
      }
    for (unsigned int i = 0; i < dim; ++i)
      {
        data_component_interpretation.push_back(
          DataComponentInterpretation::component_is_part_of_vector);
      }
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);

    data_out.add_data_vector(solution_interior,
                             solution_interior_names,
                             DataOut<dim>::type_automatic,
                             data_component_interpretation);
    data_out.build_patches(fe_trial_interior.degree);
    std::ofstream output("solution_planewave_square-" + std::to_string(cycle) +
                         ".vtu");
    data_out.write_vtu(output);

    DataOutFaces<dim> data_out_faces(false);
    data_out_faces.attach_dof_handler(dof_handler_trial_skeleton);

    std::vector<std::string> solution_skeleton_names;
    solution_skeleton_names.emplace_back("velocity_hat_real");
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
                              std::to_string(cycle) + ".vtu");
    data_out_faces.write_vtu(output_face);
  }

  // @sect3{DPGHelmholtz::calculate_L2_error}
  // In this function, we compute the $L_2$ error of each component of the
  // numerical solution, namely the real and imaginary parts of the velocity and
  // pressure, for both the interior and skeleton unknowns. Because we want to
  // have errors for the skeleton components, we cannot use the
  // VectorTools::integrate_difference function as it does not have a mechanism
  // to avoid visiting faces twice (i.e., counting the error on each cell
  // sharing the face). Therefore, we will perform the computation "by hand" for
  // both interior and skeleton solutionsby looping over all cells and faces,
  // interpolating the numerical solution at the corresponding quadrature
  // points, and comparing it to the analytical solution. To this end, FEValues
  // and FEFaceValues objects are instantiated in the same manner as during
  // assembly, appropriate volume and face quadrature rules are used, and
  // containers are allocated to store the interpolated interior and skeleton
  // fields at quadrature points. Separate accumulators are maintained for the
  // $L_2$ errors of each solution component. Finally, analytical solution
  // objects are initialized to provide the exact reference values used in the
  // error evaluation.

  template <int dim>
  void DPGHelmholtz<dim>::calculate_L2_error()
  {
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

    double L2_error_p_real     = 0;
    double L2_error_p_imag     = 0;
    double L2_error_p_hat_real = 0;
    double L2_error_p_hat_imag = 0;
    double L2_error_u_real     = 0;
    double L2_error_u_imag     = 0;
    double L2_error_u_hat_real = 0;
    double L2_error_u_hat_imag = 0;

    std::vector<Tensor<1, dim>> local_u_real(n_q_points);
    std::vector<Tensor<1, dim>> local_u_imag(n_q_points);
    std::vector<double>         local_p_real(n_q_points);
    std::vector<double>         local_p_imag(n_q_points);
    std::vector<double>         local_u_hat_real(n_face_q_points);
    std::vector<double>         local_u_hat_imag(n_face_q_points);
    std::vector<double>         local_p_hat_real(n_face_q_points);
    std::vector<double>         local_p_hat_imag(n_face_q_points);

    const AnalyticalSolution_p_real<dim> analytical_solution_p_real(wavenumber,
                                                                    theta);
    const AnalyticalSolution_p_imag<dim> analytical_solution_p_imag(wavenumber,
                                                                    theta);
    const AnalyticalSolution_u_real<dim> analytical_solution_u_real(wavenumber,
                                                                    theta);
    const AnalyticalSolution_u_imag<dim> analytical_solution_u_imag(wavenumber,
                                                                    theta);

    // We compute the $L_2$ error by looping over all active cells of the mesh
    // and evaluating both the interior and skeleton contributions. For each
    // cell, the interior velocity and pressure (real and imaginary parts) are
    // first interpolated at volume quadrature points using FEValues, and their
    // squared differences with the corresponding analytical solutions are
    // accumulated using the Jacobian-weighted quadrature weights. The skeleton
    // error is then computed by looping over the faces of the same cell and
    // interpolating the trace unknowns at face quadrature points using
    // FEFaceValues. To avoid double-counting interior faces shared by two
    // cells, we use a similar idea than the one used to define the flux
    // orientation during the assembly : each face is integrated only once by
    // retaining the contribution from the cell with the smallest active cell
    // index (boundary faces are always included). An additional important
    // detail concerns the error computation for the velocity trace variable.
    // Indeed, for the normal flux trace variable
    // $\hat{u}_n$, the analytical velocity is projected onto the outward normal
    // at each quadrature point, and the error is computed using the only the
    // magnitude (absolute value) of both numerical and analytical quantities.
    // This choice removes spurious sign changes induced by face-normal
    // orientation conventions, which may differ between neighboring cells and
    // are not physically meaningful for error estimation. Finally, the
    // accumulated errors are printed to the terminal, and stored in the
    // <code>error_table</code> for post-processing.
    for (const auto &cell : dof_handler_trial_interior.active_cell_iterators())
      {
        fe_values_trial_interior.reinit(cell);

        fe_values_trial_interior[extractor_u_real].get_function_values(
          solution_interior, local_u_real);
        fe_values_trial_interior[extractor_u_imag].get_function_values(
          solution_interior, local_u_imag);
        fe_values_trial_interior[extractor_p_real].get_function_values(
          solution_interior, local_p_real);
        fe_values_trial_interior[extractor_p_imag].get_function_values(
          solution_interior, local_p_imag);

        const auto &quadrature_points =
          fe_values_trial_interior.get_quadrature_points();

        for (const unsigned int q_index :
             fe_values_trial_interior.quadrature_point_indices())
          {
            const double JxW      = fe_values_trial_interior.JxW(q_index);
            const auto  &position = quadrature_points[q_index];

            L2_error_u_real += (local_u_real[q_index] -
                                analytical_solution_u_real.value(position))
                                 .norm_square() *
                               JxW;
            L2_error_u_imag += (local_u_imag[q_index] -
                                analytical_solution_u_imag.value(position))
                                 .norm_square() *
                               JxW;

            L2_error_p_real +=
              std::pow((local_p_real[q_index] -
                        analytical_solution_p_real.value(position, 0)),
                       2) *
              JxW;

            L2_error_p_imag +=
              std::pow((local_p_imag[q_index] -
                        analytical_solution_p_imag.value(position, 0)),
                       2) *
              JxW;
          }

        const typename DoFHandler<dim>::active_cell_iterator cell_skeleton =
          cell->as_dof_handler_iterator(dof_handler_trial_skeleton);

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

            const auto &face_quadrature_points =
              fe_values_trial_skeleton.get_quadrature_points();

            for (const unsigned int &q_index :
                 fe_values_trial_skeleton.quadrature_point_indices())
              {
                const double JxW      = fe_values_trial_skeleton.JxW(q_index);
                const auto  &position = face_quadrature_points[q_index];
                const Tensor<1, dim> normal =
                  fe_values_trial_skeleton.normal_vector(q_index);

                const unsigned int neighbor_cell_id =
                  face->at_boundary() ?
                    std::numeric_limits<unsigned int>::max() :
                    cell->neighbor(face_no)->active_cell_index();
                if (neighbor_cell_id < cell->active_cell_index())
                  {
                    continue;
                  }

                double u_hat_n_analytical_real =
                  normal * analytical_solution_u_real.value(position);
                double u_hat_n_analytical_imag =
                  normal * analytical_solution_u_imag.value(position);

                L2_error_u_hat_real += std::pow(abs(local_u_hat_real[q_index]) -
                                                  abs(u_hat_n_analytical_real),
                                                2) *
                                       JxW;
                L2_error_u_hat_imag += std::pow(abs(local_u_hat_imag[q_index]) -
                                                  abs(u_hat_n_analytical_imag),
                                                2) *
                                       JxW;

                L2_error_p_hat_real +=
                  std::pow((local_p_hat_real[q_index] -
                            analytical_solution_p_real.value(position, 0)),
                           2) *
                  JxW;
                L2_error_p_hat_imag +=
                  std::pow((local_p_hat_imag[q_index] -
                            analytical_solution_p_imag.value(position, 0)),
                           2) *
                  JxW;
              }
          }
      }

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

  // @sect3{DPGHelmholtz::refine_grid}
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

  // @sect3{DPGHelmholtz::run}
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
} // namespace Step100

// @sect3{The <code>main</code> function}

// This is the main function of the program. It creates an instance of the
// <code>DPGHelmholtz</code> class and calls its run method. It defines the
// necessary variables for our 2D DPG Helmholtz, i.e., the degree $p$ of the
// trial space <code>degree</code>, the degree difference $\Delta p$ between the
// test and trial spaces <code>delta_degree</code>, the wavenumber $k$
// <code>wavenumber</code>, and the angle of incidence of the plane wave
// $\theta$ <code>theta</code> in radians.
int main()
{
  const unsigned int dim = 2;

  try
    {
      const int    degree       = 2;
      const int    delta_degree = 1;
      const double wavenumber   = 20 * pi;
      const double theta        = pi / 4.;

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

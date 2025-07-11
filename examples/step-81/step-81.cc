/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2022 - 2025 by the deal.II authors
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
 * Authors: Manaswinee Bezbaruah, Matthias Maier, Texas A&M University, 2021.
 */

// @sect3{Include files}

// The set of include files is quite standard. The most notable include is
// the fe/fe_nedelec_sz.h file which allows us to use the FE_NedelecSZ elements.
// This is an implementation of the $H^{curl}$ conforming Nédélec Elements
// that resolves the sign conflict issues that arise from parametrization
// (for details we refer to the documentation of the FE_NedelecSZ element).

#include <deal.II/base/function.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nedelec_sz.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>


#include <fstream>
#include <iostream>
#include <memory>


// @sect3{Class Template Declarations}
//
// We begin our actual implementation by declaring all classes with their
// data structures and methods upfront.

namespace Step81
{
  using namespace dealii;
  using namespace std::complex_literals;

  // @sect4{Parameters Class}

  // The Parameters class inherits ParameterAcceptor, and instantiates all the
  // coefficients in our variational equations.
  // These coefficients are passed through ParameterAcceptor and are editable
  // through a .prm file.
  // More explanation on the use and inheritance from the ParameterAcceptor
  // can be found in step-60.

  // epsilon is the Electric Permittivity coefficient and it is a rank 2 tensor.
  // Depending on the material, we assign the i^th diagonal element of the
  // tensor to the material epsilon value (one of the private epsilon_1_ or
  // epsilon_2_ variables).
  //
  // mu_inv  is the inverse of the Magnetic Permiability coefficient and it is
  // a complex number.

  // sigma is the Surface Conductivity coefficient between material left and
  // material right and it is a rank 2 tensor. It is only changed if we are at
  // the interface between two materials. If we are at an interface, we assign
  // the i^th diagonal element of the tensor to the private sigma_ value.

  // J_a is the strength and orientation of the dipole. As mentioned in the
  // rescaling,
  // @f[
  // \mathbf{J}_a = J_0 \mathbf{e}_i\delta(x-a)
  // @f]
  // It is a rank 1 tensor
  // that depends on the private dipole_position, dipole_radius,
  // dipole_strength, dipole_orientation variables.

  template <int dim>
  class Parameters : public ParameterAcceptor
  {
  public:
    Parameters();

    using rank0_type = std::complex<double>;

    using rank1_type = Tensor<1, dim, std::complex<double>>;

    using rank2_type = Tensor<2, dim, rank0_type>;

    using curl_type = Tensor<1, dim == 2 ? 1 : dim, rank0_type>;

  public:
    rank2_type epsilon(const Point<dim> &x, types::material_id material);

    std::complex<double> mu_inv(const Point<dim>  &x,
                                types::material_id material);

    rank2_type sigma(const Point<dim>  &x,
                     types::material_id left,
                     types::material_id right);

    rank1_type J_a(const Point<dim> &point, types::material_id id);

  private:
    rank2_type           epsilon_1;
    rank2_type           epsilon_2;
    std::complex<double> mu_inv_1;
    std::complex<double> mu_inv_2;
    rank2_type           sigma_tensor;

    double         dipole_radius;
    Point<dim>     dipole_position;
    Tensor<1, dim> dipole_orientation;
    rank0_type     dipole_strength;
  };


  template <int dim>
  Parameters<dim>::Parameters()
    : ParameterAcceptor("Parameters")
  {
    epsilon_1[0][0] = 1.;
    epsilon_1[1][1] = 1.;
    add_parameter("material 1 epsilon",
                  epsilon_1,
                  "relative permittivity of material 1");

    epsilon_2[0][0] = 1.;
    epsilon_2[1][1] = 1.;
    add_parameter("material 2 epsilon",
                  epsilon_2,
                  "relative permittivity of material 2");

    mu_inv_1 = 1.;
    add_parameter("material 1 mu_inv",
                  mu_inv_1,
                  "inverse of relative permeability of material 1");

    mu_inv_2 = 1.;
    add_parameter("material 2 mu_inv",
                  mu_inv_2,
                  "inverse of relative permeability of material 2");

    sigma_tensor[0][0] = std::complex<double>(0.001, 0.2);
    sigma_tensor[1][1] = std::complex<double>(0.001, 0.2);
    add_parameter("sigma",
                  sigma_tensor,
                  "surface conductivity between material 1 and material 2");

    dipole_radius = 0.3;
    add_parameter("dipole radius", dipole_radius, "radius of the dipole");

    dipole_position = Point<dim>(0., 0.8);
    add_parameter("dipole position", dipole_position, "position of the dipole");

    dipole_orientation = Tensor<1, dim>{{0., 1.}};
    add_parameter("dipole orientation",
                  dipole_orientation,
                  "orientation of the dipole");

    dipole_strength = 1.;
    add_parameter("dipole strength", dipole_strength, "strength of the dipole");
  }

  template <int dim>
  typename Parameters<dim>::rank2_type
  Parameters<dim>::epsilon(const Point<dim> & /*x*/,
                           types::material_id material)
  {
    return (material == 1 ? epsilon_1 : epsilon_2);
  }

  template <int dim>
  std::complex<double> Parameters<dim>::mu_inv(const Point<dim> & /*x*/,
                                               types::material_id material)
  {
    return (material == 1 ? mu_inv_1 : mu_inv_2);
  }

  template <int dim>
  typename Parameters<dim>::rank2_type
  Parameters<dim>::sigma(const Point<dim> & /*x*/,
                         types::material_id left,
                         types::material_id right)
  {
    return (left == right ? rank2_type() : sigma_tensor);
  }

  template <int dim>
  typename Parameters<dim>::rank1_type
  Parameters<dim>::J_a(const Point<dim> &point, types::material_id /*id*/)
  {
    rank1_type J_a;
    const auto distance = (dipole_position - point).norm() / dipole_radius;
    if (distance > 1.)
      return J_a;
    double scale = std::cos(distance * numbers::PI / 2.) *
                   std::cos(distance * numbers::PI / 2.) /
                   (numbers::PI / 2. - 2. / numbers::PI) / dipole_radius /
                   dipole_radius;
    J_a = dipole_strength * dipole_orientation * scale;
    return J_a;
  }

  // @sect4{PerfectlyMatchedLayer Class}
  // The PerfectlyMatchedLayer class inherits ParameterAcceptor as well. It
  // implements the transformation matrices used to modify the permittivity
  // and permeability tensors supplied from the Parameters class. The
  // actual transformation of the material tensors will be done in the
  // assembly loop. The radii and the strength of the PML is specified, and
  // the coefficients will be modified using transformation matrices within
  // the PML region. The radii and strength of the PML are editable through
  // a .prm file. The rotation function $T_{exer}$ is the same as
  // introduced in the perfectly matched layer section of the introduction.
  // Similarly, the matrices A, B and C are defined as follows
  // @f[
  // A = T_{e_xe_r}^{-1}
  // \text{diag}\left(\frac{1}{\bar{d}^2},\frac{1}{d\bar{d}}\right)T_{e_xe_r},\qquad
  // B = T_{e_xe_r}^{-1} \text{diag}\left(d,\bar{d}\right)T_{e_xe_r},\qquad
  // C = T_{e_xe_r}^{-1} \text{diag}\left(\frac{1}{\bar{d}},\frac{1}{d}\right)
  // T_{e_xe_r}.\qquad
  // @f]

  template <int dim>
  class PerfectlyMatchedLayer : public ParameterAcceptor
  {
  public:
    static_assert(dim == 2,
                  "The perfectly matched layer is only implemented in 2d.");

    Parameters<dim> parameters;

    using rank1_type = Tensor<1, dim, std::complex<double>>;

    using rank2_type = Tensor<2, dim, std::complex<double>>;

    PerfectlyMatchedLayer();

    std::complex<double> d(const Point<dim> point);

    std::complex<double> d_bar(const Point<dim> point);


    rank2_type rotation(std::complex<double> d_1,
                        std::complex<double> d_2,
                        Point<dim>           point);

    rank2_type a_matrix(const Point<dim> point);

    rank2_type b_matrix(const Point<dim> point);

    rank2_type c_matrix(const Point<dim> point);

  private:
    double inner_radius;
    double outer_radius;
    double strength;
  };


  template <int dim>
  PerfectlyMatchedLayer<dim>::PerfectlyMatchedLayer()
    : ParameterAcceptor("PerfectlyMatchedLayer")
  {
    inner_radius = 12.;
    add_parameter("inner radius",
                  inner_radius,
                  "inner radius of the PML shell");
    outer_radius = 20.;
    add_parameter("outer radius",
                  outer_radius,
                  "outer radius of the PML shell");
    strength = 8.;
    add_parameter("strength", strength, "strength of the PML");
  }


  template <int dim>
  typename std::complex<double>
  PerfectlyMatchedLayer<dim>::d(const Point<dim> point)
  {
    const auto radius = point.norm();
    if (radius > inner_radius)
      {
        const double s =
          strength * ((radius - inner_radius) * (radius - inner_radius)) /
          ((outer_radius - inner_radius) * (outer_radius - inner_radius));
        return {1.0, s};
      }
    else
      {
        return 1.0;
      }
  }


  template <int dim>
  typename std::complex<double>
  PerfectlyMatchedLayer<dim>::d_bar(const Point<dim> point)
  {
    const auto radius = point.norm();
    if (radius > inner_radius)
      {
        const double s_bar =
          strength / 3. *
          ((radius - inner_radius) * (radius - inner_radius) *
           (radius - inner_radius)) /
          (radius * (outer_radius - inner_radius) *
           (outer_radius - inner_radius));
        return {1.0, s_bar};
      }
    else
      {
        return 1.0;
      }
  }


  template <int dim>
  typename PerfectlyMatchedLayer<dim>::rank2_type
  PerfectlyMatchedLayer<dim>::rotation(std::complex<double> d_1,
                                       std::complex<double> d_2,
                                       Point<dim>           point)
  {
    rank2_type result;
    result[0][0] = point[0] * point[0] * d_1 + point[1] * point[1] * d_2;
    result[0][1] = point[0] * point[1] * (d_1 - d_2);
    result[1][0] = point[0] * point[1] * (d_1 - d_2);
    result[1][1] = point[1] * point[1] * d_1 + point[0] * point[0] * d_2;
    return result;
  }


  template <int dim>
  typename PerfectlyMatchedLayer<dim>::rank2_type
  PerfectlyMatchedLayer<dim>::a_matrix(const Point<dim> point)
  {
    const auto d     = this->d(point);
    const auto d_bar = this->d_bar(point);
    return invert(rotation(d * d, d * d_bar, point)) *
           rotation(d * d, d * d_bar, point);
  }


  template <int dim>
  typename PerfectlyMatchedLayer<dim>::rank2_type
  PerfectlyMatchedLayer<dim>::b_matrix(const Point<dim> point)
  {
    const auto d     = this->d(point);
    const auto d_bar = this->d_bar(point);
    return invert(rotation(d, d_bar, point)) * rotation(d, d_bar, point);
  }


  template <int dim>
  typename PerfectlyMatchedLayer<dim>::rank2_type
  PerfectlyMatchedLayer<dim>::c_matrix(const Point<dim> point)
  {
    const auto d     = this->d(point);
    const auto d_bar = this->d_bar(point);
    return invert(rotation(1. / d_bar, 1. / d, point)) *
           rotation(1. / d_bar, 1. / d, point);
  }


  // @sect4{Maxwell Class}
  // At this point we are ready to declare all the major building blocks of
  // the finite element program which consists of the usual setup and
  // assembly routines. Most of the structure has already been introduced
  // in previous tutorial programs. The Maxwell class also holds private
  // instances of the Parameters and PerfectlyMatchedLayers classes
  // introduced above. The default values of these parameters are set to
  // show us a standing wave with absorbing boundary conditions and a PML.

  template <int dim>
  class Maxwell : public ParameterAcceptor
  {
  public:
    Maxwell();
    void run();

  private:
    /* run time parameters */
    double       scaling;
    unsigned int refinements;
    unsigned int fe_order;
    unsigned int quadrature_order;
    bool         absorbing_boundary;

    void parse_parameters_callback();
    void make_grid();
    void setup_system();
    void assemble_system();
    void solve();
    void output_results();

    Parameters<dim>            parameters;
    PerfectlyMatchedLayer<dim> perfectly_matched_layer;

    Triangulation<dim> triangulation;
    DoFHandler<dim>    dof_handler;

    std::unique_ptr<FiniteElement<dim>> fe;

    AffineConstraints<double> constraints;
    SparsityPattern           sparsity_pattern;
    SparseMatrix<double>      system_matrix;
    Vector<double>            solution;
    Vector<double>            system_rhs;
  };

  // @sect3{Class Template Definitions and Implementation}
  //
  // @sect4{The Constructor}
  // The Constructor simply consists of default initialization a number of
  // discretization parameters (such as the domain size, mesh refinement,
  // and the order of finite elements and quadrature) and declaring a
  // corresponding entry via ParameterAcceptor::add_parameter(). All of
  // these can be modified by editing the .prm file. Absorbing boundary
  // conditions can be controlled with the absorbing_boundary boolean. If
  // absorbing boundary conditions are disabled we simply enforce
  // homogeneous Dirichlet conditions on the tangential component of the
  // electric field. In the context of time-harmonic Maxwell's equations
  // these are also known as perfectly conducting boundary conditions.

  template <int dim>
  Maxwell<dim>::Maxwell()
    : ParameterAcceptor("Maxwell")
    , dof_handler(triangulation)
  {
    ParameterAcceptor::parse_parameters_call_back.connect(
      [&]() { parse_parameters_callback(); });

    scaling = 20;
    add_parameter("scaling", scaling, "scale of the hypercube geometry");

    refinements = 8;
    add_parameter("refinements",
                  refinements,
                  "number of refinements of the geometry");

    fe_order = 0;
    add_parameter("fe order", fe_order, "order of the finite element space");

    quadrature_order = 1;
    add_parameter("quadrature order",
                  quadrature_order,
                  "order of the quadrature");

    absorbing_boundary = true;
    add_parameter("absorbing boundary condition",
                  absorbing_boundary,
                  "use absorbing boundary conditions?");
  }


  template <int dim>
  void Maxwell<dim>::parse_parameters_callback()
  {
    fe = std::make_unique<FESystem<dim>>(FE_NedelecSZ<dim>(fe_order), 2);
  }

  // The Maxwell::make_grid() routine creates the mesh for the
  // computational domain which in our case is a scaled square domain.
  // Additionally, a material interface is introduced by setting the
  // material id of the upper half ($y>0$) to 1 and of the lower half
  // ($y<0$) of the computational domain to 2.
  // We are using a block decomposition into real and imaginary matrices
  // for the solution matrices. More details on this are available
  // under the Results section.

  template <int dim>
  void Maxwell<dim>::make_grid()
  {
    GridGenerator::hyper_cube(triangulation, -scaling, scaling);
    triangulation.refine_global(refinements);

    if (absorbing_boundary)
      {
        for (auto &face : triangulation.active_face_iterators())
          if (face->at_boundary())
            face->set_boundary_id(1);
      };

    for (auto &cell : triangulation.active_cell_iterators())
      if (cell->center()[1] > 0.)
        cell->set_material_id(1);
      else
        cell->set_material_id(2);


    std::cout << "Number of active cells: " << triangulation.n_active_cells()
              << std::endl;
  }

  // The Maxwell::setup_system() routine follows the usual routine of
  // enumerating all the degrees of freedom and setting up the matrix and
  // vector objects to hold the system data. Enumerating is done by using
  // DoFHandler::distribute_dofs().

  template <int dim>
  void Maxwell<dim>::setup_system()
  {
    dof_handler.distribute_dofs(*fe);
    std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    constraints.clear();

    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler,
      0, /* real part */
      Functions::ZeroFunction<dim>(2 * dim),
      0, /* boundary id */
      constraints);
    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler,
      dim, /* imaginary part */
      Functions::ZeroFunction<dim>(2 * dim),
      0, /* boundary id */
      constraints);

    constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /* keep_constrained_dofs = */ true);
    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
  }

  // This is a helper function that takes the tangential component of a tensor.
  template <int dim>
  DEAL_II_ALWAYS_INLINE inline Tensor<1, dim, std::complex<double>>
  tangential_part(const Tensor<1, dim, std::complex<double>> &tensor,
                  const Tensor<1, dim>                       &normal)
  {
    auto result = tensor;
    result[0]   = normal[1] * (tensor[0] * normal[1] - tensor[1] * normal[0]);
    result[1]   = -normal[0] * (tensor[0] * normal[1] - tensor[1] * normal[0]);
    return result;
  }


  // Assemble the stiffness matrix and the right-hand side:
  //\f{align*}{
  // A_{ij} = \int_\Omega (\mu_r^{-1}\nabla \times \varphi_j) \cdot
  // (\nabla\times\bar{\varphi}_i)\text{d}x
  // - \int_\Omega \varepsilon_r\varphi_j \cdot \bar{\varphi}_i\text{d}x
  // - i\int_\Sigma (\sigma_r^{\Sigma}(\varphi_j)_T) \cdot
  // (\bar{\varphi}_i)_T\text{do}x
  // - i\int_{\partial\Omega} (\sqrt{\mu_r^{-1}\varepsilon}(\varphi_j)_T) \cdot
  // (\nabla\times(\bar{\varphi}_i)_T)\text{d}x, \f} \f{align}{
  //  F_i = i\int_\Omega J_a \cdot \bar{\varphi_i}\text{d}x - \int_\Omega
  //  \mu_r^{-1} \cdot (\nabla \times \bar{\varphi_i}) \text{d}x.
  // \f}
  // In addition, we will be modifying the coefficients if the position of the
  // cell is within the PML region.

  template <int dim>
  void Maxwell<dim>::assemble_system()
  {
    const QGauss<dim>     quadrature_formula(quadrature_order);
    const QGauss<dim - 1> face_quadrature_formula(quadrature_order);

    FEValues<dim, dim>     fe_values(*fe,
                                 quadrature_formula,
                                 update_values | update_gradients |
                                   update_quadrature_points |
                                   update_JxW_values);
    FEFaceValues<dim, dim> fe_face_values(*fe,
                                          face_quadrature_formula,
                                          update_values | update_gradients |
                                            update_quadrature_points |
                                            update_normal_vectors |
                                            update_JxW_values);

    const unsigned int dofs_per_cell = fe->dofs_per_cell;

    const unsigned int n_q_points      = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // Next, let us assemble on the interior of the domain on the left hand
    // side. So we are computing
    // \f{align*}{
    //   \int_\Omega (\mu_r^{-1}\nabla \times \varphi_i) \cdot
    //   (\nabla\times\bar{\varphi}_j)\text{d}x
    //   -
    //   \int_\Omega \varepsilon_r\varphi_i \cdot \bar{\varphi}_j\text{d}x
    // \f}
    // and
    // \f{align}{
    //   i\int_\Omega J_a \cdot \bar{\varphi_i}\text{d}x
    //   - \int_\Omega \mu_r^{-1} \cdot (\nabla \times \bar{\varphi_i})
    //   \text{d}x.
    // \f}
    // In doing so, we need test functions $\varphi_i$ and $\varphi_j$, and the
    // curl of these test variables. We must be careful with the signs of the
    // imaginary parts of these complex test variables. Moreover, we have a
    // conditional that changes the parameters if the cell is in the PML region.
    const FEValuesExtractors::Vector real_part(0);
    const FEValuesExtractors::Vector imag_part(dim);
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);

        cell_matrix = 0.;
        cell_rhs    = 0.;

        cell->get_dof_indices(local_dof_indices);
        const auto id = cell->material_id();

        const auto &quadrature_points = fe_values.get_quadrature_points();

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {
            const Point<dim> &position = quadrature_points[q_point];

            auto       mu_inv  = parameters.mu_inv(position, id);
            auto       epsilon = parameters.epsilon(position, id);
            const auto J_a     = parameters.J_a(position, id);

            const auto A = perfectly_matched_layer.a_matrix(position);
            const auto B = perfectly_matched_layer.b_matrix(position);
            const auto d = perfectly_matched_layer.d(position);

            mu_inv  = mu_inv / d;
            epsilon = invert(A) * epsilon * invert(B);

            for (const auto i : fe_values.dof_indices())
              {
                constexpr std::complex<double> imag{0., 1.};

                const auto phi_i =
                  fe_values[real_part].value(i, q_point) -
                  imag * fe_values[imag_part].value(i, q_point);
                const auto curl_phi_i =
                  fe_values[real_part].curl(i, q_point) -
                  imag * fe_values[imag_part].curl(i, q_point);

                const auto rhs_value =
                  (imag * scalar_product(J_a, phi_i)) * fe_values.JxW(q_point);
                cell_rhs(i) += rhs_value.real();

                for (const auto j : fe_values.dof_indices())
                  {
                    const auto phi_j =
                      fe_values[real_part].value(j, q_point) +
                      imag * fe_values[imag_part].value(j, q_point);
                    const auto curl_phi_j =
                      fe_values[real_part].curl(j, q_point) +
                      imag * fe_values[imag_part].curl(j, q_point);

                    const auto temp =
                      (scalar_product(mu_inv * curl_phi_j, curl_phi_i) -
                       scalar_product(epsilon * phi_j, phi_i)) *
                      fe_values.JxW(q_point);
                    cell_matrix(i, j) += temp.real();
                  }
              }
          }

        // Now we assemble the face and the boundary. The following loops will
        // assemble
        // \f{align*}{
        // - i\int_\Sigma (\sigma_r^{\Sigma}(\varphi_i)_T) \cdot
        // (\bar{\varphi}_j)_T\text{do}x \f} and \f{align}{
        //  - i\int_{\partial\Omega} (\sqrt{\mu_r^{-1}\varepsilon}(\varphi_i)_T)
        //  \cdot (\nabla\times(\bar{\varphi}_j)_T)\text{d}x,
        // \f}
        // respectively. The test variables and the PML are implemented
        // similarly as the domain.
        //
        // If we are at the domain boundary $\partial\Omega$ and absorbing
        // boundary conditions are set (<code>id == 1</code>) we assemble
        // the corresponding boundary term:
        //
        const FEValuesExtractors::Vector real_part(0);
        const FEValuesExtractors::Vector imag_part(dim);
        for (const auto &face : cell->face_iterators())
          {
            if (face->at_boundary())
              {
                const auto id = face->boundary_id();
                if (id == 1)
                  {
                    fe_face_values.reinit(cell, face);

                    for (unsigned int q_point = 0; q_point < n_face_q_points;
                         ++q_point)
                      {
                        const auto &position = quadrature_points[q_point];

                        auto mu_inv  = parameters.mu_inv(position, id);
                        auto epsilon = parameters.epsilon(position, id);

                        const auto A =
                          perfectly_matched_layer.a_matrix(position);
                        const auto B =
                          perfectly_matched_layer.b_matrix(position);
                        const auto d = perfectly_matched_layer.d(position);

                        mu_inv  = mu_inv / d;
                        epsilon = invert(A) * epsilon * invert(B);

                        const auto normal =
                          fe_face_values.normal_vector(q_point);

                        for (const auto i : fe_face_values.dof_indices())
                          {
                            constexpr std::complex<double> imag{0., 1.};

                            const auto phi_i =
                              fe_face_values[real_part].value(i, q_point) -
                              imag *
                                fe_face_values[imag_part].value(i, q_point);
                            const auto phi_i_T = tangential_part(phi_i, normal);

                            for (const auto j : fe_face_values.dof_indices())
                              {
                                const auto phi_j =
                                  fe_face_values[real_part].value(j, q_point) +
                                  imag *
                                    fe_face_values[imag_part].value(j, q_point);
                                const auto phi_j_T =
                                  tangential_part(phi_j, normal) *
                                  fe_face_values.JxW(q_point);

                                const auto prod      = mu_inv * epsilon;
                                const auto sqrt_prod = prod;

                                const auto temp =
                                  -imag * scalar_product((sqrt_prod * phi_j_T),
                                                         phi_i_T);
                                cell_matrix(i, j) += temp.real();
                              } /* j */
                          }     /* i */
                      }         /* q_point */
                  }
              }
            else
              {
                // We are on an interior face:
                const auto face_index = cell->face_iterator_to_index(face);

                const auto id1 = cell->material_id();
                const auto id2 = cell->neighbor(face_index)->material_id();

                if (id1 == id2)
                  continue; /* skip this face */

                fe_face_values.reinit(cell, face);

                for (unsigned int q_point = 0; q_point < n_face_q_points;
                     ++q_point)
                  {
                    const auto &position = quadrature_points[q_point];

                    auto sigma = parameters.sigma(position, id1, id2);

                    const auto B = perfectly_matched_layer.b_matrix(position);
                    const auto C = perfectly_matched_layer.c_matrix(position);
                    sigma        = invert(C) * sigma * invert(B);

                    const auto normal = fe_face_values.normal_vector(q_point);

                    for (const auto i : fe_face_values.dof_indices())
                      {
                        constexpr std::complex<double> imag{0., 1.};

                        const auto phi_i =
                          fe_face_values[real_part].value(i, q_point) -
                          imag * fe_face_values[imag_part].value(i, q_point);
                        const auto phi_i_T = tangential_part(phi_i, normal);

                        for (const auto j : fe_face_values.dof_indices())
                          {
                            const auto phi_j =
                              fe_face_values[real_part].value(j, q_point) +
                              imag *
                                fe_face_values[imag_part].value(j, q_point);
                            const auto phi_j_T = tangential_part(phi_j, normal);

                            const auto temp =
                              -imag *
                              scalar_product((sigma * phi_j_T), phi_i_T) *
                              fe_face_values.JxW(q_point);
                            cell_matrix(i, j) += temp.real();
                          } /* j */
                      }     /* i */
                  }         /* q_point */
              }
          }

        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
  }

  // We use a direct solver from the SparseDirectUMFPACK to solve the system
  template <int dim>
  void Maxwell<dim>::solve()
  {
    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);
    A_direct.vmult(solution, system_rhs);
  }

  // The output is written into a vtk file with 4 components
  template <int dim>
  void Maxwell<dim>::output_results()
  {
    DataOut<2> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution,
                             {"real_Ex", "real_Ey", "imag_Ex", "imag_Ey"});
    data_out.build_patches();
    const std::string filename = "solution.vtk";
    std::ofstream     output(filename);
    data_out.write_vtk(output);
    std::cout << "Output written to " << filename << std::endl;
  }


  template <int dim>
  void Maxwell<dim>::run()
  {
    make_grid();
    setup_system();
    assemble_system();
    solve();
    output_results();
  }

} // namespace Step81

// The following main function calls the class step-81(), initializes the
// ParameterAcceptor, and calls the run() function.

int main()
{
  try
    {
      using namespace dealii;

      Step81::Maxwell<2> maxwell_2d;
      ParameterAcceptor::initialize("parameters.prm");
      maxwell_2d.run();
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
